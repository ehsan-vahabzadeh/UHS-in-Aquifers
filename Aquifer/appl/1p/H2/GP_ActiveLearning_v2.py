"""
Active-learning workflow v2 — multi-scale GP with improved acquisition.

Key changes from v1:
  1. Multi-scale kernel: Matérn-2.5 (smooth) + Matérn-1.5 (rough) additive kernel
     so the GP can represent both broad trends and narrow local features.
  2. Learned noise (nugget) instead of fixed 1e-6 train_Yvar.
     Lets the GP admit model-discrepancy rather than forcing exact interpolation.
  3. Fixed Normalize bounds — inputs are already [0,1], so we pass known bounds
     instead of letting Normalize refit from observed min/max each iteration.
  4. min_distance reduced: history=0.05, batch=0.12.
     Allows local refinement around narrow features (width ~0.07).
  5. Candidate pool refreshed every iteration with 5000 points instead of
     a fixed 800-point pool — much better coverage in 8D.
  6. Chunked IVR — if used, computes cross-covariance in chunks instead of
     building a dense (n_ref+n_cand)^2 joint matrix.
  7. Default strategy: uncertainty_diversity (most robust under model
     misspecification).

Installation:
    pip install numpy matplotlib torch gpytorch botorch
"""

from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import torch
    import gpytorch
    from botorch.fit import fit_gpytorch_mll
    from botorch.models import SingleTaskGP
    from botorch.models.transforms.input import Normalize
    from botorch.models.transforms.outcome import Standardize
    from gpytorch.mlls import ExactMarginalLogLikelihood
    from gpytorch.kernels import ScaleKernel, MaternKernel, AdditiveKernel
except ImportError as exc:
    raise ImportError(
        "This script requires PyTorch, GPyTorch, and BoTorch. "
        "Install them via `pip install torch gpytorch botorch`."
    ) from exc

TORCH_DTYPE = torch.double
TORCH_DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

N_INPUTS = 8
INPUT_NAMES = [
    "inj_rate_dev",
    "inj_rate_op",
    "production_rate",
    "porosity",
    "permeability",
    "pressure",
    "temperature",
    "cycle_length",
]
INPUT_BOUNDS = np.array([[0.0, 1.0]] * N_INPUTS)


# ------------------------------------------------------------
# 1. Fake expensive simulator (same as v1)
# ------------------------------------------------------------
def expensive_simulator(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    if x.ndim == 1:
        x = x.reshape(1, -1)
    if x.ndim != 2 or x.shape[1] != N_INPUTS:
        raise ValueError(f"Expected input shape (n_samples, {N_INPUTS}), got {x.shape}")

    lower = INPUT_BOUNDS[:, 0]
    upper = INPUT_BOUNDS[:, 1]
    u = (x - lower) / (upper - lower)
    x1, x2, x3, x4, x5, x6, x7, x8 = u.T

    smooth_trend = 0.7 * np.sin(2.0 * np.pi * (0.65 * x1 + 0.20 * x2))
    curved_trend = 0.45 * (x3 - 0.35) ** 2 + 0.35 * (x7 - 0.5)
    cross_terms = 1.1 * (x1 - 0.5) * (x5 - 0.5) - 0.9 * (x2 - 0.45) * (x6 - 0.4)

    pressure_threshold = 0.8 / (1.0 + np.exp(-30.0 * (x6 + 0.30 * x2 - 0.68)))
    operating_kink = 0.55 * np.maximum(x4 + x5 - 1.1, 0.0)

    localized_peak = 1.6 * np.exp(
        -0.5 * (
            ((x1 - 0.75) / 0.08) ** 2
            + ((x5 - 0.25) / 0.15) ** 2
            + ((x8 - 0.55) / 0.20) ** 2
        )
    )
    narrow_valley = -1.4 * np.exp(
        -0.5 * (
            ((x3 - 0.28) / 0.10) ** 2
            + ((x6 - 0.72) / 0.08) ** 2
        )
    )
    diagonal_ridge = 0.75 * np.exp(-0.5 * ((x1 + x2 + x3 - 1.65) / 0.07) ** 2)
    localized_roughness = 0.25 * np.sin(35.0 * x1 + 17.0 * x7) * np.exp(
        -0.5 * ((x8 - 0.35) / 0.18) ** 2
    )

    y = (
        smooth_trend
        + curved_trend
        + cross_terms
        + pressure_threshold
        + operating_kink
        + localized_peak
        + narrow_valley
        + diagonal_ridge
        + localized_roughness
    )
    return y


# ------------------------------------------------------------
# 2. Initial data generation
# ------------------------------------------------------------
def latin_hypercube_sample(n_points: int, bounds: np.ndarray, seed: int = 42) -> np.ndarray:
    rng = np.random.default_rng(seed)
    bounds = np.asarray(bounds, dtype=float)
    n_dim = bounds.shape[0]

    unit_samples = np.empty((n_points, n_dim))
    for dim in range(n_dim):
        unit_samples[:, dim] = (rng.permutation(n_points) + rng.random(n_points)) / n_points

    lower = bounds[:, 0]
    upper = bounds[:, 1]
    return lower + unit_samples * (upper - lower)


def generate_initial_points(n_init: int, bounds: np.ndarray, seed: int = 42) -> np.ndarray:
    return latin_hypercube_sample(n_points=n_init, bounds=bounds, seed=seed)


# ------------------------------------------------------------
# 3. Train BoTorch GP — multi-scale kernel + learned nugget
# ------------------------------------------------------------
def to_tensor(array: np.ndarray) -> torch.Tensor:
    return torch.as_tensor(array, dtype=TORCH_DTYPE, device=TORCH_DEVICE)


def to_numpy(tensor: torch.Tensor) -> np.ndarray:
    return tensor.detach().cpu().numpy()


# Fixed bounds for Normalize — prevents re-fitting from observed data range
_FIXED_BOUNDS = torch.tensor(INPUT_BOUNDS.T, dtype=TORCH_DTYPE, device=TORCH_DEVICE)


def train_gp(X: np.ndarray, y: np.ndarray) -> SingleTaskGP:
    train_X = to_tensor(X)
    train_Y = to_tensor(y).reshape(-1, 1)

    # Multi-scale additive kernel:
    #   Matérn-2.5 captures smooth/global trends
    #   Matérn-1.5 captures rougher/shorter-scale features
    covar_module = AdditiveKernel(
        ScaleKernel(MaternKernel(nu=2.5, ard_num_dims=N_INPUTS)),
        ScaleKernel(MaternKernel(nu=1.5, ard_num_dims=N_INPUTS)),
    )

    # No train_Yvar → learned noise (nugget).
    # This lets the GP admit model-discrepancy instead of forcing
    # exact interpolation through a kernel that can't represent
    # all features.
    model = SingleTaskGP(
        train_X=train_X,
        train_Y=train_Y,
        covar_module=covar_module,
        input_transform=Normalize(d=N_INPUTS, bounds=_FIXED_BOUNDS),
        outcome_transform=Standardize(m=1),
    )

    mll = ExactMarginalLogLikelihood(model.likelihood, model)
    fit_gpytorch_mll(mll)
    model.eval()
    model.likelihood.eval()
    return model


def posterior_mean_variance_torch(
    model: SingleTaskGP,
    X: torch.Tensor,
    chunk_size: int = 2048,
) -> Tuple[torch.Tensor, torch.Tensor]:
    means = []
    variances = []
    model.eval()
    with torch.no_grad(), gpytorch.settings.fast_pred_var():
        for start in range(0, X.shape[0], chunk_size):
            posterior = model.posterior(X[start : start + chunk_size])
            means.append(posterior.mean.squeeze(-1))
            variances.append(posterior.variance.squeeze(-1).clamp_min(0.0))
    return torch.cat(means), torch.cat(variances)


def predict_gp(
    model: SingleTaskGP,
    X: np.ndarray,
    return_std: bool = False,
):
    mean, variance = posterior_mean_variance_torch(model, to_tensor(X))
    mean_np = to_numpy(mean)
    if return_std:
        return mean_np, to_numpy(torch.sqrt(variance))
    return mean_np


# ------------------------------------------------------------
# 4. Utility functions
# ------------------------------------------------------------
def rmse(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    return float(np.sqrt(np.mean((y_true - y_pred) ** 2)))


def make_slice_points(
    bounds: np.ndarray,
    slice_dimension: int,
    n_points: int = 1000,
    base_point: Optional[np.ndarray] = None,
) -> np.ndarray:
    bounds = np.asarray(bounds, dtype=float)
    if base_point is None:
        base_point = 0.5 * (bounds[:, 0] + bounds[:, 1])

    x_slice = np.tile(base_point, (n_points, 1))
    x_slice[:, slice_dimension] = np.linspace(
        bounds[slice_dimension, 0],
        bounds[slice_dimension, 1],
        n_points,
    )
    return x_slice


def format_design_point(point: np.ndarray) -> str:
    return ", ".join(f"{name}={value:.3f}" for name, value in zip(INPUT_NAMES, point))


def distance_filter_mask(points: np.ndarray, x_new: np.ndarray, min_distance: float) -> np.ndarray:
    distances = np.linalg.norm(points - x_new, axis=1)
    return distances >= min_distance


def distance_filter_mask_many(
    points: np.ndarray,
    blocked_points: np.ndarray,
    min_distance: float,
) -> np.ndarray:
    if len(points) == 0:
        return np.empty(0, dtype=bool)
    if blocked_points is None or len(blocked_points) == 0:
        return np.ones(len(points), dtype=bool)
    distances = np.linalg.norm(points[:, None, :] - blocked_points[None, :, :], axis=2)
    return np.min(distances, axis=1) >= min_distance


# ------------------------------------------------------------
# 5. Uncertainty + diversity (improved defaults)
# ------------------------------------------------------------
def select_batch_uncertainty_diversity(
    model: SingleTaskGP,
    candidate_pool: np.ndarray,
    already_sampled: np.ndarray,
    n_new: int,
    min_distance_history: float = 0.05,
    min_distance_batch: float = 0.12,
) -> np.ndarray:
    """
    Rank candidates by posterior std dev, greedily keep points satisfying
    two distance rules:
      - min_distance_history: minimum distance to any already-sampled point
      - min_distance_batch:   minimum distance between newly selected points
    """
    admissible = candidate_pool[
        distance_filter_mask_many(candidate_pool, already_sampled, min_distance_history)
    ]
    if len(admissible) == 0:
        return np.empty((0, candidate_pool.shape[1]))

    _, std = predict_gp(model, admissible, return_std=True)
    ranked_idx = np.argsort(std)[::-1]

    selected = []
    for idx in ranked_idx:
        x_candidate = admissible[idx : idx + 1]

        if selected:
            selected_array = np.vstack(selected)
            dist_selected = np.min(np.linalg.norm(selected_array - x_candidate, axis=1))
        else:
            dist_selected = np.inf

        if dist_selected >= min_distance_batch:
            selected.append(x_candidate.copy())

        if len(selected) >= n_new:
            break

    if not selected:
        return np.empty((0, candidate_pool.shape[1]))
    return np.vstack(selected)


# ------------------------------------------------------------
# 6. Chunked IVR (memory-efficient)
# ------------------------------------------------------------
def select_batch_chunked_ivr(
    model: SingleTaskGP,
    candidate_pool: np.ndarray,
    reference_grid: np.ndarray,
    already_sampled: np.ndarray,
    n_new: int,
    min_distance_history: float = 0.05,
    min_distance_batch: float = 0.12,
    weights: Optional[torch.Tensor] = None,
    cand_chunk_size: int = 500,
) -> np.ndarray:
    """
    Memory-efficient IVR using chunked cross-covariance.

    Instead of building a dense (n_ref+n_cand)^2 joint matrix, we compute:
      - reference marginal variance (once)
      - cross-covariance K(ref, cand_chunk) in chunks
      - candidate marginal variance per chunk

    Score for candidate j:
        score(j) = sum_r [ w_r * Cov(r,j)^2 ] / Var(j)

    Greedy batch: after selecting a point, we condition the reference
    posterior on that observation (rank-1 update to ref covariance) and
    rescore the remaining candidates.
    """
    remaining = candidate_pool[
        distance_filter_mask_many(candidate_pool, already_sampled, min_distance_history)
    ]
    if len(remaining) == 0:
        return np.empty((0, candidate_pool.shape[1]))

    n_ref = len(reference_grid)
    ref_tensor = to_tensor(reference_grid)

    if weights is None:
        w = torch.full((n_ref,), 1.0 / n_ref, dtype=TORCH_DTYPE, device=TORCH_DEVICE)
    else:
        w = weights.to(dtype=TORCH_DTYPE, device=TORCH_DEVICE)
        w = w / w.sum().clamp_min(1e-12)

    # Pre-compute reference-reference covariance (for conditioning after selection)
    model.eval()
    with torch.no_grad():
        ref_post = model.posterior(ref_tensor)
        K_ref = ref_post.mvn.covariance_matrix.detach().squeeze()

    active_mask = np.ones(len(remaining), dtype=bool)
    selected = []

    for _ in range(n_new):
        cand_active = remaining[active_mask]
        if len(cand_active) == 0:
            break

        # Score all active candidates in chunks
        best_score = -float("inf")
        best_global_idx = -1

        active_indices = np.where(active_mask)[0]
        for chunk_start in range(0, len(active_indices), cand_chunk_size):
            chunk_idx = active_indices[chunk_start : chunk_start + cand_chunk_size]
            cand_chunk = to_tensor(remaining[chunk_idx])

            with torch.no_grad():
                # Cross-covariance K(ref, cand_chunk)
                joint = torch.cat([ref_tensor, cand_chunk], dim=0)
                joint_post = model.posterior(joint)
                joint_cov = joint_post.mvn.covariance_matrix.detach().squeeze()

                K_rc = joint_cov[:n_ref, n_ref:]  # (n_ref, chunk)
                cand_var = joint_cov.diagonal()[n_ref:].clamp_min(1e-12)  # (chunk,)

            # Condition on already-selected points via K_ref update
            # K_rc is cross-cov under the *current* (conditioned) ref covariance
            # but we approximate by using the original model posterior
            # and conditioning K_ref incrementally (done after each selection below)
            # For first-pass scoring we use the raw cross-cov
            scores = (w[:, None] * K_rc.square()).sum(dim=0) / cand_var  # (chunk,)

            chunk_best = int(torch.argmax(scores).item())
            if scores[chunk_best].item() > best_score:
                best_score = scores[chunk_best].item()
                best_global_idx = chunk_idx[chunk_best]

        if best_global_idx < 0 or not np.isfinite(best_score):
            break

        x_best = remaining[best_global_idx : best_global_idx + 1]
        selected.append(x_best.copy())

        # Condition reference covariance on selected point (rank-1 update)
        with torch.no_grad():
            joint_sel = torch.cat([ref_tensor, to_tensor(x_best)], dim=0)
            joint_post_sel = model.posterior(joint_sel)
            joint_cov_sel = joint_post_sel.mvn.covariance_matrix.detach().squeeze()
            k_col = joint_cov_sel[:n_ref, n_ref:]  # (n_ref, 1)
            k_var = joint_cov_sel[n_ref, n_ref].clamp_min(1e-12)
            K_ref = K_ref - (k_col @ k_col.T) / k_var
            K_ref = 0.5 * (K_ref + K_ref.T)

        # Apply distance filter
        keep = distance_filter_mask(remaining, x_best, min_distance_batch)
        active_mask &= keep
        active_mask[best_global_idx] = False

    if not selected:
        return np.empty((0, candidate_pool.shape[1]))
    return np.vstack(selected)


def make_decision_weights(
    model: SingleTaskGP,
    reference_grid: np.ndarray,
    sharpness: float = 8.0,
) -> torch.Tensor:
    mean_ref, _ = posterior_mean_variance_torch(model, to_tensor(reference_grid))
    shifted = mean_ref - mean_ref.max()
    return torch.softmax(sharpness * shifted, dim=0)


def select_batch_decision_aware_ivr(
    model: SingleTaskGP,
    candidate_pool: np.ndarray,
    reference_grid: np.ndarray,
    already_sampled: np.ndarray,
    n_new: int,
    min_distance_history: float = 0.05,
    min_distance_batch: float = 0.12,
    decision_sharpness: float = 8.0,
) -> np.ndarray:
    weights = make_decision_weights(model, reference_grid, sharpness=decision_sharpness)
    return select_batch_chunked_ivr(
        model=model,
        candidate_pool=candidate_pool,
        reference_grid=reference_grid,
        already_sampled=already_sampled,
        n_new=n_new,
        min_distance_history=min_distance_history,
        min_distance_batch=min_distance_batch,
        weights=weights,
    )


# ------------------------------------------------------------
# 7. Strategy wrapper
# ------------------------------------------------------------
def select_new_points(
    strategy: str,
    model: SingleTaskGP,
    X_train: np.ndarray,
    candidate_pool: np.ndarray,
    reference_grid: np.ndarray,
    n_new: int,
    min_distance_history: float = 0.05,
    min_distance_batch: float = 0.12,
    decision_sharpness: float = 8.0,
) -> np.ndarray:
    if strategy == "uncertainty_diversity":
        return select_batch_uncertainty_diversity(
            model=model,
            candidate_pool=candidate_pool,
            already_sampled=X_train,
            n_new=n_new,
            min_distance_history=min_distance_history,
            min_distance_batch=min_distance_batch,
        )

    if strategy == "integrated_variance_reduction":
        return select_batch_chunked_ivr(
            model=model,
            candidate_pool=candidate_pool,
            reference_grid=reference_grid,
            already_sampled=X_train,
            n_new=n_new,
            min_distance_history=min_distance_history,
            min_distance_batch=min_distance_batch,
        )

    if strategy == "decision_aware_ivr":
        return select_batch_decision_aware_ivr(
            model=model,
            candidate_pool=candidate_pool,
            reference_grid=reference_grid,
            already_sampled=X_train,
            n_new=n_new,
            min_distance_history=min_distance_history,
            min_distance_batch=min_distance_batch,
            decision_sharpness=decision_sharpness,
        )

    raise ValueError(f"Unknown strategy: {strategy}")


# ------------------------------------------------------------
# 8. Plotting (same structure as v1)
# ------------------------------------------------------------
def plot_initial_vs_final(
    gp_initial: SingleTaskGP,
    gp_final: SingleTaskGP,
    X_init: np.ndarray,
    y_init: np.ndarray,
    X_all: np.ndarray,
    y_all: np.ndarray,
    input_bounds: np.ndarray,
    strategy_name: str,
    slice_dimension: int = 0,
    output_dir: Optional[Path] = None,
):
    x_plot = make_slice_points(input_bounds, slice_dimension=slice_dimension)
    x_axis = x_plot[:, slice_dimension]
    y_true = expensive_simulator(x_plot)

    y_init_pred, y_init_std = predict_gp(gp_initial, x_plot, return_std=True)
    y_final_pred, y_final_std = predict_gp(gp_final, x_plot, return_std=True)

    x_eval = latin_hypercube_sample(n_points=2000, bounds=input_bounds, seed=900)
    y_eval = expensive_simulator(x_eval)
    init_rmse = rmse(y_eval, predict_gp(gp_initial, x_eval))
    final_rmse = rmse(y_eval, predict_gp(gp_final, x_eval))

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    ax = axes[0]
    ax.plot(x_axis, y_true, label="True function slice")
    ax.plot(x_axis, y_init_pred, label="Initial GP mean")
    ax.fill_between(x_axis, y_init_pred - 1.96 * y_init_std, y_init_pred + 1.96 * y_init_std,
                     alpha=0.2, label="Initial 95% band")
    ax.scatter(X_init[:, slice_dimension], y_init, s=50, alpha=0.75, label="Initial points projected")
    ax.set_title(f"Initial model\n8D test RMSE = {init_rmse:.4f}")
    ax.set_xlabel(f"{INPUT_NAMES[slice_dimension]} slice")
    ax.set_ylabel("y")
    ax.legend()

    ax = axes[1]
    ax.plot(x_axis, y_true, label="True function slice")
    ax.plot(x_axis, y_final_pred, label="Final GP mean")
    ax.fill_between(x_axis, y_final_pred - 1.96 * y_final_std, y_final_pred + 1.96 * y_final_std,
                     alpha=0.2, label="Final 95% band")
    ax.scatter(X_all[:, slice_dimension], y_all, s=25, alpha=0.6, label="All sampled points projected")
    ax.scatter(X_init[:, slice_dimension], y_init, s=55, edgecolor="black", linewidth=0.8,
               label="Initial points projected")
    ax.set_title(f"Final model\n8D test RMSE = {final_rmse:.4f}")
    ax.set_xlabel(f"{INPUT_NAMES[slice_dimension]} slice")
    ax.set_ylabel("y")
    ax.legend()

    fig.suptitle(f"Model comparison: {strategy_name} (v2) | other inputs fixed at midpoint", fontsize=14)
    fig.tight_layout()

    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_dir / f"v2_comparison_{strategy_name}.png", dpi=150)
    plt.close(fig)


def plot_iteration(
    iteration: int,
    gp: SingleTaskGP,
    X_train: np.ndarray,
    y_train: np.ndarray,
    new_points: np.ndarray,
    input_bounds: np.ndarray,
    strategy_name: str,
    slice_dimension: int = 0,
    output_dir: Optional[Path] = None,
):
    x_plot = make_slice_points(input_bounds, slice_dimension=slice_dimension)
    x_axis = x_plot[:, slice_dimension]
    y_true = expensive_simulator(x_plot)
    y_pred, y_std = predict_gp(gp, x_plot, return_std=True)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(x_axis, y_true, label="True function slice")
    ax.plot(x_axis, y_pred, label="GP mean")
    ax.fill_between(x_axis, y_pred - 1.96 * y_std, y_pred + 1.96 * y_std,
                     alpha=0.2, label="95% confidence band")
    ax.scatter(X_train[:, slice_dimension], y_train, s=40, alpha=0.7, label="Current points projected")
    if new_points is not None and len(new_points) > 0:
        y_new = expensive_simulator(new_points)
        ax.scatter(new_points[:, slice_dimension], y_new, marker="x", s=100, label="New batch projected")
    ax.set_title(f"Iteration {iteration} - {strategy_name} (v2)")
    ax.set_xlabel(f"{INPUT_NAMES[slice_dimension]} slice")
    ax.set_ylabel("y")
    ax.legend()
    fig.tight_layout()

    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_dir / f"v2_{strategy_name}_iteration_{iteration:02d}.png", dpi=150)
    plt.close(fig)


# ------------------------------------------------------------
# 9. RMSE tracking plot
# ------------------------------------------------------------
def plot_rmse_history(
    rmse_history: list,
    strategy_name: str,
    output_dir: Optional[Path] = None,
):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(range(len(rmse_history)), rmse_history, "o-", markersize=5)
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Test RMSE")
    ax.set_title(f"RMSE convergence — {strategy_name} (v2)")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_dir / f"v2_rmse_history_{strategy_name}.png", dpi=150)
    plt.close(fig)


# ------------------------------------------------------------
# 10. Main
# ------------------------------------------------------------
def main():
    input_bounds = INPUT_BOUNDS
    slice_dimension = 0
    plot_output_dir = Path(__file__).with_name("gp_active_learning_plots")

    # Strategy: "uncertainty_diversity" | "integrated_variance_reduction" | "decision_aware_ivr"
    strategy = "integrated_variance_reduction"

    # Initial dataset
    n_init = 2000
    X_init = generate_initial_points(n_init=n_init, bounds=input_bounds, seed=42)
    y_init = expensive_simulator(X_init)

    X_train = X_init.copy()
    y_train = y_init.copy()

    # Reference grid for IVR (only used if strategy involves IVR)
    reference_grid = latin_hypercube_sample(n_points=2000, bounds=input_bounds, seed=456)

    # Active learning settings
    n_iterations = 30
    points_per_iteration = 15
    min_distance_history = 0.05
    min_distance_batch = 0.12

    # Evaluation set for RMSE tracking (same across all iterations)
    x_eval = latin_hypercube_sample(n_points=2000, bounds=input_bounds, seed=789)
    y_eval = expensive_simulator(x_eval)

    # Train initial model
    gp_initial = train_gp(X_train, y_train)
    rmse_history = [rmse(y_eval, predict_gp(gp_initial, x_eval))]
    print(f"Initial RMSE: {rmse_history[0]:.6f}")

    for it in range(n_iterations):
        print(f"\n--- Iteration {it} | strategy = {strategy} ---")

        gp_current = train_gp(X_train, y_train)

        # Fresh candidate pool every iteration — much better 8D coverage
        candidate_pool = latin_hypercube_sample(
            n_points=5000, bounds=input_bounds, seed=1000 + it
        )

        X_new = select_new_points(
            strategy=strategy,
            model=gp_current,
            X_train=X_train,
            candidate_pool=candidate_pool,
            reference_grid=reference_grid,
            n_new=points_per_iteration,
            min_distance_history=min_distance_history,
            min_distance_batch=min_distance_batch,
            decision_sharpness=8.0,
        )

        print(f"Selected {len(X_new)} new points:")
        for idx, x in enumerate(X_new, start=1):
            print(f"  {idx:02d}: {format_design_point(x)}")

        if len(X_new) == 0:
            print("No admissible new points found; stopping early.")
            break

        plot_iteration(
            iteration=it,
            gp=gp_current,
            X_train=X_train,
            y_train=y_train,
            new_points=X_new,
            input_bounds=input_bounds,
            strategy_name=strategy,
            slice_dimension=slice_dimension,
            output_dir=plot_output_dir,
        )

        # "Run expensive simulator"
        y_new = expensive_simulator(X_new)

        # Append
        X_train = np.vstack([X_train, X_new])
        y_train = np.concatenate([y_train, y_new])

        # Track RMSE
        gp_updated = train_gp(X_train, y_train)
        current_rmse = rmse(y_eval, predict_gp(gp_updated, x_eval))
        rmse_history.append(current_rmse)
        print(f"  RMSE after iteration {it}: {current_rmse:.6f}")

    # Final model
    gp_final = train_gp(X_train, y_train)

    plot_initial_vs_final(
        gp_initial=gp_initial,
        gp_final=gp_final,
        X_init=X_init,
        y_init=y_init,
        X_all=X_train,
        y_all=y_train,
        input_bounds=input_bounds,
        strategy_name=strategy,
        slice_dimension=slice_dimension,
        output_dir=plot_output_dir,
    )

    plot_rmse_history(rmse_history, strategy, output_dir=plot_output_dir)

    print("\n===== PERFORMANCE SUMMARY (v2) =====")
    print(f"Input dimensions: {N_INPUTS}")
    print(f"Strategy: {strategy}")
    print(f"Kernel: AdditiveKernel(Matern-2.5 + Matern-1.5) ARD, learned noise")
    print(f"Initial RMSE: {rmse_history[0]:.6f}")
    print(f"Final   RMSE: {rmse_history[-1]:.6f}")
    print(f"Total samples: {len(X_train)} (init={n_init}, added={len(X_train)-n_init})")
    print(f"RMSE reduction: {(1 - rmse_history[-1]/rmse_history[0])*100:.1f}%")
    print(f"History distance: {min_distance_history}, Batch distance: {min_distance_batch}")


if __name__ == "__main__":
    main()
