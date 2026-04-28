"""
Toy active-learning workflow with an 8-input, single-output simulator.

Installation assumptions:
    pip install numpy matplotlib torch gpytorch botorch

This version uses BoTorch / GPyTorch instead of scikit-learn. The important
change for IVR-style selection is that we no longer refit a temporary GP for
every candidate point. The default fallback IVR score uses the fitted GP's
posterior covariance to estimate variance reduction over a finite reference
grid.
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
except ImportError as exc:
    raise ImportError(
        "This script requires PyTorch, GPyTorch, and BoTorch. "
        "Install them in the Python environment used to run this file, e.g. "
        "`pip install torch gpytorch botorch`."
    ) from exc

try:
    from botorch.acquisition.active_learning import qNegIntegratedPosteriorVariance
except Exception:
    qNegIntegratedPosteriorVariance = None


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
# 1. Fake expensive simulator
# ------------------------------------------------------------
def expensive_simulator(x: np.ndarray) -> np.ndarray:
    """
    Toy stand-in for a costly simulator with eight normalized inputs.
    Input shape: (n_samples, 8)
    Output shape: (n_samples,)
    """
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
    localized_roughness = 0.25 * np.sin(35.0 * x1 + 17.0 * x7) * np.exp(-0.5 * ((x8 - 0.35) / 0.18) ** 2)

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
# 3. Train BoTorch GP surrogate
# ------------------------------------------------------------
def to_tensor(array: np.ndarray) -> torch.Tensor:
    return torch.as_tensor(array, dtype=TORCH_DTYPE, device=TORCH_DEVICE)


def to_numpy(tensor: torch.Tensor) -> np.ndarray:
    return tensor.detach().cpu().numpy()


from gpytorch.kernels import ScaleKernel, MaternKernel

def train_gp(X: np.ndarray, y: np.ndarray) -> SingleTaskGP:
    train_X = to_tensor(X)
    train_Y = to_tensor(y).reshape(-1, 1)
    train_Yvar = torch.full_like(train_Y, 1e-6)

    covar_module = ScaleKernel(
        MaternKernel(
            nu=2.5,
            ard_num_dims=N_INPUTS,
        )
    )

    model = SingleTaskGP(
        train_X=train_X,
        train_Y=train_Y,
        train_Yvar=train_Yvar,
        covar_module=covar_module,
        input_transform=Normalize(d=N_INPUTS),
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
            posterior = model.posterior(X[start:start + chunk_size])
            means.append(posterior.mean.squeeze(-1))
            variances.append(posterior.variance.squeeze(-1).clamp_min(0.0))
    return torch.cat(means), torch.cat(variances)


def predict_gp(
    model: SingleTaskGP,
    X: np.ndarray,
    return_std: bool = False,
) -> np.ndarray | Tuple[np.ndarray, np.ndarray]:
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
    return ", ".join(
        f"{name}={value:.3f}"
        for name, value in zip(INPUT_NAMES, point)
    )


def make_decision_weights(
    model: SingleTaskGP,
    reference_grid: np.ndarray,
    sharpness: float = 8.0,
) -> torch.Tensor:
    """
    Weight reference locations where the current surrogate predicts larger y.

    This keeps the original decision_sharpness idea: larger sharpness puts more
    IVR weight on regions currently predicted to be more important.
    """
    mean_ref, _ = posterior_mean_variance_torch(model, to_tensor(reference_grid))
    shifted = mean_ref - mean_ref.max()
    return torch.softmax(sharpness * shifted, dim=0)


def distance_filter_mask(points: np.ndarray, x_new: np.ndarray, min_distance: float) -> np.ndarray:
    """
    Keep only points at least min_distance away from x_new.
    """
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
# 5. Option 1: uncertainty + diversity
# ------------------------------------------------------------
def select_batch_uncertainty_diversity(
    model: SingleTaskGP,
    candidate_pool: np.ndarray,
    already_sampled: np.ndarray,
    n_new: int,
    min_distance: float = 0.15,
) -> np.ndarray:
    """
    Fast baseline: rank candidate-pool points by posterior standard deviation,
    then greedily keep points that satisfy the diversity distance rule.
    """
    admissible = candidate_pool[
        distance_filter_mask_many(candidate_pool, already_sampled, min_distance)
    ]
    if len(admissible) == 0:
        return np.empty((0, candidate_pool.shape[1]))

    _, std = predict_gp(model, admissible, return_std=True)
    ranked_idx = np.argsort(std)[::-1]

    selected = []
    for idx in ranked_idx:
        x_candidate = admissible[idx:idx + 1]

        if selected:
            selected_array = np.vstack(selected)
            dist_selected = np.min(np.linalg.norm(selected_array - x_candidate, axis=1))
        else:
            dist_selected = np.inf

        if dist_selected >= min_distance:
            selected.append(x_candidate.copy())

        if len(selected) >= n_new:
            break

    if not selected:
        return np.empty((0, candidate_pool.shape[1]))

    return np.vstack(selected)


# ------------------------------------------------------------
# 6. Option 2 and 3: IVR-style selection
# ------------------------------------------------------------
def evaluate_acquisition_on_pool(
    acquisition,
    candidates: torch.Tensor,
    chunk_size: int = 256,
) -> torch.Tensor:
    """
    Finite-pool acquisition optimization.

    The candidate pool is the set of points where we are allowed to run the
    simulator next. Instead of optimizing continuously over the full 8D domain,
    we evaluate the acquisition on this finite pool in chunks and take argmax.
    """
    scores = []
    with torch.no_grad():
        for start in range(0, candidates.shape[0], chunk_size):
            X = candidates[start:start + chunk_size].unsqueeze(1)
            scores.append(acquisition(X).reshape(-1).detach())
    return torch.cat(scores)


def make_native_ipv_acquisition(model: SingleTaskGP, reference_grid: np.ndarray):
    if qNegIntegratedPosteriorVariance is None:
        raise RuntimeError("qNegIntegratedPosteriorVariance is not available in this BoTorch version.")

    reference_tensor = to_tensor(reference_grid)
    try:
        return qNegIntegratedPosteriorVariance(model=model, mc_points=reference_tensor)
    except TypeError:
        return qNegIntegratedPosteriorVariance(model, reference_tensor)


def select_batch_native_integrated_variance_reduction(
    model: SingleTaskGP,
    candidate_pool: np.ndarray,
    reference_grid: np.ndarray,
    already_sampled: np.ndarray,
    n_new: int,
    min_distance: float = 0.15,
) -> np.ndarray:
    """
    Native BoTorch IPV path when qNegIntegratedPosteriorVariance is available.

    This still performs greedy q-batch selection from a finite pool so we can
    preserve the diversity rule and avoid unconstrained continuous optimization.
    """
    remaining = candidate_pool[
        distance_filter_mask_many(candidate_pool, already_sampled, min_distance)
    ]
    selected = []

    for _ in range(n_new):
        if len(remaining) == 0:
            break

        acquisition = make_native_ipv_acquisition(model, reference_grid)
        if selected and hasattr(acquisition, "set_X_pending"):
            acquisition.set_X_pending(to_tensor(np.vstack(selected)))

        scores = evaluate_acquisition_on_pool(acquisition, to_tensor(remaining))
        best_idx = int(torch.argmax(scores).item())
        x_best = remaining[best_idx:best_idx + 1]
        selected.append(x_best.copy())

        keep_mask = distance_filter_mask(remaining, x_best, min_distance)
        remaining = remaining[keep_mask]

    if not selected:
        return np.empty((0, candidate_pool.shape[1]))

    return np.vstack(selected)


def select_batch_posterior_covariance_ivr(
    model: SingleTaskGP,
    candidate_pool: np.ndarray,
    reference_grid: np.ndarray,
    already_sampled: np.ndarray,
    n_new: int,
    min_distance: float = 0.15,
    weights: Optional[torch.Tensor] = None,
) -> np.ndarray:
    """
    Greedy finite-pool IVR approximation using one fitted BoTorch GP.

    Candidate pool:
        possible next simulator runs.
    Reference grid:
        locations where we care about reducing surrogate uncertainty.

    For a fitted GP, the variance reduction from observing candidate x at a
    reference point r is approximately cov(r, x)^2 / var(x). We compute this
    for all candidates at once from the BoTorch posterior covariance, then
    greedily condition the finite covariance matrix after each selected point.
    No temporary GP refits are done inside the candidate loop.
    """
    remaining = candidate_pool[
        distance_filter_mask_many(candidate_pool, already_sampled, min_distance)
    ]
    if len(remaining) == 0:
        return np.empty((0, candidate_pool.shape[1]))

    n_ref = len(reference_grid)
    reference_tensor = to_tensor(reference_grid)
    candidate_tensor = to_tensor(remaining)
    joint_points = torch.cat([reference_tensor, candidate_tensor], dim=0)

    model.eval()
    with torch.no_grad():
        posterior = model.posterior(joint_points)
        covariance = posterior.mvn.covariance_matrix.detach()
        if covariance.ndim > 2:
            covariance = covariance.squeeze()
        covariance = covariance.to(dtype=TORCH_DTYPE, device=TORCH_DEVICE)

    if weights is None:
        weights = torch.full(
            (n_ref,),
            1.0 / n_ref,
            dtype=TORCH_DTYPE,
            device=TORCH_DEVICE,
        )
    else:
        weights = weights.to(dtype=TORCH_DTYPE, device=TORCH_DEVICE)
        weights = weights / weights.sum().clamp_min(1e-12)

    active = torch.ones(len(remaining), dtype=torch.bool, device=TORCH_DEVICE)
    selected = []

    for _ in range(n_new):
        ref_candidate_cov = covariance[:n_ref, n_ref:]
        candidate_var = covariance.diagonal()[n_ref:].clamp_min(1e-12)

        scores = (weights[:, None] * ref_candidate_cov.square()).sum(dim=0) / candidate_var
        scores = torch.where(active, scores, torch.full_like(scores, -torch.inf))

        best_local_idx = int(torch.argmax(scores).item())
        if not torch.isfinite(scores[best_local_idx]):
            break

        x_best = remaining[best_local_idx:best_local_idx + 1]
        selected.append(x_best.copy())

        # Condition the finite posterior covariance on the selected pending
        # point with fixed hyperparameters. This is the key speedup over the
        # old sklearn version: no GP hyperparameter optimization is repeated.
        best_global_idx = n_ref + best_local_idx
        denom = covariance[best_global_idx, best_global_idx].clamp_min(1e-12)
        column = covariance[:, best_global_idx:best_global_idx + 1]
        covariance = covariance - (column @ column.T) / denom
        covariance = 0.5 * (covariance + covariance.T)

        keep_np = distance_filter_mask(remaining, x_best, min_distance)
        active &= torch.as_tensor(keep_np, dtype=torch.bool, device=TORCH_DEVICE)
        active[best_local_idx] = False

    if not selected:
        return np.empty((0, candidate_pool.shape[1]))

    return np.vstack(selected)


def select_batch_integrated_variance_reduction(
    model: SingleTaskGP,
    candidate_pool: np.ndarray,
    reference_grid: np.ndarray,
    already_sampled: np.ndarray,
    n_new: int,
    min_distance: float = 0.15,
) -> np.ndarray:
    """
    Unweighted IVR.

    Prefer BoTorch's qNegIntegratedPosteriorVariance when the installed version
    provides it. If that API is missing or incompatible, fall back to the
    covariance-based finite-pool IVR approximation above.
    """
    if qNegIntegratedPosteriorVariance is not None:
        try:
            return select_batch_native_integrated_variance_reduction(
                model=model,
                candidate_pool=candidate_pool,
                reference_grid=reference_grid,
                already_sampled=already_sampled,
                n_new=n_new,
                min_distance=min_distance,
            )
        except Exception as exc:
            print(
                "Native qNegIntegratedPosteriorVariance failed; "
                f"using covariance IVR fallback instead. Reason: {exc}"
            )

    return select_batch_posterior_covariance_ivr(
        model=model,
        candidate_pool=candidate_pool,
        reference_grid=reference_grid,
        already_sampled=already_sampled,
        n_new=n_new,
        min_distance=min_distance,
    )


def select_batch_decision_aware_ivr(
    model: SingleTaskGP,
    candidate_pool: np.ndarray,
    reference_grid: np.ndarray,
    already_sampled: np.ndarray,
    n_new: int,
    min_distance: float = 0.15,
    decision_sharpness: float = 8.0,
) -> np.ndarray:
    """
    Weighted IVR.

    Native IPV does not expose arbitrary decision weights cleanly across
    BoTorch versions, so this uses the covariance IVR approximation with
    weights based on the model's current predicted mean over the reference grid.
    """
    weights = make_decision_weights(
        model=model,
        reference_grid=reference_grid,
        sharpness=decision_sharpness,
    )
    return select_batch_posterior_covariance_ivr(
        model=model,
        candidate_pool=candidate_pool,
        reference_grid=reference_grid,
        already_sampled=already_sampled,
        n_new=n_new,
        min_distance=min_distance,
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
    min_distance: float = 0.15,
    decision_sharpness: float = 8.0,
) -> np.ndarray:
    if strategy == "uncertainty_diversity":
        return select_batch_uncertainty_diversity(
            model=model,
            candidate_pool=candidate_pool,
            already_sampled=X_train,
            n_new=n_new,
            min_distance=min_distance,
        )

    if strategy == "integrated_variance_reduction":
        return select_batch_integrated_variance_reduction(
            model=model,
            candidate_pool=candidate_pool,
            reference_grid=reference_grid,
            already_sampled=X_train,
            n_new=n_new,
            min_distance=min_distance,
        )

    if strategy == "decision_aware_ivr":
        return select_batch_decision_aware_ivr(
            model=model,
            candidate_pool=candidate_pool,
            reference_grid=reference_grid,
            already_sampled=X_train,
            n_new=n_new,
            min_distance=min_distance,
            decision_sharpness=decision_sharpness,
        )

    raise ValueError(f"Unknown strategy: {strategy}")


# ------------------------------------------------------------
# 8. Plot comparison: initial vs final
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
    ax.fill_between(
        x_axis,
        y_init_pred - 1.96 * y_init_std,
        y_init_pred + 1.96 * y_init_std,
        alpha=0.2,
        label="Initial 95% band",
    )
    ax.scatter(
        X_init[:, slice_dimension],
        y_init,
        s=50,
        alpha=0.75,
        label="Initial points projected",
    )
    ax.set_title(f"Initial model\n8D test RMSE = {init_rmse:.4f}")
    ax.set_xlabel(f"{INPUT_NAMES[slice_dimension]} slice")
    ax.set_ylabel("y")
    ax.legend()

    ax = axes[1]
    ax.plot(x_axis, y_true, label="True function slice")
    ax.plot(x_axis, y_final_pred, label="Final GP mean")
    ax.fill_between(
        x_axis,
        y_final_pred - 1.96 * y_final_std,
        y_final_pred + 1.96 * y_final_std,
        alpha=0.2,
        label="Final 95% band",
    )
    ax.scatter(
        X_all[:, slice_dimension],
        y_all,
        s=25,
        alpha=0.6,
        label="All sampled points projected",
    )
    ax.scatter(
        X_init[:, slice_dimension],
        y_init,
        s=55,
        edgecolor="black",
        linewidth=0.8,
        label="Initial points projected",
    )
    ax.set_title(f"Final model\n8D test RMSE = {final_rmse:.4f}")
    ax.set_xlabel(f"{INPUT_NAMES[slice_dimension]} slice")
    ax.set_ylabel("y")
    ax.legend()

    fig.suptitle(
        f"Model comparison: {strategy_name} | other inputs fixed at midpoint",
        fontsize=14,
    )
    fig.tight_layout()

    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_dir / f"comparison_{strategy_name}.png", dpi=150)

    plt.show()
    plt.close(fig)


# ------------------------------------------------------------
# 9. Optional plot for selected points per iteration
# ------------------------------------------------------------
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
    ax.fill_between(
        x_axis,
        y_pred - 1.96 * y_std,
        y_pred + 1.96 * y_std,
        alpha=0.2,
        label="95% confidence band",
    )
    ax.scatter(
        X_train[:, slice_dimension],
        y_train,
        s=40,
        alpha=0.7,
        label="Current sampled points projected",
    )
    if new_points is not None and len(new_points) > 0:
        y_new = expensive_simulator(new_points)
        ax.scatter(
            new_points[:, slice_dimension],
            y_new,
            marker="x",
            s=100,
            label="Selected new batch projected",
        )

    ax.set_title(f"Iteration {iteration} - {strategy_name}")
    ax.set_xlabel(f"{INPUT_NAMES[slice_dimension]} slice")
    ax.set_ylabel("y")
    ax.legend()
    fig.tight_layout()

    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_dir / f"{strategy_name}_iteration_{iteration:02d}.png", dpi=150)

    plt.close(fig)


# ------------------------------------------------------------
# 10. Main
# ------------------------------------------------------------
def main():
    input_bounds = INPUT_BOUNDS
    slice_dimension = 0
    plot_output_dir = Path(__file__).with_name("gp_active_learning_plots")

    # Choose strategy here:
    # "uncertainty_diversity"
    # "integrated_variance_reduction"
    # "decision_aware_ivr"
    strategy = "integrated_variance_reduction"

    # Initial dataset
    n_init = 500
    X_init = generate_initial_points(n_init=n_init, bounds=input_bounds, seed=42)
    y_init = expensive_simulator(X_init)

    X_train = X_init.copy()
    y_train = y_init.copy()

    # Candidate pool = points where the next expensive simulator run is allowed.
    candidate_pool = latin_hypercube_sample(n_points=800, bounds=input_bounds, seed=123)
     
    # Reference grid = points where integrated uncertainty reduction is measured.
    reference_grid = latin_hypercube_sample(n_points=2000, bounds=input_bounds, seed=456)

    # Active learning settings
    n_iterations = 30
    points_per_iteration = 5
    min_distance = 0.30

    # Train and keep the initial model for comparison.
    gp_initial = train_gp(X_train, y_train)

    for it in range(n_iterations):
        print(f"\n--- Iteration {it} | strategy = {strategy} ---")

        gp_current = train_gp(X_train, y_train)

        X_new = select_new_points(
            strategy=strategy,
            model=gp_current,
            X_train=X_train,
            candidate_pool=candidate_pool,
            reference_grid=reference_grid,
            n_new=points_per_iteration,
            min_distance=min_distance,
            decision_sharpness=8.0,
        )

        print("Selected new points:")
        for idx, x in enumerate(X_new, start=1):
            print(f"  {idx:02d}: {format_design_point(x)}")

        if len(X_new) == 0:
            print("No admissible new points found; stopping early.")
            break

        if len(X_new) < points_per_iteration:
            print(
                f"Warning: requested {points_per_iteration} points, "
                f"but only found {len(X_new)} that satisfied the constraints."
            )

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

        # "Run expensive simulator" on the selected batch.
        y_new = expensive_simulator(X_new)

        # Append new data and retrain on the next loop.
        X_train = np.vstack([X_train, X_new])
        y_train = np.concatenate([y_train, y_new])

        # Remove sampled / nearby points from candidate pool.
        keep_mask = np.ones(len(candidate_pool), dtype=bool)
        for x_sel in X_new:
            keep_mask &= distance_filter_mask(candidate_pool, x_sel.reshape(1, -1), min_distance)
        candidate_pool = candidate_pool[keep_mask]

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

    # Final RMSE against a large independent LHS evaluation set.
    x_eval = latin_hypercube_sample(n_points=2000, bounds=input_bounds, seed=789)
    y_true = expensive_simulator(x_eval)
    y_pred_init = predict_gp(gp_initial, x_eval)
    y_pred_final = predict_gp(gp_final, x_eval)

    print("\n===== PERFORMANCE SUMMARY =====")
    print(f"Input dimensions: {N_INPUTS}")
    print(f"Initial RMSE: {rmse(y_true, y_pred_init):.6f}")
    print(f"Final   RMSE: {rmse(y_true, y_pred_final):.6f}")
    print(f"Total samples used: {len(X_train)}")
    print(f"Added samples: {len(X_train) - len(X_init)}")


if __name__ == "__main__":
    main()
