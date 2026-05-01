"""
Active-learning workflow v3 — realistic toy function + multi-scale GP.

Changes from v2:
  - Simplified expensive_simulator that is learnable with 500 pts in 8D
  - Features have width > 0.3 (resolvable with ~2 pts/axis)
  - Still nonlinear: smooth trends, cross-terms, moderate sigmoid, curvature
  - Output bounded roughly [0, 1] to mimic Recovery Factor
  - Demonstrates that active learning genuinely improves the surrogate
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
        "Requires PyTorch, GPyTorch, and BoTorch. "
        "Install via `pip install torch gpytorch botorch`."
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
# 1. Realistic toy simulator (learnable in 8D with ~500 pts)
# ------------------------------------------------------------
def expensive_simulator(x: np.ndarray) -> np.ndarray:
    """
    Simplified 8-input simulator mimicking Recovery Factor.
    
    Output is roughly in [0, 1].
    All features have characteristic width > 0.3, resolvable
    with the ~2 pts/axis effective resolution of 500 points in 8D.
    
    Components:
      - baseline + smooth sinusoidal trends (period ~1)
      - quadratic curvature
      - cross-terms (2-way interactions)
      - moderate sigmoid transition (steepness=6, width ~0.3)
      - gentle localized bump (width 0.3)
    """
    x = np.asarray(x, dtype=float)
    if x.ndim == 1:
        x = x.reshape(1, -1)
    if x.ndim != 2 or x.shape[1] != N_INPUTS:
        raise ValueError(f"Expected (n, {N_INPUTS}), got {x.shape}")

    x1, x2, x3, x4, x5, x6, x7, x8 = x.T

    # Baseline around 0.5
    baseline = 0.45

    # Smooth sinusoidal trends (period ~1, resolvable)
    trend = (
        0.15 * np.sin(2.0 * np.pi * x1)
        + 0.10 * np.sin(2.0 * np.pi * x2)
        + 0.08 * np.cos(2.0 * np.pi * (x7 - 0.3))
    )

    # Quadratic curvature — porosity and permeability effects
    curvature = (
        -0.25 * (x4 - 0.55) ** 2   # porosity sweet spot
        + 0.18 * (x5 - 0.3)         # permeability benefit
        - 0.15 * (x3 - 0.5) ** 2    # production rate penalty at extremes
    )

    # Cross-terms — 2-way and 3-way interactions
    interactions = (
        0.22 * (x1 - 0.5) * (x5 - 0.5)     # inj_rate × permeability
        - 0.18 * (x2 - 0.5) * (x6 - 0.5)    # inj_rate_op × pressure
        + 0.14 * (x3 - 0.5) * (x8 - 0.5)    # production × cycle_length
        + 0.12 * (x1 - 0.5) * (x4 - 0.5) * (x6 - 0.5)  # 3-way interaction
    )

    # Moderate sigmoid — pressure threshold (width ~0.15, steepness=10)
    # Represents a phase behavior / breakthrough transition
    pressure_effect = 0.20 / (1.0 + np.exp(-10.0 * (x6 + 0.3 * x2 - 0.65)))

    # Localized bump (width 0.20, needs some data nearby to resolve)
    # Represents an optimal operating window
    optimal_window = 0.18 * np.exp(
        -0.5 * (
            ((x1 - 0.65) / 0.20) ** 2
            + ((x5 - 0.35) / 0.20) ** 2
        )
    )

    # Temperature effect — mild nonlinear
    temp_effect = 0.10 * x7 * (1.0 - x7)

    # Cycle length interaction — longer cycles reduce RF slightly
    cycle_effect = -0.08 * x8 * (x3 - 0.5)

    y = baseline + trend + curvature + interactions + pressure_effect + optimal_window + temp_effect + cycle_effect

    # Soft clamp to [0, 1] to mimic RF
    y = np.clip(y, 0.0, 1.0)

    return y


# ------------------------------------------------------------
# 2. Sampling
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


# ------------------------------------------------------------
# 3. GP training — multi-scale kernel + learned nugget
# ------------------------------------------------------------
def to_tensor(a: np.ndarray) -> torch.Tensor:
    return torch.as_tensor(a, dtype=TORCH_DTYPE, device=TORCH_DEVICE)

def to_numpy(t: torch.Tensor) -> np.ndarray:
    return t.detach().cpu().numpy()

_FIXED_BOUNDS = torch.tensor(INPUT_BOUNDS.T, dtype=TORCH_DTYPE, device=TORCH_DEVICE)


def train_gp(X: np.ndarray, y: np.ndarray) -> SingleTaskGP:
    train_X = to_tensor(X)
    train_Y = to_tensor(y).reshape(-1, 1)

    covar_module = AdditiveKernel(
        ScaleKernel(MaternKernel(nu=2.5, ard_num_dims=N_INPUTS)),
        ScaleKernel(MaternKernel(nu=1.5, ard_num_dims=N_INPUTS)),
    )

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


def posterior_mean_variance(model, X_tensor, chunk_size=2048):
    means, variances = [], []
    model.eval()
    with torch.no_grad(), gpytorch.settings.fast_pred_var():
        for s in range(0, X_tensor.shape[0], chunk_size):
            post = model.posterior(X_tensor[s:s + chunk_size])
            means.append(post.mean.squeeze(-1))
            variances.append(post.variance.squeeze(-1).clamp_min(0.0))
    return torch.cat(means), torch.cat(variances)


def predict_gp(model, X, return_std=False):
    mean, var = posterior_mean_variance(model, to_tensor(X))
    m = to_numpy(mean)
    if return_std:
        return m, to_numpy(torch.sqrt(var))
    return m


# ------------------------------------------------------------
# 4. Utilities
# ------------------------------------------------------------
def rmse(y_true, y_pred):
    return float(np.sqrt(np.mean((y_true - y_pred) ** 2)))


def make_slice_points(bounds, dim, n=1000, base=None):
    bounds = np.asarray(bounds, dtype=float)
    if base is None:
        base = 0.5 * (bounds[:, 0] + bounds[:, 1])
    pts = np.tile(base, (n, 1))
    pts[:, dim] = np.linspace(bounds[dim, 0], bounds[dim, 1], n)
    return pts


def format_point(pt):
    return ", ".join(f"{n}={v:.3f}" for n, v in zip(INPUT_NAMES, pt))


def dist_filter(points, x_new, min_d):
    return np.linalg.norm(points - x_new, axis=1) >= min_d


def dist_filter_many(points, blocked, min_d):
    if len(points) == 0:
        return np.empty(0, dtype=bool)
    if blocked is None or len(blocked) == 0:
        return np.ones(len(points), dtype=bool)
    dists = np.linalg.norm(points[:, None, :] - blocked[None, :, :], axis=2)
    return np.min(dists, axis=1) >= min_d


# ------------------------------------------------------------
# 5. Acquisition: uncertainty + diversity
# ------------------------------------------------------------
def select_batch_uncertainty_diversity(model, candidates, sampled, n_new,
                                       min_d_hist=0.05, min_d_batch=0.12):
    admissible = candidates[dist_filter_many(candidates, sampled, min_d_hist)]
    if len(admissible) == 0:
        return np.empty((0, candidates.shape[1]))

    _, std = predict_gp(model, admissible, return_std=True)
    ranked = np.argsort(std)[::-1]

    selected = []
    for idx in ranked:
        xc = admissible[idx:idx + 1]
        if selected:
            d = np.min(np.linalg.norm(np.vstack(selected) - xc, axis=1))
        else:
            d = np.inf
        if d >= min_d_batch:
            selected.append(xc.copy())
        if len(selected) >= n_new:
            break

    return np.vstack(selected) if selected else np.empty((0, candidates.shape[1]))


# ------------------------------------------------------------
# 6. Acquisition: chunked IVR
# ------------------------------------------------------------
def select_batch_chunked_ivr(model, candidates, ref_grid, sampled, n_new,
                              min_d_hist=0.05, min_d_batch=0.12,
                              weights=None, chunk=500):
    remaining = candidates[dist_filter_many(candidates, sampled, min_d_hist)]
    if len(remaining) == 0:
        return np.empty((0, candidates.shape[1]))

    n_ref = len(ref_grid)
    ref_t = to_tensor(ref_grid)

    if weights is None:
        w = torch.full((n_ref,), 1.0 / n_ref, dtype=TORCH_DTYPE, device=TORCH_DEVICE)
    else:
        w = weights.to(dtype=TORCH_DTYPE, device=TORCH_DEVICE)
        w = w / w.sum().clamp_min(1e-12)

    active = np.ones(len(remaining), dtype=bool)
    selected = []

    for _ in range(n_new):
        active_idx = np.where(active)[0]
        if len(active_idx) == 0:
            break

        best_score, best_gi = -float("inf"), -1

        for cs in range(0, len(active_idx), chunk):
            ci = active_idx[cs:cs + chunk]
            ct = to_tensor(remaining[ci])
            joint = torch.cat([ref_t, ct], dim=0)

            with torch.no_grad():
                post = model.posterior(joint)
                cov = post.mvn.covariance_matrix.detach().squeeze()
                K_rc = cov[:n_ref, n_ref:]
                c_var = cov.diagonal()[n_ref:].clamp_min(1e-12)

            scores = (w[:, None] * K_rc.square()).sum(dim=0) / c_var
            bi = int(torch.argmax(scores).item())
            if scores[bi].item() > best_score:
                best_score = scores[bi].item()
                best_gi = ci[bi]

        if best_gi < 0 or not np.isfinite(best_score):
            break

        xb = remaining[best_gi:best_gi + 1]
        selected.append(xb.copy())
        active &= dist_filter(remaining, xb, min_d_batch)
        active[best_gi] = False

    return np.vstack(selected) if selected else np.empty((0, candidates.shape[1]))


# ------------------------------------------------------------
# 7. Strategy wrapper
# ------------------------------------------------------------
def select_new_points(strategy, model, X_train, candidates, ref_grid, n_new,
                      min_d_hist=0.05, min_d_batch=0.12):
    if strategy == "uncertainty_diversity":
        return select_batch_uncertainty_diversity(
            model, candidates, X_train, n_new, min_d_hist, min_d_batch)
    if strategy == "integrated_variance_reduction":
        return select_batch_chunked_ivr(
            model, candidates, ref_grid, X_train, n_new, min_d_hist, min_d_batch)
    raise ValueError(f"Unknown strategy: {strategy}")


# ------------------------------------------------------------
# 8. Plotting
# ------------------------------------------------------------
def plot_initial_vs_final(gp_init, gp_final, X_init, y_init, X_all, y_all,
                          bounds, strategy, dim=0, outdir=None):
    xp = make_slice_points(bounds, dim)
    xa = xp[:, dim]
    yt = expensive_simulator(xp)

    yi_pred, yi_std = predict_gp(gp_init, xp, return_std=True)
    yf_pred, yf_std = predict_gp(gp_final, xp, return_std=True)

    xe = latin_hypercube_sample(2000, bounds, seed=900)
    ye = expensive_simulator(xe)
    ri = rmse(ye, predict_gp(gp_init, xe))
    rf = rmse(ye, predict_gp(gp_final, xe))

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    ax = axes[0]
    ax.plot(xa, yt, label="True function")
    ax.plot(xa, yi_pred, label="GP mean")
    ax.fill_between(xa, yi_pred - 1.96 * yi_std, yi_pred + 1.96 * yi_std,
                     alpha=0.2, label="95% band")
    ax.scatter(X_init[:, dim], y_init, s=30, alpha=0.5, label="Init points")
    ax.set_title(f"Initial model — RMSE = {ri:.4f}")
    ax.set_xlabel(INPUT_NAMES[dim])
    ax.set_ylabel("RF")
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.plot(xa, yt, label="True function")
    ax.plot(xa, yf_pred, label="GP mean")
    ax.fill_between(xa, yf_pred - 1.96 * yf_std, yf_pred + 1.96 * yf_std,
                     alpha=0.2, label="95% band")
    ax.scatter(X_all[:, dim], y_all, s=15, alpha=0.4, label="All points")
    ax.scatter(X_init[:, dim], y_init, s=40, edgecolor="k", linewidth=0.6,
               label="Init points")
    ax.set_title(f"Final model — RMSE = {rf:.4f}")
    ax.set_xlabel(INPUT_NAMES[dim])
    ax.set_ylabel("RF")
    ax.legend(fontsize=8)

    fig.suptitle(f"v3: {strategy} | RMSE {ri:.4f} → {rf:.4f} "
                 f"({(1-rf/ri)*100:.1f}% reduction)", fontsize=13)
    fig.tight_layout()
    if outdir:
        outdir.mkdir(parents=True, exist_ok=True)
        fig.savefig(outdir / f"v3_comparison_{strategy}.png", dpi=150)
    plt.close(fig)


def plot_iteration(it, gp, X_train, y_train, X_new, bounds, strategy, dim=0, outdir=None):
    xp = make_slice_points(bounds, dim)
    xa = xp[:, dim]
    yt = expensive_simulator(xp)
    yp, ys = predict_gp(gp, xp, return_std=True)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(xa, yt, label="True")
    ax.plot(xa, yp, label="GP mean")
    ax.fill_between(xa, yp - 1.96 * ys, yp + 1.96 * ys, alpha=0.2, label="95% band")
    ax.scatter(X_train[:, dim], y_train, s=30, alpha=0.5, label="Data")
    if X_new is not None and len(X_new) > 0:
        ax.scatter(X_new[:, dim], expensive_simulator(X_new),
                   marker="x", s=100, c="red", label="New batch")
    ax.set_title(f"Iter {it} — {strategy} (v3)")
    ax.set_xlabel(INPUT_NAMES[dim])
    ax.set_ylabel("RF")
    ax.legend(fontsize=8)
    fig.tight_layout()
    if outdir:
        outdir.mkdir(parents=True, exist_ok=True)
        fig.savefig(outdir / f"v3_{strategy}_iter_{it:02d}.png", dpi=150)
    plt.close(fig)


def plot_rmse_history(hist, strategy, outdir=None):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(range(len(hist)), hist, "o-", markersize=5)
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Test RMSE")
    ax.set_title(f"RMSE convergence — {strategy} (v3)")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    if outdir:
        outdir.mkdir(parents=True, exist_ok=True)
        fig.savefig(outdir / f"v3_rmse_{strategy}.png", dpi=150)
    plt.close(fig)


# ------------------------------------------------------------
# 9. Diagnostics
# ------------------------------------------------------------
def print_gp_diagnostics(model):
    print("\n  GP kernel diagnostics:")
    for i, k in enumerate(model.covar_module.kernels):
        ls = k.base_kernel.lengthscale.detach().cpu().numpy().flatten()
        os_ = k.outputscale.detach().cpu().item()
        print(f"    Kernel {i} (Matern nu={k.base_kernel.nu}): outputscale={os_:.4f}")
        for j, l in enumerate(ls):
            print(f"      {INPUT_NAMES[j]:20s}: lengthscale={l:.4f}")
    noise = model.likelihood.noise.detach().cpu().item()
    print(f"    Learned noise: {noise:.2e}")


# ------------------------------------------------------------
# 10. Main
# ------------------------------------------------------------
def main():
    bounds = INPUT_BOUNDS
    dim = 0
    outdir = Path(__file__).with_name("gp_active_learning_plots")

    strategy = "uncertainty_diversity"

    # Initial dataset
    n_init = 200
    X_init = latin_hypercube_sample(n_init, bounds, seed=42)
    y_init = expensive_simulator(X_init)

    X_train = X_init.copy()
    y_train = y_init.copy()

    # Reference grid (for IVR)
    ref_grid = latin_hypercube_sample(2000, bounds, seed=456)

    # Fixed eval set
    x_eval = latin_hypercube_sample(2000, bounds, seed=789)
    y_eval = expensive_simulator(x_eval)

    # Settings
    n_iter = 30
    batch_size = 10
    min_d_hist = 0.05
    min_d_batch = 0.10

    # Function statistics
    print(f"Function range: [{y_eval.min():.3f}, {y_eval.max():.3f}]")
    print(f"Function mean:  {y_eval.mean():.3f}")
    print(f"Function std:   {y_eval.std():.3f}")

    # Train initial GP
    gp_init = train_gp(X_train, y_train)
    rmse_hist = [rmse(y_eval, predict_gp(gp_init, x_eval))]
    print(f"\nInitial RMSE: {rmse_hist[0]:.6f}")
    print_gp_diagnostics(gp_init)

    for it in range(n_iter):
        print(f"\n--- Iteration {it} | {strategy} ---")

        gp = train_gp(X_train, y_train)

        # Fresh candidate pool each iteration
        candidates = latin_hypercube_sample(5000, bounds, seed=1000 + it)

        X_new = select_new_points(
            strategy, gp, X_train, candidates, ref_grid,
            n_new=batch_size, min_d_hist=min_d_hist, min_d_batch=min_d_batch,
        )

        if len(X_new) == 0:
            print("No admissible points; stopping.")
            break

        print(f"  Selected {len(X_new)} points")

        plot_iteration(it, gp, X_train, y_train, X_new, bounds, strategy,
                       dim=dim, outdir=outdir)

        y_new = expensive_simulator(X_new)
        X_train = np.vstack([X_train, X_new])
        y_train = np.concatenate([y_train, y_new])

        cur_rmse = rmse(y_eval, predict_gp(train_gp(X_train, y_train), x_eval))
        rmse_hist.append(cur_rmse)
        print(f"  RMSE: {cur_rmse:.6f}")

    gp_final = train_gp(X_train, y_train)

    plot_initial_vs_final(gp_init, gp_final, X_init, y_init, X_all=X_train,
                          y_all=y_train, bounds=bounds, strategy=strategy,
                          dim=dim, outdir=outdir)
    plot_rmse_history(rmse_hist, strategy, outdir=outdir)
    print_gp_diagnostics(gp_final)

    # Multi-slice comparison
    fig, axes = plt.subplots(2, 4, figsize=(18, 8))
    for i, ax in enumerate(axes.flat):
        xp = make_slice_points(bounds, i)
        xa = xp[:, i]
        yt = expensive_simulator(xp)
        yp, ys = predict_gp(gp_final, xp, return_std=True)
        ax.plot(xa, yt, label="True")
        ax.plot(xa, yp, label="GP")
        ax.fill_between(xa, yp - 1.96 * ys, yp + 1.96 * ys, alpha=0.2)
        ax.set_title(INPUT_NAMES[i], fontsize=10)
        ax.set_ylabel("RF")
        if i == 0:
            ax.legend(fontsize=7)
    fig.suptitle(f"v3 final: all 8 input slices | RMSE = {rmse_hist[-1]:.4f}", fontsize=13)
    fig.tight_layout()
    if outdir:
        outdir.mkdir(parents=True, exist_ok=True)
        fig.savefig(outdir / "v3_all_slices.png", dpi=150)
    plt.close(fig)

    print(f"\n{'='*50}")
    print(f"PERFORMANCE SUMMARY (v3)")
    print(f"{'='*50}")
    print(f"Strategy:       {strategy}")
    print(f"Initial points: {n_init}")
    print(f"Added points:   {len(X_train) - n_init}")
    print(f"Total points:   {len(X_train)}")
    print(f"Initial RMSE:   {rmse_hist[0]:.6f}")
    print(f"Final RMSE:     {rmse_hist[-1]:.6f}")
    print(f"Reduction:      {(1 - rmse_hist[-1]/rmse_hist[0])*100:.1f}%")


if __name__ == "__main__":
    main()
