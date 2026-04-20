import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Optional, Tuple, List

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel


# ------------------------------------------------------------
# 1. Fake expensive simulator
# ------------------------------------------------------------
def expensive_simulator(x: np.ndarray) -> np.ndarray:
    """
    Toy stand-in for a costly simulator.
    Input shape: (n_samples, 1)
    Output shape: (n_samples,)
    """
    x = np.asarray(x).reshape(-1, 1)
    y = np.sin(1.5 * x[:, 0]) + 0.2 * np.cos(4.0 * x[:, 0]) + 0.05 * (x[:, 0] - 5.0) ** 2
    return y


# ------------------------------------------------------------
# 2. Initial data generation
# ------------------------------------------------------------
def generate_initial_points(n_init: int, x_min: float, x_max: float, seed: int = 42) -> np.ndarray:
    rng = np.random.default_rng(seed)
    x_init = rng.uniform(x_min, x_max, size=(n_init, 1))
    x_init = np.sort(x_init, axis=0)
    return x_init


# ------------------------------------------------------------
# 3. Train GP surrogate
# ------------------------------------------------------------
def train_gp(X: np.ndarray, y: np.ndarray) -> GaussianProcessRegressor:
    kernel = (
        ConstantKernel(1.0, (1e-3, 1e3))
        * RBF(length_scale=1.0, length_scale_bounds=(1e-2, 1e2))
        + WhiteKernel(noise_level=1e-6, noise_level_bounds=(1e-10, 1e-2))
    )

    gp = GaussianProcessRegressor(
        kernel=kernel,
        alpha=0.0,
        normalize_y=True,
        n_restarts_optimizer=5,
        random_state=0,
    )
    gp.fit(X, y)
    return gp


# ------------------------------------------------------------
# 4. Utility functions
# ------------------------------------------------------------
def rmse(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    return float(np.sqrt(np.mean((y_true - y_pred) ** 2)))


def make_decision_weights(
    gp: GaussianProcessRegressor,
    reference_grid: np.ndarray,
    sharpness: float = 8.0,
) -> np.ndarray:
    """
    Create weights for decision-aware active learning.

    Higher weights are assigned to regions where the current GP predicts
    larger values. This is just a toy proxy for 'regions we care about more'.

    In your real problem, this could be replaced by weights based on:
      - predicted recovery factor above some threshold
      - feasible operating region
      - low-cost / high-value scenarios
      - proximity to current Pareto front, etc.
    """
    mean_ref = gp.predict(reference_grid)
    shifted = mean_ref - np.max(mean_ref)
    weights = np.exp(sharpness * shifted)
    weights = weights / np.sum(weights)
    return weights


def distance_filter_mask(points: np.ndarray, x_new: np.ndarray, min_distance: float) -> np.ndarray:
    """
    Keep only points at least min_distance away from x_new.
    """
    distances = np.linalg.norm(points - x_new, axis=1)
    return distances >= min_distance


# ------------------------------------------------------------
# 5. Option 1: uncertainty + diversity
# ------------------------------------------------------------
def select_batch_uncertainty_diversity(
    gp: GaussianProcessRegressor,
    candidate_pool: np.ndarray,
    already_sampled: np.ndarray,
    n_new: int,
    min_distance: float = 0.15,
) -> np.ndarray:
    """
    Batch selection using highest predictive std + spacing rule.
    """
    _, std = gp.predict(candidate_pool, return_std=True)
    ranked_idx = np.argsort(std)[::-1]

    selected = []
    for idx in ranked_idx:
        x_candidate = candidate_pool[idx:idx + 1]

        dist_existing = np.min(np.linalg.norm(already_sampled - x_candidate, axis=1))

        if selected:
            selected_array = np.vstack(selected)
            dist_selected = np.min(np.linalg.norm(selected_array - x_candidate, axis=1))
        else:
            dist_selected = np.inf

        if dist_existing >= min_distance and dist_selected >= min_distance:
            selected.append(x_candidate.copy())

        if len(selected) >= n_new:
            break

    if not selected:
        return np.empty((0, candidate_pool.shape[1]))

    return np.vstack(selected)


# ------------------------------------------------------------
# 6. Option 2 and 3: greedy fantasized batch selection
# ------------------------------------------------------------
def select_batch_greedy_variance_reduction(
    X_train: np.ndarray,
    y_train: np.ndarray,
    candidate_pool: np.ndarray,
    reference_grid: np.ndarray,
    n_new: int,
    min_distance: float = 0.15,
    strategy: str = "integrated_variance_reduction",
    decision_sharpness: float = 8.0,
) -> np.ndarray:
    """
    Greedy batch active learning using fantasized updates.

    strategy:
      - "integrated_variance_reduction"
      - "decision_aware_ivr"

    How it works:
      1. Fit GP to current data.
      2. For each candidate:
         - fantasy-label it with GP mean
         - refit GP
         - compute how much variance is reduced on reference_grid
      3. Pick the best candidate.
      4. Add it as a pending/fantasized point.
      5. Repeat until batch is full.

    This is a greedy approximation, not an exact joint batch optimizer.
    """
    X_aug = X_train.copy()
    y_aug = y_train.copy()
    remaining = candidate_pool.copy()
    selected = []

    for _ in range(n_new):
        if len(remaining) == 0:
            break

        gp_current = train_gp(X_aug, y_aug)
        _, std_before = gp_current.predict(reference_grid, return_std=True)
        var_before = std_before ** 2

        if strategy == "integrated_variance_reduction":
            weights = np.ones(len(reference_grid)) / len(reference_grid)
        elif strategy == "decision_aware_ivr":
            weights = make_decision_weights(
                gp=gp_current,
                reference_grid=reference_grid,
                sharpness=decision_sharpness,
            )
        else:
            raise ValueError(f"Unknown strategy: {strategy}")

        best_score = -np.inf
        best_idx = None
        best_fantasy_y = None

        for idx in range(len(remaining)):
            x_cand = remaining[idx:idx + 1]

            # Respect spacing against current augmented set
            dist_existing = np.min(np.linalg.norm(X_aug - x_cand, axis=1))
            if dist_existing < min_distance:
                continue

            # Fantasy label from current GP mean
            y_fantasy = gp_current.predict(x_cand)[0]

            # Temporary augmented set
            X_tmp = np.vstack([X_aug, x_cand])
            y_tmp = np.concatenate([y_aug, [y_fantasy]])

            gp_tmp = train_gp(X_tmp, y_tmp)
            _, std_after = gp_tmp.predict(reference_grid, return_std=True)
            var_after = std_after ** 2

            # Score = weighted average variance reduction
            score = np.sum(weights * (var_before - var_after))

            if score > best_score:
                best_score = score
                best_idx = idx
                best_fantasy_y = y_fantasy

        if best_idx is None:
            break

        x_best = remaining[best_idx:best_idx + 1]
        selected.append(x_best.copy())

        # Add fantasized point so the next batch point accounts for pending points
        X_aug = np.vstack([X_aug, x_best])
        y_aug = np.concatenate([y_aug, [best_fantasy_y]])

        # Remove nearby points from candidate pool
        keep_mask = distance_filter_mask(remaining, x_best, min_distance)
        remaining = remaining[keep_mask]

    if not selected:
        return np.empty((0, candidate_pool.shape[1]))

    return np.vstack(selected)


# ------------------------------------------------------------
# 7. Strategy wrapper
# ------------------------------------------------------------
def select_new_points(
    strategy: str,
    X_train: np.ndarray,
    y_train: np.ndarray,
    candidate_pool: np.ndarray,
    reference_grid: np.ndarray,
    n_new: int,
    min_distance: float = 0.15,
    decision_sharpness: float = 8.0,
) -> np.ndarray:
    gp = train_gp(X_train, y_train)

    if strategy == "uncertainty_diversity":
        return select_batch_uncertainty_diversity(
            gp=gp,
            candidate_pool=candidate_pool,
            already_sampled=X_train,
            n_new=n_new,
            min_distance=min_distance,
        )

    elif strategy in ["integrated_variance_reduction", "decision_aware_ivr"]:
        return select_batch_greedy_variance_reduction(
            X_train=X_train,
            y_train=y_train,
            candidate_pool=candidate_pool,
            reference_grid=reference_grid,
            n_new=n_new,
            min_distance=min_distance,
            strategy=strategy,
            decision_sharpness=decision_sharpness,
        )

    else:
        raise ValueError(f"Unknown strategy: {strategy}")


# ------------------------------------------------------------
# 8. Plot comparison: initial vs final
# ------------------------------------------------------------
def plot_initial_vs_final(
    gp_initial: GaussianProcessRegressor,
    gp_final: GaussianProcessRegressor,
    X_init: np.ndarray,
    y_init: np.ndarray,
    X_all: np.ndarray,
    y_all: np.ndarray,
    x_min: float,
    x_max: float,
    strategy_name: str,
    output_dir: Optional[Path] = None,
):
    x_plot = np.linspace(x_min, x_max, 1000).reshape(-1, 1)
    y_true = expensive_simulator(x_plot)

    y_init_pred, y_init_std = gp_initial.predict(x_plot, return_std=True)
    y_final_pred, y_final_std = gp_final.predict(x_plot, return_std=True)

    init_rmse = rmse(y_true, y_init_pred)
    final_rmse = rmse(y_true, y_final_pred)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    # Initial model
    ax = axes[0]
    ax.plot(x_plot[:, 0], y_true, label="True function")
    ax.plot(x_plot[:, 0], y_init_pred, label="Initial GP mean")
    ax.fill_between(
        x_plot[:, 0],
        y_init_pred - 1.96 * y_init_std,
        y_init_pred + 1.96 * y_init_std,
        alpha=0.2,
        label="Initial 95% band",
    )
    ax.scatter(X_init[:, 0], y_init, s=50, label="Initial points")
    ax.set_title(f"Initial model\nRMSE = {init_rmse:.4f}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.legend()

    # Final model
    ax = axes[1]
    ax.plot(x_plot[:, 0], y_true, label="True function")
    ax.plot(x_plot[:, 0], y_final_pred, label="Final GP mean")
    ax.fill_between(
        x_plot[:, 0],
        y_final_pred - 1.96 * y_final_std,
        y_final_pred + 1.96 * y_final_std,
        alpha=0.2,
        label="Final 95% band",
    )
    ax.scatter(X_all[:, 0], y_all, s=25, label="All sampled points")
    ax.scatter(
        X_init[:, 0],
        y_init,
        s=55,
        edgecolor="black",
        linewidth=0.8,
        label="Initial points",
    )
    ax.set_title(f"Final model\nRMSE = {final_rmse:.4f}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.legend()

    fig.suptitle(f"Model comparison: {strategy_name}", fontsize=14)
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
    gp: GaussianProcessRegressor,
    X_train: np.ndarray,
    y_train: np.ndarray,
    new_points: np.ndarray,
    x_min: float,
    x_max: float,
    strategy_name: str,
    output_dir: Optional[Path] = None,
):
    x_plot = np.linspace(x_min, x_max, 1000).reshape(-1, 1)
    y_true = expensive_simulator(x_plot)
    y_pred, y_std = gp.predict(x_plot, return_std=True)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(x_plot[:, 0], y_true, label="True function")
    ax.plot(x_plot[:, 0], y_pred, label="GP mean")
    ax.fill_between(
        x_plot[:, 0],
        y_pred - 1.96 * y_std,
        y_pred + 1.96 * y_std,
        alpha=0.2,
        label="95% confidence band",
    )
    ax.scatter(X_train[:, 0], y_train, s=40, label="Current sampled points")
    if new_points is not None and len(new_points) > 0:
        y_new = expensive_simulator(new_points)
        ax.scatter(new_points[:, 0], y_new, marker="x", s=100, label="Selected new batch")

    ax.set_title(f"Iteration {iteration} - {strategy_name}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.legend()
    fig.tight_layout()

    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_dir / f"{strategy_name}_iteration_{iteration:02d}.png", dpi=150)

    # plt.show()
    # plt.close(fig)


# ------------------------------------------------------------
# 10. Main
# ------------------------------------------------------------
def main():
    x_min, x_max = 0.0, 10.0
    plot_output_dir = Path(__file__).with_name("gp_active_learning_plots")

    # Choose strategy here:
    # "uncertainty_diversity"
    # "integrated_variance_reduction"
    # "decision_aware_ivr"
    strategy = "integrated_variance_reduction"

    # Initial dataset
    n_init = 6
    X_init = generate_initial_points(n_init=n_init, x_min=x_min, x_max=x_max, seed=42)
    y_init = expensive_simulator(X_init)

    X_train = X_init.copy()
    y_train = y_init.copy()

    # Candidate pool: where the algorithm is allowed to choose next simulations
    candidate_pool = np.linspace(x_min, x_max, 400).reshape(-1, 1)

    # Reference grid: where global uncertainty reduction is evaluated
    reference_grid = np.linspace(x_min, x_max, 500).reshape(-1, 1)

    # Active learning settings
    n_iterations = 4
    points_per_iteration = 20
    min_distance = 0.20

    # Train and keep the initial model for comparison
    gp_initial = train_gp(X_train, y_train)

    for it in range(n_iterations):
        print(f"\n--- Iteration {it} | strategy = {strategy} ---")

        gp_current = train_gp(X_train, y_train)

        X_new = select_new_points(
            strategy=strategy,
            X_train=X_train,
            y_train=y_train,
            candidate_pool=candidate_pool,
            reference_grid=reference_grid,
            n_new=points_per_iteration,
            min_distance=min_distance,
            decision_sharpness=8.0,
        )

        print("Selected new points:")
        for x in X_new[:, 0]:
            print(f"  x = {x:.4f}")

        if len(X_new) == 0:
            print("No admissible new points found; stopping early.")
            break

        if len(X_new) < points_per_iteration:
            print(
                f"Warning: requested {points_per_iteration} points, "
                f"but only found {len(X_new)} that satisfied the constraints."
            )

        # Optional per-iteration plot
        plot_iteration(
            iteration=it,
            gp=gp_current,
            X_train=X_train,
            y_train=y_train,
            new_points=X_new,
            x_min=x_min,
            x_max=x_max,
            strategy_name=strategy,
            output_dir=plot_output_dir,
        )

        # "Run expensive simulator" on the whole batch
        y_new = expensive_simulator(X_new)

        # Append new data
        X_train = np.vstack([X_train, X_new])
        y_train = np.concatenate([y_train, y_new])

        # Remove sampled / nearby points from candidate pool
        keep_mask = np.ones(len(candidate_pool), dtype=bool)
        for x_sel in X_new:
            keep_mask &= distance_filter_mask(candidate_pool, x_sel.reshape(1, -1), min_distance)
        candidate_pool = candidate_pool[keep_mask]

    # Final model
    gp_final = train_gp(X_train, y_train)

    # Final comparison plot: base pool vs enriched pool
    plot_initial_vs_final(
        gp_initial=gp_initial,
        gp_final=gp_final,
        X_init=X_init,
        y_init=y_init,
        X_all=X_train,
        y_all=y_train,
        x_min=x_min,
        x_max=x_max,
        strategy_name=strategy,
        output_dir=plot_output_dir,
    )

    # Print final RMSE numbers
    x_eval = np.linspace(x_min, x_max, 1000).reshape(-1, 1)
    y_true = expensive_simulator(x_eval)
    y_pred_init = gp_initial.predict(x_eval)
    y_pred_final = gp_final.predict(x_eval)

    print("\n===== PERFORMANCE SUMMARY =====")
    print(f"Initial RMSE: {rmse(y_true, y_pred_init):.6f}")
    print(f"Final   RMSE: {rmse(y_true, y_pred_final):.6f}")
    print(f"Total samples used: {len(X_train)}")
    print(f"Added samples: {len(X_train) - len(X_init)}")


if __name__ == "__main__":
    main()