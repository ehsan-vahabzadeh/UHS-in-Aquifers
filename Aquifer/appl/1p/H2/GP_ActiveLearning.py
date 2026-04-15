import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Optional

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
    """
    Fit a Gaussian Process model.
    """
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
# 4. Select new points by uncertainty
# ------------------------------------------------------------
def select_new_points_by_uncertainty(
    gp: GaussianProcessRegressor,
    candidate_pool: np.ndarray,
    already_sampled: np.ndarray,
    n_new: int,
    min_distance: float = 0.15,
) -> np.ndarray:
    """
    Select new points where GP predictive std is highest.

    We also avoid selecting points too close to already-sampled points,
    otherwise the model will waste queries on nearly identical inputs.
    """
    _, std = gp.predict(candidate_pool, return_std=True)

    # Sort candidate indices from highest uncertainty to lowest
    ranked_idx = np.argsort(std)[::-1]

    selected = []
    for idx in ranked_idx:
        x_candidate = candidate_pool[idx:idx + 1]

        # Distance to existing data
        dist_existing = np.min(np.abs(already_sampled - x_candidate))

        # Distance to already selected new points in this round
        if selected:
            selected_array = np.array(selected).reshape(-1, 1)
            dist_selected = np.min(np.abs(selected_array - x_candidate))
        else:
            dist_selected = np.inf

        if dist_existing >= min_distance and dist_selected >= min_distance:
            selected.append(float(x_candidate[0, 0]))

        if len(selected) >= n_new:
            break

    return np.array(selected).reshape(-1, 1)


# ------------------------------------------------------------
# 5. Plot current state
# ------------------------------------------------------------
def plot_iteration(
    iteration: int,
    gp: GaussianProcessRegressor,
    X_train: np.ndarray,
    y_train: np.ndarray,
    candidate_pool: np.ndarray,
    new_points: np.ndarray,
    x_min: float,
    x_max: float,
    output_dir: Optional[Path] = None,
    show_plot: bool = False,
    keep_open: bool = False,
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
    ax.scatter(X_train[:, 0], y_train, s=50, label="Sampled points")
    if new_points is not None and len(new_points) > 0:
        y_new = expensive_simulator(new_points)
        ax.scatter(new_points[:, 0], y_new, marker="x", s=100, label="Newly selected points")

    ax.set_title(f"Iteration {iteration}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.legend()
    fig.tight_layout()

    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_dir / f"iteration_{iteration:02d}.png", dpi=150)

    if show_plot:
        plt.show(block=keep_open)
        if not keep_open:
            plt.pause(0.1)

    plt.close(fig)


# ------------------------------------------------------------
# 6. Main active learning loop
# ------------------------------------------------------------
def main():
    x_min, x_max = 0.0, 10.0
    plot_output_dir = Path(__file__).with_name("gp_active_learning_plots")
    show_plots = False  # Switch to True if you want interactive windows during the run.

    # Initial dataset
    n_init = 6
    X_train = generate_initial_points(n_init=n_init, x_min=x_min, x_max=x_max, seed=42)
    y_train = expensive_simulator(X_train)

    # Candidate pool where we ask: "where should we sample next?"
    candidate_pool = np.linspace(x_min, x_max, 400).reshape(-1, 1)

    # Active learning settings
    n_iterations = 5
    points_per_iteration = 2

    for it in range(n_iterations):
        print(f"\n--- Iteration {it} ---")

        # Train surrogate
        gp = train_gp(X_train, y_train)

        # Select new points based on uncertainty
        X_new = select_new_points_by_uncertainty(
            gp=gp,
            candidate_pool=candidate_pool,
            already_sampled=X_train,
            n_new=points_per_iteration,
            min_distance=0.20,
        )

        print("Selected new points:")
        for x in X_new[:, 0]:
            print(f"  x = {x:.4f}")

        if len(X_new) == 0:
            print("No admissible new points found; stopping active learning early.")
            break

        if len(X_new) < points_per_iteration:
            print(
                f"Warning: requested {points_per_iteration} point(s), "
                f"but only found {len(X_new)} candidate(s) that satisfy min_distance."
            )

        # Plot before adding the new points
        # plot_iteration(
        #     iteration=it,
        #     gp=gp,
        #     X_train=X_train,
        #     y_train=y_train,
        #     candidate_pool=candidate_pool,
        #     new_points=X_new,
        #     x_min=x_min,
        #     x_max=x_max,
        #     output_dir=plot_output_dir,
        #     show_plot=show_plots,
        # )

        # "Run simulation" at new points
        y_new = expensive_simulator(X_new)

        # Append to training set
        X_train = np.vstack([X_train, X_new])
        y_train = np.concatenate([y_train, y_new])

    # Final model after all iterations
    gp_final = train_gp(X_train, y_train)
    # plot_iteration(
    #     iteration=n_iterations,
    #     gp=gp_final,
    #     X_train=X_train,
    #     y_train=y_train,
    #     candidate_pool=candidate_pool,
    #     new_points=None,
    #     x_min=x_min,
    #     x_max=x_max,
    #     output_dir=plot_output_dir,
    #     show_plot=show_plots,
    #     keep_open=show_plots,
    # )


if __name__ == "__main__":
    main()
