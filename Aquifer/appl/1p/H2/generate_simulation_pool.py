#!/usr/bin/env python3
"""
Generate a 500-case simulation input pool from real reservoir data.

The reservoir-property part is selected from real rows, not sampled
independently. This keeps permeability, porosity, pressure, and temperature on
the observed data manifold while still giving sparse permeability regimes a
chance to appear in the final pool.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import qmc
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT_CSV = SCRIPT_DIR / "matched_reservoir_parameters.csv"
DEFAULT_OUTPUT_CSV = SCRIPT_DIR / "simulation_input_pool_500.csv"
DEFAULT_PLOT_DIR = SCRIPT_DIR / "simulation_pool_diagnostics"

REQUIRED_COLUMN_ALIASES = {
    "permeability_md": ("Permeability [mD]", "permeability_md", "permeability", "perm_md"),
    "porosity": ("Porosity [-]", "porosity"),
    "pressure_mpa": ("Pore Pressure [MPa]", "pressure_mpa", "pore_pressure_mpa", "pressure"),
    "temperature_c": ("Formation Temp [C]", "temperature_c", "formation_temp_c", "temperature"),
}

OPTIONAL_COLUMN_ALIASES = {
    "field_name": ("Field Name", "field_name", "description", "reservoir_name"),
    "pore_volume": ("Pore Volume", "pore_volume", "porevolume"),
    "number_of_wells": ("Number of Wells", "number_of_wells", "totalwells"),
    "latitude": ("Latitude", "latitude", "lat"),
    "longitude": ("Longitude", "longitude", "lon"),
}

PROPERTY_COLUMNS = ["permeability_md", "porosity", "pressure_mpa", "temperature_c"]
CG_TYPES = ["H2", "N2", "CO2", "CH4"]


@dataclass(frozen=True)
class Config:
    input_csv: Path
    output_csv: Path
    plot_dir: Path
    n_cases: int
    n_representatives: int | None
    random_seed: int
    n_permeability_strata: int
    density_power: float
    flow_min: float
    flow_max: float
    cycle_min: float
    cycle_max: float
    cg_ratio_min: float
    cg_ratio_max: float


def resolve_column(df: pd.DataFrame, aliases: tuple[str, ...], required: bool = True) -> str | None:
    lower_lookup = {column.casefold(): column for column in df.columns}
    for alias in aliases:
        if alias.casefold() in lower_lookup:
            return lower_lookup[alias.casefold()]
    if required:
        raise ValueError(f"Missing required column. Tried: {', '.join(aliases)}")
    return None


def standardize_columns(raw: pd.DataFrame) -> pd.DataFrame:
    df = pd.DataFrame(index=raw.index)

    for output_column, aliases in REQUIRED_COLUMN_ALIASES.items():
        df[output_column] = raw[resolve_column(raw, aliases, required=True)]

    for output_column, aliases in OPTIONAL_COLUMN_ALIASES.items():
        source = resolve_column(raw, aliases, required=False)
        if source is not None:
            df[output_column] = raw[source]

    df["source_row_index"] = raw.index
    return df


def clean_reservoir_data(raw: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    df = standardize_columns(raw)
    raw_rows = len(df)

    for column in PROPERTY_COLUMNS:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    if "number_of_wells" in df:
        df["number_of_wells"] = pd.to_numeric(df["number_of_wells"], errors="coerce")
    if "pore_volume" in df:
        df["pore_volume"] = pd.to_numeric(df["pore_volume"], errors="coerce")

    before_missing = len(df)
    df = df.dropna(subset=PROPERTY_COLUMNS).copy()
    dropped_missing = before_missing - len(df)

    sensible_mask = (
        (df["permeability_md"] > 0.0)
        & (df["porosity"] > 0.0)
        & (df["porosity"] < 1.0)
        & (df["pressure_mpa"] > 0.0)
        & (df["temperature_c"] > 0.0)
    )
    before_sensible = len(df)
    df = df.loc[sensible_mask].copy()
    dropped_unphysical = before_sensible - len(df)

    before_property_duplicates = len(df)
    df = df.drop_duplicates(subset=PROPERTY_COLUMNS, keep="first").copy()
    dropped_property_duplicates = before_property_duplicates - len(df)

    dropped_name_duplicates = 0
    if "field_name" in df:
        normalized_name = df["field_name"].astype(str).str.strip().str.casefold()
        before_name_duplicates = len(df)
        df = df.loc[~normalized_name.duplicated()].copy()
        dropped_name_duplicates = before_name_duplicates - len(df)

    df = df.reset_index(drop=True)
    summary = {
        "raw_rows": raw_rows,
        "dropped_missing_required": dropped_missing,
        "dropped_unphysical": dropped_unphysical,
        "dropped_duplicate_properties": dropped_property_duplicates,
        "dropped_duplicate_names": dropped_name_duplicates,
        "clean_rows": len(df),
    }
    return df, summary


def choose_representative_count(n_rows: int, n_cases: int, requested: int | None) -> int:
    if requested is not None:
        if requested <= 0:
            raise ValueError("--representatives must be positive")
        return min(requested, n_rows, n_cases)

    # A compact but still diverse default: roughly 3*sqrt(n), clipped to the
    # requested 20-50 practical range. For ~150 reservoirs this gives ~37 reps.
    auto_count = int(round(3.0 * np.sqrt(n_rows)))
    auto_count = max(20, min(50, auto_count))
    return min(auto_count, n_rows, n_cases)


def add_sampling_features(df: pd.DataFrame, n_strata: int) -> pd.DataFrame:
    df = df.copy()
    df["log10_permeability"] = np.log10(df["permeability_md"])
    requested_strata = min(n_strata, df["log10_permeability"].nunique())
    df["permeability_stratum"] = pd.qcut(
        df["log10_permeability"],
        q=requested_strata,
        labels=False,
        duplicates="drop",
    ).astype(int)
    return df


def feature_matrix(df: pd.DataFrame) -> tuple[np.ndarray, StandardScaler]:
    features = df[["log10_permeability", "porosity", "pressure_mpa", "temperature_c"]].to_numpy()
    scaler = StandardScaler()
    return scaler.fit_transform(features), scaler


def density_balancing_weights(df: pd.DataFrame, density_power: float) -> np.ndarray:
    # KMeans normally follows the densest part of the data. These weights give
    # sparse permeability strata more pull during representative selection.
    # density_power=0 leaves row weights unchanged; density_power=1 gives each
    # permeability stratum about equal total influence.
    counts = df["permeability_stratum"].value_counts().to_dict()
    weights = np.array([counts[stratum] ** (-density_power) for stratum in df["permeability_stratum"]])
    return weights / weights.mean()


def choose_nearest_unused(
    distances: np.ndarray,
    candidate_indices: np.ndarray,
    used_indices: set[int],
) -> int | None:
    for local_idx in np.argsort(distances):
        candidate = int(candidate_indices[local_idx])
        if candidate not in used_indices:
            return candidate
    return None


def fill_by_farthest_point(
    standardized_features: np.ndarray,
    selected_indices: list[int],
    target_count: int,
) -> list[int]:
    selected = list(selected_indices)
    used = set(selected)

    while len(selected) < target_count:
        if selected:
            distances = cdist(standardized_features, standardized_features[selected]).min(axis=1)
        else:
            center = standardized_features.mean(axis=0, keepdims=True)
            distances = cdist(standardized_features, center).ravel()

        distances[list(used)] = -np.inf
        next_idx = int(np.argmax(distances))
        if not np.isfinite(distances[next_idx]):
            break
        selected.append(next_idx)
        used.add(next_idx)

    return selected


def select_representatives(
    df: pd.DataFrame,
    n_representatives: int,
    random_seed: int,
    density_power: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    standardized_features, _ = feature_matrix(df)
    sample_weights = density_balancing_weights(df, density_power)

    kmeans = KMeans(n_clusters=n_representatives, n_init=50, random_state=random_seed)
    cluster_labels = kmeans.fit_predict(standardized_features, sample_weight=sample_weights)

    selected_indices: list[int] = []
    used_indices: set[int] = set()

    for cluster_id in range(n_representatives):
        members = np.flatnonzero(cluster_labels == cluster_id)
        if len(members) == 0:
            continue
        distances = cdist(standardized_features[members], kmeans.cluster_centers_[cluster_id].reshape(1, -1)).ravel()
        selected = choose_nearest_unused(distances, members, used_indices)
        if selected is not None:
            selected_indices.append(selected)
            used_indices.add(selected)

    selected_indices = fill_by_farthest_point(
        standardized_features,
        selected_indices,
        n_representatives,
    )

    representatives = df.iloc[selected_indices].copy().reset_index(drop=True)
    representatives["representative_id"] = np.arange(len(representatives))
    representatives["reservoir_cluster"] = cluster_labels[selected_indices]

    cluster_info = (
        pd.DataFrame({"reservoir_cluster": cluster_labels, "row_count": 1})
        .groupby("reservoir_cluster", as_index=False)["row_count"]
        .sum()
    )

    cluster_to_size = cluster_info.set_index("reservoir_cluster")["row_count"].to_dict()
    representatives["cluster_size"] = representatives["reservoir_cluster"].map(cluster_to_size).astype(int)
    return representatives, cluster_info


def allocate_cases(
    representatives: pd.DataFrame,
    n_cases: int,
    density_power: float,
) -> np.ndarray:
    n_representatives = len(representatives)
    if n_cases < n_representatives:
        raise ValueError("n_cases must be at least the number of representatives")

    # Every representative tuple gets a guaranteed minimum. The remainder is
    # assigned by sqrt-like density weighting: dense regimes receive more cases,
    # but much less than they would under raw-frequency sampling.
    min_per_rep = max(2, min(5, n_cases // (4 * n_representatives)))
    base_total = min_per_rep * n_representatives
    if base_total > n_cases:
        min_per_rep = 1
        base_total = n_representatives

    counts = np.full(n_representatives, min_per_rep, dtype=int)
    remaining = n_cases - base_total
    if remaining <= 0:
        counts[: n_cases - counts.sum()] += 1
        return counts

    weights = representatives["cluster_size"].to_numpy(dtype=float) ** density_power
    weights = weights / weights.sum()
    fractional = remaining * weights
    extra = np.floor(fractional).astype(int)
    counts += extra

    shortfall = n_cases - counts.sum()
    if shortfall > 0:
        order = np.argsort(-(fractional - extra))
        counts[order[:shortfall]] += 1

    return counts


def latin_hypercube_operational_samples(config: Config, rng: np.random.Generator) -> pd.DataFrame:
    sampler = qmc.LatinHypercube(d=3, seed=config.random_seed)
    unit = sampler.random(config.n_cases)

    # Flow rate spans more than two orders of magnitude, so sample it uniformly
    # in log-space for better coverage of both small and large rates.
    log_flow_min = np.log10(config.flow_min)
    log_flow_max = np.log10(config.flow_max)
    flow_rate = 10 ** (log_flow_min + unit[:, 0] * (log_flow_max - log_flow_min))
    cycle_length = config.cycle_min + unit[:, 1] * (config.cycle_max - config.cycle_min)
    cg_ratio = config.cg_ratio_min + unit[:, 2] * (config.cg_ratio_max - config.cg_ratio_min)

    cg_types = np.resize(np.array(CG_TYPES), config.n_cases)
    rng.shuffle(cg_types)

    return pd.DataFrame(
        {
            "flow_rate_sm3_day": flow_rate,
            "cycle_length_days": cycle_length,
            "cg_ratio": cg_ratio,
            "cg_type": cg_types,
        }
    )


def build_simulation_pool(
    representatives: pd.DataFrame,
    allocation: np.ndarray,
    config: Config,
) -> pd.DataFrame:
    rng = np.random.default_rng(config.random_seed)
    rep_indices = np.repeat(np.arange(len(representatives)), allocation)
    rng.shuffle(rep_indices)

    operational = latin_hypercube_operational_samples(config, rng)
    assigned = representatives.iloc[rep_indices].reset_index(drop=True)

    pool = pd.concat(
        [
            operational,
            assigned[
                [
                    "permeability_md",
                    "porosity",
                    "pressure_mpa",
                    "temperature_c",
                    "permeability_stratum",
                    "reservoir_cluster",
                    "representative_id",
                    "source_row_index",
                ]
            ].rename(columns={"permeability_stratum": "reservoir_stratum"}),
        ],
        axis=1,
    )

    if "field_name" in assigned:
        pool["source_field_name"] = assigned["field_name"].to_numpy()

    return pool


def empirical_cdf(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    values = np.sort(values)
    probabilities = np.arange(1, len(values) + 1) / len(values)
    return values, probabilities


def diagnostic_histograms(real: pd.DataFrame, pool: pd.DataFrame, plot_dir: Path) -> None:
    columns = [
        ("permeability_md", "Permeability [mD]", True),
        ("porosity", "Porosity [-]", False),
        ("pressure_mpa", "Pore Pressure [MPa]", False),
        ("temperature_c", "Formation Temp [C]", False),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    for ax, (column, label, log_x) in zip(axes.ravel(), columns):
        if log_x:
            low = min(real[column].min(), pool[column].min())
            high = max(real[column].max(), pool[column].max())
            bins = np.logspace(np.log10(low), np.log10(high), 35)
            ax.set_xscale("log")
        else:
            bins = 25

        ax.hist(real[column], bins=bins, density=True, alpha=0.45, label="real")
        ax.hist(pool[column], bins=bins, density=True, alpha=0.45, label="simulation pool")
        ax.set_title(label)
        ax.set_ylabel("Density")
        ax.grid(True, alpha=0.25)

    axes[0, 0].legend()
    fig.suptitle("Reservoir Property Distribution Comparison")
    fig.tight_layout()
    fig.savefig(plot_dir / "diagnostic_histograms.png", dpi=200)
    plt.close(fig)


def diagnostic_ecdf(real: pd.DataFrame, pool: pd.DataFrame, plot_dir: Path) -> None:
    columns = [
        ("permeability_md", "Permeability [mD]", True),
        ("porosity", "Porosity [-]", False),
        ("pressure_mpa", "Pore Pressure [MPa]", False),
        ("temperature_c", "Formation Temp [C]", False),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    for ax, (column, label, log_x) in zip(axes.ravel(), columns):
        real_x, real_y = empirical_cdf(real[column].to_numpy())
        pool_x, pool_y = empirical_cdf(pool[column].to_numpy())
        ax.plot(real_x, real_y, label="real")
        ax.plot(pool_x, pool_y, label="simulation pool")
        if log_x:
            ax.set_xscale("log")
        ax.set_title(label)
        ax.set_ylabel("Empirical CDF")
        ax.grid(True, alpha=0.25)

    axes[0, 0].legend()
    fig.suptitle("Empirical CDF Comparison")
    fig.tight_layout()
    fig.savefig(plot_dir / "diagnostic_ecdf.png", dpi=200)
    plt.close(fig)


def diagnostic_scatter(
    real: pd.DataFrame,
    representatives: pd.DataFrame,
    allocation: np.ndarray,
    x_column: str,
    y_column: str,
    x_label: str,
    y_label: str,
    title: str,
    output_path: Path,
    y_log: bool = False,
) -> None:
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(real[x_column], real[y_column], s=24, c="0.72", alpha=0.55, label="real reservoirs")
    sizes = 25 + 7 * allocation
    scatter = ax.scatter(
        representatives[x_column],
        representatives[y_column],
        s=sizes,
        c=representatives["permeability_stratum"],
        cmap="viridis",
        edgecolor="black",
        linewidth=0.6,
        label="selected representatives",
    )
    if y_log:
        ax.set_yscale("log")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label("Permeability stratum")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def diagnostic_stratum_coverage(real: pd.DataFrame, pool: pd.DataFrame, plot_dir: Path) -> None:
    real_counts = real["permeability_stratum"].value_counts().sort_index()
    pool_counts = pool["reservoir_stratum"].value_counts().sort_index()
    coverage = pd.DataFrame({"real_cleaned": real_counts, "simulation_pool": pool_counts}).fillna(0).astype(int)
    coverage["real_fraction"] = coverage["real_cleaned"] / coverage["real_cleaned"].sum()
    coverage["pool_fraction"] = coverage["simulation_pool"] / coverage["simulation_pool"].sum()
    coverage.to_csv(plot_dir / "stratum_coverage.csv", index_label="permeability_stratum")

    ax = coverage[["real_fraction", "pool_fraction"]].plot(kind="bar", figsize=(8, 5))
    ax.set_xlabel("Permeability stratum")
    ax.set_ylabel("Fraction")
    ax.set_title("Permeability Stratum Coverage")
    ax.grid(axis="y", alpha=0.25)
    fig = ax.get_figure()
    fig.tight_layout()
    fig.savefig(plot_dir / "diagnostic_stratum_coverage.png", dpi=200)
    plt.close(fig)


def make_diagnostic_plots(
    real: pd.DataFrame,
    representatives: pd.DataFrame,
    allocation: np.ndarray,
    pool: pd.DataFrame,
    plot_dir: Path,
) -> None:
    plot_dir.mkdir(parents=True, exist_ok=True)
    diagnostic_histograms(real, pool, plot_dir)
    diagnostic_ecdf(real, pool, plot_dir)
    diagnostic_scatter(
        real,
        representatives,
        allocation,
        x_column="porosity",
        y_column="permeability_md",
        x_label="Porosity [-]",
        y_label="Permeability [mD]",
        title="Permeability vs Porosity",
        output_path=plot_dir / "diagnostic_permeability_vs_porosity.png",
        y_log=True,
    )
    diagnostic_scatter(
        real,
        representatives,
        allocation,
        x_column="temperature_c",
        y_column="pressure_mpa",
        x_label="Formation Temp [C]",
        y_label="Pore Pressure [MPa]",
        title="Pore Pressure vs Formation Temperature",
        output_path=plot_dir / "diagnostic_pressure_vs_temperature.png",
    )
    diagnostic_stratum_coverage(real, pool, plot_dir)


def print_summary(
    cleaning_summary: dict[str, int],
    representatives: pd.DataFrame,
    allocation: np.ndarray,
    pool: pd.DataFrame,
    cluster_info: pd.DataFrame,
    config: Config,
) -> None:
    print("\n=== Simulation Pool Summary ===")
    print(f"Input CSV: {config.input_csv}")
    print(f"Raw reservoir rows: {cleaning_summary['raw_rows']}")
    print(f"Rows after cleaning: {cleaning_summary['clean_rows']}")
    print(f"  dropped missing required values: {cleaning_summary['dropped_missing_required']}")
    print(f"  dropped physically invalid rows: {cleaning_summary['dropped_unphysical']}")
    print(f"  dropped duplicate property rows: {cleaning_summary['dropped_duplicate_properties']}")
    print(f"  dropped duplicate reservoir names: {cleaning_summary['dropped_duplicate_names']}")
    print(f"Representative reservoir tuples: {len(representatives)}")
    print(f"Final simulation cases: {len(pool)}")
    print(
        "Case allocation per representative: "
        f"min={allocation.min()}, median={np.median(allocation):.0f}, max={allocation.max()}"
    )
    print("\nCG type counts:")
    print(pool["cg_type"].value_counts().sort_index().to_string())
    print("\nCluster size summary:")
    print(cluster_info["row_count"].describe().to_string())
    print(f"\nWrote pool CSV: {config.output_csv}")
    print(f"Wrote diagnostics to: {config.plot_dir}")


def parse_args() -> Config:
    parser = argparse.ArgumentParser(
        description="Generate a density-aware 500-case reservoir simulation input pool."
    )
    parser.add_argument("--input-csv", type=Path, default=DEFAULT_INPUT_CSV)
    parser.add_argument("--output-csv", type=Path, default=DEFAULT_OUTPUT_CSV)
    parser.add_argument("--plot-dir", type=Path, default=DEFAULT_PLOT_DIR)
    parser.add_argument("--n-cases", type=int, default=500)
    parser.add_argument(
        "--representatives",
        type=int,
        default=None,
        help="Number of real reservoir-property tuples to select. Default: automatic 20-50.",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--permeability-strata", type=int, default=5)
    parser.add_argument(
        "--density-power",
        type=float,
        default=0.5,
        help=(
            "Hybrid power used in both stages: representative selection uses inverse "
            "permeability-stratum size**power to protect sparse regimes, while case "
            "allocation uses cluster size**power to keep some density awareness. Default: 0.5."
        ),
    )
    parser.add_argument("--flow-min", type=float, default=10_000.0)
    parser.add_argument("--flow-max", type=float, default=1_500_000.0)
    parser.add_argument("--cycle-min", type=float, default=14.0)
    parser.add_argument("--cycle-max", type=float, default=365.0)
    parser.add_argument("--cg-ratio-min", type=float, default=0.0)
    parser.add_argument("--cg-ratio-max", type=float, default=4.0)
    args = parser.parse_args()

    return Config(
        input_csv=args.input_csv,
        output_csv=args.output_csv,
        plot_dir=args.plot_dir,
        n_cases=args.n_cases,
        n_representatives=args.representatives,
        random_seed=args.seed,
        n_permeability_strata=args.permeability_strata,
        density_power=args.density_power,
        flow_min=args.flow_min,
        flow_max=args.flow_max,
        cycle_min=args.cycle_min,
        cycle_max=args.cycle_max,
        cg_ratio_min=args.cg_ratio_min,
        cg_ratio_max=args.cg_ratio_max,
    )


def main() -> int:
    config = parse_args()
    raw = pd.read_csv(config.input_csv)
    clean, cleaning_summary = clean_reservoir_data(raw)
    if clean.empty:
        raise ValueError("No valid reservoir rows remain after cleaning")

    clean = add_sampling_features(clean, config.n_permeability_strata)
    n_representatives = choose_representative_count(
        n_rows=len(clean),
        n_cases=config.n_cases,
        requested=config.n_representatives,
    )

    representatives, cluster_info = select_representatives(
        clean,
        n_representatives=n_representatives,
        random_seed=config.random_seed,
        density_power=config.density_power,
    )
    allocation = allocate_cases(
        representatives,
        n_cases=config.n_cases,
        density_power=config.density_power,
    )
    pool = build_simulation_pool(representatives, allocation, config)

    config.output_csv.parent.mkdir(parents=True, exist_ok=True)
    pool.to_csv(config.output_csv, index=False)
    make_diagnostic_plots(clean, representatives, allocation, pool, config.plot_dir)
    print_summary(cleaning_summary, representatives, allocation, pool, cluster_info, config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
