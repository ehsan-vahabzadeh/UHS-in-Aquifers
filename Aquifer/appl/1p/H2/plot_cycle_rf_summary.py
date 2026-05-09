#!/usr/bin/env python3
"""
Plot end-of-cycle cumulative RF values from an aggregated simulation summary.

Default usage:
    python3 plot_cycle_rf_summary.py

Explicit CSV:
    python3 plot_cycle_rf_summary.py pool_051_150_summary.csv

The script saves two plots:
    cycle_rf_summary_plot.png
    end_rf_vs_permeability.png
    input_effects_vs_end_rf.png
    failed_flow_vs_permeability.png
    water_purity_rf_summary.png
"""

import argparse
import csv
import statistics
from pathlib import Path


DEFAULT_CSV = Path(__file__).resolve().parent / "pool_051_150_summary.csv"
RF_COLUMNS = [f"cycle{i:02d}_end_rf" for i in range(1, 11)]
WATER_CUT_COLUMNS = [f"cycle{i:02d}_end_water_cut" for i in range(1, 11)]
PURITY_COLUMNS = [f"cycle{i:02d}_end_h2_purity" for i in range(1, 11)]
DEFAULT_FEASIBLE_STATUSES = ("success", "failed_numerical_stall")
PERMEABILITY_COLUMNS = ("PermeabilityMd", "permeability_md", "ReferencePermeability")
FLOW_RATE_COLUMNS = ("FlowRateSm3Day", "flow_rate_sm3_day")
INPUT_EFFECT_COLUMNS = [
    ("FlowRateSm3Day", "Flow rate [Sm3/day]", True),
    ("CycleLengthDays", "Cycle length [days]", False),
    ("CGRatio", "Cushion-gas ratio [-]", False),
    ("PermeabilityMd", "Permeability [mD]", True),
    ("Porosity", "Porosity [-]", False),
    ("PressureMPa", "Pressure [MPa]", False),
    ("TemperatureC", "Temperature [C]", False),
]
GAS_COLORS = {
    "CH4": "#1f77b4",
    "H2": "#2ca02c",
    "N2": "#9467bd",
    "CO2": "#d62728",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot cycle-end RF curves for feasible simulations."
    )
    parser.add_argument(
        "summary_csv",
        nargs="?",
        type=Path,
        default=DEFAULT_CSV,
        help=f"Aggregated summary CSV. Default: {DEFAULT_CSV}",
    )
    parser.add_argument(
        "--status",
        action="append",
        default=None,
        help=(
            "Summary status treated as feasible. Can be repeated. "
            "Default: success and failed_numerical_stall"
        ),
    )
    parser.add_argument(
        "--end-cycle",
        type=int,
        default=10,
        help="Cycle used for the permeability scatter plot. Default: 10",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for saved plots. Default: <summary_csv_folder>/cycle_rf_plots.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Also show the figures interactively after saving.",
    )
    return parser.parse_args()


def to_float(value):
    try:
        if value is None or str(value).strip() == "":
            return None
        return float(value)
    except ValueError:
        return None


def load_summary_rows(summary_csv):
    with open(summary_csv, newline="") as f:
        return list(csv.DictReader(f))


def load_feasible_rows(summary_csv, feasible_statuses):
    rows = load_summary_rows(summary_csv)
    feasible_statuses = {status.lower() for status in feasible_statuses}

    feasible = []
    for row in rows:
        status = (row.get("summary_status") or row.get("status") or "").strip().lower()
        if status not in feasible_statuses:
            continue

        rf_values = [to_float(row.get(col)) for col in RF_COLUMNS]
        if all(value is None for value in rf_values):
            continue

        feasible.append((row, rf_values))

    return feasible, len(rows)


def format_statuses(statuses):
    return ", ".join(statuses)


def get_status(row):
    return (row.get("summary_status") or row.get("status") or "").strip()


def get_gas_type(row):
    return (row.get("CushionGasType") or row.get("cg_type") or "").strip()


def get_first_float(row, columns):
    for column in columns:
        value = to_float(row.get(column))
        if value is not None:
            return value
    return None


def get_permeability(row):
    return get_first_float(row, PERMEABILITY_COLUMNS)


def get_flow_rate(row):
    return get_first_float(row, FLOW_RATE_COLUMNS)


def get_end_rf(rf_values, end_cycle):
    if 1 <= end_cycle <= len(rf_values) and rf_values[end_cycle - 1] is not None:
        return rf_values[end_cycle - 1]

    for value in reversed(rf_values):
        if value is not None:
            return value
    return None


def cycle_values(row, columns):
    return [to_float(row.get(column)) for column in columns]


def require_matplotlib():
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit(
            "matplotlib is required to create the RF plots. Install/load matplotlib "
            "in this Python environment and run the script again."
        ) from exc
    return plt


def add_cg_legend(ax, title="Cushion gas"):
    handles = []
    for gas_type, color in GAS_COLORS.items():
        handles.append(
            ax.scatter([], [], s=55, color=color, alpha=0.8, edgecolor="black", linewidth=0.35, label=gas_type)
        )
    return ax.legend(handles=handles, title=title, frameon=True)


def plot_cycle_rf(plt, feasible_rows, total_rows, summary_csv, feasible_statuses, output_dir):
    cycles = list(range(1, 11))
    fig, ax = plt.subplots(figsize=(10, 6))

    gas_types_seen = set()
    rf_matrix = []

    for row, rf_values in feasible_rows:
        gas_type = get_gas_type(row)
        color = GAS_COLORS.get(gas_type, "#7f7f7f")
        label = gas_type if gas_type and gas_type not in gas_types_seen else None
        if gas_type:
            gas_types_seen.add(gas_type)

        ax.plot(
            cycles,
            rf_values,
            color=color,
            alpha=0.35,
            linewidth=1.2,
            marker="o",
            markersize=3,
            label=label,
        )
        rf_matrix.append(rf_values)

    if rf_matrix:
        mean_rf = []
        median_rf = []
        for cycle_index in range(len(cycles)):
            values = [
                rf_values[cycle_index]
                for rf_values in rf_matrix
                if rf_values[cycle_index] is not None
            ]
            mean_rf.append(sum(values) / len(values) if values else None)
            median_rf.append(statistics.median(values) if values else None)

        ax.plot(cycles, mean_rf, color="black", linewidth=2.4, marker="o", label="Mean")
        ax.plot(cycles, median_rf, color="black", linewidth=2.0, linestyle="--", label="Median")

    ax.set_title(
        "End-of-Cycle Cumulative Recovery Factor\n"
        f"{len(feasible_rows)} usable simulations ({format_statuses(feasible_statuses)}) from {total_rows} rows"
    )
    ax.set_xlabel("Working gas cycle")
    ax.set_ylabel("Cumulative RF [-]")
    ax.set_xticks(cycles)
    ax.grid(True, alpha=0.3)
    ax.legend(title="Cushion gas / summary", frameon=True)
    fig.tight_layout()

    output_path = output_dir / "cycle_rf_summary_plot.png"
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Saved: {output_path}")


def plot_end_rf_vs_permeability(plt, feasible_rows, end_cycle, output_dir):
    fig, ax = plt.subplots(figsize=(9, 6))

    gas_types_seen = set()
    plotted = 0

    for row, rf_values in feasible_rows:
        permeability = get_permeability(row)
        end_rf = get_end_rf(rf_values, end_cycle)
        if permeability is None or end_rf is None:
            continue

        gas_type = get_gas_type(row)
        color = GAS_COLORS.get(gas_type, "#7f7f7f")
        label = gas_type if gas_type and gas_type not in gas_types_seen else None
        if gas_type:
            gas_types_seen.add(gas_type)

        ax.scatter(
            permeability,
            end_rf,
            s=55,
            color=color,
            alpha=0.78,
            edgecolor="black",
            linewidth=0.35,
            label=label,
        )
        plotted += 1

    ax.set_title(f"Cycle {end_cycle} Cumulative RF vs Reservoir Permeability")
    ax.set_xlabel("Reservoir permeability [mD]")
    ax.set_ylabel(f"Cycle {end_cycle} cumulative RF [-]")
    ax.grid(True, alpha=0.3)
    ax.legend(title="Cushion gas", frameon=True)
    fig.tight_layout()

    output_path = output_dir / "end_rf_vs_permeability.png"
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Saved: {output_path}")
    print(f"Scatter points plotted: {plotted}")


def plot_input_effects(plt, feasible_rows, end_cycle, output_dir):
    fig, axes = plt.subplots(2, 4, figsize=(17, 8), sharey=True)
    axes = list(axes.ravel())

    for ax, (column, label, use_log_x) in zip(axes, INPUT_EFFECT_COLUMNS):
        plotted_values = []
        plotted = 0

        for row, rf_values in feasible_rows:
            x_value = to_float(row.get(column))
            end_rf = get_end_rf(rf_values, end_cycle)
            if x_value is None or end_rf is None:
                continue

            gas_type = get_gas_type(row)
            ax.scatter(
                x_value,
                end_rf,
                s=38,
                color=GAS_COLORS.get(gas_type, "#7f7f7f"),
                alpha=0.72,
                edgecolor="black",
                linewidth=0.25,
            )
            plotted_values.append(x_value)
            plotted += 1

        if use_log_x and plotted_values and all(value > 0 for value in plotted_values):
            ax.set_xscale("log")

        ax.set_title(label)
        ax.set_xlabel(label)
        ax.grid(True, alpha=0.3)
        if plotted == 0:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)

    for ax in axes[::4]:
        ax.set_ylabel(f"Cycle {end_cycle} cumulative RF [-]")

    unused_axes = axes[len(INPUT_EFFECT_COLUMNS):]
    for ax in unused_axes:
        ax.axis("off")

    add_cg_legend(unused_axes[0] if unused_axes else axes[-1])
    fig.suptitle(f"Input Effects on Cycle {end_cycle} Cumulative RF", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    output_path = output_dir / "input_effects_vs_end_rf.png"
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Saved: {output_path}")


def plot_water_purity_rf_summary(plt, feasible_rows, output_dir):
    cycles = list(range(1, 11))
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    water_ax, purity_ax, compare_ax = axes

    gas_types_seen = set()
    water_matrix = []
    purity_matrix = []

    for row, rf_values in feasible_rows:
        gas_type = get_gas_type(row)
        color = GAS_COLORS.get(gas_type, "#7f7f7f")
        label = gas_type if gas_type and gas_type not in gas_types_seen else None
        if gas_type:
            gas_types_seen.add(gas_type)

        water_values = cycle_values(row, WATER_CUT_COLUMNS)
        purity_values = cycle_values(row, PURITY_COLUMNS)

        water_percent = [value * 100 if value is not None else None for value in water_values]
        purity_percent = [value * 100 if value is not None else None for value in purity_values]

        water_ax.plot(
            cycles,
            water_percent,
            color=color,
            alpha=0.32,
            linewidth=1.1,
            marker="o",
            markersize=2.8,
            label=label,
        )
        purity_ax.plot(
            cycles,
            purity_percent,
            color=color,
            alpha=0.32,
            linewidth=1.1,
            marker="o",
            markersize=2.8,
            label=label,
        )

        water_matrix.append(water_values)
        purity_matrix.append(purity_values)

        for rf_value, purity_value in zip(rf_values, purity_values):
            if rf_value is None or purity_value is None:
                continue
            compare_ax.scatter(
                rf_value,
                purity_value * 100,
                s=30,
                color=color,
                alpha=0.45,
                edgecolor="black",
                linewidth=0.18,
            )

    for ax, matrix, scale, label in [
        (water_ax, water_matrix, 100.0, "Mean water cut"),
        (purity_ax, purity_matrix, 100.0, "Mean purity"),
    ]:
        mean_values = []
        median_values = []
        for cycle_index in range(len(cycles)):
            values = [
                values_by_cycle[cycle_index]
                for values_by_cycle in matrix
                if values_by_cycle[cycle_index] is not None
            ]
            mean_values.append((sum(values) / len(values)) * scale if values else None)
            median_values.append(statistics.median(values) * scale if values else None)

        ax.plot(cycles, mean_values, color="black", linewidth=2.3, marker="o", label=label)
        ax.plot(cycles, median_values, color="black", linewidth=1.9, linestyle="--", label="Median")

    water_ax.set_title("Cycle-End Water Cut")
    water_ax.set_xlabel("Working gas cycle")
    water_ax.set_ylabel("Water cut [%]")
    water_ax.set_xticks(cycles)
    water_ax.grid(True, alpha=0.3)

    purity_ax.set_title("Cycle-End H2 Purity")
    purity_ax.set_xlabel("Working gas cycle")
    purity_ax.set_ylabel("H2 purity [%]")
    purity_ax.set_xticks(cycles)
    purity_ax.grid(True, alpha=0.3)

    compare_ax.set_title("H2 Purity vs Cumulative RF")
    compare_ax.set_xlabel("Cycle-end cumulative RF [-]")
    compare_ax.set_ylabel("Cycle-end H2 purity [%]")
    compare_ax.grid(True, alpha=0.3)

    water_ax.legend(title="Cushion gas / summary", frameon=True, fontsize=8)
    fig.suptitle("Water Cut, H2 Purity, and RF Relationship", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))

    output_path = output_dir / "water_purity_rf_summary.png"
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Saved: {output_path}")


def failure_label(row):
    reason = (row.get("simulation_failure_reason") or "").strip()
    if reason:
        return reason
    return get_status(row)


def plot_failed_flow_vs_permeability(plt, all_rows, feasible_statuses, output_dir):
    from matplotlib.lines import Line2D

    feasible_statuses = {status.lower() for status in feasible_statuses}
    failed_rows = [
        row for row in all_rows
        if get_status(row).lower().startswith("failed")
        and get_status(row).lower() not in feasible_statuses
    ]

    marker_by_reason = {
        "inj_pressure_exceeded": "o",
        "failed_numerical_stall": "X",
        "failed": "^",
    }

    fig, ax = plt.subplots(figsize=(9, 6))
    plotted = 0
    reasons_seen = set()

    for row in failed_rows:
        flow_rate = get_flow_rate(row)
        permeability = get_permeability(row)
        if flow_rate is None or permeability is None:
            continue

        gas_type = get_gas_type(row)
        reason = failure_label(row)
        marker = marker_by_reason.get(reason, "s")
        reasons_seen.add(reason)

        ax.scatter(
            flow_rate,
            permeability,
            s=70,
            marker=marker,
            color=GAS_COLORS.get(gas_type, "#7f7f7f"),
            alpha=0.78,
            edgecolor="black",
            linewidth=0.45,
        )
        plotted += 1

    ax.set_title("Failed Runs: Flow Rate vs Reservoir Permeability")
    ax.set_xlabel("Flow rate [Sm3/day]")
    ax.set_ylabel("Reservoir permeability [mD]")
    ax.grid(True, alpha=0.3)

    if plotted:
        ax.set_xscale("log")
        ax.set_yscale("log")

    cg_legend = add_cg_legend(ax)
    ax.add_artist(cg_legend)

    reason_handles = [
        Line2D(
            [0], [0],
            marker=marker_by_reason.get(reason, "s"),
            color="none",
            markerfacecolor="#bdbdbd",
            markeredgecolor="black",
            markersize=8,
            linestyle="None",
            label=reason,
        )
        for reason in sorted(reasons_seen)
    ]
    if reason_handles:
        ax.legend(handles=reason_handles, title="Failure type", frameon=True, loc="lower right")

    fig.tight_layout()

    output_path = output_dir / "failed_flow_vs_permeability.png"
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Saved: {output_path}")
    print(f"Failed-run points plotted: {plotted}")


def main():
    args = parse_args()
    feasible_statuses = args.status if args.status is not None else list(DEFAULT_FEASIBLE_STATUSES)
    all_rows = load_summary_rows(args.summary_csv)
    feasible_rows, total_rows = load_feasible_rows(args.summary_csv, feasible_statuses)
    if not feasible_rows:
        raise SystemExit(
            f"No rows with summary_status in {feasible_statuses!r} and RF values found in {args.summary_csv}"
        )

    output_dir = args.output_dir if args.output_dir is not None else args.summary_csv.resolve().parent / "cycle_rf_plots"
    output_dir.mkdir(parents=True, exist_ok=True)

    plt = require_matplotlib()
    print(f"Loaded: {args.summary_csv}")
    print(f"Feasible simulations plotted: {len(feasible_rows)} / {total_rows}")
    print(f"Feasible statuses: {format_statuses(feasible_statuses)}")

    plot_cycle_rf(plt, feasible_rows, total_rows, args.summary_csv, feasible_statuses, output_dir)
    plot_end_rf_vs_permeability(plt, feasible_rows, args.end_cycle, output_dir)
    plot_input_effects(plt, feasible_rows, args.end_cycle, output_dir)
    plot_water_purity_rf_summary(plt, feasible_rows, output_dir)
    plot_failed_flow_vs_permeability(plt, all_rows, feasible_statuses, output_dir)

    if args.show:
        plt.show()
    else:
        plt.close("all")


if __name__ == "__main__":
    main()
