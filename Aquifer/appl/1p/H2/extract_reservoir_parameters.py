#!/usr/bin/env python3
"""
Extract selected reservoir parameters from Reservoir_Data CSV exports.

The script reads every CSV file in Reservoir_Data, keeps the requested columns,
sorts the output rows alphabetically by field name, and writes a compact CSV.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib.pyplot as plt


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_DATA_DIR = SCRIPT_DIR / "Reservoir_Data"
DEFAULT_OUTPUT = SCRIPT_DIR / "matched_reservoir_parameters.csv"


# Output column -> possible source column names.
# The data export uses "ml" for most-likely values. The "m1" aliases are kept
# because those names are easy to read as "m1" in screenshots and copied notes.
PARAMETER_MAP: dict[str, tuple[str, ...]] = {
    "CO2 Endpoint RelPerm": ("co2endpoint_relpermeabilityml", "co2endpoint_relpermeabilitym1"),
    "Field Name": ("description",),
    "Formation Temp [C]": ("formationtempml", "formationtempm1"),
    "Latitude": ("lat",),
    "Longitude": ("lon",),
    "Number of Wells": ("totalwells",),
    "Permeability [mD]": ("storagepermeabilityml", "storagepermeabilitym1"),
    "Pore Pressure [MPa]": ("porepressureml", "porepressurem1"),
    "Pore Volume": ("porevolume",),
    "Porosity [-]": ("porosityml", "porositym1"),
    "Salinity [ppm]": ("salinity",),
    "Swr [-]": ("irreducible_water_saturationml", "irreducible_water_saturationm1"),
    "Temp Gradient": ("temperaturegradient",),
}


OUTPUT_COLUMNS = ["Field Name"] + sorted(
    column for column in PARAMETER_MAP if column != "Field Name"
)

REQUIRED_COLUMNS = (
    "Number of Wells",
    "Permeability [mD]",
    "Pore Pressure [MPa]",
    "Pore Volume",
    "Porosity [-]",
    "Formation Temp [C]",
)

PLOT_FILES = {
    "permeability_porosity": "permeability_vs_porosity.png",
    "pressure_temperature": "pressure_vs_temperature.png",
    "density": "reservoir_value_density.png",
}

PERMEABILITY_HIST_BINS = 60


def clean_value(value: object) -> str:
    if value is None:
        return ""
    value = str(value).strip()
    return "" if value == "-" else value


def first_matching_value(row: dict[str, str], source_columns: Iterable[str]) -> str:
    for column in source_columns:
        value = clean_value(row.get(column))
        if value:
            return value
    return ""


def to_float(value: str) -> float | None:
    value = clean_value(value)
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def has_required_values(row: dict[str, str]) -> bool:
    return all(to_float(row.get(column, "")) is not None for column in REQUIRED_COLUMNS)


def normalized_value(value: str) -> object:
    value = clean_value(value)
    numeric_value = to_float(value)
    if numeric_value is not None:
        return numeric_value
    return " ".join(value.casefold().split())


def duplicate_name_key(row: dict[str, str]) -> object:
    return normalized_value(row.get("Field Name", ""))


def duplicate_property_key(row: dict[str, str]) -> tuple[object, ...]:
    key = []
    for column in OUTPUT_COLUMNS:
        if column == "Field Name":
            continue
        value = clean_value(row.get(column, ""))
        key.append(normalized_value(value))
    return tuple(key)


def remove_duplicate_rows(rows: Sequence[dict[str, str]]) -> tuple[list[dict[str, str]], int]:
    unique_rows: list[dict[str, str]] = []
    seen_names: set[object] = set()
    seen_properties: set[tuple[object, ...]] = set()

    for row in rows:
        name_key = duplicate_name_key(row)
        property_key = duplicate_property_key(row)
        if name_key in seen_names or property_key in seen_properties:
            continue
        seen_names.add(name_key)
        seen_properties.add(property_key)
        unique_rows.append(row)

    return unique_rows, len(rows) - len(unique_rows)


def iter_input_csv_files(data_dir: Path, output_path: Path) -> list[Path]:
    output_path = output_path.resolve()
    return [
        path
        for path in sorted(data_dir.glob("*.csv"))
        if path.resolve() != output_path
    ]


def extract_rows(data_dir: Path, output_path: Path) -> tuple[list[dict[str, str]], list[str], int, int]:
    rows: list[dict[str, str]] = []
    missing_columns: set[str] = set()
    skipped_incomplete = 0

    for csv_path in iter_input_csv_files(data_dir, output_path):
        with csv_path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            fieldnames = set(reader.fieldnames or [])

            for output_column, source_columns in PARAMETER_MAP.items():
                if not any(source in fieldnames for source in source_columns):
                    missing_columns.add(
                        f"{csv_path.name}: {output_column} "
                        f"({', '.join(source_columns)})"
                    )

            for input_row in reader:
                if not any(clean_value(value) for value in input_row.values()):
                    continue

                output_row = {
                    output_column: first_matching_value(input_row, source_columns)
                    for output_column, source_columns in PARAMETER_MAP.items()
                }
                if not has_required_values(output_row):
                    skipped_incomplete += 1
                    continue
                rows.append(output_row)

    rows, skipped_duplicates = remove_duplicate_rows(rows)
    rows.sort(key=lambda row: (row["Field Name"].casefold(), row["Latitude"], row["Longitude"]))
    return rows, sorted(missing_columns), skipped_incomplete, skipped_duplicates


def write_output(rows: list[dict[str, str]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def values_for_plot(rows: Sequence[dict[str, str]], x_column: str, y_column: str) -> tuple[list[float], list[float]]:
    x_values: list[float] = []
    y_values: list[float] = []

    for row in rows:
        x_value = to_float(row.get(x_column, ""))
        y_value = to_float(row.get(y_column, ""))
        if x_value is None or y_value is None:
            continue
        x_values.append(x_value)
        y_values.append(y_value)

    return x_values, y_values


def log_bins(values: Sequence[float], number_of_bins: int) -> list[float]:
    positive_values = [value for value in values if value > 0.0]
    low = min(positive_values)
    high = max(positive_values)
    log_low = math.log10(low)
    log_high = math.log10(high)
    step = (log_high - log_low) / number_of_bins
    return [10 ** (log_low + i * step) for i in range(number_of_bins + 1)]


def plot_scatter(
    rows: Sequence[dict[str, str]],
    x_column: str,
    y_column: str,
    title: str,
    output_path: Path,
    y_log_scale: bool = False,
) -> None:
    x_values, y_values = values_for_plot(rows, x_column, y_column)
    if not x_values:
        raise RuntimeError(f"No plottable values for {y_column} vs {x_column}")

    plt.figure(figsize=(7, 5))
    plt.scatter(x_values, y_values, alpha=0.75)
    if y_log_scale:
        plt.yscale("log")
    plt.xlabel(x_column)
    plt.ylabel(y_column)
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def plot_density(rows: Sequence[dict[str, str]], output_path: Path) -> None:
    columns = [
        "Number of Wells",
        "Permeability [mD]",
        "Porosity [-]",
        "Pore Pressure [MPa]",
        "Formation Temp [C]",
        "Pore Volume",
    ]

    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    for ax, column in zip(axes.ravel(), columns):
        values = [to_float(row[column]) for row in rows]
        values = [value for value in values if value is not None]
        if column == "Permeability [mD]":
            values = [value for value in values if value > 0.0]
            ax.hist(values, bins=log_bins(values, PERMEABILITY_HIST_BINS), edgecolor="black", alpha=0.75)
            ax.set_xscale("log")
        elif column == "Pore Volume":
            values = [value for value in values if value > 0.0]
            ax.hist(values, bins=30, edgecolor="black", alpha=0.75)
            ax.set_xscale("log")
        else:
            ax.hist(values, bins=20, edgecolor="black", alpha=0.75)
        ax.set_title(column)
        ax.set_ylabel("Count")
        ax.grid(True, alpha=0.25)

    fig.suptitle("Reservoir Value Density", fontsize=14)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def write_plots(rows: Sequence[dict[str, str]], plot_dir: Path) -> list[Path]:
    plot_paths = [
        plot_dir / PLOT_FILES["permeability_porosity"],
        plot_dir / PLOT_FILES["pressure_temperature"],
        plot_dir / PLOT_FILES["density"],
    ]
    plot_scatter(
        rows,
        x_column="Porosity [-]",
        y_column="Permeability [mD]",
        title="Permeability vs Porosity",
        output_path=plot_paths[0],
        y_log_scale=True,
    )
    plot_scatter(
        rows,
        x_column="Formation Temp [C]",
        y_column="Pore Pressure [MPa]",
        title="Pore Pressure vs Formation Temperature",
        output_path=plot_paths[1],
    )
    plot_density(rows, plot_paths[2])
    return plot_paths


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract selected reservoir parameters into one alphabetically sorted CSV."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=DEFAULT_DATA_DIR,
        help=f"Folder containing reservoir CSV exports. Default: {DEFAULT_DATA_DIR}",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output CSV path. Default: {DEFAULT_OUTPUT}",
    )
    parser.add_argument(
        "--plot-dir",
        type=Path,
        default=None,
        help="Folder for generated PNG plots. Default: output CSV folder.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    data_dir = args.data_dir.resolve()
    output_path = args.output.resolve()
    plot_dir = args.plot_dir.resolve() if args.plot_dir is not None else output_path.parent

    if not data_dir.is_dir():
        print(f"Reservoir data folder not found: {data_dir}", file=sys.stderr)
        return 1

    rows, missing_columns, skipped_incomplete, skipped_duplicates = extract_rows(data_dir, output_path)
    if not rows:
        print(f"No reservoir rows found in {data_dir}", file=sys.stderr)
        return 1

    write_output(rows, output_path)
    plot_paths = write_plots(rows, plot_dir)

    print(f"Wrote {len(rows)} rows to {output_path}")
    print(f"Skipped {skipped_incomplete} incomplete row(s)")
    print(f"Skipped {skipped_duplicates} duplicate row(s)")
    for plot_path in plot_paths:
        print(f"Wrote plot to {plot_path.resolve()}")
    if missing_columns:
        print("Warning: some expected source columns were missing:", file=sys.stderr)
        for item in missing_columns:
            print(f"  - {item}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
