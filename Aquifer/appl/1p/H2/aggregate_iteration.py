import argparse
import csv
import json
import math
import re
from pathlib import Path


DEFAULT_MAX_CYCLES = 10


def read_csv_rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def read_json(path):
    with open(path) as f:
        return json.load(f)


def write_json(path, value):
    with open(path, "w") as f:
        json.dump(value, f, indent=2)


def scalar_to_csv(value):
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float, str)):
        return value
    return json.dumps(value, sort_keys=True)


def sanitize_column_name(value):
    value = re.sub(r"[^0-9A-Za-z]+", "_", str(value)).strip("_")
    return value.lower()


def get_nested(data, *keys, default=None):
    current = data
    for key in keys:
        if not isinstance(current, dict) or key not in current:
            return default
        current = current[key]
    return current


def list_value_at(values, index):
    if not isinstance(values, list):
        return ""
    if index is None or index < 0 or index >= len(values):
        return ""
    return values[index]


def last_value(values):
    if not isinstance(values, list) or not values:
        return ""
    return values[-1]


def as_int(value):
    try:
        if value is None:
            return None
        if isinstance(value, bool):
            return None
        if isinstance(value, float) and not math.isfinite(value):
            return None
        return int(value)
    except (TypeError, ValueError):
        return None


def add_final_list_value(outputs, column, values):
    value = last_value(values)
    if value != "":
        outputs[column] = value


def add_final_component_values(outputs, prefix, component_data):
    if not isinstance(component_data, dict):
        return

    for component, values in component_data.items():
        column = f"{prefix}_{sanitize_column_name(component)}"
        add_final_list_value(outputs, column, values)


def find_case_json(summary):
    case_dir = Path(summary.get("case_dir", ""))
    if not case_dir.is_dir():
        return None

    case_name = summary.get("name", "")
    if case_name:
        preferred = case_dir / f"{case_name}.json"
        if preferred.exists():
            return preferred

    for path in sorted(case_dir.glob("*.json")):
        if path.name != "case_summary.json":
            return path

    return None


def completed_cycle_indices(cycle_numbers, summary_status):
    latest_index_by_cycle = {}
    observed_cycles = set()

    for index, cycle_number in enumerate(cycle_numbers or []):
        cycle = as_int(cycle_number)
        if cycle is None or cycle < 0:
            continue
        observed_cycles.add(cycle)
        latest_index_by_cycle[cycle] = index

    completed_cycles = {
        cycle for cycle in observed_cycles
        if cycle + 1 in observed_cycles
    }

    if summary_status == "success" and observed_cycles:
        completed_cycles.add(max(observed_cycles))

    return {
        cycle: latest_index_by_cycle[cycle]
        for cycle in completed_cycles
        if cycle in latest_index_by_cycle
    }


def extract_cycle_end_outputs(case_data, summary_status, max_cycles):
    outputs = {}

    cycle_numbers = get_nested(case_data, "recoveryFactor", "cycleNumber", default=[])
    per_cycle_rf = get_nested(case_data, "recoveryFactor", "perCycle", default=[])
    cumulative_rf = get_nested(case_data, "recoveryFactor", "cumulative", default=[])
    cycle_h2_injected = get_nested(case_data, "recoveryFactor", "cycleH2Injected", default=[])
    cycle_h2_produced = get_nested(case_data, "recoveryFactor", "cycleH2Produced", default=[])
    cumulative_h2_injected = get_nested(case_data, "recoveryFactor", "cumulativeH2Injected", default=[])
    cumulative_h2_produced = get_nested(case_data, "recoveryFactor", "cumulativeH2Produced", default=[])
    purity = get_nested(case_data, "wellPurity", "H2MoleFraction", default=[])
    water_cut = get_nested(case_data, "wellWaterCut", default=[])
    times = get_nested(case_data, "time", default=[])

    cycle_indices = completed_cycle_indices(cycle_numbers, summary_status)
    for cycle in range(max_cycles):
        index = cycle_indices.get(cycle)
        column_prefix = f"cycle{cycle + 1:02d}_end"
        outputs[f"{column_prefix}_time"] = list_value_at(times, index)
        outputs[f"{column_prefix}_rf"] = list_value_at(cumulative_rf, index)
        outputs[f"{column_prefix}_cycle_rf"] = list_value_at(per_cycle_rf, index)
        outputs[f"{column_prefix}_water_cut"] = list_value_at(water_cut, index)
        outputs[f"{column_prefix}_h2_purity"] = list_value_at(purity, index)
        outputs[f"{column_prefix}_h2_injected"] = list_value_at(cycle_h2_injected, index)
        outputs[f"{column_prefix}_h2_produced"] = list_value_at(cycle_h2_produced, index)
        outputs[f"{column_prefix}_cumulative_h2_injected"] = list_value_at(cumulative_h2_injected, index)
        outputs[f"{column_prefix}_cumulative_h2_produced"] = list_value_at(cumulative_h2_produced, index)

    completed_count = len([v for v in cycle_indices if 0 <= v < max_cycles])
    outputs["completed_working_cycles"] = completed_count
    return outputs


def extract_case_outputs(summary, max_cycles):
    outputs = {}
    case_json_path = find_case_json(summary)
    if case_json_path is None:
        outputs["case_json_file"] = ""
        return outputs

    outputs["case_json_file"] = str(case_json_path)
    case_data = read_json(case_json_path)
    summary_status = summary.get("status", "")

    outputs["simulation_run_status"] = case_data.get("runStatus", "")
    outputs["simulation_failure_reason"] = case_data.get("failureReason", "")
    outputs["simulation_failure_time"] = case_data.get("failureTime", "")
    outputs["simulation_pressure_init_ref"] = case_data.get("pressureInitRef", "")
    outputs["simulation_p_inj_monitored"] = case_data.get("pInjMonitored", "")
    outputs["simulation_p_prod_monitored"] = case_data.get("pProdMonitored", "")
    outputs["simulation_suggested_dt_at_failure"] = case_data.get("suggestedTimeStepSize", "")
    outputs["simulation_last_accepted_dt_at_failure"] = case_data.get("lastAcceptedTimeStepSize", "")
    outputs["simulation_min_dt_limit"] = case_data.get("minTimeStepSize", "")

    add_final_list_value(outputs, "final_time", get_nested(case_data, "time", default=[]))
    add_final_list_value(outputs, "final_average_reservoir_pressure", get_nested(case_data, "averageReservoirPressure", default=[]))
    add_final_list_value(outputs, "final_average_reservoir_pressure_gas_bearing", get_nested(case_data, "averageReservoirPressureGasBearing", default=[]))
    add_final_list_value(outputs, "final_mass_balance_error_absolute", get_nested(case_data, "massBalanceError", "absolute", default=[]))
    add_final_list_value(outputs, "final_mass_balance_error_relative", get_nested(case_data, "massBalanceError", "relative", default=[]))
    add_final_list_value(outputs, "final_rf", get_nested(case_data, "recoveryFactor", "cumulative", default=[]))
    add_final_list_value(outputs, "final_cycle_rf", get_nested(case_data, "recoveryFactor", "perCycle", default=[]))
    add_final_list_value(outputs, "final_cycle_number", get_nested(case_data, "recoveryFactor", "cycleNumber", default=[]))
    add_final_list_value(outputs, "final_h2_purity", get_nested(case_data, "wellPurity", "H2MoleFraction", default=[]))
    add_final_list_value(outputs, "final_water_cut", get_nested(case_data, "wellWaterCut", default=[]))
    add_final_list_value(outputs, "final_cycle_h2_injected", get_nested(case_data, "recoveryFactor", "cycleH2Injected", default=[]))
    add_final_list_value(outputs, "final_cycle_h2_produced", get_nested(case_data, "recoveryFactor", "cycleH2Produced", default=[]))
    add_final_list_value(outputs, "final_cumulative_h2_injected", get_nested(case_data, "recoveryFactor", "cumulativeH2Injected", default=[]))
    add_final_list_value(outputs, "final_cumulative_h2_produced", get_nested(case_data, "recoveryFactor", "cumulativeH2Produced", default=[]))

    add_final_component_values(outputs, "final_inventory", get_nested(case_data, "inventory", default={}))
    add_final_component_values(outputs, "final_material_balance_error", get_nested(case_data, "materialBalanceError", default={}))
    add_final_component_values(outputs, "final_injection_value", get_nested(case_data, "InjectionValues", default={}))
    add_final_component_values(outputs, "final_production_value", get_nested(case_data, "ProductionValues", default={}))
    add_final_component_values(outputs, "final_injection_value_dt", get_nested(case_data, "InjectionValues_dt", default={}))
    add_final_component_values(outputs, "final_production_value_dt", get_nested(case_data, "ProductionValues_dt", default={}))

    outputs.update(extract_cycle_end_outputs(case_data, summary_status, max_cycles))
    return outputs


def summary_columns(summaries):
    preferred = [
        "case_id",
        "iter_id",
        "chunk_id",
        "name",
        "status",
        "case_dir",
        "returncode",
        "exception",
        "start_time",
        "end_time",
    ]
    present = {
        key for summary in summaries
        for key in summary.keys()
        if key not in {"inputs", "watchdog", "watchdog_config", "command"}
    }
    return [key for key in preferred if key in present] + sorted(present - set(preferred))


def output_columns(max_cycles, rows):
    preferred = [
        "case_json_file",
        "simulation_run_status",
        "simulation_failure_reason",
        "simulation_failure_time",
        "simulation_pressure_init_ref",
        "simulation_p_inj_monitored",
        "simulation_p_prod_monitored",
        "simulation_suggested_dt_at_failure",
        "simulation_last_accepted_dt_at_failure",
        "simulation_min_dt_limit",
        "final_time",
        "final_average_reservoir_pressure",
        "final_average_reservoir_pressure_gas_bearing",
        "final_mass_balance_error_absolute",
        "final_mass_balance_error_relative",
        "final_rf",
        "final_cycle_rf",
        "final_cycle_number",
        "final_water_cut",
        "final_h2_purity",
        "final_cycle_h2_injected",
        "final_cycle_h2_produced",
        "final_cumulative_h2_injected",
        "final_cumulative_h2_produced",
        "completed_working_cycles",
    ]

    cycle_columns = []
    for cycle in range(max_cycles):
        prefix = f"cycle{cycle + 1:02d}_end"
        cycle_columns.extend([
            f"{prefix}_time",
            f"{prefix}_rf",
            f"{prefix}_cycle_rf",
            f"{prefix}_water_cut",
            f"{prefix}_h2_purity",
            f"{prefix}_h2_injected",
            f"{prefix}_h2_produced",
            f"{prefix}_cumulative_h2_injected",
            f"{prefix}_cumulative_h2_produced",
        ])

    present = {key for row in rows for key in row.keys()}
    ordered = preferred + cycle_columns
    return ordered + sorted(present - set(ordered))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--iter-id", required=True)
    parser.add_argument(
        "--max-cycles",
        type=int,
        default=DEFAULT_MAX_CYCLES,
        help="Number of working-gas cycle end columns to write.",
    )
    args = parser.parse_args()

    manifest_rows = read_csv_rows(args.manifest)
    manifest_by_case_id = {str(row["case_id"]): row for row in manifest_rows}

    iter_dir = Path("cases") / args.iter_id
    summaries = []
    for summary_path in iter_dir.rglob("case_summary.json"):
        summaries.append(read_json(summary_path))

    summaries = sorted(summaries, key=lambda x: int(x["case_id"]))
    enriched_rows = []
    for summary in summaries:
        case_id = str(summary["case_id"])
        manifest_values = manifest_by_case_id.get(case_id, {})
        row = {}
        row.update(manifest_values)
        row.update({f"summary_{key}": value for key, value in summary.items()
                    if key not in {"inputs", "watchdog", "watchdog_config", "command"}})
        row.update(extract_case_outputs(summary, args.max_cycles))
        enriched_rows.append(row)

    out_csv = Path("results") / f"{args.iter_id}_summary.csv"
    out_json = Path("results") / f"{args.iter_id}_summary.json"
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    manifest_columns = list(manifest_rows[0].keys()) if manifest_rows else []
    summary_fieldnames = [f"summary_{column}" for column in summary_columns(summaries)]
    result_fieldnames = output_columns(args.max_cycles, enriched_rows)
    fieldnames = []
    for column in manifest_columns + summary_fieldnames + result_fieldnames:
        if column not in fieldnames:
            fieldnames.append(column)

    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in enriched_rows:
            writer.writerow({column: scalar_to_csv(row.get(column, "")) for column in fieldnames})

    write_json(out_json, enriched_rows)

    expected = len(manifest_rows)
    found = len(summaries)
    print(f"Manifest rows expected: {expected}")
    print(f"Case summaries found: {found}")
    print(f"Wrote CSV summary with simulation outputs: {out_csv}")
    print(f"Wrote JSON summary with simulation outputs: {out_json}")

    if found != expected:
        missing = sorted(
            set(row["case_id"] for row in manifest_rows)
            - set(str(summary["case_id"]) for summary in summaries)
        )
        print(f"WARNING: missing case summaries for case_id values: {', '.join(missing)}")

    print(f"Aggregated {found} case summaries for {args.iter_id}")


if __name__ == "__main__":
    main()
