import argparse
import csv
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path


# -----------------------------
# User settings
# -----------------------------
SCRIPT_DIR = Path(__file__).resolve().parent

EXECUTABLE = "appl_1p2pnc_box_H2"
POOL_CSV = SCRIPT_DIR / "simulation_input_pool_500.csv"

MANIFEST_DIR = SCRIPT_DIR / "manifests"
CASES_DIR = SCRIPT_DIR / "cases"
RESULTS_DIR = SCRIPT_DIR / "results"
LOGS_DIR = SCRIPT_DIR / "logs"

BATCH_ID = "pool_top30"
TOTAL_SIMULATIONS = 30
POOL_ROWS_TO_READ = 40
N_JOBS = 5
CASES_PER_JOB = 6

# The pool stores operational flow in standard m3/day. The DuMuX boundary
# condition currently expects mol/(m2 s). This factor matches the conversion
# used by the older runscript cases in this folder.
FLOW_RATE_SM3_DAY_TO_MOL_M2_S = 3.77625e-5
PRESSURE_MPA_TO_PA = 1.0e6
CELSIUS_TO_KELVIN = 273.15


@dataclass(frozen=True)
class ControllerConfig:
    pool_csv: Path
    batch_id: str
    total_simulations: int
    pool_rows_to_read: int
    n_jobs: int
    cases_per_job: int
    executable: str
    submit: bool
    wait: bool


def ensure_dirs():
    for directory in [MANIFEST_DIR, CASES_DIR, RESULTS_DIR, LOGS_DIR]:
        directory.mkdir(parents=True, exist_ok=True)


def executable_path(executable: str) -> Path:
    path = Path(executable)
    if path.is_absolute():
        return path
    return SCRIPT_DIR / path


def require_float(row: dict, column: str, row_number: int) -> float:
    value = row.get(column, "")
    if value == "":
        raise ValueError(f"Missing required column '{column}' in pool row {row_number}")
    return float(value)


def require_text(row: dict, column: str, row_number: int) -> str:
    value = row.get(column, "")
    if value == "":
        raise ValueError(f"Missing required column '{column}' in pool row {row_number}")
    return str(value)


def make_case_name(case: dict) -> str:
    return (
        f"case{case['case_id']:03d}"
        f"_CG{case['cg_type']}"
        f"_Q{case['flow_rate_sm3_day']:.0f}"
        f"_CL{case['cycle_length_days']:.1f}"
        f"_CGR{case['cg_ratio']:.2f}"
        f"_K{case['permeability_md']:.3g}"
        f"_Poro{case['porosity']:.3f}"
    ).replace("+", "").replace(" ", "")


def read_pool_cases(config: ControllerConfig) -> list[dict]:
    with open(config.pool_csv, newline="") as f:
        pool_rows = list(csv.DictReader(f))

    candidate_rows = pool_rows[: config.pool_rows_to_read]
    if len(candidate_rows) < config.total_simulations:
        raise ValueError(
            f"Need {config.total_simulations} simulations, but only found "
            f"{len(candidate_rows)} rows after reading the top {config.pool_rows_to_read} pool rows."
        )

    cases = []
    for case_id, row in enumerate(candidate_rows[: config.total_simulations]):
        flow_rate_sm3_day = require_float(row, "flow_rate_sm3_day", case_id)
        pressure_mpa = require_float(row, "pressure_mpa", case_id)
        temperature_c = require_float(row, "temperature_c", case_id)

        case = {
            "case_id": case_id,
            "pool_row_index": case_id,
            "flow_rate_sm3_day": flow_rate_sm3_day,
            "flow_rate_mol_m2_s": flow_rate_sm3_day * FLOW_RATE_SM3_DAY_TO_MOL_M2_S,
            "cycle_length_days": require_float(row, "cycle_length_days", case_id),
            "cg_ratio": require_float(row, "cg_ratio", case_id),
            "cg_type": require_text(row, "cg_type", case_id),
            "permeability_md": require_float(row, "permeability_md", case_id),
            "porosity": require_float(row, "porosity", case_id),
            "pressure_mpa": pressure_mpa,
            "pressure_pa": pressure_mpa * PRESSURE_MPA_TO_PA,
            "temperature_c": temperature_c,
            "temperature_k": temperature_c + CELSIUS_TO_KELVIN,
            "reservoir_stratum": row.get("reservoir_stratum", ""),
            "reservoir_cluster": row.get("reservoir_cluster", ""),
            "representative_id": row.get("representative_id", ""),
            "source_row_index": row.get("source_row_index", ""),
            "source_field_name": row.get("source_field_name", ""),
        }
        cases.append(derive_simulation_fields(case))

    return cases


def derive_simulation_fields(case: dict) -> dict:
    """
    Build the DuMuX inputs from one row of simulation_input_pool_500.csv.

    The reservoir tuple is kept together from the generated pool. Pressure is
    converted MPa -> Pa, temperature is converted C -> K, and flow rate is
    converted from standard m3/day to the mol/(m2 s) boundary rate used here.
    """
    cycle_length = float(case["cycle_length_days"])
    injection_op = cycle_length / 2.0
    extraction_op = cycle_length / 2.0
    injection_dev = float(case["cg_ratio"]) * injection_op

    # Ten operational cycles after the development injection. Time is in days.
    tend_days = injection_dev + (injection_op + extraction_op) * 10.0
    max_dt_seconds = (tend_days * 86400.0) / 2200.0
    boundary_rate = abs(float(case["flow_rate_mol_m2_s"]))

    return {
        "case_id": case["case_id"],
        "pool_row_index": case["pool_row_index"],
        "name": make_case_name(case),
        "InitialTemperature": case["temperature_k"],
        "PressureGWC": case["pressure_pa"],
        "ReferencePorosity": case["porosity"],
        "ReferencePermeability": case["permeability_md"],
        "CushionGasType": case["cg_type"],
        "InjectionRateDev": -boundary_rate,
        "InjectionRateOp": -boundary_rate,
        "ProductionRate": boundary_rate,
        "InjectionDurationDev": injection_dev,
        "InjectionDurationOp": injection_op,
        "ExtractionDurationOp": extraction_op,
        "TEnd": tend_days,
        "MaxTimeStepSize": max_dt_seconds,
        "CGRatio": case["cg_ratio"],
        "CycleLengthDays": cycle_length,
        "FlowRateSm3Day": case["flow_rate_sm3_day"],
        "FlowRateMolM2S": boundary_rate,
        "PermeabilityMd": case["permeability_md"],
        "Porosity": case["porosity"],
        "PressureMPa": case["pressure_mpa"],
        "TemperatureC": case["temperature_c"],
        "reservoir_stratum": case["reservoir_stratum"],
        "reservoir_cluster": case["reservoir_cluster"],
        "representative_id": case["representative_id"],
        "source_row_index": case["source_row_index"],
        "source_field_name": case["source_field_name"],
    }


def write_batch_manifest(config: ControllerConfig, cases: list[dict]) -> Path:
    manifest_path = MANIFEST_DIR / f"{config.batch_id}.csv"

    fieldnames = [
        "case_id",
        "chunk_id",
        "pool_row_index",
        "name",
        "InitialTemperature",
        "PressureGWC",
        "ReferencePorosity",
        "ReferencePermeability",
        "CushionGasType",
        "InjectionRateDev",
        "InjectionRateOp",
        "ProductionRate",
        "InjectionDurationDev",
        "InjectionDurationOp",
        "ExtractionDurationOp",
        "TEnd",
        "MaxTimeStepSize",
        "CGRatio",
        "CycleLengthDays",
        "FlowRateSm3Day",
        "FlowRateMolM2S",
        "PermeabilityMd",
        "Porosity",
        "PressureMPa",
        "TemperatureC",
        "reservoir_stratum",
        "reservoir_cluster",
        "representative_id",
        "source_row_index",
        "source_field_name",
    ]

    rows = []
    for local_idx, case in enumerate(cases):
        row = case.copy()
        row["chunk_id"] = local_idx // config.cases_per_job
        rows.append(row)

    with open(manifest_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    return manifest_path


def submit_batch(config: ControllerConfig, manifest_path: Path) -> tuple[str, str]:
    array_range = f"0-{config.n_jobs - 1}"
    array_cmd = [
        "sbatch",
        "--parsable",
        f"--array={array_range}",
        str(SCRIPT_DIR / "run_chunk_array.sh"),
        str(manifest_path.resolve()),
        config.batch_id,
        str(executable_path(config.executable).resolve()),
    ]
    array_jobid = subprocess.check_output(array_cmd, cwd=SCRIPT_DIR, text=True).strip()

    agg_cmd = [
        "sbatch",
        "--parsable",
        f"--dependency=afterany:{array_jobid}",
        str(SCRIPT_DIR / "aggregate_iteration.sh"),
        str(manifest_path.resolve()),
        config.batch_id,
    ]
    agg_jobid = subprocess.check_output(agg_cmd, cwd=SCRIPT_DIR, text=True).strip()

    print(f"[{config.batch_id}] array job id = {array_jobid} for tasks {array_range}")
    print(f"[{config.batch_id}] agg   job id = {agg_jobid}")
    return array_jobid, agg_jobid


def wait_for_job(jobid: str, poll_seconds: int = 20):
    while True:
        sq = subprocess.run(
            ["squeue", "-j", jobid, "-h"],
            cwd=SCRIPT_DIR,
            capture_output=True,
            text=True,
        )
        if sq.returncode == 0 and sq.stdout.strip() == "":
            break
        time.sleep(poll_seconds)


def parse_args() -> ControllerConfig:
    parser = argparse.ArgumentParser(
        description="Submit the first 30 simulation-pool cases as 5 SLURM jobs with 6 cases per job."
    )
    parser.add_argument("--pool-csv", type=Path, default=POOL_CSV)
    parser.add_argument("--batch-id", default=BATCH_ID)
    parser.add_argument("--total-simulations", type=int, default=TOTAL_SIMULATIONS)
    parser.add_argument("--pool-head", type=int, default=POOL_ROWS_TO_READ)
    parser.add_argument("--jobs", type=int, default=N_JOBS)
    parser.add_argument("--cases-per-job", type=int, default=CASES_PER_JOB)
    parser.add_argument("--executable", default=EXECUTABLE)
    parser.add_argument(
        "--manifest-only",
        action="store_true",
        help="Only write the manifest; do not call sbatch.",
    )
    parser.add_argument(
        "--wait",
        action="store_true",
        help="Wait until the aggregation job leaves the SLURM queue.",
    )
    args = parser.parse_args()

    return ControllerConfig(
        pool_csv=args.pool_csv,
        batch_id=args.batch_id,
        total_simulations=args.total_simulations,
        pool_rows_to_read=args.pool_head,
        n_jobs=args.jobs,
        cases_per_job=args.cases_per_job,
        executable=args.executable,
        submit=not args.manifest_only,
        wait=args.wait,
    )


def validate_config(config: ControllerConfig):
    expected_total = config.n_jobs * config.cases_per_job
    if config.total_simulations != expected_total:
        raise ValueError(
            f"total_simulations must equal jobs * cases_per_job. "
            f"Got {config.total_simulations}, expected {expected_total}."
        )
    if config.pool_rows_to_read < config.total_simulations:
        raise ValueError("pool-head must be at least total-simulations")


def main():
    config = parse_args()
    validate_config(config)
    ensure_dirs()

    cases = read_pool_cases(config)
    manifest_path = write_batch_manifest(config, cases)

    print(f"Wrote manifest: {manifest_path}")
    print(
        f"Prepared {len(cases)} simulations from the top {config.pool_rows_to_read} "
        f"rows of {config.pool_csv.name}: {config.n_jobs} jobs x "
        f"{config.cases_per_job} simulations/job."
    )

    if not config.submit:
        print("Manifest-only mode: not submitting SLURM jobs.")
        return

    _, agg_jobid = submit_batch(config, manifest_path)
    if config.wait:
        wait_for_job(agg_jobid)
        print(f"[{config.batch_id}] aggregation finished.")


if __name__ == "__main__":
    main()
