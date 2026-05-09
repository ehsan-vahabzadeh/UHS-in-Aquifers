import argparse
import csv
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path


# -----------------------------
# User settings
# -----------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
EXECUTABLE_NAME = "appl_1p2pnc_box_H2"


def default_build_dir() -> Path:
    # If this controller has been copied by CMake into build-cmake/appl/1p/H2,
    # use that folder directly. Otherwise fall back to the matching build folder
    # for the source tree.
    if (SCRIPT_DIR / EXECUTABLE_NAME).exists():
        return SCRIPT_DIR
    return SCRIPT_DIR.parents[2] / "build-cmake" / "appl" / "1p" / "H2"


BUILD_DIR = default_build_dir()
BUILD_EXECUTABLE = BUILD_DIR / EXECUTABLE_NAME

EXECUTABLE = str(BUILD_EXECUTABLE)
MPI_RUNNER = "/opt/software/RI/apps/OpenMPI/4.1.8-GCC-13.2.0/bin/mpirun"
POOL_CSV = SCRIPT_DIR / "simulation_input_pool_500.csv"

RUN_ROOT = BUILD_DIR
MANIFEST_DIR = RUN_ROOT / "manifests"
CASES_DIR = RUN_ROOT / "cases"
RESULTS_DIR = RUN_ROOT / "results"
LOGS_DIR = RUN_ROOT / "logs"

BATCH_ID = "pool_151_350"
TOTAL_SIMULATIONS = 200
POOL_START_ROW = 151
POOL_ROWS_TO_READ = 350
N_JOBS = 5
CASES_PER_JOB = 40

# The pool stores operational flow in standard m3/day. The DuMuX boundary
# condition currently expects mol/(m2 s). This factor matches the conversion
# used by the older runscript cases in this folder.
FLOW_RATE_SM3_DAY_TO_MOL_M2_S = 3.77625e-5
PRESSURE_MPA_TO_PA = 1.0e6
CELSIUS_TO_KELVIN = 273.15
MAX_TIME_STEP_SIZE_SECONDS = 43200.0


@dataclass(frozen=True)
class ControllerConfig:
    pool_csv: Path
    batch_id: str
    total_simulations: int
    pool_start_row: int
    pool_rows_to_read: int
    n_jobs: int
    cases_per_job: int
    executable: str
    mpi_runner: str
    submit: bool
    wait: bool
    candidates_csv: Path | None = None


def ensure_dirs():
    for directory in [RUN_ROOT, MANIFEST_DIR, CASES_DIR, RESULTS_DIR, LOGS_DIR]:
        directory.mkdir(parents=True, exist_ok=True)


def executable_path(executable: str) -> Path:
    path = Path(executable)
    if path.is_absolute():
        return path

    candidates = [
        SCRIPT_DIR / path,
        BUILD_DIR / path.name,
        RUN_ROOT / path,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


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
    if config.candidates_csv is not None:
        with open(config.candidates_csv, newline="") as f:
            candidate_rows = list(csv.DictReader(f))
        if not candidate_rows:
            raise ValueError(f"Candidates CSV is empty: {config.candidates_csv}")
        # Each row is a physical case from the pool. We trust pool_row_index
        # if present; otherwise fall back to the row order in the file.
        indices = []
        for local_idx, row in enumerate(candidate_rows):
            raw = row.get("pool_row_index", "")
            indices.append(int(raw) if raw not in ("", None) else local_idx)
    else:
        with open(config.pool_csv, newline="") as f:
            pool_rows = list(csv.DictReader(f))
        start_index = config.pool_start_row - 1
        end_index = start_index + config.total_simulations
        candidate_rows = pool_rows[start_index:end_index]
        if len(candidate_rows) < config.total_simulations:
            raise ValueError(
                f"Need {config.total_simulations} simulations, but only found "
                f"{len(candidate_rows)} pool rows from row {config.pool_start_row} "
                f"through row {config.pool_start_row + config.total_simulations - 1}."
            )
        indices = [start_index + i for i in range(len(candidate_rows))]

    cases = []
    for local_idx, row in enumerate(candidate_rows):
        pool_row_index = indices[local_idx]
        flow_rate_sm3_day = require_float(row, "flow_rate_sm3_day", pool_row_index)
        pressure_mpa = require_float(row, "pressure_mpa", pool_row_index)
        temperature_c = require_float(row, "temperature_c", pool_row_index)

        case = {
            "case_id": pool_row_index,
            "pool_row_index": pool_row_index,
            "flow_rate_sm3_day": flow_rate_sm3_day,
            "flow_rate_mol_m2_s": flow_rate_sm3_day * FLOW_RATE_SM3_DAY_TO_MOL_M2_S,
            "cycle_length_days": require_float(row, "cycle_length_days", pool_row_index),
            "cg_ratio": require_float(row, "cg_ratio", pool_row_index),
            "cg_type": require_text(row, "cg_type", pool_row_index),
            "permeability_md": require_float(row, "permeability_md", pool_row_index),
            "porosity": require_float(row, "porosity", pool_row_index),
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
    max_dt_seconds = min((tend_days * 86400.0) / 2200.0, MAX_TIME_STEP_SIZE_SECONDS)
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


def run_sbatch(command: list[str]) -> str:
    result = subprocess.run(command, cwd=RUN_ROOT, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            "sbatch failed before the job was accepted.\n"
            f"Command: {' '.join(command)}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )
    return result.stdout.strip()


def submit_batch(config: ControllerConfig, manifest_path: Path) -> tuple[str, str]:
    array_range = f"0-{config.n_jobs - 1}"
    array_stdout = LOGS_DIR / f"{config.batch_id}_chunk_%A_%a.out"
    array_stderr = LOGS_DIR / f"{config.batch_id}_chunk_%A_%a.err"
    agg_stdout = LOGS_DIR / f"{config.batch_id}_aggregate_%j.out"
    agg_stderr = LOGS_DIR / f"{config.batch_id}_aggregate_%j.err"
    exported_paths = (
        f"ALL,DUMUX_RUN_ROOT={RUN_ROOT},"
        f"DUMUX_SCRIPT_DIR={SCRIPT_DIR},"
        f"DUMUX_MPI_RUNNER={config.mpi_runner}"
    )

    array_cmd = [
        "sbatch",
        "--parsable",
        f"--chdir={RUN_ROOT}",
        f"--export={exported_paths}",
        f"--output={array_stdout}",
        f"--error={array_stderr}",
        f"--array={array_range}",
        str(SCRIPT_DIR / "run_chunk_array.sh"),
        str(manifest_path.resolve()),
        config.batch_id,
        str(executable_path(config.executable).resolve()),
    ]
    array_jobid = run_sbatch(array_cmd)

    agg_cmd = [
        "sbatch",
        "--parsable",
        f"--chdir={RUN_ROOT}",
        f"--export={exported_paths}",
        f"--output={agg_stdout}",
        f"--error={agg_stderr}",
        f"--dependency=afterany:{array_jobid}",
        str(SCRIPT_DIR / "aggregate_iteration.sh"),
        str(manifest_path.resolve()),
        config.batch_id,
    ]
    agg_jobid = run_sbatch(agg_cmd)

    print(f"[{config.batch_id}] array job id = {array_jobid} for tasks {array_range}")
    print(f"[{config.batch_id}] agg   job id = {agg_jobid}")
    print(f"[{config.batch_id}] chunk stdout pattern: {array_stdout}")
    print(f"[{config.batch_id}] chunk stderr pattern: {array_stderr}")
    print(f"[{config.batch_id}] aggregate stdout pattern: {agg_stdout}")
    print(f"[{config.batch_id}] aggregate stderr pattern: {agg_stderr}")
    return array_jobid, agg_jobid


def wait_for_job(jobid: str, poll_seconds: int = 20):
    while True:
        sq = subprocess.run(
            ["squeue", "-j", jobid, "-h"],
            cwd=RUN_ROOT,
            capture_output=True,
            text=True,
        )
        if sq.returncode == 0 and sq.stdout.strip() == "":
            break
        time.sleep(poll_seconds)


def parse_args() -> ControllerConfig:
    parser = argparse.ArgumentParser(
        description="Submit a selected range of simulation-pool cases as SLURM array jobs."
    )
    parser.add_argument("--pool-csv", type=Path, default=POOL_CSV)
    parser.add_argument(
        "--candidates-csv",
        type=Path,
        default=None,
        help=(
            "Optional CSV produced by active_learning_driver.py. "
            "When set, the controller runs exactly the cases listed here and "
            "ignores --pool-start-row / --pool-head / --total-simulations counts."
        ),
    )
    parser.add_argument("--batch-id", default=BATCH_ID)
    parser.add_argument("--total-simulations", type=int, default=TOTAL_SIMULATIONS)
    parser.add_argument(
        "--pool-start-row",
        type=int,
        default=POOL_START_ROW,
        help="1-based first row in the simulation pool to run.",
    )
    parser.add_argument("--pool-head", type=int, default=POOL_ROWS_TO_READ)
    parser.add_argument("--jobs", type=int, default=N_JOBS)
    parser.add_argument("--cases-per-job", type=int, default=CASES_PER_JOB)
    parser.add_argument("--executable", default=EXECUTABLE)
    parser.add_argument("--mpi-runner", default=MPI_RUNNER)
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

    total_simulations = args.total_simulations
    if args.candidates_csv is not None:
        with open(args.candidates_csv, newline="") as f:
            n_rows = sum(1 for _ in csv.DictReader(f))
        total_simulations = n_rows

    return ControllerConfig(
        pool_csv=args.pool_csv,
        batch_id=args.batch_id,
        total_simulations=total_simulations,
        pool_start_row=args.pool_start_row,
        pool_rows_to_read=args.pool_head,
        n_jobs=args.jobs,
        cases_per_job=args.cases_per_job,
        executable=args.executable,
        mpi_runner=args.mpi_runner,
        submit=not args.manifest_only,
        wait=args.wait,
        candidates_csv=args.candidates_csv,
    )


def validate_config(config: ControllerConfig):
    expected_total = config.n_jobs * config.cases_per_job
    if config.total_simulations != expected_total:
        raise ValueError(
            f"total_simulations must equal jobs * cases_per_job. "
            f"Got {config.total_simulations}, expected {expected_total}."
        )
    if config.candidates_csv is not None:
        # Pool-range checks do not apply when an explicit candidates list is used.
        return
    if config.pool_start_row < 1:
        raise ValueError("pool-start-row is 1-based and must be at least 1")
    last_required_row = config.pool_start_row + config.total_simulations - 1
    if config.pool_rows_to_read < last_required_row:
        raise ValueError(
            f"pool-head must be at least {last_required_row} for "
            f"{config.total_simulations} simulations starting at pool row "
            f"{config.pool_start_row}"
        )


def preflight_check(config: ControllerConfig, manifest_path: Path):
    errors = []
    resolved_executable = executable_path(config.executable)

    if config.candidates_csv is None:
        if not config.pool_csv.exists():
            errors.append(f"Pool CSV not found: {config.pool_csv}")
    else:
        if not config.candidates_csv.exists():
            errors.append(f"Candidates CSV not found: {config.candidates_csv}")
    if not RUN_ROOT.exists():
        errors.append(f"Build/run folder not found: {RUN_ROOT}")
    if not manifest_path.exists():
        errors.append(f"Manifest not found: {manifest_path}")
    if config.submit:
        if shutil.which("sbatch") is None:
            errors.append("Could not find 'sbatch' in PATH. Run this on the cluster login node.")
        if not resolved_executable.exists():
            errors.append(
                "Executable not found. Build it first or pass --executable explicitly.\n"
                f"Resolved executable path: {resolved_executable}"
            )
        if not (SCRIPT_DIR / "run_chunk_array.sh").exists():
            errors.append(f"Missing run script: {SCRIPT_DIR / 'run_chunk_array.sh'}")
        if not (SCRIPT_DIR / "run_chunk.py").exists():
            errors.append(f"Missing run_chunk.py: {SCRIPT_DIR / 'run_chunk.py'}")
        if not (SCRIPT_DIR / "aggregate_iteration.sh").exists():
            errors.append(f"Missing aggregate script: {SCRIPT_DIR / 'aggregate_iteration.sh'}")
        if not (SCRIPT_DIR / "aggregate_iteration.py").exists():
            errors.append(f"Missing aggregate_iteration.py: {SCRIPT_DIR / 'aggregate_iteration.py'}")
        mpi_runner = Path(config.mpi_runner)
        if mpi_runner.parent != Path("."):
            if not mpi_runner.exists():
                errors.append(f"MPI runner not found: {mpi_runner}")
        elif shutil.which(config.mpi_runner) is None:
            errors.append(f"Could not find MPI runner in PATH: {config.mpi_runner}")

    if errors:
        raise RuntimeError("Preflight check failed:\n- " + "\n- ".join(errors))

    print(f"Resolved executable: {resolved_executable}")
    print(f"MPI runner: {config.mpi_runner}")
    print(f"Run root: {RUN_ROOT}")
    print(f"Logs directory: {LOGS_DIR}")
    print(f"Case directories will be under: {CASES_DIR / config.batch_id}")


def main():
    config = parse_args()
    validate_config(config)
    ensure_dirs()

    cases = read_pool_cases(config)
    manifest_path = write_batch_manifest(config, cases)
    preflight_check(config, manifest_path)

    print(f"Wrote manifest: {manifest_path}")
    if config.candidates_csv is not None:
        print(
            f"Prepared {len(cases)} simulations from candidates CSV "
            f"{config.candidates_csv.name}: {config.n_jobs} jobs x "
            f"{config.cases_per_job} simulations/job."
        )
    else:
        print(
            f"Prepared {len(cases)} simulations from pool rows "
            f"{config.pool_start_row}-{config.pool_start_row + len(cases) - 1} "
            f"of {config.pool_csv.name}: {config.n_jobs} jobs x "
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
    try:
        main()
    except Exception as exc:
        raise SystemExit(f"ERROR: {exc}")
