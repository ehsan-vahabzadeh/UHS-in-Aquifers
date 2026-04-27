import csv
import json
import math
import subprocess
import time
from pathlib import Path

# -----------------------------
# User settings
# -----------------------------
EXECUTABLE = "bin/appl_1pnc_box_CH4"   # adjust this
TEMPLATE_DIR = Path("templates")
MANIFEST_DIR = Path("manifests")
CASES_DIR = Path("cases")
RESULTS_DIR = Path("results")
LOGS_DIR = Path("logs")

N_OUTER_ITERATIONS = 2
CASES_PER_ITERATION = 6
CHUNKS_PER_ITERATION = 3
CASES_PER_CHUNK = 2   # because 6 / 3 = 2


def ensure_dirs():
    for d in [TEMPLATE_DIR, MANIFEST_DIR, CASES_DIR, RESULTS_DIR, LOGS_DIR]:
        d.mkdir(parents=True, exist_ok=True)


def make_case_name(case):
    # Folder/name slug based on selected inputs
    return (
        f"CG{case['cg_type']}"
        f"_Q{case['flow_rate']:.2f}"
        f"_CL{case['cycle_length_days']:.1f}"
        f"_CGR{case['cg_ratio']:.2f}"
        f"_K{case['perm']:.3e}"
        f"_Poro{case['poro']:.4f}"
        f"_P{case['pressure_gwc']:.3e}"
        f"_T{case['temp']:.2f}"
    ).replace("+", "").replace(" ", "")


def derive_simulation_fields(case):
    """
    Build the DuMuX inputs from the 8 selected variables.
    Assumptions for this simple prototype:
      - one development injection period
      - one operational injection + one operational extraction
      - no idle periods
      - cycle length = injection_op + extraction_op
      - split cycle length equally between injection and extraction
      - CG ratio controls development injection duration:
            InjectionDurationDev = cg_ratio * InjectionDurationOp
    """
    cl = float(case["cycle_length_days"])
    inj_op = cl / 2.0
    ext_op = cl / 2.0
    inj_dev = float(case["cg_ratio"]) * inj_op

    tend_days = inj_dev + inj_op + ext_op

    # First approximation only: not a guarantee of exactly 2200 outputs
    max_dt_seconds = (tend_days * 86400.0) / 2200.0

    return {
        "name": make_case_name(case),
        "InitialTemperature": case["temp"],
        "PressureGWC": case["pressure_gwc"],
        "ReferencePorosity": case["poro"],
        "ReferencePermeability": case["perm"],
        "CushionGasType": case["cg_type"],
        "InjectionRateDev": -abs(case["flow_rate"]),
        "InjectionRateOp": -abs(case["flow_rate"]),
        "ProductionRate": abs(case["flow_rate"]),
        "InjectionDurationDev": inj_dev,
        "InjectionDurationOp": inj_op,
        "ExtractionDurationOp": ext_op,
        "TEnd": tend_days,
        "MaxTimeStepSize": max_dt_seconds,
        "CGRatio": case["cg_ratio"],
        "CycleLengthDays": cl,
    }


def build_demo_cases():
    """
    12 imaginary cases. Replace later with pool-selected or AL-selected cases.
    Units must match what your code expects.
    """
    raw = [
        {"flow_rate": 20, "cycle_length_days": 60,  "cg_ratio": 0.5, "cg_type": "N2",  "perm": 1e-14, "poro": 0.10, "pressure_gwc": 1.8e7, "temp": 330.0},
        {"flow_rate": 25, "cycle_length_days": 90,  "cg_ratio": 1.0, "cg_type": "N2",  "perm": 2e-14, "poro": 0.12, "pressure_gwc": 2.0e7, "temp": 335.0},
        {"flow_rate": 30, "cycle_length_days": 120, "cg_ratio": 1.5, "cg_type": "CH4", "perm": 4e-14, "poro": 0.14, "pressure_gwc": 2.2e7, "temp": 340.0},
        {"flow_rate": 35, "cycle_length_days": 180, "cg_ratio": 0.8, "cg_type": "H2",  "perm": 8e-15, "poro": 0.09, "pressure_gwc": 2.5e7, "temp": 345.0},
        {"flow_rate": 40, "cycle_length_days": 240, "cg_ratio": 2.0, "cg_type": "CO2", "perm": 6e-14, "poro": 0.16, "pressure_gwc": 2.8e7, "temp": 350.0},
        {"flow_rate": 22, "cycle_length_days": 75,  "cg_ratio": 1.2, "cg_type": "N2",  "perm": 3e-14, "poro": 0.11, "pressure_gwc": 1.9e7, "temp": 332.0},

        {"flow_rate": 28, "cycle_length_days": 110, "cg_ratio": 0.6, "cg_type": "CH4", "perm": 1.5e-14, "poro": 0.13, "pressure_gwc": 2.1e7, "temp": 337.0},
        {"flow_rate": 32, "cycle_length_days": 150, "cg_ratio": 1.8, "cg_type": "H2",  "perm": 5e-14,   "poro": 0.18, "pressure_gwc": 3.0e7, "temp": 355.0},
        {"flow_rate": 18, "cycle_length_days": 45,  "cg_ratio": 0.4, "cg_type": "N2",  "perm": 7e-15,   "poro": 0.08, "pressure_gwc": 1.7e7, "temp": 328.0},
        {"flow_rate": 26, "cycle_length_days": 95,  "cg_ratio": 1.1, "cg_type": "CO2", "perm": 2.5e-14, "poro": 0.15, "pressure_gwc": 2.4e7, "temp": 342.0},
        {"flow_rate": 38, "cycle_length_days": 210, "cg_ratio": 1.6, "cg_type": "CH4", "perm": 4.5e-14, "poro": 0.17, "pressure_gwc": 2.9e7, "temp": 348.0},
        {"flow_rate": 24, "cycle_length_days": 80,  "cg_ratio": 0.9, "cg_type": "H2",  "perm": 1.2e-14, "poro": 0.105,"pressure_gwc": 2.05e7,"temp": 334.0},
    ]

    cases = []
    for i, c in enumerate(raw):
        sim = derive_simulation_fields(c)
        sim["case_id"] = i
        cases.append(sim)
    return cases


def write_iteration_manifest(iter_id, cases_for_iter):
    manifest_path = MANIFEST_DIR / f"{iter_id}.csv"

    fieldnames = [
        "case_id", "chunk_id", "name",
        "InitialTemperature", "PressureGWC",
        "ReferencePorosity", "ReferencePermeability",
        "CushionGasType",
        "InjectionRateDev", "InjectionRateOp", "ProductionRate",
        "InjectionDurationDev", "InjectionDurationOp", "ExtractionDurationOp",
        "TEnd", "MaxTimeStepSize",
        "CGRatio", "CycleLengthDays",
    ]

    rows = []
    for local_idx, case in enumerate(cases_for_iter):
        row = case.copy()
        row["chunk_id"] = local_idx // CASES_PER_CHUNK   # 0,0,1,1,2,2
        rows.append(row)

    with open(manifest_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    return manifest_path


def submit_iteration(iter_id, manifest_path):
    # Submit 3 parallel jobs (array tasks 0,1,2), each with 24 cores
    array_cmd = [
        "sbatch",
        "--parsable",
        "--array=0-2",
        "scripts/run_chunk_array.sh",
        str(manifest_path),
        iter_id,
        EXECUTABLE,
    ]
    array_jobid = subprocess.check_output(array_cmd, text=True).strip()

    # Aggregation runs after array completes
    agg_cmd = [
        "sbatch",
        "--parsable",
        f"--dependency=afterany:{array_jobid}",
        "scripts/aggregate_iteration.sh",
        str(manifest_path),
        iter_id,
    ]
    agg_jobid = subprocess.check_output(agg_cmd, text=True).strip()

    print(f"[{iter_id}] array job id = {array_jobid}")
    print(f"[{iter_id}] agg   job id = {agg_jobid}")
    return array_jobid, agg_jobid


def wait_for_job(jobid, poll_seconds=20):
    """
    Simple polling. On many systems this is fine for a prototype.
    You can replace later with a better sacct parser.
    """
    while True:
        sq = subprocess.run(
            ["squeue", "-j", jobid, "-h"],
            capture_output=True, text=True
        )
        if sq.returncode == 0 and sq.stdout.strip() == "":
            # Not in queue anymore, assume finished one way or another
            break
        time.sleep(poll_seconds)


def main():
    ensure_dirs()
    all_cases = build_demo_cases()

    for outer_iter in range(N_OUTER_ITERATIONS):
        iter_id = f"iter_{outer_iter:03d}"

        start = outer_iter * CASES_PER_ITERATION
        end = start + CASES_PER_ITERATION
        cases_for_iter = all_cases[start:end]

        manifest_path = write_iteration_manifest(iter_id, cases_for_iter)
        _, agg_jobid = submit_iteration(iter_id, manifest_path)

        # Wait until this iteration is aggregated before next outer iteration
        wait_for_job(agg_jobid)

        print(f"[{iter_id}] aggregation finished. Moving to next iteration.")

    print("All demo iterations finished.")


if __name__ == "__main__":
    main()