import argparse
import csv
import json
import shutil
import subprocess
from pathlib import Path


def read_chunk_rows(manifest_path, chunk_id):
    with open(manifest_path, newline="") as f:
        rows = list(csv.DictReader(f))
    return [r for r in rows if int(r["chunk_id"]) == chunk_id]


def write_params_file(case_dir: Path, row: dict):
    """
    Creates a case-local params.input based on your template structure.
    """
    text = f"""[TimeLoop]
DtInitial = 10
TEnd = {row["TEnd"]}
MaxTimeStepSize = {row["MaxTimeStepSize"]}

[Grid]
LowerLeft = 0.2 0
UpperRight = 1000.2 60
Cells = 200 20

[Problem]
Name = {row["name"]}
EnableGravity = true
UseNitscheTypeBc = true
DispersionMode = 1
InitialTemperature = {row["InitialTemperature"]}

[BoundaryConditions]
CyclesDev = 1
InjectionDurationDev = {row["InjectionDurationDev"]}
IdleDurationDev = 0
InjectionDurationOp = {row["InjectionDurationOp"]}
ExtractionDurationOp = {row["ExtractionDurationOp"]}
IdleDurationOp = 0.0
CushionGasType = {row["CushionGasType"]}
InjectionRateDev = {row["InjectionRateDev"]}
InjectionRateOp = {row["InjectionRateOp"]}
ProductionRate = {row["ProductionRate"]}
Well_Height = 10
HydrogenInjectionConcentration = 1

[Initialization]
PressureGWC = {row["PressureGWC"]}
DepthGWC = -10

[SpatialParams]
ReferencePorosity = {row["ReferencePorosity"]}
ReferencePermeability = {row["ReferencePermeability"]}
Material.Swr = 0.15
Material.Snr = 0.0
Material.BrooksCoreyPcEntry = 5e3
Material.BrooksCoreyLambda = 0.767133
Material.n_w = 8
Material.n_nw = 2.5

[Newton]
MaxRelativeShift = 1e-6
MaxSteps = 9
MaxTimeStepDivisions = 20
UseLineSearch = true

[Assembly]
Multithreading = false

[Vtk]
AddVelocity = true

[Safety]
EnablePressureCutoff = true
InjectionPressureMultiplier = 2
MinProductionPressure = 1e5
"""
    (case_dir / "params.input").write_text(text)


def run_one_case(row, iter_id, chunk_id, ntasks, executable):
    case_id = int(row["case_id"])
    case_dir = Path("cases") / iter_id / f"chunk_{chunk_id}" / f"case_{case_id:03d}_{row['name']}"
    case_dir.mkdir(parents=True, exist_ok=True)

    write_params_file(case_dir, row)

    cmd = [
        "srun", "-n", str(ntasks),
        str(Path(executable).resolve()),
        str((case_dir / "params.input").resolve())
    ]

    summary = {
        "case_id": case_id,
        "iter_id": iter_id,
        "chunk_id": chunk_id,
        "name": row["name"],
        "status": "failed",
        "inputs": row,
    }

    try:
        result = subprocess.run(cmd, cwd=case_dir, capture_output=True, text=True)
        (case_dir / "stdout.txt").write_text(result.stdout)
        (case_dir / "stderr.txt").write_text(result.stderr)
        summary["returncode"] = result.returncode

        if result.returncode == 0:
            # If your code writes json files, keep them local to this case folder.
            json_files = sorted(case_dir.glob("*.json"))
            summary["json_outputs"] = [p.name for p in json_files]
            summary["status"] = "success"
        else:
            summary["status"] = "failed"

    except Exception as e:
        summary["status"] = "failed"
        summary["exception"] = str(e)

    with open(case_dir / "case_summary.json", "w") as f:
        json.dump(summary, f, indent=2)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--iter-id", required=True)
    parser.add_argument("--chunk-id", type=int, required=True)
    parser.add_argument("--ntasks", type=int, required=True)
    parser.add_argument("--executable", required=True)
    args = parser.parse_args()

    rows = read_chunk_rows(args.manifest, args.chunk_id)

    # Each SLURM array task runs all manifest rows for its chunk sequentially.
    # With the default controller settings that is 6 simulations per job.
    for row in rows:
        run_one_case(
            row=row,
            iter_id=args.iter_id,
            chunk_id=args.chunk_id,
            ntasks=args.ntasks,
            executable=args.executable,
        )


if __name__ == "__main__":
    main()
