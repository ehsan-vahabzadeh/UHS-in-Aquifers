import argparse
import csv
import json
import os
import signal
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path


NUMERICAL_FAILURE_MARKERS = (
    "defect=-nan",
    "is infinite or nan",
    "matrix is singular",
    "ilu failed",
    "solverabort",
    "numericalproblem",
    "linear solver did not converge",
)

MONITORED_OUTPUT_EXTENSIONS = {
    ".json",
    ".pvd",
    ".pvtu",
    ".vtu",
    ".vtm",
    ".vts",
}


def now_iso():
    return datetime.now().isoformat(timespec="seconds")


def env_float(name, default):
    value = os.environ.get(name)
    if value is None or value == "":
        return default
    return float(value)


def read_new_text(path, offset):
    try:
        size = path.stat().st_size
        if size < offset:
            offset = 0
        with open(path, "rb") as f:
            f.seek(offset)
            text = f.read().decode(errors="replace")
        return text, size
    except OSError:
        return "", offset


def read_chunk_rows(manifest_path, chunk_id):
    with open(manifest_path, newline="") as f:
        rows = list(csv.DictReader(f))
    return [r for r in rows if int(r["chunk_id"]) == chunk_id]


def resolve_mpi_runner(mpi_runner: str) -> str:
    path = Path(mpi_runner)
    if path.parent != Path("."):
        return str(path.resolve())

    resolved = shutil.which(mpi_runner)
    return resolved if resolved is not None else mpi_runner


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
Cells = 300 30

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
InjectionPressureMultiplier = 1.5
MinProductionPressure = 1e5
EnableTimeStepCutoff = true
MinTimeStepSize = 1e-4
MinTimeStepCheckTime = 86400
MinTimeStepConsecutiveSteps = 5
"""
    (case_dir / "params.input").write_text(text)


def monitored_case_files(case_dir, stdout_path, stderr_path):
    files = []
    for path in (stdout_path, stderr_path):
        if path.exists():
            files.append(path)

    for path in case_dir.iterdir():
        if not path.is_file():
            continue
        if path.name == "case_summary.json":
            continue
        if path.suffix.lower() in MONITORED_OUTPUT_EXTENSIONS:
            files.append(path)

    return sorted(set(files))


def file_activity_signature(files):
    signature = []
    for path in files:
        try:
            stat = path.stat()
        except FileNotFoundError:
            continue
        signature.append((path.name, stat.st_size, stat.st_mtime_ns))
    return tuple(signature)


def latest_json_time(case_dir):
    best_count = 0
    best_value = None
    for path in case_dir.glob("*.json"):
        if path.name == "case_summary.json":
            continue
        try:
            with open(path) as f:
                data = json.load(f)
        except (json.JSONDecodeError, OSError):
            continue

        values = data.get("time")
        if isinstance(values, list) and values:
            count = len(values)
            value = values[-1]
            if count > best_count:
                best_count = count
                best_value = value

    return best_count, best_value


def progress_snapshot(case_dir, stdout_path, stderr_path):
    files = monitored_case_files(case_dir, stdout_path, stderr_path)
    time_count, time_value = latest_json_time(case_dir)
    return {
        "file_signature": file_activity_signature(files),
        "time_count": time_count,
        "time_value": time_value,
    }


def snapshot_progress_reason(previous, current):
    if previous is None:
        return "watchdog_started"
    if current["time_count"] > previous["time_count"]:
        return "simulation_time_advanced"
    if current["time_value"] != previous["time_value"]:
        return "simulation_time_changed"
    if current["file_signature"] != previous["file_signature"]:
        return "output_file_updated"
    return None


def find_numerical_failure_marker(text):
    text = text.lower()
    for marker in NUMERICAL_FAILURE_MARKERS:
        if marker in text:
            return marker
    return None


def write_summary(summary_path, summary):
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)


def read_existing_summary_status(summary_path):
    if not summary_path.exists():
        return ""
    try:
        with open(summary_path) as f:
            return json.load(f).get("status", "")
    except (OSError, json.JSONDecodeError):
        return ""


def terminate_process_group(process, grace_seconds=30):
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return process.poll()

    try:
        return process.wait(timeout=grace_seconds)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        return process.wait()


def vtk_merge_script_path():
    script_dir = Path(os.environ.get("DUMUX_SCRIPT_DIR", Path(__file__).resolve().parent)).resolve()
    candidates = [
        script_dir / "vtk-merge-multi.py",
        Path(__file__).resolve().parent / "vtk-merge-multi.py",
        Path.cwd() / "vtk-merge-multi.py",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def run_vtk_merge(case_dir, simulation_name, ntasks, stdout_path, stderr_path):
    merge_script = vtk_merge_script_path()
    started_at = now_iso()
    result = {
        "enabled": ntasks > 1,
        "status": "skipped",
        "started_at": started_at,
        "ended_at": started_at,
        "script": str(merge_script),
        "input_dir": str(case_dir.resolve()),
        "num_parts": ntasks,
        "simulation_name": simulation_name,
        "returncode": None,
    }

    if ntasks <= 1:
        result["reason"] = "single_process_run"
        return result

    if not merge_script.exists():
        result["status"] = "failed"
        result["reason"] = "merge_script_missing"
        result["ended_at"] = now_iso()
        with open(stderr_path, "a", buffering=1) as stderr_file:
            stderr_file.write(f"\nWarning: VTK merge script not found at {merge_script}\n")
        return result

    cmd = [
        sys.executable,
        str(merge_script),
        "--input-dir", str(case_dir.resolve()),
        "--num-parts", str(ntasks),
        "--simulation-name", simulation_name,
    ]
    result["command"] = " ".join(cmd)

    with open(stdout_path, "a", buffering=1) as stdout_file, open(stderr_path, "a", buffering=1) as stderr_file:
        stdout_file.write(f"\n=== Post-processing VTK merge for {simulation_name} using {ntasks} parts ===\n")
        stdout_file.flush()
        completed = subprocess.run(
            cmd,
            cwd=case_dir,
            stdout=stdout_file,
            stderr=stderr_file,
            text=True,
        )

    result["returncode"] = completed.returncode
    result["ended_at"] = now_iso()
    if completed.returncode == 0:
        result["status"] = "success"
    else:
        result["status"] = "failed"
        result["reason"] = "merge_command_failed"

    return result


def wait_with_watchdog(
    process,
    case_dir,
    stdout_path,
    stderr_path,
    summary_path,
    summary,
    stall_timeout_seconds,
    startup_grace_seconds,
    watchdog_check_seconds,
    hard_timeout_seconds,
    numerical_stall_timeout_seconds,
):
    start_monotonic = time.monotonic()
    last_progress_monotonic = start_monotonic
    last_time_progress_monotonic = start_monotonic
    last_progress_time = now_iso()
    last_progress_reason = "process_started"
    last_snapshot = progress_snapshot(case_dir, stdout_path, stderr_path)
    last_time_count = last_snapshot["time_count"]
    last_time_value = last_snapshot["time_value"]
    last_failure_marker = None
    failure_marker_time = None
    failure_marker_monotonic = None
    output_offsets = {
        stdout_path: 0,
        stderr_path: 0,
    }

    while True:
        returncode = process.poll()
        current_snapshot = progress_snapshot(case_dir, stdout_path, stderr_path)
        progress_reason = snapshot_progress_reason(last_snapshot, current_snapshot)
        current_monotonic = time.monotonic()

        if progress_reason:
            last_snapshot = current_snapshot
            last_progress_monotonic = current_monotonic
            last_progress_time = now_iso()
            last_progress_reason = progress_reason

        time_advanced = (
            current_snapshot["time_count"] > last_time_count
            or current_snapshot["time_value"] != last_time_value
        )
        if time_advanced:
            last_time_count = current_snapshot["time_count"]
            last_time_value = current_snapshot["time_value"]
            last_time_progress_monotonic = current_monotonic
            failure_marker_time = None
            failure_marker_monotonic = None

        new_output = []
        for output_path in (stdout_path, stderr_path):
            text, new_offset = read_new_text(output_path, output_offsets[output_path])
            output_offsets[output_path] = new_offset
            new_output.append(text)

        marker = find_numerical_failure_marker("\n".join(new_output))
        if marker:
            last_failure_marker = marker
            if failure_marker_monotonic is None:
                failure_marker_time = now_iso()
                failure_marker_monotonic = current_monotonic

        elapsed_seconds = current_monotonic - start_monotonic
        idle_seconds = current_monotonic - last_progress_monotonic
        time_idle_seconds = current_monotonic - last_time_progress_monotonic
        numerical_failure_idle_seconds = (
            None if failure_marker_monotonic is None
            else current_monotonic - failure_marker_monotonic
        )
        summary["watchdog"] = {
            "enabled": True,
            "elapsed_seconds": round(elapsed_seconds, 1),
            "idle_seconds": round(idle_seconds, 1),
            "time_idle_seconds": round(time_idle_seconds, 1),
            "startup_grace_seconds": startup_grace_seconds,
            "stall_timeout_seconds": stall_timeout_seconds,
            "hard_timeout_seconds": hard_timeout_seconds,
            "numerical_stall_timeout_seconds": numerical_stall_timeout_seconds,
            "last_progress_time": last_progress_time,
            "last_progress_reason": last_progress_reason,
            "latest_json_time_count": current_snapshot["time_count"],
            "latest_json_time_value": current_snapshot["time_value"],
            "last_numerical_failure_marker": last_failure_marker,
            "numerical_failure_since": failure_marker_time,
            "numerical_failure_idle_seconds": (
                None if numerical_failure_idle_seconds is None
                else round(numerical_failure_idle_seconds, 1)
            ),
        }
        write_summary(summary_path, summary)

        if returncode is not None:
            return returncode, None

        if hard_timeout_seconds > 0 and elapsed_seconds >= hard_timeout_seconds:
            summary["watchdog"]["termination_reason"] = "hard_timeout"
            write_summary(summary_path, summary)
            terminate_process_group(process)
            return process.returncode, "failed_timeout"

        past_grace_period = elapsed_seconds >= startup_grace_seconds
        stalled = idle_seconds >= stall_timeout_seconds
        numerical_stalled = (
            numerical_stall_timeout_seconds > 0
            and numerical_failure_idle_seconds is not None
            and numerical_failure_idle_seconds >= numerical_stall_timeout_seconds
        )
        if past_grace_period and numerical_stalled:
            summary["watchdog"]["termination_reason"] = "failed_numerical_stall"
            write_summary(summary_path, summary)
            terminate_process_group(process)
            return process.returncode, "failed_numerical_stall"

        if past_grace_period and stalled:
            status = "failed_numerical_stall" if last_failure_marker else "failed_stalled"
            summary["watchdog"]["termination_reason"] = status
            write_summary(summary_path, summary)
            terminate_process_group(process)
            return process.returncode, status

        time.sleep(watchdog_check_seconds)


def run_one_case(
    row,
    iter_id,
    chunk_id,
    ntasks,
    executable,
    mpi_runner,
    stall_timeout_seconds,
    startup_grace_seconds,
    watchdog_check_seconds,
    hard_timeout_seconds,
    numerical_stall_timeout_seconds,
    skip_existing_case_dirs,
):
    case_id = int(row["case_id"])
    case_dir = Path("cases") / iter_id / f"chunk_{chunk_id}" / f"case_{case_id:03d}_{row['name']}"
    executable_path = Path(executable).resolve()
    mpi_runner_cmd = resolve_mpi_runner(mpi_runner)

    cmd = [
        mpi_runner_cmd, "-np", str(ntasks),
        str(executable_path),
        str((case_dir / "params.input").resolve())
    ]

    summary = {
        "case_id": case_id,
        "iter_id": iter_id,
        "chunk_id": chunk_id,
        "name": row["name"],
        "case_dir": str(case_dir.resolve()),
        "command": " ".join(cmd),
        "executable": str(executable_path),
        "mpi_runner": mpi_runner_cmd,
        "status": "failed",
        "inputs": row,
        "watchdog_config": {
            "stall_timeout_seconds": stall_timeout_seconds,
            "startup_grace_seconds": startup_grace_seconds,
            "watchdog_check_seconds": watchdog_check_seconds,
            "hard_timeout_seconds": hard_timeout_seconds,
            "numerical_stall_timeout_seconds": numerical_stall_timeout_seconds,
        },
    }

    stdout_path = case_dir / "stdout.txt"
    stderr_path = case_dir / "stderr.txt"
    summary_path = case_dir / "case_summary.json"

    if skip_existing_case_dirs and case_dir.exists():
        existing_status = read_existing_summary_status(summary_path)
        print(
            f"[chunk {chunk_id}] skipping case {case_id}: case folder already exists"
            + (f" (existing status: {existing_status})" if existing_status else ""),
            flush=True,
        )
        if not summary_path.exists():
            summary["status"] = "skipped_existing_case_dir"
            summary["simulation_status"] = "skipped"
            summary["skip_reason"] = "case_directory_already_exists"
            summary["skip_time"] = now_iso()
            json_files = sorted(case_dir.glob("*.json"))
            summary["json_outputs"] = [p.name for p in json_files]
            write_summary(summary_path, summary)
        return

    case_dir.mkdir(parents=True, exist_ok=True)
    write_params_file(case_dir, row)
    (case_dir / "runner_command.txt").write_text(" ".join(cmd) + "\n")

    print(f"[chunk {chunk_id}] starting case {case_id}: {row['name']}", flush=True)
    print(f"[chunk {chunk_id}] command: {' '.join(cmd)}", flush=True)

    try:
        if not executable_path.exists():
            message = f"Executable not found: {executable_path}\n"
            stdout_path.write_text("")
            stderr_path.write_text(message)
            summary["status"] = "failed_precheck"
            summary["exception"] = message.strip()
            print(f"[chunk {chunk_id}] case {case_id} failed before mpirun: {message.strip()}", flush=True)
            return

        mpi_runner_path = Path(mpi_runner_cmd)
        if mpi_runner_path.parent != Path("."):
            runner_missing = not mpi_runner_path.exists()
        else:
            runner_missing = shutil.which(mpi_runner_cmd) is None

        if runner_missing:
            message = f"Could not find MPI runner: {mpi_runner_cmd}\n"
            stdout_path.write_text("")
            stderr_path.write_text(message)
            summary["status"] = "failed_precheck"
            summary["exception"] = message.strip()
            print(f"[chunk {chunk_id}] case {case_id} failed before simulation: {message.strip()}", flush=True)
            return

        summary["status"] = "running"
        summary["start_time"] = now_iso()
        write_summary(summary_path, summary)

        env = os.environ.copy()
        env["PYTHONUNBUFFERED"] = "1"
        with open(stdout_path, "w", buffering=1) as stdout_file, open(stderr_path, "w", buffering=1) as stderr_file:
            process = subprocess.Popen(
                cmd,
                cwd=case_dir,
                stdout=stdout_file,
                stderr=stderr_file,
                text=True,
                env=env,
                start_new_session=True,
            )
            summary["pid"] = process.pid
            write_summary(summary_path, summary)
            returncode, watchdog_status = wait_with_watchdog(
                process=process,
                case_dir=case_dir,
                stdout_path=stdout_path,
                stderr_path=stderr_path,
                summary_path=summary_path,
                summary=summary,
                stall_timeout_seconds=stall_timeout_seconds,
                startup_grace_seconds=startup_grace_seconds,
                watchdog_check_seconds=watchdog_check_seconds,
                hard_timeout_seconds=hard_timeout_seconds,
                numerical_stall_timeout_seconds=numerical_stall_timeout_seconds,
            )

        summary["returncode"] = returncode
        summary["simulation_end_time"] = now_iso()

        if watchdog_status:
            summary["status"] = watchdog_status
        elif returncode == 0:
            summary["simulation_status"] = "success"
            summary["merge"] = run_vtk_merge(
                case_dir=case_dir,
                simulation_name=row["name"],
                ntasks=ntasks,
                stdout_path=stdout_path,
                stderr_path=stderr_path,
            )
            json_files = sorted(case_dir.glob("*.json"))
            summary["json_outputs"] = [p.name for p in json_files]
            summary["status"] = "success" if summary["merge"]["status"] in {"success", "skipped"} else "failed_merge"
        else:
            summary["simulation_status"] = "failed"
            summary["status"] = "failed"

        summary["end_time"] = now_iso()

    except Exception as e:
        summary["status"] = "failed"
        summary["exception"] = str(e)

    finally:
        write_summary(summary_path, summary)
        print(f"[chunk {chunk_id}] finished case {case_id}: {summary['status']}", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--iter-id", required=True)
    parser.add_argument("--chunk-id", type=int, required=True)
    parser.add_argument("--ntasks", type=int, required=True)
    parser.add_argument("--mpi-runner", default="mpirun")
    parser.add_argument("--executable", required=True)
    parser.add_argument(
        "--stall-timeout-min",
        type=float,
        default=env_float("DUMUX_STALL_TIMEOUT_MIN", 45.0),
        help="Stop a running case after this many minutes with no monitored progress.",
    )
    parser.add_argument(
        "--startup-grace-min",
        type=float,
        default=env_float("DUMUX_STARTUP_GRACE_MIN", 20.0),
        help="Minimum minutes before the no-progress watchdog can stop a case.",
    )
    parser.add_argument(
        "--watchdog-check-sec",
        type=float,
        default=env_float("DUMUX_WATCHDOG_CHECK_SEC", 600.0),
        help="Seconds between progress checks.",
    )
    parser.add_argument(
        "--case-timeout-hours",
        type=float,
        default=env_float("DUMUX_CASE_TIMEOUT_HOURS", 0.0),
        help="Optional hard wall-clock timeout per case. Use 0 to disable.",
    )
    parser.add_argument(
        "--numerical-stall-timeout-min",
        type=float,
        default=env_float("DUMUX_NUMERICAL_STALL_TIMEOUT_MIN", 30.0),
        help=(
            "After a known numerical failure marker appears, stop the case if "
            "simulation time does not advance for this many minutes. Use 0 to disable."
        ),
    )
    parser.add_argument(
        "--rerun-existing",
        action="store_true",
        help="Rerun cases even when their case directory already exists.",
    )
    args = parser.parse_args()
    stall_timeout_seconds = max(1.0, args.stall_timeout_min * 60.0)
    startup_grace_seconds = max(0.0, args.startup_grace_min * 60.0)
    watchdog_check_seconds = max(1.0, args.watchdog_check_sec)
    hard_timeout_seconds = max(0.0, args.case_timeout_hours * 3600.0)
    numerical_stall_timeout_seconds = max(0.0, args.numerical_stall_timeout_min * 60.0)

    rows = read_chunk_rows(args.manifest, args.chunk_id)
    print("=== run_chunk.py startup ===", flush=True)
    print(f"time: {now_iso()}", flush=True)
    print(f"manifest: {Path(args.manifest).resolve()}", flush=True)
    print(f"iter_id: {args.iter_id}", flush=True)
    print(f"chunk_id: {args.chunk_id}", flush=True)
    print(f"ntasks: {args.ntasks}", flush=True)
    print(f"mpi_runner: {resolve_mpi_runner(args.mpi_runner)}", flush=True)
    print(f"executable: {Path(args.executable).resolve()}", flush=True)
    print(f"rows in this chunk: {len(rows)}", flush=True)
    print(f"stall_timeout_min: {args.stall_timeout_min}", flush=True)
    print(f"startup_grace_min: {args.startup_grace_min}", flush=True)
    print(f"watchdog_check_sec: {args.watchdog_check_sec}", flush=True)
    print(f"case_timeout_hours: {args.case_timeout_hours}", flush=True)
    print(f"numerical_stall_timeout_min: {args.numerical_stall_timeout_min}", flush=True)
    print(f"skip_existing_case_dirs: {not args.rerun_existing}", flush=True)

    if not rows:
        raise SystemExit(f"No manifest rows found for chunk_id={args.chunk_id}")

    # Each SLURM array task runs all manifest rows for its chunk sequentially.
    for row in rows:
        run_one_case(
            row=row,
            iter_id=args.iter_id,
            chunk_id=args.chunk_id,
            ntasks=args.ntasks,
            mpi_runner=args.mpi_runner,
            executable=args.executable,
            stall_timeout_seconds=stall_timeout_seconds,
            startup_grace_seconds=startup_grace_seconds,
            watchdog_check_seconds=watchdog_check_seconds,
            hard_timeout_seconds=hard_timeout_seconds,
            numerical_stall_timeout_seconds=numerical_stall_timeout_seconds,
            skip_existing_case_dirs=not args.rerun_existing,
        )


if __name__ == "__main__":
    main()
