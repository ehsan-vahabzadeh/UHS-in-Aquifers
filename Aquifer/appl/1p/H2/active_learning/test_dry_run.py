"""
Synthetic dry-run for the active-learning pipeline (no DuMuX required).

Generates:
  - a fake aggregate_iteration.py-style summary CSV (with cycleXX_end_rf columns,
    a few failed rows, and the manifest input columns)
  - a small candidate-pool CSV (snake_case columns matching simulation_input_pool_500.csv)

Then runs the active_learning_driver as a subprocess and checks the output.

Usage:
    python3 active_learning/test_dry_run.py
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
DRIVER = HERE / "active_learning_driver.py"
CG_TYPES = ["CH4", "H2", "N2", "CO2"]


def fake_summary(n_success: int = 25, n_failed: int = 5, max_cycles: int = 10,
                 seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n_success + n_failed):
        is_fail = i >= n_success
        cg = rng.choice(CG_TYPES)
        flow = float(rng.uniform(5e4, 1.5e6))
        cycle_len = float(rng.uniform(60, 400))
        cgr = float(rng.uniform(2.0, 4.0))
        perm = float(rng.choice([50, 200, 500, 1000, 1600]))
        poro = float(rng.choice([0.15, 0.18, 0.22, 0.27]))
        press = float(rng.uniform(10, 35))
        temp = float(rng.uniform(50, 110))

        row = {
            "case_id": i,
            "pool_row_index": i,
            "summary_status": "failed" if is_fail else "success",
            "FlowRateSm3Day": flow,
            "CycleLengthDays": cycle_len,
            "CGRatio": cgr,
            "PermeabilityMd": perm,
            "Porosity": poro,
            "PressureMPa": press,
            "TemperatureC": temp,
            "CushionGasType": cg,
        }
        # Cycle RF: smooth function of inputs + cycle index, only for success rows
        if not is_fail:
            base = 0.05 + 0.4 * (poro - 0.15) / 0.15 + 0.1 * np.tanh((press - 20) / 10)
            for c in range(1, max_cycles + 1):
                rf = base + 0.04 * c + 0.02 * rng.standard_normal()
                row[f"cycle{c:02d}_end_rf"] = float(np.clip(rf, 0.0, 1.0))
        rows.append(row)
    return pd.DataFrame(rows)


def fake_pool(n: int = 60, seed: int = 1) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(n):
        rows.append({
            "flow_rate_sm3_day": float(rng.uniform(5e4, 1.5e6)),
            "cycle_length_days": float(rng.uniform(60, 400)),
            "cg_ratio": float(rng.uniform(2.0, 4.0)),
            "cg_type": str(rng.choice(CG_TYPES)),
            "permeability_md": float(rng.choice([50, 200, 500, 1000, 1600])),
            "porosity": float(rng.choice([0.15, 0.18, 0.22, 0.27])),
            "pressure_mpa": float(rng.uniform(10, 35)),
            "temperature_c": float(rng.uniform(50, 110)),
        })
    return pd.DataFrame(rows)


def main():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        summary_csv = td / "summary.csv"
        pool_csv = td / "pool.csv"
        out_csv = td / "selected.csv"
        diag_json = td / "diag.json"

        fake_summary().to_csv(summary_csv, index=False)
        fake_pool().to_csv(pool_csv, index=False)

        cmd = [
            sys.executable, str(DRIVER),
            "--training-csv", str(summary_csv),
            "--candidate-pool", str(pool_csv),
            "--max-cycles", "10",
            "--batch-size", "8",
            "--output", str(out_csv),
            "--diagnostics-json", str(diag_json),
        ]
        print("Running:", " ".join(cmd))
        r = subprocess.run(cmd, capture_output=True, text=True)
        print(r.stdout)
        if r.returncode != 0:
            print(r.stderr, file=sys.stderr)
            raise SystemExit(f"driver failed (rc={r.returncode})")

        out = pd.read_csv(out_csv)
        print(f"\nselected file: {out_csv}  rows={len(out)}  cols={list(out.columns)}")
        assert len(out) == 8, "batch size mismatch"
        for col in ("flow_rate_sm3_day", "cg_type", "score", "pool_row_index"):
            assert col in out.columns, f"missing {col}"
        print("OK -- dry-run passed.")


if __name__ == "__main__":
    main()
