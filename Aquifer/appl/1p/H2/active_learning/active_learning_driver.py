"""
CLI driver for the DuMuX <-> GP active-learning coupling.

Example:
    python3 active_learning_driver.py \
        --training-csv ../results/pool_top70_summary.csv \
        --candidate-pool ../simulation_input_pool_500.csv \
        --target-prefix cycle \
        --max-cycles 10 \
        --batch-size 20 \
        --output ../manifests/active_learning_next.csv

Output CSV: rows from the candidate pool (untouched columns + score) ready to be
fed into controller.py via --candidates-csv.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

# Allow running as a script from this folder.
HERE = Path(__file__).resolve().parent
if str(HERE.parent) not in sys.path:
    sys.path.insert(0, str(HERE.parent))

from active_learning.pipeline import (   # noqa: E402
    CONTINUOUS_FEATURES,
    FeaturePrep,
    _normalize_candidate_df,
    diagnose_failures,
    exclude_already_run,
    load_summary,
    score_candidates,
    train_surrogate,
    wide_to_long,
)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--training-csv", type=Path, required=True,
                   help="Aggregated DuMuX summary CSV (output of aggregate_iteration.py).")
    p.add_argument("--candidate-pool", type=Path, required=True,
                   help="Simulation input pool CSV (e.g. simulation_input_pool_500.csv).")
    p.add_argument("--target-prefix", default="cycle",
                   help="Column prefix used for cycle-end RF columns. Default: 'cycle'")
    p.add_argument("--max-cycles", type=int, default=10)
    p.add_argument("--batch-size", type=int, default=20)
    p.add_argument("--aggregate", default="mean", choices=["mean", "max", "sum"],
                   help="How to aggregate per-cycle acquisition into a per-case score.")
    p.add_argument("--output", type=Path, required=True,
                   help="CSV path for the recommended next batch.")
    p.add_argument("--diagnostics-json", type=Path, default=None,
                   help="Optional path to dump failure diagnostics as JSON.")
    p.add_argument("--no-success-filter", action="store_true",
                   help="Use all rows for training (default keeps only summary_status==success).")
    return p.parse_args()


def main():
    args = parse_args()

    print(f"[1/6] Loading training summary: {args.training_csv}")
    summary = load_summary(args.training_csv)
    print(f"    completed simulations loaded (rows): {len(summary)}")

    print(f"[2/6] Failure diagnostics")
    diag = diagnose_failures(summary)
    print(f"    success={diag['n_success']}  failed={diag['n_failed']}  total={diag['n_total']}")
    if args.diagnostics_json:
        args.diagnostics_json.parent.mkdir(parents=True, exist_ok=True)
        args.diagnostics_json.write_text(json.dumps(diag, indent=2))
        print(f"    wrote diagnostics: {args.diagnostics_json}")

    print(f"[3/6] Wide -> long target conversion (cumulative RF)")
    long_df = wide_to_long(
        summary,
        max_cycles=args.max_cycles,
        target_prefix=args.target_prefix,
        require_success=not args.no_success_filter,
    )
    print(f"    long-format training rows: {len(long_df)}")
    print(f"    target used: cumulative_rf  (one scalar; cycle_no is an input)")
    if len(long_df) < 5:
        raise SystemExit("ERROR: not enough training rows; cannot fit GP.")

    print(f"[4/6] Loading candidate pool: {args.candidate_pool}")
    pool_raw = pd.read_csv(args.candidate_pool, low_memory=False)
    pool = _normalize_candidate_df(pool_raw)
    print(f"    candidate pool rows: {len(pool)}")

    pool_left, n_excluded = exclude_already_run(pool, summary)
    print(f"    already-run candidates excluded: {n_excluded}")
    print(f"    remaining candidates: {len(pool_left)}")
    if pool_left.empty:
        raise SystemExit("ERROR: no candidates left after exclusion.")

    print(f"[5/6] Fitting GP surrogate on 12-D inputs")
    prep = FeaturePrep().fit(long_df, candidate_df=pool_left, max_cycles=args.max_cycles)
    X = prep.transform(long_df)
    y = long_df["cumulative_rf"].to_numpy(dtype=float)
    print(f"    training matrix: X={X.shape}, y={y.shape}")
    model = train_surrogate(X, y)

    print(f"[6/6] Scoring candidates (acquisition: predictive_std, "
          f"aggregate over cycles: {args.aggregate})")
    scored = score_candidates(model, prep, pool_left,
                              max_cycles=args.max_cycles,
                              aggregate=args.aggregate)
    selected = scored.head(args.batch_size).copy()
    print(f"    selected candidates: {len(selected)}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    # Re-attach all original pool columns (keeps source_row_index, source_field_name, etc.).
    keep_cols = list(pool_raw.columns)
    extras = [c for c in selected.columns
              if c not in keep_cols and c not in {"pool_row_index"}]
    out = pool_raw.iloc[selected["pool_row_index"].astype(int).values].copy().reset_index(drop=True)
    out["pool_row_index"] = selected["pool_row_index"].values
    for c in extras:
        out[c] = selected[c].values
    out.to_csv(args.output, index=False)
    print(f"    wrote: {args.output}")

    print("\nDone.")


if __name__ == "__main__":
    main()
