# active_learning/

First practical coupling between the DuMuX simulation workflow in this folder
and a Gaussian-process surrogate that predicts **cumulative recovery factor**
at a given operational cycle number.

Layout
------
- `pipeline.py` — data loader, wide-to-long converter, feature preparation
  (8 continuous + 4 one-hot = **12-D** GP input), GP training & prediction,
  candidate scoring, failure diagnostics.
- `active_learning_driver.py` — CLI that ties everything together and writes
  a "next batch" CSV that controller.py can consume.
- `test_dry_run.py` — synthetic test (no DuMuX needed) verifying the pipeline.

Inputs to the GP
----------------
8 continuous (normalized to [0, 1] using min/max from training **and** candidates):

    flow_rate_sm3_day, cycle_length_days, cg_ratio,
    permeability_md, porosity, pressure_mpa, temperature_c, cycle_no

4 one-hot (cushion-gas type):

    is_CH4, is_H2, is_N2, is_CO2

Single scalar target: `cumulative_rf` (long-format expansion of
`cycle01_end_rf … cycle10_end_rf`).

CLI
---
```bash
python3 active_learning/active_learning_driver.py \
    --training-csv results/pool_top70_summary.csv \
    --candidate-pool simulation_input_pool_500.csv \
    --target-prefix cycle \
    --max-cycles 10 \
    --batch-size 20 \
    --aggregate mean \
    --output manifests/active_learning_next.csv
```

The driver prints diagnostics for every step:
completed simulations loaded, long-format training rows, target,
candidate pool size, already-run candidates excluded,
acquisition aggregation method, selected candidates path.

Failed-simulation strategy (recommendation)
-------------------------------------------
For now: **strategy 4 — diagnostics only**. The driver always reports per-feature
ranges of failed vs. successful runs and the cushion-gas-type counts, optionally
dumped via `--diagnostics-json`. Once we have ≥ a few dozen failures we can
graduate to strategy 1 (rule-based feasibility filter built from those ranges)
or 2 (a small classifier). Both are cheap to bolt on top of the same loader.

Controller integration
----------------------
`controller.py` gained one optional flag:

```bash
python3 controller.py \
    --candidates-csv manifests/active_learning_next.csv \
    --batch-id al_round_01 \
    --jobs 4 --cases-per-job 5 \
    --manifest-only          # drop this flag to actually submit SLURM jobs
```

When `--candidates-csv` is set the controller ignores
`--pool-start-row / --pool-head / --total-simulations` and runs exactly the
cases listed in the file (preserving each case's `pool_row_index`).
Original fixed-range usage is unchanged.

Dry-run test
------------
```bash
python3 active_learning/test_dry_run.py
```
fabricates a small summary CSV (with successes + failures + cycle RF columns)
and a candidate pool, runs the driver end-to-end, and verifies the output
schema. Requires `torch`, `gpytorch`, `botorch`, `pandas`.
