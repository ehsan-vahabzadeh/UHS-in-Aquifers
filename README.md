<div align="center">

# Mixing in Underground Hydrogen Storage

### High-fidelity DuMuX simulation, Gaussian-process surrogates, and active learning in one automated campaign

[![DuMuX](https://img.shields.io/badge/simulation-DuMuX-00629B?style=flat-square)](https://dumux.org/)
[![Python](https://img.shields.io/badge/automation-Python-3776AB?style=flat-square&logo=python&logoColor=white)](https://www.python.org/)
[![BoTorch](https://img.shields.io/badge/surrogate-BoTorch-8A2BE2?style=flat-square)](https://botorch.org/)
[![Slurm](https://img.shields.io/badge/HPC-Slurm-44A833?style=flat-square)](https://slurm.schedmd.com/)

**Simulate → aggregate → learn → select → resimulate**

</div>

This repository is a computational framework for studying gas mixing and hydrogen recovery during cyclic underground hydrogen storage (UHS). It connects a compositional porous-media model in **DuMuX** to a **Gaussian Process (GP)** surrogate and an **uncertainty-driven active-learning loop**. The result is an artifact-driven workflow that can prepare cases, dispatch them on Slurm, monitor failures, collect cycle-level outputs, retrain the surrogate, and recommend the next most informative simulations.

The central idea is simple: high-fidelity reservoir simulations are expensive, so the campaign should not spend equal compute everywhere. The GP learns from completed cases and quantifies what it does not yet know; active learning then directs the next simulation batch toward those uncertain regions.

<p align="center">
  <img src="appl/1p/H2/assets/H2Fraction.gif" width="760" alt="Hydrogen mole-fraction evolution during underground storage simulation">
</p>

## The automated feedback loop

```mermaid
flowchart LR
    A[Real reservoir data] --> B[Candidate-pool design]
    B --> C[Case manifest]
    C --> D[DuMuX simulations]
    D --> E[Watchdog + case JSON]
    E --> F[Automated aggregation]
    F --> G[GP surrogate]
    G --> H[Uncertainty scoring]
    H --> I[Next informative batch]
    I --> C

    style D fill:#00629B,color:#fff,stroke:#00466f
    style G fill:#6f42c1,color:#fff,stroke:#4d2d87
    style H fill:#8a2be2,color:#fff,stroke:#5d1d99
    style I fill:#2f8f46,color:#fff,stroke:#206332
```

Each stage produces a stable file consumed by the next stage: pool CSVs, batch manifests, case-level JSON, aggregated summaries, and ranked candidate CSVs. This keeps the workflow transparent, restartable, and easy to run interactively or automate in a larger campaign script.

## What is being simulated?

The main adaptive campaign lives in [`appl/1p/H2`](appl/1p/H2). Its high-fidelity model represents cyclic hydrogen injection and production in a saline aquifer using an axisymmetric reservoir domain.

| Model element | Implementation |
| --- | --- |
| Flow model | Transient, compositional two-phase porous-media flow (`TwoPNC`) |
| Phases | Aqueous and gas |
| Components | H₂O, H₂, CH₄, CO₂, and N₂ |
| Cushion-gas options | H₂, CH₄, CO₂, or N₂ |
| Geometry | 2D axisymmetric radial section, representing a cylindrical reservoir |
| Spatial discretization | Box finite-volume method for the primary campaign; TPFA target also provided |
| Transport | Advection, molecular diffusion, and compositional dispersion |
| Porous-medium relations | Brooks–Corey capillary pressure/relative permeability and Millington–Quirk effective diffusivity |
| Operations | Cushion-gas development followed by ten H₂ injection/production cycles |
| Key responses | Cumulative and per-cycle recovery factor, H₂ purity, water cut, pressure, component inventory, and material-balance error |

The campaign varies both operational decisions and reservoir conditions:

- flow rate, cycle length, and cushion-gas ratio;
- cushion-gas identity;
- permeability and porosity;
- reservoir pressure and temperature.

The simulation code also records explicit pressure-limit and time-step-collapse failures. These outcomes are retained by the pipeline, rather than silently discarded, so failed regions can be diagnosed and later incorporated into feasibility-aware sampling.

## Gaussian Process surrogate

The practical surrogate is implemented with **BoTorch** and **GPyTorch** in [`active_learning/pipeline.py`](appl/1p/H2/active_learning/pipeline.py). Completed wide-format simulation summaries are expanded into one training observation per case and operational cycle.

The GP input is 12-dimensional:

| Features | Encoding |
| --- | --- |
| Flow rate, cycle length, cushion-gas ratio | Min–max-normalized continuous variables |
| Permeability, porosity, pressure, temperature | Min–max-normalized continuous variables |
| Cycle number | Normalized continuous variable |
| CH₄, H₂, N₂, CO₂ cushion gas | Four one-hot variables |

The scalar target is cumulative hydrogen recovery factor at a given cycle. The model uses a scaled **Matérn-5/2 kernel with automatic relevance determination (ARD)** and a standardized outcome:

$$
y(\mathbf{x}) \sim \mathcal{GP}\!\left(m(\mathbf{x}),\,k_{\text{Mat\'ern-5/2}}(\mathbf{x},\mathbf{x}')\right).
$$

For every unrun physical case, the surrogate predicts a mean and standard deviation at each cycle. The practical acquisition score is predictive standard deviation, aggregated across cycles by the mean, maximum, or sum:

$$
a(\mathbf{x}) = \operatorname{Agg}_{c=1}^{C}\left[\sigma(\mathbf{x},c)\right].
$$

The highest-scoring cases become the next simulation batch. The experimental [`GP_ActiveLearning*.py`](appl/1p/H2/GP_ActiveLearning_v3.py) studies additionally compare uncertainty-only and diversity-aware strategies.

<p align="center">
  <img src="appl/1p/H2/gp_active_learning_plots/v3_comparison_uncertainty_diversity.png" width="760" alt="Comparison of Gaussian-process active-learning strategies">
</p>

## Why active learning?

A fixed grid or a very large random design can spend substantial compute on redundant cases. Active learning replaces that static design with an adaptive one:

1. Run a modest initial design.
2. Fit the GP to successful, cycle-resolved recovery factors.
3. Evaluate all remaining candidates with the surrogate.
4. Select the cases with the greatest predictive uncertainty.
5. Run and aggregate only that batch.
6. Add the new evidence and repeat.

This targets the simulation budget where it can most improve the surrogate. The candidate pool is itself designed carefully: reservoir properties remain as observed tuples from real data, density-aware clustering preserves sparse permeability regimes, and Latin hypercube sampling covers the operational variables. Low-permeability reservoirs receive a configurable flow-rate cap to avoid obviously infeasible combinations.

## What the automation handles

| Stage | Script | Automated responsibility |
| --- | --- | --- |
| Pool construction | [`generate_simulation_pool.py`](appl/1p/H2/generate_simulation_pool.py) | Cleans reservoir data, selects representative tuples, samples operations, and writes coverage diagnostics |
| Campaign setup | [`controller.py`](appl/1p/H2/controller.py) | Converts units, derives cycling times, creates deterministic case names, chunks work, and writes a manifest |
| HPC dispatch | [`run_chunk_array.sh`](appl/1p/H2/run_chunk_array.sh) | Launches manifest chunks as a Slurm array with MPI |
| Case execution | [`run_chunk.py`](appl/1p/H2/run_chunk.py) | Creates isolated inputs, skips existing cases, monitors progress, terminates stalled runs, and records status |
| Aggregation | [`aggregate_iteration.py`](appl/1p/H2/aggregate_iteration.py) | Joins manifest inputs, runner status, simulation JSON, cycle metrics, and failure metadata into CSV/JSON |
| Surrogate update | [`active_learning_driver.py`](appl/1p/H2/active_learning/active_learning_driver.py) | Fits the GP, excludes completed candidates, ranks the remainder, and writes the next-batch CSV |

When `controller.py` submits a campaign, it also submits aggregation with an `afterany` dependency on the simulation array. Aggregation therefore runs even when some cases fail, preserving the evidence needed for campaign diagnostics. Each active-learning round is then launched from the generated next-batch CSV; no manual editing of individual simulation input files is needed.

## Quick start

### 1. Requirements

The simulation stack is based on **DUNE 2.9** and **DuMuX 3.8**. A reproduction installer with pinned revisions is provided as [`install_vahabzadeh2026a.py`](install_vahabzadeh2026a.py). Building and running the full model additionally requires a C++ compiler, CMake, MPI, and OpenMP; Slurm is required only for the cluster workflow.

The active-learning tools require Python 3.10 or newer and the following packages:

```bash
python3 -m pip install numpy pandas scipy scikit-learn matplotlib \
  torch gpytorch botorch
```

For a configured DUNE/DuMuX workspace, configure and build this module with:

```bash
cmake -S . -B build-cmake
cmake --build build-cmake --target appl_1p2pnc_box_H2 -j
```

Cluster-specific settings such as the Slurm partition, MPI path, task count, and library paths are intentionally visible in [`run_chunk_array.sh`](appl/1p/H2/run_chunk_array.sh) and should be adapted to the target system.

### 2. Validate the ML workflow without DuMuX

The synthetic dry run creates successful and failed cases, trains the GP, scores a candidate pool, and validates the selected-batch schema:

```bash
python3 appl/1p/H2/active_learning/test_dry_run.py
```

### 3. Generate a candidate pool

An example 500-case pool is included. To regenerate it from the matched reservoir table and reproduce its diagnostic plots:

```bash
python3 appl/1p/H2/generate_simulation_pool.py
```

Useful outputs are:

- `simulation_input_pool_500.csv` — the campaign's finite candidate set;
- `simulation_pool_diagnostics/` — distributions, empirical CDFs, property cross-plots, and stratum coverage.

### 4. Launch an initial simulation batch

First inspect the generated manifest without submitting work:

```bash
python3 appl/1p/H2/controller.py \
  --pool-csv appl/1p/H2/simulation_input_pool_500.csv \
  --batch-id bootstrap_01 \
  --pool-start-row 1 \
  --total-simulations 20 \
  --jobs 4 \
  --cases-per-job 5 \
  --manifest-only
```

On the configured Slurm login node, remove `--manifest-only` to submit the array. The controller automatically schedules the dependent aggregation job. With the default source-tree layout, runtime artifacts are written beneath `build-cmake/appl/1p/H2/`:

```text
manifests/bootstrap_01.csv
cases/bootstrap_01/<case-name>/...
logs/bootstrap_01_*.{out,err}
results/bootstrap_01_summary.{csv,json}
```

The controller requires `total simulations = jobs × cases per job`; this same rule applies to active-learning batches.

### 5. Train the GP and select the next batch

```bash
python3 appl/1p/H2/active_learning/active_learning_driver.py \
  --training-csv build-cmake/appl/1p/H2/results/bootstrap_01_summary.csv \
  --candidate-pool appl/1p/H2/simulation_input_pool_500.csv \
  --max-cycles 10 \
  --batch-size 20 \
  --aggregate mean \
  --diagnostics-json build-cmake/appl/1p/H2/results/bootstrap_01_failures.json \
  --output build-cmake/appl/1p/H2/manifests/al_round_01_candidates.csv
```

The driver uses successful simulations for GP training by default, reports the failed cases separately, excludes every pool row that already appears in the aggregated summary, and preserves the original candidate columns in its output.

### 6. Run the selected cases

```bash
python3 appl/1p/H2/controller.py \
  --candidates-csv build-cmake/appl/1p/H2/manifests/al_round_01_candidates.csv \
  --batch-id al_round_01 \
  --jobs 4 \
  --cases-per-job 5
```

After aggregation, combine the completed summaries used for training, rerun the active-learning driver, and submit the newly selected manifest. That repeating file-to-file handoff is the campaign loop.

## Repository map

```text
.
├── appl/
│   ├── 1p/H2/                       # adaptive UHS simulation campaign
│   │   ├── active_learning/         # practical GP training and selection pipeline
│   │   ├── controller.py            # manifests and Slurm orchestration
│   │   ├── run_chunk.py             # robust case runner and watchdog
│   │   ├── aggregate_iteration.py   # cycle-level result aggregation
│   │   ├── generate_simulation_pool.py
│   │   ├── main.cc                  # DuMuX executable entry point
│   │   ├── problem.hh               # cycling, wells, metrics, and safety limits
│   │   ├── fluidsystems/mixture.hh  # H2O–H2–CH4–CO2–N2 fluid system
│   │   └── params.input             # example physical/numerical configuration
│   ├── H2, CH4, CO2, N2/            # additional compositional scenarios
│   ├── Cartesian/                   # Cartesian-domain variants
│   └── Test1, Test2/                # biogeochemical/test configurations
├── dumux/                           # project-specific DuMuX extensions
├── PostFun/                         # broader analysis, ANN, sensitivity, and optimisation tools
├── install_vahabzadeh2026a.py       # pinned DUNE/DuMuX reproduction installer
├── CMakeLists.txt
└── dune.module
```

For lower-level details of the adaptive campaign, see the [`appl/1p/H2` workflow notes](appl/1p/H2/README.md) and the [`active_learning` technical notes](appl/1p/H2/active_learning/README.md).

## Reproducibility and failure handling

Simulation campaigns fail in more interesting ways than ordinary batch jobs. This workflow makes those states explicit:

- every case runs in its own directory with a generated `params.input`;
- manifest rows retain their original `pool_row_index`, providing stable identity across iterations;
- existing case directories are not overwritten by default;
- a watchdog monitors output and simulated-time progress;
- numerical-stall markers, wall-clock timeouts, pressure cutoffs, and collapsed time steps receive distinct statuses;
- Slurm logs, commands, timestamps, return codes, and merge outcomes are stored with each case;
- aggregation tolerates partial campaign failure and reports missing summaries;
- failure distributions can be exported as JSON alongside active-learning results.

These choices make it possible to audit why a candidate entered the training set, why another failed, and exactly which manifest produced each result.

## Current scope

The production active-learning driver currently uses **predictive uncertainty** over a finite candidate pool. Diversity-aware selection exists in the experimental GP scripts but is not yet the default in the operational driver. Failed simulations are diagnosed and excluded from GP regression; they are not yet used to train a feasibility classifier. The round-to-round interfaces are automated and composable, while advancing to the next round remains an explicit operator action so diagnostics can be reviewed before more HPC work is submitted.

That boundary is deliberate: the expensive parts run unattended, but scientific control remains visible.
