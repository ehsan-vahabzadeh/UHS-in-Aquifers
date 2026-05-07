# Gaussian Process Active Learning for Automated Simulation Design

This project demonstrates an end-to-end workflow for managing expensive simulation campaigns with Python automation, Gaussian Process surrogate modelling, and active learning. The domain example is subsurface hydrogen storage simulation, but the main contribution is the computational workflow: ingesting completed simulation results, training a surrogate model, selecting high-value new candidates, generating batch manifests, running simulations, and aggregating outputs for the next iteration.

The workflow is designed for technical audiences interested in scientific software, simulation engineering, data pipelines, machine learning for expensive models, and HPC batch automation.

## Motivation

High-fidelity simulations can be slow and expensive. Running thousands of cases with brute-force sampling often wastes compute on points that add little information.

At the same time, surrogate models need representative training data. If the training set is poorly chosen, the surrogate may be fast but unreliable.

This project addresses that problem with an adaptive loop:

- Use completed simulations as training data.
- Train a Gaussian Process surrogate model.
- Use predictive uncertainty to identify where the model needs more information.
- Select the next batch of simulations from a candidate pool.
- Run only the most useful new cases.
- Aggregate the results and repeat.

The goal is practical: reduce unnecessary simulation cost while improving coverage of important input regions.

## What This Project Does

The workflow connects simulation outputs, surrogate modelling, candidate selection, and batch execution:

1. Read completed simulation summaries.
2. Clean and reshape outputs into a training dataset.
3. Train or update a Gaussian Process surrogate model.
4. Estimate predictive uncertainty over a candidate pool.
5. Score and select high-value candidate simulations.
6. Write a next-batch CSV or simulation manifest.
7. Submit/run simulation batches.
8. Aggregate simulation outputs into a summary file.
9. Repeat the loop with the expanded dataset.

## Workflow Architecture

```text
Completed simulations
        |
        v
Aggregated summary CSV
        |
        v
Training dataset
        |
        v
Gaussian Process surrogate model
        |
        v
Uncertainty / acquisition scoring
        |
        v
Selected new candidates
        |
        v
Simulation batch manifest
        |
        v
HPC / Slurm execution
        |
        v
Aggregated outputs
        |
        +---- repeat
```

## Key Features

### Gaussian Process Surrogate Modelling

The active-learning pipeline trains a Gaussian Process model on completed simulation results. The current practical target is cumulative recovery factor at a given cycle number, represented as a single scalar output.

The model input combines continuous simulation parameters with encoded categorical variables:

| Feature group | Description |
| --- | --- |
| Continuous inputs | Flow rate, cycle length, cushion-gas ratio, permeability, porosity, pressure, temperature, cycle number |
| Categorical input | Cushion-gas type, encoded as one-hot columns |
| Target | Cumulative recovery factor for a cycle |

### Uncertainty-Aware Active Learning

Gaussian Processes provide both a prediction and an uncertainty estimate. The workflow uses this uncertainty to identify candidate simulations that are expected to improve the surrogate model most.

This supports adaptive sampling instead of blindly generating large batches of simulations.

### Candidate-Pool Generation

The project includes scripts for preparing and diagnosing a simulation input pool. Candidate pools can be generated from reservoir/property data, cleaned, normalized, and checked for coverage before simulations are selected.

### Diversity-Aware Sampling

The active-learning prototypes include uncertainty-based and diversity-aware selection logic. This helps avoid choosing a batch of candidates that are all clustered in the same region of the input space.

### Batch Workflow Orchestration

Python orchestration scripts generate simulation manifests, split cases into chunks, and coordinate simulation batches. The workflow is built to support repeated campaign rounds rather than one-off manual runs.

### Simulation Output Aggregation

Post-processing scripts collect case-level JSON outputs and runner summaries into consolidated CSV files. These aggregated files are used as training data for the next surrogate-model update.

### HPC / Slurm-Ready Execution

The project includes Slurm array-job scripts for running simulation chunks on an HPC system. Each chunk can process multiple cases sequentially, while the array structure allows several chunks to run in parallel.

### Diagnostics and Model Tracking

The workflow includes diagnostic outputs for candidate-pool coverage, active-learning behavior, and surrogate-model performance. These diagnostics are useful for checking whether the adaptive loop is selecting meaningful new cases.

## Repository Structure

```text
H2/
|-- active_learning/
|   |-- active_learning_driver.py   # CLI for GP training and candidate selection
|   |-- pipeline.py                 # Data loading, feature prep, GP model, scoring
|   |-- test_dry_run.py             # Synthetic end-to-end pipeline test
|   `-- README.md                   # Technical notes for active-learning coupling
|
|-- controller.py                   # Generates manifests and submits batch jobs
|-- run_chunk.py                    # Runs one simulation chunk and monitors progress
|-- run_chunk_array.sh              # Slurm array entry point for simulation chunks
|-- aggregate_iteration.py          # Aggregates completed case outputs
|-- aggregate_iteration.sh          # Slurm entry point for aggregation
|
|-- generate_simulation_pool.py     # Candidate-pool generation
|-- extract_reservoir_parameters.py # Input data extraction and preparation
|-- simulation_input_pool_500.csv   # Example candidate pool
|-- matched_reservoir_parameters.csv
|
|-- GP_ActiveLearning*.py           # Prototype active-learning experiments
|-- gp_active_learning_plots/       # Surrogate and acquisition diagnostics
|-- simulation_pool_diagnostics/    # Candidate-pool diagnostic plots
|
|-- manifests/                      # Generated simulation manifests
|-- cases/                          # Case-level simulation folders
|-- results/                        # Aggregated summary outputs
|-- logs/                           # Batch-job logs
|
|-- main.cc                         # Simulation executable entry point
|-- problem.hh                      # Simulation problem definition and JSON outputs
|-- params.input                    # Example simulation parameters
|-- vtk-merge-multi.py              # VTK post-processing utility
`-- CMakeLists.txt                  # Build definition for the simulation executable
```

## Methods

### Gaussian Process Regression

Gaussian Process Regression is a probabilistic surrogate modelling method. It learns a mapping from simulation inputs to outputs while also estimating uncertainty in regions where training data are sparse or inconsistent.

In this project, the GP surrogate is used to approximate expensive simulation outputs so the workflow can reason about which new simulations are likely to be most informative.

### Active Learning

Active learning means the next training points are selected intentionally instead of randomly or by a fixed grid. Here, the workflow evaluates a candidate pool and prioritizes cases where the surrogate model has high uncertainty or where a new sample would improve coverage.

### Candidate Selection

Candidate points are scored using model-based acquisition logic. The current practical workflow evaluates candidate cases across cycle numbers, aggregates their acquisition scores, and writes the highest-value candidates to a next-batch CSV.

### Why This Is Useful

For expensive simulations, the cost of a bad sampling strategy is high. Adaptive selection helps:

- avoid redundant simulation runs,
- improve surrogate quality with fewer cases,
- focus computation on uncertain or under-sampled regions,
- keep the simulation campaign organized and reproducible.

## Current Status

The project currently includes both prototype and practical workflow components:

- Prototype Gaussian Process active-learning experiments.
- A practical active-learning pipeline connected to real aggregated simulation outputs.
- Candidate-pool generation and diagnostic tooling.
- Batch simulation management through Python and Slurm scripts.
- Aggregation of case-level outputs into campaign-level CSV summaries.
- Initial coupling between surrogate training and future simulation selection.

The workflow is under active development. The current direction is to strengthen the coupling between the active-learning loop and the automated simulation submission pipeline.

## Example Usage

These commands illustrate the intended workflow. Paths and configuration values should be adjusted for the local build folder, HPC environment, and available simulation results.

### 1. Generate or Inspect a Candidate Pool

```bash
python generate_simulation_pool.py
```

Expected output examples:

- `simulation_input_pool_500.csv`
- diagnostic plots in `simulation_pool_diagnostics/`

### 2. Aggregate Completed Simulation Outputs

```bash
python aggregate_iteration.py \
  --manifest manifests/pool_051_150.csv \
  --iter-id pool_051_150
```

Expected output example:

- `results/pool_051_150_summary.csv`

### 3. Train the Surrogate and Select the Next Batch

```bash
python active_learning/active_learning_driver.py \
  --training-csv results/pool_051_150_summary.csv \
  --candidate-pool simulation_input_pool_500.csv \
  --target-prefix cycle \
  --max-cycles 10 \
  --batch-size 20 \
  --aggregate mean \
  --output manifests/active_learning_next.csv
```

Expected output example:

- `manifests/active_learning_next.csv`
- optional diagnostics JSON if `--diagnostics-json` is provided

### 4. Generate Simulation Manifests or Submit the Next Batch

Manifest-only mode:

```bash
python controller.py \
  --candidates-csv manifests/active_learning_next.csv \
  --batch-id al_round_01 \
  --jobs 4 \
  --cases-per-job 5 \
  --manifest-only
```

Submit mode:

```bash
python controller.py \
  --candidates-csv manifests/active_learning_next.csv \
  --batch-id al_round_01 \
  --jobs 4 \
  --cases-per-job 5
```

### 5. Run a Dry-Run Test of the Active-Learning Pipeline

```bash
python active_learning/test_dry_run.py
```

This creates synthetic data and verifies the core data-loading, training, scoring, and candidate-output path without requiring the simulation executable.

## Expected Outputs

| Output | Purpose |
| --- | --- |
| `training_dataset.csv` or long-format training data | Surrogate-model training table |
| `simulation_input_pool_500.csv` | Candidate pool for potential future simulations |
| `selected_next_batch.csv` or `active_learning_next.csv` | Recommended simulations for the next round |
| `manifests/*.csv` | Batch manifest consumed by the runner |
| `results/*_summary.csv` | Aggregated campaign-level outputs |
| `case_summary.json` | Per-case runner status and metadata |
| Diagnostic plots | Candidate-pool coverage, model behavior, and active-learning performance |

## Skills Demonstrated

This project demonstrates practical skills relevant to computational engineering, scientific software, data engineering, and applied machine learning:

- Python automation for simulation campaigns.
- Gaussian Process surrogate modelling.
- Active learning for expensive black-box functions.
- Uncertainty-aware candidate selection.
- Data ingestion, cleaning, reshaping, and aggregation.
- Candidate-pool generation and diagnostics.
- HPC/Slurm workflow orchestration.
- Batch execution, restart handling, and output management.
- Reproducible simulation campaign design.
- Integration of compiled simulation code with Python ML workflows.

## Future Development

Planned and natural extensions include:

- Tighter coupling between active-learning selection and simulation submission.
- More advanced acquisition functions that balance uncertainty, diversity, and objective value.
- Batch and asynchronous active learning for HPC environments.
- Improved feasibility filtering for failed or unstable simulation regions.
- Automated stopping criteria based on surrogate convergence.
- Richer diagnostics for model performance, uncertainty calibration, and input-space coverage.
- Scalable surrogate models for larger training sets and higher-dimensional candidate pools.

## Notes

The underlying simulation model is used here as an example of an expensive computational workflow. The same architecture can be adapted to other domains where simulations are costly, outputs must be aggregated carefully, and new experiments should be chosen adaptively.
