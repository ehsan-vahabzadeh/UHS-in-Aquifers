#!/bin/bash
#SBATCH --job-name=agg_iter
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

RUN_ROOT="${DUMUX_RUN_ROOT:-$PWD}"
SCRIPT_DIR="${DUMUX_SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

cd "$RUN_ROOT"
mkdir -p logs results

MANIFEST=$1
ITER_ID=$2

echo "=== aggregate_iteration.sh startup ==="
date
hostname
pwd
echo "RUN_ROOT=$RUN_ROOT"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-}"
echo "MANIFEST=$MANIFEST"
echo "ITER_ID=$ITER_ID"
echo "SCRIPT_DIR=$SCRIPT_DIR"
echo "python3=$(command -v python3)"

if [[ ! -f "$SCRIPT_DIR/aggregate_iteration.py" ]]; then
    echo "ERROR: cannot find aggregate_iteration.py at: $SCRIPT_DIR/aggregate_iteration.py" >&2
    echo "DUMUX_RUN_ROOT=$RUN_ROOT" >&2
    echo "DUMUX_SCRIPT_DIR=$SCRIPT_DIR" >&2
    exit 2
fi

python3 -u "$SCRIPT_DIR/aggregate_iteration.py" \
    --manifest "$MANIFEST" \
    --iter-id "$ITER_ID"
