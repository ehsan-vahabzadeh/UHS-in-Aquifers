#!/bin/bash
#SBATCH --job-name=agg_iter
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

MANIFEST=$1
ITER_ID=$2

python3 aggregate_iteration.py \
    --manifest "$MANIFEST" \
    --iter-id "$ITER_ID"
