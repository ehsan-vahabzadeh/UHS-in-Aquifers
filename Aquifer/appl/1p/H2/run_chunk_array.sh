#!/bin/bash
#SBATCH --job-name=dumux_chunk
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

MANIFEST=$1
ITER_ID=$2
EXECUTABLE=$3

python3 scripts/run_chunk.py \
    --manifest "$MANIFEST" \
    --iter-id "$ITER_ID" \
    --chunk-id "$SLURM_ARRAY_TASK_ID" \
    --ntasks "$SLURM_NTASKS" \
    --executable "$EXECUTABLE"