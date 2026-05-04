#!/bin/bash
#SBATCH --job-name=dumux_chunk
#SBATCH -p multicore
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1

#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

RUN_ROOT="${DUMUX_RUN_ROOT:-$PWD}"
MPI_RUNNER="${DUMUX_MPI_RUNNER:-mpirun}"
if [[ "$MPI_RUNNER" == */* ]]; then
    export PATH="$(dirname "$MPI_RUNNER"):$PATH"
fi

# Keep the batch job independent of fragile module-name spelling. These are
# the library paths reported by ldd/module list on CSF for this executable.
export LD_LIBRARY_PATH="/opt/software/RI/apps/GMP/6.2.1-GCCcore-12.3.0/lib:${LD_LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="/opt/software/RI/apps/OpenMPI/4.1.8-GCC-13.2.0/lib:${LD_LIBRARY_PATH:-}"

cd "$RUN_ROOT"
mkdir -p logs cases results

MANIFEST=$1
ITER_ID=$2
EXECUTABLE=$3

echo "=== run_chunk_array.sh startup ==="
date
hostname
pwd
echo "RUN_ROOT=$RUN_ROOT"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-}"
echo "SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID:-}"
echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-}"
echo "SLURM_NTASKS=${SLURM_NTASKS:-}"
echo "MANIFEST=$MANIFEST"
echo "ITER_ID=$ITER_ID"
echo "EXECUTABLE=$EXECUTABLE"
echo "MPI_RUNNER=$MPI_RUNNER"
echo "python3=$(command -v python3)"
echo "mpirun_in_path=$(command -v mpirun || true)"
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

if [[ ! -f "$RUN_ROOT/run_chunk.py" ]]; then
    echo "ERROR: cannot find run_chunk.py at: $RUN_ROOT/run_chunk.py" >&2
    echo "DUMUX_RUN_ROOT=$RUN_ROOT" >&2
    exit 2
fi

python3 -u "$RUN_ROOT/run_chunk.py" \
    --manifest "$MANIFEST" \
    --iter-id "$ITER_ID" \
    --chunk-id "$SLURM_ARRAY_TASK_ID" \
    --ntasks "$SLURM_NTASKS" \
    --executable "$EXECUTABLE"
