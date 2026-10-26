#!/bin/bash
#SBATCH --job-name=smoke_metrics
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Usage:
#   cd 02_integration_benchmark/utils   # submit FROM this folder
#   export DATA_DIR=<data root>
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_smoke_test.slurm

set -euo pipefail

# Node, time and job ids for traceability.
hostname; date
echo "Job name: ${SLURM_JOB_NAME:-}"; echo "Job ID: ${SLURM_JOB_ID:-}"

# Same env as 01_3 (rpy2/anndata2ri + R).
# set -u off while loading conda: its activation scripts use unset variables.
set +u
module load Miniconda3/202411
eval "$(/soft/system/software/Miniconda3/202411/bin/conda shell.bash hook)"
conda activate benchmark-py-r
set -u

# DATA_DIR is required.
: "${DATA_DIR:?set DATA_DIR (sbatch --export=ALL,DATA_DIR=... or export it first)}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

# Copy the fixture first (scp smoke_test/smoke_fixture.h5ad <cluster>:<repo>/02_integration_benchmark/utils/smoke_test/);
# rebuilding it here needs far more memory and time.

which python3
cd "$SLURM_SUBMIT_DIR"
python3 smoke_test_metrics.py
