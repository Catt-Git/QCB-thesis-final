#!/usr/bin/env bash
#SBATCH --job-name=drvi
#SBATCH --output=logs/drvi_%j.out
#SBATCH --error=logs/drvi_%j.err
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#
# 02_2 integration: DRVI on the cluster (CPU, GPU was not present in my case), a single job wrapping run_drvi.py
# once the latent size has been chosen in the notebook. `long` partition: CPU
# training is slow.
#
#   DRVI_ENV   conda env with scvi-tools + drvi-py   [default: benchmark-py-r]
#
# Usage:
#   cd 02_integration_benchmark/02_2_integration && mkdir -p logs   # submit FROM this folder
#   export DATA_DIR=<data root>
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_drvi.slurm   # n_latent 128
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_drvi.slurm --n-latent 64   # another size
#
# Arguments after the script path go to run_drvi.py. Input upload, if needed:
#   utils/sync_to_cluster.sh --what inputs --method drvi

set -euo pipefail
mkdir -p logs

# Node, time and job ids for traceability.
hostname; date
echo "Job name: ${SLURM_JOB_NAME:-}"; echo "Job ID: ${SLURM_JOB_ID:-}"

: "${DATA_DIR:?set DATA_DIR (sbatch --export=ALL,DATA_DIR=... or export it first)}"
DRVI_ENV="${DRVI_ENV:-benchmark-py-r}"

# SLURM runs a spool copy, so locate the repo via INTEGRATION_DIR or the submit dir.
SCRIPT_DIR="${INTEGRATION_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
RUNNER="$SCRIPT_DIR/run_drvi.py"
[ -f "$RUNNER" ] || { echo "run_drvi.py not found: $RUNNER (set INTEGRATION_DIR, or sbatch from 02_2_integration/)" >&2; exit 1; }

# Limit torch to the granted cores.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

# set -u off for the conda block: base's openjdk activate script reads unset JAVA_HOME.
set +u
module load Miniconda3/202411
eval "$(/soft/system/software/Miniconda3/202411/bin/conda shell.bash hook)"
conda activate "$DRVI_ENV"
set -u

which python3
echo "env: $DRVI_ENV   OMP_NUM_THREADS=$OMP_NUM_THREADS   args: ${*:-<none>}"

# Explicit branch: an empty "$@" trips set -u on older bash.
if [ $# -gt 0 ]; then python3 "$RUNNER" "$@"; else python3 "$RUNNER"; fi

date
