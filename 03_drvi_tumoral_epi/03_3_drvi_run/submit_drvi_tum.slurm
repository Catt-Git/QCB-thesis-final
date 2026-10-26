#!/usr/bin/env bash
#SBATCH --job-name=drvi_tum
#SBATCH --output=logs/drvi_tum_%j.out
#SBATCH --error=logs/drvi_tum_%j.err
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#
# 03_3 submit_drvi_tum: run_drvi_tum.py on the cluster (CPU, partition long, no GPU).
#
# Arguments after the script path reach run_drvi_tum.py unchanged. DRVI_ENV is the conda env
# with scvi-tools + drvi-py [default: benchmark-py-r].
#
# Inputs to copy up before submitting: 03_tum/shiao_tum_hvg_2k_nomt.h5ad, 03_tum/shiao_tum.h5ad
# (for the downstream object) and utils/cell_set.py (imported by the runner). To bring back:
# embed_drvi_tum_64_nomt.h5ad (03_4 on) and shiao_tum_drvi_tum_64_nomt.h5ad (03_5, 03_6).
#
# Usage:
#   cd 03_drvi_tumoral_epi/03_3_drvi_run && mkdir -p logs
#   export DATA_DIR=<data root>
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_drvi_tum.slurm
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_drvi_tum.slurm --overwrite

set -euo pipefail
mkdir -p logs

hostname; date
echo "Job name: ${SLURM_JOB_NAME:-}"; echo "Job ID: ${SLURM_JOB_ID:-}"
# Log the node, timestamp and job identifiers into the .out file for traceability.

: "${DATA_DIR:?set DATA_DIR (sbatch --export=ALL,DATA_DIR=... or export it first)}"
DRVI_ENV="${DRVI_ENV:-benchmark-py-r}"

# SLURM copies this file to the node's spool directory before running it, so
# "$0"/BASH_SOURCE point there and NOT at the repo - they cannot locate run_drvi_tum.py.
# SLURM_SUBMIT_DIR (where sbatch was called from) covers the documented usage above;
# DRVI_RUN_DIR overrides it.
SCRIPT_DIR="${DRVI_RUN_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
RUNNER="$SCRIPT_DIR/run_drvi_tum.py"
[ -f "$RUNNER" ] || { echo "run_drvi_tum.py not found: $RUNNER (set DRVI_RUN_DIR, or sbatch from 03_3_drvi_run/)" >&2; exit 1; }

# The runner imports utils/cell_set.py by relative path. Checked here so a
# half-copied repo fails in the first second rather than after the model is loaded.
CS_MODULE="$(dirname "$SCRIPT_DIR")/utils/cell_set.py"
[ -f "$CS_MODULE" ] || { echo "cell_set.py not found: $CS_MODULE (copy utils/ up as well)" >&2; exit 1; }

# torch reads this at import; keep it on the cores the scheduler granted rather than the
# node's full count.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

set +u
module load Miniconda3/202411
eval "$(/soft/system/software/Miniconda3/202411/bin/conda shell.bash hook)"
conda activate "$DRVI_ENV"
set -u

which python3
echo "env: $DRVI_ENV   OMP_NUM_THREADS=$OMP_NUM_THREADS   args: ${*:-<none>}"

if [ $# -gt 0 ]; then python3 "$RUNNER" "$@"; else python3 "$RUNNER"; fi

date
