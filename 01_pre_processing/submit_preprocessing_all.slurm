#!/usr/bin/env bash
#SBATCH --job-name=preprocessing_all
#SBATCH --output=logs/preprocessing_all_%j.out
#SBATCH --error=logs/preprocessing_all_%j.err
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=470G
#
# 01 pre-processing: the whole 01_3 -> 01_5 chain as one job. Memory is sized on scran
# (619k cells, same 470G as submit_scran_norm.slurm).
#
# Only sets up the environment and calls preprocessing_all.sh, which holds the step list
# and the resume logic: resubmitting after a failure skips the finished steps.
#
# Usage:
#   cd 01_pre_processing                # submit from this folder (relative log paths)
#   export DATA_DIR=<data root>
#   bash preprocessing_all.sh --slurm   # usual route, through the runner
#   bash preprocessing_all.sh --slurm --force
#   bash preprocessing_all.sh --slurm reduce cluster
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_preprocessing_all.slurm   # direct use
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_preprocessing_all.slurm --force reduce cluster
#
# Variables (only DATA_DIR is required):
#   DATA_DIR     directory of the .h5ad chain
#   PREPROC_DIR  the repo's 01_pre_processing/      [default: SLURM_SUBMIT_DIR]
#   CC_GENES     cell cycle gene list               [default: $DATA_DIR/regev_lab_cell_cycle_genes.txt]
#   PREPROC_ENV  conda env                          [default: benchmark-py-r]
#   FIG_DIR      scran_norm.py figures              [default: figures/01_3_normalization]

set -euo pipefail

hostname; date
echo "Job name: ${SLURM_JOB_NAME:-}"; echo "Job ID: ${SLURM_JOB_ID:-}"

: "${DATA_DIR:?set DATA_DIR (sbatch --export=ALL,DATA_DIR=... or use ./preprocessing_all.sh --slurm)}"

# The job runs from a SLURM spool copy, so the repo comes from PREPROC_DIR (set by the
# runner) or from the folder sbatch was called from.
SCRIPT_DIR="${PREPROC_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
RUNNER="$SCRIPT_DIR/preprocessing_all.sh"
[ -x "$RUNNER" ] || { echo "runner not found or not executable: $RUNNER (set PREPROC_DIR, or sbatch from 01_pre_processing/)" >&2; exit 1; }

PREPROC_ENV="${PREPROC_ENV:-benchmark-py-r}"
export CC_GENES="${CC_GENES:-$DATA_DIR/regev_lab_cell_cycle_genes.txt}"
# Threads limited to the cores granted by SLURM.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

# set -u off while loading conda: its activation scripts use unset variables.
set +u
module load Miniconda3/202411
eval "$(/soft/system/software/Miniconda3/202411/bin/conda shell.bash hook)"
conda activate "$PREPROC_ENV"
set -u

which python3
echo "env: $PREPROC_ENV   OMP_NUM_THREADS=$OMP_NUM_THREADS   runner args: ${*:-<none>}"

# No arguments -> call the runner with none (an empty "$@" fails under set -u on old bash).
if [ $# -gt 0 ]; then "$RUNNER" "$@"; else "$RUNNER"; fi

date
