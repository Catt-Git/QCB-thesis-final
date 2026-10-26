#!/usr/bin/env bash
#SBATCH --job-name=03_driver
#SBATCH --partition=long
#SBATCH --nodes=1
#
# 03 submit_driver: one phase-03 driver (infercnv_all.sh, subsetting_all.sh, ...) as one job.
# Job name, CPUs, memory and logs are set by the driver's --slurm (utils/slurm.sh); this file
# only sets up conda and runs the driver without --slurm, so the driver's own resume logic
# decides what still has to run.
#
# Usage:
#   export DATA_DIR=<data root>
#   bash <driver>.sh --slurm [driver args]      # usual route
#   sbatch --cpus-per-task=8 --mem=64G --export=ALL,DATA_DIR=$DATA_DIR \
#     utils/submit_driver.slurm <absolute path of driver>.sh [driver args]   # direct use
#
#   JOB_ENV   conda env activated before the driver [default: benchmark-py-r]

set -euo pipefail

hostname; date
echo "Job name: ${SLURM_JOB_NAME:-}"; echo "Job ID: ${SLURM_JOB_ID:-}"

: "${DATA_DIR:?set DATA_DIR (sbatch --export=ALL,DATA_DIR=... or use the driver's --slurm)}"
DRIVER="${1:?usage: sbatch submit_driver.slurm <driver.sh> [driver args]}"
shift
[ -f "$DRIVER" ] || { echo "driver not found: $DRIVER (pass its absolute path)" >&2; exit 1; }

JOB_ENV="${JOB_ENV:-benchmark-py-r}"
# Threads limited to the cores granted by SLURM.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
# Unbuffered, so an OOM kill does not lose the progress log; headless node.
export PYTHONUNBUFFERED=1
export MPLBACKEND=Agg

# set -u off for the conda block: base's openjdk activate script reads unset JAVA_HOME.
set +u
module load Miniconda3/202411
eval "$(/soft/system/software/Miniconda3/202411/bin/conda shell.bash hook)"
conda activate "$JOB_ENV"
set -u

which python3
echo "env: $JOB_ENV   OMP_NUM_THREADS=$OMP_NUM_THREADS   driver: $DRIVER   args: ${*:-<none>}"

# No arguments -> call the driver with none (an empty "$@" fails under set -u on old bash).
if [ $# -gt 0 ]; then bash "$DRIVER" "$@"; else bash "$DRIVER"; fi

date
