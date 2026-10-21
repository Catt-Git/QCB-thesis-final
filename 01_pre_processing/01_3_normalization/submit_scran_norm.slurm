#!/bin/bash
#SBATCH --job-name=scran_norm
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=470G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Usage:
#   cd 01_pre_processing/01_3_normalization   # submit FROM this folder
#   sbatch submit_scran_norm.slurm
#   sbatch --export=ALL,DATA_DIR=<data root> submit_scran_norm.slurm   # another data root

set -euo pipefail

# Log the node, timestamp, and job identifiers into the .out file for traceability.
hostname; date
echo "Job name: ${SLURM_JOB_NAME:-}"; echo "Job ID: ${SLURM_JOB_ID:-}"

# Load Miniconda and activate benchmark-py-r (needs rpy2/anndata2ri + R scran/scib stack, see environments/benchmark-py-r.yml).
# set -u off while loading conda: its activation scripts use unset variables.
set +u
module load Miniconda3/202411
eval "$(/soft/system/software/Miniconda3/202411/bin/conda shell.bash hook)"
conda activate benchmark-py-r
set -u

# Data root of the thesis run, overridable with 'sbatch --export=ALL,DATA_DIR=...'.
# FIG_DIR is not set, so figures go to the repo's 01_pre_processing/figures/01_3_normalization.
export DATA_DIR="${DATA_DIR:-/users/genomics/albertoc/Tesi/hopes_and_dreams/datasets}"
# Threads limited to the cores granted by SLURM.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

which python3
python3 scran_norm.py
