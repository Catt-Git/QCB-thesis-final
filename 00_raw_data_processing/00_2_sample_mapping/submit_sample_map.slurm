#!/bin/bash
#SBATCH --job-name=sample_map
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# submit_sample_map.slurm
# SLURM wrapper for step 2: build_sample_map.sh, which turns ena_full_metadata.tsv into the SRR -> patient/timepoint/fraction map.
# Pure text processing on a 266-row table.
# No DATA_DIR: input and outputs are small TSVs that stay IN the repo, next to build_sample_map.sh.

# Usage:
#   cd 00_raw_data_processing/00_2_sample_mapping   # submit FROM this folder, after step 1
#   sbatch submit_sample_map.slurm

# Exit on error, undefined variables, and pipe failures, so a parsing failure is reported as FAILED rather than producing a silently truncated map.
set -euo pipefail

# Log the node, timestamp, and job identifiers into the .out file for traceability.
hostname; date
echo "Job name: $SLURM_JOB_NAME"; echo "Job ID: $SLURM_JOB_ID"

awk --version 2>/dev/null | head -1 | grep -qi 'gnu awk' \
  || { echo "[error] awk on this node is not GNU awk; build_sample_map.sh needs the 3-argument match(str, re, arr), a gawk extension. Load a gawk module or activate an environment providing it." >&2; exit 1; }

# Build the map. Submit with 'sbatch submit_sample_map.slurm' FROM the 00_2_sample_mapping folder: SLURM copies this file to a spool directory, so "$0" cannot be used to locate build_sample_map.sh; SLURM_SUBMIT_DIR (the directory sbatch was called from) is the reliable anchor.
cd "$SLURM_SUBMIT_DIR"
bash build_sample_map.sh

date
