#!/bin/bash
#SBATCH --job-name=fetch_metadata
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# submit_metadata.slurm
# SLURM wrapper for step 1: fetch_metadata.sh, which queries the ENA API for the whole BioProject PRJNA1032700.
# Seconds of work and a few hundred KB of output, so the allocation is minimal (1 core, 2G); it exists so that the entire phase 00 is reproducible as a chain of sbatch jobs rather than a mix of jobs and by-hand commands.
# No DATA_DIR: the output (ena_full_metadata.tsv) is a small TSV that stays IN the repo, next to fetch_metadata.sh.
#
# Usage:
#   cd 00_raw_data_processing/00_1_metadata   # submit FROM this folder
#   sbatch submit_metadata.slurm

# Exit on error, undefined variables, and pipe failures, so a failed ENA query is reported as FAILED by SLURM instead of leaving a truncated TSV behind an exit code 0.
set -euo pipefail

# Log the node, timestamp, and job identifiers into the .out file for traceability.
hostname; date
echo "Job name: $SLURM_JOB_NAME"; echo "Job ID: $SLURM_JOB_ID"

# Run the fetch. Submit with 'sbatch submit_metadata.slurm' FROM the 00_1_metadata folder: SLURM copies this file to a spool directory, so "$0" does not point at the repo and cannot be used to locate fetch_metadata.sh; SLURM_SUBMIT_DIR (the directory sbatch was called from) is the reliable anchor.
cd "$SLURM_SUBMIT_DIR"
bash fetch_metadata.sh

date
