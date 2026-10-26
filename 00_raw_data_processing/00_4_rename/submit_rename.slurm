#!/bin/bash
#SBATCH --job-name=rename_fastq
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# submit_rename.slurm
# SLURM wrapper for step 4: rename_for_cellranger.sh, which copies the SRR FASTQ into the naming scheme Cell Ranger requires (<sample>_S1_L001_R[12]_001.fastq.gz) and writes unique_samples.txt.
# Pure I/O, single-threaded: rsync copies ~2.4 TB, so this is hours on the storage and needs the 'long' partition, but almost no CPU or RAM. It is a copy, not a move (the originals stay in fastq_raw), so the data root must hold ~4.8 TB at this point.
#
# Usage:
#   cd 00_raw_data_processing/00_4_rename   # submit FROM this folder, after the integrity check
#   export DATA_DIR=<data root>
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_rename.slurm
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR --dependency=afterok:<verify_jobid> submit_rename.slurm   # chain it to step 3b

# Exit on error, undefined variables, and pipe failures: a copy that fails halfway must be reported as FAILED, since a missing R1/R2 pair would only surface much later as a Cell Ranger error.
set -euo pipefail

# Log the node, timestamp, and job identifiers into the .out file for traceability.
hostname; date
echo "Job name: $SLURM_JOB_NAME"; echo "Job ID: $SLURM_JOB_ID"

# Data root, required: pass it with 'sbatch --export=ALL,DATA_DIR=...'. rename_for_cellranger.sh reads $DATA_DIR/fastq_raw and writes $DATA_DIR/renamed_fastq.
: "${DATA_DIR:?set DATA_DIR (sbatch --export=ALL,DATA_DIR=... or export it first)}"

# Run the renaming copy. Submit with 'sbatch submit_rename.slurm' FROM the 00_4_rename folder: SLURM copies this file to a spool directory, so "$0" cannot be used to locate rename_for_cellranger.sh; SLURM_SUBMIT_DIR (the directory sbatch was called from) is the reliable anchor.
cd "$SLURM_SUBMIT_DIR"
bash rename_for_cellranger.sh

# Report the sample count that Cell Ranger will iterate over (expected: 149), so the number is in the log before step 5 is submitted.
echo "[ok] renamed FASTQ in $DATA_DIR/renamed_fastq ($(wc -l < "$DATA_DIR/renamed_fastq/unique_samples.txt") samples)"
date
