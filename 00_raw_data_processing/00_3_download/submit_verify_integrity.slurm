#!/bin/bash
#SBATCH --job-name=verify_integrity
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# submit_verify_integrity.slurm
# SLURM wrapper for step 3b: verify_integrity.sh, the content check that runs once the whole download array has finished.
# Sized for I/O, not CPU: 'gzip -t' decompresses all 298 files (~2.4 TB) serially, which is hours of streaming at a few hundred MB/s, hence the 'long' partition. 2 cores and 8G are plenty (gzip is single-threaded and streams).
# Run it AFTER the array is complete, not per task: it is a full-dataset audit, and the download tasks are checked individually only for file presence.
#
# Usage:
#   cd 00_raw_data_processing/00_3_download   # submit FROM this folder
#   export DATA_DIR=<data root>
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_verify_integrity.slurm
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR --dependency=afterok:<download_array_jobid> submit_verify_integrity.slurm   # chain it to the download array: starts when every task is done

# Exit on undefined variables and pipe failures; -e is intentionally omitted here as in verify_integrity.sh itself, so a single corrupted file does not abort the audit before the remaining files are checked.
set -uo pipefail

# Log the node, timestamp, and job identifiers into the .out file for traceability.
hostname; date
echo "Job name: $SLURM_JOB_NAME"; echo "Job ID: $SLURM_JOB_ID"

# Data root, required: pass it with 'sbatch --export=ALL,DATA_DIR=...'.
: "${DATA_DIR:?set DATA_DIR (sbatch --export=ALL,DATA_DIR=... or export it first)}"

# Submit with 'sbatch submit_verify_integrity.slurm' FROM the 00_3_download folder: SLURM copies this file to a spool directory, so "$0" cannot be used to locate verify_integrity.sh; SLURM_SUBMIT_DIR (the directory sbatch was called from) is the reliable anchor.
cd "$SLURM_SUBMIT_DIR"

# Run the three checks (gzip integrity, R1/R2 read lengths, 5' TSO on R2) and keep a copy of the output next to the log files.
REPORT="$SLURM_SUBMIT_DIR/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.report"
bash verify_integrity.sh | tee "$REPORT"

# verify_integrity.sh reports corruption on stdout but always exits 0, so the wrapper turns 'CORRUPTED' lines into a non-zero exit: that is what makes the job show up as FAILED in 'sacct' and what lets a --dependency=afterok chain (rename, Cell Ranger) refuse to start on a broken dataset.
corrupted=$(grep -c 'CORRUPTED' "$REPORT")
if [ "$corrupted" -gt 0 ]; then
    echo "[error] $corrupted corrupted file(s), listed in $REPORT: delete both FASTQ of each of those runs and launch rawprocessing_all.sh again, which downloads only the missing runs (README, Repairing a single run)." >&2
    exit 1
fi
echo "[ok] no corrupted file; report: $REPORT"
date
