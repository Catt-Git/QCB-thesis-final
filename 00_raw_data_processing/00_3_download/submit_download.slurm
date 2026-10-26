#!/bin/bash
#SBATCH --job-name=download_fastq
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=1-149%4
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

# submit_download.slurm
# SLURM wrapper for the FASTQ download (~2.4 TB), which runs on the cluster and not locally.
# One array task = one line of sample_map_gex_final.tsv = one SRR, so the 149 runs are downloaded as 149 independent jobs instead of one multi-day serial loop: a run that fails (bad transfer, wall clock) takes down only its own task, and the failed indices can be resubmitted on their own with 'sbatch --export=ALL,DATA_DIR=$DATA_DIR --array=<idx>[,<idx>...] submit_download.slurm'.
# The '%4' throttle keeps at most 4 tasks downloading at once.

# Usage:
#   cd 00_raw_data_processing/00_3_download      # submit FROM this folder
#   export DATA_DIR=<data root>
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_download.slurm                 # all 149 runs, 4 at a time
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR --array=44,77 submit_download.slurm   # only those lines of the sample map (re-run of failed/corrupted runs)
#
# The line of a given SRR, i.e. its array index:
#   grep -n SRR26540980 ../00_2_sample_mapping/sample_map_gex_final.tsv | cut -d: -f1
# Logs (one pair per task) are written next to this script as download_fastq_<jobid>_<taskid>.{out,err}.

# Exit on error, undefined variables, and pipe failures, so a task that fails to download is reported as FAILED by SLURM rather than exiting 0 with nothing on disk.
set -euo pipefail

# Log the node, timestamp, and job identifiers into the .out file for traceability.
hostname; date
echo "Job name: $SLURM_JOB_NAME"; echo "Array job: $SLURM_ARRAY_JOB_ID  Task: $SLURM_ARRAY_TASK_ID"

# Data root (must hold ~2.4TB), required: pass it with 'sbatch --export=ALL,DATA_DIR=...'.
# DOWNLOAD_ENV is the conda environment providing the SRA Toolkit (prefetch/fasterq-dump): on the cluster that is benchmark-py-r, built from environments/benchmark-py-r.yml, which declares sra-tools.
: "${DATA_DIR:?set DATA_DIR (sbatch --export=ALL,DATA_DIR=... or export it first)}"
DOWNLOAD_ENV="${DOWNLOAD_ENV:-benchmark-py-r}"

# Load Miniconda and activate the environment holding the SRA Toolkit.
module load Miniconda3/202411
eval "$(/soft/system/software/Miniconda3/202411/bin/conda shell.bash hook)"
conda activate "$DOWNLOAD_ENV"

# Sanity check: confirm the SRA Toolkit comes from the activated environment before starting a multi-hour transfer
# (pigz is optional: download_fastq.sh falls back to gzip without it).
which prefetch fasterq-dump

# Submit with 'sbatch submit_download.slurm' FROM the 00_3_download folder: SLURM copies this file to a spool directory, so "$0" does not point at the repo and cannot be used to locate download_fastq.sh; SLURM_SUBMIT_DIR (the directory sbatch was called from) is the reliable anchor.
cd "$SLURM_SUBMIT_DIR"

# Map the array index to a run: line N of the sample map (the file has no header, one SRR per line in column 1). An index past the end of the file yields an empty SRR and fails immediately, instead of silently downloading nothing.
SAMPLE_MAP="${SAMPLE_MAP:-$SLURM_SUBMIT_DIR/../00_2_sample_mapping/sample_map_gex_final.tsv}"
SRR=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n {print $1}' "$SAMPLE_MAP")
[ -n "$SRR" ] || { echo "[error] no SRR on line $SLURM_ARRAY_TASK_ID of $SAMPLE_MAP" >&2; exit 1; }
echo ">>> [$SLURM_ARRAY_TASK_ID] $SRR"

# download_fastq.sh loops over a sample map, so the single run of this task is handed to it as a one-line map on a node-local temp file. This reuses the download logic verbatim (prefetch --max-size 100G, fasterq-dump --include-technical, compression) rather than duplicating it here, which is what keeps the local and cluster paths identical.
SHARD="$(mktemp)"
trap 'rm -f "$SHARD"' EXIT
printf '%s\n' "$SRR" > "$SHARD"

# Scratch for prefetch/fasterq-dump: the intermediate .sra and the uncompressed FASTQ of a single run are tens of GB, far more than a node's /tmp usually holds, so point the temp space at the data volume that was sized for this step.
export TMPDIR="$DATA_DIR/fastq_raw"
mkdir -p "$TMPDIR"

# Run the download for this SRR alone. THREADS follows the allocation (fasterq-dump uses it, and pigz when available), so changing --cpus-per-task above is enough to rescale.
DATA_DIR="$DATA_DIR" SAMPLE_MAP="$SHARD" THREADS="${SLURM_CPUS_PER_TASK:-8}" \
    bash download_fastq.sh

# Presence-and-size check so a task that produced nothing (or only R1) is marked FAILED and is easy to spot in 'sacct'. Content checks stay where they belong: verify_integrity.sh (gzip -t, read lengths, TSO) after the whole array has finished.
FASTQ_DIR="$DATA_DIR/fastq_raw"
for r in 1 2; do
    f="$FASTQ_DIR/${SRR}_${r}.fastq.gz"
    [ -s "$f" ] || { echo "[error] missing or empty after download: $f" >&2; exit 1; }
done
echo "[ok] $SRR -> $FASTQ_DIR/${SRR}_[12].fastq.gz"
date
