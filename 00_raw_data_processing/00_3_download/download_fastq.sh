#!/bin/bash

# download_fastq.sh
# Download the FASTQ files of the human GEX libraries via the SRA Toolkit.
# Total file size is about 2.4 TB! :(
# Why not direct ENA (wget/aria2c): the data were sequenced 151x151bp (non-standard) and released recently; ENA had not processed 
# the runs correctly and EXCLUDED the barcode read (R1) as a "technical read", so bulk ENA downloads were missing the barcode (or corrupted/truncated). 
# The reliable method is prefetch + fasterq-dump with --include-technical.

# INPUT: ../00_2_sample_mapping/sample_map_gex_final.tsv (col 1 = SRR)
# OUTPUT: <SRR>_1.fastq.gz (barcode+UMI, 26bp effective) and <SRR>_2.fastq.gz (cDNA) in $DATA_DIR/fastq_raw
# REQUIREMENTS: sra-tools (prefetch, fasterq-dump) and ~2.4TB free under $DATA_DIR; pigz is used when available (faster), gzip otherwise
# Usage:
#   export DATA_DIR=<data root>   # needs room for ~2.4TB
#   bash download_fastq.sh
#
# On the cluster (how the thesis run was done) this is not launched directly but through 'sbatch submit_download.slurm',
# which turns the loop below into a SLURM array of one task per SRR and calls this script with SAMPLE_MAP pointing at a
# one-line map. Running it by hand still works and downloads the whole list serially.

# Exit on error, undefined variables, and pipe failures, essential for detecting data corruption or misconfiguration early in the pipeline.
set -euo pipefail

# Sample map providing the SRR list (a small TSV in the repo, resolved from this script's location so any cwd works), destination for the FASTQ, and number of threads for fasterq-dump/pigz.
# DATA_DIR is required, not defaulted: ~2.4TB must not land inside the repo clone by accident, and the right location is site-specific.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
: "${DATA_DIR:?not set. Export DATA_DIR=<data root with room for ~2.4TB>; the FASTQ are written to \$DATA_DIR/fastq_raw. See README.}"
SAMPLE_MAP="${SAMPLE_MAP:-$SCRIPT_DIR/../00_2_sample_mapping/sample_map_gex_final.tsv}"
FASTQ_DIR="${FASTQ_DIR:-$DATA_DIR/fastq_raw}"
THREADS="${THREADS:-8}"

# Compressor: pigz (parallel) when the node has it, plain gzip otherwise. Both write the same .fastq.gz format;
# pigz only makes the compression faster, so it is not a requirement.
if command -v pigz >/dev/null 2>&1; then GZIP_CMD=(pigz -p "$THREADS"); else GZIP_CMD=(gzip); fi
echo "compressor: ${GZIP_CMD[0]}"

# Work inside the data folder: prefetch writes its .sra cache into the cwd and fasterq-dump then reads it from there, so the download runs where the space is rather than wherever the script was launched from.
mkdir -p "$FASTQ_DIR"
cd "$FASTQ_DIR"

# For each SRR: skip if already downloaded; prefetch with --max-size 100G (the 20G default skips the 37-40GB files); extract with --include-technical (critical: forces the barcode read R1 to be written, otherwise dropped) and --split-files; then compress.
cut -f1 "$SAMPLE_MAP" | while read -r srr; do
    [ -z "$srr" ] && continue
    if [ -f "${srr}_2.fastq.gz" ]; then
        echo "[SKIP] already present: $srr"
        continue
    fi
    echo "=== Downloading $srr ==="
    prefetch --max-size 100G "$srr"
    fasterq-dump --include-technical --split-files -p -e "$THREADS" -O . "$srr"
    "${GZIP_CMD[@]}" "${srr}_1.fastq" "${srr}_2.fastq"
    echo "=== $srr done ==="
done

# Remind the user to validate the downloaded files before renaming.
# There should be 2 FASTQ files per SRR: 149 runs in the map, so 298 FASTQ files in total.
echo "Download complete. Run verify_integrity.sh for the checks."
