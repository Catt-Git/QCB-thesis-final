#!/bin/bash

# verify_integrity.sh
# Post-download checks on the FASTQ files.
# 1) gzip integrity of every FASTQ.
# 2) R1 and R2 read lengths, to confirm the expected structure and to establish which file is the barcode
# (R1) and which is the cDNA (R2): R1 (barcode+UMI+padding) is 101bp for most samples and 151bp for 2 samples 
# (16bp barcode + 10bp UMI = 26bp effective, handled with --r1-length=26); R2 (cDNA) is the long biological 
# read and starts with the 5'v1 TSO AAGCAGTGGTATCAAC.

# Usage:
#   export DATA_DIR=<data root>
#   bash verify_integrity.sh                      # checks $DATA_DIR/fastq_raw
#   FASTQ_DIR=<folder> bash verify_integrity.sh   # checks a different folder
#
# On the cluster (how the thesis run was done) this is launched through 'sbatch submit_verify_integrity.slurm',
# which sets DATA_DIR, keeps a copy of the report and turns any CORRUPTED line into a non-zero exit
# (this script always exits 0 by design, see below).

# Exit on undefined variables and pipe failures;
# note -e is intentionally omitted so a single corrupted file does not abort the whole check loop.
set -uo pipefail

# Move into the folder holding the FASTQ written by download_fastq.sh, so the globs below match the data rather than whatever the cwd happened to be.
: "${DATA_DIR:?not set. Export DATA_DIR=<data root>; the FASTQ are expected in \$DATA_DIR/fastq_raw. See README.}"
FASTQ_DIR="${FASTQ_DIR:-$DATA_DIR/fastq_raw}"
cd "$FASTQ_DIR"

# Test each gzip archive; report OK or flag the file for re-download.
echo "### 1) gzip integrity ###"
for f in *.fastq.gz; do
    echo -n "$f: "
    gzip -t "$f" && echo "OK" || echo "CORRUPTED - delete both FASTQ of this run and launch rawprocessing_all.sh again (README, Repairing a single run)"
done

# Print the unique sequence lengths per file to confirm R1 (26/101/151bp) vs R2 (long cDNA).
echo
echo "### 2) Read lengths (first 100 reads per file) ###"
for f in *_1.fastq.gz *_2.fastq.gz; do
    [ -e "$f" ] || continue
    echo -n "$f  len: "
    zcat "$f" | head -400 | awk 'NR%4==2{print length($0)}' | sort -u | tr '\n' ','
    echo
done

# Check the first 16bp of R2 against the expected 10x Chromium 5'v1 TSO to confirm the chemistry.
echo
echo "### 3) Confirm 5'v1 TSO at the start of R2 (expected: AAGCAGTGGTATCAAC) ###"
for f in *_2.fastq.gz; do
    [ -e "$f" ] || continue
    echo -n "$f  R2 first prefix: "
    zcat "$f" | awk 'NR%4==2{print substr($0,1,16); exit}'
done
