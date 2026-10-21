#!/bin/bash

# check_download_completeness.sh
# Verify that EVERY expected human GEX run was actually downloaded (expected vs present).
# This is a completeness check (presence/count), complementary to verify_integrity.sh which checks file CONTENT. It catches runs that failed to download silently, before they surface later at Cell Ranger or as a short unique_samples.txt.
# For each SRR in the sample map it checks that BOTH <SRR>_1.fastq.gz and <SRR>_2.fastq.gz exist.

# INPUT: sample map (col 1 = SRR) and the directory holding the downloaded FASTQ files
# Usage:
#   export DATA_DIR=<data root>
#   bash check_download_completeness.sh                            # paths derived from the repo and DATA_DIR
#   bash check_download_completeness.sh <SAMPLE_MAP> <FASTQ_DIR>   # explicit paths, DATA_DIR not needed
#
# Exit code: 0 if all present, 1 if any run is missing (so it can be chained with &&).

# Exit on undefined variables and pipe failures; -e omitted so the loop reports ALL missing runs, not just the first.
set -uo pipefail

# Sample map (a small TSV that lives IN the repo, so it is resolved from this script's location and works from any cwd) and the FASTQ directory (data lives outside the repo, hence DATA_DIR). DATA_DIR is required only when the paths are not passed explicitly as arguments.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ $# -lt 2 ]; then
    : "${DATA_DIR:?not set. Export DATA_DIR=<data root, holding fastq_raw/>, or pass both paths as arguments. See README.}"
fi
SAMPLE_MAP="${1:-$SCRIPT_DIR/../00_2_sample_mapping/sample_map_gex_final.tsv}"
FASTQ_DIR="${2:-${DATA_DIR:-}/fastq_raw}"

# Counters and list of runs found to be incomplete.
expected=0
missing=0
missing_list=()

# Walk the expected SRR list and flag any run missing either FASTQ file.
while IFS=$'\t' read -r srr _rest; do
    [ -z "$srr" ] && continue
    expected=$((expected + 1))
    if [ ! -f "$FASTQ_DIR/${srr}_1.fastq.gz" ] || [ ! -f "$FASTQ_DIR/${srr}_2.fastq.gz" ]; then
        missing=$((missing + 1))
        missing_list+=("$srr")
    fi
done < "$SAMPLE_MAP"

# Report the expected/present/missing counts.
present=$((expected - missing))
echo "Expected runs (from $(basename "$SAMPLE_MAP")): $expected"
echo "Complete (both R1 and R2 present):            $present"
echo "Missing/incomplete:                           $missing"

# List the incomplete runs and fail (exit 1) if any; otherwise confirm success.
if [ "$missing" -gt 0 ]; then
    echo "Incomplete runs (delete whatever FASTQ is left of them, then launch rawprocessing_all.sh again: it downloads only the missing runs):"
    printf '  %s\n' "${missing_list[@]}"
    exit 1
fi
echo "All expected runs downloaded."