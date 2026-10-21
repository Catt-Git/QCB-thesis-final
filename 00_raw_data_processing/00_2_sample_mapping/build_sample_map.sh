#!/bin/bash

# build_sample_map.sh
# Filter human GEX libraries only (exclude TCR/BCR and mouse samples) and build the SRR -> patient / timepoint / sub-number / fraction map.
# Two DISTINCT naming formats: original library_name hNN[timepoint][subnum]_[P|N] (e.g. h19C_P, h03B2_P) 
# vs Cell Ranger sample_id PNN_[timepoint][subnum]_[P|N] (e.g. P19_C_P, P03_B2_P). 'P' = CD45+ (immune), 'N' = CD45- (non_immune). 
# Patient 03 is split into sub-portions A1/A2/B1/B2/C1/C2 (biopsy divided), but will be merged during the analysis.

# INPUT: ../00_1_metadata/ena_full_metadata.tsv
# OUTPUT: sample_map_gex_clean.tsv (all rows, incl. any UNPARSED) and sample_map_gex_final.tsv (correctly parsed rows
# with FASTQ on ENA only: every later step reads this one)

# Usage:
#   bash build_sample_map.sh
#   META=<path to ena_full_metadata.tsv> bash build_sample_map.sh   # override the default input
#
# On the cluster (how the thesis run was done) this is launched through 'sbatch submit_sample_map.slurm'.

# Exit on error, undefined variables, and pipe failures, essential for detecting data corruption or misconfiguration early in the pipeline.
set -euo pipefail

# Input metadata table and output filename for the cleaned sample map. Both are small TSVs that live IN the repo
# (no DATA_DIR needed) and are resolved from this script's location, so it can be run from any working directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
META="${META:-$SCRIPT_DIR/../00_1_metadata/ena_full_metadata.tsv}"
OUT="$SCRIPT_DIR/sample_map_gex_clean.tsv"

# Keep only run_accession (col 1) and library_name (col 7), then drop TCR/BCR libraries (_TCR suffix) 
# and mouse samples (library_name starting with 'm' + letter A-F).
tail -n +2 "$META" | awk -F'\t' '{print $1"\t"$7}' \
  | grep -vP '_TCR$' \
  | grep -vP '^SRR[0-9]+\tm[A-F]' \
  > "$SCRIPT_DIR/srr_gex_only.tsv"

# Parse each library_name into patient/timepoint/subnum/fraction and rebuild the Cell Ranger sample_id; 
# rows that do not match the expected pattern are flagged UNPARSED for manual review (in our case the 8 UNPARSED entries were mouse E0771 samples).
awk -F'\t' '{
    srr = $1
    origname = $2
    if (match(origname, /^h([0-9]+)([A-Z])([0-9]?)_(P|N)$/, arr)) {
        sample_name = "P" arr[1] "_" arr[2] arr[3] "_" arr[4]
        print srr"\t"arr[1]"\t"arr[2]"\t"arr[3]"\t"arr[4]"\t"sample_name
    } else {
        print srr"\tUNPARSED\tUNPARSED\tUNPARSED\tUNPARSED\tUNPARSED_"origname
    }
}' "$SCRIPT_DIR/srr_gex_only.tsv" > "$OUT"

# Surface any UNPARSED rows so they can be inspected before proceeding.
echo "UNPARSED samples (to check manually):"
grep UNPARSED "$OUT" || echo "  None."

# Runs with no FASTQ on ENA (empty fastq_ftp, column 10 of the metadata) have nothing to download: list them here
# and keep them out of the final map, so no later step ever sees them. In this BioProject that is SRR26541168 (P30_C_P).
EMPTY_RUNS="$(tail -n +2 "$META" | awk -F'\t' '$10 == "" {print $1}')"
echo "Runs with no FASTQ on ENA (excluded from the final map):"
awk -F'\t' -v empty="$EMPTY_RUNS" 'BEGIN { n = split(empty, a, "\n"); for (i = 1; i <= n; i++) e[a[i]] = 1 }
    ($1 in e) { print "  " $1 " (" $6 ")"; found = 1 } END { if (!found) print "  None." }' "$OUT"

# Write the final map, parsed rows with FASTQ only (columns: 1=srr 2=patient 3=timepoint 4=subnum 5=fraction 6=sample_id), and report the count.
grep -v UNPARSED "$OUT" \
  | awk -F'\t' -v empty="$EMPTY_RUNS" 'BEGIN { n = split(empty, a, "\n"); for (i = 1; i <= n; i++) e[a[i]] = 1 } !($1 in e)' \
  > "$SCRIPT_DIR/sample_map_gex_final.tsv"
echo "Valid human GEX samples: $(wc -l < "$SCRIPT_DIR/sample_map_gex_final.tsv")"
