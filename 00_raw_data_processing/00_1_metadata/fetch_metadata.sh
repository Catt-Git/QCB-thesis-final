#!/bin/bash

# fetch_metadata.sh
# Recover ENA metadata for the entire BioProject PRJNA1032700 (Shiao et al., Cancer Cell 2024; DOI 10.1016/j.ccell.2023.12.012).
# Useful to distinguish human GEX libraries from TCR/BCR and mouse samples, to reconstruct the SRR -> patient/timepoint/fraction (P/N) 
# association from 'library_name' and to record the expected number of FASTQ files per run (2 per run; 
# 0 for SRR26541168, which is empty on ENA and is dropped from the sample map at step 2).

# Provenance: ENA API filereport (public endpoint, result=read_run).
# OUTPUT: ena_full_metadata.tsv (one row per SRR, incl. header; fields: run_accession,sample_accession,sample_title,
# library_strategy,library_source,library_selection,library_name,experiment_title,fastq_md5,fastq_ftp)
# Note on fastq_md5: these MD5s refer to the files AS SERVED BY ENA. Since the pipeline downloads via prefetch + fasterq-dump
# (which rebuilds FASTQ from the .sra, not the ENA file), local MD5s will NOT match these. 
# Use fastq_md5 only as a source of truth for the expected file count per run.

# Usage:
#   bash fetch_metadata.sh   # no arguments/env vars; writes ena_full_metadata.tsv next to this script
#
# On the cluster (how the thesis run was done) this is launched through 'sbatch submit_metadata.slurm',
# which only needs a compute node with outbound internet access.

# Exit on error, undefined variables, and pipe failures, essential for detecting data corruption early in the pipeline.
set -euo pipefail

# BioProject accession and output filename for the metadata table. 
# The output is a small TSV kept IN the repo (no DATA_DIR needed) and written next to this script, so it can be run from any working directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT="PRJNA1032700"
OUT="$SCRIPT_DIR/ena_full_metadata.tsv"

# Query the ENA API for the requested fields (incl. fastq_md5) and save the full metadata table as TSV (limit=0 returns all runs). Retry up to 5 times, 10s apart, to absorb sporadic 400/5xx gateway responses from the ENA API.
wget --tries=5 --waitretry=10 --retry-connrefused -O "$OUT" \
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${PROJECT}&result=read_run&fields=run_accession,sample_accession,sample_title,library_strategy,library_source,library_selection,library_name,experiment_title,fastq_md5,fastq_ftp&format=tsv&limit=0"

# Report how many rows were retrieved, as a quick sanity check on the download (expected: 267 = 266 runs + header).
echo "Total rows (incl. header): $(wc -l < "$OUT")"
