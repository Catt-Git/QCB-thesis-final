#!/bin/bash

# run_cellranger_batch.sh
# Run 'cellranger count' on all samples in parallel (5 jobs x 16 cores = 80), inside a single SLURM allocation (see submit_cellranger.slurm).
# Key parameters (verbatim, for Methods): reference refdata-cellranger-GRCh38-3.0.0 (same as the paper); --r1-length=26 (R1 is 101/151bp = 16bp barcode + 10bp UMI + padding; without it Cell Ranger cannot detect the barcode);
#--chemistry=SC5P-R2 (10x Chromium 5'v1, confirmed by the TSO AAGCAGTGGTATCAAC at the start of R2); --localcores=16 --localmem=64; Cell Ranger 4.0.0 (via 'module load cellranger/4.0.0').

# INPUT: renamed_fastq/ and renamed_fastq/unique_samples.txt
# OUTPUT: cellranger_out/<sample>/outs/filtered_feature_bc_matrix.h5

# Usage:
#   export DATA_DIR=<data root>
#   export REF=<path to refdata-cellranger-GRCh38-3.0.0>
#   bash run_cellranger_batch.sh
#
# Normally launched via 'sbatch submit_cellranger.slurm', not directly.

# Exit on undefined variables and pipe failures; note -e is intentionally omitted so one failing sample does not abort the whole batch (failures are logged instead).
set -uo pipefail

# Data root, input FASTQ folder, reference transcriptome, output directory, and number of samples to process concurrently (5 x 16 cores = 80 on a 90-core node).
# DATA_DIR and REF are deliberately REQUIRED rather than defaulted: the ~2.4TB of FASTQ live outside the repo and their location is site-specific. Export them before running (or pass them to submit_cellranger.slurm with sbatch --export).
: "${DATA_DIR:?not set. Export DATA_DIR=<data root> and REF=<path to refdata-cellranger-GRCh38-3.0.0>. See README.}"
: "${REF:?not set. Export REF=<path to refdata-cellranger-GRCh38-3.0.0>. See README.}"
FASTQ_DIR="${FASTQ_DIR:-$DATA_DIR/renamed_fastq}"
OUTDIR="${OUTDIR:-$DATA_DIR/cellranger_out}"
PARALLEL_JOBS="${PARALLEL_JOBS:-5}"

# Create and move into the output directory (Cell Ranger writes into the current working directory).
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Define the per-sample worker: skip a sample whose filtered matrix already exists, so re-running the batch only processes
# what is missing; otherwise run cellranger count with the fixed parameters and log the sample name to failed_samples.txt
# on non-zero exit. Note: Cell Ranger refuses to restart with an existing --id, so delete folders of failed runs before relaunching.
run_one_sample() {
    sample="$1"
    if [ -f "$OUTDIR/$sample/outs/filtered_feature_bc_matrix.h5" ]; then
        echo "=== [have] $sample: filtered matrix already there, skipping ==="
        return 0
    fi
    echo "=== Starting $sample ==="
    cellranger count --id="${sample}" \
      --transcriptome="$REF" \
      --fastqs="$FASTQ_DIR" \
      --sample="${sample}" \
      --r1-length=26 \
      --chemistry=SC5P-R2 \
      --localcores=16 --localmem=64

    if [ $? -ne 0 ]; then
        echo "$sample" >> "$OUTDIR/failed_samples.txt"
    fi
}

# failed_samples.txt lists the failures of THIS batch: start it empty, so failures of an earlier run do not linger.
rm -f "$OUTDIR/failed_samples.txt"

# Export the function and variables so they are visible to the subshells spawned by xargs.
export -f run_one_sample
export FASTQ_DIR REF OUTDIR

# Feed the sample list to xargs and run PARALLEL_JOBS workers concurrently, one per sample.
cat "$FASTQ_DIR/unique_samples.txt" \
  | xargs -P "$PARALLEL_JOBS" -I{} bash -c 'run_one_sample "$@"' _ {}

# Point the user to the failure log after the batch completes.
echo "Done. Check failed_samples.txt for any errors."
