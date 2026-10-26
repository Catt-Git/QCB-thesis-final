#!/bin/bash
#SBATCH --job-name=cell_ranger
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=90
#SBATCH --mem=470G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# submit_cellranger.slurm
# SLURM wrapper for step 5: allocates the largest CPU-only node (90 cores, 470G) and launches the Cell Ranger batch (run_cellranger_batch.sh, 5 samples x 16 cores in parallel). No --gres=gpu (the cluster is CPU-only).
#
# Usage:
#   cd 00_raw_data_processing/00_5_cellranger   # submit FROM this folder, after the renaming
#   export DATA_DIR=<data root>
#   export REF=<path to refdata-cellranger-GRCh38-3.0.0>
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR,REF=$REF submit_cellranger.slurm
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR,REF=$REF --dependency=afterok:<rename_jobid> submit_cellranger.slurm   # chain it to step 4

# Log the node, timestamp, and job identifiers into the .out file for traceability.
hostname; date
echo "Job name: $SLURM_JOB_NAME"; echo "Job ID: $SLURM_JOB_ID"

# Load Miniconda and load the Cell Ranger 4.0.0 module.
module load Miniconda3/202411
eval "$(/soft/system/software/Miniconda3/202411/bin/conda shell.bash hook)"
module load cellranger/4.0.0

# Site-specific paths, required: DATA_DIR is the data root (expects renamed_fastq/ inside it, writes cellranger_out/ there), REF is the reference transcriptome. Pass them with 'sbatch --export=ALL,DATA_DIR=...,REF=...'.
: "${DATA_DIR:?set DATA_DIR (sbatch --export=ALL,DATA_DIR=...,REF=...)}"
: "${REF:?set REF to refdata-cellranger-GRCh38-3.0.0 (sbatch --export=ALL,DATA_DIR=...,REF=...)}"

# Run the parallel Cell Ranger batch within this allocation. Submit with 'sbatch submit_cellranger.slurm' FROM the 00_5_cellranger folder: SLURM copies this file to a spool directory, so "$0" does not point at the repo and cannot be used to locate run_cellranger_batch.sh; SLURM_SUBMIT_DIR (the directory sbatch was called from) is the reliable anchor.
cd "$SLURM_SUBMIT_DIR"
bash run_cellranger_batch.sh

# run_cellranger_batch.sh exits 0 even when a sample fails (it logs it and moves on), so check here that every sample has
# its filtered matrix: a missing one fails the job, and the h5ad step chained with afterok does not start on an incomplete set.
bash ../utils/check_cellranger_completeness.sh || exit 1
