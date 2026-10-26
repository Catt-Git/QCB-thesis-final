#!/bin/bash
#SBATCH --job-name=h5ad_concat
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# submit_h5ad_concat.slurm
# SLURM wrapper for step 6: runs build_requested_metadata.py (rebuilds the small sample_metadata_final.tsv from the sample map) and then build_combined_h5ad.py, which merges the 149 Cell Ranger .h5 outputs into a single AnnData.
# The two scripts are chained here rather than submitted separately because the metadata builder is a sub-second text transformation of a repo TSV: a job of its own would spend more time queueing than running, and keeping them together guarantees the metadata read by the merge is the one derived from the sample map actually present in this clone.
#
# Usage:
#   cd 00_raw_data_processing/00_6_build_h5ad   # submit FROM this folder, after Cell Ranger
#   export DATA_DIR=<data root>
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_h5ad_concat.slurm
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR --dependency=afterok:<cellranger_jobid> submit_h5ad_concat.slurm   # chain it to step 5

# Log the node, timestamp, and job identifiers into the .out file for traceability.
hostname; date
echo "Job name: $SLURM_JOB_NAME"; echo "Job ID: $SLURM_JOB_ID"

# Load Miniconda and activate the benchmark-py-r env (cellranger module NOT needed here, the merge only reads the .h5 files).
module load Miniconda3/202411
eval "$(/soft/system/software/Miniconda3/202411/bin/conda shell.bash hook)"
conda activate benchmark-py-r

# Sanity check: confirm python3 comes from the benchmark-py-r env.
which python3

# Site-specific data root, required: pass it with 'sbatch --export=ALL,DATA_DIR=...'.
: "${DATA_DIR:?set DATA_DIR (sbatch --export=ALL,DATA_DIR=... or export it first)}"

# Submit with 'sbatch submit_h5ad_concat.slurm' FROM the 00_6_build_h5ad folder: SLURM copies this file to a spool directory, so "$0" cannot be used to locate build_combined_h5ad.py; SLURM_SUBMIT_DIR (the directory sbatch was called from) is the reliable anchor.
cd "$SLURM_SUBMIT_DIR"

# Rebuild the per-sample metadata table (srr, sample_id, batch, cohort, treatment, fraction, dataset_origin, response) from the sample map. Both files are small TSVs kept IN the repo, so this touches no data under DATA_DIR; it runs here to keep the whole phase on the cluster and to make sure the table the merge is about to read matches the current sample map.
python3 build_requested_metadata.py \
    --sample-map ../00_2_sample_mapping/sample_map_gex_final.tsv \
    --out        sample_metadata_final.tsv

# Run the merge with python3 (NOT bash: bash would parse the Python syntax as shell commands and raise spurious errors, e.g. ImageMagick 'import').
# The Cell Ranger outputs, the sample list and the resulting .h5ad all live under DATA_DIR (outside the repo); only sample_metadata_final.tsv is a repo file, written next to this script by build_requested_metadata.py.
python3 build_combined_h5ad.py \
    --cellranger-dir "$DATA_DIR/cellranger_out" \
    --samples-file   "$DATA_DIR/renamed_fastq/unique_samples.txt" \
    --metadata-tsv   sample_metadata_final.tsv \
    --out            "$DATA_DIR/all_samples_combined.h5ad"
