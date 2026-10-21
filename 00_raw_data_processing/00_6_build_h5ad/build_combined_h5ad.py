#!/usr/bin/env python3

# build_combined_h5ad.py
# Load filtered_feature_bc_matrix.h5 for each sample in unique_samples.txt, attach the per-sample metadata (batch, cohort, treatment, fraction, dataset_origin, response) from sample_metadata_final.tsv, and combine everything into a single AnnData.

# Object contract (for Methods): 
# - input filtered_feature_bc_matrix.h5 (Cell Ranger 4.0.0 cell-calling already applied), so counts differ from the paper BY DESIGN
# - var_names = gene symbols made unique with var_names_make_unique() (CellTypist Cells_Adult_Breast.pkl works on symbols)
# - Ensembl IDs kept in .var['gene_ids'] via merge='same'
# - raw counts stored in .layers['counts'] before any downstream normalisation
# - barcodes made unique per sample as f"{sample}_{bc}", so concat uses index_unique=None.

# Usage:
#   export DATA_DIR=<data root>
#   python3 build_combined_h5ad.py \
#       --cellranger-dir "$DATA_DIR/cellranger_out" \
#       --samples-file   "$DATA_DIR/renamed_fastq/unique_samples.txt" \
#       --metadata-tsv   sample_metadata_final.tsv \
#       --out            "$DATA_DIR/all_samples_combined.h5ad"
#
# On the cluster (how the thesis run was done) this is launched through 'sbatch submit_h5ad_concat.slurm',
# which fills the four paths from DATA_DIR and asks for the ~150G the concatenation needs.

import argparse
import scanpy as sc
import anndata as ad
import pandas as pd
from pathlib import Path


def main():
    # Parse input/output paths
    ap = argparse.ArgumentParser()
    ap.add_argument("--cellranger-dir", default="cellranger_out")
    ap.add_argument("--samples-file", default="unique_samples.txt")
    ap.add_argument("--metadata-tsv", default="sample_metadata_final.tsv")
    ap.add_argument("--out", default="all_samples_combined.h5ad")
    args = ap.parse_args()

    # Print library versions to record them for Methods (concat/merge behaviour is version-dependent)
    print(f"anndata {ad.__version__} | scanpy {sc.__version__}")

    # Load the metadata table, indexed by sample_id for per-sample lookup
    cellranger_dir = Path(args.cellranger_dir)
    meta_df = pd.read_csv(args.metadata_tsv, sep="\t").set_index("sample_id")

    # Read the list of samples to combine
    with open(args.samples_file) as f:
        samples = [line.strip() for line in f if line.strip()]
    print(f"Samples to process: {len(samples)}")

    # Use a list: passing a dict to ad.concat makes anndata treat the
    # keys as labels, and the exact behaviour is version-dependent; a list is
    # explicit and version-robust. Trackers record skipped / unmatched samples
    adatas = []
    missing, unmatched_meta = [], []

    for sample in samples:
        # Locate each sample's filtered matrix; skip and record samples without output.
        h5_path = cellranger_dir / sample / "outs" / "filtered_feature_bc_matrix.h5"
        if not h5_path.exists():
            print(f"[SKIP] Missing output: {sample}")
            missing.append(sample)
            continue

        # Read the 10x matrix and make gene symbols unique (kept as var_names for CellTypist).
        adata = sc.read_10x_h5(h5_path)
        adata.var_names_make_unique()

        # Preserve the raw integer counts in a dedicated layer BEFORE any downstream
        # normalisation overwrites .X. Required by scVI/scANVI/DRVI (need raw counts)
        # and by several scib metrics.
        adata.layers["counts"] = adata.X.copy()

        # Prefix barcodes with the sample name so they are unique across samples before the merge.
        adata.obs_names = [f"{sample}_{bc}" for bc in adata.obs_names]

        # Inject the per-sample metadata into .obs; for samples missing from the table,
        # still set .obs['sample'] so the column exists after concat.
        if sample in meta_df.index:
            row = meta_df.loc[sample]
            adata.obs["sample"] = sample
            adata.obs["batch"] = row["batch"]
            adata.obs["cohort"] = row["cohort"]
            adata.obs["treatment"] = row["treatment"]
            adata.obs["fraction"] = row["fraction"]
            adata.obs["dataset_origin"] = row["dataset_origin"]
            adata.obs["response"] = row["response"]
        else:
            print(f"[WARN] No metadata found for: {sample}")
            unmatched_meta.append(sample)
            adata.obs["sample"] = sample

        # Keep per-sample logging for traceability of the run.
        adatas.append(adata)
        print(f"[OK] {sample}: {adata.n_obs} cells, {adata.n_vars} genes")

    # Report load outcome and any problems before merging.
    print("\n--- Summary ---")
    print(f"Samples loaded successfully: {len(adatas)}")
    if missing:
        print(f"Missing samples (h5 not found): {missing}")
    if unmatched_meta:
        print(f"Samples without clinical metadata: {unmatched_meta}")

    # Concatenate all samples: outer join on genes; index_unique=None because barcodes
    # are already unique (prefixed above); merge='same' propagates .var columns shared
    # across samples (notably .var['gene_ids'], the Ensembl IDs) into the combined object.
    print("\nMerging into a single AnnData...")
    combined = ad.concat(adatas, join="outer", index_unique=None, merge="same")

    # Write the combined object and report its final dimensions.
    combined.write_h5ad(args.out, compression="gzip")
    print(f"\nFinal total: {combined.n_obs} cells, {combined.n_vars} genes")
    print(f"Saved to: {args.out}")


if __name__ == "__main__":
    main()