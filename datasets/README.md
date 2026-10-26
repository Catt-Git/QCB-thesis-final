# datasets (`DATA_DIR`)

Where every heavy object of the thesis lives. The scripts never hardcode this path: they read
`DATA_DIR` from the environment and exit with an explanatory message if it is unset.

```bash
export DATA_DIR=<data root> # (e.g. `/home/catt/QCB-thesis-final/datasets`)
```

The clone itself can sit anywhere: no script contains the repo's path (they self-locate via
`SCRIPT_DIR` / `SLURM_SUBMIT_DIR`), so `DATA_DIR` is the only thing to set after `git clone`.

The Cell Ranger reference (`REF`, ~11 GB, `submit_cellranger.slurm`) stays outside `DATA_DIR`:
it is a shared external input, not something a phase of this pipeline writes. The exact build used
here is 10x Genomics' prebuilt human reference `refdata-cellranger-GRCh38-3.0.0` (the same as the
paper), downloadable with
`curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz` and then
untarred: `REF` points at the resulting directory.


**Almost nothing here is in git.** `.gitignore` excludes `/datasets/*` with four exceptions:
- this README,
- `Cells_Adult_Breast.pkl` (the CellTypist model),
- `regev_lab_cell_cycle_genes.txt` (the Tirosh/Regev signature),
- `signatures/` (the curated gene lists scored in phase 03).
Everything else is rebuilt by re-running the phase that writes it.

## What each phase writes here

### 00_raw_data_processing

| Path | Step |
|---|---|
| `fastq_raw/` | 00_3, the FASTQ as downloaded (~2.4 TB) |
| `renamed_fastq/` | 00_4, Cell Ranger naming + `unique_samples.txt` |
| `cellranger_out/<sample>/outs/filtered_feature_bc_matrix.h5` | 00_5, what 00_6 concatenates |
| `cellranger_out/<sample>/outs/raw_feature_bc_matrix.h5` | 00_5, the unfiltered matrix |
| `all_samples_combined.h5ad` | 00_6, the concatenated raw object |

### 01_pre_processing

| File | Step |
|---|---|
| `all_samples_combined_scrublet_unfiltered.h5ad` | 01_2, checkpoint after Scrublet, before filtering |
| `all_samples_combined_scrublet.h5ad` | 01_2, QC + Scrublet + filters |
| `all_samples_combined_scrublet_norm.h5ad` | 01_3, scran log-normalization |
| `all_samples_combined_scrublet_norm_cc.h5ad` | 01_4, + cell cycle |
| `..._cc_annotated.h5ad` | 01_4, + CellTypist (`fraction` recoded in place) |
| `..._cc_annotated_reduced.h5ad` | 01_5, + HVG / PCA / neighbours / UMAP |
| `shiao_hvg_2k_unintegrated_list.csv` | 01_5, the 2,000 HVG symbols |
| **`shiao.h5ad`** | 01_5, + leiden - **the definitive unintegrated object** (619,693 x 30,869) |

Auxiliary inputs read from here: `regev_lab_cell_cycle_genes.txt` (via `CC_GENES`) and
`Cells_Adult_Breast.pkl`.

### 02_integration_benchmark

| Path | Step |
|---|---|
| `shiao_hvg_2k.h5ad` / `shiao_hvg_2k_scaled.h5ad` | 02_1, the two Python inputs |
| `shiao_hvg_2k.rds` / `shiao_hvg_2k_scaled.rds` | 02_1, the same two for the R methods |
| `shiao_hvg_2k_unintegrated_list.rds` | 02_1, the HVG symbols for the Seurat anchors |
| `02_integration/<run_id>.h5ad` | 02_2, the integrated objects |
| `02_embeddings/<run_id>.npy` | 02_2, latent spaces kept for the figures |
| `02_scvi/`, `02_scanvi/`, `02_scgen/`, `02_drvi/` | 02_2, one flat directory per trained method |
| `02_metrics/shiao/metrics/<scaling>/hvg/<method>_<type>.csv` | 02_4, per-run metric CSVs |
| `02_metrics_merged.csv` | 02_4, the merged table behind the summary figure |

The trained methods keep their model out of `02_integration/` on purpose, so that directory holds
exactly what the metrics read. See "Trained models" in the phase-02 README.

The R methods also write a `02_integration/<run_id>.rds` intermediate, deleted once converted
(`--keep-rds` to keep it): nothing downstream reads it and it is the largest artefact of the phase.

### 03_drvi_tumoral_epi

`03_tum/` holds everything phase 03 writes. The DRVI run is `drvi_tum_64_nomt`.

| Path | Step |
|---|---|
| `03_tum/gene_order_hg38_gencode_v27.txt` | 03_1, downloaded once from the inferCNV authors |
| `03_tum/cohort_census.csv` | 03_1, the cells that went into each per-patient run |
| `03_tum/input/<cohort>/` | 03_1, `counts.mtx` + `genes.tsv` + `barcodes.tsv` + `annotations.tsv` |
| `03_tum/work/<cohort>/` | 03_1, inferCNV's working directory - **deleted on success**, 1-3 GB while it exists |
| `03_tum/summary/<cohort>_cnv.csv` | 03_1, one row per cell: group, subcluster, `cnv_score`, `cnv_corr`, `chr1`..`chr22` |
| `03_tum/sensitivity/` | 03_1, the resolution and reference sensitivity runs |
| **`03_tum/cnv_status.csv`** | 03_1, the call - one row per cell of `shiao.h5ad` |
| **`03_tum/cell_annotation.csv`** | 03_1, the CellTypist re-run - the table 03_2 joins |
| `03_tum/shiao_tum_raw.h5ad`, `_norm`, `_norm_cc` | 03_2, subset, scran, cell cycle |
| **`03_tum/shiao_tum.h5ad`** | 03_2, + HVGs without `MT-` / PCA / UMAP - the definitive object of the phase |
| `03_tum/shiao_tum_hvg_2k_nomt.h5ad` (+ `_list.csv`) | 03_2, the HVG panel alone: the DRVI input |
| `03_tum/model_drvi_tum_64_nomt.pt`, `embed_drvi_tum_64_nomt.h5ad` | 03_3, the trained model and the latent space |
| `03_tum/shiao_tum_drvi_tum_64_nomt.h5ad` | 03_3, the definitive object with the latent space in `obsm['X_drvi']` |
| `03_tum/cytotrace2_drvi_tum_64_nomt.csv` | 03_5, per-cell CytoTRACE2 scores |
| `03_tum/signature_scores_<collection>_drvi_tum_64_nomt.csv` | 03_6, per-cell signature scores |
| `03_tum/random_signature_scores_*`, `drvi_state_argmax*_drvi_tum_64_nomt.csv.gz` | 03_8, cached null draws and max(DR) assignments |
