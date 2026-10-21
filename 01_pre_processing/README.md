# 01_pre_processing

Second phase of the thesis: from the raw-count object of `00_raw_data_processing` to `shiao.h5ad`,
a quality-controlled, normalized, annotated and dimensionally-reduced unintegrated AnnData, the
input of the integration benchmark (phase 02).

## Requirements

- Conda environment `benchmark-py-r` (`environments/benchmark-py-r.yml`): scanpy, scib, CellTypist,
  and rpy2 + R scran for the normalization.
- `$DATA_DIR/all_samples_combined.h5ad` from phase 00, plus `regev_lab_cell_cycle_genes.txt` and
  `Cells_Adult_Breast.pkl` (both versioned in `datasets/`).
- ~470G of RAM for scran on 619k cells. In the thesis only this step ran on the cluster; everything
  else ran locally. The `module load` lines of the `.slurm` files are specific to the thesis cluster.

## Quick start

The three notebooks (01_1, 01_2, 01_6) are run by hand; the scripts in between are chained by
`preprocessing_all.sh`:

```bash
cd 01_pre_processing
export DATA_DIR=<data root>
# 1. run by hand: 01_1_load_and_check/load_and_check.ipynb, then 01_2_qc_scrublet_filtering/qc_scrublet_filter.ipynb
bash preprocessing_all.sh --dry-run   # what is already done, what would run
bash preprocessing_all.sh --slurm     # steps 3-8 as one job on the cluster
bash preprocessing_all.sh             # or locally, in sequence
# 2. run by hand: 01_6_visualization/visualization_unintegrated.ipynb
```

`preprocessing_all.sh` skips the steps whose output already exists, so running it again after a
failure resumes where it stopped (`--force` re-runs everything). Step names (`norm`, `cc`, `annot`,
`fraction`, `reduce`, `cluster`) run only those steps. Logs go to `logs/`. A single script can also
be run on its own, e.g. `python3 01_4_cc_and_annotation/cell_cycle_score.py`, or
`cd 01_3_normalization && sbatch submit_scran_norm.slurm` for scran alone.

## Execution order

| # | Folder | Script | What it does | Where |
|---|--------|--------|--------------|-------|
| 1 | `01_1_load_and_check` | `load_and_check.ipynb` | Sanity-check the concatenated object: raw-count `.X`/`counts`, unique gene symbols, metadata fields; `cohort × response/treatment/fraction`, NaN checks. No output file. | local (notebook) |
| 2 | `01_2_qc_scrublet_filtering` | `qc_scrublet_filter.ipynb` | QC metrics (mt/ribo), per-biopsy Scrublet, cell filter (doublets, top-1% UMI, `pct_mt<10`, `n_genes>100`), keep genes in ≥3 cells. QC violins/scatters. | local (notebook) |
| 3 | `01_3_normalization` | `scran_norm.py` | scran normalization via scib (`quickCluster → computeSumFactors → logNormCounts`); adds `size_factors`, clears `.raw`, casts `.X` to float32. | HPC (SLURM) or local |
| 3 | `01_3_normalization` | `submit_scran_norm.slurm` | SLURM wrapper for the scran step alone. | HPC (SLURM) |
| 4 | `01_4_cc_and_annotation` | `cell_cycle_score.py` | Tirosh/Regev cell-cycle scoring → `S_score`, `G2M_score`, `phase`. | local (terminal) |
| 5 | `01_4_cc_and_annotation` | `celltypist_annotation.py` | CellTypist annotation (`Cells_Adult_Breast.pkl`, majority voting) → `cell_type`, `celltypist_predicted`. | local (terminal) |
| 6 | `01_4_cc_and_annotation` | `fraction_reassignment.py` | Reassign `fraction` CD45+/CD45- → `imm`/`non_imm` from the CellTypist lineage (in place); `dataset_origin` left untouched. | local (terminal) |
| 7 | `01_5_scib_pp` | `scib_reduce_data.py` | Batch-aware HVG (`cohort`, 2000) + PCA(50) + neighbors + UMAP; also writes the selected HVG list. | local (terminal) |
| 8 | `01_5_scib_pp` | `scib_clustering.py` | Leiden optimal-resolution sweep (0.1-1.0) vs `cell_type` NMI → `shiao.h5ad`. | local (terminal) |
| 9 | `01_6_visualization` | `visualization_unintegrated.ipynb` | Figures on the final unintegrated object: cell counts and composition, sort-vs-lineage check, UMAPs, Leiden sweep. Read-only. | local (notebook) |

Step 6 must run after step 5 (it needs `cell_type`) and before step 7, so the recoded `fraction`
ends up in `shiao.h5ad`.

## Data layout

Each step reads the previous `.h5ad` and writes the next one, all under `$DATA_DIR`:

    $DATA_DIR/
    ├── all_samples_combined.h5ad                                      # input, from 00_6
    ├── all_samples_combined_scrublet_unfiltered.h5ad                  # 01_2  checkpoint after Scrublet, before filtering
    ├── all_samples_combined_scrublet.h5ad                             # 01_2  QC + Scrublet + filters
    ├── all_samples_combined_scrublet_norm.h5ad                        # 01_3  scran log-norm
    ├── all_samples_combined_scrublet_norm_cc.h5ad                     # 01_4  + cell cycle
    ├── all_samples_combined_scrublet_norm_cc_annotated.h5ad           # 01_4  + CellTypist, fraction recoded in place
    ├── all_samples_combined_scrublet_norm_cc_annotated_reduced.h5ad   # 01_5  + HVG/PCA/neighbors/UMAP
    ├── shiao_hvg_2k_unintegrated_list.csv                             # 01_5  selected HVGs
    └── shiao.h5ad                                                     # 01_5  + Leiden, final object

The repo keeps only the figures, in `figures/<step>/`.

## Key parameters

- **QC & filtering (01_2):** Scrublet per biopsy (`expected_doublet_rate=0.10`), doublet if
  `doublet_score > 0.25` (fixed cutoff, Scrublet's automatic one was unreliable on small biopsies);
  top 1% of cells by UMI count per biopsy removed as presumed doublets; keep `pct_counts_mt < 10`
  and `n_genes_by_counts > 100`; keep genes detected in at least 3 cells (zero-count genes break scran).
- **Normalization (01_3):** `scib.preprocessing.normalize` (`min_mean=0.1`, leiden `quickCluster`),
  i.e. counts divided by the scran size factors, then log1p. `.X` is in log space from here on.
- **Cell cycle (01_4):** `sc.tl.score_genes_cell_cycle` with the Tirosh/Regev list (43 S + 54 G2M),
  `random_state=0`.
- **Annotation (01_4):** CellTypist `Cells_Adult_Breast.pkl` (Kumar et al. 2023, 58 labels),
  `majority_voting=True`, run on a temporary CP10K+log1p copy because that is the normalization
  CellTypist was trained on. The `Lymph-*` labels are lymphatic endothelial cells, so non-immune.
- **HVG + reduction (01_5):** `scib.preprocessing.reduce_data`, 2000 HVGs selected per patient
  (`batch_key='cohort'`) and merged, PCA 50 components, neighbors + UMAP.
- **Clustering (01_5):** `scib.clustering.cluster_optimal_resolution`, leiden at 0.1-1.0, kept the
  resolution with the best NMI vs `cell_type`. All per-resolution columns are kept for the 01_6 plots.

## Expected output

`shiao.h5ad` should contain (numbers from the thesis object, built with an earlier gene filter
with `min_cells=1`; the current code with `min_cells=3` gives 619,692 cells × 29,206 genes):

- **619,693 cells × 30,869 genes**, 34 patients, 48 of the 58 CellTypist labels.
- `.X` = scran log-normalized (do not log1p again); `.layers['counts']` = raw counts;
  gene symbols as `var_names`, Ensembl IDs in `.var['gene_ids']`.
- `dataset_origin` = CD45 sort, unchanged from phase 00: 489,577 immune + 130,116 non-immune.
- `fraction` = immune class from CellTypist: 443,083 `imm` / 176,610 `non_imm`. The two disagree
  for 10.68% of CD45+ cells (annotated non-immune) and 4.44% of CD45- cells (annotated immune).
- Integration keys for phase 02: `batch_key = 'cohort'`, `label_key = 'cell_type'`.
