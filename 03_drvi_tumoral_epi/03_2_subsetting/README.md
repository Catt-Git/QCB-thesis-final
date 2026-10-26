# 03_2_subsetting

From `shiao.h5ad` to the object the rest of phase 03 runs on: the cells 03_1 found aneuploid,
with the phase-01 pre-processing re-run on them. The subset is defined by the inferCNV call,
not by a label; the CellTypist labels come along as metadata. scran, HVGs, PCA, neighbours and
UMAP are all recomputed. The cells are not clustered here: the clustering of these cells is
done on the DRVI latent space, in 03_3.

## Requirements

- Conda environment `benchmark-py-r` (scanpy, scib, rpy2 + R scran).
- `$DATA_DIR/shiao.h5ad`, and `cnv_status.csv` + `cell_annotation.csv` in
  `$DATA_DIR/03_tum/` from 03_1. `$DATA_DIR/regev_lab_cell_cycle_genes.txt` for the cell cycle.

## Quick start

```bash
cd 03_drvi_tumoral_epi/03_2_subsetting
export DATA_DIR=<data root>
# 1. run by hand: subset_and_qc.ipynb
bash subsetting_all.sh --dry-run # what is already done, what would run
bash subsetting_all.sh # steps 3-5, in sequence (--slurm: as one job on the cluster)
# 2. run by hand: visualization_tum.ipynb
```

`subsetting_all.sh` skips the steps whose output exists, so it resumes after a failure
(`--force` re-runs everything); step names (`norm`, `cc`, `reduce`) run only those. Logs go to `logs/`. ~42k cells: everything runs locally, `--slurm` is there for when the data
live on the cluster.

## Execution order

| # | File | What it does |
|---|------|--------------|
| 1 | `subset_and_qc.ipynb` | Join the 03_1 call and re-annotation onto `shiao.h5ad`; keep `cnv_status == 'malignant'`; write `compartment`, the post-CNV `cell_type` and the pre-CNV `cell_type_01_4`; restore `.X` = raw counts and drop every phase-01 derived slot; recompute QC and re-apply the cell thresholds; gene filter; drop cohorts under 100 cells |
| 2 | `subsetting_all.sh` | Driver for steps 3-5 |
| 3 | `scran_norm_tum.py` | scran size factors re-estimated on the subset, log-normalisation |
| 4 | `cell_cycle_score_tum.py` | Tirosh/Regev re-scoring → `S_score`, `G2M_score`, `phase` |
| 5 | `reduce_data_tum.py` | Batch-aware HVGs (`cohort`, 2,000, `MT-` genes replaced) + PCA(50) + neighbours + UMAP → `shiao_tum.h5ad` and the DRVI input of 03_3 |
| 6 | `visualization_tum.ipynb` | Figures on the final object, read-only |

[`../utils/cell_set.py`](../utils/cell_set.py) holds every file name and threshold of the phase (`MIN_GENES = 100`,
`MAX_PCT_MT = 10`, gene filter `min_cells = 3`, `MIN_CELLS_PER_COHORT = 100`) and the run id
`drvi_tum_64_nomt` that 03_3 - 03_9 use.

## Key choices

- **Cohorts with ≥ 100 malignant cells.** DRVI fits a dispersion per gene and per cohort, so a
  cohort needs enough cells. Patient10 (199 cells) is kept; the next cohort has 94. Cohorts missing a timepoint are kept
  (`DROP_INCOMPLETE_COHORTS = False`).
- **No `MT-` genes in the DRVI panel.** On this object 13 of the 2,000 HVGs are mitochondrial
  (`MT-ND1`, `MT-CO1`, `MT-CYB`, ...): they track the mitochondrial fraction of a cell's RNA,
  i.e. dissociation stress, not a state. `reduce_data_tum.py` replaces them with the next genes
  of the same batch-aware ranking, reproducing scib's selection and checking it gene for gene
  against scib's own call; PCA, neighbours and UMAP run on this panel. The metallothioneins
  (`MT1*`, `MT2A`) and `MTRNR2L12` are nuclear genes and stay.
- **Scrublet is not re-run**: of the 4,326 cells flagged as doublets in `shiao.h5ad`, 9 are in
  this subset.

## Data layout

    $DATA_DIR/03_tum/
    ├── shiao_tum_raw.h5ad                     # subset_and_qc.ipynb, raw counts
    ├── shiao_tum_norm.h5ad                    # + scran
    ├── shiao_tum_norm_cc.h5ad                 # + cell cycle
    ├── shiao_tum.h5ad                         # + HVGs (no MT-) / PCA / neighbours / UMAP, the definitive object
    └── shiao_tum_hvg_2k_nomt_list.csv, _hvg_2k_nomt.h5ad # the panel alone: the DRVI input

Figures go to `../figures/03_2_*/`, the table of the visualization notebook to
`../tables/03_2_visualization/`.

## Expected output

`shiao_tum.h5ad`: 41,943 cells × 25,197 genes, 20 cohorts; `.X` scran log-normalised,
`.layers['counts']` raw counts; `cell_type` = `malignant`, `cell_type_01_4` the pre-CNV label.
