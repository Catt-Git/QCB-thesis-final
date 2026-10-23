# 03_1_infercnv

First step of phase 03: decide which epithelial cells are malignant. The CellTypist model of 01_4
(`Cells_Adult_Breast.pkl`) is trained on normal breast and has no malignant class, so tumour cells
carry normal epithelial labels. inferCNV infers copy-number changes from expression against a
diploid reference; the aneuploid epithelium gets the single label `malignant`, and the 01_4
annotation is re-run on the remaining cells.

## Requirements

> [!IMPORTANT]
> This step uses **two conda environments**. inferCNV runs in its own environment `infercnv-r`
> (`environments/infercnv-r.yml`: R 4.5, `bioconductor-infercnv`, JAGS), to be created once.
> The Python scripts and the notebooks run in `benchmark-py-r`. The drivers (`infercnv_all.sh`,
> `sensitivity_resolution.sh`, `utils/sensitivity_reference.sh`) switch environment per stage with
> `conda run`; when a script is run by hand, `run_infercnv.R` and `utils/sensitivity_reference.R`
> need `infercnv-r`, everything else `benchmark-py-r`.

- `$DATA_DIR/shiao.h5ad` from phase 01 (raw counts in `.layers['counts']`, `celltypist_predicted`
  in `.obs`) and `$DATA_DIR/Cells_Adult_Breast.pkl`.
- Internet access once, to download the inferCNV gene ordering file (hg38, GENCODE v27).
- 1-3 GB of disk per cohort while inferCNV runs: its working directory is deleted after each cohort.

## Quick start

```bash
cd 03_drvi_tumoral_epi/03_1_infercnv
export DATA_DIR=<data root>
conda env create -f ../../environments/infercnv-r.yml # once
bash infercnv_all.sh --dry-run # what is already done, what would run
bash infercnv_all.sh --threads 12 # inputs + inferCNV on every cohort
# run by hand, in benchmark-py-r: call_malignant.ipynb
conda activate benchmark-py-r
python3 recelltypist_nonmalignant.py
```

`infercnv_all.sh` skips the cohorts whose summary already exists, so running it again after a
failure resumes where it stopped. `--force` re-runs everything, and is also needed after changing
the inputs (otherwise the prepare stage keeps the existing ones). Logs go to `logs/`.

## Execution order

| # | Script | What it does | Environment |
|---|--------|--------------|-------------|
| 1 | `prepare_infercnv_input.py` | Per-cohort counts and cell groups; downloads the gene ordering file once | `benchmark-py-r` |
| 2 | `run_infercnv.R` | inferCNV on one cohort → per-cell table + heatmap | `infercnv-r` |
| 1-2 | `infercnv_all.sh` | Driver: step 1, then step 2 on every prepared cohort | both |
| 3 | `call_malignant.ipynb` | The call: per-subcluster verdict, cohort gate → `cnv_status.csv` | `benchmark-py-r` |
| 4 | `recelltypist_nonmalignant.py` | CellTypist re-run without the malignant cells → `cell_annotation.csv` | `benchmark-py-r` |
| S1 | `sensitivity_resolution.sh` + `.py` | Same call at 4 leiden resolutions on 4 cohorts | both |
| - | `utils/sensitivity_reference.sh` + `.py` + `.R` | Same call with 4 reference configurations on 2 cohorts (not in the thesis) | both |

## Design

- **Reference: immune cells, per patient.** The dataset has no normal tissue (BASE / PD1 / RTPD1
  are three timepoints of the same tumour), so the diploid baseline comes from the immune cells of
  each tumour: T/NK and myeloid, up to 1,000 each, as two separate groups so that a gene simply
  higher in myeloid than in T cells is not read as a gain. One run per patient, because the
  residual is taken against the reference cells of the run. B/plasma cells are left out: their
  immunoglobulin loci would add structure on chr2, chr14 and chr22.
- **Stromal null.** Up to 2,000 fibroblasts, endothelial and perivascular cells, drawn stratified
  over the three, go through inferCNV as observations. The thresholds are quantiles of this null,
  per cohort.
- **Specificity control.** 500 immune cells are held out of the reference and run as observations.
  The reference cells themselves are not a check: they sit near zero by construction.
- **Two axes.** `cnv_score` is the mean squared residual against the mean of the reference;
  `cnv_corr` is the correlation with the mean residual profile of the top 5% epithelial cells of
  the cohort. Noise raises the first but not the second. Each cell's position on its cohort's null
  (`nullpos`) is the lower of its two stromal percentiles.
- **Verdict per subcluster.** inferCNV partitions each run with leiden (resolution 0.005), so the
  subclusters never mix patients. A subcluster is malignant when its median `nullpos` is ≥ 0.90,
  normal when ≤ 0.75, borderline in between; under 20 cells, its cells fall back to the per-cell
  rule (`nullpos` > 0.95).
- **Cohort gate.** A cohort needs at least 10% of its epithelium above the null (per cell or in
  malignant subclusters); otherwise its epithelium is `not_tested`.
- **No HMM.** `run_infercnv.R --hmm` turns on inferCNV's i6 HMM, which takes hours per cohort and
  is not needed for the call.

## The values of `cnv_status`

| value | cells | meaning |
|---|---|---|
| `malignant` | epithelium in a malignant subcluster of a cohort that passed the gate | aneuploid |
| `borderline` | epithelium where subcluster and per-cell rule disagree | evidence both ways |
| `non_malignant` | everything else a run covered | diploid as far as this test can tell |
| `not_tested` | B/plasma cells, cells not sampled into a run, epithelium of gated-out cohorts | no evidence either way |

`not_tested` is never merged with `non_malignant`. 03_2 keeps only `malignant` cells.
`recelltypist_nonmalignant.py` still gives every non-malignant cell a normal-breast label, with
`cnv_status` next to it, because a smaller input would change CellTypist's voting neighbourhoods.

## Sensitivity checks

- **Resolution** (`sensitivity_resolution.*`, 0.005-0.05 on Patient16, 43, 64, 52): the
  per-patient verdict never changes. At 0.005, 98-100% of the epithelium sits in subclusters of at
  least 20 cells, against 24-89% at 0.05, hence the operating value. It is close to inferCNV's
  `auto` value for the largest cohorts (~0.003-0.006 for 8,000-9,000 epithelial cells).
- **Reference** (`utils/sensitivity_reference.*`, 4 configurations on Patient16 and Patient52,
  not reported in the thesis): the AUC
  separating epithelium from the stromal null stays within 0.972-0.975 on Patient16 and at chance
  (0.43-0.49) on Patient52; the held-out immune cells are never called.

Tables in `../tables/03_1_infercnv/sensitivity/`; the production `input/` and `summary/` are never
touched.

## Data layout

    $DATA_DIR/03_tum/
    ├── gene_order_hg38_gencode_v27.txt # 03_1 downloaded once
    ├── cohort_census.csv # 03_1 cells per group in each run
    ├── input/<cohort>/ # 03_1 counts.mtx, genes.tsv, barcodes.tsv, annotations.tsv
    ├── summary/<cohort>_cnv.csv # 03_1 one row per cell: group, subcluster, cnv_score, cnv_corr, chr1..chr22
    ├── cnv_status.csv # 03_1 the call, one row per cell of shiao.h5ad
    ├── cell_annotation.csv # 03_1 the re-annotation, one row per cell of shiao.h5ad
    └── sensitivity/ # 03_1 sensitivity runs

The call is metadata, so it is written as per-cell tables, not as a new `.h5ad`; 03_2 joins it onto
the object it builds. The repo keeps the figures (`../figures/03_1_infercnv/`, one heatmap per
cohort plus `malignant_per_patient.png` from the notebook) and the tables (`../tables/03_1_infercnv/`).

## Expected output

`cnv_status.csv`, 619,693 rows. Epithelial cells by status:

- 42,291 `malignant`, 2,749 `borderline`, 16,839 `non_malignant`, 13,573 `not_tested`.
- 19 cohorts have at least 200 malignant cells, for 41,744 cells: the tumour subset of 03_2.

## Notes

- **`cnv_score` is not comparable with published values.** It is a mean of squares over the whole
  genome on inferCNV's linear scale, with denoising: a median of ~0.0028 on the malignant cells is a
  ~5.3% RMS deviation, against ~0.0007 (~2.7%) for the T/NK reference. Only contrasts within a run are meaningful;
  `cnv_corr` has no scale.
- **Subcluster names are per run**: `epi_s1` exists in every cohort, so the key is always
  (cohort, subcluster).
- **One cohort is never run**: Patient06 has 14 epithelial cells, below the minimum of 50.
