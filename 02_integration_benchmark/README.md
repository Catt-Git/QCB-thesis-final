# 02_integration_benchmark

Third phase of the thesis: the technical benchmark. Ten batch-correction methods are run on the
unintegrated object from phase 01, each output is scored with the 12 scIB metrics
(Luecken et al. 2022) against the unintegrated data, and the results are collected into one
summary table.

Integration keys: `batch_key = 'cohort'` (34 patients), `label_key = 'cell_type'` (48 CellTypist labels).

## Requirements

- `$DATA_DIR/shiao.h5ad` and `$DATA_DIR/shiao_hvg_2k_unintegrated_list.csv` from phase 01.
- Conda environments (`environments/`):
  - `benchmark-py-r` for everything except scGen. It must be activated, not just its `python`
    called, otherwise rpy2 cannot find R.
  - `scgen-py` for scGen only: its pinned stack conflicts with the main one. On the cluster scGen
    uses `$SCGEN_ENV`, an equivalent environment.
- A GPU for the integrations (scVI, scANVI, scGen, DRVI; tested on a 4 GB T1000) and a SLURM
  cluster for the metrics. Both steps can run either way, only much slower. The `module load`
  lines of the `.slurm` files are specific to the thesis cluster.

## Smoke tests: run them first

The full benchmark is the longest computation in the repository: each integration takes from
minutes to hours on 620k cells, and so does each of the 22 metrics jobs, so an error found late
costs days. The smoke tests run the same code on a 5,252-cell fixture in a few minutes and tell
whether environments and scripts are ready for the real run. They check that everything runs and
writes the expected outputs, not the biology.

```bash
cd 02_integration_benchmark
export DATA_DIR=<data root>
conda activate benchmark-py-r

# 1. metrics stack: scib, rpy2 and all 12 metrics (builds smoke_fixture.h5ad on its first run)
python3 utils/smoke_test_metrics.py # on the cluster: cd utils && sbatch submit_smoke_test.slurm

# 2. every integration method on the fixture
S=utils/smoke_test
python3 utils/make_smoke_input.py
Rscript 02_1_prepare/h5ad_to_rds.R -i $S/smoke_hvg.h5ad -o $S/smoke_hvg.rds
Rscript 02_1_prepare/hvg_csv_to_rds.R -i $S/smoke_hvg_list.csv \
    -o $S/smoke_hvg_list.rds --n-expected 1944
O=$S/smoke_out
for m in bbknn harmony scvi scanvi scanorama; do
  python3 02_2_integration/run_integration.py -m $m \
      -i $S/smoke_hvg.h5ad -o $O/$m.h5ad --emb-out $O/$m.npy --model-dir $O/models
done
for m in fastmnn seurat_cca seurat_rpca; do
  Rscript 02_2_integration/run_integration.R -m $m \
      -i $S/smoke_hvg.rds -o $O/$m.rds -v $S/smoke_hvg_list.rds
done
python3 02_2_integration/rds_to_h5ad.py -i $O/fastmnn.rds -o $O/fastmnn.h5ad --types full,embed --emb-out $O/fastmnn.npy
python3 02_2_integration/rds_to_h5ad.py -i $O/seurat_cca.rds -o $O/seurat_cca.h5ad --types full
python3 02_2_integration/rds_to_h5ad.py -i $O/seurat_rpca.rds -o $O/seurat_rpca.h5ad --types full
conda run -n scgen-py python3 02_2_integration/run_scgen.py -i $S/smoke_hvg.h5ad -o $O/scgen.h5ad --model-dir $O/models

# 3. the whole metrics pipeline on those smoke integrations (fails if an expected metric has no value)
bash utils/smoke_test_metrics_pipeline.sh
```

Everything the smoke tests write goes to `utils/smoke_test/` (not in git, disposable; set
`SMOKE_DIR` to put it elsewhere). `DATA_DIR` is only read, to build the fixture.

## Quick start

```bash
cd 02_integration_benchmark
export DATA_DIR=<data root>

# 1. inputs
python3 02_1_prepare/subset_hvg.py
python3 02_1_prepare/scale_batch.py
Rscript 02_1_prepare/h5ad_to_rds.R -i $DATA_DIR/shiao_hvg_2k.h5ad -o $DATA_DIR/shiao_hvg_2k.rds
Rscript 02_1_prepare/h5ad_to_rds.R -i $DATA_DIR/shiao_hvg_2k_scaled.h5ad -o $DATA_DIR/shiao_hvg_2k_scaled.rds
Rscript 02_1_prepare/hvg_csv_to_rds.R -i $DATA_DIR/shiao_hvg_2k_unintegrated_list.csv \
    -o $DATA_DIR/shiao_hvg_2k_unintegrated_list.rds

# 2. integrations, locally on the GPU
bash 02_2_integration/run_all.sh # every grid row except DRVI (--slurm: on the cluster)
python3 02_2_integration/run_drvi.py # DRVI, n_latent 128 (or the notebook shiao_drvi_128.ipynb)

# 3. UMAP QC panels, one folder per run in figures/
bash 02_3_plot_method_umap/plot_all.sh

# 4. metrics, on the cluster
export CLUSTER_HOST=<user>@<hostname>
bash utils/sync_to_cluster.sh # upload the integrated objects and their references
#   then, on the cluster:
#   bash 02_4_metrics/run_all_metrics.sh --slurm

# 5. summary table, after copying $DATA_DIR/02_metrics/ back from the cluster
python3 02_4_metrics/merge_metrics.py -o $DATA_DIR/02_metrics_merged.csv -r $DATA_DIR/02_metrics \
    --glob "$DATA_DIR/02_metrics/shiao/metrics/*/hvg/*.csv"
Rscript 02_4_metrics/make_summary_table.R -i $DATA_DIR/02_metrics_merged.csv
```

The three drivers (`run_all.sh`, `plot_all.sh`, `run_all_metrics.sh`) all read `benchmark_grid.tsv`
and share the same options: `--dry-run` previews, `--method`, `--scaling` or run ids select part
of the grid, `--slurm` submits a SLURM array instead of running locally, `--force` redoes what is
already done. They skip whatever is already on disk, so running one again after a crash resumes
where it stopped.

## Execution order

| # | Folder | Script | What it does | Where |
|---|--------|--------|--------------|-------|
| 1 | `02_1_prepare` | `subset_hvg.py` | Subset `shiao.h5ad` to the 2,000 HVGs → `shiao_hvg_2k.h5ad` | local |
| 2 | `02_1_prepare` | `scale_batch.py` | Per-batch z-scoring, then a new PCA and UMAP on the scaled matrix → `shiao_hvg_2k_scaled.h5ad` | local |
| 3 | `02_1_prepare` | `h5ad_to_rds.R`, `hvg_csv_to_rds.R` | Seurat `.rds` of both variants and of the HVG list, for the R methods | local |
| 4 | `02_2_integration` | `run_all.sh` | Driver: every grid row except DRVI, locally or as a SLURM array (`submit_integration.slurm`) | local or HPC |
| 4a | `02_2_integration` | `run_integration.py` + `integration_methods.py` | Python methods: BBKNN, Scanorama, Harmony, scVI, scANVI | local or HPC |
| 4b | `02_2_integration` | `run_integration.R` + `integration_methods.R`, then `rds_to_h5ad.py` | R methods: fastMNN, Seurat CCA, Seurat RPCA, converted back to `.h5ad` | local or HPC |
| 4c | `02_2_integration` | `run_scgen.py` | scGen, in its own environment | local or HPC |
| 4d | `02_2_integration` | `shiao_drvi_128.ipynb` or `run_drvi.py` | DRVI: the notebook, or the same run headless (`submit_drvi.slurm` on the cluster) | local or HPC |
| 5 | `02_3_plot_method_umap` | `plot_all.sh` → `plot_methods_umaps.py` | Driver: five UMAP QC panels per run, a visual check before the metrics | local or HPC |
| 6 | `utils` | `sync_to_cluster.sh` | Upload the files the grid names to the cluster's `DATA_DIR`, same layout | local |
| 7 | `02_4_metrics` | `run_all_metrics.sh` → `check_integrations.py` + `metrics.py` | Driver: pre-flight check and the 12 metrics, one SLURM task per (run, output type) (`submit_metrics.slurm`) | HPC or local |
| 8 | `02_4_metrics` | `merge_metrics.py`, `make_summary_table.R` | Merge the per-run CSVs; draw the summary table | local |

Shared modules in `utils/`: `metrics_shared.py` (identical preparation for every run),
`scib_compat.py` (numpy/pandas shims for scib 1.1.7), `h5ad_compat.py` (`.h5ad` that anndata 0.10
can read). `02_2_integration/model_paths.py` defines where the trained models are saved.

## The run grid

`benchmark_grid.tsv` is the single source of truth: one row per run, 18 rows (17 integrations plus
the unintegrated baseline). Its columns name the method, the dispatcher (`language`), the conda
environment, the scaling variant, the input, the output, the scIB output type(s), the unintegrated
reference the metrics score against and the `--hvgs` value. A row with two output types is scored
twice, so the 18 rows give 22 metrics jobs.

| Method | Language | Output type(s) | Scaled variant |
|---|---|---|---|
| BBKNN | Python | knn | yes |
| Scanorama | Python | full + embed | yes |
| Harmony | Python (`harmony-pytorch`) | embed | yes |
| fastMNN | R (`batchelor`) | embed + full | yes |
| Seurat CCA | R | full | yes |
| Seurat RPCA | R | full | yes |
| scGen | Python (`scgen-py`) | full | yes |
| scVI | Python | embed | no |
| scANVI | Python | embed | no |
| DRVI | Python (notebook) | embed | no (`n_latent = 128`) |

scVI, scANVI and DRVI read the raw counts and ignore `.X`, so a scaled variant would repeat the
unscaled run. The unintegrated baseline is the "before" row of the summary table: it is the input
object scored against itself, so some of its metrics are 0 or 1 by construction.

## DRVI

DRVI is outside `run_all.sh` because its latent size was chosen by looking at the trained model.
The notebook `shiao_drvi_128.ipynb` and `run_drvi.py` do the same run with the same outputs;
`submit_drvi.slurm` runs the script on the cluster (CPU, partition `long`).

**Vanished dimensions.** The DRVI paper calls a latent dimension vanished when its maximum absolute
value is below 1 and evaluates the embedding without them, while the metrics here read all 128
dimensions. `02_3_plot_method_umap/plot_drvi_pruned_umap.py --report-only` checks the effect: the
63 vanished dimensions hold 8.6e-06 of the latent variance and the 15-nearest-neighbour graph,
on which every scIB metric is built, is 99.997% identical with and without them. The scores are
therefore the same as on the pruned embedding. The pruned UMAP looks different only because a UMAP
layout changes from run to run.

## Metrics

| Batch correction (4) | Bio conservation (8) |
|---|---|
| PCR batch | NMI cluster/label |
| batch ASW | ARI cluster/label |
| graph iLISI | cell type ASW |
| graph connectivity | isolated label F1 |
| | isolated label ASW |
| | graph cLISI |
| | cell cycle conservation |
| | HVG conservation |

Trajectory conservation, in the original paper, is not used here. Not every metric exists for every
output type: `knn` output (BBKNN) supports only the graph metrics, and `embed` output has no HVG
conservation, so the summary table has empty cells by construction. PCR batch is not computed on
the scaled runs: per-batch scaling already removed the batch means from their reference, so the
score would be 0 by construction.

Overall score = 0.6 × bio conservation + 0.4 × batch correction, on min-max scaled metrics, so a
score is comparable only within one table. The summary table is drawn with the official scIB
plotting code, vendored in `utils/`.

## Expected output

`figures/02_4_metrics/shiao_summary_scores.csv` and `<timestamp>_shiao_summary_metrics.pdf/png`:

- Top of the ranking by overall score: scGen* 0.75, scANVI* 0.60, Harmony 0.57 (unscaled),
  Scanorama (embed, scaled) 0.53, then Seurat RPCA 0.51, Scanorama (embed) 0.48 and DRVI 0.47
  (unscaled); the unintegrated baseline scores 0.44.
- scGen and scANVI (starred) receive the cell-type labels in training, the same labels the
  bio-conservation metrics use, so they rank first by construction.
- Every method except BBKNN and Scanorama (embed) scores lower on the scaled variant than on the
  unscaled one.

`figures/<run_id>/` holds the five UMAP QC panels of each run, plus the latent-space figures for DRVI.

## Data layout

    $DATA_DIR/
    ├── shiao.h5ad, shiao_hvg_2k_unintegrated_list.csv # input, from phase 01
    ├── shiao_hvg_2k.h5ad / .rds # 02_1 HVG subset
    ├── shiao_hvg_2k_scaled.h5ad / .rds # 02_1 HVG subset, z-scored per batch
    ├── shiao_hvg_2k_unintegrated_list.rds # 02_1 HVG list for the Seurat anchors
    ├── 02_integration/<run_id>.h5ad # 02_2 integrated objects (02_3 adds their UMAP)
    ├── 02_embeddings/<run_id>.npy # 02_2 latent spaces
    ├── 02_scvi/, 02_scanvi/, 02_scgen/, 02_drvi/ # 02_2 trained models (+ DRVI latent space)
    ├── 02_metrics/shiao/metrics/<scaling>/hvg/*.csv # 02_4 per-run metrics
    └── 02_metrics_merged.csv # 02_4 input of the summary table

## Notes for reproducibility

- **Trained models are reused.** scVI, scANVI and scGen save their model after training and reload
  it on a re-run. To train again, delete the `.pt`: `--force` alone re-runs the integration with
  the saved model.
- **The reference must match the scaling variant**, for the metrics and the UMAP panels alike:
  per-batch scaling is a weak batch correction in itself, so a scaled run is compared with the
  scaled input.
- **`--hvgs 0`**: the references are already restricted to the HVGs, and letting scib select them
  again fails on the scaled (negative) values.
- **Cell order must not change**: the metrics compare objects cell by cell, so every dispatcher and
  `check_integrations.py` assert it.
- **`scale_batch.py` replaces `scib.preprocessing.scale_batch`**, which no longer works with
  anndata 0.13: same operation (`sc.pp.scale` per batch, no clipping), and it recomputes the PCA
  and the UMAP on the scaled matrix.
- **Every script imports `utils/scib_compat.py` before `scib`**: without it scib 1.1.7 fails
  halfway through the metrics on current numpy/pandas.
- **Seurat uses the legacy anchor API** (`FindIntegrationAnchors` + `IntegrateData`) so that it
  returns a corrected matrix (`full`), with the three largest patients as reference: full pairwise
  anchoring is not tractable on 34 batches and 620k cells.
- **scGen decodes in chunks** (`--decode-chunk`, 16,384 cells) to fit on a 4 GB GPU.
- **Leiden resolution sweep up to 1.0** for NMI/ARI, as in phase 01: it is the most expensive part
  of each metrics job.
