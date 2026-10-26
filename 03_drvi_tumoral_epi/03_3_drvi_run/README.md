# 03_3_drvi_run

DRVI on the malignant cells: 64 latent dimensions on the 2,000 batch-aware HVGs of 03_2
without the `MT-` genes. Run id **`drvi_tum_64_nomt`**.

## Requirements

- Conda environment `benchmark-py-r` locally (GPU); on the cluster an env with scvi-tools and
  drvi-py 0.2.7 (`DRVI_ENV`, default `benchmark-py-r`), CPU only.
- `$DATA_DIR/03_tum/shiao_tum_hvg_2k_nomt.h5ad` and `shiao_tum.h5ad` from 03_2, both required:
  the first is the training input, the second becomes the downstream object.
- Network for the enrichment section of the notebook (Enrichr, g:Profiler, OmniPath).

## Quick start

Locally, the notebook alone is enough: run `drvi_tum_64.ipynb` top to bottom (from
`03_3_drvi_run/`, with `DATA_DIR` exported). It writes everything 03_4 - 03_9 read, the Leiden
columns included; `OVERWRITE = True` retrains.

Without a kernel, or on the cluster, the runner does the same training, embedding and figures:

```bash
cd 03_drvi_tumoral_epi/03_3_drvi_run
export DATA_DIR=<data root>
python3 run_drvi_tum.py # train, embed, figures
python3 run_drvi_tum.py --overwrite # retrain and rewrite everything
# optional afterwards: drvi_tum_64.ipynb for SMI and enrichment (OVERWRITE = False reads the model from disk)

# on the cluster
mkdir -p logs && sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_drvi_tum.slurm
```

Existing outputs are reused, so a crash after training does not cost the training. On the
cluster, `shiao_tum_hvg_2k_nomt.h5ad`, `shiao_tum.h5ad` and `utils/cell_set.py` have
to be copied up first; `embed_drvi_tum_64_nomt.h5ad` and `shiao_tum_drvi_tum_64_nomt.h5ad`
have to come back (03_4 on, 03_5 and 03_6).

## The files

| File | What it does | Where |
|---|---|---|
| `drvi_tum_64.ipynb` | Trains (or reloads), embeds, clusters (Leiden), writes the embedding and the downstream object, draws every figure, screens the dimensions with SMI and names them by their genes (Enrichr, g:Profiler, decoupler) | local (notebook) |
| `run_drvi_tum.py` | The same training, embedding and figures without a kernel (no SMI or enrichment) | local / cluster |
| `submit_drvi_tum.slurm` | SLURM job for the runner: CPU, partition `long`, 16 CPUs, 32G | cluster |

## Key parameters

- Encoder `[128, 128]` / decoder `[128, 128]`, `dispersion='gene-batch'`, `batch_key='cohort'`,
  `n_latent = 64`, `SEED = 123`, 400 full epochs. The gene-batch dispersion is why 03_2 drops
  cohorts under 100 cells: a dispersion per gene and per batch cannot be fitted on 14 cells.
- **No early stopping**: scvi-tools raises the KL weight over 400 epochs, and that pressure is
  what makes unused dimensions vanish. Vanished dimensions are dropped downstream (03_6 on).
- **No `labels_key`**: the counts layer and `cohort` are the whole input, so no annotation can
  bias the latent space; the groupings below are only read back afterwards, as captions.

## Clustering on the latent space

The vanished dimensions are dropped as soon as the dimension statistics flag them: the
neighbours graph and the UMAP are built on the live dimensions only (`obsm['X_live']`), with
the same k and seed as the Leiden, so the clusters are drawn on the graph they were cut from.

After the integration the cells are clustered with Leiden on the live (non-vanished) DRVI
dimensions: k = 15 neighbours, resolutions 0.25 - 2.0, each kept as `drvi_leiden_<res>` in the
embedding. Clusters under 100 cells (`LEIDEN_MIN_CELLS`) are merged into the cluster that holds
most of their kNN weight: they are a handful of cells cut off a larger group and add noise to
every per-cluster reading. A small cluster with no edge to the rest of the graph (19 cells here)
is kept as it is.

Resolution **1.5** (35 clusters) is the one the figures and 03_4 - 03_9 read, as `drvi_leiden`.
`umap_drvi_leiden_all_resolutions_drvi_tum_64_nomt.png` draws every resolution on the DRVI UMAP,
`drvi_leiden_clusters_drvi_tum_64_nomt.csv` the clusters per resolution, before and after the
merge. The constants live in `utils/cell_set.py` (`LEIDEN_*`); every random step of 03_3 uses
seed 123 (`SEED` in `run_drvi_tum.py`). The settings the columns were made with are stored in
the embedding (`uns['drvi_leiden_params']`): changing one makes the runner and the notebook
recompute the columns.

## Grouping the cells

`cell_type` is `malignant` for every cell, so the figures group cells by `drvi_leiden`, the
Leiden clusters of the DRVI space at resolution 1.5: the UMAP panels, the heatmaps and the SMI
screening. The clusters are cut from the DRVI space itself, so a dimension that owns a cluster
is expected rather than a finding (every heatmap is also drawn against `cohort`). The reading that owes
nothing to any annotation is the OOD / IND interpretability score, which names a dimension by
its genes: it is what 03_6 and 03_7 use.

`cnv_score` and `cnv_corr` from 03_1 are also drawn on the latent UMAP: a gradient along a
dominant dimension would mean DRVI spent its capacity on CNV burden rather than on states.

## Outputs

    $DATA_DIR/03_tum/
    ├── model_drvi_tum_64_nomt.pt              # the trained model, one flat file
    ├── embed_drvi_tum_64_nomt.h5ad            # latent space, dimension stats (vanished), OOD/IND scores, Leiden columns: read by 03_4 on
    └── shiao_tum_drvi_tum_64_nomt.h5ad        # shiao_tum.h5ad + obsm['X_drvi'] + Leiden columns, for 03_5 and 03_6

Figures go to `../figures/03_3_drvi_tum_64_nomt/`, tables to `../tables/03_3_drvi_tum_64_nomt/`.
The enrichment section writes `factor_processes_drvi_tum_64_nomt_top<N>.csv`: per dimension
direction, its SMI match, the top hit of each tool and its top genes.

The cell cycle will enrich somewhere: a real axis, and also the confound 03_6 measures.
