# 03_6_cell_first — Route A

Prior knowledge defines the state, DRVI is the coordinate system. The unit is the cell: do the
cells prior knowledge calls stem-like, immune-evasive or hybrid-EMT sit at a distinct position
along a latent dimension? The only route that assigns cells to states.

```bash
export DATA_DIR=<data root>
bash cell_first_all.sh                               # scie, emt and gavish, resuming
bash cell_first_all.sh --force                       # re-run everything (re-scores too)
bash cell_first_all.sh --dry-run                     # print what would run
```

The wrapper runs the collections one after the other (each loads `shiao_tum.h5ad` whole);
`--slurm` submits it as one job. One script at a time:

```bash
python3 cell_first_tum.py                            # scie, the default
python3 cell_first_tum.py --collection emt           # the EMT lists
python3 cell_first_tum.py --collection gavish        # the metaprograms: scores only, no target
python3 cell_first_tum.py --collection gavish --all-metaprograms   # all 41, not just the TNBC 22
python3 cell_first_tum.py --high-q 0.80 --low-q 0.20 # a stricter target region
python3 cell_first_tum.py --overwrite                # re-score instead of reusing the csv
```

Needs 03_4 (the `.gmt`) and 03_3 (the embedding); 03_5 is optional. Scoring is on the
unintegrated all-genes object (`shiao_tum.h5ad`, 41,943 × 25,197): a signature cut down to the
HVGs is no longer that signature.

| File | What it does |
|---|---|
| `cell_first_all.sh` | runs `cell_first_tum.py` on scie, emt and gavish, then 03_4's `signature_composition_tum.py` on the three (it reads these scores) |
| `cell_first_tum.py` | scores, confounder checks, target region, dimensions × readouts |

## Standardisation within `cohort` only

The groupings available inside the malignant compartment (`cell_type_01_4`, `drvi_leiden`) are
states of the same tumor: standardising within them would remove the contrast being measured
(`Lumsec-prol`, proliferating, is 53% of the cells). Leiden is reported as a covariate of the
target set (`quadrant_per_drvi_leiden`) instead.

## How it fails, and the checks

| failure | the check |
|---|---|
| "stem-high" can just mean "cycling" | the target region is recomputed inside G1 alone and compared by Jaccard, for every collection |
| immune evasion is an absence of signal, which shallow sequencing mimics | depth of the immunogenic-low group vs the rest: Mann-Whitney + AUROC |
| a high mesenchymal score can be ambient RNA or a doublet | `doublet_score` of the mesenchymal-high group, AUROC and Spearman |
| absolute scores differ across patients | standardisation within `cohort` |
| confirmatory by construction | none inside Route A: that is why Route B exists |

Being aneuploid rules out an actual fibroblast in the object, not ambient RNA.

## The target region

For `scie` the region is defined once per stemness readout, and the stability of the cell set
across definitions is itself a result (`quadrant_stability`, Jaccard); a consensus cell is called
by a majority of definitions, at least 2. The EMT target is the hybrid/partial-EMT list alone
(`EMT_B_HYBRID` ≥ q0.75). The co-expression corner (`EMT_B_EPITHELIAL` and `EMT_B_MESENCHYMAL`
both ≥ q0.75) is drawn beside it and compared with it (`quadrant_stability_emt`), but does not
vote: the two share few cells (Jaccard 0.08).

## What comes out

Tables in `../tables/03_6_cell_first/`, every name ending in `_<collection>_drvi_tum_64_nomt.csv`.
Embedding-independent:

```
confounders_*                   readouts vs depth, mito, S, G2M
quadrant_stability_*            Jaccard across definitions (emt: hybrid vs co-expression)
quadrant_per_patient_*          target size per patient
quadrant_per_drvi_leiden_*      target share per Leiden cluster
confounder_checks_*             cell cycle (G1), depth, doublet checks on the target
$DATA_DIR/03_tum/signature_scores_<collection>_drvi_tum_64_nomt.csv   per cell, raw + z
```

Embedding-dependent:

```
dim_signature_spearman_*        dimensions × readouts (input of 03_8)
dim_target_effect_size_*        AUROC and SMD of the target per dimension
dimension_row_order_*           the row order of every later heatmap
```

On `gavish` (no target region) only `confounders`, `signature_scores`,
`dim_signature_spearman` and `dimension_row_order` are written; the script prints what it skips.

Figures in `../figures/03_6_cell_first/`: `dim_signature_heatmap_*` (every collection), the
target plane with its cutoffs (`stemness_immunogenicity_plane_scie`, `emt_target_plane_emt`),
`stemness_immunogenicity_plane_mean_scie` (the six scie planes in one panel: x the mean of
the six stemness z-scores, red the cells called by at least 4 of 6 definitions, stricter than
the consensus) and `quadrant_stability_scie`. In the heatmaps, one row per dimension: the direction is the sign of ρ,
written on the colorbar, and is what makes Route A joinable with Route B in 03_8. Vanished
dimensions are pruned (`PRUNE_VANISHED` in `utils/signature_common.py`).

## On `drvi_tum_64_nomt`

49 live dimensions. Read per direction as 03_8 joins them, |ρ| ≥ 0.30 is reached by 16
directions (14 dimensions) on `scie`, 1 on `emt` and 34 (30 dimensions) on `gavish_tnbc`.

- **scie target**: 2,704 consensus cells (6.4%). The stemness definitions agree only partly
  (Jaccard 0.06 - 0.61, `LIM_STEM` the outlier at ≤ 0.22). The target is mostly cycling: 14.7%
  of it is in G1 against 47.0% of the rest, and recomputed inside G1 it keeps Jaccard 0.35.
  The immunogenic-low cells are shallower (median 1,128 vs 1,488 detected genes, AUROC 0.56).
- **emt target**: 10,486 hybrid-high cells. Not a cycle effect (G1 Jaccard 0.96) and barely
  a doublet one (AUROC 0.56, ρ 0.22 with `doublet_score`).

### The derived EMT axis

`EMT_SCORE_B = z(EMT_B_MESENCHYMAL) − z(EMT_B_EPITHELIAL)` is a coordinate along two programs,
not a fourth one: it is in `dim_signature_spearman_emt` for 03_8, which tests it by magnitude,
but not drawn in the heatmap.
