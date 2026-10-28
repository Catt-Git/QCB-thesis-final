# 03_8_convergence — Route C

Where the two routes agree: the main result of the phase.

```bash
export DATA_DIR=<data root>
bash convergence_all.sh                        # 3 collections, null, SMI, naming, UMAPs, cluster DE/DR, resuming
bash convergence_all.sh --force                # re-run everything
bash convergence_all.sh --dry-run              # print what would run
```

One script at a time:

```bash
python3 convergence_tum.py                     # scie, the default
python3 convergence_tum.py --collection emt    # the EMT lists
python3 convergence_tum.py --collection gavish # the metaprograms: no target axis
```

Reads the tables of 03_6 and 03_7 (`../tables/03_6_cell_first/`, `../tables/03_7_factor_first/`)
and, for Route B per direction, 03_7's full ORA (`$DATA_DIR/03_tum/factor_first_top200_*_pruned.tsv`):
the signed matrix keeps only the stronger side of a pair. Always one collection, so a run only
joins a collection with itself.

| File | What it does |
|---|---|
| `convergence_all.sh` | runs the whole step: convergence on scie, emt, gavish, then null, SMI, naming, UMAPs, cluster DE and DR enrichment |
| `convergence_tum.py` | joins Route A and Route B, one verdict per dimension-direction |
| `route_a_null_tum.py` | the null the ρ bar is calibrated on |
| `signature_presence_tum.py` | library for the null: matched random gene sets (not a step) |
| `leiden_smi_tum.py` | which DR direction owns which DRVI-Leiden cluster (SMI) |
| `dr_naming_tum.py` | which directions get a name: `programs_named` (per program) and `dr_names` (per DR) |
| `dr_umap_tum.py` | each name on the UMAP: program score and dimension side by side |
| `leiden_states_tum.py` | what the DRVI-Leiden clusters contain (DE) and which directions they carry (DR enrichment) |

## Why agreement is the criterion

The two routes fail in different ways: Route A by coincidence among correlated per-cell scores,
Route B by gene-set overlap with no cellular counterpart. A dimension is unlikely to pass both
for the wrong reason, and nothing is promoted on a single route.

| `verdict` | | reading |
|---|---|---|
| `factor_only` | B, not A | the axis carries the genes, Route A did not clear the bar. `A_vs_null` says which case: `below` the noise floor (no coherent cells, possibly patient or technical) or `above` it (a real but weak association) |
| `cell_only` | A, not B | the cells separate, but the program is not encoded cleanly on one axis |
| `convergent` | A and B, same family | a name |
| `both_different_family` | A and B, different family | a signal, not a name; mostly on `gavish_tnbc`, the only collection with ten families |
| `neither` | | |

## How the routes are joined

Route A is per dimension, so its direction is the sign of ρ (ρ > 0 on DR 7 is `DR 7+`); Route B
is already per direction. One row per live dimension and direction, the dimensions taken from
03_6's `dimension_row_order`: 2 × 49 = 98 rows on `drvi_tum_64_nomt`.

A readout is flagged as confounded when its raw score correlates with depth or cycle above 0.30
(`DEPTH_FLAG`, `CYCLE_FLAG`, conventions). A convergent dimension with a flag is a lead, not a
result.

## What comes out

```
../tables/03_8_convergence/convergence_*.csv    one row per dimension-direction: both routes, the flags
../tables/03_8_convergence/target_axes_*.csv    the axes that pass the collection's criteria, if any
../figures/03_8_convergence/convergence_scatter_*.png   Route A vs Route B, coloured by verdict
../figures/03_8_convergence/routes_side_by_side_*.png   the 03_6 and 03_7 heatmaps, same row order
```

`target_axes` answers each collection's question ("a stem-like, immune-evasive state?" for
`scie`, "the hybrid E/M cells?" for `emt`) and is written only when an axis passes both
criteria: on `drvi_tum_64_nomt` none does, for either. `gavish` has no target, so it writes no
`target_axes` and its `convergence` table lacks `A_auroc_target_this_side` and
`A_standardised_mean_difference`.

## Naming: `dr_naming_tum.py`

A direction is named when Route A ρ ≥ 0.30 and Route B FDR < 0.05 on the same signature family:
a yes or a no. The magnitude of ρ is reported (`A_rho`, `max_rho`), never used to grade.

0.30 (`ROUTE_A_RHO_MIN` in `utils/signature_common.py`) against the per-direction null: 4.5% of
null directions reach it on `scie` (null p97.5 = 0.37), 3.5% on `gavish_tnbc`, 0.2% on `emt`,
i.e. ~4, ~3 and ~0 of 98 directions by chance. On `drvi_tum_64_nomt`: **31 directions on 29
dimensions** are named; one is a depth axis (`DR 1+`, `MP12_EMT_1`), the rest give **22 programs,
38 bars on 30 directions of 28 dimensions** (8 directions named by two collections); 31 = these
30 + the depth axis `DR 1+`. Writes
`dr_naming_consensus_*` and `program_dimensions_consensus_*` (the bars, read by `dr_umap_tum.py`
and 03_9) in `../tables/03_8_convergence/`, and three figures:
`programs_named_*` (per program: its directions and max ρ, plus in orange the 7 directions tied to
a DRVI-Leiden cluster by SMI), `programs_named_no_rho_*` (the same without ρ and SMI scores) and `dr_names_*`
(per DR: the name of its + and − side, 31 of 97 live directions named). The orange part also
goes to `states_vs_names_consensus_*`: none of the 7 directions carries a name. `--no-states` drops it,
`--all-states` draws every SMI match ≥ 0.7 instead of `STATE_PICKS`.

## The picture of a name — `dr_umap_tum.py`

```bash
python3 dr_umap_tum.py                       # all 38 pairs
python3 dr_umap_tum.py --only "DR 2+,DR 15-" # a few
python3 dr_umap_tum.py --score raw           # not the within-cohort z
```

One figure per bar, in `../figures/03_8_convergence/dr_umap/`, on 03_3's UMAP: left the within-cohort z of the program (purple), right
the coordinate on the dimension (red for `+`, blue for `-`, negated for `-`). Both are one-sided:
grey at or below zero, a ramp above it clipped at p99. The pairs are read from
`program_dimensions_consensus_drvi_tum_64_nomt.csv`, ρ and FDR from the 03_8 tables; only the
Spearman of the two plotted vectors is computed here.

A bar is labelled with Route B's best signature, while it only asks for the same family: 10 of
the 38 bars are labelled with a signature that is not Route A's best, and on two the label's
own ρ is under the bar (`ESC_ASSOU` on `DR 6+`, `MP13_EMT_2` on `DR 28+`).

## The null of the ρ bar — `route_a_null_tum.py`

Each gene set is replaced by size- and expression-matched random sets, scored and z-scored as in
03_6; the best set per direction gives the null of `A_rho`.

```bash
python3 route_a_null_tum.py                  # all three collections
python3 route_a_null_tum.py --collection emt
```

Random scores are cached in `$DATA_DIR/03_tum/`. Output:
`../tables/03_8_convergence/route_a_null_*.csv`.

## The DRVI-Leiden clusters — `leiden_smi_tum.py`

Which DR direction owns which cluster: SMI of every live direction against 03_3's Leiden
partition at 1.5, threshold 0.7, as in 03_3 (directions clipped at zero, vanished ones
dropped). Tables only, the figures are 03_3's: `smi_matrix_drvi_leiden_1.5_*` and
`smi_matches_drvi_leiden_1.5_*` (7 pairs, all one to one).

`STATE_PICKS` (`dr_naming_tum.py`) holds those 7 pairs: `dr_naming_tum.py` lays them over the
names, `leiden_states_tum.py` characterises them. Descriptive and circular: the clusters are
cut from the space they are described in.

## What the clusters are — `leiden_states_tum.py`

```bash
python3 leiden_states_tum.py --states-de     # the STATE_PICKS clusters: DE, composition, ORA
python3 leiden_states_tum.py --states-dr     # DR enrichment per cluster + summary (after --states-de)
```

`--states-de`: Wilcoxon one-vs-rest on all genes expressed in ≥ 1% of cells (every cluster
tested, the 7 `STATE_PICKS` characterised), then, per picked cluster: are the top 50 markers up
within each patient, the composition (patient, treatment, response, phase, depth, mito,
doublet, CNV), the overlap of the markers with the decoder's top 200 genes of each direction
(hypergeometric), and ORA against scie, emt and gavish_tnbc. Figures in
`../figures/03_8_convergence/state_de_drvi_leiden_1.5/`: `dotplot_markers_*` (top 6 markers per
picked cluster) and `composition_*` (patient and treatment shares per cluster, all cells for
scale). The decoder overlap is a table only (`state_marker_vs_decoder_*`).

`--states-dr`: every cluster, the AUROC of each live direction (cluster vs rest, also within
patient), what each direction is (03_8 name, Route B hits, 03_3 GO term, top genes), and one
row per cluster with `CLUSTER_READING` (34 clusters, read by hand from these outputs, draft of
2026-10-27). Figures in `state_dr_enrichment_drvi_leiden_1.5/`: `composition_all_clusters_*`
and `cluster_factors_barplot_*`; the AUROC is a table only (`state_dr_enrichment_*`), the
cluster × direction heatmap is 03_3's.

`STATE_PICKS` and `CLUSTER_READING` belong to one partition: re-curate them after a re-run of
03_3. The same DE and DR enrichment on the max(DR) groups is `../03_8b_argmax/`.
