# 03_7_factor_first — Route B

The dimension defines itself; prior knowledge only names it. The unit is the gene: which
program does each dimension-direction encode? The discovery route, and the check on Route A's
confirmation bias. It relies on DRVI's additive decoder, which gives every dimension a readable
gene-level footprint.

```bash
export DATA_DIR=<data root>
bash factor_first_all.sh                         # scie, emt and gavish, resuming
bash factor_first_all.sh --force                 # re-run everything
bash factor_first_all.sh --dry-run               # print what would run
```

Run after `../03_6_cell_first/cell_first_all.sh`: the wrapper stops if a collection's row order
is missing. One script at a time:

```bash
python3 factor_first_tum.py                      # scie, the default
python3 factor_first_tum.py --collection emt     # the EMT lists
python3 factor_first_tum.py --collection gavish  # the 22 TNBC-relevant metaprograms
python3 factor_first_tum.py --collection gavish --all-metaprograms  # all 41 instead
python3 factor_first_tum.py --n-top-genes 500    # a deeper list
```

Needs 03_4 (the `.gmt`), 03_3 (the embedding) and 03_6 (the row order). No GPU, no model: the
genes × direction table comes from `embed.varm`.

With `cell_type` constant, naming a dimension by its genes is the reading that owes nothing to a
borrowed vocabulary; the groupings of 03_3 and 03_6 are its cross-check.

## The method

- **Both directions are tested separately**: DRVI can encode two programs on the two sides of
  one axis.
- **Offline hypergeometric tests against the HVG background** (`gp.enrich`), not Enrichr's
  all-human-genes universe. DRVI trained on the 2,000 HVGs of 03_2, so each signature is tested
  in its HVG-restricted form (03_4's `n_in_hvg_background`).
- **BH once across every direction × gene-set pair** of one collection; pairs with no overlap
  are added back as p = 1. The overlapping lists of 03_4 are still not independent tests.

Vanished directions are pruned: 31 of 128 on `drvi_tum_64_nomt`, 97 tested. Route A keeps
`DR 48`, while this step drops `DR 48+`.

## What comes out

```
../tables/03_7_factor_first/dim_geneset_signed_significance_*.csv   the matrix 03_8 joins
../tables/03_7_factor_first/factor_first_significant_*.csv          every significant pair
$DATA_DIR/03_tum/factor_first_top<N>_<collection>_drvi_tum_64_nomt_pruned.tsv every pair tested
$DATA_DIR/03_tum/factor_first_top<N>_genes_drvi_tum_64_nomt_pruned.tsv        the top-gene lists
```

On `drvi_tum_64_nomt` (top 200, FDR < 0.05):

| collection | significant pairs | directions | dimensions | gene sets hit |
|---|---|---|---|---|
| `scie` | 119 | 53 of 97 | 45 | 10 of 10 |
| `emt` | 3 | 3 (`DR 11+`, `DR 17+`, `DR 26+`) | 3 | `EMT_B_MESENCHYMAL` only |
| `gavish_tnbc` | 368 | 88 of 97 | 49 | 21 of 22 (not `MP4_CHROMATIN`) |

Figures in `../figures/03_7_factor_first/`. The signed matrix uses Route A's row order; the sign is the direction carrying the enrichment.
The generic naming of the directions (GO terms, TFs) is 03_3's `factor_processes_*` table,
not repeated here.

**How it fails:** enrichment says nothing about cells. A dimension can be enriched for a list
while driven by a few cells, one patient or a stress response sharing the genes, and the top-gene
cut favours short sharp signatures over long diluted ones. Hence 03_8.
