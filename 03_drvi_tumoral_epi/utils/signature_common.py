"""Shared helpers for the signature-interpretation steps of the malignant phase (03_4 - 03_9).

Paths, table/figure writers and the interpretability-score
accessor for the run drvi_tum_64_nomt. The signature registries are in sig_collections.py.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

import cell_set as CS  # noqa: E402

# The run

N_LATENT = CS.n_latent()
RUN_ID = CS.run_id()

SEED = 0
SCORE_KEY = "OOD_combined"          # 03_3's interpretability scores

# Vanished dimensions (directions in 03_7) are always dropped, as DRVI does.
PRUNE_VANISHED = True

# Paths

PROJECT_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = CS.data_dir()

TUM_DIR = CS.tum_dir()                        # $DATA_DIR/03_tum, the heavy objects
# Read-only gene lists, one symbol per line.
SIG_DIR = Path(os.environ.get("SIGNATURE_DIR", DATA_DIR / "signatures"))

PHASE_DIR = PROJECT_DIR / "03_drvi_tumoral_epi"

# Phase-level tables/ and figures/, since later steps read earlier ones' tables.
# $OUT_ROOT moves both out of the repo, for a run on a different object.
OUT_ROOT = Path(os.environ.get("OUT_ROOT", "").strip() or PHASE_DIR).expanduser().resolve()
TABLE_DIR = OUT_ROOT / "tables"

FULL_H5AD = CS.path(".h5ad")                  # 03_2: all genes, log-normalised
# The panel DRVI was trained on: the Route B ORA background.
HVG_H5AD = CS.hvg_path(".h5ad")
EMBED_H5AD = TUM_DIR / f"embed_{RUN_ID}.h5ad" # 03_3: latent space + scores

# Groupings: cell_type is constant here, so tables group by Leiden, then by the
# pre-CNV label (a landmark only). Neither is an input to DRVI or Route A.
LABEL_KEY = CS.LABEL_KEY              # constant 'malignant'
GROUP_KEY = CS.PRIOR_LABEL_KEY        # 'cell_type_01_4'
LEIDEN_KEY = CS.LEIDEN_KEY            # 'drvi_leiden'
BATCH_KEY = CS.BATCH_KEY


def attach_leiden(obs: pd.DataFrame) -> pd.DataFrame:
    """Add LEIDEN_KEY to `obs` from the embedding, by cell name (shiao_tum.h5ad does not have it)."""
    if LEIDEN_KEY not in obs and EMBED_H5AD.exists():
        import h5py
        from anndata.io import read_elem
        with h5py.File(EMBED_H5AD, "r") as f:
            emb = read_elem(f["obs"])
        if LEIDEN_KEY in emb:
            obs[LEIDEN_KEY] = emb[LEIDEN_KEY].reindex(obs.index.astype(str)).values
    return obs


def grouping_keys(obs) -> list[str]:
    """The non-constant grouping keys of `obs`, best first."""
    out = []
    for key in (LEIDEN_KEY, GROUP_KEY, LABEL_KEY):
        if key in obs and obs[key].astype(str).nunique() > 1:
            out.append(key)
    return out

# Outputs carry their collection in the file name, so collections never overwrite each other.
CYTOTRACE_CSV = TUM_DIR / f"cytotrace2_{RUN_ID}.csv"   # per cell, 03_5


def gmt_path(coll) -> Path:
    """The collection as used (mapped genes only); depends on the object, not on the run."""
    return table_dir("03_4_signatures") / f"signatures_{coll.name}_{CS.compartment()}.gmt"


def scores_csv(coll) -> Path:
    """Per-cell raw + within-stratum z scores, 03_6. Heavy, so it lives outside the repo."""
    return TUM_DIR / f"signature_scores_{coll.name}_{RUN_ID}.csv"


def enrichment_tsv(coll, n_top: int) -> Path:
    """The full ORA output of 03_7, before any significance filter."""
    return TUM_DIR / f"factor_first_top{n_top}_{coll.name}_{RUN_ID}_pruned.tsv"


def factor_processes_csv(top_n: int = 100) -> Path:
    """03_3's per-direction naming (Enrichr and g:Profiler GO terms, decoupler TF), top_n genes."""
    return TABLE_DIR / f"03_3_{RUN_ID}" / f"factor_processes_{RUN_ID}_top{top_n}.csv"


def top_genes_tsv(n_top: int) -> Path:
    """The top genes per dimension-direction; not collection-scoped, every collection uses it."""
    return TUM_DIR / f"factor_first_top{n_top}_genes_{RUN_ID}_pruned.tsv"

# Floors shared by every collection

# Route A bar on |rho|, used by 03_6 and 03_8; calibrated by 03_8/route_a_null_tum.py
# (p97.5 of a direction's best random list: 0.374 scie, 0.202 emt, 0.322 gavish_tnbc on drvi_tum_64_nomt).
ROUTE_A_RHO_MIN = 0.30

MIN_SIGNATURE_GENES = 10      # below this a signature is skipped and reported
MIN_MAPPED_FRACTION = 0.60    # below this the step stops: low coverage means NOT MEASURED


# Figure sizing, so wide collections (gavish) stay readable; scie and emt are unaffected.

MAX_FIG_IN = 26.0
ANNOT_MAX_CELLS = 400    # above ~20 x 20 the printed numbers are smaller than the cells


def fig_span(n: int, per_item: float, base: float, cap: float = MAX_FIG_IN) -> float:
    """Inches for `n` rows or columns at `per_item` each, plus `base`, capped at `cap`."""
    return min(base + per_item * n, cap)


def annotate_cells(n_rows: int, n_cols: int) -> bool:
    """Whether a heatmap that size can carry its numbers legibly. Under 400 cells, yes."""
    return n_rows * n_cols <= ANNOT_MAX_CELLS


def axis_palette(n: int):
    """`n` distinguishable colours; 'deep' recycles past ten, so 'husl' above that."""
    import seaborn as _sns
    return _sns.color_palette("deep" if n <= 10 else "husl", n)

# The target region and the vote over its definitions, shared so later steps rebuild
# the same consensus cell set.


@dataclass(frozen=True)
class Cutoffs:
    """Quantiles defining 'high', 'low' and the middle band (fixed z would give uneven groups)."""

    high_q: float = 0.75
    low_q: float = 0.25
    mid_lo_q: float = 0.40
    mid_hi_q: float = 0.60

    @classmethod
    def from_args(cls, args) -> "Cutoffs":
        return cls(args.high_q, args.low_q, args.mid_lo_q, args.mid_hi_q)


def add_cutoff_arguments(parser) -> None:
    """The four quantile flags of 03_6, shared so every step defines the same region."""
    d = Cutoffs()
    parser.add_argument("--high-q", type=float, default=d.high_q,
                        help=f"quantile of the within-stratum z-score above which a cell is "
                             f"'high' (default {d.high_q})")
    parser.add_argument("--low-q", type=float, default=d.low_q,
                        help=f"quantile below which a cell is 'low' (default {d.low_q})")
    parser.add_argument("--mid-lo-q", type=float, default=d.mid_lo_q,
                        help=f"lower edge of the 'mid' band, unused by every collection now "
                             f"(default {d.mid_lo_q})")
    parser.add_argument("--mid-hi-q", type=float, default=d.mid_hi_q,
                        help=f"upper edge of the 'mid' band (default {d.mid_hi_q})")


def in_region(v: pd.Series, rule: str, cut: Cutoffs) -> pd.Series:
    """Whether each cell is at the high end, the low end, or in the middle band ("mid", unused now)."""
    if rule == "high":
        return v >= v.quantile(cut.high_q)
    if rule == "low":
        return v <= v.quantile(cut.low_q)
    if rule == "mid":
        return (v >= v.quantile(cut.mid_lo_q)) & (v <= v.quantile(cut.mid_hi_q))
    raise ValueError(f"unknown region rule {rule!r}")


def define_target(z: pd.DataFrame, plane, cut: Cutoffs) -> pd.Series:
    """The target cell set of one plane: its readouts inside their respective regions."""
    called = in_region(z[f"z_{plane.x}"], plane.x_rule, cut)
    if plane.y:
        called &= in_region(z[f"z_{plane.y}"], plane.y_rule, cut)
    return called


def consensus_vote(qdf: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    """(votes per cell, called by >= ceil(n/2) definitions, never fewer than 2 when n >= 2)."""
    n_defs = qdf.sum(axis=1)
    return n_defs, n_defs >= max(min(2, qdf.shape[1]), int(np.ceil(qdf.shape[1] / 2)))


# Signature IO


def read_signature_file(path: Path) -> list[str]:
    """One gene symbol per line, de-duplicated, order preserved (handles CRLF and BOM)."""
    seen, genes = set(), []
    with open(path, encoding="utf-8-sig") as fh:
        for line in fh:
            g = line.strip().strip('"').strip()
            if g and g not in seen:
                seen.add(g)
                genes.append(g)
    return genes


def load_signatures(coll, sig_dir: Path = SIG_DIR) -> dict[str, list[str]]:
    """One collection, as {name: [genes]}, straight from the text files, in registry order."""
    out = {}
    for name, fname in coll.files.items():
        path = sig_dir / fname
        if not path.exists():
            raise FileNotFoundError(f"signature file missing: {path}")
        genes = read_signature_file(path)
        if not genes:
            raise ValueError(f"signature file is empty: {path}")
        out[name] = genes
    return out


def write_gmt(sets: dict[str, list[str]], path: Path, descriptions: dict[str, str]) -> None:
    """GMT: name <tab> description (provenance) <tab> gene <tab> gene ..."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for name, genes in sets.items():
            fh.write("\t".join([name, descriptions.get(name, ""), *genes]) + "\n")


def read_gmt(path: Path) -> dict[str, list[str]]:
    """Back from GMT, dropping the description column."""
    sets = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 2:
                sets[parts[0]] = parts[2:]
    return sets


# Table and figure writers


# One folder per step: a table lives with the step that writes it; collection and run are
# in the file name. Matched by name prefix, so readers find it without knowing the writer.
TABLE_STEPS = {
    "03_4_signatures": ("coverage", "jaccard", "shared_genes", "emt_vs_gavish_pairs",
                        "signature_concentration", "signature_gene_contribution"),
    "03_5_cytotrace2": ("cytotrace2_by_",),
    "03_6_cell_first": ("confounder", "quadrant_", "dim_signature_spearman",
                        "dim_target_effect_size", "dimension_row_order"),
    "03_7_factor_first": ("dim_geneset_signed_significance", "factor_first_"),
    "03_8_convergence": ("convergence", "target_axes", "route_a_null", "dr_naming",
                         "program_dimensions", "states_vs_names"),
    "03_9_treatment": ("treatment_",),
}


def table_step(name: str) -> str:
    """The step folder of table `name` (longest matching prefix in TABLE_STEPS)."""
    hits = [(len(p), step) for step, prefixes in TABLE_STEPS.items()
            for p in prefixes if name.startswith(p)]
    if not hits:
        raise KeyError(f"table {name!r} has no step in TABLE_STEPS")
    return max(hits)[1]


def table_dir(step: str) -> Path:
    """tables/<step>/, created on demand."""
    d = TABLE_DIR / step
    d.mkdir(parents=True, exist_ok=True)
    return d


def table_path(name: str, coll, run_id: str | None = None) -> Path:
    """`coll` is a collection or its name; `run_id` defaults to the module's current one."""
    run_id = run_id or RUN_ID
    cname = coll if isinstance(coll, str) else coll.name
    return table_dir(table_step(name)) / f"{name}_{cname}_{run_id}.csv"


def write_table(df: pd.DataFrame, name: str, coll, index: bool = True,
                run_id: str | None = None) -> Path:
    """Write a result table under a one-line provenance comment (read back with comment='#')."""
    run_id = run_id or RUN_ID
    path = table_path(name, coll, run_id)
    pruning = "vanished dimensions PRUNED"
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"# {name} | collection {coll.name} ({coll.title}) "
                 f"| DRVI run {run_id} "
                 f"| {pruning} "
                 f"| 03_4 - 03_8 signature interpretation\n")
        df.to_csv(fh, index=index)
    print(f"[table] {path}  ({df.shape[0]} x {df.shape[1]})")
    return path


def read_table(name: str, coll, run_id: str | None = None) -> pd.DataFrame:
    """Read back a table this stage wrote, dropping the provenance comment line."""
    return pd.read_csv(table_path(name, coll, run_id), comment="#", index_col=0)


def fig_dir(step: str) -> Path:
    """figures/<step>/, `step` being e.g. '03_6_cell_first' or '03_8_convergence/dr_umap'."""
    d = OUT_ROOT / "figures" / step
    d.mkdir(parents=True, exist_ok=True)
    return d


def savefig(name: str, step: str, coll, fig=None, dpi: int = 300, run_id: str | None = None):
    """Save a figure into fig_dir(step), collection and run id appended to the name."""
    import matplotlib.pyplot as plt

    fig = plt.gcf() if fig is None else fig
    run_id = run_id or RUN_ID
    path = fig_dir(step) / f"{name}_{coll.name}_{run_id}.png"
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    print(f"[fig] {path}")
    return path


def dim_slug(dim: str) -> str:
    """'DR 7+' -> 'DR_07_pos'. Same convention as 03_3, so the two steps sort alike."""
    number, sign = dim[3:-1], dim[-1]
    return f"DR_{int(number):02d}_{'pos' if sign == '+' else 'neg'}"


def dim_sort_key(dim: str) -> tuple[int, int]:
    """Sort 'DR 7+' / 'DR 7-' numerically, positive direction first."""
    return int(dim[3:-1]), 0 if dim.endswith("+") else 1


# DRVI interpretability scores


def interpretability_scores(embed, gene_names, key: str = SCORE_KEY,
                            hide_vanished: bool = PRUNE_VANISHED):
    """Genes x dimension-directions scores from `embed.varm`, as model.get_interpretability_scores.

    Rebuilt from the embedding alone (no model). Vanished directions are dropped one by one,
    a finer cut than analysis_dimensions.
    """
    effect = np.concatenate([embed.varm[f"{key}_positive"], embed.varm[f"{key}_negative"]])

    info = (
        pd.concat([embed.var.assign(direction="+"), embed.var.assign(direction="-")])
        .assign(title=lambda df: df["title"] + df["direction"])
        .reset_index(drop=True)
    )
    info["keep"] = ~np.where(
        info["direction"] == "+",
        info["vanished_positive_direction"],
        info["vanished_negative_direction"],
    ) if hide_vanished else True

    return (
        pd.DataFrame(effect, columns=gene_names, index=info["title"])
        .loc[info.query("keep == True").sort_values(["order", "direction"])["title"]]
        .T
    )


def analysis_dimensions(embed) -> list[str]:
    """The non-vanished dimensions, in the embedding's order; the only place Route A drops one."""
    if PRUNE_VANISHED and "vanished" in embed.var:
        return embed.var.loc[~embed.var["vanished"].astype(bool), "title"].tolist()
    return embed.var["title"].tolist()


def n_vanished(embed) -> int:
    """How many dimensions DRVI flagged as vanished; 0 if the embedding has no such column."""
    if "vanished" not in embed.var:
        return 0
    return int(embed.var["vanished"].astype(bool).sum())


def banner(title: str) -> None:
    print(f"\n{'=' * 78}\n{title}\n{'=' * 78}")
