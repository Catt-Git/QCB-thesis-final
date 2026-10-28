"""03_9: the sample, the readouts and the hypotheses of treatment_programs_tum.py.

One joint DRVI run, not one per timepoint: the model corrects the patient (batch_key='cohort'),
so the treatment effect stays in one shared coordinate system. The unit is the biopsy
(cohort x treatment), paired within the patient: 42,000 cells would be pseudoreplication.
The programs are the bars of 03_8's programs_named, not re-chosen here.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

STEP_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(STEP_DIR.parent / "utils"))
import signature_common as C          # noqa: E402
import sig_collections as SC          # noqa: E402

STEP = "03_9_treatment"

# Outputs belong to the run, not to a collection, like 03_8's barplot.
CONSENSUS = SimpleNamespace(name="consensus",
                            title="all three collections, read together")

TREATMENT_KEY = C.CS.TREATMENT_KEY    # 'treatment'
RESPONSE_KEY = "response"
BATCH_KEY = C.BATCH_KEY               # 'cohort' = the patient

# Shiao et al. 2024: before treatment, after one cycle of pembrolizumab, after pembrolizumab + RT.
TREATMENTS = ("BASE", "PD1", "RTPD1")
RESPONSES = ("NR", "R1", "R2")

# (later, earlier); RTPD1 vs PD1 isolates the radiotherapy.
CONTRASTS = (("PD1", "BASE"), ("RTPD1", "BASE"), ("RTPD1", "PD1"))


def contrast_name(later: str, earlier: str) -> str:
    return f"{later}_vs_{earlier}"


CONTRAST_NAMES = tuple(contrast_name(*c) for c in CONTRASTS)

# Biopsies under this many malignant cells are not summarised (--min-cells).
MIN_CELLS_PER_SAMPLE = 20

# Fewer pairs: not tested (with 5 pairs the smallest two-sided Wilcoxon p is 0.0625).
MIN_PAIRS = 6

# Categorical slots 1-3 of the reference palette, fixed per response group.
RESPONSE_COLORS = {"NR": "#2a78d6", "R1": "#eb6834", "R2": "#1baf7a"}

# Diverging blue <-> red with a grey midpoint, for paired deltas.
DIVERGING = ("#184f95", "#6da7ec", "#f0efec", "#e98584", "#b02a29")

# The hypotheses fixed BEFORE looking at the data
# One (program, contrast, expected sign) per row, in their own BH family; everything else
# is exploratory.
HYPOTHESES = pd.DataFrame([
    ("HALLMARK_IFNG", "PD1_vs_BASE", +1, "IFN-gamma response after anti-PD1"),
    ("MP17_INTERFERON_MHC_II_1", "PD1_vs_BASE", +1, "IFN / MHC-II after anti-PD1"),
    ("HALLMARK_IFNA", "RTPD1_vs_PD1", +1, "type-I IFN after RT"),
    ("MP19_EPITHELIAL_SENESCENCE", "RTPD1_vs_PD1", +1, "senescence after RT"),
    ("MP1_CELL_CYCLE_G2_M", "RTPD1_vs_PD1", -1, "G2/M down after RT"),
    ("MP2_CELL_CYCLE_G1_S", "RTPD1_vs_PD1", -1, "G1/S down after RT"),
], columns=["program", "contrast", "expected_sign", "rationale"])

# The embedding and the programs

def load_embedding():
    """03_3's embedding and the cells x dimensions frame, columns titled 'DR 1' ... (DR n is not column n)."""
    import anndata as ad

    if not C.EMBED_H5AD.exists():
        sys.exit(f"missing {C.EMBED_H5AD}: run 03_3_drvi_run/run_drvi_tum.py first")
    adata = ad.read_h5ad(C.EMBED_H5AD)
    for key in (BATCH_KEY, TREATMENT_KEY, RESPONSE_KEY):
        if key not in adata.obs:
            sys.exit(f"{C.EMBED_H5AD.name} has no obs['{key}']")
    L = pd.DataFrame(np.asarray(adata.X), index=adata.obs_names,
                     columns=adata.var["title"].values)
    print(f"[read] {C.EMBED_H5AD.name}: {adata.n_obs} cells x {adata.n_vars} dimensions")
    return adata, L


def read_programs() -> pd.DataFrame:
    """The bars of programs_named: one row per program, its directions as a list."""
    try:
        prog = C.read_table("program_dimensions", CONSENSUS)
    except FileNotFoundError:
        sys.exit(f"missing {C.table_path('program_dimensions', CONSENSUS)}: "
                 f"run 03_8_convergence/dr_naming_tum.py first")
    prog = prog.copy()
    prog["directions"] = [[d.strip() for d in str(s).split(",") if d.strip()]
                          for s in prog["dimensions"]]
    prog["collections"] = [[c.strip() for c in str(s).split(",") if c.strip()]
                           for s in prog["collections"]]
    print(f"[read] {len(prog)} programs on "
          f"{len({d for ds in prog['directions'] for d in ds})} directions "
          f"({C.table_path('program_dimensions', CONSENSUS).parent.name}/)")
    return prog


def direction_parts(direction: str) -> tuple[str, float]:
    """'DR 18-' -> ('DR 18', -1.0)."""
    m = re.fullmatch(r"\s*(DR \d+)\s*([+-])\s*", direction)
    if not m:
        raise ValueError(f"not a dimension-direction: {direction!r}")
    return m.group(1), (1.0 if m.group(2) == "+" else -1.0)


def oriented_directions(L: pd.DataFrame, directions) -> pd.DataFrame:
    """Cells x directions: the coordinate times the sign of the direction, z-scored over cells.

    Oriented so higher means more of the program; z-scored so paired deltas share one unit.
    """
    out = {}
    for d in directions:
        dim, sign = direction_parts(d)
        if dim not in L.columns:
            sys.exit(f"{d}: no dimension titled {dim!r} in the embedding")
        v = sign * L[dim].to_numpy(float)
        sd = v.std()
        out[d] = (v - v.mean()) / (sd if sd > 0 else 1.0)
    return pd.DataFrame(out, index=L.index)


def program_composites(Z: pd.DataFrame, prog: pd.DataFrame) -> pd.DataFrame:
    """Cells x programs: the mean of the program's oriented, z-scored directions (one readout each)."""
    return pd.DataFrame({p: Z[r["directions"]].mean(axis=1) for p, r in prog.iterrows()},
                        index=Z.index)

# Samples

def sample_ids(obs: pd.DataFrame) -> pd.Series:
    """`Patient01__BASE`: one biopsy. Double underscore, safe as an R column name."""
    return obs[BATCH_KEY].astype(str) + "__" + obs[TREATMENT_KEY].astype(str)


def sample_table(obs: pd.DataFrame, min_cells: int) -> pd.DataFrame:
    """One row per biopsy: patient, timepoint, response, malignant cells, kept or not."""
    s = (obs.assign(sample=sample_ids(obs))
         .groupby("sample", observed=True)
         .agg(cohort=(BATCH_KEY, "first"), treatment=(TREATMENT_KEY, "first"),
              response=(RESPONSE_KEY, "first"), n_cells=(BATCH_KEY, "size")))
    for key in ("cohort", "treatment", "response"):
        s[key] = s[key].astype(str)
    unknown = set(s["treatment"]) - set(TREATMENTS)
    if unknown:
        sys.exit(f"unexpected treatment level(s) {sorted(unknown)}; expected {TREATMENTS}")
    s["kept"] = s["n_cells"] >= min_cells
    return s.sort_values(["response", "cohort", "treatment"])


def count_matrix(samples: pd.DataFrame) -> pd.DataFrame:
    """Patients x timepoints, malignant cells per biopsy, response as the first column."""
    m = samples.pivot_table(index="cohort", columns="treatment", values="n_cells",
                            aggfunc="sum", fill_value=0)
    m = m.reindex(columns=list(TREATMENTS), fill_value=0)
    resp = samples.groupby("cohort")["response"].first()
    m.insert(0, "response", resp.loc[m.index])
    m = m.rename_axis("cohort").reset_index()
    return m.sort_values(["response", "cohort"]).set_index("cohort")


def report_samples(samples: pd.DataFrame, min_cells: int) -> None:
    kept = samples[samples["kept"]]
    dropped = samples[~samples["kept"]]
    print(f"{len(samples)} biopsies, {len(kept)} with >= {min_cells} malignant cells "
          f"({kept['n_cells'].sum()} of {samples['n_cells'].sum()} cells)")
    if len(dropped):
        print(f"dropped ({len(dropped)}): " + ", ".join(
            f"{r.cohort} {r.treatment} ({r.n_cells})" for r in dropped.itertuples()))
    per_patient = kept.groupby("cohort")["treatment"].apply(lambda t: "+".join(
        x for x in TREATMENTS if x in set(t)))
    print("timepoints per patient after the floor: " + ", ".join(
        f"{k} {v}" for k, v in per_patient.value_counts().items()))


def bh(p) -> np.ndarray:
    """Benjamini-Hochberg over the non-missing entries; NaN stays NaN."""
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.sum():
        from statsmodels.stats.multitest import multipletests
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def diverging_cmap():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("treatment_diverging", DIVERGING)


def collection_by_slug(slug: str):
    """A collection by the slug its outputs are written under - `gavish_tnbc` included."""
    for coll in (*SC.COLLECTIONS.values(), *SC.GAVISH_VARIANTS.values()):
        if coll.name == slug:
            return coll
    raise KeyError(f"no collection writes under the slug {slug!r}")

