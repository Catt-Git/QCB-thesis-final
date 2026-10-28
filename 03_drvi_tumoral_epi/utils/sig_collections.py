"""The signature collections this phase interprets, and everything that differs between them.

03_4 - 03_8 run one procedure (Route A on the cells, Route B on the decoder, Route C where
they converge) on each collection, selected with --collection:

  * scie   - stemness x immunogenicity: is there a stem-like AND immune-evasive state?
  * emt    - which cells sit in the hybrid, partial-EMT state?
  * gavish - Gavish 2023 pan-cancer metaprograms, a vocabulary to name the dimensions
             (22 TNBC-relevant by default, all 41 with --all-metaprograms).

Outputs and the BH correction are kept per collection, so collections never affect each other.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable


# The pieces a collection is made of

@dataclass(frozen=True)
class Signature:
    """One gene list. `file` is the name in $DATA_DIR/signatures/; `provenance` goes into the
    .gmt description field.

    `primary=False` marks a robustness check or a negative control rather than an independent test.
    """

    name: str
    file: str
    axis: str
    provenance: str
    primary: bool = True


@dataclass(frozen=True)
class Derived:
    """A readout computed as z(plus) - z(minus), after the within-stratum standardisation.

    A per-cell contrast, not a gene set, so it has no Route B counterpart.
    """

    name: str
    plus: str
    minus: str
    axis: str
    description: str


@dataclass(frozen=True)
class Plane:
    """One panel of the cell-first plane: one or two readouts and the rule that defines the target.

    `rule` is "high", "low" or "mid" (the middle band; no collection uses it now). With `y` None
    the definition is one readout alone. `in_consensus=False` is drawn and compared with the
    target, but never votes on it.
    """

    label: str
    x: str
    x_rule: str
    y: str | None = None
    y_rule: str | None = None
    in_consensus: bool = True

    @property
    def readouts(self) -> list[str]:
        return [self.x] + ([self.y] if self.y else [])

    def describe(self) -> str:
        two = f" x {self.y} {self.y_rule}" if self.y else ""
        return f"{self.x} {self.x_rule}{two}" + ("" if self.in_consensus else "  (shown, no vote)")


@dataclass(frozen=True)
class Criterion:
    """One half of the 03_8 target-axis test.

    Strongest Route A rho over `axis` or `names`, oriented by `sign` (+1, -1, or 0 for magnitude only).
    """

    label: str
    sign: int
    axis: str | None = None
    names: tuple[str, ...] = ()


@dataclass(frozen=True)
class Collection:
    name: str                       # the slug in every path and filename
    title: str                      # what it is called in headings
    question: str                   # the biological question, printed by every step
    signatures: tuple[Signature, ...]
    axes: tuple[str, ...]           # column / block order in every heatmap
    # Target region: a collection asking whether a named state exists (scie, emt) defines
    # the fields below; a vocabulary (gavish) leaves them at their defaults (has_target False).
    planes: Callable[["Collection", list[str]], list[Plane]] | None = None
    plane_figure: str | None = None  # figure basename for the cell-first plane
    target_label: str | None = None  # what the target region is called in prose
    risks: tuple[str, ...] = ()      # named checks 03_6 must run, see its header
    criteria: tuple[Criterion, ...] = ()  # the 03_8 target-axis test
    derived: tuple[Derived, ...] = ()
    extra_readouts: dict[str, str] = field(default_factory=dict)  # name -> axis, joined at runtime
    depth_risk_readout: str | None = None      # the readout whose LOW group could just be shallow
    ambient_risk_axis: str | None = None       # the axis whose HIGH group could just be ambient
    # (collection, claimed signature, raw 03_6 rho) -> extra flag names. Raw, not 03_8's
    # side-oriented a_rho, which is always positive on a claiming row.
    extra_flags: Callable[["Collection", str, float], list[str]] = lambda c, n, r: []

    # Sorts names inside an axis block for figures only; tables keep registry order.
    figure_key: Callable[[str], object] | None = None

    # accessors

    @property
    def has_target(self) -> bool:
        """Whether this collection defines a target region of the cell-first plane."""
        return self.planes is not None

    @property
    def short_title(self) -> str:
        """The name in the title's last parentheses (SCIE, EMT, ...), for figure titles."""
        if self.title.endswith(")") and "(" in self.title:
            return self.title[self.title.rindex("(") + 1:-1]
        return self.title

    @property
    def names(self) -> list[str]:
        return [s.name for s in self.signatures]

    @property
    def axis_of(self) -> dict[str, str]:
        """Every readout the collection can produce, including the derived and joined ones."""
        m = {s.name: s.axis for s in self.signatures}
        m.update({d.name: d.axis for d in self.derived})
        m.update(self.extra_readouts)
        return m

    @property
    def provenance(self) -> dict[str, str]:
        return {s.name: s.provenance for s in self.signatures}

    @property
    def files(self) -> dict[str, str]:
        return {s.name: s.file for s in self.signatures}

    def by_axis(self, axis: str) -> list[str]:
        """Every readout on one axis, registry order, derived and joined ones included."""
        return [n for n, a in self.axis_of.items() if a == axis]

    def primary_names(self) -> list[str]:
        return [s.name for s in self.signatures if s.primary]

    def order(self, present: list[str], for_figure: bool = False) -> list[str]:
        """Readouts ordered by axis, then by registry order (or `figure_key` for figures)."""
        rank = {n: i for i, n in enumerate(self.names + [d.name for d in self.derived]
                                          + list(self.extra_readouts))}
        axis_rank = {a: i for i, a in enumerate(self.axes)}
        key = self.figure_key if (for_figure and self.figure_key is not None) else None
        return sorted([n for n in present if n in self.axis_of],
                      key=lambda n: (axis_rank.get(self.axis_of[n], len(self.axes)),
                                     key(n) if key is not None else 0,
                                     rank.get(n, len(rank))))

    def block_edges(self, ordered: list[str]) -> list[int]:
        """Positions where the axis changes, i.e. where a heatmap needs a separator line."""
        edges, prev = [], None
        for i, n in enumerate(ordered):
            a = self.axis_of.get(n)
            if prev is not None and a != prev:
                edges.append(i)
            prev = a
        return edges


# scie - stemness x immunogenicity
#
# No stemness consensus: the lab's intersection of the ES lists captured only proliferation.

PRIMARY_IMMUNE = "IMMUNOGENIC_CONSENSUS"

_SCIE_SIGNATURES = (
    Signature("HALLMARK_IFNA",         "HALLMARK_IFNA.txt",      "immune",   "Interferon Alpha response, MSigDB Hallmark"),
    Signature("HALLMARK_IFNG",         "HALLMARK_IFNG.txt",      "immune",   "Interferon Gamma response, MSigDB Hallmark"),
    Signature("ISDS",                  "ISDS.txt",               "immune",   "IFN-Stem Cell-Down signature, PMC5481166"),
    Signature("KEGG_APM",              "KEGG_APM.txt",           "immune",   "Antigen Presentation Machinery, KEGG"),
    Signature("IMMUNOGENIC_CONSENSUS", "Immune consensus.txt",   "immune",   "Curated by the lab from the KEGG APM signature, retaining only immunogenic genes and immunomodulators; primary immune readout"),
    Signature("BENPORATH_ES1",         "BENPORATH_ES1.txt",      "stemness", "MSigDB BENPORATH_ES_1"),
    Signature("ESC_ASSOU",             "ESC_ASSOU.txt",          "stemness", "PMC1906587, Table S3"),
    Signature("ESC_WONG",              "ESC_WONG.txt",           "stemness", "MSigDB WONG_EMBRYONIC_STEM_CELL_CORE"),
    # EMP (PMID 29784918) was dropped: 92% of its score came from RPSA, FN1 and MFAP2, and its
    # progenitor markers are below the noise floor; it was never significant on either route.
    Signature("LIM_STEM",              "LIM_STEM.txt",           "stemness", "MSigDB LIM_MAMMARY_STEM_CELL_UP"),
    Signature("FMASC",                 "fMaSC.txt",              "stemness", "Fetal mammary stem cells, PMC3277444, Supplementary Table 2"),
)


def _scie_planes(coll: "Collection", available: list[str]) -> list[Plane]:
    """One panel per stemness readout, all against the primary immune readout.

    No stemness list is primary, so the stability across these planes is itself a result.
    """
    stem = coll.order([n for n in coll.by_axis("stemness") if n in available], for_figure=True)
    return [Plane(label=s, x=s, x_rule="high", y=PRIMARY_IMMUNE, y_rule="low") for s in stem]


def _scie_figure_key(name: str) -> int:
    # CytoTRACE2 (list-free) first, LIM_STEM (mammary lineage, unlike the ES lists) last.
    return {"CytoTRACE2": 0, "LIM_STEM": 2}.get(name, 1)


def _scie_flags(coll: "Collection", claimed: str, raw_rho: float) -> list[str]:
    # Immune evasion is the absence of a signal, which shallow sequencing mimics. raw_rho is
    # 03_6's unoriented rho, where the sign is the direction.
    if coll.axis_of.get(claimed) == "immune" and raw_rho < 0:
        return ["immune_low_is_absence_of_signal"]
    return []


SCIE = Collection(
    name="scie",
    title="stemness x immunogenicity (SCIE)",
    question="is there a malignant state that is stem-like AND immune-evasive?",
    signatures=_SCIE_SIGNATURES,
    axes=("immune", "stemness"),
    planes=_scie_planes,
    plane_figure="stemness_immunogenicity_plane",
    target_label="stem-high / immunogenic-low",
    risks=("cell_cycle", "depth"),
    criteria=(
        Criterion("stem_rho", sign=+1, axis="stemness"),
        # evasion = LOW immunogenicity
        Criterion("immunogenic_low_rho", sign=-1, names=(PRIMARY_IMMUNE,)),
    ),
    # CytoTRACE2 from 03_5, joined by 03_6 when available: the only list-free stemness readout.
    extra_readouts={"CytoTRACE2": "stemness"},
    depth_risk_readout=PRIMARY_IMMUNE,
    extra_flags=_scie_flags,
    figure_key=_scie_figure_key,
)


# emt - the epithelial-to-mesenchymal axis and its hybrid state
# The collaborator's list B only; the hybrid state is called by the hybrid list.

_EMT_SIGNATURES = (
    Signature("EMT_B_EPITHELIAL",  "EMT_B_EPITHELIAL.txt",  "epithelial",  "Collaborator list B (validated + non-validated), epithelial markers"),
    Signature("EMT_B_HYBRID",      "EMT_B_HYBRID.txt",      "hybrid",      "Collaborator list B, hybrid/partial-EMT markers"),
    Signature("EMT_B_MESENCHYMAL", "EMT_B_MESENCHYMAL.txt", "mesenchymal", "Collaborator list B, mesenchymal markers"),
)

_EMT_DERIVED = (
    Derived(
        name="EMT_SCORE_B",
        plus="EMT_B_MESENCHYMAL",
        minus="EMT_B_EPITHELIAL",
        axis="emt_score",
        description="z(EMT_B_MESENCHYMAL) - z(EMT_B_EPITHELIAL), the E-to-M position of a cell",
    ),
)


def _emt_planes(coll: "Collection", available: list[str]) -> list[Plane]:
    """The target is the hybrid list high; E-high AND M-high co-expression is shown beside it."""
    planes = []
    if "EMT_B_HYBRID" in available:
        planes.append(Plane(label="hybrid list", x="EMT_B_HYBRID", x_rule="high"))
    e, m = "EMT_B_EPITHELIAL", "EMT_B_MESENCHYMAL"
    if e in available and m in available and planes:
        planes.append(Plane(label="E x M co-expression", x=e, x_rule="high", y=m, y_rule="high",
                            in_consensus=False))
    return planes


def _emt_flags(coll: "Collection", claimed: str, raw_rho: float) -> list[str]:
    # High VIM / FN1 / SPARC can be fibroblast ambient RNA or a doublet: flagged on every
    # mesenchymal claim, and 03_6 runs the doublet check.
    if coll.axis_of.get(claimed) == "mesenchymal":
        return ["mesenchymal_may_be_ambient_or_doublet"]
    return []


EMT = Collection(
    name="emt",
    title="epithelial-mesenchymal transition (EMT)",
    question="which cells sit in the hybrid, partial-EMT state?",
    signatures=_EMT_SIGNATURES,
    axes=("epithelial", "hybrid", "mesenchymal", "emt_score"),
    planes=_emt_planes,
    plane_figure="emt_target_plane",
    target_label="hybrid state: the hybrid/partial-EMT list high",
    risks=("cell_cycle", "ambient"),
    criteria=(
        Criterion("hybrid_rho", sign=+1, axis="hybrid"),
        # Magnitude only: an E-to-M axis in either orientation.
        Criterion("emt_axis_rho", sign=0, axis="emt_score"),
    ),
    derived=_EMT_DERIVED,
    ambient_risk_axis="mesenchymal",
    extra_flags=_emt_flags,
)


# gavish - Gavish et al. 2023 pan-cancer metaprograms, a vocabulary with no target region
# Families are ours (Route C's same_family); lineage metaprograms are negative controls (primary=False).

_GAVISH_PAPER = "Gavish et al. 2023, Nature 618:598-606, pan-cancer malignant metaprogram"

# Families that cannot be present here: their metaprograms are negative controls.
GAVISH_CONTROL_AXES = ("neural", "other_lineage")


def _mp(stem: str, axis: str, label: str, note: str = "", name: str | None = None,
        source: str = "the published gene list as exported to datasets/GAVISH.csv") -> Signature:
    """One metaprogram; `stem` is the filename written by utils/gavish_extraction.py.

    `name` is only set for MP3, whose file keeps the CSV typo (CYLCE).
    """
    # "as exported": eight GAVISH.csv columns carry 48 or 49 genes, not 50.
    provenance = f"{_GAVISH_PAPER} MP{stem.split('_')[0][2:]} ({label}), {source}"
    if note:
        provenance += f"; {note}"
    return Signature(name=name or stem, file=f"GAVISH_metaprograms/{stem}.txt", axis=axis,
                     provenance=provenance, primary=axis not in GAVISH_CONTROL_AXES)


_CONTROL_NOTE = ("lineage absent from a breast epithelial compartment by construction - kept "
                 "as an internal negative control, its expected correlation is zero")

_GAVISH_SIGNATURES = (
    _mp("MP1_CELL_CYCLE_G2_M",           "cell_cycle",         "Cell Cycle - G2/M",
        source="the published gene list, taken from the authors' MP_list.RDS rather than from "
               "GAVISH.csv, which does not carry it - see the note above"),
    _mp("MP2_CELL_CYCLE_G1_S",           "cell_cycle",         "Cell Cycle - G1/S"),
    _mp("MP3_CELL_CYLCE_HMG_RICH",       "cell_cycle",         "Cell Cycle - HMG-rich",
        name="MP3_CELL_CYCLE_HMG_RICH"),
    _mp("MP4_CHROMATIN",                 "chromatin",          "Chromatin"),
    _mp("MP5_STRESS",                    "stress",             "Stress"),
    _mp("MP6_HYPOXIA",                   "stress",             "Hypoxia"),
    _mp("MP7_STRESS_IN_VITRO",           "stress",             "Stress (in vitro)",
        note="derived from cultured cells, so a high score on fresh tissue is a dissociation "
             "readout before it is a biological one"),
    _mp("MP8_PROTEASOMAL_DEGRADATION",   "proteostasis",       "Proteasomal degradation"),
    _mp("MP9_UNFOLDED_PROTEIN_RESPONSE", "proteostasis",       "Unfolded protein response"),
    _mp("MP10_PROTEIN_MATURATION",       "proteostasis",       "Protein maturation"),
    _mp("MP11_TRANSLATION_INITIATION",   "proteostasis",       "Translation initiation",
        note="ribosomal-protein heavy, so it tracks library complexity as readily as a "
             "program - read it next to rho_n_genes_by_counts in the confounder table"),
    _mp("MP12_EMT_1",                    "emt",                "EMT-I"),
    _mp("MP13_EMT_2",                    "emt",                "EMT-II"),
    _mp("MP14_EMT_3",                    "emt",                "EMT-III"),
    _mp("MP15_EMT_4",                    "emt",                "EMT-IV"),
    _mp("MP16_MES_GLIOMA",               "emt",                "MES (glioma)",
        note="the glioma mesenchymal program, on the EMT axis because it is the "
             "mesenchymal one, not because gliomas are expected here"),
    _mp("MP17_INTERFERON_MHC_II_1",      "immune",             "Interferon/MHC-II (I)"),
    _mp("MP18_INTERFERON_MHC_II_2",      "immune",             "Interferon/MHC-II (II)"),
    _mp("MP19_EPITHELIAL_SENESCENCE",    "senescence",         "Epithelial Senescence"),
    _mp("MP20_MYC",                      "metabolic",          "MYC"),
    _mp("MP21_RESPIRATION",              "metabolic",          "Respiration"),
    _mp("MP38_GLUTATHIONE",              "detoxification",     "Glutathione"),
    _mp("MP39_METAL_RESPONSE",           "detoxification",     "Metal-response"),
    _mp("MP22_SECRETED_1",               "secreted",           "Secreted-I"),
    _mp("MP23_SECRETED_2",               "secreted",           "Secreted-II"),
    _mp("MP24_CILIA",                    "cilia",              "Cilia"),
    _mp("MP25_ASTROCYTES",               "neural",             "Astrocytes",           _CONTROL_NOTE),
    _mp("MP26_NPC_GLIOMA",               "neural",             "NPC (glioma)",         _CONTROL_NOTE),
    _mp("MP27_OLIGO_PROGENITOR",         "neural",             "Oligo Progenitor",     _CONTROL_NOTE),
    _mp("MP28_OLIGO_NORMAL",             "neural",             "Oligo normal",         _CONTROL_NOTE),
    _mp("MP29_NPC_OPC",                  "neural",             "NPC/OPC",              _CONTROL_NOTE),
    _mp("MP30_PDAC_CLASSICAL",           "epithelial_lineage", "PDAC-classical",
        note="a carcinoma epithelial-identity program from another organ: not a control, "
             "an out-of-tissue epithelial reference"),
    _mp("MP31_ALVEOLAR",                 "epithelial_lineage", "Alveolar",
        note="lung epithelial identity; an out-of-tissue epithelial reference, as MP30"),
    _mp("MP40_PDAC_RELATED",             "epithelial_lineage", "PDAC-related",
        note="an out-of-tissue epithelial reference, as MP30"),
    _mp("MP32_SKIN_PIGMENTATION",        "other_lineage",      "Skin-pigmentation",    _CONTROL_NOTE),
    _mp("MP33_RBCS",                     "other_lineage",      "RBCs",                 _CONTROL_NOTE),
    _mp("MP34_PLATELET_ACTIVATION",      "other_lineage",      "Platelet activation",  _CONTROL_NOTE),
    _mp("MP35_HEMATO_RELATED_1",         "other_lineage",      "Hemato-related-I",     _CONTROL_NOTE),
    _mp("MP36_IG",                       "other_lineage",      "IG",
        note="immunoglobulin; the ambient-RNA readout of this dataset as much as a control - "
             "plasma-cell transcripts are the classic soup, so read it first"),
    _mp("MP37_HEMATO_RELATED_2",         "other_lineage",      "Hemato-related-II",    _CONTROL_NOTE),
    _mp("MP41_UNASSIGNED",               "unassigned",         "Unassigned",
        note="the paper's own residual metaprogram: it has no interpretation to lend, and it "
             "is here so that a dimension matching nothing else can match it instead"),
)

# Column order of every gavish heatmap, shared by both variants.
_GAVISH_AXES = ("cell_cycle", "chromatin", "stress", "proteostasis", "emt", "immune",
                "senescence", "metabolic", "detoxification", "secreted", "cilia", "neural",
                "epithelial_lineage", "other_lineage", "unassigned")


# The TNBC-relevant subset, the default: 21 states reported in TNBC plus MP7 (dissociation control).
# MP38 and MP11 fell below a matched random-list null; --all-metaprograms keeps the controls.

_GAVISH_TNBC_MPS = (
    # states a TNBC malignant cell can be in
    1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 39,
    # in-vitro stress: separates a stress dimension from a dissociation one
    7,
)


def _mp_number(sig: Signature) -> int:
    """The MP number, from the filename stem: MP3_CELL_CYLCE_HMG_RICH.txt -> 3."""
    return int(sig.file.rsplit("/", 1)[-1].split("_")[0][2:])


_GAVISH_BY_NUMBER = {_mp_number(s): s for s in _GAVISH_SIGNATURES}

# A typo in the list above would otherwise give a quietly shorter collection.
_unknown = sorted(set(_GAVISH_TNBC_MPS) - set(_GAVISH_BY_NUMBER))
if _unknown:
    raise ValueError(f"_GAVISH_TNBC_MPS names metaprograms that are not in the registry: "
                     f"{', '.join(f'MP{n}' for n in _unknown)}")

# Registry (family) order.
_GAVISH_TNBC_SIGNATURES = tuple(s for s in _GAVISH_SIGNATURES
                                if _mp_number(s) in set(_GAVISH_TNBC_MPS))


def _gavish_flags(coll: "Collection", claimed: str, raw_rho: float) -> list[str]:
    """Flags that would make a metaprogram match mean something other than it says."""
    axis = coll.axis_of.get(claimed)
    flags = []
    if axis == "emt":
        # Same risk as the emt collection: ambient fibroblast RNA or a doublet.
        flags.append("mesenchymal_may_be_ambient_or_doublet")
    if axis in GAVISH_CONTROL_AXES:
        # A best match on a lineage that cannot be here bounds every other match.
        flags.append("lineage_control_claim_bounds_the_rest")
    return flags


def _gavish_collection(name: str, scope: str, signatures: tuple[Signature, ...]) -> Collection:
    """One of the two metaprogram collections; they differ only in the lists they carry."""
    return Collection(
        name=name,
        title=f"pan-cancer malignant metaprograms (Gavish, {scope})",
        question=("which recurrent pan-cancer metaprogram, if any, does each latent dimension "
                  "carry?"),
        signatures=signatures,
        # Only the families this variant has, in the shared order.
        axes=tuple(a for a in _GAVISH_AXES if any(s.axis == a for s in signatures)),
        extra_flags=_gavish_flags,
    )


GAVISH = _gavish_collection("gavish", "all 41", _GAVISH_SIGNATURES)
GAVISH_TNBC = _gavish_collection("gavish_tnbc", "TNBC-relevant", _GAVISH_TNBC_SIGNATURES)


COLLECTIONS = {c.name: c for c in (SCIE, EMT, GAVISH)}
DEFAULT_COLLECTION = "scie"

# gavish_tnbc is not a --collection value: --all-metaprograms widens gavish to all 41.
GAVISH_VARIANTS = {False: GAVISH_TNBC, True: GAVISH}


def get(name: str) -> Collection:
    if name not in COLLECTIONS:
        raise KeyError(f"unknown collection {name!r}; have: {', '.join(COLLECTIONS)}")
    return COLLECTIONS[name]


def resolve(args) -> Collection:
    """The collection a step runs on: --collection, narrowed by --all-metaprograms."""
    coll = get(args.collection)
    if coll is GAVISH:
        return GAVISH_VARIANTS[bool(getattr(args, "all_metaprograms", False))]
    return coll


def add_argument(parser) -> None:
    """The --collection and --all-metaprograms flags, shared by every step."""
    parser.add_argument("--collection", choices=sorted(COLLECTIONS), default=DEFAULT_COLLECTION,
                        help=f"which signature collection to interpret (default {DEFAULT_COLLECTION})")
    parser.add_argument("--all-metaprograms", action="store_true",
                        help="with --collection gavish, score all 41 metaprograms instead of "
                             f"the {len(_GAVISH_TNBC_SIGNATURES)} TNBC-relevant ones; outputs "
                             "are written under the 'gavish' slug rather than 'gavish_tnbc'. "
                             "No effect on the other collections")
