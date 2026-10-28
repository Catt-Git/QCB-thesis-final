#!/usr/bin/env bash
#
# 03_6 cell-first (Route A): every collection, in order:
#
#   1. cell_first_tum.py --collection scie    -> scores + confounder / target / row-order tables
#   2. cell_first_tum.py --collection emt     -> same, on the EMT lists
#   3. cell_first_tum.py --collection gavish  -> scores, spearman and row order only (no target)
#   4-6. ../03_4_signatures/signature_composition_tum.py on scie, emt, gavish
#                                             -> signature_concentration / gene_contribution (03_4)
#
# Needs 03_4 (the .gmt) and 03_3 (the embedding); 03_5 is optional (CytoTRACE2 joins scie only
# if its csv exists). The collections run one after the other: each loads shiao_tum.h5ad whole.
# Steps with an existing output are skipped, so a re-run resumes; --force re-runs everything
# (and re-scores, --overwrite). --slurm submits it as one job (utils/slurm.sh).
#
# Usage:
#   export DATA_DIR=<data root>
#   ./cell_first_all.sh                       # every step, resuming
#   ./cell_first_all.sh --slurm               # same chain as one job on the cluster
#   ./cell_first_all.sh --force               # re-run everything, overwriting
#   ./cell_first_all.sh --dry-run             # print what would run, do nothing
#   ./cell_first_all.sh --all-metaprograms    # gavish on all 41 metaprograms
#   ./cell_first_all.sh emt gavish            # only the named step(s), in file order
#
# Step names: scie, emt, gavish, composition_scie, composition_emt, composition_gavish.
# Logs go to 03_6_cell_first/logs/cell_first_all_<timestamp>.log and to the terminal (and to
# cell_first_all_<jobid>.out with --slurm).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PHASE_DIR="$(dirname "$SCRIPT_DIR")"
LOG_DIR="$SCRIPT_DIR/logs"
PYTHON="${PYTHON:-python3}"

STEP_NAMES=(scie emt gavish composition_scie composition_emt composition_gavish)

ARGS=("$@")
FORCE=0; DRY_RUN=0; USE_SLURM=0; ALL_METAPROGRAMS=0
STEP_FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --force|-f) FORCE=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    --slurm) USE_SLURM=1; shift ;;
    --all-metaprograms) ALL_METAPROGRAMS=1; shift ;;
    -h|--help) sed -n '3,/^$/p' "${BASH_SOURCE[0]}"; exit 0 ;;
    -*) echo "unknown option: $1" >&2; exit 1 ;;
    *) STEP_FILTER+=("$1"); shift ;;
  esac
done

# Checked after the options, so --help works without an environment.
: "${DATA_DIR:?set DATA_DIR to the directory holding the datasets (outside the repo)}"

in_list() { local x="$1"; shift; local e; for e in "$@"; do [ "$e" = "$x" ] && return 0; done; return 1; }

for requested in "${STEP_FILTER[@]:-}"; do
  [ -n "$requested" ] || continue
  in_list "$requested" "${STEP_NAMES[@]}" || {
    echo "unknown step '$requested'; valid steps: ${STEP_NAMES[*]}" >&2; exit 1; }
done

# Before resolve_run.py: the submitting shell need not have the Python env.
if [ "$USE_SLURM" -eq 1 ]; then
  . "$PHASE_DIR/utils/slurm.sh"
  submit_driver cell_first_all 8 32G "$LOG_DIR" "$DRY_RUN" "$SCRIPT_DIR/cell_first_all.sh" \
    ${ARGS[@]+"${ARGS[@]}"}
  exit 0
fi

# Run id and output root from Python, so the resume check looks where the scripts write.
read -r RUN_ID COMPARTMENT OUT_ROOT_DIR < <(cd "$PHASE_DIR" && "$PYTHON" utils/resolve_run.py)
TABLE_ROOT="$OUT_ROOT_DIR/tables"
if [ "$OUT_ROOT_DIR" != "$PHASE_DIR" ]; then LOG_DIR="$OUT_ROOT_DIR/logs"; fi

# Mirrors sig_collections.resolve(): gavish writes under gavish_tnbc unless --all-metaprograms.
GAVISH_SLUG="gavish_tnbc"
[ "$ALL_METAPROGRAMS" -eq 1 ] && GAVISH_SLUG="gavish"

STEP_DIR="$TABLE_ROOT/03_6_cell_first"
STEP_OUTPUTS=(
  "$STEP_DIR/dimension_row_order_scie_${RUN_ID}.csv"
  "$STEP_DIR/dimension_row_order_emt_${RUN_ID}.csv"
  "$STEP_DIR/dimension_row_order_${GAVISH_SLUG}_${RUN_ID}.csv"
  "$TABLE_ROOT/03_4_signatures/signature_concentration_scie_${RUN_ID}.csv"
  "$TABLE_ROOT/03_4_signatures/signature_concentration_emt_${RUN_ID}.csv"
  "$TABLE_ROOT/03_4_signatures/signature_concentration_${GAVISH_SLUG}_${RUN_ID}.csv"
)
CYTOTRACE_CSV="$DATA_DIR/03_tum/cytotrace2_${RUN_ID}.csv"

mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/cell_first_all_$(date +%Y%m%d_%H%M%S).log"
[ "$DRY_RUN" -eq 0 ] && exec > >(tee -a "$LOG_FILE") 2>&1

echo "03_6 cell-first"
echo "run id    : $RUN_ID  (compartment $COMPARTMENT)"
echo "DATA_DIR  : $DATA_DIR"
echo "gavish    : $GAVISH_SLUG"
if [ -e "$CYTOTRACE_CSV" ]; then
  echo "cytotrace : $CYTOTRACE_CSV"
else
  echo "cytotrace : absent, scie runs without CytoTRACE2 (run 03_5 first to include it)"
fi
[ "$DRY_RUN" -eq 0 ] && echo "log       : $LOG_FILE"
echo

n_run=0; n_have=0
for i in "${!STEP_NAMES[@]}"; do
  name="${STEP_NAMES[$i]}"
  output="${STEP_OUTPUTS[$i]}"

  if [ ${#STEP_FILTER[@]} -gt 0 ] && ! in_list "$name" "${STEP_FILTER[@]}"; then
    continue
  fi

  if [ "$FORCE" -eq 0 ] && [ -e "$output" ]; then
    echo "[have] $name  -> $output"
    n_have=$((n_have + 1)); continue
  fi

  # The composition steps are 03_4's script, run here because they read 03_6's scores.
  coll="${name#composition_}"
  step_dir="$SCRIPT_DIR"
  if [ "$coll" != "$name" ]; then
    cmd=("$PYTHON" signature_composition_tum.py --collection "$coll")
    step_dir="$PHASE_DIR/03_4_signatures"
  else
    cmd=("$PYTHON" cell_first_tum.py --collection "$name")
    if [ "$FORCE" -eq 1 ]; then cmd+=(--overwrite); fi
  fi
  # An `if`, not `cond && ...`, which would trip set -e.
  if [ "$ALL_METAPROGRAMS" -eq 1 ] && [ "$coll" = "gavish" ]; then cmd+=(--all-metaprograms); fi

  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry ] $name  -> ${cmd[*]}"
    n_run=$((n_run + 1)); continue
  fi

  echo "[run ] $name  -> ${cmd[*]}"
  ( cd "$step_dir" && "${cmd[@]}" )
  # Stop before the next step reads a stale file.
  [ -e "$output" ] || { echo "[error] $name finished but $output is missing" >&2; exit 1; }
  echo "[done] $name"
  echo
  n_run=$((n_run + 1))
done

echo "done: $n_run step(s) run, $n_have already done"
