#!/usr/bin/env bash
#
# 03_4 signatures: every collection's .gmt, coverage and Jaccard, then the EMT check, in order:
#
#   1. build_signatures_tum.py --collection scie    -> signatures_scie_tum.gmt + tables
#   2. build_signatures_tum.py --collection emt     -> signatures_emt_tum.gmt + tables
#   3. build_signatures_tum.py --collection gavish  -> signatures_gavish_tnbc_tum.gmt + tables
#   4. emt_vs_gavish_tum.py                         -> emt_vs_gavish_* tables (under emt)
#
# Light (the objects are read backed), so it runs locally; --slurm submits it as one job
# (utils/slurm.sh). Steps with an existing output are skipped, so a re-run resumes;
# --force re-runs everything.
#
# Usage:
#   export DATA_DIR=<data root>
#   ./signatures_all.sh                       # every step, resuming
#   ./signatures_all.sh --slurm               # same chain as one job on the cluster
#   ./signatures_all.sh --force               # re-run everything, overwriting
#   ./signatures_all.sh --dry-run             # print what would run, do nothing
#   ./signatures_all.sh --all-metaprograms    # gavish and emtgavish on all 41 metaprograms
#   ./signatures_all.sh --allow-low-coverage  # report lists under the mapping floor, do not stop
#   ./signatures_all.sh emt emtgavish         # only the named step(s), in file order
#
# Step names: scie, emt, gavish, emtgavish.
# Logs go to 03_4_signatures/logs/signatures_all_<timestamp>.log and to the terminal (and to
# signatures_all_<jobid>.out with --slurm).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PHASE_DIR="$(dirname "$SCRIPT_DIR")"
LOG_DIR="$SCRIPT_DIR/logs"
PYTHON="${PYTHON:-python3}"

STEP_NAMES=(scie emt gavish emtgavish)

ARGS=("$@")
FORCE=0; DRY_RUN=0; USE_SLURM=0; ALL_METAPROGRAMS=0; ALLOW_LOW=0
STEP_FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --force|-f) FORCE=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    --slurm) USE_SLURM=1; shift ;;
    --all-metaprograms) ALL_METAPROGRAMS=1; shift ;;
    --allow-low-coverage) ALLOW_LOW=1; shift ;;
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
  submit_driver signatures_all 2 16G "$LOG_DIR" "$DRY_RUN" "$SCRIPT_DIR/signatures_all.sh" \
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

STEP_DIR="$TABLE_ROOT/03_4_signatures"
STEP_OUTPUTS=(
  "$STEP_DIR/signatures_scie_${COMPARTMENT}.gmt"
  "$STEP_DIR/signatures_emt_${COMPARTMENT}.gmt"
  "$STEP_DIR/signatures_${GAVISH_SLUG}_${COMPARTMENT}.gmt"
  "$STEP_DIR/emt_vs_gavish_pairs_emt_${RUN_ID}.csv"
)

mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/signatures_all_$(date +%Y%m%d_%H%M%S).log"
[ "$DRY_RUN" -eq 0 ] && exec > >(tee -a "$LOG_FILE") 2>&1

echo "03_4 signatures"
echo "run id    : $RUN_ID  (compartment $COMPARTMENT)"
echo "DATA_DIR  : $DATA_DIR"
echo "gavish    : $GAVISH_SLUG"
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

  if [ "$name" = "emtgavish" ]; then
    cmd=("$PYTHON" emt_vs_gavish_tum.py)
  else
    cmd=("$PYTHON" build_signatures_tum.py --collection "$name")
    if [ "$ALLOW_LOW" -eq 1 ]; then cmd+=(--allow-low-coverage); fi
  fi
  # An `if`, not `cond && ...`, which would trip set -e.
  if [ "$ALL_METAPROGRAMS" -eq 1 ] && [ "$name" != "scie" ] && [ "$name" != "emt" ]; then
    cmd+=(--all-metaprograms)
  fi

  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry ] $name  -> ${cmd[*]}"
    n_run=$((n_run + 1)); continue
  fi

  echo "[run ] $name  -> ${cmd[*]}"
  ( cd "$SCRIPT_DIR" && "${cmd[@]}" )
  # Stop before the next step reads a stale file.
  [ -e "$output" ] || { echo "[error] $name finished but $output is missing" >&2; exit 1; }
  echo "[done] $name"
  echo
  n_run=$((n_run + 1))
done

echo "done: $n_run step(s) run, $n_have already done"
