#!/usr/bin/env bash
#
# 03_8 convergence (Route C): every collection, then the steps that read all three, in order:
#
#   1. convergence_tum.py --collection scie    -> convergence_scie (+ target_axes, figures)
#   2. convergence_tum.py --collection emt     -> same, on the EMT lists
#   3. convergence_tum.py --collection gavish  -> same, on the TNBC metaprograms (no target)
#   4. route_a_null_tum.py                     -> route_a_null_* (the null of the rho bar)
#   5. leiden_smi_tum.py                       -> SMI of the live directions vs the Leiden clusters
#   6. dr_naming_tum.py                        -> dr_naming / program_dimensions (consensus), the
#                                                 STATE_PICKS clusters overlaid in orange
#   7. dr_umap_tum.py                          -> one UMAP per named bar
#   8. leiden_states_tum.py --states-de        -> the STATE_PICKS clusters: DE, composition, ORA
#   9. leiden_states_tum.py --states-dr        -> DR enrichment per cluster (reads 8)
#
# Steps 8-9 read STATE_PICKS (dr_naming_tum.py) and CLUSTER_READING (leiden_states_tum.py),
# curated by hand on one Leiden partition: re-curate them after a re-run of 03_3.
# Needs 03_6 and 03_7 for all three collections (03_7's ORA tsv in $DATA_DIR/03_tum included)
# and 03_3's Leiden partitions. The max(DR) groups are 03_8b, run after this.
# Steps with an existing output are skipped, so a re-run resumes; --force re-runs everything.
# --slurm submits it as one job (utils/slurm.sh).
#
# Usage:
#   export DATA_DIR=<data root>
#   ./convergence_all.sh                      # steps 1-9, resuming
#   ./convergence_all.sh --slurm              # same chain as one job on the cluster
#   ./convergence_all.sh --force              # re-run everything, overwriting
#   ./convergence_all.sh --dry-run            # print what would run, do nothing
#   ./convergence_all.sh naming umap          # only the named step(s), in file order
#
# Step names: scie, emt, gavish, null, smi, naming, umap, states_de, states_dr.
# Logs go to 03_8_convergence/logs/convergence_all_<timestamp>.log and to the terminal (and to
# convergence_all_<jobid>.out with --slurm).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PHASE_DIR="$(dirname "$SCRIPT_DIR")"
LOG_DIR="$SCRIPT_DIR/logs"
PYTHON="${PYTHON:-python3}"

STEP_NAMES=(scie emt gavish null smi naming umap states_de states_dr)

ARGS=("$@")
FORCE=0; DRY_RUN=0; USE_SLURM=0
STEP_FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --force|-f) FORCE=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    --slurm) USE_SLURM=1; shift ;;
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
  submit_driver convergence_all 8 32G "$LOG_DIR" "$DRY_RUN" "$SCRIPT_DIR/convergence_all.sh" \
    ${ARGS[@]+"${ARGS[@]}"}
  exit 0
fi

# Run id and output root from Python, so the resume check looks where the scripts write.
read -r RUN_ID COMPARTMENT OUT_ROOT_DIR < <(cd "$PHASE_DIR" && "$PYTHON" utils/resolve_run.py)
RES="$(cd "$PHASE_DIR/utils" && "$PYTHON" -c 'import cell_set as CS; print(f"{CS.LEIDEN_RESOLUTION:g}")')"
TABLE_ROOT="$OUT_ROOT_DIR/tables/03_8_convergence"
FIG_ROOT="$OUT_ROOT_DIR/figures/03_8_convergence"
if [ "$OUT_ROOT_DIR" != "$PHASE_DIR" ]; then LOG_DIR="$OUT_ROOT_DIR/logs"; fi

# The output that marks each step done, in STEP_NAMES order (a directory: done if not empty).
STEP_OUTPUTS=(
  "$TABLE_ROOT/convergence_scie_${RUN_ID}.csv"
  "$TABLE_ROOT/convergence_emt_${RUN_ID}.csv"
  "$TABLE_ROOT/convergence_gavish_tnbc_${RUN_ID}.csv"
  "$TABLE_ROOT/route_a_null_gavish_tnbc_${RUN_ID}.csv"
  "$TABLE_ROOT/smi_matches_drvi_leiden_${RES}_${RUN_ID}.csv"
  "$TABLE_ROOT/program_dimensions_consensus_${RUN_ID}.csv"
  "$FIG_ROOT/dr_umap"
  "$TABLE_ROOT/state_composition_drvi_leiden_${RES}_${RUN_ID}.csv"
  "$TABLE_ROOT/state_cluster_summary_drvi_leiden_${RES}_${RUN_ID}.csv"
)

is_done() { if [ -d "$1" ]; then [ -n "$(ls -A "$1")" ]; else [ -e "$1" ]; fi; }

mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/convergence_all_$(date +%Y%m%d_%H%M%S).log"
[ "$DRY_RUN" -eq 0 ] && exec > >(tee -a "$LOG_FILE") 2>&1

echo "03_8 convergence"
echo "run id    : $RUN_ID  (compartment $COMPARTMENT, Leiden $RES)"
echo "DATA_DIR  : $DATA_DIR"
[ "$DRY_RUN" -eq 0 ] && echo "log       : $LOG_FILE"
echo

n_run=0; n_have=0
for i in "${!STEP_NAMES[@]}"; do
  name="${STEP_NAMES[$i]}"
  output="${STEP_OUTPUTS[$i]}"

  if [ ${#STEP_FILTER[@]} -gt 0 ]; then
    in_list "$name" "${STEP_FILTER[@]}" || continue
  fi

  if [ "$FORCE" -eq 0 ] && is_done "$output"; then
    echo "[have] $name  -> $output"
    n_have=$((n_have + 1)); continue
  fi

  case "$name" in
    scie|emt|gavish) cmd=("$PYTHON" convergence_tum.py --collection "$name") ;;
    null)         cmd=("$PYTHON" route_a_null_tum.py) ;;
    smi)          cmd=("$PYTHON" leiden_smi_tum.py) ;;
    naming)       cmd=("$PYTHON" dr_naming_tum.py) ;;
    umap)         cmd=("$PYTHON" dr_umap_tum.py) ;;
    states_de)    cmd=("$PYTHON" leiden_states_tum.py --states-de) ;;
    states_dr)    cmd=("$PYTHON" leiden_states_tum.py --states-dr) ;;
  esac

  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry ] $name  -> ${cmd[*]}"
    n_run=$((n_run + 1)); continue
  fi

  echo "[run ] $name  -> ${cmd[*]}"
  ( cd "$SCRIPT_DIR" && "${cmd[@]}" )
  # Stop before the next step reads a stale file.
  is_done "$output" || { echo "[error] $name finished but $output is missing" >&2; exit 1; }
  echo "[done] $name"
  echo
  n_run=$((n_run + 1))
done

echo "done: $n_run step(s) run, $n_have already done"
