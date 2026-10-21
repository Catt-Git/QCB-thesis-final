#!/usr/bin/env bash
#
# 02_3 method plots: run plot_methods_umaps.py for every integrated grid row,
# into figures/<run_id>/. Locally (in sequence, failures listed at the end) or with
# --slurm (submit_plots.slurm, knn rows submitted with more memory).
#
# Plot type: embed if the run has one, else full, else knn. Each run is compared
# with the grid's `reference` (scaled runs with the scaled reference, whose UMAP
# comes from scale_batch.py).
#
# Usage:
#   export DATA_DIR=<data root>
#   bash plot_all.sh                            # every integrated run
#   bash plot_all.sh harmony_unscaled ...       # only the named run_id(s)
#   bash plot_all.sh --scaling unscaled         # only the unscaled half
#   bash plot_all.sh --dry-run                  # print what would run, do nothing
#   bash plot_all.sh --force                    # redraw runs whose panels already exist
#   bash plot_all.sh --force --recompute-umap   # redraw *and* lay the UMAPs out again
#   bash plot_all.sh --slurm --force            # the whole grid as an array job
#
# The integrated UMAP is cached in each .h5ad (--recompute-umap discards it). Runs
# with all five panels are skipped unless --force. On the cluster, figures land in
# the repo and must be rsynced back.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PHASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
GRID="$PHASE_DIR/benchmark_grid.tsv"
FIG_DIR="$PHASE_DIR/figures"          # the phase's one figures/, not a local one
PLOTTER="$SCRIPT_DIR/plot_methods_umaps.py"

# DATA_DIR is required (cluster_env.sh).
. "$SCRIPT_DIR/../utils/cluster_env.sh"
resolve_data_dir || exit 1

ENV="${PLOT_ENV:-benchmark-py-r}"       # same env locally and in the .slurm file

# SLURM memory: BBKNN's graph (8.7e8 nonzeros) needs far more than the rest.
MAX_CONCURRENT="${PLOT_MAX_CONCURRENT:-5}"
MEM_STD="${PLOT_MEM:-96G}"
MEM_KNN="${PLOT_MEM_KNN:-300G}"

# node02 kills jobs at launch (2026-07-30); clear PLOT_EXCLUDE once fixed.
EXCLUDE_NODES="${PLOT_EXCLUDE-node02}"

[ -f "$GRID" ] || { echo "grid not found: $GRID" >&2; exit 1; }

# conda activation scripts read unbound variables, so disable set -u around them.
conda_guarded() { set +u; conda "$@"; set -u; }

# Abort if an option has no value (shift 2 would otherwise fail silently).
need_value() {
  [ $# -ge 2 ] || { echo "$1 needs a value" >&2; exit 1; }
}

# Parse the command line: --scaling, --dry-run, --force, plus optional run_id filters.
SCALING_FILTER=""
DRY_RUN=0
FORCE=0
RECOMPUTE_UMAP=0
USE_SLURM=0
FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --scaling) need_value "$@"; SCALING_FILTER="$2"; shift 2 ;;
    --scaling=*) SCALING_FILTER="${1#*=}"; shift ;;
    --force|-f) FORCE=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    --recompute-umap) RECOMPUTE_UMAP=1; shift ;;
    --slurm) USE_SLURM=1; shift ;;
    # Unknown option: a typo, not a run id.
    -*) echo "unknown option '$1'" >&2; exit 1 ;;
    *) FILTER+=("$1"); shift ;;
  esac
done
if [ -n "$SCALING_FILTER" ] && [ "$SCALING_FILTER" != "scaled" ] && [ "$SCALING_FILTER" != "unscaled" ]; then
  echo "--scaling must be 'scaled' or 'unscaled', got '$SCALING_FILTER'" >&2
  exit 1
fi

# Rows matching all filters; baseline rows never (output == reference; see 01_6).
want() {
  local run_id="$1" scaling="$2" language="$3" f
  [ "$language" = "baseline" ] && return 1
  if [ -n "$SCALING_FILTER" ] && [ "$scaling" != "$SCALING_FILTER" ]; then
    return 1
  fi
  [ ${#FILTER[@]} -eq 0 ] && return 0
  for f in "${FILTER[@]}"; do [ "$f" = "$run_id" ] && return 0; done
  return 1
}

# The five panels plot_methods_umaps.py writes, in the order it writes them.
PANELS=(cohort_integrated cohort_int_vs_unint celltype_integrated
        celltype_int_vs_unint cohort_celltype_integrated)

# Done only if all five panels exist (a crash can leave a partial set).
already_done() {
  local dir="$1" m="$2" p
  [ "$FORCE" -eq 0 ] || return 1
  for p in "${PANELS[@]}"; do [ -f "$dir/${m}_${p}.png" ] || return 1; done
  return 0
}

# Best QC view for a run's output type(s): prefer the corrected embedding.
plot_type() {
  case ",$1," in
    *,embed,*) echo embed ;;
    *,full,*)  echo full ;;
    *,knn,*)   echo knn ;;
    *)         echo "$1" ;;
  esac
}

# Check both inputs: the reference is synced separately and may be missing.
missing_inputs() {
  local int_path="$1" ref_path="$2" missing=""
  [ -f "$int_path" ] || missing="$int_path"
  [ -f "$ref_path" ] || missing="${missing:+$missing, }$ref_path"
  echo "$missing"
}

# SLURM mode: submit one array task per matching run.

# Submit one group of array indices (groups exist only to set --mem).
submit_group() {
  local label="$1" mem="$2"; shift 2
  [ $# -gt 0 ] || return 0

  local spec; spec="$(IFS=,; echo "$*")%${MAX_CONCURRENT}"
  local export_vars="ALL,DATA_DIR=$DATA_DIR,PLOT_DIR=$SCRIPT_DIR"
  [ "$RECOMPUTE_UMAP" -eq 1 ] && export_vars="$export_vars,PLOT_RECOMPUTE_UMAP=1"
  # logs/ must exist before submission: SLURM opens --output before the script runs.
  local opts=(--chdir="$SCRIPT_DIR" --export="$export_vars" --array="$spec" --mem="$mem")
  [ -n "$EXCLUDE_NODES" ] && opts+=(--exclude="$EXCLUDE_NODES")

  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry] $label: sbatch ${opts[*]} $SCRIPT_DIR/submit_plots.slurm"
    return 0
  fi
  mkdir -p "$SCRIPT_DIR/logs"
  local j; j="$(sbatch --parsable "${opts[@]}" "$SCRIPT_DIR/submit_plots.slurm")"
  echo "submitted $label array: job $j   (--mem=$mem, tasks $spec)"
}

run_slurm() {
  # Array index = grid data row; done or input-less rows are left out of the spec.
  local indices=() indices_knn=() total=0 n_have=0 n_skip=0 idx=0
  local run_id method language env scaling input output types reference hvgs
  {
    read -r _header
    while IFS=$'\t' read -r run_id method language env scaling input output types reference hvgs; do
      [ -n "$run_id" ] || continue
      idx=$((idx + 1)); total=$((total + 1))
      want "$run_id" "$scaling" "$language" || continue

      if already_done "$FIG_DIR/$run_id" "$method"; then
        echo "[have] $run_id: already plotted, not submitting"
        n_have=$((n_have + 1)); continue
      fi
      local missing; missing="$(missing_inputs "$DATA_DIR/$output" "$DATA_DIR/$reference")"
      if [ -n "$missing" ]; then
        echo "[skip] $run_id: input not here yet ($missing)"
        n_skip=$((n_skip + 1)); continue
      fi
      if [ "$(plot_type "$types")" = "knn" ]; then indices_knn+=("$idx"); else indices+=("$idx"); fi
    done
  } < "$GRID"

  local n_sel=$(( ${#indices[@]} + ${#indices_knn[@]} ))
  if [ "$n_sel" -eq 0 ]; then
    if [ "$n_have" -gt 0 ]; then
      echo "nothing to submit: all $n_have matching run(s) are already plotted (--force to redraw)"
      return 0
    fi
    echo "no grid row matches the given filters" >&2; exit 1
  fi

  echo "selected ${n_sel}/${total} run(s), ${n_have} already plotted, ${n_skip} not integrated yet, ${MAX_CONCURRENT} at a time"
  echo "  figures -> $FIG_DIR/<run_id>${EXCLUDE_NODES:+   excluding: $EXCLUDE_NODES}"
  submit_group "std" "$MEM_STD" ${indices[@]+"${indices[@]}"}
  submit_group "knn" "$MEM_KNN" ${indices_knn[@]+"${indices_knn[@]}"}
  echo "the panels land in the repo, not in DATA_DIR; bring them back with:"
  echo "  rsync -av <cluster>:$FIG_DIR/ ./figures/"
}

# Local mode: plot every matching run here, in sequence.

run_local() {
  if [ "$DRY_RUN" -eq 0 ]; then
    set +u; eval "$(conda shell.bash hook)"; set -u
    conda_guarded activate "$ENV"
  fi

  local n_done=0 n_skip=0 n_have=0 n_match=0
  local failed=()
  local run_id method language env scaling input output types reference hvgs
  {
    read -r _header
    while IFS=$'\t' read -r run_id method language env scaling input output types reference hvgs; do
      [ -n "$run_id" ] || continue
      if ! want "$run_id" "$scaling" "$language"; then continue; fi
      n_match=$((n_match + 1))

      local int_path="$DATA_DIR/$output" ref_path="$DATA_DIR/$reference"
      local out_dir="$FIG_DIR/$run_id" ptype; ptype="$(plot_type "$types")"

      # Check done before inputs, so finished runs are not reported as missing.
      if already_done "$out_dir" "$method"; then
        echo "[have] $run_id: already plotted, skipping ($out_dir)"
        n_have=$((n_have + 1)); continue
      fi

      local missing; missing="$(missing_inputs "$int_path" "$ref_path")"
      if [ -n "$missing" ]; then
        echo "[skip] $run_id: input missing ($missing)"
        n_skip=$((n_skip + 1)); continue
      fi

      local umap_args=()
      [ "$RECOMPUTE_UMAP" -eq 1 ] && umap_args+=(--recompute-umap)

      if [ "$DRY_RUN" -eq 1 ]; then
        echo "[dry] $run_id  (--type $ptype)  -> $out_dir"
        echo "      python $PLOTTER -m $method --type $ptype -i $int_path -u $ref_path -o $out_dir ${umap_args[*]}"
        n_done=$((n_done + 1)); continue
      fi

      echo "==================================================================="
      echo ">>> $run_id  ($method, --type $ptype)"
      echo "==================================================================="
      # A failure does not stop the batch; </dev/null keeps the grid on stdin.
      if python "$PLOTTER" -m "$method" --type "$ptype" \
           -i "$int_path" -u "$ref_path" -o "$out_dir" "${umap_args[@]}" </dev/null; then
        n_done=$((n_done + 1)); echo "[ok] $run_id -> $out_dir"
      else
        failed+=("$run_id"); echo "[fail] $run_id: see the traceback above, continuing" >&2
      fi
    done
  } < "$GRID"

  [ "$DRY_RUN" -eq 0 ] && conda_guarded deactivate

  if [ "$n_match" -eq 0 ]; then
    echo "no grid row matches the given filters" >&2; exit 1
  fi

  echo "==================================================================="
  echo "done: $n_done run(s) plotted, $n_have already plotted, $n_skip skipped, ${#failed[@]} failed"
  if [ "$n_have" -gt 0 ] && [ "$FORCE" -eq 0 ]; then
    echo "      (--force to redraw the $n_have existing figure set(s))"
  fi
  if [ ${#failed[@]} -gt 0 ]; then
    echo "      failed: ${failed[*]}"
    return 1
  fi
  return 0
}


if [ "$USE_SLURM" -eq 1 ]; then run_slurm; else run_local; fi
