#!/usr/bin/env bash
#
# 02_4 metrics: score every integrated grid row with the 12 scib metrics, one
# task per (run, output type). The Unintegrated baseline is scored against itself.
#
#   default (local): check_integrations.py + metrics.py in sequence, then merge.
#   --slurm        : submit_metrics.slurm, 5 at a time, knn tasks with more memory;
#                    merge + summary are run by hand afterwards (command printed).
#
# CSVs go to $DATA_DIR/02_metrics/<task>/metrics/<scaling>/hvg/<method>_<type>.csv.
#
# Usage:
#   export DATA_DIR=<data root>
#   bash run_all_metrics.sh                      # local, whole grid
#   bash run_all_metrics.sh --slurm              # submit the whole grid (5 at a time)
#   bash run_all_metrics.sh --method harmony     # one method, both scalings
#   bash run_all_metrics.sh --slurm --method harmony
#   bash run_all_metrics.sh --scaling unscaled   # the unscaled half
#   bash run_all_metrics.sh bbknn_unscaled ...   # named run_id(s)
#   bash run_all_metrics.sh --dry-run            # print what would run, do nothing
#   bash run_all_metrics.sh --force              # re-score tasks whose CSV exists
#
# Tasks whose CSV exists are skipped unless --force. Filters combine.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GRID="$SCRIPT_DIR/../benchmark_grid.tsv"

# DATA_DIR is required (cluster_env.sh).
. "$SCRIPT_DIR/../utils/cluster_env.sh"
resolve_data_dir || exit 1
TASK="shiao"
ENV="${METRICS_ENV:-benchmark-py-r}"    # same env locally and in the .slurm files
BATCH_KEY="cohort"; LABEL_KEY="cell_type"; ORGANISM="human"
ROOT="$DATA_DIR/02_metrics"
MERGED="$DATA_DIR/02_metrics_merged.csv"

# SLURM sizing: 4 CPUs per task, five at a time.
MAX_CONCURRENT=5

# Memory per output type: knn (BBKNN, 8.7e8 nonzeros) peaked at 197 GiB in leiden.
MEM_STD="${METRICS_MEM:-96G}"
MEM_KNN="${METRICS_MEM_KNN:-300G}"

# node02 kills jobs at launch (2026-07-30); clear METRICS_EXCLUDE once fixed.
EXCLUDE_NODES="${METRICS_EXCLUDE-node02}"

# Leiden flavor for NMI/ARI; must be the same for every method.
CLUSTER_FLAVOR="${CLUSTER_FLAVOR:-igraph}"

[ -f "$GRID" ] || { echo "grid not found: $GRID" >&2; exit 1; }

conda_guarded() { set +u; conda "$@"; set -u; }
in_list() { local x="$1"; shift; local e; for e in "$@"; do [ "$e" = "$x" ] && return 0; done; return 1; }

# Abort if an option has no value (shift 2 would otherwise fail silently).
need_value() {
  [ $# -ge 2 ] || { echo "$1 needs a value" >&2; exit 1; }
}

# Abort unless every given value appears in column `col` of the grid.
check_in_grid() {
  local what="$1" col="$2"; shift 2
  [ $# -gt 0 ] || return 0
  local valid; valid="$(awk -F'\t' -v c="$col" 'NR>1 && $1!="" { print $c }' "$GRID" | sort -u)"
  local v
  for v in "$@"; do
    if ! grep -qxF -- "$v" <<< "$valid"; then
      { echo "unknown $what '$v'; the grid has:"; sed 's/^/  /' <<< "$valid"; } >&2
      exit 1
    fi
  done
}

# Parse the command line.
USE_SLURM=0; NO_CHECK=0; DRY_RUN=0; FORCE=0
SCALING_FILTER=""
METHOD_FILTER=(); RUNID_FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --slurm) USE_SLURM=1; shift ;;
    --no-check) NO_CHECK=1; shift ;;
    --force|-f) FORCE=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    --scaling) need_value "$@"; SCALING_FILTER="$2"; shift 2 ;;
    --scaling=*) SCALING_FILTER="${1#*=}"; shift ;;
    --method) need_value "$@"; METHOD_FILTER+=("$2"); shift 2 ;;
    --method=*) METHOD_FILTER+=("${1#*=}"); shift ;;
    # Unknown option: a typo, not a run id.
    -*) echo "unknown option '$1'" >&2; exit 1 ;;
    *) RUNID_FILTER+=("$1"); shift ;;
  esac
done
if [ -n "$SCALING_FILTER" ] && [ "$SCALING_FILTER" != "scaled" ] && [ "$SCALING_FILTER" != "unscaled" ]; then
  echo "--scaling must be 'scaled' or 'unscaled', got '$SCALING_FILTER'" >&2; exit 1
fi

# Reject run ids / methods not in the grid instead of silently matching nothing.
check_in_grid "run id" 1 "${RUNID_FILTER[@]+"${RUNID_FILTER[@]}"}"
check_in_grid "method" 2 "${METHOD_FILTER[@]+"${METHOD_FILTER[@]}"}"

# A grid row (by run_id/scaling/method) is wanted if it matches every filter.
want_row() {
  local run_id="$1" scaling="$2" method="$3"
  [ -n "$SCALING_FILTER" ] && [ "$scaling" != "$SCALING_FILTER" ] && return 1
  [ ${#METHOD_FILTER[@]} -gt 0 ] && ! in_list "$method" "${METHOD_FILTER[@]}" && return 1
  [ ${#RUNID_FILTER[@]} -gt 0 ] && ! in_list "$run_id" "${RUNID_FILTER[@]}" && return 1
  return 0
}

# CSV path of a (run, type) pair, shared by local and SLURM mode.
csv_path() { echo "$ROOT/$TASK/metrics/$2/hvg/${1}_${3}.csv"; }

# metrics.py writes the CSV only at the end, so its existence means done.
already_done() {
  local csv="$1"
  [ "$FORCE" -eq 0 ] || return 1
  [ -f "$csv" ] || return 1
  return 0
}

# SLURM mode: resolve the matching array indices and submit the arrays.

# Submit one group of array indices (groups exist only to set --mem).
submit_group() {
  local label="$1" mem="$2"; shift 2
  [ $# -gt 0 ] || return 0

  local spec; spec="$(IFS=,; echo "$*")%${MAX_CONCURRENT}"
  # logs/ must exist before submission: SLURM opens --output before the script runs.
  local opts=(--chdir="$SCRIPT_DIR"
              --export="ALL,DATA_DIR=$DATA_DIR,METRICS_DIR=$SCRIPT_DIR,CLUSTER_FLAVOR=$CLUSTER_FLAVOR"
              --array="$spec" --mem="$mem")
  [ -n "$EXCLUDE_NODES" ] && opts+=(--exclude="$EXCLUDE_NODES")

  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry] $label: sbatch ${opts[*]} $SCRIPT_DIR/submit_metrics.slurm"
    return 0
  fi

  mkdir -p "$SCRIPT_DIR/logs"
  local m; m="$(sbatch --parsable "${opts[@]}" "$SCRIPT_DIR/submit_metrics.slurm")"
  echo "submitted $label metrics array: job $m   (--mem=$mem, tasks $spec)"
}

run_slurm() {
  # Expand the grid as submit_metrics.slurm does; scored tasks are left out of the
  # spec, and knn tasks go in their own group.
  local indices=() indices_knn=() total=0 n_have=0 n_skip=0
  while IFS=$'\t' read -r idx run_id method scaling type output reference; do
    total=$((total + 1))
    want_row "$run_id" "$scaling" "$method" || continue
    if already_done "$(csv_path "$method" "$scaling" "$type")"; then
      echo "[have] $run_id $type: already scored, not submitting"
      n_have=$((n_have + 1)); continue
    fi
    # Do not queue tasks whose inputs (integrated or reference) are not synced yet.
    local missing=""
    [ -f "$DATA_DIR/$output" ]    || missing="$DATA_DIR/$output"
    [ -f "$DATA_DIR/$reference" ] || missing="${missing:+$missing, }$DATA_DIR/$reference"
    if [ -n "$missing" ]; then
      echo "[skip] $run_id $type: input not here yet ($missing)"
      n_skip=$((n_skip + 1)); continue
    fi
    if [ "$type" = "knn" ]; then indices_knn+=("$idx"); else indices+=("$idx"); fi
  done < <(awk -F'\t' 'NR>1 && $1!="" { n=split($8, ts, ","); for (i=1;i<=n;i++){ idx++; gsub(/[ \t]/,"",ts[i]); print idx"\t"$1"\t"$2"\t"$5"\t"ts[i]"\t"$7"\t"$9 } }' "$GRID")

  local n_sel=$(( ${#indices[@]} + ${#indices_knn[@]} ))
  if [ "$n_sel" -eq 0 ]; then
    if [ "$n_have" -gt 0 ]; then
      echo "nothing to submit: all $n_have matching task(s) are already scored (--force to re-run)"
      return 0
    fi
    echo "no array task matches the given filters" >&2; exit 1
  fi

  echo "selected ${n_sel}/${total} array task(s), ${n_have} already scored, ${n_skip} not integrated yet, ${MAX_CONCURRENT} at a time"
  echo "  clustering: leiden ($CLUSTER_FLAVOR)${EXCLUDE_NODES:+   excluding: $EXCLUDE_NODES}"
  submit_group "full/embed" "$MEM_STD" ${indices[@]+"${indices[@]}"}
  submit_group "knn       " "$MEM_KNN" ${indices_knn[@]+"${indices_knn[@]}"}

  echo "when the arrays finish, merge + plot locally:"
  echo "  python $SCRIPT_DIR/merge_metrics.py -o $MERGED -r $ROOT --glob '$ROOT/$TASK/metrics/*/hvg/*.csv'"
  echo "  Rscript $SCRIPT_DIR/make_summary_table.R -i $MERGED"
}

# Local mode: run every matching (run, type) here, then merge.

run_local() {
  if [ "$DRY_RUN" -eq 0 ]; then
    set +u; eval "$(conda shell.bash hook)"; set -u
    conda_guarded activate "$ENV"
  fi

  local n_done=0 n_skip=0 n_have=0
  {
    read -r _header
    while IFS=$'\t' read -r run_id method language env scaling input output types reference hvgs; do
      [ -n "$run_id" ] || continue
      if ! want_row "$run_id" "$scaling" "$method"; then n_skip=$((n_skip + 1)); continue; fi

      local int_path="$DATA_DIR/$output" ref_path="$DATA_DIR/$reference"
      local out_dir="$ROOT/$TASK/metrics/$scaling/hvg"

      # Per (run, type), so a half-scored row resumes on the missing half.
      IFS=',' read -ra tps <<< "$types"
      local todo=() t_clean
      for t in "${tps[@]}"; do
        t_clean="$(echo "$t" | tr -d '[:space:]')"
        if already_done "$(csv_path "$method" "$scaling" "$t_clean")"; then
          echo "[have] $run_id $t_clean: already scored ($out_dir/${method}_${t_clean}.csv)"
          n_have=$((n_have + 1))
        else
          todo+=("$t_clean")
        fi
      done
      [ ${#todo[@]} -eq 0 ] && continue

      # Check done before inputs, so scored tasks are not reported as missing.
      if [ ! -f "$int_path" ]; then
        echo "[skip] $run_id: not integrated yet ($int_path)"; n_skip=$((n_skip + 1)); continue
      fi

      if [ "$DRY_RUN" -eq 1 ]; then
        echo "[dry] $run_id  types=${todo[*]}  ref=$ref_path"
        for t in "${todo[@]}"; do echo "      $method $t -> $out_dir/${method}_${t}.csv"; done
        n_done=$((n_done + 1)); continue
      fi

      mkdir -p "$out_dir"
      echo "==================================================================="
      echo ">>> $run_id  ($method, types=$types)"
      echo "==================================================================="
      if [ "$NO_CHECK" -eq 0 ]; then
        python "$SCRIPT_DIR/check_integrations.py" -i "$int_path" -u "$ref_path" \
          --types "$types" -b "$BATCH_KEY" -l "$LABEL_KEY"
      fi
      for t in "${todo[@]}"; do
        local csv="$out_dir/${method}_${t}.csv"
        echo "--- metrics: $method $t ---"
        python "$SCRIPT_DIR/metrics.py" -u "$ref_path" -i "$int_path" -o "$csv" \
          -m "$method" --type "$t" --scaling "$scaling" -b "$BATCH_KEY" -l "$LABEL_KEY" \
          --organism "$ORGANISM" --hvgs "$hvgs" --cluster-flavor "$CLUSTER_FLAVOR" -v
      done
      n_done=$((n_done + 1)); echo "[ok] $run_id"
    done
  } < "$GRID"

  if [ "$DRY_RUN" -eq 0 ]; then
    echo ">>> merging per-run CSVs"
    python "$SCRIPT_DIR/merge_metrics.py" -o "$MERGED" -r "$ROOT" \
      --glob "$ROOT/$TASK/metrics/*/hvg/*.csv" || echo "[warn] merge produced nothing"
    conda_guarded deactivate
  fi
  echo "==================================================================="
  echo "done: $n_done run(s) scored, $n_have task(s) already scored, $n_skip skipped   ->  $MERGED"
  [ "$n_have" -gt 0 ] && [ "$FORCE" -eq 0 ] && \
    echo "      (--force to re-score the $n_have existing task(s))"
  return 0
}


if [ "$USE_SLURM" -eq 1 ]; then run_slurm; else run_local; fi
