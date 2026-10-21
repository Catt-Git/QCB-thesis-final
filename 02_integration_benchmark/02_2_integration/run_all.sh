#!/usr/bin/env bash
#
# 02_2 integration: run every benchmark_grid.tsv row, locally (in sequence) or
# with --slurm (submit_integration.slurm, max 3 concurrent tasks).
#
# Per row, by language:
#   python (benchmark-py-r)  run_integration.py                    -> <run_id>.h5ad
#   python (scgen-py)        run_scgen.py                          -> <run_id>.h5ad
#   R      (benchmark-py-r)  run_integration.R + rds_to_h5ad.py    -> <run_id>.h5ad
#                            (.rds deleted afterwards unless --keep-rds)
#   notebook / baseline      skipped (DRVI: run_drvi.py; Unintegrated: nothing to do)
#
# Embeddings also go to 02_embeddings/<run_id>.npy. scVI/scANVI/scGen checkpoint
# to 02_<method>/ and reload on re-run; delete the .pt to retrain.
#
# Usage:
#   export DATA_DIR=<data root>
#   bash run_all.sh                            # local, every row
#   bash run_all.sh --slurm                    # submit every row (3 at a time)
#   bash run_all.sh harmony_unscaled ...       # only the named run_id(s)
#   bash run_all.sh --method scanorama         # one method, both scalings
#   bash run_all.sh --scaling unscaled         # every unscaled run (half the grid)
#   bash run_all.sh --scaling scaled --slurm   # every scaled run, on the cluster
#   bash run_all.sh --dry-run                  # print what would run, do nothing
#   bash run_all.sh --force                    # re-run and overwrite finished runs
#   bash run_all.sh --keep-rds                 # keep the R intermediates
#
# Finished rows (output exists) are skipped, so a re-run resumes; delete outputs
# truncated by a crash first. Filters combine.

set -euo pipefail

: "${DATA_DIR:?set DATA_DIR to the directory holding the prepared inputs}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GRID="$SCRIPT_DIR/../benchmark_grid.tsv"
HVG_RDS="$DATA_DIR/shiao_hvg_2k_unintegrated_list.rds"
MAX_CONCURRENT=3                        # SLURM: at most 3 nodes in `normal`

# Reference batches for Seurat CCA/RPCA (not the grid's `reference` column).
# Empty = run_integration.R default; e.g. SEURAT_REFERENCE=Patient53,Patient16
SEURAT_REF_ARG=()
if [ -n "${SEURAT_REFERENCE:-}" ]; then SEURAT_REF_ARG=(--reference "$SEURAT_REFERENCE"); fi

[ -f "$GRID" ] || { echo "grid not found: $GRID" >&2; exit 1; }

# conda activation scripts read unbound variables, so disable set -u around them.
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
USE_SLURM=0; DRY_RUN=0; KEEP_RDS=0; FORCE=0
SCALING_FILTER=""
METHOD_FILTER=(); RUNID_FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --slurm) USE_SLURM=1; shift ;;
    --keep-rds) KEEP_RDS=1; shift ;;
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

# A row is wanted if it matches every active filter (scaling, method, run_id).
want_row() {
  local run_id="$1" scaling="$2" method="$3"
  [ -n "$SCALING_FILTER" ] && [ "$scaling" != "$SCALING_FILTER" ] && return 1
  [ ${#METHOD_FILTER[@]} -gt 0 ] && ! in_list "$method" "${METHOD_FILTER[@]}" && return 1
  [ ${#RUNID_FILTER[@]} -gt 0 ] && ! in_list "$run_id" "${RUNID_FILTER[@]}" && return 1
  return 0
}

# Existence only: a truncated .h5ad counts as done (check_integrations.py finds it).
already_done() { [ "$FORCE" -eq 0 ] && [ -f "$1" ]; }

# Whether the run has an embedding worth exporting (embed appears in types).
has_embed() { case ",$1," in *,embed,*) return 0 ;; *) return 1 ;; esac }

# Sets model_arg: checkpoint dir for trained methods (see model_paths.py), else empty.
set_model_arg() {
  model_arg=()
  case "$1" in
    scvi|scanvi|scgen) model_arg=(--model-dir "$DATA_DIR/02_${1}") ;;
  esac
}

# Remove an integrated .rds once it has been converted, unless --keep-rds.
drop_rds() {
  local p="$1"
  if [ "$KEEP_RDS" -eq 1 ]; then echo "[keep] $p"; return 0; fi
  [ -f "$p" ] || return 0
  echo "[clean] removing intermediate $(du -h "$p" | cut -f1)  $p"
  rm -f "$p"
}

# SLURM mode: submit submit_integration.slurm for the matching grid rows.

run_slurm() {
  # One array task per grid row; notebook, baseline and finished rows are left
  # out (an instantly-exiting task still costs a throttle slot).
  local indices=() total=0 n_have=0
  while IFS=$'\t' read -r idx run_id method scaling language output; do
    total=$((total + 1))
    { [ "$language" = "notebook" ] || [ "$language" = "baseline" ]; } && continue
    want_row "$run_id" "$scaling" "$method" || continue
    if already_done "$DATA_DIR/$output"; then
      echo "[have] $run_id: already integrated, not submitting ($DATA_DIR/$output)"
      n_have=$((n_have + 1)); continue
    fi
    indices+=("$idx")
  done < <(awk -F'\t' 'NR>1 && $1!="" { idx++; print idx"\t"$1"\t"$2"\t"$5"\t"$3"\t"$7 }' "$GRID")

  if [ ${#indices[@]} -eq 0 ]; then
    if [ "$n_have" -gt 0 ]; then
      echo "nothing to submit: all $n_have matching row(s) are already integrated (--force to re-run)"
      return 0
    fi
    echo "no grid row matches the given filters (or all are notebooks)" >&2; exit 1
  fi

  local spec; spec="$(IFS=,; echo "${indices[*]}")%${MAX_CONCURRENT}"
  echo "selected ${#indices[@]}/${total} row(s), ${n_have} already integrated, throttled to ${MAX_CONCURRENT}: $spec"

  # SLURM runs a spool copy of the script, so pass the repo path.
  local exports="ALL,DATA_DIR=$DATA_DIR,KEEP_RDS=$KEEP_RDS,INTEGRATION_DIR=$SCRIPT_DIR"
  # logs/ must exist before submission: SLURM opens --output before the script runs.
  local opts=(--chdir="$SCRIPT_DIR" --export="$exports" --array="$spec")
  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry] sbatch ${opts[*]} $SCRIPT_DIR/submit_integration.slurm"
    return 0
  fi
  mkdir -p "$SCRIPT_DIR/logs"
  local j; j="$(sbatch --parsable "${opts[@]}" "$SCRIPT_DIR/submit_integration.slurm")"
  echo "submitted integration array: job $j"
  echo "DRVI is not in the array: pick the latent size in shiao_drvi_128.ipynb, then run that size"
  echo "here (python run_drvi.py --n-latent N) or on the cluster (sbatch submit_drvi.slurm --n-latent N)."
  echo "Then score with 02_4_metrics/run_all_metrics.sh."
}

# Local mode: run every matching integration here, in sequence.

run_local() {
  # conda activate inside a non-interactive shell needs the hook sourced first.
  [ "$DRY_RUN" -eq 0 ] && { set +u; eval "$(conda shell.bash hook)"; set -u; }

  local n_done=0 n_skip=0 n_have=0
  {
    read -r _header
    while IFS=$'\t' read -r run_id method language env scaling input output types ref_h5ad hvgs; do
      [ -n "$run_id" ] || continue
      if ! want_row "$run_id" "$scaling" "$method"; then n_skip=$((n_skip + 1)); continue; fi

      # Before already_done(): a baseline output is the prepared input itself,
      # and --force must never overwrite it.
      if [ "$language" = "baseline" ]; then
        echo "[skip] $run_id is the Unintegrated baseline; nothing to integrate (02_4 scores it)"
        n_skip=$((n_skip + 1)); continue
      fi
      if [ "$language" = "notebook" ]; then
        echo "[skip] $run_id is a notebook (DRVI); see run_drvi.py"; n_skip=$((n_skip + 1)); continue
      fi

      local in_path="$DATA_DIR/$input" out_path="$DATA_DIR/$output"
      local emb_path="$DATA_DIR/02_embeddings/${run_id}.npy"

      if already_done "$out_path"; then
        echo "[have] $run_id: already integrated, skipping ($out_path)"
        n_have=$((n_have + 1)); continue
      fi

      if [ "$DRY_RUN" -eq 1 ]; then
        echo "[dry] $run_id  ($method, $language, $env)  -> $out_path"
        n_done=$((n_done + 1)); continue
      fi
      mkdir -p "$(dirname "$out_path")"

      echo "==================================================================="
      echo ">>> $run_id  ($method, $language, $env)"
      echo "==================================================================="

      case "$language" in
        python)
          conda_guarded activate "$env"
          local model_arg; set_model_arg "$method"
          if [ "$method" = "scgen" ]; then
            python "$SCRIPT_DIR/run_scgen.py" -i "$in_path" -o "$out_path" "${model_arg[@]}"
          else
            local emb_arg=(); has_embed "$types" && emb_arg=(--emb-out "$emb_path")
            python "$SCRIPT_DIR/run_integration.py" -m "$method" -i "$in_path" -o "$out_path" \
              "${emb_arg[@]}" "${model_arg[@]}"
          fi
          conda_guarded deactivate
          ;;
        R)
          conda_guarded activate "$env"
          local rds_out="$DATA_DIR/02_integration/${run_id}.rds"
          Rscript "$SCRIPT_DIR/run_integration.R" -m "$method" -i "$in_path" -o "$rds_out" \
            -v "$HVG_RDS" "${SEURAT_REF_ARG[@]}"
          local emb_arg=(); has_embed "$types" && emb_arg=(--emb-out "$emb_path")
          python "$SCRIPT_DIR/rds_to_h5ad.py" -i "$rds_out" -o "$out_path" --types "$types" "${emb_arg[@]}"
          # Unused downstream and huge (15 GB for fastMNN); safe once the .h5ad passed checks.
          drop_rds "$rds_out"
          conda_guarded deactivate
          ;;
        *)
          echo "[error] unknown language '$language' for $run_id" >&2; exit 1 ;;
      esac
      n_done=$((n_done + 1)); echo "[ok] $run_id -> $out_path"
    done
  } < "$GRID"

  echo "==================================================================="
  if [ "$DRY_RUN" -eq 1 ]; then
    echo "dry run: $n_done run(s) would run, $n_have already integrated, $n_skip filtered out"
    return 0
  fi
  echo "done: $n_done run(s) executed, $n_have already integrated, $n_skip filtered out"
  [ "$n_have" -gt 0 ] && [ "$FORCE" -eq 0 ] && \
    echo "      (--force to re-run and overwrite the $n_have existing output(s))"
  return 0
}


if [ "$USE_SLURM" -eq 1 ]; then run_slurm; else run_local; fi
