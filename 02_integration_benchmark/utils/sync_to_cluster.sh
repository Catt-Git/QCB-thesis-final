#!/usr/bin/env bash
#
# 02 utils: upload the grid's files (integrated locally) to the same relative
# paths under the cluster's DATA_DIR (cluster_env.sh), for the cluster-side metrics.
#
#   --what metrics (default) : integrated objects + their references (for 02_4)
#   --what inputs            : 02_1 prepared inputs + HVG lists (for 02_2)
#   --what all               : both
#
# Usage:
#   export DATA_DIR=<data root>
#   export CLUSTER_HOST=albertoc@hpc.example.it                # user@host, required
#   bash sync_to_cluster.sh --dry-run                          # list what would be sent
#   bash sync_to_cluster.sh                                    # the whole grid
#   bash sync_to_cluster.sh --scaling unscaled                 # the unscaled half
#   bash sync_to_cluster.sh --method harmony                   # one method
#   bash sync_to_cluster.sh bbknn_unscaled drvi_unscaled_128   # named run_id(s)
#   bash sync_to_cluster.sh --what inputs                      # the 02_1 inputs instead
#   bash sync_to_cluster.sh --scp                              # scp instead of rsync
#
# Files already on the cluster with the same size are skipped unless --force.
# One shared SSH connection, so you authenticate once. Filters combine.

set -euo pipefail

: "${DATA_DIR:?set DATA_DIR to the local directory holding the objects to upload}"
: "${CLUSTER_HOST:?set CLUSTER_HOST to user@hostname of the cluster}"

UTILS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GRID="$UTILS_DIR/../benchmark_grid.tsv"

# Remote root shared with the cluster-side wrappers.
. "$UTILS_DIR/cluster_env.sh"
REMOTE="$CLUSTER_DATA_DIR"

# Not in the grid (only the R methods read them), but part of `inputs`.
HVG_LISTS=(shiao_hvg_2k_unintegrated_list.csv shiao_hvg_2k_unintegrated_list.rds)

[ -f "$GRID" ] || { echo "grid not found: $GRID" >&2; exit 1; }

in_list() { local x="$1"; shift; local e; for e in "$@"; do [ "$e" = "$x" ] && return 0; done; return 1; }
human() { numfmt --to=iec --suffix=B "$1" 2>/dev/null || echo "$1 bytes"; }

# Parse the command line.
WHAT="metrics"; USE_SCP=0; DRY_RUN=0; FORCE=0
SCALING_FILTER=""
METHOD_FILTER=(); RUNID_FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --what) WHAT="${2:-}"; shift 2 ;;
    --what=*) WHAT="${1#*=}"; shift ;;
    --scp) USE_SCP=1; shift ;;
    --force|-f) FORCE=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    --scaling) SCALING_FILTER="${2:-}"; shift 2 ;;
    --scaling=*) SCALING_FILTER="${1#*=}"; shift ;;
    --method) METHOD_FILTER+=("${2:-}"); shift 2 ;;
    --method=*) METHOD_FILTER+=("${1#*=}"); shift ;;
    -h|--help) sed -n '2,22p' "${BASH_SOURCE[0]}"; exit 0 ;;
    -*) echo "unknown option: $1" >&2; exit 1 ;;
    *) RUNID_FILTER+=("$1"); shift ;;
  esac
done
case "$WHAT" in
  metrics|inputs|all) ;;
  *) echo "--what must be 'metrics', 'inputs' or 'all', got '$WHAT'" >&2; exit 1 ;;
esac
if [ -n "$SCALING_FILTER" ] && [ "$SCALING_FILTER" != "scaled" ] && [ "$SCALING_FILTER" != "unscaled" ]; then
  echo "--scaling must be 'scaled' or 'unscaled', got '$SCALING_FILTER'" >&2; exit 1
fi

# Reject typos up front, so a misspelled run id cannot silently select nothing.
mapfile -t ALL_RUNIDS < <(awk -F'\t' 'NR>1 && $1!="" { print $1 }' "$GRID")
mapfile -t ALL_METHODS < <(awk -F'\t' 'NR>1 && $1!="" { print $2 }' "$GRID" | sort -u)
for r in ${RUNID_FILTER+"${RUNID_FILTER[@]}"}; do
  in_list "$r" "${ALL_RUNIDS[@]}" || { echo "unknown run id: $r" >&2; exit 1; }
done
for m in ${METHOD_FILTER+"${METHOD_FILTER[@]}"}; do
  in_list "$m" "${ALL_METHODS[@]}" || { echo "unknown method: $m" >&2; exit 1; }
done

want_row() {
  local run_id="$1" scaling="$2" method="$3"
  [ -n "$SCALING_FILTER" ] && [ "$scaling" != "$SCALING_FILTER" ] && return 1
  [ ${#METHOD_FILTER[@]} -gt 0 ] && ! in_list "$method" "${METHOD_FILTER[@]}" && return 1
  [ ${#RUNID_FILTER[@]} -gt 0 ] && ! in_list "$run_id" "${RUNID_FILTER[@]}" && return 1
  return 0
}

# Paths to send, relative to DATA_DIR as in the grid.

WANTED=()
add_path() { in_list "$1" ${WANTED+"${WANTED[@]}"} || WANTED+=("$1"); }

{
  read -r _header
  while IFS=$'\t' read -r run_id method language env scaling input output types reference hvgs; do
    [ -n "$run_id" ] || continue
    want_row "$run_id" "$scaling" "$method" || continue
    case "$WHAT" in
      metrics) add_path "$output"; add_path "$reference" ;;
      inputs)  add_path "$input";  add_path "$reference" ;;
      all)     add_path "$output"; add_path "$input"; add_path "$reference" ;;
    esac
  done
} < "$GRID"

if [ "$WHAT" != "metrics" ]; then
  for f in "${HVG_LISTS[@]}"; do add_path "$f"; done
fi

[ ${#WANTED[@]} -gt 0 ] || { echo "no run matches the given filters" >&2; exit 1; }

# Skip files missing locally and, unless --force, those already on the cluster.

SSH_CTL="${TMPDIR:-/tmp}/sync_to_cluster_$$.sock"
SSH_OPTS=(-o ControlPath="$SSH_CTL")
cleanup() { ssh -O exit "${SSH_OPTS[@]}" "$CLUSTER_HOST" 2>/dev/null || true; rm -f "$SSH_CTL"; }
trap cleanup EXIT
echo ">>> opening the SSH connection to $CLUSTER_HOST (authenticate once)"
ssh -fN -o ControlMaster=yes -o ControlPersist=8h "${SSH_OPTS[@]}" "$CLUSTER_HOST"

declare -A REMOTE_SIZE=()
if [ "$FORCE" -eq 0 ]; then
  # One round trip: "<size> <path>" for each remote file that exists.
  while read -r size rel; do
    [ -n "${rel:-}" ] && REMOTE_SIZE["$rel"]="$size"
  done < <(printf '%s\n' "${WANTED[@]}" \
    | ssh "${SSH_OPTS[@]}" "$CLUSTER_HOST" \
        "cd '$REMOTE' 2>/dev/null && xargs -d '\n' -r stat -c '%s %n' 2>/dev/null" || true)
fi

SEND=(); total=0; n_have=0; n_missing=0
for rel in "${WANTED[@]}"; do
  local_path="$DATA_DIR/$rel"
  if [ ! -f "$local_path" ]; then
    echo "[skip] $rel: not produced locally yet"; n_missing=$((n_missing + 1)); continue
  fi
  size="$(stat -c %s "$local_path")"
  if [ "${REMOTE_SIZE[$rel]:-}" = "$size" ]; then
    echo "[have] $rel: already on the cluster ($(human "$size"))"; n_have=$((n_have + 1)); continue
  fi
  SEND+=("$rel"); total=$((total + size))
done

if [ ${#SEND[@]} -eq 0 ]; then
  echo "nothing to send: $n_have file(s) already there, $n_missing missing locally"
  exit 0
fi

echo "==================================================================="
echo "sending ${#SEND[@]} file(s), $(human "$total") -> $CLUSTER_HOST:$REMOTE"
printf '  %s\n' "${SEND[@]}"
echo "==================================================================="

if [ "$DRY_RUN" -eq 1 ]; then
  echo "[dry] nothing transferred"
  exit 0
fi

if [ "$USE_SCP" -eq 1 ]; then
  # scp cannot create the remote tree, so do it first, then send file by file.
  mapfile -t subdirs < <(for rel in "${SEND[@]}"; do d="$(dirname "$rel")"; [ "$d" != "." ] && echo "$d"; done | sort -u)
  ssh "${SSH_OPTS[@]}" "$CLUSTER_HOST" \
    "mkdir -p '$REMOTE' $(printf "'%s' " ${subdirs+"${subdirs[@]/#/$REMOTE/}"})"
  for rel in "${SEND[@]}"; do
    echo "--- $rel"
    scp "${SSH_OPTS[@]}" -p "$DATA_DIR/$rel" "$CLUSTER_HOST:$REMOTE/$rel"
  done
else
  # --files-from keeps relative paths; --partial resumes interrupted files.
  ssh "${SSH_OPTS[@]}" "$CLUSTER_HOST" "mkdir -p '$REMOTE'"
  printf '%s\n' "${SEND[@]}" | rsync -av --partial --progress \
    -e "ssh ${SSH_OPTS[*]}" --files-from=- "$DATA_DIR/" "$CLUSTER_HOST:$REMOTE/"
fi

echo "==================================================================="
echo "done: ${#SEND[@]} sent, $n_have already there, $n_missing missing locally"
echo "      -> $CLUSTER_HOST:$REMOTE"
[ "$WHAT" = "metrics" ] && cat <<EOF
next, on the cluster:
  export DATA_DIR=$REMOTE
  cd <repo>/02_integration_benchmark/02_4_metrics && mkdir -p logs
  bash run_all_metrics.sh --slurm${SCALING_FILTER:+ --scaling $SCALING_FILTER} --dry-run
EOF
exit 0
