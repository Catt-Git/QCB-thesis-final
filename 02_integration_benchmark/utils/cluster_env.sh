#!/usr/bin/env bash
#
# 02 utils: the cluster data root, shared by sync_to_cluster.sh (upload) and the
# cluster-side wrappers (stale-copy warning), so the two cannot drift apart.
#
# Sourced, never executed:
#   . "$SCRIPT_DIR/../utils/cluster_env.sh"
#
#   CLUSTER_HOST      user@hostname; only sync_to_cluster.sh needs it
#   CLUSTER_DATA_DIR  the data root on the cluster

CLUSTER_DATA_DIR="${CLUSTER_DATA_DIR:-/users/genomics/albertoc/Tesi/hopes_and_dreams/datasets}"

# DATA_DIR is required; warns if it differs from an existing cluster root.
resolve_data_dir() {
  if [ -z "${DATA_DIR:-}" ]; then
    echo "set DATA_DIR to the directory holding the datasets" >&2
    return 1
  fi
  if [ "${DATA_DIR%/}" != "${CLUSTER_DATA_DIR%/}" ] && [ -d "$CLUSTER_DATA_DIR" ]; then
    echo "[warn] DATA_DIR=$DATA_DIR, but the cluster data root is" >&2
    echo "[warn]   $CLUSTER_DATA_DIR" >&2
    echo "[warn] and it exists on this machine. sync_to_cluster.sh fills the latter," >&2
    echo "[warn] so one of the two is stale - check which before trusting the scores." >&2
  fi
  return 0
}
