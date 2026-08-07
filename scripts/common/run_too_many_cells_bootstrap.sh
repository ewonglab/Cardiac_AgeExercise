#!/usr/bin/env bash
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Build a bootstrap TooManyCells hierarchy for the configured single-cell input.
# Steps: 1) Stage the input; 2) launch the containerised bootstrap workflow; 3) save the hierarchy and diagnostic outputs.
# =======================================
# Build one bootstrap TooManyCells tree after its matrix and labels are written.
set -euo pipefail

: "${TOO_MANY_CELLS_IMAGE:?Set TOO_MANY_CELLS_IMAGE in config/.env before running this script.}"

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <bootstrap-replicate-directory>" >&2
  exit 2
fi

cd "$1"
mkdir -p tree
singularity run "${TOO_MANY_CELLS_IMAGE}" make-tree \
  --matrix-path *data.csv \
  --labels-file *label.csv \
  --normalization NoneNorm \
  --draw-collection PieChart \
  --output tree \
  > tree/cluster_clusters.csv
