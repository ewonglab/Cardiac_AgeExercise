#!/usr/bin/env bash
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Build the shared TooManyCells hierarchies for the configured macrophage input.
# Steps: 1) Check configured inputs; 2) launch the full, cropped and group hierarchy builds; 3) write downstream hierarchy outputs.
# =======================================
# Generate the shared TooManyCells hierarchies required by the macrophage analyses.
set -euo pipefail

: "${HEART_SNRNA_ROOT:?Set HEART_SNRNA_ROOT in config/.env before running this script.}"
: "${TOO_MANY_CELLS_IMAGE:?Set TOO_MANY_CELLS_IMAGE in config/.env before running this script.}"

analysis_dir="${HEART_SNRNA_ROOT}/analysis/2023-03-15.Macro"

mkdir -p \
  "${analysis_dir}/too-many-cells/Macrophages_Full" \
  "${analysis_dir}/too-many-cells/Macrophages_full_Crop" \
  "${analysis_dir}/too-many-cells/Macrophages_Full_group"

singularity run -B "${analysis_dir}:${analysis_dir}" "${TOO_MANY_CELLS_IMAGE}" make-tree \
  --matrix-path "${analysis_dir}/Macrophages_Full_data.csv" \
  --labels-file "${analysis_dir}/Macrophage_Full_Cluster_labels.csv" \
  --normalization NoneNorm \
  --draw-collection PieChart \
  --output "${analysis_dir}/too-many-cells/Macrophages_Full" \
  --draw-scale-saturation 10 \
  --draw-node-number \
  > "${analysis_dir}/too-many-cells/Macrophages_Full/cluster_clusters.csv"

singularity run -B "${analysis_dir}:${analysis_dir}" "${TOO_MANY_CELLS_IMAGE}" make-tree \
  --labels-file "${analysis_dir}/Macrophage_Full_Cluster_labels.csv" \
  --normalization NoneNorm \
  --draw-collection PieChart \
  --output "${analysis_dir}/too-many-cells/Macrophages_full_Crop" \
  --prior "${analysis_dir}/too-many-cells/Macrophages_Full" \
  --smart-cutoff 0.5 \
  --min-size 100 \
  --draw-scale-saturation 10 \
  --draw-node-number \
  > "${analysis_dir}/too-many-cells/Macrophages_full_Crop/cluster_clusters.csv"

singularity run -B "${analysis_dir}:${analysis_dir}" "${TOO_MANY_CELLS_IMAGE}" make-tree \
  --labels-file "${analysis_dir}/Macrophage_Full_group_labels.csv" \
  --normalization NoneNorm \
  --draw-collection PieChart \
  --output "${analysis_dir}/too-many-cells/Macrophages_Full_group" \
  --prior "${analysis_dir}/too-many-cells/Macrophages_Full" \
  --smart-cutoff 0.5 \
  --min-size 100 \
  --draw-scale-saturation 10 \
  --draw-node-number \
  > "${analysis_dir}/too-many-cells/Macrophages_Full_group/cluster_clusters.csv"
