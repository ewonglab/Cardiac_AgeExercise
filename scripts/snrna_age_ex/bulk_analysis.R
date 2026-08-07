# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Run per-animal, cell-type pseudobulk DESeq2 robustness analyses.
# Steps: 1) Aggregate raw counts by mouse; 2) fit sex-adjusted DESeq2 contrasts; 3) export per-cell-type and combined results.
# =======================================
# Source reformatted from snRNA_AgeEx/Bulk_analysis(1).ipynb.
# Notebook outputs and metadata are intentionally excluded.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))
#BULK


library(Seurat)
library(SingleCellExperiment)
library(MAST)
library(future)

options(mc.cores =28)


setwd(publication_path("HEART_SNRNA_ROOT", "data/2022-10-17.SharedObject"))
Mouse = readRDS(file = "2025-09-23_Full_RNA_DateAnnotated.rds")


Validindv =  c("MAE2_1_M","MAE2_1_F","MAS1_4_F","MAS1_4_M","MAS3_2_M","MAS3_2_F",
 "YE2_1_F","YE2_1_M","YE4_3_F","YE4_3_M","YS1_4_M","YS3_2_M","YS3_2_F")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(readr)
})

options(mc.cores = 28)

# -----------------------------
# Paths / settings
# -----------------------------


outdir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-06-23.DESeq2_pseudobulk_by_celltype")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

freq_expressed <- 0.2
FCTHRESHOLD <- log2(1.5)
min_cells_per_mouse <- 20
min_total_pseudobulk_count <- 10

Validindv <- c(
  "MAE2_1_M", "MAE2_1_F",
  "MAS1_4_F", "MAS1_4_M",
  "MAS3_2_M", "MAS3_2_F",
  "YE2_1_F",  "YE2_1_M",
  "YE4_3_F",  "YE4_3_M",
  "YS1_4_M",
  "YS3_2_M",  "YS3_2_F"
)

# Contrasts matching your MAST script logic
contrast_defs <- list(
  MAE_vs_MAS = c(case = "MAE", ctrl = "MAS"),  # Aged exercise vs aged sedentary
  YE_vs_YS   = c(case = "YE",  ctrl = "YS"),   # Young exercise vs young sedentary
  MAE_vs_YE  = c(case = "MAE", ctrl = "YE"),   # Aged exercise vs young exercise
  MAS_vs_YS  = c(case = "MAS", ctrl = "YS")    # Aged sedentary vs young sedentary
)

# Map shorthand to age_condition labels
group_alias <- c(
  MAE = "Aged_Ex",
  MAS = "Aged_Sed",
  YE  = "Young_Ex",
  YS  = "Young_Sed"
)

# -----------------------------
# Load object
# -----------------------------
setwd(outdir )


RNA <- subset(Mouse, subset = mouseID %in% Validindv)

meta <- RNA[[]]

stopifnot("mouseID" %in% colnames(meta))
stopifnot("CellType" %in% colnames(meta))
stopifnot("age_condition" %in% colnames(meta))
stopifnot("sex" %in% colnames(meta))

# -----------------------------
# Get raw counts
# -----------------------------
get_counts <- function(seurat_obj, assay = "RNA") {
  if ("Assay5" %in% class(seurat_obj[[assay]])) {
    return(GetAssayData(seurat_obj, assay = assay, layer = "counts"))
  } else {
    return(seurat_obj[[assay]]@counts)
  }
}

counts_all <- get_counts(RNA, assay = "RNA")

# -----------------------------
# Cell type sets
# -----------------------------
celltype_sets <- setNames(
  as.list(sort(unique(meta$CellType))),
  sort(unique(meta$CellType))
)

# Combined macrophage population matching your MAST example
if (all(c("M1", "M2", "MacFib") %in% meta$CellType)) {
  celltype_sets[["RealMacrophage"]] <- c("M1", "M2", "MacFib")
}

# -----------------------------
# Helper functions
# -----------------------------
resolve_group <- function(x, available_groups) {
  if (x %in% available_groups) return(x)
  if (x %in% names(group_alias) && group_alias[[x]] %in% available_groups) return(group_alias[[x]])
  return(NA_character_)
}

make_pseudobulk <- function(counts_mat, cell_meta) {
  sample_id <- factor(cell_meta$mouseID)

  design_mat <- sparse.model.matrix(~ 0 + sample_id)
  colnames(design_mat) <- levels(sample_id)

  pb_counts <- counts_mat %*% design_mat
  pb_counts <- as.matrix(pb_counts)
  pb_counts <- round(pb_counts)
  storage.mode(pb_counts) <- "integer"

  sample_meta <- cell_meta %>%
    as.data.frame() %>%
    group_by(mouseID) %>%
    summarise(
      group = first(age_condition),
      sex = first(sex),
      n_cells = n(),
      .groups = "drop"
    ) %>%
    as.data.frame()

  rownames(sample_meta) <- sample_meta$mouseID

  pb_counts <- pb_counts[, rownames(sample_meta), drop = FALSE]

  list(
    counts = pb_counts,
    coldata = sample_meta
  )
}

run_deseq_celltype <- function(cell_name, cell_values) {
  message("========== ", cell_name, " ==========")

  cells_use <- rownames(meta)[meta$CellType %in% cell_values]

  if (length(cells_use) == 0) {
    message("Skipping ", cell_name, ": no cells found")
    return(NULL)
  }

  cell_meta <- meta[cells_use, , drop = FALSE]
  counts_cell <- counts_all[, cells_use, drop = FALSE]

  # -----------------------------
  # Gene filtering matching MAST logic
  # Keep genes expressed in >20% of cells in this cell type/set
  # Equivalent logic to freq(scaRaw) > freq_expressed
  # -----------------------------
  pct_expressed <- Matrix::rowMeans(counts_cell > 0)
  keep_genes <- pct_expressed > freq_expressed

  counts_cell <- counts_cell[keep_genes, , drop = FALSE]
  pct_expressed <- pct_expressed[keep_genes]

  if (nrow(counts_cell) == 0) {
    message("Skipping ", cell_name, ": no genes pass expression-frequency filter")
    return(NULL)
  }

  # -----------------------------
  # Pseudobulk aggregation
  # -----------------------------
  pb <- make_pseudobulk(counts_cell, cell_meta)

  pb_counts <- pb$counts
  sample_meta <- pb$coldata

  # Filter low-cell pseudobulk samples
  keep_samples <- sample_meta$n_cells >= min_cells_per_mouse
  sample_meta <- sample_meta[keep_samples, , drop = FALSE]
  pb_counts <- pb_counts[, rownames(sample_meta), drop = FALSE]

  if (ncol(pb_counts) < 4) {
    #message("Skipping ", cell_name, ": fewer than 4 usable pseudobulk samples")
    #return(NULL)
  }

  # Remove genes with near-zero total pseudobulk counts
  keep_count <- rowSums(pb_counts) >= min_total_pseudobulk_count
  pb_counts <- pb_counts[keep_count, , drop = FALSE]
  pct_expressed <- pct_expressed[rownames(pb_counts)]

  if (nrow(pb_counts) == 0) {
    message("Skipping ", cell_name, ": no genes left after pseudobulk count filter")
    return(NULL)
  }

  sample_meta$group <- factor(sample_meta$group)
  sample_meta$sex <- factor(sample_meta$sex)

  available_groups <- levels(sample_meta$group)

  # Save pseudobulk counts and metadata
  write.csv(
    pb_counts,
    file = file.path(outdir, paste0(cell_name, "_pseudobulk_counts.csv"))
  )

  write.csv(
    sample_meta,
    file = file.path(outdir, paste0(cell_name, "_pseudobulk_metadata.csv")),
    row.names = TRUE
  )

  # -----------------------------
  # DESeq2 model
  # -----------------------------
  dds <- DESeqDataSetFromMatrix(
    countData = pb_counts,
    colData = sample_meta,
    design = ~ sex + group
  )

  dds <- DESeq(dds, quiet = TRUE)

  all_res <- list()

  for (contrast_name in names(contrast_defs)) {
    case_raw <- contrast_defs[[contrast_name]][["case"]]
    ctrl_raw <- contrast_defs[[contrast_name]][["ctrl"]]

    case_group <- resolve_group(case_raw, available_groups)
    ctrl_group <- resolve_group(ctrl_raw, available_groups)

    if (is.na(case_group) || is.na(ctrl_group)) {
      message("Skipping ", cell_name, " ", contrast_name, ": missing group")
      next
    }

    n_case <- sum(sample_meta$group == case_group)
    n_ctrl <- sum(sample_meta$group == ctrl_group)

    if (n_case < 2 || n_ctrl < 2) {
      #message(
      #  "Skipping ", cell_name, " ", contrast_name,
      #  ": too few samples: ", case_group, " n=", n_case,
      #  ", ", ctrl_group, " n=", n_ctrl
      #)
      #next
    }

    res <- results(
      dds,
      contrast = c("group", case_group, ctrl_group),
      alpha = 0.1
    )

    res_df <- as.data.frame(res) %>%
      rownames_to_column("gene") %>%
      mutate(
        celltype = cell_name,
        contrast = contrast_name,
        case = case_group,
        ctrl = ctrl_group,
        n_case = n_case,
        n_ctrl = n_ctrl,
        pct_expressed = pct_expressed[gene],
        significant_FDR_0.1 = !is.na(padj) & padj < 0.1,
        significant_FDR_0.1_logFC = !is.na(padj) & padj < 0.1 & abs(log2FoldChange) > FCTHRESHOLD
      ) %>%
      arrange(padj)

    out_csv <- file.path(outdir, paste0("DESeq2_", cell_name, "_", contrast_name, ".csv"))
    write.csv(res_df, out_csv, row.names = FALSE)

    all_res[[contrast_name]] <- res_df
  }

  if (length(all_res) == 0) {
    return(NULL)
  }

  bind_rows(all_res)
}

# -----------------------------
# Run across cell types
# -----------------------------
all_results <- list()

for (cell_name in names(celltype_sets)) {
  res <- tryCatch(
    run_deseq_celltype(cell_name, celltype_sets[[cell_name]]),
    error = function(e) {
      message("FAILED: ", cell_name)
      message(e$message)
      return(NULL)
    }
  )

  if (!is.null(res)) {
    all_results[[cell_name]] <- res
  }
}

combined_results <- bind_rows(all_results)

write.csv(
  combined_results,
  file = file.path(outdir, "DESeq2_all_celltypes_all_contrasts.csv"),
  row.names = FALSE
)

# -----------------------------
# DEG summary table
# -----------------------------
deg_summary <- combined_results %>%
  group_by(celltype, contrast) %>%
  summarise(
    n_genes_tested = n(),
    n_FDR_0.1 = sum(significant_FDR_0.1, na.rm = TRUE),
    n_FDR_0.1_logFC = sum(significant_FDR_0.1_logFC, na.rm = TRUE),
    n_up_FDR_0.1_logFC = sum(significant_FDR_0.1_logFC & log2FoldChange > 0, na.rm = TRUE),
    n_down_FDR_0.1_logFC = sum(significant_FDR_0.1_logFC & log2FoldChange < 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(celltype, contrast)

write.csv(
  deg_summary,
  file = file.path(outdir, "DESeq2_DEG_summary.csv"),
  row.names = FALSE
)

print(deg_summary)

outdir
