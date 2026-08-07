# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Test AAV dnAP-1 transcriptional effects in female nuclei with MAST.
# Steps: 1) Load the processed object; 2) fit cell-type-specific models; 3) export differential-expression results and plots.
# =======================================
# Source reformatted from dnAP-1_snRNA/MAST_All_female.ipynb.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

#MAST Sample Removed

setwd(publication_path("DNAP1_SNRNA_ROOT", "Process/"))

setwd(publication_path("DNAP1_SNRNA_ROOT", "Process/"))
seu = readRDS("Celltype_annotated_qc_filtered.rds")

suppressPackageStartupMessages({
  library(Seurat)
  library(MAST)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(data.table)
})

# ============================================================
# USER SETTINGS
# ============================================================

obj <- seu

# New output directory so this does not overwrite your previous run
out_dir <- "MAST_DEG_by_celltype_CDR_batch_removed_male_sample"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

assay_use <- "RNA"

# Metadata column candidates
celltype_col_candidates  <- c("CellType")
condition_col_candidates <- c("Group")
batch_col_candidates     <- c("Orig_Pool")

# Used only for removing the male sample
sample_col_candidates <- c(
  "Sample_ID", "sample_id", "Sample", "sample",
  "orig.ident", "orig_ident",
  "Mouse_ID", "MouseID", "Animal_ID", "AnimalID",
  "Donor", "donor", "Library", "library"
)

sex_col_candidates <- c("Sex", "sex", "Gender", "gender")

# Set this to the exact sample ID if needed.
# If left as "AUTO", the code will remove the only sample whose metadata says male.
male_sample_to_remove <- "AUTO"

# Main comparison
ref_group  <- "AAV9-cTnT"
test_group <- "AAV9-AP-1_Dosage_2"

min_pct_expr <- 0.20
min_cells_per_group <- 20

fdr_cutoff <- 0.05
log2fc_cutoff <- 0.25

# ============================================================
# HELPER FUNCTIONS
# ============================================================

pick_col <- function(meta, candidates, required = TRUE, label = "column") {
  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) > 0) return(hit[1])

  if (required) {
    stop(
      "Could not find ", label, ". Tried: ",
      paste(candidates, collapse = ", "),
      "\nAvailable metadata columns:\n",
      paste(colnames(meta), collapse = ", ")
    )
  }

  NA_character_
}

safe_name <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(x == "", "unknown", x)
}

get_assay_mat <- function(seurat_obj, assay, data_type = c("data", "counts")) {
  data_type <- match.arg(data_type)

  mat <- tryCatch(
    {
      GetAssayData(seurat_obj, assay = assay, layer = data_type)
    },
    error = function(e1) {
      GetAssayData(seurat_obj, assay = assay, slot = data_type)
    }
  )

  mat
}

extract_mast_results <- function(zlm_fit, coef_name) {
  s <- summary(zlm_fit, doLRT = coef_name)
  dt <- as.data.table(s$datatable)

  hurdle <- dt[
    contrast == coef_name & component == "H",
    .(gene = primerid, mast_p = `Pr(>Chisq)`)
  ]

  logfc <- dt[
    contrast == coef_name & component == "logFC",
    .(
      gene = primerid,
      mast_coef_logFC = coef,
      mast_ci_low = ci.lo,
      mast_ci_high = ci.hi
    )
  ]

  continuous <- dt[
    contrast == coef_name & component == "C",
    .(
      gene = primerid,
      mast_continuous_coef = coef,
      mast_continuous_p = `Pr(>Chisq)`
    )
  ]

  discrete <- dt[
    contrast == coef_name & component == "D",
    .(
      gene = primerid,
      mast_discrete_coef = coef,
      mast_discrete_p = `Pr(>Chisq)`
    )
  ]

  res <- hurdle %>%
    left_join(logfc, by = "gene") %>%
    left_join(continuous, by = "gene") %>%
    left_join(discrete, by = "gene") %>%
    mutate(mast_FDR = p.adjust(mast_p, method = "BH"))

  res
}

detect_male_sample_to_remove <- function(meta, sample_col, sex_col) {
  sample_sex_summary <- meta %>%
    as.data.frame() %>%
    mutate(
      sample_id = as.character(.data[[sample_col]]),
      sex_value = tolower(as.character(.data[[sex_col]]))
    ) %>%
    group_by(sample_id) %>%
    summarise(
      n_cells = n(),
      sex_values = paste(sort(unique(sex_value)), collapse = ";"),
      n_male_cells = sum(sex_value %in% c("m", "male")),
      n_female_cells = sum(sex_value %in% c("f", "female")),
      .groups = "drop"
    ) %>%
    mutate(is_male_sample = n_male_cells > 0 & n_female_cells == 0)

  write_csv(
    sample_sex_summary,
    file.path(out_dir, "sample_sex_summary_before_removal.csv")
  )

  male_samples <- sample_sex_summary %>%
    filter(is_male_sample) %>%
    pull(sample_id)

  if (length(male_samples) != 1) {
    stop(
      "AUTO removal did not find exactly one male sample.\n",
      "Male-like samples found: ", paste(male_samples, collapse = ", "), "\n",
      "Set male_sample_to_remove manually to the exact sample ID.\n",
      "A summary was written to: ",
      file.path(out_dir, "sample_sex_summary_before_removal.csv")
    )
  }

  male_samples
}

run_mast_one_celltype <- function(
  ct,
  obj,
  meta,
  data_mat,
  counts_mat,
  celltype_col,
  condition_col,
  batch_col,
  ref_group,
  test_group,
  min_pct_expr,
  min_cells_per_group,
  out_dir
) {
  message("\n============================================================")
  message("Running cell type: ", ct)
  message("============================================================")

  cells_use <- rownames(meta)[
    as.character(meta[[celltype_col]]) == ct &
      as.character(meta[[condition_col]]) %in% c(ref_group, test_group)
  ]

  if (length(cells_use) == 0) {
    message("Skipping ", ct, ": no cells found.")
    return(NULL)
  }

  meta_ct <- meta[cells_use, , drop = FALSE]

  n_ref <- sum(as.character(meta_ct[[condition_col]]) == ref_group)
  n_test <- sum(as.character(meta_ct[[condition_col]]) == test_group)

  if (n_ref < min_cells_per_group || n_test < min_cells_per_group) {
    message(
      "Skipping ", ct,
      ": too few cells. ",
      ref_group, " n=", n_ref, ", ",
      test_group, " n=", n_test
    )
    return(NULL)
  }

  counts_ct <- counts_mat[, cells_use, drop = FALSE]
  data_ct   <- data_mat[, cells_use, drop = FALSE]

  ref_cells  <- rownames(meta_ct)[as.character(meta_ct[[condition_col]]) == ref_group]
  test_cells <- rownames(meta_ct)[as.character(meta_ct[[condition_col]]) == test_group]

  pct_ref  <- Matrix::rowMeans(counts_ct[, ref_cells, drop = FALSE] > 0)
  pct_test <- Matrix::rowMeans(counts_ct[, test_cells, drop = FALSE] > 0)

  keep_genes <- names(which(pct_ref >= min_pct_expr | pct_test >= min_pct_expr))

  if (length(keep_genes) < 10) {
    message("Skipping ", ct, ": fewer than 10 genes pass expression cutoff.")
    return(NULL)
  }

  data_filt <- data_ct[keep_genes, , drop = FALSE]

  meta_ct$condition_mast_original <- as.character(meta_ct[[condition_col]])
  meta_ct$condition_mast <- factor(
    make.names(meta_ct$condition_mast_original),
    levels = make.names(c(ref_group, test_group))
  )

  meta_ct$batch_mast_original <- as.character(meta_ct[[batch_col]])
  meta_ct$batch_mast_original[is.na(meta_ct$batch_mast_original) | meta_ct$batch_mast_original == ""] <- "Unknown"
  meta_ct$batch_mast <- factor(make.names(meta_ct$batch_mast_original))

  meta_ct$CDR_raw <- Matrix::colSums(counts_mat[, rownames(meta_ct), drop = FALSE] > 0)
  meta_ct$CDR <- as.numeric(scale(meta_ct$CDR_raw))

  if (anyNA(meta_ct$CDR)) {
    meta_ct$CDR <- meta_ct$CDR_raw
  }

  formula_use <- if (nlevels(meta_ct$batch_mast) > 1) {
    ~ condition_mast + CDR + batch_mast
  } else {
    message("Only one batch present for ", ct, "; fitting condition + CDR only.")
    ~ condition_mast + CDR
  }

  exprs_mat <- as.matrix(data_filt)

  fdata <- data.frame(
    primerid = rownames(exprs_mat),
    row.names = rownames(exprs_mat),
    stringsAsFactors = FALSE
  )

  cdata <- meta_ct %>%
    as.data.frame() %>%
    mutate(wellKey = rownames(meta_ct)) %>%
    as.data.frame()

  rownames(cdata) <- cdata$wellKey

  sca <- FromMatrix(
    exprsArray = exprs_mat,
    cData = cdata,
    fData = fdata
  )

  zlm_fit <- zlm(
    formula = formula_use,
    sca = sca,
    method = "bayesglm",
    ebayes = TRUE
  )

  coef_names <- colnames(zlm_fit@coefD)
  coef_target <- paste0("condition_mast", make.names(test_group))

  if (!coef_target %in% coef_names) {
    condition_hits <- grep("^condition_mast", coef_names, value = TRUE)
    stop(
      "Could not find expected MAST condition coefficient for ", ct, ".\n",
      "Expected: ", coef_target, "\n",
      "Available coefficients: ", paste(coef_names, collapse = ", "), "\n",
      "Condition-like coefficients: ", paste(condition_hits, collapse = ", ")
    )
  }

  mast_res <- extract_mast_results(zlm_fit, coef_target)

  avg_ref_norm  <- Matrix::rowMeans(expm1(data_filt[, ref_cells, drop = FALSE]))
  avg_test_norm <- Matrix::rowMeans(expm1(data_filt[, test_cells, drop = FALSE]))

  avg_log2FC <- log2((avg_test_norm + 1e-9) / (avg_ref_norm + 1e-9))

  gene_stats <- tibble(
    gene = keep_genes,
    celltype = ct,
    comparison = paste0(test_group, "_vs_", ref_group),
    n_cells_ref = n_ref,
    n_cells_test = n_test,
    pct_ref = as.numeric(pct_ref[keep_genes]),
    pct_test = as.numeric(pct_test[keep_genes]),
    pct_max = pmax(pct_ref[keep_genes], pct_test[keep_genes]),
    avg_expr_ref = as.numeric(avg_ref_norm[keep_genes]),
    avg_expr_test = as.numeric(avg_test_norm[keep_genes]),
    avg_log2FC = as.numeric(avg_log2FC[keep_genes])
  )

  final_res <- gene_stats %>%
    left_join(mast_res, by = "gene") %>%
    arrange(mast_FDR, mast_p, desc(abs(avg_log2FC)))

  ct_safe <- safe_name(ct)

  write_csv(
    final_res,
    file.path(out_dir, paste0("MAST_", ct_safe, "_", safe_name(test_group), "_vs_", safe_name(ref_group), "_all_genes.csv"))
  )

  write_csv(
    final_res %>%
      filter(!is.na(mast_FDR), mast_FDR < fdr_cutoff, abs(avg_log2FC) >= log2fc_cutoff),
    file.path(out_dir, paste0("MAST_", ct_safe, "_", safe_name(test_group), "_vs_", safe_name(ref_group), "_sig_FDR", fdr_cutoff, "_log2FC", log2fc_cutoff, ".csv"))
  )

  message(
    "Finished ", ct,
    ": tested ", nrow(final_res), " genes; ",
    sum(final_res$mast_FDR < fdr_cutoff, na.rm = TRUE), " FDR < ", fdr_cutoff
  )

  final_res
}

# ============================================================
# PREPARE METADATA
# ============================================================

DefaultAssay(obj) <- assay_use

meta <- obj@meta.data
meta$cell_barcode <- rownames(meta)

celltype_col  <- pick_col(meta, celltype_col_candidates, required = TRUE, label = "cell type column")
condition_col <- pick_col(meta, condition_col_candidates, required = TRUE, label = "condition column")
batch_col     <- pick_col(meta, batch_col_candidates, required = TRUE, label = "batch column")

sample_col <- pick_col(meta, sample_col_candidates, required = FALSE, label = "sample column")
sex_col    <- pick_col(meta, sex_col_candidates, required = FALSE, label = "sex column")

message("Using cell type column: ", celltype_col)
message("Using condition column: ", condition_col)
message("Using batch column: ", batch_col)

if (is.na(sample_col)) {
  stop(
    "Could not find a sample column for removing the male sample.\n",
    "Tried: ", paste(sample_col_candidates, collapse = ", "), "\n",
    "Available metadata columns:\n",
    paste(colnames(meta), collapse = ", ")
  )
}

message("Using sample column for removal: ", sample_col)

if (identical(male_sample_to_remove, "AUTO")) {
  if (is.na(sex_col)) {
    stop(
      "male_sample_to_remove is set to AUTO, but no sex column was found.\n",
      "Either add the correct sex column to sex_col_candidates or set male_sample_to_remove manually."
    )
  }

  message("Using sex column for AUTO male detection: ", sex_col)

  male_sample_to_remove <- detect_male_sample_to_remove(
    meta = meta,
    sample_col = sample_col,
    sex_col = sex_col
  )
}

message("Removing male sample: ", male_sample_to_remove)

cells_remove <- rownames(meta)[as.character(meta[[sample_col]]) == male_sample_to_remove]
cells_keep <- setdiff(rownames(meta), cells_remove)

if (length(cells_remove) == 0) {
  stop(
    "No cells matched male_sample_to_remove = ", male_sample_to_remove,
    " in column ", sample_col
  )
}

removal_summary <- tibble(
  removed_sample = male_sample_to_remove,
  sample_col = sample_col,
  sex_col = ifelse(is.na(sex_col), NA_character_, sex_col),
  n_cells_before = nrow(meta),
  n_cells_removed = length(cells_remove),
  n_cells_after = length(cells_keep)
)

write_csv(
  removal_summary,
  file.path(out_dir, "removed_male_sample_summary.csv")
)



obj <- subset(obj, cells = cells_keep)

# Refresh metadata after filtering
meta <- obj@meta.data
meta$cell_barcode <- rownames(meta)

message("Cells before removal: ", removal_summary$n_cells_before)
message("Cells removed: ", removal_summary$n_cells_removed)
message("Cells after removal: ", removal_summary$n_cells_after)

# ============================================================
# CHECK GROUPS AFTER REMOVAL
# ============================================================

if (!all(c(ref_group, test_group) %in% as.character(meta[[condition_col]]))) {
  stop(
    "The requested groups were not both found after removing the male sample.\n",
    "Requested ref_group: ", ref_group, "\n",
    "Requested test_group: ", test_group, "\n",
    "Available groups:\n",
    paste(sort(unique(as.character(meta[[condition_col]]))), collapse = ", ")
  )
}

group_summary_after_removal <- meta %>%
  as.data.frame() %>%
  count(.data[[condition_col]], .data[[batch_col]], name = "n_cells") %>%
  rename(condition = 1, batch = 2)

# ============================================================
# PREPARE MATRICES AFTER REMOVAL
# ============================================================

data_mat <- get_assay_mat(obj, assay = assay_use, data_type = "data")
counts_mat <- get_assay_mat(obj, assay = assay_use, data_type = "counts")

common_cells <- intersect(colnames(data_mat), rownames(meta))
common_cells <- intersect(common_cells, colnames(counts_mat))

data_mat <- data_mat[, common_cells, drop = FALSE]
counts_mat <- counts_mat[, common_cells, drop = FALSE]
meta <- meta[common_cells, , drop = FALSE]

celltypes <- sort(unique(as.character(meta[[celltype_col]])))
celltypes <- celltypes[!is.na(celltypes) & celltypes != ""]

message("Cell types to test:")
print(celltypes)

# ============================================================
# RUN MAST ACROSS CELL TYPES
# ============================================================

all_results <- lapply(celltypes, function(ct) {
  tryCatch(
    run_mast_one_celltype(
      ct = ct,
      obj = obj,
      meta = meta,
      data_mat = data_mat,
      counts_mat = counts_mat,
      celltype_col = celltype_col,
      condition_col = condition_col,
      batch_col = batch_col,
      ref_group = ref_group,
      test_group = test_group,
      min_pct_expr = min_pct_expr,
      min_cells_per_group = min_cells_per_group,
      out_dir = out_dir
    ),
    error = function(e) {
      message("FAILED for cell type: ", ct)
      message("Error: ", e$message)

      writeLines(
        e$message,
        con = file.path(out_dir, paste0("FAILED_", safe_name(ct), ".txt"))
      )

      NULL
    }
  )
})

combined_results <- bind_rows(all_results)

write_csv(
  combined_results,
  file.path(out_dir, paste0("MAST_ALL_CELLTYPES_", safe_name(test_group), "_vs_", safe_name(ref_group), "_all_genes.csv"))
)

combined_sig <- combined_results %>%
  filter(!is.na(mast_FDR), mast_FDR < fdr_cutoff, abs(avg_log2FC) >= log2fc_cutoff) %>%
  arrange(celltype, mast_FDR, mast_p, desc(abs(avg_log2FC)))

write_csv(
  combined_sig,
  file.path(out_dir, paste0("MAST_ALL_CELLTYPES_", safe_name(test_group), "_vs_", safe_name(ref_group), "_sig_FDR", fdr_cutoff, "_log2FC", log2fc_cutoff, ".csv"))
)

summary_table <- combined_results %>%
  group_by(celltype) %>%
  summarise(
    n_genes_tested = n(),
    n_FDR_0_05 = sum(mast_FDR < 0.05, na.rm = TRUE),
    n_FDR_0_10 = sum(mast_FDR < 0.10, na.rm = TRUE),
    n_sig_FDR_0_05_log2FC_0_25 = sum(mast_FDR < 0.05 & abs(avg_log2FC) >= 0.25, na.rm = TRUE),
    n_up_FDR_0_05_log2FC_0_25 = sum(mast_FDR < 0.05 & avg_log2FC >= 0.25, na.rm = TRUE),
    n_down_FDR_0_05_log2FC_0_25 = sum(mast_FDR < 0.05 & avg_log2FC <= -0.25, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sig_FDR_0_05_log2FC_0_25))

write_csv(
  summary_table,
  file.path(out_dir, paste0("MAST_ALL_CELLTYPES_", safe_name(test_group), "_vs_", safe_name(ref_group), "_summary.csv"))
)

print(summary_table)

message("\nDone.")
message("Removed male sample: ", male_sample_to_remove)
message("Output directory: ", normalizePath(out_dir))

DimPlot(seu)

FeaturePlot(seu, "Alb")



DotPlot(seu, "Alb", group.by = "seurat_clusters")

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# ============================================================
# PATHS
# ============================================================

base_dir <- publication_path("DNAP1_SNRNA_ROOT", "Process")

original_dir <- file.path(base_dir, "MAST_DEG_by_celltype_CDR_batch")
removed_dir  <- file.path(base_dir, "MAST_DEG_by_celltype_CDR_batch_removed_male_sample")

compare_dir <- file.path(base_dir, "MAST_DEG_original_vs_removed_male_sample_comparison")
dir.create(compare_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# SETTINGS
# ============================================================

fdr_cutoff <- 0.05
log2fc_cutoff <- 0.25

RNA_pallate = c(
  '#A53D99', '#F99349', '#87CEEB', '#FDDA62', '#808080',
  '#BEBADA', "#A81B1B", "#101369", "#606692", "#3A85A8",
  "#42977E", "#591f0a", "#3b7818", "#B6D66A", "#1CE382",
  "#BF862B", "#CC6A6F", "#E086B5", "#5F1296", "#5F1296",
  "#B784A7"
)

RNA_order = c(
  "Macrophage", "M1", "M2", "M_Ex", "HSPC", "Dendritic",
  "Cardiomyocyte", "Coronary_EC", "Endocardial_EC", "Lymphatic_EC",
  "Fibroblast", "Mesothelial", "B_cell", "T_cell", "Pericytes",
  "Smooth_Muscle", "Cardiac_Neuronal", "Unknown_1", "Unknown_3",
  "ATAC_unknown1", "NK"
)

names(RNA_pallate) <- RNA_order

class_cols <- c(
  "reliable_sig_both"              = "#42977E",
  "sig_both_weaker_after_removal"  = "#B6D66A",
  "likely_male_sample_driven"      = "#A81B1B",
  "original_only_male_sensitive"   = "#F99349",
  "removed_only_unmasked"          = "#3A85A8",
  "discordant_sig_both"            = "#FDDA62",
  "not_sig_either"                 = "#BEBADA"
)

theme_mytheme <- theme(
  plot.background = element_blank(),
  plot.title = element_text(size = 16, face = "bold"),
  plot.subtitle = element_text(size = 12),
  axis.text = element_text(size = 11, color = "black"),
  axis.title = element_text(size = 13, color = "black"),
  axis.line = element_line(color = "black", linewidth = 0.4),
  axis.ticks = element_line(color = "black", linewidth = 0.4),
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 11),
  legend.position = "right",
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  strip.text = element_text(size = 10, face = "bold"),
  plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
)

# ============================================================
# FUNCTIONS
# ============================================================

find_combined_file <- function(dir_path) {
  hits <- list.files(
    dir_path,
    pattern = "^MAST_ALL_CELLTYPES_.*_all_genes\\.csv$",
    full.names = TRUE
  )

  if (length(hits) == 0) {
    stop("No combined all_genes MAST file found in: ", dir_path)
  }

  if (length(hits) > 1) {
    message("Multiple combined files found in ", dir_path)
    message("Using newest file: ", basename(hits[which.max(file.info(hits)$mtime)]))
    hits <- hits[which.max(file.info(hits)$mtime)]
  }

  hits
}

read_mast <- function(file, label) {
  dat <- fread(file) %>%
    as_tibble()

  required_cols <- c("gene", "celltype", "avg_log2FC", "mast_p", "mast_FDR")

  missing_cols <- setdiff(required_cols, colnames(dat))
  if (length(missing_cols) > 0) {
    stop(
      "Missing columns in ", file, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  dat %>%
    mutate(
      run = label,
      gene = as.character(gene),
      celltype = as.character(celltype),
      avg_log2FC = as.numeric(avg_log2FC),
      mast_p = as.numeric(mast_p),
      mast_FDR = as.numeric(mast_FDR),
      sig = !is.na(mast_FDR) &
        mast_FDR < fdr_cutoff &
        !is.na(avg_log2FC) &
        abs(avg_log2FC) >= log2fc_cutoff,
      direction = case_when(
        sig & avg_log2FC > 0 ~ "Up",
        sig & avg_log2FC < 0 ~ "Down",
        TRUE ~ "Not_sig"
      )
    )
}

safe_neglog10 <- function(x) {
  -log10(pmax(as.numeric(x), 1e-300))
}

save_pdf <- function(plot, filename, width = 8.27, height = 8.27) {
  ggsave(
    filename = file.path(compare_dir, filename),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    useDingbats = FALSE
  )
}

# ============================================================
# READ INPUTS
# ============================================================

original_file <- find_combined_file(original_dir)
removed_file  <- find_combined_file(removed_dir)

message("Original file: ", original_file)
message("Removed-male file: ", removed_file)

original <- read_mast(original_file, "Original")
removed  <- read_mast(removed_file, "Removed_male_sample")

writeLines(
  c(
    paste0("Original file: ", original_file),
    paste0("Removed male sample file: ", removed_file),
    paste0("FDR cutoff: ", fdr_cutoff),
    paste0("log2FC cutoff: ", log2fc_cutoff)
  ),
  con = file.path(compare_dir, "comparison_inputs.txt")
)

# ============================================================
# JOIN ORIGINAL VS REMOVED
# ============================================================

comp <- original %>%
  dplyr::select(
    gene,
    celltype,
    comparison,
    avg_log2FC_original = avg_log2FC,
    mast_p_original = mast_p,
    mast_FDR_original = mast_FDR,
    pct_ref_original = any_of("pct_ref"),
    pct_test_original = any_of("pct_test"),
    n_cells_ref_original = any_of("n_cells_ref"),
    n_cells_test_original = any_of("n_cells_test"),
    sig_original = sig,
    direction_original = direction
  ) %>%
  full_join(
    removed %>%
      dplyr::select(
        gene,
        celltype,
        comparison,
        avg_log2FC_removed = avg_log2FC,
        mast_p_removed = mast_p,
        mast_FDR_removed = mast_FDR,
        pct_ref_removed = any_of("pct_ref"),
        pct_test_removed = any_of("pct_test"),
        n_cells_ref_removed = any_of("n_cells_ref"),
        n_cells_test_removed = any_of("n_cells_test"),
        sig_removed = sig,
        direction_removed = direction
      ),
    by = c("gene", "celltype", "comparison")
  ) %>%
  mutate(
    sig_original = replace_na(sig_original, FALSE),
    sig_removed = replace_na(sig_removed, FALSE),

    neglog10_FDR_original = safe_neglog10(mast_FDR_original),
    neglog10_FDR_removed = safe_neglog10(mast_FDR_removed),

    same_direction = sign(avg_log2FC_original) == sign(avg_log2FC_removed),
    same_direction = ifelse(is.na(same_direction), FALSE, same_direction),

    abs_log2FC_original = abs(avg_log2FC_original),
    abs_log2FC_removed = abs(avg_log2FC_removed),

    log2FC_retention = abs_log2FC_removed / pmax(abs_log2FC_original, 1e-9),
    abs_delta_log2FC = abs(avg_log2FC_removed - avg_log2FC_original),

    evidence_drop_after_removal = neglog10_FDR_original - neglog10_FDR_removed,
    abs_log2FC_drop_after_removal = abs_log2FC_original - abs_log2FC_removed,

    male_sensitivity_score = case_when(
      sig_original & !sig_removed ~ evidence_drop_after_removal + abs_log2FC_drop_after_removal,
      TRUE ~ NA_real_
    ),

    reliability_class = case_when(
      sig_original & sig_removed & !same_direction ~ "discordant_sig_both",

      sig_original & sig_removed & same_direction & log2FC_retention >= 0.50 ~
        "reliable_sig_both",

      sig_original & sig_removed & same_direction & log2FC_retention < 0.50 ~
        "sig_both_weaker_after_removal",

      sig_original & !sig_removed & log2FC_retention < 0.50 ~
        "likely_male_sample_driven",

      sig_original & !sig_removed ~
        "original_only_male_sensitive",

      !sig_original & sig_removed ~
        "removed_only_unmasked",

      TRUE ~ "not_sig_either"
    )
  )

comp <- comp %>%
  mutate(
    celltype = factor(celltype, levels = intersect(RNA_order, unique(celltype)) %>% c(setdiff(unique(celltype), .))),
    reliability_class = factor(reliability_class, levels = names(class_cols))
  )

write_csv(
  comp,
  file.path(compare_dir, "MAST_original_vs_removed_male_all_gene_comparison.csv")
)

# ============================================================
# CELLTYPE-LEVEL SIMILARITY METRICS
# ============================================================

similarity_by_celltype <- comp %>%
  group_by(celltype) %>%
  summarise(
    n_genes_compared = n(),

    n_sig_original = sum(sig_original, na.rm = TRUE),
    n_sig_removed = sum(sig_removed, na.rm = TRUE),
    n_sig_both = sum(sig_original & sig_removed, na.rm = TRUE),
    n_sig_both_same_direction = sum(sig_original & sig_removed & same_direction, na.rm = TRUE),
    n_original_only = sum(sig_original & !sig_removed, na.rm = TRUE),
    n_removed_only = sum(!sig_original & sig_removed, na.rm = TRUE),
    n_discordant_sig_both = sum(sig_original & sig_removed & !same_direction, na.rm = TRUE),

    jaccard_DEG = n_sig_both / pmax(sum(sig_original | sig_removed, na.rm = TRUE), 1),
    overlap_fraction_original_retained = n_sig_both / pmax(n_sig_original, 1),
    overlap_fraction_removed_also_original = n_sig_both / pmax(n_sig_removed, 1),

    spearman_log2FC_all = suppressWarnings(cor(avg_log2FC_original, avg_log2FC_removed, method = "spearman", use = "pairwise.complete.obs")),
    pearson_log2FC_all = suppressWarnings(cor(avg_log2FC_original, avg_log2FC_removed, method = "pearson", use = "pairwise.complete.obs")),

    median_abs_delta_log2FC_all = median(abs_delta_log2FC, na.rm = TRUE),
    median_log2FC_retention_original_sig = median(log2FC_retention[sig_original], na.rm = TRUE),

    n_reliable_sig_both = sum(reliability_class == "reliable_sig_both", na.rm = TRUE),
    n_likely_male_sample_driven = sum(reliability_class == "likely_male_sample_driven", na.rm = TRUE),
    n_original_only_male_sensitive = sum(reliability_class == "original_only_male_sensitive", na.rm = TRUE),
    n_removed_only_unmasked = sum(reliability_class == "removed_only_unmasked", na.rm = TRUE),

    .groups = "drop"
  ) %>%
  arrange(desc(n_sig_original), celltype)

write_csv(
  similarity_by_celltype,
  file.path(compare_dir, "MAST_original_vs_removed_male_similarity_by_celltype.csv")
)

class_summary <- comp %>%
  count(celltype, reliability_class, name = "n_genes") %>%
  group_by(celltype) %>%
  mutate(prop_genes = n_genes / sum(n_genes)) %>%
  ungroup()

write_csv(
  class_summary,
  file.path(compare_dir, "MAST_original_vs_removed_male_reliability_class_summary.csv")
)

# ============================================================
# GENE LISTS
# ============================================================

reliable_genes <- comp %>%
  filter(reliability_class == "reliable_sig_both") %>%
  arrange(celltype, mast_FDR_removed, mast_FDR_original, desc(abs_log2FC_removed))

male_sensitive_genes <- comp %>%
  filter(reliability_class %in% c("likely_male_sample_driven", "original_only_male_sensitive")) %>%
  arrange(celltype, desc(male_sensitivity_score), mast_FDR_original)

removed_only_genes <- comp %>%
  filter(reliability_class == "removed_only_unmasked") %>%
  arrange(celltype, mast_FDR_removed, desc(abs_log2FC_removed))

discordant_genes <- comp %>%
  filter(reliability_class == "discordant_sig_both") %>%
  arrange(celltype, desc(abs_delta_log2FC))

write_csv(reliable_genes, file.path(compare_dir, "genes_reliable_sig_in_both_same_direction.csv"))
write_csv(male_sensitive_genes, file.path(compare_dir, "genes_male_sample_sensitive_original_only.csv"))
write_csv(removed_only_genes, file.path(compare_dir, "genes_removed_male_sample_only_unmasked.csv"))
write_csv(discordant_genes, file.path(compare_dir, "genes_discordant_direction_after_removal.csv"))

# ============================================================
# PLOTS
# ============================================================

deg_counts_long <- bind_rows(
  original %>%
    filter(sig) %>%
    mutate(run = "Original") %>%
    count(celltype, run, direction, name = "n_DEG"),
  removed %>%
    filter(sig) %>%
    mutate(run = "Removed male sample") %>%
    count(celltype, run, direction, name = "n_DEG")
) %>%
  mutate(
    signed_n_DEG = ifelse(direction == "Down", -n_DEG, n_DEG),
    celltype = factor(celltype, levels = levels(comp$celltype)),
    run = factor(run, levels = c("Original", "Removed male sample"))
  )

p_deg_counts <- ggplot(deg_counts_long, aes(x = celltype, y = signed_n_DEG, fill = run)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  coord_flip() +
  scale_fill_manual(values = c("Original" = "#808080", "Removed male sample" = "#A53D99")) +
  labs(
    title = "DEG counts before and after removing the male sample",
    subtitle = paste0("Significant = FDR < ", fdr_cutoff, " and |log2FC| >= ", log2fc_cutoff, ". Up positive, down negative."),
    x = NULL,
    y = "Number of DEGs"
  ) +
  theme_mytheme

save_pdf(
  p_deg_counts,
  "plot_01_DEG_counts_original_vs_removed_male.pdf",
  width = 8.27,
  height = max(6, 0.35 * length(unique(deg_counts_long$celltype)) + 2)
)

p_jaccard <- similarity_by_celltype %>%
  mutate(celltype = factor(celltype, levels = levels(comp$celltype))) %>%
  ggplot(aes(x = celltype, y = jaccard_DEG, fill = celltype)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = RNA_pallate, na.value = "grey70", guide = "none") +
  labs(
    title = "DEG overlap similarity by cell type",
    subtitle = "Higher Jaccard means the same genes are called significant before and after removing the male sample.",
    x = NULL,
    y = "Jaccard similarity of significant DEGs"
  ) +
  theme_mytheme

save_pdf(
  p_jaccard,
  "plot_02_Jaccard_similarity_by_celltype.pdf",
  width = 8.27,
  height = max(6, 0.35 * nrow(similarity_by_celltype) + 2)
)

p_lfc_scatter <- comp %>%
  filter(!is.na(avg_log2FC_original), !is.na(avg_log2FC_removed)) %>%
  mutate(plot_alpha = ifelse(reliability_class == "not_sig_either", 0.15, 0.75)) %>%
  ggplot(aes(x = avg_log2FC_original, y = avg_log2FC_removed, color = reliability_class)) +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.35) +
  geom_point(aes(alpha = plot_alpha), size = 0.7) +
  facet_wrap(~celltype, scales = "free") +
  scale_color_manual(values = class_cols, drop = FALSE) +
  scale_alpha_identity() +
  labs(
    title = "log2FC stability after removing the male sample",
    subtitle = "Genes near the dashed line are stable. Original-only genes that collapse after removal are likely male-sample-sensitive.",
    x = "Original MAST avg_log2FC",
    y = "Removed-male MAST avg_log2FC",
    color = "Reliability class"
  ) +
  theme_mytheme

n_facets <- length(unique(comp$celltype))
save_pdf(
  p_lfc_scatter,
  "plot_03_log2FC_scatter_by_celltype.pdf",
  width = 11.69,
  height = max(8.27, ceiling(n_facets / 3) * 2.6)
)

p_fdr_scatter <- comp %>%
  filter(sig_original | sig_removed) %>%
  ggplot(aes(x = neglog10_FDR_original, y = neglog10_FDR_removed, color = reliability_class)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.35) +
  geom_point(alpha = 0.8, size = 1.1) +
  facet_wrap(~celltype, scales = "free") +
  scale_color_manual(values = class_cols, drop = FALSE) +
  labs(
    title = "Significance stability after removing the male sample",
    subtitle = "Points below the dashed line lost significance after removal.",
    x = "-log10(FDR), original",
    y = "-log10(FDR), removed male sample",
    color = "Reliability class"
  ) +
  theme_mytheme

save_pdf(
  p_fdr_scatter,
  "plot_04_FDR_scatter_significant_in_either_run.pdf",
  width = 11.69,
  height = max(8.27, ceiling(n_facets / 3) * 2.6)
)

p_class_heatmap <- class_summary %>%
  mutate(celltype = factor(celltype, levels = levels(comp$celltype))) %>%
  ggplot(aes(x = reliability_class, y = celltype, fill = n_genes)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(low = "white", high = "#A53D99") +
  labs(
    title = "Reliability class counts by cell type",
    x = NULL,
    y = NULL,
    fill = "Genes"
  ) +
  theme_mytheme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

save_pdf(
  p_class_heatmap,
  "plot_05_reliability_class_heatmap.pdf",
  width = 11.69,
  height = max(6, 0.35 * length(unique(class_summary$celltype)) + 2)
)

top_male_sensitive <- male_sensitive_genes %>%
  filter(!is.na(male_sensitivity_score)) %>%
  group_by(celltype) %>%
  slice_max(order_by = male_sensitivity_score, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    gene_label = paste(celltype, gene, sep = " | "),
    gene_label = fct_reorder(gene_label, male_sensitivity_score)
  )

if (nrow(top_male_sensitive) > 0) {
  p_top_male_sensitive <- ggplot(
    top_male_sensitive,
    aes(x = male_sensitivity_score, y = gene_label, color = reliability_class)
  ) +
    geom_segment(aes(x = 0, xend = male_sensitivity_score, yend = gene_label), linewidth = 0.35) +
    geom_point(size = 1.8) +
    scale_color_manual(values = class_cols, drop = FALSE) +
    labs(
      title = "Top male-sample-sensitive genes",
      subtitle = "Higher score means stronger loss of significance and/or effect size after removing the male sample.",
      x = "Male sensitivity score",
      y = NULL,
      color = "Reliability class"
    ) +
    theme_mytheme

  save_pdf(
    p_top_male_sensitive,
    "plot_06_top_male_sample_sensitive_genes.pdf",
    width = 8.27,
    height = max(6, 0.18 * nrow(top_male_sensitive) + 2)
  )
}

p_retention <- comp %>%
  filter(sig_original) %>%
  mutate(celltype = factor(celltype, levels = levels(comp$celltype))) %>%
  ggplot(aes(x = celltype, y = pmin(log2FC_retention, 2), fill = celltype)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.35) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.35) +
  geom_hline(yintercept = 0.5, linetype = "dotted", linewidth = 0.35) +
  coord_flip() +
  scale_fill_manual(values = RNA_pallate, na.value = "grey70", guide = "none") +
  labs(
    title = "Effect-size retention for genes significant in the original run",
    subtitle = "Values near 1 are stable. Values below 0.5 suggest the original effect weakened strongly after male-sample removal.",
    x = NULL,
    y = "Absolute log2FC retention, capped at 2"
  ) +
  theme_mytheme

save_pdf(
  p_retention,
  "plot_07_log2FC_retention_original_significant_genes.pdf",
  width = 8.27,
  height = max(6, 0.35 * length(unique(comp$celltype)) + 2)
)

# ============================================================
# SIMPLE TEXT SUMMARY
# ============================================================

overall_summary <- comp %>%
  summarise(
    n_genes_total = n(),
    n_sig_original = sum(sig_original, na.rm = TRUE),
    n_sig_removed = sum(sig_removed, na.rm = TRUE),
    n_sig_both = sum(sig_original & sig_removed, na.rm = TRUE),
    n_reliable_sig_both = sum(reliability_class == "reliable_sig_both", na.rm = TRUE),
    n_likely_male_sample_driven = sum(reliability_class == "likely_male_sample_driven", na.rm = TRUE),
    n_original_only_male_sensitive = sum(reliability_class == "original_only_male_sensitive", na.rm = TRUE),
    n_removed_only_unmasked = sum(reliability_class == "removed_only_unmasked", na.rm = TRUE),
    n_discordant_sig_both = sum(reliability_class == "discordant_sig_both", na.rm = TRUE),
    jaccard_overall = n_sig_both / pmax(sum(sig_original | sig_removed, na.rm = TRUE), 1),
    original_DEG_retention_fraction = n_sig_both / pmax(n_sig_original, 1)
  )

write_csv(
  overall_summary,
  file.path(compare_dir, "MAST_original_vs_removed_male_overall_summary.csv")
)

print(overall_summary)
print(similarity_by_celltype)

message("\nDone.")
message("Comparison output directory: ", normalizePath(compare_dir))

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggplot2)
})

# ============================================================
# PATHS
# ============================================================

female_deg_dir <- publication_path("DNAP1_SNRNA_ROOT", "Process/MAST_DEG_by_celltype_CDR_batch_removed_male_sample")

age_ex_rds_file <- publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed/All_DEG_Age_Ex.rds")

out_dir <- file.path(
  female_deg_dir,
  "Female_only_AAV_vs_Age_Ex_DEG_directional_overlap"
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# SETTINGS
# ============================================================

female_fdr_cutoff <- 0.05
female_log2fc_cutoff <- 0.25

age_ex_fdr_cutoff <- 0.05
age_ex_log2fc_cutoff <- 0.25

# These assume:
#   log2FoldChange.Age > 0 = higher in aged vs young
#   log2FoldChange.Ex  > 0 = higher in exercise vs sedentary
reverse_age_direction <- FALSE
reverse_exercise_direction <- FALSE

RNA_pallate <- c(
  '#A53D99', '#F99349', '#87CEEB', '#FDDA62', '#808080',
  '#BEBADA', "#A81B1B", "#101369", "#606692", "#3A85A8",
  "#42977E", "#591f0a", "#3b7818", "#B6D66A", "#1CE382",
  "#BF862B", "#CC6A6F", "#E086B5", "#5F1296", "#5F1296",
  "#B784A7"
)

RNA_order <- c(
  "Macrophage", "M1", "M2", "M_Ex", "HSPC", "Dendritic",
  "Cardiomyocyte", "Coronary_EC", "Endocardial_EC", "Lymphatic_EC",
  "Fibroblast", "Mesothelial", "B_cell", "T_cell", "Pericytes",
  "Smooth_Muscle", "Cardiac_Neuronal", "Unknown_1", "Unknown_3",
  "ATAC_unknown1", "NK"
)

names(RNA_pallate) <- RNA_order

relationship_cols <- c(
  "AAV_same_as_ageing" = "#A81B1B",
  "AAV_opposite_to_ageing" = "#101369",
  "AAV_same_as_exercise" = "#A020F0",
  "AAV_opposite_to_exercise" = "#377EB8"
)

theme_mytheme <- theme(
  plot.background = element_blank(),
  plot.title = element_text(size = 15, face = "bold"),
  plot.subtitle = element_text(size = 11),
  axis.text = element_text(size = 10, color = "black"),
  axis.title = element_text(size = 12, color = "black"),
  axis.line = element_line(color = "black", linewidth = 0.4),
  axis.ticks = element_line(color = "black", linewidth = 0.4),
  legend.text = element_text(size = 9),
  legend.title = element_text(size = 10),
  legend.position = "right",
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  strip.text = element_text(size = 9, face = "bold"),
  plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
)

# ============================================================
# FIND FEMALE-ONLY MAST FILE
# ============================================================

female_mast_file <- list.files(
  female_deg_dir,
  pattern = "^MAST_ALL_CELLTYPES_.*_all_genes\\.csv$",
  full.names = TRUE
)

if (length(female_mast_file) == 0) {
  stop("No combined female-only MAST all_genes file found in: ", female_deg_dir)
}

if (length(female_mast_file) > 1) {
  female_mast_file <- female_mast_file[which.max(file.info(female_mast_file)$mtime)]
}

message("Using female-only MAST file:")
message(female_mast_file)

# ============================================================
# READ FEMALE-ONLY AAV DEGs
# ============================================================

female_mast <- fread(female_mast_file) %>%
  as_tibble() %>%
  mutate(
    gene = as.character(gene),
    celltype = as.character(celltype),
    avg_log2FC = as.numeric(avg_log2FC),
    mast_p = as.numeric(mast_p),
    mast_FDR = as.numeric(mast_FDR)
  )

female_universe <- female_mast %>%
  filter(!is.na(gene), gene != "", !is.na(celltype), celltype != "") %>%
  distinct(gene, celltype)

female_aav_deg <- female_mast %>%
  filter(
    !is.na(mast_FDR),
    mast_FDR < female_fdr_cutoff,
    !is.na(avg_log2FC),
    abs(avg_log2FC) >= female_log2fc_cutoff
  ) %>%
  mutate(
    aav_log2FC = avg_log2FC,
    aav_direction = case_when(
      aav_log2FC > 0 ~ "Up",
      aav_log2FC < 0 ~ "Down",
      TRUE ~ "Zero"
    )
  ) %>%
  filter(aav_direction %in% c("Up", "Down")) %>%
  dplyr::select(
    gene,
    celltype,
    aav_log2FC,
    aav_direction,
    mast_p,
    mast_FDR,
    pct_ref,
    pct_test,
    n_cells_ref,
    n_cells_test
  ) %>%
  distinct(gene, celltype, .keep_all = TRUE)

write_csv(
  female_aav_deg,
  file.path(out_dir, "female_only_AAV_significant_DEGs.csv")
)

write_csv(
  female_universe,
  file.path(out_dir, "female_only_AAV_tested_gene_universe.csv")
)

message("Female-only AAV DEGs retained: ", nrow(female_aav_deg))

# ============================================================
# READ AGE / EXERCISE DEG RDS
# ============================================================

Age_Ex_DEG <- readRDS(age_ex_rds_file) %>%
  as_tibble()

required_cols <- c(
  "gene",
  "CellType",
  "log2FoldChange.Ex",
  "pvalue.Ex",
  "padj.Ex",
  "log2FoldChange.Age",
  "pvalue.Age",
  "padj.Age"
)

missing_cols <- setdiff(required_cols, colnames(Age_Ex_DEG))

if (length(missing_cols) > 0) {
  stop(
    "Missing expected columns in Age_Ex_DEG:\n",
    paste(missing_cols, collapse = ", "),
    "\n\nAvailable columns:\n",
    paste(colnames(Age_Ex_DEG), collapse = ", ")
  )
}

age_multiplier <- ifelse(reverse_age_direction, -1, 1)
exercise_multiplier <- ifelse(reverse_exercise_direction, -1, 1)

age_deg <- Age_Ex_DEG %>%
  transmute(
    gene = as.character(gene),
    celltype = as.character(CellType),
    contrast = "Ageing",
    deg_log2FC = as.numeric(log2FoldChange.Age) * age_multiplier,
    deg_pvalue = as.numeric(pvalue.Age),
    deg_FDR = as.numeric(padj.Age)
  ) %>%
  filter(
    !is.na(gene),
    gene != "",
    !is.na(celltype),
    celltype != "",
    !is.na(deg_log2FC),
    !is.na(deg_FDR),
    deg_FDR < age_ex_fdr_cutoff,
    abs(deg_log2FC) >= age_ex_log2fc_cutoff
  ) %>%
  mutate(
    deg_direction = case_when(
      deg_log2FC > 0 ~ "Up",
      deg_log2FC < 0 ~ "Down",
      TRUE ~ "Zero"
    )
  ) %>%
  filter(deg_direction %in% c("Up", "Down")) %>%
  distinct(gene, celltype, contrast, .keep_all = TRUE)

exercise_deg <- Age_Ex_DEG %>%
  transmute(
    gene = as.character(gene),
    celltype = as.character(CellType),
    contrast = "Exercise",
    deg_log2FC = as.numeric(log2FoldChange.Ex) * exercise_multiplier,
    deg_pvalue = as.numeric(pvalue.Ex),
    deg_FDR = as.numeric(padj.Ex)
  ) %>%
  filter(
    !is.na(gene),
    gene != "",
    !is.na(celltype),
    celltype != "",
    !is.na(deg_log2FC),
    !is.na(deg_FDR),
    deg_FDR < age_ex_fdr_cutoff,
    abs(deg_log2FC) >= age_ex_log2fc_cutoff
  ) %>%
  mutate(
    deg_direction = case_when(
      deg_log2FC > 0 ~ "Up",
      deg_log2FC < 0 ~ "Down",
      TRUE ~ "Zero"
    )
  ) %>%
  filter(deg_direction %in% c("Up", "Down")) %>%
  distinct(gene, celltype, contrast, .keep_all = TRUE)

external_deg <- bind_rows(age_deg, exercise_deg)

write_csv(
  external_deg,
  file.path(out_dir, "Age_Ex_DEG_standardised_significant_DEGs.csv")
)

message("Ageing DEGs retained: ", nrow(age_deg))
message("Exercise DEGs retained: ", nrow(exercise_deg))

# ============================================================
# DIRECTIONAL OVERLAP
# ============================================================

overlap_long <- female_aav_deg %>%
  inner_join(
    external_deg,
    by = c("gene", "celltype")
  ) %>%
  mutate(
    same_direction = aav_direction == deg_direction,
    biological_relationship = case_when(
      contrast == "Ageing" & same_direction ~ "AAV_same_as_ageing",
      contrast == "Ageing" & !same_direction ~ "AAV_opposite_to_ageing",
      contrast == "Exercise" & same_direction ~ "AAV_same_as_exercise",
      contrast == "Exercise" & !same_direction ~ "AAV_opposite_to_exercise"
    )
  ) %>%
  arrange(celltype, contrast, biological_relationship, gene)

write_csv(
  overlap_long,
  file.path(out_dir, "female_only_AAV_vs_Age_Ex_directional_overlap_long.csv")
)

overlap_summary <- overlap_long %>%
  group_by(celltype, contrast, biological_relationship, aav_direction, deg_direction) %>%
  summarise(
    n_overlap_genes = n_distinct(gene),
    mean_aav_log2FC = mean(aav_log2FC, na.rm = TRUE),
    mean_age_ex_log2FC = mean(deg_log2FC, na.rm = TRUE),
    genes = paste(sort(unique(gene)), collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(celltype, contrast, desc(n_overlap_genes))

write_csv(
  overlap_summary,
  file.path(out_dir, "female_only_AAV_vs_Age_Ex_directional_overlap_summary.csv")
)

overlap_wide <- overlap_long %>%
  dplyr::select(
    gene,
    celltype,
    aav_log2FC,
    aav_direction,
    mast_FDR,
    contrast,
    deg_log2FC,
    deg_FDR,
    deg_direction,
    biological_relationship
  ) %>%
  pivot_wider(
    names_from = contrast,
    values_from = c(deg_log2FC, deg_FDR, deg_direction, biological_relationship),
    names_sep = "__"
  ) %>%
  mutate(
    age_relationship = biological_relationship__Ageing,
    exercise_relationship = biological_relationship__Exercise,
    combined_class = case_when(
      age_relationship == "AAV_same_as_ageing" &
        exercise_relationship == "AAV_opposite_to_exercise" ~
        "AAV_matches_ageing_and_opposes_exercise",

      age_relationship == "AAV_opposite_to_ageing" &
        exercise_relationship == "AAV_same_as_exercise" ~
        "AAV_opposes_ageing_and_matches_exercise",

      age_relationship == "AAV_same_as_ageing" &
        exercise_relationship == "AAV_same_as_exercise" ~
        "AAV_matches_both_ageing_and_exercise",

      age_relationship == "AAV_opposite_to_ageing" &
        exercise_relationship == "AAV_opposite_to_exercise" ~
        "AAV_opposes_both_ageing_and_exercise",

      !is.na(age_relationship) & is.na(exercise_relationship) ~ age_relationship,
      is.na(age_relationship) & !is.na(exercise_relationship) ~ exercise_relationship,
      TRUE ~ "unclear_or_no_overlap"
    )
  ) %>%
  arrange(celltype, combined_class, gene)

write_csv(
  overlap_wide,
  file.path(out_dir, "female_only_AAV_vs_Age_Ex_directional_overlap_wide.csv")
)

combined_class_summary <- overlap_wide %>%
  count(celltype, combined_class, name = "n_genes") %>%
  group_by(celltype) %>%
  mutate(prop_genes = n_genes / sum(n_genes)) %>%
  ungroup() %>%
  arrange(celltype, desc(n_genes))

write_csv(
  combined_class_summary,
  file.path(out_dir, "female_only_AAV_vs_Age_Ex_combined_class_summary.csv")
)

# ============================================================
# PRIORITY GENE TABLES
# ============================================================

priority_age_like_exercise_opposed <- overlap_wide %>%
  filter(combined_class == "AAV_matches_ageing_and_opposes_exercise") %>%
  mutate(
    priority_score =
      abs(aav_log2FC) +
      abs(deg_log2FC__Ageing) +
      abs(deg_log2FC__Exercise)
  ) %>%
  arrange(celltype, desc(priority_score))

priority_exercise_like_age_opposed <- overlap_wide %>%
  filter(combined_class == "AAV_opposes_ageing_and_matches_exercise") %>%
  mutate(
    priority_score =
      abs(aav_log2FC) +
      abs(deg_log2FC__Ageing) +
      abs(deg_log2FC__Exercise)
  ) %>%
  arrange(celltype, desc(priority_score))

write_csv(
  priority_age_like_exercise_opposed,
  file.path(out_dir, "priority_genes_AAV_matches_ageing_and_opposes_exercise.csv")
)

write_csv(
  priority_exercise_like_age_opposed,
  file.path(out_dir, "priority_genes_AAV_opposes_ageing_and_matches_exercise.csv")
)

# ============================================================
# PLOTS
# ============================================================

celltype_levels <- c(
  intersect(RNA_order, unique(c(female_aav_deg$celltype, external_deg$celltype))),
  setdiff(unique(c(female_aav_deg$celltype, external_deg$celltype)), RNA_order)
)

overlap_summary_plot <- overlap_summary %>%
  mutate(
    celltype = factor(celltype, levels = celltype_levels),
    biological_relationship = factor(biological_relationship, levels = names(relationship_cols))
  )

if (nrow(overlap_summary_plot) > 0) {
  p_overlap <- ggplot(
    overlap_summary_plot,
    aes(x = celltype, y = n_overlap_genes, fill = biological_relationship)
  ) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    coord_flip() +
    facet_wrap(~contrast) +
    scale_fill_manual(values = relationship_cols, drop = FALSE) +
    labs(
      title = "Female-only AAV DEGs overlap with ageing and exercise DEGs",
      subtitle = "Direction-aware overlap by cell type",
      x = NULL,
      y = "Number of overlapping genes",
      fill = "Relationship"
    ) +
    theme_mytheme

  ggsave(
    file.path(out_dir, "plot_01_directional_overlap_counts.pdf"),
    p_overlap,
    width = 11.69,
    height = max(6, 0.35 * length(celltype_levels) + 2),
    units = "in",
    useDingbats = FALSE
  )
}

if (nrow(overlap_long) > 0) {
  scatter_dat <- overlap_long %>%
    mutate(
      celltype = factor(celltype, levels = celltype_levels),
      biological_relationship = factor(biological_relationship, levels = names(relationship_cols))
    )

  p_scatter <- ggplot(
    scatter_dat,
    aes(x = deg_log2FC, y = aav_log2FC, color = biological_relationship)
  ) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_point(alpha = 0.75, size = 1.1) +
    facet_grid(contrast ~ celltype, scales = "free") +
    scale_color_manual(values = relationship_cols, drop = FALSE) +
    labs(
      title = "AAV direction versus ageing and exercise direction",
      x = "Ageing or exercise log2FC",
      y = "Female-only AAV log2FC",
      color = "Relationship"
    ) +
    theme_mytheme +
    theme(
      strip.text.x = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(
    file.path(out_dir, "plot_02_AAV_vs_Age_Ex_log2FC_scatter.pdf"),
    p_scatter,
    width = max(11.69, 2.1 * length(celltype_levels)),
    height = 8.27,
    units = "in",
    useDingbats = FALSE,
    limitsize = FALSE
  )
}

if (nrow(combined_class_summary) > 0) {
  p_class <- combined_class_summary %>%
    mutate(celltype = factor(celltype, levels = celltype_levels)) %>%
    ggplot(aes(x = celltype, y = n_genes, fill = combined_class)) +
    geom_col(width = 0.75) +
    coord_flip() +
    labs(
      title = "Combined ageing and exercise directional classes",
      subtitle = "Most age-like class: AAV_matches_ageing_and_opposes_exercise",
      x = NULL,
      y = "Number of genes",
      fill = "Combined class"
    ) +
    theme_mytheme

  ggsave(
    file.path(out_dir, "plot_03_combined_directional_classes.pdf"),
    p_class,
    width = 11.69,
    height = max(6, 0.35 * length(celltype_levels) + 2),
    units = "in",
    useDingbats = FALSE
  )
}

# ============================================================
# PRINT SUMMARY
# ============================================================

message("\nDone.")
message("Output directory:")
message(normalizePath(out_dir))

message("\nTop biological class if AAV/AP-1 effect is ageing-like:")
print(
  priority_age_like_exercise_opposed %>%
    dplyr::select(
      celltype,
      gene,
      aav_log2FC,
      deg_log2FC__Ageing,
      deg_log2FC__Exercise,
      priority_score
    ) %>%
    head(50)
)

message("\nTop biological class if AAV/AP-1 effect is exercise-like:")
print(
  priority_exercise_like_age_opposed %>%
    dplyr::select(
      celltype,
      gene,
      aav_log2FC,
      deg_log2FC__Ageing,
      deg_log2FC__Exercise,
      priority_score
    ) %>%
    head(50)
)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggplot2)
  library(grid)
})

if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  stop("Please install VennDiagram first: install.packages('VennDiagram')")
}

library(VennDiagram)

# ============================================================
# PATHS
# ============================================================

female_deg_dir <- publication_path("DNAP1_SNRNA_ROOT", "Process/MAST_DEG_by_celltype_CDR_batch_removed_male_sample")
age_ex_rds_file <- publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed/All_DEG_Age_Ex.rds")

out_dir <- file.path(
  female_deg_dir,
  "CM_Fib_Venn_directional_overlap"
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# SETTINGS
# ============================================================

selected_celltypes <- c("Cardiomyocyte", "Fibroblast")

female_fdr_cutoff <- 0.05
female_log2fc_cutoff <- 0.25

age_ex_fdr_cutoff <- 0.05
age_ex_log2fc_cutoff <- 0.25

# These assume:
#   log2FoldChange.Age > 0 = higher in aged vs young
#   log2FoldChange.Ex  > 0 = higher in exercise vs sedentary
reverse_age_direction <- FALSE
reverse_exercise_direction <- FALSE

# ============================================================
# READ FEMALE-ONLY AAV RESULTS
# ============================================================

female_mast_file <- list.files(
  female_deg_dir,
  pattern = "^MAST_ALL_CELLTYPES_.*_all_genes\\.csv$",
  full.names = TRUE
)

if (length(female_mast_file) == 0) {
  stop("No combined female-only MAST file found in: ", female_deg_dir)
}

if (length(female_mast_file) > 1) {
  female_mast_file <- female_mast_file[which.max(file.info(female_mast_file)$mtime)]
}

female_mast <- fread(female_mast_file) %>%
  as_tibble() %>%
  mutate(
    gene = as.character(gene),
    celltype = as.character(celltype),
    avg_log2FC = as.numeric(avg_log2FC),
    mast_FDR = as.numeric(mast_FDR)
  ) %>%
  filter(celltype %in% selected_celltypes)

female_universe <- female_mast %>%
  filter(!is.na(gene), gene != "") %>%
  distinct(celltype, gene)

female_aav_deg <- female_mast %>%
  filter(
    !is.na(mast_FDR),
    mast_FDR < female_fdr_cutoff,
    !is.na(avg_log2FC),
    abs(avg_log2FC) >= female_log2fc_cutoff
  ) %>%
  mutate(
    aav_log2FC = avg_log2FC,
    aav_direction = ifelse(aav_log2FC > 0, "Up", "Down")
  ) %>%
  dplyr::select(celltype, gene, aav_log2FC, aav_direction, mast_FDR) %>%
  distinct(celltype, gene, .keep_all = TRUE)

# ============================================================
# READ AGE / EXERCISE RESULTS
# ============================================================

Age_Ex_DEG <- readRDS(age_ex_rds_file) %>%
  as_tibble()

required_cols <- c(
  "gene", "CellType",
  "log2FoldChange.Ex", "padj.Ex",
  "log2FoldChange.Age", "padj.Age"
)

missing_cols <- setdiff(required_cols, colnames(Age_Ex_DEG))
if (length(missing_cols) > 0) {
  stop(
    "Missing expected columns in Age_Ex_DEG: ",
    paste(missing_cols, collapse = ", ")
  )
}

age_multiplier <- ifelse(reverse_age_direction, -1, 1)
exercise_multiplier <- ifelse(reverse_exercise_direction, -1, 1)

age_deg <- Age_Ex_DEG %>%
  transmute(
    gene = as.character(gene),
    celltype = as.character(CellType),
    contrast = "Ageing",
    deg_log2FC = as.numeric(log2FoldChange.Age) * age_multiplier,
    deg_FDR = as.numeric(padj.Age)
  ) %>%
  filter(
    celltype %in% selected_celltypes,
    !is.na(gene), gene != "",
    !is.na(deg_log2FC), !is.na(deg_FDR),
    deg_FDR < age_ex_fdr_cutoff,
    abs(deg_log2FC) >= age_ex_log2fc_cutoff
  ) %>%
  mutate(
    deg_direction = ifelse(deg_log2FC > 0, "Up", "Down")
  ) %>%
  distinct(celltype, gene, contrast, .keep_all = TRUE)

exercise_deg <- Age_Ex_DEG %>%
  transmute(
    gene = as.character(gene),
    celltype = as.character(CellType),
    contrast = "Exercise",
    deg_log2FC = as.numeric(log2FoldChange.Ex) * exercise_multiplier,
    deg_FDR = as.numeric(padj.Ex)
  ) %>%
  filter(
    celltype %in% selected_celltypes,
    !is.na(gene), gene != "",
    !is.na(deg_log2FC), !is.na(deg_FDR),
    deg_FDR < age_ex_fdr_cutoff,
    abs(deg_log2FC) >= age_ex_log2fc_cutoff
  ) %>%
  mutate(
    deg_direction = ifelse(deg_log2FC > 0, "Up", "Down")
  ) %>%
  distinct(celltype, gene, contrast, .keep_all = TRUE)

external_deg <- bind_rows(age_deg, exercise_deg)

# ============================================================
# OVERLAP STATS
# ============================================================

get_overlap_stats <- function(celltype_use, contrast_use, aav_dir_use, deg_dir_use) {
  universe_genes <- female_universe %>%
    filter(celltype == celltype_use) %>%
    pull(gene) %>%
    unique()

  aav_genes <- female_aav_deg %>%
    filter(celltype == celltype_use, aav_direction == aav_dir_use) %>%
    pull(gene) %>%
    unique()

  deg_genes <- external_deg %>%
    filter(celltype == celltype_use, contrast == contrast_use, deg_direction == deg_dir_use) %>%
    pull(gene) %>%
    unique()

  aav_genes <- intersect(aav_genes, universe_genes)
  deg_genes <- intersect(deg_genes, universe_genes)

  overlap_genes <- intersect(aav_genes, deg_genes)

  a <- length(overlap_genes)
  b <- length(setdiff(aav_genes, deg_genes))
  c <- length(setdiff(deg_genes, aav_genes))
  d <- length(setdiff(universe_genes, union(aav_genes, deg_genes)))

  if (length(universe_genes) == 0 || length(aav_genes) == 0 || length(deg_genes) == 0) {
    return(list(
      stats = tibble(
        celltype = celltype_use,
        contrast = contrast_use,
        aav_direction = aav_dir_use,
        deg_direction = deg_dir_use,
        overlap_type = ifelse(aav_dir_use == deg_dir_use, "same_direction", "opposite_direction"),
        n_universe = length(universe_genes),
        n_aav = length(aav_genes),
        n_deg = length(deg_genes),
        n_overlap = a,
        odds_ratio = NA_real_,
        p_value = NA_real_,
        jaccard = NA_real_,
        overlap_genes = paste(sort(overlap_genes), collapse = ";")
      ),
      aav_genes = aav_genes,
      deg_genes = deg_genes
    ))
  }

  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")

  list(
    stats = tibble(
      celltype = celltype_use,
      contrast = contrast_use,
      aav_direction = aav_dir_use,
      deg_direction = deg_dir_use,
      overlap_type = ifelse(aav_dir_use == deg_dir_use, "same_direction", "opposite_direction"),
      n_universe = length(universe_genes),
      n_aav = length(aav_genes),
      n_deg = length(deg_genes),
      n_overlap = a,
      odds_ratio = unname(ft$estimate),
      p_value = ft$p.value,
      jaccard = a / max(length(union(aav_genes, deg_genes)), 1),
      overlap_genes = paste(sort(overlap_genes), collapse = ";")
    ),
    aav_genes = aav_genes,
    deg_genes = deg_genes
  )
}

jobs <- expand_grid(
  celltype = selected_celltypes,
  contrast = c("Ageing", "Exercise"),
  aav_direction = c("Up", "Down"),
  deg_direction = c("Up", "Down")
)

overlap_res <- pmap(
  jobs,
  ~ get_overlap_stats(..1, ..2, ..3, ..4)
)

stats_tbl <- bind_rows(lapply(overlap_res, `[[`, "stats")) %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(celltype, contrast, overlap_type, aav_direction, deg_direction)

write_csv(
  stats_tbl,
  file.path(out_dir, "directional_overlap_stats.csv")
)

# ============================================================
# VENN DRAWING
# ============================================================

make_venn_grob <- function(aav_genes, deg_genes, title_text, stats_row) {
  grid::grobTree(
    VennDiagram::draw.pairwise.venn(
      area1 = length(aav_genes),
      area2 = length(deg_genes),
      cross.area = length(intersect(aav_genes, deg_genes)),
      category = c("Female-only AAV", "Age/Ex DEG"),
      fill = c("#A53D99", "#3A85A8"),
      alpha = c(0.5, 0.5),
      cex = 1.4,
      cat.cex = 1.1,
      cat.pos = c(-20, 20),
      cat.dist = c(0.05, 0.05),
      scaled = FALSE,
      ind = FALSE
    ),
    textGrob(
      title_text,
      x = 0.5, y = 0.97,
      gp = gpar(fontsize = 13, fontface = "bold")
    ),
    textGrob(
      paste0(
        "Overlap = ", stats_row$n_overlap,
        "    OR = ", ifelse(is.na(stats_row$odds_ratio), "NA", round(stats_row$odds_ratio, 2)),
        "    P = ", ifelse(is.na(stats_row$p_value), "NA", signif(stats_row$p_value, 3)),
        "    FDR = ", ifelse(is.na(stats_row$FDR), "NA", signif(stats_row$FDR, 3)),
        "    Jaccard = ", ifelse(is.na(stats_row$jaccard), "NA", round(stats_row$jaccard, 3))
      ),
      x = 0.5, y = 0.08,
      gp = gpar(fontsize = 10)
    )
  )
}

save_venn_page <- function(celltype_use, contrast_use, mode = c("same", "opposite")) {
  mode <- match.arg(mode)

  if (mode == "same") {
    combos <- tribble(
      ~aav_direction, ~deg_direction, ~label,
      "Up",   "Up",   "Up / Up",
      "Down", "Down", "Down / Down"
    )
    out_file <- file.path(
      out_dir,
      paste0("same_direction_", celltype_use, "_", contrast_use, ".pdf")
    )
  } else {
    combos <- tribble(
      ~aav_direction, ~deg_direction, ~label,
      "Up",   "Down", "AAV Up / DEG Down",
      "Down", "Up",   "AAV Down / DEG Up"
    )
    out_file <- file.path(
      out_dir,
      paste0("opposite_direction_", celltype_use, "_", contrast_use, ".pdf")
    )
  }

  pdf(out_file, width = 10, height = 5)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 2)))

  for (i in seq_len(nrow(combos))) {
    res <- get_overlap_stats(
      celltype_use = celltype_use,
      contrast_use = contrast_use,
      aav_dir_use = combos$aav_direction[i],
      deg_dir_use = combos$deg_direction[i]
    )

    stats_row <- stats_tbl %>%
      filter(
        celltype == celltype_use,
        contrast == contrast_use,
        aav_direction == combos$aav_direction[i],
        deg_direction == combos$deg_direction[i]
      )

    venn_grob <- make_venn_grob(
      aav_genes = res$aav_genes,
      deg_genes = res$deg_genes,
      title_text = paste(celltype_use, "|", contrast_use, "|", combos$label[i]),
      stats_row = stats_row
    )

    print(
      venn_grob,
      vp = viewport(layout.pos.row = 1, layout.pos.col = i)
    )
  }

  dev.off()
}

for (ct in selected_celltypes) {
  for (co in c("Ageing", "Exercise")) {
    save_venn_page(ct, co, mode = "same")
    save_venn_page(ct, co, mode = "opposite")
  }
}

# ============================================================
# SIMPLE SUMMARY TABLES
# ============================================================

same_direction_summary <- stats_tbl %>%
  filter(overlap_type == "same_direction") %>%
  arrange(celltype, contrast, desc(n_overlap))

opposite_direction_summary <- stats_tbl %>%
  filter(overlap_type == "opposite_direction") %>%
  arrange(celltype, contrast, desc(n_overlap))

write_csv(
  same_direction_summary,
  file.path(out_dir, "same_direction_overlap_summary.csv")
)

write_csv(
  opposite_direction_summary,
  file.path(out_dir, "opposite_direction_overlap_summary.csv")
)

print(same_direction_summary)
print(opposite_direction_summary)

message("\nDone.")
message("Output directory: ", normalizePath(out_dir))

# ============================================================
# VENN DRAWING - ggplot version, no VennDiagram dependency
# ============================================================

circle_df <- function(x0, y0, r, n = 300) {
  theta <- seq(0, 2 * pi, length.out = n)
  tibble(
    x = x0 + r * cos(theta),
    y = y0 + r * sin(theta)
  )
}

make_simple_venn_plot <- function(celltype_use, contrast_use, aav_dir_use, deg_dir_use) {
  res <- get_overlap_stats(
    celltype_use = celltype_use,
    contrast_use = contrast_use,
    aav_dir_use = aav_dir_use,
    deg_dir_use = deg_dir_use
  )

  stats_row <- stats_tbl %>%
    filter(
      celltype == celltype_use,
      contrast == contrast_use,
      aav_direction == aav_dir_use,
      deg_direction == deg_dir_use
    ) %>%
    slice(1)

  aav_genes <- res$aav_genes
  deg_genes <- res$deg_genes

  n_aav <- length(aav_genes)
  n_deg <- length(deg_genes)
  n_overlap <- length(intersect(aav_genes, deg_genes))
  n_aav_only <- length(setdiff(aav_genes, deg_genes))
  n_deg_only <- length(setdiff(deg_genes, aav_genes))

  circle1 <- circle_df(-0.75, 0, 1)
  circle2 <- circle_df(0.75, 0, 1)

  subtitle_text <- paste0(
    "OR = ", ifelse(is.na(stats_row$odds_ratio), "NA", round(stats_row$odds_ratio, 2)),
    " | P = ", ifelse(is.na(stats_row$p_value), "NA", signif(stats_row$p_value, 3)),
    " | FDR = ", ifelse(is.na(stats_row$FDR), "NA", signif(stats_row$FDR, 3)),
    " | Jaccard = ", ifelse(is.na(stats_row$jaccard), "NA", round(stats_row$jaccard, 3))
  )

  ggplot() +
    geom_polygon(
      data = circle1,
      aes(x = x, y = y),
      fill = "#A53D99",
      alpha = 0.35,
      color = "black",
      linewidth = 0.5
    ) +
    geom_polygon(
      data = circle2,
      aes(x = x, y = y),
      fill = "#3A85A8",
      alpha = 0.35,
      color = "black",
      linewidth = 0.5
    ) +
    annotate(
      "text",
      x = -1.15,
      y = 0,
      label = n_aav_only,
      size = 7,
      fontface = "bold"
    ) +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = n_overlap,
      size = 7,
      fontface = "bold"
    ) +
    annotate(
      "text",
      x = 1.15,
      y = 0,
      label = n_deg_only,
      size = 7,
      fontface = "bold"
    ) +
    annotate(
      "text",
      x = -0.75,
      y = -1.25,
      label = paste0("Female-only AAV\nn = ", n_aav),
      size = 4
    ) +
    annotate(
      "text",
      x = 0.75,
      y = -1.25,
      label = paste0(contrast_use, " DEG\nn = ", n_deg),
      size = 4
    ) +
    coord_fixed(xlim = c(-2.1, 2.1), ylim = c(-1.6, 1.5), expand = FALSE) +
    labs(
      title = paste0(celltype_use, " | ", contrast_use, " | AAV ", aav_dir_use, " vs DEG ", deg_dir_use),
      subtitle = subtitle_text
    ) +
    theme_void() +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    )
}

safe_file <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

venn_jobs <- expand_grid(
  celltype = selected_celltypes,
  contrast = c("Ageing", "Exercise"),
  aav_direction = c("Up", "Down"),
  deg_direction = c("Up", "Down")
) %>%
  mutate(
    overlap_type = ifelse(aav_direction == deg_direction, "same_direction", "opposite_direction")
  )

for (i in seq_len(nrow(venn_jobs))) {
  ct <- venn_jobs$celltype[i]
  co <- venn_jobs$contrast[i]
  ad <- venn_jobs$aav_direction[i]
  dd <- venn_jobs$deg_direction[i]
  ot <- venn_jobs$overlap_type[i]

  p <- make_simple_venn_plot(
    celltype_use = ct,
    contrast_use = co,
    aav_dir_use = ad,
    deg_dir_use = dd
  )

  out_file <- file.path(
    out_dir,
    paste0(
      "venn_",
      safe_file(ot), "_",
      safe_file(ct), "_",
      safe_file(co), "_AAV_",
      safe_file(ad), "_DEG_",
      safe_file(dd),
      ".pdf"
    )
  )

  ggsave(
    filename = out_file,
    plot = p,
    width = 6.5,
    height = 5,
    units = "in",
    useDingbats = FALSE
  )
}

# ============================================================
# SUMMARY TABLES
# ============================================================

same_direction_summary <- stats_tbl %>%
  filter(overlap_type == "same_direction") %>%
  arrange(celltype, contrast, desc(n_overlap))

opposite_direction_summary <- stats_tbl %>%
  filter(overlap_type == "opposite_direction") %>%
  arrange(celltype, contrast, desc(n_overlap))

write_csv(
  same_direction_summary,
  file.path(out_dir, "same_direction_overlap_summary.csv")
)

write_csv(
  opposite_direction_summary,
  file.path(out_dir, "opposite_direction_overlap_summary.csv")
)

print(same_direction_summary)
print(opposite_direction_summary)

message("\nDone.")
message("Output directory: ", normalizePath(out_dir))
