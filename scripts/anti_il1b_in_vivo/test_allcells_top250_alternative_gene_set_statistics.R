#!/usr/bin/env Rscript
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Test anti-IL-1β responses of all-cell age and exercise top-250 gene sets.
# Steps: 1) Load the edgeR object and gene-set values; 2) run CAMERA and exact permutation tests; 3) export statistics and figure-ready plots.
# =======================================
# Source reformatted from anti-il1code/test_allcells_top250_alternative_gene_set_statistics.R.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(grid)
})

project_dir <- publication_path("H3K27AC_PROJECT_ROOT")
base_dir <- file.path(
  project_dir, "results/antiil1b_age_exercise_agent_explore",
  "edgeR_allcells_CMremoved_topN_boxplots"
)
out_dir <- file.path(
  project_dir, "results/antiil1b_age_exercise_agent_explore",
  "edgeR_allcells_alternative_gene_set_tests"
)
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

anti_file <- publication_path("ANTI_IL1B_RNA_EDGER_OBJECT")
set_file <- file.path(base_dir, "top250_gene_set_antiIL1b_values.csv")

theme_pub <- function(base_size = 6) {
  theme_classic(base_size = base_size, base_family = "ArialMT") +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1.1, margin = margin(b = 1)),
      plot.subtitle = element_text(size = base_size - 0.2, color = "grey25", margin = margin(b = 1)),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.2, color = "black"),
      axis.line = element_line(linewidth = 0.22, color = "black"),
      axis.ticks = element_line(linewidth = 0.18, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.3),
      plot.margin = margin(2, 2, 2, 2)
    )
}

save_plot <- function(p, stub, width, height) {
  qpdf <- function(filename, width, height, ...) {
    quartz(type = "pdf", file = filename, width = width, height = height, family = "ArialMT")
  }
  ggsave(file.path(plot_dir, paste0(stub, ".pdf")), p, width = width, height = height,
         units = "in", device = qpdf, bg = "white")
  ggsave(file.path(plot_dir, paste0(stub, ".png")), p, width = width, height = height,
         units = "in", dpi = 600, bg = "white")
}

fmt_p <- function(p) {
  ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
         formatC(p, format = "f", digits = 3))
}

p_stars <- function(p) {
  fifelse(
    p < 1e-4, "****",
    fifelse(p < 1e-3, "***",
            fifelse(p < 0.01, "**",
                    fifelse(p < 0.05, "*", "ns")))
  )
}

obj <- readRDS(anti_file)
design <- obj$fit$design
contrast_col <- which(colnames(design) == "condition_deOld_anti-IL-1b")
if (length(contrast_col) != 1) stop("Could not identify anti-IL1b contrast column.")

sets <- fread(set_file)[analysis == "All cells"]
sets[, set := factor(
  as.character(set),
  levels = c("Age up", "Age down", "Exercise up", "Exercise down")
)]

gene_map <- unique(as.data.table(obj$de)[, .(gene_id, gene_symbol)], by = "gene_id")
gene_symbols <- gene_map$gene_symbol[match(rownames(obj$y$counts), gene_map$gene_id)]
voom_obj <- voom(obj$y, design, plot = FALSE)

gene_sets <- split(sets$gene, sets$set)
indices <- lapply(gene_sets, function(g) which(gene_symbols %in% unique(g)))

camera_est <- as.data.table(camera(
  voom_obj, indices, design, contrast = contrast_col,
  inter.gene.cor = NA, sort = FALSE
), keep.rownames = "set")
camera_est[, method := "CAMERA, set-specific correlation"]

camera_default <- as.data.table(camera(
  voom_obj, indices, design, contrast = contrast_col,
  inter.gene.cor = 0.01, sort = FALSE
), keep.rownames = "set")
camera_default[, `:=`(
  Correlation = 0.01,
  method = "CAMERA, fixed correlation 0.01"
)]

fry_sets <- as.data.table(fry(
  voom_obj, indices, design, contrast = contrast_col, sort = "none"
), keep.rownames = "set")
fry_sets[, method := "FRY self-contained rotation approximation"]

set_test_results <- rbindlist(list(
  camera_est[, .(set, method, NGenes, Correlation, Direction, PValue, FDR)],
  camera_default[, .(set, method, NGenes, Correlation, Direction, PValue, FDR)],
  fry_sets[, .(
    set, method, NGenes, Correlation = NA_real_, Direction, PValue, FDR,
    PValue_Mixed = PValue.Mixed, FDR_Mixed = FDR.Mixed
  )]
), use.names = TRUE, fill = TRUE)
fwrite(set_test_results, file.path(out_dir, "allcells_top250_alternative_gene_set_tests.csv"))

make_weighted_signature <- function(comp_name) {
  up_genes <- sets[comparison == comp_name & direction == "Up", unique(gene)]
  down_genes <- sets[comparison == comp_name & direction == "Down", unique(gene)]
  up_idx <- which(gene_symbols %in% up_genes)
  down_idx <- which(gene_symbols %in% down_genes)
  data.frame(
    Gene = c(up_idx, down_idx),
    Weight = c(rep(1, length(up_idx)), rep(-1, length(down_idx)))
  )
}

signature_indices <- list(
  Age = make_weighted_signature("Age"),
  Exercise = make_weighted_signature("Exercise")
)

set.seed(20260731)
mroast_signature <- as.data.table(mroast(
  voom_obj, signature_indices, design, contrast = contrast_col,
  nrot = 9999, sort = "none"
), keep.rownames = "signature")
mroast_signature[, method := "MROAST weighted signature, 9999 rotations"]

fry_signature <- as.data.table(fry(
  voom_obj, signature_indices, design, contrast = contrast_col, sort = "none"
), keep.rownames = "signature")
fry_signature[, method := "FRY weighted signature"]

signature_rotation_results <- rbindlist(list(
  mroast_signature,
  fry_signature
), use.names = TRUE, fill = TRUE)
fwrite(
  signature_rotation_results,
  file.path(out_dir, "allcells_age_exercise_weighted_signature_rotation_tests.csv")
)

logcpm <- cpm(obj$y, log = TRUE, prior.count = 2)
z_logcpm <- t(scale(t(logcpm)))
groups <- as.character(obj$y$samples$group)
anti_group <- groups == "Old_anti-IL-1b"

sample_scores <- rbindlist(lapply(names(signature_indices), function(sig) {
  idx <- signature_indices[[sig]]
  up_idx <- idx$Gene[idx$Weight == 1]
  down_idx <- idx$Gene[idx$Weight == -1]
  score <- colMeans(z_logcpm[up_idx, , drop = FALSE]) -
    colMeans(z_logcpm[down_idx, , drop = FALSE])
  data.table(
    signature = sig,
    sample = colnames(logcpm),
    group = fifelse(anti_group, "anti-IL1b", "IgG"),
    score = score
  )
}))

exact_permutation <- function(score, treated) {
  observed <- mean(score[treated]) - mean(score[!treated])
  combinations <- combn(length(score), sum(treated))
  null_delta <- apply(combinations, 2, function(i) {
    mean(score[i]) - mean(score[-i])
  })
  p <- sum(abs(null_delta) >= abs(observed)) / length(null_delta)
  list(observed = observed, p = p, n_permutations = length(null_delta))
}

sample_score_stats <- sample_scores[, {
  treated <- group == "anti-IL1b"
  perm <- exact_permutation(score, treated)
  .(
    n_anti = sum(treated),
    n_IgG = sum(!treated),
    mean_anti = mean(score[treated]),
    mean_IgG = mean(score[!treated]),
    delta_anti_minus_IgG = perm$observed,
    exact_permutation_p_two_sided = perm$p,
    n_exact_permutations = perm$n_permutations,
    welch_t_p_two_sided = t.test(score[treated], score[!treated])$p.value,
    wilcoxon_p_two_sided = wilcox.test(
      score[treated], score[!treated], alternative = "two.sided", exact = TRUE
    )$p.value
  )
}, by = signature]
sample_score_stats[, FDR_exact_permutation := p.adjust(
  exact_permutation_p_two_sided, method = "BH"
)]
fwrite(sample_scores, file.path(out_dir, "allcells_sample_level_signature_scores.csv"))
fwrite(sample_score_stats, file.path(out_dir, "allcells_sample_level_signature_tests.csv"))

camera_plot_stats <- camera_est[, .(
  set = factor(set, levels = levels(sets$set)),
  p = PValue,
  label = paste0("P=", fmt_p(PValue)),
  y = 0.50
)]
fill_values <- c(
  "Age up" = "#C94F52",
  "Age down" = "#3D6F9F",
  "Exercise up" = "#C94F52",
  "Exercise down" = "#3D6F9F"
)

p_box <- ggplot(sets, aes(set, anti_logFC, fill = set)) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.24) +
  geom_boxplot(width = 0.42, outlier.shape = NA, linewidth = 0.28) +
  geom_text(
    data = camera_plot_stats,
    aes(set, y, label = label),
    inherit.aes = FALSE, family = "ArialMT", size = 1.7
  ) +
  scale_fill_manual(values = fill_values, guide = "none") +
  scale_x_discrete(labels = c(
    "Age up" = "Age up", "Age down" = "Age down",
    "Exercise up" = "Ex up", "Exercise down" = "Ex down"
  )) +
  coord_cartesian(ylim = c(-0.58, 0.58), clip = "off") +
  labs(
    title = NULL,
    subtitle = "CAMERA with set-specific gene correlation",
    x = NULL,
    y = "anti-IL1b bulk log2FC"
  ) +
  theme_pub(6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_box, "allcells_top250_boxplot_camera_estimated_correlation", 2.8, 2.5)

camera_default_plot_stats <- camera_default[, .(
  set = factor(set, levels = levels(sets$set)),
  p = PValue,
  label = paste0(p_stars(PValue), "\nP=", fmt_p(PValue)),
  y = 0.50
)]

p_box_default <- ggplot(sets, aes(set, anti_logFC, fill = set)) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.24) +
  geom_boxplot(width = 0.42, outlier.shape = NA, linewidth = 0.28) +
  geom_text(
    data = camera_default_plot_stats,
    aes(set, y, label = label),
    inherit.aes = FALSE, family = "ArialMT", size = 1.7, lineheight = 0.8
  ) +
  scale_fill_manual(values = fill_values, guide = "none") +
  scale_x_discrete(labels = c(
    "Age up" = "Age up", "Age down" = "Age down",
    "Exercise up" = "Ex up", "Exercise down" = "Ex down"
  )) +
  coord_cartesian(ylim = c(-0.58, 0.58), clip = "off") +
  labs(
    title = NULL,
    subtitle = "Standard CAMERA (inter-gene correlation = 0.01)",
    x = NULL,
    y = "anti-IL1b bulk log2FC"
  ) +
  theme_pub(6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_box_default, "allcells_top250_boxplot_camera_standard", 2.8, 2.5)

score_ann <- sample_score_stats[, .(
  label = paste0("Exact P=", fmt_p(exact_permutation_p_two_sided)),
  y = max(sample_scores[signature == .BY$signature, score]) + 0.12
), by = signature]

p_score <- ggplot(sample_scores, aes(group, score, fill = group)) +
  geom_hline(yintercept = 0, color = "grey72", linewidth = 0.2) +
  geom_boxplot(width = 0.42, outlier.shape = NA, linewidth = 0.25, alpha = 0.35) +
  geom_point(
    position = position_jitter(width = 0.07, height = 0),
    size = 1.4, shape = 21, color = "black", stroke = 0.2
  ) +
  geom_text(
    data = score_ann,
    aes(x = 1.5, y = y, label = label),
    inherit.aes = FALSE, family = "ArialMT", size = 1.7
  ) +
  facet_wrap(~signature, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c("IgG" = "grey72", "anti-IL1b" = "#2F7D7E"), guide = "none") +
  labs(
    title = "Sample-level age and exercise signature scores",
    subtitle = "Mean z-score(up genes) minus mean z-score(down genes); exact label permutation test",
    x = NULL,
    y = "Signature score"
  ) +
  theme_pub(6)
save_plot(p_score, "allcells_sample_level_signature_scores", 4.3, 2.25)

summary_lines <- c(
  "# Alternative statistics for all-cell edgeR top-250 gene sets",
  "",
  "The original one-sample Wilcoxon tests use genes as observations and do not account for inter-gene correlation.",
  "CAMERA is a competitive test using the 11 anti-IL1b bulk RNA samples and accounts for inter-gene correlation.",
  "MROAST/FRY are self-contained sample-level rotation tests.",
  "The exact permutation test compares sample-level weighted signature scores across all 462 possible assignments of six of eleven samples to anti-IL1b.",
  "",
  "## Per-set tests",
  capture.output(print(set_test_results)),
  "",
  "## Weighted signature rotation tests",
  capture.output(print(signature_rotation_results)),
  "",
  "## Sample-level signature tests",
  capture.output(print(sample_score_stats))
)
writeLines(summary_lines, file.path(out_dir, "analysis_summary.md"))

message("Wrote alternative gene-set tests to: ", out_dir)
