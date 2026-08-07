#!/usr/bin/env Rscript
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Produce the supplementary promoter H3K27ac response plot across anti-IL-1β RNA-response quintiles.
# Steps: 1) Load the promoter RNA/H3K27ac overlap table; 2) calculate median-scaled H3K27ac summaries and tests by RNA quintile; 3) export source tables and publication-ready PNG/PDF plots.
# =======================================
# Source reformatted from anti-il1code/plot_all_gene_promoter_h3k27ac_rna_quintile_boxplots.R.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

project_dir <- publication_path("H3K27AC_PROJECT_ROOT")
overlap_dir <- file.path(project_dir, "results", "h3k27ac_rna_promoter_overlap")
plot_dir <- file.path(overlap_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(overlap_dir, "promoter_all_RNA_genes_with_H3K27ac_window_annotation.csv")
dt <- fread(input_file)
dt <- dt[h3_peak_n > 0 & is.finite(logFC) & is.finite(h3_mean_logFC)]

global_median_h3_log2fc <- median(dt$h3_mean_logFC)
global_median_h3_fc <- 2^global_median_h3_log2fc
dt[, h3_median_scaled_logFC := h3_mean_logFC - global_median_h3_log2fc]

setorder(dt, logFC)
dt[, rna_percentile := 100 * (.I - 0.5) / .N]
dt[, rna_quintile_num := pmin(5L, floor(rna_percentile / 20) + 1L)]
dt[, rna_quintile := factor(rna_quintile_num, levels = 1:5, labels = paste0("Q", 1:5))]

quintile_summary <- dt[, .(
  genes = .N,
  median_rna_log2fc = median(logFC),
  median_rna_fc = 2^median(logFC),
  median_scaled_h3_log2fc = median(h3_median_scaled_logFC),
  q25_scaled_h3_log2fc = quantile(h3_median_scaled_logFC, 0.25),
  q75_scaled_h3_log2fc = quantile(h3_median_scaled_logFC, 0.75),
  raw_median_h3_log2fc = median(h3_mean_logFC),
  raw_median_h3_fc = 2^median(h3_mean_logFC)
), by = .(rna_quintile_num, rna_quintile)]
setorder(quintile_summary, rna_quintile_num)

quintile_labels <- c("Q1\nbottom 20%\nRNA down", "Q2", "Q3", "Q4", "Q5\ntop 20%\nRNA up")

wilcox_pairs <- rbindlist(lapply(combn(levels(dt$rna_quintile), 2, simplify = FALSE), function(pair) {
  x <- dt[rna_quintile == pair[1], h3_median_scaled_logFC]
  y <- dt[rna_quintile == pair[2], h3_median_scaled_logFC]
  wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
  data.table(
    group_1 = pair[1],
    group_2 = pair[2],
    n_1 = length(x),
    n_2 = length(y),
    median_1 = median(x),
    median_2 = median(y),
    median_difference_group2_minus_group1 = median(y) - median(x),
    wilcoxon_pvalue = wt$p.value
  )
}))
wilcox_pairs[, wilcoxon_fdr := p.adjust(wilcoxon_pvalue, method = "BH")]

q1_q5 <- wilcox_pairs[group_1 == "Q1" & group_2 == "Q5"]
spearman_test <- suppressWarnings(cor.test(dt$logFC, dt$h3_median_scaled_logFC, method = "spearman", exact = FALSE))

stats_table <- data.table(
  analysis = c(
    "Global median promoter H3K27ac log2FC used for scaling",
    "Global median promoter H3K27ac FC used for scaling",
    "All genes Spearman RNA log2FC vs median-scaled promoter H3K27ac log2FC",
    "Bottom RNA quintile vs top RNA quintile Wilcoxon",
    "Bottom RNA quintile median-scaled promoter H3K27ac log2FC",
    "Top RNA quintile median-scaled promoter H3K27ac log2FC",
    "Top minus bottom quintile median-scaled promoter H3K27ac log2FC"
  ),
  value = c(
    sprintf("%.4f", global_median_h3_log2fc),
    sprintf("%.3fx", global_median_h3_fc),
    paste0("rho=", signif(unname(spearman_test$estimate), 3), "; P=", signif(spearman_test$p.value, 3), "; n=", nrow(dt)),
    paste0("P=", signif(q1_q5$wilcoxon_pvalue, 3), "; FDR=", signif(q1_q5$wilcoxon_fdr, 3)),
    sprintf("%.4f", q1_q5$median_1),
    sprintf("%.4f", q1_q5$median_2),
    sprintf("%.4f", q1_q5$median_difference_group2_minus_group1)
  )
)

fwrite(quintile_summary, file.path(overlap_dir, "all_gene_RNA_quintile_promoter_H3K27ac_median_scaled_boxplot_summary.csv"))
fwrite(wilcox_pairs, file.path(overlap_dir, "all_gene_RNA_quintile_promoter_H3K27ac_median_scaled_boxplot_pairwise_wilcoxon.csv"))
fwrite(stats_table, file.path(overlap_dir, "all_gene_RNA_quintile_promoter_H3K27ac_median_scaled_boxplot_stats.csv"))

quintile_cols <- c(
  Q1 = "#2166AC",
  Q2 = "#92C5DE",
  Q3 = "#E6E1DA",
  Q4 = "#F4A582",
  Q5 = "#B2182B"
)

theme_box <- function(base_size = 8.6) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = "grey12"),
      axis.text = element_text(color = "grey18"),
      axis.title = element_text(color = "grey12"),
      axis.text.x = element_text(size = base_size - 1.0, lineheight = 0.9),
      plot.title = element_text(face = "bold", color = "grey10", size = base_size + 2.3),
      plot.subtitle = element_text(color = "grey25", size = base_size - 0.3),
      legend.position = "none",
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(6, 8, 6, 6)
    )
}

annotation_y <- 0.68
annotation_tick <- 0.035
annotation_label <- paste0(
  "Q1 vs Q5 Wilcoxon P=",
  scientific(q1_q5$wilcoxon_pvalue, digits = 2)
)

plot <- ggplot(
  dt,
  aes(
    x = rna_quintile_num,
    y = h3_median_scaled_logFC,
    group = rna_quintile,
    fill = rna_quintile
  )
) +
  geom_hline(yintercept = 0, color = "grey58", linewidth = 0.42) +
  geom_boxplot(
    width = 0.24,
    alpha = 0.86,
    color = "grey18",
    linewidth = 0.36,
    outlier.shape = NA
  ) +
  annotate("segment", x = 1, xend = 5, y = annotation_y, yend = annotation_y, linewidth = 0.36, color = "grey15") +
  annotate("segment", x = 1, xend = 1, y = annotation_y, yend = annotation_y - annotation_tick, linewidth = 0.36, color = "grey15") +
  annotate("segment", x = 5, xend = 5, y = annotation_y, yend = annotation_y - annotation_tick, linewidth = 0.36, color = "grey15") +
  annotate("text", x = 3, y = annotation_y + 0.045, label = annotation_label, size = 2.9, color = "grey15") +
  scale_x_continuous(
    breaks = quintile_summary$rna_quintile_num,
    labels = quintile_labels,
    expand = expansion(mult = c(0.035, 0.035))
  ) +
  scale_y_continuous(
    breaks = seq(-1, 1, 0.2),
    labels = label_number(accuracy = 0.01)
  ) +
  scale_fill_manual(values = quintile_cols) +
  coord_cartesian(ylim = c(-0.68, 0.78), clip = "on") +
  labs(
    title = "Promoter H3K27ac",
    subtitle = paste0(
      "RNA quintiles; centred to global promoter median (",
      sprintf("%.3fx", global_median_h3_fc),
      ")"
    ),
    x = "RNA fold-change quintile",
    y = "Median-scaled H3K27ac log2FC"
  ) +
  theme_box()

png_file <- file.path(plot_dir, "all_gene_promoter_H3K27ac_RNA_quintile_boxplots_H3_median_scaled_wilcoxon.png")
pdf_file <- file.path(plot_dir, "all_gene_promoter_H3K27ac_RNA_quintile_boxplots_H3_median_scaled_wilcoxon.pdf")

ggsave(png_file, plot, width = 4.6, height = 3.35, dpi = 450, bg = "white")
ggsave(pdf_file, plot, width = 4.6, height = 3.35, bg = "white")

cat("Wrote RNA quintile median-scaled promoter H3K27ac boxplot:\n", png_file, "\n", pdf_file, "\n", sep = "")
