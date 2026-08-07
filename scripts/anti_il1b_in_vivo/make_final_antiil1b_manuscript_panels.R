#!/usr/bin/env Rscript
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Assemble the supplied anti-IL-1β manuscript-panel set from final analysis tables.
# Steps: 1) Load verified RNA, H3K27ac, pathway and AP-1 summaries; 2) render individual panels; 3) export source tables and a review composite.
# =======================================
# Source reformatted from anti-il1code/make_final_antiil1b_manuscript_panels.R.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})

project_dir <- publication_path("H3K27AC_PROJECT_ROOT")
out_dir <- file.path(project_dir, "results", "final_antiil1b_manuscript_figure")
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "source_tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

files <- list(
  rna_object = publication_path("ANTI_IL1B_RNA_EDGER_OBJECT"),
  h3_samples = file.path(
    project_dir, "results", "cfos_h3k27ac_edger_tobias",
    "sample_normalization_used.csv"
  ),
  anti_hallmark = file.path(
    project_dir, "results", "rnaseq_hallmark_gsea_logfc_rank",
    "MSigDB_Hallmark_logFC_GSEA_results.csv"
  ),
  top250_values = file.path(
    project_dir, "results", "antiil1b_age_exercise_agent_explore",
    "edgeR_allcells_CMremoved_topN_boxplots",
    "top250_gene_set_antiIL1b_values.csv"
  ),
  camera = file.path(
    project_dir, "results", "antiil1b_age_exercise_agent_explore",
    "edgeR_allcells_alternative_gene_set_tests",
    "allcells_top250_alternative_gene_set_tests.csv"
  ),
  hallmark_concordance = file.path(
    project_dir, "results", "antiil1b_age_exercise_agent_explore",
    "edgeR_GSEA_all_terms_q020_two_tests",
    "Hallmark_q005_in_at_least_two_tests_direction_blocks.csv"
  ),
  ap1_effects = file.path(
    project_dir, "results", "antiil1b_adjusted_ap1_rna_h3k27ac_comparison",
    "adjusted_AP1_effects_RNA_and_H3K27ac.csv"
  ),
  promoter_overlap = file.path(
    project_dir, "results", "h3k27ac_rna_promoter_overlap",
    "promoter_all_RNA_genes_with_H3K27ac_window_annotation.csv"
  ),
  dynamic_ap1 = file.path(
    project_dir, "results", "cfos_h3k27ac_edger_tobias",
    "tobias_ap1_forest_adjusted_celltype",
    paste0(
      "pooled_no_M1_M2_MacFib_HSPC_TOBIAS_AP1_vs_stable_",
      "unadjusted_directional_lower_stats.csv"
    )
  )
)

missing_files <- unlist(files)[!file.exists(unlist(files))]
if (length(missing_files) > 0) {
  stop("Missing required input(s):\n", paste(missing_files, collapse = "\n"))
}

assert_columns <- function(x, required, source_name) {
  missing <- setdiff(required, names(x))
  if (length(missing) > 0) {
    stop(
      source_name, " is missing required column(s): ",
      paste(missing, collapse = ", ")
    )
  }
}

assert_equal <- function(actual, expected, message) {
  if (!identical(actual, expected)) {
    stop(
      message, "\nExpected: ", paste(expected, collapse = ", "),
      "\nObserved: ", paste(actual, collapse = ", ")
    )
  }
}

sig_stars <- function(x) {
  fcase(
    is.na(x), "",
    x < 0.001, "***",
    x < 0.01, "**",
    x < 0.05, "*",
    default = "ns"
  )
}

theme_manuscript <- function(base_size = 6) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(family = "Arial", color = "black"),
      axis.text = element_text(size = base_size - 0.2, color = "black"),
      axis.title = element_text(size = base_size, color = "black"),
      axis.line = element_line(linewidth = 0.18, color = "black"),
      axis.ticks = element_line(linewidth = 0.16, color = "black"),
      plot.title = element_text(
        size = base_size + 0.2, face = "bold", margin = margin(b = 1)
      ),
      plot.subtitle = element_text(
        size = base_size - 0.5, color = "grey25", margin = margin(b = 1)
      ),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      legend.title = element_text(size = base_size - 0.3),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.width = unit(2.6, "mm"),
      legend.key.height = unit(2.2, "mm"),
      panel.spacing = unit(1.2, "mm"),
      plot.margin = margin(2, 2, 2, 2)
    )
}

save_panel <- function(plot, stub, width, height) {
  pdf_device <- function(filename, width, height, ...) {
    if (capabilities("aqua")) {
      quartz(
        type = "pdf", file = filename, width = width, height = height,
        family = "Arial"
      )
    } else {
      cairo_pdf(filename, width = width, height = height, family = "Arial")
    }
  }
  ggsave(
    file.path(plot_dir, paste0(stub, ".pdf")), plot,
    width = width, height = height, units = "in",
    device = pdf_device, bg = "white"
  )
  ggsave(
    file.path(plot_dir, paste0(stub, ".png")), plot,
    width = width, height = height, units = "in",
    dpi = 600, bg = "white"
  )
}

blue <- "#3D6F9F"
red <- "#C94F52"
orange <- "#E35D35"
grey <- "#8F8F8F"
pvalue_cols <- c("grey90", "#F2CE73", "#E56B4A", "#781C24")

message("Auditing assay sample sizes")
rna_object <- readRDS(files$rna_object)
rna_samples <- as.data.table(rna_object$y$samples, keep.rownames = "sample")
assert_columns(rna_samples, c("group", "sample"), "RNA sample table")
rna_counts <- rna_samples[, .N, by = group]
rna_igg_n <- rna_counts[group == "Old_IgG", N]
rna_anti_n <- rna_counts[group == "Old_anti-IL-1b", N]
if (length(rna_igg_n) != 1 || length(rna_anti_n) != 1) {
  stop("Could not identify exactly one IgG and one anti-IL-1beta RNA group.")
}

h3_samples <- fread(files$h3_samples)
assert_columns(h3_samples, c("sample_id", "cond", "used_in_edgeR"), "H3K27ac sample table")
h3_samples <- h3_samples[used_in_edgeR == TRUE]
h3_counts <- h3_samples[, .N, by = cond]
h3_igg_n <- h3_counts[cond == "HeartOld_IgG", N]
h3_anti_n <- h3_counts[cond == "HeartOld_antiIL1b", N]
if (length(h3_igg_n) != 1 || length(h3_anti_n) != 1) {
  stop("Could not identify exactly one IgG and one anti-IL-1beta H3K27ac group.")
}

sample_audit <- data.table(
  assay = c("Bulk RNA-seq", "H3K27ac CUT&Tag"),
  IgG_n = c(rna_igg_n, h3_igg_n),
  anti_IL1b_n = c(rna_anti_n, h3_anti_n)
)
fwrite(sample_audit, file.path(table_dir, "panel_A_assay_sample_sizes.csv"))

message("Panel A: study design")
panel_a <- ggplot() +
  annotate("segment", x = 0.9, xend = 5.2, y = 2.25, yend = 2.25,
           linewidth = 0.22, color = "grey55") +
  annotate("segment", x = 1.0, xend = 1.0, y = 2.12, yend = 2.38,
           linewidth = 0.22) +
  annotate("segment", x = 5.0, xend = 5.0, y = 2.12, yend = 2.38,
           linewidth = 0.22) +
  annotate("rect", xmin = 1.0, xmax = 5.0, ymin = 2.55, ymax = 2.82,
           fill = "grey83", color = NA) +
  annotate("rect", xmin = 1.0, xmax = 5.0, ymin = 1.65, ymax = 1.92,
           fill = "grey45", color = NA) +
  annotate("text", x = 3.0, y = 2.685, label = "IgG", size = 2.0,
           family = "Arial") +
  annotate("text", x = 3.0, y = 1.785, label = "anti-IL-1β", size = 2.0,
           family = "Arial", color = "white") +
  annotate("text", x = 1.0, y = 3.13, label = "Baseline echo", size = 1.8,
           family = "Arial") +
  annotate("text", x = 5.0, y = 3.13, label = "Endpoint echo", size = 1.8,
           family = "Arial") +
  annotate("text", x = 1.0, y = 1.28, label = "18.5 months", size = 1.7,
           family = "Arial") +
  annotate("text", x = 3.0, y = 1.28, label = "5 weeks", size = 1.7,
           family = "Arial") +
  annotate("text", x = 5.0, y = 1.28, label = "20 months", size = 1.7,
           family = "Arial") +
  annotate("segment", x = 5.35, xend = 5.8, y = 2.25, yend = 2.25,
           arrow = arrow(length = unit(1.2, "mm")), linewidth = 0.20) +
  annotate(
    "text", x = 6.0, y = 2.75,
    label = paste0("Bulk RNA-seq\nIgG n=", rna_igg_n, "; anti-IL-1β n=", rna_anti_n),
    hjust = 0, size = 1.75, lineheight = 0.9, family = "Arial"
  ) +
  annotate(
    "text", x = 6.0, y = 1.78,
    label = paste0(
      "H3K27ac CUT&Tag\nIgG n=", h3_igg_n, "; anti-IL-1β n=", h3_anti_n
    ),
    hjust = 0, size = 1.75, lineheight = 0.9, family = "Arial"
  ) +
  annotate("text", x = 3.0, y = 3.72, label = "In vivo anti-IL-1β",
           size = 2.1, fontface = "bold.italic", family = "Arial") +
  coord_cartesian(xlim = c(0.5, 8.0), ylim = c(0.9, 3.9), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(2, 2, 2, 2))
save_panel(panel_a, "panel_A_study_design", 3.15, 1.65)

message("Panel B: anti-IL-1beta Hallmark GSEA")
hallmark <- fread(files$anti_hallmark)
if (!("FDR" %in% names(hallmark)) && "p.adjust" %in% names(hallmark)) {
  setnames(hallmark, "p.adjust", "FDR")
}
assert_columns(
  hallmark,
  c("ID", "NES", "pvalue", "FDR", "setSize"),
  "Anti-IL-1beta Hallmark GSEA table"
)
hallmark_ids <- c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_ADIPOGENESIS",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT"
)
hallmark_labels <- c(
  HALLMARK_OXIDATIVE_PHOSPHORYLATION = "OxPhos",
  HALLMARK_ADIPOGENESIS = "Adipogenesis",
  HALLMARK_TNFA_SIGNALING_VIA_NFKB = "TNFA/NFKB",
  HALLMARK_IL6_JAK_STAT3_SIGNALING = "IL6/JAK/STAT3",
  HALLMARK_COMPLEMENT = "Complement",
  HALLMARK_INTERFERON_ALPHA_RESPONSE = "IFN-alpha",
  HALLMARK_INTERFERON_GAMMA_RESPONSE = "IFN-gamma",
  HALLMARK_INFLAMMATORY_RESPONSE = "Inflammatory",
  HALLMARK_E2F_TARGETS = "E2F",
  HALLMARK_G2M_CHECKPOINT = "G2M"
)
hallmark_plot <- hallmark[ID %in% hallmark_ids]
assert_equal(
  sort(hallmark_plot$ID), sort(hallmark_ids),
  "The selected Hallmark panel terms are incomplete or duplicated."
)
hallmark_plot[, term := unname(hallmark_labels[ID])]
hallmark_plot[, term := factor(term, levels = hallmark_plot[order(NES), term])]
hallmark_plot[, neg_log10_FDR := pmin(-log10(pmax(FDR, 1e-300)), 10)]
fwrite(
  hallmark_plot[, .(ID, term, setSize, NES, pvalue, FDR)],
  file.path(table_dir, "panel_B_hallmark_GSEA.csv")
)

panel_b <- ggplot(hallmark_plot, aes(NES, term, fill = neg_log10_FDR)) +
  geom_vline(xintercept = 0, color = "grey45", linewidth = 0.17) +
  geom_col(width = 0.42) +
  scale_fill_gradientn(
    colors = pvalue_cols, limits = c(0, 10), breaks = seq(0, 10, 2),
    oob = squish, name = expression(-log[10](q))
  ) +
  scale_x_continuous(
    limits = c(-2.5, 1.85), breaks = c(-2, -1, 0, 1),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_discrete(expand = expansion(add = 0.08)) +
  guides(fill = guide_colorbar(barheight = unit(11, "mm"), barwidth = unit(1.4, "mm"))) +
  labs(x = "anti-IL-1β NES", y = NULL) +
  theme_manuscript(5.5) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "right",
    plot.margin = margin(1, 1, 1, 1)
  )
save_panel(panel_b, "panel_B_antiIL1b_hallmark_GSEA", 2.25, 2.05)
save_panel(panel_b, "panel_B_antiIL1b_hallmark_GSEA_FDR_colored", 2.25, 2.05)

message("Panel C: age and exercise gene-set boxplots")
top250 <- fread(files$top250_values)[analysis == "All cells"]
assert_columns(
  top250,
  c("analysis", "comparison", "direction", "gene", "anti_logFC", "set"),
  "Top-250 age/exercise gene-set table"
)
set_levels <- c("Age up", "Age down", "Exercise up", "Exercise down")
top250[, set := factor(as.character(set), levels = set_levels)]
set_counts <- top250[, .N, by = set]
if (nrow(set_counts) != 4 || any(set_counts$N != 250)) {
  stop("Panel C requires exactly 250 genes in each of four all-cell gene sets.")
}

camera <- fread(files$camera)
assert_columns(
  camera, c("set", "method", "NGenes", "PValue", "FDR"),
  "CAMERA results"
)
camera_plot <- camera[
  method == "CAMERA, fixed correlation 0.01" & set %in% set_levels
]
assert_equal(
  sort(camera_plot$set), sort(set_levels),
  "Panel C is missing one or more fixed-correlation CAMERA tests."
)
camera_plot[, set := factor(set, levels = set_levels)]
camera_plot[, stars := sig_stars(FDR)]
camera_plot[, y := 0.50]

fwrite(
  top250[, .(
    set, gene, source_logFC, source_PValue, source_FDR,
    anti_logFC, anti_PValue, anti_FDR
  )],
  file.path(table_dir, "panel_C_top250_gene_values.csv")
)
fwrite(
  camera_plot[, .(set, NGenes, Correlation, Direction, PValue, FDR, stars)],
  file.path(table_dir, "panel_C_CAMERA_tests.csv")
)

panel_c <- ggplot(top250, aes(set, anti_logFC, fill = set)) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.20) +
  geom_boxplot(width = 0.40, outlier.shape = NA, linewidth = 0.24) +
  geom_text(
    data = camera_plot, aes(set, y, label = stars),
    inherit.aes = FALSE, family = "Arial", size = 2.0
  ) +
  scale_fill_manual(
    values = c(
      "Age up" = red, "Age down" = blue,
      "Exercise up" = red, "Exercise down" = blue
    ),
    guide = "none"
  ) +
  scale_x_discrete(labels = c(
    "Age up" = "Age up", "Age down" = "Age down",
    "Exercise up" = "Ex up", "Exercise down" = "Ex down"
  )) +
  scale_y_continuous(
    breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  coord_cartesian(ylim = c(-0.58, 0.58), clip = "on") +
  labs(x = NULL, y = "anti-IL-1β bulk log2FC") +
  theme_manuscript(5.5) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(1, 1, 1, 1)
  )
save_panel(panel_c, "panel_C_age_exercise_top250_CAMERA", 2.0, 2.05)

message("Panel D: Hallmark concordance across anti-IL-1beta, age and exercise")
concordance <- fread(files$hallmark_concordance)
assert_columns(
  concordance,
  c(
    "term", "term_display", "direction_group", "n_tests_q005",
    "NES_anti_IL1b", "NES_edgeR_age", "NES_edgeR_exercise",
    "q_family_anti_IL1b", "q_family_edgeR_age",
    "q_family_edgeR_exercise", "display_order"
  ),
  "Hallmark concordance table"
)
desired_group <- "Anti-IL1b and exercise concordant;\nage NES opposite"
concordance <- concordance[direction_group == desired_group]
if (nrow(concordance) != 8 || any(concordance$n_tests_q005 < 2)) {
  stop(
    "Panel D must contain the eight desired-direction Hallmarks with ",
    "Hallmark-family FDR < 0.05 in at least two comparisons."
  )
}
concordance_long <- melt(
  concordance,
  id.vars = c("term", "term_display", "display_order"),
  measure.vars = patterns(
    NES = "^NES_",
    FDR = "^q_family_"
  ),
  variable.name = "comparison_index"
)
comparison_levels <- c("anti_IL1b", "edgeR_age", "edgeR_exercise")
concordance_long[, comparison := factor(
  comparison_levels[comparison_index],
  levels = comparison_levels,
  labels = c("anti-IL-1β", "age", "exercise")
)]
concordance_long[, term_display := factor(
  term_display,
  levels = rev(concordance[order(display_order), term_display])
)]
concordance_long[, FDR_border := fcase(
  FDR < 0.05, "q < 0.05",
  FDR < 0.10, "q < 0.10",
  FDR < 0.20, "q < 0.20",
  default = NA_character_
)]
concordance_long[, FDR_border := factor(
  FDR_border,
  levels = c("q < 0.05", "q < 0.10", "q < 0.20")
)]
fdr_legend <- concordance_long[rep(1L, 3L)]
fdr_legend[, FDR_border := factor(
  c("q < 0.05", "q < 0.10", "q < 0.20"),
  levels = c("q < 0.05", "q < 0.10", "q < 0.20")
)]
concordance_long[, NES_clip := pmax(pmin(NES, 2.5), -2.5)]
fwrite(
  concordance_long[, .(
    term, term_display, display_order, comparison, NES, FDR, FDR_border
  )],
  file.path(table_dir, "panel_D_hallmark_concordance.csv")
)

panel_d <- ggplot(concordance_long, aes(comparison, term_display)) +
  geom_tile(
    aes(fill = NES_clip), color = "grey88", linewidth = 0.12,
    width = 0.90, height = 0.90
  ) +
  geom_tile(
    data = concordance_long[!is.na(FDR_border)],
    aes(color = FDR_border), fill = NA, linewidth = 0.38,
    width = 0.90, height = 0.90, show.legend = FALSE
  ) +
  geom_point(
    data = fdr_legend,
    aes(color = FDR_border),
    shape = 22, fill = NA, size = 0, stroke = 0,
    show.legend = TRUE
  ) +
  geom_text(
    aes(label = sprintf("%.2f", NES)),
    family = "Arial", size = 1.45
  ) +
  scale_fill_gradient2(
    low = "#315AA6", mid = "white", high = "#C93F42",
    midpoint = 0, limits = c(-2.5, 2.5),
    breaks = c(-2, -1, 0, 1, 2), name = "NES"
  ) +
  scale_color_manual(
    values = c(
      "q < 0.05" = "black",
      "q < 0.10" = "grey35",
      "q < 0.20" = "grey65"
    ),
    breaks = c("q < 0.05", "q < 0.10", "q < 0.20"),
    name = "FDR",
    drop = FALSE
  ) +
  guides(
    fill = guide_colorbar(
      order = 1, barheight = unit(10, "mm"), barwidth = unit(1.5, "mm")
    ),
    color = guide_legend(
      order = 2,
      override.aes = list(
        shape = 22,
        size = 3,
        stroke = 0.45,
        fill = NA,
        color = c("black", "grey35", "grey65")
      )
    )
  ) +
  labs(
    title = "anti-IL-1β aligns with exercise;\nage trends opposite",
    x = NULL, y = NULL
  ) +
  theme_manuscript(5.2) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    plot.margin = margin(1, 1, 1, 1)
  )
save_panel(panel_d, "panel_D_hallmark_age_exercise_concordance", 2.8, 2.25)

message("Panel E: adjusted AP-1 effects in RNA and H3K27ac")
ap1_effects <- fread(files$ap1_effects)
assert_columns(
  ap1_effects,
  c(
    "cell_type", "n_total", "n_ap1", "n_background", "beta", "se",
    "ci_low", "ci_high", "pvalue", "FDR", "assay", "model"
  ),
  "Adjusted AP-1 effect table"
)
expected_cells <- c(
  "AllMacrophage", "Fibroblast", "Cardiomyocyte", "Coronary_EC",
  "Lymphatic_EC", "Pericytes", "Endocardial_EC"
)
expected_assays <- c("RNA\n(all expressed genes)", "H3K27ac\n(all tested peaks)")
ap1_effects <- ap1_effects[
  cell_type %in% expected_cells & assay %in% expected_assays
]
if (nrow(ap1_effects) != length(expected_cells) * length(expected_assays)) {
  stop("Panel E requires one RNA and one H3K27ac estimate for each cell type.")
}
cell_labels <- c(
  AllMacrophage = "Macrophage",
  Fibroblast = "Fibroblast",
  Cardiomyocyte = "Cardiomyocyte",
  Coronary_EC = "Coronary EC",
  Lymphatic_EC = "Lymphatic EC",
  Pericytes = "Pericytes",
  Endocardial_EC = "Endocardial EC"
)
rna_order <- ap1_effects[assay == expected_assays[1]][
  order(-abs(beta)), cell_type
]
ap1_effects[, cell_display := factor(
  unname(cell_labels[cell_type]),
  levels = unname(cell_labels[rna_order])
)]
ap1_effects[, assay := factor(assay, levels = expected_assays)]
ap1_effects[, stars := fifelse(FDR < 0.05, sig_stars(FDR), "")]
ap1_effects[, star_y := ci_low - 0.07 * diff(range(c(ci_low, ci_high, 0))),
            by = assay]
fwrite(
  ap1_effects[, .(
    assay, cell_type, cell_display, n_total, n_ap1, n_background,
    beta, se, ci_low, ci_high, pvalue, FDR, model
  )],
  file.path(table_dir, "panel_E_adjusted_AP1_effects.csv")
)

panel_e <- ggplot(ap1_effects, aes(cell_display, beta)) +
  geom_hline(yintercept = 0, color = "grey45", linewidth = 0.18) +
  geom_col(width = 0.64, fill = "white", color = "grey20", linewidth = 0.18) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.15, linewidth = 0.20
  ) +
  geom_text(
    aes(y = star_y, label = stars),
    family = "Arial", size = 1.65, vjust = 1
  ) +
  facet_wrap(~assay, nrow = 1, scales = "free_y") +
  scale_y_continuous(
    labels = label_number(accuracy = 0.001),
    expand = expansion(mult = c(0.18, 0.04))
  ) +
  labs(
    x = "TOBIAS AP-1 cell-type link",
    y = "Adjusted AP-1 effect on\nanti-IL-1β vs IgG log2FC"
  ) +
  theme_manuscript(5.2) +
  theme(
    axis.text.x = element_text(angle = 42, hjust = 1),
    plot.margin = margin(1, 1, 1, 1)
  )
save_panel(panel_e, "panel_E_adjusted_AP1_effect_RNA_H3K27ac", 4.6, 2.05)

message("Panel F: promoter H3K27ac across RNA response quintiles")
promoter <- fread(files$promoter_overlap)
assert_columns(
  promoter, c("h3_peak_n", "logFC", "h3_mean_logFC"),
  "Promoter RNA-H3K27ac overlap table"
)
promoter <- promoter[
  h3_peak_n > 0 & is.finite(logFC) & is.finite(h3_mean_logFC)
]
global_h3_median <- median(promoter$h3_mean_logFC)
promoter[, h3_median_scaled_logFC := h3_mean_logFC - global_h3_median]
setorder(promoter, logFC)
promoter[, rna_quintile_num := pmin(5L, floor((.I - 1L) * 5 / .N) + 1L)]
promoter[, rna_quintile := factor(
  rna_quintile_num, levels = 1:5, labels = paste0("Q", 1:5)
)]
quintile_summary <- promoter[, .(
  genes = .N,
  min_RNA_log2FC = min(logFC),
  max_RNA_log2FC = max(logFC),
  median_RNA_log2FC = median(logFC),
  median_scaled_H3K27ac_log2FC = median(h3_median_scaled_logFC),
  H3K27ac_Q1 = quantile(h3_median_scaled_logFC, 0.25),
  H3K27ac_Q3 = quantile(h3_median_scaled_logFC, 0.75)
), by = .(rna_quintile_num, rna_quintile)]

quintile_tests <- rbindlist(lapply(
  combn(levels(promoter$rna_quintile), 2, simplify = FALSE),
  function(pair) {
    wt <- suppressWarnings(wilcox.test(
      promoter[rna_quintile == pair[1], h3_median_scaled_logFC],
      promoter[rna_quintile == pair[2], h3_median_scaled_logFC],
      alternative = "two.sided", exact = FALSE
    ))
    data.table(group_1 = pair[1], group_2 = pair[2], pvalue = wt$p.value)
  }
))
quintile_tests[, FDR := p.adjust(pvalue, method = "BH")]
q1_q5 <- quintile_tests[group_1 == "Q1" & group_2 == "Q5"]
if (nrow(q1_q5) != 1) stop("Could not identify the Q1 versus Q5 test.")
fwrite(quintile_summary, file.path(table_dir, "panel_F_RNA_quintile_summary.csv"))
fwrite(quintile_tests, file.path(table_dir, "panel_F_RNA_quintile_tests.csv"))

panel_f <- ggplot(
  promoter,
  aes(rna_quintile_num, h3_median_scaled_logFC, group = rna_quintile,
      fill = rna_quintile)
) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.18) +
  geom_boxplot(
    width = 0.42, outlier.shape = NA, linewidth = 0.22, alpha = 0.90
  ) +
  annotate("segment", x = 1, xend = 5, y = 0.68, yend = 0.68,
           linewidth = 0.20) +
  annotate("segment", x = 1, xend = 1, y = 0.68, yend = 0.64,
           linewidth = 0.20) +
  annotate("segment", x = 5, xend = 5, y = 0.68, yend = 0.64,
           linewidth = 0.20) +
  annotate(
    "text", x = 3, y = 0.73, label = sig_stars(q1_q5$FDR),
    family = "Arial", size = 1.9
  ) +
  scale_fill_manual(
    values = c(
      Q1 = "#2166AC", Q2 = "#92C5DE", Q3 = "#E6E1DA",
      Q4 = "#F4A582", Q5 = "#B2182B"
    ),
    guide = "none"
  ) +
  scale_x_continuous(breaks = 1:5, labels = paste0("Q", 1:5)) +
  scale_y_continuous(
    breaks = seq(-0.6, 0.6, 0.2),
    labels = label_number(accuracy = 0.01)
  ) +
  coord_cartesian(ylim = c(-0.68, 0.78), clip = "on") +
  labs(
    x = "RNA fold-change quintile",
    y = "H3K27ac log2FC\n(median-scaled)"
  ) +
  theme_manuscript(5.2) +
  theme(plot.margin = margin(1, 1, 1, 1))
save_panel(panel_f, "panel_F_promoter_H3K27ac_by_RNA_quintile", 2.05, 2.05)

message("Panel G: dynamic versus stable AP-1 sites")
dynamic <- fread(files$dynamic_ap1)
assert_columns(
  dynamic,
  c(
    "group", "n_group", "n_stable", "effect", "ci_low", "ci_high",
    "two_sided_p", "two_sided_FDR", "group_label"
  ),
  "Dynamic-versus-stable AP-1 table"
)
dynamic <- dynamic[group %in% c("Age gain", "Age loss", "Ex gain", "Ex loss")]
if (nrow(dynamic) != 4) stop("Panel G requires exactly four dynamic AP-1 classes.")
dynamic[, group_label := factor(
  group_label,
  levels = rev(c("Age Up", "Age Down", "Ex Up", "Ex Down"))
)]
dynamic[, significant := two_sided_FDR < 0.05]
dynamic[, stars := fifelse(significant, sig_stars(two_sided_FDR), "")]
dynamic[, star_x := ci_high + 0.0015]
fwrite(
  dynamic[, .(
    group, group_label, n_group, n_stable, median_group, median_stable,
    effect, ci_low, ci_high, two_sided_p, two_sided_FDR, stars
  )],
  file.path(table_dir, "panel_G_dynamic_vs_stable_AP1_two_sided.csv")
)

panel_g <- ggplot(dynamic, aes(effect, group_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55",
             linewidth = 0.18) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high, color = significant),
    orientation = "y", width = 0.12, linewidth = 0.25
  ) +
  geom_point(aes(color = significant), size = 1.15) +
  geom_text(
    aes(x = star_x, label = stars),
    hjust = 0, family = "Arial", size = 1.75
  ) +
  scale_color_manual(values = c(`TRUE` = orange, `FALSE` = grey), guide = "none") +
  scale_x_continuous(
    limits = c(-0.024, 0.016), breaks = c(-0.02, -0.01, 0, 0.01),
    labels = label_number(accuracy = 0.01),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "H3K27ac log2FC shift\nDynamic vs stable AP-1 sites",
    y = NULL
  ) +
  theme_manuscript(5.2) +
  theme(plot.margin = margin(1, 1, 1, 1))
save_panel(panel_g, "panel_G_dynamic_vs_stable_AP1_two_sided_FDR", 2.15, 2.05)

message("Saving compact review composite")
review <- (
  panel_a + panel_b + panel_c + panel_d +
    plot_layout(widths = c(1.35, 0.95, 0.85, 1.20))
) / (
  panel_e + panel_f + panel_g +
    plot_layout(widths = c(1.75, 0.85, 0.90))
) +
  plot_layout(heights = c(1, 1))
save_panel(review, "final_antiIL1b_manuscript_panels_review", 11.0, 4.8)

audit <- data.table(
  item = c(
    "RNA cohort",
    "H3K27ac cohort",
    "Panel C significance annotation",
    "Panel D inclusion rule",
    "Panel E coefficient tests",
    "Panel F Q1-Q5 test",
    "Panel G significance annotation"
  ),
  verified_definition = c(
    paste0("IgG n=", rna_igg_n, "; anti-IL-1beta n=", rna_anti_n),
    paste0("IgG n=", h3_igg_n, "; anti-IL-1beta n=", h3_anti_n),
    "BH FDR from standard CAMERA with fixed inter-gene correlation 0.01",
    paste0(
      "Hallmark-family FDR < 0.05 in at least two comparisons, followed by ",
      "anti-IL-1beta/exercise same sign and age opposite sign"
    ),
    "Two-sided linear-model coefficient tests; BH FDR across seven displayed cell types within assay",
    "Two-sided Wilcoxon rank-sum; BH FDR across all ten quintile pairs",
    "Two-sided Wilcoxon rank-sum versus stable sites; BH FDR across four dynamic classes"
  )
)
fwrite(audit, file.path(out_dir, "final_figure_analysis_audit.csv"))

writeLines(
  c(
    "Final anti-IL-1beta manuscript panels completed.",
    paste0(
      "Run: Rscript ",
      file.path(project_dir, "analysis", "make_final_antiil1b_manuscript_panels.R")
    ),
    paste0("Output directory: ", out_dir),
    paste0("Standalone PDF/PNG panels: ", plot_dir),
    paste0("Exact plotted source tables: ", table_dir),
    paste0("Results text: ", file.path(out_dir, "manuscript_results.md")),
    paste0("Methods text: ", file.path(out_dir, "manuscript_methods.md")),
    paste0("Figure legend: ", file.path(out_dir, "figure_legend.md")),
    paste0("RNA samples: IgG n=", rna_igg_n, ", anti-IL-1beta n=", rna_anti_n),
    paste0("H3K27ac samples: IgG n=", h3_igg_n, ", anti-IL-1beta n=", h3_anti_n),
    paste0(
      "Panel G uses two-sided BH-adjusted tests. Ex-loss AP-1 is nominally ",
      "significant but does not meet FDR < 0.05."
    )
  ),
  file.path(out_dir, "README.txt")
)

message("Done: ", out_dir)
