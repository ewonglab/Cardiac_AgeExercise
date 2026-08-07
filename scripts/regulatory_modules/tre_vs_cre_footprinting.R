#!/usr/bin/env Rscript
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Compare TRE and CRE footprint trajectories across aging and exercise.
# Steps: 1) Load TOBIAS annotations; 2) classify motif architecture; 3) summarise and export the footprint comparison figures.
# =======================================
# Source reformatted from emily_plots/tre_vs_cre_forest_v2.R.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))
# ==============================================================================
# TRE vs CRE footprint: multiple visualization approaches
# ==============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(patchwork)
})

setwd(publication_path("PEAK2PATHWAY_ROOT"))
out_dir <- publication_path("PAPER_OUTPUT_ROOT", "regulatory_modules")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_pub <- theme_bw(base_size = 11) +
    theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey30"),
        axis.title = element_text(size = 10)
    )

col_tre <- "#984EA3"
col_cre <- "#FF7F00"
col_aging <- "#E41A1C"
col_exercise <- "#377EB8"

# ==============================================================================
# Load and prepare
# ==============================================================================
cat("Loading TOBIAS data...\n")
d <- as.data.table(readRDS("Tobias_df_Classified_Full.rds"))

tre_motifs <- grep("^Fos|^Jun", unique(d$motif), value = TRUE)
cre_motifs <- grep("^Creb[0-9]|^Crem_|^Atf[124567]", unique(d$motif), value = TRUE)

d_tre <- d[motif %in% tre_motifs]
d_cre <- d[motif %in% cre_motifs]
d_tre[, motif_class := "TRE (FOS/JUN)"]
d_cre[, motif_class := "CRE (CREB/ATF)"]
d_both <- rbind(d_tre, d_cre)
d_both[, aging_delta := MAS_score - YS_score]
d_both[, exercise_delta := MAE_score - MAS_score]

exclude_ct <- c("AllMacrophage", "MacFib", "HSPC")
d_both <- d_both[!CellType %in% exclude_ct]

ct_labels <- c(
    Fibroblast = "Fibroblast", Cardiomyocyte = "Cardiomyocyte",
    M1 = "M1 Mac.", M2 = "M2 Mac.", Pericytes = "Pericytes",
    Coronary_EC = "Coronary EC", Endocardial_EC = "Endocardial EC",
    Lymphatic_EC = "Lymphatic EC"
)
d_both[, ct_label := ct_labels[CellType]]

stats <- d_both[, .(
    mean_aging = mean(aging_delta),
    se_aging = sd(aging_delta) / sqrt(.N),
    mean_exercise = mean(exercise_delta),
    se_exercise = sd(exercise_delta) / sqrt(.N),
    n = .N
), by = .(CellType, ct_label, motif_class)]

# Wilcoxon p-values
cell_types <- unique(d_both$CellType)
pvals <- rbindlist(lapply(cell_types, function(ct) {
    tre <- d_both[CellType == ct & motif_class == "TRE (FOS/JUN)"]
    cre <- d_both[CellType == ct & motif_class == "CRE (CREB/ATF)"]
    data.table(
        CellType = ct,
        ct_label = ct_labels[ct],
        aging_p = wilcox.test(tre$aging_delta, cre$aging_delta)$p.value,
        exercise_p = wilcox.test(tre$exercise_delta, cre$exercise_delta)$p.value
    )
}))

fmt_p <- function(p) {
    ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
}
pvals[, aging_star := fmt_p(aging_p)]
pvals[, exercise_star := fmt_p(exercise_p)]

# Order by TRE aging effect
tre_order <- stats[motif_class == "TRE (FOS/JUN)"][order(-mean_aging), ct_label]

# ==============================================================================
# STYLE 1: Paired lollipop — TRE and CRE side by side, connected
# ==============================================================================
cat("Style 1: Paired lollipop...\n")

# Long format: aging + exercise in one plot
stats_long <- melt(stats,
    id.vars = c("CellType", "ct_label", "motif_class", "n"),
    measure.vars = list(mean = c("mean_aging", "mean_exercise"),
                        se = c("se_aging", "se_exercise")),
    variable.name = "contrast_idx")
stats_long[, contrast := fifelse(contrast_idx == 1,
    "Aging (Aged sed. - Young sed.)", "Exercise (Aged ex. - Aged sed.)")]
stats_long[, lo := mean - 1.96 * se]
stats_long[, hi := mean + 1.96 * se]

stats_long[, ct_label := factor(ct_label, levels = rev(tre_order))]

# Connector segments between TRE and CRE for each cell type x contrast
wide_for_seg <- dcast(stats_long, ct_label + contrast ~ motif_class, value.var = "mean")
setnames(wide_for_seg, c("TRE (FOS/JUN)", "CRE (CREB/ATF)"), c("tre_mean", "cre_mean"))

p1 <- ggplot(stats_long, aes(x = mean, y = ct_label, fill = motif_class)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_segment(data = wide_for_seg,
                 aes(x = tre_mean, xend = cre_mean, y = ct_label, yend = ct_label),
                 color = "grey60", linewidth = 0.5, inherit.aes = FALSE) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, linewidth = 0.35,
                   position = position_dodge(width = 0)) +
    geom_point(size = 3, shape = 21, stroke = 0.4, color = "grey20") +
    facet_wrap(~contrast, ncol = 2) +
    scale_fill_manual(values = c("TRE (FOS/JUN)" = col_tre,
                                  "CRE (CREB/ATF)" = col_cre), name = "") +
    labs(x = expression(Mean~Delta*" TOBIAS score"),
         y = "",
         title = "TRE vs CRE footprint changes with aging and exercise",
         subtitle = "Points = mean +/- 95% CI | Grey segments connect TRE-CRE pairs") +
    theme_pub +
    theme(legend.position = "bottom",
          panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 10))

ggsave(file.path(out_dir, "fig_tre_cre_paired_lollipop.pdf"), p1, width = 10, height = 5)
cat("  Written: fig_tre_cre_paired_lollipop.pdf\n")

# ==============================================================================
# STYLE 2: Diverging bar — TRE-CRE difference per cell type
# ==============================================================================
cat("Style 2: Diverging bar (TRE - CRE difference)...\n")

# Compute difference: TRE mean - CRE mean per cell type
diff_dt <- dcast(stats, CellType + ct_label ~ motif_class,
                 value.var = c("mean_aging", "mean_exercise"))
setnames(diff_dt,
    c("mean_aging_TRE (FOS/JUN)", "mean_aging_CRE (CREB/ATF)",
      "mean_exercise_TRE (FOS/JUN)", "mean_exercise_CRE (CREB/ATF)"),
    c("tre_aging", "cre_aging", "tre_exer", "cre_exer"))
diff_dt[, aging_diff := tre_aging - cre_aging]
diff_dt[, exercise_diff := tre_exer - cre_exer]

diff_long <- melt(diff_dt, id.vars = c("CellType", "ct_label"),
    measure.vars = c("aging_diff", "exercise_diff"),
    variable.name = "contrast", value.name = "diff")
diff_long[, contrast := fifelse(contrast == "aging_diff",
    "Aging", "Exercise")]
diff_long[, contrast := factor(contrast, levels = c("Aging", "Exercise"))]

# Merge p-values
pvals_long <- melt(pvals, id.vars = c("CellType", "ct_label"),
    measure.vars = c("aging_star", "exercise_star"),
    variable.name = "contrast_v", value.name = "star")
pvals_long[, contrast := fifelse(contrast_v == "aging_star", "Aging", "Exercise")]
pvals_long[, contrast := factor(contrast, levels = c("Aging", "Exercise"))]
diff_long <- merge(diff_long, pvals_long[, .(ct_label, contrast, star)],
                   by = c("ct_label", "contrast"))

diff_long[, ct_label := factor(ct_label, levels = rev(tre_order))]

p2 <- ggplot(diff_long, aes(x = diff, y = ct_label, fill = contrast)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_col(width = 0.6, position = position_dodge(width = 0.7)) +
    geom_text(aes(label = star,
                  x = fifelse(diff >= 0, diff + 0.003, diff - 0.003),
                  hjust = fifelse(diff >= 0, 0, 1)),
              size = 3.5, fontface = "bold", color = "grey30",
              position = position_dodge(width = 0.7)) +
    scale_fill_manual(values = c("Aging" = col_aging, "Exercise" = col_exercise),
                      name = "") +
    labs(x = expression(Delta*" TOBIAS score difference (TRE minus CRE)"),
         y = "",
         title = "TRE-CRE divergence across cell types",
         subtitle = "Positive = TRE gains more (or loses less) than CRE | Wilcoxon TRE vs CRE") +
    theme_pub +
    theme(legend.position = "bottom",
          panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 10))

ggsave(file.path(out_dir, "fig_tre_cre_diverging_bar.pdf"), p2, width = 7, height = 5)
cat("  Written: fig_tre_cre_diverging_bar.pdf\n")

# ==============================================================================
# STYLE 3: Scatter — aging delta vs exercise delta, TRE vs CRE per cell type
# ==============================================================================
cat("Style 3: Scatter (aging vs exercise delta)...\n")

stats[, ct_label_s := factor(ct_label, levels = tre_order)]

p3 <- ggplot(stats, aes(x = mean_aging, y = mean_exercise,
                          fill = motif_class, label = ct_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    # Quadrant shading
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
             fill = "#984EA3", alpha = 0.04) +
    annotate("text", x = 0.03, y = -0.06, label = "Aging gain +\nExercise reversal",
             size = 2.8, color = col_tre, fontface = "italic", hjust = 0) +
    # Connect TRE-CRE pairs with arrows
    geom_segment(data = dcast(stats, ct_label ~ motif_class,
                              value.var = c("mean_aging", "mean_exercise")),
                 aes(x = `mean_aging_CRE (CREB/ATF)`,
                     xend = `mean_aging_TRE (FOS/JUN)`,
                     y = `mean_exercise_CRE (CREB/ATF)`,
                     yend = `mean_exercise_TRE (FOS/JUN)`),
                 arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
                 color = "grey50", linewidth = 0.4, inherit.aes = FALSE) +
    geom_errorbarh(aes(xmin = mean_aging - 1.96 * se_aging,
                       xmax = mean_aging + 1.96 * se_aging),
                   height = 0, linewidth = 0.3, color = "grey50") +
    geom_errorbar(aes(ymin = mean_exercise - 1.96 * se_exercise,
                      ymax = mean_exercise + 1.96 * se_exercise),
                  width = 0, linewidth = 0.3, color = "grey50") +
    geom_point(size = 3.5, shape = 21, stroke = 0.4, color = "grey20") +
    ggrepel::geom_text_repel(size = 2.8, max.overlaps = 20, seed = 42,
                              segment.color = "grey70", segment.size = 0.3,
                              min.segment.length = 0.2) +
    scale_fill_manual(values = c("TRE (FOS/JUN)" = col_tre,
                                  "CRE (CREB/ATF)" = col_cre), name = "") +
    labs(x = expression(Mean~aging~Delta~(Aged~sed.~-~Young~sed.)),
         y = expression(Mean~exercise~Delta~(Aged~ex.~-~Aged~sed.)),
         title = "Aging vs exercise footprint changes: TRE vs CRE",
         subtitle = "Arrows point from CRE to TRE for each cell type | Mean +/- 95% CI") +
    theme_pub +
    theme(legend.position = "bottom")

ggsave(file.path(out_dir, "fig_tre_cre_aging_vs_exercise_scatter.pdf"),
       p3, width = 7, height = 6)
cat("  Written: fig_tre_cre_aging_vs_exercise_scatter.pdf\n")

# ==============================================================================
# STYLE 4: Heatmap — cell type x (motif class x contrast)
# ==============================================================================
cat("Style 4: Heatmap...\n")

heat_dt <- rbind(
    stats[, .(ct_label, motif_class, value = mean_aging,
              panel = "Aging")],
    stats[, .(ct_label, motif_class, value = mean_exercise,
              panel = "Exercise")]
)
heat_dt[, panel := factor(panel, levels = c("Aging", "Exercise"))]
heat_dt[, ct_label := factor(ct_label, levels = rev(tre_order))]

# Add significance
pvals_heat <- rbind(
    pvals[, .(ct_label, panel = "Aging", star = aging_star)],
    pvals[, .(ct_label, panel = "Exercise", star = exercise_star)]
)
pvals_heat[, panel := factor(panel, levels = c("Aging", "Exercise"))]

# Position star at midpoint between TRE and CRE
star_pos <- heat_dt[, .(mid = mean(value)), by = .(ct_label, panel)]
star_pos <- merge(star_pos, pvals_heat, by = c("ct_label", "panel"))

p4 <- ggplot(heat_dt, aes(x = motif_class, y = ct_label, fill = value)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%+.3f", value)),
              size = 3, color = "white", fontface = "bold") +
    facet_wrap(~panel, ncol = 2) +
    scale_fill_gradient2(low = col_exercise, mid = "grey90", high = col_aging,
                         midpoint = 0, name = expression(Mean~Delta~score),
                         guide = guide_colorbar(barwidth = 8, barheight = 0.5,
                                                title.position = "top")) +
    labs(x = "", y = "",
         title = "TRE vs CRE footprint changes: aging and exercise",
         subtitle = "Mean TOBIAS score delta per cell type") +
    theme_pub +
    theme(legend.position = "bottom",
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 9, face = "bold"),
          axis.text.y = element_text(size = 10))

ggsave(file.path(out_dir, "fig_tre_cre_heatmap.pdf"), p4, width = 8, height = 5)
cat("  Written: fig_tre_cre_heatmap.pdf\n")

cat("\nDone. All 4 styles written.\n")
