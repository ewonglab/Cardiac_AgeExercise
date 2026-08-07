#!/usr/bin/env Rscript
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Compare anti-IL-1β H3K27ac shifts at dynamic and stable AP-1 site classes.
# Steps: 1) Load AP-1 class membership; 2) run directional Wilcoxon comparisons; 3) export the statistics table and forest plot.
# =======================================
# Source reformatted from anti-il1code/plot_pooled_ap1_directional_excluding_mac_subsets_hspc.R.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

project_dir <- publication_path("H3K27AC_PROJECT_ROOT")
out_dir <- file.path(
  project_dir,
  "results",
  "cfos_h3k27ac_edger_tobias",
  "tobias_ap1_forest_adjusted_celltype"
)
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

membership_file <- file.path(out_dir, "celltype_TOBIAS_AP1_classes_H3K27ac_logFC_plot_input.csv")

excluded_cell_types <- c("M1", "M2", "MacFib", "HSPC")
groups <- c("Stable", "Age gain", "Age loss", "Ex gain", "Ex loss")
dynamic_groups <- setdiff(groups, "Stable")
group_labels <- c(
  "Stable" = "Stable",
  "Age gain" = "Age Up",
  "Age loss" = "Age Down",
  "Ex gain" = "Ex Up",
  "Ex loss" = "Ex Down"
)

fdr_stars <- function(x) {
  fifelse(
    is.na(x), "",
    fifelse(x < 0.001, "***",
      fifelse(x < 0.01, "**",
        fifelse(x < 0.05, "*", "")
      )
    )
  )
}

theme_pub <- function(base_size = 7.5) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(size = base_size - 1, color = "grey30"),
      axis.text = element_text(size = base_size - 1, color = "black"),
      axis.title = element_text(size = base_size),
      axis.line = element_line(linewidth = 0.18),
      axis.ticks = element_line(linewidth = 0.16),
      legend.position = "none",
      plot.margin = margin(4, 16, 4, 4)
    )
}

safe_shift <- function(dt, group_name) {
  x <- dt[group == group_name, logFC]
  y <- dt[group == "Stable", logFC]
  wt_two <- suppressWarnings(wilcox.test(
    x,
    y,
    alternative = "two.sided",
    exact = FALSE,
    conf.int = TRUE
  ))
  wt_lower <- suppressWarnings(wilcox.test(
    x,
    y,
    alternative = "less",
    exact = FALSE
  ))
  data.table(
    group = group_name,
    n_group = length(x),
    n_stable = length(y),
    median_group = median(x, na.rm = TRUE),
    median_stable = median(y, na.rm = TRUE),
    median_delta = median(x, na.rm = TRUE) - median(y, na.rm = TRUE),
    effect = unname(wt_two$estimate),
    ci_low = unname(wt_two$conf.int[1]),
    ci_high = unname(wt_two$conf.int[2]),
    two_sided_p = wt_two$p.value,
    one_sided_lower_p = wt_lower$p.value
  )
}

message("Reading cell-type AP-1 class membership")
membership <- fread(membership_file)
pooled <- unique(
  membership[
    !cell_type %in% excluded_cell_types & group %in% groups & is.finite(logFC),
    .(peak_id, group = as.character(group), logFC)
  ]
)
pooled[, group := factor(group, levels = groups)]

message("Running unadjusted directional Wilcoxon tests")
stats <- rbindlist(lapply(dynamic_groups, function(g) safe_shift(pooled, g)))
stats[, `:=`(
  one_sided_lower_FDR = p.adjust(one_sided_lower_p, method = "BH"),
  two_sided_FDR = p.adjust(two_sided_p, method = "BH")
)]
stats[, `:=`(
  stars = fdr_stars(one_sided_lower_FDR),
  group_label = factor(group_labels[group], levels = rev(group_labels[dynamic_groups]))
)]
fwrite(
  stats,
  file.path(out_dir, "pooled_no_M1_M2_MacFib_HSPC_TOBIAS_AP1_vs_stable_unadjusted_directional_lower_stats.csv")
)

counts <- pooled[, .(peak_n = uniqueN(peak_id)), by = group]
counts[, group_label := group_labels[as.character(group)]]
fwrite(counts, file.path(out_dir, "pooled_no_M1_M2_MacFib_HSPC_TOBIAS_AP1_class_peak_counts.csv"))

plot_dt <- copy(stats)
plot_dt[, `:=`(
  significant_lower_FDR_0.05 = one_sided_lower_FDR < 0.05,
  star_x = ci_high + 0.0022
)]
x_lim <- c(
  min(plot_dt$ci_low, na.rm = TRUE) - 0.004,
  max(plot_dt$star_x, na.rm = TRUE) + 0.004
)

p <- ggplot(plot_dt, aes(effect, group_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.25) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high, color = significant_lower_FDR_0.05),
    orientation = "y",
    width = 0.13,
    linewidth = 0.36
  ) +
  geom_point(aes(color = significant_lower_FDR_0.05), size = 2.05, shape = 16) +
  geom_text(aes(x = star_x, label = stars), hjust = 0, size = 3.1, color = "grey15") +
  scale_color_manual(values = c(`TRUE` = "#d95d39", `FALSE` = "#8d8d8d"), guide = "none") +
  coord_cartesian(xlim = x_lim, clip = "off") +
  labs(
    title = "Pooled AP-1 classes versus stable sites",
    subtitle = "Excludes M1, M2, Mac/Fib and HSPC; one-sided Wilcoxon FDR tests lower response",
    x = "H3K27ac log2FC shift vs stable AP-1 sites",
    y = NULL
  ) +
  theme_pub()

ggsave(
  file.path(plot_dir, "pooled_no_M1_M2_MacFib_HSPC_TOBIAS_AP1_vs_stable_unadjusted_directional_lower_forest_stars_grey_nonsig.png"),
  p,
  width = 4.75,
  height = 2.45,
  dpi = 450
)
ggsave(
  file.path(plot_dir, "pooled_no_M1_M2_MacFib_HSPC_TOBIAS_AP1_vs_stable_unadjusted_directional_lower_forest_stars_grey_nonsig.pdf"),
  p,
  width = 4.75,
  height = 2.45,
  useDingbats = FALSE
)

message("Wrote HSPC-excluded plot and stats to ", normalizePath(out_dir))
print(stats[, .(
  group,
  n_group,
  n_stable,
  effect = round(effect, 4),
  two_sided_p = signif(two_sided_p, 3),
  two_sided_FDR = signif(two_sided_FDR, 3),
  one_sided_lower_p = signif(one_sided_lower_p, 3),
  one_sided_lower_FDR = signif(one_sided_lower_FDR, 3),
  stars
)])
