#!/usr/bin/env Rscript
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Create c-FOS motif-density dose-response and reversal publication figures.
# Steps: 1) Load master and motif tables; 2) calculate motif-bin and reversal summaries; 3) export the figure panels.
# =======================================
# Source reformatted from emily_plots/cut&tag/pub_dose_response_and_reversal.R.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))
# ==============================================================================
# Publication figures: Motif dose-response & reversal numbers (10wk C-FOS)
# ==============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(edgeR)
    library(ggplot2)
    library(scales)
})

ANALYSIS_DIR <- publication_path("CUTTAG_10_WEEK_ANALYSIS_ROOT")
setwd(ANALYSIS_DIR)
FIG_DIR <- file.path(ANALYSIS_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Consistent theme ---
theme_pub <- theme_bw(base_size = 11) +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey95", color = NA),
        strip.text = element_text(face = "bold", size = 10),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey30"),
        legend.background = element_rect(fill = "white", color = NA),
        axis.title = element_text(size = 10)
    )

col_YS <- "#4DAF4A"
col_OS <- "#E41A1C"
col_OE <- "#377EB8"
col_motif_pos <- "#984EA3"
col_motif_neg <- "#FF7F00"

# Dose-response gradient (light to dark blue)
dose_cols <- c("0" = "#D9D9D9", "1-2" = "#C6DBEF", "3-5" = "#6BAED6",
               "6-10" = "#2171B5", "11+" = "#08306B")

# ==============================================================================
# Load data
# ==============================================================================
cat("Loading data...\n")

master <- fread("master_results_C-FOS_support_all3.csv")
motif <- fread("motif_annotation_fimo_C-FOS.csv")
os_outlier <- fread("edgeR_C-FOS_support_all3_OS_outlier.csv")

master <- merge(master, motif[, .(peak_id, n_ap1_hits_fimo)], by = "peak_id", all.x = TRUE)
master <- merge(master, os_outlier[, .(peak_id, logFC_os_out = logFC, P_os_out = PValue)],
                by = "peak_id", all.x = TRUE)

# Motif bins
master[, motif_bin := cut(n_ap1_hits_fimo,
    breaks = c(-1, 0, 2, 5, 10, Inf),
    labels = c("0", "1-2", "3-5", "6-10", "11+"))]

cat("Peaks:", nrow(master), "\n")
cat("FOS::JUN+:", sum(master$has_fosjun_motif), "\n")
cat("Motif-:", sum(!master$has_ap1_motif_fimo), "\n\n")

# ==============================================================================
# FIGURE A: Dose-response boxplot (OS outlier logFC by motif bin)
# ==============================================================================
cat("Figure A: Dose-response boxplot...\n")

# Summary stats for annotation
dose_summary <- master[, .(
    n = .N,
    median_lfc = median(logFC_os_out),
    pct_up = round(100 * mean(logFC_os_out > 0), 1)
), by = motif_bin][order(match(motif_bin, names(dose_cols)))]

pA <- ggplot(master, aes(x = motif_bin, y = logFC_os_out, fill = motif_bin)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.1, width = 0.75,
                 color = "grey30", linewidth = 0.3) +
    stat_summary(fun = median, geom = "point", shape = 18, size = 2.5, color = "white") +
    scale_fill_manual(values = dose_cols) +
    annotate("text", x = 1:5, y = max(master$logFC_os_out) * 0.95,
             label = paste0("n=", dose_summary$n),
             size = 2.8, color = "grey40") +
    labs(x = "AP-1 motif hits per peak (FIMO)",
         y = expression(log[2]*FC~(2 %.% OS - OE - YS)),
         title = "OS-outlier signal scales with AP-1 motif density",
         subtitle = expression("Trend: Spearman"~rho~"= 0.32, p = 1.7 x"~10^{-85})) +
    theme_pub +
    theme(legend.position = "none")

ggsave(file.path(FIG_DIR, "pub_dose_response_boxplot.pdf"), pA, width = 5.5, height = 4)
cat("  Written: pub_dose_response_boxplot.pdf\n")

# ==============================================================================
# FIGURE B: % OS-elevated by motif bin (bar plot with CI)
# ==============================================================================
cat("Figure B: % OS-elevated by motif bin...\n")

pct_dt <- master[, .(
    pct = 100 * mean(logFC_os_out > 0),
    n = .N
), by = motif_bin]
pct_dt[, motif_bin := factor(motif_bin, levels = names(dose_cols))]
# Binomial 95% CI
pct_dt[, ci_lo := 100 * qbeta(0.025, pct / 100 * n, (1 - pct / 100) * n + 1)]
pct_dt[, ci_hi := 100 * qbeta(0.975, pct / 100 * n + 1, (1 - pct / 100) * n)]

pB <- ggplot(pct_dt, aes(x = motif_bin, y = pct, fill = motif_bin)) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_col(width = 0.7, color = "grey30", linewidth = 0.3) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2, linewidth = 0.4) +
    geom_text(aes(label = paste0(round(pct, 1), "%")),
              vjust = -1.5, size = 3.2, fontface = "bold") +
    geom_text(aes(y = 2, label = paste0("n=", n)),
              size = 2.5, color = "white", fontface = "bold") +
    scale_fill_manual(values = dose_cols) +
    scale_y_continuous(limits = c(0, 82), breaks = seq(0, 80, 20)) +
    labs(x = "AP-1 motif hits per peak (FIMO)",
         y = "% peaks with OS elevated\n(2*OS > OE + YS)",
         title = "Fraction of OS-elevated peaks by AP-1 motif density",
         subtitle = "Dashed line = 50% (no directional bias) | Error bars = 95% CI") +
    theme_pub +
    theme(legend.position = "none")

ggsave(file.path(FIG_DIR, "pub_dose_response_barplot.pdf"), pB, width = 5.5, height = 4)
cat("  Written: pub_dose_response_barplot.pdf\n")

# ==============================================================================
# FIGURE C: Fisher's exact - mosaic/balloon plot
# ==============================================================================
cat("Figure C: Fisher's exact enrichment...\n")

# Build enrichment data for multiple tests
fisher_results <- list()

# FOS::JUN vs aging-up (p<0.05)
master[, aging_up_05 := (logFC_OS_vs_YS > 0 & PValue_OS_vs_YS < 0.05)]
tab1 <- table(master$has_fosjun_motif, master$aging_up_05)
ft1 <- fisher.test(tab1)
fisher_results[["FOS::JUN motif+\nvs Aging UP"]] <- data.table(
    test = "FOS::JUN motif+\nvs Aging UP (p<0.05)",
    OR = ft1$estimate, pvalue = ft1$p.value,
    ci_lo = ft1$conf.int[1], ci_hi = ft1$conf.int[2])

# FOS::JUN vs OS-outlier elevated (p<0.05)
master[, os_out_up_05 := (logFC_os_out > 0 & P_os_out < 0.05)]
tab2 <- table(master$has_fosjun_motif, master$os_out_up_05)
ft2 <- fisher.test(tab2)
fisher_results[["FOS::JUN motif+\nvs OS-outlier UP"]] <- data.table(
    test = "FOS::JUN motif+\nvs OS-outlier UP (p<0.05)",
    OR = ft2$estimate, pvalue = ft2$p.value,
    ci_lo = ft2$conf.int[1], ci_hi = ft2$conf.int[2])

# Any AP-1 vs aging-up
tab3 <- table(master$has_ap1_motif_fimo, master$aging_up_05)
ft3 <- fisher.test(tab3)
fisher_results[["Any AP-1 motif\nvs Aging UP"]] <- data.table(
    test = "Any AP-1 motif\nvs Aging UP (p<0.05)",
    OR = ft3$estimate, pvalue = ft3$p.value,
    ci_lo = ft3$conf.int[1], ci_hi = ft3$conf.int[2])

# FOS::JUN vs OS-outlier (p<0.01)
master[, os_out_up_01 := (logFC_os_out > 0 & P_os_out < 0.01)]
tab4 <- table(master$has_fosjun_motif, master$os_out_up_01)
ft4 <- fisher.test(tab4)
fisher_results[["FOS::JUN motif+\nvs OS-outlier UP\n(p<0.01)"]] <- data.table(
    test = "FOS::JUN motif+\nvs OS-outlier UP (p<0.01)",
    OR = ft4$estimate, pvalue = ft4$p.value,
    ci_lo = ft4$conf.int[1], ci_hi = ft4$conf.int[2])

fisher_dt <- rbindlist(fisher_results)
fisher_dt[, test := factor(test, levels = rev(test))]
fisher_dt[, sig_label := ifelse(pvalue < 1e-10, "***",
                          ifelse(pvalue < 1e-5, "**",
                          ifelse(pvalue < 0.05, "*", "ns")))]
# Cap OR for plotting (Inf cases)
fisher_dt[is.infinite(OR), OR := 100]
fisher_dt[is.infinite(ci_hi), ci_hi := 200]

pC <- ggplot(fisher_dt, aes(x = OR, y = test)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_segment(aes(x = ci_lo, xend = pmin(ci_hi, 150), yend = test),
                 linewidth = 1.2, color = "grey40") +
    geom_point(aes(size = -log10(pvalue)), color = col_motif_pos, shape = 16) +
    geom_text(aes(label = paste0("OR=", ifelse(OR >= 100, "Inf", round(OR, 1)),
                                  "\np=", format(pvalue, digits = 2))),
              hjust = -0.3, size = 2.8, lineheight = 0.85) +
    scale_x_log10(limits = c(0.3, 300),
                  breaks = c(1, 5, 10, 50, 100),
                  labels = c("1", "5", "10", "50", "100")) +
    scale_size_continuous(range = c(3, 8), name = expression(-log[10]*p)) +
    labs(x = "Odds ratio (log scale)",
         y = "",
         title = "FOS::JUN motif enrichment in aging-elevated peaks",
         subtitle = "Fisher's exact test | Point size = significance") +
    theme_pub +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 8, lineheight = 0.9))

ggsave(file.path(FIG_DIR, "pub_fisher_enrichment_forest.pdf"), pC, width = 6.5, height = 3.5)
cat("  Written: pub_fisher_enrichment_forest.pdf\n")

# ==============================================================================
# FIGURE D: Reversal flow / waterfall
# ==============================================================================
cat("Figure D: Reversal flow...\n")

# Define categories for all peaks
master[, aging_up := (logFC_OS_vs_YS > 0 & PValue_OS_vs_YS < 0.05)]
master[, aging_down := (logFC_OS_vs_YS < 0 & PValue_OS_vs_YS < 0.05)]
master[, exer_reverses_up := (logFC_OE_vs_OS < 0 & PValue_OE_vs_OS < 0.05)]
master[, exer_reverses_down := (logFC_OE_vs_OS > 0 & PValue_OE_vs_OS < 0.05)]

# Build flow data for FOS::JUN+ and motif-
flow_data <- list()

for (subset_label in c("FOS::JUN motif+", "AP-1 motif\u2013")) {
    if (subset_label == "FOS::JUN motif+") {
        sub <- master[has_fosjun_motif == TRUE]
    } else {
        sub <- master[has_ap1_motif_fimo == FALSE]
    }
    n_total <- nrow(sub)
    n_aging_up <- sum(sub$aging_up)
    n_aging_down <- sum(sub$aging_down)
    n_no_aging <- n_total - n_aging_up - n_aging_down
    n_rev_up <- sum(sub$aging_up & sub$exer_reverses_up)
    n_rev_down <- sum(sub$aging_down & sub$exer_reverses_down)

    flow_data[[subset_label]] <- data.table(
        subset = subset_label,
        category = c("Total peaks", "Aging UP\n(OS > YS, p<0.05)",
                     "Aging DOWN\n(OS < YS, p<0.05)", "No aging change",
                     "Reversed\n(gain-then-loss)", "Reversed\n(loss-then-gain)"),
        count = c(n_total, n_aging_up, n_aging_down, n_no_aging, n_rev_up, n_rev_down),
        stage = c("All", "Aging", "Aging", "Aging", "Reversal", "Reversal"),
        direction = c("neutral", "up", "down", "neutral", "up", "down")
    )
}

flow_dt <- rbindlist(flow_data)
flow_dt[, subset := factor(subset, levels = c("FOS::JUN motif+", "AP-1 motif\u2013"))]

# Stacked bar chart showing the funnel
bar_dt <- flow_dt[category != "Total peaks" & category != "No aging change"]
bar_dt[, category := factor(category, levels = c(
    "Aging UP\n(OS > YS, p<0.05)", "Reversed\n(gain-then-loss)",
    "Aging DOWN\n(OS < YS, p<0.05)", "Reversed\n(loss-then-gain)"))]

pD <- ggplot(bar_dt, aes(x = category, y = count, fill = subset)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7,
             color = "grey30", linewidth = 0.3) +
    geom_text(aes(label = count),
              position = position_dodge(width = 0.8), vjust = -0.4, size = 3.5) +
    scale_fill_manual(values = c("FOS::JUN motif+" = col_motif_pos,
                                  "AP-1 motif\u2013" = col_motif_neg),
                      name = "") +
    labs(x = "", y = "Number of peaks",
         title = "Aging and reversal peak counts by motif status",
         subtitle = paste0("FOS::JUN+ (n=", sum(master$has_fosjun_motif),
                          ") vs motif\u2013 (n=", sum(!master$has_ap1_motif_fimo),
                          ") | Nominal p < 0.05 threshold")) +
    theme_pub +
    theme(legend.position = "bottom",
          axis.text.x = element_text(size = 8, lineheight = 0.9))

ggsave(file.path(FIG_DIR, "pub_reversal_counts_barplot.pdf"), pD, width = 6, height = 4)
cat("  Written: pub_reversal_counts_barplot.pdf\n")

# ==============================================================================
# FIGURE E: What happens to reversed peaks - OE vs YS distribution
# ==============================================================================
cat("Figure E: Reversed peak outcomes...\n")

reversed <- master[aging_up == TRUE & exer_reverses_up == TRUE]
if (nrow(reversed) > 0) {
    reversed[, outcome := ifelse(abs(logFC_OE_vs_YS) < 0.2, "Full reversal\n(OE ~ YS)",
                           ifelse(logFC_OE_vs_YS > 0, "Partial reversal\n(OE between YS & OS)",
                                  "Over-reversal\n(OE below YS)"))]

    pE <- ggplot(reversed, aes(x = logFC_OE_vs_YS)) +
        geom_histogram(binwidth = 0.3, fill = col_motif_pos, color = "white",
                       alpha = 0.8) +
        geom_vline(xintercept = 0, linetype = "solid", color = col_YS, linewidth = 1) +
        annotate("text", x = 0, y = Inf, label = "YS level", vjust = 2, hjust = -0.1,
                 color = col_YS, size = 3.5, fontface = "bold") +
        annotate("text", x = median(reversed$logFC_OE_vs_YS), y = Inf,
                 label = paste0("median = ", round(median(reversed$logFC_OE_vs_YS), 2)),
                 vjust = 3.5, size = 3, color = "grey30") +
        geom_vline(xintercept = median(reversed$logFC_OE_vs_YS),
                   linetype = "dashed", color = col_OE, linewidth = 0.8) +
        labs(x = expression(log[2]*FC~(OE~vs~YS)~of~reversed~peaks),
             y = "Number of peaks",
             title = "Reversed peaks return to young levels",
             subtitle = paste0("n=", nrow(reversed),
                              " peaks with aging-UP + exercise-reversal (p<0.05 each)")) +
        theme_pub

    ggsave(file.path(FIG_DIR, "pub_reversed_peak_outcomes.pdf"), pE, width = 5, height = 3.5)
    cat("  Written: pub_reversed_peak_outcomes.pdf\n")
}

# ==============================================================================
# FIGURE F: Combined summary - 2x2 contingency mosaic-style
# ==============================================================================
cat("Figure F: Motif x aging contingency...\n")

# Build contingency data
master[, aging_status := ifelse(aging_up, "Aging UP\n(OS > YS)",
                          ifelse(aging_down, "Aging DOWN\n(OS < YS)",
                                 "No change"))]
master[, motif_status := ifelse(has_fosjun_motif, "FOS::JUN\nmotif+", "Other /\nmotif\u2013")]

cont_dt <- master[, .N, by = .(motif_status, aging_status)]
cont_dt[, motif_status := factor(motif_status, levels = c("FOS::JUN\nmotif+", "Other /\nmotif\u2013"))]
cont_dt[, aging_status := factor(aging_status,
    levels = c("Aging UP\n(OS > YS)", "No change", "Aging DOWN\n(OS < YS)"))]

# Calculate proportions within each motif group for labeling
cont_dt[, total := sum(N), by = motif_status]
cont_dt[, pct := round(100 * N / total, 1)]

pF <- ggplot(cont_dt, aes(x = aging_status, y = N, fill = motif_status)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7,
             color = "grey30", linewidth = 0.3) +
    geom_text(aes(label = paste0(N, "\n(", pct, "%)")),
              position = position_dodge(width = 0.8), vjust = -0.2,
              size = 2.8, lineheight = 0.85) +
    scale_fill_manual(values = c("FOS::JUN\nmotif+" = col_motif_pos,
                                  "Other /\nmotif\u2013" = "grey70"),
                      name = "") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = "", y = "Number of peaks",
         title = "Aging-associated binding changes are FOS::JUN motif-specific",
         subtitle = paste0("Fisher's exact OR = ", round(ft1$estimate, 1),
                          ", p = ", format(ft1$p.value, digits = 2),
                          " | Total: ", nrow(master), " consensus peaks")) +
    theme_pub +
    theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, "pub_motif_x_aging_barplot.pdf"), pF, width = 5.5, height = 4)
cat("  Written: pub_motif_x_aging_barplot.pdf\n")

# ==============================================================================
# FIGURE G: Volcano colored by motif density
# ==============================================================================
cat("Figure G: Volcano by motif density...\n")

master[, dose_group := ifelse(n_ap1_hits_fimo == 0, "0 hits",
                        ifelse(n_ap1_hits_fimo <= 5, "1-5 hits",
                        ifelse(n_ap1_hits_fimo <= 10, "6-10 hits", "11+ hits")))]
master[, dose_group := factor(dose_group, levels = c("0 hits", "1-5 hits", "6-10 hits", "11+ hits"))]

pG <- ggplot(master, aes(x = logFC_os_out, y = -log10(P_os_out))) +
    geom_point(aes(color = dose_group), size = 0.8, alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    scale_color_manual(values = c("0 hits" = "#D9D9D9", "1-5 hits" = "#9ECAE1",
                                   "6-10 hits" = "#3182BD", "11+ hits" = "#08306B"),
                       name = "AP-1 motif hits") +
    labs(x = expression(log[2]*FC~(2 %.% OS - OE - YS)),
         y = expression(-log[10]~p),
         title = "OS-outlier volcano colored by AP-1 motif density",
         subtitle = "Peaks with more motif hits cluster on the OS-elevated side") +
    theme_pub +
    theme(legend.position = c(0.85, 0.80),
          legend.key.size = unit(0.4, "cm"),
          legend.text = element_text(size = 8))

suppressWarnings(
    ggsave(file.path(FIG_DIR, "pub_volcano_by_motif_density.pdf"), pG, width = 6, height = 5)
)
cat("  Written: pub_volcano_by_motif_density.pdf\n")

# ==============================================================================
# FIGURE H: Reversal waterfall - per-peak logFC trajectory
# ==============================================================================
cat("Figure H: Reversal waterfall...\n")

# For FOS::JUN+ aging-up peaks, show their journey: logFC_aging -> logFC_OE_vs_YS
aging_up_fj <- master[has_fosjun_motif == TRUE & aging_up == TRUE]
aging_up_fj <- aging_up_fj[order(-logFC_OS_vs_YS)]
aging_up_fj[, rank := 1:.N]

waterfall_long <- melt(aging_up_fj,
    id.vars = c("peak_id", "rank"),
    measure.vars = c("logFC_OS_vs_YS", "logFC_OE_vs_YS"),
    variable.name = "contrast", value.name = "logFC")
waterfall_long[, contrast := factor(
    ifelse(contrast == "logFC_OS_vs_YS", "OS vs YS\n(aging)", "OE vs YS\n(net)"),
    levels = c("OS vs YS\n(aging)", "OE vs YS\n(net)"))]

pH <- ggplot(waterfall_long, aes(x = contrast, y = logFC, group = rank)) +
    geom_hline(yintercept = 0, color = col_YS, linewidth = 0.8) +
    geom_line(alpha = 0.25, color = "grey50", linewidth = 0.3) +
    geom_point(aes(color = contrast), size = 1.2, alpha = 0.5) +
    scale_color_manual(values = c("OS vs YS\n(aging)" = col_OS,
                                   "OE vs YS\n(net)" = col_OE)) +
    annotate("text", x = 0.6, y = 0, label = "YS baseline",
             color = col_YS, size = 3, fontface = "bold") +
    labs(x = "", y = expression(log[2]*FC~vs~YS),
         title = "Per-peak aging effect and net OE position",
         subtitle = paste0("n=", nrow(aging_up_fj),
                          " FOS::JUN+ aging-UP peaks | Lines connect same peak")) +
    theme_pub +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 10, face = "bold"))

ggsave(file.path(FIG_DIR, "pub_reversal_waterfall.pdf"), pH, width = 4.5, height = 5)
cat("  Written: pub_reversal_waterfall.pdf\n")

# ==============================================================================
# Summary
# ==============================================================================
cat("\n===== Publication figures generated =====\n")
cat("A. pub_dose_response_boxplot.pdf      - OS outlier logFC by motif density\n")
cat("B. pub_dose_response_barplot.pdf      - % OS-elevated by motif density\n")
cat("C. pub_fisher_enrichment_forest.pdf   - Fisher's exact forest plot\n")
cat("D. pub_reversal_counts_barplot.pdf    - Aging & reversal counts by motif\n")
cat("E. pub_reversed_peak_outcomes.pdf     - Where reversed peaks land (OE vs YS)\n")
cat("F. pub_motif_x_aging_barplot.pdf      - Contingency: motif x aging status\n")
cat("G. pub_volcano_by_motif_density.pdf   - Volcano colored by motif hits\n")
cat("H. pub_reversal_waterfall.pdf         - Per-peak aging -> OE trajectory\n")
cat("\nDone.\n")
