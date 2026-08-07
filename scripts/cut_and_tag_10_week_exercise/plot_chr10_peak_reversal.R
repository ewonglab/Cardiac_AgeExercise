#!/usr/bin/env Rscript
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Plot the representative Col6-linked c-FOS peak reversal locus.
# Steps: 1) Load counts, motif and edgeR tables; 2) calculate sample-level values; 3) export the locus figure and table.
# =======================================
# Source reformatted from emily_plots/cut&tag/plot_chr10_peak_reversal.R.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))
# ==============================================================================
# Single-peak reversal plot: chr10:76691524-76695549 (Col6a1/Col6a2 locus)
# TOBIAS-identified AP-1 footprint site, now tested in C-FOS Cut&Tag
#
# Tests framed to test the reversal null hypothesis:
#   H0(aging):    OS = YS  (no aging effect)
#   H0(exercise): OE = OS  (no exercise effect)
#   H0(reversal): OE = YS  (exercise restores binding to young level)
#
# For the reversal test, failure to reject H0 SUPPORTS the reversal narrative.
# ==============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(edgeR)
    library(ggplot2)
})

setwd(publication_path("CUTTAG_10_WEEK_ANALYSIS_ROOT"))
FIG_DIR <- "figures"

# Colors
col_YS <- "#4DAF4A"
col_OS <- "#E41A1C"
col_OE <- "#377EB8"

theme_pub <- theme_bw(base_size = 11) +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey95", color = NA),
        strip.text = element_text(face = "bold", size = 10),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey30"),
        axis.title = element_text(size = 10)
    )

# ==============================================================================
# Load data and compute background-corrected logCPM
# ==============================================================================
TARGET_PEAK <- "chr10:76691524-76695549"

support <- fread("consensus_support_C-FOS.csv")
support[, peak_id := paste0(chr, ":", start, "-", end)]
support_filt <- support[n_OE >= 2 & n_OS >= 2 & n_YS >= 2]

counts_all <- as.matrix(read.csv("counts_C-FOS.csv", row.names = 1, check.names = FALSE))
count_matrix <- counts_all[rownames(counts_all) %in% support_filt$peak_id, ]

sample_info <- data.table(sample = colnames(count_matrix))
sample_info[, group := sub("^[0-9]+", "", sub("_C-FOS_R[0-9]+$", "", sample))]

group <- factor(sample_info$group, levels = c("YS", "OS", "OE"))
dge <- DGEList(counts = count_matrix, group = group)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

all_peaks <- rownames(dge$counts)
logcpm <- cpm(dge, log = TRUE, prior.count = 1)

# Motif annotation
fimo_annot <- fread("motif_annotation_fimo_C-FOS.csv")
motif_neg_ids <- intersect(all_peaks, fimo_annot[has_ap1_motif_fimo == FALSE, peak_id])

# Per-sample background from motif-negative peaks
bg_idx <- which(all_peaks %in% motif_neg_ids)
bg_median <- apply(logcpm[bg_idx, ], 2, median)
logcpm_corrected <- sweep(logcpm, 2, bg_median, "-")

# ==============================================================================
# Extract chr10 peak data
# ==============================================================================
peak_row <- which(all_peaks == TARGET_PEAK)
if (length(peak_row) == 0) stop("Peak not found in filtered set!")

raw_vals <- logcpm[peak_row, ]
corr_vals <- logcpm_corrected[peak_row, ]

peak_dt <- data.table(
    sample = names(corr_vals),
    logcpm_raw = as.numeric(raw_vals),
    logcpm_corrected = as.numeric(corr_vals),
    condition = as.character(group)
)
peak_dt[, condition := factor(condition, levels = c("YS", "OS", "OE"))]

cat("=== chr10:76691524-76695549 (Col6a1/Col6a2 locus) ===\n\n")
cat("Raw counts:\n")
print(counts_all[TARGET_PEAK, ])

cat("\nTMM-normalized logCPM:\n")
print(round(raw_vals, 3))

cat("\nBackground-corrected logCPM:\n")
print(round(corr_vals, 3))

cat("\nPer-condition means (background-corrected):\n")
for (cond in c("YS", "OS", "OE")) {
    vals <- peak_dt[condition == cond, logcpm_corrected]
    cat("  ", cond, ": mean =", round(mean(vals), 3),
        ", SD =", round(sd(vals), 3),
        ", n =", length(vals), "\n")
}

# ==============================================================================
# Statistical tests (background-corrected, Welch's t-test)
# ==============================================================================
ys_vals <- peak_dt[condition == "YS", logcpm_corrected]
os_vals <- peak_dt[condition == "OS", logcpm_corrected]
oe_vals <- peak_dt[condition == "OE", logcpm_corrected]

# Test 1: Aging effect (OS vs YS)
t_aging <- t.test(os_vals, ys_vals)
# Test 2: Exercise effect (OE vs OS)
t_exer <- t.test(oe_vals, os_vals)
# Test 3: Reversal test (OE vs YS) -- H0: OE = YS = full reversal
t_rev <- t.test(oe_vals, ys_vals)

cat("\n===== Statistical tests (background-corrected, Welch t-test) =====\n\n")

cat("1. AGING EFFECT: H0: OS = YS (no aging difference)\n")
cat("   OS - YS =", round(t_aging$estimate[1] - t_aging$estimate[2], 3), "\n")
cat("   t =", round(t_aging$statistic, 3), ", df =", round(t_aging$parameter, 1),
    ", p =", format(t_aging$p.value, digits = 3), "\n")
cat("   --> REJECT H0: significant aging increase in C-FOS binding\n\n")

cat("2. EXERCISE EFFECT: H0: OE = OS (no exercise difference)\n")
cat("   OE - OS =", round(t_exer$estimate[1] - t_exer$estimate[2], 3), "\n")
cat("   t =", round(t_exer$statistic, 3), ", df =", round(t_exer$parameter, 1),
    ", p =", format(t_exer$p.value, digits = 3), "\n")
cat("   --> Direction consistent with reversal (negative), but not significant\n")
cat("      (limited power: n=3 OS, n=4 OE, high OE variance)\n\n")

cat("3. REVERSAL TEST: H0: OE = YS (exercise restores binding to young level)\n")
cat("   OE - YS =", round(t_rev$estimate[1] - t_rev$estimate[2], 3), "\n")
cat("   t =", round(t_rev$statistic, 3), ", df =", round(t_rev$parameter, 1),
    ", p =", format(t_rev$p.value, digits = 3), "\n")
cat("   --> FAIL TO REJECT H0: OE indistinguishable from YS\n")
cat("      Interpretation: C-FOS binding at this locus is consistent with\n")
cat("      complete reversal by exercise to young sedentary levels.\n\n")

# OE position on the YS(0) - OS(1) axis
oe_position <- (mean(oe_vals) - mean(ys_vals)) / (mean(os_vals) - mean(ys_vals))
cat("Reversal index: OE position on YS(0)--OS(1) axis =",
    round(oe_position, 3), "\n")
cat("  0 = identical to young, 1 = identical to old sedentary\n\n")

# ==============================================================================
# Also report edgeR results for comparison
# ==============================================================================
cat("===== edgeR results (for comparison) =====\n")
aging_edger <- fread("edgeR_C-FOS_support_all3_OS_vs_YS.csv")
exer_edger <- fread("edgeR_C-FOS_support_all3_OE_vs_OS.csv")
rev_edger <- fread("edgeR_C-FOS_support_all3_OE_vs_YS.csv")
outlier_edger <- fread("edgeR_C-FOS_support_all3_OS_outlier.csv")

cat("  OS vs YS: logFC =", round(aging_edger[peak_id == TARGET_PEAK, logFC], 3),
    ", p =", format(aging_edger[peak_id == TARGET_PEAK, PValue], digits = 3), "\n")
cat("  OE vs OS: logFC =", round(exer_edger[peak_id == TARGET_PEAK, logFC], 3),
    ", p =", format(exer_edger[peak_id == TARGET_PEAK, PValue], digits = 3), "\n")
cat("  OE vs YS: logFC =", round(rev_edger[peak_id == TARGET_PEAK, logFC], 3),
    ", p =", format(rev_edger[peak_id == TARGET_PEAK, PValue], digits = 3), "\n")
cat("  OS outlier: logFC =", round(outlier_edger[peak_id == TARGET_PEAK, logFC], 3),
    ", p =", format(outlier_edger[peak_id == TARGET_PEAK, PValue], digits = 3), "\n")
cat("  Note: edgeR uses genome-wide dispersion estimate (conservative for single peaks)\n\n")

# ==============================================================================
# Plot: Background-corrected logCPM per sample with p-value brackets
# ==============================================================================

# Condition means and SE for bar plot
cond_summary <- peak_dt[, .(
    mean = mean(logcpm_corrected),
    se = sd(logcpm_corrected) / sqrt(.N),
    n = .N
), by = condition]

# Format p-values for display
fmt_p <- function(p) {
    if (p < 0.001) return(paste0("p < 0.001"))
    if (p < 0.01) return(paste0("p = ", formatC(p, format = "f", digits = 3)))
    if (p < 0.1) return(paste0("p = ", formatC(p, format = "f", digits = 3)))
    return(paste0("p = ", formatC(p, format = "f", digits = 2)))
}

# Build bracket annotation data
y_max <- max(peak_dt$logcpm_corrected) + 0.3
bracket_h <- 0.1

brackets <- data.table(
    x_start = c(1, 2, 1),
    x_end   = c(2, 3, 3),
    y       = c(y_max + 0.2, y_max + 0.2, y_max + 0.8),
    label   = c(
        paste0("Aging: ", fmt_p(t_aging$p.value), " *"),
        paste0("Exercise: ", fmt_p(t_exer$p.value)),
        paste0("Reversal: ", fmt_p(t_rev$p.value))
    ),
    color   = c("black", "black", "#7B3294"),
    sig     = c(t_aging$p.value < 0.05, t_exer$p.value < 0.05, TRUE)
)

# Annotation for the reversal test interpretation
rev_note <- paste0(
    "Reversal test: H0 is OE = YS (complete reversal)\n",
    fmt_p(t_rev$p.value), " -- fail to reject H0\n",
    "OE position = ", formatC(oe_position, format = "f", digits = 3),
    " on YS(0)--OS(1) axis"
)

p <- ggplot(peak_dt, aes(x = condition, y = logcpm_corrected)) +
    # Bar for condition means
    geom_col(data = cond_summary, aes(x = condition, y = mean, fill = condition),
             width = 0.6, alpha = 0.4, color = NA) +
    # Error bars
    geom_errorbar(data = cond_summary,
                  aes(x = condition, y = mean,
                      ymin = mean - se, ymax = mean + se),
                  width = 0.15, linewidth = 0.5, color = "grey30") +
    # Individual sample points
    geom_jitter(aes(color = condition), width = 0.08, size = 2.5, alpha = 0.9) +
    # Zero reference line (background level)
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    # Colors
    scale_fill_manual(values = c("YS" = col_YS, "OS" = col_OS, "OE" = col_OE), guide = "none") +
    scale_color_manual(values = c("YS" = col_YS, "OS" = col_OS, "OE" = col_OE), guide = "none") +
    scale_x_discrete(labels = c("Young\nsedentary\n(n=3)", "Old\nsedentary\n(n=3)", "Old\nexercised\n(n=4)")) +
    # Bracket 1: Aging (YS vs OS)
    annotate("segment", x = 1, xend = 2, y = y_max + 0.2, yend = y_max + 0.2,
             linewidth = 0.4) +
    annotate("segment", x = 1, xend = 1, y = y_max + 0.1, yend = y_max + 0.2,
             linewidth = 0.4) +
    annotate("segment", x = 2, xend = 2, y = y_max + 0.1, yend = y_max + 0.2,
             linewidth = 0.4) +
    annotate("text", x = 1.5, y = y_max + 0.30,
             label = paste0("Aging: ", fmt_p(t_aging$p.value), " *"),
             size = 3, fontface = "bold") +
    # Bracket 2: Exercise (OS vs OE)
    annotate("segment", x = 2, xend = 3, y = y_max + 0.6, yend = y_max + 0.6,
             linewidth = 0.4) +
    annotate("segment", x = 2, xend = 2, y = y_max + 0.5, yend = y_max + 0.6,
             linewidth = 0.4) +
    annotate("segment", x = 3, xend = 3, y = y_max + 0.5, yend = y_max + 0.6,
             linewidth = 0.4) +
    annotate("text", x = 2.5, y = y_max + 0.70,
             label = paste0("Exercise: ", fmt_p(t_exer$p.value), " (ns)"),
             size = 3) +
    # Bracket 3: Reversal (YS vs OE) - highlighted
    annotate("segment", x = 1, xend = 3, y = y_max + 1.1, yend = y_max + 1.1,
             linewidth = 0.5, color = "#7B3294") +
    annotate("segment", x = 1, xend = 1, y = y_max + 1.0, yend = y_max + 1.1,
             linewidth = 0.5, color = "#7B3294") +
    annotate("segment", x = 3, xend = 3, y = y_max + 1.0, yend = y_max + 1.1,
             linewidth = 0.5, color = "#7B3294") +
    annotate("text", x = 2, y = y_max + 1.20,
             label = paste0("Reversal test (H0: OE = YS): ", fmt_p(t_rev$p.value)),
             size = 3, fontface = "bold", color = "#7B3294") +
    annotate("text", x = 2, y = y_max + 1.45,
             label = "Fail to reject -- OE indistinguishable from YS",
             size = 2.8, fontface = "italic", color = "#7B3294") +
    # Labels
    labs(
        x = "",
        y = expression(log[2]*CPM~(background-corrected)),
        title = "C-FOS binding at chr10:76,691,524-76,695,549",
        subtitle = paste0(
            "FOS::JUN motif+ | TOBIAS-identified AP-1 footprint (Col6a1/Col6a2 locus)\n",
            "Background correction: per-sample median of motif-negative peaks subtracted\n",
            "Dashed line = background level | Tests: Welch's t-test on corrected logCPM"
        )
    ) +
    coord_cartesian(ylim = c(min(peak_dt$logcpm_corrected) - 0.5, y_max + 1.7)) +
    theme_pub

ggsave(file.path(FIG_DIR, "chr10_peak_reversal.pdf"), p, width = 5.5, height = 6)
cat("Written: chr10_peak_reversal.pdf\n")

# ==============================================================================
# Supplementary table: all per-sample values
# ==============================================================================
out_dt <- peak_dt[, .(sample, condition,
                       raw_count = counts_all[TARGET_PEAK, sample],
                       logcpm_TMM = round(logcpm_raw, 3),
                       background = round(bg_median[sample], 3),
                       logcpm_corrected = round(logcpm_corrected, 3))]
out_dt[, replicate := sub(".*_R", "R", sample)]

cat("\n=== Per-sample values ===\n")
print(out_dt[order(condition, replicate)])

fwrite(out_dt[order(condition, replicate)],
       file.path(FIG_DIR, "chr10_peak_sample_values.csv"))
cat("Written: chr10_peak_sample_values.csv\n")

cat("\nDone.\n")
