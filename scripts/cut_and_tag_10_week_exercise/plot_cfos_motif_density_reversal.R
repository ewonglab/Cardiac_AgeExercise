#!/usr/bin/env Rscript
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Plot c-FOS aging and exercise effects stratified by FOS::JUN motif status.
# Steps: 1) Load consensus counts and motif annotations; 2) rebuild edgeR contrasts; 3) export the two-panel density figure.
# =======================================
# Source reformatted from emily_plots/cut&tag/pub_density_reversal.R.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))
# ==============================================================================
# Two-panel density: logFC distributions by FOS::JUN motif status
#   Panel A: Aging effect (OS vs YS)
#   Panel B: Exercise effect (OE vs OS)
# ==============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(edgeR)
    library(ggplot2)
    library(patchwork)
})

ANALYSIS_DIR <- publication_path("CUTTAG_10_WEEK_ANALYSIS_ROOT")
setwd(ANALYSIS_DIR)
FIG_DIR <- file.path(ANALYSIS_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Theme & colors ---
theme_pub <- theme_bw(base_size = 11) +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", color = NA),
        strip.text = element_text(face = "bold", size = 10),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey30"),
        axis.title = element_text(size = 10)
    )

col_motif_pos <- "#984EA3"  # purple
col_motif_neg <- "#FF7F00"  # orange

# ==============================================================================
# Rebuild DGE (same pipeline as publication_figures.R)
# ==============================================================================
cat("Setting up data...\n")

support <- fread("consensus_support_C-FOS.csv")
support[, peak_id := paste0(chr, ":", start, "-", end)]
support_filt <- support[n_OE >= 2 & n_OS >= 2 & n_YS >= 2]

counts_all <- as.matrix(read.csv("counts_C-FOS.csv", row.names = 1, check.names = FALSE))
count_matrix <- counts_all[rownames(counts_all) %in% support_filt$peak_id, ]

sample_info <- data.table(sample = colnames(count_matrix))
sample_info[, condition := sub("^([0-9]+[A-Z]+)_.*", "\\1", sample)]
sample_info[, group := sub("^[0-9]+", "", condition)]

group <- factor(sample_info$group, levels = c("YS", "OS", "OE"))
dge <- DGEList(counts = count_matrix, group = group)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

all_peaks <- rownames(dge$counts)

# Motif annotation
fimo_annot <- fread("motif_annotation_fimo_C-FOS.csv")
fj_col <- grep("FOS.*JUN|MA0099", names(fimo_annot), value = TRUE)[1]
fosjun_ids <- intersect(all_peaks, fimo_annot[get(fj_col) == TRUE, peak_id])
motif_neg_ids <- intersect(all_peaks, fimo_annot[has_ap1_motif_fimo == FALSE, peak_id])

# ==============================================================================
# Get logFC for the two contrasts
# ==============================================================================

# Aging: OS vs YS
con_aging <- makeContrasts(OS - YS, levels = design)
res_aging <- glmQLFTest(fit, contrast = con_aging)
tt_aging  <- as.data.table(topTags(res_aging, n = Inf, sort.by = "none")$table,
                            keep.rownames = "peak_id")

# Exercise: OE vs OS
con_exer <- makeContrasts(OE - OS, levels = design)
res_exer <- glmQLFTest(fit, contrast = con_exer)
tt_exer  <- as.data.table(topTags(res_exer, n = Inf, sort.by = "none")$table,
                            keep.rownames = "peak_id")

# Assign motif groups
assign_motif <- function(dt) {
    dt[, motif_group := fifelse(
        peak_id %in% fosjun_ids, "FOS::JUN motif+",
        fifelse(peak_id %in% motif_neg_ids, "No AP-1 motif", NA_character_)
    )]
    dt[!is.na(motif_group)]
}

# Keep only FOS::JUN+ peaks for both contrasts
aging_fj <- tt_aging[peak_id %in% fosjun_ids]
exer_fj  <- tt_exer[peak_id %in% fosjun_ids]

n_fj <- nrow(aging_fj)
cat("\nFOS::JUN+ peaks:", n_fj, "\n")

# Combine into one long table for superimposed density
overlay_dt <- rbind(
    aging_fj[, .(peak_id, logFC, contrast = "Aging (Aged sed. vs Young sed.)")],
    exer_fj[, .(peak_id, logFC, contrast = "Exercise (Aged ex. vs Aged sed.)")]
)

# Medians
overlay_med <- overlay_dt[, .(med = median(logFC)), by = contrast]
cat("\nMedian logFC at FOS::JUN+ peaks:\n")
print(overlay_med)

# Paired Wilcoxon: are the two distributions mirror images?
# (same peaks, aging logFC vs exercise logFC)
merged <- merge(aging_fj[, .(peak_id, aging_logFC = logFC)],
                exer_fj[, .(peak_id, exer_logFC = logFC)],
                by = "peak_id")
cor_test <- cor.test(merged$aging_logFC, merged$exer_logFC, method = "spearman")
cat("\nSpearman correlation (aging vs exercise logFC):",
    round(cor_test$estimate, 3), ", p =", format(cor_test$p.value, digits = 3), "\n")

# Colors for the two contrasts
col_aging <- "#E41A1C"   # red  (aging = gain)
col_exercise <- "#377EB8" # blue (exercise = reversal)

# ==============================================================================
# Single-panel superimposed density
# ==============================================================================

p1 <- ggplot(overlay_dt, aes(x = logFC, fill = contrast)) +
    geom_density(alpha = 0.45, color = NA) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_vline(data = overlay_med, aes(xintercept = med, color = contrast),
               linetype = "solid", linewidth = 0.8, show.legend = FALSE) +
    scale_fill_manual(values = c(
        "Aging (Aged sed. vs Young sed.)" = col_aging,
        "Exercise (Aged ex. vs Aged sed.)" = col_exercise),
        name = "") +
    scale_color_manual(values = c(
        "Aging (Aged sed. vs Young sed.)" = col_aging,
        "Exercise (Aged ex. vs Aged sed.)" = col_exercise)) +
    labs(x = expression(log[2]*fold~change),
         y = "Density",
         title = "Exercise opposes aging at FOS::JUN motif sites",
         subtitle = paste0(
             "FOS::JUN+ peaks (n=", n_fj,
             ") | Solid lines = median logFC | ",
             "Spearman r = ", round(cor_test$estimate, 2)
         )) +
    coord_cartesian(xlim = c(-3, 3)) +
    theme_pub +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 9))

ggsave(file.path(FIG_DIR, "pub_density_aging_vs_exercise_fosjun.pdf"),
       p1, width = 6, height = 4.5)
cat("\nWritten: pub_density_aging_vs_exercise_fosjun.pdf\n")

cat("\nDone.\n")
