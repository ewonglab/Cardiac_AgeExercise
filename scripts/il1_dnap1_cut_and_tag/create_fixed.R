# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Compile IL-1 and dnAP-1 CUT&Tag differential analyses and figures.
# Steps: 1) Load processed peak/count results; 2) calculate contrasts and annotations; 3) export summaries and publication plots.
# =======================================
# Source reformatted from IL1_dnAP1_CUT&Tag/Create-Fixed.ipynb.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

#!/usr/bin/env Rscript
# CM-specific analysis of AAV DN-AP1 Cut&Tag C-FOS differential binding
# Intersects edgeR results with CM-specific TOBIAS footprint peaks
# Tests whether DN-AP1 preferentially reduces binding at CM peaks

library(data.table)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ggplot2)
library(patchwork)
library(dplyr)
# ── Paths ──────────────────────────────────────────────────────────────────────
# When sourced, use script location; when run interactively, use working dir
base_dir <- publication_path("PAPER_OUTPUT_ROOT", "Paper_Tables/AAV_PLot/")
# Parse mark + normalization from command line
# Usage: Rscript cm_specific_analysis.R [mark] [norm]
#   mark: C-FOS, C-JUN, P-C-JUN (default C-FOS)
#   norm: tmm or spikein (default tmm)
cat("Loading TOBIAS data (CM only)...\n")
tobias_file <- file.path(publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Output/Tobias_df_Classified_Full.rds"))
tobias <- readRDS(tobias_file)
tobias <- tobias %>% select(c("peak_chr", "peak_start", "peak_end",
                                         "name", "category", "CellType",
                                         "in_MAS", "in_MAE", "in_YS",
                                         "MAS_score", "MAE_score", "YS_score"))

tobias_cm <- tobias[tobias$CellType == "Cardiomyocyte",]
setDT(tobias_cm)
rm(tobias); gc()
for (mark in c("C-FOS", "C-JUN", "P-C-JUN")){
  for (norm in c("spikein","tmm")){
mark_dir <- tolower(gsub("-", "", mark))  # C-FOS -> cfos, C-JUN -> cjun
edger_subdir <- if (norm == "spikein") {
  file.path(mark, "males_spikein")
} else {
  file.path(mark, "males")
}
edger_file <- file.path(publication_path("AAV_CUTTAG_EDGER_DIR"),edger_subdir,
                        paste0(mark, ".edgeR.AAV9_AP_1_vs_AAV9_Control.full.csv"))

out_dir    <- file.path(base_dir, paste0("cm_specific_results_", mark, "_", norm))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("Mark: %s | Normalization: %s\n", mark, norm))
cat(sprintf("edgeR file: %s\n", edger_file))
cat("Loading edgeR results...\n")
edger <- fread(edger_file)
setnames(edger, c("logFC", "logCPM", "F", "PValue", "FDR",
                   "peak_id", "chr", "start", "end"))

# ── Load TOBIAS CM peaks ──────────────────────────────────────────────────────



# Unique CM ATAC peaks (collapse across motifs)
cm_peaks_all <- unique(tobias_cm[, .(peak_chr, peak_start, peak_end)])
cat(sprintf("  Unique CM ATAC peaks: %d\n", nrow(cm_peaks_all)))

# AP-1 family TFs
ap1_tfs <- c("Fos", "Fosb", "Fosl1", "Fosl2",
             "Jun", "Junb", "Jund",
             "Batf", "Batf3",
             "Atf2", "Atf3", "Atf4", "Atf7",
             "Bach1", "Bach2",
             "Maf", "Mafa", "Mafb", "Maff", "Mafg", "Mafk")

# CM peaks with AP-1 motifs
tobias_cm_ap1 <- tobias_cm[name %in% ap1_tfs]
cm_peaks_ap1 <- unique(tobias_cm_ap1[, .(peak_chr, peak_start, peak_end)])
cat(sprintf("  CM ATAC peaks with AP-1 motifs: %d\n", nrow(cm_peaks_ap1)))

# AP-1 category breakdown in CM
cat("\nAP-1 binding categories in CM peaks:\n")
print(tobias_cm_ap1[, .N, by = category][order(-N)])

# ── Build GRanges ─────────────────────────────────────────────────────────────
edger_gr <- GRanges(seqnames = edger$chr,
                    ranges = IRanges(start = edger$start, end = edger$end))

cm_all_gr <- GRanges(seqnames = cm_peaks_all$peak_chr,
                     ranges = IRanges(start = cm_peaks_all$peak_start,
                                      end = cm_peaks_all$peak_end))

cm_ap1_gr <- GRanges(seqnames = cm_peaks_ap1$peak_chr,
                     ranges = IRanges(start = cm_peaks_ap1$peak_start,
                                      end = cm_peaks_ap1$peak_end))

# ── Overlap: which edgeR peaks overlap CM ATAC peaks? ─────────────────────────
hits_all  <- overlapsAny(edger_gr, cm_all_gr)
hits_ap1  <- overlapsAny(edger_gr, cm_ap1_gr)

edger[, is_cm_peak := hits_all]
edger[, is_cm_ap1_peak := hits_ap1]

cat(sprintf("\n── Overlap summary ──\n"))
cat(sprintf("  edgeR peaks overlapping any CM ATAC peak:     %d / %d (%.1f%%)\n",
            sum(hits_all), nrow(edger), 100 * mean(hits_all)))
cat(sprintf("  edgeR peaks overlapping CM AP-1 ATAC peak:    %d / %d (%.1f%%)\n",
            sum(hits_ap1), nrow(edger), 100 * mean(hits_ap1)))

# ── Test 1: Directional shift — CM peaks vs non-CM peaks ─────────────────────
cat("\n── Test 1: logFC comparison (CM vs non-CM peaks) ──\n")
for (label in c("is_cm_peak", "is_cm_ap1_peak")) {
  cm_logfc    <- edger[get(label) == TRUE, logFC]
  noncm_logfc <- edger[get(label) == FALSE, logFC]

  cat(sprintf("\n  %s:\n", label))
  cat(sprintf("    CM peaks (n=%d):     mean logFC = %.4f, median = %.4f\n",
              length(cm_logfc), mean(cm_logfc), median(cm_logfc)))
  cat(sprintf("    Non-CM peaks (n=%d): mean logFC = %.4f, median = %.4f\n",
              length(noncm_logfc), mean(noncm_logfc), median(noncm_logfc)))

  # Wilcoxon test: are CM peaks shifted toward negative logFC?
  wt <- wilcox.test(cm_logfc, noncm_logfc, alternative = "less")
  cat(sprintf("    Wilcoxon (CM < non-CM): P = %.4g\n", wt$p.value))

  # Also proportion negative
  prop_cm    <- mean(cm_logfc < 0)
  prop_noncm <- mean(noncm_logfc < 0)
  pt <- prop.test(c(sum(cm_logfc < 0), sum(noncm_logfc < 0)),
                  c(length(cm_logfc), length(noncm_logfc)),
                  alternative = "greater")
  cat(sprintf("    Prop negative: CM=%.3f, non-CM=%.3f, prop.test P=%.4g\n",
              prop_cm, prop_noncm, pt$p.value))
}

# ── Test 2: Age-gained AP-1 sites — most expected to decrease ────────────────
cat("\n── Test 2: Age-gained AP-1 sites in CM ──\n")
# Categories indicating aging-specific binding
# "In Aged Sed Only" = bound only in aged sedentary
# "Lost YoungSed" = lost from young but present in aged conditions
age_cats <- c("In Aged Sed Only", "Lost YoungSed")
tobias_age_gained <- tobias_cm_ap1[category %in% age_cats]
# Also include peaks where TOBIAS score increased with aging
tobias_age_increased <- tobias_cm_ap1[MAS_score > YS_score & MAS_score > 0]
age_increased_peaks <- unique(tobias_age_increased[, .(peak_chr, peak_start, peak_end)])
age_peaks <- unique(tobias_age_gained[, .(peak_chr, peak_start, peak_end)])
cat(sprintf("  CM peaks with age-gained AP-1 motifs: %d\n", nrow(age_peaks)))

if (nrow(age_peaks) > 0) {
  age_gr <- GRanges(seqnames = age_peaks$peak_chr,
                    ranges = IRanges(start = age_peaks$peak_start,
                                    end = age_peaks$peak_end))
  hits_age <- overlapsAny(edger_gr, age_gr)
  edger[, is_age_gained_ap1 := hits_age]

  age_logfc <- edger[is_age_gained_ap1 == TRUE, logFC]
  other_logfc <- edger[is_age_gained_ap1 == FALSE, logFC]

  cat(sprintf("  edgeR peaks at age-gained AP-1 sites: %d\n", sum(hits_age)))
  cat(sprintf("  Mean logFC (age-gained): %.4f\n", mean(age_logfc)))
  cat(sprintf("  Mean logFC (other):      %.4f\n", mean(other_logfc)))

  wt <- wilcox.test(age_logfc, other_logfc, alternative = "less")
  cat(sprintf("  Wilcoxon (age-gained < other): P = %.4g\n", wt$p.value))
}

# Also test peaks where TOBIAS score increased with aging (continuous)
cat(sprintf("  CM peaks with age-increased AP-1 score: %d\n",
            nrow(age_increased_peaks)))
if (nrow(age_increased_peaks) > 0) {
  age_inc_gr <- GRanges(seqnames = age_increased_peaks$peak_chr,
                        ranges = IRanges(start = age_increased_peaks$peak_start,
                                         end = age_increased_peaks$peak_end))
  hits_inc <- overlapsAny(edger_gr, age_inc_gr)
  edger[, is_age_increased_ap1 := hits_inc]

  inc_logfc <- edger[is_age_increased_ap1 == TRUE, logFC]
  dec_logfc <- edger[is_age_increased_ap1 == FALSE, logFC]
  cat(sprintf("  edgeR peaks at age-increased AP-1 sites: %d\n", sum(hits_inc)))
  cat(sprintf("  Mean logFC (age-increased): %.4f\n", mean(inc_logfc)))
  cat(sprintf("  Mean logFC (other):         %.4f\n", mean(dec_logfc)))
  wt2 <- wilcox.test(inc_logfc, dec_logfc, alternative = "less")
  cat(sprintf("  Wilcoxon (age-increased < other): P = %.4g\n", wt2$p.value))
}

# ── Test 3: Gene-specific analysis ──────────────────────────────────────────
cat("\n── Test 3: Gene-specific loci ──\n")

# Get gene coordinates from TxDb
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes_of_interest <- c("Cdkn1a", "Dusp1", "Dusp4", "Spry2", "Grb2",
                        "Braf", "Itga5", "Yap1", "Gadd45a",
                        "Fos", "Jun", "Fosb", "Fosl2", "Junb", "Jund")

# Map gene symbols to entrez IDs
symbol_to_entrez <- AnnotationDbi::select(org.Mm.eg.db, keys = genes_of_interest,
                           columns = "ENTREZID", keytype = "SYMBOL")
entrez_ids <- na.omit(symbol_to_entrez$ENTREZID)

# Get gene ranges
all_genes <- genes(txdb)
gene_loci <- all_genes[names(all_genes) %in% entrez_ids]

# Extend by 200kb flanking (to capture distal regulatory elements)
gene_windows <- resize(gene_loci, width = width(gene_loci) + 400000,
                       fix = "center")

# For each gene, find overlapping edgeR peaks
results_list <- list()

for (i in seq_along(genes_of_interest)) {
  sym <- genes_of_interest[i]
  eid <- symbol_to_entrez$ENTREZID[symbol_to_entrez$SYMBOL == sym]
  if (is.na(eid) || length(eid) == 0) {
    cat(sprintf("  %s: gene not found in annotation\n", sym))
    next
  }

  if (!(eid %in% names(gene_loci))) {
    cat(sprintf("  %s: gene not in TxDb\n", sym))
    next
  }

  win <- gene_windows[names(gene_windows) == eid]
  hits_gene <- overlapsAny(edger_gr, win)
  gene_peaks <- edger[hits_gene]

  # Also check which of these overlap CM AP-1 peaks
  gene_peaks_cm <- gene_peaks[is_cm_ap1_peak == TRUE]

  n_total <- nrow(gene_peaks)
  n_cm    <- nrow(gene_peaks_cm)
  n_neg   <- sum(gene_peaks$logFC < 0)
  mean_lfc <- mean(gene_peaks$logFC)
  mean_lfc_cm <- if (n_cm > 0) mean(gene_peaks_cm$logFC) else NA

  cat(sprintf("  %s: %d peaks (200kb window), %d CM-AP1, mean logFC=%.3f",
              sym, n_total, n_cm, mean_lfc))
  if (!is.na(mean_lfc_cm)) {
    cat(sprintf(", CM-AP1 mean logFC=%.3f", mean_lfc_cm))
  }
  cat(sprintf(", %d/%d negative (%.0f%%)\n", n_neg, n_total,
              100 * n_neg / max(n_total, 1)))

  if (n_cm > 0) {
    results_list[[sym]] <- gene_peaks_cm[, .(gene = sym, peak_id, chr, start,
                                              end, logFC, logCPM, PValue, FDR)]
  }
}

# Save gene-level CM-AP1 peak details
if (length(results_list) > 0) {
  gene_results <- rbindlist(results_list)
  fwrite(gene_results,
         file.path(out_dir, "gene_specific_cm_ap1_peaks.csv"))
  cat(sprintf("\n  Saved %d CM-AP1 peaks near genes of interest\n",
              nrow(gene_results)))
}

# ── Test 4: TOBIAS score-weighted analysis ───────────────────────────────────
# For CM-AP1 peaks, check if age-increased TOBIAS score predicts more
# negative logFC in the DN experiment
cat("\n── Test 4: TOBIAS aging delta vs Cut&Tag logFC ──\n")

# Get per-peak mean AP-1 aging delta from TOBIAS
tobias_cm_ap1[, aging_delta := MAS_score - YS_score]
peak_aging <- tobias_cm_ap1[, .(mean_aging_delta = mean(aging_delta),
                                 mean_MAS = mean(MAS_score),
                                 mean_YS = mean(YS_score),
                                 n_motifs = .N),
                             by = .(peak_chr, peak_start, peak_end)]

# Build GRanges for matching
peak_aging_gr <- GRanges(seqnames = peak_aging$peak_chr,
                         ranges = IRanges(start = peak_aging$peak_start,
                                          end = peak_aging$peak_end))

# Find overlaps with edgeR peaks
olaps <- findOverlaps(edger_gr, peak_aging_gr)

if (length(olaps) > 0) {
  matched <- data.table(
    edger_idx = queryHits(olaps),
    tobias_idx = subjectHits(olaps)
  )
  matched[, logFC := edger$logFC[edger_idx]]
  matched[, aging_delta := peak_aging$mean_aging_delta[tobias_idx]]
  matched[, mean_MAS := peak_aging$mean_MAS[tobias_idx]]

  # Correlation: peaks that gained AP-1 with aging should lose it with DN
  cor_test <- cor.test(matched$aging_delta, matched$logFC,
                       method = "spearman")
  cat(sprintf("  Matched peak pairs: %d\n", nrow(matched)))
  cat(sprintf("  Spearman correlation (aging_delta vs DN logFC): rho=%.4f, P=%.4g\n",
              cor_test$estimate, cor_test$p.value))
  cat("  (Negative rho = peaks that gained AP-1 with aging lose it with DN)\n")

  # Split by aging direction
  gained <- matched[aging_delta > 0]
  lost   <- matched[aging_delta < 0]
  cat(sprintf("\n  Age-gained AP-1 sites (delta>0): n=%d, mean DN logFC=%.4f\n",
              nrow(gained), mean(gained$logFC)))
  cat(sprintf("  Age-lost AP-1 sites (delta<0):   n=%d, mean DN logFC=%.4f\n",
              nrow(lost), mean(lost$logFC)))

  if (nrow(gained) > 0 & nrow(lost) > 0) {
    wt <- wilcox.test(gained$logFC, lost$logFC, alternative = "less")
    cat(sprintf("  Wilcoxon (gained < lost): P = %.4g\n", wt$p.value))
  }

  # Save matched data
  fwrite(matched, file.path(out_dir, "tobias_aging_vs_dn_logfc.csv"))
}

# ── Plots ────────────────────────────────────────────────────────────────────
cat("\n── Generating plots ──\n")

# Plot 1: logFC density by peak type
p1 <- ggplot(edger, aes(x = logFC, fill = is_cm_ap1_peak)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "grey60"),
                    labels = c("TRUE" = "CM AP-1", "FALSE" = "Other"),
                    name = "Peak type") +
  labs(title = paste0(mark, " logFC: CM AP-1 peaks vs other peaks"),
       subtitle = "AAV9 DN-AP1 vs Control (males)",
       x = "logFC (DN - Control)", y = "Density") +
  theme_bw(base_size = 12)

# Plot 2: Cumulative distribution
edger_plot <- edger[, .(logFC, type = fifelse(is_cm_ap1_peak,
                                               "CM AP-1 peak", "Other peak"))]
p2 <- ggplot(edger_plot, aes(x = logFC, colour = type)) +
  stat_ecdf(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = c("CM AP-1 peak" = "firebrick",
                                  "Other peak" = "grey40")) +
  labs(title = "Cumulative distribution of logFC",
       x = "logFC (DN - Control)", y = "Cumulative proportion",
       colour = "Peak type") +
  theme_bw(base_size = 12)

# Plot 3: TOBIAS aging delta vs DN logFC (if available)
if (exists("matched") && nrow(matched) > 0) {
  p3 <- ggplot(matched, aes(x = aging_delta, y = logFC)) +
    geom_point(alpha = 0.1, size = 0.5) +
    geom_smooth(method = "lm", colour = "firebrick") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(title = "TOBIAS aging delta vs DN Cut&Tag logFC",
         subtitle = "CM AP-1 peaks",
         x = "TOBIAS aging delta (MAS - YS score)",
         y = "Cut&Tag logFC (DN - Control)") +
    theme_bw(base_size = 12)

  combined <- (p1 | p2) / p3 + plot_annotation(
    title = paste0("AAV DN-AP1: CM-specific ", mark, " binding analysis")
  )
} else {
  combined <- p1 | p2
}

ggsave(file.path(out_dir, "cm_specific_cfos_analysis.pdf"),
       combined, width = 14, height = 10)

# ── Plot 4: Gene-level barplot ───────────────────────────────────────────────
if (length(results_list) > 0) {
  gene_summary <- gene_results[, .(mean_logFC = mean(logFC),
                                    se_logFC = sd(logFC) / sqrt(.N),
                                    n_peaks = .N,
                                    prop_negative = mean(logFC < 0)),
                                by = gene]
  gene_summary[, gene := factor(gene, levels = gene_summary[order(mean_logFC), gene])]

  p4 <- ggplot(gene_summary, aes(x = gene, y = mean_logFC)) +
    geom_col(aes(fill = mean_logFC < 0)) +
    geom_errorbar(aes(ymin = mean_logFC - se_logFC,
                      ymax = mean_logFC + se_logFC), width = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_text(aes(label = sprintf("n=%d", n_peaks)),
              vjust = -0.5, size = 3) +
    scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "coral"),
                      guide = "none") +
    labs(title = "Mean logFC at CM AP-1 peaks near genes of interest",
         subtitle = paste0(mark, " Cut&Tag: AAV9 DN-AP1 vs Control (males)"),
         x = NULL, y = "Mean logFC (DN - Control)") +
    coord_flip() +
    theme_bw(base_size = 12)

  ggsave(file.path(out_dir, "gene_specific_logfc.pdf"),
         p4, width = 8, height = 6)
}

# ── Save full annotated edgeR table ──────────────────────────────────────────
fwrite(edger, file.path(out_dir, "edgeR_cm_annotated.csv"))

cat(sprintf("\n── Done. Results saved to %s ──\n", out_dir))
}
}










suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(ggplot2)
})

# Combined from:
# - Create.html: CM AP-1 TOBIAS binding calls and Cut&Tag peak overlap
# - PLot.html: 200 kb peak-to-gene mapping, mark-specific background, Fisher test

basedir <- publication_path("DNAP1_SNRNA_ROOT", "Process/")
plotdir <- publication_path("PAPER_OUTPUT_ROOT", "Paper_Tables/AAV_PLot/")
tobias_file <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Output/Tobias_df_Classified_Full.rds")
edger_base <- publication_path("AAV_CUTTAG_EDGER_DIR")

window_bp <- 200000
norm <- "spikein"

# Match the example panel. Add "C-FOS" here if you want the C-FOS facet too.
marks <- c("C-JUN", "P-C-JUN")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

category_map <- c(
  "In Aged Sed Only" = "Age Up, Ex Down",
  "Lost YoungSed" = "Age Up, Ex maintained",
  "Lost AgeSed" = "Age Down, Ex Up",
  "In Young Sedentary Only" = "Age Down, Ex maintained"
)

category_levels <- unname(category_map)

ap1_tfs <- c(
  "Fos", "Fosb", "Fosl1", "Fosl2",
  "Jun", "Junb", "Jund",
  "Batf", "Batf3",
  "Atf2", "Atf3", "Atf4", "Atf7",
  "Bach1", "Bach2",
  "Maf", "Mafa", "Mafb", "Maff", "Mafg", "Mafk"
)

mast_file <- file.path(
  basedir,
  "MAST_DEG_by_celltype_CDR_batch_removed_male_sample",
  "MAST_ALL_CELLTYPES_AAV9_AP_1_Dosage_2_vs_AAV9_cTnT_all_genes.csv"
)

mast <- read.csv(mast_file)

cm <- mast %>%
  filter(celltype == "Cardiomyocyte", !is.na(gene), gene != "")

cm_all <- unique(cm$gene)

cm_deg_down <- cm %>%
  filter(mast_FDR < 0.05, mast_coef_logFC < 0) %>%
  pull(gene) %>%
  unique()

cm_deg_up <- cm %>%
  filter(mast_FDR < 0.05, mast_coef_logFC > 0) %>%
  pull(gene) %>%
  unique()

message("CM MAST-tested genes: ", length(cm_all))
message("CM UP DEGs: ", length(cm_deg_up))
message("CM DOWN DEGs: ", length(cm_deg_down))

read_tobias_ap1 <- function(tobias_file) {
  tobias <- readRDS(tobias_file) %>%
    dplyr::select(
      peak_chr, peak_start, peak_end,
      name, category, CellType,
      in_MAS, in_MAE, in_YS,
      MAS_score, MAE_score, YS_score
    ) %>%
    filter(CellType == "Cardiomyocyte", name %in% ap1_tfs)

  as.data.table(tobias)
}

tobias_cm_ap1 <- read_tobias_ap1(tobias_file)

tobias_all_peaks <- unique(tobias_cm_ap1[, .(peak_chr, peak_start, peak_end)])

tobias_category_peaks <- tobias_cm_ap1[
  category %in% names(category_map),
  .(
    peak_chr,
    peak_start,
    peak_end,
    tobias_category = category,
    category_label = unname(category_map[category])
  )
] %>%
  unique()

message("CM AP-1 TOBIAS peaks: ", nrow(tobias_all_peaks))
message("CM AP-1 TOBIAS category peaks used:")
print(tobias_category_peaks[, .N, by = .(tobias_category, category_label)][order(category_label)])

make_granges <- function(df, chr_col, start_col, end_col) {
  GRanges(
    seqnames = df[[chr_col]],
    ranges = IRanges(
      start = as.integer(df[[start_col]]),
      end = as.integer(df[[end_col]])
    )
  )
}

tobias_all_gr <- make_granges(tobias_all_peaks, "peak_chr", "peak_start", "peak_end")
tobias_cat_gr <- make_granges(tobias_category_peaks, "peak_chr", "peak_start", "peak_end")
mcols(tobias_cat_gr)$category_label <- tobias_category_peaks$category_label
mcols(tobias_cat_gr)$tobias_category <- tobias_category_peaks$tobias_category

edger_file_for_mark <- function(mark, norm = "spikein") {
  edger_subdir <- if (norm == "spikein") {
    file.path(mark, "males_spikein")
  } else {
    file.path(mark, "males")
  }

  file.path(
    edger_base,
    edger_subdir,
    paste0(mark, ".edgeR.AAV9_AP_1_vs_AAV9_Control.full.csv")
  )
}

read_edger <- function(edger_file) {
  e <- fread(edger_file)

  if (!all(c("chr", "start", "end") %in% names(e))) {
    if (ncol(e) < 9) {
      stop("edgeR file does not contain chr/start/end columns and has fewer than 9 columns: ", edger_file)
    }
    setnames(
      e,
      old = names(e)[seq_len(9)],
      new = c("logFC", "logCPM", "F", "PValue", "FDR", "peak_id", "chr", "start", "end")
    )
  }

  e[, peak_row_id := .I]
  e
}

annotate_peak_table_to_genes <- function(peak_table, window_bp) {
  if (nrow(peak_table) == 0) {
    return(data.frame())
  }

  gr <- make_granges(peak_table, "chr", "start", "end")

  metadata_cols <- setdiff(colnames(peak_table), c("chr", "start", "end"))
  for (col in metadata_cols) {
    mcols(gr)[[col]] <- peak_table[[col]]
  }

  anno <- as.data.frame(
    annotatePeak(
      gr,
      TxDb = txdb,
      annoDb = "org.Mm.eg.db",
      level = "gene"
    )
  )

  anno %>%
    mutate(gene = SYMBOL) %>%
    filter(!is.na(gene), gene != "") %>%
    filter(!is.na(distanceToTSS), abs(distanceToTSS) <= window_bp)
}

map_mark_to_tobias_gene_sets <- function(mark, norm, window_bp) {
  edger_file <- edger_file_for_mark(mark, norm)
  message("\nReading ", mark, " ", norm, ": ", edger_file)

  e <- read_edger(edger_file)
  edger_gr <- make_granges(e, "chr", "start", "end")

  hits_all <- findOverlaps(edger_gr, tobias_all_gr)
  bg_peak_ids <- unique(queryHits(hits_all))

  bg_peaks <- e[bg_peak_ids]
  bg_gene_map <- annotate_peak_table_to_genes(bg_peaks, window_bp)

  bg_mark <- intersect(unique(bg_gene_map$gene), cm_all)

  hits_cat <- findOverlaps(edger_gr, tobias_cat_gr)
  cat_hit_table <- data.table(
    peak_row_id = queryHits(hits_cat),
    category_label = mcols(tobias_cat_gr)$category_label[subjectHits(hits_cat)],
    tobias_category = mcols(tobias_cat_gr)$tobias_category[subjectHits(hits_cat)]
  ) %>%
    unique(by = c("peak_row_id", "category_label"))

  cat_peaks <- merge(e, cat_hit_table, by = "peak_row_id")
  cat_gene_map <- annotate_peak_table_to_genes(cat_peaks, window_bp)

  category_gene_sets <- lapply(category_levels, function(category_label) {
    if (!"category_label" %in% colnames(cat_gene_map)) {
      return(character(0))
    }

    intersect(
      unique(cat_gene_map$gene[cat_gene_map$category_label == category_label]),
      bg_mark
    )
  })
  names(category_gene_sets) <- category_levels

  list(
    bg_mark = bg_mark,
    category_gene_sets = category_gene_sets,
    bg_gene_map = bg_gene_map,
    cat_gene_map = cat_gene_map
  )
}

fisher_enrich <- function(category_genes, deg_genes, bg) {
  bg <- unique(bg)
  category_genes <- intersect(unique(category_genes), bg)
  deg_genes <- intersect(unique(deg_genes), bg)

  a <- length(intersect(category_genes, deg_genes))
  b <- length(setdiff(category_genes, deg_genes))
  c <- length(setdiff(deg_genes, category_genes))
  d <- length(setdiff(bg, union(category_genes, deg_genes)))

  out_base <- data.frame(
    n_category_genes = length(category_genes),
    n_deg_genes_in_background = length(deg_genes),
    overlap = a,
    background_n = length(bg)
  )

  # Continuity-corrected log2 OR and CI keep the bar plot finite even with zeros.
  log2_or <- log2(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)))
  log2_or_se <- sqrt(
    1 / (a + 0.5) +
      1 / (b + 0.5) +
      1 / (c + 0.5) +
      1 / (d + 0.5)
  ) / log(2)

  out_base$log2_OR <- log2_or
  out_base$log2_OR_low <- log2_or - 1.96 * log2_or_se
  out_base$log2_OR_high <- log2_or + 1.96 * log2_or_se

  if (length(bg) == 0 || length(category_genes) == 0 || length(deg_genes) == 0) {
    out_base$OR <- NA_real_
    out_base$pval <- NA_real_
    return(out_base)
  }

  mat <- matrix(
    c(a, b, c, d),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      category_status = c("In_TOBIAS_category", "Other_TOBIAS_AP1_peak_genes"),
      DEG_status = c("DEG", "Not_DEG")
    )
  )

  ft <- fisher.test(mat)
  out_base$OR <- as.numeric(ft$estimate)
  out_base$pval <- ft$p.value
  out_base
}

mark_maps <- lapply(marks, map_mark_to_tobias_gene_sets, norm = norm, window_bp = window_bp)
names(mark_maps) <- marks

res <- list()

for (mk in marks) {
  bg_mark <- mark_maps[[mk]]$bg_mark
  category_gene_sets <- mark_maps[[mk]]$category_gene_sets

  message(mk, " mark-specific 200 kb background genes: ", length(bg_mark))

  for (category_label in category_levels) {
    for (dir in c("DOWN", "UP")) {
      deg <- if (dir == "DOWN") cm_deg_down else cm_deg_up

      res[[length(res) + 1]] <- cbind(
        mark = mk,
        category_label = category_label,
        DEG_direction = dir,
        fisher_enrich(
          category_genes = category_gene_sets[[category_label]],
          deg_genes = deg,
          bg = bg_mark
        )
      )
    }
  }
}

fisher_df <- bind_rows(res) %>%
  mutate(
    padj = p.adjust(pval, method = "BH"),
    category_label = factor(category_label, levels = category_levels),
    DEG_direction = factor(DEG_direction, levels = c("DOWN", "UP")),
    mark = factor(mark, levels = marks),
    mark_label = recode(
      as.character(mark),
      "C-FOS" = "C-Fos",
      "C-JUN" = "C-Jun",
      "P-C-JUN" = "P-C-Jun"
    ),
    deg_label = recode(
      as.character(DEG_direction),
      "DOWN" = "CM DEG down",
      "UP" = "CM DEG up"
    ),
    sig_star = case_when(
      !is.na(padj) & padj < 0.001 ~ "***",
      !is.na(padj) & padj < 0.01 ~ "**",
      !is.na(padj) & padj < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

out_csv <- file.path(plotdir, "Fig5b_TOBIAS_binding_category_200kb_enrichment.csv")
write.csv(fisher_df, out_csv, row.names = FALSE)

yr <- max(abs(c(fisher_df$log2_OR_low, fisher_df$log2_OR_high)), na.rm = TRUE)

plot_df <- fisher_df %>%
  mutate(
    star_y = if_else(log2_OR >= 0, log2_OR_high + 0.05 * yr, log2_OR_low - 0.05 * yr)
  )

p <- ggplot(plot_df, aes(category_label, log2_OR)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "black") +
  geom_col(width = 0.65, fill = "grey85", colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = log2_OR_low, ymax = log2_OR_high), width = 0.18, linewidth = 0.3) +
  geom_text(aes(y = star_y, label = sig_star), size = 2.7, fontface = "bold") +
  facet_grid(deg_label ~ mark_label) +
  coord_cartesian(ylim = c(-yr * 1.15, yr * 1.15), clip = "off") +
  scale_x_discrete(labels = function(x) gsub(", ", "\n", x, fixed = TRUE)) +
  labs(
    x = NULL,
    y = expression(log[2]~"(odds ratio)")
  ) +
  theme_classic(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "#ded8cf", colour = NA),
    strip.text = element_text(face = "bold", size = 9),
    panel.spacing.x = unit(1.5, "mm"),
    panel.spacing.y = unit(2.0, "mm"),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.3),
    legend.position = "none",
    plot.margin = margin(4, 8, 4, 4)
  )

out_pdf <- file.path(plotdir, "Fig5b_TOBIAS_binding_category_200kb_barplot.pdf")
ggsave(out_pdf, p, width = 6.2, height = 4.2)

message("Saved: ", out_csv)
message("Saved: ", out_pdf)

fisher_df

DT::datatable(fisher_df, filter = 'top')

library(dplyr)
library(ggplot2)

dodge <- position_dodge(width = 0.75)

yr <- max(abs(c(fisher_df$log2_OR_low, fisher_df$log2_OR_high)), na.rm = TRUE)

plot_df <- fisher_df %>%
  mutate(
    deg_label = recode(as.character(DEG_direction),
      DOWN = "Down",
      UP = "Up"
    ),
    star_y = if_else(log2_OR >= 0,
      log2_OR_high + 0.05 * yr,
      log2_OR_low - 0.05 * yr
    )
  )

p <- ggplot(plot_df, aes(category_label, log2_OR, fill = deg_label)) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  geom_col(position = dodge, width = 0.65, colour = "black", linewidth = 0.25) +
  geom_errorbar(
    aes(ymin = log2_OR_low, ymax = log2_OR_high),
    position = dodge,
    width = 0.18,
    linewidth = 0.3
  ) +
  geom_text(
    aes(y = star_y, label = sig_star),
    position = dodge,
    size = 2.7,
    fontface = "bold"
  ) +
  facet_wrap(~ mark_label, nrow = 1) +
  scale_fill_manual(values = c("Down" = "light grey", "Up" = "dark grey")) +
  scale_x_discrete(labels = function(x) gsub(", ", "\n", x, fixed = TRUE)) +
  coord_cartesian(ylim = c(-yr * 1.15, yr * 1.15), clip = "off") +
  labs(x = NULL, y = expression(log[2]~"(odds ratio)"), fill = "CM DEG") +
  theme_classic(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "#ded8cf", colour = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1)
  )

p
                   
     ggsave(publication_path("PAPER_OUTPUT_ROOT", "Paper_Tables/AAV_PLot/Fig5b_TOBIAS_binding_category_200kb_barplot_flat.pdf"), height = 2, width = 5)

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(ggplot2)
})

basedir <- publication_path("DNAP1_SNRNA_ROOT", "Process/")
plotdir <- publication_path("PAPER_OUTPUT_ROOT", "Paper_Tables/AAV_PLot/")
tobias_file <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Output/Tobias_df_Classified_Full.rds")
edger_base <- publication_path("AAV_CUTTAG_EDGER_DIR")

window_bp <- 200000
norm <- "spikein"
use_dn_binding_fdr <- FALSE
dn_binding_fdr_cutoff <- 0.05

marks <- c("C-JUN", "P-C-JUN")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

category_map <- c(
  "In Aged Sed Only" = "Age Up, Ex Down",
  "Lost YoungSed" = "Age Up, Ex maintained",
  "Lost AgeSed" = "Age Down, Ex Up",
  "In Young Sedentary Only" = "Age Down, Ex maintained"
)

category_levels <- unname(category_map)
dn_category_levels <- c("DN-AP1 down", "DN-AP1 up")
plot_category_levels <- c(category_levels, dn_category_levels)

ap1_tfs <- c(
  "Fos", "Fosb", "Fosl1", "Fosl2",
  "Jun", "Junb", "Jund",
  "Batf", "Batf3",
  "Atf2", "Atf3", "Atf4", "Atf7",
  "Bach1", "Bach2",
  "Maf", "Mafa", "Mafb", "Maff", "Mafg", "Mafk"
)

mast_file <- file.path(
  basedir,
  "MAST_DEG_by_celltype_CDR_batch_removed_male_sample",
  "MAST_ALL_CELLTYPES_AAV9_AP_1_Dosage_2_vs_AAV9_cTnT_all_genes.csv"
)

mast <- read.csv(mast_file)

cm <- mast %>%
  filter(celltype == "Cardiomyocyte", !is.na(gene), gene != "")

cm_all <- unique(cm$gene)

cm_deg_down <- cm %>%
  filter(mast_FDR < 0.05, mast_coef_logFC < 0) %>%
  pull(gene) %>%
  unique()

cm_deg_up <- cm %>%
  filter(mast_FDR < 0.05, mast_coef_logFC > 0) %>%
  pull(gene) %>%
  unique()

read_tobias_ap1 <- function(tobias_file) {
  tobias <- readRDS(tobias_file) %>%
    dplyr::select(
      peak_chr, peak_start, peak_end,
      name, category, CellType,
      in_MAS, in_MAE, in_YS,
      MAS_score, MAE_score, YS_score
    ) %>%
    filter(CellType == "Cardiomyocyte", name %in% ap1_tfs)

  as.data.table(tobias)
}

tobias_cm_ap1 <- read_tobias_ap1(tobias_file)

tobias_all_peaks <- unique(tobias_cm_ap1[, .(peak_chr, peak_start, peak_end)])

tobias_category_peaks <- tobias_cm_ap1[
  category %in% names(category_map),
  .(
    peak_chr,
    peak_start,
    peak_end,
    tobias_category = category,
    category_label = unname(category_map[category])
  )
] %>%
  unique()

make_granges <- function(df, chr_col, start_col, end_col) {
  GRanges(
    seqnames = df[[chr_col]],
    ranges = IRanges(
      start = as.integer(df[[start_col]]),
      end = as.integer(df[[end_col]])
    )
  )
}

tobias_all_gr <- make_granges(tobias_all_peaks, "peak_chr", "peak_start", "peak_end")
tobias_cat_gr <- make_granges(tobias_category_peaks, "peak_chr", "peak_start", "peak_end")

mcols(tobias_cat_gr)$category_label <- tobias_category_peaks$category_label
mcols(tobias_cat_gr)$tobias_category <- tobias_category_peaks$tobias_category

edger_file_for_mark <- function(mark, norm = "spikein") {
  edger_subdir <- if (norm == "spikein") {
    file.path(mark, "males_spikein")
  } else {
    file.path(mark, "males")
  }

  file.path(
    edger_base,
    edger_subdir,
    paste0(mark, ".edgeR.AAV9_AP_1_vs_AAV9_Control.full.csv")
  )
}

read_edger <- function(edger_file) {
  e <- fread(edger_file)

  if (!all(c("chr", "start", "end") %in% names(e))) {
    setnames(
      e,
      old = names(e)[seq_len(9)],
      new = c("logFC", "logCPM", "F", "PValue", "FDR", "peak_id", "chr", "start", "end")
    )
  }

  e[, peak_row_id := .I]
  e
}

annotate_peak_table_to_genes <- function(peak_table, window_bp) {
  if (nrow(peak_table) == 0) {
    return(data.frame())
  }

  gr <- make_granges(peak_table, "chr", "start", "end")

  metadata_cols <- setdiff(colnames(peak_table), c("chr", "start", "end"))
  for (col in metadata_cols) {
    mcols(gr)[[col]] <- peak_table[[col]]
  }

  anno <- as.data.frame(
    annotatePeak(
      gr,
      TxDb = txdb,
      annoDb = "org.Mm.eg.db",
      level = "gene"
    )
  )

  anno %>%
    mutate(gene = SYMBOL) %>%
    filter(!is.na(gene), gene != "") %>%
    filter(!is.na(distanceToTSS), abs(distanceToTSS) <= window_bp)
}

map_mark_to_tobias_gene_sets <- function(mark, norm, window_bp) {
  edger_file <- edger_file_for_mark(mark, norm)
  message("\nReading ", mark, " ", norm, ": ", edger_file)

  e <- read_edger(edger_file)
  edger_gr <- make_granges(e, "chr", "start", "end")

  hits_all <- findOverlaps(edger_gr, tobias_all_gr)
  bg_peak_ids <- unique(queryHits(hits_all))

  bg_peaks <- e[bg_peak_ids]
  bg_gene_map <- annotate_peak_table_to_genes(bg_peaks, window_bp)
  bg_mark <- intersect(unique(bg_gene_map$gene), cm_all)

  hits_cat <- findOverlaps(edger_gr, tobias_cat_gr)

  cat_hit_table <- data.table(
    peak_row_id = queryHits(hits_cat),
    category_label = mcols(tobias_cat_gr)$category_label[subjectHits(hits_cat)],
    tobias_category = mcols(tobias_cat_gr)$tobias_category[subjectHits(hits_cat)]
  ) %>%
    unique(by = c("peak_row_id", "category_label"))

  cat_peaks <- merge(e, cat_hit_table, by = "peak_row_id")
  cat_gene_map <- annotate_peak_table_to_genes(cat_peaks, window_bp)

  tobias_category_gene_sets <- lapply(category_levels, function(category_label) {
    if (!"category_label" %in% colnames(cat_gene_map)) {
      return(character(0))
    }

    intersect(
      unique(cat_gene_map$gene[cat_gene_map$category_label == category_label]),
      bg_mark
    )
  })

  names(tobias_category_gene_sets) <- category_levels

  dn_gene_map <- bg_gene_map %>%
    filter(!is.na(logFC))

  if (use_dn_binding_fdr && "FDR" %in% colnames(dn_gene_map)) {
    dn_gene_map <- dn_gene_map %>%
      filter(!is.na(FDR), FDR < dn_binding_fdr_cutoff)
  }

  dn_category_gene_sets <- list(
    "DN-AP1 down" = intersect(unique(dn_gene_map$gene[dn_gene_map$logFC < 0]), bg_mark),
    "DN-AP1 up" = intersect(unique(dn_gene_map$gene[dn_gene_map$logFC > 0]), bg_mark)
  )

  category_gene_sets <- c(tobias_category_gene_sets, dn_category_gene_sets)

  list(
    bg_mark = bg_mark,
    category_gene_sets = category_gene_sets,
    bg_gene_map = bg_gene_map,
    cat_gene_map = cat_gene_map
  )
}

fisher_enrich <- function(category_genes, deg_genes, bg) {
  bg <- unique(bg)
  category_genes <- intersect(unique(category_genes), bg)
  deg_genes <- intersect(unique(deg_genes), bg)

  a <- length(intersect(category_genes, deg_genes))
  b <- length(setdiff(category_genes, deg_genes))
  c <- length(setdiff(deg_genes, category_genes))
  d <- length(setdiff(bg, union(category_genes, deg_genes)))

  log2_or <- log2(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)))
  log2_or_se <- sqrt(
    1 / (a + 0.5) +
      1 / (b + 0.5) +
      1 / (c + 0.5) +
      1 / (d + 0.5)
  ) / log(2)

  out <- data.frame(
    n_category_genes = length(category_genes),
    n_deg_genes_in_background = length(deg_genes),
    overlap = a,
    background_n = length(bg),
    log2_OR = log2_or,
    log2_OR_low = log2_or - 1.96 * log2_or_se,
    log2_OR_high = log2_or + 1.96 * log2_or_se
  )

  if (length(bg) == 0 || length(category_genes) == 0 || length(deg_genes) == 0) {
    out$OR <- NA_real_
    out$pval <- NA_real_
    return(out)
  }

  mat <- matrix(
    c(a, b, c, d),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      category_status = c("In_category", "Other_background_genes"),
      DEG_status = c("DEG", "Not_DEG")
    )
  )

  ft <- fisher.test(mat)

  out$OR <- as.numeric(ft$estimate)
  out$pval <- ft$p.value
  out
}

mark_maps <- lapply(marks, map_mark_to_tobias_gene_sets, norm = norm, window_bp = window_bp)
names(mark_maps) <- marks

res <- list()

for (mk in marks) {
  bg_mark <- mark_maps[[mk]]$bg_mark
  category_gene_sets <- mark_maps[[mk]]$category_gene_sets

  message(mk, " mark-specific 200 kb background genes: ", length(bg_mark))

  for (category_label in plot_category_levels) {
    for (dir in c("DOWN", "UP")) {
      deg <- if (dir == "DOWN") cm_deg_down else cm_deg_up

      res[[length(res) + 1]] <- cbind(
        mark = mk,
        category_label = category_label,
        DEG_direction = dir,
        fisher_enrich(
          category_genes = category_gene_sets[[category_label]],
          deg_genes = deg,
          bg = bg_mark
        )
      )
    }
  }
}

fisher_df <- bind_rows(res) %>%
  mutate(
    padj = p.adjust(pval, method = "BH"),
    category_label = factor(category_label, levels = plot_category_levels),
    DEG_direction = factor(DEG_direction, levels = c("DOWN", "UP")),
    mark = factor(mark, levels = marks),
    mark_label = recode(
      as.character(mark),
      "C-FOS" = "C-Fos",
      "C-JUN" = "C-Jun",
      "P-C-JUN" = "P-C-Jun"
    ),
    deg_label = recode(
      as.character(DEG_direction),
      "DOWN" = "Down",
      "UP" = "Up"
    ),
    sig_star = case_when(
      !is.na(padj) & padj < 0.001 ~ "***",
      !is.na(padj) & padj < 0.01 ~ "**",
      !is.na(padj) & padj < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

out_csv <- file.path(plotdir, "Fig5b_TOBIAS_plus_DNAP1_category_200kb_enrichment.csv")
write.csv(fisher_df, out_csv, row.names = FALSE)

yr <- max(abs(c(fisher_df$log2_OR_low, fisher_df$log2_OR_high)), na.rm = TRUE)
dodge <- position_dodge(width = 0.75)

plot_df <- fisher_df %>%
  mutate(
    star_y = if_else(
      log2_OR >= 0,
      log2_OR_high + 0.05 * yr,
      log2_OR_low - 0.05 * yr
    )
  )

p <- ggplot(plot_df, aes(category_label, log2_OR, fill = deg_label)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "black") +
  geom_col(position = dodge, width = 0.65, colour = "black", linewidth = 0.25) +
  geom_errorbar(
    aes(ymin = log2_OR_low, ymax = log2_OR_high),
    position = dodge,
    width = 0.18,
    linewidth = 0.3
  ) +
  geom_text(
    aes(y = star_y, label = sig_star),
    position = dodge,
    size = 2.7,
    fontface = "bold"
  ) +
  facet_wrap(~ mark_label, nrow = 1) +
  scale_fill_manual(values = c(
    "Down" = "grey25",
    "Up" = "grey80"
  )) +
  coord_cartesian(ylim = c(-yr * 1.15, yr * 1.15), clip = "off") +
  scale_x_discrete(labels = function(x) gsub(", ", "\n", x, fixed = TRUE)) +
  labs(
    x = NULL,
    y = expression(log[2]~"(odds ratio)"),
    fill = "CM DEG"
  ) +
  theme_classic(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "#ded8cf", colour = NA),
    strip.text = element_text(face = "bold", size = 9),
    panel.spacing.x = unit(1.5, "mm"),
    panel.spacing.y = unit(2.0, "mm"),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.3),
    legend.position = "right",
    plot.margin = margin(4, 8, 4, 4)
  )

out_pdf <- file.path(plotdir, "Fig5b_TOBIAS_plus_DNAP1_category_200kb_barplot.pdf")
ggsave(out_pdf, p, width = 7.4, height = 3.2)

message("Saved: ", out_csv)
message("Saved: ", out_pdf)

fisher_df















#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(ggplot2)
})

plotdir <- publication_path("PAPER_OUTPUT_ROOT", "Paper_Tables/AAV_PLot/")

aav9_mast_file <- file.path(
  publication_path("DNAP1_SNRNA_ROOT", "Process"),
  "MAST_DEG_by_celltype_CDR_batch_removed_male_sample",
  "MAST_ALL_CELLTYPES_AAV9_AP_1_Dosage_2_vs_AAV9_cTnT_all_genes.csv"
)

age_ex_deg_file <- publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed/All_DEG_Age_Ex.rds")
tobias_file <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Output/Tobias_df_Classified_Full.rds")

aav9_cuttag_file <- file.path(
  publication_path("AAV_CUTTAG_EDGER_DIR"),
  "C-FOS",
  "males_spikein",
  "C-FOS.edgeR.AAV9_AP_1_vs_AAV9_Control.full.csv"
)

window_bp <- 200000
aav9_fdr_cutoff <- 0.05
age_ex_padj_cutoff <- 0.05
age_ex_abs_logfc_cutoff <- 0.5
celltype_keep <- "Cardiomyocyte"
aav9_null_seed <- 20260618
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

ap1_tfs <- c("Fos", "Fosb", "Fosl1", "Fosl2", "Jun", "Junb", "Jund")
aav9_deg_levels <- c("AAV9 Down", "AAV9 Null Down-size", "AAV9 Up", "AAV9 Any DEG")
aav9_plot_levels <- c("AAV9 Down", "AAV9 Null Down-size")
age_ex_levels <- c("Age Down", "Age Up", "Age Any DEG", "Exercise Down", "Exercise Up", "Exercise Any DEG")
age_ex_plot_levels <- c("Age Down", "Age Up", "Exercise Down", "Exercise Up")
aav9_plot_labels <- c("AAV9 Down" = "AAV9 Inhibited", "AAV9 Null Down-size" = "AAV9 Null (inhibited-size)")
age_ex_plot_labels <- c("Age Down" = "Age Inhibited", "Age Up" = "Age Up", "Exercise Down" = "Exercise Inhibited", "Exercise Up" = "Exercise Up")
aav9_deg_cols <- c(
  "AAV9 Down" = "grey25",
  "AAV9 Null Down-size" = "#4C78A8",
  "AAV9 Up" = "grey80",
  "AAV9 Any DEG" = "grey55"
)

std_chr <- function(x) {
  x <- as.character(x)
  ifelse(grepl("^chr", x), x, paste0("chr", x))
}

peak_id <- function(df, chr_col, start_col, end_col) {
  paste(std_chr(df[[chr_col]]), as.integer(df[[start_col]]), as.integer(df[[end_col]]), sep = ":")
}

gr_from <- function(df, chr_col, start_col, end_col) {
  GRanges(
    seqnames = std_chr(df[[chr_col]]),
    ranges = IRanges(start = as.integer(df[[start_col]]), end = as.integer(df[[end_col]]))
  )
}

unique_peak_table <- function(df, chr_col = "peak_chr", start_col = "peak_start", end_col = "peak_end") {
  df <- as.data.frame(df)
  if (nrow(df) == 0) {
    return(data.frame(peak_chr = character(0), peak_start = integer(0), peak_end = integer(0), peak_id = character(0)))
  }

  out <- data.frame(
    peak_chr = std_chr(df[[chr_col]]),
    peak_start = as.integer(df[[start_col]]),
    peak_end = as.integer(df[[end_col]])
  )
  out$peak_id <- peak_id(out, "peak_chr", "peak_start", "peak_end")
  out[!duplicated(out$peak_id), ]
}

safe_label <- function(x) {
  out <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_|_$", "", out)
}

p_stars <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

print_fisher_table <- function(mat, title) {
  cat("\n", title, "\n", sep = "")
  cat(paste(rep("=", nchar(title)), collapse = ""), "\n", sep = "")
  print(mat)
}

read_aav9_degs <- function(file) {
  if (!file.exists(file)) {
    stop("AAV9 MAST file was not found:\n  ", file)
  }

  mast <- fread(file)
  needed <- c("gene", "celltype", "mast_coef_logFC", "mast_FDR")
  missing <- setdiff(needed, colnames(mast))
  if (length(missing) > 0) {
    stop("AAV9 MAST file missing columns: ", paste(missing, collapse = ", "))
  }

  mast <- mast[celltype == celltype_keep & !is.na(gene) & gene != ""]
  mast <- mast[!is.na(mast_coef_logFC) & !is.na(mast_FDR)]
  mast <- mast[!duplicated(gene)]

  bg_genes <- unique(mast$gene)
  aav9_down_genes <- mast[mast_FDR < aav9_fdr_cutoff & mast_coef_logFC < 0, unique(gene)]
  set.seed(aav9_null_seed)
  aav9_null_down_size_genes <- sample(bg_genes, length(aav9_down_genes), replace = FALSE)

  deg_sets <- list(
    "AAV9 Down" = aav9_down_genes,
    "AAV9 Null Down-size" = aav9_null_down_size_genes,
    "AAV9 Up" = mast[mast_FDR < aav9_fdr_cutoff & mast_coef_logFC > 0, unique(gene)],
    "AAV9 Any DEG" = mast[mast_FDR < aav9_fdr_cutoff, unique(gene)]
  )

  message("AAV9 CM MAST-tested genes: ", length(bg_genes))
  message("AAV9 Down DEGs: ", length(deg_sets[["AAV9 Down"]]))
  message("AAV9 Null Down-size genes: ", length(deg_sets[["AAV9 Null Down-size"]]), " sampled from tested genes with seed ", aav9_null_seed)
  message("AAV9 Up DEGs: ", length(deg_sets[["AAV9 Up"]]))
  message("AAV9 Any DEGs: ", length(deg_sets[["AAV9 Any DEG"]]))

  list(bg_genes = bg_genes, deg_sets = deg_sets, mast = mast)
}

read_age_ex_degs <- function(file) {
  if (!file.exists(file)) {
    stop("Age/exercise MAST RDS file was not found:\n  ", file)
  }

  mast <- as.data.table(readRDS(file))
  needed <- c("gene", "CellType", "log2FoldChange.Ex", "padj.Ex", "log2FoldChange.Age", "padj.Age")
  missing <- setdiff(needed, colnames(mast))
  if (length(missing) > 0) {
    stop("Age/exercise MAST object missing columns: ", paste(missing, collapse = ", "))
  }

  mast[, gene := as.character(gene)]
  mast[, CellType := trimws(as.character(CellType))]
  mast <- mast[CellType == celltype_keep & !is.na(gene) & gene != ""]
  mast <- mast[!duplicated(gene)]

  age <- mast[!is.na(log2FoldChange.Age) & !is.na(padj.Age)]
  exercise <- mast[!is.na(log2FoldChange.Ex) & !is.na(padj.Ex)]

  out <- list(
    "Age Down" = list(
      contrast = "Age",
      bg_genes = unique(age$gene),
      genes = age[padj.Age < age_ex_padj_cutoff & log2FoldChange.Age < -age_ex_abs_logfc_cutoff, unique(gene)]
    ),
    "Age Up" = list(
      contrast = "Age",
      bg_genes = unique(age$gene),
      genes = age[padj.Age < age_ex_padj_cutoff & log2FoldChange.Age > age_ex_abs_logfc_cutoff, unique(gene)]
    ),
    "Age Any DEG" = list(
      contrast = "Age",
      bg_genes = unique(age$gene),
      genes = age[padj.Age < age_ex_padj_cutoff & abs(log2FoldChange.Age) > age_ex_abs_logfc_cutoff, unique(gene)]
    ),
    "Exercise Down" = list(
      contrast = "Exercise",
      bg_genes = unique(exercise$gene),
      genes = exercise[padj.Ex < age_ex_padj_cutoff & log2FoldChange.Ex < -age_ex_abs_logfc_cutoff, unique(gene)]
    ),
    "Exercise Up" = list(
      contrast = "Exercise",
      bg_genes = unique(exercise$gene),
      genes = exercise[padj.Ex < age_ex_padj_cutoff & log2FoldChange.Ex > age_ex_abs_logfc_cutoff, unique(gene)]
    ),
    "Exercise Any DEG" = list(
      contrast = "Exercise",
      bg_genes = unique(exercise$gene),
      genes = exercise[padj.Ex < age_ex_padj_cutoff & abs(log2FoldChange.Ex) > age_ex_abs_logfc_cutoff, unique(gene)]
    )
  )

  message("\nAge/exercise CM DEG sets:")
  for (nm in names(out)) {
    message("  ", nm, ": ", length(out[[nm]]$genes), " genes; background ", length(out[[nm]]$bg_genes))
  }

  out
}

annotate_peak_gene_map <- function(peaks) {
  if (nrow(peaks) == 0) {
    return(data.frame(peak_id = character(0), gene = character(0), distanceToTSS = numeric(0)))
  }

  gr <- gr_from(peaks, "peak_chr", "peak_start", "peak_end")
  mcols(gr)$peak_id <- peaks$peak_id

  anno <- as.data.frame(
    annotatePeak(
      gr,
      TxDb = txdb,
      annoDb = "org.Mm.eg.db",
      level = "gene"
    )
  )

  anno %>%
    mutate(gene = SYMBOL) %>%
    filter(!is.na(gene), gene != "", !is.na(distanceToTSS)) %>%
    dplyr::select(peak_id, gene, distanceToTSS) %>%
    distinct()
}

genes_within <- function(gene_map, window_bp, background_genes) {
  gene_map %>%
    filter(abs(distanceToTSS) <= window_bp) %>%
    pull(gene) %>%
    unique() %>%
    intersect(background_genes)
}

table_feature_vs_aav9_deg <- function(aav9_deg_genes, bg_genes, feature_genes, feature_label) {
  aav9_deg_genes <- intersect(unique(aav9_deg_genes), bg_genes)
  feature_genes <- intersect(unique(feature_genes), bg_genes)
  non_feature_genes <- setdiff(bg_genes, feature_genes)

  matrix(
    c(
      length(intersect(feature_genes, aav9_deg_genes)),
      length(setdiff(feature_genes, aav9_deg_genes)),
      length(intersect(non_feature_genes, aav9_deg_genes)),
      length(setdiff(non_feature_genes, aav9_deg_genes))
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      Feature_status = c(paste0("In_", safe_label(feature_label)), paste0("Not_in_", safe_label(feature_label))),
      AAV9_DEG_status = c("AAV9_DEG", "Not_AAV9_DEG")
    )
  )
}

run_fisher <- function(mat, plot_name, aav9_deg_set, feature, window_bp = NA_integer_, background_n = NA_integer_, feature_n = NA_integer_, aav9_deg_n = NA_integer_) {
  if (any(rowSums(mat) == 0) || any(colSums(mat) == 0)) {
    or <- ci_low <- ci_high <- pval <- NA_real_
  } else {
    ft <- fisher.test(mat, alternative = "two.sided", conf.int = TRUE)
    or <- as.numeric(ft$estimate)
    ci_low <- ft$conf.int[[1]]
    ci_high <- ft$conf.int[[2]]
    pval <- ft$p.value
  }

  data.frame(
    plot = plot_name,
    aav9_deg_set = aav9_deg_set,
    feature = feature,
    window_bp = window_bp,
    window_kb = ifelse(is.na(window_bp), NA_real_, window_bp / 1000),
    background_n = background_n,
    feature_n = feature_n,
    aav9_deg_n = aav9_deg_n,
    OR = or,
    ci_low = ci_low,
    ci_high = ci_high,
    pval = pval,
    p_star = p_stars(pval),
    fisher_alternative = "two.sided",
    row1_col1 = mat[1, 1],
    row1_col2 = mat[1, 2],
    row2_col1 = mat[2, 1],
    row2_col2 = mat[2, 2],
    row1_name = rownames(mat)[1],
    row2_name = rownames(mat)[2],
    col1_name = colnames(mat)[1],
    col2_name = colnames(mat)[2],
    stringsAsFactors = FALSE
  )
}

prep_log_or_plot <- function(df) {
  finite_log_values <- c(
    log(df$OR[is.finite(df$OR) & df$OR > 0]),
    log(df$ci_low[is.finite(df$ci_low) & df$ci_low > 0]),
    log(df$ci_high[is.finite(df$ci_high) & df$ci_high > 0])
  )
  finite_log_values <- finite_log_values[is.finite(finite_log_values)]
  log_cap <- if (length(finite_log_values) == 0) 1 else max(abs(finite_log_values), 1, na.rm = TRUE) * 1.15

  out <- df %>%
    mutate(
      logOR_plot = case_when(
        is.na(OR) ~ NA_real_,
        OR <= 0 ~ -log_cap,
        is.finite(OR) ~ log(OR),
        TRUE ~ log_cap
      ),
      log_ci_high_plot = case_when(
        is.na(ci_high) ~ NA_real_,
        ci_high <= 0 ~ -log_cap,
        is.finite(ci_high) ~ log(ci_high),
        TRUE ~ log_cap
      ),
      log_ci_low_plot = case_when(
        is.na(ci_low) ~ NA_real_,
        ci_low <= 0 ~ -log_cap,
        is.finite(ci_low) ~ log(ci_low),
        TRUE ~ log_cap
      ),
      aav9_deg_set = factor(aav9_deg_set, levels = aav9_deg_levels),
      feature = factor(feature, levels = c("Any AP-1 TOBIAS motif", "Any AAV9 C-FOS CUT&Tag peak", age_ex_levels))
    )

  out$star_y <- pmax(out$logOR_plot, out$log_ci_high_plot, na.rm = TRUE) + (0.08 * log_cap)
  out$star_y[!is.finite(out$star_y)] <- NA_real_
  out
}

save_plots_one_row <- function(plots, out_pdf, width = 16, height = 3.8) {
  grDevices::pdf(out_pdf, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)

  grid::grid.newpage()
  pushViewport <- grid::pushViewport
  popViewport <- grid::popViewport
  viewport <- grid::viewport
  grid.layout <- grid::grid.layout

  pushViewport(viewport(layout = grid.layout(1, length(plots))))
  for (i in seq_along(plots)) {
    print(plots[[i]], vp = viewport(layout.pos.row = 1, layout.pos.col = i))
  }
  popViewport()
}

plot_aav9_deg_sets <- function(df, title) {
  plot_df <- df %>%
    filter(aav9_deg_set %in% aav9_plot_levels) %>%
    prep_log_or_plot() %>%
    mutate(aav9_deg_set = factor(as.character(aav9_deg_set), levels = aav9_plot_levels))

  p <- ggplot(plot_df, aes(x = aav9_deg_set, y = logOR_plot, fill = aav9_deg_set)) +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "black") +
    geom_col(width = 0.62, colour = "black", linewidth = 0.25, na.rm = TRUE) +
    geom_errorbar(aes(ymin = log_ci_low_plot, ymax = log_ci_high_plot), width = 0.18, linewidth = 0.3, na.rm = TRUE) +
    geom_text(aes(y = star_y, label = p_star), vjust = 0, size = 3.5, na.rm = TRUE) +
    scale_fill_manual(values = aav9_deg_cols, guide = "none") +
    scale_x_discrete(labels = aav9_plot_labels) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.18))) +
    labs(
      title = title,
      subtitle = paste0("AAV9 inhibited DEGs versus same-size random null from AAV9-tested genes; seed ", aav9_null_seed),
      x = NULL,
      y = "log odds ratio"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1)
    )

  p
}

plot_age_ex_overlap <- function(df) {
  plot_df <- df %>%
    filter(feature %in% age_ex_plot_levels, aav9_deg_set %in% aav9_plot_levels) %>%
    prep_log_or_plot() %>%
    mutate(
      aav9_deg_set = factor(as.character(aav9_deg_set), levels = aav9_plot_levels),
      feature = factor(as.character(feature), levels = age_ex_plot_levels)
    )

  p <- ggplot(plot_df, aes(x = feature, y = logOR_plot, fill = aav9_deg_set)) +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "black") +
    geom_col(position = position_dodge(width = 0.72), width = 0.62, colour = "black", linewidth = 0.25, na.rm = TRUE) +
    geom_errorbar(
      aes(ymin = log_ci_low_plot, ymax = log_ci_high_plot, group = aav9_deg_set),
      position = position_dodge(width = 0.72),
      width = 0.18,
      linewidth = 0.3,
      na.rm = TRUE
    ) +
    geom_text(
      aes(y = star_y, label = p_star, group = aav9_deg_set),
      position = position_dodge(width = 0.72),
      vjust = 0,
      size = 3.2,
      na.rm = TRUE
    ) +
    scale_fill_manual(values = aav9_deg_cols, breaks = aav9_plot_levels, labels = aav9_plot_labels, drop = FALSE) +
    scale_x_discrete(labels = age_ex_plot_labels) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.18))) +
    labs(
      title = "AAV9 inhibited/null genes overlapping CM age and exercise DEGs",
      subtitle = paste0("Age/exercise DEGs: padj < ", age_ex_padj_cutoff, ", |logFC| > ", age_ex_abs_logfc_cutoff),
      x = NULL,
      y = "log odds ratio",
      fill = "AAV9 gene set"
    ) +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1)
    )

  p
}

aav9 <- read_aav9_degs(aav9_mast_file)
age_ex <- read_age_ex_degs(age_ex_deg_file)

message("\nReading TOBIAS AP-1 motifs: ", tobias_file)
tobias <- readRDS(tobias_file) %>%
  dplyr::select(peak_chr, peak_start, peak_end, name, CellType) %>%
  filter(CellType == celltype_keep, name %in% ap1_tfs)
setDT(tobias)
ap1_peaks <- unique_peak_table(tobias)
message("CM AP-1 TOBIAS motif peaks: ", nrow(ap1_peaks))

message("\nReading all AAV9 C-FOS CUT&Tag peaks: ", aav9_cuttag_file)
cuttag <- fread(aav9_cuttag_file)
needed_cuttag <- c("chr", "start", "end")
missing_cuttag <- setdiff(needed_cuttag, colnames(cuttag))
if (length(missing_cuttag) > 0) {
  stop("AAV9 CUT&Tag file missing columns: ", paste(missing_cuttag, collapse = ", "))
}
cuttag_peaks <- unique_peak_table(cuttag, chr_col = "chr", start_col = "start", end_col = "end")
message("All AAV9 C-FOS CUT&Tag peaks: ", nrow(cuttag_peaks))

ap1_genes <- genes_within(annotate_peak_gene_map(ap1_peaks), window_bp, aav9$bg_genes)
cuttag_genes <- genes_within(annotate_peak_gene_map(cuttag_peaks), window_bp, aav9$bg_genes)

message("\nGenes within ", window_bp / 1000, " kb of any AP-1 TOBIAS motif: ", length(ap1_genes))
message("Genes within ", window_bp / 1000, " kb of any AAV9 C-FOS CUT&Tag peak: ", length(cuttag_genes))

near_ap1_df <- bind_rows(lapply(names(aav9$deg_sets), function(aav9_deg_set) {
  mat <- table_feature_vs_aav9_deg(aav9$deg_sets[[aav9_deg_set]], aav9$bg_genes, ap1_genes, "Any AP-1 TOBIAS motif")
  print_fisher_table(mat, paste0("Fisher table: ", aav9_deg_set, " near any AP-1 TOBIAS motif"))
  run_fisher(
    mat,
    plot_name = "AAV9 DEG near AP-1 TOBIAS motifs",
    aav9_deg_set = aav9_deg_set,
    feature = "Any AP-1 TOBIAS motif",
    window_bp = window_bp,
    background_n = length(aav9$bg_genes),
    feature_n = length(ap1_genes),
    aav9_deg_n = length(intersect(aav9$deg_sets[[aav9_deg_set]], aav9$bg_genes))
  )
}))

near_cuttag_df <- bind_rows(lapply(names(aav9$deg_sets), function(aav9_deg_set) {
  mat <- table_feature_vs_aav9_deg(aav9$deg_sets[[aav9_deg_set]], aav9$bg_genes, cuttag_genes, "Any AAV9 C-FOS CUT&Tag peak")
  print_fisher_table(mat, paste0("Fisher table: ", aav9_deg_set, " near any AAV9 C-FOS CUT&Tag peak"))
  run_fisher(
    mat,
    plot_name = "AAV9 DEG near AAV9 C-FOS CUT&Tag peaks",
    aav9_deg_set = aav9_deg_set,
    feature = "Any AAV9 C-FOS CUT&Tag peak",
    window_bp = window_bp,
    background_n = length(aav9$bg_genes),
    feature_n = length(cuttag_genes),
    aav9_deg_n = length(intersect(aav9$deg_sets[[aav9_deg_set]], aav9$bg_genes))
  )
}))

overlap_df <- bind_rows(unlist(lapply(names(age_ex), function(age_ex_set) {
  lapply(names(aav9$deg_sets), function(aav9_deg_set) {
    bg <- intersect(aav9$bg_genes, age_ex[[age_ex_set]]$bg_genes)
    feature_genes <- intersect(age_ex[[age_ex_set]]$genes, bg)
    aav9_genes <- intersect(aav9$deg_sets[[aav9_deg_set]], bg)

    mat <- table_feature_vs_aav9_deg(aav9_genes, bg, feature_genes, age_ex_set)
    print_fisher_table(mat, paste0("Fisher table: ", aav9_deg_set, " overlapping ", age_ex_set))
    run_fisher(
      mat,
      plot_name = "AAV9 DEG overlapping age/exercise DEGs",
      aav9_deg_set = aav9_deg_set,
      feature = age_ex_set,
      background_n = length(bg),
      feature_n = length(feature_genes),
      aav9_deg_n = length(aav9_genes)
    )
  })
}), recursive = FALSE))

results_df <- bind_rows(near_ap1_df, near_cuttag_df, overlap_df) %>%
  mutate(padj = p.adjust(pval, method = "BH"), padj_star = p_stars(padj))

out_csv <- file.path(plotdir, "AAV9_DEG_three_enrichment_fisher_tables.csv")
write.csv(results_df, out_csv, row.names = FALSE)

p_ap1 <- plot_aav9_deg_sets(
  near_ap1_df,
  title = "AAV9 inhibited/null genes near AP-1 TOBIAS motifs"
)

p_cuttag <- plot_aav9_deg_sets(
  near_cuttag_df,
  title = "AAV9 inhibited/null genes near AAV9 C-FOS CUT&Tag peaks"
)

p_overlap <- plot_age_ex_overlap(overlap_df)

out_pdf <- file.path(plotdir, "AAV9_inhibited_enrichment_plots_single_row.pdf")
save_plots_one_row(list(p_ap1, p_cuttag, p_overlap), out_pdf, width = 17.5, height = 3.8)

message("Saved: ", out_csv)
message("Saved: ", out_pdf)

p_ap1
p_cuttag
p_overlap

save_plots_one_row(list(p_ap1, p_cuttag, p_overlap), out_pdf, width = 8, height = 3.8)

save_plots_one_row <- function(plots, out_pdf, width = 17.5, height = 3.8) {
  plot_grobs <- lapply(plots, ggplot2::ggplotGrob)

  common_widths <- do.call(
    grid::unit.pmax,
    lapply(plot_grobs, function(x) x$widths)
  )

  plot_grobs <- lapply(plot_grobs, function(x) {
    x$widths <- common_widths
    x
  })

  grDevices::pdf(out_pdf, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(
      nrow = 1,
      ncol = length(plot_grobs),
      widths = grid::unit(rep(1, length(plot_grobs)), "null")
    )
  ))

  for (i in seq_along(plot_grobs)) {
    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = i))
    grid::grid.draw(plot_grobs[[i]])
    grid::popViewport()
  }

  grid::popViewport()
}
      out_pdf <- file.path(plotdir, "AAV9_inhibited_enrichment_plots_single_row.pdf")

save_plots_one_row(
  list(p_ap1, p_cuttag, p_overlap),
  out_pdf,
  width = 9,
  height = 3.8
)

message("Saved: ", out_pdf)



#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(ggplot2)
})

plotdir <- publication_path("PAPER_OUTPUT_ROOT", "Paper_Tables/AAV_PLot/")

aav9_mast_file <- file.path(
  publication_path("DNAP1_SNRNA_ROOT", "Process"),
  "MAST_DEG_by_celltype_CDR_batch_removed_male_sample",
  "MAST_ALL_CELLTYPES_AAV9_AP_1_Dosage_2_vs_AAV9_cTnT_all_genes.csv"
)

age_ex_deg_file <- publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed/All_DEG_Age_Ex.rds")
tobias_file <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Output/Tobias_df_Classified_Full.rds")

aav9_cuttag_file <- file.path(
  publication_path("AAV_CUTTAG_EDGER_DIR"),
  "C-FOS",
  "males_spikein",
  "C-FOS.edgeR.AAV9_AP_1_vs_AAV9_Control.full.csv"
)

window_bp <- 50000
aav9_fdr_cutoff <- 0.05
age_ex_padj_cutoff <- 0.05
age_ex_abs_logfc_cutoff <- 0.5
celltype_keep <- "Cardiomyocyte"
aav9_null_seed <- 20260618
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

ap1_tfs <- c("Fos", "Fosb", "Fosl1", "Fosl2", "Jun", "Junb", "Jund")
aav9_deg_levels <- c("AAV9 Down", "AAV9 Null Down-size", "AAV9 Up", "AAV9 Any DEG")
aav9_plot_levels <- c("AAV9 Down", "AAV9 Null Down-size")
age_ex_levels <- c("Age Down", "Age Up", "Age Any DEG", "Exercise Down", "Exercise Up", "Exercise Any DEG")
age_ex_plot_levels <- c("Age Down", "Age Up", "Exercise Down", "Exercise Up")
aav9_plot_labels <- c("AAV9 Down" = "AAV9 Inhibited", "AAV9 Null Down-size" = "AAV9 Null (inhibited-size)")
age_ex_plot_labels <- c("Age Down" = "Age Inhibited", "Age Up" = "Age Up", "Exercise Down" = "Exercise Inhibited", "Exercise Up" = "Exercise Up")
aav9_deg_cols <- c(
  "AAV9 Down" = "grey25",
  "AAV9 Null Down-size" = "#4C78A8",
  "AAV9 Up" = "grey80",
  "AAV9 Any DEG" = "grey55"
)

std_chr <- function(x) {
  x <- as.character(x)
  ifelse(grepl("^chr", x), x, paste0("chr", x))
}

peak_id <- function(df, chr_col, start_col, end_col) {
  paste(std_chr(df[[chr_col]]), as.integer(df[[start_col]]), as.integer(df[[end_col]]), sep = ":")
}

gr_from <- function(df, chr_col, start_col, end_col) {
  GRanges(
    seqnames = std_chr(df[[chr_col]]),
    ranges = IRanges(start = as.integer(df[[start_col]]), end = as.integer(df[[end_col]]))
  )
}

unique_peak_table <- function(df, chr_col = "peak_chr", start_col = "peak_start", end_col = "peak_end") {
  df <- as.data.frame(df)
  if (nrow(df) == 0) {
    return(data.frame(peak_chr = character(0), peak_start = integer(0), peak_end = integer(0), peak_id = character(0)))
  }

  out <- data.frame(
    peak_chr = std_chr(df[[chr_col]]),
    peak_start = as.integer(df[[start_col]]),
    peak_end = as.integer(df[[end_col]])
  )
  out$peak_id <- peak_id(out, "peak_chr", "peak_start", "peak_end")
  out[!duplicated(out$peak_id), ]
}

safe_label <- function(x) {
  out <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_|_$", "", out)
}

p_stars <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

print_fisher_table <- function(mat, title) {
  cat("\n", title, "\n", sep = "")
  cat(paste(rep("=", nchar(title)), collapse = ""), "\n", sep = "")
  print(mat)
}

read_aav9_degs <- function(file) {
  if (!file.exists(file)) {
    stop("AAV9 MAST file was not found:\n  ", file)
  }

  mast <- fread(file)
  needed <- c("gene", "celltype", "mast_coef_logFC", "mast_FDR")
  missing <- setdiff(needed, colnames(mast))
  if (length(missing) > 0) {
    stop("AAV9 MAST file missing columns: ", paste(missing, collapse = ", "))
  }

  mast <- mast[celltype == celltype_keep & !is.na(gene) & gene != ""]
  mast <- mast[!is.na(mast_coef_logFC) & !is.na(mast_FDR)]
  mast <- mast[!duplicated(gene)]

  bg_genes <- unique(mast$gene)
  aav9_down_genes <- mast[mast_FDR < aav9_fdr_cutoff & mast_coef_logFC < 0, unique(gene)]
  set.seed(aav9_null_seed)
  aav9_null_down_size_genes <- sample(bg_genes, length(aav9_down_genes), replace = FALSE)

  deg_sets <- list(
    "AAV9 Down" = aav9_down_genes,
    "AAV9 Null Down-size" = aav9_null_down_size_genes,
    "AAV9 Up" = mast[mast_FDR < aav9_fdr_cutoff & mast_coef_logFC > 0, unique(gene)],
    "AAV9 Any DEG" = mast[mast_FDR < aav9_fdr_cutoff, unique(gene)]
  )

  message("AAV9 CM MAST-tested genes: ", length(bg_genes))
  message("AAV9 Down DEGs: ", length(deg_sets[["AAV9 Down"]]))
  message("AAV9 Null Down-size genes: ", length(deg_sets[["AAV9 Null Down-size"]]), " sampled from tested genes with seed ", aav9_null_seed)
  message("AAV9 Up DEGs: ", length(deg_sets[["AAV9 Up"]]))
  message("AAV9 Any DEGs: ", length(deg_sets[["AAV9 Any DEG"]]))

  list(bg_genes = bg_genes, deg_sets = deg_sets, mast = mast)
}

read_age_ex_degs <- function(file) {
  if (!file.exists(file)) {
    stop("Age/exercise MAST RDS file was not found:\n  ", file)
  }

  mast <- as.data.table(readRDS(file))
  needed <- c("gene", "CellType", "log2FoldChange.Ex", "padj.Ex", "log2FoldChange.Age", "padj.Age")
  missing <- setdiff(needed, colnames(mast))
  if (length(missing) > 0) {
    stop("Age/exercise MAST object missing columns: ", paste(missing, collapse = ", "))
  }

  mast[, gene := as.character(gene)]
  mast[, CellType := trimws(as.character(CellType))]
  mast <- mast[CellType == celltype_keep & !is.na(gene) & gene != ""]
  mast <- mast[!duplicated(gene)]

  age <- mast[!is.na(log2FoldChange.Age) & !is.na(padj.Age)]
  exercise <- mast[!is.na(log2FoldChange.Ex) & !is.na(padj.Ex)]

  out <- list(
    "Age Down" = list(
      contrast = "Age",
      bg_genes = unique(age$gene),
      genes = age[padj.Age < age_ex_padj_cutoff & log2FoldChange.Age < -age_ex_abs_logfc_cutoff, unique(gene)]
    ),
    "Age Up" = list(
      contrast = "Age",
      bg_genes = unique(age$gene),
      genes = age[padj.Age < age_ex_padj_cutoff & log2FoldChange.Age > age_ex_abs_logfc_cutoff, unique(gene)]
    ),
    "Age Any DEG" = list(
      contrast = "Age",
      bg_genes = unique(age$gene),
      genes = age[padj.Age < age_ex_padj_cutoff & abs(log2FoldChange.Age) > age_ex_abs_logfc_cutoff, unique(gene)]
    ),
    "Exercise Down" = list(
      contrast = "Exercise",
      bg_genes = unique(exercise$gene),
      genes = exercise[padj.Ex < age_ex_padj_cutoff & log2FoldChange.Ex < -age_ex_abs_logfc_cutoff, unique(gene)]
    ),
    "Exercise Up" = list(
      contrast = "Exercise",
      bg_genes = unique(exercise$gene),
      genes = exercise[padj.Ex < age_ex_padj_cutoff & log2FoldChange.Ex > age_ex_abs_logfc_cutoff, unique(gene)]
    ),
    "Exercise Any DEG" = list(
      contrast = "Exercise",
      bg_genes = unique(exercise$gene),
      genes = exercise[padj.Ex < age_ex_padj_cutoff & abs(log2FoldChange.Ex) > age_ex_abs_logfc_cutoff, unique(gene)]
    )
  )

  message("\nAge/exercise CM DEG sets:")
  for (nm in names(out)) {
    message("  ", nm, ": ", length(out[[nm]]$genes), " genes; background ", length(out[[nm]]$bg_genes))
  }

  out
}

annotate_peak_gene_map <- function(peaks) {
  if (nrow(peaks) == 0) {
    return(data.frame(peak_id = character(0), gene = character(0), distanceToTSS = numeric(0)))
  }

  gr <- gr_from(peaks, "peak_chr", "peak_start", "peak_end")
  mcols(gr)$peak_id <- peaks$peak_id

  anno <- as.data.frame(
    annotatePeak(
      gr,
      TxDb = txdb,
      annoDb = "org.Mm.eg.db",
      level = "gene"
    )
  )

  anno %>%
    mutate(gene = SYMBOL) %>%
    filter(!is.na(gene), gene != "", !is.na(distanceToTSS)) %>%
    dplyr::select(peak_id, gene, distanceToTSS) %>%
    distinct()
}

genes_within <- function(gene_map, window_bp, background_genes) {
  gene_map %>%
    filter(abs(distanceToTSS) <= window_bp) %>%
    pull(gene) %>%
    unique() %>%
    intersect(background_genes)
}

table_feature_vs_aav9_deg <- function(aav9_deg_genes, bg_genes, feature_genes, feature_label) {
  aav9_deg_genes <- intersect(unique(aav9_deg_genes), bg_genes)
  feature_genes <- intersect(unique(feature_genes), bg_genes)
  non_feature_genes <- setdiff(bg_genes, feature_genes)

  matrix(
    c(
      length(intersect(feature_genes, aav9_deg_genes)),
      length(setdiff(feature_genes, aav9_deg_genes)),
      length(intersect(non_feature_genes, aav9_deg_genes)),
      length(setdiff(non_feature_genes, aav9_deg_genes))
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      Feature_status = c(paste0("In_", safe_label(feature_label)), paste0("Not_in_", safe_label(feature_label))),
      AAV9_DEG_status = c("AAV9_DEG", "Not_AAV9_DEG")
    )
  )
}

run_fisher <- function(mat, plot_name, aav9_deg_set, feature, window_bp = NA_integer_, background_n = NA_integer_, feature_n = NA_integer_, aav9_deg_n = NA_integer_) {
  if (any(rowSums(mat) == 0) || any(colSums(mat) == 0)) {
    or <- ci_low <- ci_high <- pval <- NA_real_
  } else {
    ft <- fisher.test(mat, alternative = "two.sided", conf.int = TRUE)
    or <- as.numeric(ft$estimate)
    ci_low <- ft$conf.int[[1]]
    ci_high <- ft$conf.int[[2]]
    pval <- ft$p.value
  }

  data.frame(
    plot = plot_name,
    aav9_deg_set = aav9_deg_set,
    feature = feature,
    window_bp = window_bp,
    window_kb = ifelse(is.na(window_bp), NA_real_, window_bp / 1000),
    background_n = background_n,
    feature_n = feature_n,
    aav9_deg_n = aav9_deg_n,
    OR = or,
    ci_low = ci_low,
    ci_high = ci_high,
    pval = pval,
    p_star = p_stars(pval),
    fisher_alternative = "two.sided",
    row1_col1 = mat[1, 1],
    row1_col2 = mat[1, 2],
    row2_col1 = mat[2, 1],
    row2_col2 = mat[2, 2],
    row1_name = rownames(mat)[1],
    row2_name = rownames(mat)[2],
    col1_name = colnames(mat)[1],
    col2_name = colnames(mat)[2],
    stringsAsFactors = FALSE
  )
}

prep_log_or_plot <- function(df) {
  finite_log_values <- c(
    log(df$OR[is.finite(df$OR) & df$OR > 0]),
    log(df$ci_low[is.finite(df$ci_low) & df$ci_low > 0]),
    log(df$ci_high[is.finite(df$ci_high) & df$ci_high > 0])
  )
  finite_log_values <- finite_log_values[is.finite(finite_log_values)]
  log_cap <- if (length(finite_log_values) == 0) 1 else max(abs(finite_log_values), 1, na.rm = TRUE) * 1.15

  out <- df %>%
    mutate(
      logOR_plot = case_when(
        is.na(OR) ~ NA_real_,
        OR <= 0 ~ -log_cap,
        is.finite(OR) ~ log(OR),
        TRUE ~ log_cap
      ),
      log_ci_high_plot = case_when(
        is.na(ci_high) ~ NA_real_,
        ci_high <= 0 ~ -log_cap,
        is.finite(ci_high) ~ log(ci_high),
        TRUE ~ log_cap
      ),
      log_ci_low_plot = case_when(
        is.na(ci_low) ~ NA_real_,
        ci_low <= 0 ~ -log_cap,
        is.finite(ci_low) ~ log(ci_low),
        TRUE ~ log_cap
      ),
      aav9_deg_set = factor(aav9_deg_set, levels = aav9_deg_levels),
      feature = factor(feature, levels = c("Any AP-1 TOBIAS motif", "Any AAV9 C-FOS CUT&Tag peak", age_ex_levels))
    )

  out$star_y <- pmax(out$logOR_plot, out$log_ci_high_plot, na.rm = TRUE) + (0.08 * log_cap)
  out$star_y[!is.finite(out$star_y)] <- NA_real_
  out
}

save_plots_one_row <- function(plots, out_pdf, width = 16, height = 3.8) {
  grDevices::pdf(out_pdf, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)

  grid::grid.newpage()
  pushViewport <- grid::pushViewport
  popViewport <- grid::popViewport
  viewport <- grid::viewport
  grid.layout <- grid::grid.layout

  pushViewport(viewport(layout = grid.layout(1, length(plots))))
  for (i in seq_along(plots)) {
    print(plots[[i]], vp = viewport(layout.pos.row = 1, layout.pos.col = i))
  }
  popViewport()
}

plot_aav9_deg_sets <- function(df, title) {
  plot_df <- df %>%
    filter(aav9_deg_set %in% aav9_plot_levels) %>%
    prep_log_or_plot() %>%
    mutate(aav9_deg_set = factor(as.character(aav9_deg_set), levels = aav9_plot_levels))

  p <- ggplot(plot_df, aes(x = aav9_deg_set, y = logOR_plot, fill = aav9_deg_set)) +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "black") +
    geom_col(width = 0.62, colour = "black", linewidth = 0.25, na.rm = TRUE) +
    geom_errorbar(aes(ymin = log_ci_low_plot, ymax = log_ci_high_plot), width = 0.18, linewidth = 0.3, na.rm = TRUE) +
    geom_text(aes(y = star_y, label = p_star), vjust = 0, size = 3.5, na.rm = TRUE) +
    scale_fill_manual(values = aav9_deg_cols, guide = "none") +
    scale_x_discrete(labels = aav9_plot_labels) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.18))) +
    labs(
      title = title,
      subtitle = paste0("AAV9 inhibited DEGs versus same-size random null from AAV9-tested genes; seed ", aav9_null_seed),
      x = NULL,
      y = "log odds ratio"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1)
    )

  p
}

plot_age_ex_overlap <- function(df) {
  plot_df <- df %>%
    filter(feature %in% age_ex_plot_levels, aav9_deg_set %in% aav9_plot_levels) %>%
    prep_log_or_plot() %>%
    mutate(
      aav9_deg_set = factor(as.character(aav9_deg_set), levels = aav9_plot_levels),
      feature = factor(as.character(feature), levels = age_ex_plot_levels)
    )

  p <- ggplot(plot_df, aes(x = feature, y = logOR_plot, fill = aav9_deg_set)) +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "black") +
    geom_col(position = position_dodge(width = 0.72), width = 0.62, colour = "black", linewidth = 0.25, na.rm = TRUE) +
    geom_errorbar(
      aes(ymin = log_ci_low_plot, ymax = log_ci_high_plot, group = aav9_deg_set),
      position = position_dodge(width = 0.72),
      width = 0.18,
      linewidth = 0.3,
      na.rm = TRUE
    ) +
    geom_text(
      aes(y = star_y, label = p_star, group = aav9_deg_set),
      position = position_dodge(width = 0.72),
      vjust = 0,
      size = 3.2,
      na.rm = TRUE
    ) +
    scale_fill_manual(values = aav9_deg_cols, breaks = aav9_plot_levels, labels = aav9_plot_labels, drop = FALSE) +
    scale_x_discrete(labels = age_ex_plot_labels) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.18))) +
    labs(
      title = "AAV9 inhibited/null genes overlapping CM age and exercise DEGs",
      subtitle = paste0("Age/exercise DEGs: padj < ", age_ex_padj_cutoff, ", |logFC| > ", age_ex_abs_logfc_cutoff),
      x = NULL,
      y = "log odds ratio",
      fill = "AAV9 gene set"
    ) +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1)
    )

  p
}

aav9 <- read_aav9_degs(aav9_mast_file)
age_ex <- read_age_ex_degs(age_ex_deg_file)

message("\nReading TOBIAS AP-1 motifs: ", tobias_file)
tobias <- readRDS(tobias_file) %>%
  dplyr::select(peak_chr, peak_start, peak_end, name, CellType) %>%
  filter(CellType == celltype_keep, name %in% ap1_tfs)
setDT(tobias)
ap1_peaks <- unique_peak_table(tobias)
message("CM AP-1 TOBIAS motif peaks: ", nrow(ap1_peaks))

message("\nReading all AAV9 C-FOS CUT&Tag peaks: ", aav9_cuttag_file)
cuttag <- fread(aav9_cuttag_file)
needed_cuttag <- c("chr", "start", "end")
missing_cuttag <- setdiff(needed_cuttag, colnames(cuttag))
if (length(missing_cuttag) > 0) {
  stop("AAV9 CUT&Tag file missing columns: ", paste(missing_cuttag, collapse = ", "))
}
cuttag_peaks <- unique_peak_table(cuttag, chr_col = "chr", start_col = "start", end_col = "end")
message("All AAV9 C-FOS CUT&Tag peaks: ", nrow(cuttag_peaks))

ap1_genes <- genes_within(annotate_peak_gene_map(ap1_peaks), window_bp, aav9$bg_genes)
cuttag_genes <- genes_within(annotate_peak_gene_map(cuttag_peaks), window_bp, aav9$bg_genes)

message("\nGenes within ", window_bp / 1000, " kb of any AP-1 TOBIAS motif: ", length(ap1_genes))
message("Genes within ", window_bp / 1000, " kb of any AAV9 C-FOS CUT&Tag peak: ", length(cuttag_genes))

near_ap1_df <- bind_rows(lapply(names(aav9$deg_sets), function(aav9_deg_set) {
  mat <- table_feature_vs_aav9_deg(aav9$deg_sets[[aav9_deg_set]], aav9$bg_genes, ap1_genes, "Any AP-1 TOBIAS motif")
  print_fisher_table(mat, paste0("Fisher table: ", aav9_deg_set, " near any AP-1 TOBIAS motif"))
  run_fisher(
    mat,
    plot_name = "AAV9 DEG near AP-1 TOBIAS motifs",
    aav9_deg_set = aav9_deg_set,
    feature = "Any AP-1 TOBIAS motif",
    window_bp = window_bp,
    background_n = length(aav9$bg_genes),
    feature_n = length(ap1_genes),
    aav9_deg_n = length(intersect(aav9$deg_sets[[aav9_deg_set]], aav9$bg_genes))
  )
}))

near_cuttag_df <- bind_rows(lapply(names(aav9$deg_sets), function(aav9_deg_set) {
  mat <- table_feature_vs_aav9_deg(aav9$deg_sets[[aav9_deg_set]], aav9$bg_genes, cuttag_genes, "Any AAV9 C-FOS CUT&Tag peak")
  print_fisher_table(mat, paste0("Fisher table: ", aav9_deg_set, " near any AAV9 C-FOS CUT&Tag peak"))
  run_fisher(
    mat,
    plot_name = "AAV9 DEG near AAV9 C-FOS CUT&Tag peaks",
    aav9_deg_set = aav9_deg_set,
    feature = "Any AAV9 C-FOS CUT&Tag peak",
    window_bp = window_bp,
    background_n = length(aav9$bg_genes),
    feature_n = length(cuttag_genes),
    aav9_deg_n = length(intersect(aav9$deg_sets[[aav9_deg_set]], aav9$bg_genes))
  )
}))

overlap_df <- bind_rows(unlist(lapply(names(age_ex), function(age_ex_set) {
  lapply(names(aav9$deg_sets), function(aav9_deg_set) {
    bg <- intersect(aav9$bg_genes, age_ex[[age_ex_set]]$bg_genes)
    feature_genes <- intersect(age_ex[[age_ex_set]]$genes, bg)
    aav9_genes <- intersect(aav9$deg_sets[[aav9_deg_set]], bg)

    mat <- table_feature_vs_aav9_deg(aav9_genes, bg, feature_genes, age_ex_set)
    print_fisher_table(mat, paste0("Fisher table: ", aav9_deg_set, " overlapping ", age_ex_set))
    run_fisher(
      mat,
      plot_name = "AAV9 DEG overlapping age/exercise DEGs",
      aav9_deg_set = aav9_deg_set,
      feature = age_ex_set,
      background_n = length(bg),
      feature_n = length(feature_genes),
      aav9_deg_n = length(aav9_genes)
    )
  })
}), recursive = FALSE))

results_df <- bind_rows(near_ap1_df, near_cuttag_df, overlap_df) %>%
  mutate(padj = p.adjust(pval, method = "BH"), padj_star = p_stars(padj))

out_csv <- file.path(plotdir, "AAV9_DEG_three_enrichment_fisher_tables.csv")
write.csv(results_df, out_csv, row.names = FALSE)

p_ap1 <- plot_aav9_deg_sets(
  near_ap1_df,
  title = "AAV9 inhibited/null genes near AP-1 TOBIAS motifs"
)

p_cuttag <- plot_aav9_deg_sets(
  near_cuttag_df,
  title = "AAV9 inhibited/null genes near AAV9 C-FOS CUT&Tag peaks"
)

p_overlap <- plot_age_ex_overlap(overlap_df)

out_pdf <- file.path(plotdir, "AAV9_inhibited_enrichment_plots_single_row.pdf")
save_plots_one_row(list(p_ap1, p_cuttag, p_overlap), out_pdf, width = 17.5, height = 3.8)

message("Saved: ", out_csv)
message("Saved: ", out_pdf)

p_ap1
p_cuttag
p_overlap
save_plots_one_row <- function(plots, out_pdf, width = 17.5, height = 3.8) {
  plot_grobs <- lapply(plots, ggplot2::ggplotGrob)

  common_widths <- do.call(
    grid::unit.pmax,
    lapply(plot_grobs, function(x) x$widths)
  )

  plot_grobs <- lapply(plot_grobs, function(x) {
    x$widths <- common_widths
    x
  })

  grDevices::pdf(out_pdf, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(
      nrow = 1,
      ncol = length(plot_grobs),
      widths = grid::unit(rep(1, length(plot_grobs)), "null")
    )
  ))

  for (i in seq_along(plot_grobs)) {
    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = i))
    grid::grid.draw(plot_grobs[[i]])
    grid::popViewport()
  }

  grid::popViewport()
}
      out_pdf <- file.path(plotdir, "AAV9_inhibited_enrichment_plots_single_row_20kb.pdf")

save_plots_one_row(
  list(p_ap1, p_cuttag, p_overlap),
  out_pdf,
  width = 9,
  height = 3.8
)

message("Saved: ", out_pdf)

