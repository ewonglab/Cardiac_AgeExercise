#!/usr/bin/env Rscript
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Compare 10-week c-FOS CUT&Tag reversal with TOBIAS AP-1 footprints and pathway links.
# Steps: 1) Load CUT&Tag and TOBIAS results; 2) compute genomic overlaps and enrichment; 3) export overlap summaries and figures.
# =======================================
# Source reformatted from emily_plots/cut&tag/contains_pathway_plot.R.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))
# ==============================================================================
# venn_cuttag_vs_tobias_ap1.R
# Overlap analysis: C-FOS CUT&Tag peaks vs fibroblast AP-1 TOBIAS footprinting
# Includes reversed-vs-reversed overlap, Venn diagrams, Fisher's exact tests,
# and pathway enrichment of reversed CUT&Tag sites
# ==============================================================================

library(data.table)
library(GenomicRanges)
library(VennDiagram)
library(grid)
library(ggplot2)

cuttag_dir <- publication_path("CUTTAG_10_WEEK_ANALYSIS_ROOT")
p2p_dir    <- publication_path("PEAK2PATHWAY_ROOT")

fig_dir <- file.path(cuttag_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# ==============================================================================
# 1. Load data
# ==============================================================================

cat("=" , strrep("=", 69), "\n")
cat("REVERSAL DEFINITIONS\n")
cat(strrep("=", 70), "\n\n")

cat("CUT&Tag reversal:\n")
cat("  Peaks where C-FOS binding INCREASED with aging (logFC_OS_vs_YS > 0)\n")
cat("  and DECREASED with exercise (logFC_OE_vs_OS < 0).\n")
cat("  'Directional': both directions present regardless of P-value\n")
cat("  'One-sig':     + at least one comparison P < 0.05\n")
cat("  'Both-sig':    + both comparisons P < 0.05\n\n")

cat("Footprint reversal:\n")
cat("  AP-1 sites that are BOUND in aged sedentary (in_MAS) but NOT in\n")
cat("  young sedentary (not in_YS) [= age-gained], AND then LOST with\n")
cat("  exercise (not in_MAE) [= exercise-reversed].\n")
cat("  i.e., in_MAS & !in_YS & !in_MAE\n\n")

cat("Loading CUT&Tag consensus peaks...\n")
peaks_bed <- fread(file.path(cuttag_dir, "consensus_peaks_C-FOS_support_all3.bed"),
                   col.names = c("chr", "start", "end", "peak_id", "score", "strand"))
cat("  C-FOS CUT&Tag peaks:", nrow(peaks_bed), "\n")

# Load reversal info and merge
reversal <- fread(file.path(cuttag_dir, "reversal_C-FOS_support_all3_all_peaks.csv"))
peaks_bed <- merge(peaks_bed, reversal[, .(peak_id, logFC_aging, PValue_aging,
                                            logFC_exer, PValue_exer, reversal_index)],
                   by = "peak_id", all.x = TRUE)

# Define reversed CUT&Tag sites
peaks_bed$reversed_directional <- with(peaks_bed,
  !is.na(logFC_aging) & logFC_aging > 0 & logFC_exer < 0)
peaks_bed$reversed_one_sig <- with(peaks_bed,
  !is.na(logFC_aging) & logFC_aging > 0 & logFC_exer < 0 &
  (PValue_aging < 0.05 | PValue_exer < 0.05))
peaks_bed$reversed_nominal <- with(peaks_bed,
  !is.na(logFC_aging) & logFC_aging > 0 & PValue_aging < 0.05 &
  logFC_exer < 0 & PValue_exer < 0.05)

cat("\nCUT&Tag reversal counts:\n")
cat("  Directional (aging UP & exer DOWN):     ", sum(peaks_bed$reversed_directional), "\n")
cat("  One-sig (+ at least one P<0.05):        ", sum(peaks_bed$reversed_one_sig), "\n")
cat("  Both-sig (aging P<0.05 & exer P<0.05):  ", sum(peaks_bed$reversed_nominal), "\n")

# Load TOBIAS footprinting
cat("\nLoading TOBIAS fibroblast AP-1 footprinting sites...\n")
tob <- readRDS(file.path(p2p_dir, "Tobias_df_Classified_Full.rds"))
ap1_members <- c("Fos", "Fosb", "Fosl1", "Fosl2", "Jun", "Junb", "Jund")
fib_ap1 <- tob[tob$CellType == "Fibroblast" & tob$name %in% ap1_members, ]
cat("  Fibroblast AP-1 sites (all):", nrow(fib_ap1), "\n")

fib_ap1_unique <- unique(fib_ap1[, c("seqnames", "start", "end")])
cat("  Unique genomic locations:", nrow(fib_ap1_unique), "\n")

# Age-gained: bound in aged-sed but not young-sed
fib_ag <- fib_ap1[fib_ap1$in_MAS & !fib_ap1$in_YS, ]
fib_ag_unique <- unique(fib_ag[, c("seqnames", "start", "end")])
cat("  Age-gained (in_MAS & !in_YS):", nrow(fib_ag), "(", nrow(fib_ag_unique), "unique loci)\n")

# Reversed footprint: age-gained AND lost with exercise
fib_rev <- fib_ap1[fib_ap1$in_MAS & !fib_ap1$in_YS & !fib_ap1$in_MAE, ]
fib_rev_unique <- unique(fib_rev[, c("seqnames", "start", "end")])
cat("  Reversed (in_MAS & !in_YS & !in_MAE):", nrow(fib_rev), "(", nrow(fib_rev_unique), "unique loci)\n")

# ==============================================================================
# 2. Build GRanges
# ==============================================================================

cuttag_gr <- GRanges(seqnames = peaks_bed$chr,
                      ranges = IRanges(start = peaks_bed$start, end = peaks_bed$end),
                      peak_id = peaks_bed$peak_id)

fp_all_gr <- GRanges(seqnames = fib_ap1_unique$seqnames,
                      ranges = IRanges(start = fib_ap1_unique$start,
                                       end = fib_ap1_unique$end))
fp_ag_gr <- GRanges(seqnames = fib_ag_unique$seqnames,
                     ranges = IRanges(start = fib_ag_unique$start,
                                      end = fib_ag_unique$end))
fp_rev_gr <- GRanges(seqnames = fib_rev_unique$seqnames,
                      ranges = IRanges(start = fib_rev_unique$start,
                                       end = fib_rev_unique$end))

cuttag_rev_gr <- cuttag_gr[cuttag_gr$peak_id %in%
  peaks_bed$peak_id[peaks_bed$reversed_directional]]
cuttag_rev_sig_gr <- cuttag_gr[cuttag_gr$peak_id %in%
  peaks_bed$peak_id[peaks_bed$reversed_one_sig]]

# ==============================================================================
# 3. All overlaps
# ==============================================================================

overlap_report <- function(name1, gr1, name2, gr2) {
  hits <- findOverlaps(gr1, gr2)
  n1 <- length(unique(queryHits(hits)))
  n2 <- length(unique(subjectHits(hits)))
  cat(sprintf("  %s in %s: %d/%d (%.1f%%)\n", name1, name2, n1, length(gr1),
              100 * n1 / length(gr1)))
  cat(sprintf("  %s in %s: %d/%d (%.1f%%)\n", name2, name1, n2, length(gr2),
              100 * n2 / length(gr2)))
  list(hits = hits, n_q = n1, n_s = n2)
}

cat("\n", strrep("=", 70), "\n")
cat("OVERLAP ANALYSIS\n")
cat(strrep("=", 70), "\n\n")

cat("--- All C-FOS CUT&Tag vs All Fib AP-1 ---\n")
ov1 <- overlap_report("CUT&Tag", cuttag_gr, "FP", fp_all_gr)

cat("\n--- All C-FOS CUT&Tag vs Age-Gained Fib AP-1 ---\n")
ov2 <- overlap_report("CUT&Tag", cuttag_gr, "AgeGained_FP", fp_ag_gr)

cat("\n--- Reversed CUT&Tag vs All Fib AP-1 ---\n")
ov3 <- overlap_report("Rev_CUT&Tag", cuttag_rev_gr, "FP", fp_all_gr)

cat("\n--- Reversed CUT&Tag vs Age-Gained Fib AP-1 ---\n")
ov4 <- overlap_report("Rev_CUT&Tag", cuttag_rev_gr, "AgeGained_FP", fp_ag_gr)

cat("\n--- Reversed CUT&Tag vs Reversed Fib AP-1 ---\n")
ov5 <- overlap_report("Rev_CUT&Tag", cuttag_rev_gr, "Rev_FP", fp_rev_gr)

cat("\n--- All CUT&Tag vs Reversed Fib AP-1 ---\n")
ov6 <- overlap_report("CUT&Tag", cuttag_gr, "Rev_FP", fp_rev_gr)

cat("\n--- One-sig Reversed CUT&Tag vs Reversed Fib AP-1 ---\n")
ov7 <- overlap_report("RevSig_CUT&Tag", cuttag_rev_sig_gr, "Rev_FP", fp_rev_gr)

# ==============================================================================
# 4. Fisher's exact tests
# ==============================================================================

cat("\n", strrep("=", 70), "\n")
cat("FISHER'S EXACT TESTS\n")
cat(strrep("=", 70), "\n\n")

peaks_bed$overlaps_fp_all <- peaks_bed$peak_id %in%
  cuttag_gr$peak_id[unique(queryHits(ov1$hits))]
peaks_bed$overlaps_fp_ag <- peaks_bed$peak_id %in%
  cuttag_gr$peak_id[unique(queryHits(ov2$hits))]
peaks_bed$overlaps_fp_rev <- peaks_bed$peak_id %in%
  cuttag_gr$peak_id[unique(queryHits(ov6$hits))]

run_fisher <- function(label, pred_col, resp_col) {
  tab <- table(factor(peaks_bed[[pred_col]], levels = c(FALSE, TRUE)),
               factor(peaks_bed[[resp_col]], levels = c(FALSE, TRUE)))
  dimnames(tab) <- list(pred_col = c("FALSE", "TRUE"),
                         resp_col = c("FALSE", "TRUE"))
  ft <- fisher.test(tab)
  cat(label, "\n")
  print(tab)
  cat(sprintf("  OR = %.3f, 95%% CI [%.3f, %.3f], P = %s\n\n",
              ft$estimate, ft$conf.int[1], ft$conf.int[2],
              format.pval(ft$p.value, digits = 3)))
  ft
}

ft1 <- run_fisher("Test 1: Dir-reversed CUT&Tag enriched for ALL AP-1 FP overlap?",
                   "reversed_directional", "overlaps_fp_all")
ft2 <- run_fisher("Test 2: Dir-reversed CUT&Tag enriched for age-gained AP-1 overlap?",
                   "reversed_directional", "overlaps_fp_ag")
ft3 <- run_fisher("Test 3: Dir-reversed CUT&Tag enriched for REVERSED AP-1 FP overlap?",
                   "reversed_directional", "overlaps_fp_rev")
ft4 <- run_fisher("Test 4: One-sig reversed CUT&Tag enriched for ALL AP-1 FP overlap?",
                   "reversed_one_sig", "overlaps_fp_all")
ft5 <- run_fisher("Test 5: One-sig reversed CUT&Tag enriched for REVERSED AP-1 FP overlap?",
                   "reversed_one_sig", "overlaps_fp_rev")

# ==============================================================================
# 5. Pathway enrichment of reversed CUT&Tag peaks
# ==============================================================================

cat("\n", strrep("=", 70), "\n")
cat("PATHWAY ENRICHMENT: REVERSED CUT&TAG PEAKS\n")
cat(strrep("=", 70), "\n\n")

# Gene modules from peak2pathway (01_peak2module_enrichment.R)
gene_modules <- list(
  Angiogenesis = c(
    "Vegfa", "Vegfb", "Vegfc", "Flt1", "Kdr", "Flt4", "Angpt1", "Angpt2",
    "Tek", "Pecam1", "Cdh5", "Nos3", "Eng", "Dll4", "Notch1", "Ednra",
    "Ednrb", "Edn1", "Ace", "Ace2", "Hif1a", "Epas1", "Nrp1", "Nrp2",
    "Pdgfb", "Pdgfrb", "Pgf"
  ),
  Inflammation = c(
    "Tnf", "Il6", "Il1b", "Il1a", "Il18", "Ccl2", "Ccl3", "Ccl4", "Ccl5",
    "Ccl7", "Cxcl1", "Cxcl2", "Cxcl10", "Cxcl12", "Ccr2", "Cxcr4",
    "Nfkb1", "Nfkb2", "Rela", "Ikbkb", "Nlrp3", "Pycard", "Casp1",
    "Ptgs2", "Nos2", "Mmp9", "Il10", "Tgfb1", "Tlr2", "Tlr4"
  ),
  ECM_Remodeling = c(
    "Col1a1", "Col1a2", "Col3a1", "Col4a1", "Col6a1", "Col6a3", "Fn1",
    "Lama2", "Lama4", "Mmp2", "Mmp3", "Mmp9", "Mmp14", "Timp1", "Timp2",
    "Timp3", "Lox", "Loxl2", "Dcn", "Bgn", "Vcan", "Postn", "Thbs1",
    "Thbs2", "Sparc", "Tnc", "Tgfb1", "Tgfb2", "Tgfb3", "Tgfbr1",
    "Tgfbr2", "Adamts2"
  ),
  Oxidative_Stress = c(
    "Sod1", "Sod2", "Sod3", "Cat", "Gpx1", "Gpx4", "Gsr", "Gclc", "Gclm",
    "Gstm1", "Txn1", "Txnrd1", "Prdx1", "Prdx2", "Nfe2l2", "Keap1",
    "Hmox1", "Nqo1", "Cybb", "Cyba", "Nox4", "Hif1a", "Epas1", "Egln1",
    "Egln3", "Bnip3", "Gadd45a", "Gadd45b", "Ldha", "Pdk1"
  ),
  Mitochondrial = c(
    "Ppargc1a", "Ppargc1b", "Tfam", "Nrf1", "Ndufs1", "Ndufv1", "Sdha",
    "Sdhb", "Uqcrc1", "Cox5a", "Cox7a1", "Atp5f1a", "Dnm1l", "Fis1",
    "Mff", "Mfn1", "Mfn2", "Opa1", "Pink1", "Prkn", "Bnip3l", "Fundc1",
    "Hspd1", "Hspe1", "Lonp1", "Clpp", "Cpt1a", "Cpt1b", "Cpt2",
    "Acadl", "Acadm", "Cs", "Idh2"
  ),
  Mechanotransduction = c(
    "Itgb1", "Itgb3", "Itga5", "Itga7", "Itgav", "Vcl", "Tln1", "Pxn",
    "Ptk2", "Yap1", "Wwtr1", "Tead1", "Tead4", "Lats1", "Lats2",
    "Piezo1", "Piezo2", "Trpv4", "Rock1", "Rock2", "Rhoa", "Dmd",
    "Dag1", "Sgcd", "Ttn", "Cd44", "Sdc1", "Sdc4", "Actn1", "Flna"
  ),
  IL1_Pathway = c(
    "Il1a", "Il1b", "Il18", "Il33", "Il36a", "Il36b", "Il36g", "Il36rn",
    "Il1r1", "Il1r2", "Il1rap", "Il1rl1", "Il1rl2", "Il18r1", "Il18rap",
    "Il1rn", "Myd88", "Irak1", "Irak2", "Irak4", "Traf6", "Map3k7",
    "Mapk14", "Mapk8", "Nfkb1", "Rela", "Fos", "Jun", "Junb", "Nlrp3",
    "Pycard", "Casp1", "Casp4", "Ptgs2", "Ccl2", "Cxcl1", "Cxcl2",
    "Mmp3", "Mmp9", "Icam1", "Vcam1", "Sele"
  ),
  Cardiac_Hypertrophy = c(
    "Nppa", "Nppb", "Myh7", "Acta1", "Mybpc3", "Gata4", "Nkx2-5",
    "Mef2a", "Mef2c", "Hand2", "Tbx5", "Myocd", "Srf", "Ankrd1", "Ctgf",
    "Tgfb1"
  ),
  Senescence = c(
    "Cdkn1a", "Cdkn2a", "Tp53", "Rb1", "Serpine1", "Igfbp3", "Igfbp5",
    "Gadd45a", "Gadd45b", "Bmi1", "Hmga2", "Trp53bp1", "Pml", "Ccnd1"
  ),
  Calcium_Handling = c(
    "Ryr2", "Atp2a2", "Pln", "Casq2", "Cacna1c", "Slc8a1", "Camk2d",
    "Calm1", "Calm2", "Calm3", "Trdn", "Jph2", "Hrc", "Mybpc3"
  ),
  Fibrosis_Signaling = c(
    "Smad2", "Smad3", "Smad4", "Smad7", "Ctgf", "Pdgfa", "Pdgfb",
    "Pdgfra", "Pdgfrb", "Wnt3a", "Wnt5a", "Fzd1", "Fzd2", "Axin2",
    "Dkk1", "Wisp1"
  )
)

# Annotate peaks to nearest gene
cat("Annotating peaks to nearest gene...\n")
tss <- readRDS(file.path(p2p_dir, "mm10_refGene_tss.rds"))

if (is(tss, "GRanges")) {
  tss_gr <- tss
  gene_col <- colnames(mcols(tss_gr))[1]
} else {
  chr_col <- intersect(c("chrom", "chr", "seqnames"), colnames(tss))[1]
  start_col <- intersect(c("txStart", "start", "tss"), colnames(tss))[1]
  end_col <- intersect(c("txEnd", "end"), colnames(tss))[1]
  gene_col <- intersect(c("name2", "gene_name", "symbol", "gene"), colnames(tss))[1]
  strand_col <- intersect(c("strand"), colnames(tss))[1]
  tss_pos <- ifelse(tss[[strand_col]] == "+", tss[[start_col]], tss[[end_col]])
  tss_gr <- GRanges(seqnames = tss[[chr_col]],
                     ranges = IRanges(start = tss_pos, width = 1))
  mcols(tss_gr)$gene_name <- tss[[gene_col]]
  gene_col <- "gene_name"
}

peak_gr <- GRanges(seqnames = peaks_bed$chr,
                    ranges = IRanges(start = peaks_bed$start, end = peaks_bed$end))
nearest_idx <- nearest(peak_gr, tss_gr)
peaks_bed$nearest_gene <- mcols(tss_gr)[[gene_col]][nearest_idx]
peaks_bed$dist_to_tss <- distance(peak_gr, tss_gr[nearest_idx])

# Classify peaks by module (within 200kb)
all_module_genes <- unique(unlist(gene_modules))
peaks_bed$in_any_module <- peaks_bed$nearest_gene %in% all_module_genes &
  peaks_bed$dist_to_tss <= 200000

for (mod_name in names(gene_modules)) {
  peaks_bed[[paste0("in_", mod_name)]] <- peaks_bed$nearest_gene %in% gene_modules[[mod_name]] &
    peaks_bed$dist_to_tss <= 200000
}

# Fisher's exact test: is each module enriched among reversed peaks?
cat("\nPathway enrichment (Fisher's exact): reversed vs non-reversed CUT&Tag peaks\n")
cat("(Peaks within 200kb of TSS, assigned to nearest gene)\n\n")

cat(sprintf("%-22s  %5s  %5s  %8s  %8s  %8s\n",
            "Module", "Rev", "Total", "Rev%", "OR", "P"))
cat(strrep("-", 75), "\n")

pathway_results <- list()

for (mod_name in names(gene_modules)) {
  col <- paste0("in_", mod_name)

  # 2x2: reversed_directional x in_module
  tab <- table(
    reversed = factor(peaks_bed$reversed_directional, levels = c(FALSE, TRUE)),
    in_module = factor(peaks_bed[[col]], levels = c(FALSE, TRUE))
  )

  ft <- fisher.test(tab)

  n_rev_in_mod <- sum(peaks_bed$reversed_directional & peaks_bed[[col]])
  n_total_in_mod <- sum(peaks_bed[[col]])
  pct_rev <- if (n_total_in_mod > 0) round(100 * n_rev_in_mod / n_total_in_mod, 1) else 0

  cat(sprintf("%-22s  %5d  %5d  %7.1f%%  %8.2f  %s\n",
              mod_name, n_rev_in_mod, n_total_in_mod, pct_rev,
              ft$estimate, format.pval(ft$p.value, digits = 3)))

  pathway_results[[mod_name]] <- data.frame(
    module = mod_name,
    n_reversed_in_module = n_rev_in_mod,
    n_total_in_module = n_total_in_mod,
    pct_reversed = pct_rev,
    n_reversed_not_module = sum(peaks_bed$reversed_directional & !peaks_bed[[col]]),
    n_not_reversed_in_module = sum(!peaks_bed$reversed_directional & peaks_bed[[col]]),
    OR = round(ft$estimate, 3),
    CI_low = round(ft$conf.int[1], 3),
    CI_high = round(ft$conf.int[2], 3),
    P = ft$p.value,
    stringsAsFactors = FALSE
  )
}

pathway_df <- do.call(rbind, pathway_results)
pathway_df$FDR <- p.adjust(pathway_df$P, method = "BH")
pathway_df <- pathway_df[order(pathway_df$P), ]

cat("\n\nWith FDR correction:\n")
cat(sprintf("%-22s  %8s  %8s  %8s  %10s\n", "Module", "OR", "P", "FDR", "Sig"))
cat(strrep("-", 65), "\n")
for (i in seq_len(nrow(pathway_df))) {
  sig <- ifelse(pathway_df$FDR[i] < 0.001, "***",
         ifelse(pathway_df$FDR[i] < 0.01, "**",
         ifelse(pathway_df$FDR[i] < 0.05, "*",
         ifelse(pathway_df$FDR[i] < 0.1, ".", ""))))
  cat(sprintf("%-22s  %8.2f  %8s  %8s  %10s\n",
              pathway_df$module[i],
              pathway_df$OR[i],
              format.pval(pathway_df$P[i], digits = 3),
              format.pval(pathway_df$FDR[i], digits = 3),
              sig))
}

write.csv(pathway_df, file.path(cuttag_dir, "reversed_cuttag_pathway_enrichment.csv"),
          row.names = FALSE)

# Also check: what genes are the reversed peaks near?
cat("\n\nTop genes near reversed CUT&Tag peaks (directional, within 200kb):\n")
rev_genes <- peaks_bed[peaks_bed$reversed_directional & peaks_bed$dist_to_tss <= 200000,
                        .(n_peaks = .N, mean_RI = mean(reversal_index, na.rm = TRUE)),
                        by = nearest_gene]
rev_genes <- rev_genes[order(-n_peaks)]

# Flag module membership
rev_genes$modules <- sapply(rev_genes$nearest_gene, function(g) {
  mods <- names(gene_modules)[sapply(gene_modules, function(m) g %in% m)]
  if (length(mods) == 0) return("")
  paste(mods, collapse = ";")
})

cat(sprintf("%-15s  %5s  %8s  %s\n", "Gene", "Peaks", "Mean_RI", "Modules"))
cat(strrep("-", 60), "\n")
for (i in seq_len(min(30, nrow(rev_genes)))) {
  cat(sprintf("%-15s  %5d  %8.3f  %s\n",
              rev_genes$nearest_gene[i], rev_genes$n_peaks[i],
              rev_genes$mean_RI[i], rev_genes$modules[i]))
}

# ==============================================================================
# 6. Venn diagrams
# ==============================================================================

cat("\n\nGenerating Venn diagrams...\n")

# Individual Venns
make_venn <- function(file, area1, area2, cross, cat1, cat2, fill, title) {
  pdf(file, width = 8, height = 6)
  grid.newpage()
  draw.pairwise.venn(
    area1 = area1, area2 = area2, cross.area = cross,
    category = c(cat1, cat2), fill = fill, alpha = 0.5,
    cat.cex = 1.2, cex = 1.5, cat.pos = c(-30, 30), cat.dist = c(0.05, 0.05)
  )
  grid.text(title, y = 0.95, gp = gpar(fontsize = 13, fontface = "bold"))
  dev.off()
}

make_venn(file.path(fig_dir, "venn_cuttag_vs_tobias_ap1_all.pdf"),
          length(cuttag_gr), length(fp_all_gr), ov1$n_q,
          paste0("C-FOS CUT&Tag\n(", format(length(cuttag_gr), big.mark=","), ")"),
          paste0("Fib AP-1 FP\n(", format(length(fp_all_gr), big.mark=","), ")"),
          c("#E41A1C", "#377EB8"),
          "All C-FOS CUT&Tag vs All Fib AP-1 Footprints")

make_venn(file.path(fig_dir, "venn_cuttag_vs_tobias_ap1_agegained.pdf"),
          length(cuttag_gr), length(fp_ag_gr), ov2$n_q,
          paste0("C-FOS CUT&Tag\n(", format(length(cuttag_gr), big.mark=","), ")"),
          paste0("Fib Age-Gained AP-1\n(", format(length(fp_ag_gr), big.mark=","), ")"),
          c("#E41A1C", "#FF7F00"),
          "All C-FOS CUT&Tag vs Age-Gained Fib AP-1")

make_venn(file.path(fig_dir, "venn_reversed_cuttag_vs_reversed_fp.pdf"),
          length(cuttag_rev_gr), length(fp_rev_gr), ov5$n_q,
          paste0("Reversed C-FOS\n(", format(length(cuttag_rev_gr), big.mark=","), ")"),
          paste0("Reversed Fib AP-1\n(", format(length(fp_rev_gr), big.mark=","), ")"),
          c("#4DAF4A", "#984EA3"),
          "Reversed C-FOS CUT&Tag vs Reversed Fib AP-1 Footprints")

# --- Combined 6-panel figure ---
make_combined <- function(outfile, dev_fun, ...) {
  dev_fun(outfile, ...)
  pushViewport(viewport(layout = grid.layout(2, 3)))

  draw_panel <- function(row, col, label, a1, a2, cross, c1, c2, fill, title) {
    pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
    grid.text(label, x = 0.02, y = 0.98, gp = gpar(fontsize = 16, fontface = "bold"),
              just = c("left", "top"))
    vd <- draw.pairwise.venn(
      area1 = a1, area2 = a2, cross.area = cross,
      category = c(c1, c2), fill = fill, alpha = 0.5,
      cat.cex = 0.8, cex = 1.0, cat.pos = c(-25, 25), cat.dist = c(0.06, 0.06),
      ind = FALSE
    )
    grid.draw(vd)
    grid.text(title, y = 0.93, gp = gpar(fontsize = 9, fontface = "bold"))
    popViewport()
  }

  # A: All vs All
  draw_panel(1, 1, "A", length(cuttag_gr), length(fp_all_gr), ov1$n_q,
             paste0("CUT&Tag\n(", format(length(cuttag_gr), big.mark=","), ")"),
             paste0("Fib AP-1\n(", format(length(fp_all_gr), big.mark=","), ")"),
             c("#E41A1C", "#377EB8"), "All C-FOS vs All Fib AP-1")

  # B: All vs Age-gained
  draw_panel(1, 2, "B", length(cuttag_gr), length(fp_ag_gr), ov2$n_q,
             paste0("CUT&Tag\n(", format(length(cuttag_gr), big.mark=","), ")"),
             paste0("Age-Gained\n(", format(length(fp_ag_gr), big.mark=","), ")"),
             c("#E41A1C", "#FF7F00"), "All C-FOS vs Age-Gained AP-1")

  # C: All vs Reversed FP
  draw_panel(1, 3, "C", length(cuttag_gr), length(fp_rev_gr), ov6$n_q,
             paste0("CUT&Tag\n(", format(length(cuttag_gr), big.mark=","), ")"),
             paste0("Rev FP\n(", format(length(fp_rev_gr), big.mark=","), ")"),
             c("#E41A1C", "#984EA3"), "All C-FOS vs Reversed AP-1")

  # D: Reversed CUT&Tag vs All FP
  draw_panel(2, 1, "D", length(cuttag_rev_gr), length(fp_all_gr), ov3$n_q,
             paste0("Rev CUT&Tag\n(", format(length(cuttag_rev_gr), big.mark=","), ")"),
             paste0("Fib AP-1\n(", format(length(fp_all_gr), big.mark=","), ")"),
             c("#4DAF4A", "#377EB8"), "Reversed C-FOS vs All Fib AP-1")

  # E: Reversed CUT&Tag vs Reversed FP
  draw_panel(2, 2, "E", length(cuttag_rev_gr), length(fp_rev_gr), ov5$n_q,
             paste0("Rev CUT&Tag\n(", format(length(cuttag_rev_gr), big.mark=","), ")"),
             paste0("Rev FP\n(", format(length(fp_rev_gr), big.mark=","), ")"),
             c("#4DAF4A", "#984EA3"), "Reversed CUT&Tag vs Reversed FP")

  # F: Stats summary
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3))
  grid.text("F", x = 0.02, y = 0.98, gp = gpar(fontsize = 16, fontface = "bold"),
            just = c("left", "top"))
  grid.text("Fisher's Exact Tests", y = 0.88,
            gp = gpar(fontsize = 11, fontface = "bold"))

  stats_text <- paste0(
    "1. Rev CUT&Tag vs ALL FP:\n",
    "   OR=", round(ft1$estimate, 1), " P=", format.pval(ft1$p.value, 3), "\n\n",
    "2. Rev CUT&Tag vs age-gained FP:\n",
    "   OR=", round(ft2$estimate, 1), " P=", format.pval(ft2$p.value, 3), "\n\n",
    "3. Rev CUT&Tag vs REVERSED FP:\n",
    "   OR=", round(ft3$estimate, 1), " P=", format.pval(ft3$p.value, 3), "\n\n",
    "4. RevSig CUT&Tag vs ALL FP:\n",
    "   OR=", round(ft4$estimate, 1), " P=", format.pval(ft4$p.value, 3), "\n\n",
    "5. RevSig CUT&Tag vs REVERSED FP:\n",
    "   OR=", round(ft5$estimate, 1), " P=", format.pval(ft5$p.value, 3)
  )
  grid.text(stats_text, y = 0.40, gp = gpar(fontsize = 8, fontfamily = "mono"),
            just = "center")
  popViewport()

  dev.off()
}

make_combined(file.path(fig_dir, "venn_cuttag_vs_tobias_ap1_combined.png"),
              png, width = 18, height = 12, units = "in", res = 300)
make_combined(file.path(fig_dir, "venn_cuttag_vs_tobias_ap1_combined.pdf"),
              pdf, width = 18, height = 12)

# ==============================================================================
# 7. Pathway enrichment barplot
# ==============================================================================

cat("Generating pathway enrichment barplot...\n")

pathway_df$sig_label <- ifelse(pathway_df$FDR < 0.001, "***",
                        ifelse(pathway_df$FDR < 0.01, "**",
                        ifelse(pathway_df$FDR < 0.05, "*",
                        ifelse(pathway_df$FDR < 0.1, ".", ""))))

p_enrich <- ggplot(pathway_df,
                    aes(x = reorder(module, -log10(P)), y = -log10(P))) +
  geom_col(aes(fill = OR > 1), width = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  geom_text(aes(label = sprintf("OR=%.1f\nn=%d", OR, n_reversed_in_module)),
            hjust = -0.1, size = 3) +
  scale_fill_manual(values = c("TRUE" = "#D73027", "FALSE" = "#4575B4"),
                     labels = c("Depleted", "Enriched"), name = "Reversed peaks") +
  coord_flip() +
  labs(x = NULL, y = expression(-log[10](P)),
       title = "Pathway enrichment among reversed C-FOS CUT&Tag peaks",
       subtitle = "Fisher's exact test: reversed vs non-reversed peaks near pathway genes (200kb)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(fig_dir, "pathway_enrichment_reversed_cuttag.pdf"),
       p_enrich, width = 10, height = 6, bg = "white")
ggsave(file.path(fig_dir, "pathway_enrichment_reversed_cuttag.png"),
       p_enrich, width = 10, height = 6, dpi = 300, bg = "white")

# ==============================================================================
# 8. Summary CSV
# ==============================================================================

summary_df <- data.frame(
  comparison = c(
    "All_cFOS_vs_All_FibAP1",
    "All_cFOS_vs_AgeGained_FibAP1",
    "All_cFOS_vs_Reversed_FibAP1",
    "Reversed_cFOS_vs_All_FibAP1",
    "Reversed_cFOS_vs_AgeGained_FibAP1",
    "Reversed_cFOS_vs_Reversed_FibAP1",
    "ReversedSig_cFOS_vs_Reversed_FibAP1"
  ),
  set1_name = c("C-FOS CUT&Tag", "C-FOS CUT&Tag", "C-FOS CUT&Tag",
                 "Rev C-FOS", "Rev C-FOS", "Rev C-FOS", "RevSig C-FOS"),
  set1_total = c(length(cuttag_gr), length(cuttag_gr), length(cuttag_gr),
                  length(cuttag_rev_gr), length(cuttag_rev_gr), length(cuttag_rev_gr),
                  length(cuttag_rev_sig_gr)),
  set2_name = c("Fib AP-1 All", "Fib AP-1 AgeGained", "Fib AP-1 Reversed",
                 "Fib AP-1 All", "Fib AP-1 AgeGained", "Fib AP-1 Reversed",
                 "Fib AP-1 Reversed"),
  set2_total = c(length(fp_all_gr), length(fp_ag_gr), length(fp_rev_gr),
                  length(fp_all_gr), length(fp_ag_gr), length(fp_rev_gr),
                  length(fp_rev_gr)),
  set1_overlap = c(ov1$n_q, ov2$n_q, ov6$n_q, ov3$n_q, ov4$n_q, ov5$n_q, ov7$n_q),
  set1_pct = round(100 * c(ov1$n_q/length(cuttag_gr), ov2$n_q/length(cuttag_gr),
                             ov6$n_q/length(cuttag_gr),
                             ov3$n_q/length(cuttag_rev_gr), ov4$n_q/length(cuttag_rev_gr),
                             ov5$n_q/length(cuttag_rev_gr),
                             ov7$n_q/length(cuttag_rev_sig_gr)), 1),
  set2_overlap = c(ov1$n_s, ov2$n_s, ov6$n_s, ov3$n_s, ov4$n_s, ov5$n_s, ov7$n_s),
  set2_pct = round(100 * c(ov1$n_s/length(fp_all_gr), ov2$n_s/length(fp_ag_gr),
                             ov6$n_s/length(fp_rev_gr),
                             ov3$n_s/length(fp_all_gr), ov4$n_s/length(fp_ag_gr),
                             ov5$n_s/length(fp_rev_gr),
                             ov7$n_s/length(fp_rev_gr)), 1),
  stringsAsFactors = FALSE
)

write.csv(summary_df, file.path(cuttag_dir, "cross_overlap_cuttag_tobias_summary.csv"),
          row.names = FALSE)

cat("\n", strrep("=", 70), "\n")
cat("FINAL SUMMARY TABLE\n")
cat(strrep("=", 70), "\n\n")
print(summary_df, row.names = FALSE)

cat("\n\nAll output files:\n")
cat("  figures/venn_cuttag_vs_tobias_ap1_combined.png/pdf\n")
cat("  figures/venn_cuttag_vs_tobias_ap1_all.pdf\n")
cat("  figures/venn_cuttag_vs_tobias_ap1_agegained.pdf\n")
cat("  figures/venn_reversed_cuttag_vs_reversed_fp.pdf\n")
cat("  figures/pathway_enrichment_reversed_cuttag.png/pdf\n")
cat("  cross_overlap_cuttag_tobias_summary.csv\n")
cat("  reversed_cuttag_pathway_enrichment.csv\n")
cat("\nDone.\n")
