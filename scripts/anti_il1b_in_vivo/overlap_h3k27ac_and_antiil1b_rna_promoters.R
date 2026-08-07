#!/usr/bin/env Rscript
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Link anti-IL-1β RNA responses to promoter and nearby H3K27ac peaks.
# Steps: 1) Load RNA and H3K27ac differential results; 2) map peaks to gene windows; 3) export association tests, annotations and diagnostic plots.
# =======================================
# Source reformatted from anti-il1code/overlap_h3k27ac_and_antiil1b_rna_promoters.R.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
})

project_dir <- publication_path("H3K27AC_PROJECT_ROOT")
rna_file <- file.path(project_dir, "results", "rnaseq_gsea_kegg", "edgeR_DE_with_Entrez_mapping.csv")
h3_file <- file.path(project_dir, "results", "cfos_h3k27ac_edger_tobias", "edgeR_H3K27ac_tested_peaks.csv")

out_dir <- file.path(project_dir, "results", "h3k27ac_rna_promoter_overlap")
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

promoter_upstream_bp <- 2000L
promoter_downstream_bp <- 500L
tss_window_bp <- 100000L
rna_nominal_p_cutoff <- 0.05
rna_abs_logfc_cutoff <- 0.5
h3_nominal_p_cutoff <- 0.05
h3_abs_logfc_cutoff <- 0.5

collapse_unique <- function(x) {
  x <- sort(unique(na.omit(as.character(x))))
  x <- x[nzchar(x) & x != "NA"]
  if (length(x) == 0) {
    return(NA_character_)
  }
  paste(x, collapse = ";")
}

fill_missing <- function(dt, cols, value) {
  for (col in cols) {
    if (!col %in% names(dt)) {
      set(dt, j = col, value = value)
    } else {
      idx <- which(is.na(dt[[col]]))
      if (length(idx) > 0) {
        set(dt, idx, col, value)
      }
    }
  }
  invisible(dt)
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "fg", flag = "#"))
}

fisher_feature <- function(dt, feature_col, feature_label, window_label) {
  universe <- dt[window_mapped == TRUE]
  feature_value <- universe[[feature_col]]
  if (is.logical(feature_value)) {
    universe[, feature_present_tmp := feature_value == TRUE]
  } else {
    universe[, feature_present_tmp := as.numeric(feature_value) > 0]
  }
  foreground <- universe[rna_nominal_altered == TRUE]
  background <- universe[rna_nominal_altered != TRUE]

  a <- foreground[feature_present_tmp == TRUE, .N]
  b <- foreground[feature_present_tmp != TRUE, .N]
  c <- background[feature_present_tmp == TRUE, .N]
  d <- background[feature_present_tmp != TRUE, .N]

  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
  data.table(
    window = window_label,
    feature = feature_label,
    foreground = "RNA nominal P<0.05 and abs(log2FC)>=0.5",
    foreground_genes = a + b,
    foreground_feature_genes = a,
    foreground_feature_percent = ifelse(a + b > 0, 100 * a / (a + b), NA_real_),
    background_genes = c + d,
    background_feature_genes = c,
    background_feature_percent = ifelse(c + d > 0, 100 * c / (c + d), NA_real_),
    odds_ratio = unname(ft$estimate),
    pvalue = ft$p.value
  )
}

run_cor <- function(dt, y_col, label, window_label, subset_label) {
  use <- dt[
    h3_peak_n > 0 &
      is.finite(logFC) &
      is.finite(get(y_col))
  ]
  if (nrow(use) < 4 || uniqueN(use$logFC) < 2 || uniqueN(use[[y_col]]) < 2) {
    return(data.table(
      window = window_label,
      subset = subset_label,
      h3_metric = label,
      genes = nrow(use),
      spearman_rho = NA_real_,
      pvalue = NA_real_
    ))
  }
  ct <- suppressWarnings(cor.test(use$logFC, use[[y_col]], method = "spearman", exact = FALSE))
  data.table(
    window = window_label,
    subset = subset_label,
    h3_metric = label,
    genes = nrow(use),
    spearman_rho = unname(ct$estimate),
    pvalue = ct$p.value
  )
}

direction_test <- function(dt, window_label, class_filter = NULL, label = "All altered RNA genes") {
  use <- dt[
    rna_nominal_altered == TRUE &
      h3_peak_n > 0 &
      is.finite(h3_best_logFC)
  ]
  if (!is.null(class_filter)) {
    use <- use[rna_nominal_class %in% class_filter]
  }
  use <- use[sign(logFC) != 0 & sign(h3_best_logFC) != 0]
  concordant <- use[sign(logFC) == sign(h3_best_logFC), .N]
  discordant <- use[sign(logFC) != sign(h3_best_logFC), .N]
  total <- concordant + discordant
  p <- if (total > 0) binom.test(concordant, total, p = 0.5)$p.value else NA_real_
  data.table(
    window = window_label,
    subset = label,
    genes = total,
    concordant_best_peak_genes = concordant,
    discordant_best_peak_genes = discordant,
    concordant_percent = ifelse(total > 0, 100 * concordant / total, NA_real_),
    binomial_pvalue = p
  )
}

message("Reading anti-IL1B RNA-seq DE table")
rna <- fread(rna_file)
required_rna <- c("gene_symbol", "logFC", "logCPM", "PValue", "FDR")
missing_rna <- setdiff(required_rna, names(rna))
if (length(missing_rna) > 0) {
  stop("RNA table is missing required column(s): ", paste(missing_rna, collapse = ", "))
}
rna <- rna[!is.na(gene_symbol) & gene_symbol != "" & gene_symbol != "NA"]
rna[, `:=`(
  gene_symbol = as.character(gene_symbol),
  gene_id = as.character(gene_id),
  logFC = as.numeric(logFC),
  logCPM = as.numeric(logCPM),
  PValue = as.numeric(PValue),
  FDR = as.numeric(FDR)
)]
rna[, abs_logFC_for_sort := abs(logFC)]
setorder(rna, PValue, -abs_logFC_for_sort)
rna <- rna[!duplicated(gene_symbol)]
rna[, abs_logFC_for_sort := NULL]
rna[, rna_direction := fcase(
  logFC > 0, "up",
  logFC < 0, "down",
  default = "flat"
)]
rna[, rna_nominal_p05 := is.finite(PValue) & PValue < rna_nominal_p_cutoff]
rna[, rna_nominal_altered := rna_nominal_p05 == TRUE & abs(logFC) >= rna_abs_logfc_cutoff]
rna[, rna_nominal_class := fcase(
  rna_nominal_altered == TRUE & logFC > 0, "RNA nominal up",
  rna_nominal_altered == TRUE & logFC < 0, "RNA nominal down",
  default = "RNA not altered"
)]
rna[, rna_nominal_class := factor(
  rna_nominal_class,
  levels = c("RNA not altered", "RNA nominal down", "RNA nominal up")
)]

message("Reading H3K27ac edgeR peak table")
h3 <- fread(h3_file)
required_h3 <- c("peak_id", "chr", "start", "end", "logFC", "PValue", "FDR", "tested")
missing_h3 <- setdiff(required_h3, names(h3))
if (length(missing_h3) > 0) {
  stop("H3K27ac table is missing required column(s): ", paste(missing_h3, collapse = ", "))
}
h3[, `:=`(
  peak_id = as.character(peak_id),
  chr = as.character(chr),
  start = as.integer(start),
  end = as.integer(end),
  logFC = as.numeric(logFC),
  PValue = as.numeric(PValue),
  FDR = as.numeric(FDR),
  tested = as.logical(tested)
)]
h3 <- h3[tested == TRUE]
h3[, h3_peak_row := .I]
h3[, h3_direction := fcase(
  logFC > 0, "up",
  logFC < 0, "down",
  default = "flat"
)]
h3[, h3_nominal_p05 := is.finite(PValue) & PValue < h3_nominal_p_cutoff]
h3[, h3_nominal_p05_abs := h3_nominal_p05 == TRUE & abs(logFC) >= h3_abs_logfc_cutoff]
h3[, h3_fdr_005 := is.finite(FDR) & FDR < 0.05]

h3_gr <- GRanges(
  seqnames = h3$chr,
  ranges = IRanges(start = h3$start + 1L, end = h3$end),
  h3_peak_row = h3$h3_peak_row
)
h3_gr <- keepStandardChromosomes(h3_gr, pruning.mode = "coarse")

message("Building mm10 transcript TSS windows")
tx_by_gene <- transcriptsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by = "gene")
tx <- unlist(tx_by_gene, use.names = FALSE)
entrez_ids <- rep(names(tx_by_gene), elementNROWS(tx_by_gene))
tss_pos <- ifelse(as.character(strand(tx)) == "+", start(tx), end(tx))
symbol_map <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys = unique(as.character(entrez_ids)),
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

make_tss_windows <- function(upstream, downstream, label) {
  gr <- promoters(tx, upstream = upstream, downstream = downstream)
  mcols(gr)$entrez_id <- as.character(entrez_ids)
  mcols(gr)$gene_symbol <- unname(symbol_map[as.character(entrez_ids)])
  mcols(gr)$tss_position <- as.integer(tss_pos)
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  gr <- trim(gr)
  gr <- gr[
    width(gr) > 0 &
      !is.na(mcols(gr)$gene_symbol) &
      mcols(gr)$gene_symbol != ""
  ]
  mcols(gr)$window_id <- paste0(
    label,
    "|",
    mcols(gr)$gene_symbol,
    "|",
    mcols(gr)$entrez_id,
    "|",
    as.character(seqnames(gr)),
    ":",
    start(gr),
    "-",
    end(gr),
    ":",
    as.character(strand(gr))
  )
  gr
}

window_gene_summary <- function(gr, label) {
  unique(data.table(
    window = label,
    gene_symbol = as.character(mcols(gr)$gene_symbol),
    entrez_id = as.character(mcols(gr)$entrez_id),
    window_id = as.character(mcols(gr)$window_id)
  ))[, .(
    window_mapped = TRUE,
    window_entrez_ids = collapse_unique(entrez_id),
    transcript_tss_window_n = uniqueN(window_id)
  ), by = .(window, gene_symbol)]
}

summarise_window <- function(gr, label) {
  message("Intersecting H3K27ac peaks with ", label, " windows")
  hits <- findOverlaps(gr, h3_gr, ignore.strand = TRUE)
  if (length(hits) == 0) {
    gene_peak <- data.table(
      window = character(),
      gene_symbol = character(),
      peak_id = character()
    )
  } else {
    window_idx <- queryHits(hits)
    peak_idx <- mcols(h3_gr)$h3_peak_row[subjectHits(hits)]
    peak_mid <- floor((h3$start[peak_idx] + 1L + h3$end[peak_idx]) / 2)
    hit_dt <- data.table(
      window = label,
      gene_symbol = as.character(mcols(gr)$gene_symbol[window_idx]),
      entrez_id = as.character(mcols(gr)$entrez_id[window_idx]),
      window_id = as.character(mcols(gr)$window_id[window_idx]),
      window_chr = as.character(seqnames(gr)[window_idx]),
      window_start = start(gr)[window_idx],
      window_end = end(gr)[window_idx],
      window_strand = as.character(strand(gr)[window_idx]),
      tss_position = as.integer(mcols(gr)$tss_position[window_idx]),
      peak_id = h3$peak_id[peak_idx],
      peak_chr = h3$chr[peak_idx],
      peak_start_0based = h3$start[peak_idx],
      peak_end = h3$end[peak_idx],
      peak_midpoint_1based = peak_mid,
      distance_to_tss_bp = abs(peak_mid - as.integer(mcols(gr)$tss_position[window_idx])),
      h3_logFC = h3$logFC[peak_idx],
      h3_logCPM = h3$logCPM[peak_idx],
      h3_PValue = h3$PValue[peak_idx],
      h3_FDR = h3$FDR[peak_idx],
      h3_nominal_p05 = h3$h3_nominal_p05[peak_idx],
      h3_nominal_p05_abs = h3$h3_nominal_p05_abs[peak_idx],
      h3_fdr_005 = h3$h3_fdr_005[peak_idx],
      h3_direction = h3$h3_direction[peak_idx]
    )

    gene_peak <- hit_dt[, .(
      window_entrez_ids = collapse_unique(entrez_id),
      window_ids = collapse_unique(window_id),
      window_chr = collapse_unique(window_chr),
      peak_chr = peak_chr[1],
      peak_start_0based = peak_start_0based[1],
      peak_end = peak_end[1],
      peak_midpoint_1based = peak_midpoint_1based[1],
      min_distance_to_tss_bp = min(distance_to_tss_bp, na.rm = TRUE),
      h3_logFC = h3_logFC[1],
      h3_logCPM = h3_logCPM[1],
      h3_PValue = h3_PValue[1],
      h3_FDR = h3_FDR[1],
      h3_nominal_p05 = h3_nominal_p05[1],
      h3_nominal_p05_abs = h3_nominal_p05_abs[1],
      h3_fdr_005 = h3_fdr_005[1],
      h3_direction = h3_direction[1]
    ), by = .(window, gene_symbol, peak_id)]
  }

  if (nrow(gene_peak) > 0) {
    gene_peak[, abs_h3_logFC_for_sort := abs(h3_logFC)]
    setorder(gene_peak, gene_symbol, h3_PValue, -abs_h3_logFC_for_sort)
    gene_peak[, abs_h3_logFC_for_sort := NULL]
    best <- gene_peak[, .SD[1], by = gene_symbol]
    best <- best[, .(
      gene_symbol,
      h3_best_peak_id = peak_id,
      h3_best_peak_chr = peak_chr,
      h3_best_peak_start_0based = peak_start_0based,
      h3_best_peak_end = peak_end,
      h3_best_peak_min_distance_to_tss_bp = min_distance_to_tss_bp,
      h3_best_logFC = h3_logFC,
      h3_best_logCPM = h3_logCPM,
      h3_best_PValue = h3_PValue,
      h3_best_FDR = h3_FDR,
      h3_best_nominal_p05 = h3_nominal_p05,
      h3_best_nominal_p05_abs = h3_nominal_p05_abs,
      h3_best_fdr_005 = h3_fdr_005
    )]

    summary <- gene_peak[, .(
      h3_peak_n = .N,
      h3_nominal_p05_peak_n = sum(h3_nominal_p05 == TRUE),
      h3_nominal_p05_abs_peak_n = sum(h3_nominal_p05_abs == TRUE),
      h3_fdr_005_peak_n = sum(h3_fdr_005 == TRUE),
      h3_any_nominal_up = any(h3_nominal_p05 == TRUE & h3_logFC > 0),
      h3_any_nominal_down = any(h3_nominal_p05 == TRUE & h3_logFC < 0),
      h3_any_nominal_abs_up = any(h3_nominal_p05_abs == TRUE & h3_logFC > 0),
      h3_any_nominal_abs_down = any(h3_nominal_p05_abs == TRUE & h3_logFC < 0),
      h3_any_fdr_up = any(h3_fdr_005 == TRUE & h3_logFC > 0),
      h3_any_fdr_down = any(h3_fdr_005 == TRUE & h3_logFC < 0),
      h3_mean_logFC = mean(h3_logFC, na.rm = TRUE),
      h3_median_logFC = median(h3_logFC, na.rm = TRUE),
      h3_min_PValue = min(h3_PValue, na.rm = TRUE),
      h3_min_FDR = min(h3_FDR, na.rm = TRUE),
      h3_max_abs_logFC = max(abs(h3_logFC), na.rm = TRUE),
      h3_min_distance_to_tss_bp = min(min_distance_to_tss_bp, na.rm = TRUE),
      h3_median_distance_to_tss_bp = median(min_distance_to_tss_bp, na.rm = TRUE)
    ), by = gene_symbol]
    summary <- merge(summary, best, by = "gene_symbol", all.x = TRUE, sort = FALSE)
  } else {
    summary <- data.table(gene_symbol = character())
  }

  gene_map <- window_gene_summary(gr, label)
  annot <- merge(rna, gene_map, by = "gene_symbol", all.x = TRUE, sort = FALSE)
  annot <- merge(annot, summary, by = "gene_symbol", all.x = TRUE, sort = FALSE)
  annot[, window := label]
  fill_missing(annot, "window_mapped", FALSE)
  fill_missing(
    annot,
    c("h3_peak_n", "h3_nominal_p05_peak_n", "h3_nominal_p05_abs_peak_n", "h3_fdr_005_peak_n"),
    0L
  )
  fill_missing(
    annot,
    c(
      "h3_any_nominal_up", "h3_any_nominal_down",
      "h3_any_nominal_abs_up", "h3_any_nominal_abs_down",
      "h3_any_fdr_up", "h3_any_fdr_down"
    ),
    FALSE
  )
  annot[, h3_best_same_direction := h3_peak_n > 0 & sign(logFC) == sign(h3_best_logFC)]
  annot[, h3_mean_same_direction := h3_peak_n > 0 & sign(logFC) == sign(h3_mean_logFC)]
  annot[, h3_nominal_same_direction := fcase(
    rna_direction == "up" & h3_any_nominal_up == TRUE, TRUE,
    rna_direction == "down" & h3_any_nominal_down == TRUE, TRUE,
    default = FALSE
  )]
  annot[, h3_nominal_opposite_direction := fcase(
    rna_direction == "up" & h3_any_nominal_down == TRUE, TRUE,
    rna_direction == "down" & h3_any_nominal_up == TRUE, TRUE,
    default = FALSE
  )]
  annot[, h3_nominal_abs_same_direction := fcase(
    rna_direction == "up" & h3_any_nominal_abs_up == TRUE, TRUE,
    rna_direction == "down" & h3_any_nominal_abs_down == TRUE, TRUE,
    default = FALSE
  )]
  annot[, h3_nominal_abs_opposite_direction := fcase(
    rna_direction == "up" & h3_any_nominal_abs_down == TRUE, TRUE,
    rna_direction == "down" & h3_any_nominal_abs_up == TRUE, TRUE,
    default = FALSE
  )]

  status_summary <- annot[window_mapped == TRUE, .(
    genes = .N,
    genes_with_h3_peak = sum(h3_peak_n > 0),
    percent_with_h3_peak = 100 * mean(h3_peak_n > 0),
    genes_with_h3_nominal_peak = sum(h3_nominal_p05_peak_n > 0),
    percent_with_h3_nominal_peak = 100 * mean(h3_nominal_p05_peak_n > 0),
    median_h3_best_logFC = median(h3_best_logFC, na.rm = TRUE),
    median_h3_mean_logFC = median(h3_mean_logFC, na.rm = TRUE),
    best_peak_concordant_genes = sum(h3_peak_n > 0 & h3_best_same_direction == TRUE),
    best_peak_discordant_genes = sum(h3_peak_n > 0 & h3_best_same_direction != TRUE),
    best_peak_concordant_percent = ifelse(
      sum(h3_peak_n > 0) > 0,
      100 * sum(h3_peak_n > 0 & h3_best_same_direction == TRUE) / sum(h3_peak_n > 0),
      NA_real_
    )
  ), by = rna_nominal_class]
  status_summary[, window := label]
  setcolorder(status_summary, c("window", setdiff(names(status_summary), "window")))

  enrichment <- rbindlist(list(
    fisher_feature(annot, "h3_peak_n", "Any H3K27ac peak in window", label),
    fisher_feature(annot, "h3_nominal_p05_peak_n", "Any nominal P<0.05 H3K27ac peak in window", label),
    fisher_feature(annot, "h3_nominal_p05_abs_peak_n", "Any nominal P<0.05 and abs(log2FC)>=0.5 H3K27ac peak in window", label),
    fisher_feature(annot, "h3_fdr_005_peak_n", "Any FDR<0.05 H3K27ac peak in window", label)
  ), use.names = TRUE, fill = TRUE)

  correlations <- rbindlist(list(
    run_cor(annot, "h3_best_logFC", "Best promoter/nearby H3K27ac peak log2FC", label, "All genes with H3K27ac peak"),
    run_cor(annot[rna_nominal_p05 == TRUE], "h3_best_logFC", "Best promoter/nearby H3K27ac peak log2FC", label, "RNA nominal P<0.05"),
    run_cor(annot[rna_nominal_altered == TRUE], "h3_best_logFC", "Best promoter/nearby H3K27ac peak log2FC", label, "RNA nominal P<0.05 and abs(log2FC)>=0.5"),
    run_cor(annot, "h3_mean_logFC", "Mean H3K27ac peak log2FC", label, "All genes with H3K27ac peak"),
    run_cor(annot[rna_nominal_p05 == TRUE], "h3_mean_logFC", "Mean H3K27ac peak log2FC", label, "RNA nominal P<0.05"),
    run_cor(annot[rna_nominal_altered == TRUE], "h3_mean_logFC", "Mean H3K27ac peak log2FC", label, "RNA nominal P<0.05 and abs(log2FC)>=0.5")
  ), use.names = TRUE, fill = TRUE)

  best_direction <- rbindlist(list(
    direction_test(annot, label),
    direction_test(annot, label, "RNA nominal up", "RNA nominal up"),
    direction_test(annot, label, "RNA nominal down", "RNA nominal down")
  ), use.names = TRUE, fill = TRUE)

  peak_direction <- data.table()
  if (nrow(gene_peak) > 0) {
    peak_direction <- merge(
      gene_peak,
      rna[, .(gene_symbol, rna_logFC = logFC, rna_PValue = PValue, rna_FDR = FDR, rna_nominal_altered, rna_nominal_class)],
      by = "gene_symbol",
      all.x = TRUE,
      sort = FALSE
    )
    peak_direction <- peak_direction[
      rna_nominal_altered == TRUE &
        h3_nominal_p05 == TRUE &
        sign(rna_logFC) != 0 &
        sign(h3_logFC) != 0
    ]
    peak_direction[, pair_direction := fifelse(sign(rna_logFC) == sign(h3_logFC), "same", "opposite")]
    peak_direction <- peak_direction[, .(
      window = label,
      nominal_h3_gene_peak_pairs = .N,
      genes = uniqueN(gene_symbol),
      same_direction_pairs = sum(pair_direction == "same"),
      opposite_direction_pairs = sum(pair_direction == "opposite"),
      same_direction_percent = 100 * mean(pair_direction == "same")
    ), by = rna_nominal_class]
  }

  candidates <- annot[
    rna_nominal_altered == TRUE & h3_peak_n > 0,
    .(
      window,
      gene_symbol,
      gene_id,
      rna_logFC = logFC,
      rna_PValue = PValue,
      rna_FDR = FDR,
      rna_nominal_class,
      h3_peak_n,
      h3_nominal_p05_peak_n,
      h3_nominal_p05_abs_peak_n,
      h3_fdr_005_peak_n,
      h3_best_same_direction,
      h3_nominal_same_direction,
      h3_nominal_opposite_direction,
      h3_best_peak_id,
      h3_best_peak_chr,
      h3_best_peak_start_0based,
      h3_best_peak_end,
      h3_best_peak_min_distance_to_tss_bp,
      h3_best_logFC,
      h3_best_PValue,
      h3_best_FDR,
      h3_mean_logFC,
      h3_median_logFC,
      h3_min_PValue,
      h3_min_FDR
    )
  ]
  setorder(candidates, -h3_best_same_direction, h3_best_PValue, rna_PValue)

  mapping_summary <- data.table(
    window = label,
    metric = c(
      "RNA genes in DE table",
      "RNA genes mapped to mm10 transcript window",
      "RNA nominal P<0.05 genes",
      "RNA nominal P<0.05 and abs(log2FC)>=0.5 genes",
      "RNA altered genes with any H3K27ac peak in window",
      "RNA altered genes with any nominal H3K27ac peak in window",
      "RNA altered genes with same-direction best H3K27ac peak",
      "RNA altered genes with same-direction nominal H3K27ac peak",
      "Unique H3K27ac peaks overlapping any RNA gene window",
      "H3K27ac gene-peak pairs in RNA gene windows",
      "Unique FDR<0.05 H3K27ac peaks in RNA gene windows"
    ),
    value = c(
      nrow(rna),
      annot[window_mapped == TRUE, .N],
      annot[rna_nominal_p05 == TRUE, .N],
      annot[rna_nominal_altered == TRUE, .N],
      annot[rna_nominal_altered == TRUE & h3_peak_n > 0, .N],
      annot[rna_nominal_altered == TRUE & h3_nominal_p05_peak_n > 0, .N],
      annot[rna_nominal_altered == TRUE & h3_peak_n > 0 & h3_best_same_direction == TRUE, .N],
      annot[rna_nominal_altered == TRUE & h3_nominal_same_direction == TRUE, .N],
      ifelse(nrow(gene_peak) > 0, uniqueN(gene_peak$peak_id), 0L),
      nrow(gene_peak),
      ifelse(nrow(gene_peak) > 0, uniqueN(gene_peak[h3_fdr_005 == TRUE, peak_id]), 0L)
    )
  )

  fwrite(gene_peak, file.path(out_dir, paste0(label, "_H3K27ac_peak_gene_pairs.csv")))
  fwrite(annot, file.path(out_dir, paste0(label, "_all_RNA_genes_with_H3K27ac_window_annotation.csv")))
  fwrite(
    annot[rna_nominal_altered == TRUE],
    file.path(out_dir, paste0(label, "_RNA_nominal_altered_genes_with_H3K27ac_window_annotation.csv"))
  )
  fwrite(candidates, file.path(out_dir, paste0(label, "_RNA_H3K27ac_directional_candidate_genes.csv")))
  fwrite(status_summary, file.path(out_dir, paste0(label, "_RNA_class_H3K27ac_summary.csv")))
  fwrite(enrichment, file.path(out_dir, paste0(label, "_RNA_altered_H3K27ac_enrichment_tests.csv")))
  fwrite(correlations, file.path(out_dir, paste0(label, "_RNA_H3K27ac_correlation_tests.csv")))
  fwrite(best_direction, file.path(out_dir, paste0(label, "_RNA_H3K27ac_best_peak_direction_tests.csv")))
  fwrite(peak_direction, file.path(out_dir, paste0(label, "_RNA_H3K27ac_nominal_peak_direction_summary.csv")))
  fwrite(mapping_summary, file.path(out_dir, paste0(label, "_mapping_summary.csv")))

  plot_dt <- annot[h3_peak_n > 0 & is.finite(h3_best_logFC)]
  if (nrow(plot_dt) > 0) {
    scatter <- ggplot(plot_dt, aes(x = logFC, y = h3_best_logFC, color = rna_nominal_class)) +
      geom_hline(yintercept = 0, linewidth = 0.3, color = "grey55") +
      geom_vline(xintercept = 0, linewidth = 0.3, color = "grey55") +
      geom_point(alpha = 0.45, size = 1.4) +
      scale_color_manual(
        values = c("RNA not altered" = "grey70", "RNA nominal down" = "#2C7BB6", "RNA nominal up" = "#D7191C"),
        drop = FALSE
      ) +
      theme_classic(base_size = 10) +
      labs(
        x = "RNA log2FC: anti-IL1B vs IgG",
        y = paste0(label, " best H3K27ac peak log2FC"),
        color = NULL
      )
    ggsave(
      file.path(plot_dir, paste0(label, "_RNA_logFC_vs_best_H3K27ac_logFC.png")),
      scatter,
      width = 5.8,
      height = 4.4,
      dpi = 300
    )
    ggsave(
      file.path(plot_dir, paste0(label, "_RNA_logFC_vs_best_H3K27ac_logFC.pdf")),
      scatter,
      width = 5.8,
      height = 4.4
    )

    altered_plot <- plot_dt[rna_nominal_altered == TRUE]
    if (nrow(altered_plot) > 0) {
      dir_plot_dt <- altered_plot[, .(
        genes = .N
      ), by = .(
        rna_nominal_class,
        direction = fifelse(h3_best_same_direction == TRUE, "Same direction", "Opposite direction")
      )]
      dir_plot <- ggplot(dir_plot_dt, aes(x = rna_nominal_class, y = genes, fill = direction)) +
        geom_col(width = 0.72, color = "white", linewidth = 0.2) +
        scale_fill_manual(values = c("Same direction" = "#1A9850", "Opposite direction" = "#B2182B")) +
        theme_classic(base_size = 10) +
        labs(x = NULL, y = "RNA-altered genes with H3K27ac peak", fill = NULL)
      ggsave(
        file.path(plot_dir, paste0(label, "_RNA_altered_best_H3K27ac_direction_bar.png")),
        dir_plot,
        width = 4.4,
        height = 3.6,
        dpi = 300
      )
      ggsave(
        file.path(plot_dir, paste0(label, "_RNA_altered_best_H3K27ac_direction_bar.pdf")),
        dir_plot,
        width = 4.4,
        height = 3.6
      )
    }

    box <- ggplot(plot_dt, aes(x = rna_nominal_class, y = h3_best_logFC, fill = rna_nominal_class)) +
      geom_hline(yintercept = 0, linewidth = 0.3, color = "grey55") +
      geom_boxplot(outlier.shape = NA, width = 0.65, alpha = 0.85) +
      geom_point(position = position_jitter(width = 0.15, height = 0), alpha = 0.18, size = 0.7) +
      scale_fill_manual(
        values = c("RNA not altered" = "grey80", "RNA nominal down" = "#74ADD1", "RNA nominal up" = "#F46D43"),
        drop = FALSE
      ) +
      theme_classic(base_size = 10) +
      theme(legend.position = "none") +
      labs(x = NULL, y = paste0(label, " best H3K27ac peak log2FC"))
    ggsave(
      file.path(plot_dir, paste0(label, "_best_H3K27ac_logFC_by_RNA_class.png")),
      box,
      width = 5.2,
      height = 3.8,
      dpi = 300
    )
    ggsave(
      file.path(plot_dir, paste0(label, "_best_H3K27ac_logFC_by_RNA_class.pdf")),
      box,
      width = 5.2,
      height = 3.8
    )
  }

  list(
    annot = annot,
    gene_peak = gene_peak,
    status_summary = status_summary,
    enrichment = enrichment,
    correlations = correlations,
    best_direction = best_direction,
    peak_direction = peak_direction,
    mapping_summary = mapping_summary
  )
}

prom_gr <- make_tss_windows(promoter_upstream_bp, promoter_downstream_bp, "promoter")
tss100_gr <- make_tss_windows(tss_window_bp, tss_window_bp + 1L, "TSS100kb")

promoter_res <- summarise_window(prom_gr, "promoter")
tss100_res <- summarise_window(tss100_gr, "TSS100kb")

combined_mapping <- rbindlist(list(promoter_res$mapping_summary, tss100_res$mapping_summary), use.names = TRUE, fill = TRUE)
combined_status <- rbindlist(list(promoter_res$status_summary, tss100_res$status_summary), use.names = TRUE, fill = TRUE)
combined_enrichment <- rbindlist(list(promoter_res$enrichment, tss100_res$enrichment), use.names = TRUE, fill = TRUE)
combined_cor <- rbindlist(list(promoter_res$correlations, tss100_res$correlations), use.names = TRUE, fill = TRUE)
combined_best_dir <- rbindlist(list(promoter_res$best_direction, tss100_res$best_direction), use.names = TRUE, fill = TRUE)
combined_peak_dir <- rbindlist(list(promoter_res$peak_direction, tss100_res$peak_direction), use.names = TRUE, fill = TRUE)

fwrite(combined_mapping, file.path(out_dir, "combined_mapping_summary.csv"))
fwrite(combined_status, file.path(out_dir, "combined_RNA_class_H3K27ac_summary.csv"))
fwrite(combined_enrichment, file.path(out_dir, "combined_RNA_altered_H3K27ac_enrichment_tests.csv"))
fwrite(combined_cor, file.path(out_dir, "combined_RNA_H3K27ac_correlation_tests.csv"))
fwrite(combined_best_dir, file.path(out_dir, "combined_RNA_H3K27ac_best_peak_direction_tests.csv"))
fwrite(combined_peak_dir, file.path(out_dir, "combined_RNA_H3K27ac_nominal_peak_direction_summary.csv"))

promoter_key <- promoter_res$best_direction[subset == "All altered RNA genes"]
tss100_key <- tss100_res$best_direction[subset == "All altered RNA genes"]
promoter_cor <- promoter_res$correlations[
  subset == "RNA nominal P<0.05 and abs(log2FC)>=0.5" &
    h3_metric == "Best promoter/nearby H3K27ac peak log2FC"
]
tss100_cor <- tss100_res$correlations[
  subset == "RNA nominal P<0.05 and abs(log2FC)>=0.5" &
    h3_metric == "Best promoter/nearby H3K27ac peak log2FC"
]

summary_lines <- c(
  "# Anti-IL1B RNA vs H3K27ac overlap",
  "",
  paste0("- RNA contrast: anti-IL1B vs IgG; altered RNA genes were defined as P < ", rna_nominal_p_cutoff, " and abs(log2FC) >= ", rna_abs_logfc_cutoff, "."),
  paste0("- H3K27ac contrast: anti-IL1B vs IgG; promoter was TSS -", promoter_upstream_bp, " bp to +", promoter_downstream_bp, " bp."),
  paste0("- Broader regulatory window was TSS +/-", tss_window_bp / 1000, " kb."),
  paste0("- H3K27ac nominal peaks used P < ", h3_nominal_p_cutoff, "; stricter nominal peaks additionally used abs(log2FC) >= ", h3_abs_logfc_cutoff, "."),
  "",
  "## Directional support",
  "",
  paste0(
    "- Promoter best-peak direction among altered RNA genes: ",
    promoter_key$concordant_best_peak_genes, "/", promoter_key$genes,
    " concordant (", fmt_num(promoter_key$concordant_percent), "%; binomial P = ",
    fmt_num(promoter_key$binomial_pvalue, 4), ")."
  ),
  paste0(
    "- TSS +/-100 kb best-peak direction among altered RNA genes: ",
    tss100_key$concordant_best_peak_genes, "/", tss100_key$genes,
    " concordant (", fmt_num(tss100_key$concordant_percent), "%; binomial P = ",
    fmt_num(tss100_key$binomial_pvalue, 4), ")."
  ),
  paste0(
    "- Promoter RNA-H3K27ac Spearman correlation among altered RNA genes: rho = ",
    fmt_num(promoter_cor$spearman_rho), ", P = ", fmt_num(promoter_cor$pvalue, 4), "."
  ),
  paste0(
    "- TSS +/-100 kb RNA-H3K27ac Spearman correlation among altered RNA genes: rho = ",
    fmt_num(tss100_cor$spearman_rho), ", P = ", fmt_num(tss100_cor$pvalue, 4), "."
  ),
  "",
  "## Main outputs",
  "",
  "- `combined_mapping_summary.csv`",
  "- `combined_RNA_class_H3K27ac_summary.csv`",
  "- `combined_RNA_altered_H3K27ac_enrichment_tests.csv`",
  "- `combined_RNA_H3K27ac_correlation_tests.csv`",
  "- `combined_RNA_H3K27ac_best_peak_direction_tests.csv`",
  "- `promoter_RNA_H3K27ac_directional_candidate_genes.csv`",
  "- `TSS100kb_RNA_H3K27ac_directional_candidate_genes.csv`",
  "- `plots/promoter_RNA_logFC_vs_best_H3K27ac_logFC.png`",
  "- `plots/TSS100kb_RNA_logFC_vs_best_H3K27ac_logFC.png`"
)
writeLines(summary_lines, file.path(out_dir, "antiIL1B_RNA_H3K27ac_overlap_summary.md"))

cat("Wrote anti-IL1B RNA/H3K27ac overlap analysis to ", out_dir, "\n", sep = "")
