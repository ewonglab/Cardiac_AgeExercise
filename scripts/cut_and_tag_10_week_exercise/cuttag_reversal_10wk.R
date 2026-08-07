# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Analyse 10-week c-FOS CUT&Tag peaks and reversal across aging and exercise.
# Steps: 1) Build consensus peaks and counts; 2) fit edgeR contrasts; 3) annotate motif/reversal classes and write summary outputs.
# =======================================
# Source reformatted from 10_week_exercsie_Cut&Tag/cuttag_reversal_10wk.R.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

#!/usr/bin/env Rscript
# ==============================================================================
# Cut&Tag 10wk: AP-1 (C-FOS) reversal analysis
#
# Steps:
#   1. QC peak counts and sample exclusion
#   2. Consensus peaks via igraph clustering of reciprocal-overlap pairs
#   3. H3K27ac consensus + TF-at-H3K27ac subset
#   4. Fragment counting with featureCounts (Rsubread)
#   5. edgeR differential binding
#   6. Reversal analysis (intersection + continuous index)
#   7. AP-1 motif annotation + reversal in motif+ subset
#   8. Master results table and summary figures
#
# Requires: run_consensus_and_counts.sh to have been run first (bedtools steps)
# ==============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(igraph)
    library(edgeR)
    library(ggplot2)
    library(Rsubread)
})

# --- Configuration -----------------------------------------------------------
EXPERIMENT <- "10wk"
TF_MARK    <- "C-FOS"
PREFIX     <- "1"  # sample prefix (1OE, 1OS, 1YS)

ANALYSIS_DIR <- publication_path("CUTTAG_10_WEEK_ANALYSIS_ROOT")
setwd(ANALYSIS_DIR)

PEAK_DIR <- publication_path("CUTTAG_10_WEEK_ROOT", "post_analysis_20260217_v1/deliverables_10wk/KEY_RESULTS/02_Peaks/called_peaks_nfcore/seacr")

BAM_DIR <- publication_path("CUTTAG_10_WEEK_ROOT", "results_mm10/02_alignment/bowtie2/target/markdup")

MOTIF_DIR <- publication_path("CUTTAG_10_WEEK_ROOT", "post_analysis_20260217_v1/deliverables_10wk/KEY_RESULTS/05_Motif_Locations/ap1_hits_parsed")

FIG_DIR <- file.path(ANALYSIS_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

cat("Working directory:", ANALYSIS_DIR, "\n")
cat("Peak directory:", PEAK_DIR, "\n")

# ==============================================================================
# STEP 1: QC -Peak count table and sample exclusion
# ==============================================================================
cat("\n===== STEP 1: QC peak counts =====\n")

bed_files <- list.files(PEAK_DIR, pattern = "\\.seacr\\.peaks\\.stringent\\.bed$",
                        full.names = TRUE)

qc <- data.table(
    file = basename(bed_files),
    sample = sub("\\.seacr\\.peaks\\.stringent\\.bed$", "", basename(bed_files))
)

# Count peaks per file
qc[, n_peaks := sapply(bed_files, function(f) {
    nrow(fread(f, header = FALSE, select = 1:3))
})]

# Parse sample metadata
qc[, mark := sub("^[12][A-Z]+_", "", sub("_R[0-9]+$", "", sample))]
qc[, condition := sub("^([0-9]+[A-Z]+)_.*", "\\1", sample)]
qc[, replicate := sub(".*_R", "R", sample)]

fwrite(qc, "qc_peak_counts.csv")
cat("QC peak counts written to qc_peak_counts.csv\n")
print(qc[order(mark, condition, replicate)])

# Define exclusions (< 100 peaks)
excluded_samples <- qc[n_peaks < 100, sample]
cat("\nExcluded samples (< 100 peaks):\n")
if (length(excluded_samples) > 0) {
    print(qc[sample %in% excluded_samples, .(sample, n_peaks)])
} else {
    cat("  None\n")
}

# QC-passing samples by mark
qc_pass_tf      <- qc[mark == TF_MARK & n_peaks >= 100, sample]
qc_pass_h3k27ac <- qc[mark == "H3K27ac" & n_peaks >= 100, sample]

cat("\nQC-passing C-FOS samples:", length(qc_pass_tf), "\n")
cat("QC-passing H3K27ac samples:", length(qc_pass_h3k27ac), "\n")

# ==============================================================================
# STEP 2: Consensus peaks -igraph clustering of reciprocal-overlap pairs
# ==============================================================================
cat("\n===== STEP 2: Consensus peak building =====\n")

build_consensus_from_pairs <- function(mark, qc_pass_samples) {
    pairs_file <- paste0("reciprocal_pairs_", mark, ".tsv")
    pooled_file <- paste0("all_peaks_pooled_", mark, ".bed")

    if (!file.exists(pairs_file)) {
        stop("Missing ", pairs_file, " -- run run_consensus_and_counts.sh first")
    }

    # Read pooled peaks (these are the nodes)
    pooled <- fread(pooled_file, header = FALSE,
                    col.names = c("chr", "start", "end", "sample"))

    # Create unique peak ID
    pooled[, peak_id := paste0(chr, ":", start, "-", end, "|", sample)]

    # Read reciprocal overlap pairs
    pairs <- fread(pairs_file, header = FALSE,
                   col.names = c("chr_a", "start_a", "end_a", "sample_a",
                                 "chr_b", "start_b", "end_b", "sample_b"))

    pairs[, id_a := paste0(chr_a, ":", start_a, "-", end_a, "|", sample_a)]
    pairs[, id_b := paste0(chr_b, ":", start_b, "-", end_b, "|", sample_b)]

    # Remove self-overlaps
    pairs <- pairs[id_a != id_b]

    cat("  Building graph for", mark, ":",
        nrow(pooled), "peaks,", nrow(pairs), "edges\n")

    # Build graph and find connected components
    g <- graph_from_data_frame(pairs[, .(id_a, id_b)], directed = FALSE,
                               vertices = pooled$peak_id)
    components <- components(g)

    pooled[, cluster := components$membership[peak_id]]

    # For each cluster: compute consensus region (union) and count support
    consensus <- pooled[, .(
        chr   = chr[1],
        start = min(start),
        end   = max(end),
        n_peaks = .N,
        samples = paste(sort(unique(sample)), collapse = ",")
    ), by = cluster]

    # Count replicates per condition
    conditions <- c(paste0(PREFIX, "OE"), paste0(PREFIX, "OS"),
                    paste0(PREFIX, "YS"))

    for (cond in conditions) {
        col_name <- paste0("n_", sub(PREFIX, "", cond))
        consensus[, (col_name) := sapply(strsplit(samples, ","), function(s) {
            sum(grepl(paste0("^", cond, "_"), s))
        })]
    }

    # Support filter: >= 2 replicates in at least one condition
    consensus[, max_support := pmax(n_OE, n_OS, n_YS)]
    consensus_filt <- consensus[max_support >= 2]

    cat("  Clusters before filter:", nrow(consensus), "\n")
    cat("  Clusters after support filter (>=2 in any condition):",
        nrow(consensus_filt), "\n")
    cat("  Support per condition (median): OE=", median(consensus_filt$n_OE),
        " OS=", median(consensus_filt$n_OS),
        " YS=", median(consensus_filt$n_YS), "\n")

    # Create BED output
    consensus_bed <- consensus_filt[order(chr, start),
        .(chr, start, end,
          name = paste0(chr, ":", start, "-", end),
          score = n_peaks,
          strand = ".")]

    bed_out <- paste0("consensus_peaks_", mark, ".bed")
    fwrite(consensus_bed, bed_out, sep = "\t", col.names = FALSE)
    cat("  Written:", bed_out, "(", nrow(consensus_bed), "peaks )\n")

    # Support CSV
    support_out <- paste0("consensus_support_", mark, ".csv")
    fwrite(consensus_filt[order(chr, start)], support_out)
    cat("  Written:", support_out, "\n")

    return(consensus_filt)
}

# Build TF consensus
consensus_tf <- build_consensus_from_pairs(TF_MARK, qc_pass_tf)

# Build H3K27ac consensus
consensus_h3k <- build_consensus_from_pairs("H3K27ac", qc_pass_h3k27ac)

# ==============================================================================
# STEP 2d: AP-1 peaks overlapping H3K27ac (active regulatory elements)
# ==============================================================================
cat("\n===== STEP 2d: TF at H3K27ac overlap =====\n")

# Check bedtools is available
bedtools_check <- suppressWarnings(system2("bedtools", "--version", stdout = TRUE, stderr = TRUE))
if (length(bedtools_check) == 0 || !grepl("bedtools", bedtools_check[1])) {
    cat("WARNING: bedtools not found in PATH. Using R-based overlap instead.\n")
    # R-based fallback: load both BEDs, find any overlap with >= 10% reciprocal
    tf_bed <- fread(paste0("consensus_peaks_", TF_MARK, ".bed"), header = FALSE,
                    col.names = c("chr","start","end","name","score","strand"))
    h3k_bed <- fread("consensus_peaks_H3K27ac.bed", header = FALSE,
                     col.names = c("chr","start","end","name","score","strand"))
    # Use GenomicRanges for overlap
    suppressPackageStartupMessages(library(GenomicRanges))
    tf_gr  <- GRanges(tf_bed$chr, IRanges(tf_bed$start, tf_bed$end), name = tf_bed$name)
    h3k_gr <- GRanges(h3k_bed$chr, IRanges(h3k_bed$start, h3k_bed$end))
    hits <- findOverlaps(tf_gr, h3k_gr)
    # Check 10% reciprocal overlap
    tf_w  <- width(tf_gr[queryHits(hits)])
    h3k_w <- width(h3k_gr[subjectHits(hits)])
    olap  <- pmin(end(tf_gr[queryHits(hits)]), end(h3k_gr[subjectHits(hits)])) -
             pmax(start(tf_gr[queryHits(hits)]), start(h3k_gr[subjectHits(hits)]))
    keep_hits <- (olap / tf_w >= 0.10) & (olap / h3k_w >= 0.10)
    overlap_idx <- unique(queryHits(hits)[keep_hits])
    h3k_overlap_bed <- tf_bed[overlap_idx, ]
    fwrite(h3k_overlap_bed, paste0("consensus_peaks_", TF_MARK, "_at_H3K27ac.bed"),
           sep = "\t", col.names = FALSE)
} else {
    cat("bedtools found:", bedtools_check[1], "\n")
    ret <- system2("bedtools", args = c("intersect",
        "-a", paste0("consensus_peaks_", TF_MARK, ".bed"),
        "-b", "consensus_peaks_H3K27ac.bed",
        "-f", "0.10", "-u"),
        stdout = paste0("consensus_peaks_", TF_MARK, "_at_H3K27ac.bed"))
    if (ret != 0) warning("bedtools intersect returned non-zero exit status")
    h3k_overlap_bed <- fread(paste0("consensus_peaks_", TF_MARK, "_at_H3K27ac.bed"),
                             header = FALSE, col.names = c("chr","start","end","name","score","strand"))
}
cat("TF peaks at H3K27ac:", nrow(h3k_overlap_bed), "\n")

# ==============================================================================
# STEP 3: Fragment counting -featureCounts (Rsubread)
# ==============================================================================
cat("\n===== STEP 3: featureCounts =====\n")

# Convert consensus BED to SAF format
consensus_bed <- fread(paste0("consensus_peaks_", TF_MARK, ".bed"),
                       header = FALSE,
                       col.names = c("chr","start","end","name","score","strand"))

saf <- data.frame(
    GeneID = consensus_bed$name,
    Chr    = consensus_bed$chr,
    Start  = consensus_bed$start,
    End    = consensus_bed$end,
    Strand = "."
)

# BAM file paths for QC-passing TF samples
bam_paths <- file.path(BAM_DIR, paste0(qc_pass_tf, ".target.markdup.sorted.bam"))
names(bam_paths) <- qc_pass_tf

# Check which BAMs exist
bam_exists <- file.exists(bam_paths)
if (!all(bam_exists)) {
    cat("WARNING: Missing BAM files:\n")
    cat(paste("  ", bam_paths[!bam_exists], collapse = "\n"), "\n")
    cat("Continuing with available BAMs...\n")
    bam_paths <- bam_paths[bam_exists]
}

if (length(bam_paths) == 0) {
    cat("ERROR: No BAM files found. Skipping featureCounts.\n")
    cat("Expected BAM directory:", BAM_DIR, "\n")
    cat("Will look for pre-existing count matrix...\n")

    counts_file <- paste0("counts_", TF_MARK, ".txt")
    if (!file.exists(counts_file)) {
        stop("No BAM files and no pre-existing count matrix. Cannot proceed.")
    }
} else {
    cat("Running featureCounts on", length(bam_paths), "BAM files\n")

    fc <- featureCounts(
        files             = bam_paths,
        annot.ext         = saf,
        isGTFAnnotationFile = FALSE,
        isPairedEnd       = TRUE,
        countReadPairs    = TRUE,
        nthreads          = 8
    )

    # Save count matrix
    count_matrix <- fc$counts
    colnames(count_matrix) <- names(bam_paths)
    write.csv(count_matrix, paste0("counts_", TF_MARK, ".csv"), quote = FALSE)

    # Also save featureCounts summary
    write.csv(fc$stat, paste0("counts_", TF_MARK, "_stat.csv"),
              row.names = FALSE, quote = FALSE)
    cat("Counts written to counts_", TF_MARK, ".csv\n")
}

# ==============================================================================
# STEP 4: edgeR differential binding
# ==============================================================================
cat("\n===== STEP 4: edgeR differential binding =====\n")

# Load count matrix
count_matrix <- as.matrix(read.csv(paste0("counts_", TF_MARK, ".csv"),
                                   row.names = 1, check.names = FALSE))

# Sample metadata
sample_info <- data.table(
    sample = colnames(count_matrix)
)
sample_info[, condition := sub(paste0("_", TF_MARK, "_R[0-9]+$"), "", sample)]
# Simplify condition names: 1OE -> OE, 1OS -> OS, 1YS -> YS
sample_info[, group := sub("^[0-9]+", "", condition)]

cat("Sample info:\n")
print(sample_info)

# Create DGEList
group <- factor(sample_info$group, levels = c("YS", "OS", "OE"))
dge <- DGEList(counts = count_matrix, group = group)

# Filter lowly expressed peaks
keep <- filterByExpr(dge)
cat("Peaks before filtering:", nrow(dge), "\n")
dge <- dge[keep, , keep.lib.sizes = FALSE]
cat("Peaks after filtering:", nrow(dge), "\n")

# TMM normalization
dge <- calcNormFactors(dge, method = "TMM")
cat("Norm factors:", dge$samples$norm.factors, "\n")

# Design matrix
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Estimate dispersion
dge <- estimateDisp(dge, design)

# MDS plot
pdf(file.path(FIG_DIR, paste0("MDS_", TF_MARK, ".pdf")), width = 6, height = 5)
plotMDS(dge, col = as.numeric(group), pch = 16, main = paste(EXPERIMENT, TF_MARK, "MDS"))
legend("topright", levels(group), col = 1:3, pch = 16)
dev.off()

# BCV plot
pdf(file.path(FIG_DIR, paste0("BCV_", TF_MARK, ".pdf")), width = 6, height = 5)
plotBCV(dge, main = paste(EXPERIMENT, TF_MARK, "BCV"))
dev.off()

# Fit model
fit <- glmQLFit(dge, design)

# Contrasts
con_aging <- makeContrasts(OS - YS, levels = design)
con_exer  <- makeContrasts(OE - OS, levels = design)
con_oe_ys <- makeContrasts(OE - YS, levels = design)

res_aging <- glmQLFTest(fit, contrast = con_aging)
res_exer  <- glmQLFTest(fit, contrast = con_exer)
res_oe_ys <- glmQLFTest(fit, contrast = con_oe_ys)

# Extract and save results
save_edger <- function(res, label) {
    tt <- topTags(res, n = Inf, sort.by = "none")$table
    tt$peak_id <- rownames(tt)
    outfile <- paste0("edgeR_", TF_MARK, "_", label, ".csv")
    fwrite(as.data.table(tt), outfile)
    cat(label, ": FDR<0.05 =", sum(tt$FDR < 0.05),
        "(up=", sum(tt$FDR < 0.05 & tt$logFC > 0),
        ", down=", sum(tt$FDR < 0.05 & tt$logFC < 0), ")\n")

    # Volcano plot
    pdf(file.path(FIG_DIR, paste0("volcano_", TF_MARK, "_", label, ".pdf")),
        width = 6, height = 5)
    plot(tt$logFC, -log10(tt$PValue),
         pch = 16, cex = 0.3,
         col = ifelse(tt$FDR < 0.05, "red", "grey60"),
         xlab = "logFC", ylab = "-log10(p-value)",
         main = paste(EXPERIMENT, TF_MARK, label))
    abline(v = 0, lty = 2)
    dev.off()

    return(tt)
}

tt_aging <- save_edger(res_aging, "OS_vs_YS")
tt_exer  <- save_edger(res_exer,  "OE_vs_OS")
tt_oe_ys <- save_edger(res_oe_ys, "OE_vs_YS")

# ==============================================================================
# STEP 5: Reversal analysis
# ==============================================================================
cat("\n===== STEP 5: Reversal analysis =====\n")

run_reversal <- function(peaks_to_use, peakset_label) {
    cat("\n--- Reversal analysis:", peakset_label, "---\n")

    # Subset edgeR results to peaks in this set
    aging <- as.data.table(tt_aging)[peak_id %in% peaks_to_use]
    exer  <- as.data.table(tt_exer)[peak_id %in% peaks_to_use]
    oe_ys <- as.data.table(tt_oe_ys)[peak_id %in% peaks_to_use]

    # Merge
    merged <- merge(
        aging[, .(peak_id, logFC_aging = logFC, PValue_aging = PValue, FDR_aging = FDR)],
        exer[, .(peak_id, logFC_exer = logFC, PValue_exer = PValue, FDR_exer = FDR)],
        by = "peak_id"
    )
    merged <- merge(merged,
        oe_ys[, .(peak_id, logFC_oe_ys = logFC, PValue_oe_ys = PValue, FDR_oe_ys = FDR)],
        by = "peak_id"
    )

    cat("  Peaks in set:", nrow(merged), "\n")

    # Guard: skip if no peaks
    if (nrow(merged) == 0) {
        cat("  WARNING: No peaks in this set. Skipping reversal.\n")
        return(list(results = merged, n_strict = 0, n_relaxed = 0,
                    wilcox = list(p.value = NA, conf.int = c(NA, NA)),
                    median_ri = NA))
    }

    # --- 5a. Intersection approach ---

    # Strict: aging FDR<0.05 AND exercise FDR<0.05 AND opposite direction
    merged[, reversal_strict := (FDR_aging < 0.05 & FDR_exer < 0.05 &
                                  sign(logFC_aging) != sign(logFC_exer))]
    # Subdivide: gain-then-loss vs loss-then-gain
    merged[, reversal_strict_gain_loss := (reversal_strict &
                                            logFC_aging > 0 & logFC_exer < 0)]
    merged[, reversal_strict_loss_gain := (reversal_strict &
                                            logFC_aging < 0 & logFC_exer > 0)]

    # Relaxed: aging FDR<0.05 AND exercise nominal p<0.05 AND opposite direction
    merged[, reversal_relaxed := (FDR_aging < 0.05 & PValue_exer < 0.05 &
                                   sign(logFC_aging) != sign(logFC_exer))]
    merged[, reversal_relaxed_gain_loss := (reversal_relaxed &
                                             logFC_aging > 0 & logFC_exer < 0)]
    merged[, reversal_relaxed_loss_gain := (reversal_relaxed &
                                             logFC_aging < 0 & logFC_exer > 0)]

    cat("  Strict reversal: ", sum(merged$reversal_strict),
        " (gain->loss:", sum(merged$reversal_strict_gain_loss),
        ", loss->gain:", sum(merged$reversal_strict_loss_gain), ")\n")
    cat("  Relaxed reversal:", sum(merged$reversal_relaxed),
        " (gain->loss:", sum(merged$reversal_relaxed_gain_loss),
        ", loss->gain:", sum(merged$reversal_relaxed_loss_gain), ")\n")

    # --- 5b. Continuous reversal index ---
    # For peaks with nominally significant aging effect (p < 0.05)
    aging_sig <- merged[PValue_aging < 0.05]
    cat("  Peaks with nominal aging effect (p<0.05):", nrow(aging_sig), "\n")

    if (nrow(aging_sig) > 0) {
        aging_sig[, reversal_index := 1 - logFC_oe_ys / logFC_aging]

        # Winsorize extreme values for display
        aging_sig[, ri_winsor := pmax(pmin(reversal_index, 3), -2)]

        # Wilcoxon test: is median different from 0?
        wt <- wilcox.test(aging_sig$reversal_index, mu = 0, conf.int = TRUE)
        cat("  Reversal index: median =", median(aging_sig$reversal_index, na.rm = TRUE), "\n")
        cat("  Wilcoxon p-value:", wt$p.value, "\n")
        cat("  95% CI:", wt$conf.int, "\n")

        # Merge reversal index back
        merged <- merge(merged, aging_sig[, .(peak_id, reversal_index)],
                        by = "peak_id", all.x = TRUE)
    } else {
        merged[, reversal_index := NA_real_]
        wt <- list(p.value = NA, conf.int = c(NA, NA),
                   estimate = NA)
    }

    # --- 5c. Plots ---

    # Scatter: logFC(aging) vs logFC(exercise)
    pdf(file.path(FIG_DIR, paste0("reversal_scatter_", TF_MARK, "_", peakset_label, ".pdf")),
        width = 6, height = 5)
    cor_val <- cor(merged$logFC_aging, merged$logFC_exer, use = "complete.obs")
    plot(merged$logFC_aging, merged$logFC_exer,
         pch = 16, cex = 0.3, col = "grey40",
         xlab = "logFC (OS vs YS) [aging]",
         ylab = "logFC (OE vs OS) [exercise]",
         main = paste0(EXPERIMENT, " ", TF_MARK, " -", peakset_label,
                       "\nr = ", round(cor_val, 3)))
    abline(h = 0, v = 0, lty = 2, col = "grey70")
    # Highlight reversed peaks (relaxed)
    rev_pts <- merged[reversal_relaxed == TRUE]
    if (nrow(rev_pts) > 0) {
        points(rev_pts$logFC_aging, rev_pts$logFC_exer,
               pch = 16, cex = 0.5, col = "red")
    }
    legend("topright",
           c(paste("All peaks (n=", nrow(merged), ")"),
             paste("Reversed (n=", nrow(rev_pts), ")")),
           col = c("grey40", "red"), pch = 16, cex = 0.8)
    dev.off()

    # Reversal index histogram
    if (nrow(aging_sig) > 0) {
        pdf(file.path(FIG_DIR, paste0("reversal_index_hist_", TF_MARK, "_", peakset_label, ".pdf")),
            width = 6, height = 5)
        hist(aging_sig$ri_winsor, breaks = 50,
             main = paste0(EXPERIMENT, " ", TF_MARK, " -", peakset_label,
                           "\nReversal index (nominal aging p<0.05)"),
             xlab = "Reversal index (0=none, 1=full)",
             col = "steelblue", border = "white")
        abline(v = 0, lty = 2, col = "red", lwd = 2)
        abline(v = 1, lty = 2, col = "darkgreen", lwd = 2)
        abline(v = median(aging_sig$reversal_index, na.rm = TRUE),
               lty = 1, col = "orange", lwd = 2)
        legend("topright",
               c("No reversal (0)", "Full reversal (1)",
                 paste0("Median (", round(median(aging_sig$reversal_index, na.rm = TRUE), 3), ")")),
               col = c("red", "darkgreen", "orange"), lty = c(2, 2, 1), lwd = 2, cex = 0.8)
        dev.off()
    }

    # Save reversal results
    outfile <- paste0("reversal_", TF_MARK, "_", peakset_label, ".csv")
    fwrite(merged, outfile)
    cat("  Written:", outfile, "\n")

    return(list(
        results = merged,
        n_strict = sum(merged$reversal_strict, na.rm = TRUE),
        n_relaxed = sum(merged$reversal_relaxed, na.rm = TRUE),
        wilcox = wt,
        median_ri = if (nrow(aging_sig) > 0) median(aging_sig$reversal_index, na.rm = TRUE) else NA
    ))
}

# Get peak IDs from edgeR results (after filterByExpr)
all_peak_ids <- rownames(dge$counts)

# H3K27ac overlap peak IDs
h3k_peak_ids <- h3k_overlap_bed$name

# Run reversal on all peaks
rev_all <- run_reversal(all_peak_ids, "all_peaks")

# Run reversal on H3K27ac-overlapping peaks
rev_h3k <- run_reversal(intersect(all_peak_ids, h3k_peak_ids), "at_H3K27ac")

# ==============================================================================
# STEP 6: AP-1 motif annotation
# ==============================================================================
cat("\n===== STEP 6: AP-1 motif annotation =====\n")

# Look for the master peaks FIMO hits file
motif_files <- list.files(MOTIF_DIR, pattern = "master_peaks.*parsed.*\\.gz$",
                          full.names = TRUE)
# Also try for C-FOS specific
motif_files_tf <- list.files(MOTIF_DIR,
    pattern = paste0(TF_MARK, ".*parsed.*\\.gz$"), full.names = TRUE)

cat("Available motif files:\n")
cat(paste(" ", c(motif_files, motif_files_tf)), sep = "\n")

# Use the C-FOS master peaks file if available
motif_file <- motif_files_tf[grep("master_peaks", motif_files_tf)]
if (length(motif_file) == 0) {
    motif_file <- motif_files_tf[1]
}
if (length(motif_file) == 0) {
    motif_file <- motif_files[1]
}

if (length(motif_file) > 0 && !is.na(motif_file)) {
    cat("Using motif file:", motif_file, "\n")

    motif_hits <- fread(cmd = paste("zcat", motif_file))
    cat("Motif hits loaded:", nrow(motif_hits), "\n")

    # The motif hits have peak_id, motif_start0, motif_end0
    # We need to check which NEW consensus peaks contain a motif hit

    # Get motif hit genomic coordinates
    # motif_hits columns: peak_id, window_chr, window_start0, window_end0,
    #                     offset, motif_start0, motif_end0, sequence, motif_name,
    #                     strand, motif_score
    motif_gr <- unique(motif_hits[, .(chr = window_chr,
                                       start = motif_start0,
                                       end = motif_end0)])

    # Write motif hits as BED for bedtools overlap
    motif_bed_file <- "ap1_motif_hits.bed"
    fwrite(motif_gr[, .(chr, start, end)], motif_bed_file,
           sep = "\t", col.names = FALSE)
    system2("sort", args = c("-k1,1", "-k2,2n", motif_bed_file),
            stdout = paste0(motif_bed_file, ".sorted"))
    file.rename(paste0(motif_bed_file, ".sorted"), motif_bed_file)

    # Overlap consensus peaks with motif hits (any overlap)
    motif_overlap_file <- paste0("consensus_peaks_", TF_MARK, "_motif_overlap.bed")
    system2("bedtools", args = c("intersect",
        "-a", paste0("consensus_peaks_", TF_MARK, ".bed"),
        "-b", motif_bed_file,
        "-u"),
        stdout = motif_overlap_file)

    motif_peaks <- fread(motif_overlap_file, header = FALSE,
                         col.names = c("chr","start","end","name","score","strand"))
    cat("Consensus peaks with AP-1 motif:", nrow(motif_peaks),
        "out of", nrow(consensus_bed), "\n")

    # Save motif annotation
    motif_status <- data.table(
        peak_id = consensus_bed$name,
        has_ap1_motif = consensus_bed$name %in% motif_peaks$name
    )
    fwrite(motif_status, paste0("motif_annotation_", TF_MARK, ".csv"))

    # Run reversal on motif+ subset
    motif_peak_ids <- motif_peaks$name
    rev_motif <- run_reversal(intersect(all_peak_ids, motif_peak_ids), "motif_positive")

} else {
    cat("WARNING: No motif hit file found. Skipping motif annotation.\n")
    cat("Expected files in:", MOTIF_DIR, "\n")
    motif_status <- data.table(peak_id = consensus_bed$name, has_ap1_motif = NA)
    rev_motif <- NULL
}

# ==============================================================================
# STEP 7: Master results table
# ==============================================================================
cat("\n===== STEP 7: Master results table =====\n")

# Start with edgeR results
master <- merge(
    as.data.table(tt_aging)[, .(peak_id,
        logFC_OS_vs_YS = logFC, PValue_OS_vs_YS = PValue, FDR_OS_vs_YS = FDR)],
    as.data.table(tt_exer)[, .(peak_id,
        logFC_OE_vs_OS = logFC, PValue_OE_vs_OS = PValue, FDR_OE_vs_OS = FDR)],
    by = "peak_id"
)
master <- merge(master,
    as.data.table(tt_oe_ys)[, .(peak_id,
        logFC_OE_vs_YS = logFC, PValue_OE_vs_YS = PValue, FDR_OE_vs_YS = FDR)],
    by = "peak_id"
)

# Add coordinates from peak_id
master[, c("chr", "coords") := tstrsplit(peak_id, ":", fixed = TRUE)]
master[, c("start", "end") := tstrsplit(coords, "-", fixed = TRUE, type.convert = TRUE)]
master[, coords := NULL]

# Add support info
support <- fread(paste0("consensus_support_", TF_MARK, ".csv"))
support[, peak_id := paste0(chr, ":", start, "-", end)]
master <- merge(master,
    support[, .(peak_id, n_OE, n_OS, n_YS, n_peaks)],
    by = "peak_id", all.x = TRUE
)

# Add reversal info from the all-peaks analysis
rev_results <- rev_all$results
master <- merge(master,
    rev_results[, .(peak_id, reversal_strict, reversal_relaxed, reversal_index)],
    by = "peak_id", all.x = TRUE
)

# Add motif status
master <- merge(master, motif_status, by = "peak_id", all.x = TRUE)

# Add H3K27ac overlap status
master[, at_H3K27ac := peak_id %in% h3k_peak_ids]

# Order columns
setcolorder(master, c("peak_id", "chr", "start", "end",
    "n_OE", "n_OS", "n_YS", "n_peaks",
    "logFC_OS_vs_YS", "PValue_OS_vs_YS", "FDR_OS_vs_YS",
    "logFC_OE_vs_OS", "PValue_OE_vs_OS", "FDR_OE_vs_OS",
    "logFC_OE_vs_YS", "PValue_OE_vs_YS", "FDR_OE_vs_YS",
    "reversal_strict", "reversal_relaxed", "reversal_index",
    "has_ap1_motif", "at_H3K27ac"))

fwrite(master[order(chr, start)], paste0("master_results_", TF_MARK, ".csv"))
cat("Master results:", nrow(master), "peaks\n")

# ==============================================================================
# Summary
# ==============================================================================
cat("\n")
cat("=============================================================\n")
cat("  SUMMARY: ", EXPERIMENT, " ", TF_MARK, " Reversal Analysis\n")
cat("=============================================================\n")
cat("QC-passing TF samples:       ", length(qc_pass_tf), "\n")
cat("Consensus TF peaks:          ", nrow(consensus_bed), "\n")
cat("After edgeR filterByExpr:    ", length(all_peak_ids), "\n")
cat("H3K27ac consensus peaks:     ", nrow(fread(paste0("consensus_peaks_H3K27ac.bed"), header=FALSE)), "\n")
cat("TF at H3K27ac:               ", nrow(h3k_overlap_bed), "\n")
if (exists("motif_peaks")) {
    cat("TF with AP-1 motif:          ", nrow(motif_peaks), "\n")
}
cat("\n--- Aging (OS vs YS) ---\n")
cat("  FDR<0.05: ", sum(tt_aging$FDR < 0.05),
    " (up=", sum(tt_aging$FDR < 0.05 & tt_aging$logFC > 0),
    ", down=", sum(tt_aging$FDR < 0.05 & tt_aging$logFC < 0), ")\n")
cat("\n--- Exercise (OE vs OS) ---\n")
cat("  FDR<0.05: ", sum(tt_exer$FDR < 0.05),
    " (up=", sum(tt_exer$FDR < 0.05 & tt_exer$logFC > 0),
    ", down=", sum(tt_exer$FDR < 0.05 & tt_exer$logFC < 0), ")\n")
cat("\n--- Reversal (all peaks) ---\n")
cat("  Strict:  ", rev_all$n_strict, "\n")
cat("  Relaxed: ", rev_all$n_relaxed, "\n")
cat("  Median reversal index: ", rev_all$median_ri, "\n")
cat("  Wilcoxon p:", rev_all$wilcox$p.value, "\n")
cat("\n--- Reversal (at H3K27ac) ---\n")
cat("  Strict:  ", rev_h3k$n_strict, "\n")
cat("  Relaxed: ", rev_h3k$n_relaxed, "\n")
cat("  Median reversal index: ", rev_h3k$median_ri, "\n")
if (!is.null(rev_motif)) {
    cat("\n--- Reversal (motif+) ---\n")
    cat("  Strict:  ", rev_motif$n_strict, "\n")
    cat("  Relaxed: ", rev_motif$n_relaxed, "\n")
    cat("  Median reversal index: ", rev_motif$median_ri, "\n")
}
cat("\n=============================================================\n")
cat("Done! Output files in:", ANALYSIS_DIR, "\n")

sessionInfo()
