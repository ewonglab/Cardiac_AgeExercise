#!/usr/bin/env Rscript
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Calculate GC content for supplied NMF-module peak assignments.
# Steps: 1) Load module peaks; 2) obtain mm10 sequence GC percentages; 3) aggregate module summaries for downstream figures.
# =======================================
# Source reformatted from emily_plots/GC/08b_module_gc_content.R.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))
# ==============================================================================
# 08b_module_gc_content.R
# Reproduces outputs_tobias/module_gc_content.csv from chromatin_modules.csv.
# For each NMF module's top-500 peak set: compute per-peak GC% via
# BSgenome.Mmusculus.UCSC.mm10, then aggregate to module-level summaries.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(BSgenome.Mmusculus.UCSC.mm10)
})

# --- Paths ---
default_base <- publication_path("SIGNALCRE_ROOT")
args <- commandArgs(trailingOnly = TRUE)
base <- if (length(args) >= 1) args[1] else default_base
modules_path <- file.path(base, "outputs_pooled", "chromatin_modules.csv")
out_path     <- file.path(base, "outputs_tobias", "module_gc_content.csv")

stopifnot(file.exists(modules_path))

# --- Load module -> peak assignments (top-500 peaks per module) ---
mods <- fread(modules_path)
message("Modules loaded: ", uniqueN(mods$module_id),
        " across ", uniqueN(mods$cell_type), " cell types; ",
        nrow(mods), " module x peak rows")

# --- Parse peak coordinate (chrN-START-END) -> GRanges ---
peak_parts <- tstrsplit(mods$peak, "-", fixed = TRUE)
mods[, `:=`(
  peak_chr   = peak_parts[[1]],
  peak_start = as.integer(peak_parts[[2]]),
  peak_end   = as.integer(peak_parts[[3]])
)]
mods[, peak_width := peak_end - peak_start]

# Unique peaks (a peak can belong to multiple modules)
unique_peaks <- unique(mods[, .(peak, peak_chr, peak_start, peak_end, peak_width)])
message("Unique peaks: ", nrow(unique_peaks))

# Restrict to canonical chromosomes present in BSgenome
valid_chr <- seqnames(BSgenome.Mmusculus.UCSC.mm10)
unique_peaks <- unique_peaks[peak_chr %in% as.character(valid_chr)]

# --- GC% from BSgenome ---
peak_gr <- GRanges(seqnames = unique_peaks$peak_chr,
                   ranges = IRanges(unique_peaks$peak_start, unique_peaks$peak_end))
peak_seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, peak_gr)
gc_freq    <- letterFrequency(peak_seqs, letters = c("G", "C"))
total_freq <- letterFrequency(peak_seqs, letters = c("A", "C", "G", "T"))
unique_peaks[, gc_pct := rowSums(gc_freq) / rowSums(total_freq) * 100]

message("GC range: ", round(min(unique_peaks$gc_pct), 1), "% - ",
        round(max(unique_peaks$gc_pct), 1), "% | mean: ",
        round(mean(unique_peaks$gc_pct), 1), "%")

# --- Aggregate per module ---
mods <- merge(mods, unique_peaks[, .(peak, gc_pct)], by = "peak", all.x = TRUE)

module_gc <- mods[, .(
  mean_gc     = mean(gc_pct, na.rm = TRUE),
  median_gc   = median(gc_pct, na.rm = TRUE),
  weighted_gc = weighted.mean(gc_pct, peak_loading, na.rm = TRUE),
  sd_gc       = sd(gc_pct, na.rm = TRUE),
  n_peaks     = .N,
  mean_width  = mean(peak_width, na.rm = TRUE)
), by = .(module_id, module_label, cell_type)]

# KLF-family flag (label contains a KLF-family TF)
module_gc[, is_klf := grepl("(^|_)(Klf[0-9]+|Sp[1-9])($|_)", module_label)]

# Sort by mean_gc (desc) to match committed file's row order
setorder(module_gc, -mean_gc)

# --- Write ---
fwrite(module_gc, out_path)
message("Wrote: ", out_path, " (", nrow(module_gc), " modules)")
