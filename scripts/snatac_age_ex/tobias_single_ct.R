# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Run cell-type-specific TOBIAS footprinting and differential binding.
# Steps: 1) Locate peaks/BAMs and motif resources; 2) execute TOBIAS steps; 3) write footprint and motif outputs.
# =======================================
# Source reformatted from snATAC_AgeEx/tobias_single_ct.R.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

library(tidyverse)
library(stringr)
library(parallel)

# --- 0. Read Command Line Argument ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Error: No cell type specified.")
target_ct <- args[1]

message(paste(">>> STARTING PIPELINE FOR:", target_ct, "<<<"))

# --- 1. Configuration ---
base_dir     <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Output/")
peak_dir     <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_RefinedPeaks/")
bam_dir      <- publication_path("HEART_SNRNA_ROOT", "analysis/2024-02-07.ChromBp/pooled/bams/")
blacklist_file <- publication_path("TOBIAS_RESOURCES_DIR", "mm10-blacklist.v2.bed")
genome         <- publication_path("MM10_RESOURCES_DIR", "mm10.fa")
motifs         <- publication_path("TOBIAS_RESOURCES_DIR", "all_motifs.pfm")

# Resource usage: 4 conditions * 4 cores = 16 cores total
cores_parallel_jobs <- 4
cores_tobias_internal <- 4 

if(!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)

# --- 2. File Discovery ---
peak_pattern <- paste0("^", target_ct, "_consensus.bed$")
peak_file <- list.files(peak_dir, pattern = peak_pattern, full.names = TRUE)

if (length(peak_file) != 1) stop(paste("Error: Peak file missing for", target_ct))

all_bams <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE, recursive = TRUE)
all_bams <- all_bams[grep('cellID', all_bams, invert = TRUE)]
bam_matches <- grep(paste0("/", target_ct, "_"), all_bams, value = TRUE)

if (length(bam_matches) == 0) stop(paste("No BAM files found for", target_ct))

work_df <- tibble(bam_path = bam_matches) %>%
  mutate(
    filename = basename(bam_path),
    condition = str_extract(filename, "(MAE|MAS|YE|YS)"),
    check_ct = str_remove(str_extract(filename, "^(.*?)_(?=MAE|MAS|YE|YS)"), "_$")
  ) %>%
  filter(check_ct == target_ct & !is.na(condition))

# --- 3. Execution ---
# A. ATACorrect & Footprinting
mclapply(seq_len(nrow(work_df)), function(i) {
  row <- work_df[i, ]
  out_subdir <- file.path(base_dir, paste0(target_ct, "_", row$condition))
  if(!dir.exists(out_subdir)) dir.create(out_subdir, recursive = TRUE)
  
  # Correct
  cmd_correct <- paste('TOBIAS ATACorrect', '--bam', row$bam_path, '--genome', genome,
                       '--peaks', peak_file, '--blacklist', blacklist_file,
                       '--outdir', out_subdir, '--cores', cores_tobias_internal)
  system(cmd_correct)
  
  # Footprint
  bam_base <- tools::file_path_sans_ext(basename(row$bam_path))
  cmd_footprint <- paste('TOBIAS FootprintScores', 
                         '--signal', file.path(out_subdir, paste0(bam_base, "_corrected.bw")),
                         '--regions', peak_file,
                         '--output', file.path(out_subdir, paste0(target_ct, "_", row$condition, "_footprints.bw")),
                         '--cores', cores_tobias_internal)
  system(cmd_footprint)
}, mc.cores = cores_parallel_jobs)

# B. BINDetect
get_bw <- function(c) file.path(base_dir, paste0(target_ct, "_", c), paste0(target_ct, "_", c, "_footprints.bw"))
signals <- c(get_bw("MAS"), get_bw("MAE"), get_bw("YS"), get_bw("YE"))
labels  <- c("MAS", "MAE", "YS", "YE")
valid <- which(file.exists(signals))

if (length(valid) >= 2) {
  cmd_bin <- paste('TOBIAS BINDetect', '--motifs', motifs,
                   '--signals', paste(signals[valid], collapse=" "),
                   '--genome', genome, '--peaks', peak_file,
                   '--outdir', file.path(base_dir, 'BINDetect', paste0(target_ct, "_comp_all")),
                   '--cond_names', paste(labels[valid], collapse=" "),
                   '--cores', 16) # Use all 16 allocated cores
  system(cmd_bin)
}
