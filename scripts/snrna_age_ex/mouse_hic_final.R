# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Link AP-1 regulatory elements to genes using Cicero and mouse Hi-C context.
# Steps: 1) Load ATAC/Hi-C resources; 2) calculate co-accessibility and loop overlaps; 3) export enhancer-gene linkage summaries.
# =======================================
# Source reformatted from snRNA_AgeEx/Mouse_HiC_Final.ipynb.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

# Redo_Find Result

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# Load and Setup
# ============================================================

TOBIAS_Peak_level <- read.csv(
  publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Output/Tobias_df_Classified_Full.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

outdir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-05-11.HiC_Cm/AP1_TOBIAS_site_sets")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Helper for BED writing
write_bed <- function(dt, file) {
  if (nrow(dt) == 0) return(invisible(NULL))
  bed <- dt[, .(chr, start, end, name = paste(site_id, cell, sep = "|"), score = 0L, strand = ".")]
  fwrite(unique(bed), file, sep = "\t", col.names = FALSE)
}

# ============================================================
# Processing
# ============================================================

dt <- as.data.table(TOBIAS_Peak_level)
dt[, cell := CellType]
dt[, chr := as.character(seqnames)]
dt[, start := as.integer(start)]
dt[, end := as.integer(end)]

# Keep valid standard intervals
standard_chroms <- paste0("chr", c(1:19, "X", "Y"))
dt <- dt[chr %in% standard_chroms & !is.na(start) & end > start]

# Restrict to Fos in Cardiomyocytes
dt_fos_cm <- dt[motif == "Fos_Homo_sapiens_DBD_1.00_M08069_2.00" & CellType == "Cardiomyocyte"]
if (nrow(dt_fos_cm) == 0) stop("No data retained.")

dt_fos_cm[, site_id := paste(chr, start, end, sep = "_")]

# Determine binding status per row
dt_fos_cm[, `:=`(
  bound_YS  = !is.na(in_YS) & (in_YS %in% c("TRUE", "T", "1", "YES", "Y", "BOUND") | in_YS == TRUE),
  bound_MAS = !is.na(in_MAS) & (in_MAS %in% c("TRUE", "T", "1", "YES", "Y", "BOUND") | in_MAS == TRUE),
  bound_MAE = !is.na(in_MAE) & (in_MAE %in% c("TRUE", "T", "1", "YES", "Y", "BOUND") | in_MAE == TRUE)
)]

# Collapse to unique sites (any bound row = bound site)
site_classes <- dt_fos_cm[, .(
  cell,
  in_YS_any = any(bound_YS, na.rm = TRUE),
  in_MAS_any = any(bound_MAS, na.rm = TRUE),
  in_MAE_any = any(bound_MAE, na.rm = TRUE)
), by = .(chr, start, end, site_id)]

# ============================================================
# Define Classes Post-Collapse
# ============================================================

# Age Gained/Lost
site_classes[, age_gained := !in_YS_any & in_MAS_any]
site_classes[, age_lost   := in_YS_any & !in_MAS_any]

# Ex Gained/Lost
site_classes[, ex_gained  := !in_MAS_any & in_MAE_any]
site_classes[, ex_lost    := in_MAS_any & !in_MAE_any]

# Stable Across All
site_classes[, stable_all := in_YS_any & in_MAS_any & in_MAE_any]

# ============================================================
# Write Outputs
# ============================================================

write_bed(site_classes[age_gained == TRUE], file.path(outdir, "AP1_FOS_age_gained.bed"))
write_bed(site_classes[age_lost == TRUE],   file.path(outdir, "AP1_FOS_age_lost.bed"))
write_bed(site_classes[ex_gained == TRUE],  file.path(outdir, "AP1_FOS_ex_gained.bed"))
write_bed(site_classes[ex_lost == TRUE],    file.path(outdir, "AP1_FOS_ex_lost.bed"))
write_bed(site_classes[stable_all == TRUE], file.path(outdir, "AP1_FOS_stable_across_all.bed"))

# Save main table
fwrite(site_classes, file.path(outdir, "AP1_FOS_site_classes.tsv"), sep = "\t")

cat("Processing complete. Files saved to:", outdir, "\n")



file.path(outdir)

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
})

# ============================================================
# Build Hi-C loop anchor and hotspot BEDs from HiCCUPS BEDPE
# ============================================================

base_dir <- publication_path("HEART_SNRNA_ROOT", "2026-05-11.HiC_Cm/")

# EDIT this if your 5 kb HiCCUPS output folder has a different name
loops_file <- file.path(base_dir, "ctrl_hiccups_10kb_GPU/merged_loops.bedpe")

# This matches the path your enrichment script expected
outdir <- file.path(base_dir, "AP1_HiC_hotspot_overlap")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

standard_chroms <- paste0("chr", c(1:19, "X"))

# ============================================================
# Read loops
# ============================================================

if (!file.exists(loops_file)) {
  stop("Cannot find loops_file: ", loops_file)
}

loops <- fread(loops_file, header = TRUE, fill = TRUE)

# Fix first column if named "#chr1"
setnames(loops, old = names(loops)[1], new = sub("^#", "", names(loops)[1]))

required_cols <- c(
  "chr1", "x1", "x2",
  "chr2", "y1", "y2",
  "observed",
  "expectedBL",
  "expectedDonut",
  "expectedH",
  "expectedV"
)

missing_cols <- setdiff(required_cols, names(loops))

if (length(missing_cols) > 0) {
  stop("Missing columns from HiCCUPS BEDPE: ", paste(missing_cols, collapse = ", "))
}

loops <- loops[
  chr1 %in% standard_chroms &
    chr2 %in% standard_chroms
]

loops[, `:=`(
  x1 = as.integer(x1),
  x2 = as.integer(x2),
  y1 = as.integer(y1),
  y2 = as.integer(y2),
  observed = as.numeric(observed),
  expectedBL = as.numeric(expectedBL),
  expectedDonut = as.numeric(expectedDonut),
  expectedH = as.numeric(expectedH),
  expectedV = as.numeric(expectedV)
)]

loops <- loops[
  !is.na(x1) & !is.na(x2) &
    !is.na(y1) & !is.na(y2) &
    x2 > x1 &
    y2 > y1
]

# ============================================================
# Define loop strength
# Conservative O/E score across all local HiCCUPS background masks
# ============================================================

eps <- 1e-6

loops[, oeBL := observed / pmax(expectedBL, eps)]
loops[, oeDonut := observed / pmax(expectedDonut, eps)]
loops[, oeH := observed / pmax(expectedH, eps)]
loops[, oeV := observed / pmax(expectedV, eps)]

loops[, loop_strength := pmin(
  log2(oeBL),
  log2(oeDonut),
  log2(oeH),
  log2(oeV)
)]

# Keep positively enriched loops only
#loops <- loops[loop_strength > 0]
loops[, loop_strength := matrixStats::rowMedians(
    as.matrix(.SD)
), .SDcols = c("oeBL", "oeDonut", "oeH", "oeV")]

# Keep loops where the median enrichment is positive
loops <- loops[loop_strength > 0]

#loops[, n_positive := (oeBL > 1) + (oeDonut > 1) + (oeH > 1) + (oeV > 1)]

# Keep loops that are enriched in at least 3 out of 4 models
#loops <- loops[n_positive >= 3]

# Define strength as the mean for the remaining set
cols <- c("oeBL", "oeDonut", "oeH", "oeV")
loops[, loop_strength := rowMeans(log2(.SD)), .SDcols = cols]
loops[, loop_id := paste0("loop_", seq_len(.N))]

cat("Loops retained:", nrow(loops), "\n")

# ============================================================
# Convert loops to anchors
# ============================================================

anchor1 <- loops[, .(
  chr = chr1,
  start = x1,
  end = x2,
  loop_id,
  loop_strength,
  observed,
  expectedDonut,
  numCollapsed
)]

anchor2 <- loops[, .(
  chr = chr2,
  start = y1,
  end = y2,
  loop_id,
  loop_strength,
  observed,
  expectedDonut,
  numCollapsed
)]

anchors <- rbind(anchor1, anchor2, use.names = TRUE, fill = TRUE)
anchors <- anchors[end > start]

anchor_gr <- GRanges(
  seqnames = anchors$chr,
  ranges = IRanges(start = anchors$start + 1L, end = anchors$end)
)

mcols(anchor_gr)$loop_id <- anchors$loop_id
mcols(anchor_gr)$loop_strength <- anchors$loop_strength
mcols(anchor_gr)$observed <- anchors$observed
mcols(anchor_gr)$expectedDonut <- anchors$expectedDonut
mcols(anchor_gr)$numCollapsed <- anchors$numCollapsed

merged_anchor_gr <- reduce(anchor_gr, with.revmap = TRUE, ignore.strand = TRUE)
revmap <- mcols(merged_anchor_gr)$revmap


setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2026-05-11.HiC_Cm/AP1_TOBIAS_site_sets"))



all_loop_anchors <- anchor_gr




library(GenomicRanges)
library(data.table)
outdir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-05-11.HiC_Cm/AP1_TOBIAS_site_sets")
# 1. Load your BED files
stable_bed <- fread(file.path(outdir, "AP1_FOS_stable_across_all.bed"))
gained_bed <- fread(file.path(outdir, "AP1_FOS_age_gained.bed"))

# 2. Convert BEDs to GRanges
# We define a function for consistency
bed_to_gr <- function(dt) {
  GRanges(seqnames = dt$V1, ranges = IRanges(start = dt$V2, end = dt$V3))
}

gr_stable <- bed_to_gr(stable_bed)
gr_gained <- bed_to_gr(gained_bed)

# # 3. Flatten loops into a single GRanges object (or check both anchors)
# # Usually, we check if either anchor overlaps the motif

all_loop_anchors <- anchor_gr
# 1. Ensure you are using a unique set of loops (anchors)
# If your loops are stored as two sets of coordinates, 
# combine them to create a representative set of loop boundaries.
all_loops <- all_loop_anchors
total_loops <- length(all_loops)

# 2. Get hit/miss vectors for both sets
overlaps_stable <- countOverlaps(all_loops, gr_stable) > 0
overlaps_gained <- countOverlaps(all_loops, gr_gained) > 0

# 3. Create the contingency table
# Stable row: [Overlap Yes, Overlap No]
# Gained row: [Overlap Yes, Overlap No]

A <- sum(overlaps_stable)
C <- total_loops - A
B <- sum(overlaps_gained)
D <- total_loops - B

contingency_table <- matrix(
  c(A, C, B, D), 
  nrow = 2, 
  byrow = TRUE,
  dimnames = list(
    Set = c("Stable", "Gained"),
    Overlap = c("Yes", "No")
  )
)

# 4. Perform Fisher's Exact Test
fisher_result <- fisher.test(contingency_table, alternative = "two.sided")

# 5. Output results
print("Contingency Table:")
print(contingency_table)

print("Fisher's Exact Test Results:")
print(fisher_result)

library(GenomicRanges)
library(data.table)
library(ggplot2)

outdir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-05-11.HiC_Cm/AP1_TOBIAS_site_sets")

# Helper: Load BED as GRanges
load_bed_as_gr <- function(filename) {
  dt <- fread(file.path(outdir, filename))
  GRanges(seqnames = dt$V1, ranges = IRanges(start = dt$V2, end = dt$V3))
}

# Define your files
files <- c(
  stable = "AP1_FOS_stable_across_all.bed",
  age_gained = "AP1_FOS_age_gained.bed",
  age_lost = "AP1_FOS_age_lost.bed",
  ex_gained = "AP1_FOS_ex_gained.bed",
  ex_lost = "AP1_FOS_ex_lost.bed"
)

# Load all sets into a named list
motif_sets <- lapply(files, load_bed_as_gr)

# Calculate overlap boolean vectors for all loops
# (Assumes 'all_loop_anchors' is already defined in your environment)
overlap_list <- lapply(motif_sets, function(gr) countOverlaps(all_loop_anchors, gr) > 0)
total_loops <- length(all_loop_anchors)

# Define pairs to test
comparisons <- list(
  c("age_gained", "stable"),
  c("age_lost", "stable"),
  c("ex_gained", "stable"),
  c("ex_lost", "stable"),
  c("age_gained", "age_lost"),
  c("ex_gained", "ex_lost")
)

# Perform Fisher's Test for each pair
results <- rbindlist(lapply(comparisons, function(pair) {
  setA <- overlap_list[[pair[1]]]
  setB <- overlap_list[[pair[2]]]
  
  # Contingency Table: [SetA_Yes, SetA_No] vs [SetB_Yes, SetB_No]
  # Note: A loop can technically overlap both, but we are testing the 
  # association of the loop set with Motif A vs Motif B.
  
  mat <- matrix(c(
    sum(setA), total_loops - sum(setA),
    sum(setB), total_loops - sum(setB)
  ), nrow = 2, byrow = TRUE)
  
  test <- fisher.test(mat)
  
  data.table(
    Comparison = paste(pair[1], "vs", pair[2]),
    OR = test$estimate,
    CI_low = test$conf.int[1],
    CI_high = test$conf.int[2],
    P_val = test$p.value
  )
}))

print(results)

# 1. Combine dynamic sets: TRUE if any dynamic site is hit
overlap_list$any_dynamic <- (
  overlap_list$age_gained | 
  overlap_list$age_lost   | 
  overlap_list$ex_gained  | 
  overlap_list$ex_lost
)

# 2. Add the new comparison to your list
comparisons <- list(
  c("age_gained", "stable"),
  c("age_lost", "stable"),
  c("ex_gained", "stable"),
  c("ex_lost", "stable"),
  c("any_dynamic", "stable"), # <--- The requested addition
  c("age_gained", "age_lost"),
  c("ex_gained", "ex_lost")
)

# 3. Now the rest of your loop will automatically calculate this result
results <- rbindlist(lapply(comparisons, function(pair) {
  setA <- overlap_list[[pair[1]]]
  setB <- overlap_list[[pair[2]]]
  
  mat <- matrix(c(
    sum(setA), total_loops - sum(setA),
    sum(setB), total_loops - sum(setB)
  ), nrow = 2, byrow = TRUE)
  
  test <- fisher.test(mat)
  
  data.table(
    Comparison = paste(pair[1], "vs", pair[2]),
    OR = test$estimate,
    CI_low = test$conf.int[1],
    CI_high = test$conf.int[2],
    P_val = test$p.value
  )
}))
# 1. Apply multiple testing correction (Benjamini-Hochberg)
results[, FDR := p.adjust(P_val, method = "BH")]

# 2. Define significance category
results[, Sig_Level := ifelse(FDR < 0.001, "FDR < 0.001", "NS")]

# 3. Factor for plotting order (ensures consistent legend)
results[, Sig_Level := factor(Sig_Level, levels = c("FDR < 0.001", "NS"))]

sessionInfo()

# 1. Transform the data to log2 scale
results[, log2OR := log2(OR)]
results[, log2CI_low := log2(CI_low)]
results[, log2CI_high := log2(CI_high)]

# 2. Updated Plot
ggplot(results, aes(x = log2OR, y = Comparison, color = Sig_Level)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = log2CI_low, xmax = log2CI_high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  scale_color_manual(
    values = c("FDR < 0.001" = "#E41A1C", "NS" = "#999999"),
    name = "Significance"
  ) +
  theme_classic() +
  labs(
    title = "Enrichment Log-Odds (FDR Corrected)",
    subtitle = "Values > 0 are enriched; values < 0 are depleted",
    x = "log2(Odds Ratio)",
    y = "Comparison Pair"
  ) +
  theme(legend.position = "bottom")

ggsave(publication_path("PAPER_OUTPUT_ROOT", "Paper_Tables/HiC_mouse/MouseHIC.pdf"))


