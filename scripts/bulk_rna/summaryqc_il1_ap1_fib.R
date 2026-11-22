# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Summarise IL-1/AP-1 fibroblast bulk-RNA quality control and response plots.
# Steps: 1) Load count/QC inputs; 2) calculate library and expression summaries; 3) export publication plots and tables.
# =======================================
# Source reformatted from Bulk_RNA/SummaryQC-Il1_AP1_Fib(2).ipynb.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

# Summary table QC

#!/usr/bin/env Rscript
# rnaseq_qc_summary.R
# Extract key QC metrics from nf-core RNA-seq pipeline MultiQC output files.
#
# Usage:
#   Rscript rnaseq_qc_summary.R <results_dir> [output.tsv]
#
# Arguments:
#   results_dir  Path to the nf-core rnaseq results directory
#                MultiQC data files are expected under:
#                  <results_dir>/multiqc/star_salmon/multiqc_data/  (nf-core default)
#                Falls back to: multiqc/multiqc_data/, multiqc_data/, or results_dir itself
#   output.tsv   Output file path (default: qc_summary.tsv)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})
setwd(publication_path("BULK_RNA_IL1_ROOT"))
# ── Arguments ─────────────────────────────────────────────────────────────────
args <-publication_path("BULK_RNA_IL1_ROOT", "multiqc/star_salmon/multiqc_report_data")
if (length(args) < 1) {
  cat("Usage: Rscript rnaseq_qc_summary.R <results_dir> [output.tsv]\n")
  quit(status = 1)
}
results_dir <- normalizePath(args[1], mustWork = FALSE)
outfile     <- if (length(args) >= 2) args[2] else "qc_summary.tsv"

# ── Locate MultiQC data directory ─────────────────────────────────────────────
find_multiqc_dir <- function(base) {
  candidates <- c(
    file.path(base, "multiqc", "star_salmon", "multiqc_data"),
    file.path(base, "multiqc", "multiqc_data"),
    file.path(base, "multiqc_data"),
    base
  )
  for (p in candidates) {
    if (dir.exists(p) && file.exists(file.path(p, "multiqc_star.txt"))) return(p)
  }
  stop(paste("Cannot locate multiqc_data directory with multiqc_star.txt under:", base))
}

multiqc_dir <- find_multiqc_dir(results_dir)
message("Using MultiQC data from: ", multiqc_dir)

read_mqc <- function(filename) {
  path <- file.path(multiqc_dir, filename)
  if (!file.exists(path)) stop(paste("Required file not found:", path))
  read_tsv(path, show_col_types = FALSE)
}

# ── 1. STAR: canonical sample list, reads passing QC, unique mapping % ────────
star <- read_mqc("multiqc_star.txt")
samples <- star$Sample

star_metrics <- star %>%
  transmute(
    Sample,
    Reads_passing_QC_M          = total_reads / 1e6,
    Aligned_Uniquely_mapped_pct = uniquely_mapped_percent,
    star_total_reads            = total_reads   # kept for downstream ratio calculations
  )
head(star_metrics)

# ── 2. General stats: FastQC raw/trimmed per-read rows ────────────────────────
# Per-read rows are named <sample>_1 (R1) and <sample>_2 (R2)
gen <- read_mqc("multiqc_general_stats.txt")

read_names <- unlist(lapply(samples, function(s) paste0(s, c("_1", "_2"))))
gen_reads  <- gen %>%
  filter(Sample %in% read_names) %>%
  mutate(
    parent_sample = str_replace(Sample, "_[12]$", ""),
    read_num      = str_extract(Sample, "[12]$")
  )
gen_r1 <- filter(gen_reads, read_num == "1")
gen_r2 <- filter(gen_reads, read_num == "2")

col_raw_total  <- "FastQC (raw)_mqc-generalstats-fastqc_raw-total_sequences"
col_raw_dup    <- "FastQC (raw)_mqc-generalstats-fastqc_raw-percent_duplicates"
col_raw_gc     <- "FastQC (raw)_mqc-generalstats-fastqc_raw-percent_gc"
col_trim_total <- "FastQC (trimmed)_mqc-generalstats-fastqc_trimmed-total_sequences"

# Total raw reads: R1 count * 2 (paired-end), in millions
total_raw_M <- gen_r1 %>%
  transmute(
    Sample            = parent_sample,
    Total_raw_reads_M = .data[[col_raw_total]] * 2 / 1e6
  )

# Duplication %: average of R1 and R2
dup <- inner_join(
  gen_r1 %>% transmute(Sample = parent_sample, d1 = .data[[col_raw_dup]]),
  gen_r2 %>% transmute(Sample = parent_sample, d2 = .data[[col_raw_dup]]),
  by = "Sample"
) %>% transmute(Sample, FastQC_duplication_pct = (d1 + d2) / 2)

# GC content: average of R1 and R2
gc_pct <- inner_join(
  gen_r1 %>% transmute(Sample = parent_sample, g1 = .data[[col_raw_gc]]),
  gen_r2 %>% transmute(Sample = parent_sample, g2 = .data[[col_raw_gc]]),
  by = "Sample"
) %>% transmute(Sample, GC_content = (g1 + g2) / 2)

# Filter pass %: trimmed reads / raw reads * 100 (from R1 row; R1 == R2 here)
filter_pass <- gen_r1 %>%
  transmute(
    Sample          = parent_sample,
    Filter_pass_pct = .data[[col_trim_total]] / .data[[col_raw_total]] * 100
  )
                            head(filter_pass)

# ── 3. Cutadapt: adapter trimming % ──────────────────────────────────────────
cutadapt <- read_mqc("multiqc_cutadapt.txt") %>%
  filter(Sample %in% read_names) %>%
  mutate(
    parent_sample = str_replace(Sample, "_[12]$", ""),
    adapter_pct   = r_with_adapters / r_processed * 100
  )

adapter <- cutadapt %>%
  group_by(Sample = parent_sample) %>%
  summarise(Adapter_pct = mean(adapter_pct), .groups = "drop")

# ── 4. Q30 rate from per-sequence quality scores plot ─────────────────────────
# File structure: alternating header row (quality bins) / data row (counts)
# Both R1 and R2 samples are interleaved with their respective headers.
parse_q30_file <- function(path) {
  lines <- readLines(path)
  out   <- vector("list", length(lines) %/% 2)
  k     <- 0L
  i     <- 1L
  while (i + 1 <= length(lines)) {
    hdr    <- str_split(lines[i],     "\t")[[1]]
    dat    <- str_split(lines[i + 1], "\t")[[1]]
    quals  <- suppressWarnings(as.numeric(hdr[-1]))   # first field is empty
    counts <- suppressWarnings(as.numeric(dat[-1]))   # first field is sample name
    sname  <- dat[1]
    total  <- sum(counts, na.rm = TRUE)
    q30    <- if (total > 0) sum(counts[!is.na(quals) & quals >= 30], na.rm = TRUE) / total * 100
              else NA_real_
    k      <- k + 1L
    out[[k]] <- list(sample_read = sname, q30 = q30)
    i <- i + 2L
  }
  tibble(
    sample_read = sapply(out[seq_len(k)], `[[`, "sample_read"),
    q30         = sapply(out[seq_len(k)], `[[`, "q30")
  )
}

q30_path <- file.path(multiqc_dir, "mqc_fastqc_per_sequence_quality_scores_plot_1.txt")
if (!file.exists(q30_path)) {
  warning("Q30 quality scores file not found; Q30_rate set to NA.")
  q30 <- tibble(Sample = samples, Q30_rate = NA_real_)
} else {
  q30 <- parse_q30_file(q30_path) %>%
    filter(sample_read %in% read_names) %>%
    mutate(Sample = str_replace(sample_read, "_[12]$", "")) %>%
    group_by(Sample) %>%
    summarise(Q30_rate = mean(q30, na.rm = TRUE), .groups = "drop")
}

# ── 5. Featurecounts biotype: assigned reads, MT biotype % ───────────────────
biotype_raw <- read_mqc("multiqc_featurecounts_biotype_plot.txt")
biotype_num <- biotype_raw %>%
  select(-Sample) %>%
  mutate(across(everything(), as.numeric))

mt_cols <- intersect(c("Mt_tRNA", "Mt_rRNA"), names(biotype_num))
if (length(mt_cols) == 0) {
  warning("Mt_tRNA / Mt_rRNA columns not found in featurecounts biotype file.")
}

biotype_metrics <- biotype_raw %>%
  transmute(
    Sample,
    total_assigned            = rowSums(biotype_num, na.rm = TRUE),
    mt_biotype                = if (length(mt_cols) > 0)
                                  rowSums(biotype_num[, mt_cols, drop = FALSE], na.rm = TRUE)
                                else 0,
    FeatureCount_Assigned_reads_M = total_assigned / 1e6,
    Mt_pct_biotype            = mt_biotype / total_assigned * 100
  ) %>%
  select(Sample, FeatureCount_Assigned_reads_M, Mt_pct_biotype)

# Gene_Assigned_pct = featurecount assigned reads / STAR total reads * 100
gene_assigned_pct <- biotype_metrics %>%
  left_join(star_metrics %>% select(Sample, star_total_reads), by = "Sample") %>%
  transmute(
    Sample,
    FeatureCount_Assigned_reads_M,
    Mt_pct_biotype,
    Gene_Assigned_pct = FeatureCount_Assigned_reads_M * 1e6 / star_total_reads * 100
  )

# ── 6. Samtools idxstats: MT % of genome-mapped reads ────────────────────────
# Cell format: [mapped_reads, chr_length]  — extract first number per cell.
idxstats_path <- file.path(multiqc_dir, "multiqc_samtools_idxstats.txt")
if (!file.exists(idxstats_path)) {
  warning("multiqc_samtools_idxstats.txt not found; Mt_pct_genome set to NA.")
  mt_genome <- tibble(Sample = samples, Mt_pct_genome = NA_real_)
} else {
  idx_lines  <- readLines(idxstats_path)
  idx_header <- str_split(idx_lines[1], "\t")[[1]]
  chroms     <- idx_header[-1]   # chromosome names (first field is "Sample")
  mt_chr     <- intersect(c("MT", "chrM", "chrMT"), chroms)[1]

  if (is.na(mt_chr)) {
    warning("MT chromosome column not found in idxstats; Mt_pct_genome set to NA.")
    mt_genome <- tibble(Sample = samples, Mt_pct_genome = NA_real_)
  } else {
    mt_genome <- lapply(idx_lines[-1], function(line) {
      parts  <- str_split(line, "\t")[[1]]
      sname  <- parts[1]
      values <- parts[-1]
      names(values) <- chroms
      # Extract the mapped-read count (first integer inside brackets)
      mapped <- suppressWarnings(as.numeric(str_extract(values, "(?<=\\[)\\d+")))
      names(mapped) <- chroms
      mt_m  <- mapped[[mt_chr]]
      tot_m <- sum(mapped[chroms != "*"], na.rm = TRUE)
      tibble(
        Sample        = sname,
        Mt_pct_genome = if (!is.na(tot_m) && tot_m > 0) mt_m / tot_m * 100 else NA_real_
      )
    }) %>% bind_rows()
  }
}

# ── 7. Genes detected from salmon merged gene counts ─────────────────────────
results_dir <- normalizePath(publication_path("BULK_RNA_IL1_ROOT"), mustWork = FALSE)
gene_counts_path <- file.path(results_dir, "star_salmon", "salmon.merged.gene_counts.tsv")

if (file.exists(gene_counts_path)) {
  gc_dat      <- read_tsv(gene_counts_path, show_col_types = FALSE)
  sample_cols <- intersect(names(gc_dat), samples)
  if (length(sample_cols) > 0) {
    genes_detected <- gc_dat %>%
      select(all_of(sample_cols)) %>%
      summarise(across(everything(), ~ sum(. > 0, na.rm = TRUE))) %>%
      pivot_longer(everything(), names_to = "Sample", values_to = "Genes_detected")
  } else {
    message("Warning: no sample columns matched in gene counts file; Genes_detected set to NA.")
    genes_detected <- tibble(Sample = samples, Genes_detected = NA_integer_)
  }
} else {
  message("Note: gene counts file not found at\n  ", gene_counts_path,
          "\n  Genes_detected will be NA (fill in after copying the file).")
  genes_detected <- tibble(Sample = samples, Genes_detected = NA_integer_)
}

# ── 8. Derive Group from sample name ─────────────────────────────────────────
# Group = prefix before the trailing _<digits>  (e.g. IgG_1 → IgG, il1b_3 → il1b)
group_df <- tibble(
  Sample     = samples,
  Experiment = NA_character_,
  Group      = str_replace(samples, "_\\d+$", "")
)
group_df

genes_detected

mt_genome

# ── 9. Assemble final table ───────────────────────────────────────────────────
summary_tbl <- group_df %>%
  left_join(total_raw_M,                                by = "Sample") %>%
  left_join(dup,                                        by = "Sample") %>%
  left_join(q30,                                        by = "Sample") %>%
  left_join(gc_pct,                                     by = "Sample") %>%
  left_join(adapter,                                    by = "Sample") %>%
  left_join(filter_pass,                                by = "Sample") %>%
  left_join(star_metrics %>% select(-star_total_reads), by = "Sample") %>%
  left_join(gene_assigned_pct,                          by = "Sample") %>%
  left_join(genes_detected,                             by = "Sample") %>%
  left_join(mt_genome,                                  by = "Sample") %>%
  select(
    Sample, Experiment, Group,
    Total_raw_reads_M, FastQC_duplication_pct, Q30_rate, GC_content,
    Adapter_pct, Filter_pass_pct,
    Reads_passing_QC_M, Aligned_Uniquely_mapped_pct,
    FeatureCount_Assigned_reads_M, Gene_Assigned_pct, Genes_detected,
    Mt_pct_biotype, Mt_pct_genome
  ) %>%
  mutate(across(where(is.double), ~ round(., 3))) %>%
  arrange(Sample)

write_tsv(summary_tbl, outfile, na = "NA")
message(sprintf(
  "Saved QC summary: %d samples x %d columns → %s",
  nrow(summary_tbl), ncol(summary_tbl), outfile
))

summary_tbl

# Run analysis 



# Load required libraries
library(edgeR)
library(tidyverse)
setwd(publication_path("BULK_RNA_IL1_ROOT", "star_salmon"))
# 1. Read and format the count data
counts_data <- read_tsv("salmon.merged.gene_counts.tsv")

# nf-core RNA-seq outputs decimal expected counts from Salmon. edgeR requires integers.
# We also handle potential duplicate gene_names by summing their counts.
count_mat <- counts_data %>%
  dplyr::select(-gene_id) %>% 
  group_by(gene_name) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("gene_name") %>%
  as.matrix() %>%
  round() # edgeR requires integer counts

# 2. Extract metadata from column names
# Exclude the BioBanker samples from the main DE analysis
samples <- colnames(count_mat)
samples <- samples[!grepl("BioBanker", samples)]
count_mat <- count_mat[, samples]

# Create grouping factor (e.g., extracts "SP_Veh_3h" from "SP_Veh_3h_R1")
groups <- factor(sub("_R[1-3]$", "", samples))

# 3. Setup edgeR DGEList object
y <- DGEList(counts = count_mat, group = groups)

# Filter lowly expressed genes
keep <- filterByExpr(y, group = groups)
y <- y[keep, , keep.lib.sizes = FALSE]

# Calculate TMM normalization factors
y <- calcNormFactors(y)

# 4. Define Design Matrix
# A ~ 0 + group design creates a column for each group, making contrasts straightforward
design <- model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)

# Estimate dispersion
y <- estimateDisp(y, design)

# Fit the Quasi-Likelihood model (recommended for bulk RNA-seq)
fit <- glmQLFit(y, design)

# 5. Define Contrasts
# T5224 is the AP-1 inhibitor. Relevant contrasts test the stimulus, the inhibitor alone, 
# and the inhibitor's ability to block the stimulus.
my_contrasts <- makeContrasts(
  
  # --- 3 Hour Timepoint ---
  # Effect of IL-1 stimulation
  IL1_vs_Veh_3h = SP_IL1_3h - SP_Veh_3h,
  # Baseline effect of AP-1 inhibitor
  T5224_vs_Veh_3h = SP_T5224_3h - SP_Veh_3h,
  # Rescue effect: How does pre-treatment alter the IL-1 response?
  T5224preIL1_vs_IL1_3h = SP_T5224pre_IL1_3h - SP_IL1_3h,
  T5224preIL1_vs_Veh_3h = SP_T5224pre_IL1_3h - SP_Veh_3h,
  T5224preIL1_vs_Veh_6h = SP_T5224pre_IL1_6h - SP_Veh_6h,
  # --- 6 Hour Timepoint ---
  IL1_vs_Veh_6h = SP_IL1_6h - SP_Veh_6h,
  T5224_vs_Veh_6h = SP_T5224_6h - SP_Veh_6h,
  T5224preIL1_vs_IL1_6h = SP_T5224pre_IL1_6h - SP_IL1_6h,
  
  # --- Time/Treatment Interactions (Optional) ---
  # Is the effect of IL-1 significantly different at 6h compared to 3h?
  #Interaction_IL1_Time = (SP_IL1_6h - SP_Veh_6h) - (SP_IL1_3h - SP_Veh_3h),
    # Temporal Contrasts
  #Time_Baseline = SP_Veh_6h - SP_Veh_3h,
  #Time_IL1_Trajectory = SP_IL1_6h - SP_IL1_3h,
  #Time_T5224_Trajectory = SP_T5224_6h - SP_T5224_3h,
  
  # Interaction Contrasts
  #Interaction_Rescue_Time = (SP_T5224pre_IL1_6h - SP_IL1_6h) - (SP_T5224pre_IL1_3h - SP_IL1_3h),

  levels = design

)


norm_counts_log2 <- cpm(y, log = TRUE, prior.count = 1)

# 2. Convert to a data frame and fix the row names
norm_counts_df <- as.data.frame(norm_counts_log2) %>%
  rownames_to_column("gene_name")

# 3. Write to CSV or TSV
write_csv(norm_counts_df, "normalized_counts_log2_cpm.csv")

colnames(norm_counts_df)
getwd()

# Load required libraries
library(edgeR)
library(tidyverse)

# ----------------------------------------------------------------------------
# 1. Prepare Metadata for Plotting
# ----------------------------------------------------------------------------
# We extract the condition and timepoint directly from your group names
# to use for grouping and colouring in ggplot.
sample_info <- data.frame(
  sample_id = colnames(y),
  group = y$samples$group,
  lib_size = y$samples$lib.size
) %>%
  mutate(
    # Extract the base condition (e.g., SP_Veh, SP_T5224pre_IL1)
    condition = sub("_[36]h$", "", group),
    # Extract the timepoint
    time = str_extract(group, "[36]h"),
    # Count how many genes have at least 1 read in each sample
    detected_genes = colSums(y$counts > 0)
  ) %>%
  dplyr::select(sample_id, group, condition, time, lib_size, detected_genes)

# ----------------------------------------------------------------------------
# 2. Total Mapped Reads (Library Size)
# ----------------------------------------------------------------------------
# Shows if any samples failed sequencing or mapping.
p_reads <- ggplot(sample_info, aes(x = sample_id, y = lib_size / 1e6, fill = condition)) +
  geom_col() +
  facet_grid(~ time, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(
    title = "Total Mapped Reads per Sample", 
    y = "Millions of Reads", 
    x = ""
  )
p_reads
#ggsave("QC_Total_Reads.png", p_reads, width = 10, height = 6)

# ----------------------------------------------------------------------------
# 3. Number of Detected Genes
# ----------------------------------------------------------------------------
# Shows library complexity. A drop here indicates RNA degradation.
p_genes <- ggplot(sample_info, aes(x = sample_id, y = detected_genes, fill = condition)) +
  geom_col() +
  facet_grid(~ time, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(
    title = "Detected Genes (>0 counts)", 
    y = "Number of Genes", 
    x = ""
  )
p_genes
#ggsave("QC_Detected_Genes.png", p_genes, width = 10, height = 6)

# ----------------------------------------------------------------------------
# 4. Multi-Dimensional Scaling (MDS) Plot
# ----------------------------------------------------------------------------
# The most important plot. Shows how samples cluster based on leading logFC.
# We calculate coordinates using edgeR, then plot with ggplot2.
mds_calc <- plotMDS(y, plot = FALSE)
mds_data <- data.frame(
  Dim1 = mds_calc$x,
  Dim2 = mds_calc$y,
  sample_info
)

p_mds <- ggplot(mds_data, aes(x = Dim1, y = Dim2, colour = condition, shape = time)) +
  geom_point(size = 5, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "MDS Plot (Sample Clustering)", 
    x = "Leading logFC Dim 1", 
    y = "Leading logFC Dim 2"
  )
p_mds
#ggsave("QC_MDS_Plot.png", p_mds, width = 8, height = 6)

# ----------------------------------------------------------------------------
# 5. Normalised Expression Distributions (Boxplots)
# ----------------------------------------------------------------------------
# Ensures TMM normalisation successfully aligned the sample distributions.
logCPM <- cpm(y, log = TRUE)
logCPM_long <- as.data.frame(logCPM) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "logCPM") %>%
  left_join(sample_info, by = "sample_id")

p_dist <- ggplot(logCPM_long, aes(x = sample_id, y = logCPM, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  facet_grid(~ time, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(
    title = "Normalised Log2 CPM Distributions", 
    y = "Log2 CPM", 
    x = ""
  )
p_dist
#ggsave("QC_Normalised_Distributions.png", p_dist, width = 10, height = 6)

# ----------------------------------------------------------------------------
# 6. Biological Coefficient of Variation (BCV)
# ----------------------------------------------------------------------------
# This is best generated using the base R edgeR function. 
# It plots the biological variation against gene abundance.
#png("QC_BCV_Plot.png", width = 800, height = 600, res = 100)
plotBCV(y, main = "Biological Coefficient of Variation")
#dev.off()

# 3. Loop through each contrast, run the test, save the CSV, and extract the summary
# Initialize an empty list to store the decideTests results

library(edgeR)
library(tidyverse)

# 1. Create a dedicated directory to keep your outputs organised
output_dir <- "DGE_Results"
dir.create(output_dir, showWarnings = FALSE)
dt_list <- list()
# 2. Extract the names of all the contrasts defined in 'my_contrasts'
contrast_names <- colnames(my_contrasts)
# 3. Loop through each contrast, run the test, save the CSV, and extract the summary
# Initialize an empty list to store the summaries
dt_summary_list <- list()

for (c_name in contrast_names) {
  
  # Run the Quasi-Likelihood F-test for the specific contrast
  qlf <- glmQLFTest(fit, contrast = my_contrasts[, c_name])
  
  # Extract and save the full DGE table
  res_table <- topTags(qlf, n = Inf)$table %>%
    rownames_to_column("gene_name")
  file_path <- paste0(output_dir, "/DGE_", c_name, ".csv")
  write_csv(res_table, file_path)
  
  # Extract the significance calls for this test
  dt <- decideTests(qlf, p.value = 0.05, adjust.method = "BH")
  
  # Run the correct edgeR/limma summary immediately and store it
  dt_summary_list[[c_name]] <- summary(dt)
}

# 4. Generate a master summary of significant genes
# Bind the summarized columns together
summary_matrix <- do.call(cbind, dt_summary_list)
colnames(summary_matrix) <- contrast_names

# Convert it into a clean, readable dataframe
summary_table <- as.data.frame(summary_matrix) %>%
  rownames_to_column("Direction") %>%
  mutate(
    Direction = case_when(
      Direction == "-1" ~ "Downregulated",
      Direction == "0"  ~ "Not Significant",
      Direction == "1"  ~ "Upregulated"
    )
  )

# Save the summary table
write_csv(summary_table, paste0(output_dir, "/Master_Significance_Summary.csv"))

# Print the summary to the console for immediate review
print(summary_table)

# Load required libraries
library(clusterProfiler)
library(org.Mm.eg.db) # Mouse annotation database
library(tidyverse)
library(enrichplot)
setwd(publication_path("BULK_RNA_IL1_ROOT", "star_salmon"))
# 1. Define your background (universe)
all_genes <- rownames(y) 
universe_entrez <- bitr(all_genes, fromType = "SYMBOL", 
                        toType = "ENTREZID", 
                        OrgDb = org.Mm.eg.db)

# 2. Setup a directory for the split enrichment outputs
enrich_dir <- "Enrichment_Split_Results"
dir.create(enrich_dir, showWarnings = FALSE)

# 3. Get the list of your significant result files
result_files <- list.files("DGE_Results", pattern = "^DGE_.*\\.csv$", full.names = TRUE)

# 4. Loop through each contrast and run the split enrichment
for (file in result_files) {
  
  c_name <- str_remove(basename(file), "^DGE_") %>% str_remove("\\.csv$")
  if (file.exists(paste0(enrich_dir, "/Split_KEGG_plot_", c_name, ".png"))){next}
  res <- read_csv(file, show_col_types = FALSE)
  
  # Separate genes into Upregulated (logFC > 1) and Downregulated (logFC < -1)
  up_genes <- res %>% filter(FDR < 0.05 & logFC > 1) %>% pull(gene_name)
  down_genes <- res %>% filter(FDR < 0.05 & logFC < -1) %>% pull(gene_name)
  
  # Initialize an empty list to hold the ID clusters
  gene_clusters <- list()
  
  # Convert and add Upregulated genes if any exist
  if (length(up_genes) > 0) {
    up_entrez <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
    gene_clusters[["Upregulated"]] <- up_entrez$ENTREZID
  }
  
  # Convert and add Downregulated genes if any exist
  if (length(down_genes) > 0) {
    down_entrez <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
    gene_clusters[["Downregulated"]] <- down_entrez$ENTREZID
  }
  
  # If there are no significant genes in either direction, skip to the next file
  if (length(gene_clusters) == 0) next
  
  # ---------------------------------------------------------
  # GO Biological Process - Split Comparison
  # ---------------------------------------------------------
  comp_go <- compareCluster(geneCluster   = gene_clusters,
                            fun           = "enrichGO",
                            universe      = universe_entrez$ENTREZID,
                            OrgDb         = org.Mm.eg.db,
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.2,
                            readable      = TRUE)
  
  if (!is.null(comp_go) && nrow(comp_go) > 0) {
    write_csv(as.data.frame(comp_go), paste0(enrich_dir, "/Split_GO_BP_", c_name, ".csv"))
    
    # showCategory = 10 will pull the top 10 terms for Up and top 10 for Down
    p_go <- dotplot(comp_go, showCategory = 10, title = paste("GO BP (Split):", c_name))
    ggsave(paste0(enrich_dir, "/Split_GO_BP_plot_", c_name, ".png"), p_go, width = 11, height = 8)
  }
  
  # ---------------------------------------------------------
  # KEGG Pathway - Split Comparison
  # ---------------------------------------------------------
  comp_kegg <- compareCluster(geneCluster   = gene_clusters,
                              fun           = "enrichKEGG",
                              universe      = universe_entrez$ENTREZID,
                              organism      = "mmu", # Mouse organism code
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05)
  
  # Translate Entrez IDs back to Gene Symbols
  if (!is.null(comp_kegg)) {
      comp_kegg <- setReadable(comp_kegg, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
  }
  
  if (!is.null(comp_kegg) && nrow(comp_kegg) > 0) {
    write_csv(as.data.frame(comp_kegg), paste0(enrich_dir, "/Split_KEGG_", c_name, ".csv"))
    
    p_kegg <- dotplot(comp_kegg, showCategory = 10, title = paste("KEGG (Split):", c_name))
    ggsave(paste0(enrich_dir, "/Split_KEGG_plot_", c_name, ".png"), p_kegg, width = 11, height = 8)
  }
}

print("Split enrichment analysis for mm10 complete. Check the 'Enrichment_Split_Results' directory.")

# Load required libraries
library(clusterProfiler)
library(org.Mm.eg.db) # Mouse annotation database
library(tidyverse)
library(enrichplot)

# 1. Setup a directory for GSEA outputs
gsea_dir <- "GSEA_Results"
dir.create(gsea_dir, showWarnings = FALSE)

# 2. Get the list of your result files from the previous edgeR step
result_files <- list.files("DGE_Results", pattern = "^DGE_.*\\.csv$", full.names = TRUE)

# 3. Loop through each contrast and run GSEA
for (file in result_files) {
  
  c_name <- str_remove(basename(file), "^DGE_") %>% str_remove("\\.csv$")
    if (file.exists(paste0(gsea_dir, "/GSEA_KEGG_ridge_", c_name, ".png"))){next}
  res <- read_csv(file, show_col_types = FALSE)
  
  # ---------------------------------------------------------
  # Prepare the Ranked Gene List
  # ---------------------------------------------------------
  # GSEA requires a named numeric vector sorted in descending order.
  # We must also handle duplicate gene symbols by keeping the one with the highest absolute logFC.
  res_clean <- res %>%
    filter(!is.na(logFC) & !is.na(gene_name)) %>%
    arrange(desc(abs(logFC))) %>%
    distinct(gene_name, .keep_all = TRUE)
  
  # Convert symbols to Entrez IDs for KEGG compatibility
  entrez_map <- bitr(res_clean$gene_name, 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Mm.eg.db)
  
  # Merge the Entrez IDs back onto our results
  res_mapped <- res_clean %>%
    inner_join(entrez_map, by = c("gene_name" = "SYMBOL"))
  
  # Create the final ranked vector: Values are logFC, names are Entrez IDs
  gene_list <- res_mapped$logFC
  names(gene_list) <- res_mapped$ENTREZID
  
  # Sort strictly in descending order (critical for GSEA)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # ---------------------------------------------------------
  # Gene Ontology (GO) - Biological Process GSEA
  # ---------------------------------------------------------
  gse_go <- gseGO(geneList      = gene_list,
                  OrgDb         = org.Mm.eg.db,
                  ont           = "BP",
                  minGSSize     = 15,
                  maxGSSize     = 500,
                  pvalueCutoff  = 0.2,
                  pAdjustMethod = "BH",
                  verbose       = FALSE)
  
  if (!is.null(gse_go) && nrow(gse_go) > 0) {
    # Translate Entrez IDs back to Gene Symbols in the core enrichment column
    gse_go <- setReadable(gse_go, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
    write_csv(as.data.frame(gse_go), paste0(gsea_dir, "/GSEA_GO_BP_", c_name, ".csv"))
    
    # Ridge plots are excellent for visualising the distribution of core enriched genes
    p_go_ridge <- ridgeplot(gse_go, showCategory = 15) + 
      ggtitle(paste("GSEA GO BP:", c_name))
    ggsave(paste0(gsea_dir, "/GSEA_GO_BP_ridge_", c_name, ".png"), p_go_ridge, width = 12, height = 8)
  }
  
  # ---------------------------------------------------------
  # KEGG Pathway GSEA
  # ---------------------------------------------------------
  gse_kegg <- gseKEGG(geneList      = gene_list,
                      organism      = "mmu",
                      minGSSize     = 15,
                      maxGSSize     = 500,
                      pvalueCutoff  = 0.2,
                      pAdjustMethod = "BH",
                      verbose       = FALSE)
  
  if (!is.null(gse_kegg) && nrow(gse_kegg) > 0) {
    gse_kegg <- setReadable(gse_kegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
    write_csv(as.data.frame(gse_kegg), paste0(gsea_dir, "/GSEA_KEGG_", c_name, ".csv"))
    
    p_kegg_ridge <- ridgeplot(gse_kegg, showCategory = 15) + 
      ggtitle(paste("GSEA KEGG:", c_name))
    ggsave(paste0(gsea_dir, "/GSEA_KEGG_ridge_", c_name, ".png"), p_kegg_ridge, width = 12, height = 8)
  }
}

print("GSEA complete. Check the 'GSEA_Results' directory.")



result_files

# Load required libraries
library(tidyverse)

# 1. Define the AP-1 complex members
# Note: Mouse gene symbols use title case (e.g., Jun, not JUN)
ap1_genes <- c("Jun", "Junb", "Jund", "Fos", "Fosb", "Fosl1", "Fosl2", 
               "Atf2", "Atf3", "Atf4", "Atf7", "Batf", "Batf2", "Batf3", "Jdp2")

# 2. Locate and load the DGE results
result_files <- list.files("DGE_Results", pattern = "^DGE_.*\\.csv$", full.names = TRUE)
result_files = result_files[!grepl(result_files, pattern = "Trajectory|Interaction|Time")]
# Read all files into a single dataframe and label them by contrast
ap1_data <- map_dfr(result_files, function(file) {
  contrast_name <- str_remove(basename(file), "^DGE_") %>% str_remove("\\.csv$")
  
  read_csv(file, show_col_types = FALSE) %>%
    mutate(Contrast = contrast_name)
})

# 3. Filter for AP-1 genes and format the data
plot_data <- ap1_data %>%
  filter(gene_name %in% ap1_genes) %>%
  mutate(
    # Add an asterisk if the gene is significantly differentially expressed
    Significance = if_else(FDR < 0.05, "*", ""),
    
    # Cap extreme logFC values for better color scaling (optional but recommended)
    # This prevents one massive outlier from washing out the color gradient
    logFC_capped = case_when(
      logFC > 3 ~ 3,
      logFC < -3 ~ -3,
      TRUE ~ logFC
    )
  ) %>%
  # Filter out AP-1 genes that were completely unexpressed across all samples
  filter(!is.na(logFC)) 

# 4. Generate the Bubble Plot
p_ap1 <- ggplot(plot_data, aes(x = Contrast, y = gene_name)) +
  # Map bubble size to Expression (logCPM) and color to Fold Change (logFC)
  geom_point(aes(size = logCPM, color = logFC_capped)) +
  # Add the significance asterisks directly on top of the bubbles
  geom_text(aes(label = Significance), vjust = 0.75, size = 6, color = "black") +
  # Use a diverging color palette: Blue for down, White for neutral, Red for up
  scale_color_gradient2(
    low = "dodgerblue3", mid = "white", high = "firebrick3", midpoint = 0,
    name = "Log2 Fold Change",
    limits = c(-3, 3), oob = scales::squish # squish keeps capped values at the extreme colors
  ) +
  scale_size_continuous(range = c(2, 10), name = "Expression (logCPM)") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 10),
    axis.text.y = element_text(face = "italic", size = 11), # Gene names should be italicized
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  ) +
  labs(
    title = "Expression and Regulation of AP-1 Complex Elements",
    subtitle = "Size = Expression Level; Color = Fold Change; * = FDR < 0.05",
    x = "",
    y = ""
  )

# 5. Save and view the plot
ggsave("AP1_Elements_Expression_Plot.png", p_ap1, width = 12, height = 8, dpi = 300)
print(p_ap1)

# Load required libraries
library(tidyverse)

# 1. Define the AP-1 complex members
# Note: Mouse gene symbols use title case (e.g., Jun, not JUN)
IL1_targets_mouse <- c(
  # Cytokines / inflammatory mediators
  "Il6", "Il1a", "Il1b", "Il11", "Lif", "Tnf", "Csf2", "Csf3",

  # Chemokines
  "Cxcl1", "Cxcl2", "Cxcl3", "Cxcl5", "Cxcl10",
  "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Ccl7", "Ccl20",

  # Prostaglandin / lipid inflammatory response
  "Ptgs2", "Ptges", "Ptges2", "Pla2g4a",

  # NF-kB / IL-1 feedback and signalling response
  "Nfkbia", "Nfkbiz", "Nfkb2", "Relb",
  "Tnfaip3", "Tnfaip2", 

  # MAPK / AP-1 immediate early response
  "Egr1", "Egr2",
  "Dusp1", "Dusp2", "Dusp5",

  # Adhesion / inflammatory activation
  "Icam1", "Vcam1",

  # Matrix remodelling / fibroblast inflammatory output
  "Mmp3", "Mmp9", "Mmp13", "Mmp14",
  "Timp1", "Has1", "Has2",

  # Inflammasome / innate immune response
  "Nos2", "Nlrp3", "Saa3"
)

IL1_targets_mouse <- unique(IL1_targets_mouse)
# 2. Locate and load the DGE results
result_files <- list.files("DGE_Results", pattern = "^DGE_.*\\.csv$", full.names = TRUE)
result_files = result_files[!grepl(result_files, pattern = "Trajectory|Interaction|Time")]
# Read all files into a single dataframe and label them by contrast
ap1_data <- map_dfr(result_files, function(file) {
  contrast_name <- str_remove(basename(file), "^DGE_") %>% str_remove("\\.csv$")
  
  read_csv(file, show_col_types = FALSE) %>%
    mutate(Contrast = contrast_name)
})

# 3. Filter for AP-1 genes and format the data
plot_data <- ap1_data %>%
  filter(gene_name %in% IL1_targets_mouse) %>%
  mutate(
    # Add an asterisk if the gene is significantly differentially expressed
    Significance = if_else(FDR < 0.05, "*", ""),
    
    # Cap extreme logFC values for better color scaling (optional but recommended)
    # This prevents one massive outlier from washing out the color gradient
    logFC_capped = case_when(
      logFC > 3 ~ 3,
      logFC < -3 ~ -3,
      TRUE ~ logFC
    )
  ) %>%
  # Filter out AP-1 genes that were completely unexpressed across all samples
  filter(!is.na(logFC)) 

# 4. Generate the Bubble Plot
p_ap1 <- ggplot(plot_data, aes(x = Contrast, y = gene_name)) +
  # Map bubble size to Expression (logCPM) and color to Fold Change (logFC)
  geom_point(aes(size = logCPM, color = logFC_capped)) +
  # Add the significance asterisks directly on top of the bubbles
  geom_text(aes(label = Significance), vjust = 0.75, size = 6, color = "black") +
  # Use a diverging color palette: Blue for down, White for neutral, Red for up
  scale_color_gradient2(
    low = "dodgerblue3", mid = "white", high = "firebrick3", midpoint = 0,
    name = "Log2 Fold Change",
    limits = c(-3, 3), oob = scales::squish # squish keeps capped values at the extreme colors
  ) +
  scale_size_continuous(range = c(2, 10), name = "Expression (logCPM)") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 10),
    axis.text.y = element_text(face = "italic", size = 11), # Gene names should be italicized
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  ) +
  labs(
    title = "Expression and Regulation of IL-1 Complex Elements",
    subtitle = "Size = Expression Level; Color = Fold Change; * = FDR < 0.05",
    x = "",
    y = ""
  )

# 5. Save and view the plot
ggsave("IL-1_Elements_Expression_Plot.pdf", p_ap1, width = 12, height = 10)
print(p_ap1)

getwd()

# Load required libraries
library(edgeR)
library(tidyverse)

# 1. Define the AP-1 complex members (Mouse symbols)
ap1_genes <- c("Jun", "Junb", "Jund", "Fos", "Fosb", "Fosl1", "Fosl2", 
               "Atf2", "Atf3", "Atf4", "Atf7", "Batf", "Batf2", "Batf3", "Jdp2")

# 2. Extract and format the logCPM data from your DGEList ('y')
# Ensure 'y' is the filtered and normalised object from the earlier steps
logCPM <- cpm(y, log = TRUE)

# Filter for the AP-1 genes that are actually expressed in your dataset
ap1_cpm <- logCPM[rownames(logCPM) %in% ap1_genes, ]

# 3. Calculate Z-scores (row-wise scaling)
# This standardises each gene so its mean is 0 and standard deviation is 1
ap1_scaled <- t(scale(t(ap1_cpm)))

# 4. Prepare data for ggplot
heatmap_data <- as.data.frame(ap1_scaled) %>%
  rownames_to_column("gene_name") %>%
  pivot_longer(-gene_name, names_to = "sample_id", values_to = "z_score") %>%
  # Extract the group from the sample name to control the plotting order
  mutate(
    group = sub("_R[1-3]$", "", sample_id),
    time = str_extract(group, "[36]h"),
    condition = sub("_[36]h$", "", group)
  )

# 5. Order the factors for a logical layout on the plot
# We want Vehicle -> IL1 -> T5224 -> T5224+IL1
condition_order <- c("SP_Veh", "SP_IL1", "SP_T5224", "SP_T5224pre_IL1")
heatmap_data$condition <- factor(heatmap_data$condition, levels = condition_order)

# Create a combined factor for the x-axis to group replicates together neatly
heatmap_data <- heatmap_data %>%
  arrange(time, condition, sample_id) %>%
  mutate(sample_id = factor(sample_id, levels = unique(sample_id)))

# 6. Generate the Heatmap
p_heatmap <- ggplot(heatmap_data, aes(x = sample_id, y = gene_name, fill = z_score)) +
  geom_tile(color = "white", size = 0.5) + # Add white borders for crispness
  # Use a diverging Red/Blue palette centered on 0
  scale_fill_gradient2(
    low = "dodgerblue3", mid = "white", high = "firebrick3", midpoint = 0,
    name = "Relative Expression\n(Z-score)",
    limits = c(-3, 3), oob = scales::squish
  ) +
  # Split the plot by timepoint
  facet_grid(~ time, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9, face = "bold"),
    axis.text.y = element_text(face = "italic", size = 11),
    strip.text = element_text(size = 14, face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 16)
  ) +
  labs(
    title = "AP-1 Complex Expression Across All Samples",
    subtitle = "Data row-scaled (Z-score). Red = Higher than average, Blue = Lower than average.",
    x = "",
    y = ""
  )

# 7. Save and display
ggsave("AP1_Sample_Expression_Heatmap.png", p_heatmap, width = 14, height = 8, dpi = 300)
print(p_heatmap)

# Load required libraries
library(tidyverse)

# ============================================================================
# PLOT 1: Unified GSEA Dotplot (Directional Enrichment)
# ============================================================================

# 1. Locate and load all GSEA GO BP results
gsea_files <- list.files("GSEA_Results", pattern = "^GSEA_GO_BP_.*\\.csv$", full.names = TRUE)

gsea_data <- map_dfr(gsea_files, function(file) {
  # Extract contrast name
  c_name <- str_remove(basename(file), "^GSEA_GO_BP_") %>% str_remove("\\.csv$")
  read_csv(file, show_col_types = FALSE) %>% mutate(Contrast = c_name)
})

# 2. Extract the top 5 most significant terms per contrast
top_gsea_terms <- gsea_data %>%
  group_by(Contrast) %>%
  slice_min(order_by = p.adjust, n =7, with_ties = FALSE) %>%
  pull(Description) %>%
  unique()

# 3. Format data for plotting
plot_gsea <- gsea_data %>%
  filter(Description %in% top_gsea_terms) %>%
  mutate(
    # Order y-axis by overall best significance
    Description = fct_reorder(Description, p.adjust, .fun = min, .desc = TRUE),
    # Convert FDR to -log10 scale so smaller p-values make bigger dots
    log_FDR = -log10(p.adjust)
  )

# 4. Generate GSEA Plot
p_gsea <- ggplot(plot_gsea, aes(x = Contrast, y = Description)) +
  geom_point(aes(size = log_FDR, color = NES)) +
  # Diverging palette: Blue for negative NES (Down), Red for positive NES (Up)
  scale_color_gradient2(
    low = "dodgerblue3", mid = "white", high = "firebrick3", midpoint = 0, 
    name = "NES (Enrichment)"
  ) +
  scale_size_continuous(range = c(3, 8), name = "-log10(FDR)") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 14)
  ) +
  labs(
    title = "Master GSEA: Top GO Biological Processes",
    subtitle = "Colour indicates direction of pathway shift (Normalised Enrichment Score)",
    x = "", y = ""
  )

ggsave("Master_GSEA_GO_Plot.png", p_gsea, width = 13, height = 9, dpi = 300)


# ============================================================================
# PLOT 2: Unified Split GO Dotplot (Upregulated vs Downregulated)
# ============================================================================

# 1. Locate and load all Split GO BP results
split_files <- list.files("Enrichment_Split_Results", pattern = "^Split_GO_BP_.*\\.csv$", full.names = TRUE)

split_data <- map_dfr(split_files, function(file) {
  c_name <- str_remove(basename(file), "^Split_GO_BP_") %>% str_remove("\\.csv$")
  read_csv(file, show_col_types = FALSE) %>% mutate(Contrast = c_name)
})

# 2. Parse the GeneRatio column from character ("15/100") to numeric (0.15)
# 2 & 3. Parse the GeneRatio column safely
# We use sub() to extract the numbers before and after the slash
split_data <- split_data %>%
  mutate(
    # sub("/.*", "") removes the slash and everything after it
    num = as.numeric(sub("/.*", "", GeneRatio)),
    # sub(".*/", "") removes the slash and everything before it
    den = as.numeric(sub(".*/", "", GeneRatio)),
    # Calculate the decimal
    GeneRatio_Num = num / den
  ) %>%
  dplyr::select(-num, -den)
# 3. Extract the top 4 most significant terms per contrast AND per direction
# We use n = 4 here because splitting by direction doubles the amount of terms retrieved
top_split_terms <- split_data %>%
  group_by(Contrast, Cluster) %>%
  slice_min(order_by = p.adjust, n = 4, with_ties = FALSE) %>%
  pull(Description) %>%
  unique()

# 4. Format data for plotting
plot_split <- split_data %>%
  filter(Description %in% top_split_terms) %>%
  mutate(Description = fct_reorder(Description, p.adjust, .fun = min, .desc = TRUE))

# 5. Generate Split Plot
p_split <- ggplot(plot_split, aes(x = Contrast, y = Description)) +
  geom_point(aes(size = GeneRatio_Num, color = p.adjust)) +
  scale_color_viridis_c(guide = guide_colorbar(reverse = TRUE), name = "FDR (p.adjust)") +
  scale_size_continuous(range = c(2, 8), name = "Gene Ratio") +
  # Split the plot into two distinct panels based on direction
  facet_wrap(~ Cluster) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 11),
    strip.text = element_text(face = "bold", size = 13), # Makes the Up/Down labels pop
    strip.background = element_rect(fill = "grey90"),
    plot.title = element_text(face = "bold", size = 14)
  ) +
  labs(
    title = "Master Split GO: Dysregulated Biological Processes",
    subtitle = "Pathways are separated by their driving fold-change direction",
    x = "", y = ""
  )

ggsave("Master_Split_GO_Plot.pdf", p_split, width = 16, height = 10, dpi = 300)

# Print to console/viewer
print(p_gsea)
print(p_split)

# Load required libraries
library(tidyverse)
library(ggVennDiagram)
library(UpSetR)

# 1. Locate the DGE results
result_files <- list.files("DGE_Results", pattern = "^DGE_.*\\.csv$", full.names = TRUE)

# 2. Extract lists of significant genes for each contrast
# We will create a named list where each element is a vector of significant gene names (FDR < 0.05)
sig_gene_lists <- list()

for (file in result_files) {
  # Clean up the contrast name for the plot labels
  c_name <- str_remove(basename(file), "^DGE_") %>% str_remove("\\.csv$")
  
  # Read data and pull only the significant genes
  res <- read_csv(file, show_col_types = FALSE)
  sig_genes <- res %>%
    filter(FDR < 0.05) %>%
    pull(gene_name) # Extracts as a vector
  
  # Store in our list
  sig_gene_lists[[c_name]] <- sig_genes
}

# ============================================================================
# PLOT 1: Targeted Venn Diagram (e.g., The 3-Hour Rescue Effect)
# ============================================================================
# Select only the specific contrasts you want in the Venn (max 4 recommended)
venn_contrasts <- c("IL1_vs_Veh_3h", "T5224_vs_Veh_3h", "T5224preIL1_vs_IL1_3h")

# Subset our master list
venn_list <- sig_gene_lists[venn_contrasts]

# Generate the Venn Diagram
p_venn <- ggVennDiagram(venn_list, 
                        label_alpha = 0, # Removes the background box behind numbers
                        set_size = 4.5) +
  scale_fill_gradient(low = "white", high = "firebrick3", name = "Gene Count") +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5)) +
  labs(title = "Overlap of Significant Genes at 3 Hours")

ggsave("Venn_Diagram_3h_Rescue.png", p_venn, width = 8, height = 8, bg = "white")
print(p_venn)

venn_contrasts <- c("IL1_vs_Veh_6h", "T5224_vs_Veh_6h", "T5224preIL1_vs_IL1_6h")

# Subset our master list
venn_list <- sig_gene_lists[venn_contrasts]

# Generate the Venn Diagram
p_venn <- ggVennDiagram(venn_list, 
                        label_alpha = 0, # Removes the background box behind numbers
                        set_size = 4.5) +
  scale_fill_gradient(low = "white", high = "firebrick3", name = "Gene Count") +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5)) +
  labs(title = "Overlap of Significant Genes at 6 Hours")

ggsave("Venn_Diagram_6h_Rescue.png", p_venn, width = 8, height = 8, bg = "white")
print(p_venn)

venn_contrasts <- c("IL1_vs_Veh_6h", "T5224_vs_Veh_6h","IL1_vs_Veh_3h", "T5224_vs_Veh_3h")

# Subset our master list
venn_list <- sig_gene_lists[venn_contrasts]

# Generate the Venn Diagram
p_venn <- ggVennDiagram(venn_list, 
                        label_alpha = 0, # Removes the background box behind numbers
                        set_size = 4.5) +
  scale_fill_gradient(low = "white", high = "firebrick3", name = "Gene Count") +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5)) +
  labs(title = "Overlap of Significant Genes")

ggsave("Venn_Diagram_Basic.png", p_venn, width = 8, height = 8, bg = "white")
print(p_venn)
# ============================================================================
# PLOT 2: Global UpSet Plot (All Contrasts)
# ============================================================================
# UpSet plots are brilliant for showing intersections across many sets. 
# UpSetR requires a specific matrix format, but 'fromList' converts our named list automatically.

# Save the plot directly to a PNG (UpSetR uses base R plotting, not ggplot2)
pdf("UpSet_Plot_All_Contrasts.pdf", width = 35, height = 10,)

upset(fromList(sig_gene_lists), 
      nsets = length(sig_gene_lists), # Include all your contrasts
      nintersects = 80,               # Show the top 20 largest intersections
      order.by = "freq",              # Order the bars by the size of the intersection
      keep.order = TRUE,
      mainbar.y.label = "Intersection Size (Number of Genes)", 
      sets.x.label = "Total Significant Genes per Contrast", 
      text.scale = c(1.5, 1.5, 1.2, 1.2, 1.5, 1.3), # Scales various text elements
      point.size = 3.5, 
      line.size = 1,
      shade.color = "grey80",
      matrix.color = "dodgerblue4",
      main.bar.color = "grey30")

dev.off()

# Load required libraries
library(tidyverse)

# ============================================================================
# PART 1: Find and Plot ECM-Related GSEA Terms
# ============================================================================

# 1. Load all GSEA results
gsea_files <- list.files("GSEA_Results", pattern = "^GSEA_GO_BP_.*\\.csv$", full.names = TRUE)
gsea_files = gsea_files[!grepl(gsea_files, pattern = "Trajectory|Interaction|Time")]
gsea_data <- map_dfr(gsea_files, function(file) {
  c_name <- str_remove(basename(file), "^GSEA_GO_BP_") %>% str_remove("\\.csv$")
  read_csv(file, show_col_types = FALSE) %>% mutate(Contrast = c_name)
})

# 2. Filter for ECM and Collagen related terms using regular expressions
# ignore_case = TRUE ensures we catch "Extracellular", "extracellular", etc.
ecm_gsea <- gsea_data %>%
  filter(str_detect(Description, regex("extracellular matrix|collagen", ignore_case = TRUE)))

# 3. Format data for the GSEA Dotplot
plot_ecm_terms <- ecm_gsea %>%
  mutate(
    # Order the y-axis by the average NES so upregulated and downregulated terms group together
    Description = fct_reorder(Description, NES, .fun = mean),
    log_FDR = -log10(p.adjust)
  )

# 4. Generate the ECM GSEA Plot
p_ecm_terms <- ggplot(plot_ecm_terms, aes(x = Contrast, y = Description)) +
  geom_point(aes(size = log_FDR, color = NES)) +
  scale_color_gradient2(
    low = "dodgerblue3", mid = "white", high = "firebrick3", midpoint = 0, 
    name = "NES (Enrichment)"
  ) +
  scale_size_continuous(range = c(3, 8), name = "-log10(FDR)") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 14)
  ) +
  labs(
    title = "GSEA: Extracellular Matrix (ECM) Dynamics",
    subtitle = "Pathway-level shifts across all contrasts",
    x = "", y = ""
  )

ggsave("ECM_GSEA_Terms_Plot.png", p_ecm_terms, width = 12, height = 7, dpi = 300)


# ============================================================================
# PART 2: Extract and Plot the Contributing (Leading Edge) Genes
# ============================================================================

# 1. Extract the genes driving the significant ECM pathways
# GSEA stores these in the 'core_enrichment' column separated by slashes (e.g., "Col1a1/Mmp9/Fn1")
ecm_drivers <- ecm_gsea %>%
  filter(p.adjust < 0.05) %>%
  # separate_rows splits the slash-separated string into individual rows for counting
  separate_rows(core_enrichment, sep = "/") %>%
  # Count how many times each gene appears in a significant ECM pathway
  count(core_enrichment, sort = TRUE) 

# Take the top 30 most frequent ECM drivers to keep the plot readable
top_ecm_genes <- ecm_drivers %>%
  slice_head(n = 30) %>%
  pull(core_enrichment)

# 2. Locate and load the DGE results to get the fold-changes for these specific genes
dge_files <- list.files("DGE_Results", pattern = "^DGE_.*\\.csv$", full.names = TRUE)

ecm_dge_data <- map_dfr(dge_files, function(file) {
  c_name <- str_remove(basename(file), "^DGE_") %>% str_remove("\\.csv$")
  read_csv(file, show_col_types = FALSE) %>% mutate(Contrast = c_name)
})

# 3. Filter the DGE master list for our top ECM genes
plot_ecm_genes <- ecm_dge_data %>%
  filter(gene_name %in% top_ecm_genes) %>%
  mutate(
    Significance = if_else(FDR < 0.05, "*", ""),
    logFC_capped = case_when(
      logFC > 3 ~ 3,
      logFC < -3 ~ -3,
      TRUE ~ logFC
    )
  ) %>%
  filter(!is.na(logFC))

# 4. Generate the Contributing Genes Bubble Plot
p_ecm_genes <- ggplot(plot_ecm_genes, aes(x = Contrast, y = gene_name)) +
  geom_point(aes(size = logCPM, color = logFC_capped)) +
  geom_text(aes(label = Significance), vjust = 0.75, size = 5, color = "black") +
  scale_color_gradient2(
    low = "dodgerblue3", mid = "white", high = "firebrick3", midpoint = 0,
    name = "Log2 Fold Change",
    limits = c(-3, 3), oob = scales::squish
  ) +
  scale_size_continuous(range = c(2, 8), name = "Expression (logCPM)") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 10),
    axis.text.y = element_text(face = "italic", size = 10),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  ) +
  labs(
    title = "Expression & Regulation of Key ECM Contributing Genes",
    subtitle = "Top 30 leading-edge genes driving the ECM GSEA pathways",
    x = "", y = ""
  )

ggsave("ECM_Contributing_Genes_Plot.png", p_ecm_genes, width = 12, height = 10, dpi = 300)

# Display in console
print(p_ecm_terms)
print(p_ecm_genes)

# Load required libraries
library(tidyverse)

# ============================================================================
# PART 1: Find and Plot ECM-Related GSEA Terms
# ============================================================================

# 1. Load all GSEA results
gsea_files <- list.files("GSEA_Results", pattern = "^GSEA_GO_BP_.*\\.csv$", full.names = TRUE)
gsea_files = gsea_files[!grepl(gsea_files, pattern = "Trajectory|Interaction|Time")]
gsea_data <- map_dfr(gsea_files, function(file) {
  c_name <- str_remove(basename(file), "^GSEA_GO_BP_") %>% str_remove("\\.csv$")
  read_csv(file, show_col_types = FALSE) %>% mutate(Contrast = c_name)
})

# 2. Filter for ECM and Collagen related terms using regular expressions
# ignore_case = TRUE ensures we catch "Extracellular", "extracellular", etc.
ecm_gsea <- gsea_data %>%
  filter(str_detect(Description, regex("Immune response", ignore_case = TRUE)))

# 3. Format data for the GSEA Dotplot
plot_ecm_terms <- ecm_gsea %>%
  mutate(
    # Order the y-axis by the average NES so upregulated and downregulated terms group together
    Description = fct_reorder(Description, NES, .fun = mean),
    log_FDR = -log10(p.adjust)
  )
p_ecm_terms = plot_ecm_terms[nchar(as.character(plot_ecm_terms$Description)) < 30,]
# 4. Generate the ECM GSEA Plot
p_ecm_terms <- ggplot(plot_ecm_terms, aes(x = Contrast, y = Description)) +
  geom_point(aes(size = log_FDR, color = NES)) +
  scale_color_gradient2(
    low = "dodgerblue3", mid = "white", high = "firebrick3", midpoint = 0, 
    name = "NES (Enrichment)"
  ) +
  scale_size_continuous(range = c(3, 8), name = "-log10(FDR)") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 14)
  ) +
  labs(
    title = "GSEA: Extracellular Matrix (ECM) Dynamics",
    subtitle = "Pathway-level shifts across all contrasts",
    x = "", y = ""
  )

ggsave("Immune_GSEA_Terms_Plot.png", p_ecm_terms, width = 12, height = 7, dpi = 300)


# ============================================================================
# PART 2: Extract and Plot the Contributing (Leading Edge) Genes
# ============================================================================

# 1. Extract the genes driving the significant ECM pathways
# GSEA stores these in the 'core_enrichment' column separated by slashes (e.g., "Col1a1/Mmp9/Fn1")
ecm_drivers <- ecm_gsea %>%
  filter(p.adjust < 0.05) %>%
  # separate_rows splits the slash-separated string into individual rows for counting
  separate_rows(core_enrichment, sep = "/") %>%
  # Count how many times each gene appears in a significant ECM pathway
  count(core_enrichment, sort = TRUE) 

# Take the top 30 most frequent ECM drivers to keep the plot readable
top_ecm_genes <- ecm_drivers %>%
  slice_head(n = 30) %>%
  pull(core_enrichment)

# 2. Locate and load the DGE results to get the fold-changes for these specific genes
dge_files <- list.files("DGE_Results", pattern = "^DGE_.*\\.csv$", full.names = TRUE)

ecm_dge_data <- map_dfr(dge_files, function(file) {
  c_name <- str_remove(basename(file), "^DGE_") %>% str_remove("\\.csv$")
  read_csv(file, show_col_types = FALSE) %>% mutate(Contrast = c_name)
})

# 3. Filter the DGE master list for our top ECM genes
plot_ecm_genes <- ecm_dge_data %>%
  filter(gene_name %in% top_ecm_genes) %>%
  mutate(
    Significance = if_else(FDR < 0.05, "*", ""),
    logFC_capped = case_when(
      logFC > 3 ~ 3,
      logFC < -3 ~ -3,
      TRUE ~ logFC
    )
  ) %>%
  filter(!is.na(logFC))

# 4. Generate the Contributing Genes Bubble Plot
p_ecm_genes <- ggplot(plot_ecm_genes, aes(x = Contrast, y = gene_name)) +
  geom_point(aes(size = logCPM, color = logFC_capped)) +
  geom_text(aes(label = Significance), vjust = 0.75, size = 5, color = "black") +
  scale_color_gradient2(
    low = "dodgerblue3", mid = "white", high = "firebrick3", midpoint = 0,
    name = "Log2 Fold Change",
    limits = c(-3, 3), oob = scales::squish
  ) +
  scale_size_continuous(range = c(2, 8), name = "Expression (logCPM)") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 10),
    axis.text.y = element_text(face = "italic", size = 10),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  ) +
  labs(
    title = "Expression & Regulation of Key ECM Contributing Genes",
    subtitle = "Top 30 leading-edge genes driving the ECM GSEA pathways",
    x = "", y = ""
  )

ggsave("ECM_Contributing_Genes_Plot.png", p_ecm_genes, width = 12, height = 10, dpi = 300)

# Display in console
print(p_ecm_terms)
print(p_ecm_genes)



getwd()

gsea_files

DT::datatable(gsea_data %>% filter(ID == "GO:0030198") , filter = 'top')

# Load required libraries
library(tidyverse)

# 1. Load the results for the primary stimulus and the rescue at 3h
res_il1 <- read_csv("DGE_Results/DGE_IL1_vs_Veh_3h.csv", show_col_types = FALSE)
res_rescue <- read_csv("DGE_Results/DGE_T5224preIL1_vs_IL1_3h.csv", show_col_types = FALSE)

# 2. Merge and rename columns to keep track of the specific contrasts
# We use an inner_join to ensure we only look at genes present in both tables
overlap_df <- res_il1 %>%
  dplyr::select(gene_name, logFC_IL1 = logFC, FDR_IL1 = FDR) %>%
  inner_join(
    res_rescue %>% dplyr::select(gene_name, logFC_Rescue = logFC, FDR_Rescue = FDR),
    by = "gene_name"
  )

# 3. Categorise the genes based on their directional shifts
# We are using an absolute logFC > 1 and FDR < 0.05 to ensure we are tracking strong biological drivers
classified_genes <- overlap_df %>%
  mutate(
    Status = case_when(
      
      # --- AP-1 Dependent (Successfully Rescued) ---
      # Upregulated by IL-1, but significantly knocked back down by T5224
       FDR_IL1 < 0.05 & logFC_IL1 > 0.5 & FDR_Rescue < 0.05 & logFC_Rescue < -0.5 ~ "Dependent_Up",
      FDR_IL1 < 0.05 & logFC_IL1 < -0.5 & FDR_Rescue < 0.05 & logFC_Rescue > 0.5 ~ "Dependent_Down",
      FDR_IL1 < 0.05 & logFC_IL1 > 0.5 & (FDR_Rescue >= 0.05 | logFC_Rescue >= 0) ~ "Escaping_Up",
      FDR_IL1 < 0.05 & logFC_IL1 < -0.5 & (FDR_Rescue >= 0.05 | logFC_Rescue <= 0) ~ "Escaping_Down",
      
      # Genes that didn't respond strongly to the primary IL-1 stimulus
      TRUE ~ "Not_Primary_Driver"
    )
  )

# 4. View the counts for each category in the console
print(table(classified_genes$Status))

# 5. Extract and save the biologically relevant lists
target_genes <- classified_genes %>% filter(Status != "Not_Primary_Driver")
write_csv(target_genes, "AP1_Target_Categorisation_3h.csv")

# 6. Generate a Quadrant Scatter Plot
# This maps the IL-1 fold-change on the X-axis against the Rescue fold-change on the Y-axis.
p_quadrant <- ggplot(target_genes, aes(x = logFC_IL1, y = logFC_Rescue, colour = Status)) +
  geom_point(alpha = 0.7, size = 2) +
  # Add reference lines for no change
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
  # Customise colours to make the dependent vs escaping clear
  scale_colour_manual(values = c(
    "Dependent_Up" = "firebrick3",
    "Dependent_Down" = "dodgerblue3",
    "Escaping_Up" = "lightpink",
    "Escaping_Down" = "lightblue"
  )) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "IL-1 Response vs AP-1 Rescue (3 Hours)",
    subtitle = "Dependent genes reverse direction; Escaping genes resist the inhibitor.",
    x = "Log2 Fold Change (IL-1 vs Vehicle)",
    y = "Log2 Fold Change (T5224+IL1 vs IL-1)",
    colour = "Gene Classification"
  )

# Save the plot
ggsave("AP1_Rescue_Quadrant_Plot_3h.png", p_quadrant, width = 10, height = 8, dpi = 300)

# Display the plot
print(p_quadrant)

# Load required libraries
library(tidyverse)

# ============================================================================
# PART 1: The 6-Hour Analysis
# ============================================================================

# 1. Load the 6h results
res_il1_6h <- read_csv("DGE_Results/DGE_IL1_vs_Veh_6h.csv", show_col_types = FALSE)
res_rescue_6h <- read_csv("DGE_Results/DGE_T5224preIL1_vs_IL1_6h.csv", show_col_types = FALSE)

# 2. Merge and categorise the 6h genes
classified_6h <- res_il1_6h %>%
  dplyr::select(gene_name, logFC_IL1 = logFC, FDR_IL1 = FDR) %>%
  inner_join(
    res_rescue_6h %>% dplyr::select(gene_name, logFC_Rescue = logFC, FDR_Rescue = FDR),
    by = "gene_name"
  ) %>%
  mutate(
    Timepoint = "6 Hours",
    Status = case_when(
      FDR_IL1 < 0.05 & logFC_IL1 > 0.5 & FDR_Rescue < 0.05 & logFC_Rescue < -0.5 ~ "Dependent_Up",
      FDR_IL1 < 0.05 & logFC_IL1 < -0.5 & FDR_Rescue < 0.05 & logFC_Rescue > 0.5 ~ "Dependent_Down",
      FDR_IL1 < 0.05 & logFC_IL1 > 0.5 & (FDR_Rescue >= 0.05 | logFC_Rescue >= 0) ~ "Escaping_Up",
      FDR_IL1 < 0.05 & logFC_IL1 < -0.5 & (FDR_Rescue >= 0.05 | logFC_Rescue <= 0) ~ "Escaping_Down",
      TRUE ~ "Not_Primary_Driver"
    )
  )

# Extract biologically relevant 6h targets
target_genes_6h <- classified_6h %>% filter(Status != "Not_Primary_Driver")
write_csv(target_genes_6h, "AP1_Target_Categorisation_6h.csv")

# Print 6h summary to console
cat("\n--- 6 Hour Categorisation Summary ---\n")
print(table(target_genes_6h$Status))

# ============================================================================
# PART 2: The Combined 3-Hour and 6-Hour Analysis
# ============================================================================

# 1. Recreate the 3h categorisation (assuming the environment was cleared)
res_il1_3h <- read_csv("DGE_Results/DGE_IL1_vs_Veh_3h.csv", show_col_types = FALSE)
res_rescue_3h <- read_csv("DGE_Results/DGE_T5224preIL1_vs_IL1_3h.csv", show_col_types = FALSE)

target_genes_3h <- res_il1_3h %>%
  dplyr::select(gene_name, logFC_IL1 = logFC, FDR_IL1 = FDR) %>%
  inner_join(
    res_rescue_3h %>% dplyr::select(gene_name, logFC_Rescue = logFC, FDR_Rescue = FDR),
    by = "gene_name"
  ) %>%
  mutate(
    Timepoint = "3 Hours",
    Status = case_when(
      FDR_IL1 < 0.05 & logFC_IL1 > 0.5 & FDR_Rescue < 0.05 & logFC_Rescue < -0.5 ~ "Dependent_Up",
      FDR_IL1 < 0.05 & logFC_IL1 < -0.5 & FDR_Rescue < 0.05 & logFC_Rescue > 0.5 ~ "Dependent_Down",
      FDR_IL1 < 0.05 & logFC_IL1 > 0.5 & (FDR_Rescue >= 0.05 | logFC_Rescue >= 0) ~ "Escaping_Up",
      FDR_IL1 < 0.05 & logFC_IL1 < -0.5 & (FDR_Rescue >= 0.05 | logFC_Rescue <= 0) ~ "Escaping_Down",
      TRUE ~ "Not_Primary_Driver"
    )
  ) %>%
  filter(Status != "Not_Primary_Driver")

# 2. Bind the 3h and 6h target lists together
combined_targets <- bind_rows(target_genes_3h, target_genes_6h) %>%
  # Ensure factors are ordered chronologically for plotting
  mutate(Timepoint = factor(Timepoint, levels = c("3 Hours", "6 Hours")))

# 3. Generate the Side-by-Side Faceted Plot
p_combined_quadrant <- ggplot(combined_targets, aes(x = logFC_IL1, y = logFC_Rescue, colour = Status)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
  scale_colour_manual(values = c(
    "Dependent_Up" = "firebrick3",
    "Dependent_Down" = "dodgerblue3",
    "Escaping_Up" = "lightpink",
    "Escaping_Down" = "lightblue"
  )) +
  facet_wrap(~ Timepoint, ncol = 2) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    strip.text = element_text(face = "bold", size = 13),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "IL-1 Response vs AP-1 Rescue Over Time",
    subtitle = "Tracking Dependent and Escaping gene populations at 3h and 6h",
    x = "Log2 Fold Change (IL-1 vs Vehicle)",
    y = "Log2 Fold Change (T5224+IL1 vs IL-1)",
    colour = "Gene Classification"
  )

ggsave("AP1_Rescue_Combined_Quadrant_Plot.png", p_combined_quadrant, width = 12, height = 7, dpi = 300)
print(p_combined_quadrant)

# ============================================================================
# PART 3: Identify Core Sustained AP-1 Targets (CORRECTED)
# ============================================================================

# Find genes that are successfully rescued (Dependent) at BOTH timepoints
temporal_tracking <- combined_targets %>%
  dplyr::select(gene_name, Timepoint, Status) %>%
  # Safely rename the timepoint strings BEFORE pivoting to avoid column name errors
  mutate(Timepoint = case_when(
    Timepoint == "3 Hours" ~ "Status_3h",
    Timepoint == "6 Hours" ~ "Status_6h"
  )) %>%
  pivot_wider(names_from = Timepoint, values_from = Status)

# Filter for genes strictly dependent at both timepoints
# We use replace_na() to safely handle genes that are only present at one timepoint
core_dependent <- temporal_tracking %>%
  filter(
    str_detect(replace_na(Status_3h, "None"), "Dependent") & 
    str_detect(replace_na(Status_6h, "None"), "Dependent")
  )

write_csv(core_dependent, "Core_Sustained_AP1_Dependent_Genes.csv")

cat("\n--- Temporal Overlap Summary ---\n")
cat("Total genes responding strongly to IL-1 in at least one timepoint:", nrow(temporal_tracking), "\n")
cat("Number of Core Sustained Targets (Dependent at BOTH 3h and 6h):", nrow(core_dependent), "\n")

DT::datatable(temporal_tracking, filter = 'top')

# Load required libraries
library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse)
library(enrichplot)

# 1. Filter the combined targets for ONLY the successfully rescued (Dependent) genes
# We use the 'combined_targets' dataframe generated in the previous step
dependent_genes <- combined_targets %>%
  filter(str_detect(Status, "Dependent")) %>%
  # Create a clean grouping name for the plot (e.g., "3h_Dependent_Up")
  mutate(
    Clean_Time = if_else(Timepoint == "3 Hours", "3h", "6h"),
    Group = paste(Clean_Time, Status, sep = "_")
  )

# 2. Map Gene Symbols to Entrez IDs for clusterProfiler
entrez_map <- bitr(unique(dependent_genes$gene_name), 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Mm.eg.db)

# Join the Entrez IDs back to our groups
dependent_mapped <- dependent_genes %>%
  inner_join(entrez_map, by = c("gene_name" = "SYMBOL"))

# 3. Create the named list of gene clusters required by compareCluster
# The split() function automatically breaks the Entrez IDs into separate vectors based on the 'Group' column
gene_list <- split(dependent_mapped$ENTREZID, dependent_mapped$Group)

# (Optional but recommended) Define your background universe again
# universe_entrez <- bitr(rownames(y), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# 4. Run the comparative GO Biological Process enrichment
comp_go_dependent <- compareCluster(
  geneCluster   = gene_list,
  fun           = "enrichGO",
  # universe      = universe_entrez$ENTREZID, # Uncomment if 'y' is loaded in your environment
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE # Translates Entrez IDs back to Gene Symbols in the CSV
)

# 5. Plot and Save the Results
if (!is.null(comp_go_dependent) && nrow(comp_go_dependent) > 0) {
  
  # Export the full data table
  write_csv(as.data.frame(comp_go_dependent), "GO_AP1_Dependent_Genes_3h_6h.csv")
  
  # Generate a comparative dotplot
  p_dep_go <- dotplot(comp_go_dependent, showCategory = 5) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11),
      axis.text.y = element_text(size = 11),
      plot.title = element_text(face = "bold", size = 14)
    ) +
    labs(
      title = "GO Terms: AP-1 Dependent (Rescued) Genes",
      subtitle = "Comparing pathways successfully rescued by T5224 at 3h vs 6h",
      x = ""
    )
  
  ggsave("GO_AP1_Dependent_Plot.png", p_dep_go, width = 12, height = 8, dpi = 300)
  print(p_dep_go)
  
} else {
  print("No significant GO terms were found for the dependent gene sets at these thresholds.")
}

# Load required libraries
library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse)
library(enrichplot)

# 1. Filter the combined targets for ONLY the Escaping genes
# We use the 'combined_targets' dataframe generated in the earlier combined analysis script
escaping_genes <- combined_targets %>%
  filter(str_detect(Status, "Escaping")) %>%
  mutate(
    Clean_Time = if_else(Timepoint == "3 Hours", "3h", "6h"),
    Group = paste(Clean_Time, Status, sep = "_")
  )

# 2. Map Gene Symbols to Entrez IDs for clusterProfiler
entrez_map <- bitr(unique(escaping_genes$gene_name), 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Mm.eg.db)

# Join the Entrez IDs back to our groups
escaping_mapped <- escaping_genes %>%
  inner_join(entrez_map, by = c("gene_name" = "SYMBOL"))

# 3. Create the named list of gene clusters required by compareCluster
gene_list <- split(escaping_mapped$ENTREZID, escaping_mapped$Group)

# 4. Run the comparative GO Biological Process enrichment
comp_go_escaping <- compareCluster(
  geneCluster   = gene_list,
  fun           = "enrichGO",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE 
)

# 5. Plot and Save the Results
if (!is.null(comp_go_escaping) && nrow(comp_go_escaping) > 0) {
  
  # Export the full data table
  write_csv(as.data.frame(comp_go_escaping), "GO_AP1_Escaping_Genes_3h_6h.csv")
  
  # Generate a comparative dotplot
  p_esc_go <- dotplot(comp_go_escaping, showCategory = 5) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11),
      axis.text.y = element_text(size = 11),
      plot.title = element_text(face = "bold", size = 14)
    ) +
    labs(
      title = "GO Terms: AP-1 Escaping (Independent) Genes",
      subtitle = "Pathways driven by IL-1 that resist T5224 inhibition at 3h vs 6h",
      x = ""
    )
  
  ggsave("GO_AP1_Escaping_Plot.png", p_esc_go, width = 12, height = 8, dpi = 300)
  print(p_esc_go)
  
} else {
  print("No significant GO terms were found for the escaping gene sets at these thresholds.")
}

# Load the extra package for statistical plotting
library(ggpubr)
library(tidyverse)

# ============================================================================
# PART 2: Merge Annotation, Calculate 'n', and Plot with Stats
# ============================================================================

# 1. Load the 3h results for all three relevant contrasts
res_il1 <- read_csv("DGE_Results/DGE_IL1_vs_Veh_3h.csv", show_col_types = FALSE)
res_rescue <- read_csv("DGE_Results/DGE_T5224preIL1_vs_IL1_3h.csv", show_col_types = FALSE)
res_t5224 <- read_csv("DGE_Results/DGE_T5224_vs_Veh_3h.csv", show_col_types = FALSE)

# 2. Merge data and append the base motif annotation
motif_df <- res_il1 %>%
  dplyr::select(gene_name, logFC_IL1 = logFC, FDR_IL1 = FDR) %>%
  inner_join(
    res_rescue %>% dplyr::select(gene_name, logFC_Rescue = logFC, FDR_Rescue = FDR),
    by = "gene_name"
  ) %>%
  inner_join(
    res_t5224 %>% dplyr::select(gene_name, logFC_T5224 = logFC, FDR_T5224 = FDR),
    by = "gene_name"
  ) %>%
  mutate(
    Base_Motif_Status = if_else(gene_name %in% ap1_target_genes, 
                                "AP-1 Motif Present", 
                                "No AP-1 Motif")
  )

# 3. Filter for significant genes and dynamically calculate the 'n'
filtered_df <- motif_df %>%
  filter(FDR_IL1 < 0.05 | FDR_T5224 < 0.05)

# Count how many genes fall into each category
counts_df <- filtered_df %>%
  count(Base_Motif_Status) %>%
  # Create a new label that includes the exact 'n'
  mutate(Motif_Label = paste0(Base_Motif_Status, "\n(n = ", n, ")"))

# 4. Format the data for plotting
plot_data <- filtered_df %>%
  # Bring in our new labels with the 'n' attached
  inner_join(counts_df, by = "Base_Motif_Status") %>%
  pivot_longer(
    cols = c(logFC_IL1, logFC_Rescue, logFC_T5224),
    names_to = "Contrast",
    values_to = "logFC"
  ) %>%
  mutate(
    Contrast = case_when(
      Contrast == "logFC_T5224" ~ "T5224 Baseline (Inhibitor vs Veh)",
      Contrast == "logFC_IL1" ~ "IL-1 Response (Stimulus)", 
      Contrast == "logFC_Rescue" ~ "T5224 Rescue (Inhibitor + IL-1)"
    ),
    Contrast = factor(Contrast, levels = c(
      "T5224 Baseline (Inhibitor vs Veh)",
      "IL-1 Response (Stimulus)",
      "T5224 Rescue (Inhibitor + IL-1)"
    ))
  )

# 5. Extract exact strings for colour mapping
label_present <- counts_df$Motif_Label[counts_df$Base_Motif_Status == "AP-1 Motif Present"]
label_missing <- counts_df$Motif_Label[counts_df$Base_Motif_Status == "No AP-1 Motif"]

# 6. Generate the overlaid Violin/Boxplot with Stats
p_motif <- ggplot(plot_data, aes(x = Motif_Label, y = logFC, fill = Motif_Label)) +
  geom_violin(alpha = 0.5, trim = FALSE, colour = NA) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA, colour = "black") +
  facet_wrap(~ Contrast, scales = "free_y") +
  
  # Apply colours to our dynamic labels
  scale_fill_manual(values = setNames(c("darkorange", "grey60"), 
                                      c(label_present, label_missing))) +
  
  coord_cartesian(ylim = c(-5, 6)) + # Slightly expanded upper y-limit for the p-value text
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey90"),
    axis.text.x = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  
  # Add the Statistical Test (Wilcoxon Rank-Sum / Mann-Whitney)
  stat_compare_means(
    method = "wilcox.test", 
    label = "p.format",      # Shows exact p-value (e.g., "p = 0.002"). Use "p.signif" if you just want asterisks (***)
    label.y = 5.2,           # The y-axis height where the p-value will be printed
    size = 4.5,
    hjust = 0.5
  ) +
  
  labs(
    title = "Impact of AP-1 Promoter Motifs on Expression Dynamics",
    subtitle = "Comparing significant responders (FDR < 0.05) across all 3h contrasts",
    x = "",
    y = "Log2 Fold Change"
  )

# Save and Display
ggsave("AP1_Motif_Distribution_Plot_All_Contrasts_Stats.png", p_motif, width = 14, height = 7, dpi = 300)
print(p_motif)
# Load the extra package for statistical plotting
library(ggpubr)
library(tidyverse)

# ============================================================================
# PART 2: Merge Annotation, Calculate 'n', and Plot with Stats
# ============================================================================

# 1. Load the 3h results for all three relevant contrasts
res_il1 <- read_csv("DGE_Results/DGE_IL1_vs_Veh_3h.csv", show_col_types = FALSE)
res_rescue <- read_csv("DGE_Results/DGE_T5224preIL1_vs_IL1_3h.csv", show_col_types = FALSE)
res_t5224 <- read_csv("DGE_Results/DGE_T5224_vs_Veh_3h.csv", show_col_types = FALSE)

# 2. Merge data and append the base motif annotation
motif_df <- res_il1 %>%
  dplyr::select(gene_name, logFC_IL1 = logFC, FDR_IL1 = FDR) %>%
  inner_join(
    res_rescue %>% dplyr::select(gene_name, logFC_Rescue = logFC, FDR_Rescue = FDR),
    by = "gene_name"
  ) %>%
  inner_join(
    res_t5224 %>% dplyr::select(gene_name, logFC_T5224 = logFC, FDR_T5224 = FDR),
    by = "gene_name"
  ) %>%
  mutate(
    Base_Motif_Status = if_else(gene_name %in% ap1_target_genes, 
                                "AP-1 Motif Present", 
                                "No AP-1 Motif")
  )

# 3. Filter for significant genes and dynamically calculate the 'n'
filtered_df <- motif_df #%>%
  #filter(FDR_IL1 < 0.05 | FDR_T5224 < 0.05)

# Count how many genes fall into each category
counts_df <- filtered_df %>%
  count(Base_Motif_Status) %>%
  # Create a new label that includes the exact 'n'
  mutate(Motif_Label = paste0(Base_Motif_Status, "\n(n = ", n, ")"))

# 4. Format the data for plotting
plot_data <- filtered_df %>%
  # Bring in our new labels with the 'n' attached
  inner_join(counts_df, by = "Base_Motif_Status") %>%
  pivot_longer(
    cols = c(logFC_IL1, logFC_Rescue, logFC_T5224),
    names_to = "Contrast",
    values_to = "logFC"
  ) %>%
  mutate(
    Contrast = case_when(
      Contrast == "logFC_T5224" ~ "T5224 Baseline (Inhibitor vs Veh)",
      Contrast == "logFC_IL1" ~ "IL-1 Response (Stimulus)", 
      Contrast == "logFC_Rescue" ~ "T5224 Rescue (Inhibitor + IL-1)"
    ),
    Contrast = factor(Contrast, levels = c(
      "T5224 Baseline (Inhibitor vs Veh)",
      "IL-1 Response (Stimulus)",
      "T5224 Rescue (Inhibitor + IL-1)"
    ))
  )

# 5. Extract exact strings for colour mapping
label_present <- counts_df$Motif_Label[counts_df$Base_Motif_Status == "AP-1 Motif Present"]
label_missing <- counts_df$Motif_Label[counts_df$Base_Motif_Status == "No AP-1 Motif"]

# 6. Generate the overlaid Violin/Boxplot with Stats
p_motif <- ggplot(plot_data, aes(x = Motif_Label, y = logFC, fill = Motif_Label)) +
  geom_violin(alpha = 0.5, trim = FALSE, colour = NA) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA, colour = "black") +
  facet_wrap(~ Contrast, scales = "free_y") +
  
  # Apply colours to our dynamic labels
  scale_fill_manual(values = setNames(c("darkorange", "grey60"), 
                                      c(label_present, label_missing))) +
  
  coord_cartesian(ylim = c(-5, 6)) + # Slightly expanded upper y-limit for the p-value text
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey90"),
    axis.text.x = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  
  # Add the Statistical Test (Wilcoxon Rank-Sum / Mann-Whitney)
  stat_compare_means(
    method = "wilcox.test", 
    label = "p.format",      # Shows exact p-value (e.g., "p = 0.002"). Use "p.signif" if you just want asterisks (***)
    label.y = 5.2,           # The y-axis height where the p-value will be printed
    size = 4.5,
    hjust = 0.5
  ) +
  
  labs(
    title = "Impact of AP-1 Promoter Motifs on Expression Dynamics",
    subtitle = "Comparing responders across all 3h contrasts",
    x = "",
    y = "Log2 Fold Change"
  )

# Save and Display
ggsave("AP1_Motif_Distribution_Plot_All_Contrasts_Stats.png", p_motif, width = 14, height = 7, dpi = 300)
print(p_motif)

# Load the extra package for statistical plotting
library(ggpubr)
library(tidyverse)

# ============================================================================
# PART 2: Merge Annotation, Calculate 'n', and Plot with Stats
# ============================================================================

# 1. Load the 3h results for all three relevant contrasts
res_il1 <- read_csv("DGE_Results/DGE_IL1_vs_Veh_6h.csv", show_col_types = FALSE)
res_rescue <- read_csv("DGE_Results/DGE_T5224preIL1_vs_IL1_6h.csv", show_col_types = FALSE)
res_t5224 <- read_csv("DGE_Results/DGE_T5224_vs_Veh_6h.csv", show_col_types = FALSE)

# 2. Merge data and append the base motif annotation
motif_df <- res_il1 %>%
  dplyr::select(gene_name, logFC_IL1 = logFC, FDR_IL1 = FDR) %>%
  inner_join(
    res_rescue %>% dplyr::select(gene_name, logFC_Rescue = logFC, FDR_Rescue = FDR),
    by = "gene_name"
  ) %>%
  inner_join(
    res_t5224 %>% dplyr::select(gene_name, logFC_T5224 = logFC, FDR_T5224 = FDR),
    by = "gene_name"
  ) %>%
  mutate(
    Base_Motif_Status = if_else(gene_name %in% ap1_target_genes, 
                                "AP-1 Motif Present", 
                                "No AP-1 Motif")
  )

# 3. Filter for significant genes and dynamically calculate the 'n'
filtered_df <- motif_df %>%
  filter(FDR_IL1 < 0.05 | FDR_T5224 < 0.05)

# Count how many genes fall into each category
counts_df <- filtered_df %>%
  count(Base_Motif_Status) %>%
  # Create a new label that includes the exact 'n'
  mutate(Motif_Label = paste0(Base_Motif_Status, "\n(n = ", n, ")"))

# 4. Format the data for plotting
plot_data <- filtered_df %>%
  # Bring in our new labels with the 'n' attached
  inner_join(counts_df, by = "Base_Motif_Status") %>%
  pivot_longer(
    cols = c(logFC_IL1, logFC_Rescue, logFC_T5224),
    names_to = "Contrast",
    values_to = "logFC"
  ) %>%
  mutate(
    Contrast = case_when(
      Contrast == "logFC_T5224" ~ "T5224 Baseline (Inhibitor vs Veh)",
      Contrast == "logFC_IL1" ~ "IL-1 Response (Stimulus)", 
      Contrast == "logFC_Rescue" ~ "T5224 Rescue (Inhibitor + IL-1)"
    ),
    Contrast = factor(Contrast, levels = c(
      "T5224 Baseline (Inhibitor vs Veh)",
      "IL-1 Response (Stimulus)",
      "T5224 Rescue (Inhibitor + IL-1)"
    ))
  )

# 5. Extract exact strings for colour mapping
label_present <- counts_df$Motif_Label[counts_df$Base_Motif_Status == "AP-1 Motif Present"]
label_missing <- counts_df$Motif_Label[counts_df$Base_Motif_Status == "No AP-1 Motif"]

# 6. Generate the overlaid Violin/Boxplot with Stats
p_motif <- ggplot(plot_data, aes(x = Motif_Label, y = logFC, fill = Motif_Label)) +
  geom_violin(alpha = 0.5, trim = FALSE, colour = NA) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA, colour = "black") +
  facet_wrap(~ Contrast, scales = "free_y") +
  
  # Apply colours to our dynamic labels
  scale_fill_manual(values = setNames(c("darkorange", "grey60"), 
                                      c(label_present, label_missing))) +
  
  coord_cartesian(ylim = c(-5, 6)) + # Slightly expanded upper y-limit for the p-value text
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey90"),
    axis.text.x = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  
  # Add the Statistical Test (Wilcoxon Rank-Sum / Mann-Whitney)
  stat_compare_means(
    method = "wilcox.test", 
    label = "p.format",      # Shows exact p-value (e.g., "p = 0.002"). Use "p.signif" if you just want asterisks (***)
    label.y = 5.2,           # The y-axis height where the p-value will be printed
    size = 4.5,
    hjust = 0.5
  ) +
  
  labs(
    title = "Impact of AP-1 Promoter Motifs on Expression Dynamics",
    subtitle = "Comparing significant responders (FDR < 0.05) across all 6h contrasts",
    x = "",
    y = "Log2 Fold Change"
  )

# Save and Display
ggsave("AP1_Motif_Distribution_Plot_All_Contrasts_Stats.png", p_motif, width = 14, height = 7, dpi = 300)
print(p_motif)
# Load the extra package for statistical plotting
library(ggpubr)
library(tidyverse)

# ============================================================================
# PART 2: Merge Annotation, Calculate 'n', and Plot with Stats
# ============================================================================

# 1. Load the 3h results for all three relevant contrasts
res_il1 <- read_csv("DGE_Results/DGE_IL1_vs_Veh_6h.csv", show_col_types = FALSE)
res_rescue <- read_csv("DGE_Results/DGE_T5224preIL1_vs_IL1_6h.csv", show_col_types = FALSE)
res_t5224 <- read_csv("DGE_Results/DGE_T5224_vs_Veh_6h.csv", show_col_types = FALSE)

# 2. Merge data and append the base motif annotation
motif_df <- res_il1 %>%
  dplyr::select(gene_name, logFC_IL1 = logFC, FDR_IL1 = FDR) %>%
  inner_join(
    res_rescue %>% dplyr::select(gene_name, logFC_Rescue = logFC, FDR_Rescue = FDR),
    by = "gene_name"
  ) %>%
  inner_join(
    res_t5224 %>% dplyr::select(gene_name, logFC_T5224 = logFC, FDR_T5224 = FDR),
    by = "gene_name"
  ) %>%
  mutate(
    Base_Motif_Status = if_else(gene_name %in% ap1_target_genes, 
                                "AP-1 Motif Present", 
                                "No AP-1 Motif")
  )

# 3. Filter for significant genes and dynamically calculate the 'n'
filtered_df <- motif_df #%>%
  #filter(FDR_IL1 < 0.05 | FDR_T5224 < 0.05)

# Count how many genes fall into each category
counts_df <- filtered_df %>%
  count(Base_Motif_Status) %>%
  # Create a new label that includes the exact 'n'
  mutate(Motif_Label = paste0(Base_Motif_Status, "\n(n = ", n, ")"))

# 4. Format the data for plotting
plot_data <- filtered_df %>%
  # Bring in our new labels with the 'n' attached
  inner_join(counts_df, by = "Base_Motif_Status") %>%
  pivot_longer(
    cols = c(logFC_IL1, logFC_Rescue, logFC_T5224),
    names_to = "Contrast",
    values_to = "logFC"
  ) %>%
  mutate(
    Contrast = case_when(
      Contrast == "logFC_T5224" ~ "T5224 Baseline (Inhibitor vs Veh)",
      Contrast == "logFC_IL1" ~ "IL-1 Response (Stimulus)", 
      Contrast == "logFC_Rescue" ~ "T5224 Rescue (Inhibitor + IL-1)"
    ),
    Contrast = factor(Contrast, levels = c(
      "T5224 Baseline (Inhibitor vs Veh)",
      "IL-1 Response (Stimulus)",
      "T5224 Rescue (Inhibitor + IL-1)"
    ))
  )

# 5. Extract exact strings for colour mapping
label_present <- counts_df$Motif_Label[counts_df$Base_Motif_Status == "AP-1 Motif Present"]
label_missing <- counts_df$Motif_Label[counts_df$Base_Motif_Status == "No AP-1 Motif"]

# 6. Generate the overlaid Violin/Boxplot with Stats
p_motif <- ggplot(plot_data, aes(x = Motif_Label, y = logFC, fill = Motif_Label)) +
  geom_violin(alpha = 0.5, trim = FALSE, colour = NA) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA, colour = "black") +
  facet_wrap(~ Contrast, scales = "free_y") +
  
  # Apply colours to our dynamic labels
  scale_fill_manual(values = setNames(c("darkorange", "grey60"), 
                                      c(label_present, label_missing))) +
  
  coord_cartesian(ylim = c(-5, 6)) + # Slightly expanded upper y-limit for the p-value text
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey90"),
    axis.text.x = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  
  # Add the Statistical Test (Wilcoxon Rank-Sum / Mann-Whitney)
  stat_compare_means(
    method = "wilcox.test", 
    label = "p.format",      # Shows exact p-value (e.g., "p = 0.002"). Use "p.signif" if you just want asterisks (***)
    label.y = 5.2,           # The y-axis height where the p-value will be printed
    size = 4.5,
    hjust = 0.5
  ) +
  
  labs(
    title = "Impact of AP-1 Promoter Motifs on Expression Dynamics",
    subtitle = "Comparing  responders across all 6h contrasts",
    x = "",
    y = "Log2 Fold Change"
  )

# Save and Display
ggsave("AP1_Motif_Distribution_Plot_All_Contrasts_Stats.png", p_motif, width = 14, height = 7, dpi = 300)
print(p_motif)

# head(All_TOBIAS_BINDING)  # Removed: All_TOBIAS_BINDING is first created in "PART 2" below, so this call failed when the script was run top to bottom.

# Load required libraries
library(tidyverse)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(clusterProfiler)

# ============================================================================
# PART 1: Define mm10 Promoters
# ============================================================================
print("Defining mm10 promoter regions...")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
# 2000bp upstream to 500bp downstream of the Transcription Start Site
promoters_gr <- promoters(genes(txdb), upstream = 2000, downstream = 500)

# ============================================================================
# PART 2: Load and Filter TOBIAS scATAC-seq Data
# ============================================================================
print("Loading TOBIAS RDS file...")
All_TOBIAS_BINDING <- readRDS(publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Output/Tobias_df_Classified_Full.rds"))

# Apply your cell type standardisation
All_TOBIAS_BINDING$cell <- All_TOBIAS_BINDING$CellType

# Filter for active AP-1 footprints specifically in Fibroblasts
fibro_ap1_footprints <- All_TOBIAS_BINDING %>%
  # Target the fibroblast population
  filter(grepl("Fibroblast", cell, ignore.case = TRUE)) %>%
  # Filter for AP-1 complex members (Update 'TF_name' to match your column name if needed)
  filter(grepl("Jun|Fos|Batf|Atf", name, ignore.case = TRUE)) 

print(paste("Found", nrow(fibro_ap1_footprints), "active AP-1 footprints in fibroblasts."))

# ============================================================================
# PART 3: Overlap Footprints with Promoters
# ============================================================================
# Convert the TOBIAS dataframe into a GenomicRanges object
# (Update 'chr', 'start', 'end' if your TOBIAS columns use TFBS_chr, TFBS_start, etc.)
tobias_gr <- GRanges(
  seqnames = fibro_ap1_footprints$seqnames,
  ranges = IRanges(start = fibro_ap1_footprints$start, end = fibro_ap1_footprints$end)
)

# Find which footprints physically overlap with our defined promoters
overlaps <- findOverlaps(tobias_gr, promoters_gr)

# Extract the Entrez IDs of the genes whose promoters had an overlapping footprint
entrez_with_empirical_motif <- unique(names(promoters_gr)[subjectHits(overlaps)])

# Convert Entrez IDs to Gene Symbols to match the RNA-seq data
symbol_map <- bitr(entrez_with_empirical_motif, 
                   fromType = "ENTREZID", 
                   toType = "SYMBOL", 
                   OrgDb = org.Mm.eg.db)

empirical_ap1_targets <- unique(symbol_map$SYMBOL)

print(paste("Identified", length(empirical_ap1_targets), "genes with an empirical AP-1 footprint in their promoter."))

# ============================================================================
# PART 4: Merge with DGE Data and Plot
# ============================================================================

# Load the extra package for statistical plotting
library(ggpubr)
library(tidyverse)

# ============================================================================
# PART 4: Merge ATAC Annotation, Calculate 'n', and Plot with Stats
# ============================================================================

# 1. Load the 3h results for all three relevant contrasts
res_il1 <- read_csv("DGE_Results/DGE_IL1_vs_Veh_3h.csv", show_col_types = FALSE)
res_rescue <- read_csv("DGE_Results/DGE_T5224preIL1_vs_IL1_3h.csv", show_col_types = FALSE)
res_t5224 <- read_csv("DGE_Results/DGE_T5224_vs_Veh_3h.csv", show_col_types = FALSE)

# 2. Merge data and append the base empirical footprint annotation
footprint_df <- res_il1 %>%
  dplyr::select(gene_name, logFC_IL1 = logFC, FDR_IL1 = FDR) %>%
  inner_join(
    res_rescue %>% dplyr::select(gene_name, logFC_Rescue = logFC, FDR_Rescue = FDR),
    by = "gene_name"
  ) %>%
  inner_join(
    res_t5224 %>% dplyr::select(gene_name, logFC_T5224 = logFC, FDR_T5224 = FDR),
    by = "gene_name"
  ) %>%
  mutate(
    Base_Motif_Status = if_else(gene_name %in% empirical_ap1_targets, 
                                "Bound AP-1 Footprint (ATAC)", 
                                "No AP-1 Footprint")
  )

# 3. Filter for significant genes and dynamically calculate the 'n'
filtered_df <- footprint_df %>%
  filter(FDR_IL1 < 0.05 | FDR_T5224 < 0.05)

# Count how many genes fall into each category
counts_df <- filtered_df %>%
  count(Base_Motif_Status) %>%
  # Create a new label that includes the exact 'n'
  mutate(Motif_Label = paste0(Base_Motif_Status, "\n(n = ", n, ")"))

# 4. Format the data for plotting
plot_data <- filtered_df %>%
  # Bring in our new labels with the 'n' attached
  inner_join(counts_df, by = "Base_Motif_Status") %>%
  pivot_longer(
    cols = c(logFC_IL1, logFC_Rescue, logFC_T5224),
    names_to = "Contrast",
    values_to = "logFC"
  ) %>%
  mutate(
    Contrast = case_when(
      Contrast == "logFC_T5224" ~ "T5224 Baseline (Inhibitor vs Veh)",
      Contrast == "logFC_IL1" ~ "IL-1 Response (Stimulus)", 
      Contrast == "logFC_Rescue" ~ "T5224 Rescue (Inhibitor + IL-1)"
    ),
    # Enforce a logical plotting order from left to right
    Contrast = factor(Contrast, levels = c(
      "T5224 Baseline (Inhibitor vs Veh)",
      "IL-1 Response (Stimulus)",
      "T5224 Rescue (Inhibitor + IL-1)"
    ))
  )

# 5. Extract exact strings for colour mapping
label_bound <- counts_df$Motif_Label[counts_df$Base_Motif_Status == "Bound AP-1 Footprint (ATAC)"]
label_unbound <- counts_df$Motif_Label[counts_df$Base_Motif_Status == "No AP-1 Footprint"]

# 6. Generate the overlaid Violin/Boxplot with Stats
p_footprint <- ggplot(plot_data, aes(x = Motif_Label, y = logFC, fill = Motif_Label)) +
  geom_violin(alpha = 0.5, trim = FALSE, colour = NA) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA, colour = "black") +
  # Now facet across all 3 contrasts
  facet_wrap(~ Contrast, scales = "free_y") +
  
  # Apply our purple/grey colours to the dynamically generated labels
  scale_fill_manual(values = setNames(c("purple3", "grey60"), 
                                      c(label_bound, label_unbound))) +
  
  # Slightly expand upper y-limit to ensure the p-value text fits
  coord_cartesian(ylim = c(-5, 6.5)) + 
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey90"),
    axis.text.x = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  
  # Add the Statistical Test
  stat_compare_means(
    method = "wilcox.test", 
    label = "p.format",      
    label.y = 5.5,           # Anchor height for the p-value text
    size = 4.5,
    hjust = 0.5
  ) +
  
  labs(
    title = "Impact of Empirical AP-1 Footprints on Gene Expression",
    subtitle = "Integrating scATAC-seq TOBIAS binding with bulk RNA-seq dynamics (3h)",
    x = "",
    y = "Log2 Fold Change"
  )

# Save and Display (Width increased to 14 for 3 panels)
ggsave("Empirical_AP1_Footprint_Distribution_Stats.png", p_footprint, width = 14, height = 7, dpi = 300)
print(p_footprint)

# ============================================================================
# PART 3: Overlap Footprints with Promoters
# ============================================================================
# Convert the TOBIAS dataframe into a GenomicRanges object
# (Update 'chr', 'start', 'end' if your TOBIAS columns use TFBS_chr, TFBS_start, etc.)
tobias_gr <- GRanges(
  seqnames = fibro_ap1_footprints$seqnames,
  ranges = IRanges(start = fibro_ap1_footprints$start, end = fibro_ap1_footprints$end)
)

# Find which footprints physically overlap with our defined promoters
overlaps <- findOverlaps(tobias_gr, promoters_gr)

# Extract the Entrez IDs of the genes whose promoters had an overlapping footprint
entrez_with_empirical_motif <- unique(names(promoters_gr)[subjectHits(overlaps)])

# Convert Entrez IDs to Gene Symbols to match the RNA-seq data
symbol_map <- bitr(entrez_with_empirical_motif, 
                   fromType = "ENTREZID", 
                   toType = "SYMBOL", 
                   OrgDb = org.Mm.eg.db)

empirical_ap1_targets <- unique(symbol_map$SYMBOL)

print(paste("Identified", length(empirical_ap1_targets), "genes with an empirical AP-1 footprint in their promoter."))

# ============================================================================
# PART 4: Merge with DGE Data and Plot
# ============================================================================

# Load the extra package for statistical plotting
library(ggpubr)
library(tidyverse)

# ============================================================================
# PART 4: Merge ATAC Annotation, Calculate 'n', and Plot with Stats
# ============================================================================

# 1. Load the 3h results for all three relevant contrasts
res_il1 <- read_csv("DGE_Results/DGE_IL1_vs_Veh_3h.csv", show_col_types = FALSE)
res_rescue <- read_csv("DGE_Results/DGE_T5224preIL1_vs_IL1_3h.csv", show_col_types = FALSE)
res_t5224 <- read_csv("DGE_Results/DGE_T5224_vs_Veh_3h.csv", show_col_types = FALSE)

# 2. Merge data and append the base empirical footprint annotation
footprint_df <- res_il1 %>%
  dplyr::select(gene_name, logFC_IL1 = logFC, FDR_IL1 = FDR) %>%
  inner_join(
    res_rescue %>% dplyr::select(gene_name, logFC_Rescue = logFC, FDR_Rescue = FDR),
    by = "gene_name"
  ) %>%
  inner_join(
    res_t5224 %>% dplyr::select(gene_name, logFC_T5224 = logFC, FDR_T5224 = FDR),
    by = "gene_name"
  ) %>%
  mutate(
    Base_Motif_Status = if_else(gene_name %in% empirical_ap1_targets, 
                                "Bound AP-1 Footprint (ATAC)", 
                                "No AP-1 Footprint")
  )

# 3. Filter for significant genes and dynamically calculate the 'n'
filtered_df <- footprint_df #%>%
  #filter(FDR_IL1 < 0.05 | FDR_T5224 < 0.05)

# Count how many genes fall into each category
counts_df <- filtered_df %>%
  count(Base_Motif_Status) %>%
  # Create a new label that includes the exact 'n'
  mutate(Motif_Label = paste0(Base_Motif_Status, "\n(n = ", n, ")"))

# 4. Format the data for plotting
plot_data <- filtered_df %>%
  # Bring in our new labels with the 'n' attached
  inner_join(counts_df, by = "Base_Motif_Status") %>%
  pivot_longer(
    cols = c(logFC_IL1, logFC_Rescue, logFC_T5224),
    names_to = "Contrast",
    values_to = "logFC"
  ) %>%
  mutate(
    Contrast = case_when(
      Contrast == "logFC_T5224" ~ "T5224 Baseline (Inhibitor vs Veh)",
      Contrast == "logFC_IL1" ~ "IL-1 Response (Stimulus)", 
      Contrast == "logFC_Rescue" ~ "T5224 Rescue (Inhibitor + IL-1)"
    ),
    # Enforce a logical plotting order from left to right
    Contrast = factor(Contrast, levels = c(
      "T5224 Baseline (Inhibitor vs Veh)",
      "IL-1 Response (Stimulus)",
      "T5224 Rescue (Inhibitor + IL-1)"
    ))
  )

# 5. Extract exact strings for colour mapping
label_bound <- counts_df$Motif_Label[counts_df$Base_Motif_Status == "Bound AP-1 Footprint (ATAC)"]
label_unbound <- counts_df$Motif_Label[counts_df$Base_Motif_Status == "No AP-1 Footprint"]

# 6. Generate the overlaid Violin/Boxplot with Stats
p_footprint <- ggplot(plot_data, aes(x = Motif_Label, y = logFC, fill = Motif_Label)) +
  geom_violin(alpha = 0.5, trim = FALSE, colour = NA) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA, colour = "black") +
  # Now facet across all 3 contrasts
  facet_wrap(~ Contrast, scales = "free_y") +
  
  # Apply our purple/grey colours to the dynamically generated labels
  scale_fill_manual(values = setNames(c("purple3", "grey60"), 
                                      c(label_bound, label_unbound))) +
  
  # Slightly expand upper y-limit to ensure the p-value text fits
  coord_cartesian(ylim = c(-5, 6.5)) + 
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey90"),
    axis.text.x = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  
  # Add the Statistical Test
  stat_compare_means(
    method = "wilcox.test", 
    label = "p.format",      
    label.y = 5.5,           # Anchor height for the p-value text
    size = 4.5,
    hjust = 0.5
  ) +
  
  labs(
    title = "Impact of Empirical AP-1 Footprints on Gene Expression",
    subtitle = "Integrating scATAC-seq TOBIAS binding with bulk RNA-seq dynamics (3h)",
    x = "",
    y = "Log2 Fold Change"
  )

# Save and Display (Width increased to 14 for 3 panels)
ggsave("Empirical_AP1_Footprint_Distribution_Stats.png", p_footprint, width = 14, height = 7, dpi = 300)
print(p_footprint)

# Load required libraries
library(ggpubr)
library(tidyverse)

# ============================================================================
# PART 1: Load Data and Annotate Footprints
# ============================================================================

# 1. Load ONLY the T5224 vs Vehicle baseline contrast (3h)
res_t5224 <- read_csv("DGE_Results/DGE_T5224_vs_Veh_3h.csv", show_col_types = FALSE)

# 2. Annotate the genome with the TOBIAS ATAC footprint data
anno_t5224 <- res_t5224 %>%
  # Remove any genes with NA logFC to ensure clean math
  filter(!is.na(logFC)) %>%
  mutate(
    Base_Motif_Status = if_else(gene_name %in% empirical_ap1_targets, 
                                "Bound AP-1 Footprint (ATAC)", 
                                "No AP-1 Footprint")
  )

# ============================================================================
# PART 2 & 3: The Bulletproof Stats Fix
# ============================================================================

# 1. Create the two data subsets
df_all <- anno_t5224 %>% mutate(Scope = "All Detected Genes")
df_sig <- anno_t5224 %>% filter(FDR < 0.05) %>% mutate(Scope = "Significant Genes Only")

# 2. Combine and calculate 'n' for the facet headers
# By moving the counts here, the x-axis remains perfectly uniform for ggpubr
plot_data <- bind_rows(df_all, df_sig) %>%
  # Lock the logical plotting order
  mutate(Scope = factor(Scope, levels = c("All Detected Genes", "Significant Genes Only"))) %>%
  arrange(Scope) %>%
  group_by(Scope) %>%
  mutate(
    # Count the instances within each facet
    n_bound = sum(Base_Motif_Status == "Bound AP-1 Footprint (ATAC)"),
    n_unbound = sum(Base_Motif_Status == "No AP-1 Footprint"),
    # Create a dynamic facet title containing the numbers
    Facet_Label = paste0(Scope, "\n(Bound n = ", n_bound, " | Unbound n = ", n_unbound, ")")
  ) %>%
  ungroup() %>%
  # fct_inorder ensures the facets plot left-to-right exactly as we arranged them
  mutate(Facet_Label = fct_inorder(Facet_Label))

# 3. Define the single, universal comparison
my_comparisons <- list(c("Bound AP-1 Footprint (ATAC)", "No AP-1 Footprint"))

# 4. Generate the Plot
p_validation <- ggplot(plot_data, aes(x = Base_Motif_Status, y = logFC, fill = Base_Motif_Status)) +
  geom_violin(alpha = 0.5, trim = FALSE, colour = NA) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA, colour = "black") +
  
  # Facet by our new dynamic header
  facet_wrap(~ Facet_Label) + 
  
  # Use static colours since the x-axis categories are now uniform
  scale_fill_manual(values = c("Bound AP-1 Footprint (ATAC)" = "purple3", 
                               "No AP-1 Footprint" = "grey60")) +
  
  coord_cartesian(ylim = c(-4, 4.5)) + 
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey90"),
    axis.text.x = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  
  # Add Wilcoxon Rank-Sum Test
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test", 
    label = "p.format",      
    label.y = 3.8,           
    size = 4.5
  ) +
  
  labs(
    title = "Validation of T5224 Inhibitor Mechanism of Action",
    subtitle = "Comparing basal transcription shifts at physically bound AP-1 motifs (T5224 vs Vehicle)",
    x = "",
    y = "Log2 Fold Change (T5224 vs Vehicle)"
  )

# Save and Display
ggsave("T5224_Validation_Global_vs_Sig_Fixed.png", p_validation, width = 12, height = 7, dpi = 300)
print(p_validation)

# Load required libraries
library(ggpubr)
library(tidyverse)

# ============================================================================
# PART 1: Load Data and Annotate Footprints
# ============================================================================

# 1. Load ONLY the T5224 vs Vehicle baseline contrast (3h)
res_t5224 <- read_csv("DGE_Results/DGE_T5224_vs_Veh_6h.csv", show_col_types = FALSE)

# 2. Annotate the genome with the TOBIAS ATAC footprint data
anno_t5224 <- res_t5224 %>%
  # Remove any genes with NA logFC to ensure clean math
  filter(!is.na(logFC)) %>%
  mutate(
    Base_Motif_Status = if_else(gene_name %in% empirical_ap1_targets, 
                                "Bound AP-1 Footprint (ATAC)", 
                                "No AP-1 Footprint")
  )

# ============================================================================
# PART 2 & 3: The Bulletproof Stats Fix
# ============================================================================

# 1. Create the two data subsets
df_all <- anno_t5224 %>% mutate(Scope = "All Detected Genes")
df_sig <- anno_t5224 %>% filter(FDR < 0.05) %>% mutate(Scope = "Significant Genes Only")

# 2. Combine and calculate 'n' for the facet headers
# By moving the counts here, the x-axis remains perfectly uniform for ggpubr
plot_data <- bind_rows(df_all, df_sig) %>%
  # Lock the logical plotting order
  mutate(Scope = factor(Scope, levels = c("All Detected Genes", "Significant Genes Only"))) %>%
  arrange(Scope) %>%
  group_by(Scope) %>%
  mutate(
    # Count the instances within each facet
    n_bound = sum(Base_Motif_Status == "Bound AP-1 Footprint (ATAC)"),
    n_unbound = sum(Base_Motif_Status == "No AP-1 Footprint"),
    # Create a dynamic facet title containing the numbers
    Facet_Label = paste0(Scope, "\n(Bound n = ", n_bound, " | Unbound n = ", n_unbound, ")")
  ) %>%
  ungroup() %>%
  # fct_inorder ensures the facets plot left-to-right exactly as we arranged them
  mutate(Facet_Label = fct_inorder(Facet_Label))

# 3. Define the single, universal comparison
my_comparisons <- list(c("Bound AP-1 Footprint (ATAC)", "No AP-1 Footprint"))

# 4. Generate the Plot
p_validation <- ggplot(plot_data, aes(x = Base_Motif_Status, y = logFC, fill = Base_Motif_Status)) +
  geom_violin(alpha = 0.5, trim = FALSE, colour = NA) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA, colour = "black") +
  
  # Facet by our new dynamic header
  facet_wrap(~ Facet_Label) + 
  
  # Use static colours since the x-axis categories are now uniform
  scale_fill_manual(values = c("Bound AP-1 Footprint (ATAC)" = "purple3", 
                               "No AP-1 Footprint" = "grey60")) +
  
  coord_cartesian(ylim = c(-4, 4.5)) + 
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey90"),
    axis.text.x = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  
  # Add Wilcoxon Rank-Sum Test
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test", 
    label = "p.format",      
    label.y = 3.8,           
    size = 4.5
  ) +
  
  labs(
    title = "Validation of T5224 Inhibitor Mechanism of Action",
    subtitle = "Comparing basal transcription shifts at physically bound AP-1 motifs (T5224 vs Vehicle)",
    x = "",
    y = "Log2 Fold Change (T5224 vs Vehicle)"
  )

# Save and Display
ggsave("T5224_Validation_Global_vs_Sig_Fixed.png", p_validation, width = 12, height = 7, dpi = 300)
print(p_validation)

# Load required libraries
library(tidyverse)

# 1. Load the baseline inhibitor data (assuming it's not already in memory)
res_t5224 <- read_csv("DGE_Results/DGE_T5224_vs_Veh_3h.csv", show_col_types = FALSE)

# 2. Categorise the data for the contingency table
fisher_df <- res_t5224 %>%
  filter(!is.na(logFC) & !is.na(FDR)) %>%
  mutate(
    # Bucket 1: The biological response we care about
    Response = case_when(
      FDR < 0.05 & logFC < 0 ~ "Significantly Downregulated",
      TRUE ~ "Other (Unchanged or Up)"
    ),
    # Bucket 2: The physical binding data (from our earlier TOBIAS steps)
    Motif = if_else(gene_name %in% empirical_ap1_targets, 
                    "Bound Footprint", 
                    "No Footprint")
  ) %>%
  # Convert to factors to strictly control the matrix orientation for the Odds Ratio
  mutate(
    Response = factor(Response, levels = c("Significantly Downregulated", "Other (Unchanged or Up)")),
    Motif = factor(Motif, levels = c("Bound Footprint", "No Footprint"))
  )

# 3. Create the Contingency Table
contingency_table <- table(fisher_df$Motif, fisher_df$Response)

cat("\n--- Contingency Table ---\n")
print(contingency_table)

# 4. Run the Fisher's Exact Test
fisher_res <- fisher.test(contingency_table)

cat("\n--- Fisher's Exact Test Results ---\n")
print(fisher_res)

# 5. Prepare data for a proportional bar chart
plot_data <- fisher_df %>%
  count(Response, Motif) %>%
  group_by(Response) %>%
  # Calculate percentages for the labels
  mutate(
    Proportion = n / sum(n),
    Label = paste0(round(Proportion * 100, 1), "%\n(n=", n, ")")
  )

# 6. Generate the Stacked Bar Plot
p_fisher <- ggplot(plot_data, aes(x = Response, y = Proportion, fill = Motif)) +
  # position = "fill" forces the bars to 100% height for direct visual comparison
  geom_bar(stat = "identity", position = "fill", colour = "black", width = 0.6) +
  # Add the text labels inside the bars
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), 
            colour = "white", fontface = "bold", size = 4.5) +
  scale_fill_manual(values = c("Bound Footprint" = "purple3", "No Footprint" = "grey60")) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 11),
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  # Convert y-axis to percentage
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Enrichment of AP-1 Footprints in Downregulated Genes",
    subtitle = sprintf("Fisher's Exact Test: p = %.2e | Odds Ratio = %.2f", 
                       fisher_res$p.value, fisher_res$estimate),
    x = "",
    y = "Proportion of Genes"
  )

# Save and Display
ggsave("Fisher_Test_AP1_Enrichment.png", p_fisher, width = 9, height = 7, dpi = 300)
print(p_fisher)

# Load required libraries
library(tidyverse)

# 1. Load the baseline inhibitor data (assuming it's not already in memory)
res_t5224 <- read_csv("DGE_Results/DGE_T5224_vs_Veh_6h.csv", show_col_types = FALSE)

# 2. Categorise the data for the contingency table
fisher_df <- res_t5224 %>%
  filter(!is.na(logFC) & !is.na(FDR)) %>%
  mutate(
    # Bucket 1: The biological response we care about
    Response = case_when(
      FDR < 0.05 & logFC < 0 ~ "Significantly Downregulated",
      TRUE ~ "Other (Unchanged or Up)"
    ),
    # Bucket 2: The physical binding data (from our earlier TOBIAS steps)
    Motif = if_else(gene_name %in% empirical_ap1_targets, 
                    "Bound Footprint", 
                    "No Footprint")
  ) %>%
  # Convert to factors to strictly control the matrix orientation for the Odds Ratio
  mutate(
    Response = factor(Response, levels = c("Significantly Downregulated", "Other (Unchanged or Up)")),
    Motif = factor(Motif, levels = c("Bound Footprint", "No Footprint"))
  )

# 3. Create the Contingency Table
contingency_table <- table(fisher_df$Motif, fisher_df$Response)

cat("\n--- Contingency Table ---\n")
print(contingency_table)

# 4. Run the Fisher's Exact Test
fisher_res <- fisher.test(contingency_table)

cat("\n--- Fisher's Exact Test Results ---\n")
print(fisher_res)

# 5. Prepare data for a proportional bar chart
plot_data <- fisher_df %>%
  count(Response, Motif) %>%
  group_by(Response) %>%
  # Calculate percentages for the labels
  mutate(
    Proportion = n / sum(n),
    Label = paste0(round(Proportion * 100, 1), "%\n(n=", n, ")")
  )

# 6. Generate the Stacked Bar Plot
p_fisher <- ggplot(plot_data, aes(x = Response, y = Proportion, fill = Motif)) +
  # position = "fill" forces the bars to 100% height for direct visual comparison
  geom_bar(stat = "identity", position = "fill", colour = "black", width = 0.6) +
  # Add the text labels inside the bars
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), 
            colour = "white", fontface = "bold", size = 4.5) +
  scale_fill_manual(values = c("Bound Footprint" = "purple3", "No Footprint" = "grey60")) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 11),
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  # Convert y-axis to percentage
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Enrichment of AP-1 Footprints in Downregulated Genes",
    subtitle = sprintf("Fisher's Exact Test: p = %.2e | Odds Ratio = %.2f", 
                       fisher_res$p.value, fisher_res$estimate),
    x = "",
    y = "Proportion of Genes"
  )

# Save and Display
ggsave("Fisher_Test_AP1_Enrichment.png", p_fisher, width = 9, height = 7, dpi = 300)
print(p_fisher)

library(tidyverse)

# 1. Define the Feedback Loop Validation Genes
feedback_genes <- c(
  "Jun", "Fos",
  "Dusp1", "Dusp2", "Dusp8", "Dusp10", 
  "Il1r1", "Il1rap", "Myd88", "Mapk8", "Mapk9",
  "Nfkbia", "Tnfaip3"
)

# 2. Load the 3h results (Selecting BOTH logFC and FDR this time)
res_t5224 <- read_csv("DGE_Results/DGE_T5224_vs_Veh_3h.csv", show_col_types = FALSE) %>%
  dplyr::select(gene_name, logFC_T5224 = logFC, FDR_T5224 = FDR)

res_il1 <- read_csv("DGE_Results/DGE_IL1_vs_Veh_3h.csv", show_col_types = FALSE) %>%
  dplyr::select(gene_name, logFC_IL1 = logFC, FDR_IL1 = FDR)

res_rescue <- read_csv("DGE_Results/DGE_T5224preIL1_vs_IL1_3h.csv", show_col_types = FALSE) %>%
  dplyr::select(gene_name, logFC_Rescue = logFC, FDR_Rescue = FDR)

# 3. Merge, Categorise, and Pivot
feedback_df <- res_t5224 %>%
  inner_join(res_il1, by = "gene_name") %>%
  inner_join(res_rescue, by = "gene_name") %>%
  filter(gene_name %in% feedback_genes) %>%
  mutate(
    Category = case_when(
      gene_name %in% c("Jun", "Fos") ~ "1. AP-1 Core Factors",
      grepl("Dusp", gene_name) ~ "2. MAPK Phosphatases (Brakes)",
      grepl("Il1|Myd|Mapk", gene_name) ~ "3. Upstream Amplifiers",
      grepl("Nfk", gene_name) | gene_name == "Tnfaip3" ~ "4. Crosstalk Regulators",
      TRUE ~ "Other"
    )
  ) %>%
  # Advanced pivot: This smoothly handles BOTH logFC and FDR columns at once
  pivot_longer(
    cols = -c(gene_name, Category),
    names_to = c(".value", "Contrast"),
    names_sep = "_"
  ) %>%
  # Translate the FDR math into standard scientific asterisks
  mutate(
    Significance = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE        ~ "" # Leave blank if not significant
    ),
    Contrast = case_when(
      Contrast == "T5224" ~ "1. Baseline\n(T5224 vs Veh)",
      Contrast == "IL1" ~ "2. Stimulus\n(IL-1 vs Veh)", 
      Contrast == "Rescue" ~ "3. Rescue\n(T5224+IL1 vs IL-1)"
    )
  )

# 4. Create the Validation Dot Plot with Asterisks
p_feedback_sig <- ggplot(feedback_df, aes(x = Contrast, y = gene_name)) +
  # The coloured dots
  geom_point(aes(fill = logFC), shape = 21, size = 9, colour = "grey30", stroke = 1) +
  
  # The significance asterisks overlaid directly on the dots
  # vjust = 0.75 pushes the asterisk slightly down so it visually centres perfectly inside the dot
  geom_text(aes(label = Significance), vjust = 0.75, size = 6, colour = "black", fontface = "bold") +
  
  scale_fill_gradient2(
    low = "dodgerblue3", mid = "white", high = "firebrick3", midpoint = 0,
    name = "Log2 Fold Change\n(Red = Up, Blue = Down)"
  ) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(face = "bold", size = 11),
    axis.text.y = element_text(face = "italic", size = 11),
    strip.text.y = element_text(face = "bold", angle = 0),
    strip.background = element_rect(fill = "grey90"),
    plot.title = element_text(face = "bold", size = 15),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Validating the Compensatory Feedback Loop",
    subtitle = "Asterisks denote statistical significance (*p<0.05, **p<0.01, ***p<0.001)",
    x = "", y = ""
  )

ggsave("Feedback_Loop_Validation_With_Stats.png", p_feedback_sig, width = 11, height = 8, dpi = 300)
print(p_feedback_sig)

library(tidyverse)

# 1. Define the exploration targets
exploration_genes <- c(
  # Substitutes & Repressors
  "Junb", "Jund", "Fosb", "Fosl1", "Fosl2", "Atf2", "Atf3", "Atf4", "Batf", "Jdp2",
  # Parallel Regulators
  "Rela", "Nfkb1", "Nfkb2", "Relb", "Cebpb", "Cebpd", "Stat3", "Il6st",
  # Core Kinases
  "Map3k1", "Map3k5", "Map3k7", "Map2k4", "Map2k7",
  # Autocrine Cytokines
  "Il6", "Lif", "Cxcl1", "Cxcl2", "Ccl2", "Tnf"
)

# 2. Merge your existing dataframes
# Assuming res_t5224, res_il1, and res_rescue are still in your environment
exploration_df <- res_t5224 %>%
  dplyr::select(gene_name, logFC_T5224 = logFC, FDR_T5224 = FDR) %>%
  inner_join(
    res_il1 %>% dplyr::select(gene_name, logFC_IL1 = logFC, FDR_IL1 = FDR),
    by = "gene_name"
  ) %>%
  inner_join(
    res_rescue %>% dplyr::select(gene_name, logFC_Rescue = logFC, FDR_Rescue = FDR),
    by = "gene_name"
  ) %>%
  filter(gene_name %in% exploration_genes) %>%
  # Arrange by the baseline T5224 response to find the strongest drivers
  arrange(desc(logFC_T5224))

# 3. View the table in RStudio's viewer for easy scrolling
View(exploration_df)

# 4. Save to CSV for a closer look in Excel if needed
write_csv(exploration_df, "Extended_Pathway_Exploration.csv")

suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
})

base_dir <- publication_path("PAPER_OUTPUT_ROOT", "IL1_RNASeq")

out_dir <- file.path(base_dir, "peak_gene_DGE_enrichment")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fdr_cutoff <- 0.05
abs_log2fc_cutoff <- 0

dge_files <- list.files(base_dir, pattern = "^DGE_.*\\.csv$", full.names = TRUE)

bed_files <- list.files(
  base_dir,
  pattern = "(age-gained|background)\\.bed$",
  full.names = TRUE
)

bed_tbl <- tibble(file = bed_files) %>%
  mutate(
    bed_file = basename(file),
    peak_class = case_when(
      str_detect(bed_file, "age-gained") ~ "age_gained",
      str_detect(bed_file, "background") ~ "background",
      TRUE ~ NA_character_
    ),
    peak_group = bed_file %>%
      str_remove("_age-gained\\.bed$") %>%
      str_remove("_background\\.bed$")
  ) %>%
  filter(!is.na(peak_class))

bed_pairs <- bed_tbl %>%
  dplyr::select(peak_group, peak_class, file) %>%
  pivot_wider(names_from = peak_class, values_from = file) %>%
  filter(!is.na(age_gained), !is.na(background))

stopifnot(nrow(bed_pairs) > 0)
stopifnot(length(dge_files) > 0)

genes_gr <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)

genes_gr <- keepStandardChromosomes(genes_gr, pruning.mode = "coarse")
seqlevelsStyle(genes_gr) <- "UCSC"

genes_gr$ENTREZID <- names(genes_gr)

genes_gr$SYMBOL <- mapIds(
  org.Mm.eg.db,
  keys = genes_gr$ENTREZID,
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

genes_gr$ENSEMBL <- mapIds(
  org.Mm.eg.db,
  keys = genes_gr$ENTREZID,
  keytype = "ENTREZID",
  column = "ENSEMBL",
  multiVals = "first"
)

genes_gr <- genes_gr[!is.na(genes_gr$SYMBOL)]

import_bed_clean <- function(path) {
  gr <- rtracklayer::import(path)

  if (!inherits(gr, "GRanges")) {
    stop("Could not import as GRanges: ", path)
  }

  seqlevelsStyle(gr) <- "UCSC"
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  gr <- gr[seqnames(gr) %in% seqlevels(genes_gr)]

  gr
}

nearest_gene_symbols <- function(peaks, genes) {
  idx <- nearest(peaks, genes, ignore.strand = TRUE)
  unique(genes$SYMBOL[idx[!is.na(idx)]])
}

within_100kb_gene_symbols <- function(peaks, genes, flank_bp = 100000) {
  peaks_100kb <- GRanges(
    seqnames = seqnames(peaks),
    ranges = IRanges(
      start = pmax(1, start(peaks) - flank_bp),
      end = end(peaks) + flank_bp
    ),
    strand = "*"
  )

  seqlevelsStyle(peaks_100kb) <- "UCSC"

  hits <- findOverlaps(peaks_100kb, genes, ignore.strand = TRUE)

  unique(genes$SYMBOL[subjectHits(hits)])
}

guess_col <- function(df, candidates, fallback_regex = NULL) {
  nms <- names(df)
  nms_clean <- tolower(nms)

  candidates_clean <- tolower(candidates)

  exact_hit <- which(nms_clean %in% candidates_clean)
  if (length(exact_hit) > 0) {
    return(nms[exact_hit[1]])
  }

  if (!is.null(fallback_regex)) {
    regex_hit <- which(str_detect(nms_clean, fallback_regex))
    if (length(regex_hit) > 0) {
      return(nms[regex_hit[1]])
    }
  }

  NA_character_
}

map_gene_ids_to_symbol <- function(x) {
  x <- as.character(x)
  x <- str_replace(x, "\\.\\d+$", "")

  prop_ensmusg <- mean(str_detect(x, "^ENSMUSG"), na.rm = TRUE)
  prop_numeric <- mean(str_detect(x, "^[0-9]+$"), na.rm = TRUE)

  if (is.nan(prop_ensmusg)) prop_ensmusg <- 0
  if (is.nan(prop_numeric)) prop_numeric <- 0

  if (prop_ensmusg > 0.5) {
    mapped <- mapIds(
      org.Mm.eg.db,
      keys = x,
      keytype = "ENSEMBL",
      column = "SYMBOL",
      multiVals = "first"
    )
    return(unname(mapped))
  }

  if (prop_numeric > 0.5) {
    mapped <- mapIds(
      org.Mm.eg.db,
      keys = x,
      keytype = "ENTREZID",
      column = "SYMBOL",
      multiVals = "first"
    )
    return(unname(mapped))
  }

  x
}

read_dge_file <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)

  gene_col <- guess_col(
    df,
    candidates = c(
      "gene", "genes", "symbol", "gene_symbol", "gene_name",
      "external_gene_name", "feature", "features", "X", "...1"
    ),
    fallback_regex = "gene|symbol|feature"
  )

  p_col <- guess_col(
    df,
    candidates = c(
      "padj", "p_val_adj", "p_val_adj.", "FDR", "adj.P.Val",
      "qval", "q_value", "qvalue", "fdr_bh"
    ),
    fallback_regex = "padj|adj.*p|p.*adj|fdr|qval|q_value"
  )

  lfc_col <- guess_col(
    df,
    candidates = c(
      "log2FoldChange", "logFC", "avg_log2FC", "avg_logFC",
      "log2FC", "estimate", "coef"
    ),
    fallback_regex = "log2.*fold|logfc|log2fc|avg_log"
  )

  if (is.na(gene_col)) {
    stop("Could not identify gene column in: ", basename(path))
  }

  if (is.na(p_col)) {
    stop("Could not identify adjusted P/FDR column in: ", basename(path))
  }

  if (is.na(lfc_col)) {
    stop("Could not identify logFC/log2FC column in: ", basename(path))
  }

  out <- df %>%
    transmute(
      raw_gene = as.character(.data[[gene_col]]),
      gene = map_gene_ids_to_symbol(raw_gene),
      padj = suppressWarnings(as.numeric(.data[[p_col]])),
      log2FC = suppressWarnings(as.numeric(.data[[lfc_col]]))
    ) %>%
    filter(!is.na(gene), gene != "", !is.na(padj), !is.na(log2FC)) %>%
    dplyr::distinct(gene, .keep_all = TRUE)

  attr(out, "gene_col") <- gene_col
  attr(out, "p_col") <- p_col
  attr(out, "lfc_col") <- lfc_col

  out
}

make_dge_sets <- function(dge_df) {
  sig <- dge_df %>%
    filter(padj <= fdr_cutoff, abs(log2FC) >= abs_log2fc_cutoff)

  list(
    all = sig$gene,
    up = sig %>% filter(log2FC > 0) %>% pull(gene),
    down = sig %>% filter(log2FC < 0) %>% pull(gene)
  )
}

run_fisher <- function(age_genes, background_genes, dge_genes, measured_genes) {
  age_genes <- unique(age_genes)
  background_genes <- unique(background_genes)
  dge_genes <- unique(dge_genes)
  measured_genes <- unique(measured_genes)

  universe <- intersect(unique(c(background_genes, age_genes)), measured_genes)

  in_age <- universe %in% age_genes
  in_dge <- universe %in% dge_genes

  mat <- matrix(
    c(
      sum(in_age & in_dge),
      sum(in_age & !in_dge),
      sum(!in_age & in_dge),
      sum(!in_age & !in_dge)
    ),
    nrow = 2,
    byrow = TRUE
  )

  rownames(mat) <- c("age_gained_linked", "background_only_or_not_age_gained")
  colnames(mat) <- c("DGE", "not_DGE")

  if (sum(mat) == 0 || any(rowSums(mat) == 0) || any(colSums(mat) == 0)) {
    return(tibble(
      universe_n = length(universe),
      age_linked_n = sum(in_age),
      dge_n = sum(in_dge),
      overlap_n = sum(in_age & in_dge),
      expected_overlap_n = NA_real_,
      odds_ratio = NA_real_,
      p_value = NA_real_,
      overlap_genes = paste(sort(intersect(age_genes, dge_genes)), collapse = ";")
    ))
  }

  ft <- fisher.test(mat, alternative = "greater")

  tibble(
    universe_n = length(universe),
    age_linked_n = sum(in_age),
    dge_n = sum(in_dge),
    overlap_n = sum(in_age & in_dge),
    expected_overlap_n = sum(in_age) * sum(in_dge) / length(universe),
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value,
    overlap_genes = paste(sort(intersect(age_genes, dge_genes)), collapse = ";")
  )
}

linked_gene_tbl <- map_dfr(seq_len(nrow(bed_pairs)), function(i) {
  peak_group <- bed_pairs$peak_group[i]

  age_peaks <- import_bed_clean(bed_pairs$age_gained[i])
  bg_peaks <- import_bed_clean(bed_pairs$background[i])

  age_nearest <- nearest_gene_symbols(age_peaks, genes_gr)
  bg_nearest <- nearest_gene_symbols(bg_peaks, genes_gr)

  age_100kb <- within_100kb_gene_symbols(age_peaks, genes_gr, flank_bp = 100000)
  bg_100kb <- within_100kb_gene_symbols(bg_peaks, genes_gr, flank_bp = 100000)

  bind_rows(
    tibble(
      peak_group = peak_group,
      peak_class = "age_gained",
      link_method = "nearest",
      gene = age_nearest
    ),
    tibble(
      peak_group = peak_group,
      peak_class = "background",
      link_method = "nearest",
      gene = bg_nearest
    ),
    tibble(
      peak_group = peak_group,
      peak_class = "age_gained",
      link_method = "within_100kb",
      gene = age_100kb
    ),
    tibble(
      peak_group = peak_group,
      peak_class = "background",
      link_method = "within_100kb",
      gene = bg_100kb
    )
  )
}) %>%
  filter(!is.na(gene), gene != "") %>%
  dplyr::distinct()

readr::write_csv(
  linked_gene_tbl,
  file.path(out_dir, "linked_genes_by_peak_set.csv")
)

dge_tbl <- tibble(dge_file = dge_files) %>%
  mutate(
    dge_name = basename(dge_file) %>%
      str_remove("^DGE_") %>%
      str_remove("\\.csv$")
  )

dge_data <- setNames(
  map(dge_tbl$dge_file, read_dge_file),
  dge_tbl$dge_name
)

enrichment_results <- map_dfr(names(dge_data), function(dge_name) {
  dge_df <- dge_data[[dge_name]]
  dge_sets <- make_dge_sets(dge_df)

  map_dfr(unique(linked_gene_tbl$peak_group), function(pg) {
    map_dfr(unique(linked_gene_tbl$link_method), function(method) {
      age_genes <- linked_gene_tbl %>%
        filter(
          peak_group == pg,
          peak_class == "age_gained",
          link_method == method
        ) %>%
        pull(gene) %>%
        unique()

      background_genes <- linked_gene_tbl %>%
        filter(
          peak_group == pg,
          peak_class == "background",
          link_method == method
        ) %>%
        pull(gene) %>%
        unique()

      map_dfr(names(dge_sets), function(direction) {
        run_fisher(
          age_genes = age_genes,
          background_genes = background_genes,
          dge_genes = dge_sets[[direction]],
          measured_genes = dge_df$gene
        ) %>%
          mutate(
            peak_group = pg,
            link_method = method,
            dge_name = dge_name,
            dge_direction = direction,
            .before = 1
          )
      })
    })
  })
}) %>%
  mutate(
    padj_BH_global = p.adjust(p_value, method = "BH"),
    padj_BH_by_link_method = ave(
      p_value,
      link_method,
      FUN = function(x) p.adjust(x, method = "BH")
    ),
    enrichment_direction = case_when(
      is.na(odds_ratio) ~ NA_character_,
      odds_ratio > 1 ~ "enriched",
      odds_ratio < 1 ~ "depleted",
      TRUE ~ "neutral"
    )
  ) %>%
  arrange(padj_BH_global, p_value)

readr::write_csv(
  enrichment_results,
  file.path(out_dir, "DGE_enrichment_age_gained_vs_background.csv")
)

sig_results <- enrichment_results %>%
  filter(!is.na(padj_BH_global), padj_BH_global < 0.1) %>%
  arrange(padj_BH_global, p_value)

readr::write_csv(
  sig_results,
  file.path(out_dir, "DGE_enrichment_significant_FDR10.csv")
)

plot_df <- enrichment_results %>%
  mutate(
    neglog10_padj = -log10(padj_BH_global),
    label = paste(peak_group, link_method, sep = " | ")
  ) %>%
  filter(!is.na(neglog10_padj), is.finite(neglog10_padj))

p <- ggplot(
  plot_df,
  aes(
    x = dge_name,
    y = label,
    size = overlap_n,
    fill = log2(odds_ratio)
  )
) +
  geom_point(shape = 21, colour = "black") +
  facet_wrap(~ dge_direction, nrow = 1) +
  scale_fill_gradient2(
    low = "#377EB8",
    mid = "white",
    high = "#E41A1C",
    midpoint = 0,
    na.value = "grey80"
  ) +
  labs(
    x = "DGE comparison",
    y = "Peak set | gene-linking method",
    size = "Overlap genes",
    fill = "log2 odds ratio",
    title = "Enrichment of DGE genes near age-gained AP-1/FOS/JUN sites",
    subtitle = "Fisher test: age-gained linked genes vs matched background linked genes"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(out_dir, "DGE_enrichment_dotplot.pdf"),
  plot = p,
  width = 13,
  height = 7
)

message("Done.")
message("Main output:")
message(file.path(out_dir, "DGE_enrichment_age_gained_vs_background.csv"))
message("Significant FDR<0.1 output:")
message(file.path(out_dir, "DGE_enrichment_significant_FDR10.csv"))
message("Linked genes output:")
message(file.path(out_dir, "linked_genes_by_peak_set.csv"))

print(
  enrichment_results %>%
    dplyr::select(
      peak_group,
      link_method,
      dge_name,
      dge_direction,
      universe_n,
      age_linked_n,
      dge_n,
      overlap_n,
      expected_overlap_n,
      odds_ratio,
      p_value,
      padj_BH_global
    ) %>%
    arrange(padj_BH_global, p_value) %>%
    head(30)
)

plot_df <- enrichment_results %>%
  mutate(
    neglog10_padj = -log10(padj_BH_global),
    label = paste(peak_group, link_method, sep = " | "),
    sig_label = case_when(
      padj_BH_global < 0.001 ~ "***",
      padj_BH_global < 0.01  ~ "**",
      padj_BH_global < 0.05  ~ "*",
      padj_BH_global < 0.1   ~ ".",
      TRUE ~ ""
    )
  ) %>%
  filter(!is.na(neglog10_padj), is.finite(neglog10_padj))

p <- ggplot(
  plot_df,
  aes(
    x = dge_name,
    y = label,
    size = overlap_n,
    fill = log2(odds_ratio)
  )
) +
  geom_point(shape = 21, colour = "black") +
  geom_text(
    aes(label = sig_label),
    size = 4,
    colour = "black",
    vjust = 0.4
  ) +
  facet_wrap(~ dge_direction, nrow = 1) +
  scale_fill_gradient2(
    low = "#377EB8",
    mid = "white",
    high = "#E41A1C",
    midpoint = 0,
    na.value = "grey80"
  ) +
  labs(
    x = "DGE comparison",
    y = "Peak set | gene-linking method",
    size = "Overlap genes",
    fill = "log2 odds ratio",
    title = "Enrichment of DGE genes near age-gained AP-1/FOS/JUN sites",
    subtitle = "Fisher test: age-gained linked genes vs matched background linked genes; * FDR<0.05, . FDR<0.1"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(out_dir, "DGE_enrichment_dotplot_FDR_annotated.pdf"),
  plot = p,
  width = 13,
  height = 7
)

select      <- dplyr::select
filter      <- dplyr::filter
mutate      <- dplyr::mutate
transmute   <- dplyr::transmute
summarise   <- dplyr::summarise
arrange     <- dplyr::arrange
distinct    <- dplyr::distinct
count       <- dplyr::count
rename      <- dplyr::rename
left_join   <- dplyr::left_join
inner_join  <- dplyr::inner_join
full_join   <- dplyr::full_join
anti_join   <- dplyr::anti_join
group_by    <- dplyr::group_by
ungroup     <- dplyr::ungroup
rowwise     <- dplyr::rowwise
case_when   <- dplyr::case_when
pull        <- dplyr::pull
slice_head  <- dplyr::slice_head
slice_max   <- dplyr::slice_max

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(forcats)
})

# -----------------------------
# User settings
# -----------------------------

base_dir <- publication_path("PAPER_OUTPUT_ROOT", "IL1_RNASeq")

out_dir <- file.path(base_dir, "DGE_overlap_AP1_dependency")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fdr_cutoff <- 0.05
abs_log2fc_cutoff <- 0
top_label_n <- 20
top_heatmap_n <- 60

# -----------------------------
# Input files
# -----------------------------

dge_paths <- list.files(base_dir, pattern = "^DGE_.*\\.csv$", full.names = TRUE)

if (length(dge_paths) == 0) {
  stop("No DGE_*.csv files found in: ", base_dir)
}

# Optional linked-gene file from the previous BED enrichment script
linked_gene_file <- file.path(
  base_dir,
  "peak_gene_DGE_enrichment",
  "linked_genes_by_peak_set.csv"
)

has_linked_gene_file <- file.exists(linked_gene_file)

# -----------------------------
# Helper functions
# -----------------------------

guess_col <- function(df, candidates, fallback_regex = NULL) {
  nms <- names(df)
  nms_clean <- tolower(nms)
  candidates_clean <- tolower(candidates)

  exact_hit <- which(nms_clean %in% candidates_clean)
  if (length(exact_hit) > 0) {
    return(nms[exact_hit[1]])
  }

  if (!is.null(fallback_regex)) {
    regex_hit <- which(stringr::str_detect(nms_clean, fallback_regex))
    if (length(regex_hit) > 0) {
      return(nms[regex_hit[1]])
    }
  }

  NA_character_
}

read_dge_file <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)

  gene_col <- guess_col(
    df,
    candidates = c(
      "gene", "genes", "symbol", "gene_symbol", "gene_name",
      "external_gene_name", "feature", "features", "X", "...1"
    ),
    fallback_regex = "gene|symbol|feature"
  )

  p_col <- guess_col(
    df,
    candidates = c(
      "padj", "p_val_adj", "p_val_adj.", "FDR", "adj.P.Val",
      "qval", "q_value", "qvalue", "fdr_bh"
    ),
    fallback_regex = "padj|adj.*p|p.*adj|fdr|qval|q_value"
  )

  lfc_col <- guess_col(
    df,
    candidates = c(
      "log2FoldChange", "logFC", "avg_log2FC", "avg_logFC",
      "log2FC", "estimate", "coef"
    ),
    fallback_regex = "log2.*fold|logfc|log2fc|avg_log"
  )

  if (is.na(gene_col)) {
    stop("Could not identify gene column in: ", basename(path))
  }

  if (is.na(p_col)) {
    stop("Could not identify adjusted P/FDR column in: ", basename(path))
  }

  if (is.na(lfc_col)) {
    stop("Could not identify logFC/log2FC column in: ", basename(path))
  }

  contrast <- basename(path) |>
    stringr::str_remove("^DGE_") |>
    stringr::str_remove("\\.csv$")

  df |>
    dplyr::transmute(
      gene = as.character(.data[[gene_col]]),
      log2FC = suppressWarnings(as.numeric(.data[[lfc_col]])),
      padj = suppressWarnings(as.numeric(.data[[p_col]])),
      contrast = contrast,
      gene_col_used = gene_col,
      lfc_col_used = lfc_col,
      p_col_used = p_col
    ) |>
    dplyr::mutate(
      gene = stringr::str_replace(gene, "\\.\\d+$", ""),
      sig = !is.na(padj) &
        padj <= fdr_cutoff &
        !is.na(log2FC) &
        abs(log2FC) >= abs_log2fc_cutoff,
      direction = dplyr::case_when(
        sig & log2FC > 0 ~ "up",
        sig & log2FC < 0 ~ "down",
        TRUE ~ "ns"
      )
    ) |>
    dplyr::filter(!is.na(gene), gene != "", !is.na(log2FC), !is.na(padj)) |>
    dplyr::distinct(gene, contrast, .keep_all = TRUE)
}

make_intersection_table <- function(binary_membership_tbl, gene_col = "gene") {
  set_cols <- setdiff(names(binary_membership_tbl), gene_col)

  if (length(set_cols) == 0) {
    return(tibble::tibble(
      intersection = character(),
      n_sets = integer(),
      n_genes = integer()
    ))
  }

  binary_membership_tbl |>
    dplyr::mutate(
      n_sets = rowSums(dplyr::across(dplyr::all_of(set_cols))),
      intersection = apply(
        as.data.frame(dplyr::across(dplyr::all_of(set_cols))),
        1,
        function(x) {
          paste(set_cols[as.logical(x)], collapse = " & ")
        }
      )
    ) |>
    dplyr::filter(n_sets > 0) |>
    dplyr::count(intersection, n_sets, name = "n_genes") |>
    dplyr::arrange(dplyr::desc(n_sets), dplyr::desc(n_genes))
}

safe_first <- function(x) {
  if (length(x) == 0) {
    return(NA_real_)
  }
  x[1]
}

known_il1_targets <- c(
  "Il6", "Il1a", "Il1b", "Il11", "Lif", "Tnf", "Csf2", "Csf3",
  "Cxcl1", "Cxcl2", "Cxcl3", "Cxcl5", "Cxcl10",
  "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Ccl7", "Ccl20",
  "Ptgs2", "Ptges", "Ptges2", "Pla2g4a",
  "Nfkbia", "Nfkbiz", "Nfkb2", "Relb",
  "Tnfaip3", "Tnfaip2", "Zc3h12a", "Ier3",
  "Socs3", "Map3k8",
  "Fos", "Fosb", "Jun", "Junb", "Jund", "Egr1", "Egr2",
  "Dusp1", "Dusp2", "Dusp5",
  "Icam1", "Vcam1", "Sele", "Selp",
  "Mmp3", "Mmp9", "Mmp13", "Mmp14",
  "Timp1", "Has1", "Has2",
  "Nos2", "Nlrp3", "Saa3"
)

# -----------------------------
# Read DGE files
# -----------------------------

dge_long <- purrr::map_dfr(dge_paths, read_dge_file)

readr::write_csv(
  dge_long,
  file.path(out_dir, "all_DGE_long_format.csv")
)

contrast_check <- dge_long |>
  dplyr::distinct(contrast, gene_col_used, lfc_col_used, p_col_used) |>
  dplyr::arrange(contrast)

readr::write_csv(
  contrast_check,
  file.path(out_dir, "DGE_column_detection_check.csv")
)

message("Detected DGE columns:")
print(contrast_check)

required_contrasts <- c(
  "IL1_vs_Veh_3h",
  "IL1_vs_Veh_6h",
  "T5224preIL1_vs_IL1_3h",
  "T5224preIL1_vs_IL1_6h"
)

missing_required <- setdiff(required_contrasts, unique(dge_long$contrast))

if (length(missing_required) > 0) {
  warning("Missing expected contrasts: ", paste(missing_required, collapse = ", "))
}

# -----------------------------
# DEG sets
# -----------------------------

dge_sets <- dge_long |>
  dplyr::filter(sig) |>
  dplyr::mutate(set_name = paste(contrast, direction, sep = "__")) |>
  dplyr::select(gene, set_name) |>
  dplyr::distinct()

set_sizes <- dge_sets |>
  dplyr::count(set_name, name = "n_genes") |>
  dplyr::arrange(dplyr::desc(n_genes))

readr::write_csv(
  set_sizes,
  file.path(out_dir, "DEG_set_sizes.csv")
)

p_set_sizes <- ggplot2::ggplot(
  set_sizes,
  ggplot2::aes(x = reorder(set_name, n_genes), y = n_genes)
) +
  ggplot2::geom_col(fill = "grey40") +
  ggplot2::coord_flip() +
  ggplot2::labs(
    x = NULL,
    y = "Significant genes",
    title = "DGE set sizes",
    subtitle = paste0("FDR <= ", fdr_cutoff, ", abs(log2FC) >= ", abs_log2fc_cutoff)
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

ggplot2::ggsave(
  file.path(out_dir, "DGE_set_sizes.pdf"),
  p_set_sizes,
  width = 8,
  height = 6
)

# -----------------------------
# Focused overlap between IL1 and T5224preIL1_vs_IL1 DEG sets
# -----------------------------

focus_sets <- c(
  "IL1_vs_Veh_3h__up",
  "IL1_vs_Veh_3h__down",
  "IL1_vs_Veh_6h__up",
  "IL1_vs_Veh_6h__down",
  "T5224preIL1_vs_IL1_3h__up",
  "T5224preIL1_vs_IL1_3h__down",
  "T5224preIL1_vs_IL1_6h__up",
  "T5224preIL1_vs_IL1_6h__down"
)

focus_sets <- intersect(focus_sets, unique(dge_sets$set_name))

membership <- dge_sets |>
  dplyr::filter(set_name %in% focus_sets) |>
  dplyr::mutate(value = TRUE) |>
  tidyr::pivot_wider(
    names_from = set_name,
    values_from = value,
    values_fill = FALSE
  )

intersection_tbl <- make_intersection_table(membership, gene_col = "gene")

readr::write_csv(
  intersection_tbl,
  file.path(out_dir, "focused_DEG_set_intersections.csv")
)

p_intersections <- intersection_tbl |>
  dplyr::filter(n_sets >= 2) |>
  dplyr::slice_max(n_genes, n = 30) |>
  dplyr::mutate(intersection = forcats::fct_reorder(intersection, n_genes)) |>
  ggplot2::ggplot(ggplot2::aes(x = intersection, y = n_genes, fill = n_sets)) +
  ggplot2::geom_col(colour = "black", linewidth = 0.2) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_gradient(low = "grey80", high = "grey20") +
  ggplot2::labs(
    x = NULL,
    y = "Genes in intersection",
    fill = "Number of sets",
    title = "Major overlaps among IL-1 and T5224 DGE sets",
    subtitle = "Focused on IL1_vs_Veh and T5224preIL1_vs_IL1 contrasts"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

ggplot2::ggsave(
  file.path(out_dir, "focused_DEG_set_intersections.pdf"),
  p_intersections,
  width = 11,
  height = 8
)

# -----------------------------
# Pairwise Jaccard overlap between all DEG sets
# -----------------------------

set_membership_binary <- dge_sets |>
  dplyr::mutate(value = 1) |>
  tidyr::pivot_wider(
    names_from = set_name,
    values_from = value,
    values_fill = 0
  )

all_set_names <- setdiff(names(set_membership_binary), "gene")

pairwise_jaccard <- tidyr::expand_grid(
  set_a = all_set_names,
  set_b = all_set_names
) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    n_a = sum(set_membership_binary[[set_a]] == 1),
    n_b = sum(set_membership_binary[[set_b]] == 1),
    overlap = sum(set_membership_binary[[set_a]] == 1 & set_membership_binary[[set_b]] == 1),
    union = sum(set_membership_binary[[set_a]] == 1 | set_membership_binary[[set_b]] == 1),
    jaccard = ifelse(union > 0, overlap / union, NA_real_)
  ) |>
  dplyr::ungroup()

readr::write_csv(
  pairwise_jaccard,
  file.path(out_dir, "pairwise_DEG_set_jaccard.csv")
)

p_jaccard <- ggplot2::ggplot(
  pairwise_jaccard,
  ggplot2::aes(x = set_a, y = set_b, fill = jaccard)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
  ggplot2::scale_fill_gradient(low = "white", high = "black", na.value = "grey90") +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    fill = "Jaccard",
    title = "Pairwise overlap between DGE sets"
  ) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    panel.grid = ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(out_dir, "pairwise_DEG_set_jaccard_heatmap.pdf"),
  p_jaccard,
  width = 10,
  height = 9
)

# -----------------------------
# AP-1 dependency classification
# -----------------------------

dge_wide <- dge_long |>
  dplyr::select(gene, contrast, log2FC, padj, sig, direction) |>
  tidyr::pivot_wider(
    names_from = contrast,
    values_from = c(log2FC, padj, sig, direction),
    names_sep = "__"
  )

has_col <- function(x) {
  x %in% names(dge_wide)
}

classify_timepoint <- function(df, timepoint) {
  il1_lfc <- paste0("log2FC__IL1_vs_Veh_", timepoint)
  il1_padj <- paste0("padj__IL1_vs_Veh_", timepoint)
  il1_sig <- paste0("sig__IL1_vs_Veh_", timepoint)

  t5224_il1_lfc <- paste0("log2FC__T5224preIL1_vs_IL1_", timepoint)
  t5224_il1_padj <- paste0("padj__T5224preIL1_vs_IL1_", timepoint)
  t5224_il1_sig <- paste0("sig__T5224preIL1_vs_IL1_", timepoint)

  t5224_alone_lfc <- paste0("log2FC__T5224_vs_Veh_", timepoint)
  t5224_alone_padj <- paste0("padj__T5224_vs_Veh_", timepoint)
  t5224_alone_sig <- paste0("sig__T5224_vs_Veh_", timepoint)

  needed <- c(
    il1_lfc,
    il1_padj,
    il1_sig,
    t5224_il1_lfc,
    t5224_il1_padj,
    t5224_il1_sig
  )

  if (!all(needed %in% names(df))) {
    warning(
      "Skipping ",
      timepoint,
      ". Missing: ",
      paste(setdiff(needed, names(df)), collapse = ", ")
    )
    return(tibble::tibble())
  }

  df |>
    dplyr::transmute(
      gene = gene,
      timepoint = timepoint,
      IL1_log2FC = .data[[il1_lfc]],
      IL1_padj = .data[[il1_padj]],
      IL1_sig = .data[[il1_sig]],
      T5224preIL1_vs_IL1_log2FC = .data[[t5224_il1_lfc]],
      T5224preIL1_vs_IL1_padj = .data[[t5224_il1_padj]],
      T5224preIL1_vs_IL1_sig = .data[[t5224_il1_sig]],
      T5224_vs_Veh_log2FC = if (has_col(t5224_alone_lfc)) .data[[t5224_alone_lfc]] else NA_real_,
      T5224_vs_Veh_padj = if (has_col(t5224_alone_padj)) .data[[t5224_alone_padj]] else NA_real_,
      T5224_vs_Veh_sig = if (has_col(t5224_alone_sig)) .data[[t5224_alone_sig]] else NA
    ) |>
    dplyr::mutate(
      IL1_direction = dplyr::case_when(
        IL1_sig & IL1_log2FC > 0 ~ "IL1_up",
        IL1_sig & IL1_log2FC < 0 ~ "IL1_down",
        TRUE ~ "IL1_ns"
      ),
      T5224_effect_under_IL1 = dplyr::case_when(
        T5224preIL1_vs_IL1_sig & T5224preIL1_vs_IL1_log2FC > 0 ~ "T5224_increases_vs_IL1",
        T5224preIL1_vs_IL1_sig & T5224preIL1_vs_IL1_log2FC < 0 ~ "T5224_decreases_vs_IL1",
        TRUE ~ "T5224_no_effect_vs_IL1"
      ),
      AP1_dependency_class = dplyr::case_when(
        IL1_sig & IL1_log2FC > 0 &
          T5224preIL1_vs_IL1_sig & T5224preIL1_vs_IL1_log2FC < 0 ~
          "IL1_induced_AP1_dependent",

        IL1_sig & IL1_log2FC < 0 &
          T5224preIL1_vs_IL1_sig & T5224preIL1_vs_IL1_log2FC > 0 ~
          "IL1_repressed_AP1_dependent",

        IL1_sig & IL1_log2FC > 0 &
          T5224preIL1_vs_IL1_sig & T5224preIL1_vs_IL1_log2FC > 0 ~
          "IL1_induced_T5224_enhanced",

        IL1_sig & IL1_log2FC < 0 &
          T5224preIL1_vs_IL1_sig & T5224preIL1_vs_IL1_log2FC < 0 ~
          "IL1_repressed_T5224_enhanced",

        IL1_sig & !T5224preIL1_vs_IL1_sig ~
          "IL1_changed_AP1_independent_or_not_detected",

        TRUE ~ "not_IL1_DEG"
      ),
      AP1_dependency_score = dplyr::case_when(
        IL1_log2FC > 0 ~ IL1_log2FC - T5224preIL1_vs_IL1_log2FC,
        IL1_log2FC < 0 ~ abs(IL1_log2FC) + T5224preIL1_vs_IL1_log2FC,
        TRUE ~ NA_real_
      ),
      T5224_alone_changed = !is.na(T5224_vs_Veh_sig) & T5224_vs_Veh_sig,
      strict_IL1_specific_AP1_dependent = AP1_dependency_class %in% c(
        "IL1_induced_AP1_dependent",
        "IL1_repressed_AP1_dependent"
      ) & !T5224_alone_changed
    )
}

ap1_tbl <- dplyr::bind_rows(
  classify_timepoint(dge_wide, "3h"),
  classify_timepoint(dge_wide, "6h")
) |>
  dplyr::arrange(
    timepoint,
    dplyr::desc(AP1_dependency_class == "IL1_induced_AP1_dependent"),
    dplyr::desc(AP1_dependency_score)
  )

readr::write_csv(
  ap1_tbl,
  file.path(out_dir, "AP1_dependency_all_genes_by_timepoint.csv")
)

ap1_dependent_genes <- ap1_tbl |>
  dplyr::filter(AP1_dependency_class %in% c(
    "IL1_induced_AP1_dependent",
    "IL1_repressed_AP1_dependent"
  )) |>
  dplyr::arrange(timepoint, AP1_dependency_class, dplyr::desc(AP1_dependency_score))

readr::write_csv(
  ap1_dependent_genes,
  file.path(out_dir, "AP1_dependent_IL1_genes.csv")
)

ap1_induced_dependent_genes <- ap1_tbl |>
  dplyr::filter(AP1_dependency_class == "IL1_induced_AP1_dependent") |>
  dplyr::arrange(timepoint, dplyr::desc(AP1_dependency_score))

readr::write_csv(
  ap1_induced_dependent_genes,
  file.path(out_dir, "IL1_induced_AP1_dependent_genes.csv")
)

ap1_strict <- ap1_dependent_genes |>
  dplyr::filter(strict_IL1_specific_AP1_dependent)

readr::write_csv(
  ap1_strict,
  file.path(out_dir, "strict_IL1_specific_AP1_dependent_genes_excluding_T5224_alone_DEGs.csv")
)

class_counts <- ap1_tbl |>
  dplyr::count(timepoint, AP1_dependency_class, name = "n_genes") |>
  dplyr::arrange(timepoint, dplyr::desc(n_genes))

readr::write_csv(
  class_counts,
  file.path(out_dir, "AP1_dependency_class_counts.csv")
)

p_class_counts <- ggplot2::ggplot(
  class_counts,
  ggplot2::aes(x = AP1_dependency_class, y = n_genes, fill = timepoint)
) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), colour = "black", linewidth = 0.2) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    x = NULL,
    y = "Genes",
    fill = "Timepoint",
    title = "Classification of IL-1 and AP-1 inhibitor responses",
    subtitle = "AP-1-dependent induced genes are IL1 up and T5224preIL1_vs_IL1 down"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

ggplot2::ggsave(
  file.path(out_dir, "AP1_dependency_class_counts.pdf"),
  p_class_counts,
  width = 9,
  height = 6
)

# -----------------------------
# Scatter: IL1 response vs T5224 effect under IL1
# -----------------------------

plot_scatter_df <- ap1_tbl |>
  dplyr::filter(!is.na(IL1_log2FC), !is.na(T5224preIL1_vs_IL1_log2FC)) |>
  dplyr::mutate(
    plot_class = dplyr::case_when(
      AP1_dependency_class == "IL1_induced_AP1_dependent" ~ "IL1 induced, AP1 dependent",
      AP1_dependency_class == "IL1_repressed_AP1_dependent" ~ "IL1 repressed, AP1 dependent",
      AP1_dependency_class == "IL1_induced_T5224_enhanced" ~ "IL1 induced, T5224 enhanced",
      AP1_dependency_class == "IL1_repressed_T5224_enhanced" ~ "IL1 repressed, T5224 enhanced",
      AP1_dependency_class == "IL1_changed_AP1_independent_or_not_detected" ~ "IL1 DEG, no AP1 dependence",
      TRUE ~ "Other"
    ),
    label_gene = AP1_dependency_class == "IL1_induced_AP1_dependent"
  )

label_df <- plot_scatter_df |>
  dplyr::filter(label_gene) |>
  dplyr::group_by(timepoint) |>
  dplyr::slice_max(AP1_dependency_score, n = top_label_n, with_ties = FALSE) |>
  dplyr::ungroup()

p_scatter <- ggplot2::ggplot(
  plot_scatter_df,
  ggplot2::aes(x = IL1_log2FC, y = T5224preIL1_vs_IL1_log2FC)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_point(ggplot2::aes(colour = plot_class), alpha = 0.65, size = 1.7) +
  ggplot2::geom_text(
    data = label_df,
    ggplot2::aes(label = gene),
    size = 2.5,
    check_overlap = TRUE,
    vjust = -0.5
  ) +
  ggplot2::facet_wrap(~ timepoint) +
  ggplot2::scale_colour_manual(
    values = c(
      "IL1 induced, AP1 dependent" = "#A81B1B",
      "IL1 repressed, AP1 dependent" = "#101369",
      "IL1 induced, T5224 enhanced" = "#F99349",
      "IL1 repressed, T5224 enhanced" = "#3A85A8",
      "IL1 DEG, no AP1 dependence" = "grey45",
      "Other" = "grey85"
    )
  ) +
  ggplot2::labs(
    x = "IL1 vs Vehicle log2FC",
    y = "T5224preIL1 vs IL1 log2FC",
    colour = NULL,
    title = "IL-1 response versus AP-1 inhibitor effect",
    subtitle = "Key quadrant: IL1 log2FC > 0 and T5224preIL1_vs_IL1 log2FC < 0"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  file.path(out_dir, "IL1_vs_T5224_AP1_dependency_scatter.pdf"),
  p_scatter,
  width = 10,
  height = 7
)

# -----------------------------
# AP1-dependent gene overlap across 3h and 6h
# -----------------------------

ap1_membership_summary <- ap1_tbl |>
  dplyr::filter(AP1_dependency_class %in% c(
    "IL1_induced_AP1_dependent",
    "IL1_repressed_AP1_dependent"
  )) |>
  dplyr::select(gene, timepoint, AP1_dependency_class, AP1_dependency_score) |>
  dplyr::mutate(value = TRUE) |>
  tidyr::unite("set", timepoint, AP1_dependency_class, sep = "__") |>
  dplyr::select(gene, set, value) |>
  dplyr::distinct() |>
  tidyr::pivot_wider(
    names_from = set,
    values_from = value,
    values_fill = FALSE
  )

ap1_intersections <- make_intersection_table(ap1_membership_summary, gene_col = "gene")

readr::write_csv(
  ap1_intersections,
  file.path(out_dir, "AP1_dependent_gene_intersections_3h_6h.csv")
)

p_ap1_intersections <- ap1_intersections |>
  dplyr::mutate(intersection = forcats::fct_reorder(intersection, n_genes)) |>
  ggplot2::ggplot(ggplot2::aes(x = intersection, y = n_genes, fill = n_sets)) +
  ggplot2::geom_col(colour = "black", linewidth = 0.2) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_gradient(low = "grey75", high = "grey20") +
  ggplot2::labs(
    x = NULL,
    y = "Genes",
    fill = "Number of AP1-dependent sets",
    title = "Overlap of AP-1-dependent IL-1 gene classes across 3 h and 6 h"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

ggplot2::ggsave(
  file.path(out_dir, "AP1_dependent_gene_intersections_3h_6h.pdf"),
  p_ap1_intersections,
  width = 10,
  height = 5
)

# -----------------------------
# Heatmap of top IL1-induced AP1-dependent genes
# -----------------------------

top_genes_for_heatmap <- ap1_tbl |>
  dplyr::filter(AP1_dependency_class == "IL1_induced_AP1_dependent") |>
  dplyr::group_by(gene) |>
  dplyr::summarise(
    max_dependency_score = max(AP1_dependency_score, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(max_dependency_score)) |>
  dplyr::slice_head(n = top_heatmap_n) |>
  dplyr::pull(gene)

if (length(top_genes_for_heatmap) > 0) {
  heatmap_df <- dge_long |>
    dplyr::filter(gene %in% top_genes_for_heatmap) |>
    dplyr::mutate(
      contrast = factor(
        contrast,
        levels = c(
          "IL1_vs_Veh_3h",
          "T5224preIL1_vs_IL1_3h",
          "T5224preIL1_vs_Veh_3h",
          "T5224_vs_Veh_3h",
          "IL1_vs_Veh_6h",
          "T5224preIL1_vs_IL1_6h",
          "T5224preIL1_vs_Veh_6h",
          "T5224_vs_Veh_6h"
        )
      ),
      gene = factor(gene, levels = rev(top_genes_for_heatmap)),
      sig_label = ifelse(sig, "*", "")
    ) |>
    dplyr::filter(!is.na(contrast))

  p_heatmap <- ggplot2::ggplot(
    heatmap_df,
    ggplot2::aes(x = contrast, y = gene, fill = log2FC)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = sig_label), size = 2.5) +
    ggplot2::scale_fill_gradient2(
      low = "#377EB8",
      mid = "white",
      high = "#E41A1C",
      midpoint = 0,
      na.value = "grey90"
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      fill = "log2FC",
      title = "Top IL-1-induced AP-1-dependent genes",
      subtitle = "* indicates significant in that contrast"
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )

  ggplot2::ggsave(
    file.path(out_dir, "top_IL1_induced_AP1_dependent_genes_heatmap.pdf"),
    p_heatmap,
    width = 9,
    height = max(5, length(top_genes_for_heatmap) * 0.16)
  )
}

# -----------------------------
# Known IL1 target vector overlay
# -----------------------------

known_il1_ap1_hits <- ap1_tbl |>
  dplyr::filter(gene %in% known_il1_targets) |>
  dplyr::arrange(
    timepoint,
    dplyr::desc(AP1_dependency_class == "IL1_induced_AP1_dependent"),
    dplyr::desc(AP1_dependency_score)
  )

readr::write_csv(
  known_il1_ap1_hits,
  file.path(out_dir, "known_IL1_targets_AP1_dependency_status.csv")
)

p_known_scatter <- plot_scatter_df |>
  dplyr::mutate(is_known_il1_target = gene %in% known_il1_targets) |>
  ggplot2::ggplot(ggplot2::aes(x = IL1_log2FC, y = T5224preIL1_vs_IL1_log2FC)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_point(ggplot2::aes(colour = is_known_il1_target), alpha = 0.7, size = 1.8) +
  ggplot2::geom_text(
    data = plot_scatter_df |> dplyr::filter(gene %in% known_il1_targets),
    ggplot2::aes(label = gene),
    size = 2.5,
    check_overlap = TRUE,
    vjust = -0.5
  ) +
  ggplot2::facet_wrap(~ timepoint) +
  ggplot2::scale_colour_manual(
    values = c("TRUE" = "#A81B1B", "FALSE" = "grey80"),
    labels = c("FALSE" = "Other genes", "TRUE" = "Known IL-1 target vector")
  ) +
  ggplot2::labs(
    x = "IL1 vs Vehicle log2FC",
    y = "T5224preIL1 vs IL1 log2FC",
    colour = NULL,
    title = "Known IL-1 target genes in the AP-1 dependency space"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  file.path(out_dir, "known_IL1_targets_AP1_dependency_scatter.pdf"),
  p_known_scatter,
  width = 10,
  height = 7
)

# -----------------------------
# Ranked summary table
# -----------------------------

summary_table <- ap1_tbl |>
  dplyr::filter(AP1_dependency_class == "IL1_induced_AP1_dependent") |>
  dplyr::group_by(gene) |>
  dplyr::summarise(
    timepoints = paste(sort(unique(timepoint)), collapse = ";"),
    n_timepoints = dplyr::n_distinct(timepoint),
    max_AP1_dependency_score = max(AP1_dependency_score, na.rm = TRUE),

    IL1_log2FC_3h = safe_first(IL1_log2FC[timepoint == "3h"]),
    IL1_padj_3h = safe_first(IL1_padj[timepoint == "3h"]),
    T5224preIL1_vs_IL1_log2FC_3h = safe_first(T5224preIL1_vs_IL1_log2FC[timepoint == "3h"]),
    T5224preIL1_vs_IL1_padj_3h = safe_first(T5224preIL1_vs_IL1_padj[timepoint == "3h"]),

    IL1_log2FC_6h = safe_first(IL1_log2FC[timepoint == "6h"]),
    IL1_padj_6h = safe_first(IL1_padj[timepoint == "6h"]),
    T5224preIL1_vs_IL1_log2FC_6h = safe_first(T5224preIL1_vs_IL1_log2FC[timepoint == "6h"]),
    T5224preIL1_vs_IL1_padj_6h = safe_first(T5224preIL1_vs_IL1_padj[timepoint == "6h"]),

    known_IL1_target = gene %in% known_il1_targets,
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(n_timepoints), dplyr::desc(max_AP1_dependency_score))

readr::write_csv(
  summary_table,
  file.path(out_dir, "ranked_IL1_induced_AP1_dependent_gene_summary.csv")
)

# -----------------------------
# Optional: annotate AP1-dependent genes with age-gained/background BED-linked genes
# -----------------------------

if (has_linked_gene_file) {
  linked_gene_tbl <- readr::read_csv(linked_gene_file, show_col_types = FALSE)

  linked_gene_summary <- linked_gene_tbl |>
    dplyr::mutate(
      linked_label = paste(peak_group, peak_class, link_method, sep = "__")
    ) |>
    dplyr::select(gene, linked_label) |>
    dplyr::distinct() |>
    dplyr::mutate(value = TRUE) |>
    tidyr::pivot_wider(
      names_from = linked_label,
      values_from = value,
      values_fill = FALSE
    )

  summary_table_with_links <- summary_table |>
    dplyr::left_join(linked_gene_summary, by = "gene")

  readr::write_csv(
    summary_table_with_links,
    file.path(out_dir, "ranked_IL1_induced_AP1_dependent_gene_summary_with_peak_links.csv")
  )

  ap1_tbl_with_links <- ap1_tbl |>
    dplyr::left_join(linked_gene_summary, by = "gene")

  readr::write_csv(
    ap1_tbl_with_links,
    file.path(out_dir, "AP1_dependency_all_genes_by_timepoint_with_peak_links.csv")
  )

  link_cols <- setdiff(names(linked_gene_summary), "gene")

  if (length(link_cols) > 0) {
    link_count_tbl <- summary_table_with_links |>
      dplyr::mutate(
        n_peak_link_classes = rowSums(
          dplyr::across(dplyr::all_of(link_cols), ~ tidyr::replace_na(.x, FALSE))
        )
      ) |>
      dplyr::arrange(dplyr::desc(n_peak_link_classes), dplyr::desc(max_AP1_dependency_score))

    readr::write_csv(
      link_count_tbl,
      file.path(out_dir, "ranked_IL1_induced_AP1_dependent_gene_summary_peak_link_counts.csv")
    )
  }
}

# -----------------------------
# Print summary
# -----------------------------

message("Done.")
message("Output directory:")
message(out_dir)

message("\nMost useful files:")
message(file.path(out_dir, "ranked_IL1_induced_AP1_dependent_gene_summary.csv"))
message(file.path(out_dir, "IL1_induced_AP1_dependent_genes.csv"))
message(file.path(out_dir, "AP1_dependent_IL1_genes.csv"))
message(file.path(out_dir, "strict_IL1_specific_AP1_dependent_genes_excluding_T5224_alone_DEGs.csv"))
message(file.path(out_dir, "focused_DEG_set_intersections.csv"))

if (has_linked_gene_file) {
  message(file.path(out_dir, "ranked_IL1_induced_AP1_dependent_gene_summary_with_peak_links.csv"))
}

message("\nMost useful plots:")
message(file.path(out_dir, "IL1_vs_T5224_AP1_dependency_scatter.pdf"))
message(file.path(out_dir, "AP1_dependency_class_counts.pdf"))
message(file.path(out_dir, "focused_DEG_set_intersections.pdf"))
message(file.path(out_dir, "pairwise_DEG_set_jaccard_heatmap.pdf"))
message(file.path(out_dir, "known_IL1_targets_AP1_dependency_scatter.pdf"))

if (length(top_genes_for_heatmap) > 0) {
  message(file.path(out_dir, "top_IL1_induced_AP1_dependent_genes_heatmap.pdf"))
}

message("\nTop IL-1-induced AP-1-dependent genes:")

1+1

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(forcats)
})

# ============================================================
# Fibroblast-only overlap between:
#   1. IL-1/AP-1 dependency RNA-seq results
#   2. Fibroblast age/exercise DEG effects
#
# Key AP-1-dependent definition:
#   IL1_vs_Veh log2FC > 0 and significant
#   T5224preIL1_vs_IL1 log2FC < 0 and significant
#
# Key biological overlap:
#   IL-1-induced AP-1-dependent genes
#   that are age-up and exercise-down in fibroblasts
# ============================================================

# -----------------------------
# Paths
# -----------------------------

base_dir <- publication_path("PAPER_OUTPUT_ROOT", "IL1_RNASeq")

ap1_dir <- file.path(base_dir, "DGE_overlap_AP1_dependency")

ap1_file <- file.path(
  ap1_dir,
  "AP1_dependency_all_genes_by_timepoint.csv"
)

age_ex_rds <- file.path(base_dir, "All_DEG_Age_Ex.rds")

# Allows the same script to work if testing locally with uploaded file
if (!file.exists(age_ex_rds) && file.exists("/mnt/data/All_DEG_Age_Ex.rds")) {
  age_ex_rds <- "/mnt/data/All_DEG_Age_Ex.rds"
}

out_dir <- file.path(ap1_dir, "Fibroblast_AgeEx_overlap")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(ap1_file)) {
  stop("Missing AP-1 dependency file. Run the AP-1 dependency script first: ", ap1_file)
}

if (!file.exists(age_ex_rds)) {
  stop("Missing All_DEG_Age_Ex.rds. Expected at: ", age_ex_rds)
}

# -----------------------------
# Thresholds
# -----------------------------

age_ex_fdr_cutoff <- 0.05
age_ex_abs_log2fc_cutoff <- 0

top_label_n <- 25
top_heatmap_n <- 60

# -----------------------------
# Colours
# -----------------------------

RNA_pallate <- c(
  '#A53D99',
  '#F99349',
  '#87CEEB',
  '#FDDA62',
  '#808080',
  '#BEBADA',
  "#A81B1B",
  "#101369",
  "#606692",
  "#3A85A8",
  "#42977E",
  "#591f0a",
  "#3b7818",
  "#B6D66A",
  "#1CE382",
  "#BF862B",
  "#CC6A6F",
  "#E086B5",
  "#5F1296",
  "#5F1296",
  "#B784A7"
)

RNA_order <- c(
  "Macrophage",
  "M1",
  "M2",
  "M_Ex",
  "HSPC",
  "Dendritic",
  "Cardiomyocyte",
  "Coronary_EC",
  "Endocardial_EC",
  "Lymphatic_EC",
  "Fibroblast",
  "Mesothelial",
  "B_cell",
  "T_cell",
  "Pericytes",
  "Smooth_Muscle",
  "Cardiac_Neuronal",
  "Unknown_1",
  "Unknown_3",
  "ATAC_unknown1",
  "NK"
)

names(RNA_pallate) <- RNA_order

status_cols <- c(
  "Age up + Exercise down" = "#A81B1B",
  "Age down + Exercise up" = "#377EB8",
  "Reversed by exercise" = "#A53D99",
  "Age up only" = "#F99349",
  "Exercise down only" = "#42977E",
  "Other fibroblast DEG" = "grey45",
  "Not fibroblast Age/Exercise DEG" = "grey82"
)

ap1_class_cols <- c(
  "IL1_induced_AP1_dependent" = "#A81B1B",
  "IL1_repressed_AP1_dependent" = "#101369",
  "IL1_changed_AP1_independent_or_not_detected" = "grey55",
  "not_IL1_DEG" = "grey85"
)

known_il1_targets <- c(
  "Il6", "Il1a", "Il1b", "Il11", "Lif", "Tnf", "Csf2", "Csf3",
  "Cxcl1", "Cxcl2", "Cxcl3", "Cxcl5", "Cxcl10",
  "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Ccl7", "Ccl20",
  "Ptgs2", "Ptges", "Ptges2", "Pla2g4a",
  "Nfkbia", "Nfkbiz", "Nfkb2", "Relb",
  "Tnfaip3", "Tnfaip2", "Zc3h12a", "Ier3",
  "Socs3", "Map3k8",
  "Fos", "Fosb", "Jun", "Junb", "Jund", "Egr1", "Egr2",
  "Dusp1", "Dusp2", "Dusp5",
  "Icam1", "Vcam1", "Sele", "Selp",
  "Mmp3", "Mmp9", "Mmp13", "Mmp14",
  "Timp1", "Has1", "Has2",
  "Nos2", "Nlrp3", "Saa3"
)

# -----------------------------
# Helper functions
# -----------------------------

guess_col <- function(df, candidates, fallback_regex = NULL, exclude_regex = NULL) {
  nms <- names(df)
  nms_clean <- tolower(nms)
  candidates_clean <- tolower(candidates)

  exact_hit <- which(nms_clean %in% candidates_clean)

  if (!is.null(exclude_regex) && length(exact_hit) > 0) {
    exact_hit <- exact_hit[!stringr::str_detect(nms_clean[exact_hit], exclude_regex)]
  }

  if (length(exact_hit) > 0) {
    return(nms[exact_hit[1]])
  }

  if (!is.null(fallback_regex)) {
    regex_hit <- which(stringr::str_detect(nms_clean, fallback_regex))

    if (!is.null(exclude_regex) && length(regex_hit) > 0) {
      regex_hit <- regex_hit[!stringr::str_detect(nms_clean[regex_hit], exclude_regex)]
    }

    if (length(regex_hit) > 0) {
      return(nms[regex_hit[1]])
    }
  }

  NA_character_
}

flatten_rds_to_tibble <- function(obj) {
  if (is.data.frame(obj)) {
    df <- as.data.frame(obj)

    if (
      !"gene" %in% tolower(names(df)) &&
        !all(rownames(df) == as.character(seq_len(nrow(df))))
    ) {
      df <- tibble::rownames_to_column(df, var = "rowname_gene")
    }

    return(tibble::as_tibble(df))
  }

  if (is.list(obj)) {
    df_list <- purrr::imap(obj, function(x, nm) {
      if (!is.data.frame(x)) {
        return(NULL)
      }

      df <- as.data.frame(x)

      if (
        !"gene" %in% tolower(names(df)) &&
          !all(rownames(df) == as.character(seq_len(nrow(df))))
      ) {
        df <- tibble::rownames_to_column(df, var = "rowname_gene")
      }

      df$rds_list_name <- nm
      tibble::as_tibble(df)
    })

    df_list <- purrr::compact(df_list)

    if (length(df_list) == 0) {
      stop("RDS is a list, but no data.frames were found inside it.")
    }

    return(dplyr::bind_rows(df_list))
  }

  stop("Unsupported RDS object type: ", paste(class(obj), collapse = ", "))
}

safe_first_numeric <- function(x) {
  if (length(x) == 0) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(x[1]))
}

safe_first_character <- function(x) {
  if (length(x) == 0) {
    return(NA_character_)
  }
  as.character(x[1])
}

run_fisher_overlap <- function(query_genes, target_genes, universe_genes) {
  query_genes <- unique(query_genes)
  target_genes <- unique(target_genes)
  universe_genes <- unique(universe_genes)

  query_genes <- intersect(query_genes, universe_genes)
  target_genes <- intersect(target_genes, universe_genes)

  in_query <- universe_genes %in% query_genes
  in_target <- universe_genes %in% target_genes

  overlap_genes <- sort(intersect(query_genes, target_genes))

  mat <- matrix(
    c(
      sum(in_query & in_target),
      sum(in_query & !in_target),
      sum(!in_query & in_target),
      sum(!in_query & !in_target)
    ),
    nrow = 2,
    byrow = TRUE
  )

  rownames(mat) <- c("IL1_AP1_gene_set", "background_genes")
  colnames(mat) <- c("fibroblast_DEG_set", "not_fibroblast_DEG_set")

  if (
    length(universe_genes) == 0 ||
      any(rowSums(mat) == 0) ||
      any(colSums(mat) == 0)
  ) {
    return(tibble::tibble(
      universe_n = length(universe_genes),
      query_n = length(query_genes),
      target_n = length(target_genes),
      overlap_n = length(overlap_genes),
      expected_overlap_n = NA_real_,
      odds_ratio = NA_real_,
      p_value = NA_real_,
      overlap_genes = paste(overlap_genes, collapse = ";")
    ))
  }

  ft <- fisher.test(mat, alternative = "greater")

  tibble::tibble(
    universe_n = length(universe_genes),
    query_n = length(query_genes),
    target_n = length(target_genes),
    overlap_n = length(overlap_genes),
    expected_overlap_n = length(query_genes) * length(target_genes) / length(universe_genes),
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value,
    overlap_genes = paste(overlap_genes, collapse = ";")
  )
}

make_sig_label <- function(padj) {
  dplyr::case_when(
    !is.na(padj) & padj < 0.001 ~ "***",
    !is.na(padj) & padj < 0.01 ~ "**",
    !is.na(padj) & padj < 0.05 ~ "*",
    !is.na(padj) & padj < 0.1 ~ ".",
    TRUE ~ ""
  )
}

# -----------------------------
# Load AP-1 dependency table
# -----------------------------

ap1_tbl <- readr::read_csv(ap1_file, show_col_types = FALSE) |>
  dplyr::mutate(
    gene = as.character(gene),
    gene = stringr::str_replace(gene, "\\.\\d+$", ""),
    AP1_dependency_class = as.character(AP1_dependency_class),
    timepoint = as.character(timepoint),
    known_IL1_target = gene %in% known_il1_targets
  ) |>
  dplyr::filter(!is.na(gene), gene != "")

# -----------------------------
# Load and standardise fibroblast age/exercise table
# -----------------------------

age_ex_raw <- readRDS(age_ex_rds)
age_ex_df <- flatten_rds_to_tibble(age_ex_raw)

gene_col <- guess_col(
  age_ex_df,
  candidates = c(
    "gene", "genes", "Gene", "symbol", "Symbol", "gene_symbol",
    "gene_name", "external_gene_name", "rowname_gene"
  ),
  fallback_regex = "gene|symbol|rowname"
)

celltype_col <- guess_col(
  age_ex_df,
  candidates = c(
    "CellType", "celltype", "cell_type", "cluster", "annotation",
    "ident", "cell", "rds_list_name"
  ),
  fallback_regex = "cell.*type|cluster|annotation|ident|rds_list_name"
)

deg_col <- guess_col(
  age_ex_df,
  candidates = c(
    "DEG", "deg", "Deg", "DEG_class", "class", "category",
    "status", "gene_class"
  ),
  fallback_regex = "^deg$|class|category|status"
)

age_lfc_col <- guess_col(
  age_ex_df,
  candidates = c(
    "Age_log2FC", "Age_logFC", "Age_avg_log2FC",
    "log2FC_Age", "logFC_Age", "avg_log2FC_Age",
    "log2FoldChange.Age", "log2FC.Age", "logFC.Age",
    "Age.log2FC", "Age.logFC"
  ),
  fallback_regex = "age.*(log2|logfc|fc|fold)|((log2|logfc|fc|fold).*)age",
  exclude_regex = "padj|fdr|p_val|pvalue|p\\.val|qval|stat|se"
)

age_padj_col <- guess_col(
  age_ex_df,
  candidates = c(
    "Age_padj", "Age_FDR", "Age_adj_p", "Age_adj.P.Val",
    "padj_Age", "FDR_Age", "adj.P.Val_Age",
    "padj.Age", "FDR.Age", "adj.P.Val.Age"
  ),
  fallback_regex = "age.*(padj|fdr|qval|adj)|((padj|fdr|qval|adj).*)age"
)

ex_lfc_col <- guess_col(
  age_ex_df,
  candidates = c(
    "Exercise_log2FC", "Exercise_logFC", "Exercise_avg_log2FC",
    "Ex_log2FC", "Ex_logFC", "Ex_avg_log2FC",
    "log2FC_Exercise", "logFC_Exercise", "avg_log2FC_Exercise",
    "log2FC_Ex", "logFC_Ex", "avg_log2FC_Ex",
    "log2FoldChange.Ex", "log2FC.Ex", "logFC.Ex",
    "Exercise.log2FC", "Exercise.logFC",
    "Ex.log2FC", "Ex.logFC"
  ),
  fallback_regex = "(exercise|(^|[_.])ex($|[_.]))(.*)(log2|logfc|fc|fold)|((log2|logfc|fc|fold).*)(exercise|(^|[_.])ex($|[_.]))",
  exclude_regex = "padj|fdr|p_val|pvalue|p\\.val|qval|stat|se"
)

ex_padj_col <- guess_col(
  age_ex_df,
  candidates = c(
    "Exercise_padj", "Exercise_FDR", "Exercise_adj_p", "Exercise_adj.P.Val",
    "Ex_padj", "Ex_FDR", "Ex_adj_p", "Ex_adj.P.Val",
    "padj_Exercise", "FDR_Exercise", "adj.P.Val_Exercise",
    "padj_Ex", "FDR_Ex", "adj.P.Val_Ex",
    "padj.Ex", "FDR.Ex", "adj.P.Val.Ex"
  ),
  fallback_regex = "(exercise|(^|[_.])ex($|[_.]))(.*)(padj|fdr|qval|adj)|((padj|fdr|qval|adj).*)(exercise|(^|[_.])ex($|[_.]))"
)

if (is.na(gene_col)) {
  stop("Could not detect gene column in All_DEG_Age_Ex.rds. Columns are:\n", paste(names(age_ex_df), collapse = ", "))
}

if (is.na(celltype_col)) {
  stop("Could not detect cell type column in All_DEG_Age_Ex.rds. Columns are:\n", paste(names(age_ex_df), collapse = ", "))
}

if (is.na(age_lfc_col) || is.na(age_padj_col)) {
  stop(
    "Could not detect Age logFC/padj columns.\n",
    "Detected age_lfc_col: ", age_lfc_col, "\n",
    "Detected age_padj_col: ", age_padj_col, "\n",
    "Columns are:\n", paste(names(age_ex_df), collapse = ", ")
  )
}

if (is.na(ex_lfc_col) || is.na(ex_padj_col)) {
  stop(
    "Could not detect Exercise logFC/padj columns.\n",
    "Detected ex_lfc_col: ", ex_lfc_col, "\n",
    "Detected ex_padj_col: ", ex_padj_col, "\n",
    "Columns are:\n", paste(names(age_ex_df), collapse = ", ")
  )
}

message("Detected columns:")
message("  gene_col: ", gene_col)
message("  celltype_col: ", celltype_col)
message("  deg_col: ", deg_col)
message("  age_lfc_col: ", age_lfc_col)
message("  age_padj_col: ", age_padj_col)
message("  ex_lfc_col: ", ex_lfc_col)
message("  ex_padj_col: ", ex_padj_col)

fibro_age_ex <- age_ex_df |>
  dplyr::transmute(
    gene = as.character(.data[[gene_col]]),
    CellType_raw = as.character(.data[[celltype_col]]),
    DEG_label = if (!is.na(deg_col)) as.character(.data[[deg_col]]) else NA_character_,
    Age_log2FC = suppressWarnings(as.numeric(.data[[age_lfc_col]])),
    Age_padj = suppressWarnings(as.numeric(.data[[age_padj_col]])),
    Exercise_log2FC = suppressWarnings(as.numeric(.data[[ex_lfc_col]])),
    Exercise_padj = suppressWarnings(as.numeric(.data[[ex_padj_col]]))
  ) |>
  dplyr::mutate(
    gene = stringr::str_replace(gene, "\\.\\d+$", ""),
    CellType = dplyr::case_when(
      stringr::str_to_lower(CellType_raw) %in% c("fibro", "fibroblast", "fibroblasts") ~ "Fibroblast",
      stringr::str_detect(stringr::str_to_lower(CellType_raw), "^fibro") ~ "Fibroblast",
      TRUE ~ CellType_raw
    )
  ) |>
  dplyr::filter(CellType == "Fibroblast") |>
  dplyr::filter(!is.na(gene), gene != "") |>
  dplyr::mutate(
    DEG_label_lower = stringr::str_to_lower(tidyr::replace_na(DEG_label, "")),

    Age_sig = !is.na(Age_padj) &
      Age_padj <= age_ex_fdr_cutoff &
      !is.na(Age_log2FC) &
      abs(Age_log2FC) >= age_ex_abs_log2fc_cutoff,

    Exercise_sig = !is.na(Exercise_padj) &
      Exercise_padj <= age_ex_fdr_cutoff &
      !is.na(Exercise_log2FC) &
      abs(Exercise_log2FC) >= age_ex_abs_log2fc_cutoff,

    Age_direction = dplyr::case_when(
      Age_sig & Age_log2FC > 0 ~ "Age_up",
      Age_sig & Age_log2FC < 0 ~ "Age_down",
      TRUE ~ "Age_NS"
    ),

    Exercise_direction = dplyr::case_when(
      Exercise_sig & Exercise_log2FC > 0 ~ "Exercise_up",
      Exercise_sig & Exercise_log2FC < 0 ~ "Exercise_down",
      TRUE ~ "Exercise_NS"
    ),

    Age_up = Age_sig & Age_log2FC > 0,
    Age_down = Age_sig & Age_log2FC < 0,
    Exercise_up = Exercise_sig & Exercise_log2FC > 0,
    Exercise_down = Exercise_sig & Exercise_log2FC < 0,

    Age_up_Exercise_down = Age_up & Exercise_down,
    Age_down_Exercise_up = Age_down & Exercise_up,

    Reversed_by_exercise = Age_sig &
      Exercise_sig &
      sign(Age_log2FC) == -sign(Exercise_log2FC),

    Exacerbated_by_exercise = Age_sig &
      Exercise_sig &
      sign(Age_log2FC) == sign(Exercise_log2FC),

    Reversed_by_exercise = Reversed_by_exercise |
      stringr::str_detect(DEG_label_lower, "revers"),

    Exacerbated_by_exercise = Exacerbated_by_exercise |
      stringr::str_detect(DEG_label_lower, "exacerb"),

    Fibroblast_AgeEx_status = dplyr::case_when(
      Age_up_Exercise_down ~ "Age up + Exercise down",
      Age_down_Exercise_up ~ "Age down + Exercise up",
      Reversed_by_exercise ~ "Reversed by exercise",
      Age_up & !Exercise_down ~ "Age up only",
      Exercise_down & !Age_up ~ "Exercise down only",
      Age_sig | Exercise_sig ~ "Other fibroblast DEG",
      TRUE ~ "Not fibroblast Age/Exercise DEG"
    )
  ) |>
  dplyr::distinct(gene, .keep_all = TRUE)

if (nrow(fibro_age_ex) == 0) {
  stop("No fibroblast rows found after filtering. Check CellType labels in the RDS.")
}

readr::write_csv(
  fibro_age_ex,
  file.path(out_dir, "fibroblast_AgeEx_standardised.csv")
)

message("Fibroblast rows: ", nrow(fibro_age_ex))
message("Fibroblast genes: ", length(unique(fibro_age_ex$gene)))

# -----------------------------
# Build AP-1 gene summaries
# -----------------------------

ap1_induced_summary <- ap1_tbl |>
  dplyr::filter(AP1_dependency_class == "IL1_induced_AP1_dependent") |>
  dplyr::group_by(gene) |>
  dplyr::summarise(
    AP1_dependent_timepoints = paste(sort(unique(timepoint)), collapse = ";"),
    AP1_dependent_n_timepoints = dplyr::n_distinct(timepoint),
    max_AP1_dependency_score = max(AP1_dependency_score, na.rm = TRUE),

    IL1_log2FC_3h = safe_first_numeric(IL1_log2FC[timepoint == "3h"]),
    IL1_padj_3h = safe_first_numeric(IL1_padj[timepoint == "3h"]),
    T5224preIL1_vs_IL1_log2FC_3h = safe_first_numeric(T5224preIL1_vs_IL1_log2FC[timepoint == "3h"]),
    T5224preIL1_vs_IL1_padj_3h = safe_first_numeric(T5224preIL1_vs_IL1_padj[timepoint == "3h"]),

    IL1_log2FC_6h = safe_first_numeric(IL1_log2FC[timepoint == "6h"]),
    IL1_padj_6h = safe_first_numeric(IL1_padj[timepoint == "6h"]),
    T5224preIL1_vs_IL1_log2FC_6h = safe_first_numeric(T5224preIL1_vs_IL1_log2FC[timepoint == "6h"]),
    T5224preIL1_vs_IL1_padj_6h = safe_first_numeric(T5224preIL1_vs_IL1_padj[timepoint == "6h"]),

    known_IL1_target = any(known_IL1_target),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    dplyr::desc(AP1_dependent_n_timepoints),
    dplyr::desc(max_AP1_dependency_score)
  )

readr::write_csv(
  ap1_induced_summary,
  file.path(out_dir, "IL1_induced_AP1_dependent_gene_summary.csv")
)

# -----------------------------
# Combined fibroblast overlap table
# -----------------------------

fibro_ap1_overlap <- ap1_induced_summary |>
  dplyr::inner_join(fibro_age_ex, by = "gene") |>
  dplyr::mutate(
    high_confidence_model_gene = Age_up_Exercise_down | Reversed_by_exercise,
    very_high_confidence_model_gene = Age_up_Exercise_down,
    model_rank = dplyr::case_when(
      Age_up_Exercise_down ~ 1L,
      Reversed_by_exercise ~ 2L,
      Age_up ~ 3L,
      Exercise_down ~ 4L,
      TRUE ~ 5L
    )
  ) |>
  dplyr::arrange(
    model_rank,
    dplyr::desc(AP1_dependent_n_timepoints),
    dplyr::desc(max_AP1_dependency_score)
  )

readr::write_csv(
  fibro_ap1_overlap,
  file.path(out_dir, "fibroblast_IL1_induced_AP1_dependent_AgeEx_overlap_gene_table.csv")
)

high_confidence_genes <- fibro_ap1_overlap |>
  dplyr::filter(high_confidence_model_gene) |>
  dplyr::arrange(
    dplyr::desc(very_high_confidence_model_gene),
    dplyr::desc(AP1_dependent_n_timepoints),
    dplyr::desc(max_AP1_dependency_score)
  )

readr::write_csv(
  high_confidence_genes,
  file.path(out_dir, "HIGH_CONFIDENCE_fibroblast_IL1_AP1_AgeEx_overlap_genes.csv")
)

very_high_confidence_genes <- fibro_ap1_overlap |>
  dplyr::filter(very_high_confidence_model_gene) |>
  dplyr::arrange(
    dplyr::desc(AP1_dependent_n_timepoints),
    dplyr::desc(max_AP1_dependency_score)
  )

readr::write_csv(
  very_high_confidence_genes,
  file.path(out_dir, "VERY_HIGH_CONFIDENCE_AgeUp_ExerciseDown_IL1_AP1_genes.csv")
)

# -----------------------------
# Fibroblast DEG set membership
# -----------------------------

fibro_gene_sets <- fibro_age_ex |>
  dplyr::transmute(
    gene,
    Age_up,
    Age_down,
    Exercise_up,
    Exercise_down,
    Age_up_Exercise_down,
    Age_down_Exercise_up,
    Reversed_by_exercise,
    Exacerbated_by_exercise
  ) |>
  tidyr::pivot_longer(
    cols = -gene,
    names_to = "fibroblast_category",
    values_to = "in_category"
  ) |>
  dplyr::filter(in_category) |>
  dplyr::select(gene, fibroblast_category) |>
  dplyr::distinct()

readr::write_csv(
  fibro_gene_sets,
  file.path(out_dir, "fibroblast_AgeEx_gene_sets_long.csv")
)

fibro_set_sizes <- fibro_gene_sets |>
  dplyr::count(fibroblast_category, name = "n_genes") |>
  dplyr::arrange(dplyr::desc(n_genes))

readr::write_csv(
  fibro_set_sizes,
  file.path(out_dir, "fibroblast_AgeEx_gene_set_sizes.csv")
)

# -----------------------------
# Fisher enrichment:
# IL1/AP1-dependent sets vs fibroblast age/exercise sets
# -----------------------------

ap1_sets_for_testing <- dplyr::bind_rows(
  ap1_tbl |>
    dplyr::filter(AP1_dependency_class == "IL1_induced_AP1_dependent") |>
    dplyr::mutate(
      AP1_set = paste0(timepoint, "__IL1_induced_AP1_dependent"),
      AP1_set_plot = paste0(timepoint, ": IL-1 induced, AP-1 dependent")
    ) |>
    dplyr::select(gene, AP1_set, AP1_set_plot, timepoint, AP1_dependency_class),

  ap1_tbl |>
    dplyr::filter(AP1_dependency_class == "IL1_repressed_AP1_dependent") |>
    dplyr::mutate(
      AP1_set = paste0(timepoint, "__IL1_repressed_AP1_dependent"),
      AP1_set_plot = paste0(timepoint, ": IL-1 repressed, AP-1 dependent")
    ) |>
    dplyr::select(gene, AP1_set, AP1_set_plot, timepoint, AP1_dependency_class),

  ap1_tbl |>
    dplyr::filter(AP1_dependency_class == "IL1_induced_AP1_dependent") |>
    dplyr::mutate(
      AP1_set = "any_timepoint__IL1_induced_AP1_dependent",
      AP1_set_plot = "Any timepoint: IL-1 induced, AP-1 dependent",
      timepoint = "any_timepoint"
    ) |>
    dplyr::select(gene, AP1_set, AP1_set_plot, timepoint, AP1_dependency_class),

  ap1_tbl |>
    dplyr::filter(AP1_dependency_class == "IL1_repressed_AP1_dependent") |>
    dplyr::mutate(
      AP1_set = "any_timepoint__IL1_repressed_AP1_dependent",
      AP1_set_plot = "Any timepoint: IL-1 repressed, AP-1 dependent",
      timepoint = "any_timepoint"
    ) |>
    dplyr::select(gene, AP1_set, AP1_set_plot, timepoint, AP1_dependency_class)
) |>
  dplyr::distinct()

universe_genes <- intersect(
  unique(ap1_tbl$gene),
  unique(fibro_age_ex$gene)
)

test_grid <- tidyr::expand_grid(
  AP1_set = unique(ap1_sets_for_testing$AP1_set),
  fibroblast_category = c(
    "Age_up",
    "Age_down",
    "Exercise_up",
    "Exercise_down",
    "Age_up_Exercise_down",
    "Age_down_Exercise_up",
    "Reversed_by_exercise",
    "Exacerbated_by_exercise"
  )
)

overlap_enrichment <- purrr::map_dfr(seq_len(nrow(test_grid)), function(i) {
  this_ap1_set <- test_grid$AP1_set[i]
  this_category <- test_grid$fibroblast_category[i]

  query_genes <- ap1_sets_for_testing |>
    dplyr::filter(AP1_set == this_ap1_set) |>
    dplyr::pull(gene) |>
    unique()

  target_genes <- fibro_gene_sets |>
    dplyr::filter(fibroblast_category == this_category) |>
    dplyr::pull(gene) |>
    unique()

  meta <- ap1_sets_for_testing |>
    dplyr::filter(AP1_set == this_ap1_set) |>
    dplyr::distinct(AP1_set, AP1_set_plot, timepoint, AP1_dependency_class) |>
    dplyr::slice(1)

  run_fisher_overlap(
    query_genes = query_genes,
    target_genes = target_genes,
    universe_genes = universe_genes
  ) |>
    dplyr::mutate(
      AP1_set = this_ap1_set,
      AP1_set_plot = meta$AP1_set_plot,
      timepoint = meta$timepoint,
      AP1_dependency_class = meta$AP1_dependency_class,
      fibroblast_category = this_category,
      .before = 1
    )
}) |>
  dplyr::mutate(
    padj_global = p.adjust(p_value, method = "BH"),
    padj_within_AP1_set = ave(
      p_value,
      AP1_set,
      FUN = function(x) p.adjust(x, method = "BH")
    ),
    log2_odds_ratio = log2(odds_ratio),
    sig_label = make_sig_label(padj_global)
  ) |>
  dplyr::arrange(padj_global, p_value)

readr::write_csv(
  overlap_enrichment,
  file.path(out_dir, "fibroblast_IL1_AP1_vs_AgeEx_enrichment.csv")
)

overlap_genes_long <- overlap_enrichment |>
  dplyr::filter(!is.na(overlap_genes), overlap_genes != "") |>
  dplyr::select(
    AP1_set,
    AP1_set_plot,
    timepoint,
    AP1_dependency_class,
    fibroblast_category,
    overlap_genes
  ) |>
  tidyr::separate_rows(overlap_genes, sep = ";") |>
  dplyr::rename(gene = overlap_genes) |>
  dplyr::left_join(fibro_age_ex, by = "gene") |>
  dplyr::left_join(ap1_induced_summary, by = "gene") |>
  dplyr::arrange(
    fibroblast_category,
    AP1_set,
    dplyr::desc(max_AP1_dependency_score)
  )

readr::write_csv(
  overlap_genes_long,
  file.path(out_dir, "fibroblast_IL1_AP1_vs_AgeEx_overlap_genes_long.csv")
)

# ============================================================
# PLOTS
# ============================================================

# -----------------------------
# Plot 1:
# Enrichment dotplot.
# This is the main statistical plot.
# -----------------------------

plot_enrich_df <- overlap_enrichment |>
  dplyr::filter(
    AP1_dependency_class == "IL1_induced_AP1_dependent",
    fibroblast_category %in% c(
      "Age_up",
      "Age_down",
      "Exercise_up",
      "Exercise_down",
      "Age_up_Exercise_down",
      "Age_down_Exercise_up",
      "Reversed_by_exercise",
      "Exacerbated_by_exercise"
    )
  ) |>
  dplyr::mutate(
    fibroblast_category = factor(
      fibroblast_category,
      levels = c(
        "Age_up",
        "Exercise_down",
        "Age_up_Exercise_down",
        "Reversed_by_exercise",
        "Age_down",
        "Exercise_up",
        "Age_down_Exercise_up",
        "Exacerbated_by_exercise"
      ),
      labels = c(
        "Age up",
        "Exercise down",
        "Age up + Exercise down",
        "Reversed by exercise",
        "Age down",
        "Exercise up",
        "Age down + Exercise up",
        "Exacerbated by exercise"
      )
    ),
    AP1_set_plot = factor(
      AP1_set_plot,
      levels = c(
        "3h: IL-1 induced, AP-1 dependent",
        "6h: IL-1 induced, AP-1 dependent",
        "Any timepoint: IL-1 induced, AP-1 dependent"
      )
    )
  ) |>
  dplyr::filter(!is.na(AP1_set_plot), !is.na(fibroblast_category))

p_enrichment <- ggplot2::ggplot(
  plot_enrich_df,
  ggplot2::aes(
    x = fibroblast_category,
    y = AP1_set_plot,
    size = overlap_n,
    fill = log2_odds_ratio
  )
) +
  ggplot2::geom_point(shape = 21, colour = "black", alpha = 0.95) +
  ggplot2::geom_text(
    ggplot2::aes(label = sig_label),
    size = 4,
    colour = "black",
    vjust = 0.35
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#377EB8",
    mid = "white",
    high = "#E41A1C",
    midpoint = 0,
    na.value = "grey85"
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    size = "Overlap genes",
    fill = "log2 odds ratio",
    title = "Fibroblast age/exercise overlap with IL-1-induced AP-1-dependent genes",
    subtitle = "Fisher test. Universe = genes measured in both IL-1/AP-1 and fibroblast age/exercise datasets"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
    panel.grid.minor = ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(out_dir, "PLOT1_fibroblast_IL1_AP1_vs_AgeEx_enrichment_dotplot.pdf"),
  p_enrichment,
  width = 11,
  height = 4.8
)

# -----------------------------
# Plot 2:
# Fibroblast Age vs Exercise effect sizes.
# This shows whether IL1/AP1 genes fall in the reversal quadrant.
# -----------------------------

fibro_age_ex_scatter <- fibro_age_ex |>
  dplyr::left_join(ap1_induced_summary, by = "gene") |>
  dplyr::mutate(
    IL1_induced_AP1_dependent = !is.na(max_AP1_dependency_score),
    AP1_label_group = dplyr::case_when(
      IL1_induced_AP1_dependent & Age_up_Exercise_down ~ "IL1/AP1 dependent + Age up/Exercise down",
      IL1_induced_AP1_dependent & Reversed_by_exercise ~ "IL1/AP1 dependent + Reversed",
      IL1_induced_AP1_dependent ~ "IL1/AP1 dependent",
      TRUE ~ "Other fibroblast genes"
    ),
    label_gene = IL1_induced_AP1_dependent & (
      Age_up_Exercise_down | Reversed_by_exercise
    )
  )

label_age_ex <- fibro_age_ex_scatter |>
  dplyr::filter(label_gene) |>
  dplyr::arrange(
    dplyr::desc(Age_up_Exercise_down),
    dplyr::desc(max_AP1_dependency_score),
    Age_padj
  ) |>
  dplyr::slice_head(n = top_label_n)

p_age_ex_scatter <- ggplot2::ggplot(
  fibro_age_ex_scatter,
  ggplot2::aes(x = Age_log2FC, y = Exercise_log2FC)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  ggplot2::annotate(
    "rect",
    xmin = 0,
    xmax = Inf,
    ymin = -Inf,
    ymax = 0,
    alpha = 0.08,
    fill = "#A81B1B"
  ) +
  ggplot2::geom_point(
    data = fibro_age_ex_scatter |>
      dplyr::filter(AP1_label_group == "Other fibroblast genes"),
    colour = "grey82",
    size = 1.1,
    alpha = 0.45
  ) +
  ggplot2::geom_point(
    data = fibro_age_ex_scatter |>
      dplyr::filter(AP1_label_group != "Other fibroblast genes"),
    ggplot2::aes(colour = AP1_label_group, size = max_AP1_dependency_score),
    alpha = 0.85
  ) +
  ggplot2::geom_text(
    data = label_age_ex,
    ggplot2::aes(label = gene),
    size = 2.7,
    check_overlap = TRUE,
    vjust = -0.5
  ) +
  ggplot2::scale_colour_manual(
    values = c(
      "IL1/AP1 dependent + Age up/Exercise down" = "#A81B1B",
      "IL1/AP1 dependent + Reversed" = "#A53D99",
      "IL1/AP1 dependent" = "#42977E"
    )
  ) +
  ggplot2::labs(
    x = "Fibroblast age effect log2FC",
    y = "Fibroblast exercise effect log2FC",
    colour = NULL,
    size = "AP-1 dependency score",
    title = "Fibroblast age versus exercise effects",
    subtitle = "Red quadrant marks genes increased with age and decreased by exercise"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  file.path(out_dir, "PLOT2_fibroblast_Age_vs_Exercise_highlight_IL1_AP1_genes.pdf"),
  p_age_ex_scatter,
  width = 8,
  height = 7
)

# -----------------------------
# Plot 3:
# IL1 response versus AP-1 inhibitor response,
# coloured by fibroblast age/exercise status.
# This shows whether the AP-1 dependent IL1 response is also age/exercise relevant.
# -----------------------------

ap1_scatter_with_fibro <- ap1_tbl |>
  dplyr::left_join(
    fibro_age_ex |>
      dplyr::select(
        gene,
        Age_log2FC,
        Age_padj,
        Exercise_log2FC,
        Exercise_padj,
        Age_up,
        Age_down,
        Exercise_up,
        Exercise_down,
        Age_up_Exercise_down,
        Age_down_Exercise_up,
        Reversed_by_exercise,
        Exacerbated_by_exercise,
        Fibroblast_AgeEx_status
      ),
    by = "gene"
  ) |>
  dplyr::mutate(
    Fibroblast_AgeEx_status = tidyr::replace_na(
      Fibroblast_AgeEx_status,
      "Not fibroblast Age/Exercise DEG"
    ),
    Fibroblast_AgeEx_status = factor(
      Fibroblast_AgeEx_status,
      levels = names(status_cols)
    ),
    is_key_AP1_class = AP1_dependency_class == "IL1_induced_AP1_dependent",
    label_gene = is_key_AP1_class &
      Fibroblast_AgeEx_status %in% c(
        "Age up + Exercise down",
        "Reversed by exercise"
      )
  )

label_il1_scatter <- ap1_scatter_with_fibro |>
  dplyr::filter(label_gene) |>
  dplyr::group_by(timepoint) |>
  dplyr::arrange(
    dplyr::desc(Age_up_Exercise_down),
    dplyr::desc(AP1_dependency_score)
  ) |>
  dplyr::slice_head(n = top_label_n) |>
  dplyr::ungroup()

p_il1_scatter <- ggplot2::ggplot(
  ap1_scatter_with_fibro,
  ggplot2::aes(x = IL1_log2FC, y = T5224preIL1_vs_IL1_log2FC)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  ggplot2::annotate(
    "rect",
    xmin = 0,
    xmax = Inf,
    ymin = -Inf,
    ymax = 0,
    alpha = 0.08,
    fill = "#A81B1B"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(colour = Fibroblast_AgeEx_status),
    alpha = 0.65,
    size = 1.7
  ) +
  ggplot2::geom_point(
    data = ap1_scatter_with_fibro |>
      dplyr::filter(is_key_AP1_class),
    shape = 21,
    fill = NA,
    colour = "black",
    size = 2.2,
    stroke = 0.55
  ) +
  ggplot2::geom_text(
    data = label_il1_scatter,
    ggplot2::aes(label = gene),
    size = 2.6,
    check_overlap = TRUE,
    vjust = -0.5
  ) +
  ggplot2::facet_wrap(~ timepoint) +
  ggplot2::scale_colour_manual(values = status_cols, drop = FALSE) +
  ggplot2::labs(
    x = "IL1 vs Vehicle log2FC",
    y = "T5224preIL1 vs IL1 log2FC",
    colour = "Fibroblast age/exercise status",
    title = "IL-1/AP-1 dependency space coloured by fibroblast age/exercise status",
    subtitle = "Key quadrant: IL1-induced and suppressed by AP-1 inhibitor. Black rings mark IL1-induced AP1-dependent genes"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  file.path(out_dir, "PLOT3_IL1_vs_T5224_coloured_by_fibroblast_AgeEx_status.pdf"),
  p_il1_scatter,
  width = 11,
  height = 7
)

# -----------------------------
# Plot 4:
# Counts of AP1-dependent IL1 genes by fibroblast age/exercise category.
# This is easier to read than a huge upset plot.
# -----------------------------

ap1_count_by_status <- ap1_tbl |>
  dplyr::filter(AP1_dependency_class == "IL1_induced_AP1_dependent") |>
  dplyr::left_join(
    fibro_age_ex |>
      dplyr::select(gene, Fibroblast_AgeEx_status),
    by = "gene"
  ) |>
  dplyr::mutate(
    Fibroblast_AgeEx_status = tidyr::replace_na(
      Fibroblast_AgeEx_status,
      "Not fibroblast Age/Exercise DEG"
    ),
    Fibroblast_AgeEx_status = factor(
      Fibroblast_AgeEx_status,
      levels = names(status_cols)
    )
  ) |>
  dplyr::count(timepoint, Fibroblast_AgeEx_status, name = "n_genes") |>
  dplyr::group_by(timepoint) |>
  dplyr::mutate(percent = 100 * n_genes / sum(n_genes)) |>
  dplyr::ungroup()

readr::write_csv(
  ap1_count_by_status,
  file.path(out_dir, "IL1_induced_AP1_dependent_counts_by_fibroblast_AgeEx_status.csv")
)

p_counts <- ggplot2::ggplot(
  ap1_count_by_status,
  ggplot2::aes(x = Fibroblast_AgeEx_status, y = n_genes, fill = Fibroblast_AgeEx_status)
) +
  ggplot2::geom_col(colour = "black", linewidth = 0.2) +
  ggplot2::geom_text(
    ggplot2::aes(label = n_genes),
    vjust = -0.25,
    size = 3
  ) +
  ggplot2::facet_wrap(~ timepoint, nrow = 1) +
  ggplot2::scale_fill_manual(values = status_cols, drop = FALSE) +
  ggplot2::labs(
    x = NULL,
    y = "Genes",
    fill = NULL,
    title = "Fibroblast age/exercise status of IL-1-induced AP-1-dependent genes"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "none"
  )

ggplot2::ggsave(
  file.path(out_dir, "PLOT4_counts_by_fibroblast_AgeEx_status.pdf"),
  p_counts,
  width = 10,
  height = 5
)

top_heatmap_genes <- fibro_ap1_overlap %>%
  dplyr::arrange(
    model_rank,
    dplyr::desc(AP1_dependent_n_timepoints),
    dplyr::desc(max_AP1_dependency_score)
  ) %>%
  dplyr::slice_head(n = top_heatmap_n) %>%
  dplyr::pull(gene) %>%
  unique()

heatmap_ap1_base <- ap1_tbl %>%
  dplyr::filter(gene %in% top_heatmap_genes) %>%
  dplyr::select(
    gene,
    timepoint,
    IL1_log2FC,
    IL1_padj,
    T5224preIL1_vs_IL1_log2FC,
    T5224preIL1_vs_IL1_padj
  ) %>%
  dplyr::mutate(
    effect_IL1 = paste0("IL1 vs Veh ", timepoint),
    effect_T5224 = paste0("T5224preIL1 vs IL1 ", timepoint)
  )

heatmap_ap1_values <- dplyr::bind_rows(
  heatmap_ap1_base %>%
    dplyr::transmute(
      gene = gene,
      effect = effect_IL1,
      log2FC = IL1_log2FC,
      padj = IL1_padj
    ),
  heatmap_ap1_base %>%
    dplyr::transmute(
      gene = gene,
      effect = effect_T5224,
      log2FC = T5224preIL1_vs_IL1_log2FC,
      padj = T5224preIL1_vs_IL1_padj
    )
)

heatmap_fibro_base <- fibro_age_ex %>%
  dplyr::filter(gene %in% top_heatmap_genes) %>%
  dplyr::select(
    gene,
    Age_log2FC,
    Age_padj,
    Exercise_log2FC,
    Exercise_padj
  )

heatmap_fibro_values <- dplyr::bind_rows(
  heatmap_fibro_base %>%
    dplyr::transmute(
      gene = gene,
      effect = "Fibroblast Age effect",
      log2FC = Age_log2FC,
      padj = Age_padj
    ),
  heatmap_fibro_base %>%
    dplyr::transmute(
      gene = gene,
      effect = "Fibroblast Exercise effect",
      log2FC = Exercise_log2FC,
      padj = Exercise_padj
    )
)

heatmap_df <- dplyr::bind_rows(
  heatmap_ap1_values,
  heatmap_fibro_values
) %>%
  dplyr::mutate(
    gene = factor(gene, levels = rev(top_heatmap_genes)),
    effect = factor(
      effect,
      levels = c(
        "IL1 vs Veh 3h",
        "T5224preIL1 vs IL1 3h",
        "IL1 vs Veh 6h",
        "T5224preIL1 vs IL1 6h",
        "Fibroblast Age effect",
        "Fibroblast Exercise effect"
      )
    ),
    sig_label = ifelse(!is.na(padj) & padj <= 0.05, "*", "")
  ) %>%
  dplyr::filter(!is.na(effect))

p_heatmap <- ggplot2::ggplot(
  heatmap_df,
  ggplot2::aes(x = effect, y = gene, fill = log2FC)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.18) +
  ggplot2::geom_text(
    ggplot2::aes(label = sig_label),
    size = 2.4
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#377EB8",
    mid = "white",
    high = "#E41A1C",
    midpoint = 0,
    na.value = "grey90"
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    fill = "log2FC",
    title = "Top fibroblast-relevant IL-1-induced AP-1-dependent genes",
    subtitle = "Genes ranked by age-up/exercise-down status, persistence across timepoints, and AP-1 dependency score. * = FDR <= 0.05"
  ) +
  ggplot2::theme_bw(base_size = 8.5) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
    panel.grid = ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(out_dir, "PLOT5_top_fibroblast_relevant_IL1_AP1_genes_heatmap.pdf"),
  p_heatmap,
  width = 10,
  height = max(6, length(top_heatmap_genes) * 0.16)
)
# -----------------------------
# Plot 6:
# Ranked high-confidence genes.
# This is a compact gene-prioritisation plot.
# -----------------------------

rank_plot_df <- high_confidence_genes |>
  dplyr::mutate(
    gene = factor(gene, levels = rev(unique(gene))),
    Fibroblast_AgeEx_status = factor(
      Fibroblast_AgeEx_status,
      levels = names(status_cols)
    )
  ) |>
  dplyr::slice_head(n = 40)

if (nrow(rank_plot_df) > 0) {
  p_rank <- ggplot2::ggplot(
    rank_plot_df,
    ggplot2::aes(
      x = max_AP1_dependency_score,
      y = gene,
      fill = Fibroblast_AgeEx_status
    )
  ) +
    ggplot2::geom_col(colour = "black", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = status_cols, drop = FALSE) +
    ggplot2::labs(
      x = "Maximum AP-1 dependency score",
      y = NULL,
      fill = "Fibroblast age/exercise status",
      title = "Prioritised fibroblast IL-1/AP-1 genes",
      subtitle = "High-confidence genes are reversed by exercise or age-up/exercise-down"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )

  ggplot2::ggsave(
    file.path(out_dir, "PLOT6_ranked_high_confidence_fibroblast_IL1_AP1_genes.pdf"),
    p_rank,
    width = 8,
    height = max(5, nrow(rank_plot_df) * 0.18)
  )
}

# -----------------------------
# Print summaries
# -----------------------------

message("\nDone.")
message("Output directory:")
message(out_dir)

message("\nMain tables:")
message(file.path(out_dir, "fibroblast_IL1_induced_AP1_dependent_AgeEx_overlap_gene_table.csv"))
message(file.path(out_dir, "HIGH_CONFIDENCE_fibroblast_IL1_AP1_AgeEx_overlap_genes.csv"))
message(file.path(out_dir, "VERY_HIGH_CONFIDENCE_AgeUp_ExerciseDown_IL1_AP1_genes.csv"))
message(file.path(out_dir, "fibroblast_IL1_AP1_vs_AgeEx_enrichment.csv"))
message(file.path(out_dir, "fibroblast_IL1_AP1_vs_AgeEx_overlap_genes_long.csv"))

message("\nMain plots:")
message(file.path(out_dir, "PLOT1_fibroblast_IL1_AP1_vs_AgeEx_enrichment_dotplot.pdf"))
message(file.path(out_dir, "PLOT2_fibroblast_Age_vs_Exercise_highlight_IL1_AP1_genes.pdf"))
message(file.path(out_dir, "PLOT3_IL1_vs_T5224_coloured_by_fibroblast_AgeEx_status.pdf"))
message(file.path(out_dir, "PLOT4_counts_by_fibroblast_AgeEx_status.pdf"))
message(file.path(out_dir, "PLOT5_top_fibroblast_relevant_IL1_AP1_genes_heatmap.pdf"))

if (nrow(high_confidence_genes) > 0) {
  message(file.path(out_dir, "PLOT6_ranked_high_confidence_fibroblast_IL1_AP1_genes.pdf"))
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(forcats)
})

# ============================================================
# Basic overlap plots:
#   1. IL-1 vs Vehicle responsive genes overlapped with fibroblast Age/Exercise genes
#   2. AP-1 inhibitor vs Vehicle responsive genes overlapped with fibroblast Age/Exercise genes
#
# This does NOT use T5224preIL1_vs_IL1.
# ============================================================

base_dir <- publication_path("PAPER_OUTPUT_ROOT", "IL1_RNASeq")

age_ex_rds <- file.path(base_dir, "All_DEG_Age_Ex.rds")

out_dir <- file.path(base_dir, "Basic_IL1_vs_Control_Inhibitor_vs_Control_AgeEx_overlap")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dge_fdr_cutoff <- 0.05
dge_abs_log2fc_cutoff <- 0

age_ex_fdr_cutoff <- 0.05
age_ex_abs_log2fc_cutoff <- 0

top_label_n <- 30
top_heatmap_n <- 60

# -----------------------------
# Helper functions
# -----------------------------

guess_col <- function(df, candidates, fallback_regex = NULL, exclude_regex = NULL) {
  nms <- names(df)
  nms_clean <- tolower(nms)
  candidates_clean <- tolower(candidates)

  exact_hit <- which(nms_clean %in% candidates_clean)

  if (!is.null(exclude_regex) && length(exact_hit) > 0) {
    exact_hit <- exact_hit[!stringr::str_detect(nms_clean[exact_hit], exclude_regex)]
  }

  if (length(exact_hit) > 0) {
    return(nms[exact_hit[1]])
  }

  if (!is.null(fallback_regex)) {
    regex_hit <- which(stringr::str_detect(nms_clean, fallback_regex))

    if (!is.null(exclude_regex) && length(regex_hit) > 0) {
      regex_hit <- regex_hit[!stringr::str_detect(nms_clean[regex_hit], exclude_regex)]
    }

    if (length(regex_hit) > 0) {
      return(nms[regex_hit[1]])
    }
  }

  NA_character_
}

flatten_rds_to_tibble <- function(obj) {
  if (is.data.frame(obj)) {
    df <- as.data.frame(obj)

    if (
      !"gene" %in% tolower(names(df)) &&
        !all(rownames(df) == as.character(seq_len(nrow(df))))
    ) {
      df <- tibble::rownames_to_column(df, var = "rowname_gene")
    }

    return(tibble::as_tibble(df))
  }

  if (is.list(obj)) {
    df_list <- purrr::imap(obj, function(x, nm) {
      if (!is.data.frame(x)) return(NULL)

      df <- as.data.frame(x)

      if (
        !"gene" %in% tolower(names(df)) &&
          !all(rownames(df) == as.character(seq_len(nrow(df))))
      ) {
        df <- tibble::rownames_to_column(df, var = "rowname_gene")
      }

      df$rds_list_name <- nm
      tibble::as_tibble(df)
    })

    df_list <- purrr::compact(df_list)

    if (length(df_list) == 0) {
      stop("RDS is a list, but no data.frames were found inside it.")
    }

    return(dplyr::bind_rows(df_list))
  }

  stop("Unsupported RDS object type: ", paste(class(obj), collapse = ", "))
}

read_dge_file <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)

  gene_col <- guess_col(
    df,
    candidates = c(
      "gene", "genes", "symbol", "gene_symbol", "gene_name",
      "external_gene_name", "feature", "features", "X", "...1"
    ),
    fallback_regex = "gene|symbol|feature"
  )

  p_col <- guess_col(
    df,
    candidates = c(
      "padj", "p_val_adj", "p_val_adj.", "FDR", "adj.P.Val",
      "qval", "q_value", "qvalue", "fdr_bh"
    ),
    fallback_regex = "padj|adj.*p|p.*adj|fdr|qval|q_value"
  )

  lfc_col <- guess_col(
    df,
    candidates = c(
      "log2FoldChange", "logFC", "avg_log2FC", "avg_logFC",
      "log2FC", "estimate", "coef"
    ),
    fallback_regex = "log2.*fold|logfc|log2fc|avg_log",
    exclude_regex = "padj|fdr|p_val|pvalue|qval"
  )

  if (is.na(gene_col)) stop("Could not detect gene column in: ", basename(path))
  if (is.na(p_col)) stop("Could not detect FDR/padj column in: ", basename(path))
  if (is.na(lfc_col)) stop("Could not detect logFC column in: ", basename(path))

  contrast <- basename(path) %>%
    stringr::str_remove("^DGE_") %>%
    stringr::str_remove("\\.csv$")

  df %>%
    dplyr::transmute(
      gene = as.character(.data[[gene_col]]),
      contrast = contrast,
      log2FC = suppressWarnings(as.numeric(.data[[lfc_col]])),
      padj = suppressWarnings(as.numeric(.data[[p_col]])),
      gene_col_used = gene_col,
      lfc_col_used = lfc_col,
      p_col_used = p_col
    ) %>%
    dplyr::mutate(
      gene = stringr::str_replace(gene, "\\.\\d+$", ""),
      sig = !is.na(padj) &
        padj <= dge_fdr_cutoff &
        !is.na(log2FC) &
        abs(log2FC) >= dge_abs_log2fc_cutoff,
      direction = dplyr::case_when(
        sig & log2FC > 0 ~ "up",
        sig & log2FC < 0 ~ "down",
        TRUE ~ "ns"
      )
    ) %>%
    dplyr::filter(!is.na(gene), gene != "", !is.na(log2FC), !is.na(padj)) %>%
    dplyr::distinct(gene, contrast, .keep_all = TRUE)
}

run_fisher_overlap <- function(query_genes, target_genes, universe_genes) {
  universe_genes <- unique(universe_genes)
  query_genes <- unique(intersect(query_genes, universe_genes))
  target_genes <- unique(intersect(target_genes, universe_genes))

  in_query <- universe_genes %in% query_genes
  in_target <- universe_genes %in% target_genes

  overlap_genes <- sort(intersect(query_genes, target_genes))

  mat <- matrix(
    c(
      sum(in_query & in_target),
      sum(in_query & !in_target),
      sum(!in_query & in_target),
      sum(!in_query & !in_target)
    ),
    nrow = 2,
    byrow = TRUE
  )

  if (
    length(universe_genes) == 0 ||
      any(rowSums(mat) == 0) ||
      any(colSums(mat) == 0)
  ) {
    return(tibble::tibble(
      query_n = length(query_genes),
      target_n = length(target_genes),
      universe_n = length(universe_genes),
      overlap_n = length(overlap_genes),
      expected_overlap_n = NA_real_,
      odds_ratio = NA_real_,
      p_value = NA_real_,
      overlap_genes = paste(overlap_genes, collapse = ";")
    ))
  }

  ft <- fisher.test(mat, alternative = "greater")

  tibble::tibble(
    query_n = length(query_genes),
    target_n = length(target_genes),
    universe_n = length(universe_genes),
    overlap_n = length(overlap_genes),
    expected_overlap_n = length(query_genes) * length(target_genes) / length(universe_genes),
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value,
    overlap_genes = paste(overlap_genes, collapse = ";")
  )
}

sig_label <- function(padj) {
  dplyr::case_when(
    !is.na(padj) & padj < 0.001 ~ "***",
    !is.na(padj) & padj < 0.01 ~ "**",
    !is.na(padj) & padj < 0.05 ~ "*",
    !is.na(padj) & padj < 0.1 ~ ".",
    TRUE ~ ""
  )
}

# -----------------------------
# Read DGE files
# -----------------------------

dge_files <- list.files(base_dir, pattern = "^DGE_.*\\.csv$", full.names = TRUE)

if (length(dge_files) == 0) {
  stop("No DGE_*.csv files found in: ", base_dir)
}

dge_long <- purrr::map_dfr(dge_files, read_dge_file)

readr::write_csv(
  dge_long,
  file.path(out_dir, "all_DGE_long.csv")
)

readr::write_csv(
  dge_long %>%
    dplyr::distinct(contrast, gene_col_used, lfc_col_used, p_col_used) %>%
    dplyr::arrange(contrast),
  file.path(out_dir, "DGE_column_detection_check.csv")
)

# -----------------------------
# Build response gene sets
# -----------------------------
# Only these contrasts are used:
#   IL1_vs_Veh_3h
#   IL1_vs_Veh_6h
#   T5224_vs_Veh_3h
#   T5224_vs_Veh_6h
#
# T5224preIL1_vs_IL1 is deliberately excluded.

response_sets <- dge_long %>%
  dplyr::filter(sig, direction %in% c("up", "down")) %>%
  dplyr::mutate(
    response_type = dplyr::case_when(
      stringr::str_detect(contrast, "^IL1_vs_Veh") ~ "IL-1 vs Control",
      stringr::str_detect(contrast, "^T5224_vs_Veh") ~ "AP-1 inhibitor vs Control",
      TRUE ~ NA_character_
    ),
    timepoint = dplyr::case_when(
      stringr::str_detect(contrast, "3h$") ~ "3 h",
      stringr::str_detect(contrast, "6h$") ~ "6 h",
      TRUE ~ NA_character_
    ),
    response_set = paste(response_type, direction, timepoint, sep = " | ")
  ) %>%
  dplyr::filter(!is.na(response_type), !is.na(timepoint)) %>%
  dplyr::select(
    gene,
    response_type,
    timepoint,
    direction,
    response_set,
    contrast,
    log2FC,
    padj
  ) %>%
  dplyr::distinct()

response_sets_any_time <- response_sets %>%
  dplyr::group_by(gene, response_type, direction) %>%
  dplyr::summarise(
    timepoint = "Any time",
    response_set = paste(response_type[1], direction[1], "Any time", sep = " | "),
    min_padj = min(padj, na.rm = TRUE),
    max_abs_log2FC = max(abs(log2FC), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::select(gene, response_type, timepoint, direction, response_set)

response_sets_all <- dplyr::bind_rows(
  response_sets %>%
    dplyr::select(gene, response_type, timepoint, direction, response_set),
  response_sets_any_time
) %>%
  dplyr::distinct()

readr::write_csv(
  response_sets_all,
  file.path(out_dir, "IL1_vs_Control_and_Inhibitor_vs_Control_response_gene_sets.csv")
)

response_set_sizes <- response_sets_all %>%
  dplyr::count(response_type, timepoint, direction, response_set, name = "n_genes") %>%
  dplyr::arrange(response_type, timepoint, direction)

readr::write_csv(
  response_set_sizes,
  file.path(out_dir, "response_gene_set_sizes.csv")
)

# -----------------------------
# Load fibroblast Age/Exercise RDS
# -----------------------------

if (!file.exists(age_ex_rds) && file.exists("/mnt/data/All_DEG_Age_Ex.rds")) {
  age_ex_rds <- "/mnt/data/All_DEG_Age_Ex.rds"
}

if (!file.exists(age_ex_rds)) {
  stop("Cannot find All_DEG_Age_Ex.rds. Set age_ex_rds to the correct path.")
}

age_ex_raw <- readRDS(age_ex_rds)
age_ex_df <- flatten_rds_to_tibble(age_ex_raw)

gene_col <- guess_col(
  age_ex_df,
  candidates = c(
    "gene", "genes", "Gene", "symbol", "Symbol", "gene_symbol",
    "gene_name", "external_gene_name", "rowname_gene"
  ),
  fallback_regex = "gene|symbol|rowname"
)

celltype_col <- guess_col(
  age_ex_df,
  candidates = c(
    "CellType", "celltype", "cell_type", "cluster", "annotation",
    "ident", "cell", "rds_list_name"
  ),
  fallback_regex = "cell.*type|cluster|annotation|ident|rds_list_name"
)

age_lfc_col <- guess_col(
  age_ex_df,
  candidates = c(
    "Age_log2FC", "Age_logFC", "Age_avg_log2FC",
    "log2FC_Age", "logFC_Age", "avg_log2FC_Age",
    "log2FoldChange.Age", "log2FC.Age", "logFC.Age",
    "Age.log2FC", "Age.logFC"
  ),
  fallback_regex = "age.*(log2|logfc|fc|fold)|((log2|logfc|fc|fold).*)age",
  exclude_regex = "padj|fdr|p_val|pvalue|p\\.val|qval|stat|se"
)

age_padj_col <- guess_col(
  age_ex_df,
  candidates = c(
    "Age_padj", "Age_FDR", "Age_adj_p", "Age_adj.P.Val",
    "padj_Age", "FDR_Age", "adj.P.Val_Age",
    "padj.Age", "FDR.Age", "adj.P.Val.Age"
  ),
  fallback_regex = "age.*(padj|fdr|qval|adj)|((padj|fdr|qval|adj).*)age"
)

ex_lfc_col <- guess_col(
  age_ex_df,
  candidates = c(
    "Exercise_log2FC", "Exercise_logFC", "Exercise_avg_log2FC",
    "Ex_log2FC", "Ex_logFC", "Ex_avg_log2FC",
    "log2FC_Exercise", "logFC_Exercise", "avg_log2FC_Exercise",
    "log2FC_Ex", "logFC_Ex", "avg_log2FC_Ex",
    "log2FoldChange.Ex", "log2FC.Ex", "logFC.Ex",
    "Exercise.log2FC", "Exercise.logFC",
    "Ex.log2FC", "Ex.logFC"
  ),
  fallback_regex = "(exercise|(^|[_.])ex($|[_.]))(.*)(log2|logfc|fc|fold)|((log2|logfc|fc|fold).*)(exercise|(^|[_.])ex($|[_.]))",
  exclude_regex = "padj|fdr|p_val|pvalue|p\\.val|qval|stat|se"
)

ex_padj_col <- guess_col(
  age_ex_df,
  candidates = c(
    "Exercise_padj", "Exercise_FDR", "Exercise_adj_p", "Exercise_adj.P.Val",
    "Ex_padj", "Ex_FDR", "Ex_adj_p", "Ex_adj.P.Val",
    "padj_Exercise", "FDR_Exercise", "adj.P.Val_Exercise",
    "padj_Ex", "FDR_Ex", "adj.P.Val_Ex",
    "padj.Ex", "FDR.Ex", "adj.P.Val.Ex"
  ),
  fallback_regex = "(exercise|(^|[_.])ex($|[_.]))(.*)(padj|fdr|qval|adj)|((padj|fdr|qval|adj).*)(exercise|(^|[_.])ex($|[_.]))"
)

if (is.na(gene_col)) stop("Could not detect gene column in All_DEG_Age_Ex.rds")
if (is.na(celltype_col)) stop("Could not detect cell type column in All_DEG_Age_Ex.rds")
if (is.na(age_lfc_col)) stop("Could not detect age logFC column in All_DEG_Age_Ex.rds")
if (is.na(age_padj_col)) stop("Could not detect age padj/FDR column in All_DEG_Age_Ex.rds")
if (is.na(ex_lfc_col)) stop("Could not detect exercise logFC column in All_DEG_Age_Ex.rds")
if (is.na(ex_padj_col)) stop("Could not detect exercise padj/FDR column in All_DEG_Age_Ex.rds")

message("Detected Age/Exercise columns:")
message("gene_col: ", gene_col)
message("celltype_col: ", celltype_col)
message("age_lfc_col: ", age_lfc_col)
message("age_padj_col: ", age_padj_col)
message("ex_lfc_col: ", ex_lfc_col)
message("ex_padj_col: ", ex_padj_col)

fibro_age_ex <- age_ex_df %>%
  dplyr::transmute(
    gene = as.character(.data[[gene_col]]),
    CellType_raw = as.character(.data[[celltype_col]]),
    Age_log2FC = suppressWarnings(as.numeric(.data[[age_lfc_col]])),
    Age_padj = suppressWarnings(as.numeric(.data[[age_padj_col]])),
    Exercise_log2FC = suppressWarnings(as.numeric(.data[[ex_lfc_col]])),
    Exercise_padj = suppressWarnings(as.numeric(.data[[ex_padj_col]]))
  ) %>%
  dplyr::mutate(
    gene = stringr::str_replace(gene, "\\.\\d+$", ""),
    CellType = dplyr::case_when(
      stringr::str_to_lower(CellType_raw) %in% c("fibro", "fibroblast", "fibroblasts") ~ "Fibroblast",
      stringr::str_detect(stringr::str_to_lower(CellType_raw), "^fibro") ~ "Fibroblast",
      TRUE ~ CellType_raw
    )
  ) %>%
  dplyr::filter(CellType == "Fibroblast") %>%
  dplyr::filter(!is.na(gene), gene != "") %>%
  dplyr::mutate(
    Age_sig = !is.na(Age_padj) &
      Age_padj <= age_ex_fdr_cutoff &
      !is.na(Age_log2FC) &
      abs(Age_log2FC) >= age_ex_abs_log2fc_cutoff,

    Exercise_sig = !is.na(Exercise_padj) &
      Exercise_padj <= age_ex_fdr_cutoff &
      !is.na(Exercise_log2FC) &
      abs(Exercise_log2FC) >= age_ex_abs_log2fc_cutoff,

    Age_up = Age_sig & Age_log2FC > 0,
    Age_down = Age_sig & Age_log2FC < 0,
    Exercise_up = Exercise_sig & Exercise_log2FC > 0,
    Exercise_down = Exercise_sig & Exercise_log2FC < 0,

    Age_up_Exercise_down = Age_up & Exercise_down,
    Age_down_Exercise_up = Age_down & Exercise_up,

    Reversed_by_exercise = Age_sig &
      Exercise_sig &
      sign(Age_log2FC) == -sign(Exercise_log2FC),

    Exacerbated_by_exercise = Age_sig &
      Exercise_sig &
      sign(Age_log2FC) == sign(Exercise_log2FC)
  ) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

if (nrow(fibro_age_ex) == 0) {
  stop("No fibroblast rows found. Check CellType labels in All_DEG_Age_Ex.rds.")
}

readr::write_csv(
  fibro_age_ex,
  file.path(out_dir, "fibroblast_AgeEx_standardised.csv")
)

# -----------------------------
# Build fibroblast Age/Exercise sets
# -----------------------------

fibro_sets <- fibro_age_ex %>%
  dplyr::transmute(
    gene,
    Age_up,
    Age_down,
    Exercise_up,
    Exercise_down,
    Age_up_Exercise_down,
    Age_down_Exercise_up,
    Reversed_by_exercise,
    Exacerbated_by_exercise
  ) %>%
  tidyr::pivot_longer(
    cols = -gene,
    names_to = "fibroblast_set",
    values_to = "in_set"
  ) %>%
  dplyr::filter(in_set) %>%
  dplyr::select(gene, fibroblast_set) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    fibroblast_set = factor(
      fibroblast_set,
      levels = c(
        "Age_up",
        "Exercise_down",
        "Age_up_Exercise_down",
        "Reversed_by_exercise",
        "Age_down",
        "Exercise_up",
        "Age_down_Exercise_up",
        "Exacerbated_by_exercise"
      ),
      labels = c(
        "Age up",
        "Exercise down",
        "Age up + Exercise down",
        "Reversed by exercise",
        "Age down",
        "Exercise up",
        "Age down + Exercise up",
        "Exacerbated by exercise"
      )
    )
  )

readr::write_csv(
  fibro_sets,
  file.path(out_dir, "fibroblast_AgeEx_gene_sets.csv")
)

# -----------------------------
# Overlap tests
# -----------------------------

universe_genes <- intersect(
  unique(dge_long$gene),
  unique(fibro_age_ex$gene)
)

overlap_results <- tidyr::expand_grid(
  response_set = unique(response_sets_all$response_set),
  fibroblast_set = unique(fibro_sets$fibroblast_set)
) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    test = list({
      query_genes <- response_sets_all %>%
        dplyr::filter(response_set == .env$response_set) %>%
        dplyr::pull(gene) %>%
        unique()

      target_genes <- fibro_sets %>%
        dplyr::filter(fibroblast_set == .env$fibroblast_set) %>%
        dplyr::pull(gene) %>%
        unique()

      run_fisher_overlap(
        query_genes = query_genes,
        target_genes = target_genes,
        universe_genes = universe_genes
      )
    })
  ) %>%
  dplyr::ungroup() %>%
  tidyr::unnest(test) %>%
  tidyr::separate(
    response_set,
    into = c("response_type", "direction", "timepoint"),
    sep = " \\| ",
    remove = FALSE
  ) %>%
  dplyr::mutate(
    response_type = factor(
      response_type,
      levels = c("IL-1 vs Control", "AP-1 inhibitor vs Control")
    ),
    direction = factor(direction, levels = c("up", "down")),
    timepoint = factor(timepoint, levels = c("3 h", "6 h", "Any time")),
    fibroblast_set = factor(
      fibroblast_set,
      levels = c(
        "Age up",
        "Exercise down",
        "Age up + Exercise down",
        "Reversed by exercise",
        "Age down",
        "Exercise up",
        "Age down + Exercise up",
        "Exacerbated by exercise"
      )
    ),
    padj_global = p.adjust(p_value, method = "BH"),
    log2_odds_ratio = log2(odds_ratio),
    sig = sig_label(padj_global)
  ) %>%
  dplyr::arrange(padj_global, p_value)

readr::write_csv(
  overlap_results,
  file.path(out_dir, "IL1_vs_Control_Inhibitor_vs_Control_fibroblast_AgeEx_overlap_results.csv")
)

overlap_genes_long <- overlap_results %>%
  dplyr::filter(!is.na(overlap_genes), overlap_genes != "") %>%
  dplyr::select(
    response_type,
    direction,
    timepoint,
    response_set,
    fibroblast_set,
    overlap_genes
  ) %>%
  tidyr::separate_rows(overlap_genes, sep = ";") %>%
  dplyr::rename(gene = overlap_genes) %>%
  dplyr::left_join(fibro_age_ex, by = "gene") %>%
  dplyr::arrange(response_type, direction, timepoint, fibroblast_set, gene)

readr::write_csv(
  overlap_genes_long,
  file.path(out_dir, "IL1_vs_Control_Inhibitor_vs_Control_fibroblast_AgeEx_overlap_genes_long.csv")
)

# ============================================================
# Plots
# ============================================================

# -----------------------------
# Plot 1: main enrichment dotplot
# -----------------------------

p1 <- ggplot2::ggplot(
  overlap_results,
  ggplot2::aes(
    x = fibroblast_set,
    y = paste(direction, timepoint),
    size = overlap_n,
    fill = log2_odds_ratio
  )
) +
  ggplot2::geom_point(shape = 21, colour = "black", alpha = 0.9) +
  ggplot2::geom_text(
    ggplot2::aes(label = sig),
    size = 3.5,
    colour = "black",
    vjust = 0.35
  ) +
  ggplot2::facet_wrap(~ response_type, ncol = 1) +
  ggplot2::scale_fill_gradient2(
    low = "#377EB8",
    mid = "white",
    high = "#E41A1C",
    midpoint = 0,
    na.value = "grey85"
  ) +
  ggplot2::labs(
    x = "Fibroblast age/exercise set",
    y = "Response direction and timepoint",
    size = "Overlap genes",
    fill = "log2 odds ratio",
    title = "IL-1 vs Control and AP-1 inhibitor vs Control overlap with fibroblast age/exercise genes",
    subtitle = "Fisher test. *, **, *** indicate global BH-adjusted enrichment"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
    panel.grid.minor = ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(out_dir, "PLOT1_overlap_dotplot_IL1_vs_Control_Inhibitor_vs_Control_AgeEx.pdf"),
  p1,
  width = 12,
  height = 8
)

# -----------------------------
# Plot 2: model-relevant overlap counts
# -----------------------------

model_sets <- c(
  "Age up",
  "Exercise down",
  "Age up + Exercise down",
  "Reversed by exercise"
)

p2_df <- overlap_results %>%
  dplyr::filter(fibroblast_set %in% model_sets) %>%
  dplyr::mutate(
    fibroblast_set = factor(fibroblast_set, levels = model_sets)
  )

p2 <- ggplot2::ggplot(
  p2_df,
  ggplot2::aes(x = fibroblast_set, y = overlap_n, fill = response_type)
) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.8),
    colour = "black",
    linewidth = 0.2
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = overlap_n),
    position = ggplot2::position_dodge(width = 0.8),
    vjust = -0.25,
    size = 2.8
  ) +
  ggplot2::facet_grid(direction ~ timepoint) +
  ggplot2::scale_fill_manual(
    values = c(
      "IL-1 vs Control" = "#A81B1B",
      "AP-1 inhibitor vs Control" = "#A53D99"
    )
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Overlap genes",
    fill = NULL,
    title = "Overlap counts with model-relevant fibroblast age/exercise sets"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  file.path(out_dir, "PLOT2_overlap_counts_model_relevant_sets.pdf"),
  p2,
  width = 12,
  height = 7
)

# -----------------------------
# Plot 3: fibroblast Age vs Exercise scatter
# Highlight IL-1 vs Control and inhibitor vs Control responsive genes
# -----------------------------

highlight_genes <- response_sets_all %>%
  dplyr::filter(
    timepoint == "Any time",
    direction == "up",
    response_type %in% c("IL-1 vs Control", "AP-1 inhibitor vs Control")
  ) %>%
  dplyr::select(gene, response_type) %>%
  dplyr::distinct() %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    highlight_group = paste(sort(unique(response_type)), collapse = " + "),
    .groups = "drop"
  )

scatter_df <- fibro_age_ex %>%
  dplyr::left_join(highlight_genes, by = "gene") %>%
  dplyr::mutate(
    highlight_group = tidyr::replace_na(highlight_group, "Other fibroblast genes"),
    model_relevant = Age_up_Exercise_down | Reversed_by_exercise,
    label_gene = highlight_group != "Other fibroblast genes" & model_relevant
  )

label_df <- scatter_df %>%
  dplyr::filter(label_gene) %>%
  dplyr::arrange(Age_padj, Exercise_padj) %>%
  dplyr::slice_head(n = top_label_n)

p3 <- ggplot2::ggplot(
  scatter_df,
  ggplot2::aes(x = Age_log2FC, y = Exercise_log2FC)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  ggplot2::annotate(
    "rect",
    xmin = 0,
    xmax = Inf,
    ymin = -Inf,
    ymax = 0,
    alpha = 0.08,
    fill = "#A81B1B"
  ) +
  ggplot2::geom_point(
    data = scatter_df %>% dplyr::filter(highlight_group == "Other fibroblast genes"),
    colour = "grey82",
    alpha = 0.45,
    size = 1.1
  ) +
  ggplot2::geom_point(
    data = scatter_df %>% dplyr::filter(highlight_group != "Other fibroblast genes"),
    ggplot2::aes(colour = highlight_group),
    alpha = 0.85,
    size = 2
  ) +
  ggplot2::geom_text(
    data = label_df,
    ggplot2::aes(label = gene),
    size = 2.6,
    check_overlap = TRUE,
    vjust = -0.5
  ) +
  ggplot2::scale_colour_manual(
    values = c(
      "IL-1 vs Control" = "#A81B1B",
      "AP-1 inhibitor vs Control" = "#A53D99",
      "AP-1 inhibitor vs Control + IL-1 vs Control" = "#F99349"
    )
  ) +
  ggplot2::labs(
    x = "Fibroblast age effect log2FC",
    y = "Fibroblast exercise effect log2FC",
    colour = NULL,
    title = "Fibroblast age/exercise effects for IL-1 and AP-1 inhibitor responsive genes",
    subtitle = "Red shaded quadrant = age up and exercise down"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  file.path(out_dir, "PLOT3_fibroblast_Age_vs_Exercise_highlight_IL1_vs_Control_Inhibitor_vs_Control.pdf"),
  p3,
  width = 8,
  height = 7
)

# -----------------------------
# Plot 4: simple heatmap of top overlapping genes
# -----------------------------

top_overlap_genes <- overlap_genes_long %>%
  dplyr::filter(
    fibroblast_set %in% c(
      "Age up",
      "Exercise down",
      "Age up + Exercise down",
      "Reversed by exercise"
    ),
    response_type %in% c("IL-1 vs Control", "AP-1 inhibitor vs Control")
  ) %>%
  dplyr::count(gene, name = "n_overlap_categories") %>%
  dplyr::arrange(dplyr::desc(n_overlap_categories)) %>%
  dplyr::slice_head(n = top_heatmap_n) %>%
  dplyr::pull(gene)

heatmap_dge <- dge_long %>%
  dplyr::filter(
    gene %in% top_overlap_genes,
    contrast %in% c(
      "IL1_vs_Veh_3h",
      "IL1_vs_Veh_6h",
      "T5224_vs_Veh_3h",
      "T5224_vs_Veh_6h"
    )
  ) %>%
  dplyr::transmute(
    gene,
    effect = dplyr::case_when(
      contrast == "IL1_vs_Veh_3h" ~ "IL1 vs Veh 3 h",
      contrast == "IL1_vs_Veh_6h" ~ "IL1 vs Veh 6 h",
      contrast == "T5224_vs_Veh_3h" ~ "T5224 vs Veh 3 h",
      contrast == "T5224_vs_Veh_6h" ~ "T5224 vs Veh 6 h",
      TRUE ~ contrast
    ),
    log2FC,
    padj
  )

heatmap_fibro <- fibro_age_ex %>%
  dplyr::filter(gene %in% top_overlap_genes) %>%
  dplyr::select(gene, Age_log2FC, Age_padj, Exercise_log2FC, Exercise_padj)

heatmap_fibro_long <- dplyr::bind_rows(
  heatmap_fibro %>%
    dplyr::transmute(
      gene,
      effect = "Fibroblast Age effect",
      log2FC = Age_log2FC,
      padj = Age_padj
    ),
  heatmap_fibro %>%
    dplyr::transmute(
      gene,
      effect = "Fibroblast Exercise effect",
      log2FC = Exercise_log2FC,
      padj = Exercise_padj
    )
)

heatmap_df <- dplyr::bind_rows(
  heatmap_dge,
  heatmap_fibro_long
) %>%
  dplyr::mutate(
    gene = factor(gene, levels = rev(top_overlap_genes)),
    effect = factor(
      effect,
      levels = c(
        "IL1 vs Veh 3 h",
        "T5224 vs Veh 3 h",
        "IL1 vs Veh 6 h",
        "T5224 vs Veh 6 h",
        "Fibroblast Age effect",
        "Fibroblast Exercise effect"
      )
    ),
    sig_label = ifelse(!is.na(padj) & padj <= 0.05, "*", "")
  ) %>%
  dplyr::filter(!is.na(effect))

p4 <- ggplot2::ggplot(
  heatmap_df,
  ggplot2::aes(x = effect, y = gene, fill = log2FC)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
  ggplot2::geom_text(ggplot2::aes(label = sig_label), size = 2.4) +
  ggplot2::scale_fill_gradient2(
    low = "#377EB8",
    mid = "white",
    high = "#E41A1C",
    midpoint = 0,
    na.value = "grey90"
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    fill = "log2FC",
    title = "Top overlapping IL-1 and AP-1 inhibitor responsive genes",
    subtitle = "* = FDR <= 0.05"
  ) +
  ggplot2::theme_bw(base_size = 8.5) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
    panel.grid = ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(out_dir, "PLOT4_top_overlap_genes_heatmap.pdf"),
  p4,
  width = 11,
  height = max(6, length(top_overlap_genes) * 0.16)
)

# -----------------------------
# Print summary
# -----------------------------

message("Done.")
message("Output directory:")
message(out_dir)

message("\nMain tables:")
message(file.path(out_dir, "IL1_vs_Control_and_Inhibitor_vs_Control_response_gene_sets.csv"))
message(file.path(out_dir, "fibroblast_AgeEx_gene_sets.csv"))
message(file.path(out_dir, "IL1_vs_Control_Inhibitor_vs_Control_fibroblast_AgeEx_overlap_results.csv"))
message(file.path(out_dir, "IL1_vs_Control_Inhibitor_vs_Control_fibroblast_AgeEx_overlap_genes_long.csv"))

message("\nMain plots:")
message(file.path(out_dir, "PLOT1_overlap_dotplot_IL1_vs_Control_Inhibitor_vs_Control_AgeEx.pdf"))
message(file.path(out_dir, "PLOT2_overlap_counts_model_relevant_sets.pdf"))
message(file.path(out_dir, "PLOT3_fibroblast_Age_vs_Exercise_highlight_IL1_vs_Control_Inhibitor_vs_Control.pdf"))
message(file.path(out_dir, "PLOT4_top_overlap_genes_heatmap.pdf"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(forcats)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(enrichplot)
})

# ============================================================
# GO enrichment for:
#   1. IL-1 responsive genes that increase with age
#   2. IL-1 responsive genes that decrease with age
#   3. AP-1 inhibitor responsive genes that change with exercise
#
# Uses:
#   IL-1 genes        = DGE_IL1_vs_Veh_3h/6h
#   Inhibitor genes  = DGE_T5224_vs_Veh_3h/6h
#   Age/Exercise     = fibroblast rows from All_DEG_Age_Ex.rds
#
# Important:
#   This does NOT use T5224preIL1_vs_IL1.
# ============================================================

base_dir <- publication_path("PAPER_OUTPUT_ROOT", "IL1_RNASeq")

age_ex_rds <- file.path(base_dir, "All_DEG_Age_Ex.rds")

if (!file.exists(age_ex_rds) && file.exists("/mnt/data/All_DEG_Age_Ex.rds")) {
  age_ex_rds <- "/mnt/data/All_DEG_Age_Ex.rds"
}

out_dir <- file.path(base_dir, "GO_IL1_Inhibitor_Fibroblast_AgeEx")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dge_fdr_cutoff <- 0.05
dge_abs_log2fc_cutoff <- 0

age_ex_fdr_cutoff <- 0.05
age_ex_abs_log2fc_cutoff <- 0

go_ontology <- "BP"
go_p_cutoff <- 0.05
go_q_cutoff <- 0.2
min_genes_for_go <- 5
top_n_terms <- 20

# -----------------------------
# Helper functions
# -----------------------------

guess_col <- function(df, candidates, fallback_regex = NULL, exclude_regex = NULL) {
  nms <- names(df)
  nms_clean <- tolower(nms)
  candidates_clean <- tolower(candidates)

  exact_hit <- which(nms_clean %in% candidates_clean)

  if (!is.null(exclude_regex) && length(exact_hit) > 0) {
    exact_hit <- exact_hit[!stringr::str_detect(nms_clean[exact_hit], exclude_regex)]
  }

  if (length(exact_hit) > 0) {
    return(nms[exact_hit[1]])
  }

  if (!is.null(fallback_regex)) {
    regex_hit <- which(stringr::str_detect(nms_clean, fallback_regex))

    if (!is.null(exclude_regex) && length(regex_hit) > 0) {
      regex_hit <- regex_hit[!stringr::str_detect(nms_clean[regex_hit], exclude_regex)]
    }

    if (length(regex_hit) > 0) {
      return(nms[regex_hit[1]])
    }
  }

  NA_character_
}

flatten_rds_to_tibble <- function(obj) {
  if (is.data.frame(obj)) {
    df <- as.data.frame(obj)

    if (
      !"gene" %in% tolower(names(df)) &&
        !all(rownames(df) == as.character(seq_len(nrow(df))))
    ) {
      df <- tibble::rownames_to_column(df, var = "rowname_gene")
    }

    return(tibble::as_tibble(df))
  }

  if (is.list(obj)) {
    df_list <- purrr::imap(obj, function(x, nm) {
      if (!is.data.frame(x)) {
        return(NULL)
      }

      df <- as.data.frame(x)

      if (
        !"gene" %in% tolower(names(df)) &&
          !all(rownames(df) == as.character(seq_len(nrow(df))))
      ) {
        df <- tibble::rownames_to_column(df, var = "rowname_gene")
      }

      df$rds_list_name <- nm
      tibble::as_tibble(df)
    })

    df_list <- purrr::compact(df_list)

    if (length(df_list) == 0) {
      stop("RDS is a list, but no data.frames were found inside it.")
    }

    return(dplyr::bind_rows(df_list))
  }

  stop("Unsupported RDS object type: ", paste(class(obj), collapse = ", "))
}

read_dge_file <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)

  gene_col <- guess_col(
    df,
    candidates = c(
      "gene", "genes", "symbol", "gene_symbol", "gene_name",
      "external_gene_name", "feature", "features", "X", "...1"
    ),
    fallback_regex = "gene|symbol|feature"
  )

  p_col <- guess_col(
    df,
    candidates = c(
      "padj", "p_val_adj", "p_val_adj.", "FDR", "adj.P.Val",
      "qval", "q_value", "qvalue", "fdr_bh"
    ),
    fallback_regex = "padj|adj.*p|p.*adj|fdr|qval|q_value"
  )

  lfc_col <- guess_col(
    df,
    candidates = c(
      "log2FoldChange", "logFC", "avg_log2FC", "avg_logFC",
      "log2FC", "estimate", "coef"
    ),
    fallback_regex = "log2.*fold|logfc|log2fc|avg_log",
    exclude_regex = "padj|fdr|p_val|pvalue|qval"
  )

  if (is.na(gene_col)) stop("Could not detect gene column in: ", basename(path))
  if (is.na(p_col)) stop("Could not detect FDR/padj column in: ", basename(path))
  if (is.na(lfc_col)) stop("Could not detect logFC column in: ", basename(path))

  contrast <- basename(path) %>%
    stringr::str_remove("^DGE_") %>%
    stringr::str_remove("\\.csv$")

  df %>%
    dplyr::transmute(
      gene = as.character(.data[[gene_col]]),
      contrast = contrast,
      log2FC = suppressWarnings(as.numeric(.data[[lfc_col]])),
      padj = suppressWarnings(as.numeric(.data[[p_col]])),
      gene_col_used = gene_col,
      lfc_col_used = lfc_col,
      p_col_used = p_col
    ) %>%
    dplyr::mutate(
      gene = stringr::str_replace(gene, "\\.\\d+$", ""),
      sig = !is.na(padj) &
        padj <= dge_fdr_cutoff &
        !is.na(log2FC) &
        abs(log2FC) >= dge_abs_log2fc_cutoff,
      direction = dplyr::case_when(
        sig & log2FC > 0 ~ "up",
        sig & log2FC < 0 ~ "down",
        TRUE ~ "ns"
      )
    ) %>%
    dplyr::filter(!is.na(gene), gene != "", !is.na(log2FC), !is.na(padj)) %>%
    dplyr::distinct(gene, contrast, .keep_all = TRUE)
}

symbols_to_entrez <- function(symbols) {
  symbols <- unique(symbols)
  symbols <- symbols[!is.na(symbols) & symbols != ""]

  mapped <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = symbols,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )

  mapped <- mapped[!is.na(mapped)]
  unique(unname(mapped))
}

run_go <- function(gene_symbols, universe_symbols, set_name) {
  gene_entrez <- symbols_to_entrez(gene_symbols)
  universe_entrez <- symbols_to_entrez(universe_symbols)

  gene_entrez <- intersect(gene_entrez, universe_entrez)

  if (length(gene_entrez) < min_genes_for_go) {
    message("Skipping ", set_name, ": fewer than ", min_genes_for_go, " mapped genes.")
    return(NULL)
  }

  ego <- clusterProfiler::enrichGO(
    gene = gene_entrez,
    universe = universe_entrez,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = go_ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = go_p_cutoff,
    qvalueCutoff = go_q_cutoff,
    readable = TRUE
  )

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    message("No enriched GO terms for: ", set_name)
    return(NULL)
  }

  ego
}

plot_go_bar <- function(go_df, set_name, output_file, top_n = 20) {
  plot_df <- go_df %>%
    dplyr::arrange(p.adjust, pvalue) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(
      Description = stringr::str_wrap(Description, width = 45),
      Description = forcats::fct_reorder(Description, -log10(p.adjust))
    )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = -log10(p.adjust),
      y = Description,
      size = Count,
      fill = GeneRatio_numeric
    )
  ) +
    ggplot2::geom_point(shape = 21, colour = "black", alpha = 0.9) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "#A81B1B",
      na.value = "grey85"
    ) +
    ggplot2::labs(
      x = "-log10(BH adjusted P)",
      y = NULL,
      size = "Genes",
      fill = "Gene ratio",
      title = set_name
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank()
    )

  ggplot2::ggsave(
    output_file,
    p,
    width = 8,
    height = max(4.5, nrow(plot_df) * 0.28)
  )
}

parse_gene_ratio <- function(x) {
  purrr::map_dbl(x, function(z) {
    parts <- stringr::str_split(z, "/", simplify = TRUE)
    as.numeric(parts[1]) / as.numeric(parts[2])
  })
}

# -----------------------------
# Load DGE files
# -----------------------------

dge_files <- list.files(base_dir, pattern = "^DGE_.*\\.csv$", full.names = TRUE)

if (length(dge_files) == 0) {
  stop("No DGE_*.csv files found in: ", base_dir)
}

dge_long <- purrr::map_dfr(dge_files, read_dge_file)

readr::write_csv(
  dge_long,
  file.path(out_dir, "all_DGE_long.csv")
)

readr::write_csv(
  dge_long %>%
    dplyr::distinct(contrast, gene_col_used, lfc_col_used, p_col_used) %>%
    dplyr::arrange(contrast),
  file.path(out_dir, "DGE_column_detection_check.csv")
)

# -----------------------------
# Build IL-1 and inhibitor response gene sets
# -----------------------------

il1_response <- dge_long %>%
  dplyr::filter(
    contrast %in% c("IL1_vs_Veh_3h", "IL1_vs_Veh_6h"),
    sig,
    direction %in% c("up", "down")
  )

inhibitor_response <- dge_long %>%
  dplyr::filter(
    contrast %in% c("T5224_vs_Veh_3h", "T5224_vs_Veh_6h"),
    sig,
    direction %in% c("up", "down")
  )

il1_any <- unique(il1_response$gene)
il1_up <- il1_response %>% dplyr::filter(direction == "up") %>% dplyr::pull(gene) %>% unique()
il1_down <- il1_response %>% dplyr::filter(direction == "down") %>% dplyr::pull(gene) %>% unique()

inhibitor_any <- unique(inhibitor_response$gene)
inhibitor_up <- inhibitor_response %>% dplyr::filter(direction == "up") %>% dplyr::pull(gene) %>% unique()
inhibitor_down <- inhibitor_response %>% dplyr::filter(direction == "down") %>% dplyr::pull(gene) %>% unique()

# -----------------------------
# Load fibroblast age/exercise table
# -----------------------------

if (!file.exists(age_ex_rds)) {
  stop("Cannot find All_DEG_Age_Ex.rds. Set age_ex_rds to the correct path.")
}

age_ex_raw <- readRDS(age_ex_rds)
age_ex_df <- flatten_rds_to_tibble(age_ex_raw)

gene_col <- guess_col(
  age_ex_df,
  candidates = c(
    "gene", "genes", "Gene", "symbol", "Symbol", "gene_symbol",
    "gene_name", "external_gene_name", "rowname_gene"
  ),
  fallback_regex = "gene|symbol|rowname"
)

celltype_col <- guess_col(
  age_ex_df,
  candidates = c(
    "CellType", "celltype", "cell_type", "cluster", "annotation",
    "ident", "cell", "rds_list_name"
  ),
  fallback_regex = "cell.*type|cluster|annotation|ident|rds_list_name"
)

age_lfc_col <- guess_col(
  age_ex_df,
  candidates = c(
    "Age_log2FC", "Age_logFC", "Age_avg_log2FC",
    "log2FC_Age", "logFC_Age", "avg_log2FC_Age",
    "log2FoldChange.Age", "log2FC.Age", "logFC.Age",
    "Age.log2FC", "Age.logFC"
  ),
  fallback_regex = "age.*(log2|logfc|fc|fold)|((log2|logfc|fc|fold).*)age",
  exclude_regex = "padj|fdr|p_val|pvalue|p\\.val|qval|stat|se"
)

age_padj_col <- guess_col(
  age_ex_df,
  candidates = c(
    "Age_padj", "Age_FDR", "Age_adj_p", "Age_adj.P.Val",
    "padj_Age", "FDR_Age", "adj.P.Val_Age",
    "padj.Age", "FDR.Age", "adj.P.Val.Age"
  ),
  fallback_regex = "age.*(padj|fdr|qval|adj)|((padj|fdr|qval|adj).*)age"
)

ex_lfc_col <- guess_col(
  age_ex_df,
  candidates = c(
    "Exercise_log2FC", "Exercise_logFC", "Exercise_avg_log2FC",
    "Ex_log2FC", "Ex_logFC", "Ex_avg_log2FC",
    "log2FC_Exercise", "logFC_Exercise", "avg_log2FC_Exercise",
    "log2FC_Ex", "logFC_Ex", "avg_log2FC_Ex",
    "log2FoldChange.Ex", "log2FC.Ex", "logFC.Ex",
    "Exercise.log2FC", "Exercise.logFC",
    "Ex.log2FC", "Ex.logFC"
  ),
  fallback_regex = "(exercise|(^|[_.])ex($|[_.]))(.*)(log2|logfc|fc|fold)|((log2|logfc|fc|fold).*)(exercise|(^|[_.])ex($|[_.]))",
  exclude_regex = "padj|fdr|p_val|pvalue|p\\.val|qval|stat|se"
)

ex_padj_col <- guess_col(
  age_ex_df,
  candidates = c(
    "Exercise_padj", "Exercise_FDR", "Exercise_adj_p", "Exercise_adj.P.Val",
    "Ex_padj", "Ex_FDR", "Ex_adj_p", "Ex_adj.P.Val",
    "padj_Exercise", "FDR_Exercise", "adj.P.Val_Exercise",
    "padj_Ex", "FDR_Ex", "adj.P.Val_Ex",
    "padj.Ex", "FDR.Ex", "adj.P.Val.Ex"
  ),
  fallback_regex = "(exercise|(^|[_.])ex($|[_.]))(.*)(padj|fdr|qval|adj)|((padj|fdr|qval|adj).*)(exercise|(^|[_.])ex($|[_.]))"
)

if (is.na(gene_col)) stop("Could not detect gene column in All_DEG_Age_Ex.rds")
if (is.na(celltype_col)) stop("Could not detect cell type column in All_DEG_Age_Ex.rds")
if (is.na(age_lfc_col)) stop("Could not detect age logFC column in All_DEG_Age_Ex.rds")
if (is.na(age_padj_col)) stop("Could not detect age padj/FDR column in All_DEG_Age_Ex.rds")
if (is.na(ex_lfc_col)) stop("Could not detect exercise logFC column in All_DEG_Age_Ex.rds")
if (is.na(ex_padj_col)) stop("Could not detect exercise padj/FDR column in All_DEG_Age_Ex.rds")

message("Detected Age/Exercise columns:")
message("gene_col: ", gene_col)
message("celltype_col: ", celltype_col)
message("age_lfc_col: ", age_lfc_col)
message("age_padj_col: ", age_padj_col)
message("ex_lfc_col: ", ex_lfc_col)
message("ex_padj_col: ", ex_padj_col)

fibro_age_ex <- age_ex_df %>%
  dplyr::transmute(
    gene = as.character(.data[[gene_col]]),
    CellType_raw = as.character(.data[[celltype_col]]),
    Age_log2FC = suppressWarnings(as.numeric(.data[[age_lfc_col]])),
    Age_padj = suppressWarnings(as.numeric(.data[[age_padj_col]])),
    Exercise_log2FC = suppressWarnings(as.numeric(.data[[ex_lfc_col]])),
    Exercise_padj = suppressWarnings(as.numeric(.data[[ex_padj_col]]))
  ) %>%
  dplyr::mutate(
    gene = stringr::str_replace(gene, "\\.\\d+$", ""),
    CellType = dplyr::case_when(
      stringr::str_to_lower(CellType_raw) %in% c("fibro", "fibroblast", "fibroblasts") ~ "Fibroblast",
      stringr::str_detect(stringr::str_to_lower(CellType_raw), "^fibro") ~ "Fibroblast",
      TRUE ~ CellType_raw
    )
  ) %>%
  dplyr::filter(CellType == "Fibroblast") %>%
  dplyr::filter(!is.na(gene), gene != "") %>%
  dplyr::mutate(
    Age_sig = !is.na(Age_padj) &
      Age_padj <= age_ex_fdr_cutoff &
      !is.na(Age_log2FC) &
      abs(Age_log2FC) >= age_ex_abs_log2fc_cutoff,

    Exercise_sig = !is.na(Exercise_padj) &
      Exercise_padj <= age_ex_fdr_cutoff &
      !is.na(Exercise_log2FC) &
      abs(Exercise_log2FC) >= age_ex_abs_log2fc_cutoff,

    Age_up = Age_sig & Age_log2FC > 0,
    Age_down = Age_sig & Age_log2FC < 0,
    Exercise_up = Exercise_sig & Exercise_log2FC > 0,
    Exercise_down = Exercise_sig & Exercise_log2FC < 0
  ) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

if (nrow(fibro_age_ex) == 0) {
  stop("No fibroblast rows found. Check CellType labels in All_DEG_Age_Ex.rds.")
}

readr::write_csv(
  fibro_age_ex,
  file.path(out_dir, "fibroblast_AgeEx_standardised.csv")
)

# -----------------------------
# Define overlap gene sets for GO
# -----------------------------

age_up_genes <- fibro_age_ex %>% dplyr::filter(Age_up) %>% dplyr::pull(gene) %>% unique()
age_down_genes <- fibro_age_ex %>% dplyr::filter(Age_down) %>% dplyr::pull(gene) %>% unique()

exercise_changed_genes <- fibro_age_ex %>% dplyr::filter(Exercise_sig) %>% dplyr::pull(gene) %>% unique()
exercise_up_genes <- fibro_age_ex %>% dplyr::filter(Exercise_up) %>% dplyr::pull(gene) %>% unique()
exercise_down_genes <- fibro_age_ex %>% dplyr::filter(Exercise_down) %>% dplyr::pull(gene) %>% unique()

go_gene_sets <- list(
  IL1_any_Age_up = intersect(il1_any, age_up_genes),
  IL1_up_Age_up = intersect(il1_up, age_up_genes),
  IL1_down_Age_up = intersect(il1_down, age_up_genes),

  IL1_any_Age_down = intersect(il1_any, age_down_genes),
  IL1_up_Age_down = intersect(il1_up, age_down_genes),
  IL1_down_Age_down = intersect(il1_down, age_down_genes),

  Inhibitor_any_Exercise_changed = intersect(inhibitor_any, exercise_changed_genes),
  Inhibitor_up_Exercise_changed = intersect(inhibitor_up, exercise_changed_genes),
  Inhibitor_down_Exercise_changed = intersect(inhibitor_down, exercise_changed_genes),

  Inhibitor_any_Exercise_up = intersect(inhibitor_any, exercise_up_genes),
  Inhibitor_any_Exercise_down = intersect(inhibitor_any, exercise_down_genes),
  Inhibitor_up_Exercise_up = intersect(inhibitor_up, exercise_up_genes),
  Inhibitor_up_Exercise_down = intersect(inhibitor_up, exercise_down_genes),
  Inhibitor_down_Exercise_up = intersect(inhibitor_down, exercise_up_genes),
  Inhibitor_down_Exercise_down = intersect(inhibitor_down, exercise_down_genes)
)

go_input_table <- purrr::imap_dfr(
  go_gene_sets,
  function(genes, set_name) {
    tibble::tibble(
      set_name = set_name,
      gene = sort(unique(genes)),
      n_genes = length(unique(genes))
    )
  }
)

readr::write_csv(
  go_input_table,
  file.path(out_dir, "GO_input_gene_sets_long.csv")
)

go_input_sizes <- go_input_table %>%
  dplyr::distinct(set_name, n_genes) %>%
  dplyr::arrange(dplyr::desc(n_genes))

readr::write_csv(
  go_input_sizes,
  file.path(out_dir, "GO_input_gene_set_sizes.csv")
)

# Background universe:
# genes measured in both the IL-1/inhibitor DGE tables and fibroblast age/exercise table
universe_symbols <- intersect(
  unique(dge_long$gene),
  unique(fibro_age_ex$gene)
)

readr::write_csv(
  tibble::tibble(gene = sort(universe_symbols)),
  file.path(out_dir, "GO_background_universe_genes.csv")
)

# -----------------------------
# Run GO enrichment
# -----------------------------

go_results <- list()

for (set_name in names(go_gene_sets)) {
  message("Running GO for: ", set_name)

  ego <- run_go(
    gene_symbols = go_gene_sets[[set_name]],
    universe_symbols = universe_symbols,
    set_name = set_name
  )

  go_results[[set_name]] <- ego

  if (!is.null(ego)) {
    go_df <- as.data.frame(ego) %>%
      dplyr::mutate(
        set_name = set_name,
        GeneRatio_numeric = parse_gene_ratio(GeneRatio),
        .before = 1
      )

    readr::write_csv(
      go_df,
      file.path(out_dir, paste0("GO_", set_name, ".csv"))
    )

    plot_go_bar(
      go_df = go_df,
      set_name = set_name,
      output_file = file.path(out_dir, paste0("PLOT_GO_", set_name, ".pdf")),
      top_n = top_n_terms
    )
  }
}

go_results_long <- purrr::imap_dfr(
  go_results,
  function(ego, set_name) {
    if (is.null(ego)) {
      return(tibble::tibble())
    }

    as.data.frame(ego) %>%
      dplyr::mutate(
        set_name = set_name,
        GeneRatio_numeric = parse_gene_ratio(GeneRatio),
        .before = 1
      )
  }
) %>%
  dplyr::arrange(set_name, p.adjust, pvalue)

readr::write_csv(
  go_results_long,
  file.path(out_dir, "GO_enrichment_results_long.csv")
)

# -----------------------------
# Combined GO plot for the most useful sets
# -----------------------------

focused_sets <- c(
  "IL1_any_Age_up",
  "IL1_up_Age_up",
  "IL1_any_Age_down",
  "IL1_down_Age_down",
  "Inhibitor_any_Exercise_changed",
  "Inhibitor_up_Exercise_changed",
  "Inhibitor_down_Exercise_changed",
  "Inhibitor_any_Exercise_up",
  "Inhibitor_any_Exercise_down"
)

focused_go_df <- go_results_long %>%
  dplyr::filter(set_name %in% focused_sets) %>%
  dplyr::group_by(set_name) %>%
  dplyr::arrange(p.adjust, pvalue, .by_group = TRUE) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    set_name = factor(set_name, levels = focused_sets),
    Description_wrapped = stringr::str_wrap(Description, width = 42),
    Description_wrapped = forcats::fct_reorder(Description_wrapped, -log10(p.adjust))
  )

if (nrow(focused_go_df) > 0) {
  p_combined <- ggplot2::ggplot(
    focused_go_df,
    ggplot2::aes(
      x = set_name,
      y = Description_wrapped,
      size = Count,
      fill = -log10(p.adjust)
    )
  ) +
    ggplot2::geom_point(shape = 21, colour = "black", alpha = 0.9) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "#A81B1B",
      na.value = "grey85"
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      size = "Genes",
      fill = "-log10(FDR)",
      title = "GO enrichment of IL-1 and AP-1 inhibitor overlap gene sets",
      subtitle = "Universe = genes measured in both DGE and fibroblast age/exercise datasets"
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )

  ggplot2::ggsave(
    file.path(out_dir, "PLOT_combined_GO_focused_sets.pdf"),
    p_combined,
    width = 12,
    height = max(6, length(unique(focused_go_df$Description_wrapped)) * 0.22)
  )
}

# -----------------------------
# Optional compareCluster plot
# -----------------------------

compare_gene_sets <- go_gene_sets[focused_sets]
compare_gene_sets <- compare_gene_sets[purrr::map_int(compare_gene_sets, length) >= min_genes_for_go]

if (length(compare_gene_sets) >= 2) {
  compare_entrez_sets <- purrr::map(compare_gene_sets, symbols_to_entrez)
  compare_entrez_sets <- purrr::map(compare_entrez_sets, ~ intersect(.x, symbols_to_entrez(universe_symbols)))

  cc <- clusterProfiler::compareCluster(
    geneCluster = compare_entrez_sets,
    fun = "enrichGO",
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = go_ontology,
    universe = symbols_to_entrez(universe_symbols),
    pAdjustMethod = "BH",
    pvalueCutoff = go_p_cutoff,
    qvalueCutoff = go_q_cutoff,
    readable = TRUE
  )

  cc_df <- as.data.frame(cc)

  readr::write_csv(
    cc_df,
    file.path(out_dir, "GO_compareCluster_focused_sets.csv")
  )

  if (nrow(cc_df) > 0) {
    p_cc <- enrichplot::dotplot(
      cc,
      showCategory = 10,
      font.size = 8
    ) +
      ggplot2::labs(
        title = "GO compareCluster: IL-1 and inhibitor overlap gene sets"
      ) +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)
      )

    ggplot2::ggsave(
      file.path(out_dir, "PLOT_compareCluster_GO_focused_sets.pdf"),
      p_cc,
      width = 12,
      height = 8
    )
  }
}

# -----------------------------
# Print useful summary
# -----------------------------

message("Done.")
message("Output directory:")
message(out_dir)

message("\nGene set sizes:")
print(go_input_sizes)

message("\nMain outputs:")
message(file.path(out_dir, "GO_input_gene_sets_long.csv"))
message(file.path(out_dir, "GO_input_gene_set_sizes.csv"))
message(file.path(out_dir, "GO_background_universe_genes.csv"))
message(file.path(out_dir, "GO_enrichment_results_long.csv"))
message(file.path(out_dir, "PLOT_combined_GO_focused_sets.pdf"))
message(file.path(out_dir, "PLOT_compareCluster_GO_focused_sets.pdf"))

message("\nTop GO terms:")

# ============================================================
# Replace from:
#   # Run GO enrichment
# onwards with this safer version
# ============================================================

parse_gene_ratio <- function(x) {
  purrr::map_dbl(x, function(z) {
    if (is.na(z) || !stringr::str_detect(z, "/")) {
      return(NA_real_)
    }

    parts <- stringr::str_split(z, "/", simplify = TRUE)

    if (ncol(parts) < 2) {
      return(NA_real_)
    }

    as.numeric(parts[1]) / as.numeric(parts[2])
  })
}

standardise_go_df <- function(ego, set_name) {
  if (is.null(ego)) {
    return(tibble::tibble())
  }

  go_df <- as.data.frame(ego)

  if (nrow(go_df) == 0) {
    return(tibble::tibble())
  }

  required_cols <- c(
    "ID",
    "Description",
    "GeneRatio",
    "BgRatio",
    "pvalue",
    "p.adjust",
    "qvalue",
    "geneID",
    "Count"
  )

  missing_cols <- setdiff(required_cols, names(go_df))

  for (mc in missing_cols) {
    go_df[[mc]] <- NA
  }

  go_df %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      set_name = set_name,
      GeneRatio_numeric = parse_gene_ratio(GeneRatio),
      p_value = as.numeric(pvalue),
      padj = as.numeric(`p.adjust`),
      q_value = as.numeric(qvalue),
      Count = as.numeric(Count),
      .before = 1
    ) %>%
    dplyr::select(
      set_name,
      ID,
      Description,
      GeneRatio,
      GeneRatio_numeric,
      BgRatio,
      p_value,
      padj,
      q_value,
      geneID,
      Count,
      dplyr::everything()
    )
}

plot_go_bar_safe <- function(go_df, set_name, output_file, top_n = 20) {
  if (is.null(go_df) || nrow(go_df) == 0) {
    message("No GO terms to plot for: ", set_name)
    return(NULL)
  }

  plot_df <- go_df %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::arrange(padj, p_value) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(
      Description = stringr::str_wrap(Description, width = 45),
      Description = forcats::fct_reorder(Description, -log10(padj))
    )

  if (nrow(plot_df) == 0) {
    message("No plottable GO terms for: ", set_name)
    return(NULL)
  }

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = -log10(padj),
      y = Description,
      size = Count,
      fill = GeneRatio_numeric
    )
  ) +
    ggplot2::geom_point(shape = 21, colour = "black", alpha = 0.9) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "#A81B1B",
      na.value = "grey85"
    ) +
    ggplot2::labs(
      x = "-log10(BH adjusted P)",
      y = NULL,
      size = "Genes",
      fill = "Gene ratio",
      title = set_name
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank()
    )

  ggplot2::ggsave(
    output_file,
    p,
    width = 8,
    height = max(4.5, nrow(plot_df) * 0.28)
  )

  p
}

# -----------------------------
# Run GO enrichment safely
# -----------------------------

go_results <- list()
go_results_df_list <- list()

for (set_name in names(go_gene_sets)) {
  message("Running GO for: ", set_name)

  genes_this_set <- unique(go_gene_sets[[set_name]])
  genes_this_set <- genes_this_set[!is.na(genes_this_set) & genes_this_set != ""]

  message("  Input genes: ", length(genes_this_set))

  if (length(genes_this_set) < min_genes_for_go) {
    message("  Skipping ", set_name, ": fewer than ", min_genes_for_go, " genes before ID mapping.")
    next
  }

  ego <- run_go(
    gene_symbols = genes_this_set,
    universe_symbols = universe_symbols,
    set_name = set_name
  )

  go_results[[set_name]] <- ego

  go_df <- standardise_go_df(ego, set_name)

  if (nrow(go_df) == 0) {
    message("  No enriched GO terms for: ", set_name)
    next
  }

  go_results_df_list[[set_name]] <- go_df

  readr::write_csv(
    go_df,
    file.path(out_dir, paste0("GO_", set_name, ".csv"))
  )

  plot_go_bar_safe(
    go_df = go_df,
    set_name = set_name,
    output_file = file.path(out_dir, paste0("PLOT_GO_", set_name, ".pdf")),
    top_n = top_n_terms
  )
}

go_results_long <- dplyr::bind_rows(go_results_df_list)

if (nrow(go_results_long) == 0) {
  warning("No GO enrichment results were returned for any gene set.")

  go_results_long <- tibble::tibble(
    set_name = character(),
    ID = character(),
    Description = character(),
    GeneRatio = character(),
    GeneRatio_numeric = numeric(),
    BgRatio = character(),
    p_value = numeric(),
    padj = numeric(),
    q_value = numeric(),
    geneID = character(),
    Count = numeric()
  )
} else {
  go_results_long <- go_results_long %>%
    dplyr::arrange(set_name, padj, p_value)
}

readr::write_csv(
  go_results_long,
  file.path(out_dir, "GO_enrichment_results_long.csv")
)

# -----------------------------
# Combined GO plot for focused sets
# -----------------------------

focused_sets <- c(
  "IL1_any_Age_up",
  "IL1_up_Age_up",
  "IL1_down_Age_up",
  "IL1_any_Age_down",
  "IL1_up_Age_down",
  "IL1_down_Age_down",
  "Inhibitor_any_Exercise_changed",
  "Inhibitor_up_Exercise_changed",
  "Inhibitor_down_Exercise_changed",
  "Inhibitor_any_Exercise_up",
  "Inhibitor_any_Exercise_down"
)

focused_go_df <- go_results_long %>%
  dplyr::filter(set_name %in% focused_sets) %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::group_by(set_name) %>%
  dplyr::arrange(padj, p_value, .by_group = TRUE) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    set_name = factor(set_name, levels = focused_sets),
    Description_wrapped = stringr::str_wrap(Description, width = 42),
    Description_wrapped = forcats::fct_reorder(Description_wrapped, -log10(padj))
  )

if (nrow(focused_go_df) > 0) {
  p_combined <- ggplot2::ggplot(
    focused_go_df,
    ggplot2::aes(
      x = set_name,
      y = Description_wrapped,
      size = Count,
      fill = -log10(padj)
    )
  ) +
    ggplot2::geom_point(shape = 21, colour = "black", alpha = 0.9) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "#A81B1B",
      na.value = "grey85"
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      size = "Genes",
      fill = "-log10(FDR)",
      title = "GO enrichment of IL-1 and AP-1 inhibitor overlap gene sets",
      subtitle = "Universe = genes measured in both DGE and fibroblast age/exercise datasets"
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )

  ggplot2::ggsave(
    file.path(out_dir, "PLOT_combined_GO_focused_sets.pdf"),
    p_combined,
    width = 12,
    height = max(6, length(unique(focused_go_df$Description_wrapped)) * 0.22)
  )
} else {
  message("No focused GO terms available for combined plot.")
}

# -----------------------------
# Optional compareCluster plot
# -----------------------------

compare_gene_sets <- go_gene_sets[intersect(focused_sets, names(go_gene_sets))]
compare_gene_sets <- compare_gene_sets[purrr::map_int(compare_gene_sets, length) >= min_genes_for_go]

if (length(compare_gene_sets) >= 2) {
  universe_entrez <- symbols_to_entrez(universe_symbols)

  compare_entrez_sets <- purrr::map(compare_gene_sets, symbols_to_entrez)
  compare_entrez_sets <- purrr::map(compare_entrez_sets, ~ intersect(.x, universe_entrez))
  compare_entrez_sets <- compare_entrez_sets[purrr::map_int(compare_entrez_sets, length) >= min_genes_for_go]

  if (length(compare_entrez_sets) >= 2) {
    cc <- clusterProfiler::compareCluster(
      geneCluster = compare_entrez_sets,
      fun = "enrichGO",
      OrgDb = org.Mm.eg.db,
      keyType = "ENTREZID",
      ont = go_ontology,
      universe = universe_entrez,
      pAdjustMethod = "BH",
      pvalueCutoff = go_p_cutoff,
      qvalueCutoff = go_q_cutoff,
      readable = TRUE
    )

    cc_df <- as.data.frame(cc)

    readr::write_csv(
      cc_df,
      file.path(out_dir, "GO_compareCluster_focused_sets.csv")
    )

    if (nrow(cc_df) > 0) {
      p_cc <- enrichplot::dotplot(
        cc,
        showCategory = 10,
        font.size = 8
      ) +
        ggplot2::labs(
          title = "GO compareCluster: IL-1 and inhibitor overlap gene sets"
        ) +
        ggplot2::theme_bw(base_size = 9) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)
        )

      ggplot2::ggsave(
        file.path(out_dir, "PLOT_compareCluster_GO_focused_sets.pdf"),
        p_cc,
        width = 12,
        height = 8
      )
    }
  } else {
    message("Skipping compareCluster: fewer than two mapped gene sets.")
  }
} else {
  message("Skipping compareCluster: fewer than two input gene sets with enough genes.")
}

# -----------------------------
# Print summary
# -----------------------------

message("Done.")
message("Output directory:")
message(out_dir)

message("\nGene set sizes:")
print(go_input_sizes)

message("\nMain outputs:")
message(file.path(out_dir, "GO_input_gene_sets_long.csv"))
message(file.path(out_dir, "GO_input_gene_set_sizes.csv"))
message(file.path(out_dir, "GO_background_universe_genes.csv"))
message(file.path(out_dir, "GO_enrichment_results_long.csv"))
message(file.path(out_dir, "PLOT_combined_GO_focused_sets.pdf"))
message(file.path(out_dir, "PLOT_compareCluster_GO_focused_sets.pdf"))

message("\nTop GO terms:")
print(
  go_results_long %>%
    dplyr::select(
      set_name,
      ID,
      Description,
      GeneRatio,
      BgRatio,
      p_value,
      padj,
      q_value,
      Count
    ) %>%
    dplyr::arrange(padj, p_value) %>%
    dplyr::head(40)
)

# ============================================================
# Scatter plots:
#   1. IL1 vs Control log2FC on X, fibroblast Age log2FC on Y
#   2. AP-1 inhibitor vs Control log2FC on X, fibroblast Exercise log2FC on Y
#
# For both 3 h and 6 h.
#
# Requires objects from the previous script:
#   dge_long
#   fibro_age_ex
#   out_dir
# ============================================================

if (!exists("dge_long")) {
  stop("dge_long not found. Run the basic overlap script first.")
}

if (!exists("fibro_age_ex")) {
  stop("fibro_age_ex not found. Run the basic overlap script first.")
}

if (!exists("out_dir")) {
  out_dir <- file.path(base_dir, "Basic_IL1_vs_Control_Inhibitor_vs_Control_AgeEx_overlap")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
}

scatter_out_dir <- file.path(out_dir, "Correlation_scatterplots")
dir.create(scatter_out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Helper function for correlation labels
# -----------------------------

make_cor_label <- function(df, x_col, y_col, method = "pearson") {
  df_use <- df %>%
    dplyr::filter(
      !is.na(.data[[x_col]]),
      !is.na(.data[[y_col]])
    )

  if (nrow(df_use) < 3) {
    return(tibble::tibble(
      label = "n < 3",
      x = NA_real_,
      y = NA_real_
    ))
  }

  ct <- suppressWarnings(
    cor.test(
      df_use[[x_col]],
      df_use[[y_col]],
      method = method
    )
  )

  tibble::tibble(
    label = paste0(
      "r = ", signif(unname(ct$estimate), 3),
      "\nP = ", signif(ct$p.value, 3),
      "\nn = ", nrow(df_use)
    ),
    x = quantile(df_use[[x_col]], 0.05, na.rm = TRUE),
    y = quantile(df_use[[y_col]], 0.95, na.rm = TRUE)
  )
}

# -----------------------------
# 1. IL1 log2FC versus fibroblast Age log2FC
# -----------------------------

il1_age_scatter <- dge_long %>%
  dplyr::filter(
    contrast %in% c("IL1_vs_Veh_3h", "IL1_vs_Veh_6h")
  ) %>%
  dplyr::mutate(
    timepoint = dplyr::case_when(
      contrast == "IL1_vs_Veh_3h" ~ "3 h",
      contrast == "IL1_vs_Veh_6h" ~ "6 h",
      TRUE ~ NA_character_
    ),
    IL1_log2FC = log2FC,
    IL1_padj = padj,
    IL1_sig = sig
  ) %>%
  dplyr::select(
    gene,
    timepoint,
    IL1_log2FC,
    IL1_padj,
    IL1_sig
  ) %>%
  dplyr::inner_join(
    fibro_age_ex %>%
      dplyr::select(
        gene,
        Age_log2FC,
        Age_padj,
        Age_sig,
        Age_up,
        Age_down,
        Exercise_log2FC,
        Exercise_padj,
        Exercise_sig
      ),
    by = "gene"
  ) %>%
  dplyr::mutate(
    point_class = dplyr::case_when(
      IL1_sig & Age_sig & IL1_log2FC > 0 & Age_log2FC > 0 ~ "IL1 up + Age up",
      IL1_sig & Age_sig & IL1_log2FC < 0 & Age_log2FC < 0 ~ "IL1 down + Age down",
      IL1_sig & Age_sig & sign(IL1_log2FC) != sign(Age_log2FC) ~ "Opposite significant",
      IL1_sig & !Age_sig ~ "IL1 significant only",
      !IL1_sig & Age_sig ~ "Age significant only",
      TRUE ~ "Not significant"
    ),
    timepoint = factor(timepoint, levels = c("3 h", "6 h"))
  )

readr::write_csv(
  il1_age_scatter,
  file.path(scatter_out_dir, "IL1_vs_Control_log2FC_vs_fibroblast_Age_log2FC_scatter_data.csv")
)

il1_age_cor_labels <- il1_age_scatter %>%
  dplyr::group_by(timepoint) %>%
  dplyr::group_modify(~ make_cor_label(.x, "IL1_log2FC", "Age_log2FC", method = "pearson")) %>%
  dplyr::ungroup()

p_il1_age <- ggplot2::ggplot(
  il1_age_scatter,
  ggplot2::aes(x = IL1_log2FC, y = Age_log2FC)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  ggplot2::geom_point(
    ggplot2::aes(colour = point_class),
    alpha = 0.7,
    size = 1.7
  ) +
  ggplot2::geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    colour = "black",
    linewidth = 0.8
  ) +
  ggplot2::geom_text(
    data = il1_age_cor_labels,
    ggplot2::aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.3
  ) +
  ggplot2::facet_wrap(~ timepoint, nrow = 1) +
  ggplot2::scale_colour_manual(
    values = c(
      "IL1 up + Age up" = "#A81B1B",
      "IL1 down + Age down" = "#101369",
      "Opposite significant" = "#F99349",
      "IL1 significant only" = "#A53D99",
      "Age significant only" = "#42977E",
      "Not significant" = "grey78"
    )
  ) +
  ggplot2::labs(
    x = "IL1 vs Vehicle log2FC",
    y = "Fibroblast Age effect log2FC",
    colour = NULL,
    title = "IL-1 response versus fibroblast ageing effect",
    subtitle = "Each point is a gene measured in both datasets. Trend line = linear model; r = Pearson correlation."
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  file.path(scatter_out_dir, "PLOT_IL1_vs_Control_log2FC_vs_fibroblast_Age_log2FC.pdf"),
  p_il1_age,
  width = 10,
  height = 6
)

# -----------------------------
# 2. AP-1 inhibitor log2FC versus fibroblast Exercise log2FC
# -----------------------------
# Here AP-1 inhibitor = T5224_vs_Veh, not T5224preIL1_vs_IL1.

inhibitor_exercise_scatter <- dge_long %>%
  dplyr::filter(
    contrast %in% c("T5224_vs_Veh_3h", "T5224_vs_Veh_6h")
  ) %>%
  dplyr::mutate(
    timepoint = dplyr::case_when(
      contrast == "T5224_vs_Veh_3h" ~ "3 h",
      contrast == "T5224_vs_Veh_6h" ~ "6 h",
      TRUE ~ NA_character_
    ),
    Inhibitor_log2FC = log2FC,
    Inhibitor_padj = padj,
    Inhibitor_sig = sig
  ) %>%
  dplyr::select(
    gene,
    timepoint,
    Inhibitor_log2FC,
    Inhibitor_padj,
    Inhibitor_sig
  ) %>%
  dplyr::inner_join(
    fibro_age_ex %>%
      dplyr::select(
        gene,
        Exercise_log2FC,
        Exercise_padj,
        Exercise_sig,
        Exercise_up,
        Exercise_down,
        Age_log2FC,
        Age_padj,
        Age_sig
      ),
    by = "gene"
  ) %>%
  dplyr::mutate(
    point_class = dplyr::case_when(
      Inhibitor_sig & Exercise_sig & Inhibitor_log2FC > 0 & Exercise_log2FC > 0 ~ "Inhibitor up + Exercise up",
      Inhibitor_sig & Exercise_sig & Inhibitor_log2FC < 0 & Exercise_log2FC < 0 ~ "Inhibitor down + Exercise down",
      Inhibitor_sig & Exercise_sig & sign(Inhibitor_log2FC) != sign(Exercise_log2FC) ~ "Opposite significant",
      Inhibitor_sig & !Exercise_sig ~ "Inhibitor significant only",
      !Inhibitor_sig & Exercise_sig ~ "Exercise significant only",
      TRUE ~ "Not significant"
    ),
    timepoint = factor(timepoint, levels = c("3 h", "6 h"))
  )

readr::write_csv(
  inhibitor_exercise_scatter,
  file.path(scatter_out_dir, "Inhibitor_vs_Control_log2FC_vs_fibroblast_Exercise_log2FC_scatter_data.csv")
)

inhibitor_exercise_cor_labels <- inhibitor_exercise_scatter %>%
  dplyr::group_by(timepoint) %>%
  dplyr::group_modify(~ make_cor_label(.x, "Inhibitor_log2FC", "Exercise_log2FC", method = "pearson")) %>%
  dplyr::ungroup()

p_inhibitor_exercise <- ggplot2::ggplot(
  inhibitor_exercise_scatter,
  ggplot2::aes(x = Inhibitor_log2FC, y = Exercise_log2FC)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  ggplot2::geom_point(
    ggplot2::aes(colour = point_class),
    alpha = 0.7,
    size = 1.7
  ) +
  ggplot2::geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    colour = "black",
    linewidth = 0.8
  ) +
  ggplot2::geom_text(
    data = inhibitor_exercise_cor_labels,
    ggplot2::aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.3
  ) +
  ggplot2::facet_wrap(~ timepoint, nrow = 1) +
  ggplot2::scale_colour_manual(
    values = c(
      "Inhibitor up + Exercise up" = "#A81B1B",
      "Inhibitor down + Exercise down" = "#101369",
      "Opposite significant" = "#F99349",
      "Inhibitor significant only" = "#A53D99",
      "Exercise significant only" = "#42977E",
      "Not significant" = "grey78"
    )
  ) +
  ggplot2::labs(
    x = "AP-1 inhibitor vs Vehicle log2FC",
    y = "Fibroblast Exercise effect log2FC",
    colour = NULL,
    title = "AP-1 inhibitor response versus fibroblast exercise effect",
    subtitle = "AP-1 inhibitor contrast = T5224 vs Vehicle. Trend line = linear model; r = Pearson correlation."
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  file.path(scatter_out_dir, "PLOT_Inhibitor_vs_Control_log2FC_vs_fibroblast_Exercise_log2FC.pdf"),
  p_inhibitor_exercise,
  width = 10,
  height = 6
)

# -----------------------------
# 3. Combined long-format plot, if you want both in one PDF
# -----------------------------

combined_scatter <- dplyr::bind_rows(
  il1_age_scatter %>%
    dplyr::transmute(
      gene,
      comparison = "IL1 vs Control x Age",
      timepoint,
      x_log2FC = IL1_log2FC,
      y_log2FC = Age_log2FC,
      x_padj = IL1_padj,
      y_padj = Age_padj,
      x_sig = IL1_sig,
      y_sig = Age_sig,
      point_class = point_class
    ),
  inhibitor_exercise_scatter %>%
    dplyr::transmute(
      gene,
      comparison = "AP-1 inhibitor vs Control x Exercise",
      timepoint,
      x_log2FC = Inhibitor_log2FC,
      y_log2FC = Exercise_log2FC,
      x_padj = Inhibitor_padj,
      y_padj = Exercise_padj,
      x_sig = Inhibitor_sig,
      y_sig = Exercise_sig,
      point_class = point_class
    )
) %>%
  dplyr::mutate(
    comparison = factor(
      comparison,
      levels = c(
        "IL1 vs Control x Age",
        "AP-1 inhibitor vs Control x Exercise"
      )
    )
  )

combined_cor_labels <- combined_scatter %>%
  dplyr::group_by(comparison, timepoint) %>%
  dplyr::group_modify(~ make_cor_label(.x, "x_log2FC", "y_log2FC", method = "pearson")) %>%
  dplyr::ungroup()

readr::write_csv(
  combined_scatter,
  file.path(scatter_out_dir, "combined_correlation_scatter_data.csv")
)

p_combined <- ggplot2::ggplot(
  combined_scatter,
  ggplot2::aes(x = x_log2FC, y = y_log2FC)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  ggplot2::geom_point(
    alpha = 0.55,
    size = 1.4,
    colour = "grey35"
  ) +
  ggplot2::geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    colour = "#A81B1B",
    linewidth = 0.8
  ) +
  ggplot2::geom_text(
    data = combined_cor_labels,
    ggplot2::aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.1
  ) +
  ggplot2::facet_grid(comparison ~ timepoint, scales = "free") +
  ggplot2::labs(
    x = "Perturbation log2FC",
    y = "In vivo fibroblast log2FC",
    title = "Correlation between perturbation responses and in vivo fibroblast effects",
    subtitle = "Top row: IL1 vs Control compared with Age. Bottom row: AP-1 inhibitor vs Control compared with Exercise."
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(scatter_out_dir, "PLOT_combined_IL1_Age_and_Inhibitor_Exercise_correlation_scatterplots.pdf"),
  p_combined,
  width = 10,
  height = 8
)

# -----------------------------
# Print summary
# -----------------------------

message("Done.")
message("Scatter output directory:")
message(scatter_out_dir)

message("\nPlots:")
message(file.path(scatter_out_dir, "PLOT_IL1_vs_Control_log2FC_vs_fibroblast_Age_log2FC.pdf"))
message(file.path(scatter_out_dir, "PLOT_Inhibitor_vs_Control_log2FC_vs_fibroblast_Exercise_log2FC.pdf"))
message(file.path(scatter_out_dir, "PLOT_combined_IL1_Age_and_Inhibitor_Exercise_correlation_scatterplots.pdf"))

message("\nData:")
message(file.path(scatter_out_dir, "IL1_vs_Control_log2FC_vs_fibroblast_Age_log2FC_scatter_data.csv"))
message(file.path(scatter_out_dir, "Inhibitor_vs_Control_log2FC_vs_fibroblast_Exercise_log2FC_scatter_data.csv"))
message(file.path(scatter_out_dir, "combined_correlation_scatter_data.csv"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(forcats)
})

# ============================================================
# BED enrichment against IL-1 and AP-1 responsive genes
#
# Foreground BEDs:
#   aged_gained_fosjun_fg.bed
#   exer_lost_fosjun_fg.bed
#
# Background BED:
#   background_no_ap1.bed
#
# Gene-linking methods:
#   nearest gene
#   genes within +/-100 kb
#
# Response gene sets:
#   IL1_vs_Veh
#   T5224_vs_Veh
#   T5224preIL1_vs_IL1
#   IL1 induced and AP-1 inhibitor suppressed
#
# Region-overlap tests:
#   New foreground BEDs vs existing C-FOS/C-JUN BED files
#
# Fixed:
#   No rowwise() + .env$foreground_bed_set scoping bug
# ============================================================

base_dir <- publication_path("PAPER_OUTPUT_ROOT", "IL1_RNASeq")

out_dir <- file.path(base_dir, "FOSJUN_BED_IL1_AP1_enrichment")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dge_fdr_cutoff <- 0.05
dge_abs_log2fc_cutoff <- 0
window_bp <- 100000

new_bed_files <- tibble::tibble(
  bed_set = c(
    "aged_gained_fosjun_fg",
    "exer_lost_fosjun_fg",
    "background_no_ap1"
  ),
  bed_class = c(
    "foreground",
    "foreground",
    "background"
  ),
  file = file.path(
    base_dir,
    c(
      "aged_gained_fosjun_fg.bed",
      "exer_lost_fosjun_fg.bed",
      "background_no_ap1.bed"
    )
  )
)

missing_beds <- new_bed_files$file[!file.exists(new_bed_files$file)]

if (length(missing_beds) > 0) {
  stop(
    "Missing BED file(s):\n",
    paste(missing_beds, collapse = "\n")
  )
}

background_bed_set <- "background_no_ap1"

# -----------------------------
# Helper functions
# -----------------------------

guess_col <- function(df, candidates, fallback_regex = NULL, exclude_regex = NULL) {
  nms <- names(df)
  nms_clean <- tolower(nms)
  candidates_clean <- tolower(candidates)

  exact_hit <- which(nms_clean %in% candidates_clean)

  if (!is.null(exclude_regex) && length(exact_hit) > 0) {
    exact_hit <- exact_hit[!stringr::str_detect(nms_clean[exact_hit], exclude_regex)]
  }

  if (length(exact_hit) > 0) {
    return(nms[exact_hit[1]])
  }

  if (!is.null(fallback_regex)) {
    regex_hit <- which(stringr::str_detect(nms_clean, fallback_regex))

    if (!is.null(exclude_regex) && length(regex_hit) > 0) {
      regex_hit <- regex_hit[!stringr::str_detect(nms_clean[regex_hit], exclude_regex)]
    }

    if (length(regex_hit) > 0) {
      return(nms[regex_hit[1]])
    }
  }

  NA_character_
}

read_dge_file <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)

  gene_col <- guess_col(
    df,
    candidates = c(
      "gene", "genes", "symbol", "gene_symbol", "gene_name",
      "external_gene_name", "feature", "features", "X", "...1"
    ),
    fallback_regex = "gene|symbol|feature"
  )

  p_col <- guess_col(
    df,
    candidates = c(
      "padj", "p_val_adj", "p_val_adj.", "FDR", "adj.P.Val",
      "qval", "q_value", "qvalue", "fdr_bh"
    ),
    fallback_regex = "padj|adj.*p|p.*adj|fdr|qval|q_value"
  )

  lfc_col <- guess_col(
    df,
    candidates = c(
      "log2FoldChange", "logFC", "avg_log2FC", "avg_logFC",
      "log2FC", "estimate", "coef"
    ),
    fallback_regex = "log2.*fold|logfc|log2fc|avg_log",
    exclude_regex = "padj|fdr|p_val|pvalue|qval"
  )

  if (is.na(gene_col)) {
    stop("Could not detect gene column in: ", basename(path))
  }

  if (is.na(p_col)) {
    stop("Could not detect FDR/padj column in: ", basename(path))
  }

  if (is.na(lfc_col)) {
    stop("Could not detect logFC column in: ", basename(path))
  }

  contrast <- basename(path) %>%
    stringr::str_remove("^DGE_") %>%
    stringr::str_remove("\\.csv$")

  df %>%
    dplyr::transmute(
      gene = as.character(.data[[gene_col]]),
      contrast = contrast,
      log2FC = suppressWarnings(as.numeric(.data[[lfc_col]])),
      padj = suppressWarnings(as.numeric(.data[[p_col]])),
      gene_col_used = gene_col,
      lfc_col_used = lfc_col,
      p_col_used = p_col
    ) %>%
    dplyr::mutate(
      gene = stringr::str_replace(gene, "\\.\\d+$", ""),
      sig = !is.na(padj) &
        padj <= dge_fdr_cutoff &
        !is.na(log2FC) &
        abs(log2FC) >= dge_abs_log2fc_cutoff,
      direction = dplyr::case_when(
        sig & log2FC > 0 ~ "up",
        sig & log2FC < 0 ~ "down",
        TRUE ~ "ns"
      )
    ) %>%
    dplyr::filter(!is.na(gene), gene != "", !is.na(log2FC), !is.na(padj)) %>%
    dplyr::distinct(gene, contrast, .keep_all = TRUE)
}

clean_gr <- function(gr) {
  gr <- GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode = "coarse")

  suppressWarnings({
    GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
  })

  gr <- sort(gr)

  gr$region_id <- paste0(
    as.character(GenomicRanges::seqnames(gr)),
    ":",
    GenomicRanges::start(gr),
    "-",
    GenomicRanges::end(gr)
  )

  gr
}

import_bed_clean <- function(path) {
  gr <- rtracklayer::import(path)

  if (!inherits(gr, "GRanges")) {
    stop("Could not import as GRanges: ", path)
  }

  clean_gr(gr)
}

sig_label <- function(padj) {
  dplyr::case_when(
    !is.na(padj) & padj < 0.001 ~ "***",
    !is.na(padj) & padj < 0.01 ~ "**",
    !is.na(padj) & padj < 0.05 ~ "*",
    !is.na(padj) & padj < 0.1 ~ ".",
    TRUE ~ ""
  )
}

run_fg_bg_gene_fisher <- function(fg_genes, bg_genes, target_genes, measured_genes) {
  fg_genes <- unique(fg_genes)
  bg_genes <- unique(bg_genes)
  target_genes <- unique(target_genes)
  measured_genes <- unique(measured_genes)

  universe_genes <- intersect(unique(c(fg_genes, bg_genes)), measured_genes)

  fg_genes <- intersect(fg_genes, universe_genes)
  bg_genes <- intersect(bg_genes, universe_genes)
  target_genes <- intersect(target_genes, universe_genes)

  in_fg <- universe_genes %in% fg_genes
  in_target <- universe_genes %in% target_genes

  overlap_genes <- sort(intersect(fg_genes, target_genes))

  mat <- matrix(
    c(
      sum(in_fg & in_target),
      sum(in_fg & !in_target),
      sum(!in_fg & in_target),
      sum(!in_fg & !in_target)
    ),
    nrow = 2,
    byrow = TRUE
  )

  rownames(mat) <- c("foreground_linked_genes", "background_linked_genes")
  colnames(mat) <- c("response_genes", "not_response_genes")

  if (
    length(universe_genes) == 0 ||
      any(rowSums(mat) == 0) ||
      any(colSums(mat) == 0)
  ) {
    return(tibble::tibble(
      universe_n = length(universe_genes),
      foreground_gene_n = length(fg_genes),
      background_gene_n = length(bg_genes),
      response_gene_n = length(target_genes),
      overlap_n = length(overlap_genes),
      expected_overlap_n = NA_real_,
      odds_ratio = NA_real_,
      p_value = NA_real_,
      overlap_genes = paste(overlap_genes, collapse = ";")
    ))
  }

  ft <- fisher.test(mat, alternative = "greater")

  tibble::tibble(
    universe_n = length(universe_genes),
    foreground_gene_n = length(fg_genes),
    background_gene_n = length(bg_genes),
    response_gene_n = length(target_genes),
    overlap_n = length(overlap_genes),
    expected_overlap_n = length(fg_genes) * length(target_genes) / length(universe_genes),
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value,
    overlap_genes = paste(overlap_genes, collapse = ";")
  )
}

run_region_overlap_fisher <- function(fg_gr, bg_gr, ref_gr) {
  fg_has_overlap <- GenomicRanges::countOverlaps(fg_gr, ref_gr, ignore.strand = TRUE) > 0
  bg_has_overlap <- GenomicRanges::countOverlaps(bg_gr, ref_gr, ignore.strand = TRUE) > 0

  mat <- matrix(
    c(
      sum(fg_has_overlap),
      sum(!fg_has_overlap),
      sum(bg_has_overlap),
      sum(!bg_has_overlap)
    ),
    nrow = 2,
    byrow = TRUE
  )

  rownames(mat) <- c("foreground_regions", "background_regions")
  colnames(mat) <- c("overlaps_reference", "does_not_overlap_reference")

  if (any(rowSums(mat) == 0) || any(colSums(mat) == 0)) {
    return(tibble::tibble(
      foreground_region_n = length(fg_gr),
      background_region_n = length(bg_gr),
      foreground_overlap_n = sum(fg_has_overlap),
      background_overlap_n = sum(bg_has_overlap),
      foreground_overlap_fraction = mean(fg_has_overlap),
      background_overlap_fraction = mean(bg_has_overlap),
      odds_ratio = NA_real_,
      p_value = NA_real_
    ))
  }

  ft <- fisher.test(mat, alternative = "greater")

  tibble::tibble(
    foreground_region_n = length(fg_gr),
    background_region_n = length(bg_gr),
    foreground_overlap_n = sum(fg_has_overlap),
    background_overlap_n = sum(bg_has_overlap),
    foreground_overlap_fraction = mean(fg_has_overlap),
    background_overlap_fraction = mean(bg_has_overlap),
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value
  )
}

regions_to_tbl <- function(gr, prefix = "query") {
  tibble::tibble(
    !!paste0(prefix, "_chr") := as.character(GenomicRanges::seqnames(gr)),
    !!paste0(prefix, "_start") := GenomicRanges::start(gr),
    !!paste0(prefix, "_end") := GenomicRanges::end(gr),
    !!paste0(prefix, "_width") := GenomicRanges::width(gr),
    !!paste0(prefix, "_region_id") := gr$region_id
  )
}

# -----------------------------
# Gene annotation
# -----------------------------

genes_gr <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)

genes_gr <- GenomeInfoDb::keepStandardChromosomes(genes_gr, pruning.mode = "coarse")

suppressWarnings({
  GenomeInfoDb::seqlevelsStyle(genes_gr) <- "UCSC"
})

genes_gr$ENTREZID <- names(genes_gr)

genes_gr$SYMBOL <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys = genes_gr$ENTREZID,
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

genes_gr <- genes_gr[!is.na(genes_gr$SYMBOL)]
genes_gr <- sort(genes_gr)

link_nearest_genes <- function(peaks_gr, genes_gr, bed_set) {
  idx <- GenomicRanges::nearest(peaks_gr, genes_gr, ignore.strand = TRUE)

  keep <- !is.na(idx)

  if (sum(keep) == 0) {
    return(tibble::tibble())
  }

  peak_keep <- peaks_gr[keep]
  gene_keep <- genes_gr[idx[keep]]

  tibble::tibble(
    bed_set = bed_set,
    link_method = "nearest",
    region_id = peak_keep$region_id,
    peak_chr = as.character(GenomicRanges::seqnames(peak_keep)),
    peak_start = GenomicRanges::start(peak_keep),
    peak_end = GenomicRanges::end(peak_keep),
    gene = gene_keep$SYMBOL,
    gene_chr = as.character(GenomicRanges::seqnames(gene_keep)),
    gene_start = GenomicRanges::start(gene_keep),
    gene_end = GenomicRanges::end(gene_keep),
    distance_to_gene = GenomicRanges::distance(peak_keep, gene_keep)
  ) %>%
    dplyr::filter(!is.na(gene), gene != "") %>%
    dplyr::distinct()
}

link_100kb_genes <- function(peaks_gr, genes_gr, bed_set, window_bp = 100000) {
  peak_window <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(peaks_gr),
    ranges = IRanges::IRanges(
      start = pmax(1, GenomicRanges::start(peaks_gr) - window_bp),
      end = GenomicRanges::end(peaks_gr) + window_bp
    ),
    strand = "*"
  )

  suppressWarnings({
    GenomeInfoDb::seqlevelsStyle(peak_window) <- GenomeInfoDb::seqlevelsStyle(peaks_gr)[1]
  })

  GenomeInfoDb::seqlevels(peak_window) <- GenomeInfoDb::seqlevels(peaks_gr)
  peak_window$region_id <- peaks_gr$region_id

  hits <- GenomicRanges::findOverlaps(
    peak_window,
    genes_gr,
    ignore.strand = TRUE
  )

  if (length(hits) == 0) {
    return(tibble::tibble())
  }

  peak_hit <- peaks_gr[S4Vectors::queryHits(hits)]
  gene_hit <- genes_gr[S4Vectors::subjectHits(hits)]

  tibble::tibble(
    bed_set = bed_set,
    link_method = paste0("within_", window_bp / 1000, "kb"),
    region_id = peak_hit$region_id,
    peak_chr = as.character(GenomicRanges::seqnames(peak_hit)),
    peak_start = GenomicRanges::start(peak_hit),
    peak_end = GenomicRanges::end(peak_hit),
    gene = gene_hit$SYMBOL,
    gene_chr = as.character(GenomicRanges::seqnames(gene_hit)),
    gene_start = GenomicRanges::start(gene_hit),
    gene_end = GenomicRanges::end(gene_hit),
    distance_to_gene = GenomicRanges::distance(peak_hit, gene_hit)
  ) %>%
    dplyr::filter(!is.na(gene), gene != "") %>%
    dplyr::distinct()
}

# -----------------------------
# Read DGE files and make response gene sets
# -----------------------------

dge_files <- list.files(base_dir, pattern = "^DGE_.*\\.csv$", full.names = TRUE)

if (length(dge_files) == 0) {
  stop("No DGE_*.csv files found in: ", base_dir)
}

dge_long <- purrr::map_dfr(dge_files, read_dge_file)

readr::write_csv(
  dge_long,
  file.path(out_dir, "all_DGE_long.csv")
)

column_check <- dge_long %>%
  dplyr::distinct(contrast, gene_col_used, lfc_col_used, p_col_used) %>%
  dplyr::arrange(contrast)

readr::write_csv(
  column_check,
  file.path(out_dir, "DGE_column_detection_check.csv")
)

response_sets_basic <- dge_long %>%
  dplyr::filter(sig, direction %in% c("up", "down")) %>%
  dplyr::mutate(
    response_family = dplyr::case_when(
      stringr::str_detect(contrast, "^IL1_vs_Veh") ~ "IL1_response",
      stringr::str_detect(contrast, "^T5224preIL1_vs_IL1") ~ "AP1_inhibitor_response_under_IL1",
      stringr::str_detect(contrast, "^T5224_vs_Veh") ~ "AP1_inhibitor_response_alone",
      TRUE ~ NA_character_
    ),
    timepoint = dplyr::case_when(
      stringr::str_detect(contrast, "3h$") ~ "3h",
      stringr::str_detect(contrast, "6h$") ~ "6h",
      TRUE ~ NA_character_
    ),
    response_set = paste(response_family, direction, timepoint, sep = "__"),
    response_label = dplyr::case_when(
      response_family == "IL1_response" ~ paste("IL-1", direction, timepoint),
      response_family == "AP1_inhibitor_response_under_IL1" ~ paste("AP-1 inhibitor under IL-1", direction, timepoint),
      response_family == "AP1_inhibitor_response_alone" ~ paste("AP-1 inhibitor alone", direction, timepoint),
      TRUE ~ response_set
    )
  ) %>%
  dplyr::filter(!is.na(response_family), !is.na(timepoint)) %>%
  dplyr::select(
    gene,
    response_family,
    response_set,
    response_label,
    contrast,
    direction,
    timepoint,
    log2FC,
    padj
  ) %>%
  dplyr::distinct()

response_sets_any_time <- response_sets_basic %>%
  dplyr::mutate(
    timepoint = "any_time",
    response_set = paste(response_family, direction, timepoint, sep = "__"),
    response_label = dplyr::case_when(
      response_family == "IL1_response" ~ paste("IL-1", direction, "any time"),
      response_family == "AP1_inhibitor_response_under_IL1" ~ paste("AP-1 inhibitor under IL-1", direction, "any time"),
      response_family == "AP1_inhibitor_response_alone" ~ paste("AP-1 inhibitor alone", direction, "any time"),
      TRUE ~ response_set
    )
  ) %>%
  dplyr::select(
    gene,
    response_family,
    response_set,
    response_label,
    direction,
    timepoint
  ) %>%
  dplyr::distinct()

dge_wide <- dge_long %>%
  dplyr::select(gene, contrast, log2FC, padj, sig, direction) %>%
  tidyr::pivot_wider(
    names_from = contrast,
    values_from = c(log2FC, padj, sig, direction),
    names_sep = "__"
  )

make_dependency_sets_one_timepoint <- function(dge_wide, timepoint) {
  il1_lfc <- paste0("log2FC__IL1_vs_Veh_", timepoint)
  il1_sig <- paste0("sig__IL1_vs_Veh_", timepoint)
  inh_lfc <- paste0("log2FC__T5224preIL1_vs_IL1_", timepoint)
  inh_sig <- paste0("sig__T5224preIL1_vs_IL1_", timepoint)

  needed <- c(il1_lfc, il1_sig, inh_lfc, inh_sig)

  if (!all(needed %in% names(dge_wide))) {
    return(tibble::tibble())
  }

  dplyr::bind_rows(
    dge_wide %>%
      dplyr::filter(
        .data[[il1_sig]],
        .data[[il1_lfc]] > 0,
        .data[[inh_sig]],
        .data[[inh_lfc]] < 0
      ) %>%
      dplyr::transmute(
        gene,
        response_family = "IL1_induced_AP1_inhibitor_suppressed",
        response_set = paste0("IL1_induced_AP1_inhibitor_suppressed__", timepoint),
        response_label = paste("IL-1 induced and AP-1 inhibitor suppressed", timepoint),
        direction = "model",
        timepoint = timepoint
      ),
    dge_wide %>%
      dplyr::filter(
        .data[[il1_sig]],
        .data[[il1_lfc]] < 0,
        .data[[inh_sig]],
        .data[[inh_lfc]] > 0
      ) %>%
      dplyr::transmute(
        gene,
        response_family = "IL1_repressed_AP1_inhibitor_rescued",
        response_set = paste0("IL1_repressed_AP1_inhibitor_rescued__", timepoint),
        response_label = paste("IL-1 repressed and AP-1 inhibitor rescued", timepoint),
        direction = "model",
        timepoint = timepoint
      )
  ) %>%
    dplyr::distinct()
}

response_sets_model <- dplyr::bind_rows(
  make_dependency_sets_one_timepoint(dge_wide, "3h"),
  make_dependency_sets_one_timepoint(dge_wide, "6h")
)

response_sets_model_any <- response_sets_model %>%
  dplyr::mutate(
    timepoint = "any_time",
    response_set = paste0(response_family, "__any_time"),
    response_label = dplyr::case_when(
      response_family == "IL1_induced_AP1_inhibitor_suppressed" ~ "IL-1 induced and AP-1 inhibitor suppressed any time",
      response_family == "IL1_repressed_AP1_inhibitor_rescued" ~ "IL-1 repressed and AP-1 inhibitor rescued any time",
      TRUE ~ response_set
    )
  ) %>%
  dplyr::select(
    gene,
    response_family,
    response_set,
    response_label,
    direction,
    timepoint
  ) %>%
  dplyr::distinct()

response_sets_all <- dplyr::bind_rows(
  response_sets_basic %>%
    dplyr::select(gene, response_family, response_set, response_label, direction, timepoint),
  response_sets_any_time,
  response_sets_model,
  response_sets_model_any
) %>%
  dplyr::distinct()

response_sets_all <- response_sets_all %>%
  dplyr::mutate(
    response_label = factor(
      response_label,
      levels = unique(response_label)
    )
  )

readr::write_csv(
  response_sets_all,
  file.path(out_dir, "IL1_AP1_response_gene_sets.csv")
)

response_set_sizes <- response_sets_all %>%
  dplyr::count(
    response_family,
    response_set,
    response_label,
    direction,
    timepoint,
    name = "n_genes"
  ) %>%
  dplyr::arrange(response_family, timepoint, direction)

readr::write_csv(
  response_set_sizes,
  file.path(out_dir, "IL1_AP1_response_gene_set_sizes.csv")
)

# -----------------------------
# Import new BEDs and link to genes
# -----------------------------

bed_gr_list <- purrr::set_names(
  purrr::map(new_bed_files$file, import_bed_clean),
  new_bed_files$bed_set
)

linked_genes <- purrr::map_dfr(names(bed_gr_list), function(bed_set) {
  gr <- bed_gr_list[[bed_set]]

  dplyr::bind_rows(
    link_nearest_genes(gr, genes_gr, bed_set),
    link_100kb_genes(gr, genes_gr, bed_set, window_bp = window_bp)
  )
}) %>%
  dplyr::left_join(
    new_bed_files %>%
      dplyr::select(bed_set, bed_class),
    by = "bed_set"
  ) %>%
  dplyr::distinct()

readr::write_csv(
  linked_genes,
  file.path(out_dir, "new_BED_linked_genes_nearest_and_100kb.csv")
)

linked_gene_summary <- linked_genes %>%
  dplyr::count(bed_set, bed_class, link_method, name = "linked_gene_records") %>%
  dplyr::left_join(
    linked_genes %>%
      dplyr::group_by(bed_set, link_method) %>%
      dplyr::summarise(
        unique_linked_genes = dplyr::n_distinct(gene),
        .groups = "drop"
      ),
    by = c("bed_set", "link_method")
  ) %>%
  dplyr::arrange(bed_set, link_method)

readr::write_csv(
  linked_gene_summary,
  file.path(out_dir, "new_BED_linked_gene_summary.csv")
)

# -----------------------------
# Enrichment of IL-1/AP-1 response genes near new BEDs
# Fixed with pmap_dfr, no rowwise/.env scoping
# -----------------------------

foreground_sets <- new_bed_files %>%
  dplyr::filter(bed_class == "foreground") %>%
  dplyr::pull(bed_set)

link_methods <- unique(linked_genes$link_method)

measured_genes <- unique(dge_long$gene)

gene_test_grid <- tidyr::expand_grid(
  foreground_bed_set = foreground_sets,
  link_method = link_methods,
  response_set = unique(response_sets_all$response_set)
) %>%
  dplyr::left_join(
    response_sets_all %>%
      dplyr::distinct(
        response_set,
        response_family,
        response_label,
        direction,
        timepoint
      ),
    by = "response_set"
  )

gene_enrichment_results <- purrr::pmap_dfr(
  gene_test_grid,
  function(foreground_bed_set,
           link_method,
           response_set,
           response_family,
           response_label,
           direction,
           timepoint) {

    fg_genes <- linked_genes %>%
      dplyr::filter(
        bed_set == foreground_bed_set,
        link_method == !!link_method
      ) %>%
      dplyr::pull(gene) %>%
      unique()

    bg_genes <- linked_genes %>%
      dplyr::filter(
        bed_set == background_bed_set,
        link_method == !!link_method
      ) %>%
      dplyr::pull(gene) %>%
      unique()

    target_genes <- response_sets_all %>%
      dplyr::filter(response_set == !!response_set) %>%
      dplyr::pull(gene) %>%
      unique()

    run_fg_bg_gene_fisher(
      fg_genes = fg_genes,
      bg_genes = bg_genes,
      target_genes = target_genes,
      measured_genes = measured_genes
    ) %>%
      dplyr::mutate(
        foreground_bed_set = foreground_bed_set,
        link_method = link_method,
        response_set = response_set,
        response_family = response_family,
        response_label = response_label,
        direction = direction,
        timepoint = timepoint,
        .before = 1
      )
  }
) %>%
  dplyr::mutate(
    padj_global = p.adjust(p_value, method = "BH"),
    padj_by_foreground_and_link = ave(
      p_value,
      foreground_bed_set,
      link_method,
      FUN = function(x) p.adjust(x, method = "BH")
    ),
    log2_odds_ratio = log2(odds_ratio),
    sig = sig_label(padj_global)
  ) %>%
  dplyr::arrange(padj_global, p_value)

readr::write_csv(
  gene_enrichment_results,
  file.path(out_dir, "new_BED_vs_IL1_AP1_response_gene_enrichment.csv")
)

gene_overlap_long <- gene_enrichment_results %>%
  dplyr::filter(!is.na(overlap_genes), overlap_genes != "") %>%
  dplyr::select(
    foreground_bed_set,
    link_method,
    response_family,
    response_set,
    response_label,
    direction,
    timepoint,
    overlap_genes
  ) %>%
  tidyr::separate_rows(overlap_genes, sep = ";") %>%
  dplyr::rename(gene = overlap_genes) %>%
  dplyr::left_join(
    dge_long %>%
      dplyr::select(gene, contrast, log2FC, padj, sig, direction) %>%
      dplyr::rename(dge_direction = direction),
    by = "gene"
  ) %>%
  dplyr::arrange(
    foreground_bed_set,
    link_method,
    response_family,
    response_label,
    gene
  )

readr::write_csv(
  gene_overlap_long,
  file.path(out_dir, "new_BED_vs_IL1_AP1_response_overlap_genes_long.csv")
)

# -----------------------------
# Plot gene enrichment
# -----------------------------

plot_gene_enrichment <- gene_enrichment_results %>%
  dplyr::filter(
    response_family %in% c(
      "IL1_response",
      "AP1_inhibitor_response_under_IL1",
      "AP1_inhibitor_response_alone",
      "IL1_induced_AP1_inhibitor_suppressed",
      "IL1_repressed_AP1_inhibitor_rescued"
    )
  ) %>%
  dplyr::mutate(
    foreground_bed_set = factor(
      foreground_bed_set,
      levels = c("aged_gained_fosjun_fg", "exer_lost_fosjun_fg")
    ),
    link_method = factor(
      link_method,
      levels = c("nearest", paste0("within_", window_bp / 1000, "kb"))
    ),
    response_label = forcats::fct_inorder(as.character(response_label))
  )

p_gene_enrichment <- ggplot2::ggplot(
  plot_gene_enrichment,
  ggplot2::aes(
    x = response_label,
    y = paste(foreground_bed_set, link_method, sep = " | "),
    size = overlap_n,
    fill = log2_odds_ratio
  )
) +
  ggplot2::geom_point(shape = 21, colour = "black", alpha = 0.9) +
  ggplot2::geom_text(
    ggplot2::aes(label = sig),
    size = 3.2,
    vjust = 0.35,
    colour = "black"
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#377EB8",
    mid = "white",
    high = "#E41A1C",
    midpoint = 0,
    na.value = "grey85"
  ) +
  ggplot2::labs(
    x = "IL-1 / AP-1 response gene set",
    y = "BED set | gene-linking method",
    size = "Overlap genes",
    fill = "log2 odds ratio",
    title = "Enrichment of IL-1/AP-1 response genes near new FOS/JUN BED regions",
    subtitle = "Foreground BED-linked genes tested against background_no_ap1 linked genes. *, **, *** are global BH FDR"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    panel.grid.minor = ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(out_dir, "PLOT1_new_BED_IL1_AP1_response_gene_enrichment_dotplot.pdf"),
  p_gene_enrichment,
  width = 16,
  height = 6.2
)

plot_gene_enrichment_model <- gene_enrichment_results %>%
  dplyr::filter(
    response_family %in% c(
      "IL1_response",
      "AP1_inhibitor_response_under_IL1",
      "AP1_inhibitor_response_alone",
      "IL1_induced_AP1_inhibitor_suppressed"
    ),
    timepoint %in% c("3h", "6h", "any_time"),
    direction %in% c("up", "down", "model")
  ) %>%
  dplyr::mutate(
    response_label = forcats::fct_inorder(as.character(response_label)),
    foreground_bed_set = factor(
      foreground_bed_set,
      levels = c("aged_gained_fosjun_fg", "exer_lost_fosjun_fg")
    )
  )

p_gene_enrichment_model <- ggplot2::ggplot(
  plot_gene_enrichment_model,
  ggplot2::aes(
    x = response_label,
    y = foreground_bed_set,
    size = overlap_n,
    fill = log2_odds_ratio
  )
) +
  ggplot2::geom_point(shape = 21, colour = "black", alpha = 0.9) +
  ggplot2::geom_text(
    ggplot2::aes(label = sig),
    size = 3.3,
    vjust = 0.35,
    colour = "black"
  ) +
  ggplot2::facet_wrap(~ link_method, ncol = 1) +
  ggplot2::scale_fill_gradient2(
    low = "#377EB8",
    mid = "white",
    high = "#E41A1C",
    midpoint = 0,
    na.value = "grey85"
  ) +
  ggplot2::labs(
    x = "Response gene set",
    y = "Foreground BED set",
    size = "Overlap genes",
    fill = "log2 odds ratio",
    title = "Focused enrichment summary",
    subtitle = "Best signal would be enrichment for IL-1 induced and AP-1 inhibitor suppressed genes"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    panel.grid.minor = ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(out_dir, "PLOT2_focused_new_BED_IL1_AP1_response_enrichment.pdf"),
  p_gene_enrichment_model,
  width = 14,
  height = 7.5
)

# -----------------------------
# Direct region overlap with existing C-FOS/C-JUN BED files
# -----------------------------

cfos_cjun_bed_files <- list.files(
  base_dir,
  pattern = "^C-(FOS|JUN).*\\.bed$",
  full.names = TRUE
)

cfos_cjun_tbl <- tibble::tibble(file = cfos_cjun_bed_files) %>%
  dplyr::mutate(
    reference_set = basename(file) %>%
      stringr::str_remove("\\.bed$"),
    factor = dplyr::case_when(
      stringr::str_detect(reference_set, "C-FOS") ~ "C-FOS",
      stringr::str_detect(reference_set, "C-JUN") ~ "C-JUN",
      TRUE ~ "Other"
    ),
    reference_class = dplyr::case_when(
      stringr::str_detect(reference_set, "age-gained") ~ "age-gained",
      stringr::str_detect(reference_set, "background") ~ "background",
      TRUE ~ "other"
    )
  ) %>%
  dplyr::filter(factor %in% c("C-FOS", "C-JUN"))

if (nrow(cfos_cjun_tbl) == 0) {
  warning("No C-FOS or C-JUN BED files found.")
} else {
  ref_gr_list <- purrr::set_names(
    purrr::map(cfos_cjun_tbl$file, import_bed_clean),
    cfos_cjun_tbl$reference_set
  )

  bg_gr <- bed_gr_list[[background_bed_set]]

  # -----------------------------
  # Region enrichment with pmap_dfr
  # -----------------------------

  region_test_grid <- tidyr::expand_grid(
    foreground_bed_set = foreground_sets,
    reference_set = names(ref_gr_list)
  ) %>%
    dplyr::left_join(
      cfos_cjun_tbl %>%
        dplyr::select(reference_set, factor, reference_class),
      by = "reference_set"
    )

  region_overlap_enrichment <- purrr::pmap_dfr(
    region_test_grid,
    function(foreground_bed_set,
             reference_set,
             factor,
             reference_class) {

      fg_gr <- bed_gr_list[[foreground_bed_set]]
      ref_gr <- ref_gr_list[[reference_set]]
      bg_gr <- bed_gr_list[[background_bed_set]]

      run_region_overlap_fisher(
        fg_gr = fg_gr,
        bg_gr = bg_gr,
        ref_gr = ref_gr
      ) %>%
        dplyr::mutate(
          foreground_bed_set = foreground_bed_set,
          reference_set = reference_set,
          factor = factor,
          reference_class = reference_class,
          .before = 1
        )
    }
  ) %>%
    dplyr::mutate(
      padj_global = p.adjust(p_value, method = "BH"),
      log2_odds_ratio = log2(odds_ratio),
      sig = sig_label(padj_global)
    ) %>%
    dplyr::arrange(padj_global, p_value)

  readr::write_csv(
    region_overlap_enrichment,
    file.path(out_dir, "new_BED_region_overlap_enrichment_vs_CFOS_CJUN.csv")
  )

  # -----------------------------
  # Direct overlap summary with pmap_dfr
  # -----------------------------

  direct_overlap_grid <- tidyr::expand_grid(
    query_bed_set = c(foreground_sets, background_bed_set),
    reference_set = names(ref_gr_list)
  ) %>%
    dplyr::left_join(
      cfos_cjun_tbl %>%
        dplyr::select(reference_set, factor, reference_class),
      by = "reference_set"
    )

  direct_overlap_summary <- purrr::pmap_dfr(
    direct_overlap_grid,
    function(query_bed_set,
             reference_set,
             factor,
             reference_class) {

      query_gr <- bed_gr_list[[query_bed_set]]
      reference_gr <- ref_gr_list[[reference_set]]

      query_regions_overlapping_reference <- sum(
        GenomicRanges::countOverlaps(
          query_gr,
          reference_gr,
          ignore.strand = TRUE
        ) > 0
      )

      reference_regions_overlapping_query <- sum(
        GenomicRanges::countOverlaps(
          reference_gr,
          query_gr,
          ignore.strand = TRUE
        ) > 0
      )

      tibble::tibble(
        query_bed_set = query_bed_set,
        reference_set = reference_set,
        factor = factor,
        reference_class = reference_class,
        query_n = length(query_gr),
        reference_n = length(reference_gr),
        query_regions_overlapping_reference = query_regions_overlapping_reference,
        reference_regions_overlapping_query = reference_regions_overlapping_query,
        query_overlap_fraction = query_regions_overlapping_reference / length(query_gr),
        reference_overlap_fraction = reference_regions_overlapping_query / length(reference_gr)
      )
    }
  ) %>%
    dplyr::arrange(query_bed_set, factor, reference_class, reference_set)

  readr::write_csv(
    direct_overlap_summary,
    file.path(out_dir, "new_BED_direct_region_overlap_summary_vs_CFOS_CJUN.csv")
  )

  overlapping_regions_long <- purrr::map_dfr(foreground_sets, function(query_set) {
    query_gr <- bed_gr_list[[query_set]]

    purrr::map_dfr(names(ref_gr_list), function(reference_set) {
      ref_gr <- ref_gr_list[[reference_set]]

      hits <- GenomicRanges::findOverlaps(query_gr, ref_gr, ignore.strand = TRUE)

      if (length(hits) == 0) {
        return(tibble::tibble())
      }

      query_hit <- query_gr[S4Vectors::queryHits(hits)]
      ref_hit <- ref_gr[S4Vectors::subjectHits(hits)]

      query_tbl <- regions_to_tbl(query_hit, prefix = "query")
      ref_tbl <- regions_to_tbl(ref_hit, prefix = "reference")

      dplyr::bind_cols(query_tbl, ref_tbl) %>%
        dplyr::mutate(
          query_bed_set = query_set,
          reference_set = reference_set,
          .before = 1
        )
    })
  }) %>%
    dplyr::left_join(
      cfos_cjun_tbl %>%
        dplyr::select(reference_set, factor, reference_class),
      by = "reference_set"
    ) %>%
    dplyr::arrange(
      query_bed_set,
      reference_set,
      query_chr,
      query_start,
      reference_start
    )

  readr::write_csv(
    overlapping_regions_long,
    file.path(out_dir, "new_BED_overlapping_regions_long_vs_CFOS_CJUN.csv")
  )

  p_region_enrichment <- region_overlap_enrichment %>%
    dplyr::mutate(
      foreground_bed_set = factor(
        foreground_bed_set,
        levels = c("aged_gained_fosjun_fg", "exer_lost_fosjun_fg")
      ),
      reference_set = factor(reference_set, levels = unique(reference_set))
    ) %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = reference_set,
        y = foreground_bed_set,
        size = foreground_overlap_n,
        fill = log2_odds_ratio
      )
    ) +
    ggplot2::geom_point(shape = 21, colour = "black", alpha = 0.9) +
    ggplot2::geom_text(
      ggplot2::aes(label = sig),
      size = 3.2,
      vjust = 0.35,
      colour = "black"
    ) +
    ggplot2::facet_grid(reference_class ~ factor, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_gradient2(
      low = "#377EB8",
      mid = "white",
      high = "#E41A1C",
      midpoint = 0,
      na.value = "grey85"
    ) +
    ggplot2::labs(
      x = "Reference C-FOS / C-JUN BED set",
      y = "New foreground BED set",
      size = "Foreground regions overlapping reference",
      fill = "log2 odds ratio",
      title = "Direct genomic overlap of new FOS/JUN BEDs with existing C-FOS/C-JUN BEDs",
      subtitle = "Foreground overlap tested against background_no_ap1. *, **, *** are global BH FDR"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )

  ggplot2::ggsave(
    file.path(out_dir, "PLOT3_region_overlap_enrichment_vs_CFOS_CJUN.pdf"),
    p_region_enrichment,
    width = 13,
    height = 5.5
  )

  p_region_fraction <- direct_overlap_summary %>%
    dplyr::filter(query_bed_set %in% c(foreground_sets, background_bed_set)) %>%
    dplyr::mutate(
      query_bed_set = factor(
        query_bed_set,
        levels = c("aged_gained_fosjun_fg", "exer_lost_fosjun_fg", "background_no_ap1")
      ),
      reference_set = factor(reference_set, levels = unique(reference_set))
    ) %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = reference_set,
        y = query_overlap_fraction,
        fill = query_bed_set
      )
    ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.8),
      colour = "black",
      linewidth = 0.2
    ) +
    ggplot2::facet_grid(reference_class ~ factor, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_manual(
      values = c(
        "aged_gained_fosjun_fg" = "#A81B1B",
        "exer_lost_fosjun_fg" = "#A53D99",
        "background_no_ap1" = "grey60"
      )
    ) +
    ggplot2::labs(
      x = "Reference C-FOS / C-JUN BED set",
      y = "Fraction of query regions overlapping reference",
      fill = "Query BED set",
      title = "Raw fraction of new BED regions overlapping C-FOS/C-JUN BEDs"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )

  ggplot2::ggsave(
    file.path(out_dir, "PLOT4_region_overlap_fraction_vs_CFOS_CJUN.pdf"),
    p_region_fraction,
    width = 13,
    height = 6
  )
}

# -----------------------------
# Ranked tables
# -----------------------------

top_gene_enrichments <- gene_enrichment_results %>%
  dplyr::select(
    foreground_bed_set,
    link_method,
    response_family,
    response_label,
    foreground_gene_n,
    background_gene_n,
    response_gene_n,
    overlap_n,
    expected_overlap_n,
    odds_ratio,
    p_value,
    padj_global,
    overlap_genes
  ) %>%
  dplyr::arrange(padj_global, p_value)

readr::write_csv(
  top_gene_enrichments,
  file.path(out_dir, "TOP_TABLE_gene_enrichments_ranked.csv")
)

if (exists("region_overlap_enrichment")) {
  top_region_enrichments <- region_overlap_enrichment %>%
    dplyr::select(
      foreground_bed_set,
      reference_set,
      factor,
      reference_class,
      foreground_region_n,
      background_region_n,
      foreground_overlap_n,
      background_overlap_n,
      foreground_overlap_fraction,
      background_overlap_fraction,
      odds_ratio,
      p_value,
      padj_global
    ) %>%
    dplyr::arrange(padj_global, p_value)

  readr::write_csv(
    top_region_enrichments,
    file.path(out_dir, "TOP_TABLE_region_overlap_enrichments_ranked.csv")
  )
}

# -----------------------------
# Print summary
# -----------------------------

message("Done.")
message("Output directory:")
message(out_dir)

message("\nMain gene-linking outputs:")
message(file.path(out_dir, "new_BED_linked_genes_nearest_and_100kb.csv"))
message(file.path(out_dir, "new_BED_vs_IL1_AP1_response_gene_enrichment.csv"))
message(file.path(out_dir, "new_BED_vs_IL1_AP1_response_overlap_genes_long.csv"))
message(file.path(out_dir, "TOP_TABLE_gene_enrichments_ranked.csv"))

message("\nMain region-overlap outputs:")
message(file.path(out_dir, "new_BED_region_overlap_enrichment_vs_CFOS_CJUN.csv"))
message(file.path(out_dir, "new_BED_direct_region_overlap_summary_vs_CFOS_CJUN.csv"))
message(file.path(out_dir, "new_BED_overlapping_regions_long_vs_CFOS_CJUN.csv"))
message(file.path(out_dir, "TOP_TABLE_region_overlap_enrichments_ranked.csv"))

message("\nMain plots:")
message(file.path(out_dir, "PLOT1_new_BED_IL1_AP1_response_gene_enrichment_dotplot.pdf"))
message(file.path(out_dir, "PLOT2_focused_new_BED_IL1_AP1_response_enrichment.pdf"))
message(file.path(out_dir, "PLOT3_region_overlap_enrichment_vs_CFOS_CJUN.pdf"))
message(file.path(out_dir, "PLOT4_region_overlap_fraction_vs_CFOS_CJUN.pdf"))

message("\nTop gene enrichments:")
print(
  top_gene_enrichments %>%
    dplyr::select(
      foreground_bed_set,
      link_method,
      response_label,
      overlap_n,
      expected_overlap_n,
      odds_ratio,
      p_value,
      padj_global
    ) %>%
    dplyr::head(30)
)

if (exists("top_region_enrichments")) {
  message("\nTop region overlap enrichments:")
  print(
    top_region_enrichments %>%
      dplyr::select(
        foreground_bed_set,
        reference_set,
        factor,
        reference_class,
        foreground_overlap_n,
        background_overlap_n,
        foreground_overlap_fraction,
        background_overlap_fraction,
        odds_ratio,
        p_value,
        padj_global
      ) %>%
      dplyr::head(30)
  )
}

# ============================================================
# Plot expression/DGE effects for genes linked to BED regions
#
# Requires either:
#   linked_genes and dge_long objects from the BED enrichment script
#
# Or reads:
#   new_BED_linked_genes_nearest_and_100kb.csv
#   all_DGE_long.csv
#
# Outputs:
#   - linked DEG table
#   - heatmap of linked DEG log2FC across contrasts
#   - dotplot of linked DEG effects
#   - barplot of linked DEG counts by contrast/direction
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(forcats)
})

# -----------------------------
# Paths
# -----------------------------

base_dir <- publication_path("PAPER_OUTPUT_ROOT", "IL1_RNASeq")

bed_out_dir <- file.path(base_dir, "FOSJUN_BED_IL1_AP1_enrichment")

linked_deg_out_dir <- file.path(bed_out_dir, "linked_DEG_expression_plots")
dir.create(linked_deg_out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load objects if not already in memory
# -----------------------------

if (!exists("linked_genes")) {
  linked_gene_file <- file.path(
    bed_out_dir,
    "new_BED_linked_genes_nearest_and_100kb.csv"
  )

  if (!file.exists(linked_gene_file)) {
    stop("Could not find linked genes file: ", linked_gene_file)
  }

  linked_genes <- readr::read_csv(linked_gene_file, show_col_types = FALSE)
}

if (!exists("dge_long")) {
  dge_file <- file.path(
    bed_out_dir,
    "all_DGE_long.csv"
  )

  if (!file.exists(dge_file)) {
    stop("Could not find DGE long file: ", dge_file)
  }

  dge_long <- readr::read_csv(dge_file, show_col_types = FALSE)
}

# -----------------------------
# Settings
# -----------------------------

fdr_cutoff <- 0.05
top_n_genes <- 80

contrast_order <- c(
  "IL1_vs_Veh_3h",
  "T5224_vs_Veh_3h",
  "T5224preIL1_vs_IL1_3h",
  "T5224preIL1_vs_Veh_3h",
  "IL1_vs_Veh_6h",
  "T5224_vs_Veh_6h",
  "T5224preIL1_vs_IL1_6h",
  "T5224preIL1_vs_Veh_6h"
)

contrast_labels <- c(
  "IL1_vs_Veh_3h" = "IL-1 vs Veh 3 h",
  "T5224_vs_Veh_3h" = "T5224 vs Veh 3 h",
  "T5224preIL1_vs_IL1_3h" = "T5224+IL1 vs IL1 3 h",
  "T5224preIL1_vs_Veh_3h" = "T5224+IL1 vs Veh 3 h",
  "IL1_vs_Veh_6h" = "IL-1 vs Veh 6 h",
  "T5224_vs_Veh_6h" = "T5224 vs Veh 6 h",
  "T5224preIL1_vs_IL1_6h" = "T5224+IL1 vs IL1 6 h",
  "T5224preIL1_vs_Veh_6h" = "T5224+IL1 vs Veh 6 h"
)

bed_set_order <- c(
  "aged_gained_fosjun_fg",
  "exer_lost_fosjun_fg",
  "background_no_ap1"
)

bed_set_labels <- c(
  "aged_gained_fosjun_fg" = "Aged gained FOS/JUN",
  "exer_lost_fosjun_fg" = "Exercise lost FOS/JUN",
  "background_no_ap1" = "Background no AP-1"
)

link_method_order <- c(
  "nearest",
  "within_100kb"
)

# -----------------------------
# Clean and join BED-linked genes with DGE effects
# -----------------------------

dge_clean <- dge_long %>%
  dplyr::mutate(
    gene = as.character(gene),
    gene = stringr::str_replace(gene, "\\.\\d+$", ""),
    contrast = as.character(contrast),
    log2FC = as.numeric(log2FC),
    padj = as.numeric(padj),
    sig = !is.na(padj) & padj <= fdr_cutoff & !is.na(log2FC),
    dge_direction = dplyr::case_when(
      sig & log2FC > 0 ~ "up",
      sig & log2FC < 0 ~ "down",
      TRUE ~ "ns"
    ),
    contrast_label = dplyr::recode(contrast, !!!contrast_labels, .default = contrast)
  )

linked_clean <- linked_genes %>%
  dplyr::mutate(
    gene = as.character(gene),
    gene = stringr::str_replace(gene, "\\.\\d+$", ""),
    bed_set = as.character(bed_set),
    link_method = as.character(link_method),
    bed_set_label = dplyr::recode(bed_set, !!!bed_set_labels, .default = bed_set),
    link_method_label = dplyr::case_when(
      link_method == "nearest" ~ "Nearest gene",
      link_method == "within_100kb" ~ "Within 100 kb",
      TRUE ~ link_method
    )
  )

linked_gene_summary <- linked_clean %>%
  dplyr::group_by(gene, bed_set, bed_set_label, link_method, link_method_label, bed_class) %>%
  dplyr::summarise(
    n_linked_regions = dplyr::n_distinct(region_id),
    min_distance_to_gene = suppressWarnings(min(distance_to_gene, na.rm = TRUE)),
    median_distance_to_gene = suppressWarnings(median(distance_to_gene, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    min_distance_to_gene = ifelse(is.infinite(min_distance_to_gene), NA_real_, min_distance_to_gene),
    median_distance_to_gene = ifelse(is.infinite(median_distance_to_gene), NA_real_, median_distance_to_gene)
  )

linked_dge <- linked_gene_summary %>%
  dplyr::inner_join(dge_clean, by = "gene") %>%
  dplyr::mutate(
    bed_set = factor(bed_set, levels = bed_set_order),
    bed_set_label = factor(
      bed_set_label,
      levels = bed_set_labels[bed_set_order]
    ),
    link_method = factor(link_method, levels = link_method_order),
    link_method_label = factor(
      link_method_label,
      levels = c("Nearest gene", "Within 100 kb")
    ),
    contrast = factor(contrast, levels = contrast_order),
    contrast_label = factor(
      contrast_label,
      levels = contrast_labels[contrast_order]
    )
  ) %>%
  dplyr::filter(!is.na(contrast))

linked_degs_only <- linked_dge %>%
  dplyr::filter(sig)

readr::write_csv(
  linked_dge,
  file.path(linked_deg_out_dir, "BED_linked_genes_all_DGE_effects.csv")
)

readr::write_csv(
  linked_degs_only,
  file.path(linked_deg_out_dir, "BED_linked_genes_that_are_DEGs.csv")
)

# -----------------------------
# Summarise linked DEG counts
# -----------------------------

linked_deg_counts <- linked_degs_only %>%
  dplyr::count(
    bed_set,
    bed_set_label,
    link_method,
    link_method_label,
    contrast,
    contrast_label,
    dge_direction,
    name = "n_linked_DEGs"
  ) %>%
  dplyr::arrange(bed_set, link_method, contrast, dge_direction)

readr::write_csv(
  linked_deg_counts,
  file.path(linked_deg_out_dir, "BED_linked_DEG_counts_by_contrast.csv")
)

# -----------------------------
# Pick genes for plotting
# -----------------------------
# Prioritise foreground-linked genes with:
#   many significant DGE contrasts
#   strong absolute log2FC
#   many linked regions

top_genes <- linked_degs_only %>%
  dplyr::filter(bed_class == "foreground") %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    n_bed_sets = dplyr::n_distinct(bed_set),
    n_link_methods = dplyr::n_distinct(link_method),
    n_sig_contrasts = dplyr::n_distinct(contrast),
    max_abs_log2FC = max(abs(log2FC), na.rm = TRUE),
    min_padj = min(padj, na.rm = TRUE),
    total_linked_regions = sum(n_linked_regions, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_bed_sets),
    dplyr::desc(n_sig_contrasts),
    dplyr::desc(max_abs_log2FC),
    min_padj
  ) %>%
  dplyr::slice_head(n = top_n_genes) %>%
  dplyr::pull(gene)

if (length(top_genes) == 0) {
  stop("No linked DEGs found for foreground BED sets.")
}

readr::write_csv(
  tibble::tibble(gene = top_genes),
  file.path(linked_deg_out_dir, "genes_used_for_linked_DEG_plots.csv")
)

# -----------------------------
# Plot 1: heatmap of DGE log2FC for linked DEGs
# -----------------------------

heatmap_df <- linked_dge %>%
  dplyr::filter(
    gene %in% top_genes,
    bed_class == "foreground"
  ) %>%
  dplyr::mutate(
    sig_label = ifelse(sig, "*", ""),
    gene = factor(gene, levels = rev(top_genes))
  )

p_heatmap <- ggplot2::ggplot(
  heatmap_df,
  ggplot2::aes(x = contrast_label, y = gene, fill = log2FC)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.18) +
  ggplot2::geom_text(
    ggplot2::aes(label = sig_label),
    size = 2.3
  ) +
  ggplot2::facet_grid(
    link_method_label ~ bed_set_label,
    scales = "free_y",
    space = "free_y"
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#377EB8",
    mid = "white",
    high = "#E41A1C",
    midpoint = 0,
    na.value = "grey90"
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    fill = "DGE log2FC",
    title = "DGE expression effects for genes linked to FOS/JUN BED regions",
    subtitle = paste0(
      "Top ", length(top_genes),
      " foreground-linked DEG genes. * = FDR <= ",
      fdr_cutoff
    )
  ) +
  ggplot2::theme_bw(base_size = 8.5) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
    panel.grid = ggplot2::element_blank(),
    strip.text.y = ggplot2::element_text(angle = 0)
  )

ggplot2::ggsave(
  file.path(linked_deg_out_dir, "PLOT1_linked_DEG_log2FC_heatmap.pdf"),
  p_heatmap,
  width = 13,
  height = max(6, length(top_genes) * 0.16)
)

# -----------------------------
# Plot 2: dotplot of significant linked DEGs only
# -----------------------------

dotplot_df <- linked_degs_only %>%
  dplyr::filter(
    gene %in% top_genes,
    bed_class == "foreground"
  ) %>%
  dplyr::mutate(
    gene = factor(gene, levels = rev(top_genes)),
    neglog10_padj = -log10(padj)
  )

p_dot <- ggplot2::ggplot(
  dotplot_df,
  ggplot2::aes(
    x = contrast_label,
    y = gene,
    size = neglog10_padj,
    fill = log2FC
  )
) +
  ggplot2::geom_point(shape = 21, colour = "black", alpha = 0.9) +
  ggplot2::facet_grid(
    link_method_label ~ bed_set_label,
    scales = "free_y",
    space = "free_y"
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#377EB8",
    mid = "white",
    high = "#E41A1C",
    midpoint = 0,
    na.value = "grey90"
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    size = "-log10(FDR)",
    fill = "DGE log2FC",
    title = "Significant DGE effects among genes linked to FOS/JUN BED regions",
    subtitle = "Only significant DGE entries are shown"
  ) +
  ggplot2::theme_bw(base_size = 8.5) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
    panel.grid.minor = ggplot2::element_blank(),
    strip.text.y = ggplot2::element_text(angle = 0)
  )

ggplot2::ggsave(
  file.path(linked_deg_out_dir, "PLOT2_linked_significant_DEG_dotplot.pdf"),
  p_dot,
  width = 13,
  height = max(6, length(top_genes) * 0.16)
)

# -----------------------------
# Plot 3: linked DEG counts by BED set, contrast, direction
# -----------------------------

count_plot_df <- linked_deg_counts %>%
  dplyr::filter(
    bed_set %in% c("aged_gained_fosjun_fg", "exer_lost_fosjun_fg")
  ) %>%
  dplyr::mutate(
    dge_direction = factor(dge_direction, levels = c("up", "down"))
  )

p_counts <- ggplot2::ggplot(
  count_plot_df,
  ggplot2::aes(
    x = contrast_label,
    y = n_linked_DEGs,
    fill = dge_direction
  )
) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.75),
    colour = "black",
    linewidth = 0.2
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = n_linked_DEGs),
    position = ggplot2::position_dodge(width = 0.75),
    vjust = -0.25,
    size = 2.6
  ) +
  ggplot2::facet_grid(link_method_label ~ bed_set_label) +
  ggplot2::scale_fill_manual(
    values = c(
      "up" = "#A81B1B",
      "down" = "#377EB8"
    )
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Linked DEG genes",
    fill = "DGE direction",
    title = "Number of BED-linked genes that are DEGs",
    subtitle = paste0("DGE threshold: FDR <= ", fdr_cutoff)
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  file.path(linked_deg_out_dir, "PLOT3_linked_DEG_counts_by_contrast.pdf"),
  p_counts,
  width = 13,
  height = 7
)

# -----------------------------
# Plot 4: ranked linked DEG genes by strongest DGE effect
# -----------------------------

rank_df <- linked_degs_only %>%
  dplyr::filter(
    bed_class == "foreground",
    gene %in% top_genes
  ) %>%
  dplyr::group_by(gene, bed_set_label, link_method_label) %>%
  dplyr::summarise(
    strongest_abs_log2FC = max(abs(log2FC), na.rm = TRUE),
    strongest_signed_log2FC = log2FC[which.max(abs(log2FC))][1],
    best_contrast = as.character(contrast_label[which.max(abs(log2FC))][1]),
    min_padj = min(padj, na.rm = TRUE),
    n_sig_contrasts = dplyr::n_distinct(contrast),
    n_linked_regions = max(n_linked_regions, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_sig_contrasts),
    dplyr::desc(strongest_abs_log2FC),
    min_padj
  ) %>%
  dplyr::slice_head(n = 50) %>%
  dplyr::mutate(
    gene = factor(gene, levels = rev(unique(gene)))
  )

readr::write_csv(
  rank_df,
  file.path(linked_deg_out_dir, "ranked_foreground_linked_DEG_genes.csv")
)

p_rank <- ggplot2::ggplot(
  rank_df,
  ggplot2::aes(
    x = strongest_signed_log2FC,
    y = gene,
    fill = bed_set_label
  )
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_col(colour = "black", linewidth = 0.2) +
  ggplot2::facet_wrap(~ link_method_label) +
  ggplot2::scale_fill_manual(
    values = c(
      "Aged gained FOS/JUN" = "#A81B1B",
      "Exercise lost FOS/JUN" = "#A53D99"
    )
  ) +
  ggplot2::labs(
    x = "Strongest signed DGE log2FC",
    y = NULL,
    fill = "BED set",
    title = "Ranked foreground BED-linked DEG genes",
    subtitle = "Bar shows strongest signed DGE effect across IL-1/T5224 contrasts"
  ) +
  ggplot2::theme_bw(base_size = 9.5) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  file.path(linked_deg_out_dir, "PLOT4_ranked_foreground_linked_DEG_genes.pdf"),
  p_rank,
  width = 10,
  height = max(5.5, nrow(rank_df) * 0.13)
)

# -----------------------------
# Optional focused plot:
# only IL-1 and T5224preIL1_vs_IL1 contrasts
# -----------------------------

focused_contrasts <- c(
  "IL1_vs_Veh_3h",
  "T5224preIL1_vs_IL1_3h",
  "IL1_vs_Veh_6h",
  "T5224preIL1_vs_IL1_6h"
)

focused_top_genes <- linked_degs_only %>%
  dplyr::filter(
    bed_class == "foreground",
    contrast %in% focused_contrasts
  ) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    n_sig_focused_contrasts = dplyr::n_distinct(contrast),
    max_abs_log2FC = max(abs(log2FC), na.rm = TRUE),
    min_padj = min(padj, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_sig_focused_contrasts),
    dplyr::desc(max_abs_log2FC),
    min_padj
  ) %>%
  dplyr::slice_head(n = 50) %>%
  dplyr::pull(gene)

if (length(focused_top_genes) > 0) {
  focused_df <- linked_dge %>%
    dplyr::filter(
      gene %in% focused_top_genes,
      bed_class == "foreground",
      contrast %in% focused_contrasts
    ) %>%
    dplyr::mutate(
      gene = factor(gene, levels = rev(focused_top_genes)),
      sig_label = ifelse(sig, "*", "")
    )

  p_focused <- ggplot2::ggplot(
    focused_df,
    ggplot2::aes(x = contrast_label, y = gene, fill = log2FC)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.18) +
    ggplot2::geom_text(
      ggplot2::aes(label = sig_label),
      size = 2.4
    ) +
    ggplot2::facet_grid(
      link_method_label ~ bed_set_label,
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#377EB8",
      mid = "white",
      high = "#E41A1C",
      midpoint = 0,
      na.value = "grey90"
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      fill = "DGE log2FC",
      title = "Focused IL-1/AP-1 inhibitor expression effects for BED-linked genes",
      subtitle = "IL1_vs_Veh and T5224preIL1_vs_IL1 only. * = FDR <= 0.05"
    ) +
    ggplot2::theme_bw(base_size = 8.5) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      panel.grid = ggplot2::element_blank(),
      strip.text.y = ggplot2::element_text(angle = 0)
    )

  ggplot2::ggsave(
    file.path(linked_deg_out_dir, "PLOT5_focused_IL1_AP1_linked_DEG_heatmap.pdf"),
    p_focused,
    width = 11,
    height = max(6, length(focused_top_genes) * 0.16)
  )
}

# -----------------------------
# Print summary
# -----------------------------

message("Done.")
message("Output directory:")
message(linked_deg_out_dir)

message("\nMain tables:")
message(file.path(linked_deg_out_dir, "BED_linked_genes_all_DGE_effects.csv"))
message(file.path(linked_deg_out_dir, "BED_linked_genes_that_are_DEGs.csv"))
message(file.path(linked_deg_out_dir, "BED_linked_DEG_counts_by_contrast.csv"))
message(file.path(linked_deg_out_dir, "ranked_foreground_linked_DEG_genes.csv"))

message("\nMain plots:")
message(file.path(linked_deg_out_dir, "PLOT1_linked_DEG_log2FC_heatmap.pdf"))
message(file.path(linked_deg_out_dir, "PLOT2_linked_significant_DEG_dotplot.pdf"))
message(file.path(linked_deg_out_dir, "PLOT3_linked_DEG_counts_by_contrast.pdf"))
message(file.path(linked_deg_out_dir, "PLOT4_ranked_foreground_linked_DEG_genes.pdf"))
message(file.path(linked_deg_out_dir, "PLOT5_focused_IL1_AP1_linked_DEG_heatmap.pdf"))

message("\nTop linked DEG genes:")
print(
  rank_df %>%
    dplyr::select(
      gene,
      bed_set_label,
      link_method_label,
      strongest_signed_log2FC,
      best_contrast,
      min_padj,
      n_sig_contrasts,
      n_linked_regions
    ) %>%
    dplyr::head(40)
)



# Check emily 
#!/usr/bin/env Rscript
# Pathway enrichment on the full il1_RNA DEG sets (no fib-universe filter).
# Background = all 13,264 genes tested by edgeR in the bulk RNA.
# Tests GO BP, KEGG, MSigDB Hallmark, Reactome.

suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(msigdbr)
})

base <- publication_path("PAPER_OUTPUT_ROOT", "IL1_RNASeq/")
dge_dir <- file.path(base)
out_dir <- file.path(base, "enrichment")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)  # Added: fwrite() below needs this folder to exist

load_dge <- function(c, logfc_thr=0.585, p_thr=0.05) {
  d <- fread(file.path(dge_dir, sprintf("DGE_%s.csv", c)))
  d[, up := !is.na(PValue) & PValue < p_thr & logFC >  logfc_thr]
  d[, dn := !is.na(PValue) & PValue < p_thr & logFC < -logfc_thr]
  d
}

contrasts <- list(
  T5224_6h_up             = load_dge("T5224_vs_Veh_6h"),
  T5224_6h_dn             = load_dge("T5224_vs_Veh_6h")
)

get_set <- function(nm) {
  d <- contrasts[[nm]]
  dir <- sub(".*_", "", nm)
  g <- if (dir == "up") d[up==TRUE, gene_name] else d[dn==TRUE, gene_name]
  list(genes = g, universe = d$gene_name)
}

symbol_to_entrez <- function(sym) {
  x <- suppressMessages(bitr(sym, fromType="SYMBOL", toType="ENTREZID",
                               OrgDb=org.Mm.eg.db))
  x$ENTREZID
}

h_t2g <- as.data.table(msigdbr(species="Mus musculus", collection="H"))
h_t2g <- h_t2g[, .(gs_name, entrez_gene=as.character(ncbi_gene))]
reactome_t2g <- as.data.table(msigdbr(species="Mus musculus",
                                        collection="C2",
                                        subcollection="CP:REACTOME"))
reactome_t2g <- reactome_t2g[, .(gs_name, entrez_gene=as.character(ncbi_gene))]

all_results <- list()
for (nm in names(contrasts)) {
  s <- get_set(nm)
  cat(sprintf("\n=== %s : %d genes ===\n", nm, length(s$genes)))
  entrez <- symbol_to_entrez(s$genes)
  uni <- symbol_to_entrez(s$universe)
  if (length(entrez) < 5) next

  go <- enrichGO(gene=entrez, universe=uni, OrgDb=org.Mm.eg.db,
                  keyType="ENTREZID", ont="BP",
                  pvalueCutoff=0.05, qvalueCutoff=0.1, readable=TRUE)
  if (!is.null(go) && nrow(as.data.frame(go))>0) {
    df <- as.data.frame(go); df$set <- nm; df$db <- "GO_BP"
    all_results[[paste0(nm,"_GO")]] <- df
  }

  kegg <- enrichKEGG(gene=entrez, universe=uni, organism="mmu",
                      pvalueCutoff=0.05, qvalueCutoff=0.1)
  if (!is.null(kegg) && nrow(as.data.frame(kegg))>0) {
    df <- as.data.frame(setReadable(kegg, OrgDb=org.Mm.eg.db,
                                       keyType="ENTREZID"))
    df$set <- nm; df$db <- "KEGG"
    all_results[[paste0(nm,"_KEGG")]] <- df
  }

  hall <- enricher(gene=entrez, universe=uni, TERM2GENE=h_t2g,
                    pvalueCutoff=0.05, qvalueCutoff=0.1)
  if (!is.null(hall) && nrow(as.data.frame(hall))>0) {
    df <- as.data.frame(setReadable(hall, OrgDb=org.Mm.eg.db,
                                       keyType="ENTREZID"))
    df$set <- nm; df$db <- "MSigDB_Hallmark"
    all_results[[paste0(nm,"_Hallmark")]] <- df
  }

  rea <- enricher(gene=entrez, universe=uni, TERM2GENE=reactome_t2g,
                   pvalueCutoff=0.05, qvalueCutoff=0.1)
  if (!is.null(rea) && nrow(as.data.frame(rea))>0) {
    df <- as.data.frame(setReadable(rea, OrgDb=org.Mm.eg.db,
                                       keyType="ENTREZID"))
    df$set <- nm; df$db <- "Reactome"
    all_results[[paste0(nm,"_Reactome")]] <- df
  }
}

all_res <- rbindlist(all_results, fill=TRUE)
setcolorder(all_res, c("set","db","ID","Description","pvalue","p.adjust",
                        "qvalue","Count","GeneRatio","BgRatio","geneID"))
fwrite(all_res, file.path(out_dir, "pathway_enrichment_full_il1RNA.tsv"),
       sep="\t")
cat("\nDone. Rows per set & db:\n")
print(all_res[, .N, by=.(set, db)])

base_dir <- publication_path("PAPER_OUTPUT_ROOT", "IL1_RNASeq")

out_dir <- file.path(base_dir, "DGE_overlap_AP1_dependency")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fdr_cutoff <- 0.05
abs_log2fc_cutoff <- 0
top_label_n <- 20
top_heatmap_n <- 60

# -----------------------------
# Input files
# -----------------------------

dge_paths <- list.files(base_dir, pattern = "^DGE_.*\\.csv$", full.names = TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
})

# -----------------------------
# Inputs
# -----------------------------

base_dir <- publication_path("PAPER_OUTPUT_ROOT", "IL1_RNASeq")

age_ex_rds <- file.path(base_dir, "All_DEG_Age_Ex.rds")
t5224_file <- file.path(base_dir, "DGE_T5224_vs_Veh_6h.csv")

out_pdf <- file.path(base_dir, "manual_genes_Fibroblast_T5224_Age_Exercise_heatmap.pdf")
out_csv <- file.path(base_dir, "manual_genes_Fibroblast_T5224_Age_Exercise_heatmap_values.csv")
missing_csv <- file.path(base_dir, "manual_genes_missing_from_heatmap_inputs.csv")

genes <- c(
  "Adamtsl3", "Bmper", "Tbx20", "Tsc22d3", "Aldh1a1", "mt-Co3", "Bnc2",
  "Dnm3os", "Ptprj", "mt-Cytb", "Bcl2l1", "Lpar1", "mt-Co2", "Fbn1",
  "Parva", "Mid1", "Mmp14", "Arhgap18", "Zcchc24", "Slit2", "Plekhh2",
  "Pcolce", "Ccnd3", "Pknox2", "Runx1t1", "Alpk1", "Map4k3", "Wwc2",
  "Greb1l", "Zdhhc20", "Nedd4", "Hs2st1", "Acaca", "Ntn4", "Cntln",
  "Stom", "Col6a2", "Ahdc1", "Itpr1", "Myo1d", "Zmynd8", "Sp100",
  "Rnf150", "Col6a1", "Stxbp6", "Tpst1", "Efemp1", "Fmo2", "Dclk1",
  "Adamts2", "Frmd4b", "Tgfbr3", "Dennd2a", "Klhl13", "Cebpd"
)

# -----------------------------
# Helpers
# -----------------------------

flatten_rds <- function(x) {
  if (is.data.frame(x)) return(as_tibble(x))
  
  if (is.list(x)) {
    return(bind_rows(lapply(names(x), function(nm) {
      if (!is.data.frame(x[[nm]])) return(NULL)
      as_tibble(x[[nm]]) %>% mutate(rds_list_name = nm)
    })))
  }
  
  stop("All_DEG_Age_Ex.rds is not a data.frame or list of data.frames")
}

pick_gene_col <- function(df) {
  hits <- intersect(c("gene", "Gene", "genes", "symbol", "Symbol", "gene_name"), names(df))
  if (length(hits) == 0) stop("No gene column found. Columns are: ", paste(names(df), collapse = ", "))
  hits[1]
}

pick_lfc_col <- function(df) {
  hits <- intersect(c("log2FoldChange", "logFC", "log2FC", "avg_log2FC", "avg_logFC"), names(df))
  if (length(hits) == 0) stop("No logFC column found. Columns are: ", paste(names(df), collapse = ", "))
  hits[1]
}

pick_padj_col <- function(df) {
  hits <- intersect(c("padj", "FDR", "adj.P.Val", "p_val_adj", "qval"), names(df))
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

# -----------------------------
# Read fibroblast age/exercise DEG
# -----------------------------

age_ex_raw <- readRDS(age_ex_rds)
age_ex <- flatten_rds(age_ex_raw)

age_ex_gene_col <- pick_gene_col(age_ex)

age_ex <- age_ex %>%
  rename(gene = all_of(age_ex_gene_col)) %>%
  mutate(gene = str_replace(as.character(gene), "\\.\\d+$", ""))

if (!"CellType" %in% names(age_ex)) {
  if ("celltype" %in% names(age_ex)) {
    age_ex <- age_ex %>% rename(CellType = celltype)
  } else if ("rds_list_name" %in% names(age_ex)) {
    age_ex <- age_ex %>% rename(CellType = rds_list_name)
  } else {
    stop("No CellType/celltype/rds_list_name column found in All_DEG_Age_Ex.rds")
  }
}

fibro_age_ex <- age_ex %>%
  filter(CellType == "Fibroblast" | str_detect(tolower(CellType), "^fibro")) %>%
  transmute(
    gene,
    Ageing = log2FoldChange.Age,
    Ageing_padj = padj.Age,
    Exercise = log2FoldChange.Ex,
    Exercise_padj = padj.Ex
  ) %>%
  filter(gene %in% genes) %>%
  group_by(gene) %>%
  summarise(
    Ageing = mean(Ageing, na.rm = TRUE),
    Ageing_padj = min(Ageing_padj, na.rm = TRUE),
    Exercise = mean(Exercise, na.rm = TRUE),
    Exercise_padj = min(Exercise_padj, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Ageing_padj = ifelse(is.infinite(Ageing_padj), NA_real_, Ageing_padj),
    Exercise_padj = ifelse(is.infinite(Exercise_padj), NA_real_, Exercise_padj)
  )

# -----------------------------
# Read T5224 AP-1 inhibitor DEG
# -----------------------------

read_t5224 <- function(file, effect_name) {
  df <- read.csv(file)
  
  gene_col <- pick_gene_col(df)
  lfc_col <- pick_lfc_col(df)
  padj_col <- pick_padj_col(df)
  
  df %>%
    rename(gene = all_of(gene_col)) %>%
    mutate(gene = str_replace(as.character(gene), "\\.\\d+$", "")) %>%
    filter(gene %in% genes) %>%
    transmute(
      gene,
      effect = effect_name,
      log2FC = as.numeric(.data[[lfc_col]]),
      padj = if (!is.na(padj_col)) as.numeric(.data[[padj_col]]) else NA_real_
    ) %>%
    group_by(gene, effect) %>%
    summarise(
      log2FC = mean(log2FC, na.rm = TRUE),
      padj = min(padj, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(padj = ifelse(is.infinite(padj), NA_real_, padj))
}

t5224 <- read_t5224(t5224_file, "T5224 vs Veh 6h")

# -----------------------------
# Build plotting table
# -----------------------------

age_ex_long <- bind_rows(
  fibro_age_ex %>%
    transmute(gene, effect = "Fibroblast ageing", log2FC = Ageing, padj = Ageing_padj),
  fibro_age_ex %>%
    transmute(gene, effect = "Fibroblast exercise", log2FC = Exercise, padj = Exercise_padj)
)

effect_levels <- c(
  "T5224 vs Veh 6h",
  "Fibroblast ageing",
  "Fibroblast exercise"
)

heat_df <- expand_grid(
  gene = genes,
  effect = effect_levels
) %>%
  left_join(
    bind_rows(t5224, age_ex_long),
    by = c("gene", "effect")
  ) %>%
  mutate(
    gene = factor(gene, levels = rev(genes)),
    effect = factor(effect, levels = effect_levels),
    sig = case_when(
      !is.na(padj) & padj <= 0.001 ~ "***",
      !is.na(padj) & padj <= 0.01 ~ "**",
      !is.na(padj) & padj <= 0.05 ~ "*",
      TRUE ~ ""
    )
  )

missing_df <- heat_df %>%
  group_by(gene) %>%
  summarise(
    n_missing = sum(is.na(log2FC)),
    missing_effects = paste(effect[is.na(log2FC)], collapse = "; "),
    .groups = "drop"
  ) %>%
  filter(n_missing > 0)

write.csv(heat_df, out_csv, row.names = FALSE)
write.csv(missing_df, missing_csv, row.names = FALSE)

# -----------------------------
# Plot
# -----------------------------

cap <- quantile(abs(heat_df$log2FC), 0.98, na.rm = TRUE)
if (!is.finite(cap) || cap <= 0) cap <- max(abs(heat_df$log2FC), na.rm = TRUE)
if (!is.finite(cap) || cap <= 0) cap <- 1

heat_df <- heat_df %>%
  mutate(log2FC_plot = pmax(pmin(log2FC, cap), -cap))

p <- ggplot(heat_df, aes(x = effect, y = gene, fill = log2FC_plot)) +
  geom_tile(colour = "white", linewidth = 0.18) +
  geom_text(aes(label = sig), size = 2.4) +
  scale_fill_gradient2(
    low = "#377EB8",
    mid = "white",
    high = "#E41A1C",
    midpoint = 0,
    #limits = c(-cap, cap),
    na.value = "grey88"
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = paste0("log2FC"),
    title = "Fibroblast T5224, ageing, and exercise effects",
    subtitle = "* FDR <= 0.05, ** FDR <= 0.01, *** FDR <= 0.001"
  ) +
  theme_bw(base_size = 8.5) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(
  out_pdf,
  p,
  width = 6.5,
  height = max(6, length(genes) * 0.14),
  useDingbats = FALSE
)

message("PDF written to: ", out_pdf)
message("Values written to: ", out_csv)
message("Missing values written to: ", missing_csv)
