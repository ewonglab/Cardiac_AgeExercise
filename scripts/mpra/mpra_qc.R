# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Perform quality control and exploratory visualisation of MPRA measurements.
# Steps: 1) Load replicated reporter tables; 2) calculate QC summaries; 3) export diagnostic figures used to assess the library.
# =======================================
# Source reformatted from MPRA/MPRA_QC.ipynb.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

# MPRA QC PLots
MPRA = read.csv(publication_path("MPRA_ROOT", "results/counts/MPRA_Results_log2Alpha_replicated_Peak.csv"))
MPRA_fullShuffle = read.csv(publication_path("MPRA_ROOT", "results/counts/MPRA_Comparative_Results_FullShuffle.csv"))
MPRA_motifShuffle  = read.csv(publication_path("MPRA_ROOT", "results/counts/MPRA_Comparative_Results_MotifShuffle.csv"))

library(ggplot2)
library(scales)  # Added: comma / label_percent are used before the later library(scales) call
# Calculate MPRA activity scores (RNA/DNA ratios)
# Usage: Rscript calculate_activity.R <count_matrix.tsv> <samples.csv> <output_dir>
theme_mpra <- theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(
      size = 12,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 9,
      hjust = 0,
      margin = margin(b = 6)
    ),
    axis.title = element_text(
      size = 9
    ),
    axis.text = element_text(
      size = 7,
      colour = "black"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    strip.text = element_text(
      size = 8,
      face = "bold"
    ),
    strip.background = element_rect(
      fill = "grey92",
      colour = "black",
      linewidth = 0.25
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      colour = "black",
      linewidth = 0.3
    ),
    legend.title = element_text(
      size = 8,
      face = "bold"
    ),
    legend.text = element_text(
      size = 8
    ),
    plot.margin = margin(8, 8, 8, 8)
  )
suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(tidyverse)
})

# Parse command line arguments
COUNT_MATRIX=publication_path("PAPER_OUTPUT_ROOT", "Paper_Tables/MPRA/count_matrix.tsv")


count_file <- publication_path("PAPER_OUTPUT_ROOT", "Paper_Tables/MPRA/count_matrix.tsv")
sample_file <- publication_path("PAPER_OUTPUT_ROOT", "Paper_Tables/MPRA/samples.csv")
output_dir <- publication_path("PAPER_OUTPUT_ROOT", "Paper_Tables/MPRA/")

cat("=== MPRA Activity Score Calculation ===\n")
cat("Count matrix:", count_file, "\n")
cat("Sample metadata:", sample_file, "\n")
cat("Output directory:", output_dir, "\n\n")

# Read count matrix
cat("Reading count matrix...\n")
counts <- fread(count_file, header = TRUE)
element_names <- counts$element
counts <- as.matrix(counts[, -1])

rownames(counts) <- element_names
counts <- round(counts / 2)
cat("  Elements:", nrow(counts), "\n")
cat("  Samples:", ncol(counts), "\n\n")

# Read sample metadata
cat("Reading sample metadata...\n")
samples <- fread(sample_file, header = TRUE)
samples$condition <- samples$type  # DNA or RNA

cat("  DNA samples:", sum(samples$type == "DNA"), "\n")
cat("  RNA samples:", sum(samples$type == "RNA"), "\n\n")

#Replicate information 
# --- 1. Define Batch/Species Metadata ---
# RNA mapping
samples_rna <- samples[samples$type == "RNA", ]
samples_rna$species <- factor(ifelse(1:nrow(samples_rna) <= 3, "Human", "Mouse"))
samples_rna$bio_rep <- factor(case_when(
  1:nrow(samples_rna) <= 3 ~ "Human_Rep1",
  1:nrow(samples_rna) %in% 4:6 ~ "Mouse_RepA",
  1:nrow(samples_rna) %in% 7:9 ~ "Mouse_RepB"
))
samples_rna
# DNA mapping 
# (Assuming DNA samples were also sequenced in these batches; 
# if DNA is a single pool, assign it a common baseline species)
samples_dna <- samples[samples$type == "DNA", ]
samples_dna$species <- factor(ifelse(1:nrow(samples_dna) <= 3, "Human", "Mouse")) 
samples_dna$bio_rep <- "DNA_Pool"
samples_dna


dna_cols <- samples_dna %>% filter(type == "DNA" & species == "Mouse") %>% pull(sample_id)
rna_cols <- samples_rna %>% filter(type == "RNA" & species == "Mouse") %>% pull(sample_id)

dna_mat <- counts %>% 
  as.data.frame() %>% 
  dplyr::select(all_of(dna_cols)) %>% 
  as.matrix()

rna_mat <- counts %>% 
  as.data.frame() %>% 
  dplyr::select(all_of(rna_cols)) %>% 
  as.matrix()


# MPRA QC PLots
MPRA = read.csv(publication_path("MPRA_ROOT", "results/counts/MPRA_Results_log2Alpha_replicated_Peak.csv"))
MPRA_fullShuffle = read.csv(publication_path("MPRA_ROOT", "results/counts/MPRA_Comparative_Results_FullShuffle.csv"))
MPRA_motifShuffle  = read.csv(publication_path("MPRA_ROOT", "results/counts/MPRA_Comparative_Results_MotifShuffle.csv"))



head(MPRA)
head(MPRA_fullShuffle)
head(MPRA_motifShuffle)

# ------------------------------------------------------------
# 3. Add metadata matching your current setup
# ------------------------------------------------------------

samples_rna <- samples %>%
  filter(type == "RNA") %>%
  mutate(
    species = factor(ifelse(row_number() <= 3, "Human", "Mouse"),
                     levels = c("Human", "Mouse")),
    bio_rep = factor(case_when(
      row_number() <= 3 ~ "Human_Rep1",
      row_number() %in% 4:6 ~ "Mouse_RepA",
      row_number() %in% 7:9 ~ "Mouse_RepB",
      TRUE ~ "Unassigned"
    ))
  )

samples_dna <- samples %>%
  filter(type == "DNA") %>%
  mutate(
    species = factor(ifelse(row_number() <= 3, "Human", "Mouse"),
                     levels = c("Human", "Mouse")),
    bio_rep = factor("DNA_Pool")
  )

samples_annotated <- bind_rows(samples_dna, samples_rna) %>%
  mutate(
    type = factor(type, levels = c("DNA", "RNA")),
    sample_label = paste(sample_id, species, type, bio_rep, sep = " | ")
  )

if (nrow(samples_rna) != 9) {
  warning("RNA sample count is not 9. Check that the Human/Mouse/replicate assignment still matches your experiment.")
}

# ------------------------------------------------------------
# 4. Restrict everything to Mouse samples only
# ------------------------------------------------------------

samples_mouse <- samples_annotated %>%
  filter(species == "Mouse") %>%
  mutate(
    species = droplevels(species),
    bio_rep = droplevels(bio_rep),
    type = droplevels(type),
    sample_id = factor(sample_id, levels = sample_id)
  )

mouse_cols <- samples_mouse %>%
  pull(sample_id) %>%
  as.character()

counts_mouse <- counts[, mouse_cols, drop = FALSE]

cat("Mouse-only QC dataset:\n")
cat("  Elements:", nrow(counts_mouse), "\n")
cat("  Samples:", ncol(counts_mouse), "\n")
cat("  DNA samples:", sum(samples_mouse$type == "DNA"), "\n")
cat("  RNA samples:", sum(samples_mouse$type == "RNA"), "\n\n")

# Explicit mouse DNA/RNA matrices
dna_cols_mouse <- samples_mouse %>%
  filter(type == "DNA") %>%
  pull(sample_id) %>%
  as.character()

rna_cols_mouse <- samples_mouse %>%
  filter(type == "RNA") %>%
  pull(sample_id) %>%
  as.character()

dna_mat_mouse <- counts_mouse[, dna_cols_mouse, drop = FALSE]
rna_mat_mouse <- counts_mouse[, rna_cols_mouse, drop = FALSE]

# ------------------------------------------------------------
# 5. Sample-level raw count QC: Mouse only
# ------------------------------------------------------------

sample_qc <- tibble(
  sample_id = colnames(counts_mouse),
  library_size = colSums(counts_mouse, na.rm = TRUE),
  n_detected_elements = colSums(counts_mouse > 0, na.rm = TRUE),
  pct_detected_elements = 100 * colSums(counts_mouse > 0, na.rm = TRUE) / nrow(counts_mouse),
  median_count_all_elements = apply(counts_mouse, 2, median, na.rm = TRUE),
  median_count_detected_elements = apply(
    counts_mouse,
    2,
    function(x) median(x[x > 0], na.rm = TRUE)
  ),
  pct_zero_elements = 100 * colSums(counts_mouse == 0, na.rm = TRUE) / nrow(counts_mouse)
) %>%
  left_join(samples_mouse, by = "sample_id") %>%
  mutate(
    sample_id = factor(sample_id, levels = samples_mouse$sample_id)
  )

# ------------------------------------------------------------
# 6. Long-form count table for distributions: Mouse only
# ------------------------------------------------------------

count_long <- as.data.table(counts_mouse, keep.rownames = "element") %>%
  melt(
    id.vars = "element",
    variable.name = "sample_id",
    value.name = "count"
  ) %>%
  as_tibble() %>%
  left_join(samples_mouse, by = "sample_id") %>%
  mutate(
    sample_id = factor(sample_id, levels = samples_mouse$sample_id),
    log10_count = log10(count + 1)
  )

# ------------------------------------------------------------
# 7. Coverage threshold summaries: Mouse only
# ------------------------------------------------------------

coverage_thresholds <- c(0, 1, 5, 10, 25, 50, 100)

coverage_summary <- map_dfr(
  coverage_thresholds,
  function(thresh) {
    tibble(
      sample_id = colnames(counts_mouse),
      threshold = thresh,
      n_elements_at_or_above = colSums(counts_mouse >= thresh, na.rm = TRUE),
      pct_elements_at_or_above = 100 * colSums(counts_mouse >= thresh, na.rm = TRUE) / nrow(counts_mouse)
    )
  }
) %>%
  left_join(samples_mouse, by = "sample_id") %>%
  mutate(
    sample_id = factor(sample_id, levels = samples_mouse$sample_id),
    threshold_label = paste0("≥", threshold)
  )

# ------------------------------------------------------------
# 8. Sample correlation and PCA from log2 CPM: Mouse only
# ------------------------------------------------------------

library_sizes <- colSums(counts_mouse, na.rm = TRUE)
cpm_mat <- t(t(counts_mouse) / library_sizes) * 1e6
log2_cpm_mat <- log2(cpm_mat + 1)

# Remove elements with zero variance across mouse samples before PCA
element_sd <- apply(log2_cpm_mat, 1, sd, na.rm = TRUE)
pca_keep <- is.finite(element_sd) & element_sd > 0

cat("PCA elements retained:", sum(pca_keep), "of", length(pca_keep), "\n")
cat("Zero-variance / invalid elements removed:", sum(!pca_keep), "\n")

log2_cpm_mat_pca <- log2_cpm_mat[pca_keep, , drop = FALSE]

pca <- prcomp(
  t(log2_cpm_mat_pca),
  center = TRUE,
  scale. = TRUE
)

pca_df <- as_tibble(
  pca$x[, 1:2, drop = FALSE],
  rownames = "sample_id"
) %>%
  left_join(samples_mouse, by = "sample_id")

pc_var <- 100 * summary(pca)$importance[2, 1:2]

cor_mat <- cor(
  log2_cpm_mat,
  method = "spearman",
  use = "pairwise.complete.obs"
)

cor_order <- hclust(as.dist(1 - cor_mat))$order
cor_sample_order <- colnames(cor_mat)[cor_order]

cor_df <- as.data.frame(cor_mat) %>%
  rownames_to_column("sample_x") %>%
  pivot_longer(
    cols = -sample_x,
    names_to = "sample_y",
    values_to = "spearman_rho"
  ) %>%
  mutate(
    sample_x = factor(sample_x, levels = rev(cor_sample_order)),
    sample_y = factor(sample_y, levels = cor_sample_order)
  )

# ------------------------------------------------------------
# 9. DNA/RNA aggregate representation: Mouse only
# ------------------------------------------------------------

dna_rna_representation <- tibble(
  element = rownames(counts_mouse),
  DNA_total = rowSums(dna_mat_mouse, na.rm = TRUE),
  RNA_total = rowSums(rna_mat_mouse, na.rm = TRUE),
  log10_DNA = log10(DNA_total + 1),
  log10_RNA = log10(RNA_total + 1)
)

dna_rna_corr <- dna_rna_representation %>%
  summarise(
    spearman_rho = cor(
      log10_DNA,
      log10_RNA,
      method = "spearman",
      use = "pairwise.complete.obs"
    )
  ) %>%
  mutate(
    label = paste0("Spearman \u03c1 = ", round(spearman_rho, 3))
  )

write.csv(
  dna_rna_corr,
  file.path(output_dir, "MPRA_Mouse_DNA_RNA_aggregate_representation_correlation.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. Update the plots that previously faceted by species
# ------------------------------------------------------------

p_library_size <- ggplot(
  sample_qc,
  aes(x = sample_id, y = library_size, fill = type)
) +
  geom_col(colour = "black", linewidth = 0.25) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = NULL,
    y = "Total counts",
    title = "Mouse MPRA library size per sample"
  ) +
  theme_mpra +
  theme(legend.position = "top")

p_detected_elements <- ggplot(
  sample_qc,
  aes(x = sample_id, y = pct_detected_elements, fill = type)
) +
  geom_col(colour = "black", linewidth = 0.25) +
  scale_y_continuous(
    labels = label_percent(scale = 1),
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = NULL,
    y = "Elements with count > 0",
    title = "Mouse MPRA element detection rate per sample"
  ) +
  theme_mpra +
  theme(legend.position = "top")

p_count_distribution <- ggplot(
  count_long,
  aes(x = sample_id, y = log10_count, fill = type)
) +
  geom_violin(scale = "width", linewidth = 0.2, trim = TRUE) +
  geom_boxplot(
    width = 0.10,
    outlier.shape = NA,
    linewidth = 0.2,
    fill = "white"
  ) +
  labs(
    x = NULL,
    y = expression(log[10] * "(count + 1)"),
    title = "Mouse MPRA raw count distributions"
  ) +
  theme_mpra +
  theme(legend.position = "top")

p_coverage_thresholds <- coverage_summary %>%
  filter(threshold %in% c(1, 5, 10, 25, 50, 100)) %>%
  ggplot(
    aes(
      x = threshold_label,
      y = pct_elements_at_or_above,
      group = sample_id,
      colour = type
    )
  ) +
  geom_line(alpha = 0.65, linewidth = 0.45) +
  geom_point(size = 1.2) +
  scale_y_continuous(labels = label_percent(scale = 1)) +
  labs(
    x = "Read-count threshold",
    y = "Elements retained",
    title = "Mouse MPRA element retention across minimum count thresholds"
  ) +
  theme_mpra +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

p_pca <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    shape = type,
    colour = bio_rep
  )
) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text(
    aes(label = sample_id),
    size = 2.5,
    vjust = -0.8,
    show.legend = FALSE
  ) +
  labs(
    x = paste0("PC1 (", round(pc_var[1], 1), "%)"),
    y = paste0("PC2 (", round(pc_var[2], 1), "%)"),
    colour = "Mouse group",
    shape = "Molecule",
    title = "PCA of mouse MPRA log2 CPM counts"
  ) +
  theme_mpra +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

p_dna_rna_representation <- ggplot(
  dna_rna_representation,
  aes(x = log10_DNA, y = log10_RNA)
) +
  geom_point(size = 0.35, alpha = 0.20) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, colour = "black") +
  geom_text(
    data = dna_rna_corr,
    aes(
      x = -Inf,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.2,
    size = 3
  ) +
  labs(
    x = expression(log[10] * "(aggregated mouse DNA counts + 1)"),
    y = expression(log[10] * "(aggregated mouse RNA counts + 1)"),
    title = "Mouse MPRA DNA input versus RNA output representation"
  ) +
  theme_mpra +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# ============================================================
# Mouse-only MPRA QC plots and summaries
# - Restricts raw count QC to Mouse DNA/RNA samples only
# - Prints all plots
# - Saves all plots as PDFs only
# - Adds class retention plot:
#   Age Up, Age Down, Ex Up, Ex Down, Stable
#   passing all filters, defined as replicated == TRUE
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------

count_file <- publication_path("PAPER_OUTPUT_ROOT", "Paper_Tables/MPRA/count_matrix.tsv")
sample_file <- publication_path("PAPER_OUTPUT_ROOT", "Paper_Tables/MPRA/samples.csv")
output_dir <- publication_path("PAPER_OUTPUT_ROOT", "Paper_Tables/MPRA/")

mpra_results_file <- publication_path("MPRA_ROOT", "results/counts/MPRA_Results_log2Alpha_replicated_Peak.csv")
fullshuffle_file <- publication_path("MPRA_ROOT", "results/counts/MPRA_Comparative_Results_FullShuffle.csv")
motifshuffle_file <- publication_path("MPRA_ROOT", "results/counts/MPRA_Comparative_Results_MotifShuffle.csv")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Plot theme
# ------------------------------------------------------------

theme_mpra <- theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(
      size = 12,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 9,
      hjust = 0,
      margin = margin(b = 6)
    ),
    axis.title = element_text(
      size = 9
    ),
    axis.text = element_text(
      size = 7,
      colour = "black"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    strip.text = element_text(
      size = 8,
      face = "bold"
    ),
    strip.background = element_rect(
      fill = "grey92",
      colour = "black",
      linewidth = 0.25
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      colour = "black",
      linewidth = 0.3
    ),
    legend.title = element_text(
      size = 8,
      face = "bold"
    ),
    legend.text = element_text(
      size = 8
    ),
    plot.margin = margin(8, 8, 8, 8)
  )

# ------------------------------------------------------------
# 3. Read count matrix
# ------------------------------------------------------------

cat("Reading count matrix...\n")

counts_dt <- fread(count_file, header = TRUE)

element_names <- counts_dt$element
counts <- as.matrix(counts_dt[, -1])
rownames(counts) <- element_names

# Keep your existing preprocessing
counts <- round(counts / 2)

cat("  Elements:", nrow(counts), "\n")
cat("  Samples:", ncol(counts), "\n\n")

# ------------------------------------------------------------
# 4. Read sample metadata
# ------------------------------------------------------------

cat("Reading sample metadata...\n")

samples <- fread(sample_file, header = TRUE) %>%
  as_tibble()

stopifnot("sample_id" %in% colnames(samples))
stopifnot("type" %in% colnames(samples))
stopifnot(all(samples$sample_id %in% colnames(counts)))

# Match sample metadata to count matrix order
samples <- samples %>%
  mutate(sample_id = factor(sample_id, levels = colnames(counts))) %>%
  arrange(sample_id) %>%
  mutate(sample_id = as.character(sample_id))

cat("  DNA samples:", sum(samples$type == "DNA"), "\n")
cat("  RNA samples:", sum(samples$type == "RNA"), "\n\n")

# ------------------------------------------------------------
# 5. Add species and biological replicate metadata
# ------------------------------------------------------------

samples_rna <- samples %>%
  filter(type == "RNA") %>%
  mutate(
    species = factor(
      ifelse(row_number() <= 3, "Human", "Mouse"),
      levels = c("Human", "Mouse")
    ),
    bio_rep = factor(case_when(
      row_number() <= 3 ~ "Human_Rep1",
      row_number() %in% 4:6 ~ "Mouse_RepA",
      row_number() %in% 7:9 ~ "Mouse_RepB",
      TRUE ~ "Unassigned"
    ))
  )

samples_dna <- samples %>%
  filter(type == "DNA") %>%
  mutate(
    species = factor(
      ifelse(row_number() <= 3, "Human", "Mouse"),
      levels = c("Human", "Mouse")
    ),
    bio_rep = factor("DNA_Pool")
  )

samples_annotated <- bind_rows(samples_dna, samples_rna) %>%
  mutate(
    type = factor(type, levels = c("DNA", "RNA")),
    sample_label = paste(sample_id, species, type, bio_rep, sep = " | ")
  )

if (nrow(samples_rna) != 9) {
  warning(
    "RNA sample count is not 9. Check that Human/Mouse/replicate assignment still matches your experiment."
  )
}

# ------------------------------------------------------------
# 6. Restrict raw QC to Mouse samples only
# ------------------------------------------------------------

samples_mouse <- samples_annotated %>%
  filter(species == "Mouse") %>%
  mutate(
    sample_id = as.character(sample_id),
    species = droplevels(species),
    bio_rep = droplevels(bio_rep),
    type = droplevels(type)
  )

mouse_sample_levels <- samples_mouse$sample_id

mouse_cols <- samples_mouse %>%
  pull(sample_id)

counts_mouse <- counts[, mouse_cols, drop = FALSE]

cat("Mouse-only QC dataset:\n")
cat("  Elements:", nrow(counts_mouse), "\n")
cat("  Samples:", ncol(counts_mouse), "\n")
cat("  Mouse DNA samples:", sum(samples_mouse$type == "DNA"), "\n")
cat("  Mouse RNA samples:", sum(samples_mouse$type == "RNA"), "\n\n")

# Explicit Mouse DNA/RNA matrices
dna_cols_mouse <- samples_mouse %>%
  filter(type == "DNA") %>%
  pull(sample_id)

rna_cols_mouse <- samples_mouse %>%
  filter(type == "RNA") %>%
  pull(sample_id)

dna_mat_mouse <- counts_mouse[, dna_cols_mouse, drop = FALSE]
rna_mat_mouse <- counts_mouse[, rna_cols_mouse, drop = FALSE]

# ------------------------------------------------------------
# 7. Mouse sample-level raw count QC
# ------------------------------------------------------------

sample_qc <- tibble(
  sample_id = colnames(counts_mouse),
  library_size = colSums(counts_mouse, na.rm = TRUE),
  n_detected_elements = colSums(counts_mouse > 0, na.rm = TRUE),
  pct_detected_elements = 100 * colSums(counts_mouse > 0, na.rm = TRUE) / nrow(counts_mouse),
  median_count_all_elements = apply(counts_mouse, 2, median, na.rm = TRUE),
  median_count_detected_elements = apply(
    counts_mouse,
    2,
    function(x) {
      detected <- x[x > 0]
      if (length(detected) == 0) {
        return(NA_real_)
      } else {
        return(median(detected, na.rm = TRUE))
      }
    }
  ),
  pct_zero_elements = 100 * colSums(counts_mouse == 0, na.rm = TRUE) / nrow(counts_mouse)
) %>%
  left_join(samples_mouse, by = "sample_id") %>%
  mutate(
    sample_id = factor(sample_id, levels = mouse_sample_levels)
  )

# ------------------------------------------------------------
# 8. Long-form mouse count table for distributions
# ------------------------------------------------------------

count_long <- as.data.table(counts_mouse, keep.rownames = "element") %>%
  melt(
    id.vars = "element",
    variable.name = "sample_id",
    value.name = "count"
  ) %>%
  as_tibble() %>%
  left_join(samples_mouse, by = "sample_id") %>%
  mutate(
    sample_id = factor(sample_id, levels = mouse_sample_levels),
    log10_count = log10(count + 1)
  )

# ------------------------------------------------------------
# 9. Coverage threshold summaries: Mouse only
# ------------------------------------------------------------

coverage_thresholds <- c(0, 1, 5, 10, 25, 50, 100)

coverage_summary <- map_dfr(
  coverage_thresholds,
  function(thresh) {
    tibble(
      sample_id = colnames(counts_mouse),
      threshold = thresh,
      n_elements_at_or_above = colSums(counts_mouse >= thresh, na.rm = TRUE),
      pct_elements_at_or_above = 100 * colSums(counts_mouse >= thresh, na.rm = TRUE) / nrow(counts_mouse)
    )
  }
) %>%
  left_join(samples_mouse, by = "sample_id") %>%
  mutate(
    sample_id = factor(sample_id, levels = mouse_sample_levels),
    threshold_label = paste0("≥", threshold)
  )

# ------------------------------------------------------------
# 10. Mouse sample correlation and PCA from log2 CPM
# ------------------------------------------------------------

library_sizes_mouse <- colSums(counts_mouse, na.rm = TRUE)

cpm_mat_mouse <- t(t(counts_mouse) / library_sizes_mouse) * 1e6
log2_cpm_mat_mouse <- log2(cpm_mat_mouse + 1)

# Correlation matrix
cor_mat_mouse <- cor(
  log2_cpm_mat_mouse,
  method = "spearman",
  use = "pairwise.complete.obs"
)

cor_order_mouse <- hclust(as.dist(1 - cor_mat_mouse))$order
cor_sample_order_mouse <- colnames(cor_mat_mouse)[cor_order_mouse]

cor_df <- as.data.frame(cor_mat_mouse) %>%
  rownames_to_column("sample_x") %>%
  pivot_longer(
    cols = -sample_x,
    names_to = "sample_y",
    values_to = "spearman_rho"
  ) %>%
  mutate(
    sample_x = factor(sample_x, levels = rev(cor_sample_order_mouse)),
    sample_y = factor(sample_y, levels = cor_sample_order_mouse)
  )

# PCA: remove zero-variance features first
element_sd_mouse <- apply(log2_cpm_mat_mouse, 1, sd, na.rm = TRUE)

pca_keep_mouse <- is.finite(element_sd_mouse) & element_sd_mouse > 0

cat("PCA elements retained:", sum(pca_keep_mouse), "of", length(pca_keep_mouse), "\n")
cat("Zero-variance / invalid elements removed:", sum(!pca_keep_mouse), "\n\n")

if (sum(pca_keep_mouse) < 2) {
  stop("Fewer than 2 non-zero variance MPRA elements remain for PCA.")
}

log2_cpm_mat_mouse_pca <- log2_cpm_mat_mouse[pca_keep_mouse, , drop = FALSE]

pca_mouse <- prcomp(
  t(log2_cpm_mat_mouse_pca),
  center = TRUE,
  scale. = TRUE
)

pca_df <- as_tibble(
  pca_mouse$x[, 1:2, drop = FALSE],
  rownames = "sample_id"
) %>%
  left_join(samples_mouse, by = "sample_id")

pc_var_mouse <- 100 * summary(pca_mouse)$importance[2, 1:2]

# ------------------------------------------------------------
# 11. Mouse DNA/RNA aggregate representation
# ------------------------------------------------------------

dna_rna_representation <- tibble(
  element = rownames(counts_mouse),
  DNA_total = rowSums(dna_mat_mouse, na.rm = TRUE),
  RNA_total = rowSums(rna_mat_mouse, na.rm = TRUE),
  log10_DNA = log10(DNA_total + 1),
  log10_RNA = log10(RNA_total + 1)
)

dna_rna_corr <- dna_rna_representation %>%
  summarise(
    spearman_rho = cor(
      log10_DNA,
      log10_RNA,
      method = "spearman",
      use = "pairwise.complete.obs"
    )
  ) %>%
  mutate(
    label = paste0("Spearman \u03c1 = ", round(spearman_rho, 3))
  )

# ------------------------------------------------------------
# 12. Read fitted mouse MPRA result tables
# ------------------------------------------------------------

MPRA <- read.csv(mpra_results_file, stringsAsFactors = FALSE)
MPRA_fullShuffle <- read.csv(fullshuffle_file, stringsAsFactors = FALSE)
MPRA_motifShuffle <- read.csv(motifshuffle_file, stringsAsFactors = FALSE)

# ------------------------------------------------------------
# 13. MPRA fitted activity score distributions
# ------------------------------------------------------------

activity_long <- MPRA %>%
  select(
    Coords,
    category,
    Group,
    replicated,
    log2_alpha_Endo,
    log2_alpha_MotifShuff,
    log2_alpha_FullShuff
  ) %>%
  pivot_longer(
    cols = starts_with("log2_alpha_"),
    names_to = "Sequence_Type",
    values_to = "log2_alpha"
  ) %>%
  mutate(
    Sequence_Type = recode(
      Sequence_Type,
      "log2_alpha_Endo" = "Endogenous",
      "log2_alpha_MotifShuff" = "Motif shuffle",
      "log2_alpha_FullShuff" = "Full shuffle"
    ),
    Sequence_Type = factor(
      Sequence_Type,
      levels = c("Endogenous", "Motif shuffle", "Full shuffle")
    )
  )

# ------------------------------------------------------------
# 14. Shuffle comparison tables
# ------------------------------------------------------------

shuffle_scatter <- bind_rows(
  MPRA_fullShuffle %>%
    transmute(
      Coords,
      category,
      Group,
      enhancer,
      is_significant,
      Comparison = "Endogenous vs Full shuffle",
      log2_alpha_Endo,
      log2_alpha_Shuffle = log2_alpha_FullShuff,
      logFC = logFC_Endo_FullShuff,
      pval,
      fdr
    ),
  MPRA_motifShuffle %>%
    transmute(
      Coords,
      category,
      Group,
      enhancer,
      is_significant,
      Comparison = "Endogenous vs Motif shuffle",
      log2_alpha_Endo,
      log2_alpha_Shuffle = log2_alpha_MotifShuff,
      logFC = logFC_Endo_MotifShuff,
      pval,
      fdr
    )
) %>%
  mutate(
    Comparison = factor(
      Comparison,
      levels = c("Endogenous vs Motif shuffle", "Endogenous vs Full shuffle")
    ),
    is_significant = factor(
      is_significant,
      levels = c("Non-Significant", "Significant")
    )
  )

shuffle_corr <- shuffle_scatter %>%
  group_by(Comparison) %>%
  summarise(
    spearman_rho = cor(
      log2_alpha_Endo,
      log2_alpha_Shuffle,
      method = "spearman",
      use = "pairwise.complete.obs"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("Spearman \u03c1 = ", round(spearman_rho, 3))
  )

shuffle_significance_summary <- shuffle_scatter %>%
  count(Comparison, category, is_significant, name = "n_elements") %>%
  group_by(Comparison, category) %>%
  mutate(
    fraction = n_elements / sum(n_elements)
  ) %>%
  ungroup()

# ------------------------------------------------------------
# 15. Replicated element summary
# ------------------------------------------------------------

replicated_summary <- MPRA %>%
  count(category, Group, replicated, name = "n_elements") %>%
  group_by(category, Group) %>%
  mutate(
    fraction = n_elements / sum(n_elements)
  ) %>%
  ungroup()

# ------------------------------------------------------------
# 16. Sequence-class retention summary
#
# Counts how many regulatory-class sequences pass all filters.
# "Pass all filters" = replicated == TRUE in MPRA table.
#
# Contrast classes are intentionally non-exclusive:
# a sequence can be both Age Up and Ex Down, for example.
# ------------------------------------------------------------

passed_filter_coords <- MPRA %>%
  filter(replicated == TRUE) %>%
  distinct(Coords)

sequence_class_membership <- MPRA_fullShuffle %>%
  distinct(Coords, in_YS, in_MAS, in_MAE) %>%
  mutate(
    `Age Up` = !in_YS & in_MAS,
    `Age Down` = in_YS & !in_MAS,
    `Ex Up` = !in_MAS & in_MAE,
    `Ex Down` = in_MAS & !in_MAE,
    `Stable` = in_YS & in_MAS & in_MAE
  ) %>%
  pivot_longer(
    cols = c(`Age Up`, `Age Down`, `Ex Up`, `Ex Down`, `Stable`),
    names_to = "Sequence_Class",
    values_to = "In_Class"
  ) %>%
  filter(In_Class) %>%
  select(Coords, Sequence_Class)

class_filter_summary <- sequence_class_membership %>%
  mutate(
    Passed_All_Filters = Coords %in% passed_filter_coords$Coords
  ) %>%
  group_by(Sequence_Class) %>%
  summarise(
    n_total_in_class = n_distinct(Coords),
    n_passed_all_filters = n_distinct(Coords[Passed_All_Filters]),
    fraction_passed = n_passed_all_filters / n_total_in_class,
    .groups = "drop"
  ) %>%
  mutate(
    Sequence_Class = factor(
      Sequence_Class,
      levels = c("Age Up", "Age Down", "Ex Up", "Ex Down", "Stable")
    )
  ) %>%
  arrange(Sequence_Class)

print(class_filter_summary)

# ------------------------------------------------------------
# 17. Plot objects
# ------------------------------------------------------------

p_library_size <- ggplot(
  sample_qc,
  aes(x = sample_id, y = library_size, fill = type)
) +
  geom_col(
    colour = "black",
    linewidth = 0.25
  ) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = NULL,
    y = "Total counts",
    title = "Mouse MPRA library size per sample"
  ) +
  theme_mpra +
  theme(
    legend.position = "top"
  )

p_detected_elements <- ggplot(
  sample_qc,
  aes(x = sample_id, y = pct_detected_elements, fill = type)
) +
  geom_col(
    colour = "black",
    linewidth = 0.25
  ) +
  scale_y_continuous(
    labels = label_percent(scale = 1),
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = NULL,
    y = "Elements with count > 0",
    title = "Mouse MPRA element detection rate per sample"
  ) +
  theme_mpra +
  theme(
    legend.position = "top"
  )

p_count_distribution <- ggplot(
  count_long,
  aes(x = sample_id, y = log10_count, fill = type)
) +
  geom_violin(
    scale = "width",
    linewidth = 0.2,
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.10,
    outlier.shape = NA,
    linewidth = 0.2,
    fill = "white"
  ) +
  labs(
    x = NULL,
    y = expression(log[10] * "(count + 1)"),
    title = "Mouse MPRA raw count distributions"
  ) +
  theme_mpra +
  theme(
    legend.position = "top"
  )

p_coverage_thresholds <- coverage_summary %>%
  filter(threshold %in% c(1, 5, 10, 25, 50, 100)) %>%
  ggplot(
    aes(
      x = threshold_label,
      y = pct_elements_at_or_above,
      group = sample_id,
      colour = type
    )
  ) +
  geom_line(
    alpha = 0.65,
    linewidth = 0.45
  ) +
  geom_point(
    size = 1.2
  ) +
  scale_y_continuous(
    labels = label_percent(scale = 1)
  ) +
  labs(
    x = "Read-count threshold",
    y = "Elements retained",
    title = "Mouse MPRA element retention across minimum count thresholds"
  ) +
  theme_mpra +
  theme(
    legend.position = "top",
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    )
  )

p_sample_correlation <- ggplot(
  cor_df,
  aes(x = sample_y, y = sample_x, fill = spearman_rho)
) +
  geom_tile(
    colour = "white",
    linewidth = 0.2
  ) +
  geom_text(
    aes(label = sprintf("%.2f", spearman_rho)),
    size = 2.2,
    colour = "black"
  ) +
  scale_fill_gradient(
    low = "white",
    high = "grey25",
    limits = c(min(cor_df$spearman_rho, na.rm = TRUE), 1)
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Spearman\ncorrelation",
    title = "Mouse MPRA sample-to-sample correlation of log2 CPM counts"
  ) +
  theme_mpra +
  theme(
    legend.position = "right",
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.text.y = element_text(
      angle = 0
    )
  )

p_pca <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    shape = type,
    colour = bio_rep
  )
) +
  geom_point(
    size = 3,
    alpha = 0.9
  ) +
  geom_text(
    aes(label = sample_id),
    size = 2.5,
    vjust = -0.8,
    show.legend = FALSE
  ) +
  labs(
    x = paste0("PC1 (", round(pc_var_mouse[1], 1), "%)"),
    y = paste0("PC2 (", round(pc_var_mouse[2], 1), "%)"),
    colour = "Mouse group",
    shape = "Molecule",
    title = "PCA of mouse MPRA log2 CPM counts"
  ) +
  theme_mpra +
  theme(
    legend.position = "right",
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    )
  )

p_dna_rna_representation <- ggplot(
  dna_rna_representation,
  aes(x = log10_DNA, y = log10_RNA)
) +
  geom_point(
    size = 0.35,
    alpha = 0.20
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_text(
    data = dna_rna_corr,
    aes(
      x = -Inf,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.2,
    size = 3
  ) +
  labs(
    x = expression(log[10] * "(aggregated mouse DNA counts + 1)"),
    y = expression(log[10] * "(aggregated mouse RNA counts + 1)"),
    title = "Mouse MPRA DNA input versus RNA output representation"
  ) +
  theme_mpra +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    )
  )

p_activity_distribution <- ggplot(
  activity_long,
  aes(x = Sequence_Type, y = log2_alpha, fill = Sequence_Type)
) +
  geom_violin(
    trim = TRUE,
    scale = "width",
    linewidth = 0.25
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    linewidth = 0.25,
    fill = "white"
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.25,
    linetype = "dashed"
  ) +
  labs(
    x = NULL,
    y = expression(log[2] * " alpha"),
    title = "Distribution of fitted mouse MPRA activity scores"
  ) +
  theme_mpra +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    )
  )

p_replicated_fraction <- ggplot(
  replicated_summary,
  aes(x = category, y = fraction, fill = replicated)
) +
  geom_col(
    colour = "black",
    linewidth = 0.2
  ) +
  facet_wrap(
    ~ Group,
    scales = "free_x"
  ) +
  scale_y_continuous(
    labels = percent_format()
  ) +
  labs(
    x = NULL,
    y = "Fraction of elements",
    fill = "Replicated",
    title = "Fraction of mouse MPRA elements classified as replicated"
  ) +
  theme_mpra +
  theme(
    legend.position = "top",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

p_shuffle_activity_scatter <- ggplot(
  shuffle_scatter,
  aes(
    x = log2_alpha_Shuffle,
    y = log2_alpha_Endo,
    colour = is_significant
  )
) +
  geom_point(
    size = 0.45,
    alpha = 0.45
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linewidth = 0.35,
    linetype = "dashed",
    colour = "black"
  ) +
  geom_text(
    data = shuffle_corr,
    aes(
      x = -Inf,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.2,
    size = 3
  ) +
  facet_wrap(
    ~ Comparison,
    nrow = 1
  ) +
  scale_colour_manual(
    values = c(
      "Non-Significant" = "grey70",
      "Significant" = "black"
    ),
    drop = FALSE
  ) +
  labs(
    x = expression(log[2] * " alpha: shuffled sequence"),
    y = expression(log[2] * " alpha: endogenous sequence"),
    colour = NULL,
    title = "Mouse MPRA endogenous versus shuffled activity"
  ) +
  theme_mpra +
  theme(
    legend.position = "top",
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    )
  )

p_shuffle_logfc_distribution <- ggplot(
  shuffle_scatter,
  aes(x = logFC, fill = Comparison)
) +
  geom_histogram(
    bins = 70,
    colour = "black",
    linewidth = 0.15,
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.3,
    linetype = "dashed"
  ) +
  facet_wrap(
    ~ Comparison,
    nrow = 1,
    scales = "free_y"
  ) +
  labs(
    x = expression(log[2] * " fold change: Endogenous - Shuffle"),
    y = "Number of elements",
    title = "Mouse MPRA endogenous-versus-shuffle effect distributions"
  ) +
  theme_mpra +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    )
  )

p_shuffle_pvalue_hist <- ggplot(
  shuffle_scatter,
  aes(x = pval)
) +
  geom_histogram(
    bins = 50,
    fill = "grey65",
    colour = "black",
    linewidth = 0.2
  ) +
  facet_wrap(
    ~ Comparison,
    nrow = 1,
    scales = "free_y"
  ) +
  labs(
    x = "Nominal P value",
    y = "Number of elements",
    title = "Mouse MPRA P-value distributions for shuffle comparisons"
  ) +
  theme_mpra +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    )
  )

p_shuffle_significance_fraction <- ggplot(
  shuffle_significance_summary,
  aes(x = category, y = fraction, fill = is_significant)
) +
  geom_col(
    colour = "black",
    linewidth = 0.2
  ) +
  facet_wrap(
    ~ Comparison,
    nrow = 1,
    scales = "free_x"
  ) +
  scale_y_continuous(
    labels = percent_format()
  ) +
  scale_fill_manual(
    values = c(
      "Non-Significant" = "grey75",
      "Significant" = "grey20"
    ),
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "Fraction of elements",
    fill = NULL,
    title = "Fraction of significant mouse MPRA endogenous-versus-shuffle effects"
  ) +
  theme_mpra +
  theme(
    legend.position = "top",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

class_colours <- c(
  "Age Up" = "#E41A1C",
  "Age Down" = "#377EB8",
  "Ex Up" = "#A020F0",
  "Ex Down" = "#00008B",
  "Stable" = "grey60"
)

p_class_filter_retention <- ggplot(
  class_filter_summary,
  aes(
    x = Sequence_Class,
    y = n_passed_all_filters,
    fill = Sequence_Class
  )
) +
  geom_col(
    width = 0.72,
    colour = "black",
    linewidth = 0.25
  ) +
  geom_text(
    aes(
      label = paste0(
        comma(n_passed_all_filters),
        "\n",
        "(",
        percent(fraction_passed, accuracy = 0.1),
        ")"
      )
    ),
    vjust = -0.35,
    size = 3
  ) +
  scale_fill_manual(
    values = class_colours,
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.14))
  ) +
  labs(
    x = NULL,
    y = "Sequences passing all filters",
    title = "Mouse MPRA sequences retained after filtering by regulatory class",
    subtitle = "Retained = replicated == TRUE; classes are non-exclusive contrast classes"
  ) +
  theme_mpra +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    )
  )

# ------------------------------------------------------------
# 18. Print every plot
# ------------------------------------------------------------

plots_to_print <- list(
  p_library_size,
  p_detected_elements,
  p_count_distribution,
  p_coverage_thresholds,
  p_sample_correlation,
  p_pca,
  p_dna_rna_representation,
  p_activity_distribution,
  p_replicated_fraction,
  p_shuffle_activity_scatter,
  p_shuffle_logfc_distribution,
  p_shuffle_pvalue_hist,
  p_shuffle_significance_fraction,
  p_class_filter_retention
)

for (p in plots_to_print) {
  print(p)
}

# ------------------------------------------------------------
# 19. Save every plot as an individual PDF
# ------------------------------------------------------------

ggsave(
  filename = file.path(output_dir, "MPRA_Mouse_QC_01_library_sizes.pdf"),
  plot = p_library_size,
  width = 8.27,
  height = 4.5
)

ggsave(
  filename = file.path(output_dir, "MPRA_Mouse_QC_02_detected_elements.pdf"),
  plot = p_detected_elements,
  width = 8.27,
  height = 4.5
)

ggsave(
  filename = file.path(output_dir, "MPRA_Mouse_QC_03_count_distributions.pdf"),
  plot = p_count_distribution,
  width = 8.27,
  height = 4.8
)

ggsave(
  filename = file.path(output_dir, "MPRA_Mouse_QC_04_coverage_thresholds.pdf"),
  plot = p_coverage_thresholds,
  width = 8.27,
  height = 4.8
)

ggsave(
  filename = file.path(output_dir, "MPRA_Mouse_QC_05_sample_correlation_heatmap.pdf"),
  plot = p_sample_correlation,
  width = 8.27,
  height = 7.5
)

ggsave(
  filename = file.path(output_dir, "MPRA_Mouse_QC_06_PCA.pdf"),
  plot = p_pca,
  width = 6.5,
  height = 5.2
)

ggsave(
  filename = file.path(output_dir, "MPRA_Mouse_QC_07_DNA_vs_RNA_representation.pdf"),
  plot = p_dna_rna_representation,
  width = 6.5,
  height = 5.0
)

ggsave(
  filename = file.path(output_dir, "MPRA_Mouse_QC_08_activity_distributions.pdf"),
  plot = p_activity_distribution,
  width = 6.5,
  height = 4.5
)

ggsave(
  filename = file.path(output_dir, "MPRA_Mouse_QC_09_replicated_fraction.pdf"),
  plot = p_replicated_fraction,
  width = 8.27,
  height = 5.5
)

ggsave(
  filename = file.path(output_dir, "MPRA_Mouse_QC_10_endogenous_vs_shuffle_scatter.pdf"),
  plot = p_shuffle_activity_scatter,
  width = 8.27,
  height = 4.8
)

ggsave(
  filename = file.path(output_dir, "MPRA_Mouse_QC_11_shuffle_logFC_distributions.pdf"),
  plot = p_shuffle_logfc_distribution,
  width = 8.27,
  height = 4.8
)

ggsave(
  filename = file.path(output_dir, "MPRA_Mouse_QC_12_shuffle_Pvalue_histograms.pdf"),
  plot = p_shuffle_pvalue_hist,
  width = 8.27,
  height = 4.8
)

ggsave(
  filename = file.path(output_dir, "MPRA_Mouse_QC_13_shuffle_significance_fraction.pdf"),
  plot = p_shuffle_significance_fraction,
  width = 8.27,
  height = 4.8
)

ggsave(
  filename = file.path(output_dir, "MPRA_Mouse_QC_14_sequence_class_filter_retention.pdf"),
  plot = p_class_filter_retention,
  width = 7.2,
  height = 5.2
)

# ------------------------------------------------------------
# 20. Also save one multipage PDF containing all plots
# ------------------------------------------------------------

all_plots_pdf <- file.path(output_dir, "MPRA_Mouse_QC_All_Plots.pdf")

pdf(
  file = all_plots_pdf,
  width = 8.27,
  height = 11.69,
  onefile = TRUE
)

for (p in plots_to_print) {
  print(p)
}

dev.off()

cat("\nAll mouse MPRA QC plots saved as PDFs in:\n", output_dir, "\n")
cat("\nMultipage summary PDF:\n", all_plots_pdf, "\n")

# ============================================================
# Supplementary mouse MPRA QC figure
# Library size, detection rate, and raw count distributions
# arranged in one row and saved as a single PDF
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------------------
# 1. Slightly tighten individual panels for supplementary figure
# ------------------------------------------------------------

p_library_size_supp <- p_library_size +
  labs(
    title = "Library size"
  ) +
  theme_mpra +
  theme(
    plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6.5, colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )

p_detected_elements_supp <- p_detected_elements +
  labs(
    title = "Element detection rate"
  ) +
  theme_mpra +
  theme(
    plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6.5, colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )

p_count_distribution_supp <- p_count_distribution +
  labs(
    title = "Raw count distribution"
  ) +
  theme_mpra +
  theme(
    plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6.5, colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    plot.margin = margin(5, 5, 5, 5)
  )

# ------------------------------------------------------------
# 2. Combine into one-row supplementary figure
# ------------------------------------------------------------

p_mouse_mpra_basic_qc_supp <- (
  p_library_size_supp |
    p_detected_elements_supp |
    p_count_distribution_supp
) +
  plot_layout(
    ncol = 3,
    widths = c(1, 1, 1.15),
    guides = "collect"
  ) &
  theme(
    legend.position = "top"
  )

# ------------------------------------------------------------
# 3. Print figure
# ------------------------------------------------------------

print(p_mouse_mpra_basic_qc_supp)

# ------------------------------------------------------------
# 4. Save as PDF
# Fixed A4 portrait-width export, suitable for a supplementary panel
# ------------------------------------------------------------

ggsave(
  filename = file.path(output_dir, "MPRA_Mouse_QC_Supplementary_Basic_Count_QC.pdf"),
  plot = p_mouse_mpra_basic_qc_supp,
  width = 8.27,
  height = 3.45
)



suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
})

# ------------------------------------------------------------
# Define binding classes
# ------------------------------------------------------------

class_defs <- list(
  Stable = function(df) {
    df %>% filter(in_YS == TRUE, in_MAS == TRUE, in_MAE == TRUE)
  },
  Age_up = function(df) {
    df %>% filter(in_YS == FALSE, in_MAS == TRUE)
  },
  Age_down = function(df) {
    df %>% filter(in_YS == TRUE, in_MAS == FALSE)
  },
  Ex_up = function(df) {
    df %>% filter(in_MAS == FALSE, in_MAE == TRUE)
  },
  Ex_down = function(df) {
    df %>% filter(in_MAS == TRUE, in_MAE == FALSE)
  }
)

dynamic_classes <- c("Age_up", "Age_down", "Ex_up", "Ex_down")

# ------------------------------------------------------------
# 1. Wilcoxon tests:
#    log2_alpha_Endo in ACTIVE enhancers only
# ------------------------------------------------------------

stable_active <- class_defs$Stable(MPRA_fullShuffle) %>%
  filter(enhancer == "Active") %>%
  mutate(binding_class = "Stable")

wilcox_results <- map_dfr(dynamic_classes, function(cls) {
  
  dynamic_active <- class_defs[[cls]](MPRA_fullShuffle) %>%
    filter(enhancer == "Active") %>%
    mutate(binding_class = cls)
  
  test_df <- bind_rows(stable_active, dynamic_active) %>%
    filter(!is.na(log2_alpha_Endo))
  
  wt <- wilcox.test(
    log2_alpha_Endo ~ binding_class,
    data = test_df,
    exact = FALSE
  )
  
  tibble(
    comparison = paste(cls, "vs Stable"),
    dynamic_class = cls,
    n_stable = sum(test_df$binding_class == "Stable"),
    n_dynamic = sum(test_df$binding_class == cls),
    median_stable = median(test_df$log2_alpha_Endo[test_df$binding_class == "Stable"], na.rm = TRUE),
    median_dynamic = median(test_df$log2_alpha_Endo[test_df$binding_class == cls], na.rm = TRUE),
    median_difference_dynamic_minus_stable = median_dynamic - median_stable,
    p_value = wt$p.value
  )
}) %>%
  mutate(
    fdr = p.adjust(p_value, method = "BH")
  )

wilcox_results

# ------------------------------------------------------------
# 2. Fisher tests:
#    proportion of ACTIVE sequences in each class vs Stable
# ------------------------------------------------------------

stable_all <- class_defs$Stable(MPRA_fullShuffle) %>%
  mutate(binding_class = "Stable")

fisher_results <- map_dfr(dynamic_classes, function(cls) {
  
  dynamic_all <- class_defs[[cls]](MPRA_fullShuffle) %>%
    mutate(binding_class = cls)
  
  test_df <- bind_rows(stable_all, dynamic_all) %>%
    mutate(
      is_active = enhancer == "Active"
    )
  
  active_stable <- sum(test_df$binding_class == "Stable" & test_df$is_active, na.rm = TRUE)
  inactive_stable <- sum(test_df$binding_class == "Stable" & !test_df$is_active, na.rm = TRUE)
  
  active_dynamic <- sum(test_df$binding_class == cls & test_df$is_active, na.rm = TRUE)
  inactive_dynamic <- sum(test_df$binding_class == cls & !test_df$is_active, na.rm = TRUE)
  
  fisher_mat <- matrix(
    c(
      active_dynamic, inactive_dynamic,
      active_stable, inactive_stable
    ),
    nrow = 2,
    byrow = TRUE
  )
  
  rownames(fisher_mat) <- c(cls, "Stable")
  colnames(fisher_mat) <- c("Active", "Inactive")
  
  ft <- fisher.test(fisher_mat)
  
  tibble(
    comparison = paste(cls, "vs Stable"),
    dynamic_class = cls,
    active_dynamic = active_dynamic,
    inactive_dynamic = inactive_dynamic,
    active_stable = active_stable,
    inactive_stable = inactive_stable,
    prop_active_dynamic = active_dynamic / (active_dynamic + inactive_dynamic),
    prop_active_stable = active_stable / (active_stable + inactive_stable),
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value
  )
}) %>%
  mutate(
    fdr = p.adjust(p_value, method = "BH")
  )

fisher_results

# ------------------------------------------------------------
# Optional: combine both result tables for export
# ------------------------------------------------------------

wilcox_results_export <- wilcox_results %>%
  mutate(test = "Wilcoxon: log2_alpha_Endo among Active enhancers")

fisher_results_export <- fisher_results %>%
  mutate(test = "Fisher: proportion Active vs Inactive")

write.csv(
  wilcox_results_export,
  "MPRA_Wilcoxon_dynamic_vs_stable_active_enhancers.csv",
  row.names = FALSE
)

write.csv(
  fisher_results_export,
  "MPRA_Fisher_dynamic_vs_stable_active_proportions.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Optional: plot Fisher proportions
# ------------------------------------------------------------

fisher_plot_df <- fisher_results %>%
  select(dynamic_class, prop_active_dynamic, prop_active_stable) %>%
  tidyr::pivot_longer(
    cols = c(prop_active_dynamic, prop_active_stable),
    names_to = "group",
    values_to = "prop_active"
  ) %>%
  mutate(
    group = recode(
      group,
      prop_active_dynamic = "Dynamic",
      prop_active_stable = "Stable"
    ),
    dynamic_class = factor(dynamic_class, levels = dynamic_classes)
  )

ggplot(fisher_plot_df, aes(x = dynamic_class, y = prop_active, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), colour = "black", linewidth = 0.3) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = NULL,
    y = "Proportion Active",
    title = "Proportion of active MPRA sequences by binding class",
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

plot_cols <- c(
  "Stable"   = "#808080",
  "Age_up"   = "#87CEEB",
  "Age_down" = "#87CEEB",
  "Ex_up"    = "#87CEEB",
  "Ex_down"  = "#87CEEB"
)

MPRA_enhancers <- MPRA_fullShuffle %>%
  filter(enhancer == "Active")

plot_df <- bind_rows(
  MPRA_enhancers %>%
    filter(in_YS == TRUE, in_MAS == TRUE, in_MAE == TRUE) %>%
    mutate(binding_class = "Stable"),
  
  MPRA_enhancers %>%
    filter(in_YS == FALSE, in_MAS == TRUE) %>%
    mutate(binding_class = "Age_up"),
  
  MPRA_enhancers %>%
    filter(in_YS == TRUE, in_MAS == FALSE) %>%
    mutate(binding_class = "Age_down"),
  
  MPRA_enhancers %>%
    filter(in_MAS == FALSE, in_MAE == TRUE) %>%
    mutate(binding_class = "Ex_up"),
  
  MPRA_enhancers %>%
    filter(in_MAS == TRUE, in_MAE == FALSE) %>%
    mutate(binding_class = "Ex_down")
) %>%
  mutate(
    binding_class = factor(
      binding_class,
      levels = c("Stable", "Age_up", "Age_down", "Ex_up", "Ex_down")
    )
  )


p <- ggplot(plot_df, aes(x = binding_class, y = log2_alpha_Endo, fill = binding_class)) +
  geom_violin(trim = FALSE, scale = "width", colour = "black", linewidth = 0.3, width = 0.5) +
  geom_boxplot(width = 0.15, outlier.shape = NA, colour = "black", fill = "white", linewidth = 0.3) +
  geom_jitter(width = 0.12, size = 0.4, alpha = 0.1) +
  scale_fill_manual(values = plot_cols, drop = FALSE) +
  labs(
    x = NULL,
    y = "Endogenous sequence activity\nlog2 alpha",
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
  )

p
setwd(publication_path("PAPER_OUTPUT_ROOT", "Paper_Tables/SupTabble/"))
ggsave("MPRA_violin.pdf", width = 2.5, height = 2.5)

MPRA_sup = read.csv(publication_path("PAPER_OUTPUT_ROOT", "Paper_Tables/SupTabble/Supplementary_Table_MPRA_Results_with_TSS_distance.csv"))

head(MPRA_sup)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# Basic GC-content violin plot by MPRA category
# Change x = category to x = Group or x = is_stable if that is the comparison you want.

plot_df <- MPRA_sup %>%
  filter(!is.na(GC_content)) %>%
  mutate(
    category = factor(Group),
    is_stable = factor(is_stable, levels = c(1, 0), labels = c("Stable", "Dynamic")),
    Group = factor(Group,
      levels = c("Stable", "Age Up", "Age Down", "Ex Up", "Ex Down")
    )
  )


  

gc_cols <- c(
  "Stable" = "#808080",
  "Dynamic" = "#87CEEB"
)

p_gc <- ggplot(plot_df, aes(x = Group, y = GC_content)) +
  geom_violin(
    aes(fill = is_stable),
    trim = FALSE,
    scale = "width",
    colour = "black",
    linewidth = 0.3
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    colour = "black",
    fill = "white",
    linewidth = 0.3
  ) +
  geom_jitter(
    width = 0.12,
    size = 0.4,
    alpha = 0.1
  ) +
  scale_fill_manual(values = gc_cols, drop = FALSE) +
  scale_colour_manual(values = gc_cols, drop = FALSE) +
  labs(
    x = NULL,
    y = "GC content",
    title = "GC content of MPRA sequences",
    fill = NULL,
    colour = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )
library(Seurat)
p_gc +NoLegend()
setwd(publication_path("PAPER_OUTPUT_ROOT", "Paper_Tables/SupTabble/"))
ggsave("MPRA_violin_GC.pdf", width = 2.5, height = 2.5)
