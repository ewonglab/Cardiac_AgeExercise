# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Calculate MPRA reporter activity from count matrices and sample metadata.
# Steps: 1) Load DNA/RNA counts; 2) estimate activity and contrasts; 3) export analysis tables for QC and figures.
# =======================================
# Source reformatted from MPRA/MPRA_analysis.Rmd.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

#!/usr/bin/env Rscript
# Calculate MPRA activity scores (RNA/DNA ratios)
# Usage: Rscript calculate_activity.R <count_matrix.tsv> <samples.csv> <output_dir>

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(dplyr)  # Added: case_when() is used before the later library(dplyr) call
})

# Parse command line arguments
COUNT_MATRIX="${WORK_DIR}/results/counts/count_matrix.tsv"


count_file <- publication_path("MPRA_ROOT", "results/counts/count_matrix.tsv")
sample_file <- publication_path("MPRA_ROOT", "config/samples.csv")
output_dir <- publication_path("MPRA_ROOT", "results/analysis")

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

library(dplyr)

# 1. Filter samples using dplyr::select where appropriate
dna_cols <- samples %>% filter(type == "DNA") %>% pull(sample_id)
rna_cols <- samples %>% filter(type == "RNA") %>% pull(sample_id)

dna_mat <- counts %>% 
  as.data.frame() %>% 
  dplyr::select(all_of(dna_cols)) %>% 
  as.matrix()

rna_mat <- counts %>% 
  as.data.frame() %>% 
  dplyr::select(all_of(rna_cols)) %>% 
  as.matrix()

# Extract seqeunces 
All_TOBIAS_BINDING <- readRDS(publication_path("HEART_SNRNA_ROOT", "analysis/2024-11-30.Tobias_Ap1/Tobias_df_Full_table.rds"))
library(dplyr)
library(rGREAT)
Fos = All_TOBIAS_BINDING %>% filter(motif == "Fos" & cell == "Fibroblast") 

print(class(rna_mat))
ids <- rownames(rna_mat)
print(length(ids))
print(head(ids))

# 1. Extract IDs
ids <- rownames(rna_mat)

# 2. Split by the pipe character (fixed = TRUE is faster and safer here)
# This creates a list of vectors
split_list <- strsplit(ids, split = "|", fixed = TRUE)

# 3. Build the data frame
# We use sapply to grab the first and second elements of each list item
element_metadata <- data.frame(
  ElementID = ids,
  Category  = sapply(split_list, `[`, 1),
  Coords    = sapply(split_list, `[`, 2),
  stringsAsFactors = FALSE
)

# 4. Verify the result
head(element_metadata)
#head(element_metadata)
element_metadata$Coords = stringr::str_remove(element_metadata$Coords, "_2$|_3$")


#Bring in DAP information 
# Read in DAP 
DAP = read.csv(publication_path("HEART_SNRNA_ROOT", "analysis/2023-10-11.AllMacs/All_DAP.csv"))
#head(DAP)
#colnames(DAP)
#colnames(DAP) = c("ID","peak_coord","pvalue","log2FoldChange","pct.1","pct.2","p_val_adj","cell","condition")
#boxplot(DAP$log2FoldChange)
DAP$cell = ifelse(DAP$cell == "Macrophage","AllMacrophage", DAP$cell)
DAP$state = paste(DAP$peak_coord, DAP$cell)

Age = DAP[ DAP$cond == 'Aged_Sed - Young_Sed',]
Ex = DAP[ DAP$cond == 'Aged_Ex - Aged_Sed',]

#table(Age$state %in% Ex$state)

Age =Age[Age$state %in% Ex$state,]
Age =Age[order(Age$state),]
#Age = na.omit(Age)
#table(Age$cell)
Ex =Ex[Ex$state %in% Age$state,]
Ex =Ex[order(Ex$state),]
#Ex = na.omit(Ex)
#table(Ex$cell)
colnames(Ex) = paste0(colnames(Ex),"_Ex")
colnames(Age) = paste0(colnames(Age),"_Age")

table(Ex$state == Age$state)
Age_Ex_ATAC = cbind(Ex, Age)
sigLvl = 0.1

FC = 0
Joint = Age_Ex_ATAC
FC = 0
P_thre = .1
# Directions
Age_Ex_ATAC$Direction = "NSig"
Age_Ex_ATAC$Direction = ifelse( Age_Ex_ATAC$pvalue_Age < P_thre & abs(Age_Ex_ATAC$log2FoldChange_Age) > FC &  (Age_Ex_ATAC$pvalue_Ex > P_thre | abs(Age_Ex_ATAC$log2FoldChange_Ex) < FC) ,"Age DEG Age Down" ,Age_Ex_ATAC$Direction) 
Age_Ex_ATAC$Direction = ifelse( Age_Ex_ATAC$pvalue_Ex < P_thre & abs(Age_Ex_ATAC$log2FoldChange_Ex) > FC &  (Age_Ex_ATAC$pvalue_Age > P_thre | abs(Age_Ex_ATAC$log2FoldChange_Age) < FC) ,"Ex DEG Ex Down" ,Age_Ex_ATAC$Direction) 
Age_Ex_ATAC$Direction = ifelse( Age_Ex_ATAC$pvalue_Age < P_thre & abs(Age_Ex_ATAC$log2FoldChange_Age) > FC & (Age_Ex_ATAC$log2FoldChange_Age) > 0 &  (Age_Ex_ATAC$pvalue_Ex > P_thre | abs(Age_Ex_ATAC$log2FoldChange_Ex) < FC) ,"Age DEG Age Up" ,Age_Ex_ATAC$Direction) 
Age_Ex_ATAC$Direction = ifelse( Age_Ex_ATAC$pvalue_Ex < P_thre & abs(Age_Ex_ATAC$log2FoldChange_Ex) > FC & (Age_Ex_ATAC$log2FoldChange_Ex) > 0 &  (Age_Ex_ATAC$pvalue_Age > P_thre | abs(Age_Ex_ATAC$log2FoldChange_Age) < FC) ,"Ex DEG Ex Up" ,Age_Ex_ATAC$Direction) 


Age_Ex_ATAC$Direction = ifelse(Age_Ex_ATAC$pvalue_Age < P_thre & Age_Ex_ATAC$pvalue_Ex < P_thre & Age_Ex_ATAC$log2FoldChange_Age > FC & Age_Ex_ATAC$log2FoldChange_Ex > FC ,"Exercise exacerbates Age Up",Age_Ex_ATAC$Direction) 
Age_Ex_ATAC$Direction = ifelse(Age_Ex_ATAC$pvalue_Age < P_thre & Age_Ex_ATAC$pvalue_Ex < P_thre & Age_Ex_ATAC$log2FoldChange_Age < -FC & Age_Ex_ATAC$log2FoldChange_Ex < -FC ,"Exercise exacerbates Age Down",Age_Ex_ATAC$Direction)
Age_Ex_ATAC$Direction = ifelse(Age_Ex_ATAC$pvalue_Age < P_thre & Age_Ex_ATAC$pvalue_Ex < P_thre & Age_Ex_ATAC$log2FoldChange_Age < -FC & Age_Ex_ATAC$log2FoldChange_Ex > FC ,"Exercise reverses Age Down",Age_Ex_ATAC$Direction)
Age_Ex_ATAC$Direction = ifelse(Age_Ex_ATAC$pvalue_Age < P_thre & Age_Ex_ATAC$pvalue_Ex < P_thre & Age_Ex_ATAC$log2FoldChange_Age > FC & Age_Ex_ATAC$log2FoldChange_Ex < -FC ,"Exercise reverses Age Up",Age_Ex_ATAC$Direction)


Age_Ex_ATAC$AgeDirection = "NSig"
Age_Ex_ATAC$AgeDirection = ifelse( Age_Ex_ATAC$pvalue_Age < P_thre & abs(Age_Ex_ATAC$log2FoldChange_Age) > FC &  (Age_Ex_ATAC$pvalue_Ex > P_thre | abs(Age_Ex_ATAC$log2FoldChange_Ex) < FC) ,"Age DEG" ,Age_Ex_ATAC$AgeDirection) 
Age_Ex_ATAC$ExDirection = "NSig"
Age_Ex_ATAC$ExDirection = ifelse( Age_Ex_ATAC$pvalue_Ex < P_thre & abs(Age_Ex_ATAC$log2FoldChange_Ex) > FC &  (Age_Ex_ATAC$pvalue_Age > P_thre | abs(Age_Ex_ATAC$log2FoldChange_Age) < FC) ,"Ex DEG" ,Age_Ex_ATAC$ExDirection) 
df =data.frame((table(Age_Ex_ATAC$Direction , Age_Ex_ATAC$cell_Age)))

library(tidyverse)


adjust_start_site_Tob <- function(string) {
  # Use regular expression to extract chromosome, start, and end
  matches <- regmatches(string, regexec("([^-]+)-(\\d+)-(\\d+)", string))
  
  # Extract components
  chr <- matches[[1]][2]
  start <- as.numeric(matches[[1]][3])
  end <- matches[[1]][4]
  
  # Add 1 to the start site
  new_start <- start + 1
  
  # Rebuild the string
  return(paste0(chr, "-", new_start, "-", end))
}

regions_tob  = lapply( Fos$identifier, adjust_start_site_Tob)
regions_tob = unlist(regions_tob)
Fos$matched_ID = regions_tob
table(Fos$matched_ID %in% Age_Ex_ATAC$peak_coord_Age)

Age_Ex_ATAC = Age_Ex_ATAC %>% subset( cell_Age =="Fibroblast")
table(Fos$matched_ID %in% Age_Ex_ATAC$peak_coord_Age)
table(Fos$matched_ID %in% Age_Ex_ATAC$peak_coord_Ex)

Joint = merge(Fos, Age_Ex_ATAC, by.x = 'matched_ID', by.y ='peak_coord_Age', all = FALSE)

DT::datatable(element_metadata, filter = 'top')

extend_to_170 <- function(gr) {
  stopifnot(inherits(gr, "GRanges"))
  target_w <- pmax(width(gr), 170L)                    # don't shrink wider peaks
  gr_ext <- resize(gr, width = target_w, fix = "center", ignore.strand = TRUE)
  if (all(!is.na(seqlengths(gr_ext)))) gr_ext <- trim(gr_ext)
  gr_ext
}

# Example:
motif_gr = (GRanges(seqnames = Joint$motif_chr, ranges = IRanges(Joint$motif_start, Joint$motif_end)))
 real_gr170 <- extend_to_170(motif_gr)
Joint$extend = paste0(seqnames(real_gr170), ":",start(real_gr170),"-",end(real_gr170))

table(element_metadata$Coords %in% Joint$extend)
element_metadata$Coords[!element_metadata$Coords %in% Joint$extend]
element_metadata = merge(element_metadata, Joint, by.x = 'Coords', by.y = 'extend', all.x = TRUE, all.y = FALSE)

library(tidyr)
library(dplyr)
library(MPRAnalyze)

# Create the master annotation table for all elements in the library
element_metadata
# For the MpraObject, we also need the sample-level metadata (your 9 RNA replicates)
# If all 9 RNA replicates represent the same biological state (e.g. one cell type),
# our rnaAnnot remains simple.
rna_annot <- samples_rna

dna_annot <- samples_dna

table(element_metadata$Category)
table(element_metadata$category)

library(dplyr)
library(tidyr)
library(MPRAnalyze)

library(dplyr)
library(tidyr)
library(MPRAnalyze)
# Inside run_mpra_comparison:


# --- Updated Helper Function with Unique Column Names Fix ---
library(dplyr)
library(tidyr)
library(MPRAnalyze)

run_mpra_comparison <- function(rna_in, dna_in, meta, cond_a, cond_b) {
  
  message(paste0("Running Pooled Comparison (Human + Mouse): ", cond_a, " vs ", cond_b))
  
  # 1. Standard Subset and Alignment (same logic as before)
  ids_a <- meta %>% filter(Category == cond_a) %>% pull(ElementID)
  ids_b <- meta %>% filter(Category == cond_b) %>% pull(ElementID)
  
  mat_a_rna <- rna_in[ids_a, ]
  rownames(mat_a_rna) <- meta$Coords[match(ids_a, meta$ElementID)]
  mat_b_rna <- rna_in[ids_b, ]
  rownames(mat_b_rna) <- meta$Coords[match(ids_b, meta$ElementID)]
  
  mat_a_dna <- dna_in[ids_a, ]
  rownames(mat_a_dna) <- meta$Coords[match(ids_a, meta$ElementID)]
  mat_b_dna <- dna_in[ids_b, ]
  rownames(mat_b_dna) <- meta$Coords[match(ids_b, meta$ElementID)]
  
  common <- intersect(rownames(mat_a_rna), rownames(mat_b_rna))
  message(paste("Common Coords found:", length(common)))
  
    # 2. Add an extra check for DNA "health"
    # Calculate the average DNA across replicates for both conditions
    avg_dna_a <- rowMeans(mat_a_dna[common, ])
    avg_dna_b <- rowMeans(mat_b_dna[common, ])

    # Keep only those where both states have an average of at least 10 reads
    healthy_pairs <- common[avg_dna_a >= 10 & avg_dna_b >= 10]

    message(paste("Coordinates passing DNA health check:", length(healthy_pairs)))

    # Use healthy_pairs instead of common for the rest of the function
    rna_paired <- cbind(mat_a_rna[healthy_pairs, ], mat_b_rna[healthy_pairs, ])
      colnames(rna_paired) <- c(paste0(colnames(mat_a_rna), "_", cond_a), 
                                paste0(colnames(mat_b_rna), "_", cond_b))

      dna_paired <- cbind(mat_a_dna[healthy_pairs, ], mat_b_dna[healthy_pairs, ])
      colnames(dna_paired) <- c(paste0(colnames(mat_a_dna), "_", cond_a), 
                                paste0(colnames(mat_b_dna), "_", cond_b))

  # 2. CREATE COMPLEX ANNOTATIONS
  # We repeat the metadata for both 'Original' and 'Shuffle' halves of the matrix
  n_reps <- ncol(rna_in) # Should be 9
  
  # Species: 1-3 Human, 4-9 Mouse
  species_vec <- factor(c(rep("Human", 3), rep("Mouse", 6)))
  # Bio Reps: 1-3 Rep1, 4-6 RepA, 7-9 RepB
  replicate_vec <- factor(c(rep("H_Rep1", 3), rep("M_RepA", 3), rep("M_RepB", 3)))
  
  rna_annot <- data.frame(
    condition = factor(rep(c(cond_a, cond_b), each = n_reps), levels = c(cond_b, cond_a)),
    species = rep(species_vec, 2),
    batch = rep(replicate_vec, 2),
    row.names = colnames(rna_paired)
  )
  
  # For DNA, we replicate the batch info as well to match the structure
  dna_annot <- data.frame(
    condition = factor(rep(c(cond_a, cond_b), each = ncol(dna_in)), levels = c(cond_b, cond_a)),
    species = rep(species_vec, 2),
       batch = rep(replicate_vec, 2),# Assuming DNA mapping follows RNA
    row.names = colnames(dna_paired)
  )
  
  # 3. Initialize MpraObject
  obj <- MpraObject(dnaCounts = dna_paired, rnaCounts = rna_paired, 
                    dnaAnnot = dna_annot, rnaAnnot = rna_annot)
  
  obj <- estimateDepthFactors(obj, lib.factor = "condition")
  
  # 4. Fit Model with Species Adjustment
  # We include 'species' in the rnaDesign to account for the cell-line batch effect.
  # We keep 'condition' in both to test for sequence-specific activity.
  obj <- analyzeComparative(obj, 
                            dnaDesign = ~ condition + batch, 
                            rnaDesign = ~ condition + batch,
                            reducedDesign = ~ batch) # Testing 'condition' given 'species'
  
  # 5. Extract Results
  res <- tryCatch({
    testLrt(obj)
  }, error = function(e) {
    message("testLrt failed, falling back to testCondition (Wald Test)...")
    obj <- analyzeComparative(obj, 
                              dnaDesign = ~ condition + species, 
                              rnaDesign = ~ condition + species,
                              fit_SE = TRUE)
    # Testing the 'condition' effect specifically
    return(testCondition(obj, fixedEffect = "condition", term = cond_a))
  })

  res$Coords <- rownames(res)
  res$Comparison <- paste0(cond_a, "_vs_", cond_b)
  res$fdr <- p.adjust(res$pval, method = "BH")

  # Keep a simple effect-size column for plotting (same scale as logFC)
  res$estimate <- res$logFC

  # Alpha (activity) estimates for each element under each condition (ORIG vs SHUF_*)
  alpha <- getAlpha(obj, by.factor = "condition")
  log2_alpha <- log2(pmax(alpha, .Machine$double.xmin))
  res$log2_alpha_a <- log2_alpha[res$Coords, as.character(cond_a)]
  res$log2_alpha_b <- log2_alpha[res$Coords, as.character(cond_b)]

  return(res)
}
message("Run")
# --- Run the Comparisons as before ---
res_full_shuf <- run_mpra_comparison(rna_mat, dna_mat, element_metadata, 
                                     cond_a = "ORIG", cond_b = "SHUF_FULL")

res_motif_shuf <- run_mpra_comparison(rna_mat, dna_mat, element_metadata, 
                                      cond_a = "ORIG", cond_b = "MOTIF_SHUF")



# --- 3. Combine and Classify ---

all_results <- bind_rows(res_full_shuf, res_motif_shuf)

# Define "Significant Impact" at p < 0.1
all_results <- all_results %>%
  mutate(is_significant = ifelse(fdr < 0.1, "Significant", "Non-Significant"))

# Quick summary
table(all_results$Comparison, all_results$is_significant)

# Save results
saveRDS(all_results, "MPRA_Comparative_Results_All_Conditions.rds")



all_results <- bind_rows(res_full_shuf, res_motif_shuf)

# Define "Significant Impact" at p < 0.1
all_results <- all_results %>%
  mutate(is_significant = ifelse(fdr < 0.05, "Significant", "Non-Significant"))

# Quick summary
table(all_results$Comparison, all_results$is_significant)

# Save results
saveRDS(all_results, "MPRA_Comparative_Results_All_Conditions.rds")
