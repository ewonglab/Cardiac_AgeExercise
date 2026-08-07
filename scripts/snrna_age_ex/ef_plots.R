# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Create echocardiography and strain outcome plots.
# Steps: 1) Load animal-level cardiac function data; 2) fit reported comparisons; 3) export figure-ready summary plots.
# =======================================
# Source reformatted from snRNA_AgeEx/EF_plots.ipynb.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

#EF 
library(OmnipathR)
library(Seurat)
library(Seurat)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(dplyr)
library(Seurat)
library(patchwork)
library(sctransform)
library(ggplot2)
library(pheatmap)
library(Seurat)
library(ggplot2)
library(Matrix)
library(RColorBrewer)
library(dplyr)
library(scales)
library(data.table)
library(stats)
library("Nebulosa")
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(magrittr)
library(dplyr)
library(ggpmisc)



library(ggpubr)

setwd(publication_path("HEART_SNRNA_ROOT", "data/2024-05-22.EF"))

EF_data = readRDS("harmony_annotated.rds")

colnames(EF_data@meta.data)
table(EF_data$cell.type)

library(openxlsx)
EF_value = read.xlsx(publication_path("HEART_SNRNA_ROOT", "data/2024-04-20.Stravos/Ef/supplementarytable1.xlsx"))

library(readxl)
library(dplyr)

# Load the Excel file
file_path <- 'NIHMS1915768-supplement-supplementary_tables.xlsx'

# Get the names of all sheets in the Excel file
sheet_names <- excel_sheets(file_path)

sheet_names = grep(sheet_names, pattern = "table 5|table 6|table 7", value = TRUE)
sheet_names

# Initialize a list to store the dataframes
sheets_data <- list()

# Iterate through each sheet in the Excel file
for (sheet_name in sheet_names) {
  # Read the sheet into a DataFrame
  df <- read_excel(file_path, sheet = sheet_name)
    df = data.frame(df)
  #print(head(df))
  # Get the first line (first row) as the name
       colnames(df)= df[1,]
    df = df[-1,]
    df$baseMean = as.double(df$baseMean)
    df$log2FoldChange = as.double(df$log2FoldChange)
    df$pvalue = as.double(df$pvalue)
    df$padj = as.double(df$padj)
  #print(sheet_name)
  # Store the DataFrame in the list with the name from the first line
  sheets_data[[sheet_name]] <- df
}



library(stringr)

names(sheets_data) = str_replace(names(sheets_data),"supplementary table 5","RRpostVsDonor")
names(sheets_data) = str_replace(names(sheets_data),"supplementary table 6","preLVADHFvsDonor")
names(sheets_data) = str_replace(names(sheets_data),"supplementary table 7","RRpostvspreLVADHF")

names(sheets_data)
head(sheets_data[['RRpostVsDonor adipocyte']])

Adipocyte <- sheets_data[['RRpostVsDonor adipocyte']]
Adipocyte$cell <- "Adipocyte"
Adipocyte$sigpvalue <- ifelse(Adipocyte$padj < 0.05, "p < 0.05","p > 0.05")
Adipocyte$sig <- ifelse(Adipocyte$padj < 0.05 & abs(Adipocyte$log2FoldChange) > 0.58, "Significant","Not Significant")

Cardiomyocyte <- sheets_data[['RRpostVsDonor cardiomyo']]
Cardiomyocyte$cell <- "Cardiomyocyte"
Cardiomyocyte$sigpvalue <- ifelse(Cardiomyocyte$padj < 0.05, "p < 0.05","p > 0.05")
Cardiomyocyte$sig <- ifelse(Cardiomyocyte$padj < 0.05 & abs(Cardiomyocyte$log2FoldChange) > 0.58, "Significant","Not Significant")

Endocardium <- sheets_data[['RRpostVsDonor endocardi']]
Endocardium$cell <- "Endocardium"
Endocardium$sigpvalue <- ifelse(Endocardium$padj < 0.05, "p < 0.05","p > 0.05")
Endocardium$sig <- ifelse(Endocardium$padj < 0.05 & abs(Endocardium$log2FoldChange) > 0.58, "Significant","Not Significant")

Endothelium <-sheets_data[['RRpostVsDonor endotheli']]
Endothelium$cell <- "Endothelium"
Endothelium$sigpvalue <- ifelse(Endothelium$padj < 0.05, "p < 0.05","p > 0.05")
Endothelium$sig <- ifelse(Endothelium$padj < 0.05 & abs(Endothelium$log2FoldChange) > 0.58, "Significant","Not Significant")

Epicardium <- sheets_data[['RRpostVsDonor epicardiu']]
Epicardium$cell <- "Epicardium"
Epicardium$sigpvalue <- ifelse(Epicardium$padj < 0.05, "p < 0.05","p > 0.05")
Epicardium$sig <- ifelse(Epicardium$padj < 0.05 & abs(Epicardium$log2FoldChange) > 0.58, "Significant","Not Significant")

Fibroblast <- sheets_data[['RRpostVsDonor fibroblas']]
Fibroblast$cell <- "Fibroblast"
Fibroblast$sigpvalue <- ifelse(Fibroblast$padj < 0.05, "p < 0.05","p > 0.05")
Fibroblast$sig <- ifelse(Fibroblast$padj < 0.05 & abs(Fibroblast$log2FoldChange) > 0.58, "Significant","Not Significant")

Glia <- sheets_data[['RRpostVsDonor glia']]
Glia$cell <- "Glia"
Glia$sigpvalue <- ifelse(Glia$padj < 0.05, "p < 0.05","p > 0.05")
Glia$sig <- ifelse(Glia$padj < 0.05 & abs(Glia$log2FoldChange) > 0.58, "Significant","Not Significant")

Lymphatic <- sheets_data[['RRpostVsDonor lymphatic']]
Lymphatic$cell <- "Lymphatic"
Lymphatic$sigpvalue <- ifelse(Lymphatic$padj < 0.05, "p < 0.05","p > 0.05")
Lymphatic$sig <- ifelse(Lymphatic$padj < 0.05 & abs(Lymphatic$log2FoldChange) > 0.58, "Significant","Not Significant")

Mast <- sheets_data[['RRpostVsDonor mast']]
Mast$cell <- "Mast"
Mast$sigpvalue <- ifelse(Mast$padj < 0.05, "p < 0.05","p > 0.05")
Mast$sig <- ifelse(Mast$padj < 0.05 & abs(Mast$log2FoldChange) > 0.58, "Significant","Not Significant")

Myeloid <- sheets_data[['RRpostVsDonor myeloid']]
Myeloid$cell <- "Myeloid"
Myeloid$sigpvalue <- ifelse(Myeloid$padj < 0.05, "p < 0.05","p > 0.05")
Myeloid$sig <- ifelse(Myeloid$padj < 0.05 & abs(Myeloid$log2FoldChange) > 0.58, "Significant","Not Significant")

Pericyte <- sheets_data[['RRpostVsDonor pericyte']]
Pericyte$cell <- "Pericyte"
Pericyte$sigpvalue <- ifelse(Pericyte$padj < 0.05, "p < 0.05","p > 0.05")
Pericyte$sig <- ifelse(Pericyte$padj < 0.05 & abs(Pericyte$log2FoldChange) > 0.58, "Significant","Not Significant")

SMC <- sheets_data[['RRpostVsDonor smc']]
SMC$cell <- "SMC"
SMC$sigpvalue <- ifelse(SMC$padj < 0.05, "p < 0.05","p > 0.05")
SMC$sig <- ifelse(SMC$padj < 0.05 & abs(SMC$log2FoldChange) > 0.58, "Significant","Not Significant")

TNKCells <- sheets_data[['RRpostVsDonor tnk cells']]
TNKCells$cell <- "TNKCells"
TNKCells$sigpvalue <- ifelse(TNKCells$padj < 0.05, "p < 0.05","p > 0.05")
TNKCells$sig <- ifelse(TNKCells$padj < 0.05 & abs(TNKCells$log2FoldChange) > 0.58, "Significant","Not Significant")
data <- data.frame(rbind(Adipocyte, Cardiomyocyte, Endocardium, Endothelium, Epicardium, Fibroblast, Glia, Lymphatic, Mast, Myeloid, Pericyte, SMC, TNKCells))
data$cell <- factor(data$cell, levels = c("Adipocyte","Cardiomyocyte", "Endocardium", "Endothelium", "Epicardium", "Fibroblast", "Glia", "Lymphatic", "Mast", "Myeloid", "Pericyte", "SMC", "TNKCells"))

data$cluster <- ifelse(data$log2FoldChange > 0, "Rpost", "Donor")
data_rpost_donor <- filter(data, data$sig == "Significant")

Adipocyte <- sheets_data[['preLVADHFvsDonor adipocyte']]
Adipocyte$cell <- "Adipocyte"
Adipocyte$sigpvalue <- ifelse(Adipocyte$padj < 0.05, "p < 0.05","p > 0.05")
Adipocyte$sig <- ifelse(Adipocyte$padj < 0.05 & abs(Adipocyte$log2FoldChange) > 0.58, "Significant","Not Significant")

Cardiomyocyte <- sheets_data[['preLVADHFvsDonor cardiomyo']]
Cardiomyocyte$cell <- "Cardiomyocyte"
Cardiomyocyte$sigpvalue <- ifelse(Cardiomyocyte$padj < 0.05, "p < 0.05","p > 0.05")
Cardiomyocyte$sig <- ifelse(Cardiomyocyte$padj < 0.05 & abs(Cardiomyocyte$log2FoldChange) > 0.58, "Significant","Not Significant")

Endocardium <- sheets_data[['preLVADHFvsDonor endocardi']]
Endocardium$cell <- "Endocardium"
Endocardium$sigpvalue <- ifelse(Endocardium$padj < 0.05, "p < 0.05","p > 0.05")
Endocardium$sig <- ifelse(Endocardium$padj < 0.05 & abs(Endocardium$log2FoldChange) > 0.58, "Significant","Not Significant")

Endothelium <-sheets_data[['preLVADHFvsDonor endotheli']]
Endothelium$cell <- "Endothelium"
Endothelium$sigpvalue <- ifelse(Endothelium$padj < 0.05, "p < 0.05","p > 0.05")
Endothelium$sig <- ifelse(Endothelium$padj < 0.05 & abs(Endothelium$log2FoldChange) > 0.58, "Significant","Not Significant")

Epicardium <- sheets_data[['preLVADHFvsDonor epicardiu']]
Epicardium$cell <- "Epicardium"
Epicardium$sigpvalue <- ifelse(Epicardium$padj < 0.05, "p < 0.05","p > 0.05")
Epicardium$sig <- ifelse(Epicardium$padj < 0.05 & abs(Epicardium$log2FoldChange) > 0.58, "Significant","Not Significant")

Fibroblast <- sheets_data[['preLVADHFvsDonor fibroblas']]
Fibroblast$cell <- "Fibroblast"
Fibroblast$sigpvalue <- ifelse(Fibroblast$padj < 0.05, "p < 0.05","p > 0.05")
Fibroblast$sig <- ifelse(Fibroblast$padj < 0.05 & abs(Fibroblast$log2FoldChange) > 0.58, "Significant","Not Significant")

Glia <- sheets_data[['preLVADHFvsDonor glia']]
Glia$cell <- "Glia"
Glia$sigpvalue <- ifelse(Glia$padj < 0.05, "p < 0.05","p > 0.05")
Glia$sig <- ifelse(Glia$padj < 0.05 & abs(Glia$log2FoldChange) > 0.58, "Significant","Not Significant")

Lymphatic <- sheets_data[['preLVADHFvsDonor lymphatic']]
Lymphatic$cell <- "Lymphatic"
Lymphatic$sigpvalue <- ifelse(Lymphatic$padj < 0.05, "p < 0.05","p > 0.05")
Lymphatic$sig <- ifelse(Lymphatic$padj < 0.05 & abs(Lymphatic$log2FoldChange) > 0.58, "Significant","Not Significant")

Mast <- sheets_data[['preLVADHFvsDonor mast']]
Mast$cell <- "Mast"
Mast$sigpvalue <- ifelse(Mast$padj < 0.05, "p < 0.05","p > 0.05")
Mast$sig <- ifelse(Mast$padj < 0.05 & abs(Mast$log2FoldChange) > 0.58, "Significant","Not Significant")

Myeloid <- sheets_data[['preLVADHFvsDonor myeloid']]
Myeloid$cell <- "Myeloid"
Myeloid$sigpvalue <- ifelse(Myeloid$padj < 0.05, "p < 0.05","p > 0.05")
Myeloid$sig <- ifelse(Myeloid$padj < 0.05 & abs(Myeloid$log2FoldChange) > 0.58, "Significant","Not Significant")

Pericyte <- sheets_data[['preLVADHFvsDonor pericyte']]
Pericyte$cell <- "Pericyte"
Pericyte$sigpvalue <- ifelse(Pericyte$padj < 0.05, "p < 0.05","p > 0.05")
Pericyte$sig <- ifelse(Pericyte$padj < 0.05 & abs(Pericyte$log2FoldChange) > 0.58, "Significant","Not Significant")

SMC <- sheets_data[['preLVADHFvsDonor smc']]
SMC$cell <- "SMC"
SMC$sigpvalue <- ifelse(SMC$padj < 0.05, "p < 0.05","p > 0.05")
SMC$sig <- ifelse(SMC$padj < 0.05 & abs(SMC$log2FoldChange) > 0.58, "Significant","Not Significant")

TNKCells <- sheets_data[['preLVADHFvsDonor tnk cells']]
TNKCells$cell <- "TNKCells"
TNKCells$sigpvalue <- ifelse(TNKCells$padj < 0.05, "p < 0.05","p > 0.05")
TNKCells$sig <- ifelse(TNKCells$padj < 0.05 & abs(TNKCells$log2FoldChange) > 0.58, "Significant","Not Significant")
data <- data.frame(rbind(Adipocyte, Cardiomyocyte, Endocardium, Endothelium, Epicardium, Fibroblast, Glia, Lymphatic, Mast, Myeloid, Pericyte, SMC, TNKCells))
data$cell <- factor(data$cell, levels = c("Adipocyte","Cardiomyocyte", "Endocardium", "Endothelium", "Epicardium", "Fibroblast", "Glia", "Lymphatic", "Mast", "Myeloid", "Pericyte", "SMC", "TNKCells"))

data$cluster <- ifelse(data$log2FoldChange > 0, "HF", "Donor")
data_hf_donor <- filter(data, data$sig == "Significant")

Adipocyte <- sheets_data[['RRpostvspreLVADHF adipocyte']]
Adipocyte$cell <- "Adipocyte"
Adipocyte$sigpvalue <- ifelse(Adipocyte$padj < 0.05, "p < 0.05","p > 0.05")
Adipocyte$sig <- ifelse(Adipocyte$padj < 0.05 & abs(Adipocyte$log2FoldChange) > 0.58, "Significant","Not Significant")

Cardiomyocyte <- sheets_data[['RRpostvspreLVADHF cardiomyo']]
Cardiomyocyte$cell <- "Cardiomyocyte"
Cardiomyocyte$sigpvalue <- ifelse(Cardiomyocyte$padj < 0.05, "p < 0.05","p > 0.05")
Cardiomyocyte$sig <- ifelse(Cardiomyocyte$padj < 0.05 & abs(Cardiomyocyte$log2FoldChange) > 0.58, "Significant","Not Significant")

Endocardium <- sheets_data[['RRpostvspreLVADHF endocardi']]
Endocardium$cell <- "Endocardium"
Endocardium$sigpvalue <- ifelse(Endocardium$padj < 0.05, "p < 0.05","p > 0.05")
Endocardium$sig <- ifelse(Endocardium$padj < 0.05 & abs(Endocardium$log2FoldChange) > 0.58, "Significant","Not Significant")

Endothelium <-sheets_data[['RRpostvspreLVADHF endotheli']]
Endothelium$cell <- "Endothelium"
Endothelium$sigpvalue <- ifelse(Endothelium$padj < 0.05, "p < 0.05","p > 0.05")
Endothelium$sig <- ifelse(Endothelium$padj < 0.05 & abs(Endothelium$log2FoldChange) > 0.58, "Significant","Not Significant")

Epicardium <- sheets_data[['RRpostvspreLVADHF epicardiu']]
Epicardium$cell <- "Epicardium"
Epicardium$sigpvalue <- ifelse(Epicardium$padj < 0.05, "p < 0.05","p > 0.05")
Epicardium$sig <- ifelse(Epicardium$padj < 0.05 & abs(Epicardium$log2FoldChange) > 0.58, "Significant","Not Significant")

Fibroblast <- sheets_data[['RRpostvspreLVADHF fibroblas']]
Fibroblast$cell <- "Fibroblast"
Fibroblast$sigpvalue <- ifelse(Fibroblast$padj < 0.05, "p < 0.05","p > 0.05")
Fibroblast$sig <- ifelse(Fibroblast$padj < 0.05 & abs(Fibroblast$log2FoldChange) > 0.58, "Significant","Not Significant")

Glia <- sheets_data[['RRpostvspreLVADHF glia']]
Glia$cell <- "Glia"
Glia$sigpvalue <- ifelse(Glia$padj < 0.05, "p < 0.05","p > 0.05")
Glia$sig <- ifelse(Glia$padj < 0.05 & abs(Glia$log2FoldChange) > 0.58, "Significant","Not Significant")

Lymphatic <- sheets_data[['RRpostvspreLVADHF lymphatic']]
Lymphatic$cell <- "Lymphatic"
Lymphatic$sigpvalue <- ifelse(Lymphatic$padj < 0.05, "p < 0.05","p > 0.05")
Lymphatic$sig <- ifelse(Lymphatic$padj < 0.05 & abs(Lymphatic$log2FoldChange) > 0.58, "Significant","Not Significant")

Mast <- sheets_data[['RRpostvspreLVADHF mast']]
Mast$cell <- "Mast"
Mast$sigpvalue <- ifelse(Mast$padj < 0.05, "p < 0.05","p > 0.05")
Mast$sig <- ifelse(Mast$padj < 0.05 & abs(Mast$log2FoldChange) > 0.58, "Significant","Not Significant")

Myeloid <- sheets_data[['RRpostvspreLVADHF myeloid']]
Myeloid$cell <- "Myeloid"
Myeloid$sigpvalue <- ifelse(Myeloid$padj < 0.05, "p < 0.05","p > 0.05")
Myeloid$sig <- ifelse(Myeloid$padj < 0.05 & abs(Myeloid$log2FoldChange) > 0.58, "Significant","Not Significant")

Pericyte <- sheets_data[['RRpostvspreLVADHF pericyte']]
Pericyte$cell <- "Pericyte"
Pericyte$sigpvalue <- ifelse(Pericyte$padj < 0.05, "p < 0.05","p > 0.05")
Pericyte$sig <- ifelse(Pericyte$padj < 0.05 & abs(Pericyte$log2FoldChange) > 0.58, "Significant","Not Significant")

SMC <- sheets_data[['RRpostvspreLVADHF smc']]
SMC$cell <- "SMC"
SMC$sigpvalue <- ifelse(SMC$padj < 0.05, "p < 0.05","p > 0.05")
SMC$sig <- ifelse(SMC$padj < 0.05 & abs(SMC$log2FoldChange) > 0.58, "Significant","Not Significant")

TNKCells <- sheets_data[['RRpostvspreLVADHF tnk cells']]
TNKCells$cell <- "TNKCells"
TNKCells$sigpvalue <- ifelse(TNKCells$padj < 0.05, "p < 0.05","p > 0.05")
TNKCells$sig <- ifelse(TNKCells$padj < 0.05 & abs(TNKCells$log2FoldChange) > 0.58, "Significant","Not Significant")
data <- data.frame(rbind(Adipocyte, Cardiomyocyte, Endocardium, Endothelium, Epicardium, Fibroblast, Glia, Lymphatic, Mast, Myeloid, Pericyte, SMC, TNKCells))
data$cell <- factor(data$cell, levels = c("Adipocyte","Cardiomyocyte", "Endocardium", "Endothelium", "Epicardium", "Fibroblast", "Glia", "Lymphatic", "Mast", "Myeloid", "Pericyte", "SMC", "TNKCells"))

data$cluster <- ifelse(data$log2FoldChange > 0, "HF", "Rpost")
data_rpost_hf <- filter(data, data$sig == "Significant")

head(data_rpost_hf)

data_rpost_hf_RPOST <- filter(data_rpost_hf, data_rpost_hf$cluster == "Rpost")
data_hf_donor_HF <- filter(data_hf_donor, data_hf_donor$cluster == "HF")
data_rpost_donor_RPOST <- filter(data_rpost_donor, data_rpost_donor$cluster == "Rpost")

fun <- function(x) {
  if (x == "Rpost") {"Rpost"} 
  else if (x == "Donor") {"Donor"}
  else if (x == "NRpost") {"NRpost"}
  else if (x == "Rpre") {"HF"}
  else if (x == "NRpre") {"HF"}
}
EF_data$condition2 <- mapply(fun, EF_data$condition)

metadata <- EF_data@meta.data

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(EF_data$cell.type)
metadata$sample_id <- factor(EF_data$orig.ident)
#metadata$sample_id <- factor(seurat$condition2)
metadata$group_id <- factor(EF_data$condition2)
counts <- EF_data@assays$RNA@counts 

library(SingleCellExperiment)
sce <- SingleCellExperiment(assays = list(counts = counts), 
                           colData = metadata)

# Identify groups for aggregation of counts
groups <- colData(sce)[, c("cluster_id", "sample_id")]

# Named vector of cluster names
kids <- purrr::set_names(levels(sce$cluster_id))

# Total number of clusters
nk <- length(kids)

# Named vector of sample names
sids <- purrr::set_names(levels(sce$sample_id))

# Total number of samples 
ns <- length(sids)

# Generate sample level metadata

## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$sample_id))

## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$sample_id)

## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ], 
                  n_cells, row.names = NULL) %>% 
                select(-"cluster_id")



sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
sce

library(Matrix)
library(Matrix.utils)
groups <- colData(sce)[, c("cluster_id", "sample_id")]

# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 2), 
                 `[`, 1)
# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
# pb <- split.data.frame(pb, 
#                        factor(splitf)) %>%
#         lapply(function(u) 
#                 set_colnames(t(u), 
#                              stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
        lapply(function(u) 
                set_colnames(t(u), 
                             stringr::str_extract(rownames(u), "(?<=_).*")))



options(width = 100)
table(sce$cluster_id)
      table(sce$sample_id)


get_sample_ids <- function(x){
        pb[[x]] %>%
                colnames()
}

de_samples <- map(1:length(kids), get_sample_ids) %>%
        unlist()



# Get cluster IDs for each of the samples

samples_list <- map(1:length(kids), get_sample_ids)

get_cluster_ids <- function(x){
        rep(names(pb)[x], 
            each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
        unlist()



gg_df <- data.frame(cluster_id = de_cluster_ids,
                    sample_id = de_samples)


gg_df <- left_join(gg_df, ei[, c("sample_id", "group_id")]) 


metadata <- gg_df %>%
        dplyr::select(cluster_id, sample_id, group_id) 
        
metadata      


clusters <- unique(metadata$cluster_id)
clusters











# CD163
cluster = 'Myeloid'
cluster_num <- 2

# Subset the metadata to only the desired cells
cluster_metadata <- metadata[which(metadata$cluster_id == cluster), ]

# Assign the rownames of the metadata to be the sample IDs
rownames(cluster_metadata) <- cluster_metadata$sample_id

# Subset the counts to only the desired cells
counts <- pb[[cluster]]

cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
rownames(cluster_metadata) = str_replace_all(rownames(cluster_metadata),"-",".")
# Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
all(rownames(cluster_metadata) == colnames(cluster_counts))   
library(DESeq2)
rld <- rlog(as.matrix(cluster_counts))

#Genes that are Up in RRPost vs PreLvad
#Genes that are Up in RRPost vs donr
#Genes that aren't up in  RRPost vs donr


rpost_genes = "CD163"
x <- t(t(rld)[,rpost_genes])

expdata <- x
Pop1 <- rpost_genes
pops<-list(Pop1)
#Z-Scores



z_scores = t(expdata)
z_scores = data.frame(z_scores[,1])
z_scores$sample = rownames(z_scores)



EF_value2 = read.xlsx(publication_path("HEART_SNRNA_ROOT", "data/2024-04-20.Stravos/Ef/supplementarytable2.xlsx"))

colnames(z_scores) = c("z","sample")
Pre = z_scores[grepl(z_scores$sample , pattern = "pre"),]
Pre = cbind(Pre,EF_value2)
Pre$EF = Pre$Baseline.LVEF
Post = z_scores[grepl(z_scores$sample , pattern = "post"),]
Post = cbind(Post,EF_value2)

Post$EF = Post$`Post-LVAD.LVEF.(max)`


plot_df = rbind(Pre,Post)




print(ggplot(plot_df, aes(x = EF, y = z))+geom_point()+ geom_smooth(method = "lm")+
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top")+ggtitle(paste(cluster,rpost_genes)))
 















#ReadIn M1Vs M2 markers 
HeartMacro = readRDS(publication_path("HEART_SNRNA_ROOT", "data/2022-10-17.SharedObject/2023-10-12_Macrophages_RNA_TooManyAnnotated.rds"))

Idents(HeartMacro) = "CellType"
Markers = FindAllMarkers(HeartMacro ,only.pos = TRUE, min.pct = 0.1)

# Seurat v4-style call within a v5 environment
Markers <- FindAllMarkers(
  object = HeartMacro,
  only.pos = TRUE,
  logfc.threshold = 0.1,
  test.use = "wilcox",        # Standard Wilcoxon
  base.helper = NULL,         # Disables the v5 group-level pseudocount logic
  mean.fxn = rowMeans,        # Forces cell-level averaging
  fc.name = "avg_log2FC"      # Ensures Base 2 is used
)



DT::datatable(Markers, filter = 'top')



# Map Mouse to Human
reticulate::py_require("matplotlib")
reticulate::py_require("scvelo")
library(seuratTools)
 dest_symbols <- Markers$gene %>%
            tibble::enframe("gene_index", "MGI.symbol") %>%
            dplyr::left_join(human_to_mouse_homologs, by = "MGI.symbol") 


map = dest_symbols[(dest_symbols$HGNC.symbol %in% rownames(rld)),]
All_Human = merge(Markers,map, by.x = "gene",by.y ="MGI.symbol")

top20 <- All_Human %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)

sessionInfo()

theme_mytheme <- theme(
  # Specify the overall plot appearance
  plot.background = element_blank(),
  plot.title = element_text(size = 16, face = "bold"),
  plot.subtitle = element_text(size = 14),
  
  # Specify the axis appearance
  axis.text = element_text(size = 12, color = "black"),
  axis.title = element_text(size = 14, color = "black"),
  axis.line = element_line(color = "black", size = 0.5),
  axis.ticks = element_line(color = "black", size = 0.5),
  
  # Specify the legend appearance
  #legend.title = element_blank(),
  legend.text = element_text(size = 12),
  legend.position = "right",
  
  # Specify the panel/grid appearance, no grid line 
  panel.background = element_blank(),#element_rect(fill = "white", color = "black"),
  panel.grid.major = element_blank(),#element_line(color = "lightgray", size = 0.2),
  panel.grid.minor = element_blank(),#element_line(color = "lightgray", size = 0.2),
  
  # Specify the plot margins
  #plot.margin = margin(1, 1, 1, 1, "cm"),
  
  #set the ratio of the plot 
  #aspect.ratio=1
)

# M1 Markers
top20 <- All_Human %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)
cluster = 'Myeloid'
cluster_num <- 2

# Subset the metadata to only the desired cells
cluster_metadata <- metadata[which(metadata$cluster_id == cluster), ]

# Assign the rownames of the metadata to be the sample IDs
rownames(cluster_metadata) <- cluster_metadata$sample_id

# Subset the counts to only the desired cells
counts <- pb[[cluster]]

cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
rownames(cluster_metadata) = str_replace_all(rownames(cluster_metadata),"-",".")
# Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
all(rownames(cluster_metadata) == colnames(cluster_counts))   
library(DESeq2)
rld <- rlog(as.matrix(cluster_counts))

#Genes that are Up in RRPost vs PreLvad
#Genes that are Up in RRPost vs donr
#Genes that aren't up in  RRPost vs donr
rpost_genes = unlist(top20[top20$cluster == "M1",'HGNC.symbol'])
x <- t(t(rld)[,rpost_genes])

expdata <- x
Pop1 <- rpost_genes
pops<-list(Pop1)
#Z-Scores
z_scores<-NULL

for (i in 1:length(pops)) {
genes <- pops[[i]]
zz <- which(tolower(rownames(expdata)) %in% tolower(genes))
av <- numeric(ncol(expdata))

geneExp <- as.matrix(expdata[zz, ])
geneExp <- t(scale(t(geneExp)))
geneExp[is.nan(geneExp)] <- 0
z_scores <- rbind(z_scores,(av + colSums(geneExp) / length(zz)))
}
save = z_scores

z_scores = t(save)
z_scores = data.frame(z_scores[,1])
z_scores$sample = rownames(z_scores)



EF_value2 = read.xlsx(publication_path("HEART_SNRNA_ROOT", "data/2024-04-20.Stravos/Ef/supplementarytable2.xlsx"))

colnames(z_scores) = c("z","sample")
Pre = z_scores[grepl(z_scores$sample , pattern = "pre"),]
Pre = cbind(Pre,EF_value2)
Pre$EF = Pre$Baseline.LVEF
Post = z_scores[grepl(z_scores$sample , pattern = "post"),]
Post = cbind(Post,EF_value2)

Post$EF = Post$`Post-LVAD.LVEF.(max)`


plot_df = rbind(Pre,Post)




print(ggplot(plot_df, aes(x = EF, y = z))+geom_point()+ geom_smooth(method = "lm")+theme_mytheme+ ylab("M1 Score") +
    stat_cor(method = "pearson", label.x.npc = 0.7, label.y.npc = "bottom")+ggtitle(paste(cluster,"M1 markers")))
 #ggsave("M1_Markers_Ef.pdf")

M1_score = plot_df





# M2 Markers
top20 <- All_Human %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)

cluster = 'Myeloid'
cluster_num <- 2

# Subset the metadata to only the desired cells
cluster_metadata <- metadata[which(metadata$cluster_id == cluster), ]

# Assign the rownames of the metadata to be the sample IDs
rownames(cluster_metadata) <- cluster_metadata$sample_id

# Subset the counts to only the desired cells
counts <- pb[[cluster]]

cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
rownames(cluster_metadata) = str_replace_all(rownames(cluster_metadata),"-",".")
# Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
all(rownames(cluster_metadata) == colnames(cluster_counts))   
library(DESeq2)
rld <- rlog(as.matrix(cluster_counts))

#Genes that are Up in RRPost vs PreLvad
#Genes that are Up in RRPost vs donr
#Genes that aren't up in  RRPost vs donr
rpost_genes = unlist(top20[top20$cluster == "M2",'HGNC.symbol'])
x <- t(t(rld)[,rpost_genes])

expdata <- x
Pop1 <- rpost_genes
pops<-list(Pop1)
#Z-Scores
z_scores<-NULL

for (i in 1:length(pops)) {
genes <- pops[[i]]
zz <- which(tolower(rownames(expdata)) %in% tolower(genes))
av <- numeric(ncol(expdata))

geneExp <- as.matrix(expdata[zz, ])
geneExp <- t(scale(t(geneExp)))
geneExp[is.nan(geneExp)] <- 0
z_scores <- rbind(z_scores,(av + colSums(geneExp) / length(zz)))
}
save = z_scores

z_scores = t(save)
z_scores = data.frame(z_scores[,1])
z_scores$sample = rownames(z_scores)



EF_value2 = read.xlsx(publication_path("HEART_SNRNA_ROOT", "data/2024-04-20.Stravos/Ef/supplementarytable2.xlsx"))

colnames(z_scores) = c("z","sample")
Pre = z_scores[grepl(z_scores$sample , pattern = "pre"),]
Pre = cbind(Pre,EF_value2)
Pre$EF = Pre$Baseline.LVEF
Post = z_scores[grepl(z_scores$sample , pattern = "post"),]
Post = cbind(Post,EF_value2)

Post$EF = Post$`Post-LVAD.LVEF.(max)`


plot_df = rbind(Pre,Post)




print(ggplot(plot_df, aes(x = EF, y = z))+geom_point()+ geom_smooth(method = "lm")+theme_mytheme+ ylab("M2 Score") +
    stat_cor(method = "pearson", label.x.npc = 0.7, label.y.npc = "bottom")+ggtitle(paste(cluster,"M2 markers")))
   #ggsave("M2.pdf")
M2_score = plot_df

# M1 Markers
top20 <- All_Human %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)
cluster = 'Myeloid'
cluster_num <- 2

# Subset the metadata to only the desired cells
cluster_metadata <- metadata[which(metadata$cluster_id == cluster), ]

# Assign the rownames of the metadata to be the sample IDs
rownames(cluster_metadata) <- cluster_metadata$sample_id

# Subset the counts to only the desired cells
counts <- pb[[cluster]]

cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
rownames(cluster_metadata) = str_replace_all(rownames(cluster_metadata),"-",".")
# Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
all(rownames(cluster_metadata) == colnames(cluster_counts))   
library(DESeq2)
rld <- rlog(as.matrix(cluster_counts))

#Genes that are Up in RRPost vs PreLvad
#Genes that are Up in RRPost vs donr
#Genes that aren't up in  RRPost vs donr
rpost_genes = unlist(top20[top20$cluster == "M_Ex",'HGNC.symbol'])
x <- t(t(rld)[,rpost_genes])

expdata <- x
Pop1 <- rpost_genes
pops<-list(Pop1)
#Z-Scores
z_scores<-NULL

for (i in 1:length(pops)) {
genes <- pops[[i]]
zz <- which(tolower(rownames(expdata)) %in% tolower(genes))
av <- numeric(ncol(expdata))

geneExp <- as.matrix(expdata[zz, ])
geneExp <- t(scale(t(geneExp)))
geneExp[is.nan(geneExp)] <- 0
z_scores <- rbind(z_scores,(av + colSums(geneExp) / length(zz)))
}
save = z_scores

z_scores = t(save)
z_scores = data.frame(z_scores[,1])
z_scores$sample = rownames(z_scores)



EF_value2 = read.xlsx(publication_path("HEART_SNRNA_ROOT", "data/2024-04-20.Stravos/Ef/supplementarytable2.xlsx"))

colnames(z_scores) = c("z","sample")
Pre = z_scores[grepl(z_scores$sample , pattern = "pre"),]
Pre = cbind(Pre,EF_value2)
Pre$EF = Pre$Baseline.LVEF
Post = z_scores[grepl(z_scores$sample , pattern = "post"),]
Post = cbind(Post,EF_value2)

Post$EF = Post$`Post-LVAD.LVEF.(max)`


plot_df = rbind(Pre,Post)




print(ggplot(plot_df, aes(x = EF, y = z))+geom_point()+ geom_smooth(method = "lm")+theme_mytheme+ ylab("MacFib Score") +
    stat_cor(method = "pearson", label.x.npc = 0.7, label.y.npc = "bottom")+ggtitle(paste(cluster,"MacFib markers")))
 #ggsave("MacFib.pdf")
MacFib_score = plot_df



#  Human single-cell integration — macrophage signature scores Per human cell: cell_id, dataset, macrophage_state_score (Cd163hi, MHC-IIhi, MacFib), ejection_fraction (where available), Pearson/Spearman correlation coefficients. MAKE

#Individual, State Score, Ejection fraction, 
M1_score$M1_score = M1_score$z
M2_score$M2_score = M2_score$z
MacFib_score$MacFib_score = MacFib_score$z

M1_score
M2_score
MacFib_score

library(dplyr)
library(tibble)

# Convert rownames to a column, join, then convert back to rownames
combined_scores <- M1_score %>%
  rownames_to_column("ID") %>%
  left_join(rownames_to_column(M2_score, "ID"), by = "ID") %>%
  left_join(rownames_to_column(MacFib_score, "ID"), by = "ID") %>%
  column_to_rownames("ID")

colnames(combined_scores)
combined_scores = combined_scores[,c("sample",'EF','MacFib_score','M1_score','M2_score')]
write.csv(combined_scores, "EF_macScores_summary.csv")
getwd()
