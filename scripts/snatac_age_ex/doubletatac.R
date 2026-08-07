# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Detect candidate snATAC doublets with scDblFinder.
# Steps: 1) Convert the ATAC object to SingleCellExperiment; 2) run scDblFinder; 3) inspect and save doublet classifications.
# =======================================
# Source reformatted from snATAC_AgeEx/DoubletATAC.ipynb.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

#Doublet FInder

sessionInfo()

library(scDblFinder)
library(Signac)
library(Seurat)
library(BiocParallel)

setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-02-03.MACS/"))
merged = readRDS('MergedAtac_QC_annotated.rds')

merged.sce <- as.SingleCellExperiment(merged)

sce <- scDblFinder(merged.sce, clusters = 'seurat_clusters', samples = 'run', aggregateFeatures = TRUE, returnType = "full")

table(sce$scDblFinder.class)

sce.seurat <- as.Seurat(sce, counts = "counts", data = "logcounts")

DimPlot(sce.seurat, group.by = "scDblFinder.class")

DimPlot(sce.seurat, group.by = "scDblFinder.class", split.by = "seurat_clusters")
library(ggplot2)
ggsave("DoubleClusters.png", width = 49)

getwd()
write.csv(cbind(colnames(sce),sce$scDblFinder.class), "Doublets.csv")
keep = colnames(sce)[sce$scDblFinder.class == "singlet"]

merged$CellID = colnames(merged)
Filtered = subset(merged, CellID %in% keep)

saveRDS(Filtered, "ATAC_Doublet_filtered.rds")
