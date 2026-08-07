# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Generate global Cicero co-accessibility connections from the annotated snATAC-seq object.
# Steps: 1) Convert the global ATAC object to a Cicero cell-data set; 2) infer genomic co-accessibility connections; 3) save the Cicero object, connections and co-accessibility networks for downstream link analyses.
# =======================================
# Source reformatted from snATAC_AgeEx/cicero-Copy1.ipynb.
# The required source was stored across code, raw and Markdown notebook cells; only the confirmed global connection-generation path is retained here.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

library(Seurat)
library(cicero)
library(Signac)
library(SeuratWrappers)

combined = readRDS(publication_path("HEART_SNRNA_ROOT", "data/2022-10-17.SharedObject/2023-10-11_Full_ATAC_TooManyAnnotated.rds"))

DefaultAssay( combined) = "ATAC"

ATAC.cds <- as.cell_data_set(x = combined)
ATAC.cicero <- make_cicero_cds(ATAC.cds, reduced_coordinates = reducedDims(ATAC.cds)$UMAP)

setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-10-31.Cicero"))

saveRDS(ATAC.cicero,'ATAC.cicero.rds')

genome <- seqlengths(Annotation(combined))
# use chromosome 1 to save some timeb
# omit this step to run on the whole genome


# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)

# run cicero
conns <- run_cicero(ATAC.cicero, genomic_coords = genome.df, sample_num = 100)

saveRDS(conns,"conns.rds")

ccans <- generate_ccans(conns)

saveRDS(ccans,"ccans.rds")
