# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Apply recorded doublet calls to the ATAC object and rebuild embeddings.
# Steps: 1) Load annotated ATAC data and doublet labels; 2) retain singlets; 3) rerun TF-IDF/LSI and save filtered outputs.
# =======================================
# Source reformatted from snATAC_AgeEx/ATACDOUBLET_filter.ipynb.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

library(scDblFinder)
library(Signac)
library(Seurat)
library(BiocParallel)
library(ggplot2)

setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-02-03.MACS/"))
merged = readRDS('MergedAtac_QC_annotated.rds')



head(Doublets)

merged$doublets = Doublets$x

merged$doublets = Doublets$x
DimPlot(merged, group.by = "doublets")

DimPlot(merged, group.by = "doublets", split.by = "seurat_clusters")

ggsave("DoubleClusters.png", width = 49)

Filtered = subset(merged, doublets =="singlet")

Filtered <- RunTFIDF(Filtered)
Filtered <- FindTopFeatures(Filtered, min.cutoff = 20)
Filtered <- RunSVD(Filtered)

DepthCor(Filtered)

Filtered <- RunUMAP(
  object = Filtered,
  reduction = 'lsi',
  dims = 2:50
)

DimPlot(object = Filtered, label = TRUE) + NoLegend()
DimPlot(object = Filtered, label = TRUE, group.by = "doublets") + NoLegend()

DefaultAssay(Filtered) = "ATAC"
Filtered <- FindNeighbors(object = Filtered, reduction = 'lsi', dims = 2:50)
Filtered <- FindClusters(object = Filtered, verbose = TRUE, algorithm = 3)

DimPlot(object = Filtered, label = TRUE) + NoLegend()

DefaultAssay(Filtered) ="All_RNA"
print('EC_marker_feature_plot')
DotPlot(Filtered, features = c("Kdr", "Pecam1","Ly6a", "Fabp4"))


print( 'Fibroblast_marker_feature_plot')
DotPlot(Filtered, features = c("Col1a1", "Pdgfra","Ddr2", "Tcf21"))

print('CM_marker_feature_plot')
DotPlot(Filtered, features = c("Tnnt2", "Ttn","Myh6", "Tnni3"))

print('Mural_marker_feature')
DotPlot(Filtered, features = c("Pdgfrb", "Vtn","Acta2", "Rgs5"))


print('ribo_marker')
DotPlot(Filtered, features = c("Rpl7a", "Rpl26","Rps3","Rps6"))


print('leukocytes_marker')
DotPlot(Filtered, features = c("Ptprc"))


print('TCELL_marker')
DotPlot(Filtered, features = c("Cd3e","Cd4","Cd8a"  ))


print('BCELL_marker')
DotPlot(Filtered, features = c("Cd19","Cd20"  ))


print('macrophages_marker')
DotPlot(Filtered, features = c("Itgam","Itgax","Ly6c1","Cd163"  ))


print('neutrophils_marker')
DotPlot(Filtered, features = c("Ly6g"))


print('Adipocytes_marker')
DotPlot(Filtered, features = c("Gpam","Fasn","Lep"))


print('AtrialCM')
DotPlot(Filtered, features = c("Nppa","Myl7","Myl4"))


print('Endothilial')
DotPlot(Filtered, features = c("Vwf","Pecam1","Cdh5"))


print('Fibro_marker')
DotPlot(Filtered, features = c("Dcn","Gsn","Pdgfra"))


print('lymphoid_marker')
DotPlot(Filtered, features = c("Cd8a","Il7r","Cd40lg"))


print('mesothelial_marker')
DotPlot(Filtered, features = c("Msln","Wt1","Bnc1"))


print('myeloid_marker')
DotPlot(Filtered, features = c("Cd14","C1qa","Cd68"))


print('Neuronal_marker')
DotPlot(Filtered, features = c("Plp1","Nrxn1","Nrxn3"))


print('Pericytes_marker')
DotPlot(Filtered, features = c("Rgs5","Abcc9","Kcnj8"))


print( 'SmoothMuscle_marker')
DotPlot(Filtered, features = c("Myh11","Tagln","Acta2"))


print('VentricalCM')
DotPlot(Filtered, features = c("Myh7","Myl2","Fhl2"))

setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-03-01.CellbenderFiltered_MAS3_2_Adjust"))
RNA = readRDS('CellTypes_M1_M2.seurat.rds')
Immune = readRDS(publication_path("HEART_SNRNA_ROOT", "analysis/2023-03-15.Macro/Immune_clustering.rds"))

table(Immune$Too_Many_Type)
head(RNA$cellID)
head(RNA@meta.data[RNA$cellID %in% Immune$cellID,'Too_many'])

RNA$Too_many = NA
RNA@meta.data[RNA$cellID %in% Immune$cellID,'Too_many'] = Immune$Too_Many_Type
RNA@meta.data[!(RNA$cellID %in% Immune$cellID),'Too_many'] = RNA@meta.data[!(RNA$cellID %in% Immune$cellID),'cellType2']
table(RNA$Too_many)
table(Immune$Too_Many_Type)

transfer.anchors = readRDS("Doublet_anchor.rds")

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = RNA$Too_many,
  weight.reduction = Filtered[['lsi']],
  dims = 2:50
)


Filtered <- AddMetaData(object = Filtered, metadata = predicted.labels)
Filtered$predicted.cellType_Doublet = Filtered$predicted.id

saveRDS(Filtered, "DoubletFiltered_ATAC.rds")

setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-02-03.MACS/"))
plot1 <- DimPlot(
  object = RNA,
  group.by = 'Too_many',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = Filtered,
  group.by = 'predicted.cellType_Doublet',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
plot1 + plot2
ggsave(file = "Cellbender_IntegratedcellType_Doublet.png", width = 20)

setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-02-03.MACS/"))
plot1 <- DimPlot(
  object = Filtered,
  group.by = 'seurat_clusters',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot2 <- DimPlot(
  object = Filtered,
  group.by = 'predicted.cellType_Doublet',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
plot1 + plot2
ggsave(file = "Cellbender_IntegratedcellType_Doublet_clusters.png", width = 20)

setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-02-03.MACS/"))
Filtered = readRDS("DoubletFiltered_celltype_ATAC.rds")


CM = c(1,2,4,5,23)
Coronary_EC = c(7,10,16)
Endocardial_EC = c(9)
Lympathic_EC = c(14)
Fibroblast= c(0,3,11,12,13)
Myeloid = c(6,20)
Bcell = 15
Tcell = 17
Cardiac_Neuronal = 22 
Mesothelial = 18
Pericyte = 8
smooth_muscle = 19

ATAC_unknown3 = 21


Filtered[["cellType_Doublet"]] <- Filtered[["seurat_clusters"]]
Filtered[["cellType_Doublet"]]  <- ifelse(unlist(Filtered[["seurat_clusters"]]) %in% CM, "Cardiomyocyte",unlist(Filtered[["cellType_Doublet"]]))
Filtered[["cellType_Doublet"]]  <- ifelse(unlist(Filtered[["seurat_clusters"]]) %in% Coronary_EC, "Coronary_EC",unlist(Filtered[["cellType_Doublet"]]))
Filtered[["cellType_Doublet"]]  <- ifelse(unlist(Filtered[["seurat_clusters"]]) %in% Endocardial_EC, "Endocardial_EC",unlist(Filtered[["cellType_Doublet"]]))
Filtered[["cellType_Doublet"]]  <- ifelse(unlist(Filtered[["seurat_clusters"]]) %in% Lympathic_EC, "Lymphatic_EC",unlist(Filtered[["cellType_Doublet"]]))
Filtered[["cellType_Doublet"]]  <- ifelse(unlist(Filtered[["seurat_clusters"]]) %in% Mesothelial, "Mesothelial",unlist(Filtered[["cellType_Doublet"]]))
Filtered[["cellType_Doublet"]]  <- ifelse(unlist(Filtered[["seurat_clusters"]]) %in% Fibroblast, "Fibroblast",unlist(Filtered[["cellType_Doublet"]]))
Filtered[["cellType_Doublet"]]  <- ifelse(unlist(Filtered[["seurat_clusters"]]) %in% Myeloid, "Myeloid",unlist(Filtered[["cellType_Doublet"]]))
Filtered[["cellType_Doublet"]]  <- ifelse(unlist(Filtered[["seurat_clusters"]]) %in% Tcell, "T_cell",unlist(Filtered[["cellType_Doublet"]]))
Filtered[["cellType_Doublet"]]  <- ifelse(unlist(Filtered[["seurat_clusters"]]) %in% Bcell, "B_cell",unlist(Filtered[["cellType_Doublet"]]))
Filtered[["cellType_Doublet"]]  <- ifelse(unlist(Filtered[["seurat_clusters"]]) %in% Pericyte, "Pericyte",unlist(Filtered[["cellType_Doublet"]]))
Filtered[["cellType_Doublet"]]  <- ifelse(unlist(Filtered[["seurat_clusters"]]) %in% smooth_muscle, "Smooth_muscle",unlist(Filtered[["cellType_Doublet"]]))
Filtered[["cellType_Doublet"]]  <- ifelse(unlist(Filtered[["seurat_clusters"]]) %in% Cardiac_Neuronal, "Cardiac_Neuronal",unlist(Filtered[["cellType_Doublet"]]))



Filtered[["cellType_Doublet"]]  <- ifelse(unlist(Filtered[["seurat_clusters"]]) %in% ATAC_unknown3, "ATAC_unknown3",unlist(Filtered[["cellType_Doublet"]]))


saveRDS(Filtered, "DoubletFiltered_celltype2_ATAC.rds")

setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-02-03.MACS/"))
Filtered = readRDS("DoubletFiltered_celltype_ATAC.rds")

Filtered
table(Filtered$ind)



table(Filtered$sex, useNA = "always")
Filtered_sex = subset(Filtered, !is.na(sex))
table(Filtered_sex$sex, useNA = "always")

table(Filtered$sex, useNA = "always")

saveRDS(Filtered, "DoubletFiltered_Sex_celltype_ATAC.rds")

