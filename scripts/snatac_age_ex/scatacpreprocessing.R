# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Preprocess snATAC-seq libraries and construct the annotated atlas.
# Steps: 1) Load Cell Ranger peaks/matrices; 2) create Signac objects and QC embeddings; 3) transfer labels and save the combined object.
# =======================================
# Source reformatted from snATAC_AgeEx/scATACPreProcessing (1).Rmd.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

library(Signac)
library(EnsDb.Mmusculus.v79)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
library(Seurat)
library(future)
library(ggplot2)


annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"
setwd(publication_path("HEART_SNRNA_ROOT", "data/2022-11-08.ATACObject"))
saveRDS(annotations,"annotations.rds")

setwd(publication_path("HEART_SNRNA_ROOT", "data/2022-11-08.ATACObject"))
annotations = readRDS("annotations.rds")
plan(if (future::supportsMulticore()) "multicore" else "multisession", workers = 32)  # "multiprocess" is defunct in current 'future'
options(future.globals.maxSize = 100000 * 1024^2)

dir = publication_path("SCATAC_CELLRANGER_ROOT")
readPeaks = function(x){
  return( read.table(file = paste0(dir, x,"/outs/peaks.bed"),
  col.names = c("chr", "start", "end")
))
}

set = c("MAE2_1",  "MAE4_3",  "MAS1_4", "MAS3_2_v2",  "YE2_1",  "YE4_3",  "YS1_4",  "YS3_2")
peaks = list()
gr= list()
for (i in 1:length(set)){
peaks[[i]] <- readPeaks(set[i])
gr[[i]] = makeGRangesFromDataFrame(peaks[[i]])

}



combined.peaks <- reduce(x = c(gr[[1]],gr[[2]],gr[[3]],gr[[4]],gr[[5]],gr[[6]],gr[[7]],gr[[8]]))

peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks


readMeta = function(x){
  return( read.table(file = paste0(dir, x,"/outs/singlecell.csv"),
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]
)
}

meta = list()
for (i in 1:length(set)){
meta[[i]] <- readMeta(set[i])
print(summary(meta[[i]]$passed_filters))
meta[[i]] <- meta[[i]][meta[[i]]$passed_filters >500,]
print(summary(meta[[i]]$passed_filters))
}

makeFrags = function(x){
  return( CreateFragmentObject(
  path = paste0(dir, x,"/outs/fragments.tsv.gz"),
  cells = rownames(meta[[i]]))

)
}

frags = list()
for (i in 1:length(set)){
frags[[i]] <- makeFrags(set[i])

}



makeMatrix = function(i){
  return( FeatureMatrix(
  fragments = frags[[i]],
  features = combined.peaks,
  cells = rownames(meta[[i]])
))
}
Matrix = list()
for (i in 1:length(set)){
Matrix[[i]] <- makeMatrix(i)

}
setwd(publication_path("HEART_SNRNA_ROOT", "data/2022-11-08.ATACObject"))
saveRDS(Matrix, "AllMatrix.rds")



assay = list()
object = list()
for (i in 1:length(set)){
  assay[[i]] <- CreateChromatinAssay(Matrix[[i]], fragments = frags[[i]])
object[[i]] <- CreateSeuratObject(assay[[i]], assay = "ATAC", meta.data=meta[[i]])
  
  
}


setwd(publication_path("HEART_SNRNA_ROOT", "data/2022-11-08.ATACObject"))
for (i in 1:length(set)){
  Annotation(object[[i]]) <- annotations
  object[[i]] <- NucleosomeSignal(object = object[[i]])

  # compute TSS enrichment score per cell
  object[[i]] <- TSSEnrichment(object = object[[i]], fast = FALSE)
  
  # add blacklist ratio and fraction of reads in peaks
  object[[i]]$pct_reads_in_peaks <- object[[i]]$peak_region_fragments / object[[i]]$passed_filters * 100
  object[[i]]$blacklist_ratio <- object[[i]]$blacklist_region_fragments / object[[i]]$peak_region_fragments
  
  object[[i]]$high.tss <- ifelse(object[[i]]$TSS.enrichment > 2, 'High', 'Low')
  TSSPlot(object[[i]], group.by = 'high.tss') + NoLegend()
  ggsave(paste0(set[i],".TSS.png"), width = 16)
  object[[i]]$nucleosome_group <- ifelse(object[[i]]$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  FragmentHistogram(object = object[[i]], group.by = 'nucleosome_group',region = "chr1-1-10000000")
  ggsave(paste0(set[i],".Fragment.png"), width = 16)
  VlnPlot(
  object = object[[i]],
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
  )
  ggsave(paste0(set[i],".Vln.png"), width = 16)
  
  
}

setwd(publication_path("HEART_SNRNA_ROOT", "data/2022-11-08.ATACObject"))
saveRDS(object, "Allobject.rds")




setwd(publication_path("HEART_SNRNA_ROOT", "data/2022-11-08.ATACObject"))
object = readRDS("Allobject.rds")
peakFragmentsFilters= c(35000,35000,35000,30000,15000,15000,25000,30000)
  
filter = list()
for (i in 1:length(set)){
  filter[[i]] <- subset(
    x = object[[i]],
    subset = peak_region_fragments > 3000 &
    peak_region_fragments < peakFragmentsFilters[i] &
      pct_reads_in_peaks > 30 &
      blacklist_ratio < 0.05 &
      nucleosome_signal < 4 &
      TSS.enrichment > 2
  )
   VlnPlot(
  object = filter[[i]],
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
  )
  ggsave(paste0(set[i],"Filter.Vln.png"), width = 16)
}


for (i in 1:length(set)){
  filter[[i]]$dataset = set[i]
  
}



combined <- merge(
  x =  filter[[1]],
  y =  c(filter[[2]],filter[[3]],filter[[4]],filter[[5]],filter[[6]],filter[[7]],filter[[8]]),
  add.cell.ids = set
)
combined[["ATAC"]]

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)

DepthCor(combined)
ggsave("DepthCor.png")

combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')

DimPlot(combined, split.by = 'dataset', pt.size = 0.1)
ggsave("CombinedUmap_split_2:50.png", width = 20)

combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:50)
combined <- FindClusters(object = combined, verbose = TRUE, algorithm = 3)
DimPlot(object = combined, label = TRUE,split.by = 'dataset') + NoLegend()
ggsave("CombinedUmap_clusters.png", width = 20)


gene.activities <- GeneActivity(combined)


combined[['RNA']] <- CreateAssayObject(counts = gene.activities)
combined <- NormalizeData(
  object = combined,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined$nCount_RNA)
)
saveRDS(combined,"combined.RDS")




ggsave(FeaturePlot(combined, reduction = "umap", features = c("Kdr", "Pecam1","Ly6a", "Fabp4"),  ncol =2, cols = c("grey", "red")), file = 'EC_marker_feature_plot.png')
ggsave(FeaturePlot(combined, reduction = "umap", features = c("Col1a1", "Pdgfra","Ddr2", "Tcf21"),  ncol =2, cols = c("grey", "red")), file = 'Fibroblast_marker_feature_plot.png')
ggsave(FeaturePlot(combined, reduction = "umap", features = c("Tnnt2", "Ttn","Myh6", "Tnni3"),  ncol =2, cols = c("grey", "red")), file= 'CM_marker_feature_plot.png')
ggsave(FeaturePlot(combined, reduction = "umap", features = c("Pdgfrb", "Vtn","Acta2", "Rgs5"),  ncol =2, cols = c("grey", "red")), file = 'Mural_marker_feature_plot.png')
ggsave(FeaturePlot(combined, reduction = "umap", features = c("Rpl7a", "Rpl26","Rps3","Rps6"),  ncol =2, cols = c("grey", "red")), file = 'ribo_marker_feature_plot.png')
ggsave(FeaturePlot(combined, reduction = "umap", features = c("Ptprc"),  ncol =2, cols = c("grey", "red")), file = 'leukocytes_marker_feature_plot.png')
ggsave(FeaturePlot(combined, reduction = "umap", features = c("Cd3e","Cd4","Cd8a"  ),  ncol =2, cols = c("grey", "red")), file = 'TCELL_marker_feature_plot.png')
ggsave(FeaturePlot(combined, reduction = "umap", features = c("Cd19","Cd20"  ),  ncol =2, cols = c("grey", "red")), file = 'BCELL_marker_feature_plot.png')
ggsave(FeaturePlot(combined, reduction = "umap", features = c("Itgam","Itgax","Ly6c1"  ),  ncol =2, cols = c("grey", "red")), file = 'macrophages_marker_feature_plot.png')
ggsave(FeaturePlot(combined, reduction = "umap", features = c("Ly6g"),  ncol =2, cols = c("grey", "red")), file = 'neutrophils_marker_feature_plot.png')

ggsave(FeaturePlot(combined, reduction = "umap", features = c("Gpam","Fasn","Lep"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'Adipocytes_marker_feature_plot.png')
ggsave(FeaturePlot(combined, reduction = "umap", features = c("Nppa","Myl7","Myl4"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'AtrialCM_marker_feature_plot.png')
ggsave(FeaturePlot(combined, reduction = "umap", features = c("Vwf","Pecam1","Cdh5"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'Endothilial_marker_feature_plot.png')
ggsave(FeaturePlot(combined, reduction = "umap", features = c("Dcn","Gsn","Pdgfra"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'Fibro_marker_feature_plot.png')
ggsave(FeaturePlot(combined, reduction = "umap", features = c("Cd8a","Il7r","Cd40lg"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'lymphoid_marker_feature_plot.png')
ggsave(FeaturePlot(combined, reduction = "umap", features = c("Msln","Wt1","Bnc1"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'mesothelial_marker_feature_plot.png')
ggsave(FeaturePlot(combined, reduction = "umap", features = c("Cd14","C1qa","Cd68"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'myeloid_marker_feature_plot.png')
ggsave(FeaturePlot(combined, reduction = "umap", features = c("Plp1","Nrxn1","Nrxn3"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'Neuronal_marker_feature_plot.png')
ggsave(FeaturePlot(combined, reduction = "umap", features = c("Rgs5","Abcc9","Kcnj8"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'Pericytes_marker_feature_plot.png')
ggsave(FeaturePlot(combined, reduction = "umap", features = c("Myh11","Tagln","Acta2"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'SmoothMuscle_marker_feature_plot.png')
ggsave(FeaturePlot(combined, reduction = "umap", features = c("Myh7","Myl2","Fhl2"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'VentricalCM_marker_feature_plot.png')

ggsave(FeaturePlot(combined, reduction = "umap", features = c("Gpam","Fasn","Lep",
"Vwf","Pecam1","Cdh5",
"Dcn","Gsn","Pdgfra",
"Cd8a","Il7r","Cd40lg",
"Msln","Wt1","Bnc1",
"Cd14","C1qa","Cd68",
"Plp1","Nrxn1","Nrxn3",
"Rgs5","Abcc9","Kcnj8",
"Myh11","Tagln","Acta2",
"Myh7","Myl2","Fhl2"
),  ncol =9, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'All_marker_feature_plot.png', height=32, width = 48 )


ggsave(FeaturePlot(combined, reduction = "umap", features = c("Cav1","Nr2f2","Thbd","Stard10","Vsnl1","Nr2f1","Gm15543","Myh6"),  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'AtrialScienceCM_marker_feature_plot.png')
ggsave(FeaturePlot(combined, reduction = "umap", features = c("Myl2","Pln","Mpped2","Sall4","Capn6","Hey2","Myh7"),  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'VentricularScienceCM_marker_feature_plot.png')

Artery = c("Stmn2", "Fos", "Id1", "Sox17", "Klf4", "Alpl", "Junb", "Crip1", "Mecom","Hspa1a")
Capillary_arterial = c("Cxcl12", "Rbp7", "Mgll", "Ly6c1", "Aqp7", "Btnl9", "Pdgfd","Sema7a", "Clec2d","Unc5b")
Capillary = c("Kdr", "Ssh2", "AU020206", "Fam57b", "Kifc3", "Ntf3", "Tmem182", "Ebf3","Coro2b", "Fam117b")
Capillary_venous = c("Ier3", "Pltp", "Fmo2", "Fmo1", "Calcrl", "Emp1", "AU021092", "Bsg","Nfkbia", "Eln")
Lymphatic = c("Ccl21a", "Mmrn1", "Fth1", "Lyve1", "Prelp", "Fxyd6", "S100a6", "Cd9","Pard6g", "Flt4")
endocardial_EC = c("Kank1", "Gpm6a", "Smoc1", "Cdh11", "Plagl1", "Tmem108", "Npr3", "Comp","Ahsg", "Mgp")
coronary_EC = c("Sparcl1", "Mgll", "Ly6c1", "Tspan13", "Aqp1", "Cldn5", "Timp4")
aorta_EC = c("Ehd3", "Fam167b", "Esm1", "AU021092", "Ppap2a")


ggsave(FeaturePlot(combined, reduction = "umap", features =Artery,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'Artery_Endo_marker_feature_plot.png', height = 15, width = 15)
ggsave(FeaturePlot(combined, reduction = "umap", features =Capillary_arterial,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'Capillary_arterial_Endo_marker_feature_plot.png', height = 15, width = 15)
ggsave(FeaturePlot(combined, reduction = "umap", features =Capillary,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'Capillary_Endo_marker_feature_plot.png', height = 15, width = 15)
ggsave(FeaturePlot(combined, reduction = "umap", features =Capillary_venous,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'Capillary_venous_Endo_marker_feature_plot.png', height = 15, width = 15)
ggsave(FeaturePlot(combined, reduction = "umap", features =Lymphatic,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'Lymphatic_Endo_marker_feature_plot.png', height = 15, width = 15)
ggsave(FeaturePlot(combined, reduction = "umap", features =endocardial_EC,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'endocardial_EC_Endo_marker_feature_plot.png', height = 15, width = 15)
ggsave(FeaturePlot(combined, reduction = "umap", features =coronary_EC,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'coronary_EC_Endo_marker_feature_plot.png', height = 15, width = 15)
ggsave(FeaturePlot(combined, reduction = "umap", features =aorta_EC,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'aorta_EC_Endo_marker_feature_plot.png', height = 15, width = 15)



setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2022-08-18.HSP"))
Mouse = readRDS("MouseWithPredictedSexCelltypes.rds")

Validindv =  c("MAE2_1_M","MAE2_1_F","MAS1_4_F","MAS1_4_M","MAS3_2_M","MAS3_2_F",
 "YE2_1_F","YE2_1_M","YE4_3_F","YE4_3_M","YS1_4_M","YS3_2_M","YS3_2_F")
Mouse = subset(Mouse, subset = ind %in% Validindv)

setwd(publication_path("HEART_SNRNA_ROOT", "data/2022-11-08.ATACObject"))
combined = readRDS("combined.RDS")



#label transfer

DefaultAssay(combined) <- 'RNA'

transfer.anchors <- FindTransferAnchors(
  reference = Mouse,
  query = combined,
  reduction = 'cca',
   k.filter = NA,

)




predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = Mouse$cellType,
  weight.reduction = combined[['lsi']],
  dims = 2:30
)

combined <- AddMetaData(object = combined, metadata = predicted.labels)
combined$predicted.cellType = combined$predicted.id

predicted.Cluster.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = Mouse$seurat_clusters,
  weight.reduction = combined[['lsi']],
  dims = 2:30
)

combined <- AddMetaData(object = combined, metadata = predicted.Cluster.labels)
combined$predicted.cluster = combined$predicted.id
combined$predicted.cluster = as.numeric(combined$predicted.cluster)
saveRDS(combined,"combinedWithTransfer.RDS")


setwd(publication_path("HEART_SNRNA_ROOT", "data/2022-11-08.ATACObject"))
combined = readRDS("combinedWithTransfer.RDS")

plot1 <- DimPlot(
  object = Mouse,
  group.by = 'cellType',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = combined,
  group.by = 'predicted.cellType',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot1 + plot2
ggsave(file = "IntegratedCellType.png", width = 20)



plot3 <- DimPlot(
  object = combined,
   split.by = 'predicted.cellType',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
plot3
ggsave(file = "AllCellType.png", width = 40)

Idents(Mouse) = Mouse$ind
plot4 <- DimPlot(
  object = Mouse,
   split.by = 'cellType',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
plot4
ggsave(file = "RNAAllCellType.png", width = 40)

Idents(Mouse) = Mouse$ind
plot5 <- DimPlot(
  object = combined,
   group.by = 'seurat_clusters',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
plot2+plot5
ggsave(file = "Cluster_CellType.png", width = 20)


plot6 <- DimPlot(
  object = combined,
   group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
plot6
ggsave(file = "AllCluster.png")
plot7 <- DimPlot(
  object = Mouse,
  group.by = 'seurat_clusters',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
plot6+plot7
ggsave(file = "ClusterMap.png", width = 15)


setwd(publication_path("LUNG_QC_ROOT", "analysis/2021-08-20.Entropy/"))

load("clean_nodatasets_060720.RData")

setwd(publication_path("HEART_SNRNA_ROOT", "data/2022-11-08.ATACObject"))
countTable = combined@assays$RNA@counts
condition_list = combined@meta.data$seurat_clusters
other = combined@meta.data



data = countTable
  cgenesB = cgenesA[cgenesA %in% rownames(data)]
  xpairs_list = as.data.frame(matrix(unlist(strsplit(xpairs, "_")), nrow=length(xpairs), byrow=T))
  xpairsB = xpairs[xpairs_list$V1 %in% cgenesB & xpairs_list$V2 %in% cgenesB]
  system.time(pdTrain<-query_transform(expT[cgenesB, ], xpairsB))
  system.time(rf_tspAll<-sc_makeClassifier(pdTrain[xpairsB,], genes=xpairsB, groups=grps, nRand=10, ntrees=500))
  print("Done making classifier")
  system.time(kidTransAll<-query_transform(data[cgenesB,], xpairsB))
  print("A")
  system.time(crParkall<-rf_classPredict(rf_tspAll, kidTransAll, numRand=0))
  print("B")
  





combined[["SingleCellNet"]] <- rownames(crParkall)[apply(crParkall,2,which.max)]
  combined[["SingleCellNetScore"]] <- apply(crParkall, 2, max)
  
  

Idents(combined) <- combined$SingleCellNet

setwd(publication_path("HEART_SNRNA_ROOT", "data/2022-11-08.ATACObject"))

tplot <- DimPlot(combined, reduction = "umap", label = T, label.size = 2)
tplot[[1]]$layers[[1]]$aes_params$alpha = .2
ggsave(tplot+ theme(legend.text=element_text(size=3)),file = "PredictedCellType.png", width = 20)



keep = names(table(combined[["SingleCellNet"]])[table(combined[["SingleCellNet"]])> 500])



combined[["SingleCellNetPlot"]] = ifelse((unlist(combined[['SingleCellNet']])) %in% keep  ,unlist(combined[["SingleCellNet"]]) ,'small')

Idents(combined) <- combined$SingleCellNetPlot

tplot <- DimPlot(combined, reduction = "umap", label = T, label.size = 2)
tplot[[1]]$layers[[1]]$aes_params$alpha = .2
ggsave(tplot+ theme(legend.text=element_text(size=3)),file = "PredictedCellTypeSmall.png", width = 20)


keep = subset(combined, subset = SingleCellNetScore > 0.2)

Idents(keep) <- keep$SingleCellNet

tplot <- DimPlot(keep, reduction = "umap", label = T, label.size = 2)
tplot[[1]]$layers[[1]]$aes_params$alpha = .2
ggsave(tplot+ theme(legend.text=element_text(size=3)),file = "PredictedCellTypeConfident.png", width = 20)




#mito


combined$mitopercent=combined$mitochondrial/combined$total
combined$mitopercent2=combined$mitochondrial/(combined$total-combined$unmapped)




setwd(publication_path("HEART_SNRNA_ROOT", "data/2022-11-08.ATACObject"))

VlnPlot(
  object = combined, features = c('mitopercent',"mitopercent2"), split.by = "dataset")




ggsave("scATAC_mitochondrial.png")



#plot  cell proportions
#Y = proportions 
#X = cell type split by datatype 


RNA_prop = (table(Mouse$cellType, Mouse$run))
RNA_prop = prop.table(RNA_prop, margin = 1)
RNA_prop = data.frame(RNA_prop)


ATAC_prop = (table(combined$predicted.cellType, combined$dataset))
ATAC_prop = prop.table(ATAC_prop, margin = 1)
ATAC_prop = data.frame(ATAC_prop)
ATAC_prop$Var2 = as.character(ATAC_prop$Var2)
ATAC_prop[ATAC_prop$Var2 == "MAS3_2_v2",'Var2'] = "MAS3_2"

ATAC_prop$assay = "ATAC"
RNA_prop$assay = "RNA"


plotData = rbind(ATAC_prop,RNA_prop)

setwd(publication_path("HEART_SNRNA_ROOT", "data/2022-11-08.ATACObject"))
ggplot(plotData, aes(x = assay, y = Freq, fill = Var2))+ geom_col()+facet_wrap(~Var1)
ggsave("cellProportion.png")

  RNA_prop = (table(Mouse$cellType, Mouse$run))
  RNA_prop = prop.table(RNA_prop, margin = 1)
  RNA_prop = data.frame(RNA_prop)
  
  
  ATAC_prop = (table(combined$predicted.cellType, combined$dataset))
  ATAC_prop= ATAC_prop[,c(1,3:8)]
  ATAC_prop = prop.table(ATAC_prop, margin = 1)
  ATAC_prop = data.frame(ATAC_prop)
  ATAC_prop$Var2 = as.character(ATAC_prop$Var2)
  ATAC_prop[ATAC_prop$Var2 == "MAS3_2_v2",'Var2'] = "MAS3_2"
  ATAC_prop$assay = "ATAC"
  RNA_prop$assay = "RNA"
  
  
  plotData = rbind(ATAC_prop,RNA_prop)
  
  setwd(publication_path("HEART_SNRNA_ROOT", "data/2022-11-08.ATACObject"))
  ggplot(plotData, aes(x = assay, y = Freq, fill = Var2))+ geom_col()+facet_wrap(~Var1)
  ggsave("cellProportion_sameRuns.png")
  
  
  


#Transfer run label make sure each is the same sample 


setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2022-08-18.HSP"))
Mouse = readRDS("MouseWithPredictedSexCelltypes.rds")

DefaultAssay(combined) <- 'RNA'

transfer.anchors <- FindTransferAnchors(
  reference = Mouse,
  query = combined,
  reduction = 'cca',
   k.filter = NA,

)




predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = Mouse$run,
  weight.reduction = combined[['lsi']],
  dims = 2:30
)

combined <- AddMetaData(object = combined, metadata = predicted.labels)
combined$predicted.run = combined$predicted.id
  
setwd(publication_path("HEART_SNRNA_ROOT", "data/2022-11-08.ATACObject"))
saveRDS(combined,"combinedWithTransfer.RDS")
saveRDS(transfer.anchors,"transfer.anchors.RDS")

table(combined$dataset)
table(combined$predicted.run)

combined$dataset[combined$dataset == "MAS3_2_v2"] = "MAS3_2"
combined$matchedRun = combined$dataset == combined$predicted.run
sub = subset(combined, subset = dataset == "MAE4_3", invert = TRUE)
plot6 <- DimPlot(
  object = sub,
   split.by = 'matchedRun',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
plot6
ggsave(file = "AgreedRun.png", width = 20, height = 10)


data = data.frame(cbind(combined$dataset,combined$predicted.run))
data = data[data$X1 != "MAE4_3",]
data[data$X1 == "MAS3_2_v2","X1"] = "MAS3_2"

table(data$X1 == data$X2)

#FALSE  TRUE
#66242 15796

setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2022-08-18.HSP"))
Mouse = readRDS("MouseWithPredictedSexCelltypes.rds")

DefaultAssay(combined) <- 'RNA'
combined$dataset[combined$dataset == "MAS3_2_v2"] = "MAS3_2"
combined$matchedRun = combined$dataset == combined$predicted.run
sub = subset(combined, subset = dataset == "MAE4_3", invert = TRUE)
transfer.anchors <- FindTransferAnchors(
  reference = Mouse,
  query = sub,
  reduction = 'cca',
   k.filter = NA,

)




predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = Mouse$run,
  weight.reduction = combined[['lsi']],
  dims = 2:30
)

sub <- AddMetaData(object = sub, metadata = predicted.labels)
sub$predicted.run = sub$predicted.id

table(sub$predicted.run == sub$dataset)
#66242 15796


#sex split 

DefaultAssay(combined) <- 'ATAC'
gene.activities <- GeneActivity(combined,biotypes = NULL)


combined[['Expression']] <- CreateAssayObject(counts = gene.activities)
combined <- NormalizeData(
  object = combined,
  assay = 'Expression',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined$nCount_Expression)
)

saveRDS(combined, "Expression.rds")
setwd(publication_path("HEART_SNRNA_ROOT", "data/2022-11-08.ATACObject"))
combined = readRDS("Expression.rds")
DefaultAssay(combined) = "Expression"

#move over whole expression modue 

Xgene.set <-rownames(combined)[which(toupper(rownames(combined)) %in%  Xgenes)]
Ygene.set <-rownames(combined)[which(toupper(rownames(combined)) %in%  Ygenes)]

plot2 <- FeatureScatter(combined, feature1 = "Xist", feature2 = "Uty")
ggsave("sexgenes.png")

sex = FetchData(object = combined, vars = c("Uty", "Xist"))

sex$result = NA
sex[ sex$Uty == sex$Xist & sex$Uty == 0,]$result = "Invalid"
sex[ sex$Uty > 0 &  sex$Xist == 0,]$result = "MALE"
sex[ sex$Uty == 0 &  sex$Xist > 0,]$result = "FEMALE"
sex[ sex$Uty > 0 &  sex$Xist >0,]$result = "Confused"


  # X-inactivation
  # http://bioinf.wehi.edu.au/software/GenderGenes/index.html
  Xgenes<- c("ARHGAP4","STS","ARSD", "ARSL", "AVPR2", "BRS3", "S100G", "CHM",
             "CLCN4", "DDX3X","EIF1AX","EIF2S3", "GPM6B", "GRPR", "HCFC1",
             "L1CAM", "MAOA", "MYCLP1", "NAP1L3", "GPR143", "CDK16", "PLXNB3",
             "PRKX", "RBBP7", "RENBP", "RPS4X", "TRAPPC2", "SH3BGRL", "TBL1X",
             "UBA1", "KDM6A", "XG", "XIST", "ZFX", "PUDP", "PNPLA4", "USP9X",
             "KDM5C", "SMC1A", "NAA10", "OFD1", "IKBKG", "PIR", "INE2", "INE1",
             "AP1S2", "GYG2", "MED14", "RAB9A", "ITM2A", "MORF4L2", "CA5B",
             "SRPX2", "GEMIN8", "CTPS2", "CLTRN", "NLGN4X", "DUSP21", "ALG13",
             "SYAP1", "SYTL4", "FUNDC1", "GAB3", "RIBC1", "FAM9C","CA5BP1")
  # genes belonging to the male-specific region of chromosome Y (unique genes)
  # http://bioinf.wehi.edu.au/software/GenderGenes/index.html
  Ygenes<-c("AMELY", "DAZ1", "PRKY", "RBMY1A1", "RBMY1HP", "RPS4Y1", "SRY",
            "TSPY1", "UTY", "ZFY","KDM5D", "USP9Y", "DDX3Y", "PRY", "XKRY",
            "BPY2", "VCY", "CDY1", "EIF1AY", "TMSB4Y","CDY2A", "NLGN4Y",
            "PCDH11Y", "HSFY1", "TGIF2LY", "TBL1Y", "RPS4Y2", "HSFY2",
            "CDY2B", "TXLNGY","CDY1B", "DAZ3", "DAZ2", "DAZ4")
#https://www.r-bloggers.com/2017/02/an-intro-to-gaussian-mixture-modeling/
# Load the package
library(mclust)
library(Seurat)
#library(Signac)
setwd(publication_path("LUNG_QC_ROOT", "data/2021-08-23.Qingmerge/"))
dat =readRDS( file = "Merged3runs_object.rds")

# add the gene activity matrix to the Seurat object as a new assay and normalize it

dat <- NormalizeData(
  object = dat,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(dat$nCount_RNA)
)





DefaultAssay(dat) <- 'RNA'
meta=dat[[]]
#top variable genes that separate sample in PC1 and PC2 by sex
#Used top 200 genes that contribute to PC1 (|loading|>0.05)
x=GetAssayData(dat)
# build artificial genes


Xgene.set <-rownames(x)[which(toupper(rownames(x)) %in%  Xgenes)]
Ygene.set <-rownames(x)[which(toupper(rownames(x)) %in%  Ygenes)]
cm.new<-as.data.frame(matrix(rep(0, 2*ncol(x)), ncol = ncol(x),nrow = 2))
    row.names(cm.new) <- c("superX","superY")
    colnames(cm.new) <- colnames(x)
    #cm.new["Xist", ]<- x["Xist", ]
    cm.new["superX", ] <-colMeans(x[Xgene.set,])
    cm.new["superY", ] <-colMeans(x[Ygene.set,])

tcm.final <- t(cm.new)
tcm.final <- as.data.frame(tcm.final)
zero.cells <- NA
dnc <- tcm.final$superY==0 & tcm.final$superX==0
if(any(dnc)==TRUE){
  zero.cells <- row.names(tcm.final)[dnc]
  cat(length(zero.cells), "cell/s are unable to be classified due to an
      abundance of zeroes on X and Y chromosome genes\n")
}
tcm.final <- tcm.final[!dnc, ]
# the processed transposed count matrix
tcm <-tcm.final
tcm$sex=factor(meta[!dnc, ]$sex, levels=c('M', 'F'))

LungTCM = tcm
## R program to implement
# validation set approach
library(caret)
# setting seed to generate a
# reproducible random sampling
set.seed(123)
# creating training data as 80% of the dataset

# generating training dataset
# from the random_sample
training_dataset <- tcm
LungModel <- glm(sex ~., data = training_dataset, family=binomial)

  

    meta=combined[[]]
    #top variable genes that separate sample in PC1 and PC2 by sex
    #Used top 200 genes that contribute to PC1 (|loading|>0.05)
    x=GetAssayData(combined)
    # build artificial genes
    Xgene.set <-rownames(x)[which(toupper(rownames(x)) %in%  Xgenes)]
    Ygene.set <-rownames(x)[which(toupper(rownames(x)) %in%  Ygenes)]
    cm.new<-as.data.frame(matrix(rep(0, 2*ncol(x)), ncol = ncol(x),nrow = 2))
    row.names(cm.new) <- c("superX","superY")
    colnames(cm.new) <- colnames(x)
    #cm.new["Xist", ]<- x["Xist", ]
    cm.new["superX", ] <-colMeans(x[Xgene.set,])
    cm.new["superY", ] <-colMeans(x[Ygene.set,])
    
    tcm.final <- t(cm.new)
    tcm.final <- as.data.frame(tcm.final)
    zero.cells <- NA
    dnc <- tcm.final$superY==0 & tcm.final$superX==0
    if(any(dnc)==TRUE){
      zero.cells <- row.names(tcm.final)[dnc]
      cat(length(zero.cells), "cell/s are unable to be classified due to an
          abundance of zeroes on X and Y chromosome genes\n")
    }
    tcm.final <- tcm.final[!dnc, ]
    # the processed transposed count matrix
    tcm <-tcm.final
  
    HumanHeartTCM = tcm
    predictions <- predict(LungModel, HumanHeartTCM,type = "response")
    # computing model performance metrics
    
    pred_ <- as.factor(ifelse(predict(LungModel, HumanHeartTCM, type="response")>0.5,"F","M"))
    combined$sexPredict_NoXist = pred_

 setwd(publication_path("HEART_SNRNA_ROOT", "data/2022-11-08.ATACObject"))   
    png('PredictionConfidencenoXist.png')
    hist(predictions)
    abline(v = 0.5, col = "red")
    dev.off()


plot6 <- DimPlot(
  object = combined,,
  split.by = 'sexPredict_NoXist',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
ggsave(file = "sex_prediction_noXist_umap.png", width = 20, height = 10)


VlnPlot(combined, "Xist", split.by = 'sexPredict_NoXist')
ggsave(file = "sex_prediction_NoXist_Xist.png", width = 20, height = 10)

VlnPlot(combined, "Uty", split.by = 'sexPredict_NoXist')
ggsave(file = "sex_prediction_NoXist_Uty.png", width = 20, height = 10)


saveRDS(combined,"combinedWithSex.rds")




 

  setwd(publication_path("HEART_SNRNA_ROOT", "data/2022-11-08.ATACObject"))   
 combined = readRDS("combinedWithSex.rds")
Atrial = 9
Vent = c(1,2,5)
Coronary_EC = c(14,10,4)
Endocardial = c(8)
Lympathic_EC = c(15)
Fibroblast= c(3,0,11,12)
Myeloid = c(6,20)
Bcell = 13
Tcell = 16
Neuronal = 19
Mesothelial = 17 
Pericyte = 7
smooth_muscle = 18
unknown2 = 21
unknown4 = 22
unknown5 = 23

combined[["cellType"]] <- combined[["seurat_clusters"]]
combined[["cellType"]]  <- ifelse(unlist(combined[["seurat_clusters"]]) %in% Atrial, "Atrial",unlist(combined[["cellType"]]))
combined[["cellType"]]  <- ifelse(unlist(combined[["seurat_clusters"]]) %in% Vent, "Vent",unlist(combined[["cellType"]]))
combined[["cellType"]]  <- ifelse(unlist(combined[["seurat_clusters"]]) %in% Coronary_EC, "Coronary_EC",unlist(combined[["cellType"]]))
combined[["cellType"]]  <- ifelse(unlist(combined[["seurat_clusters"]]) %in% Endocardial, "Endocardial",unlist(combined[["cellType"]]))
combined[["cellType"]]  <- ifelse(unlist(combined[["seurat_clusters"]]) %in% Lympathic_EC, "Lympathic_EC",unlist(combined[["cellType"]]))
combined[["cellType"]]  <- ifelse(unlist(combined[["seurat_clusters"]]) %in% Mesothelial, "Mesothelial",unlist(combined[["cellType"]]))
combined[["cellType"]]  <- ifelse(unlist(combined[["seurat_clusters"]]) %in% Fibroblast, "Fibroblast",unlist(combined[["cellType"]]))
combined[["cellType"]]  <- ifelse(unlist(combined[["seurat_clusters"]]) %in% Myeloid, "Myeloid",unlist(combined[["cellType"]]))
combined[["cellType"]]  <- ifelse(unlist(combined[["seurat_clusters"]]) %in% Tcell, "Tcell",unlist(combined[["cellType"]]))
combined[["cellType"]]  <- ifelse(unlist(combined[["seurat_clusters"]]) %in% Bcell, "Bcell",unlist(combined[["cellType"]]))
combined[["cellType"]]  <- ifelse(unlist(combined[["seurat_clusters"]]) %in% Pericyte, "Pericyte",unlist(combined[["cellType"]]))
combined[["cellType"]]  <- ifelse(unlist(combined[["seurat_clusters"]]) %in% smooth_muscle, "smooth_muscle",unlist(combined[["cellType"]]))
combined[["cellType"]]  <- ifelse(unlist(combined[["seurat_clusters"]]) %in% Neuronal, "Neuronal",unlist(combined[["cellType"]]))


combined[["cellType"]]  <- ifelse(unlist(combined[["seurat_clusters"]]) %in% unknown2, "unknown2",unlist(combined[["cellType"]]))
combined[["cellType"]]  <- ifelse(unlist(combined[["seurat_clusters"]]) %in% unknown4, "unknown4",unlist(combined[["cellType"]]))
combined[["cellType"]]  <- ifelse(unlist(combined[["seurat_clusters"]]) %in% unknown5, "unknown5",unlist(combined[["cellType"]]))


#assign metadata
combined[["dataset"]][combined[["dataset"]] == "MAS3_2_v2"] = "MAS3_2"
combined$group = table(combined$dataset)
combined$group = (substr(combined$dataset,1,nchar(combined$dataset)-3))
table(combined$group)
combined$age = (substr(combined$group,1,nchar(combined$group)-1))
table(combined$age)

combined$conditon = (substr(combined$group,nchar(combined$group),nchar(combined$group)))
table(combined$conditon)
 combined$sex = combined$sexPredict_NoXist
combined$ind = paste(combined$dataset,combined$sex, sep = "_")


saveRDS(combined,"MergedAtac_QC_annotated.rds")
  setwd(publication_path("HEART_SNRNA_ROOT", "data/2022-11-08.ATACObject"))   
 combined = readRDS("MergedAtac_QC_annotated.rds")
