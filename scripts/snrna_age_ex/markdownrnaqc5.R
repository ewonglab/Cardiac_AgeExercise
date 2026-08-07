# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Perform primary snRNA-seq QC, integration, annotation and fixed-effects discovery analysis.
# Steps: 1) Load corrected matrices; 2) QC and integrate nuclei; 3) annotate cell types and export downstream analysis objects/results.
# =======================================
# Source reformatted from snRNA_AgeEx/MarkdownRNAQC5.Rmd.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

knitr::opts_chunk$set(echo = TRUE)
library(DropletUtils)
library(Seurat)
library(ggplot2)
library(scran)
library(dplyr)
library(knitr)
library(stringr)
knitr::opts_chunk$set(error = TRUE)



Runs = c("MAE2_1","MAS1_4","MAS3_2","YE2_1","YE4_3","YS1_4","YS3_2")


Seurat_objs = list()
#loop through for each run
for (s in 1:length(Runs)){
    Run = Runs[s]
    fname <- file.path(paste0(publication_path("HEART_SNRNA_ROOT", "data/2022-01-23.Cellbender/2ndRun/"),Run,"_cellbender_filtered.h5"))
    
    expression_matrix = Read10X_h5(fname, use.names = TRUE, unique.features = TRUE)
    
    obj = CreateSeuratObject(counts = expression_matrix)
    
    library(tibble)
    
    #annotate Mitochondria genes
    mouse_mito = as_tibble(read.delim(publication_path("MM10_RESOURCES_DIR", "Mouse.MitoCarta3.0.txt"), header = TRUE))
    mouse_mito = mouse_mito %>% dplyr::select(c(Symbol, MitoCarta2.0_FDR)) %>% dplyr::slice(1:100)
    mito.genes = as.character(mouse_mito$Symbol)
    mito.genes = mito.genes[mito.genes %in% rownames(obj[['RNA']])]
    percent.mito = Matrix::colSums(obj[['RNA']][mito.genes,]) / Matrix::colSums(obj[['RNA']])  * 100
    obj[["percent.mt"]] <- percent.mito
    print(summary(obj@meta.data$percent.mt))
    print(obj) #27665
    Seurat_objs[[s]] = obj
}
setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-03-01.CellbenderFiltered_MAS3_2_Adjust"))
#saveRDS(Seurat_objs, file="Seurat_Cellbender_objs.rds")


#filter all
library(Seurat)
library(ggplot2)
library(scran)

#Define Runs and loop through them
Runs = c("MAE2_1","MAS1_4","MAS3_2","YE2_1","YE4_3","YS1_4","YS3_2")
setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-03-01.CellbenderFiltered_MAS3_2_Adjust"))
#Seurat_objs = readRDS(file="Seurat_Cellbender_objs.rds")
#loop through each run
for (i in 1:length(Runs)){
    Run = Runs[i]
  print(Run)
  
  
  #read in emptyDrops

  nonEmpty = Seurat_objs[[i]]
  nonEmpty
  #read in doublets
  setwd(publication_path("HEART_SNRNA_ROOT", "data/2022-07-04.Scrublet"))
  doublets = read.table(file=paste0(Run,"_doublet_scores_doublet.txt"))
  head(doublets)
  doublet_cells = doublets$cell
  
  #subsetting
  
  
  #Remove Doublets From Empty Drops Object
  doublet_Noempty <- nonEmpty[,!colnames(nonEmpty) %in% doublet_cells]
  
  print("nonDoublet_nonEmpty cells")
  print(doublet_Noempty)
  
  setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-03-01.CellbenderFiltered_MAS3_2_Adjust"))
  #run QC
  Idents(doublet_Noempty) = doublet_Noempty$orig.ident
  doublet_Noempty[["Seurat_percent.mt"]] <- PercentageFeatureSet(doublet_Noempty, pattern = "^mt-")
  print(paste(Run, "Quality Control"))
  print(VlnPlot(doublet_Noempty, features = c("nFeature_RNA", "nCount_RNA", "Seurat_percent.mt"),ncol = 3, pt.size = 0.05))
  ggsave( file = paste0(Run,"vlnPlot.png"))
  doublet_Noempty$run = Run
  plot1 <- FeatureScatter(doublet_Noempty, feature1 = "nCount_RNA", feature2 = "Seurat_percent.mt")
  plot2 <- FeatureScatter(doublet_Noempty, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(paste(Run, "Count vs Features"))
  print(plot1 + plot2)
  ggsave(plot1 + plot2, file = paste0(Run,"featuresVsCountsPlot.pdf"))
  
  #Filter Using Seurat method 
  doublet_Noempty <- subset(doublet_Noempty, subset = nFeature_RNA > 500 & nFeature_RNA < 5500 & Seurat_percent.mt < 5) 
  print(paste(Run, "Quality Control Post filtering"))
  print(VlnPlot(doublet_Noempty, features = c("nFeature_RNA", "nCount_RNA", "Seurat_percent.mt"), ncol = 3, pt.size = 0.05))
  ggsave(VlnPlot(doublet_Noempty, features = c("nFeature_RNA", "nCount_RNA", "Seurat_percent.mt"), ncol = 3, pt.size = 0.05), file = paste0(Run,"PostFilter_vlnPlot.png"))
  
  #doublet_Noempty =NormalizeData(doublet_Noempty)
  setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-03-01.CellbenderFiltered_MAS3_2_Adjust"))
  #saveRDS(doublet_Noempty, file = paste0(Run,"Normalized_hashed_doublet_empty_nonemptyCells.rds"))
  
  print(paste(Run,"Filtered_Hashed_nonDoublet_nonEmpty cells"))
  print(doublet_Noempty)
  
  print("\n\n\n\n")
  
}


library(Seurat)
library(ggplot2)

 #Merge good runs together 
setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-03-01.CellbenderFiltered_MAS3_2_Adjust"))
Runs = c("MAE2_1","MAS1_4","MAS3_2","YE2_1","YE4_3","YS1_4","YS3_2")
MAE2_1 = readRDS(paste0(Runs[1],"Normalized_hashed_doublet_empty_nonemptyCells.rds"))
MAS1_4  = readRDS(paste0(Runs[2],"Normalized_hashed_doublet_empty_nonemptyCells.rds"))
MAS3_2  = readRDS(paste0(Runs[3],"Normalized_hashed_doublet_empty_nonemptyCells.rds"))
YE2_1  = readRDS(paste0(Runs[4],"Normalized_hashed_doublet_empty_nonemptyCells.rds"))
YE4_3  = readRDS(paste0(Runs[5],"Normalized_hashed_doublet_empty_nonemptyCells.rds"))
YS1_4  = readRDS(paste0(Runs[6],"Normalized_hashed_doublet_empty_nonemptyCells.rds"))
YS3_2  = readRDS(paste0(Runs[7],"Normalized_hashed_doublet_empty_nonemptyCells.rds"))
MAE2_1 

MAS1_4  
MAS3_2  
YE2_1 
YE4_3  
YS1_4  
YS3_2

#NOTE: run MAE4_3 excluded for low cell count. 
setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-03-01.CellbenderFiltered_MAS3_2_Adjust"))
run.merge <- merge(MAE2_1, y = c(MAS1_4, MAS3_2,YE2_1,YE4_3,YS1_4,YS3_2), add.cell.ids = Runs, project = "All_Heart_runs")
run.merge
run.merge =NormalizeData(run.merge)


setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-03-01.CellbenderFiltered_MAS3_2_Adjust"))

#Find variable features
run.merge <- FindVariableFeatures(run.merge, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(run.merge), 10)
plot1 <- VariableFeaturePlot(run.merge)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(paste("Variable Features"))
print(plot2)
ggsave(plot2,file = "run.merge_topFeatures.png")

all.genes <- rownames(run.merge)
run.merge <- ScaleData(run.merge, features = all.genes)    
print("scaleDone")

#Run PCA
run.merge <- RunPCA(run.merge, features = VariableFeatures(object = run.merge))
print(paste("All merged PCA"))
print(DimPlot(run.merge, reduction = "pca"))
ggsave(DimPlot(run.merge, reduction = "pca"), file = "run.merge_pca.png")
print(paste("All merged PCA group by run"))
print(DimPlot(run.merge, reduction = "pca", group.by = "run"))
ggsave(DimPlot(run.merge, reduction = "pca", group.by = "run"), file = "run.merge_run_pca.png")


#Check dimensionality
run.merge <- JackStraw(run.merge, num.replicate = 100)
run.merge <- ScoreJackStraw(run.merge, dims = 1:20)


print(JackStrawPlot(run.merge, dims = 1:20))
ggsave(JackStrawPlot(run.merge, dims = 1:20), file = "run.merge_jackstraw.png")


print(ElbowPlot(run.merge))
ggsave(ElbowPlot(run.merge),file = "run.merge_Elbowplot.png")


#Run First UMAP
run.merge <- FindNeighbors(run.merge, dims = 1:20) #note dimensions 
run.merge <- FindClusters(run.merge, resolution = 0.5)

run.merge <- RunUMAP(run.merge, dims = 1:20)

print(paste("All Runs UMAP"))
print(DimPlot(run.merge, reduction = "umap"))

ggsave(DimPlot(run.merge, reduction = "umap"),file = "run.merge_umap.png")

Idents(run.merge) <- run.merge$run
tplot <- DimPlot(run.merge, reduction = "umap")
tplot[[1]]$layers[[1]]$aes_params$alpha = .2
print(paste("All Runs colour by run"))
tplot
ggsave(tplot,file = "run.merge_run_umap.png")
Idents(run.merge) <- run.merge$seurat_clusters
print(paste("QC plots"))
print(VlnPlot(run.merge, features = c("nFeature_RNA", "nCount_RNA", "Seurat_percent.mt")))
 ggsave(VlnPlot(run.merge, features = c("nFeature_RNA", "nCount_RNA", "Seurat_percent.mt"), ncol = 3, pt.size = 0.05, split.by = "run", group.by = "seurat_clusters"), file = paste0("PostFilter_vlnPlot.png"), width = 30)
 
 #Check contamination gene quality
 print(paste("Contamination genes"))
 print(VlnPlot(run.merge, features = c("Ryr2","Dmd"), ncol = 2, pt.size = 0.05, split.by = "run"))
       
 print(paste("Contamination genes Per Cluster"))
print(VlnPlot(run.merge, features = c("Ryr2","Dmd"), ncol = 2, pt.size = 0.05, split.by = "run", group.by = "seurat_clusters"))
  ggsave(VlnPlot(run.merge, features = c("Ryr2","Dmd"), ncol = 3, pt.size = 0.05, split.by = "run", group.by = "seurat_clusters"), file = paste0("Contamination_vlnPlot.png"), width = 30)
saveRDS(run.merge, file = "Mergedruns_object.rds")



setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-03-01.CellbenderFiltered_MAS3_2_Adjust"))
Mouse = run.merge
#renormalize using defualt method
Mouse <- NormalizeData(
  object = Mouse,
  assay = 'RNA'
)


#Check Different number of dimensions effect

Mouse <- FindNeighbors(Mouse, dims = 1:50) #note dimensions 
Mouse <- FindClusters(Mouse, resolution = 0.5)

Mouse <- RunUMAP(Mouse, dims = 1:50)
print(paste("50 Dim UMAP"))
print(DimPlot(Mouse, reduction = "umap"))
ggsave(DimPlot(Mouse, reduction = "umap"),file = "Mouse_50_umap.png")


Mouse <- FindNeighbors(Mouse, dims = 1:15) #note dimensions 
Mouse <- FindClusters(Mouse, resolution = 0.5)

Mouse <- RunUMAP(Mouse, dims = 1:15)
print(paste("15 Dim UMAP"))
print(DimPlot(Mouse, reduction = "umap",group.by = "run"))
ggsave(DimPlot(Mouse, reduction = "umap",group.by = "run"),file = "Mouse_15_umap.png")


Mouse <- FindNeighbors(Mouse, dims = 1:5) #note dimensions 
Mouse <- FindClusters(Mouse, resolution = 0.5)

Mouse <- RunUMAP(Mouse, dims = 1:5)
print(paste("5 Dim UMAP"))
print(DimPlot(Mouse, reduction = "umap", group.by = "run"))
ggsave(DimPlot(Mouse, reduction = "umap", group.by = "run"),file = "Mouse_5_umap.png")


Mouse <- FindNeighbors(Mouse, dims = 1:17) #note dimensions 
Mouse <- FindClusters(Mouse, resolution = 0.5)

Mouse <- RunUMAP(Mouse, dims = 1:17)
print(paste("17 Dim UMAP"))
print(DimPlot(Mouse, reduction = "umap",group.by = "run"))
ggsave(DimPlot(Mouse, reduction = "umap",group.by = "run"),file = "Mouse_17_umap.png")



#15 works best 




#check resolutions 


Mouse <- FindNeighbors(Mouse, dims = 1:15) #note dimensions 
Mouse <- FindClusters(Mouse, resolution = 0.5)
Mouse <- RunUMAP(Mouse, dims = 1:15)
Idents(Mouse) = Mouse$seurat_clusters
print(paste("0.5 Res UMAP"))
print(DimPlot(Mouse, reduction = "umap",group.by = "seurat_clusters", label = T))
ggsave(DimPlot(Mouse, reduction = "umap",group.by = "seurat_clusters", label = T),file = "Mouse_15_0.5_umap.png")


Mouse <- FindClusters(Mouse, resolution = 0.2)
Mouse <- RunUMAP(Mouse, dims = 1:15)
Idents(Mouse) = Mouse$seurat_clusters
print(paste("0.2 Res UMAP"))
print(DimPlot(Mouse, reduction = "umap",group.by = "seurat_clusters", label = T))
ggsave(DimPlot(Mouse, reduction = "umap",group.by = "seurat_clusters", label = T),file = "Mouse_15_0.2_umap.png")


Mouse <- FindClusters(Mouse, resolution = 0.8)
Mouse <- RunUMAP(Mouse, dims = 1:15)
Idents(Mouse) = Mouse$seurat_clusters
print(paste("0.8 Res UMAP"))
print(DimPlot(Mouse, reduction = "umap",group.by = "seurat_clusters", label = T))
ggsave(DimPlot(Mouse, reduction = "umap",group.by = "seurat_clusters", label = T),file = "Mouse_15_0.8_umap.png")


Mouse <- FindClusters(Mouse, resolution = 0.3)
Mouse <- RunUMAP(Mouse, dims = 1:15)
Idents(Mouse) = Mouse$seurat_clusters
print(paste("0.3 Res UMAP"))
print(DimPlot(Mouse, reduction = "umap",group.by = "seurat_clusters", label = T))
ggsave(DimPlot(Mouse, reduction = "umap",group.by = "seurat_clusters", label = T),file = "Mouse_15_0.3_umap.png")

Mouse <- FindClusters(Mouse, resolution = 0.4)
Mouse <- RunUMAP(Mouse, dims = 1:15)
Idents(Mouse) = Mouse$seurat_clusters
print(paste("0.4 Res UMAP"))
print(DimPlot(Mouse, reduction = "umap",group.by = "seurat_clusters", label = T))
ggsave(DimPlot(Mouse, reduction = "umap",group.by = "seurat_clusters", label = T),file = "Mouse_15_0.4_umap.png")

#0.5 resolution best for cellTypes 




Mouse <- FindNeighbors(Mouse, dims = 1:15) #note dimensions 
Mouse <- FindClusters(Mouse, resolution = 0.5)
Mouse <- RunUMAP(Mouse, dims = 1:15)
Idents(Mouse) = Mouse$seurat_clusters
ggsave(DimPlot(Mouse, reduction = "umap",group.by = "seurat_clusters", label = T),file = "Mouse_clusters_umap.png")


# genes located in the X chromosome that have been reported to escape
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
  scale.factor = 10000
)


DefaultAssay(dat) <- 'RNA'
meta=dat[[]]
#top variable genes that separate sample in PC1 and PC2 by sex
#Used top 200 genes that contribute to PC1 (|loading|>0.05)
x=GetAssayData(dat)
# build artificial genes
Xgene.set <-rownames(x)[which(toupper(rownames(x)) %in%  Xgenes)]
Ygene.set <-rownames(x)[which(toupper(rownames(x)) %in%  Ygenes)]
cm.new<-as.data.frame(matrix(rep(0, 3*ncol(x)), ncol = ncol(x),nrow = 3))
row.names(cm.new) <- c('Xist',"superX","superY")
colnames(cm.new) <- colnames(x)
cm.new["Xist", ]<- x["Xist", ]
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
random_sample <- createDataPartition(tcm$sex, p = 0.8, list = FALSE)
# generating training dataset
# from the random_sample
training_dataset <- tcm[random_sample, ]
# generating testing dataset
# from rows which are not
# included in random_sample
testing_dataset <- tcm[-random_sample, ]
# Building the model
# training the model by assigning sales column
# as target variable and rest other columns
# as independent variables
#https://stackoverflow.com/questions/62437321/set-cutoff-threshold-when-predicting-in-r
LungModel <- glm(sex ~., data = training_dataset, family=binomial)
# predicting the target variable
predictions <- predict(LungModel, testing_dataset)
# computing model performance metrics
pred_ <- as.factor(ifelse(predict(LungModel, testing_dataset, type="response")>0.7,"F","M"))
# here we go!
confusionMatrix(pred_, as.factor(testing_dataset$sex))
#https://daviddalpiaz.github.io/r4sl/logistic-regression.html
library(pROC)
test_prob = predict(LungModel, newdata = testing_dataset, type = "response") #predict probability
pred_ <- as.factor(ifelse(predict(LungModel, testing_dataset, type="response")>0.5,"F","M"))
# here we go!
confusionMatrix(pred_, as.factor(testing_dataset$sex))
test_roc = roc(testing_dataset$sex ~ test_prob, plot = TRUE, print.auc = TRUE)
as.numeric(test_roc$auc)
setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2022-08-16.sex"))
pdf('gene_Lung_auc.pdf')
test_roc = roc(testing_dataset$sex ~ test_prob, plot = TRUE, print.auc = TRUE)
dev.off()

    
   
#Mouse = run.merge
setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-03-01.CellbenderFiltered_MAS3_2_Adjust"))


Mouse <- NormalizeData(
  object = Mouse,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = 10000
)


DefaultAssay(Mouse) <- 'RNA'
meta=Mouse[[]]
#top variable genes that separate sample in PC1 and PC2 by sex
#Used top 200 genes that contribute to PC1 (|loading|>0.05)
x=GetAssayData(Mouse)
# build artificial genes
Xgene.set <-rownames(x)[which(toupper(rownames(x)) %in%  Xgenes)]
Ygene.set <-rownames(x)[which(toupper(rownames(x)) %in%  Ygenes)]
cm.new<-as.data.frame(matrix(rep(0, 3*ncol(x)), ncol = ncol(x),nrow = 3))
row.names(cm.new) <- c('Xist',"superX","superY")
colnames(cm.new) <- colnames(x)
cm.new["Xist", ]<- x["Xist", ]
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

HeartTCM = tcm


predictions <- predict(LungModel, HeartTCM)
pred_ <- as.factor(ifelse(predict(LungModel, HeartTCM, type="response")>0.5,"F","M"))
table(pred_, useNA = "always")
#pred_
#    F     M
#61127 71428
returnPredict = predict(LungModel, tcm, type="response")
png('PredictionConfidence.png')
hist(returnPredict)
abline(v = 0.5, col = "red")
dev.off()

Mouse$sexPredict = pred_

Mouse$sexPredictValue = returnPredict


plot1 <- DimPlot(
  object = Mouse,
  group.by = 'sexPredict', #seurat_clusters', #"celltype",
  label = TRUE,
  repel = TRUE) + ggtitle('Predicted sex')
print(paste("UMAP with predicted sex"))
plot1
ggsave("sexLabels.png")

plot2 <- DimPlot(
  object = Mouse,
  group.by = 'sexPredict',
  split.by = 'sexPredict',#seurat_clusters', #"celltype",
  label = TRUE,
  repel = TRUE) + ggtitle('Predicted sex')
print(paste("UMAP with predicted sex split"))
plot2
ggsave("sexSplitLabels.png")

setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-03-01.CellbenderFiltered_MAS3_2_Adjust"))



Mouse$ID = paste0(Mouse$run,"_",Mouse$sexPredict)

setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-03-01.CellbenderFiltered_MAS3_2_Adjust"))
plot1 <- DimPlot(
  object = Mouse,
  group.by = 'ID', #seurat_clusters', #"celltype",
  ) + ggtitle('snRNA-seq')
print(paste("UMAP with ID"))
plot1
ggsave("IndividualMapping.umap.png")
setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-03-01.CellbenderFiltered_MAS3_2_Adjust"))
plot1 <- DimPlot(
  object = Mouse,
  group.by = 'ID',
  split.by = 'ID',#seurat_clusters', #"celltype",
  ) + ggtitle('snRNA-seq')
print(paste("UMAP split by ID"))
plot1
ggsave("IndividualMappingSplit.umap.png", height =10, width = 20)


#saveRDS(Mouse,"MouseWithPredictedSex.rds")
setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-03-01.CellbenderFiltered_MAS3_2_Adjust"))
#Mouse = readRDS("MouseWithPredictedSex.rds")




setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-03-01.CellbenderFiltered_MAS3_2_Adjust"))
#Mouse = readRDS("MouseWithPredictedSex.rds")
#renormalize using defualt method
Mouse <- NormalizeData(
  object = Mouse,
  assay = 'RNA'
)

#Assign metadata from run info

Mouse$age = ifelse(nchar(Mouse$run) == 6,"MA","Y")
#substring(Mouse$run, 2, 3)
Mouse$cond = ifelse(substring(Mouse$run, 2, 3) %in% c("AS", "S1","S3"),"S","E")

Mouse$group = paste0(Mouse$age,Mouse$cond)
setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-03-01.CellbenderFiltered_MAS3_2_Adjust"))

Mouse$cellID = colnames(Mouse)


#Generate all Feautre plots of cell markers
print(paste("Neurons Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Nrxn1", "Nrxn3","Kcnmb4"),  ncol =2, cols = c("grey", "red"))
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Nrxn1", "Nrxn3","Kcnmb4"),  ncol =2, cols = c("grey", "red")), file = 'Neurons2_marker_feature_plot_cellbender.png')

print(paste("myofibroblast Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Acta2", "Des","Pecam1"),  ncol =2, cols = c("grey", "red"))
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Acta2", "Des","Pecam1"),  ncol =2, cols = c("grey", "red")), file = 'myofibroblast_marker_feature_plot_cellbender.png')

print(paste("EC Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Kdr", "Pecam1","Ly6a", "Fabp4"),  ncol =2, cols = c("grey", "red"))
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Kdr", "Pecam1","Ly6a", "Fabp4"),  ncol =2, cols = c("grey", "red")), file = 'EC_marker_feature_plot_cellbender.png')

print(paste("Fibroblast Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Col1a1", "Pdgfra","Ddr2", "Tcf21"),  ncol =2, cols = c("grey", "red"))
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Col1a1", "Pdgfra","Ddr2", "Tcf21"),  ncol =2, cols = c("grey", "red")), file = 'Fibroblast_marker_feature_plot_cellbender.png')

print(paste("CM Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Tnnt2", "Ttn","Myh6", "Tnni3"),  ncol =2, cols = c("grey", "red"))
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Tnnt2", "Ttn","Myh6", "Tnni3"),  ncol =2, cols = c("grey", "red")), file= 'CM_marker_feature_plot_cellbender.png')

print(paste("Mural Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Pdgfrb", "Vtn","Acta2", "Rgs5"),  ncol =2, cols = c("grey", "red"))
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Pdgfrb", "Vtn","Acta2", "Rgs5"),  ncol =2, cols = c("grey", "red")), file = 'Mural_marker_feature_plot_cellbender.png')

print(paste("Ribo Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Rpl7a", "Rpl26","Rps3","Rps6"),  ncol =2, cols = c("grey", "red"))
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Rpl7a", "Rpl26","Rps3","Rps6"),  ncol =2, cols = c("grey", "red")), file = 'ribo_marker_feature_plot_cellbender.png')

print(paste("leukocytes Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Ptprc"),  ncol =2, cols = c("grey", "red"))
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Ptprc"),  ncol =2, cols = c("grey", "red")), file = 'leukocytes_marker_feature_plot_cellbender.png')

print(paste("T Cell Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Cd3e","Cd4","Cd8a"  ),  ncol =2, cols = c("grey", "red"))
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Cd3e","Cd4","Cd8a"  ),  ncol =2, cols = c("grey", "red")), file = 'TCELL_marker_feature_plot_cellbender.png')

print(paste("B Cell Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Cd19","Cd20","Ms4a1"  ),  ncol =2, cols = c("grey", "red"))
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Cd19","Cd20","Ms4a1"  ),  ncol =2, cols = c("grey", "red")), file = 'BCELL_marker_feature_plot_cellbender.png')

print(paste("Macrophages Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Itgam","Itgax","Ly6c1"  ),  ncol =2, cols = c("grey", "red"))
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Itgam","Itgax","Ly6c1"  ),  ncol =2, cols = c("grey", "red")), file = 'macrophages_marker_feature_plot_cellbender.png')

print(paste("neutrophils Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Ly6g"),  ncol =2, cols = c("grey", "red"))
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Ly6g"),  ncol =2, cols = c("grey", "red")), file = 'neutrophils_marker_feature_plot_cellbender.png')

print(paste("Harvey Atrial_Vent Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Hand2","Nr2f1","Hand1","Hey2"),  ncol =2, cols = c("grey", "red"))
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Hand2","Nr2f1","Hand1","Hey2"),  ncol =2, cols = c("grey", "red")), file = 'Harvey_Atrial_Vent_marker_feature_plot_cellbender.png')

print(paste("macrophages2 Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Cd86","Cd80","Cd163","Mhcii","Fizz1","Tlr1","Tlr8","Vegf","Ccl2","Fyn","Cd11b"),  ncol =2, cols = c("grey", "red"))
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Cd86","Cd80","Cd163","Mhcii","Fizz1","Tlr1","Tlr8","Vegf","Ccl2","Fyn","Cd11b"),  ncol =2, cols = c("grey", "red")), file = 'macrophages2_marker_feature_plot_cellbender.png')

print(paste("macrophages_4 Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c('F13a1', 'Rbpj', 'Lilrb4a', 'Slc9a9', 'Mrc1', 'Cd163'),  ncol =2, cols = c("grey", "red"))
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c('F13a1', 'Rbpj', 'Lilrb4a', 'Slc9a9', 'Mrc1', 'Cd163'),  ncol =2, cols = c("grey", "red")), file = 'macrophages_4_marker_feature_plot_cellbender.png')

print(paste("Adipocytes Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Gpam","Fasn","Lep","Hoxc9"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE)
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Gpam","Fasn","Lep","Hoxc9","Mpzl2","Ebf3","Fbox31","Tcf21"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'Adipocytes_marker_feature_plot_cellbender.png')

print(paste("AtrialCM Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Nppa","Myl7","Myl4"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE)
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Nppa","Myl7","Myl4"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'AtrialCM_marker_feature_plot_cellbender.png')

print(paste("Endothilial Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Vwf","Pecam1","Cdh5"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE)
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Vwf","Pecam1","Cdh5"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'Endothilial_marker_feature_plot_cellbender.png')


print(paste("lymphoid Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Cd8a","Il7r","Cd40lg"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE)
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Cd8a","Il7r","Cd40lg"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'lymphoid_marker_feature_plot_cellbender.png')


print(paste("mesothelial Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Msln","Wt1","Bnc1"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE)
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Msln","Wt1","Bnc1"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'mesothelial_marker_feature_plot_cellbender.png')


print(paste("myeloid Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Cd14","C1qa","Cd68"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE)
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Cd14","C1qa","Cd68"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'myeloid_marker_feature_plot_cellbender.png')


print(paste("Neuronal Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Plp1","Nrxn1","Nrxn3"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE)
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Plp1","Nrxn1","Nrxn3"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'Neuronal_marker_feature_plot_cellbender.png')


print(paste("Pericytes Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Rgs5","Abcc9","Kcnj8"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE)
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Rgs5","Abcc9","Kcnj8"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'Pericytes_marker_feature_plot_cellbender.png')


print(paste("SmoothMuscle Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Myh11","Tagln","Acta2"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE)
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Myh11","Tagln","Acta2"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'SmoothMuscle_marker_feature_plot_cellbender.png')


print(paste("VentricalCM Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Myh7","Myl2","Fhl2"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE)
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Myh7","Myl2","Fhl2"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'VentricalCM_marker_feature_plot_cellbender.png')

print(paste("hepatocyte Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Alb","Cyp2e1"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE)
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Alb","Cyp2e1"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'hepatocyte_marker_feature_plot_cellbender.png')


print(paste("proliferating Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Top2a","Tubb5","Prc1","Cdk1","Stmn1","Cenpf"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE)
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Top2a","Tubb5","Prc1","Cdk1","Stmn1","Cenpf"),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'proliferating_marker_feature_plot_cellbender.png')


print(paste("macrophage Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c('C1qb','C1qa','C1qc','Lyz2','Ctsc','Ctss'),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE)
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c('C1qb','C1qa','C1qc','Lyz2','Ctsc','Ctss'),  ncol =2, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'macrophage_marker_feature_plot_cellbender.png')


#https://www.sciencedirect.com/science/article/pii/S1534580716307523?via%3Dihub

print(paste("AtrialScienceCM_marker Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Cav1","Nr2f2","Thbd","Stard10","Vsnl1","Nr2f1","Gm15543","Myh6"),  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE)
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Cav1","Nr2f2","Thbd","Stard10","Vsnl1","Nr2f1","Gm15543","Myh6"),  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'AtrialScienceCM_marker_feature_plot_cellbender.png')

print(paste("VentricularScienceCM Markers"))
FeaturePlot(Mouse, reduction = "umap", features = c("Myl2","Pln","Mpped2","Sall4","Capn6","Hey2","Myh7"),  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE)
ggsave(FeaturePlot(Mouse, reduction = "umap", features = c("Myl2","Pln","Mpped2","Sall4","Capn6","Hey2","Myh7"),  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'VentricularScienceCM_marker_feature_plot_cellbender.png')

Granulocyte  = str_to_title(c('CSF2RB', 'CSF3R', 'IL1R2', 'IL1RN', 'IL8RB', 'IL13RA1', 'FCGRT', 'FCGR3B','CD16', 'C3AR1','FPR1', 'CD10','MME',  'EMR2', 'TM7'))
neutrophils  = str_to_title(c('CD15', 'CD11b', 'CD16',  'CD10','CD16', 'CXCR2', 'CXCR4',  'CD62L'))
Myeloids = str_to_title(c("CD11b",'CD14','EMR1','ITGAM','Ly6G','CD115','CD68','CD206','MRC1'))

print(paste("Myeloids Markers"))
FeaturePlot(Mouse, reduction = "umap", features =Myeloids,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE)
ggsave(FeaturePlot(Mouse, reduction = "umap", features =Myeloids,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'Myeloids_marker_feature_plot_cellbender.png', height = 15, width = 15)

print(paste("neutrophils  Markers"))
FeaturePlot(Mouse, reduction = "umap", features =neutrophils ,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE)
ggsave(FeaturePlot(Mouse, reduction = "umap", features =neutrophils,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'neutrophils _marker_feature_plot_cellbender.png', height = 15, width = 15)


print(paste("Granulocyte Markers"))
FeaturePlot(Mouse, reduction = "umap", features =Granulocyte,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE)
ggsave(FeaturePlot(Mouse, reduction = "umap", features =Granulocyte,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'Granulocyte_marker_feature_plot_cellbender.png', height = 15, width = 15)


#gene lists from Matt


Artery = c("Stmn2", "Fos", "Id1", "Sox17", "Klf4", "Alpl", "Junb", "Crip1", "Mecom","Hspa1a")
Capillary_arterial = c("Cxcl12", "Rbp7", "Mgll", "Ly6c1", "Aqp7", "Btnl9", "Pdgfd","Sema7a", "Clec2d","Unc5b")
Capillary = c("Kdr", "Ssh2", "AU020206", "Fam57b", "Kifc3", "Ntf3", "Tmem182", "Ebf3","Coro2b", "Fam117b")
Capillary_venous = c("Ier3", "Pltp", "Fmo2", "Fmo1", "Calcrl", "Emp1", "AU021092", "Bsg","Nfkbia", "Eln")
Lymphatic = c("Ccl21a", "Mmrn1", "Fth1", "Lyve1", "Prelp", "Fxyd6", "S100a6", "Cd9","Pard6g", "Flt4")
endocardial_EC = c("Kank1", "Gpm6a", "Smoc1", "Cdh11", "Plagl1", "Tmem108", "Npr3", "Comp","Ahsg", "Mgp")
coronary_EC = c("Sparcl1", "Mgll", "Ly6c1", "Tspan13", "Aqp1", "Cldn5", "Timp4")
aorta_EC = c("Ehd3", "Fam167b", "Esm1", "AU021092", "Ppap2a")

print(paste("Artery_Endo Markers"))
FeaturePlot(Mouse, reduction = "umap", features =Artery,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE)
ggsave(FeaturePlot(Mouse, reduction = "umap", features =Artery,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'Artery_Endo_marker_feature_plot_cellbender.png', height = 15, width = 15)

print(paste("Capillary_arterial Markers"))
FeaturePlot(Mouse, reduction = "umap", features =Capillary_arterial,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE)
ggsave(FeaturePlot(Mouse, reduction = "umap", features =Capillary_arterial,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'Capillary_arterial_Endo_marker_feature_plot_cellbender.png', height = 15, width = 15)

print(paste("Capillary_Endo Markers"))
FeaturePlot(Mouse, reduction = "umap", features =Capillary,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE)
ggsave(FeaturePlot(Mouse, reduction = "umap", features =Capillary,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'Capillary_Endo_marker_feature_plot_cellbender.png', height = 15, width = 15)

print(paste("Capillary_venous_Endo Markers"))
FeaturePlot(Mouse, reduction = "umap", features =Capillary_venous,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE)
ggsave(FeaturePlot(Mouse, reduction = "umap", features =Capillary_venous,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'Capillary_venous_Endo_marker_feature_plot_cellbender.png', height = 15, width = 15)

print(paste("Lymphatic Markers"))
FeaturePlot(Mouse, reduction = "umap", features =Lymphatic,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE)
ggsave(FeaturePlot(Mouse, reduction = "umap", features =Lymphatic,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'Lymphatic_Endo_marker_feature_plot_cellbender.png', height = 15, width = 15)

print(paste("endocardial_EC Markers"))
FeaturePlot(Mouse, reduction = "umap", features =endocardial_EC,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE)
ggsave(FeaturePlot(Mouse, reduction = "umap", features =endocardial_EC,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'endocardial_EC_Endo_marker_feature_plot_cellbender.png', height = 15, width = 15)

print(paste("coronary_EC Markers"))
FeaturePlot(Mouse, reduction = "umap", features =coronary_EC,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE)
ggsave(FeaturePlot(Mouse, reduction = "umap", features =coronary_EC,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'coronary_EC_Endo_marker_feature_plot_cellbender.png', height = 15, width = 15)

print(paste("aorta_EC Markers"))
FeaturePlot(Mouse, reduction = "umap", features =aorta_EC,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE)
ggsave(FeaturePlot(Mouse, reduction = "umap", features =aorta_EC,  ncol =3, cols = c("grey", "red"), order = TRUE, raster = FALSE), file = 'aorta_EC_Endo_marker_feature_plot_cellbender.png', height = 15, width = 15)



#add Mouse celltype labels to seurat clusters
B_cell = c(13)
Cardiac_Neuronal = c(18)
Cardiomyocyte = c(1,5,9,10)
Endocardial_EC = c(7)
Fibroblast = c(0,2)
Lymphatic_EC = c(11)
Mesothelial = c(15)
Myeloid =c(4,17)
Pericytes = c(6,19)
Smooth_Muscle = c(16)
Unknown_hepatocyte = 20
Coronary_EC = c(3,8)
T_cell = c(12)
Proliferating_cell = c(21)
Macrophage = c(14)
Unknown_1 = 22



Mouse[["cellType"]] <- Mouse[["seurat_clusters"]]
Mouse[["cellType"]]  <- ifelse(unlist(Mouse[["seurat_clusters"]]) %in% B_cell, "B_cell",unlist(Mouse[["cellType"]]))
Mouse[["cellType"]]  <- ifelse(unlist(Mouse[["seurat_clusters"]]) %in% Cardiac_Neuronal, "Cardiac_Neuronal",unlist(Mouse[["cellType"]]))
Mouse[["cellType"]]  <- ifelse(unlist(Mouse[["seurat_clusters"]]) %in% Cardiomyocyte, "Cardiomyocyte",unlist(Mouse[["cellType"]]))
Mouse[["cellType"]]  <- ifelse(unlist(Mouse[["seurat_clusters"]]) %in% Endocardial_EC, "Endocardial_EC",unlist(Mouse[["cellType"]]))
Mouse[["cellType"]]  <- ifelse(unlist(Mouse[["seurat_clusters"]]) %in% Fibroblast, "Fibroblast",unlist(Mouse[["cellType"]]))
Mouse[["cellType"]]  <- ifelse(unlist(Mouse[["seurat_clusters"]]) %in% Lymphatic_EC, "Lymphatic_EC",unlist(Mouse[["cellType"]]))
Mouse[["cellType"]]  <- ifelse(unlist(Mouse[["seurat_clusters"]]) %in% Mesothelial, "Mesothelial",unlist(Mouse[["cellType"]]))
Mouse[["cellType"]]  <- ifelse(unlist(Mouse[["seurat_clusters"]]) %in% Myeloid, "Myeloid",unlist(Mouse[["cellType"]]))
Mouse[["cellType"]]  <- ifelse(unlist(Mouse[["seurat_clusters"]]) %in% Pericytes, "Pericytes",unlist(Mouse[["cellType"]]))

Mouse[["cellType"]]  <- ifelse(unlist(Mouse[["seurat_clusters"]]) %in% Smooth_Muscle, "Smooth_Muscle",unlist(Mouse[["cellType"]]))
Mouse[["cellType"]]  <- ifelse(unlist(Mouse[["seurat_clusters"]]) %in% Unknown_hepatocyte, "Unknown_3",unlist(Mouse[["cellType"]]))
Mouse[["cellType"]]  <- ifelse(unlist(Mouse[["seurat_clusters"]]) %in% Coronary_EC, "Coronary_EC",unlist(Mouse[["cellType"]]))


Mouse[["cellType"]] <- ifelse(unlist(Mouse[["seurat_clusters"]]) %in% T_cell, "T_cell",unlist(Mouse[["cellType"]]))


Mouse[["cellType"]] <- ifelse(unlist(Mouse[["seurat_clusters"]]) %in% Macrophage, "Macrophage",unlist(Mouse[["cellType"]]))

Mouse[["cellType"]] <- ifelse(unlist(Mouse[["seurat_clusters"]]) %in% Proliferating_cell, "Unknown_2",unlist(Mouse[["cellType"]]))

Mouse[["cellType"]] <- ifelse(unlist(Mouse[["seurat_clusters"]]) %in% Unknown_1, "Unknown_1",unlist(Mouse[["cellType"]]))

table(Mouse$cellType)
Idents(Mouse) <- Mouse$cellType
#Markers = FindAllMarkers(Mouse, only.pos = TRUE, logfc.threshold = 0.5)
#saveRDS(Markers,file = "Markers_cellType.seurat.rds")
Markers = readRDS("Markers_cellType.seurat.rds")
Markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(Mouse, features = top10$gene)
ggsave("MarkerGenes_cellType.png", height = 40 , width = 40)



#Assign indiviudal Mouse IDS    


setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-03-01.CellbenderFiltered_MAS3_2_Adjust"))


Mouse$sex  = Mouse$sexPredict
Mouse$condition = Mouse$cond



Idents(Mouse) <- Mouse$cellType
tplot <- DimPlot(Mouse, reduction = "umap", label = T, label.size = 4, repel = TRUE)
tplot[[1]]$layers[[1]]$aes_params$alpha = .2
print(paste("Umap with cellTypes"))
tplot
ggsave(tplot+ theme(legend.text=element_text(size=3)),file = "AllCellType.png")

Mouse$ind = paste0(Mouse$run,'_',Mouse$sexPredict)
saveRDS(Mouse,"Seurat_cellbender_annotatedCelltype.rds")





setwd(publication_path("LUNG_QC_ROOT", "analysis/2021-08-20.Entropy/"))

load("clean_nodatasets_060720.RData")

setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-03-01.CellbenderFiltered_MAS3_2_Adjust"))
countTable = Mouse@assays$RNA@counts

condition_list = Mouse@meta.data$seurat_clusters
other = Mouse@meta.data
species = "Mouse"
study = "in House"
data = get("countTable")

Skip = TRUE 
if (Skip == FALSE){
    
    #Build QC object how it is expected
    qc = as.data.frame(matrix(nrow = ncol(data), ncol = 27))
    colnames(qc) = c("timepoint", "study", "depth", "genes", "entropy", "top5", "top5_norm", "depth_norm", "mito", "mito_all", "ribo", "cm_score", "max_score", "max_celltype", "good_cell", "data", "cellname", "other_meta", "include_dataset", "reason", "isolation", "sequencing", "mapping", "datatype", "species", "sample_type", "full_label")
    rownames(qc) = colnames(data)
    qc$study = study
    qc$timepoint = NA
    qc$depth = colSums(data)
    qc$genes = colSums(data > 0)
    qc$top5 = top_n_genes(data, 5)
    qc$top5_norm = 0
    qc$depth_norm = 0
    qc$mito = mito(data, species = species)
    qc$mito_all = mito(data, species = species, all = TRUE)
    qc$ribo = ribo(data, species = species)
    #scn_data = scn(data, qc, study, qc)
    
    #Run SCN without heatmap
    cgenesB = cgenesA[cgenesA %in% rownames(data)]
    xpairs_list = as.data.frame(matrix(unlist(strsplit(xpairs, "_")), nrow=length(xpairs), byrow=T))
    xpairsB = xpairs[xpairs_list$V1 %in% cgenesB & xpairs_list$V2 %in% cgenesB]
    system.time(pdTrain<-query_transform(expT[cgenesB, ], xpairsB))
    system.time(rf_tspAll<-sc_makeClassifier(pdTrain[xpairsB,], genes=xpairsB, groups=grps, nRand=10, ntrees=500))
    print("Done making classifier")
    system.time(kidTransAll<-query_transform(data[cgenesB,], xpairsB))
    print("A")
    system.time(crParkall<-rf_classPredict(rf_tspAll, kidTransAll, numRand=10))
    print("B")

    
    #Assign lableled Predictions
    qc$cm_score = crParkall["cardiac muscle cell",1:(dim(crParkall)[2]-10) ]
    qc$max_score= apply(crParkall[,1:(dim(crParkall)[2]-10)], 2, max)
    qc$max_celltype = rownames(crParkall[,1:(dim(crParkall)[2]-10)])[apply(crParkall[,1:(dim(crParkall)[2]-10)],2,which.max)]
    
    
    
    qc$good_cell = ((qc$top5_norm < 1.3) & (qc$depth_norm > -0.5) )
    #ggsave(file = "celltypes")
    saveRDS(qc, file = "qc.rds")
}

#Read in result
setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-03-01.CellbenderFiltered_MAS3_2_Adjust"))
qc = readRDS( file = "qc.rds")
Mouse[["SingleCellNet"]] <-  qc$max_celltype
keep = names(table(Mouse[["SingleCellNet"]])[table(Mouse[["SingleCellNet"]])> 1000])
Mouse[["SingleCellNetPlot"]] = ifelse((unlist(Mouse[['SingleCellNet']])) %in% keep  ,unlist(Mouse[["SingleCellNet"]]) ,'small')

Idents(Mouse) <- Mouse$SingleCellNetPlot


tplot <- DimPlot(Mouse, reduction = "umap", label = T, label.size = 2)
tplot[[1]]$layers[[1]]$aes_params$alpha = .2
print(paste("SingleCellNet Major predicted Cell Types"))
tplot

ggsave(file = "PredictedCellType.png")


Mouse[["SingleCellNet"]] <-  qc$max_celltype
Idents(Mouse) <- Mouse$SingleCellNet
tplot <- DimPlot(Mouse, reduction = "umap", label = T, label.size = 2)
tplot[[1]]$layers[[1]]$aes_params$alpha = .2
print(paste("SingleCellNet predicted Cell Types"))
tplot+ theme(legend.text=element_text(size=1))
ggsave(tplot+ theme(legend.text=element_text(size=3)),file = "AllPredictedCellType.png")

for (x in sort(unique(Mouse$seurat_clusters))){
  sub = subset(Mouse, seurat_clusters ==x)
  Idents(sub) <- sub$SingleCellNet
  tplot <- DimPlot(sub, reduction = "umap", label = T, label.size = 2)
  tplot[[1]]$layers[[1]]$aes_params$alpha = .2
  print(paste(x,"SingleCellNet predicted Cell Types"))
  print(tplot)
  ggsave(tplot+ theme(legend.text=element_text(size=3)),file = paste0("cluster_",x,"_PredictedCellType.png"))
    
  
  
}


library(Seurat)

setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-03-01.CellbenderFiltered_MAS3_2_Adjust"))
saveRDS(Mouse,file = "CellTypes.seurat.rds")



library(SeuratDisk)
library(seuratTools)
library(DropletUtils)
library(Seurat)
library(ggplot2)
library(scran)
setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-03-01.CellbenderFiltered_MAS3_2_Adjust"))
#Mouse = readRDS("Seurat_cellbender_annotatedCelltype.rds")


HumanHeart <- LoadH5Seurat(publication_path("HUMAN_HEART_ROOT", "data/2021-06-28.HumanHeart/global_raw.h5seurat"))
HumanHeart

HumanHeart =NormalizeData(
    HumanHeart)

new_rownames <- convert_symbols_by_species(src_genes = rownames(HumanHeart), src_species = "human")
seu_slots <- c("counts", "data")
seu = HumanHeart


  for (i in seu_slots) {
    rownames(slot(seu@assays$RNA, i)) <- new_rownames
  }

HumanMouseNames = seu


Idents(HumanMouseNames) =HumanMouseNames$cell_type


feat = rownames(Mouse)[rownames(Mouse) %in% rownames(HumanMouseNames)]


#Already done so load in object
#transfer.anchors <- FindTransferAnchors(
#  reference = HumanMouseNames,
#  query = Mouse,
#  features = feat,
#  k.filter = NA
#)


setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-03-01.CellbenderFiltered_MAS3_2_Adjust"))
#saveRDS(transfer.anchors, file='transfer.anchorsCell_type_HumanHeartNAKfilter.RData')

transfer.anchors = readRDS('transfer.anchorsCell_type_HumanHeartNAKfilter.RData')

  predicted.labels <- TransferData(
      anchorset = transfer.anchors,
      refdata = HumanMouseNames$cell_type,
      dims = 1:30,
    )

Mouse <- AddMetaData(object = Mouse, metadata = predicted.labels)

#plot the distribution of prediction score
png('predictscore_cellType_human.png')
hist(Mouse$prediction.score.max)
abline(v = 0.5, col = "red")
dev.off()


table(Mouse$prediction.score.max > 0.5) #examine the distribution of prediction scores and optionally filter out those cells with low scores.
Mouse.filtered <- subset(Mouse, subset = prediction.score.max > 0.5) #the prediction scores are high which suggest a high degree of confidence in our cell type assignments

#scATAC.filtered$predicted.id <- factor(scATAC.filtered$predicted.id, levels = levels(scRNA))  # to make the colors match
Mouse.filtered$predicted.id <- factor(Mouse.filtered$predicted.id, levels = levels(HumanMouseNames))

plot1 <- DimPlot(
  object = HumanMouseNames,
  group.by = 'cell_type', #seurat_clusters', #"celltype",
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('other snRNA-seq')



plot2 <- DimPlot(
  object = Mouse.filtered,#scATAC.filtered,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE)  + ggtitle('snRNA-seq')


print(paste("Label Trnasfred cellTypes"))
plot1 + plot2

ggsave('integrateUMAP_cellType_human.png', width=16)



#Subset to cm

sub = subset(Mouse, cellType =="Cardiomyocyte")
sub <- FindNeighbors(sub, dims = 1:15) #note dimensions 
sub <- FindClusters(sub, resolution = 0.2)
sub <- RunUMAP(sub, dims = 1:15)
Idents(sub) = Mouse$seurat_clusters

HumanSub = subset(HumanMouseNames, cell_type %in% c('Atrial_Cardiomyocyte','Ventricular_Cardiomyocyte'))
  
  
  
  
  
transfer.anchors <- FindTransferAnchors(
  reference = HumanSub,
  query = sub,
  features = feat,
  k.filter = NA
)


setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-03-01.CellbenderFiltered_MAS3_2_Adjust"))
saveRDS(transfer.anchors, file='CM_ONLY_transfer.anchorsCell_type_HumanHeartNAKfilter.RData')

transfer.anchors = readRDS('CM_ONLY_transfer.anchorsCell_type_HumanHeartNAKfilter.RData')

  predicted.labels <- TransferData(
      anchorset = transfer.anchors,
      refdata = HumanSub$cell_type,
      dims = 1:15,
    )

sub <- AddMetaData(object = sub, metadata = predicted.labels)

#plot the distribution of prediction score
png('CM_only_predictscore_cellType_human.png')
hist(sub$prediction.score.max)
abline(v = 0.5, col = "red")
dev.off()

table(sub$prediction.score.max > 0.5) #examine the distribution of prediction scores and optionally filter out those cells with low scores.
sub.filtered <- subset(sub, subset = prediction.score.max > 0.5) #the prediction scores are high which suggest a high degree of confidence in our cell type assignments

#scATAC.filtered$predicted.id <- factor(scATAC.filtered$predicted.id, levels = levels(scRNA))  # to make the colors match
sub.filtered$predicted.id <- factor(sub.filtered$predicted.id, levels = levels(HumanMouseNames))

plot1 <- DimPlot(
  object = HumanSub,
  group.by = 'cell_type', #seurat_clusters', #"celltype",
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('other snRNA-seq')



plot2 <- DimPlot(
  object = sub.filtered,#scATAC.filtered,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE)  + ggtitle('snRNA-seq')


print(paste("CM Label Transfered cellTypes"))
plot1 + plot2

ggsave('CM_integrateUMAP_cellType_human.png', width=16)


#subcluster
setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-03-01.CellbenderFiltered_MAS3_2_Adjust"))
for (x in unique(Mouse$cellType)){
  sub = subset(Mouse, cellType ==x)
  
  sub <- FindNeighbors(sub, dims = 1:15) #note dimensions 
  sub <- FindClusters(sub, resolution = 0.2)
  sub <- RunUMAP(sub, dims = 1:15)
  Idents(sub) = sub$seurat_clusters
  print(paste(x,"Reclustered"))
  print(DimPlot(sub, reduction = "umap",group.by = "seurat_clusters", label = T))
  ggsave(DimPlot(sub, reduction = "umap",group.by = "seurat_clusters", label = T),file = paste(x,"_Mouse_cluster_15_0.5_umap.png"))
  print(paste(x,"Reclustered individuals"))
  print(DimPlot(sub, reduction = "umap",group.by = "ind", label = T))
  ggsave(DimPlot(sub, reduction = "umap",group.by = "ind", label = T),file = paste0(x,"_Mouse_ind_15_0.5_umap.png"))
  
  Idents(sub) = sub$seurat_clusters
  Markers = FindAllMarkers(sub, only.pos = TRUE, logfc.threshold = 0.5)
  saveRDS(Markers,file = paste0(x,"_Markers_cellType.seurat.rds"))
  Markers %>%
      group_by(cluster) %>%
      top_n(n = 10, wt = avg_log2FC) -> top10
    print(paste(x,"Top Markers"))
    DoHeatmap(sub, features = top10$gene)
    ggsave(paste0(x,"MarkerGenes_cellType.png"), height = 40 , width = 40)
    print(paste(x,"QC"))
    print(VlnPlot(sub, features = c("nFeature_RNA", "nCount_RNA", "Seurat_percent.mt"),group.by = 'seurat_clusters', ncol = 3, pt.size = 0.05))
    ggsave(VlnPlot(sub, features = c("nFeature_RNA", "nCount_RNA", "Seurat_percent.mt"),group.by = 'seurat_clusters', ncol = 3, pt.size = 0.05), file = paste0(x,"_vlnPlot.png"))
}
    


  



setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-03-01.CellbenderFiltered_MAS3_2_Adjust"))
#Markers = FindAllMarkers(Mouse, only.pos = TRUE, logfc.threshold = 0.5)
#saveRDS(Markers,file = "Markers_cellType.seurat.rds")
Idents(Mouse) = Mouse$cellType
Markers = readRDS('Markers_cellType.seurat.rds')
Markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
print("CellTypeMarkers")
DoHeatmap(Mouse, features = top10$gene)+ theme( axis.text.x = 5)
ggsave("MarkerGenes_cellType.png", height = 40 , width = 40)

Idents(Mouse) =Mouse$seurat_clusters
#Markers = FindAllMarkers(Mouse, only.pos = TRUE, logfc.threshold = 0.5)
#saveRDS(Markers,file = "Markers_cluster.seurat.rds")
Markers = readRDS('Markers_cluster.seurat.rds')
Markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
print("Cluster Markers")
DoHeatmap(Mouse, features = top10$gene)+ theme( axis.text.x = 5)
ggsave("MarkerGenes_cluster.png", height = 40 , width = 40)
