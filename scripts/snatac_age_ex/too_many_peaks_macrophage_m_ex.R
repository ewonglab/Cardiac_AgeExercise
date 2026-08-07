# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Prepare macrophage ATAC peaks for TooManyCells hierarchy analysis.
# Steps: 1) Load immune ATAC clusters; 2) select macrophage populations; 3) export the matrix and annotations for hierarchy building.
# =======================================
# Source reformatted from snATAC_AgeEx/Too_Many_peaks_Macrophage_M_Ex.ipynb.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

library(Seurat)
library(Signac)
library(tidyr)
library(dplyr)
library(SingleCellExperiment)
library(biomaRt)
library(scran)
custom_theme <- theme(
         panel.background = element_rect(fill=NA,color=NA),
         plot.background = element_rect(fill=NA, color=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.background = element_rect(fill=NA,color = NA),
         legend.box.background = element_rect(fill=NA,color = NA),
         text=element_text(size=17),
         axis.text = element_text(size = 15),
         axis.line = element_line(size = 1, color = "black")
       )

# Set the custom theme as the default
theme_set(custom_theme)

setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-03-15.Macro"))
Immune = readRDS("Immune_clustering.rds")


DimPlot(Immune, reduction = "umap",group.by = "seurat_clusters")

Macrophages = subset(Immune, seurat_clusters %in% c(0,1,4,8))

Macrophages = FindVariableFeatures(Macrophages)
Macrophages <- RunPCA(Macrophages, features = VariableFeatures(object = Macrophages))

Macrophages <- JackStraw(Macrophages, num.replicate = 100)
Macrophages <- ScoreJackStraw(Macrophages, dims = 1:20)

setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-03-15.Macro"))
print(JackStrawPlot(Macrophages, dims = 1:20))

Macrophages <- FindNeighbors(Macrophages, dims = 2:9) #note dimensions 
Macrophages <- FindClusters(Macrophages, resolution = 0.2)
Idents(Macrophages) = Macrophages$seurat_clusters
Macrophages <- RunUMAP(Macrophages, dims = 2:9)

DimPlot(Macrophages, reduction = "umap",group.by = "seurat_clusters")

meta = data.frame(cbind(colnames(Macrophages),Macrophages$seurat_clusters))
colnames(meta) = c("item","label")
head(meta)
library(stringr)
meta$item = str_replace(meta$item,"-",".")
head(meta)
write.csv(meta, "Macrophage_Full_Cluster_labels.csv",row.names = F, quote = F)

setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-03-15.Macro/"))
meta = data.frame(cbind(colnames(Macrophages),Macrophages$group))
colnames(meta) = c("item","label")
head(meta)
library(stringr)
meta$item = str_replace(meta$item,"-",".")
head(meta)
write.csv(meta, "Macrophage_Full_group_labels.csv",row.names = F, quote = F)

write.csv(slot(Macrophages@assays[["RNA"]],'data'),"Macrophages_Full_data.csv",quote = F)

# Generate the required full, cropped and group hierarchies before reading downstream clustering information:
# bash scripts/common/run_too_many_cells.sh

clustering = read.csv(publication_path("HEART_SNRNA_ROOT", "analysis/2023-03-15.Macro/too-many-cells/Macrophages_full_Crop/cluster_clusters.csv"))
clustering$pathSave = clustering$path
head(clustering)
clustering_split <- separate(clustering, path, into = paste0("cluster", 1:7), sep = "/", fill = "right")
Macrophages$too_many_ID = str_replace(Macrophages$cellID,"-",".")
clustering_split_ordered<- clustering_split %>%
  arrange(match(cell, Macrophages$too_many_ID))
table(clustering_split_ordered$cell == Macrophages$too_many_ID)
table(clustering_split_ordered$cluster5)

#Want to split at 

#node 15 M2
# node 2 M2
# node 26 M1
# node 31 M1
# node 38 UNk

clustering_split_ordered$select = clustering_split_ordered$pathSave
clustering_split_ordered$select = ifelse(grepl(clustering_split_ordered$pathSave, pattern = "/15/"),"15",clustering_split_ordered$select)
clustering_split_ordered$select = ifelse(grepl(clustering_split_ordered$pathSave, pattern = "/2/"),"2",clustering_split_ordered$select)
clustering_split_ordered$select = ifelse(grepl(clustering_split_ordered$pathSave, pattern = "/26/"),"26",clustering_split_ordered$select)
clustering_split_ordered$select = ifelse(grepl(clustering_split_ordered$pathSave, pattern = "/31/"),"31",clustering_split_ordered$select)
clustering_split_ordered$select = ifelse(grepl(clustering_split_ordered$pathSave, pattern = "/38/"),"38",clustering_split_ordered$select)
table(clustering_split_ordered$select)

Macrophages$too_many_cluster = clustering_split_ordered$select
Macrophages$too_many_node = clustering_split_ordered$cluster

DimPlot(Macrophages, reduction = "umap",group.by = "seurat_clusters")
DimPlot(Macrophages, reduction = "umap",group.by = "too_many_cluster")
FeaturePlot(Macrophages, "Mrc1")
FeaturePlot(Macrophages, "Il1b")
FeaturePlot(Macrophages, "H2-Ab1")
VlnPlot(Macrophages, "Cd163",group.by = "seurat_clusters")
VlnPlot(Macrophages, "Cd163",group.by = "too_many_cluster")
VlnPlot(Macrophages, "Cd86",group.by = "too_many_cluster")
VlnPlot(Macrophages, "H2-Ab1",group.by = "seurat_clusters")
VlnPlot(Macrophages, "H2-Ab1",group.by = "too_many_cluster")
VlnPlot(Macrophages, "H2-Eb1",group.by = "too_many_cluster")

Macrophages$Too_cellType = Macrophages$too_many_cluster
Macrophages$Too_cellType = ifelse(Macrophages$too_many_cluster %in% c(2,15),"M2",Macrophages$Too_cellType)
Macrophages$Too_cellType = ifelse(Macrophages$too_many_cluster %in% c(26,31),"M1",Macrophages$Too_cellType)

Macrophages$cluster_cellType = Macrophages$seurat_clusters
Macrophages$cluster_cellType = ifelse(Macrophages$seurat_clusters %in% c(1),"M1",Macrophages$cluster_cellType)
Macrophages$cluster_cellType = ifelse(Macrophages$seurat_clusters %in% c(0,3),"M2",Macrophages$cluster_cellType)
table(Macrophages$cluster_cellType,Macrophages$group)

pt_ind <- table(Macrophages$Too_cellType, Macrophages$group)
pt <- as.data.frame(pt_ind)
M1_M2 = pt[pt$Var1 %in% c('M1','M2'),]
M1_M2
MAE_M1_M2 = c(296,786)
MAS_M1_M2 = c(1588,2185)
YE_M1_M2 = c(601,1384)
YS_M1_M2 = c(1308,1982)
a= chisq.test(x = rbind(MAE_M1_M2,MAS_M1_M2))
a$expected
a$observed
a
a= chisq.test(x = rbind(YE_M1_M2,YS_M1_M2))
a$expected
a$observed
a
a= chisq.test(x = rbind(MAE_M1_M2,YE_M1_M2))
a$expected
a$observed
a
a= chisq.test(x = rbind(MAS_M1_M2,YS_M1_M2))
a$expected
a$observed
a



Cluster_info = read.csv(publication_path("HEART_SNRNA_ROOT", "analysis/2023-03-15.Macro/too-many-cells/Macrophages_Full_group/node_info.csv"))

Cluster_info$compositionSave = Cluster_info$composition
Cluster_info_split <- separate(Cluster_info, composition, into = c("MAE","MAS","YE","YS"), sep = "/", fill = "right")
Cluster_info_split$MAE = str_remove(Cluster_info_split$MAE, pattern = "MAE:")
Cluster_info_split$MAE = str_remove(Cluster_info_split$MAE, pattern = "\\(.*\\)")
Cluster_info_split$YE = str_remove(Cluster_info_split$YE, pattern = "YE:")
Cluster_info_split$YE = str_remove(Cluster_info_split$YE, pattern = "\\(.*\\)")
Cluster_info_split$YS = str_remove(Cluster_info_split$YS, pattern = "YS:")
Cluster_info_split$YS = str_remove(Cluster_info_split$YS, pattern = "\\(.*\\)")
Cluster_info_split$MAS = str_remove(Cluster_info_split$MAS, pattern = "MAS:")
Cluster_info_split$MAS = str_remove(Cluster_info_split$MAS, pattern = "\\(.*\\)")
Cluster_info_split$MAS = as.numeric(Cluster_info_split$MAS)
Cluster_info_split$YS = as.numeric(Cluster_info_split$YS)
Cluster_info_split$YE = as.numeric(Cluster_info_split$YE)
Cluster_info_split$MAE = as.numeric(Cluster_info_split$MAE)
#DT::datatable(Cluster_info_split)

MAE_Large = Cluster_info_split[Cluster_info_split$MAE > 0.5,]
MAS_Large = Cluster_info_split[Cluster_info_split$MAS > 0.5,]
YE_Large = Cluster_info_split[Cluster_info_split$YE > 0.5,]
YS_Large = Cluster_info_split[Cluster_info_split$YS > 0.5,]

MAE_Large
MAS_Large
YE_Large
YS_Large

Ex_Large = Cluster_info_split[(Cluster_info_split$MAE + Cluster_info_split$YE)> 0.75,]
Sed_Large =Cluster_info_split[(Cluster_info_split$MAS + Cluster_info_split$YS)> 0.75,]
Young_Large =Cluster_info_split[(Cluster_info_split$YE + Cluster_info_split$YS)> 0.75,]
Old_Large =Cluster_info_split[(Cluster_info_split$MAE + Cluster_info_split$MAS)> 0.75,]



Ex_Large
Sed_Large
Young_Large
Old_Large

VlnPlot(Macrophages, "Etv6",group.by = "Too_cellType",split.by = "group")+ geom_boxplot(width = 0.3, position=position_dodge(0.9))
VlnPlot(Macrophages, "Hivep2",group.by = "Too_cellType",split.by = "group")+ geom_boxplot(width = 0.3, position=position_dodge(0.9))
VlnPlot(Macrophages, "Nfkb1",group.by = "Too_cellType",split.by = "group")+ geom_boxplot(width = 0.3, position=position_dodge(0.9))


for (x in unique( Macrophages$group)){
    sub = subset(Macrophages , group == x)
    a = FeaturePlot(Macrophages,"Xkr4")


      dat <- data.frame(sub@meta.data)
          dat$UMAP1 <- Embeddings(sub, "umap")[,1]
          dat$UMAP2 <- Embeddings(sub, "umap")[,2]
          dat$group <- factor(dat$group, levels=c("MAE", "MAS", "YE", "YS"))
          dat$sex <- factor(dat$sex, levels=c("M", "F"))
          dat$cond <- factor(dat$cond, levels=c("E", "S"))
          dat$age <- factor(dat$age, levels=c("MA", "Y"))
          dat$ind <- factor(dat$ind, levels=unique(dat$ind))





    pl = ggplot(dat, aes(UMAP1, UMAP2)) + 
        ylim(layer_scales(a)$y$range$range*1.3)+
        xlim(layer_scales(a)$x$range$range*1.3)+
        stat_density_2d(
          geom = "raster",
          aes(alpha = after_stat(ndensity)*10 ,fill = after_stat(ndensity)),
          contour = FALSE
        )+ scale_fill_gradientn(colors = c("yellow","orange","red","darkred")) +
         scale_alpha_continuous(range = c(0, 1), limits = c(0, 3), 
                                 guide = guide_none())+ ggtitle(x)+ custom_theme
    print(pl)
    #ggsave(paste0(x,"_Density.pdf"))
}



pt_ind <- table(Macrophages$Too_cellType, Macrophages$group)
pt_ind <- as.data.frame(pt_ind)
pt_ind$CellType <- as.character(pt_ind$Var1)
pt_ind$CellType = factor(pt_ind$CellType)



df = pt_ind
df_total <- df %>%
  group_by(Var2) %>%
  summarise(Total = sum(Freq))

df <- left_join(df, df_total, by = "Var2")

df$percent = df$Freq/df$Total*100
df

pie_chart <- ggplot(df, aes(x = "", y = percent, fill = CellType)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(Freq, " (", round(percent, 1), "%)")),
            position = position_stack(vjust = 0.5)) +
  labs(title = "Pie Chart with Counts and Percentages", fill = "Category") +
  #scale_fill_manual(values = M1_M2_Pallate) +
  theme_minimal() +
  theme(legend.position = "right")+ facet_wrap(~ Var2)

# Display the pie chart
print(pie_chart)

#ggsave("PieChartratio.pdf")

# Create a pie chart with counts and percentages

getwd()
saveRDS(Macrophages,"Macrophages_Toomany.rds")

setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-03-15.Macro"))
Macrophages = readRDS("Macrophages_Toomany.rds")

table(Macrophages$too_many_cluster)

#Add labels onto Immune population 
DimPlot(Immune, reduction = "umap",group.by = "seurat_clusters")
Immune$Too_Many_Type =Immune$seurat_clusters
Immune$Too_Many_Type = ifelse(Immune$seurat_clusters %in% c(3,6),"T_cell",Immune$Too_Many_Type)
Immune$Too_Many_Type = ifelse(Immune$seurat_clusters %in% c(2),"B_cell",Immune$Too_Many_Type)
Immune$Too_Many_Type = ifelse(Immune$seurat_clusters %in% c(7),"Dendritic",Immune$Too_Many_Type)
Immune$Too_Many_Type = ifelse(Immune$seurat_clusters %in% c(5),"HSPC",Immune$Too_Many_Type)
Immune$Too_Many_Type = ifelse(Immune$cellID %in% Macrophages$cellID[Macrophages$Too_cellType == "M1"],"M1",Immune$Too_Many_Type)
Immune$Too_Many_Type = ifelse(Immune$cellID %in% Macrophages$cellID[Macrophages$Too_cellType == "M2"],"M2",Immune$Too_Many_Type)
Immune$Too_Many_Type = ifelse(Immune$cellID %in% Macrophages$cellID[Macrophages$Too_cellType == "38"],"M_Ex",Immune$Too_Many_Type)
DimPlot(Immune, reduction = "umap",group.by = "Too_Many_Type")
DimPlot(Macrophages, reduction = "umap",group.by = "Too_cellType")
DimPlot(Macrophages, reduction = "umap",group.by = "too_many_cluster")
DimPlot(Macrophages, reduction = "umap",group.by = "too_many_node")
table(Macrophages$too_many_node)

saveRDS(Macrophages,"Macrophages_Toomany.rds")
saveRDS(Immune,"Immune_clustering.rds")

FeaturePlot(Macrophages, "Ccr2", reduction = "umap")
FeaturePlot(Macrophages, "Cd86", reduction = "umap")
FeaturePlot(Macrophages, "Cd163", reduction = "umap")
FeaturePlot(Macrophages, "Cx3cr1", reduction = "umap")
FeaturePlot(Macrophages, "H2-Ab1", reduction = "umap")
FeaturePlot(Macrophages, "Ly6c1", reduction = "umap")
FeaturePlot(Macrophages, "Lgals3", reduction = "umap")

mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))


sce <-SingleCellExperiment(assays=list(counts=Macrophages[['RNA']]@counts), colData=Macrophages@meta.data)


mouse <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

genes = as.data.frame(row.names(sce))
colnames(genes) <- "mgi_symbol"
#genes
genesAsEns <- getBM(filters= "mgi_symbol", attributes= c("ensembl_gene_id","mgi_symbol"),values=genes,mart= mouse)
#genesAsEns





colnames(genesAsEns) <- c("ensembl_gene_id", "mgi_symbol")
both <- merge(genes, genesAsEns, all.x = TRUE, sort = FALSE)
dupRemoved = distinct(both,mgi_symbol, .keep_all= TRUE)
#remove duplicate name pairs
dupRemoved$ensembl_gene_id <- ifelse(is.na(dupRemoved$ensembl_gene_id), dupRemoved$mgi_symbol, dupRemoved$ensembl_gene_id)
reorder_idx <- match(genes$mgi_symbol,dupRemoved$mgi_symbol) 
ordered = dupRemoved[reorder_idx,]
rownames(ordered) <- c(1:nrow(ordered))

rownames(sce) <- ordered$ensembl_gene_id

#run cyclone
assignments <- cyclone(sce, mm.pairs)
saveRDS(assignments,publication_path("HEART_SNRNA_ROOT", "analysis/2023-03-15.Macro/assignemnts.rds"))

assignments = readRDS(publication_path("HEART_SNRNA_ROOT", "analysis/2023-03-15.Macro/assignemnts.rds"))

cycles = as.data.frame(assignments)
#cycles
Macrophages$cycle = cycles$phases
Macrophages$G1 = as.numeric(cycles$normalized.scores.G1)
Macrophages$G2 = as.numeric(cycles$normalized.scores.G2)
Macrophages$S = as.numeric(cycles$normalized.scores.S)

DimPlot(Macrophages, reduction = "umap",group.by = "cycle")
DimPlot(Macrophages, reduction = "umap",group.by = "G1")
DimPlot(Macrophages, reduction = "umap",group.by = "G2")
DimPlot(Macrophages, reduction = "umap",group.by = "S")

setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-03-15.Macro/"))
Macrophages = readRDS("Macrophages_Toomany.rds")

round(ncol(Macrophages)/10)
i = 1

table(downsampled.obj$too_many_node)

setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-03-15.Macro/too-many-cells/Bootstrap/"))
fold = paste0("/replicate_",i)
dir.create(file.path(getwd(), fold))
setwd(file.path(getwd(), fold))
downsampled.obj <- Macrophages[, sample(colnames(Macrophages), size = round(ncol(Macrophages)), replace=T)]
write.csv(slot(downsampled.obj@assays[["RNA"]],'data'),paste0("replicate_",i,"_data.csv"),quote = F)


meta = data.frame(cbind(colnames(downsampled.obj),downsampled.obj$too_many_node))
colnames(meta) = c("item","label")
head(meta)
library(stringr)
meta$item = str_replace(meta$item,"-",".")
head(meta)
write.csv(meta, paste0("replicate_",i,"_label.csv"),row.names = F, quote = F)



# Run `bash scripts/common/run_too_many_cells_bootstrap.sh <replicate-dir>` here.

i =2
setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-03-15.Macro/too-many-cells/Bootstrap/"))
fold = paste0("/replicate_",i)
dir.create(file.path(getwd(), fold))
setwd(file.path(getwd(), fold))
downsampled.obj <- Macrophages[sample(rownames(Macrophages), size = round(nrow(Macrophages)), replace=T),]
write.csv(slot(downsampled.obj@assays[["RNA"]],'data'),paste0("replicate_",i,"_data.csv"),quote = F)


meta = data.frame(cbind(colnames(downsampled.obj),downsampled.obj$too_many_node))
colnames(meta) = c("item","label")
head(meta)
library(stringr)
meta$item = str_replace(meta$item,"-",".")
head(meta)
write.csv(meta, paste0("replicate_",i,"_label.csv"),row.names = F, quote = F)



clustering = read.csv(publication_path("HEART_SNRNA_ROOT", "analysis/2023-03-15.Macro/too-many-cells/Macrophages_full_Crop/cluster_clusters.csv"))
clustering$pathSave = clustering$path

clustering_split <- separate(clustering, path, into = paste0("cluster", 1:7), sep = "/", fill = "right")


OG_Cluster  = unique(unlist(clustering_split[,3:9]))
OG_Cluster = as.numeric(OG_Cluster)
OG_Cluster = na.omit(OG_Cluster)
OG_Cluster

#Get cells in eahc cluster

OG_Cluster_Cells = list()
for (cluster in OG_Cluster){

    rows = clustering_split %>% filter_all(any_vars(. %in% !!cluster))
    OG_Cluster_Cells[[paste0('cluster',cluster)]] = rows$cell
    
}
OG_Cluster_Cells <- OG_Cluster_Cells[order(sapply(OG_Cluster_Cells, length), decreasing = TRUE)]

lapply(OG_Cluster_Cells, length)


setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-03-15.Macro/too-many-cells/"))
      Colours = c('#A020F0','#E41A1C','#377EB8','#00008B')
blank_theme <- theme_minimal()+
  theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  )


library(stringr)
ClusteCond <- lapply(OG_Cluster_Cells, function(x) str_remove(x, "\\d.*"))

proportions = prop.table(table((ClusteCond[[1]])))
proportions
numbers = table((ClusteCond[[1]]))
sig = c()
for (x in    1:length(ClusteCond)){
    
    cluster = ClusteCond[[x]]
    name = names(ClusteCond)[x]
    NormalProp = prop.table(table((cluster))/proportions)
    print(name)
    print(table((cluster)))
    print(NormalProp)
    
    chis = (chisq.test(x =table(cluster), p =proportions))
    p.values <- chis$p.value
Signif <- symnum(p.values, corr = FALSE, na = FALSE, cutpoints = c(0, 
    0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", "NS"))

    data = data.frame(NormalProp)
    sig = c(sig,Signif)

# Create the pie chart using ggplot2
    data = data.frame(NormalProp)
    data
    
    # Create the pie chart using ggplot2
    pie_chart <- ggplot(data, aes(x = "", y = Freq, fill = Var1)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar(theta = "y") +
      labs(title = paste(name,"Ncells =", length(cluster), "\n P = ",Signif), fill = "") +
      scale_fill_manual(values = Colours) +   theme_minimal()   + blank_theme + theme(axis.text.x=element_blank())+ theme(legend.position = "none", plot.title = element_text(size = 25, hjust = 0.5)) 


    # Print the pie chart
    print(pie_chart)
    ggsave(paste0(name,"_pie.pdf"))
    
}
  

getwd()

p.values <- c(9.5e-15)
Signif <- symnum(p.values, corr = FALSE, na = FALSE, cutpoints = c(0, 
    0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))

Signif

library(ggplot2)  

Colours = c('#A020F0','#E41A1C','#377EB8','#00008B')
# Custom labels for categories
data = data.frame(NormalProp)
  data

# Create the pie chart using ggplot2
pie_chart <- ggplot(data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = paste("Ncells =", length(cluster)), fill = "") +
  scale_fill_manual(values = Colours) +   theme_minimal()   + blank_theme + theme(axis.text.x=element_blank())+ theme(legend.position = "none", plot.title = element_text(size = 25, hjust = 0.5)) 
 

# Print the pie chart
print(pie_chart)
