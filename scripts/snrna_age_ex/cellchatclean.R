# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Infer and compare cell-cell communication across experimental conditions.
# Steps: 1) Prepare CellChat objects by condition; 2) estimate interaction probabilities; 3) export inflammatory and collagen pathway summaries.
# =======================================
# Source reformatted from snRNA_AgeEx/CellChatClean.ipynb.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))


library(rtracklayer)
library(GenomicRanges)
library(dittoSeq)
library(ggsignif)
library(scales)
library(CellChat, lib.loc=publication_path("R_LIBS_DIR"))
library(patchwork)
library(rentrez)
library(Seurat)

# Define the cellTypes colouring 

RNA_pallate = c(
              
       '#F99349', #M1
         '#87CEEB',      #M2
       '#FDDA62'  ,      #MEX
        '#808080' ,       #HSPC
        '#BEBADA',        #DEND
    "#A81B1B", #CM
    "#101369", #C_EC
    "#606692" ,#E_EC
    "#3A85A8", #L_EC
    "#42977E", #Fibro
    "#591f0a", #Meso
    "#3b7818", #B_Cell
    "#B6D66A", #T_Cell
    "#1CE382", #Pericytes
    "#BF862B", #SM
    "#CC6A6F", #CN
    "#E086B5", #UN1
    "#5F1296") #UN2


RNA_order = c("M1","M2","M_Ex","HSPC","Dendritic","Cardiomyocyte","Coronary_EC",
  "Endocardial_EC","Lymphatic_EC","Fibroblast",
 "Mesothelial", "B_cell",
 "T_cell","Pericytes","Smooth_Muscle","Cardiac_Neuronal","Unknown_1","Unknown_3")

names(RNA_pallate) = RNA_order

# Redifne some Heatmap code to make work as inteded Orginal from CellChat 
library(ComplexHeatmap)
library(grid)
#' @export
netVisual_heatmap2 <- function(object, comparison = c(1,2), measure = c("count", "weight"), signaling = NULL, slot.name = c("netP", "net"), color.use = NULL, color.heatmap = c("#2166ac","#b2182b"),
                              title.name = NULL, width = NULL, height = NULL, font.size = 8, font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE,
                              sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, row.show = NULL, col.show = NULL, minCol = NULL, maxCol = NULL, row_order = NULL, col_order = NULL,
                              bar_row_colour = NULL, bar_col_colour = NULL ){
  # obj1 <- object.list[[comparison[1]]]
  # obj2 <- object.list[[comparison[2]]]
  if (!is.null(measure)) {
    measure <- match.arg(measure)
  }
  slot.name <- match.arg(slot.name)
  if (is.list(object@net[[1]])) {
    message("Do heatmap based on a merged object \n")
    
    obj1 <- object@net[[comparison[1]]][[measure]]
    obj2 <- object@net[[comparison[2]]][[measure]]
    net.diff <- obj2 - obj1
    #print(net.diff)
    if (measure == "count") {
      if (is.null(title.name)) {
        title.name = "Differential number of interactions"
      }
    } else if (measure == "weight") {
      if (is.null(title.name)) {
        title.name = "Differential interaction strength"
      }
    }
    legend.name = "Relative values"
  } else {
    message("Do heatmap based on a single object \n")
    if (!is.null(signaling)) {
      net.diff <- slot(object, slot.name)$prob[,,signaling]
      if (is.null(title.name)) {
        title.name = paste0(signaling, " signaling network")
      }
      legend.name <- "Communication Prob."
    } else if (!is.null(measure)) {
      net.diff <- object@net[[measure]]
      if (measure == "count") {
        if (is.null(title.name)) {
          title.name = "Number of interactions"
        }
      } else if (measure == "weight") {
        if (is.null(title.name)) {
          title.name = "Interaction strength"
        }
      }
      legend.name <- title.name
    }
  }

  net <- net.diff
  #print(net)

  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0

  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    if (length(idx) > 0) {
      net <- net[-idx, ]
      net <- net[, -idx]
    }
  }

  mat <- net
    
    #print(mat)
  if (is.null(color.use)) {
    color.use <- scPalette(ncol(mat))
  }
  names(color.use) <- colnames(mat)

  if (!is.null(row.show)) {
    mat <- mat[row.show, ]
  }
  if (!is.null(col.show)) {
    mat <- mat[ ,col.show]
    color.use <- color.use[col.show]
  }
  
    if(is.null(minCol)){
        minCol =min(mat) 
    }

    if(is.null(maxCol)){
        maxCol =max(mat)
    }
   #print( paste("MINMAX col =", minCol, maxCol))
  if (min(mat) < 0) {
      print('inside2')
    color.heatmap.use = colorRamp3(c(minCol, 0, maxCol), c(color.heatmap[1], "#f7f7f7", color.heatmap[2]))
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",min(mat, na.rm = T)))+1), 0, round(max(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",max(mat, na.rm = T)))+1))
    print(colorbar.break)
      # color.heatmap.use = colorRamp3(c(seq(min(mat), -(max(mat)-min(max(mat)))/9, length.out = 4), 0, seq((max(mat)-min(max(mat)))/9, max(mat), length.out = 4)), RColorBrewer::brewer.pal(n = 9, name = color.heatmap))
  } else {
    if (length(color.heatmap) == 3) {
      color.heatmap.use = colorRamp3(c(0, minCol, maxCol), color.heatmap)
    } else if (length(color.heatmap) == 2) {
      color.heatmap.use = colorRamp3(c(minCol, maxCol), color.heatmap)
    } else if (length(color.heatmap) == 1) {
      color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100)
    }
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",min(mat, na.rm = T)))+1), round(max(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",max(mat, na.rm = T)))+1))
  }
  # col_fun(as.vector(mat))
   #message(min(mat))
    #message(max(mat))
   if (!is.null(bar_col_colour)) {
     bar.col.use = bar_col_colour
  } else{
       bar.col.use = color.use
 }
    
     if (!is.null(bar_row_colour)) {
     bar.row.use = bar_row_colour
  } else{
       bar.row.use = color.use
 }
    
  df<- data.frame(group = colnames(mat)); rownames(df) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),which = "column",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), which = "row",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))

  ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), border = FALSE,gp = gpar(fill = bar.row.use, col=bar.row.use)), show_annotation_name = FALSE)
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), border = FALSE,gp = gpar(fill = bar.col.use, col=bar.col.use)), show_annotation_name = FALSE)
 
   # print(sum(abs(mat)))
   if ((sum(mat != 0) == 1)){ 
       print("single")
     color.heatmap.use = c("white", 'red')
  }else if (sum(abs(mat) > 0) == 1) {
      message("inside2")
    color.heatmap.use = c("white", color.heatmap.use)
  } else {
    mat[is.na(mat)] <- 0
  }
  #  print(mat)
    #print(mat)
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = legend.name, row_order = row_order, column_order  = col_order,
                bottom_annotation = col_annotation, left_annotation =row_annotation, top_annotation = ha2, right_annotation = ha1,
                cluster_rows = cluster.rows,cluster_columns = cluster.cols,
                row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
               # width = unit(width, "cm"), height = unit(height, "cm"),
                column_title = title.name,column_title_gp = gpar(fontsize = font.size.title),column_names_rot = 90,
                row_title = "Sources (Sender)",row_title_gp = gpar(fontsize = font.size.title),row_title_rot = 90,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                            border = NA, #at = colorbar.break,
                                            legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
  )
    
   #  draw(ht1)
  return(ht1)
}


# Account for cell changes

aggregateNet2 <- function(object, sources.use = NULL, targets.use = NULL, signaling = NULL, pairLR.use = NULL, remove.isolate = TRUE, thresh = 0.05, return.object = TRUE) {
  net <- object@net
  if (is.null(sources.use) & is.null(targets.use) & is.null(signaling) & is.null(pairLR.use)) {
    prob <- net$prob
    pval <- net$pval
    pval[prob == 0] <- 1
    prob[pval >= thresh] <- 0
    net$count <- apply(prob > 0, c(1,2), sum)
    net$weight <- apply(prob, c(1,2), sum)
    net$weight[is.na(net$weight)] <- 0
    net$count[is.na(net$count)] <- 0
  } else {
    df.net <- subsetCommunication(object, slot.name = "net",
                                  sources.use = sources.use, targets.use = targets.use,
                                  signaling = signaling,
                                  pairLR.use = pairLR.use,
                                  thresh = thresh)
    df.net$source_target <- paste(df.net$source, df.net$target, sep = ":")
    df.net2 <- df.net %>% group_by(source_target) %>% summarize(count = n(), .groups = 'drop')
    df.net3 <- df.net %>% group_by(source_target) %>% summarize(prob = sum(prob), .groups = 'drop')
    df.net2$prob <- df.net3$prob
    a <- stringr::str_split(df.net2$source_target, ":", simplify = T)
    df.net2$source <- as.character(a[, 1])
    df.net2$target <- as.character(a[, 2])
    cells.level <- levels(object@idents)
    if (remove.isolate) {
      message("Isolate cell groups without any interactions are removed. To block it, set `remove.isolate = FALSE`")
      df.net2$source <- factor(df.net2$source, levels = cells.level[cells.level %in% unique(df.net2$source)])
      df.net2$target <- factor(df.net2$target, levels = cells.level[cells.level %in% unique(df.net2$target)])
    } else {
      df.net2$source <- factor(df.net2$source, levels = cells.level)
      df.net2$target <- factor(df.net2$target, levels = cells.level)
    }

    count <- tapply(df.net2[["count"]], list(df.net2[["source"]], df.net2[["target"]]), sum)
    prob <- tapply(df.net2[["prob"]], list(df.net2[["source"]], df.net2[["target"]]), sum)
    net$count <- count
    net$weight <- prob
    net$weight[is.na(net$weight)] <- 0
    net$count[is.na(net$count)] <- 0
  }
  if (return.object) {
    object@net <- net
    return(object)
  } else {
    return(net)
  }

}

# Read in RNA Object and remove Unknown celltypes 
RNA = readRDS(publication_path("HEART_SNRNA_ROOT", "data/2022-10-17.SharedObject/2023-10-9_Full_RNA_TooManyAnnotated.rds"))
RNA = subset(RNA, CellType %in% c("Cardiomyocyte","Cardiac_Neuronal","Mesothelial",
  "M2","HSPC","Dendritic",
 "B_cell","M_Ex","Fibroblast",
"T_cell","Endocardial_EC","M1",
 "Lymphatic_EC","Pericytes","Coronary_EC","Smooth_Muscle"
))

RNA$CellTypeMacro = ifelse(RNA$CellType %in% c("M1","M2","M_Ex"),"Macrophage", RNA$CellType)

Idents(RNA) = 'CellTypeMacro'

# Define Database 
CellChatDB <- CellChatDB.mouse 
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB

# Run cell chat on each condition
for (cond in unique(RNA$age_condition)){
    seurat_object = subset(RNA, age_condition == cond)
    print(seurat_object)
    data.input <- GetAssayData(seurat_object, assay = "RNA", slot = "data") # normalized data matrix
    labels <- Idents(seurat_object)
    meta <- data.frame(age_condition = labels, row.names = names(labels)) # create a dataframe of the cell labels
    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "age_condition")

    cellchat <- addMeta(cellchat, meta = meta)
    cellchat <- setIdent(cellchat, ident.use = "age_condition") # set "labels" as default cell identity
    print(levels(cellchat@idents)) # show factor levels of the cell labels
    groupSize <- as.numeric(table(cellchat@idents))

    cellchat@DB <- CellChatDB.use
    cellchat <- subsetData(cellchat)
    future::plan(if (future::supportsMulticore()) "multicore" else "multisession", workers = 7) # do parallel ("multiprocess" is defunct in current 'future')
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    groupSize <- as.numeric(table(cellchat@idents))
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
    netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2024-02-22.CellChat_Clean_Macro"))
    saveRDS(cellchat, paste0(cond,"_Filt_cellChat.rds"))
}

# Read in the shared results and merge togetehr 
setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2024-02-22.CellChat_Clean"))
    Aged_Ex = readRDS("Aged_Ex_Filt_cellChat.rds")
    Aged_Sed = readRDS("Aged_Sed_Filt_cellChat.rds")
    Young_Ex = readRDS("Young_Ex_Filt_cellChat.rds")
    Young_Sed = readRDS("Young_Sed_Filt_cellChat.rds")
object.list <- list(Aged_Ex = Aged_Ex, Aged_Sed = Aged_Sed,Young_Ex = Young_Ex,Young_Sed = Young_Sed)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

object.list <- list(Aged_Ex = Aged_Ex, Aged_Sed = Aged_Sed,Young_Ex = Young_Ex,Young_Sed = Young_Sed)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Write out all signficant interactions 
sources.use = NULL
 targets.use = NULL
signaling = NULL
pairLR.use = NULL
thresh = 0.05
#netVisual_bubble(cellchat, sources.use = 6, targets.use = NULL,  comparison = c(2, 1), angle.x = 45, thresh = 0.1)
df.net <- subsetCommunication(cellchat, slot.name = "net",
                                  sources.use = sources.use, targets.use = targets.use,
                                  signaling = signaling,
                                  pairLR.use = pairLR.use,
                                  thresh = thresh)

df.net$Aged_Ex$Cond = "Aged_Ex"
df.net$Aged_Sed$Cond = "Aged_Sed"  # Fixed: was mislabelled "Aged_Ex", so the exported table had two sets of rows called Aged_Ex
df.net$Young_Ex$Cond = "Young_Ex"
df.net$Young_Sed$Cond = "Young_Sed"
df_all = rbind(df.net$Aged_Ex,df.net$Aged_Sed,df.net$Young_Ex,df.net$Young_Sed)

nrow(df_all)
write.csv(df_all,"2024-02-22.Significant_Filt_interactions_allCell.csv")

# Generate Initial heatmap of everthing 
RNA_pallate = (RNA_pallate)[rownames(cellchat@net$Aged_Ex$weight)]
unorder = netVisual_heatmap2(cellchat,comparison =  c("Young_Sed", "Aged_Sed"), minCol = -1, maxCol = 1, measure = "weight", cluster.rows = FALSE, cluster.cols = FALSE, color.use = RNA_pallate)
unorder

RNA_pallate = (RNA_pallate)[rownames(cellchat@net$Aged_Ex$weight)]
RNA_pallate_bar_row = (RNA_pallate)[row_order(unorder)]
RNA_pallate_bar_col = (RNA_pallate)[column_order(unorder)]

# Generate Initial heatmap of Age comparison 
Age <- netVisual_heatmap2(cellchat,comparison =  c("Young_Sed", "Aged_Sed"), minCol = -1, maxCol = 1, measure = "weight", cluster.rows = TRUE, cluster.cols = TRUE, color.use = RNA_pallate,
                          bar_row_colour = RNA_pallate_bar_row, bar_col_colour = RNA_pallate_bar_col
                         )
Age
#GetHere














# Generate Initial heatmap of Age comparison Loosing the clustering
age <- netVisual_heatmap2(cellchat,comparison =  c("Young_Sed", "Aged_Sed"), measure = "weight", minCol = -1.07, maxCol = 1.07, cluster.rows = FALSE , row_order = row_order(Age), col_order = column_order(Age),
                          bar_row_colour = RNA_pallate_bar_row, bar_col_colour = RNA_pallate_bar_col,color.use = RNA_pallate)
pdf("Age_orderd_Filt_Change_Weights_Heatmap.pdf")
age
dev.off()
age

# Generate Initial heatmap of Ex comparison 
ex <- netVisual_heatmap2(cellchat,comparison =  c("Aged_Sed", "Aged_Ex"), measure = "weight",  minCol = -1.6, maxCol = 1.6, cluster.rows = FALSE , row_order = row_order(Age), col_order = column_order(Age),
                         bar_row_colour = RNA_pallate_bar_row, bar_col_colour = RNA_pallate_bar_col,color.use = RNA_pallate)
pdf("AgeEx_orderd_Filt_Change_Weights_Heatmap.pdf")
ex
dev.off()

ex

ex

head(read.csv(publication_path("HEART_SNRNA_ROOT", "analysis/2023-04-27.Karthik/GO_term_summary_20230907_223136.csv"), header = FALSE))

# Get the genes in the interactionDB
genes = read.csv(publication_path("HEART_SNRNA_ROOT", "analysis/2023-04-27.Karthik/GO_term_summary_20230907_223136.csv"), header = FALSE)
genes=genes[-c(1:2),]
head(genes)

posGenes = genes[genes$V6 == "positive regulation of inflammatory response",]
negGenes = genes[genes$V6 == "negative regulation of inflammatory response",]
genes = genes[genes$V6 == "positive regulation of inflammatory response",]

Remove = (posGenes$V2[posGenes$V2 %in% negGenes$V2])

length(unique(genes$V2))
genes = genes$V2
genes = genes[!(genes  %in% Remove)]


genes = unique(genes)
#head(cellchat@DB)
ValidGenes = unique(c(cellchat@DB$interaction$ligand,cellchat@DB$interaction$receptor,unlist(cellchat@DB$complex),unlist(cellchat@DB$cofactor)))



# Get the genes in the interactionDB
genes = read.csv(publication_path("HEART_SNRNA_ROOT", "analysis/2023-04-27.Karthik/GO_term_summary_20230907_223136.csv"), header = FALSE)
genes=genes[-c(1:2),]
head(genes)

posGenes = genes[genes$V6 == "positive regulation of inflammatory response",]
negGenes = genes[genes$V6 == "negative regulation of inflammatory response",]
genes = genes[genes$V6 == "positive regulation of inflammatory response",]

Remove = (posGenes$V2[posGenes$V2 %in% negGenes$V2])

length(unique(genes$V2))
genes = genes$V2
genes = genes[!(genes  %in% Remove)]


genes = unique(genes)
#head(cellchat@DB)
ValidGenes = unique(c(cellchat@DB$interaction$ligand,cellchat@DB$interaction$receptor,unlist(cellchat@DB$complex),unlist(cellchat@DB$cofactor)))

complex = cellchat@DB$complex
cofactor = cellchat@DB$cofactor
length(ValidGenes)
table(genes %in% ValidGenes)
KeptGenes = genes[genes %in% ValidGenes]
ligands = KeptGenes[(KeptGenes %in% cellchat@DB$interaction$ligand)]
receptors = KeptGenes[(KeptGenes %in% cellchat@DB$interaction$receptor)]

Complex_Cofactor = KeptGenes[!(KeptGenes %in% c(ligands,receptors))]



#Check if gene is in the cofactor and keep, then grab all interactions of intrerest 
length(Complex_Cofactor)
unit = c()
for (gene in Complex_Cofactor){
    hits = complex[gene == complex$subunit_1,]
    unit = c(unit, rownames(hits))
    hits = complex[gene == complex$subunit_2,]
    unit = c(unit, rownames(hits))
    hits = complex[gene == complex$subunit_3,]
    unit = c(unit, rownames(hits))
    hits = complex[gene == complex$subunit_4,]
    unit = c(unit, rownames(hits))  
}

for (gene in Complex_Cofactor){
    for (col in 1:ncol(cofactor)){
        hits = cofactor[gene == cofactor[,col],]
        unit = c(unit, rownames(hits))
        
    }
}

table(cellchat@DB$interaction$ligand %in% c(KeptGenes,unit))
Interactions_of_intrest = cellchat@DB$interaction[cellchat@DB$interaction$ligand %in% c(KeptGenes,unit) |cellchat@DB$interaction$receptor %in% c(KeptGenes,unit),]

# Save pathways of interest 
Pathways = data.frame(unique(Interactions_of_intrest$interaction_name))

# Recreate objects using just these pathways 
Pathways = data.frame(unique(Interactions_of_intrest$interaction_name))
colnames(Pathways) = "interaction_name"
MAEPathways =Pathways[!Pathways %in% c()]
MASPathways =Pathways[!Pathways %in% c()]
YEPathways =Pathways[!Pathways %in% c()]
YSPathways =Pathways[!Pathways %in% c()]

Aged_Ex_Go = aggregateNet2(
  Aged_Ex,
  sources.use = NULL,
  targets.use = NULL,
  signaling = NULL,
  pairLR.use = MAEPathways,
  remove.isolate = FALSE,
  thresh = 0.05,
  return.object = TRUE
)

print("Aged_Ex DONE")
Aged_Sed_Go = aggregateNet2(
  Aged_Sed,
  sources.use = NULL,
  targets.use = NULL,
  signaling = NULL,
  pairLR.use = MASPathways,
  remove.isolate = FALSE,
  thresh = 0.05,
  return.object = TRUE
)
print("Aged_Sed DONE")
Young_Ex_Go = aggregateNet2(
  Young_Ex,
  sources.use = NULL,
  targets.use = NULL,
  signaling = NULL,
  pairLR.use = YEPathways,
  remove.isolate = FALSE,
  thresh = 0.05,
  return.object = TRUE
)
print("Young_Ex DONE")
Young_Sed_Go = aggregateNet2(
  Young_Sed,
  sources.use = NULL,
  targets.use = NULL,
  signaling = NULL,
  pairLR.use = YSPathways,
  remove.isolate = FALSE,
  thresh = 0.05,
  return.object = TRUE
)
print("YS DONE")


library(dplyr)  # Added: the summarise() pipelines below run before the later library(dplyr) call
df_interactions <- subsetCommunication(Aged_Ex_Go)

# 2. Count the number of interactions per condition
interaction_counts <- df_interactions %>%

  summarise(
    n_interactions = n()
  )

print(interaction_counts)
df_interactions <- subsetCommunication(Aged_Sed_Go)

# 2. Count the number of interactions per condition
interaction_counts <- df_interactions %>%

  summarise(
    n_interactions = n()
  )

print(interaction_counts)
df_interactions <- subsetCommunication(Young_Sed_Go)

# 2. Count the number of interactions per condition
interaction_counts <- df_interactions %>%

  summarise(
    n_interactions = n()
  )

print(interaction_counts)

df_interactions

object.list <- list(Aged_Ex = Aged_Ex_Go, Aged_Sed = Aged_Sed_Go,Young_Ex = Young_Ex_Go,Young_Sed = Young_Sed_Go)
cellchatGO <- mergeCellChat(object.list, add.names = names(object.list))

# Check which interactions 
netVisual_bubble(Aged_Ex_Go, remove.isolate = TRUE, targets.use = "M1", font.size = 20, thresh = 0.05)
netVisual_bubble(Aged_Sed_Go, remove.isolate = TRUE, targets.use = "M1", font.size = 20, thresh = 0.05)

genes = c(  ligands,
receptors,
Complex_Cofactor)
getwd()
write.csv(genes,"Immune_genes_cellChat.csv")
write.csv(ligands,"Immune_ligands_cellChat.csv")
write.csv(receptors,"Immune_receptors_cellChat.csv")

genes =read.csv("Immune_genes_cellChat.csv")
genes$x

# Create heatmaps of these subsets 
PosGo <- netVisual_heatmap2(cellchatGO,comparison = c("Aged_Sed", "Aged_Ex"), measure = "weight", remove.isolate = TRUE, minCol = -.1, maxCol =  0.1,
cluster.rows = TRUE,
cluster.cols = TRUE,
                         bar_row_colour = RNA_pallate_bar_row, bar_col_colour = RNA_pallate_bar_col,color.use = RNA_pallate)

PosGo

# Create heatmaps of these subsets Ex
gg1 = netVisual_heatmap2(cellchatGO,comparison = c("Aged_Sed", "Aged_Ex"), measure = "weight", remove.isolate = TRUE, minCol = -.1, maxCol =  0.1,
cluster.rows = FALSE, row_order = row_order(PosGo),
cluster.cols = FALSE,  col_order = column_order(PosGo),
                         bar_row_colour = RNA_pallate_bar_row, bar_col_colour = RNA_pallate_bar_col,color.use = RNA_pallate)
#> Do heatmap based on a merged object
gg1
getwd()
pdf("AGE_Ex_GO_Immune_Weight2.pdf")
gg1
dev.off()

# Create heatmaps of these subsets Age
gg1 = netVisual_heatmap2(cellchatGO,comparison = c("Young_Sed", "Aged_Sed"), measure = "weight", remove.isolate = TRUE, minCol = -.1, maxCol =  0.1,
cluster.rows = FALSE, row_order = row_order(PosGo),
cluster.cols = FALSE,  col_order = column_order(PosGo),
                         bar_row_colour = RNA_pallate_bar_row, bar_col_colour = RNA_pallate_bar_col,color.use = RNA_pallate)

gg1
pdf("AGE_GO_Immune_Weight2.pdf")
gg1
dev.off()

# Get full csv 

object.list <- list(MAE = Aged_Ex_Go, MAS = Aged_Sed_Go)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "MAE"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0, thresh.fc = 0, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
Ex_Up_Age <- subsetCommunication(cellchat, net = net, datasets = "MAS",pairLR.use = Pathways)
Ex_Up_Ex <- subsetCommunication(cellchat, net = net, datasets = "MAE",pairLR.use = Pathways)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
#net.down <- subsetCommunication(cellchat, net = net, datasets = "YS",ligand.logFC = -0.1, receptor.logFC = -0.1)
#gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
#gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
#netVisual_chord_gene(object.list[[2]], sources.use = 8, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#pdf("DownRegulated_Age_FibroblastSource.pdf")
#netVisual_chord_gene(object.list[[1]], sources.use = 8, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in Aged Sed vs Young Sed"))
#dev.off()

#dev.off()
#pdf("UpRegulated_Age_Immune.pdf")
#netVisual_chord_gene(object.list[[2]], slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in Aged Sed vs Young Sed"))
#dev.off()

# Get full csv 

object.list <- list(YS = Young_Sed_Go, MAS = Aged_Sed_Go)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "MAS"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.0, thresh.fc = 0.0, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
Age_Up_Age <- subsetCommunication(cellchat, net = net, datasets = "MAS",pairLR.use = Pathways)
Age_Up_Young <- subsetCommunication(cellchat, net = net, datasets = "YS",pairLR.use = Pathways)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
#net.down <- subsetCommunication(cellchat, net = net, datasets = "YS",ligand.logFC = -0.1, receptor.logFC = -0.1)
#gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
#gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
#netVisual_chord_gene(object.list[[2]], sources.use = 8, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#pdf("DownRegulated_Age_FibroblastSource.pdf")
#netVisual_chord_gene(object.list[[1]], sources.use = 8, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in Aged Sed vs Young Sed"))
#dev.off()

#dev.off()
#pdf("UpRegulated_Age_Immune.pdf")
#netVisual_chord_gene(object.list[[2]], slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in Aged Sed vs Young Sed"))
#dev.off()
# Immune_Age = net.up
# write.csv(Immune_Age, publication_path("LEGACY_ATAC_RESULTS_ROOT", "CellChat_GO_Immune_Age_SplitMacrophage.csv"))







library(dplyr)
Ex_Up_Age$ID = paste( Ex_Up_Age$source, Ex_Up_Age$target, Ex_Up_Age$interaction_name)
Ex_Up_Ex$ID = paste( Ex_Up_Ex$source, Ex_Up_Ex$target, Ex_Up_Ex$interaction_name)
Age_Up_Age$ID = paste( Age_Up_Age$source, Age_Up_Age$target, Age_Up_Age$interaction_name)
Age_Up_Young$ID = paste( Age_Up_Young$source, Age_Up_Young$target, Age_Up_Young$interaction_name)
All = merge(Ex_Up_Age,Ex_Up_Ex, by = "ID", suffixes  = c(".Ex_Up_Age",".Ex_Up_Ex"), all = TRUE)
All = merge(All,Age_Up_Age, by = "ID", suffixes = c("",".Age_Up_Age"), all = TRUE)
All = merge(All,Age_Up_Young, by = "ID", suffixes = c("",".Age_Up_Young"), all = TRUE)
All <- All %>%
  mutate(
    # 1. Update source
    source = coalesce(source, source.Age_Up_Young, source.Ex_Up_Age, source.Ex_Up_Ex),
    
    # 2. Update target
    target = coalesce(target, target.Age_Up_Young, target.Ex_Up_Age, target.Ex_Up_Ex),
    
    # 3. Update ligand
    ligand = coalesce(ligand, ligand.Age_Up_Young, ligand.Ex_Up_Age, ligand.Ex_Up_Ex),
    
    # 4. Update receptor
    receptor = coalesce(receptor, receptor.Age_Up_Young, receptor.Ex_Up_Age, receptor.Ex_Up_Ex)
  )
write.csv(All, publication_path("LEGACY_ATAC_RESULTS_ROOT", "CellChat_GO_Immune_SplitMacrophage_AllData.csv"))

write.csv(All, publication_path("LEGACY_ATAC_RESULTS_ROOT", "CellChat_GO_Immune_SplitMacrophage_AllData.csv"))

Pathways
DT::datatable( CellChatDB.use$interaction , filter = 'top')

DT::datatable( Aged_Sed@netP , filter = 'top')

options(warn=-1)
#Check some specific types of interactions 
PathsFilt = c("IL6")
for (test in PathsFilt){
    

    setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2024-02-22.CellChat_Clean/pathways"))
    
    print(test)
    Aged_Ex_Path = aggregateNet2(
      Aged_Ex,
      sources.use = NULL,
      targets.use = NULL,
      signaling = test,
      pairLR.use = NULL,
      thresh = 0.05,
      remove.isolate = FALSE,
      return.object = TRUE
    )

    print("Aged_Ex DONE")
    Aged_Sed_Path = aggregateNet2(
      Aged_Sed,
      sources.use = NULL,
      targets.use = NULL,
      signaling = test,
      pairLR.use = NULL,
      thresh = 0.05,
         remove.isolate = FALSE,
      return.object = TRUE
    )
    print("Aged_Sed DONE")
    Young_Ex_Path = aggregateNet2(
      Young_Ex,
      sources.use = NULL,
      targets.use = NULL,
      signaling = test,
      pairLR.use = NULL,
      thresh = 0.05,
         remove.isolate = FALSE,
      return.object = TRUE
    )
    print("Young_Ex DONE")
    Young_Sed_Path = aggregateNet2(
      Young_Sed,
      sources.use = NULL,
      targets.use = NULL,
      signaling = test,
      pairLR.use = NULL,
      thresh = 0.05,
         remove.isolate = FALSE,
      return.object = TRUE
    )
    print("Young_Sed DONE")

    object.list <- list(Aged_Ex = Aged_Ex_Path, Aged_Sed = Aged_Sed_Path,Young_Ex = Young_Ex_Path,Young_Sed = Young_Sed_Path)
    
    cellchat_Path <- mergeCellChat(object.list, add.names = names(object.list))
    
    enriched_interactions = c(Aged_Sed = nrow(extractEnrichedLR(
      Aged_Sed_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    )),


     Aged_Ex = nrow(extractEnrichedLR(
      Aged_Ex_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    )),

     Young_Sed = nrow(extractEnrichedLR(
      Young_Sed_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    )),

     Young_Ex= nrow(extractEnrichedLR(
      Young_Ex_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    ))

    )
    enriched_interactionsframe =data.frame(enriched_interactions)

    enriched_interactionsframe$cond = rownames(enriched_interactionsframe)

    #enriched_interactionsframe

  

    
    Age = netVisual_heatmap2(cellchat_Path,comparison =  c("Young_Sed", "Aged_Sed"), measure = "weight", signaling = test, title.name = "Age Effect")
    Ex = netVisual_heatmap2(cellchat_Path,comparison =  c("Aged_Sed", "Aged_Ex"), measure = "weight", signaling = test, title.name = "Ex Effect")
    YoungEx = netVisual_heatmap2(cellchat_Path,comparison =  c("Young_Sed", "Young_Ex"), measure = "weight", signaling = test, title.name = "Young Ex Effect")
    MinDiff = min(c(min(as.vector(Age@matrix),na.rm = TRUE),min(as.vector(Ex@matrix),na.rm = TRUE),min(as.vector(YoungEx@matrix),na.rm = TRUE)))
    MaxDiff = max(c(max(as.vector(Age@matrix),na.rm = TRUE),max(as.vector(Ex@matrix),na.rm = TRUE),max(as.vector(YoungEx@matrix),na.rm = TRUE)))

    
    print("SetDif")
    YoungSed = netVisual_heatmap2(Young_Sed_Path, measure = "weight", signaling = test, title.name = "Young Sed")
    YoungEx = netVisual_heatmap2(Young_Ex_Path, measure = "weight", signaling = test, title.name = "Young Sed")
    AgeSed = netVisual_heatmap2(Aged_Sed_Path, measure = "weight", signaling = test, title.name = "Age Sed")
    AgeEx = netVisual_heatmap2(Aged_Ex_Path, measure = "weight", signaling = test, title.name = "Age Ex")



    Max = max(c(max(as.vector(YoungSed@matrix),na.rm = TRUE),max(as.vector(AgeSed@matrix) ,na.rm = TRUE),max(as.vector(AgeEx@matrix),na.rm = TRUE),max(as.vector(YoungEx@matrix),na.rm = TRUE)))

       MinDiff   
       MaxDiff
       Max       
    
    print("SetMinMax")
    Age =  netVisual_heatmap2(cellchat_Path,comparison =  c("Young_Sed", "Aged_Sed"), measure = "weight", signaling = test, title.name = paste(test,"Age Effect"), minCol = MinDiff, maxCol = MaxDiff)

    Ex = netVisual_heatmap2(cellchat_Path,comparison =  c("Aged_Sed", "Aged_Ex"), measure = "weight", signaling = test, title.name =  paste(test,"Age Ex Effect"), minCol = MinDiff, maxCol = MaxDiff)
    YoungEx = netVisual_heatmap2(cellchat_Path,comparison =  c("Young_Sed", "Young_Ex"), measure = "weight", signaling = test,title.name =  paste(test,"Young Ex Effect"), minCol = MinDiff, maxCol = MaxDiff)
    
    
    YoungSed = netVisual_heatmap2(Young_Sed_Path, measure = "weight", signaling = test, title.name =  paste(test,"Young Sed"), minCol = 0, maxCol = Max, color.heatmap = c("#f7f7f7", "#b2182b"))
    YoungExSing = netVisual_heatmap2(Young_Ex_Path, measure = "weight", signaling = test, title.name =  paste(test,"Young Ex"), minCol = 0, maxCol = Max, color.heatmap = c("#f7f7f7", "#b2182b"))
    AgeSed = netVisual_heatmap2(Aged_Sed_Path, measure = "weight", signaling = test, title.name =  paste(test,"Age Sed"), minCol = 0, maxCol = Max, color.heatmap = c("#f7f7f7", "#b2182b"))
    AgeEx = netVisual_heatmap2(Aged_Ex_Path, measure = "weight", signaling = test, title.name =  paste(test,"Age Ex"), minCol = 0, maxCol = Max, color.heatmap = c("#f7f7f7", "#b2182b"))

    nGenes = ggplot(enriched_interactionsframe, aes(x=cond, y=enriched_interactions)) +geom_bar(stat="identity")
    plots = list(Age,YoungEx,Ex,nGenes,YoungSed,YoungExSing,AgeSed,AgeEx)

    
    print("Draw")
    pdf(paste(test,"combined_heatmaps.pdf"), width = 15, height = 10)
    pushViewport(viewport(layout = grid.layout(3, 4))) # Create a layout with 2 row and 3 columns

    # Add the first heatmap to the left column
    pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
    draw(Age)
    popViewport() # Exit the left column viewport

    
    pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 1))
    draw(YoungEx)
    popViewport() # 
    
    
    
    
    # Add the second heatmap to the right column
    pushViewport(viewport(layout.pos.col = 3, layout.pos.row = 1))
    draw(Ex)
    popViewport() # Exit the right column viewport

    print(ggplot(enriched_interactionsframe, aes(x=cond, y=enriched_interactions)) +geom_bar(stat="identity"), vp = viewport(layout.pos.row = 1, layout.pos.col = 4))
  
    
    pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2))
    draw(YoungSed)
    popViewport()
    
    pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 2))
    draw(YoungExSing)
    popViewport()

    pushViewport(viewport(layout.pos.col = 3, layout.pos.row = 2))
    draw(AgeSed)
    popViewport()

    pushViewport(viewport(layout.pos.col = 4, layout.pos.row = 2))
    draw(AgeEx)
    popViewport()

    popViewport() # Exit the main viewport
    
    
    pushViewport(viewport(layout = grid.layout(3, 1)))
    dev.off() # Close the PDF device
    
}


#Check some specific types of interactions 
PathsFilt = c("COLLAGEN","LAMININ", "TGFb")

options(warn=-1)
#Check some specific types of interactions 
PathsFilt = c("COLLAGEN","LAMININ", "TGFb")
for (test in PathsFilt){
    

    setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2024-02-22.CellChat_Clean/pathways2"))
    
    print(test)
    Aged_Ex_Path = aggregateNet2(
      Aged_Ex,
      sources.use = NULL,
      targets.use = NULL,
      signaling = test,
      pairLR.use = NULL,
      thresh = 0.05,
      remove.isolate = FALSE,
      return.object = TRUE
    )

    print("Aged_Ex DONE")
    Aged_Sed_Path = aggregateNet2(
      Aged_Sed,
      sources.use = NULL,
      targets.use = NULL,
      signaling = test,
      pairLR.use = NULL,
      thresh = 0.05,
         remove.isolate = FALSE,
      return.object = TRUE
    )
    print("Aged_Sed DONE")
    Young_Ex_Path = aggregateNet2(
      Young_Ex,
      sources.use = NULL,
      targets.use = NULL,
      signaling = test,
      pairLR.use = NULL,
      thresh = 0.05,
         remove.isolate = FALSE,
      return.object = TRUE
    )
    print("Young_Ex DONE")
    Young_Sed_Path = aggregateNet2(
      Young_Sed,
      sources.use = NULL,
      targets.use = NULL,
      signaling = test,
      pairLR.use = NULL,
      thresh = 0.05,
         remove.isolate = FALSE,
      return.object = TRUE
    )
    print("Young_Sed DONE")

    object.list <- list(Aged_Ex = Aged_Ex_Path, Aged_Sed = Aged_Sed_Path,Young_Ex = Young_Ex_Path,Young_Sed = Young_Sed_Path)
    
    cellchat_Path <- mergeCellChat(object.list, add.names = names(object.list))
    
    enriched_interactions = c(Aged_Sed = nrow(extractEnrichedLR(
      Aged_Sed_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    )),


     Aged_Ex = nrow(extractEnrichedLR(
      Aged_Ex_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    )),

     Young_Sed = nrow(extractEnrichedLR(
      Young_Sed_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    )),

     Young_Ex= nrow(extractEnrichedLR(
      Young_Ex_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    ))

    )
    enriched_interactionsframe =data.frame(enriched_interactions)

    enriched_interactionsframe$cond = rownames(enriched_interactionsframe)

    #enriched_interactionsframe

  

    
    Age = netVisual_heatmap2(cellchat_Path,comparison =  c("Young_Sed", "Aged_Sed"), measure = "weight", signaling = test, title.name = "Age Effect")
    Ex = netVisual_heatmap2(cellchat_Path,comparison =  c("Aged_Sed", "Aged_Ex"), measure = "weight", signaling = test, title.name = "Ex Effect")
    YoungEx = netVisual_heatmap2(cellchat_Path,comparison =  c("Young_Sed", "Young_Ex"), measure = "weight", signaling = test, title.name = "Young Ex Effect")
    MinDiff = min(c(min(as.vector(Age@matrix),na.rm = TRUE),min(as.vector(Ex@matrix),na.rm = TRUE),min(as.vector(YoungEx@matrix),na.rm = TRUE)))
    MaxDiff = max(c(max(as.vector(Age@matrix),na.rm = TRUE),max(as.vector(Ex@matrix),na.rm = TRUE),max(as.vector(YoungEx@matrix),na.rm = TRUE)))

    
    print("SetDif")
    YoungSed = netVisual_heatmap2(Young_Sed_Path, measure = "weight", signaling = test, title.name = "Young Sed")
    YoungEx = netVisual_heatmap2(Young_Ex_Path, measure = "weight", signaling = test, title.name = "Young Sed")
    AgeSed = netVisual_heatmap2(Aged_Sed_Path, measure = "weight", signaling = test, title.name = "Age Sed")
    AgeEx = netVisual_heatmap2(Aged_Ex_Path, measure = "weight", signaling = test, title.name = "Age Ex")



    Max = max(c(max(as.vector(YoungSed@matrix),na.rm = TRUE),max(as.vector(AgeSed@matrix) ,na.rm = TRUE),max(as.vector(AgeEx@matrix),na.rm = TRUE),max(as.vector(YoungEx@matrix),na.rm = TRUE)))

       MinDiff   
       MaxDiff
       Max       
    
    print("SetMinMax")
    Age =  netVisual_heatmap2(cellchat_Path,comparison =  c("Young_Sed", "Aged_Sed"), measure = "weight", signaling = test, title.name = paste(test,"Age Effect"), minCol = MinDiff, maxCol = MaxDiff)

    Ex = netVisual_heatmap2(cellchat_Path,comparison =  c("Aged_Sed", "Aged_Ex"), measure = "weight", signaling = test, title.name =  paste(test,"Age Ex Effect"), minCol = MinDiff, maxCol = MaxDiff)
    YoungEx = netVisual_heatmap2(cellchat_Path,comparison =  c("Young_Sed", "Young_Ex"), measure = "weight", signaling = test,title.name =  paste(test,"Young Ex Effect"), minCol = MinDiff, maxCol = MaxDiff)
    
    
    YoungSed = netVisual_heatmap2(Young_Sed_Path, measure = "weight", signaling = test, title.name =  paste(test,"Young Sed"), minCol = 0, maxCol = Max, color.heatmap = c("#f7f7f7", "#b2182b"))
    YoungExSing = netVisual_heatmap2(Young_Ex_Path, measure = "weight", signaling = test, title.name =  paste(test,"Young Ex"), minCol = 0, maxCol = Max, color.heatmap = c("#f7f7f7", "#b2182b"))
    AgeSed = netVisual_heatmap2(Aged_Sed_Path, measure = "weight", signaling = test, title.name =  paste(test,"Age Sed"), minCol = 0, maxCol = Max, color.heatmap = c("#f7f7f7", "#b2182b"))
    AgeEx = netVisual_heatmap2(Aged_Ex_Path, measure = "weight", signaling = test, title.name =  paste(test,"Age Ex"), minCol = 0, maxCol = Max, color.heatmap = c("#f7f7f7", "#b2182b"))

    nGenes = ggplot(enriched_interactionsframe, aes(x=cond, y=enriched_interactions)) +geom_bar(stat="identity")
    plots = list(Age,YoungEx,Ex,nGenes,YoungSed,YoungExSing,AgeSed,AgeEx)

    
    print("Draw")
    pdf(paste(test,"combined_heatmaps.pdf"), width = 15, height = 10)
    pushViewport(viewport(layout = grid.layout(3, 4))) # Create a layout with 2 row and 3 columns

    # Add the first heatmap to the left column
    pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
    draw(Age)
    popViewport() # Exit the left column viewport

    
    pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 1))
    draw(YoungEx)
    popViewport() # 
    
    
    
    
    # Add the second heatmap to the right column
    pushViewport(viewport(layout.pos.col = 3, layout.pos.row = 1))
    draw(Ex)
    popViewport() # Exit the right column viewport

    print(ggplot(enriched_interactionsframe, aes(x=cond, y=enriched_interactions)) +geom_bar(stat="identity"), vp = viewport(layout.pos.row = 1, layout.pos.col = 4))
  
    
    pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2))
    draw(YoungSed)
    popViewport()
    
    pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 2))
    draw(YoungExSing)
    popViewport()

    pushViewport(viewport(layout.pos.col = 3, layout.pos.row = 2))
    draw(AgeSed)
    popViewport()

    pushViewport(viewport(layout.pos.col = 4, layout.pos.row = 2))
    draw(AgeEx)
    popViewport()

    popViewport() # Exit the main viewport
    
    
    pushViewport(viewport(layout = grid.layout(3, 1)))
    dev.off() # Close the PDF device
    
}




# Do the same pathways but only create the specifci plots of interest 

PathsFilt = c("COLLAGEN","LAMININ", "TGFb")

options(warn=-1)
PathsFilt = c("COLLAGEN","LAMININ", "TGFb")
for (test in PathsFilt){
    

    setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2024-02-22.CellChat_Clean/pathways2"))
    
    print(test)
    Aged_Ex_Path = aggregateNet2(
      Aged_Ex,
      sources.use = NULL,
      targets.use = NULL,
      signaling = test,
      pairLR.use = NULL,
      thresh = 0.05,
      remove.isolate = FALSE,
      return.object = TRUE
    )

    print("Aged_Ex DONE")
    Aged_Sed_Path = aggregateNet2(
      Aged_Sed,
      sources.use = NULL,
      targets.use = NULL,
      signaling = test,
      pairLR.use = NULL,
      thresh = 0.05,
         remove.isolate = FALSE,
      return.object = TRUE
    )
    print("Aged_Sed DONE")
    Young_Ex_Path = aggregateNet2(
      Young_Ex,
      sources.use = NULL,
      targets.use = NULL,
      signaling = test,
      pairLR.use = NULL,
      thresh = 0.05,
         remove.isolate = FALSE,
      return.object = TRUE
    )
    print("Young_Ex DONE")
    Young_Sed_Path = aggregateNet2(
      Young_Sed,
      sources.use = NULL,
      targets.use = NULL,
      signaling = test,
      pairLR.use = NULL,
      thresh = 0.05,
         remove.isolate = FALSE,
      return.object = TRUE
    )
    print("Young_Sed DONE")

    object.list <- list(Aged_Ex = Aged_Ex_Path, Aged_Sed = Aged_Sed_Path,Young_Ex = Young_Ex_Path,Young_Sed = Young_Sed_Path)
    
    cellchat_Path <- mergeCellChat(object.list, add.names = names(object.list))
    
    enriched_interactions = c(Aged_Sed = nrow(extractEnrichedLR(
      Aged_Sed_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    )),


     Aged_Ex = nrow(extractEnrichedLR(
      Aged_Ex_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    )),

     Young_Sed = nrow(extractEnrichedLR(
      Young_Sed_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    )),

     Young_Ex= nrow(extractEnrichedLR(
      Young_Ex_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    ))

    )
    enriched_interactionsframe =data.frame(enriched_interactions)

    enriched_interactionsframe$cond = rownames(enriched_interactionsframe)

   # enriched_interactionsframe


 print("here")

    Age = netVisual_heatmap2(cellchat_Path,comparison =  c("Young_Sed", "Aged_Sed"), measure = "weight", signaling = test, title.name = "Age Effect",  cluster.rows = TRUE, cluster.cols = TRUE)
    Ex = netVisual_heatmap2(cellchat_Path,comparison =  c("Aged_Sed", "Aged_Ex"), measure = "weight", signaling = test, title.name = "Ex Effect",  cluster.rows = TRUE,cluster.cols = TRUE)
    YoungEx = netVisual_heatmap2(cellchat_Path,comparison =  c("Young_Sed", "Young_Ex"), measure = "weight", signaling = test, title.name = "Young Ex Effect")
    MinDiff = min(c(min(as.vector(Age@matrix),na.rm = TRUE),min(as.vector(Ex@matrix),na.rm = TRUE),min(as.vector(YoungEx@matrix),na.rm = TRUE)))
    MaxDiff = max(c(max(as.vector(Age@matrix),na.rm = TRUE),max(as.vector(Ex@matrix),na.rm = TRUE),max(as.vector(YoungEx@matrix),na.rm = TRUE)))

    
    print("SetDif")
   

    Max = max(c(max(as.vector(Age@matrix),na.rm = TRUE),max(as.vector(Ex@matrix) ,na.rm = TRUE)))

       MinDiff   
       MaxDiff
       Max       
    
    print("SetMinMax")
    Age =  netVisual_heatmap2(cellchat_Path,comparison =  c("Young_Sed", "Aged_Sed"), measure = "weight", signaling = test, title.name = paste(test,"Age Effect"), minCol = MinDiff, maxCol = MaxDiff,
                       cluster.rows = FALSE, row_order = row_order(Age),
                        cluster.cols = FALSE,  col_order = column_order(Age),
                         bar_row_colour = RNA_pallate_bar_row, bar_col_colour = RNA_pallate_bar_col,color.use = RNA_pallate)

    Ex = netVisual_heatmap2(cellchat_Path,comparison =  c("Aged_Sed", "Aged_Ex"), measure = "weight", signaling = test, title.name =  paste(test,"Age Ex Effect"), minCol = MinDiff, maxCol = MaxDiff,
                                   cluster.rows = FALSE, row_order = row_order(Age),
                                    cluster.cols = FALSE,  col_order = column_order(Age),
                         bar_row_colour = RNA_pallate_bar_row, bar_col_colour = RNA_pallate_bar_col,color.use = RNA_pallate)
     

    
    print("Draw")
    pdf(paste(test,"Age_heatmaps2.pdf"), width = 10, height = 10)
   
    draw(Age)
  
    dev.off() 
    
    
    pdf(paste(test,"AgeEx_heatmaps2.pdf"), width = 10, height = 10)
   
    draw(Ex)
  
    dev.off() # Close the PDF device
    
}


# Recreate objects using just these pathways 

Aged_Ex_Col = aggregateNet2(
  Aged_Ex,
  sources.use = NULL,
  targets.use = NULL,
  signaling = 'COLLAGEN',
  pairLR.use = NULL,
  remove.isolate = FALSE,
  thresh = 0.05,
  return.object = TRUE
)

print("Aged_Ex DONE")
Aged_Sed_Col = aggregateNet2(
  Aged_Sed,
  sources.use = NULL,
  targets.use = NULL,
  signaling = 'COLLAGEN',
  pairLR.use = NULL,
  remove.isolate = FALSE,
  thresh = 0.05,
  return.object = TRUE
)
print("Aged_Sed DONE")
Young_Ex_Col = aggregateNet2(
  Young_Ex,
  sources.use = NULL,
  targets.use = NULL,
  signaling = 'COLLAGEN',
  pairLR.use = NULL,
  remove.isolate = FALSE,
  thresh = 0.05,
  return.object = TRUE
)
print("Young_Ex DONE")
Young_Sed_Col = aggregateNet2(
  Young_Sed,
  sources.use = NULL,
  targets.use = NULL,
  signaling = "COLLAGEN",
  pairLR.use = NULL,
  remove.isolate = FALSE,
  thresh = 0.05,
  return.object = TRUE
)
print("YS DONE")

# Get full csv 

object.list <- list(MAE = Aged_Ex_Col, MAS = Aged_Sed_Col)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "MAE"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0, thresh.fc = 0, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
Ex_Up_Age <- subsetCommunication(cellchat, net = net, datasets = "MAS",signaling = "COLLAGEN" )
Ex_Up_Ex <- subsetCommunication(cellchat, net = net, datasets = "MAE",signaling = "COLLAGEN" )

# Get full csv 

object.list <- list(YS = Young_Sed_Col, MAS = Aged_Sed_Col)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "MAS"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.0, thresh.fc = 0.0, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
Age_Up_Age <- subsetCommunication(cellchat, net = net, datasets = "MAS",signaling = "COLLAGEN" )
Age_Up_Young <- subsetCommunication(cellchat, net = net, datasets = "YS",signaling = "COLLAGEN" )



library(dplyr)
Ex_Up_Age$ID = paste( Ex_Up_Age$source, Ex_Up_Age$target, Ex_Up_Age$interaction_name)
Ex_Up_Ex$ID = paste( Ex_Up_Ex$source, Ex_Up_Ex$target, Ex_Up_Ex$interaction_name)
Age_Up_Age$ID = paste( Age_Up_Age$source, Age_Up_Age$target, Age_Up_Age$interaction_name)
Age_Up_Young$ID = paste( Age_Up_Young$source, Age_Up_Young$target, Age_Up_Young$interaction_name)
All = merge(Ex_Up_Age,Ex_Up_Ex, by = "ID", suffixes  = c(".Ex_Up_Age",".Ex_Up_Ex"), all = TRUE)
All = merge(All,Age_Up_Age, by = "ID", suffixes = c("",".Age_Up_Age"), all = TRUE)
All = merge(All,Age_Up_Young, by = "ID", suffixes = c("",".Age_Up_Young"), all = TRUE)
All <- All %>%
  mutate(
    # 1. Update source
    source = coalesce(source, source.Age_Up_Young, source.Ex_Up_Age, source.Ex_Up_Ex),
    
    # 2. Update target
    target = coalesce(target, target.Age_Up_Young, target.Ex_Up_Age, target.Ex_Up_Ex),
    
    # 3. Update ligand
    ligand = coalesce(ligand, ligand.Age_Up_Young, ligand.Ex_Up_Age, ligand.Ex_Up_Ex),
    
    # 4. Update receptor
    receptor = coalesce(receptor, receptor.Age_Up_Young, receptor.Ex_Up_Age, receptor.Ex_Up_Ex)
  )
write.csv(All, publication_path("LEGACY_ATAC_RESULTS_ROOT", "CellChat_GO_COLLAGEN_SplitMacrophage_AllData.csv"))









# Check the actual interactions between each celltype 
setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2024-02-22.CellChat_Clean"))
    Aged_Ex = readRDS("Aged_Ex_Filt_cellChat.rds")
    Aged_Sed = readRDS("Aged_Sed_Filt_cellChat.rds")
    Young_Ex = readRDS("Young_Ex_Filt_cellChat.rds")
    Young_Sed = readRDS("Young_Sed_Filt_cellChat.rds")

object.list <- list(Aged_Ex = Aged_Ex, Aged_Sed = Aged_Sed,Young_Ex = Young_Ex,Young_Sed = Young_Sed)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

PathsFilt = c("COLLAGEN","LAMININ","TGFb")

groupSize <- as.numeric(table(Aged_Ex@idents))
netVisual_circle(Aged_Ex@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- Aged_Ex@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

i =8
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])

i =4
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])

pathways.show <- c("COLLAGEN") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(Aged_Ex, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(Aged_Ex, signaling = pathways.show, layout = "circle")

netVisual_aggregate(Aged_Ex, signaling = pathways.show, layout = "chord")


netVisual_heatmap(Aged_Ex, signaling = "COLLAGEN", color.heatmap = "Reds")

netVisual_heatmap(Aged_Sed, signaling = "COLLAGEN", color.heatmap = "Reds")

netVisual_heatmap(Aged_Ex, signaling = "LAMININ", color.heatmap = "Reds")

netAnalysis_contribution(Aged_Ex, signaling = "COLLAGEN")

netVisual_bubble(Aged_Ex, sources.use = "M_Ex", remove.isolate = FALSE, signaling = c("COLLAGEN"))

netVisual_bubble(Aged_Sed, sources.use = "M_Ex", remove.isolate = FALSE, signaling = c("COLLAGEN"))



# :Look at IL1r1

netVisual_heatmap(Aged_Ex, signaling = "IL1", color.heatmap = "Reds")

netVisual_heatmap(Aged_Sed, signaling = "IL1", color.heatmap = "Reds") 

netVisual_heatmap(Young_Sed, signaling = "IL1", color.heatmap = "Reds") 

netVisual_heatmap(Young_Ex, signaling = "IL1", color.heatmap = "Reds") 

netVisual_heatmap(Aged_Ex, signaling = "IL1", color.heatmap = "Reds") 

netVisual_bubble(Aged_Ex, remove.isolate = FALSE, signaling = c("IL1"), font.size = 20, thresh = 1)

netVisual_heatmap(Aged_Ex, signaling = "TGFb", color.heatmap = "Reds") 

netVisual_heatmap(Aged_Ex, signaling = "COLLAGEN", color.heatmap = "Reds",measure = "weight") 

netVisual_heatmap(Aged_Ex, signaling = "COLLAGEN", color.heatmap = "Reds",measure = "count") 

netVisual_bubble(Aged_Ex, remove.isolate = FALSE,sources.use = "M_Ex", signaling = c("TGFb"), font.size = 20, thresh = 1)

enriched_interactions
enriched_interactionsframe

extractEnrichedLR(
      Aged_Sed_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    )
extractEnrichedLR(
      Aged_Ex_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    )
extractEnrichedLR(
      Young_Sed_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    )

TgfGenes = c('Tgfb1',
'Tgfb2','Acvr1','Tgfbr1',
'Tgfb3','Tgfbr2')

library(openxlsx)
library(stringr)
library(EnhancedVolcano)
setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed"))
RandomBatch=list.files(pattern = 'MAST.*csv', recursive = TRUE, include.dirs = TRUE)
RandomBatch = RandomBatch[!grepl(RandomBatch, pattern = "Global")]
RandomBatch

All = data.frame(matrix(nrow = 0,ncol = 15))
for (i in RandomBatch) {
 
  names =i
 
     
   # print(i)
    setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed"))
    table =  read.csv(i)
      name = str_remove(i, "MAST")
      name = str_remove(name, ".csv")
      condition = str_extract(name,"Aged.*|Young.*")
      cell = str_remove(name, "_Aged.*|_Young.*")
      table$cell = cell
      table$condition = condition
   
    
      All = rbind(All,table)
  }

Tgfb_DEG = subset(All,primerid %in% TgfGenes)
Tgfb_DEG = Tgfb_DEG[!grepl(Tgfb_DEG$cell, pattern = "Unk"),]



AgeEffect = subset(Tgfb_DEG,condition %in% 'Aged_Sed-Young_Sed')
ExEffect = subset(Tgfb_DEG,condition %in% 'Aged_Ex-Aged_Sed')

setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2024-02-22.CellChat_Clean/"))

library(EnhancedVolcano)

AgeEffect$name = paste(AgeEffect$cell, AgeEffect$primerid)
ExEffect$name = paste(ExEffect$cell, ExEffect$primerid)

DT::datatable(AgeEffect)
DT::datatable(ExEffect)

library(EnhancedVolcano)
# Create the volcano for the genes we are interesetned in 
# Assuming your data frame is called `data`, with columns `coef` (for fold change), `fdr`, and `name` for labels
AgeVolc = EnhancedVolcano(AgeEffect,
                lab = AgeEffect$name,          # Column with labels (gene names or identifiers)
                x = 'coef',               # Column for x-axis (fold change)
                y = 'fdr',                # Column for y-axis (FDR or p-value)
                xlab = "Effect size (coef)", 
                ylab = "-Log10 FDR",
                xlim = c(-2.1, 2.1),          # Adjust x limits as needed
                pCutoff = 0.1,            # FDR cutoff
                FCcutoff = 0.1,           # Fold change cutoff
                pointSize = 2,            # Point size (optional)
                labSize = 3,              # Label size (optional)
                col = c("grey30", "forestgreen", "royalblue", "red2"), # Colors for points
                legendPosition = "right", # Position of legend
                legendLabSize = 12,       # Size of legend labels
                title = "Volcano Plot",
                subtitle = "Effect of Age")


library(EnhancedVolcano)

# Assuming your data frame is called `data`, with columns `coef` (for fold change), `fdr`, and `name` for labels
ExVoc = EnhancedVolcano(ExEffect,
                lab = ExEffect$name,          # Column with labels (gene names or identifiers)
                x = 'coef',               # Column for x-axis (fold change)
                y = 'fdr',                # Column for y-axis (FDR or p-value)
                xlab = "Effect size (coef)", 
                ylab = "-Log10 FDR",
               xlim = c(-2.1, 2.1),          # Adjust x limits as needed
                pCutoff = 0.1,            # FDR cutoff
                FCcutoff = 0.1,           # Fold change cutoff
                pointSize = 2,            # Point size (optional)
                labSize = 3,              # Label size (optional)
                col = c("grey30", "forestgreen", "royalblue", "red2"), # Colors for points
                legendPosition = "right", # Position of legend
                legendLabSize = 12,       # Size of legend labels
                title = "Volcano Plot",
                subtitle = "Effect of Exercsie")

library(dplyr)
library(scales)
library(ggplot2)

# Combine the data and identify opposite directions
combined_data <- AgeEffect %>%
  inner_join(ExEffect, by = c("cell", "primerid"), suffix = c("_age", "_ex")) %>%
  mutate(opposite_direction = coef_age * coef_ex < 0)  # TRUE if directions are opposite

# Plot for Age Changes with asterisk overlay
Age = ggplot(combined_data, aes(y = cell, x = primerid, size = pmin(-log10(fdr_age),10), color = coef_age)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  geom_text(data = filter(combined_data, opposite_direction), 
            aes(label = "*"), 
            vjust = -0.5, color = "black", size = 4) + 
  ggtitle("Age Changes") +
  theme_minimal()

# Plot for Ex Changes with asterisk overlay
Ex = ggplot(combined_data, aes(y = cell, x = primerid, size = pmin(-log10(fdr_ex),10), color = coef_ex)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  geom_text(data = filter(combined_data, opposite_direction), 
            aes(label = "*"), 
            vjust = -0.5, color = "black", size = 4) +
  ggtitle("Ex Changes") +
  theme_minimal()
Age + Ex
ggsave("TGFB_Gene_DEG.pdf", width = 15)

(Age + Ex) / (AgeVolc +ExVoc)
ggsave("TGFB_Gene_DEG.pdf", width = 30, height = 15)

# Check for reversal of the interactions 

setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2024-02-22.CellChat_Clean"))
    Aged_Ex = readRDS("Aged_Ex_Filt_cellChat.rds")
    Aged_Sed = readRDS("Aged_Sed_Filt_cellChat.rds")
    Young_Ex = readRDS("Young_Ex_Filt_cellChat.rds")
    Young_Sed = readRDS("Young_Sed_Filt_cellChat.rds")

object.list <- list(Aged_Ex = Aged_Ex, Aged_Sed = Aged_Sed,Young_Ex = Young_Ex,Young_Sed = Young_Sed)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

object.list <- list(MAE = Aged_Ex, MAS = Aged_Sed)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "MAE"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0, thresh.fc = 0, thresh.p = 0.1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS

# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.downEx <- subsetCommunication(cellchat, net = net, datasets = "MAS",ligand.logFC = NULL, receptor.logFC = NULL)
net.downEx = net.downEx[ net.downEx$receptor.logFC < 0 |  net.downEx$ligand.logFC < 0,]



object.list <- list(YS = Young_Sed, MAS = Aged_Sed)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "MAS"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0, thresh.fc =0, thresh.p = 0.1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.upAge <- subsetCommunication(cellchat, net = net, datasets = "MAS",ligand.logFC = NULL, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS

net.upAge = net.upAge[ net.upAge$receptor.logFC > 0 |  net.upAge$ligand.logFC > 0,]



net.upAge$all = paste(net.upAge$interaction_name,net.upAge$source,net.upAge$target)
net.downEx$all = paste(net.downEx$interaction_name,net.downEx$source,net.downEx$target)

reversed = data.frame(all = net.upAge$all[net.upAge$all %in% net.downEx$all], Path  = net.upAge$pathway_name[net.upAge$all %in% net.downEx$all], interaction_name = net.upAge$interaction_name[net.upAge$all %in% net.downEx$all], target = net.upAge$target[net.upAge$all %in% net.downEx$all], source = net.upAge$source[net.upAge$all %in% net.downEx$all])



plot_Rev = net.upAge[ net.upAge$all %in% reversed$all,]

reversed = reversed[ reversed$Path == "TGFb",] 
reversed = na.omit(reversed)
reversed_df = data.frame(reversed$interaction_name)


table(net.upAge[ net.upAge$all %in% reversed$all,]$receptor.logFC > 0 )
table(net.upAge[ net.upAge$all %in% reversed$all,]$ligand.logFC > 0 )

table(net.downEx[ net.downEx$all %in% reversed$all,]$receptor.logFC < 0 )
table(net.downEx[ net.downEx$all %in% reversed$all,]$ligand.logFC < 0 )

filt = plot_Rev[ plot_Rev$pathway_name == "TGFb",]
filt$interaction_name = unfactor(filt$interaction_name)
filt$cells = paste(filt$source, filt$target)

DT::datatable(na.omit(plot_Rev), filter = 'top')



netVisual_chord_gene(object.list[[1]], slot.name = 'net', net = plot_Rev, pairLR.use = reversed_df,lab.cex = 0.8,legend.pos.x = 0.5, legend.pos.y = 0.5,
                     small.gap = 3.5, title.name = paste0("Down-regulated signaling in Aged Sed vs Young Sed"), show.legend = T)





netVisual_chord_gene(object.list[[1]], slot.name = 'net', net = plot_Rev, pairLR.use = reversed_df,lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in Aged Sed vs Young Sed"),show.legend = T)

getwd()
dev.off()
pdf("Reversed_Tgfb_Legend2.pdf")
netVisual_chord_gene(object.list[[1]], slot.name = 'net', net = plot_Rev, pairLR.use = reversed_df,lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in Aged Sed vs Young Sed"),show.legend = T)
dev.off()


pdf("modify.pdf")
netVisual_chord_gene(object.list[[2]], slot.name = 'net', net = plot_Rev, pairLR.use = reversed_df,lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in Aged Sed vs Young Sed"),show.legend = T)
dev.off()









Tgfb = plot_Rev[ plot_Rev$pathway_name == "TGFb",]
#Tgfb= na.omit(Tgfb)
na.omit(net.downEx[ net.downEx$pathway_name == "TGFB",])
head(na.omit(Tgfb))
Tgfb = Tgfb[,c(1:9,21)]

Tgfb = Tgfb[!is.na(Tgfb$source),]
dim(Tgfb)
Joint = merge(Tgfb,net.downEx, by="all", suffixes = c("","_Ex"))
#Joint = na.omit(Joint)
Joint = merge(Joint, net.upAge , by="all", suffixes = c("","_Age"))
colnames(Joint)
dim(Joint)
#Joint = Joint[,c('interaction_name','sourcenames','targetnames','probnames','ligandnames','receptornames','ligand.logFC','ligand.pvalues',
          #       'receptor.logFC','receptor.pvalues','ligand.logFC_Age','ligand.pvalues_Age','receptor.logFC_Age','receptor.pvalues_Age'),]
dim(Joint)

DT::datatable(Joint)

Joint$Ligand_Exname = ""
Joint$Ligand_Exname = ifelse(Joint$ligand.logFC > 0 & Joint$ligand.pvalues < 0.05, "Ex Increase",Joint$Ligand_Exname)
Joint$Ligand_Exname = ifelse(Joint$ligand.logFC < 0  & Joint$ligand.pvalues < 0.05, "Ex Decrease",Joint$Ligand_Exname)

Joint$Ligand_Agename = ""
Joint$Ligand_Agename = ifelse(Joint$ligand.logFC_Age > 0  & Joint$ligand.pvalues_Age < 0.05, "Age Increase",Joint$Ligand_Agename)
Joint$Ligand_Agename = ifelse(Joint$ligand.logFC_Age < 0  & Joint$ligand.pvalues_Age < 0.05, "Age Decrease",Joint$Ligand_Agename)

Joint$Receptor_Exname = ""
Joint$Receptor_Exname = ifelse(Joint$receptor.logFC > 0  & Joint$receptor.pvalues < 0.05, "Ex Increase",Joint$Receptor_Exname)
Joint$Receptor_Exname = ifelse(Joint$receptor.logFC < 0& Joint$receptor.pvalues < 0.05, "Ex Decrease",Joint$Receptor_Exname)

Joint$Receptor_Agename = ""
Joint$Receptor_Agename = ifelse(Joint$receptor.logFC_Age > 0 & Joint$receptor.pvalues_Age < 0.05, "Age Increase",Joint$Receptor_Agename)
Joint$Receptor_Agename = ifelse(Joint$receptor.logFC_Age < 0 & Joint$receptor.pvalues_Age < 0.05, "Age Decrease",Joint$Receptor_Agename)



DT::datatable(Joint[,c(1:6,23,27,43,47)], filter = 'top')

CellChatDB <- CellChatDB.mouse 
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB

setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2024-02-22.CellChat_Clean"))
    Aged_Ex = readRDS("Aged_Ex_Filt_cellChat.rds")
    Aged_Sed = readRDS("Aged_Sed_Filt_cellChat.rds")
    Young_Ex = readRDS("Young_Ex_Filt_cellChat.rds")
    Young_Sed = readRDS("Young_Sed_Filt_cellChat.rds")

object.list <- list(Aged_Ex = Aged_Ex, Aged_Sed = Aged_Sed,Young_Ex = Young_Ex,Young_Sed = Young_Sed)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

options(warn=-1)
PathsFilt = c("IGF")
for (test in PathsFilt){
    

    setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2024-02-22.CellChat_Clean/pathways"))
    
    print(test)
    Aged_Ex_Path = aggregateNet2(
      Aged_Ex,
      sources.use = NULL,
      targets.use = NULL,
      signaling = test,
      pairLR.use = NULL,
      thresh = 0.05,
      remove.isolate = FALSE,
      return.object = TRUE
    )

    print("Aged_Ex DONE")
    Aged_Sed_Path = aggregateNet2(
      Aged_Sed,
      sources.use = NULL,
      targets.use = NULL,
      signaling = test,
      pairLR.use = NULL,
      thresh = 0.05,
         remove.isolate = FALSE,
      return.object = TRUE
    )
    print("Aged_Sed DONE")
    Young_Ex_Path = aggregateNet2(
      Young_Ex,
      sources.use = NULL,
      targets.use = NULL,
      signaling = test,
      pairLR.use = NULL,
      thresh = 0.05,
         remove.isolate = FALSE,
      return.object = TRUE
    )
    print("Young_Ex DONE")
    Young_Sed_Path = aggregateNet2(
      Young_Sed,
      sources.use = NULL,
      targets.use = NULL,
      signaling = test,
      pairLR.use = NULL,
      thresh = 0.05,
         remove.isolate = FALSE,
      return.object = TRUE
    )
    print("Young_Sed DONE")

    object.list <- list(Aged_Ex = Aged_Ex_Path, Aged_Sed = Aged_Sed_Path,Young_Ex = Young_Ex_Path,Young_Sed = Young_Sed_Path)
    
    cellchat_Path <- mergeCellChat(object.list, add.names = names(object.list))
    
    enriched_interactions = c(Aged_Sed = nrow(extractEnrichedLR(
      Aged_Sed_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    )),


     Aged_Ex = nrow(extractEnrichedLR(
      Aged_Ex_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    )),

     Young_Sed = nrow(extractEnrichedLR(
      Young_Sed_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    )),

     Young_Ex= nrow(extractEnrichedLR(
      Young_Ex_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    ))

    )
    enriched_interactionsframe =data.frame(enriched_interactions)

    enriched_interactionsframe$cond = rownames(enriched_interactionsframe)

    enriched_interactionsframe

  

    
    Age = netVisual_heatmap2(cellchat_Path,comparison =  c("Young_Sed", "Aged_Sed"), measure = "weight", signaling = test, title.name = "Age Effect")
    Ex = netVisual_heatmap2(cellchat_Path,comparison =  c("Aged_Sed", "Aged_Ex"), measure = "weight", signaling = test, title.name = "Ex Effect")
    YoungEx = netVisual_heatmap2(cellchat_Path,comparison =  c("Young_Sed", "Young_Ex"), measure = "weight", signaling = test, title.name = "Young Ex Effect")
    MinDiff = min(c(min(as.vector(Age@matrix),na.rm = TRUE),min(as.vector(Ex@matrix),na.rm = TRUE),min(as.vector(YoungEx@matrix),na.rm = TRUE)))
    MaxDiff = max(c(max(as.vector(Age@matrix),na.rm = TRUE),max(as.vector(Ex@matrix),na.rm = TRUE),max(as.vector(YoungEx@matrix),na.rm = TRUE)))

    
    print("SetDif")
    YoungSed = netVisual_heatmap2(Young_Sed_Path, measure = "weight", signaling = test, title.name = "Young Sed")
    YoungEx = netVisual_heatmap2(Young_Ex_Path, measure = "weight", signaling = test, title.name = "Young Sed")
    AgeSed = netVisual_heatmap2(Aged_Sed_Path, measure = "weight", signaling = test, title.name = "Age Sed")
    AgeEx = netVisual_heatmap2(Aged_Ex_Path, measure = "weight", signaling = test, title.name = "Age Ex")



    Max = max(c(max(as.vector(YoungSed@matrix),na.rm = TRUE),max(as.vector(AgeSed@matrix) ,na.rm = TRUE),max(as.vector(AgeEx@matrix),na.rm = TRUE),max(as.vector(YoungEx@matrix),na.rm = TRUE)))

       MinDiff   
       MaxDiff
       Max       
    
    print("SetMinMax")
    Age =  netVisual_heatmap2(cellchat_Path,comparison =  c("Young_Sed", "Aged_Sed"), measure = "weight", signaling = test, title.name = paste(test,"Age Effect"), minCol = MinDiff, maxCol = MaxDiff)

    Ex = netVisual_heatmap2(cellchat_Path,comparison =  c("Aged_Sed", "Aged_Ex"), measure = "weight", signaling = test, title.name =  paste(test,"Age Ex Effect"), minCol = MinDiff, maxCol = MaxDiff)
    YoungEx = netVisual_heatmap2(cellchat_Path,comparison =  c("Young_Sed", "Young_Ex"), measure = "weight", signaling = test,title.name =  paste(test,"Young Ex Effect"), minCol = MinDiff, maxCol = MaxDiff)
    
    
    YoungSed = netVisual_heatmap2(Young_Sed_Path, measure = "weight", signaling = test, title.name =  paste(test,"Young Sed"), minCol = 0, maxCol = Max, color.heatmap = c("#f7f7f7", "#b2182b"))
    YoungExSing = netVisual_heatmap2(Young_Ex_Path, measure = "weight", signaling = test, title.name =  paste(test,"Young Ex"), minCol = 0, maxCol = Max, color.heatmap = c("#f7f7f7", "#b2182b"))
    AgeSed = netVisual_heatmap2(Aged_Sed_Path, measure = "weight", signaling = test, title.name =  paste(test,"Age Sed"), minCol = 0, maxCol = Max, color.heatmap = c("#f7f7f7", "#b2182b"))
    AgeEx = netVisual_heatmap2(Aged_Ex_Path, measure = "weight", signaling = test, title.name =  paste(test,"Age Ex"), minCol = 0, maxCol = Max, color.heatmap = c("#f7f7f7", "#b2182b"))

    nGenes = ggplot(enriched_interactionsframe, aes(x=cond, y=enriched_interactions)) +geom_bar(stat="identity")
    plots = list(Age,YoungEx,Ex,nGenes,YoungSed,YoungExSing,AgeSed,AgeEx)

    
    print("Draw")
    pdf(paste(test,"combined_heatmaps.pdf"), width = 15, height = 10)
    pushViewport(viewport(layout = grid.layout(3, 4))) # Create a layout with 2 row and 3 columns

    # Add the first heatmap to the left column
    pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
    draw(Age)
    popViewport() # Exit the left column viewport

    
    pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 1))
    draw(YoungEx)
    popViewport() # 
    
    
    
    
    # Add the second heatmap to the right column
    pushViewport(viewport(layout.pos.col = 3, layout.pos.row = 1))
    draw(Ex)
    popViewport() # Exit the right column viewport

    print(ggplot(enriched_interactionsframe, aes(x=cond, y=enriched_interactions)) +geom_bar(stat="identity"), vp = viewport(layout.pos.row = 1, layout.pos.col = 4))
  
    
    pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2))
    draw(YoungSed)
    popViewport()
    
    pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 2))
    draw(YoungExSing)
    popViewport()

    pushViewport(viewport(layout.pos.col = 3, layout.pos.row = 2))
    draw(AgeSed)
    popViewport()

    pushViewport(viewport(layout.pos.col = 4, layout.pos.row = 2))
    draw(AgeEx)
    popViewport()

    popViewport() # Exit the main viewport
    
    
    pushViewport(viewport(layout = grid.layout(3, 1)))
    dev.off() # Close the PDF device
    
}


options(warn=-1)

for (test in PathsFilt){
    

    setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2024-02-22.CellChat_Clean/pathways"))
    
    print(test)
    Aged_Ex_Path = aggregateNet2(
      Aged_Ex,
      sources.use = NULL,
      targets.use = NULL,
      signaling = test,
      pairLR.use = NULL,
      thresh = 0.05,
      remove.isolate = FALSE,
      return.object = TRUE
    )

    print("Aged_Ex DONE")
    Aged_Sed_Path = aggregateNet2(
      Aged_Sed,
      sources.use = NULL,
      targets.use = NULL,
      signaling = test,
      pairLR.use = NULL,
      thresh = 0.05,
         remove.isolate = FALSE,
      return.object = TRUE
    )
    print("Aged_Sed DONE")
    Young_Ex_Path = aggregateNet2(
      Young_Ex,
      sources.use = NULL,
      targets.use = NULL,
      signaling = test,
      pairLR.use = NULL,
      thresh = 0.05,
         remove.isolate = FALSE,
      return.object = TRUE
    )
    print("Young_Ex DONE")
    Young_Sed_Path = aggregateNet2(
      Young_Sed,
      sources.use = NULL,
      targets.use = NULL,
      signaling = test,
      pairLR.use = NULL,
      thresh = 0.05,
         remove.isolate = FALSE,
      return.object = TRUE
    )
    print("Young_Sed DONE")

    object.list <- list(Aged_Ex = Aged_Ex_Path, Aged_Sed = Aged_Sed_Path,Young_Ex = Young_Ex_Path,Young_Sed = Young_Sed_Path)
    
    cellchat_Path <- mergeCellChat(object.list, add.names = names(object.list))
    
    enriched_interactions = c(Aged_Sed = nrow(extractEnrichedLR(
      Aged_Sed_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    )),


     Aged_Ex = nrow(extractEnrichedLR(
      Aged_Ex_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    )),

     Young_Sed = nrow(extractEnrichedLR(
      Young_Sed_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    )),

     Young_Ex= nrow(extractEnrichedLR(
      Young_Ex_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    ))

    )
    enriched_interactionsframe =data.frame(enriched_interactions)

    enriched_interactionsframe$cond = rownames(enriched_interactionsframe)

    enriched_interactionsframe


 print("here")

    Age = netVisual_heatmap2(cellchat_Path,comparison =  c("Young_Sed", "Aged_Sed"), measure = "weight", signaling = test, title.name = "Age Effect",  cluster.rows = TRUE, cluster.cols = TRUE)
    Ex = netVisual_heatmap2(cellchat_Path,comparison =  c("Aged_Sed", "Aged_Ex"), measure = "weight", signaling = test, title.name = "Ex Effect",  cluster.rows = TRUE,cluster.cols = TRUE)
    YoungEx = netVisual_heatmap2(cellchat_Path,comparison =  c("Young_Sed", "Young_Ex"), measure = "weight", signaling = test, title.name = "Young Ex Effect")
    MinDiff = min(c(min(as.vector(Age@matrix),na.rm = TRUE),min(as.vector(Ex@matrix),na.rm = TRUE),min(as.vector(YoungEx@matrix),na.rm = TRUE)))
    MaxDiff = max(c(max(as.vector(Age@matrix),na.rm = TRUE),max(as.vector(Ex@matrix),na.rm = TRUE),max(as.vector(YoungEx@matrix),na.rm = TRUE)))

    
    print("SetDif")
   

    Max = max(c(max(as.vector(Age@matrix),na.rm = TRUE),max(as.vector(Ex@matrix) ,na.rm = TRUE)))

       MinDiff   
       MaxDiff
       Max       
    
    print("SetMinMax")
    Age =  netVisual_heatmap2(cellchat_Path,comparison =  c("Young_Sed", "Aged_Sed"), measure = "weight", signaling = test, title.name = paste(test,"Age Effect"), minCol = MinDiff, maxCol = MaxDiff,
                       cluster.rows = FALSE, row_order = row_order(Age),
                        cluster.cols = FALSE,  col_order = column_order(Age),
                         bar_row_colour = RNA_pallate_bar_row, bar_col_colour = RNA_pallate_bar_col,color.use = RNA_pallate)

    Ex = netVisual_heatmap2(cellchat_Path,comparison =  c("Aged_Sed", "Aged_Ex"), measure = "weight", signaling = test, title.name =  paste(test,"Age Ex Effect"), minCol = MinDiff, maxCol = MaxDiff,
                                   cluster.rows = FALSE, row_order = row_order(Age),
                                    cluster.cols = FALSE,  col_order = column_order(Age),
                         bar_row_colour = RNA_pallate_bar_row, bar_col_colour = RNA_pallate_bar_col,color.use = RNA_pallate)
     

    
    print("Draw")
    pdf(paste(test,"Age_heatmaps2.pdf"), width = 10, height = 10)
   
    draw(Age)
  
    dev.off() 
    
    
    pdf(paste(test,"AgeEx_heatmaps2.pdf"), width = 10, height = 10)
   
    draw(Ex)
  
    dev.off() # Close the PDF device
    
}


object.list <- list(MAE = Aged_Ex, MAS = Aged_Sed)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "MAE"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0, thresh.fc = 0, thresh.p = 0.1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS

# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.downEx <- subsetCommunication(cellchat, net = net, datasets = "MAS",ligand.logFC = NULL, receptor.logFC = NULL)
net.downEx = net.downEx[ net.downEx$receptor.logFC < 0 |  net.downEx$ligand.logFC < 0,]

object.list <- list(YS = Young_Sed, MAS = Aged_Sed)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "MAS"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0, thresh.fc =0, thresh.p = 0.1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.upAge <- subsetCommunication(cellchat, net = net, datasets = "MAS",ligand.logFC = NULL, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS

net.upAge = net.upAge[ net.upAge$receptor.logFC > 0 |  net.upAge$ligand.logFC > 0,]



net.upAge$all = paste(net.upAge$interaction_name,net.upAge$source,net.upAge$target)
net.downEx$all = paste(net.downEx$interaction_name,net.downEx$source,net.downEx$target)


reversed = data.frame(all = net.upAge$all[net.upAge$all %in% net.downEx$all], Path  = net.upAge$pathway_name[net.upAge$all %in% net.downEx$all], interaction_name = net.upAge$interaction_name[net.upAge$all %in% net.downEx$all], target = net.upAge$target[net.upAge$all %in% net.downEx$all], source = net.upAge$source[net.upAge$all %in% net.downEx$all])
plot_Rev = net.upAge[ net.upAge$all %in% reversed$all,]

reversed = reversed[ reversed$Path == "IGF",] 
reversed = na.omit(reversed)
reversed_df = data.frame(reversed$interaction_name)





na.omit(plot_Rev)
DT::datatable(na.omit(plot_Rev), filter = 'top')

  pdf(paste(test,"ReversedComs.pdf"), width = 10, height = 10)
netVisual_chord_gene(object.list[[1]], slot.name = 'net', net = plot_Rev, pairLR.use = reversed_df,lab.cex = 0.8,legend.pos.x = 0.5, legend.pos.y = 0.5,
                     small.gap = 3.5, title.name = paste0("Down-regulated signaling in Aged Sed vs Young Sed"), show.legend = T)

dev.off()

netVisual_bubble(Aged_Ex, remove.isolate = TRUE, signaling = c("IGF"), font.size = 20, thresh = 0.05)
netVisual_bubble(Aged_Sed, remove.isolate = TRUE, signaling = c("IGF"), font.size = 20, thresh = 0.05)
netVisual_bubble(Young_Sed, remove.isolate = TRUE, targets.use = "Lymphatic_EC", signaling = c("IGF"), font.size = 20, thresh = 0.05)



# TGF B and collagen

# Read in the shared results and merge togetehr 
setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2024-02-22.CellChat_Clean"))
    Aged_Ex = readRDS("Aged_Ex_Filt_cellChat.rds")
    Aged_Sed = readRDS("Aged_Sed_Filt_cellChat.rds")
    Young_Ex = readRDS("Young_Ex_Filt_cellChat.rds")
    Young_Sed = readRDS("Young_Sed_Filt_cellChat.rds")
object.list <- list(Aged_Ex = Aged_Ex, Aged_Sed = Aged_Sed,Young_Ex = Young_Ex,Young_Sed = Young_Sed)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

net_age <- subsetCommunication(Aged_Sed,ligand.logFC = NULL, receptor.logFC = NULL)
DT::datatable(net_age, filter = 'top')

# Specify comparison = c(1, 2, 3, 4) to pull from all four objects in the list

# (a bare "gg3" call stood here; removed because gg3 is only assigned a few lines below)
# Specify comparison = c(1, 2, 3, 4) to pull from all four objects in the list
gg3 <- netVisual_bubble(Aged_Sed,
                        sources.use = "Cardiomyocyte", 
                        targets.use = "Fibroblast",
                        signaling = c("TGFb"),
                        angle.x = 45)
gg3


head(cellchat@meta)

# 2. Set the main idents to that column (replace 'cell_type' with your actual column name)
cellchat <- setIdent(cellchat, ident.use = "age_condition") 

# 3. Verify it worked - this should no longer be NULL
levels(cellchat@idents)



# 1. Re-merge the objects
object.list <- list(Aged_Ex = Aged_Ex, Aged_Sed = Aged_Sed, 
                   Young_Ex = Young_Ex, Young_Sed = Young_Sed)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# 2. Update the pathway probabilities for the merged object (Crucial step)
#cellchat <- computeCommunProbPathway(cellchat)

library(dplyr)
library(ggplot2)

# 1. Extract the data for your specific cell pairs and pathways
df_comm <- subsetCommunication(cellchat, 
                               sources.use = "Cardiomyocyte", 
                               targets.use = "Fibroblast",
                               signaling = c("COLLAGEN", "TGFb"))

library(dplyr)
library(ggplot2)
library(viridis)

# 1. Combine the list of data frames into one 'long' format
# We use .id = "Condition" to turn the list names (Aged_Ex, etc.) into a column
df_combined <- bind_rows(df_comm, .id = "Condition") %>%
  dplyr::select(Condition, pathway_name, interaction_name_2, prob, pval) %>%
  # Create a cleaner label for the Y-axis combining pathway and LR pair
  mutate(interaction = paste0(pathway_name, ": ", interaction_name_2))

# 2. Set the factor levels to ensure the plot follows a logical order
# (e.g., Young to Aged, Sedentary to Exercise)
df_combined$Condition <- factor(df_combined$Condition, 
                                levels = c("Young_Sed", "Young_Ex", "Aged_Sed", "Aged_Ex"))

library(dplyr)
library(ggplot2)

# 1. Add the transformed p-value column to your existing dataframe
df_plot_ready <- df_combined %>%
  dplyr::mutate(log_pval = -log10(pval))

# 2. Generate the updated plot
ggplot(df_plot_ready, aes(x = Condition, y = interaction)) +
  # Map size to the transformed p-value and colour to the probability
  geom_point(aes(colour = prob, size = 5)) +
  
  # Colour scale for Strength (prob)
  scale_colour_viridis_c(option = "plasma", name = "Strength\n(Comm. Prob)") +
  
  # Size scale for Significance (-log10(pval))
  # You can adjust the 'range' to make the dots larger or smaller overall
  
  theme_bw() +
  labs(
    title = "CM-to-Fibroblast Signaling",
    subtitle = "Size indicates significance; Colour indicates interaction strength",
    x = "Experimental Group",
    y = "Signaling Interaction"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_line(colour = "grey90"),
    legend.position = "right"
  )

library(dplyr)
library(ggplot2)

# 1. Extract the data for your specific cell pairs and pathways
df_comm <- subsetCommunication(cellchat, 
                               sources.use = "Fibroblast", 
                               targets.use = "Cardiomyocyte",
                               signaling = c("COLLAGEN", "TGFb"))

library(dplyr)
library(ggplot2)
library(viridis)

# 1. Combine the list of data frames into one 'long' format
# We use .id = "Condition" to turn the list names (Aged_Ex, etc.) into a column
df_combined <- bind_rows(df_comm, .id = "Condition") %>%
  dplyr::select(Condition, pathway_name, interaction_name_2, prob, pval) %>%
  # Create a cleaner label for the Y-axis combining pathway and LR pair
  mutate(interaction = paste0(pathway_name, ": ", interaction_name_2))

# 2. Set the factor levels to ensure the plot follows a logical order
# (e.g., Young to Aged, Sedentary to Exercise)
df_combined$Condition <- factor(df_combined$Condition, 
                                levels = c("Young_Sed", "Young_Ex", "Aged_Sed", "Aged_Ex"))

library(dplyr)
library(ggplot2)

# 1. Add the transformed p-value column to your existing dataframe
df_plot_ready <- df_combined %>%
  dplyr::mutate(log_pval = -log10(pval))

# 2. Generate the updated plot
ggplot(df_plot_ready, aes(x = Condition, y = interaction)) +
  # Map size to the transformed p-value and colour to the probability
  geom_point(aes(colour = prob, size = 5)) +
  
  # Colour scale for Strength (prob)
  scale_colour_viridis_c(option = "plasma", name = "Strength\n(Comm. Prob)") +
  
  # Size scale for Significance (-log10(pval))
  # You can adjust the 'range' to make the dots larger or smaller overal
  
  theme_bw() +
  labs(
    title = "Fibroblast-to-CM Signaling",
    subtitle = "Size indicates significance; Colour indicates interaction strength",
    x = "Experimental Group",
    y = "Signaling Interaction"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_line(colour = "grey90"),
    legend.position = "right"
  )

df_comm <- subsetCommunication(cellchat, 
                               signaling = c("COLLAGEN", "TGFb"))
df_combined <- bind_rows(df_comm, .id = "Condition") %>%
  # Create a cleaner label for the Y-axis combining pathway and LR pair
  mutate(interaction = paste0(pathway_name, ": ", interaction_name_2))

DT::datatable(df_combined, filter = 'top')



options(warn=-1)
PathsFilt = c("COLLAGEN")

test = 'COLLAGEN'
    setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2024-02-22.CellChat_Clean/pathways"))
    
    print(test)
    Aged_Ex_Path = aggregateNet2(
      Aged_Ex,
      sources.use = NULL,
      targets.use = NULL,
      signaling = test,
      pairLR.use = NULL,
      thresh = 0.05,
      remove.isolate = FALSE,
      return.object = TRUE
    )

    print("Aged_Ex DONE")
    Aged_Sed_Path = aggregateNet2(
      Aged_Sed,
      sources.use = NULL,
      targets.use = NULL,
      signaling = test,
      pairLR.use = NULL,
      thresh = 0.05,
         remove.isolate = FALSE,
      return.object = TRUE
    )
    print("Aged_Sed DONE")
    Young_Ex_Path = aggregateNet2(
      Young_Ex,
      sources.use = NULL,
      targets.use = NULL,
      signaling = test,
      pairLR.use = NULL,
      thresh = 0.05,
         remove.isolate = FALSE,
      return.object = TRUE
    )
    print("Young_Ex DONE")
    Young_Sed_Path = aggregateNet2(
      Young_Sed,
      sources.use = NULL,
      targets.use = NULL,
      signaling = test,
      pairLR.use = NULL,
      thresh = 0.05,
         remove.isolate = FALSE,
      return.object = TRUE
    )
    print("Young_Sed DONE")

    object.list <- list(Aged_Ex = Aged_Ex_Path, Aged_Sed = Aged_Sed_Path,Young_Ex = Young_Ex_Path,Young_Sed = Young_Sed_Path)
    
    cellchat_Path <- mergeCellChat(object.list, add.names = names(object.list))
    
    enriched_interactions = c(Aged_Sed = nrow(extractEnrichedLR(
      Aged_Sed_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    )),


     Aged_Ex = nrow(extractEnrichedLR(
      Aged_Ex_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    )),

     Young_Sed = nrow(extractEnrichedLR(
      Young_Sed_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    )),

     Young_Ex= nrow(extractEnrichedLR(
      Young_Ex_Path,
      test,
      geneLR.return = FALSE,
      enriched.only = TRUE,
      thresh = 0.05,
      geneInfo = NULL,
      complex_input = NULL
    ))

    )
    enriched_interactionsframe =data.frame(enriched_interactions)

    enriched_interactionsframe$cond = rownames(enriched_interactionsframe)

   # enriched_interactionsframe



df_comm <- subsetCommunication(Aged_Ex_Path)
DT::datatable(df_comm, filter ='top')

 AgeEx = netVisual_heatmap2(Aged_Ex_Path, measure = "weight", signaling = test, title.name = "Age Ex",color.heatmap = "Reds", )
AgeEx
