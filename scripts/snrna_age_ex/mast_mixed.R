# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Fit mixed-effects MAST models as a robustness analysis.
# Steps: 1) Load the annotated RNA object; 2) include animal-level structure; 3) compare effect estimates with primary fixed-effects results.
# =======================================
# Source reformatted from snRNA_AgeEx/MAST_Mixed(1).r.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

library(Seurat)
library(SingleCellExperiment)
library(MAST)
library(future)

options(mc.cores =28)


setwd(publication_path("HEART_SNRNA_ROOT", "data/2023-03-01.CellbenderFiltered_MAS3_2_Adjust"))
Mouse = readRDS(file = "CellTypes.seurat.rds")


Validindv =  c("MAE2_1_M","MAE2_1_F","MAS1_4_F","MAS1_4_M","MAS3_2_M","MAS3_2_F",
 "YE2_1_F","YE2_1_M","YE4_3_F","YE4_3_M","YS1_4_M","YS3_2_M","YS3_2_F")



RNA = subset(Mouse, subset = ID %in% Validindv)


freq_expressed <- 0.2
FCTHRESHOLD <- log2(1.5)

cell ="B_cell"
args = commandArgs(trailingOnly=TRUE)

cell = args[1] 
Allcase = c("MAE","YE","MAE","MAS")
    
Allctrl = c("MAS","YS","YE","YS")

case = Allcase[4]
ctrl = Allctrl[4]
sub = subset(RNA, cellType == cell)


#get count table 

  sce <-SingleCellExperiment(assays=list(counts=sub[['RNA']]@counts), colData=sub@meta.data)
  #normlization TPM
  library(scater)
  TPM=calculateTPM(sce, exprs_values="counts")
  logTPM<-log2(TPM+1)
  fData <- data.frame(names=rownames(logTPM))
  rownames(fData) <- rownames(logTPM)
  cData <- sub[[]] #data.frame(condition=seurat$condition, age=seurat$age, condition=seurat$condition, ID=seurat$ID, cluster=seurat$seurat_clusters, sex=seurat$sex, batch=seurat$batch)
  rownames(cData) <- colnames(logTPM)
  library(MAST)
  scaRaw<- FromMatrix(as.matrix(logTPM), cData, fData)
  expressed_genes <- freq(scaRaw) > freq_expressed
  scaRaw <- scaRaw[expressed_genes,]  
  cdr2 <- colSums(SummarizedExperiment::assay(scaRaw)>0)

#annotate metadata
SummarizedExperiment::colData(scaRaw)$ngeneson <- scale(cdr2)
SummarizedExperiment::colData(scaRaw)$group <-factor(SummarizedExperiment::colData(scaRaw)$group)
SummarizedExperiment::colData(scaRaw)$ID <-factor(SummarizedExperiment::colData(scaRaw)$ID)




#Run MAST
setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-03-09.MAST_Mixed_sex"))
#Age Sed 
#Base = YS
#con = MAS
scaRaw$group = relevel(scaRaw$group, ref="YS") 
zlmCond <- (MAST::zlm(~group+ ngeneson  +sex + (1|ID),scaRaw, method='glmer',ebayes = F,strictConvergence = FALSE))

saveRDS(zlmCond,paste0(cell,"MAS-YS.rds"))
summaryCond <- summary(zlmCond, doLRT='groupMAS') 
summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[summaryDt$contrast=='groupMAS'
                                  & summaryDt$component=='logFC', c(1,7,5,6,8)],
                        summaryDt[summaryDt$contrast=='groupMAS'
                                  & summaryDt$component=='H', c(1,4)],
                        by = 'primerid')

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdle <- stats::na.omit(as.data.frame(fcHurdle))
fcHurdle=fcHurdle[order(fcHurdle$fdr),]
head(fcHurdle)
write.csv(fcHurdle,paste0("MAST",cell,"_MAS-YS.csv"))
                 
                 
#Age E 
#Base = YE 
#con = MAE
scaRaw$group = relevel(scaRaw$group, ref="YE") 
zlmCond <- (MAST::zlm(~group+ ngeneson +sex   + (1|ID),scaRaw, method='glmer',ebayes = F,strictConvergence = FALSE))
saveRDS(zlmCond,paste0(cell,"MAE-YE.rds"))
summaryCond <- summary(zlmCond, doLRT='groupMAE') 
summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[summaryDt$contrast=='groupMAE'
                                  & summaryDt$component=='logFC', c(1,7,5,6,8)],
                        summaryDt[summaryDt$contrast=='groupMAE'
                                  & summaryDt$component=='H', c(1,4)],
                        by = 'primerid')
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdle <- stats::na.omit(as.data.frame(fcHurdle))
fcHurdle=fcHurdle[order(fcHurdle$fdr),]
head(fcHurdle)
write.csv(fcHurdle,paste0("MAST",cell,"_MAE-YE.csv"))                 


#Ex Aged 
#Base = MAS 
#con = MAE
scaRaw$group = relevel(scaRaw$group, ref="MAS") 
zlmCond <- (MAST::zlm(~group+ ngeneson +sex   + (1|ID),scaRaw, method='glmer',ebayes = F,strictConvergence = FALSE))
saveRDS(zlmCond,paste0(cell,"MAE-MAS.rds"))
summaryCond <- summary(zlmCond, doLRT='groupMAE') 
summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[summaryDt$contrast=='groupMAE'
                                  & summaryDt$component=='logFC', c(1,7,5,6,8)],
                        summaryDt[summaryDt$contrast=='groupMAE'
                                  & summaryDt$component=='H', c(1,4)],
                        by = 'primerid')
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdle <- stats::na.omit(as.data.frame(fcHurdle))
fcHurdle=fcHurdle[order(fcHurdle$fdr),]
head(fcHurdle)
write.csv(fcHurdle,paste0("MAST",cell,"_MAE-MAS.csv"))  


#Ex Young 
#Base = YS 
#con = YE
scaRaw$group = relevel(scaRaw$group, ref="YS") 
zlmCond <- (MAST::zlm(~group+ ngeneson +sex   + (1|ID),scaRaw, method='glmer',ebayes = F,strictConvergence = FALSE))
saveRDS(zlmCond,paste0(cell,"YE-YS.rds"))
summaryCond <- summary(zlmCond, doLRT='groupYE') 
summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[summaryDt$contrast=='groupYE'
                                  & summaryDt$component=='logFC', c(1,7,5,6,8)],
                        summaryDt[summaryDt$contrast=='groupYE'
                                  & summaryDt$component=='H', c(1,4)],
                        by = 'primerid')
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdle <- stats::na.omit(as.data.frame(fcHurdle))
fcHurdle=fcHurdle[order(fcHurdle$fdr),]
head(fcHurdle)
write.csv(fcHurdle,paste0("MAST",cell,"_YE-YS.csv"))  


