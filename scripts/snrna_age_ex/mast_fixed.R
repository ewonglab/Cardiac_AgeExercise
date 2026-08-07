# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Fit primary fixed-effects MAST models for aging and exercise contrasts.
# Steps: 1) Load the annotated RNA object; 2) test genes by cell type with sex/CDR covariates; 3) export differential-expression tables.
# =======================================
# Source reformatted from snRNA_AgeEx/MAST_Fixed.r.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

library(MAST)
library(SingleCellExperiment)
setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed"))



freq_expressed <- 0.1
FCTHRESHOLD <- log2(1.5)

setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed"))
Mouse = readRDS(publication_path("HEART_SNRNA_ROOT", "data/2022-10-17.SharedObject/2023-10-9_Full_RNA_TooManyAnnotated.rds"))



freq_expressed <- 0.2
FCTHRESHOLD <- log2(1.5)
args = commandArgs(trailingOnly=TRUE)

cell = args[1]
sub <- subset(Mouse, CellType == cell)



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
SummarizedExperiment::colData(scaRaw)$age_condition <-factor(SummarizedExperiment::colData(scaRaw)$age_condition)
SummarizedExperiment::colData(scaRaw)$mouseID <-factor(SummarizedExperiment::colData(scaRaw)$mouseID)




#Run MAST
setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed"))
#Age Sed
#Base = Young_Sed
#con = Aged_Sed
scaRaw$age_condition = relevel(scaRaw$age_condition, ref="Young_Sed")
zlmCond <- (MAST::zlm(~age_condition+ ngeneson   + sex ,scaRaw, method='bayesglm'))

saveRDS(zlmCond,paste0(cell,"Aged_Sed-Young_Sed.rds"))
summaryCond <- summary(zlmCond, doLRT='age_conditionAged_Sed')
summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[summaryDt$contrast=='age_conditionAged_Sed'
                                  & summaryDt$component=='logFC', c(1,7,5,6,8)],
                        summaryDt[summaryDt$contrast=='age_conditionAged_Sed'
                                  & summaryDt$component=='H', c(1,4)],
                        by = 'primerid')

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdle <- stats::na.omit(as.data.frame(fcHurdle))
fcHurdle=fcHurdle[order(fcHurdle$fdr),]
head(fcHurdle)
write.csv(fcHurdle,paste0("MAST",cell,"_Aged_Sed-Young_Sed.csv"))


#Age E
#Base = Young_Ex
#con = Aged_Ex
scaRaw$age_condition = relevel(scaRaw$age_condition, ref="Young_Ex")
zlmCond <- (MAST::zlm(~age_condition+ ngeneson   + sex ,scaRaw, method='bayesglm'))
saveRDS(zlmCond,paste0(cell,"Aged_Ex-Young_Ex.rds"))
summaryCond <- summary(zlmCond, doLRT='age_conditionAged_Ex')
summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[summaryDt$contrast=='age_conditionAged_Ex'
                                  & summaryDt$component=='logFC', c(1,7,5,6,8)],
                        summaryDt[summaryDt$contrast=='age_conditionAged_Ex'
                                  & summaryDt$component=='H', c(1,4)],
                        by = 'primerid')
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdle <- stats::na.omit(as.data.frame(fcHurdle))
fcHurdle=fcHurdle[order(fcHurdle$fdr),]
head(fcHurdle)
write.csv(fcHurdle,paste0("MAST",cell,"_Aged_Ex-Young_Ex.csv"))


#Ex Aged
#Base = Aged_Sed
#con = Aged_Ex
scaRaw$age_condition = relevel(scaRaw$age_condition, ref="Aged_Sed")
zlmCond <- (MAST::zlm(~age_condition+ ngeneson   + sex ,scaRaw, method='bayesglm'))
saveRDS(zlmCond,paste0(cell,"Aged_Ex-Aged_Sed.rds"))
summaryCond <- summary(zlmCond, doLRT='age_conditionAged_Ex')
summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[summaryDt$contrast=='age_conditionAged_Ex'
                                  & summaryDt$component=='logFC', c(1,7,5,6,8)],
                        summaryDt[summaryDt$contrast=='age_conditionAged_Ex'
                                  & summaryDt$component=='H', c(1,4)],
                        by = 'primerid')
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdle <- stats::na.omit(as.data.frame(fcHurdle))
fcHurdle=fcHurdle[order(fcHurdle$fdr),]
head(fcHurdle)
write.csv(fcHurdle,paste0("MAST",cell,"_Aged_Ex-Aged_Sed.csv"))


#Ex Young
#Base = Young_Sed
#con = Young_Ex
scaRaw$age_condition = relevel(scaRaw$age_condition, ref="Young_Sed")
zlmCond <- (MAST::zlm(~age_condition+ ngeneson   + sex ,scaRaw, method='bayesglm'))
saveRDS(zlmCond,paste0(cell,"Young_Ex-Young_Sed.rds"))
summaryCond <- summary(zlmCond, doLRT='age_conditionYoung_Ex')
summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[summaryDt$contrast=='age_conditionYoung_Ex'
                                  & summaryDt$component=='logFC', c(1,7,5,6,8)],
                        summaryDt[summaryDt$contrast=='age_conditionYoung_Ex'
                                  & summaryDt$component=='H', c(1,4)],
                        by = 'primerid')
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdle <- stats::na.omit(as.data.frame(fcHurdle))
fcHurdle=fcHurdle[order(fcHurdle$fdr),]
head(fcHurdle)
write.csv(fcHurdle,paste0("MAST",cell,"_Young_Ex-Young_Sed.csv"))





