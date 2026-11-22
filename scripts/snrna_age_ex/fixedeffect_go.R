# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Run fixed-effects MAST downstream enrichment and GO visualisation.
# Steps: 1) Load differential-expression results; 2) select significant genes; 3) perform enrichment and export figures/tables.
# =======================================
# Source reformatted from snRNA_AgeEx/FixedEffect_Go.ipynb.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

library(ggplot2)
library(clusterProfiler, lib.loc = publication_path("R_LIBS_DIR"))
library(Seurat)
library(cicero)
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(rtracklayer)
library(ggpubr)
library(dplyr)
library(tidyr)
library(stringr)
library(tidyquant)
library(ggrepel)
library(openxlsx)



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
library(ggplot2)



library(openxlsx)
library(stringr)
library(EnhancedVolcano)
setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed"))
RandomBatch=list.files(pattern = 'MAST.*csv', recursive = TRUE, include.dirs = TRUE)

i = RandomBatch[1]
FC = 0.5
P = 0.05
All = data.frame(matrix(nrow = 0,ncol = 15))
for (i in RandomBatch) {
 
  names =i
 
    
    print(i)
    setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed"))
    table =  read.csv(i)
    setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed/GO"))
      name = str_remove(i, "MAST")
      name = str_remove(name, ".csv")
      condition = str_extract(name,"Aged.*|Young.*")
      cell = str_remove(name, "_Aged.*|_Young.*")
      table$cell = cell
      table$condition = condition
    #table$P = as.numeric(table$P)
    
   
    
    
    
      All = rbind(All,table)
  }



#Age common
Sig = All[ All$fdr < 0.01,]
Age = Sig[ Sig$condition == 'Aged_Sed-Young_Sed',]



AgeIncreasing = Age[ Age$coef > 0 ,]
n_AgeIncreasing = data.frame(sort(table(AgeIncreasing$primerid), decreasing = TRUE))
AgeDecreasing = Age[ Age$coef < 0 ,]
n_AgeDecreasing = data.frame(sort(table(AgeDecreasing$primerid), decreasing = TRUE))

n_AgeIncreasing_non1 = n_AgeIncreasing[ n_AgeIncreasing$Freq !=1 ,]
ggplot(n_AgeIncreasing_non1, aes(x = as.factor(Freq))) +
  geom_bar() + scale_y_continuous( expand = c(0,0.1))+
  labs(title = "Frequency Distribution", x = "Frequency", y = "Count")
n_AgeDecreasing_non1 = n_AgeDecreasing[ n_AgeDecreasing$Freq !=1 ,]
ggplot(n_AgeDecreasing_non1, aes(x = as.factor(Freq))) +
  geom_bar() + scale_y_continuous( expand = c(0,0.1))+
  labs(title = "Frequency Distribution", x = "Frequency", y = "Count")

#Age common
Sig = All[ All$fdr < 0.01,]
Ex = Sig[ Sig$condition == 'Aged_Ex-Aged_Sed',]



ExIncreasing = Ex[ Ex$coef > 0 ,]
n_ExIncreasing = data.frame(sort(table(ExIncreasing$primerid), decreasing = TRUE))
ExDecreasing = Ex[ Ex$coef < 0 ,]
n_ExDecreasing = data.frame(sort(table(ExDecreasing$primerid), decreasing = TRUE))

n_ExIncreasing_non1 = n_ExIncreasing[ n_ExIncreasing$Freq !=1 ,]
ggplot(n_ExIncreasing_non1, aes(x = as.factor(Freq))) +
  geom_bar() + scale_y_continuous( expand = c(0,0.1))+
  labs(title = "Frequency Distribution", x = "Frequency", y = "Count")
n_ExDecreasing_non1 = n_ExDecreasing[ n_ExDecreasing$Freq !=1 ,]
ggplot(n_ExDecreasing_non1, aes(x = as.factor(Freq))) +
  geom_bar() + scale_y_continuous( expand = c(0,0.1))+
  labs(title = "Frequency Distribution", x = "Frequency", y = "Count")

#get ones common to 4 or more 
n_ExIncreasingTop = unique(n_ExIncreasing[ n_ExIncreasing$Freq >= 4 ,'Var1'])
head(n_ExIncreasingTop)
n_AgeIncreasingTop = unique(n_AgeIncreasing[ n_AgeIncreasing$Freq >= 4 ,'Var1'])
head(n_AgeIncreasingTop)
n_ExDecreasingTop = unique( n_ExDecreasing[ n_ExDecreasing$Freq >= 4 ,'Var1'])
head(n_ExDecreasingTop)
n_AgeDecreasingTop = unique(n_AgeDecreasing[ n_AgeDecreasing$Freq >= 4 ,'Var1'])
head(n_AgeDecreasingTop)

library(org.Mm.eg.db)
AnnotationDbi::species(org.Mm.eg.db)
library(clusterProfiler)
 
 ExIncreasing_Go <- enrichGO(gene          = n_ExIncreasingTop,
                            universe      = unique(All$primerid),
                            OrgDb         = org.Mm.eg.db,
                            ont           = "ALL",
                            pAdjustMethod = "fdr",
                            pvalueCutoff=0.9,keyType = 'SYMBOL')
ExIncreasing_Go

ExDecreasing_Go <- enrichGO(gene          = n_ExDecreasingTop,
                            universe      = unique(All$primerid),
                            OrgDb         = org.Mm.eg.db,
                            ont           = "ALL",
                            pAdjustMethod = "fdr",
                            pvalueCutoff=1,keyType = 'SYMBOL')
ExDecreasing_Go



AgeIncreasing_Go <- enrichGO(gene          = n_AgeIncreasingTop,
                            universe      = unique(All$primerid),
                            OrgDb         = org.Mm.eg.db,
                            ont           = "ALL",
                            pAdjustMethod = "fdr",
                            pvalueCutoff=1,keyType = 'SYMBOL')
AgeIncreasing_Go

AgeDecreasing_Go <- enrichGO(gene          = n_AgeDecreasingTop,
                            universe      = unique(All$primerid),
                            OrgDb         = org.Mm.eg.db,
                            ont           = "ALL",
                            pAdjustMethod = "fdr",
                            pvalueCutoff=1,keyType = 'SYMBOL')
AgeDecreasing_Go

dotplot(ExIncreasing_Go)

dotplot(ExDecreasing_Go)

dotplot(AgeIncreasing_Go)

dotplot(AgeDecreasing_Go)

1+1

library(openxlsx)
library(stringr)
library(EnhancedVolcano)
setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed"))
RandomBatch=list.files(pattern = 'MAST.*csv', recursive = TRUE, include.dirs = TRUE)
#RandomBatch = RandomBatch[grepl(RandomBatch, pattern = "M1|M2|Macrop")]
RandomBatch

library(org.Mm.eg.db)
AnnotationDbi::species(org.Mm.eg.db)
library(clusterProfiler)
i = RandomBatch[1]
FC = 0.2
P = 0.05
All = data.frame(matrix(nrow = 0,ncol = 15))
for (i in RandomBatch) {
 
  names =i
 
    
    print(i)
    setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed"))
    table =  read.csv(i)
    setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed/GO"))
      name = str_remove(i, "MAST")
      name = str_remove(name, ".csv")
      condition = str_extract(name,"Aged.*|Young.*")
      cell = str_remove(name, "_Aged.*|_Young.*")
      table$cell = cell
      table$condition = condition
   
    
    AnnotationDbi::species(org.Mm.eg.db)
        Order = table[order(table$coef, decreasing = TRUE ),]
    ego_supr = NULL
    ego_act = NULL
    testFr = NULL
        if(length(Order[Order$coef > FC & Order$fdr < P,"primerid"]) > 0 ){
            ego_act <- enrichGO(gene          = Order[Order$coef > FC & Order$fdr < P,"primerid"],
                            universe      = unique(table$primerid),
                            OrgDb         = org.Mm.eg.db,
                            ont           = "ALL",
                            pAdjustMethod = "fdr",
                            pvalueCutoff=0.9,keyType = 'SYMBOL')
            
        }
        if(length(Order[Order$coef <  -FC & Order$fdr < P,"primerid"]) > 0 ){
            ego_supr <- enrichGO(gene          = Order[Order$coef < - FC & Order$fdr < P,"primerid"],
                            universe      =  unique(table$primerid),
                            OrgDb         = org.Mm.eg.db,
                            ont           = "ALL",
                            pAdjustMethod = "fdr",
                            pvalueCutoff=0.9,keyType = 'SYMBOL')
            
        }
        if (!(is.null(ego_act))){
             if (nrow(ego_act@result) != 0 ){
            print(sort(ego_act@gene))
            ego_act@result$Sign = "+"
            }
            testFr = ego_act
            testFr@result = ego_act@result
            
        }
        if (!(is.null(ego_supr))){
            if (nrow(ego_supr@result) != 0 ){
            print(sort(ego_supr@gene))
            ego_supr@result$Sign = "-"
            }
            
            if(is.null(testFr)){
                testFr = ego_supr
                testFr@result = ego_supr@result
            }else{
               testFr@result = rbind(testFr@result,            ego_supr@result)
                
            }
        }
       
       
        
             
            if (!is.null(testFr)){
            if (nrow(testFr@result != 0 )){
            dotplot(testFr, showCategory=20,split="Sign", orderBy = "Sign",title  =paste(cell,condition)) + facet_grid(.~Sign)
            ggsave(paste0(cell,"_",condition, "_FC",FC,"_Fdr",P,"_GO.png"), width = 15, height = 15)
            }
            }

        saveRDS(testFr,paste0(cell,"_",condition, "_FC",FC,"_Fdr",P,"_GO.rds"))

    
    
    
      All = rbind(All,table)
  }




library(openxlsx)
library(stringr)
library(EnhancedVolcano)
setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed/GO"))
RandomBatch=list.files(pattern = '_0.0_GO.rds"', recursive = TRUE, include.dirs = TRUE)
RandomBatch

library(clusterProfiler)
RandomBatch
All = data.frame(matrix(nrow = 0,ncol = 15))
for (i in RandomBatch) {
 
  names =i

      name = str_remove(i, "_GO.rds")
      condition = str_extract(name,"Aged.*|Young.*")
      cell = str_remove(name, "_Aged.*|_Young.*")
      
    
    print(i)
    table =  readRDS(i)
    if (is.null(table)){(next)}
  table@result = table@result[ table@result$Sign == '+',] 
    if( nrow(table@result) == 0){pos = NULL}else{
        simple = simplify(
          table,
              cutoff = 0.5,
              by = "p.adjust",
              select_fun = min,
              measure = "Wang",
              semData = NULL
           )
       pos =  dotplot(simple, showCategory=20) 
    }
    
    table =  readRDS(i)
    table@result = table@result[ table@result$Sign == '-',] 
     if( nrow(table@result) == 0){neg = NULL}else{
        simple = simplify(
          table,
              cutoff = 0.5,
              by = "p.adjust",
              select_fun = min,
              measure = "Wang",
              semData = NULL
           )
       neg =  dotplot(simple, showCategory=20) 
     }
  ggarrange(pos,neg)
  ggsave(paste0("simple__0.0_GO",cell,"_",condition,".png"),width =15)
} 
    
    

setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed/GO"))
RandomBatchgsea=list.files(pattern = '_GO.rds', recursive = TRUE, include.dirs = TRUE)
All = data.frame(matrix(nrow = 0,ncol = 16))
for (i in RandomBatchgsea) {
 
  names =i
 
    
    print(i)
    table =  readRDS(i)
    if(is.null(table)){next}
    table = table@result
     if(nrow(table)==0){next}
      name = str_remove(i, "MAST")
      name = str_remove(name, "_GO.rds")
      condition = str_remove(name, ".*_")
      cell = str_remove(name, "_.*")
      table$cell = cell
      table$condition = condition
    table$name = name
    #table$P = as.numeric(table$P)
   

      All = rbind(All,table)
} 
  All = All[!grepl(pattern = "Unk", All$name),]
#map so ID Decriptuon and ontololgy match
 
 writOut =All[All$p.adjust < 0.1,]
 
 write.csv(All, "AllGO_0.0_GO_results.csv")
 

df = All
#df = df[df$Description %in% PathwaysKeep$Var1,]
library(ggplot2)
library(dplyr)
library(stringr)
library(forcats)
## from Tommy's code

  
  cellResult = df
    
    x_order = unique(cellResult$name)

  cellResult  = cellResult[cellResult$p.adjust < 0.1,]
 cellResult$GeneRatioCal =sapply(cellResult$GeneRatio, function(x) eval(parse(text = x)))

  #cellResult  = cellResult[nchar(cellResult$Description) < 100,]
  AGE  = cellResult
    if(nrow(AGE) < 1){next}
  p <- ggplot(AGE, aes(x = name, y = fct_reorder(Description, p.adjust))) + 
               geom_point(aes(size = GeneRatioCal, color = p.adjust)) +
               theme_bw(base_size = 14) +
         scale_x_discrete(limits = x_order, labels = x_order)+
    #scale_y_discrete(limits = filter$Var1)+
        scale_colour_gradient(limits=c(0, 0.10), low="red") +
        ylab(NULL) +
        ggtitle(paste("GO pathway enrichment"))+
     theme(axis.text.x = element_text(angle = 45, hjust=1, size = 30), text =element_text(size = 20),axis.text.y =element_text(size = 20))

    p+ facet_grid(.~Sign)


  ggsave(paste0("barplot/","All.png"), width = 30, height = 30)
                                 
                                 
  #dotplot_rev(AGE,showCategory = 10, by ="count")

df = All
#df = df[df$Description %in% PathwaysKeep$Var1,]
library(ggplot2)
library(dplyr)
library(stringr)
library(forcats)
## from Tommy's code

  
  cellResult = df
    
    x_order = unique(cellResult$name)

  cellResult  = cellResult[cellResult$p.adjust < 0.05,]
 cellResult$GeneRatioCal =sapply(cellResult$GeneRatio, function(x) eval(parse(text = x)))

  #cellResult  = cellResult[nchar(cellResult$Description) < 100,]
  AGE  = cellResult
    if(nrow(AGE) < 1){next}
  p <- ggplot(AGE, aes(x = name, y = fct_reorder(Description, p.adjust))) + 
               geom_point(aes(size = GeneRatioCal, color = p.adjust)) +
               theme_bw(base_size = 14) +
         scale_x_discrete(limits = x_order, labels = x_order)+
    #scale_y_discrete(limits = filter$Var1)+
        scale_colour_gradient(limits=c(0, 0.10), low="red") +
        ylab(NULL) +
        ggtitle(paste("GO pathway enrichment"))+
     theme(axis.text.x = element_text(angle = 45, hjust=1, size = 30), text =element_text(size = 20),axis.text.y =element_text(size = 20))

    p+ facet_grid(.~Sign)


  ggsave(paste0("barplot/","All.png"), width = 30, height = 30)
                                 
                                 
  #dotplot_rev(AGE,showCategory = 10, by ="count")

df = All
#df = df[df$Description %in% PathwaysKeep$Var1,]
library(ggplot2)
library(dplyr)
library(stringr)
library(forcats)
## from Tommy's code

 for (celltype in unique(cellResult$cell)){ 
  cellResult = df
    cellResult = cellResult[cellResult$cell == celltype,]
    x_order = unique(cellResult$name)

  cellResult  = cellResult[cellResult$p.adjust < 0.05,]
 cellResult$GeneRatioCal =sapply(cellResult$GeneRatio, function(x) eval(parse(text = x)))
#for each cell

                                 #cellResult  = cellResult[nchar(cellResult$Description) < 100,]
  AGE  = cellResult
    if(nrow(AGE) < 1){next}
  p <- ggplot(AGE, aes(x = name, y = fct_reorder(Description, p.adjust))) + 
               geom_point(aes(size = GeneRatioCal, color = p.adjust)) +
               theme_bw(base_size = 14) +
         scale_x_discrete(limits = x_order, labels = x_order)+
    #scale_y_discrete(limits = filter$Var1)+
        scale_colour_gradient(limits=c(0, 0.10), low="red") +
        ylab(NULL) +
        ggtitle(paste(celltype,"GO pathway enrichment"))+
     theme(axis.text.x = element_text(angle = 45, hjust=1, size = 30), text =element_text(size = 20),axis.text.y =element_text(size = 20))

    p+ facet_grid(.~Sign)


  ggsave(paste0("barplot/",celltype,".png"), width = 30, height = 30)
}                        
                                 
  #dotplot_rev(AGE,showCategory = 10, by ="count")

library(org.Mm.eg.db)
AnnotationDbi::species(org.Mm.eg.db)
library(clusterProfiler)

library(openxlsx)
library(stringr)
library(EnhancedVolcano)
setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed"))
RandomBatch=list.files(pattern = 'MAST.*csv', recursive = TRUE, include.dirs = TRUE)
RandomBatch = RandomBatch[grepl("Aged_Ex-Aged_Sed.csv|_Aged_Sed-Young_Sed.csv",RandomBatch)]
RandomBatch

library(org.Mm.eg.db)
AnnotationDbi::species(org.Mm.eg.db)
library(clusterProfiler)
i = RandomBatch[1]
FC = 0.5
P = 0.05
All = data.frame(matrix(nrow = 0,ncol = 15))
print("running")
for (i in RandomBatch) {
       print(i)
      names =i
      name = str_remove(i, "MAST")
      name = str_remove(name, ".csv")
      condition = str_extract(name,"Aged.*|Young.*")
      cell = str_remove(name, "_Aged.*|_Young.*")
    
   
    
     setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed/GSEA_1000"))
    
   if (file.exists(paste0(cell,"_",condition, "GSEA_GO.rds"))) {
       next
    }
    setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed"))
     print(i)
    table =  read.csv(i)
    setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed/GSEA_1000"))
     
      table$cell = cell
      table$condition = condition
   
    
    AnnotationDbi::species(org.Mm.eg.db)
     geneList = table[,'coef']
    names(geneList) = as.character(table$primerid)
    geneList = sort(geneList, decreasing = TRUE)
    
    ego_supr = NULL
    ego_act = NULL
    testFr = NULL
     ego_act <- gseGO(geneList     = geneList,
              OrgDb        = org.Mm.eg.db,
              ont          = "All",
              minGSSize    = 10,
              maxGSSize    = 1000,
              pvalueCutoff = 0.5,
                         keyType = "SYMBOL",
              verbose      = FALSE)
     
            
   
        
             
            if (!is.null(ego_act)){
            if (nrow(ego_act@result != 0 )){
            ego_act@result$direction = ifelse(ego_act@result$NES >0 ,"Pos","Neg")
            ego_act@result$direction <- factor(ego_act@result$direction, levels = c("Pos", "Neg"))    
            dotplot(ego_act, split = 'direction', orderBy = 'direction', showCategory=20,title  =paste(cell,condition))+ facet_wrap(. ~ direction)
            ggsave(paste0(cell,"_",condition, "_GSEA_GO.png"), width = 15, height = 15)
            }
            }

        saveRDS(ego_act,paste0(cell,"_",condition, "GSEA_GO.rds"))

    
    
    
      #All = rbind(All,table)
  }


# GO on Reversed 
#Age
library(openxlsx)
library(stringr)
library(EnhancedVolcano)
setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed"))


Files = list.files(pattern = "Aged_Sed-Young_Sed.csv")

Files = Files[!grepl(Files, pattern = "AllCells|M1M2|MacFib|M1|M2|Unk")]
Files


Result_Mixed = list()
for (file in Files){
    Cells =  str_remove(file,"_Aged_Sed-Young_Sed.csv")
    cell =  str_remove(Cells,"MAST")
     Age = read.csv(file)
    Age$cell = cell
     Result_Mixed[[cell]] = Age
    
}

All_Age =do.call(rbind, Result_Mixed)

table(All_Age$fdr < 0.05, All_Age$cell)
table(All_Age$fdr < 0.05 & abs(All_Age$coef) > 0.5)

#Age
Files = list.files(pattern = "Aged_Ex-Aged_Sed.csv")
Files
Files = Files[!grepl(Files, pattern = "AllCells|M1M2|MacFib|M1|M2|Unk")]
Files


Result_Mixed = list()
for (file in Files){
    Cells =  str_remove(file,"_Aged_Ex-Aged_Sed.csv")
    cell =  str_remove(Cells,"MAST")
    Age = read.csv(file)
    Age$cell = cell
    Result_Mixed[[cell]] = Age
    
}

All_Ex =do.call(rbind, Result_Mixed)

table(All_Ex$fdr < 0.05, All_Ex$cell)
table(All_Ex$fdr < 0.05& abs(All_Ex$coef) > 0.5)


All_Ex$state  = paste(All_Ex$primerid, All_Ex$cell)
All_Age$state  = paste(All_Age$primerid, All_Age$cell)
Joint = merge(All_Age, All_Ex, by = 'state', suffixes = c("_Age","_Ex"))


#Reversed Genes

Sig.lvl = 0.05
FC = 0.5
# Directions
AgeUp = (Joint$fdr_Age < Sig.lvl & (Joint$coef_Age) > FC)
AgeDown = (Joint$fdr_Age < Sig.lvl & (Joint$coef_Age) < -FC)

ExUp = (Joint$fdr_Ex < Sig.lvl & (Joint$coef_Ex) > FC)
ExDown = (Joint$fdr_Ex < Sig.lvl & (Joint$coef_Ex) < -FC)

Joint$Direction = "Unset"
Joint$Direction = ifelse( !AgeUp & !AgeDown &!ExUp &!ExDown ,"NSig",Joint$Direction) 
Joint$Direction = ifelse( AgeUp ,"Age DEG (Up)",Joint$Direction) 
Joint$Direction = ifelse( AgeDown ,"Age DEG (Down)",Joint$Direction) 
Joint$Direction = ifelse(AgeUp & ExDown ,"Exercise reverses (Up)",Joint$Direction)
Joint$Direction = ifelse(AgeDown & ExUp,"Exercise reverses (Down)",Joint$Direction)


#df =data.frame((table(Joint$Direction , Joint$cell_Age)))

library(tidyverse)
table(Joint$Direction, Joint$cell_Age)

DT::datatable(Joint)

 setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed/GSEA_1000"))
for (x in unique(Joint$cell_Age)){
    sub = Joint[Joint$cell_Age == x,]
    rev = sub[ sub$Direction == "Exercise reverses (Up)",]
    if(length(unique(rev$primerid_Age)) == 0){next}
ExDecreasing_Go <- enrichGO(gene          = unique(rev$primerid_Age),
                            universe      = unique(sub$primerid_Age),
                            OrgDb         = org.Mm.eg.db,
                            ont           = "ALL",
                            pAdjustMethod = "fdr",minGSSize    = 10,
              maxGSSize    = 1000,
                            pvalueCutoff=1, qvalueCutoff = 1,
                            keyType = 'SYMBOL',
                                   )
dotplot(ExDecreasing_Go)
    ggsave(paste0(x,"_Reversed_AgeUp.pdf"), width = 10)
    write.csv(ExDecreasing_Go@result,paste0(x,"_Reversed_AgeUp.csv"))
}

unique(rev$primerid_Age)

 setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed/GSEA_1000"))
for (x in unique(Joint$cell_Age)){
    sub = Joint[Joint$cell_Age == x,]
    rev = sub[ sub$Direction == "Exercise reverses (Down)",]
  if(length(unique(rev$primerid_Age)) == 0){next}  
ExDecreasing_Go <- enrichGO(gene          = unique(rev$primerid_Age),
                            universe      = unique(sub$primerid_Age),
                            OrgDb         = org.Mm.eg.db,
                            ont           = "ALL",
                            pAdjustMethod = "fdr",minGSSize    = 10,
              maxGSSize    = 1000,
                            pvalueCutoff=1, qvalueCutoff = 1,
                            keyType = 'SYMBOL',
                                   )
dotplot(ExDecreasing_Go)

    ggsave(paste0(x,"_Reversed_AgeDown.pdf"), width = 10)
    write.csv(ExDecreasing_Go@result,paste0(x,"_Reversed_AgeDown.csv"))
}



#Merge together 

library(openxlsx)
library(stringr)
library(EnhancedVolcano)
setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed/GSEA_1000"))
RandomBatch=list.files(pattern = '*GSEA_GO.rds', recursive = TRUE, include.dirs = TRUE)
RandomBatch


AllResults = data.frame(matrix(nrow = 0, ncol = 16))

for (file in RandomBatch){
    setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed/GSEA_1000"))
    load = readRDS(file)
    setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed/GSEA_1000"))
    table = load@result
    i =file
      name = str_remove(i, "GSEA_GO.rds")
      condition = str_extract(name,"Aged.*|Young.*")
      cell = str_remove(name, "_Aged.*|_Young.*")
      if (nrow(table) != 0){
           dotplot(load, split = 'direction', orderBy = 'direction', showCategory=20,title  =paste(cell,condition))+ facet_wrap(. ~ direction)
    ggsave(paste0(cell,"_",condition, "_GSEA_GO.pdf"), width = 15, height = 15)
 
        table$condition = condition
        table$cell = cell
        table$cell_condtion = paste(cell, condition)
        AllResults = rbind(AllResults,table)
    }
}

write.csv(AllResults,"GSEA_ALL_results.csv")

GSEA_Results =list.files(pattern = '*_GSEA_GO.pdf', recursive = TRUE, include.dirs = TRUE)
qpdf::pdf_combine(GSEA_Results, "All_GSEA_Plots.pdf")

DT::datatable(AllResults)

manipulate =load 
manipulate@result = AllResults
manipulate@result$cond_Cell_direction = paste(manipulate@result$cell_condtion,manipulate@result$direction)
dotplot(manipulate, split = 'cond_Cell_direction', orderBy = 'cond_Cell_direction', showCategory=5,title  =paste("ALL"))+ facet_wrap(. ~ cond_Cell_direction)
ggsave("All_GSEA.pdf", width = 40 ,height = 40)

  ego_act <- gseGO(geneList     = geneList,
              OrgDb        = org.Mm.eg.db,
              ont          = "All",
              minGSSize    = 10,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
                         keyType = "SYMBOL",
              verbose      = FALSE)
     
ego_act@result$direction = ifelse(ego_act@result$NES >0 ,"Pos","Neg")
ego_act@result$direction <- factor(ego_act@result$direction, levels = c("Pos", "Neg"))

dotplot(ego_act, split = 'direction', orderBy = 'direction')+ facet_wrap(~ direction)

head(ego_act@result)

# Get reversed for each cellType 
library(openxlsx)
library(stringr)
library(EnhancedVolcano)
setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed"))
RandomBatch=list.files(pattern = 'MAST.*csv', recursive = TRUE, include.dirs = TRUE)
#RandomBatch = RandomBatch[grepl(RandomBatch, pattern = "M1|M2|Macrop")]
RandomBatch

library(org.Mm.eg.db)
AnnotationDbi::species(org.Mm.eg.db)
library(clusterProfiler)
i = RandomBatch[1]
FC = 0.2
P = 0.05
All = data.frame(matrix(nrow = 0,ncol = 15))
for (i in RandomBatch) {
 
  names =i
 
    
    print(i)
    setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed"))
    table =  read.csv(i)
    setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed/GO"))
      name = str_remove(i, "MAST")
      name = str_remove(name, ".csv")
      condition = str_extract(name,"Aged.*|Young.*")
      cell = str_remove(name, "_Aged.*|_Young.*")
      table$cell = cell
      table$condition = condition
   
    
    AnnotationDbi::species(org.Mm.eg.db)
        Order = table[order(table$coef, decreasing = TRUE ),]
    ego_supr = NULL
    ego_act = NULL
    testFr = NULL
   
    
    
      All = rbind(All,table)
  }


table(AllResults$cell)
table(AllResults$condition)

AgeEffect = AllResults[ AllResults$condition == "Aged_Sed-Young_Sed",]
ExEffect = AllResults[ AllResults$condition == "Aged_Ex-Aged_Sed",]

head(sort(table(AgeEffect$Description),decreasing = TRUE),20)
head(sort(table(ExEffect$Description),decreasing = TRUE),20)



DT::datatable(AgeEffect)

DT::datatable(ExEffect)

library(openxlsx)
library(stringr)
library(EnhancedVolcano)
setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed"))
RandomBatch=list.files(pattern = 'MAST.*csv', recursive = TRUE, include.dirs = TRUE)
#RandomBatch = RandomBatch[grepl(RandomBatch, pattern = "_Aged_Ex-Aged_Sed|Aged_Sed-Young_Sed")]
RandomBatch = RandomBatch[!grepl(RandomBatch, pattern = "Global")]
RandomBatch

library(org.Mm.eg.db)
AnnotationDbi::species(org.Mm.eg.db)
library(clusterProfiler)
i = RandomBatch[1]
FC = 0.5
P = 0.05
All = data.frame(matrix(nrow = 0,ncol = 15))
print("running")
for (i in RandomBatch) {
       print(i)
      names =i
      name = str_remove(i, "MAST")
      name = str_remove(name, ".csv")
      condition = str_extract(name,"Aged.*|Young.*")
      cell = str_remove(name, "_Aged.*|_Young.*")
    
   
    
     setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed/GSEA_0.5"))
    
   if (file.exists(paste0(cell,"_",condition, "GSEA_GO_0.5.rds"))) {
       next
    }
    setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed"))
     print(i)
    table =  read.csv(i)
    setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed/GSEA_0.5"))
     
      table$cell = cell
      table$condition = condition
   
    
    AnnotationDbi::species(org.Mm.eg.db)
     geneList = table[,'coef']
    names(geneList) = as.character(table$primerid)
    geneList = sort(geneList, decreasing = TRUE)
    
    ego_supr = NULL
    ego_act = NULL
    testFr = NULL
     ego_act <- gseGO(geneList     = geneList,
              OrgDb        = org.Mm.eg.db,
              ont          = "All",
              minGSSize    = 10,
              maxGSSize    = 500,
              pvalueCutoff = 0.5,
                         keyType = "SYMBOL",
              verbose      = FALSE)
     
            
   
        
             
            if (!is.null(ego_act)){
            if (nrow(ego_act@result != 0 )){
            ego_act@result$direction = ifelse(ego_act@result$NES >0 ,"Pos","Neg")
            ego_act@result$direction <- factor(ego_act@result$direction, levels = c("Pos", "Neg"))    
            dotplot(ego_act, split = 'direction', orderBy = 'direction', showCategory=20,title  =paste(cell,condition))+ facet_wrap(. ~ direction)
            ggsave(paste0(cell,"_",condition, "_GSEA_GO_0.5.png"), width = 15, height = 15)
            }
            }

        saveRDS(ego_act,paste0(cell,"_",condition, "GSEA_GO_0.5.rds"))

    
    
    
      #All = rbind(All,table)
  }


library(openxlsx)
library(stringr)
library(EnhancedVolcano)
library(clusterProfiler)
setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed/GSEA_0.5"))
RandomBatch=list.files(pattern = 'GSEA_GO_0.5.rds', recursive = TRUE, include.dirs = TRUE)
#RandomBatch = RandomBatch[grepl(RandomBatch, pattern = "_Aged_Ex-Aged_Sed|Aged_Sed-Young_Sed")]
RandomBatch
AllResults = data.frame(matrix(nrow = 0, ncol = 16))

for (file in RandomBatch){
    setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed/GSEA_0.5"))
    load = readRDS(file)
    setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed/GSEA_0.5"))
    table = load@result
    i =file
      name = str_remove(i, "GSEA_GO_0.5.rds")
      condition = str_extract(name,"Aged.*|Young.*")
      cell = str_remove(name, "_Aged.*|_Young.*")
      if (nrow(table) != 0){
           dotplot(load, split = 'direction', orderBy = 'direction', showCategory=20,title  =paste(cell,condition))+ facet_wrap(. ~ direction)
    ggsave(paste0(cell,"_",condition, "_GSEA_GO.pdf"), width = 15, height = 15)
 
        table$condition = condition
        table$cell = cell
        table$cell_condtion = paste(cell, condition)
        AllResults = rbind(AllResults,table)
    }
}

AllResults$absNES = abs(AllResults$NES)

AgeEffect = AllResults[ AllResults$condition == "Aged_Sed-Young_Sed",]
ExEffect = AllResults[ AllResults$condition == "Aged_Ex-Aged_Sed",]
allterms = c(AgeEffect$Description,ExEffect$Description)
head(allterms)
head(sort(table(allterms), decreasing = TRUE),20)

head(sort(table(AgeEffect$Description),decreasing = TRUE),20)
head(sort(table(ExEffect$Description),decreasing = TRUE),20)

DT::datatable(AgeEffect, filter = 'top')



DT::datatable(ExEffect, filter = 'top')

ExSplicing =  ExEffect[ ExEffect$ID %in% c('GO:0043484','GO:0055006','GO:0043292','GO:0060537'),]
ExSplicing = na.omit(ExSplicing)
AgeSplicing =  AgeEffect[ AgeEffect$ID %in% c('GO:0043484','GO:0055006','GO:0043292','GO:0060537'),]
AgeSplicing = na.omit(AgeSplicing)


AgeSplicing= AgeSplicing[ !AgeSplicing$cell %in% c( 'M1M2',     'Unknown_3',        'M2', 'M1',          'M_Ex'),]
#
ExSplicing= ExSplicing[ !ExSplicing$cell %in% c( 'M1M2',     'Unknown_3',        'M2',   'M1',         'M_Ex'),]

levels =  c(unique(c(AgeSplicing$cell,ExSplicing$cell)))
levels
AgeSplicing$cell = factor(AgeSplicing$cell, levels =levels)
ExSplicing$cell = factor(ExSplicing$cell, levels = levels)                          

library(DOSE)
  ggplot(ExSplicing, aes_string(x='cell', y="Description", size='pvalue', color='NES')) +
        geom_point() +
        scale_color_continuous(low="blue", high="red", name = 'Enrichment', limits = c(-2.5, 2.5)) +
        scale_y_discrete() +
        ylab(NULL) + ggtitle('Exercise Effect') + theme_dose() +
        scale_size(range=c(3, 8)) + theme( axis.text.x = element_text(angle = 45, vjust = 1 , hjust = 1))


  ggplot(AgeSplicing, aes_string(x='cell', y="Description", size='pvalue', color='NES')) +
        geom_point() +
        scale_color_continuous(low="blue", high="red", name = 'Enrichment',limits = c(-2.5, 2.5)) +
        scale_y_discrete() +

        ylab(NULL) + ggtitle('Age Effect') + theme_dose() +
        scale_size(range=c(3, 8)) + theme( axis.text.x = element_text(angle = 45, vjust = 1 , hjust = 1))




dotplot_internal <- function(object, x = "NES", color = "p.adjust",
                             showCategory=10, size=NULL, split = NULL,
                             font.size=12, title = "", orderBy="x",
                             decreasing=TRUE, label_format = 30) {

    colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))
    if (x == "geneRatio" || x == "GeneRatio") {
        x <- "GeneRatio"
        if (is.null(size))
            size <- "fold_enrichment"
    }

    df <- fortify.internal(object, showCategory = showCategory, split=split)
    ## already parsed in fortify
    ## df$GeneRatio <- parse_ratio(df$GeneRatio)

    if (orderBy !=  'x' && !orderBy %in% colnames(df)) {
        message('wrong orderBy parameter; set to default `orderBy = "x"`')
        orderBy <- "x"
    }

    if (orderBy == "x") {
        df <- dplyr::mutate(df, x = eval(parse(text=x)))
    }

    label_func <- default_labeller(label_format)
    if(is.function(label_format)) {
        label_func <- label_format
    }

    idx <- order(df[[orderBy]], decreasing = decreasing)
    df$Description <- factor(df$Description,
                          levels=rev(unique(df$Description[idx])))
   # return(df)  
    ggplot(df, aes_string(x=x, y="Description", size=size, color=colorBy)) +
        geom_point() +
        scale_color_continuous(low="red", high="blue", name = color,
            guide=guide_colorbar(reverse=TRUE)) +
        scale_y_discrete(labels = label_func) +
        ylab(NULL) + ggtitle(title) + theme_dose(font.size) +
        scale_size(range=c(3, 8))
    
}
default_labeller <- function(n) {
    function(str){
        str <- gsub("_", " ", str)
        ep_str_wrap(str, n)
    }
}

ep_str_wrap <- function(string, width) {
    x <- gregexpr(' ', string)
    vapply(seq_along(x),
           FUN = function(i) {
               y <- x[[i]]
               n <- nchar(string[i])
               len <- (c(y,n) - c(0, y)) ## length + 1
               idx <- len > width
               j <- which(!idx)
               if (length(j) && max(j) == length(len)) {
                   j <- j[-length(j)]
               }
               if (length(j)) {
                   idx[j] <- len[j] + len[j+1] > width
               }
               idx <- idx[-length(idx)] ## length - 1
               start <- c(1, y[idx] + 1)
               end <- c(y[idx] - 1, n)
               words <- substring(string[i], start, end)
               paste0(words, collapse="\n")
           },
           FUN.VALUE = character(1)
    )
}
fortify.internal <- function(model, data, showCategory=5, by = "Count",
                             order=FALSE, drop=FALSE, split=NULL, ...) {
    res <- as.data.frame(model)
    res <- res[!is.na(res$Description), ]

    topN <- function(res, showCategory) {
        if ( is.numeric(showCategory) ) {
            if ( showCategory <= nrow(res) ) {
                res <- res[1:showCategory,]
            }
        } else { ## selected categories
            res <- res[res$Description %in% showCategory,]
        }
        return(res)
    }

    if (is.null(split)) {
        res <- topN(res, showCategory)
    } 

    res$Description <- factor(res$Description,
                              levels=rev(unique(res$Description)))

    return(res)
}



ExEffect$join = paste(ExEffect$cell, ExEffect$ID)
AgeEffect$join = paste(AgeEffect$cell, AgeEffect$ID)


ExEffect = ExEffect[ExEffect$join %in% AgeEffect$join,]
AgeEffect = AgeEffect[AgeEffect$join %in% ExEffect$join,]
ExEffect = ExEffect[order(ExEffect$join),]
AgeEffect = AgeEffect[order(AgeEffect$join),]
dim(ExEffect)
dim(AgeEffect)
head(ExEffect$join)
head(AgeEffect$join)

ExEffectSearch = ExEffect
colnames(ExEffectSearch) =paste0(colnames(ExEffectSearch), ".Ex")
AgeEffectSearch = AgeEffect
colnames(AgeEffectSearch) =paste0(colnames(AgeEffectSearch), ".Age")
Joined = cbind(ExEffectSearch,AgeEffectSearch)
colnames(Joined)
Joined$diff = abs(Joined$NES.Age - Joined$NES.Ex)
#Joined[,c('ID.Ex','Description.Ex','NES.Ex','pvalue.Ex','NES.Age')]
#Joined[,c('ID.Ex','Description.Ex','NES.Ex','pvalue.Ex','NES.Age','pvalue.Age','cell.Ex')]

# Selected GO terms are defined explicitly in `Keep` below.
DT::datatable(Joined[,c('diff','Description.Ex','NES.Ex','pvalue.Ex','NES.Age','pvalue.Age','cell.Ex','join.Age')], filter = 'top')


'Cardiomyocyte GO:0015078' %in% Joined$join.Age
Keep = c(
'Cardiomyocyte GO:0015078',
'AllMacrophage GO:2000648',
'Coronary_EC GO:0007623',
'Endocardial_EC GO:0060419',
'Fibroblast GO:0043292',
'Pericytes GO:0055006',
'T_cell GO:0002708',
'B_cell GO:0008380',
'Lymphatic_EC GO:0061077')
head(Joined$join.Age)
PlotDiff = Joined[ Joined$join.Age %in% Keep,]
PlotDiff

library(DOSE)
getwd()
  ggplot(PlotDiff, aes_string(x='cell.Ex', y="Description.Ex", color='NES.Ex')) +
        geom_point() +
        scale_color_continuous(low="blue", high="red", name = 'Enrichment', limits = c(-2.5, 2.5)) +
        scale_y_discrete() +
        ylab(NULL) + ggtitle('Exercise Effect') + theme_dose() +
        scale_size(range=c(3, 8)) + theme( axis.text.x = element_text(angle = 45, vjust = 1 , hjust = 1))
ggsave("IntersetingEx.pdf", width = 10)

  ggplot(PlotDiff, aes_string(x='cell.Age', y="Description.Age", color='NES.Age')) +
        geom_point() +
        scale_color_continuous(low="blue", high="red", name = 'Enrichment',limits = c(-2.5, 2.5)) +
        scale_y_discrete() +

        ylab(NULL) + ggtitle('Age Effect') + theme_dose() +
        scale_size(range=c(3, 8)) + theme( axis.text.x = element_text(angle = 45, vjust = 1 , hjust = 1))
ggsave("IntersetingAge.pdf", width = 10)

PlotDiffAge = PlotDiff[,c('cell.Age','Description.Age','NES.Age'),]
PlotDiffAge$Cond = 'ASed-YSed'
PlotDiffEx = PlotDiff[,c('cell.Ex','Description.Ex','NES.Ex'),]
PlotDiffEx$Cond = 'AEx-ASed'

colnames(PlotDiffAge) =c('Cell','Description','NES','Comparison')
colnames(PlotDiffEx) =c('Cell','Description','NES','Comparison')



# Moved up from below: these axis orders are needed by the two plots that follow (they were defined after use).
cellLevel = c('AllMacrophage','B_cell','Cardiomyocyte','Coronary_EC','Endocardial_EC','Fibroblast','Pericytes','T_cell','Lymphatic_EC')
DesLevel = c('positive regulation of stem cell proliferation',
             'RNA splicing',
             'proton transmembrane transporter activity',
             'circadian rhythm',
             'heart growth',
             'contractile fiber',
             'cardiac cell development',
             'positive regulation of lymphocyte mediated immunity','chaperone-mediated protein folding')

PlotDiffLong = rbind(PlotDiffAge,PlotDiffEx)
PlotDiffLong = PlotDiffLong[order(PlotDiffLong$NES),]
PlotDiffLong$Description = factor(PlotDiffLong$Description, levels = unique(PlotDiffLong$Description))
ggplot(PlotDiffLong, aes_string(x='NES', y="Description", color='Comparison')) +
        geom_point(size = 5) +
          scale_y_discrete(limits = DesLevel) +
        ylab(NULL) + ggtitle('Reversed Terms Effect') + theme_dose()
ggsave("IntersetingBoth_Term.pdf", width = 10)

PlotDiffLong = rbind(PlotDiffAge,PlotDiffEx)
ggplot(PlotDiffLong, aes_string(x='NES', y="Cell", color='Comparison')) +
        geom_point(size = 5) +
          scale_y_discrete( limits =cellLevel ) +

        ylab(NULL) + ggtitle('Reversed Terms Effect') + theme_dose()
ggsave("IntersetingBoth_Cell.pdf", width = 10)

unique(PlotDiffLong$Description)
# (cellLevel and DesLevel are now defined above, before the two plots that use them.)

PlotDiffLong
