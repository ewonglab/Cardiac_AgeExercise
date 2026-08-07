# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Run GREAT enrichment for shared aging/exercise ATAC changes.
# Steps: 1) Load differential accessibility results; 2) define overlapping contrasts; 3) submit/test region sets and export enrichment summaries.
# =======================================
# Source reformatted from snATAC_AgeEx/mouse_BG_GREAT.ipynb.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))



setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2025-06-28.GREAT_local_mousebg"))

DAP = read.csv(publication_path("HEART_SNRNA_ROOT", "analysis/2023-10-11.AllMacs/All_DAP.csv"))
head(DAP)
colnames(DAP)
#colnames(DAP) = c("ID","peak_coord","pvalue","log2FoldChange","pct.1","pct.2","p_val_adj","cell","condition")
boxplot(DAP$log2FoldChange)
DAP$state = paste(DAP$peak_coord, DAP$cell)

Age = DAP[ DAP$cond == 'Aged_Sed - Young_Sed',]
Ex = DAP[ DAP$cond == 'Aged_Ex - Aged_Sed',]

table(Age$state %in% Ex$state)

Age =Age[Age$state %in% Ex$state,]
Age =Age[order(Age$state),]
#Age = na.omit(Age)
#table(Age$cell)
Ex =Ex[Ex$state %in% Age$state,]
Ex =Ex[order(Ex$state),]
#Ex = na.omit(Ex)
#table(Ex$cell)
colnames(Ex) = paste0(colnames(Ex),"_Ex")
colnames(Age) = paste0(colnames(Age),"_Age")

table(Ex$state == Age$state)
Age_Ex_ATAC = cbind(Ex, Age)
sigLvl = 0.1

FC = 0
Joint = Age_Ex_ATAC

FC = 0
P_thre = .1
Sig.lvl = .1
# Directions
Age_Ex_ATAC$Direction = "NSig"
Age_Ex_ATAC$Direction = ifelse( Age_Ex_ATAC$pvalue_Age < P_thre & abs(Age_Ex_ATAC$log2FoldChange_Age) > FC &  (Age_Ex_ATAC$pvalue_Ex > P_thre | abs(Age_Ex_ATAC$log2FoldChange_Ex) < FC) ,"Age DEG" ,Age_Ex_ATAC$Direction) 
Age_Ex_ATAC$Direction = ifelse( Age_Ex_ATAC$pvalue_Ex < P_thre & abs(Age_Ex_ATAC$log2FoldChange_Ex) > FC &  (Age_Ex_ATAC$pvalue_Age > P_thre | abs(Age_Ex_ATAC$log2FoldChange_Age) < FC) ,"Ex DEG" ,Age_Ex_ATAC$Direction) 
Age_Ex_ATAC$Direction = ifelse(Age_Ex_ATAC$pvalue_Age < P_thre & Age_Ex_ATAC$pvalue_Ex < P_thre & Age_Ex_ATAC$log2FoldChange_Age > FC & Age_Ex_ATAC$log2FoldChange_Ex > FC ,"Exercise exacerbates AgeUp ",Age_Ex_ATAC$Direction) 
Age_Ex_ATAC$Direction = ifelse(Age_Ex_ATAC$pvalue_Age < P_thre & Age_Ex_ATAC$pvalue_Ex < P_thre & Age_Ex_ATAC$log2FoldChange_Age < -FC & Age_Ex_ATAC$log2FoldChange_Ex < -FC ,"Exercise exacerbates AgeDown",Age_Ex_ATAC$Direction)
Age_Ex_ATAC$Direction = ifelse(Age_Ex_ATAC$pvalue_Age < P_thre & Age_Ex_ATAC$pvalue_Ex < P_thre & Age_Ex_ATAC$log2FoldChange_Age < -FC & Age_Ex_ATAC$log2FoldChange_Ex > FC ,"Exercise reverses AgeDown",Age_Ex_ATAC$Direction)
Age_Ex_ATAC$Direction = ifelse(Age_Ex_ATAC$pvalue_Age < P_thre & Age_Ex_ATAC$pvalue_Ex < P_thre & Age_Ex_ATAC$log2FoldChange_Age > FC & Age_Ex_ATAC$log2FoldChange_Ex < -FC ,"Exercise reverses AgeUp",Age_Ex_ATAC$Direction)

Age_Ex_ATAC$AgeDirection = "NSig"
Age_Ex_ATAC$AgeDirection = ifelse( Age_Ex_ATAC$pvalue_Age < P_thre & abs(Age_Ex_ATAC$log2FoldChange_Age) > FC &  (Age_Ex_ATAC$pvalue_Ex > P_thre | abs(Age_Ex_ATAC$log2FoldChange_Ex) < FC) ,"Age DEG" ,Age_Ex_ATAC$AgeDirection) 
Age_Ex_ATAC$ExDirection = "NSig"
Age_Ex_ATAC$ExDirection = ifelse( Age_Ex_ATAC$pvalue_Ex < P_thre & abs(Age_Ex_ATAC$log2FoldChange_Ex) > FC &  (Age_Ex_ATAC$pvalue_Age > P_thre | abs(Age_Ex_ATAC$log2FoldChange_Age) < FC) ,"Ex DEG" ,Age_Ex_ATAC$ExDirection) 
df =data.frame((table(Age_Ex_ATAC$Direction , Age_Ex_ATAC$cell_Age)))

library(tidyverse)

Age_Ex_ATAC <- Age_Ex_ATAC %>%
  tidyr::separate(peak_coord_Ex, into = c("chr", "start", "end"), sep = "-", convert = TRUE)

# Create GRanges object
 #Reversed Genes
   Joint =  Age_Ex_ATAC
    # Directions,
    AgeUp = (Joint$pvalue_Age < Sig.lvl & (Joint$log2FoldChange_Age) > FC)
    AgeDown = (Joint$pvalue_Age < Sig.lvl & (Joint$log2FoldChange_Age) < -FC)
    
    ExUp = (Joint$pvalue_Ex < Sig.lvl & (Joint$log2FoldChange_Ex) > FC)
    ExDown = (Joint$pvalue_Ex < Sig.lvl & (Joint$log2FoldChange_Ex) < -FC)
    



library(rGREAT)
library(GenomicRanges)
library(dplyr)
setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2025-06-28.GREAT_local_mousebg"))
options(scipen = 999)

# Parameters
Sig.lvl <- 0.1
FC <- 0

# Change this if your cell type column has a different name
cell_type_col <- "cell_Age"

# Helper to convert 'chr:start-end' to GRanges
make_gr_safe <- function(coords_vec) {
  coords_vec <- unique(coords_vec[!is.na(coords_vec)])
  
  if(length(coords_vec) == 0) return(NULL)
  
  tmp <- do.call(rbind, strsplit(coords_vec, "[:-]"))
  
  GRanges(
    seqnames = tmp[, 1],
    ranges = IRanges(
      start = as.numeric(tmp[, 2]),
      end = as.numeric(tmp[, 3])
    )
  )
}

# Define comparisons
comparisons <- list(
  "Age Up" = function(df) {
    df %>% filter(pvalue_Age < Sig.lvl, log2FoldChange_Age > FC)
  },
  "Age Down" = function(df) {
    df %>% filter(pvalue_Age < Sig.lvl, log2FoldChange_Age < -FC)
  },
  "Ex Up" = function(df) {
    df %>% filter(pvalue_Ex < Sig.lvl, log2FoldChange_Ex > FC)
  },
  "Ex Down" = function(df) {
    df %>% filter(pvalue_Ex < Sig.lvl, log2FoldChange_Ex < -FC)
  }
)

great_results <- list()

cell_types <- unique(Joint[[cell_type_col]])

for(cell in cell_types) {
  
  message("Processing cell type: ", cell)
  
  cell_df <- Joint %>%
    filter(.data[[cell_type_col]] == cell)
  
  for(comp_name in names(comparisons)) {
    
    fg_df <- comparisons[[comp_name]](cell_df)
    
    if(nrow(fg_df) > 10) {
      
      cond_name <- paste(cell, comp_name, sep = "_")
      message("Running Genome-wide GREAT for: ", cond_name)
      
      fg_gr <- make_gr_safe(fg_df$peak_coord_Age)
      if (file.exists( paste0(cell,"_MF_",comp_name,".rds"))){next}
      if(is.null(fg_gr) || length(fg_gr) <= 10) {
        message("Skipping ", cond_name, ": not enough valid coordinates")
        next
      }
        fg_df
      bg_gr = make_gr_safe(cell_df$peak_coord_Age)
      job <- great( 
        fg_gr,
        tss_source = "mm10",
        gene_sets = c("GO:CC"),
        background = bg_gr
      )
      
      res <- getEnrichmentTable(job)
      
      res$CellType <- cell
      res$Comparison <- comp_name
      res$Condition <- cond_name
      res$N_regions <- length(fg_gr)
      
      great_results[[cond_name]] <- res
      saveRDS(great_results, paste0(cell,"_CC_",comp_name,".rds"))
        job <- great( 
        fg_gr,
        tss_source = "mm10",
        gene_sets = c("GO:BP"),
        background = bg_gr
      )
      
      res <- getEnrichmentTable(job)
      
      res$CellType <- cell
      res$Comparison <- comp_name
      res$Condition <- cond_name
      res$N_regions <- length(fg_gr)
      
      great_results[[cond_name]] <- res
      saveRDS(great_results, paste0(cell,"_BP_",comp_name,".rds"))
        job <- great( 
        fg_gr,
        tss_source = "mm10",
        gene_sets = c("GO:MF"),
        background = bg_gr
      )
      
      res <- getEnrichmentTable(job)
      
      res$CellType <- cell
      res$Comparison <- comp_name
      res$Condition <- cond_name
      res$N_regions <- length(fg_gr)
      
      great_results[[cond_name]] <- res
      saveRDS(great_results, paste0(cell,"_MF_",comp_name,".rds"))
    } else {
      message("Skipping ", cell, " ", comp_name, ": only ", nrow(fg_df), " regions")
    }
  }
}

final_great_df <- bind_rows(great_results)

write.csv(
  final_great_df,
  "Joint_celltype_GREAT_GO_BP_cellwide.csv",
  row.names = FALSE
)

saveRDS(
  final_great_df,
  "Joint_celltype_GREAT_GO_BP_cellwide.rds"
)

library(rGREAT)
library(GenomicRanges)
library(dplyr)
great

library(tidyverse)

in_dir = getwd()
files <- list.files(
  in_dir,
  pattern = ".rds",
  full.names = TRUE
)

parse_filename <- function(file) {
  fname <- basename(file)

  m <- stringr::str_match(
    fname,
    "^(.*)_(Age Down|Age Up|Ex Down|Ex Up|Exercise reverses \\(AgeDown\\)|Exercise reverses \\(AgeUp\\)|Exercise exacerbates \\(AgeDown\\)|Exercise exacerbates \\(AgeUp\\))_GO:(BP|CC|MF)_GreatJob_local\\.rds$"
  )

  tibble(
    file = fname,
    cell_type = m[, 2],
    comparison = m[, 3],
    ontology = paste0("GO:", m[, 4])
  )
}

extract_great_table <- function(file) {
  meta <- parse_filename(file)

  obj <- readRDS(file)

  tab <- tryCatch({
    if (inherits(obj, "GreatJob")) {
      # For rGREAT GreatJob objects
      rGREAT::getEnrichmentTables(obj)
    } else if (is.data.frame(obj)) {
      obj
    } else if (is.list(obj)) {
      # If the RDS already stores one or more result tables
      bind_rows(obj, .id = "table_name")
    } else {
      stop("Unsupported RDS object type")
    }
  }, error = function(e) {
    warning("Failed to extract: ", basename(file), " | ", conditionMessage(e))
    return(NULL)
  })

  if (is.null(tab)) return(NULL)

  # rGREAT::getEnrichmentTables() may return a list of tables
  if (is.list(tab) && !is.data.frame(tab)) {
    tab <- bind_rows(tab, .id = "great_category")
  }

  bind_cols(meta, as_tibble(tab))
}

all_great_results <- map_dfr(files, extract_great_table)

write_csv(
  all_great_results,
  file.path(in_dir, "all_GREAT_results_combined.csv")
)

saveRDS(
  all_great_results,
  file.path(in_dir, "all_GREAT_results_combined.rds")
)

all_great_results

# =============================================================================
# 2. USER'S RNA CODE (Generates RNA_oneDirection and RNA_otherDirection)
# =============================================================================
theme_mytheme <- theme(
  plot.background = element_blank(),
  plot.title = element_text(size = 16, face = "bold"),
  plot.subtitle = element_text(size = 14),
  axis.text = element_text(size = 12, color = "black"),
  axis.title = element_text(size = 14, color = "black"),
  axis.line = element_line(color = "black", size = 0.5),
  axis.ticks = element_line(color = "black", size = 0.5),
  legend.text = element_text(size = 12),
  legend.position = "right",
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.margin = margin(1, 1, 1, 1, "cm")
)

AllResults <- read.csv(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed/GSEA_1000/GSEA_ALL_results.csv"))
AllResults$cell <- ifelse(AllResults$cell == 'AllMacrophage', "Macrophage", AllResults$cell)
AllResults$cond_Cell_direction <- paste(AllResults$cell_condtion, AllResults$direction)

Keep <- read.csv(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed/GSEA_1000/GO_Doubles.csv"))
Keep$X...Term <- str_remove(Keep$X...Term, "'")
Keep <- Keep[is.na(Keep$Remove),]
Keep <- Keep[order(Keep$Class),]
Keep$X...Term <- stringr::str_remove(Keep$X...Term, "'")
Keep$X...Term <- factor(Keep$X...Term, levels = unique(Keep$X...Term))

LargeCell <- AllResults[!grepl(AllResults$cell, pattern = "Unk|M1|M2|HSPC|Dendritic|Neuro|Smooth|M_Ex|Meso"),]
LargeCell <- LargeCell[LargeCell$p.adjust < 0.2,]

# Panel 1 RNA Terms (Age Up / Ex Down)
KeppPlot_One <- LargeCell[LargeCell$Description %in% Keep$X...Term,]
KeppPlot_One$Description <- factor(KeppPlot_One$Description, levels = levels(Keep$X...Term))
KeppPlot_One <- KeppPlot_One[KeppPlot_One$condition %in% c("Aged_Ex-Aged_Sed", "Aged_Sed-Young_Sed"),]
KeppPlot_One$split <- paste(KeppPlot_One$condition, KeppPlot_One$direction)
KeppPlot_One <- KeppPlot_One[KeppPlot_One$split %in% c('Aged_Sed-Young_Sed Pos', 'Aged_Ex-Aged_Sed Neg'),]
RNA_oneDirection <- KeppPlot_One

# Panel 2 RNA Terms (Age Down / Ex Up)
KeppPlot_Other <- LargeCell[LargeCell$Description %in% Keep$X...Term,]
KeppPlot_Other$Description <- factor(KeppPlot_Other$Description, levels = levels(Keep$X...Term))
KeppPlot_Other <- KeppPlot_Other[KeppPlot_Other$condition %in% c("Aged_Ex-Aged_Sed", "Aged_Sed-Young_Sed"),]
KeppPlot_Other <- KeppPlot_Other[KeppPlot_Other$Description %in% names(table(KeppPlot_Other$Description))[(table(KeppPlot_Other$Description) > 1)],]
KeppPlot_Other$split <- paste(KeppPlot_Other$condition, KeppPlot_Other$direction)
KeppPlot_Other <- KeppPlot_Other[KeppPlot_Other$split %in% c('Aged_Sed-Young_Sed Neg', 'Aged_Ex-Aged_Sed Pos'),]
RNA_otherDirection <- KeppPlot_Other

# EXTRACT EXACT Y-AXIS LIMITS FROM RNA DATA
y_limits_panel1 <- as.character(unique(Keep$X...Term[Keep$X...Term %in% RNA_oneDirection$Description]))
y_limits_panel2 <- as.character(unique(Keep$X...Term[Keep$X...Term %in% RNA_otherDirection$Description]))


keep = c(
  "immune response",
  "defense response",
  "extracellular structure organization",
  "extracellular matrix organization",
  "RNA splicing",
  "RNA processing",
  "RNA metabolic process",
  "RNA binding",
  "myofibril",
  "muscle system process",
  "muscle adaptation",
  "heart contraction",
  "cardiac cell development",
  "regulation of cardiac tissue growth",
  "muscle tissue development",
  "muscle cell differentiation",
  "heart growth"
)
y_limits_panel1 = y_limits_panel1[y_limits_panel1 %in% keep]
y_limits_panel2 = y_limits_panel2[y_limits_panel2 %in% keep]
terms_panel1_lower <- tolower(y_limits_panel1)
terms_panel2_lower <- tolower(y_limits_panel2)



mouse_bg <- readRDS("all_GREAT_results_combined.rds")
head(mouse_bg)
unique(mouse_bg$CellType)

# =============================================================================
# 3. ATAC DATA PREPARATION AND PLOTTING
# =============================================================================
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)

mouse_bg <- readRDS("all_GREAT_results_combined.rds") %>% filter(Comparison %in% c("Age Down","Age Up","Ex Down","Ex Up")) %>%
  dplyr::select(Comparison, Description = description, p.adjust = p_adjust, fold_enrichment,CellType) %>%
  mutate(
    log2_FC_enrichment = log2(fold_enrichment),
    cell = CellType,
    Experiment = str_split(Comparison, " ", simplify = TRUE)[,1],
    direction = str_split(Comparison, " ", simplify = TRUE)[,2],
    split = paste(Experiment, direction) 
  )

# Rename cells

# Filter ATAC data:
# 1. Must be a term we care about
# 2. Must be significant (p < 0.2)
ATAC_PlotData <- mouse_bg %>% 
  filter(tolower(Description) %in% tolower(Keep$X...Term)) %>%
  filter(p.adjust < 0.2)

# Fix capitalization to perfectly match the RNA plots
ATAC_PlotData$Description <- Keep$X...Term[match(tolower(ATAC_PlotData$Description), tolower(Keep$X...Term))]
ATAC_PlotData$Description <- factor(ATAC_PlotData$Description, levels = levels(Keep$X...Term))
#ATAC_PlotData$cell <- factor(ATAC_PlotData$cell, levels = c("Cardiomyocyte", "Fibroblast"))
ATAC_PlotData <- ATAC_PlotData %>% 
filter(cell %in% c('Cardiomyocyte','Pericytes','Fibroblast',
                   'Macrophage','Endocardial_EC','Coronary_EC','B_cell',
                   'Lymphatic_EC','T_cell'))
ATAC_PlotData$cell = factor(ATAC_PlotData$cell, levels = c('Cardiomyocyte','Coronary_EC','Endocardial_EC','Lymphatic_EC','Fibroblast',
                   'Macrophage','B_cell','T_cell','Pericytes'))

ATAC_PlotData <- ATAC_PlotData %>% filter(Description %in% c(
  "immune response",
  "defense response",
  "extracellular structure organization",
  "extracellular matrix organization",
  "RNA splicing",
  "RNA processing",
  "RNA metabolic process",
  "RNA binding",
  "myofibril",
  "muscle system process",
  "muscle adaptation",
  "heart contraction",
  "cardiac cell development",
  "regulation of cardiac tissue growth",
  "muscle tissue development",
  "muscle cell differentiation",
  "heart growth"
))
# Find max log2 enrichment to keep dot sizes uniform
max_enrichment <- max(abs(ATAC_PlotData$log2_FC_enrichment), na.rm=TRUE)

# -----------------------------------------------------------------------------
# PANEL 1: Age Up / Ex Down
# -----------------------------------------------------------------------------
ATAC_Panel1 <- ATAC_PlotData %>% 
  filter(split %in% c('Age Up', 'Ex Down')) %>%
  filter(tolower(Description) %in% terms_panel1_lower) %>%
  mutate(split = factor(split, levels = c('Age Up', 'Ex Down')))

Age_Up_Plot <- ggplot(ATAC_Panel1, aes(x = cell, y = Description, size = abs(log2_FC_enrichment), color = p.adjust)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue", limits = c(0, 0.2)) + 
  scale_size_continuous(limits = c(0, max_enrichment)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  facet_wrap(~split) +
  theme_mytheme +
  #scale_x_discrete(drop = FALSE, limits = c("Cardiomyocyte", "Fibroblast")) +
  # Force the Y-axis to exactly match the RNA Panel 1 terms
  scale_y_discrete(limits = y_limits_panel1) +
scale_x_discrete(limits = c('Cardiomyocyte','Coronary_EC','Endocardial_EC','Lymphatic_EC','Fibroblast',
                   'Macrophage','B_cell','T_cell','Pericytes')) +
  labs(title = "ATAC: Age Up vs Ex Down", x = "", y = "", size = "Log2(Enrichment)", color = "FDR")

# -----------------------------------------------------------------------------
# PANEL 2: Age Down / Ex Up
# -----------------------------------------------------------------------------
ATAC_Panel2 <- ATAC_PlotData %>% 
  filter(split %in% c('Age Down', 'Ex Up')) %>%
  filter(tolower(Description) %in% terms_panel2_lower) %>%
  mutate(split = factor(split, levels = c('Age Down', 'Ex Up')))

Age_Down_Plot <- ggplot(ATAC_Panel2, aes(x = cell, y = Description, size = abs(log2_FC_enrichment), color = p.adjust)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue", limits = c(0, 0.2)) + 
  scale_size_continuous(limits = c(0, max_enrichment)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  facet_wrap(~split) + 
  theme_mytheme +
  #scale_x_discrete(drop = FALSE, limits = c("Cardiomyocyte", "Fibroblast")) +
  # Force the Y-axis to exactly match the RNA Panel 2 terms
  scale_y_discrete(limits = y_limits_panel2) +
scale_x_discrete(limits = c('Cardiomyocyte','Coronary_EC','Endocardial_EC','Lymphatic_EC','Fibroblast',
                   'Macrophage','B_cell','T_cell','Pericytes')) +
  labs(title = "ATAC: Age Down vs Ex Up", x = "", y = "", size = "Log2(Enrichment)", color = "FDR")

# -----------------------------------------------------------------------------
# COMBINE AND SAVE
# -----------------------------------------------------------------------------
combined <- plot_grid(
  Age_Up_Plot, Age_Down_Plot, 
  ncol = 1, 
  align = 'v', 
  rel_heights = c(10, 18) 
)
ggsave(plot = combined, "ATAC_Matched_to_RNA_Filtered_0.2_test_bg.pdf", width = 8, height = 10)
print(combined)



# =============================================================================
# 3. ATAC DATA PREPARATION AND PLOTTING
# =============================================================================
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)

mouse_bg <- readRDS("all_GREAT_results_combined.rds") %>% filter(Comparison %in% c("Age Down","Age Up","Ex Down","Ex Up")) %>%
  dplyr::select(Comparison, Description = description, p.adjust = p_adjust, fold_enrichment,CellType) %>%
  mutate(
    log2_FC_enrichment = log2(fold_enrichment),
    cell = CellType,
    Experiment = str_split(Comparison, " ", simplify = TRUE)[,1],
    direction = str_split(Comparison, " ", simplify = TRUE)[,2],
    split = paste(Experiment, direction) 
  )

# Rename cells

# Filter ATAC data:
# 1. Must be a term we care about
# 2. Must be significant (p < 0.2)
ATAC_PlotData <- mouse_bg %>% 
  filter(tolower(Description) %in% tolower(Keep$X...Term)) %>%
  filter(p.adjust < 0.2)

# Fix capitalization to perfectly match the RNA plots
ATAC_PlotData$Description <- Keep$X...Term[match(tolower(ATAC_PlotData$Description), tolower(Keep$X...Term))]
ATAC_PlotData$Description <- factor(ATAC_PlotData$Description, levels = levels(Keep$X...Term))
ATAC_PlotData <- ATAC_PlotData %>% 
filter(cell %in% c('Cardiomyocyte','Pericytes','Fibroblast',
                   'Macrophage','Endocardial_EC','Coronary_EC','B_cell',
                   'Lymphatic_EC','T_cell'))
ATAC_PlotData$cell = factor(ATAC_PlotData$cell, levels = c('Cardiomyocyte','Coronary_EC','Endocardial_EC','Lymphatic_EC','Fibroblast',
                   'Macrophage','B_cell','T_cell','Pericytes'))


ATAC_PlotData <- ATAC_PlotData %>% filter(Description %in% c(
  "immune response",
  "defense response",
  "extracellular structure organization",
  "extracellular matrix organization",
  "RNA splicing",
  "RNA processing",
  "RNA metabolic process",
  "RNA binding",
  "myofibril",
  "muscle system process",
  "muscle adaptation",
  "heart contraction",
  "cardiac cell development",
  "regulation of cardiac tissue growth",
  "muscle tissue development",
  "muscle cell differentiation",
  "heart growth"
))
# Find max log2 enrichment to keep dot sizes uniform
max_enrichment <- max(abs(ATAC_PlotData$log2_FC_enrichment), na.rm=TRUE)

# -----------------------------------------------------------------------------
# PANEL 1: Age Up / Ex Down
# -----------------------------------------------------------------------------
ATAC_Panel1 <- ATAC_PlotData %>% 
  filter(split %in% c('Age Down', 'Ex Up')) %>%
  filter(tolower(Description) %in% terms_panel1_lower) %>%
  mutate(split = factor(split, levels = c('Age Down', 'Ex Up')))

Age_Up_Plot <- ggplot(ATAC_Panel1, aes(x = cell, y = Description, size = abs(log2_FC_enrichment), color = p.adjust)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue", limits = c(0, 0.2)) + 
  scale_size_continuous(limits = c(0, max_enrichment)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  facet_wrap(~split) +
  theme_mytheme +
  #scale_x_discrete(drop = FALSE, limits = c("Cardiomyocyte", "Fibroblast")) +
  # Force the Y-axis to exactly match the RNA Panel 1 terms
  scale_y_discrete(limits = y_limits_panel1) +
scale_x_discrete(limits = c('Cardiomyocyte','Coronary_EC','Endocardial_EC','Lymphatic_EC','Fibroblast',
                   'Macrophage','B_cell','T_cell','Pericytes')) +
  labs(title = "ATAC: Age Up vs Ex Down", x = "", y = "", size = "Log2(Enrichment)", color = "FDR")

# -----------------------------------------------------------------------------
# PANEL 2: Age Down / Ex Up
# -----------------------------------------------------------------------------
ATAC_Panel2 <- ATAC_PlotData %>% 
  filter(split %in% c('Age Up', 'Ex Down')) %>%
  filter(tolower(Description) %in% terms_panel2_lower) %>%
  mutate(split = factor(split, levels = c('Age Up', 'Ex Down')))

Age_Down_Plot <- ggplot(ATAC_Panel2, aes(x = cell, y = Description, size = abs(log2_FC_enrichment), color = p.adjust)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue", limits = c(0, 0.2)) + 
  scale_size_continuous(limits = c(0, max_enrichment)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  facet_wrap(~split) + 
  theme_mytheme +
  #scale_x_discrete(drop = FALSE, limits = c("Cardiomyocyte", "Fibroblast")) +
  # Force the Y-axis to exactly match the RNA Panel 2 terms
  scale_y_discrete(limits = y_limits_panel2) +
scale_x_discrete(limits = c('Cardiomyocyte','Coronary_EC','Endocardial_EC','Lymphatic_EC','Fibroblast',
                   'Macrophage','B_cell','T_cell','Pericytes')) +
  labs(title = "ATAC: Age Down vs Ex Up", x = "", y = "", size = "Log2(Enrichment)", color = "FDR")

# -----------------------------------------------------------------------------
# COMBINE AND SAVE
# -----------------------------------------------------------------------------
combined <- plot_grid(
  Age_Up_Plot, Age_Down_Plot, 
  ncol = 1, 
  align = 'v', 
  rel_heights = c(10, 18) 
)
ggsave(plot = combined, "ATAC_Matched_to_RNA_Filtered_0.2_test_bg_rev.pdf", width = 8, height = 10)
print(combined)



Keep$X...Term
