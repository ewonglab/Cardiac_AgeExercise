# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Run KEGG enrichment for the reported transcriptional contrasts.
# Steps: 1) Load contrast gene lists; 2) calculate KEGG enrichment; 3) export pathway summaries and plots.
# =======================================
# Source reformatted from snRNA_AgeEx/KEGG.ipynb.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

# Ampk Pathway plot 



keepCell = c('M1M2',
'Pericytes','Fibroblast','Cardiomyocyte',
'HSPC','Coronary_EC','Endocardial_EC','Lymphatic_EC')

cells = c('Cardiomyocyte','Coronary_EC','Endocardial_EC','Lymphatic_EC','Fibroblast','HSPC','M1M2',
'Pericytes')

# Read in MASt results 

# Read in DEG 
library(openxlsx)
library(stringr)
library(EnhancedVolcano)
setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed"))
RandomBatch=list.files(pattern = 'MAST.*csv', recursive = TRUE, include.dirs = TRUE)
RandomBatch = RandomBatch[!grepl(RandomBatch, pattern = "Global")]
#RandomBatch

All = data.frame(matrix(nrow = 0,ncol = 15))
for (i in RandomBatch) {
 
  names =i
 
     
    #print(i)
    setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed"))
    table =  read.csv(i)
    setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-11-06.MAST_Final_Fixed/GO"))
      name = str_remove(i, "MAST")
      name = str_remove(name, ".csv")
      condition = str_extract(name,"Aged.*|Young.*")
      cell = str_remove(name, "_Aged.*|_Young.*")
      table$cell = cell
      table$condition = condition
   
    
      All = rbind(All,table)
  }
All$cell_cond = paste(All$cell,All$condition)
All$state = paste(All$primerid, All$cell)
Age = All[ All$condition == 'Aged_Sed-Young_Sed',]
Ex = All[ All$condition == 'Aged_Ex-Aged_Sed',]
Age =Age[Age$state %in% Ex$state,]
Age =Age[order(Age$state),]

Ex =Ex[Ex$state %in% Age$state,]
Ex =Ex[order(Ex$state),]

colnames(Ex) = paste0(colnames(Ex),".Ex")
colnames(Age) = paste0(colnames(Age),".Age")

Age_Ex = cbind(Ex, Age)
Age_Ex = Age_Ex[rowSums(is.na(Age_Ex)) != ncol(Age_Ex), ]
FC = 0.5
# Directions
Age_Ex$Direction = "NSig"
Age_Ex$Direction = ifelse( Age_Ex$fdr.Age < 0.05 & abs(Age_Ex$coef.Age) > FC &  (Age_Ex$fdr.Ex > 0.05 | abs(Age_Ex$coef.Ex) < FC) ,"Age DEG" ,Age_Ex$Direction) 
Age_Ex$Direction = ifelse( Age_Ex$fdr.Ex < 0.05 & abs(Age_Ex$coef.Ex) > FC &  (Age_Ex$fdr.Age > 0.05 | abs(Age_Ex$coef.Age) < FC) ,"Ex DEG" ,Age_Ex$Direction) 
Age_Ex$Direction = ifelse(Age_Ex$fdr.Age < 0.05 & Age_Ex$fdr.Ex < 0.05 & Age_Ex$coef.Age > FC & Age_Ex$coef.Ex > FC ,"Exercise exacerbates",Age_Ex$Direction) 
Age_Ex$Direction = ifelse(Age_Ex$fdr.Age < 0.05 & Age_Ex$fdr.Ex < 0.05 & Age_Ex$coef.Age < -FC & Age_Ex$coef.Ex < -FC ,"Exercise exacerbates",Age_Ex$Direction)
Age_Ex$Direction = ifelse(Age_Ex$fdr.Age < 0.05 & Age_Ex$fdr.Ex < 0.05 & Age_Ex$coef.Age < -FC & Age_Ex$coef.Ex > FC ,"Exercise reverses",Age_Ex$Direction)
Age_Ex$Direction = ifelse(Age_Ex$fdr.Age < 0.05 & Age_Ex$fdr.Ex < 0.05 & Age_Ex$coef.Age > FC & Age_Ex$coef.Ex < -FC ,"Exercise reverses",Age_Ex$Direction)

Age_Ex$AgeDirection = "NSig"
Age_Ex$AgeDirection = ifelse( Age_Ex$fdr.Age < 0.05 & abs(Age_Ex$coef.Age) > FC &  (Age_Ex$fdr.Ex > 0.05 | abs(Age_Ex$coef.Ex) < FC) ,"Age DEG" ,Age_Ex$AgeDirection) 
Age_Ex$ExDirection = "NSig"
Age_Ex$ExDirection = ifelse( Age_Ex$fdr.Ex < 0.05 & abs(Age_Ex$coef.Ex) > FC &  (Age_Ex$fdr.Age > 0.05 | abs(Age_Ex$coef.Age) < FC) ,"Ex DEG" ,Age_Ex$ExDirection) 
df =data.frame((table(Age_Ex$Direction , Age_Ex$cell.Age)))

library(tidyverse)




theme_mytheme <- theme(
  # Specify the overall plot appearance
  plot.background = element_blank(),
  plot.title = element_text(size = 8, face = "bold"),
  plot.subtitle = element_text(size = 8),
  
  # Specify the axis appearance
  axis.text = element_text(size = 8, color = "black"),
  axis.title = element_text(size = 8, color = "black"),
  axis.line = element_line(color = "black", size = 0.5),
  axis.ticks = element_line(color = "black", size = 0.5),
  
  # Specify the legend appearance
  legend.title = element_blank(),
  legend.text = element_text(size = 8),
  legend.position = "right",
  
  # Specify the panel/grid appearance, no grid line 
  panel.background = element_blank(),#element_rect(fill = "white", color = "black"),
  panel.grid.major = element_blank(),#element_line(color = "lightgray", size = 0.2),
  panel.grid.minor = element_blank(),#element_line(color = "lightgray", size = 0.2),
  
  # Specify the plot margins
  plot.margin = margin(10, 0, 0, 0, "mm"),
  
  #set the ratio of the plot 
  #aspect.ratio=1
)

# 3. Map Entrez IDs → gene symbols
library(KEGGREST)
library(org.Mm.eg.db)
library(dplyr)
library(tidyr)
library(ggplot2)
res_list = list()
pltlist = list()
# 1. Specify the KEGG pathway (replace with your pathway of interest)
ctrl_ids <- c("mmu04010","mmu04151","mmu04150","mmu04310","mmu04330",
              "mmu04340","mmu04630","mmu04064","mmu04020","mmu04350","mmu04152")

more_ctrl_ids <- c(
  "mmu04022","mmu04024","mmu04910",
  "mmu04620","mmu04621","mmu04622","mmu04660",
  "mmu04662","mmu04062","mmu04722","mmu04370"
)

ctrl_ids <- c(ctrl_ids, more_ctrl_ids)
for (x in ctrl_ids){
#x = ctrl_ids[1]
    # 2. Pull Entrez IDs for that pathway
    pathway_info <- keggGet(x)[[1]]
    gene_entries <- pathway_info$GENE
    descs      <- gene_entries[seq(2, length(gene_entries), by = 2)]
    gene_symbols <- sub(";.*", "", descs)

    length(gene_symbols)

    # 4. Subset your FC data (Age_Ex) to only those pathway genes
    #    assumes Age_Ex has columns: cell, primerid.Ex (gene), coef.Age, coef.Ex
    df_fc <- Age_Ex %>%
      filter(primerid.Ex %in% gene_symbols, cell.Age %in% keepCell,) %>%
      dplyr::select(cell.Age, primerid.Ex , coef.Age, coef.Ex)

    length(unique(df_fc$primerid.Ex))
    # 5. Pivot into long form for ggplot
    df_plot <- df_fc %>%
      pivot_longer(
        cols      = c(coef.Age, coef.Ex),
        names_to  = "Coefficient",
        values_to = "FC"
      )
    length(unique(df_plot$primerid.Ex))
    library(dplyr)
    library(ggplot2)

    # 1. Keep only observed FC values
    df_obs <- df_plot %>% 
      filter(!is.na(FC))

    genes_multi <- df_obs %>%
      group_by(primerid.Ex) %>%
      filter(n_distinct(cell.Age) >= 1) %>%
      distinct(primerid.Ex) %>%
      pull(primerid.Ex)

    # 3. Subset to those genes
    df_filtered <- df_obs %>%
      filter(primerid.Ex %in% genes_multi)

    gene_order <- df_filtered %>%
      filter(Coefficient == "coef.Age") %>%
      group_by(primerid.Ex) %>%
      summarise(meanFC = mean(FC, na.rm = TRUE), .groups = "drop") %>%
      arrange(meanFC) %>%
      pull(primerid.Ex)

    # 2) Turn primerid.Ex into an ordered factor
    df_plot2 <- df_filtered %>%
      mutate(primerid.Ex = factor(primerid.Ex, levels = gene_order))

    cells = cells
    plots <- lapply(cells, function(cell) {
      df_sub <- df_filtered %>% filter(cell.Age == cell)
       if (nrow(df_sub) == 0) {
    message("Skipping ", cell, " because no data")
    return(NULL)
      }
      gene_order <- df_sub %>% 
        filter(Coefficient == "coef.Age") %>% 
        arrange(FC) %>% 
        pull(primerid.Ex)
      df_sub <- df_sub %>% 
        mutate(primerid.Ex = factor(primerid.Ex, levels = gene_order))



     prop_opposite <- df_sub %>%
      # restrict to genes with both measurements
      pivot_wider(
        id_cols     = c(cell.Age, primerid.Ex),
        names_from  = Coefficient,
        values_from = FC
      ) %>%
      filter(!is.na(coef.Age) & !is.na(coef.Ex)) %>%
      group_by(cell.Age) %>%
      summarise(
        prop_opposite = mean(coef.Age * coef.Ex < 0),
           n_opposite = sum(coef.Age * coef.Ex < 0),
        n_genes       = n(),
          prop_total = sum(coef.Age * coef.Ex < 0)/length(unique(gene_symbols)),
          prop_total_n = length(unique(gene_symbols)),
          score_raw = sum(coef.Age * coef.Ex),
          score_raw_norm = sum(coef.Age * coef.Ex)* (n() /length(unique(gene_symbols))),
          score_raw_norm_other = sum(coef.Age * coef.Ex)/ (n()),
           score_raw_norm_Best = sum(coef.Age * coef.Ex)/length(unique(gene_symbols)),
        .groups       = "drop"
      )
    prop_opposite
 prop_opposite$path= pathway_info$NAME
    res_list[[ paste0(cell,":",pathway_info$NAME) ]] <<- prop_opposite

 df_sub$Coefficient = ifelse(df_sub$Coefficient ==  'coef.Age', "Age FC", "Ex FC")
      ggplot(df_sub, aes(x = Coefficient, y = primerid.Ex, fill = FC)) +
        geom_tile(color = "white") +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits   = c(-2.5, 2.5),        
                             na.value = "grey80", name = "Fold Change", ) +
        theme_minimal() +theme_mytheme +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text.y = element_blank(),
              strip.background = element_blank()) +
        labs(x = NULL, y = NULL, title = paste(cell, "\nscore =", signif(-prop_opposite$score_raw_norm, digits = 2)))
    })

    library(cowplot)

    # Remove legends from each plot
    plots_noleg <- lapply(plots, function(p) 
      p + theme(legend.position = "none")
    )

    # Extract a shared legend from the first plot
    shared_legend <- get_legend(
      plots[[1]] +
        theme(
          legend.position = "right",
          legend.direction = "horizontal"
        )
    )

    # Arrange plots without legends into a grid
    plots_grid <- plot_grid(plotlist = plots_noleg, nrow =1 , align = "hv")

    # Combine the grid and the shared legend on the right
    final_plot <- plot_grid(
      plots_grid,
      rel_widths = c(1, 0.2), labels = c(paste("A",pathway_info$NAME))
    )

    #final_plot
 #pltlist[[ paste0(pathway_info$NAME) ]] = final_plot
    #ggsave(paste0(publication_path("HEART_SNRNA_ROOT", "analysis/2025-04-23.GREAT_rev_enrich/"),pathway_info$NAME,"_Kegg.pdf"), width = 210, height = 70, units = 'mm')

}

library(dplyr)
library(tidyr)
library(ggplot2)
All =  do.call(rbind,res_list)
All$Path_Short = stringr::str_remove(All$path , " .*")
#res_list
    #All$cell.Age

score_df_scaled <- All %>%
  group_by(cell.Age) %>%
  mutate(
    score_scaled = (-score_raw_norm - min(-score_raw_norm, na.rm = TRUE)) /
                   (max(-score_raw_norm, na.rm = TRUE) - min(-score_raw_norm, na.rm = TRUE))
  ) %>%ungroup()
# 1. Build a matrix of scaled scores: rows = pathways, cols = cell types
mat_df <- score_df_scaled %>%
  dplyr::select(cell.Age, Path_Short, score_scaled) %>%
  pivot_wider(names_from = cell.Age, values_from = score_scaled, values_fill = 0)

# convert to matrix
mat <- as.matrix(mat_df[,-1])
rownames(mat) <- mat_df$Path_Short

# 2. Hierarchical clustering on pathways (rows) and cell types (columns)
row_dist <- dist(mat)
row_hclust <- hclust(row_dist)
row_order <- row_hclust$labels[row_hclust$order]

col_dist <- dist(t(mat))
col_hclust <- hclust(col_dist)
col_order <- col_hclust$labels[col_hclust$order]

# 3. Reorder factors in the original data frame
score_df_scaled <- score_df_scaled %>%
  mutate(
    Path_Short = factor(Path_Short, levels = row_order),
    cell.Age   = factor(cell.Age, levels = col_order)
  )

# 4. Re-plot heatmap with clustered order
Scaled_summary_clustered <- ggplot(score_df_scaled, aes(x = cell.Age, y = Path_Short, fill = score_scaled)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low      = "white",
    high     = "darkred",
    midpoint = 0,
    name     = "Cell Scaled\nReversal\nScore"
  ) +
  theme_minimal() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    axis.text.y     = element_text(size = 10),
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  labs(
    x     = "Cell Type",
    y     = "Pathway",
    title = "Clustered Pathway Reversal Scores Across Cell Types"
  )
write.csv(score_df_scaled,publication_path("PAPER_OUTPUT_ROOT", "Paper_Tables/KEGG_Attenuation.csv"))

Scaled_summary_clustered
ggsave(publication_path("PAPER_OUTPUT_ROOT", "Paper_Tables/KEGG_Attenuation/KEGG.pdf"), width = 100, height = 150, units = 'mm')

library(dplyr)
library(tidyr)
library(ggplot2)
All =  do.call(rbind,res_list)
All$Path_Short = stringr::str_remove(All$path , " .*")
#res_list
    #All$cell.Age

score_df_scaled <- All %>%
  group_by(cell.Age) %>%
  mutate(
    score_scaled = (-score_raw_norm_Best - min(-score_raw_norm_Best, na.rm = TRUE)) /
                   (max(-score_raw_norm_Best, na.rm = TRUE) - min(-score_raw_norm_Best, na.rm = TRUE))
  ) %>%
  ungroup()
# 1. Build a matrix of scaled scores: rows = pathways, cols = cell types
mat_df <- score_df_scaled %>%
  dplyr::select(cell.Age, Path_Short, score_scaled) %>%
  pivot_wider(names_from = cell.Age, values_from = score_scaled, values_fill = 0)

# convert to matrix
mat <- as.matrix(mat_df[,-1])
rownames(mat) <- mat_df$Path_Short

# 2. Hierarchical clustering on pathways (rows) and cell types (columns)
row_dist <- dist(mat)
row_hclust <- hclust(row_dist)
row_order <- row_hclust$labels[row_hclust$order]

col_dist <- dist(t(mat))
col_hclust <- hclust(col_dist)
col_order <- col_hclust$labels[col_hclust$order]

# 3. Reorder factors in the original data frame
score_df_scaled <- score_df_scaled %>%
  mutate(
    Path_Short = factor(Path_Short, levels = row_order),
    cell.Age   = factor(cell.Age, levels = col_order)
  )

# 4. Re-plot heatmap with clustered order
Scaled_summary_clustered <- ggplot(score_df_scaled, aes(x = cell.Age, y = Path_Short, fill = score_scaled)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low      = "white",
    high     = "darkred",
    midpoint = 0,
    name     = "Cell Scaled\nReversal\nScore"
  ) +
  theme_minimal() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    axis.text.y     = element_text(size = 10),
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  labs(
    x     = "Cell Type",
    y     = "Pathway",
    title = "Clustered Pathway Reversal Scores Across Cell Types"
  )
write.csv(score_df_scaled,publication_path("PAPER_OUTPUT_ROOT", "Paper_Tables/KEGG_Attenuation.csv"))

Scaled_summary_clustered
ggsave(publication_path("PAPER_OUTPUT_ROOT", "Paper_Tables/KEGG_Attenuation/KEGG2.pdf"), width = 100, height = 150, units = 'mm')
