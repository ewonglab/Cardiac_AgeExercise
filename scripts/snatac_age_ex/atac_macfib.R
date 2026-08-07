# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Compare chromatin accessibility between MacFib and conventional macrophage states.
# Steps: 1) Load the annotated ATAC object; 2) define comparison groups; 3) run marker analysis and export results.
# =======================================
# Source reformatted from snATAC_AgeEx/ATAC_MacFib.ipynb.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

# ATAC Mac Fib comparison 

library(Signac)

library(Seurat)



ATAC = readRDS(publication_path("HEART_SNRNA_ROOT", "data/2022-10-17.SharedObject/2024-07-16_Full_ATAC_MatchAnnotated.rds"))



ATAC$comparison_group <- ifelse(ATAC$CellType %in% c("M1", "M2"), "M1M2", ATAC$CellType)

# Step 2: Verify the new grouping
table(ATAC$comparison_group)
DefaultAssay(ATAC) = "ATAC"

# Step 3: Find markers between MacFib and M1M2
markers <- FindMarkers(
  object = ATAC,
  ident.1 = "MacFib",  # Compare MacFib
  ident.2 = "M1M2",    # Against pooled M1 and M2
  group.by = "comparison_group",
     test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_ATAC'# Use the new grouping
)




setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-10-5.M1M2_ATAC"))

saveRDS(markers,"markers.MacFib_M1M2.rds")

setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2023-10-5.M1M2_ATAC"))

markers = readRDS("markers.MacFib_M1M2.rds")

head(markers)

library(GenomicRanges)
top.da.peak <- rownames(markers[markers$p_val < 0.005 & markers$pct.1 > 0.2 & markers$avg_log2FC > 1, ])
length(top.da.peak)

ranges.show_Enrich <- StringToGRanges(top.da.peak)
ranges.show_BG <- StringToGRanges(rownames(markers))

library(ggplot2)
theme_mytheme <- theme(
  # Specify the overall plot appearance
  plot.background = element_blank(),
  plot.title = element_text(size = 16, face = "bold"),
  plot.subtitle = element_text(size = 14),
  
  # Specify the axis appearance
  axis.text = element_text(size = 12, color = "black"),
  axis.title = element_text(size = 14, color = "black"),
  axis.line = element_line(color = "black", size = 0.5),
  axis.ticks = element_line(color = "black", size = 0.5),
  
  # Specify the legend appearance
  legend.title = element_blank(),
  legend.text = element_text(size = 12),
  legend.position = "right",
  
  # Specify the panel/grid appearance, no grid line 
  panel.background = element_blank(),#element_rect(fill = "white", color = "black"),
  panel.grid.major = element_blank(),#element_line(color = "lightgray", size = 0.2),
  panel.grid.minor = element_blank(),#element_line(color = "lightgray", size = 0.2),
  
  # Specify the plot margins
  plot.margin = margin(1, 1, 1, 1, "cm"),
  
  #set the ratio of the plot 
  #aspect.ratio=1
)

library(rGREAT)


great_results <- great(    ranges.show_Enrich,    background  = ranges.show_BG,  tss_source =   "mm10", gene_sets = c("GO:BP") # Use mouse genome
)

great_results2 <- great(    ranges.show_Enrich,    background  = ranges.show_BG,  tss_source =   "mm10", gene_sets = c("GO:MF") # Use mouse genome
)
great_results3 <- great(    ranges.show_Enrich,    background  = ranges.show_BG,  tss_source =   "mm10", gene_sets = c("GO:CC") # Use mouse genome
)

tb = rbind(getEnrichmentTable(great_results),getEnrichmentTable(great_results2),getEnrichmentTable(great_results3))



DT::datatable(tb)

top.da.peak <- rownames(markers[markers$p_val < 0.005 & markers$pct.2 > 0.2 & markers$avg_log2FC < -1, ])
length(top.da.peak)
ranges.show_Enrich <- StringToGRanges(top.da.peak)
ranges.show_BG <- StringToGRanges(rownames(markers))

great_results_M1M2_1 <- great(    ranges.show_Enrich,    background  = ranges.show_BG,  tss_source =   "mm10", gene_sets = c("GO:BP") # Use mouse genome
)

great_results_M1M2_2 <- great(    ranges.show_Enrich,    background  = ranges.show_BG,  tss_source =   "mm10", gene_sets = c("GO:MF") # Use mouse genome
)
great_results_M1M2_3 <- great(    ranges.show_Enrich,    background  = ranges.show_BG,  tss_source =   "mm10", gene_sets = c("GO:CC") # Use mouse genome
)

tb2 = rbind(getEnrichmentTable(great_results),getEnrichmentTable(great_results2),getEnrichmentTable(great_results3))

getwd()

write.csv(tb2 ,"ATAC_MACFIB_Down.csv")

write.csv(tb ,"ATAC_MACFIB_Up.csv")

getwd()
library(ggplot2)
library(stringr)
tb2$description <- str_wrap(tb2$description, width = 40)

ggplot(tb2[1:5,], aes(x = fold_enrichment, y = description, size = observed_gene_hits, color = (p_adjust))) +
  geom_point() +
scale_size_continuous(limits = c(0, 500))+
   scale_color_gradient2(low = "red", high = "yellow",mid = 'red',    limits = c(0, 1) ) +  # Use plasma color scale for p-values
  labs(
    title = "GREAT Enrichment Results",
    x = "Fold Enrichment",
    y = "Biological Process",
    size = "Gene Ratio",
    color = "-log10(P-value)"
  ) +
  theme_mytheme + xlim(0,3)+
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14)
  )
ggsave("Macrophage_GreatTerms.pdf")

tb2[1:5,]

library(ggplot2)
library(stringr)
tb$description <- str_wrap(tb$description, width = 40)
ggplot(tb[1:5,], aes(x = fold_enrichment, y = description, size = observed_gene_hits, color = (p_adjust))) +
  geom_point() +
scale_size_continuous(limits = c(0, 500))+
  scale_color_gradient2(low = "red", high = "yellow",mid = 'red',    limits = c(0, 1) ) +  # Use plasma color scale for p-values
  labs(
    title = "GREAT Enrichment Results",
    x = "Fold Enrichment",
    y = "Biological Process",
    size = "Gene Ratio",
    color = "-log10(P-value)"
  ) +
    theme_mytheme + xlim(0,3)+
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14)
  )
ggsave("MacFib_GreatTerms.pdf")

library(ggplot2)
library(stringr)
tb$description <- str_wrap(tb$description, width = 40)
ggplot(tb[1:5,], aes(x = fold_enrichment, y = description, size = observed_gene_hits, color = (p_adjust))) +
  geom_point() +
scale_size_continuous(limits = c(0, 500))+
  scale_color_gradient2(low = "red", high = "yellow",mid = 'red',    limits = c(0, 1) ) +  # Use plasma color scale for p-values
  labs(
    title = "GREAT Enrichment Results",
    x = "Fold Enrichment",
    y = "Biological Process",
    size = "Gene Ratio",
    color = "-log10(P-value)"
  ) +
    theme_mytheme + xlim(0,3)+
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14)
  )
ggsave("MacFib_GreatTerms.pdf")

#DT::datatable(tb, filter = 'top')
Keep= c('immune response',
'vesicle-mediated transport',
'regulation of immune response',
'innate immune response',
'myeloid leukocyte differentiation')
plt = tb2[ tb2$description %in% Keep,]
plt$description = factor(plt$description, levels = rev(Keep))
library(ggplot2)
library(stringr)
tb$description <- str_wrap(tb$description, width = 40)
ggplot(plt, aes(x = fold_enrichment, y = description, size = observed_gene_hits, color = (p_adjust))) +
  geom_point() +
scale_size_continuous(limits = c(0, 500))+
  scale_color_gradient2(low = "red", high = "yellow",mid = 'red',    limits = c(0, 1) ) +  # Use plasma color scale for p-values
  labs(
    title = "GREAT Enrichment Results",
    x = "Fold Enrichment",
    y = "Biological Process",
    size = "Gene Ratio",
    color = "-log10(P-value)"
  ) +
    theme_mytheme + xlim(0,3)+
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14)
  )
ggsave("Macrophage_GreatTerms_filt.pdf")

DT::datatable(tb,filter = 'top')

#DT::datatable(tb, filter = 'top')
Keep= c('extracellular matrix structural constituent',
'external encapsulating structure',
'collagen-containing extracellular matrix',
'structural molecule activity',
'animal organ morphogenesis')
plt = tb[ tb$description %in% Keep,]
plt$description = factor(plt$description, levels = rev(Keep))
library(ggplot2)
library(stringr)
tb$description <- str_wrap(tb$description, width = 40)
ggplot(plt, aes(x = fold_enrichment, y = description, size = observed_gene_hits, color = (p_adjust))) +
  geom_point() +
scale_size_continuous(limits = c(0, 500))+
  scale_color_gradient2(low = "red", high = "yellow",mid = 'red',    limits = c(0, 1) ) +  # Use plasma color scale for p-values
  labs(
    title = "GREAT Enrichment Results",
    x = "Fold Enrichment",
    y = "Biological Process",
    size = "Gene Ratio",
    color = "-log10(P-value)"
  ) +
    theme_mytheme + xlim(0,3)+
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14)
  )

#DT::datatable(tb, filter = 'top')
Keep= c('extracellular matrix structural constituent',
'external encapsulating structure',
'collagen-containing extracellular matrix',
'structural molecule activity',
'animal organ morphogenesis')
plt = tb[ tb$description %in% Keep,]
plt$description = factor(plt$description, levels = rev(Keep))
library(ggplot2)
library(stringr)
tb$description <- str_wrap(tb$description, width = 40)
ggplot(plt, aes(x = fold_enrichment, y = description, size = observed_gene_hits, color = (p_adjust))) +
  geom_point() +
scale_size_continuous(limits = c(0, 500))+
  scale_color_gradient2(low = "red", high = "yellow",mid = 'red',    limits = c(0, 1) ) +  # Use plasma color scale for p-values
  labs(
    title = "GREAT Enrichment Results",
    x = "Fold Enrichment",
    y = "Biological Process",
    size = "Gene Ratio",
    color = "-log10(P-value)"
  ) +
    theme_mytheme + xlim(0,3)+
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14)
  )
ggsave("M_Ex_GreatTerms_filt.pdf")

DAP_GREAT = read.csv(publication_path("HEART_SNRNA_ROOT", "analysis/2023-10-11.AllMacs/GREAT_Back/All.GREAT.csv"))

table(DAP_GREAT$p_adjust_hyper < 0.1)
 DAP_GREAT$Direction  = factor(DAP_GREAT$Direction)
 DAP_GREAT$description  = factor(DAP_GREAT$description)
DT::datatable( DAP_GREAT[DAP_GREAT$p_adjust_hyper < 0.1,], filter = 'top')

DAP_GREAT$state =paste(DAP_GREAT$Direction, DAP_GREAT$cell)
Sig  = DAP_GREAT[ DAP_GREAT$p_adjust_hyper < 0.1,]
length(unique(Sig$cell))
table(Sig$Description == "MAPK cascade", Sig$Direction)
table(Sig$Description == "extracellular matrix assembly", Sig$Direction)
table(Sig$Description == "MAPK cascade")

#sort(table(Sig$Description), decreasing = TRUE)

table(Sig$Description == "muscle structure development", Sig$state)

DT::datatable(Sig[Sig$Description == "positive regulation of p38MAPK cascade",], filter = 'top')
