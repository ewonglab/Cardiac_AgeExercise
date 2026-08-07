# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Compare cardiac macrophage signatures with the supplied human reference dataset.
# Steps: 1) Load reference and study signatures; 2) map homologous markers; 3) export cross-species comparison results.
# =======================================
# Source reformatted from snRNA_AgeEx/ShuEtAl(1).ipynb.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

# Shu_Et Al

# Vaibhao is interested in getting seeing which subpop of fibroblast is enriched for  wif1  - we can also compare the subpops in our fibroblast and also macs with this data

library(DropletUtils)
library(Seurat)
library(ggplot2)
library(scran)
library(dplyr)
library(knitr)
library(stringr)



setwd(publication_path("HEART_SNRNA_ROOT", "data/2024-10-16.Shu_etal"))
run.merge_Filt = readRDS("Shu_Filtered.rds")




Shu_Etal = run.merge_Filt
Shu_Etal$CellType = Shu_Etal$cell_type

HeartMacro = readRDS(publication_path("HEART_SNRNA_ROOT", "data/2022-10-17.SharedObject/2023-10-12_Macrophages_RNA_TooManyAnnotated.rds"))

Shu_Etal$CellType = Shu_Etal$cell_type
DimPlot(Shu_Etal, group.by = "CellType")

reticulate::py_require("matplotlib")
reticulate::py_require("scvelo")
library(seuratTools)

convert_mouse_seu_to_human <- function(seu) {
    # transfer default species expression data to a species-specific assay
    seu[["mouse"]] <- seu[["RNA"]]

    new_rownames <- convert_symbols_by_species(src_genes = rownames(seu), src_species = "mouse")

    seu_slots <- c("counts", "data", "scale.data", "meta.features")

    for (i in seu_slots) {
        current_slot <- slot(seu@assays[["RNA"]], i)
        if (!(dim(current_slot)[1] == 0)) {
            rownames(slot(seu@assays[["RNA"]], i)) <- new_rownames
        }
    }

    return(seu)
}
HeartMacro = convert_mouse_seu_to_human(HeartMacro)

HeartMacro





setwd(publication_path("PAPER_OUTPUT_ROOT", "Paper_Tables/Shu_EtAl_QC/"))

# Subset To Macrophage 

Macro_Shu = subset(Shu_Etal, CellType == "Myeloid")

Macro_Shu <- FindVariableFeatures(Macro_Shu)
Macro_Shu <- ScaleData(Macro_Shu, verbose = FALSE)
Macro_Shu <- RunPCA(Macro_Shu, npcs = 30, verbose = FALSE)
Macro_Shu <- FindNeighbors(Macro_Shu, dims = 1:10) #note dimensions 
Macro_Shu <- FindClusters(Macro_Shu, resolution = 0.5)
Macro_Shu <- RunUMAP(Macro_Shu, dims = 1:10)

Macro_Shu$ShuCluster = Macro_Shu$seurat_clusters

DimPlot(Macro_Shu, reduction = "umap", group.by = "run")
 DimPlot(Macro_Shu, reduction = "umap", group.by = "CellType")
 DimPlot(Macro_Shu, reduction = "umap", group.by = "seurat_clusters", label = T)

genes_to_plot <- rev(c(
  "Mmp2", "Col1a1", "Gli2", "Scara5", "Col1a2", "Pdgfra", "Col3a1",
  "Col5a2", "Cthrc1", "Acta2", "Tead3", "Tcf21", "Postn", "Fn1",
  "Cx3cr1", "Ccr2", "Timd4", "Lyve1", "Itgam", "Cd163", "Adgre1", "HLA-A"
))
DotPlot(Macro_Shu, features = toupper(genes_to_plot), group.by ="ShuCluster")+coord_flip()+
  scale_color_gradientn(colors = c("blue", "yellow", "red")) 

genes_to_plot <- rev(c(
 "CCL3L1","CCL4A2","Il1B","RGS1","PLTP","LGMN","CD163","DCN","EREG","TOP2A"
))
DotPlot(Macro_Shu, features = toupper(genes_to_plot), group.by ="ShuCluster")+coord_flip()+
  scale_color_gradientn(colors = c("blue", "yellow", "red")) 

# Labels 
#M1 = 1,2,3
#M2 =0,4,5,6 
#Remodeling = 7
#Monocyte = 8
#Dividing = 9
# Subset To Macrophage 

# Macro_Shu = subset(Shu_Etal, CellType == "Myeloid")

# Macro_Shu <- FindVariableFeatures(Macro_Shu)
# Macro_Shu <- ScaleData(Macro_Shu, verbose = FALSE)
# Macro_Shu <- RunPCA(Macro_Shu, npcs = 30, verbose = FALSE)
# Macro_Shu <- FindNeighbors(Macro_Shu, dims = 1:10) #note dimensions 
# Macro_Shu <- FindClusters(Macro_Shu, resolution = 0.5)
# Macro_Shu <- RunUMAP(Macro_Shu, dims = 1:10)

Macro_Shu$Label = "NA"
Macro_Shu$Label = ifelse(Macro_Shu$ShuCluster %in% c(1,2,3),"M1",Macro_Shu$Label)
Macro_Shu$Label = ifelse(Macro_Shu$ShuCluster %in% c(0,4,5,6),"M2",Macro_Shu$Label)
Macro_Shu$Label = ifelse(Macro_Shu$ShuCluster %in% c(7),"Remodeling",Macro_Shu$Label)
Macro_Shu$Label = ifelse(Macro_Shu$ShuCluster %in% c(8),"Monocyte",Macro_Shu$Label)
Macro_Shu$Label = ifelse(Macro_Shu$ShuCluster %in% c(9),"Dividing",Macro_Shu$Label)
# Macro_Shu = subset(Macro_Shu, Label != "NA")
# Idents(Macro_Shu) = 'Label'
# MAcrophageMarks = FindAllMarkers(Macro_Shu)

 DimPlot(Macro_Shu, reduction = "umap", group.by = "Label", label = T)

# Map Tentitive MacFib back onto Glbal
Idents(Macro_Shu) = 'Label'

MAcrophageMarks <- FindAllMarkers(
  object = Macro_Shu,
  only.pos = TRUE,
  logfc.threshold = 0.1,
  test.use = "wilcox",        # Standard Wilcoxon
  base.helper = NULL,         # Disables the v5 group-level pseudocount logic
  mean.fxn = rowMeans,        # Forces cell-level averaging
  fc.name = "avg_log2FC"      # Ensures Base 2 is used
)

genes= MAcrophageMarks[ MAcrophageMarks$avg_log2FC > 0.5, 'gene']
length(genes)
DoHeatmap(Macro_Shu,features = genes)

getwd()

Macro_Shu$Macfib = ifelse(Macro_Shu$seurat_clusters == 7,"MacFib","Macrophage")

DimPlot(Macro_Shu, group.by = "Macfib")+NoLegend()
ggsave("MacFib_Umap_Shu.pdf")



DimPlot(Shu_Etal, reduction = "umap", group.by = "CellType", label = T)

#Mac Fib Markers 
genes_to_plot <- rev(c(
  "Mmp2", "Col1a1", "Gli2", "Scara5", "Col1a2", "Pdgfra", "Col3a1",
  "Col5a2", "Cthrc1", "Acta2", "Tead3", "Tcf21", "Postn", "Fn1",
  "Cx3cr1", "Ccr2", "Timd4", "Lyve1", "Itgam", "Cd163", "Adgre1", "HLA-A"
))
colours = c("#4575B4","#74ADD1","#ABD9E9","#E0F3F8","#FFFFBF", 
               "#FEE090","#FDAE61","#F46D43","#D73027")[c(1,1:9,9)]

FibMac = subset(Shu_Etal,CellType %in% c("M1","M2","TentitiveMacFib","VIC"))
FibMac$CellType = factor(FibMac$CellType, levels = c("VIC","TentitiveMacFib","M1","M2"))
DotPlot(FibMac, features = toupper(genes_to_plot), group.by ="CellType")+coord_flip()+
  scale_color_gradientn(colors = colours)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

DotPlot(FibMac, features = toupper(genes_to_plot), group.by ="CellType", scale = FALSE)+coord_flip()+
  scale_color_gradientn(colors = colours)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
ggsave("MacFibExression.Shu.pdf")

# MErge by itegration


Macro.integrated = readRDS(publication_path("HEART_SNRNA_ROOT", "data/2024-10-16.Shu_etal/Integrated_Macro.rds"))



DimPlot(Macro.integrated, reduction = "umap", group.by = "run", label = TRUE)
 DimPlot(Macro.integrated, reduction = "umap", group.by = "CellType")
 DimPlot(Macro.integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
 DimPlot(Macro.integrated, reduction = "umap", group.by = "ShuCluster", label = TRUE)
 DimPlot(Macro.integrated, reduction = "umap", group.by = "Label", label = TRUE)

FeaturePlot(Macro.integrated, features = c("IL1B","CCL3L1","CCL4L2","TNF"))

DotPlot(Macro.integrated, features = c("IL1B","CCL3L1","CCL4L2","TNF"), group.by = "Mouse_Human_Fact", scale = TRUE)
DotPlot(Macro.integrated, features = c("IL1B","CCL3L1","CCL4L2","TNF"), group.by = "Mouse_Human_Fact", scale = FALSE)

table(Macro.integrated$Mouse_Human)

Macro.integrated$Mouse_Human_Fact = factor(Macro.integrated$Mouse_Human, levels = c("M2_Mouse","M2_Human","M1_Mouse","M1_Human","M_Ex_Mouse","Remodeling_Human",'Monocyte_Human','Dividing_Human'))

library(ggplot2)
library(viridis)
Idents(Macro.integrated) = "Mouse_Human_Fact"
plot = DoHeatmap(Macro.integrated,features = genes)
#ggsave("IntegratedMacrophageMarkers.pdf", width = 20 , height = 20)







MAcrophageMarks %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5) %>%
    slice_head(n = 100) %>%
    ungroup() -> top10
top10 = top10[ top10$cluster %in% c("M1","M2","Remodeling"),]

test = subset(Macro.integrated,Mouse_Human_Fact %in% c('M1_Human','M1_Mouse','M2_Human','M2_Mouse','M_Ex_Mouse','Remodeling_Human'))
test = ScaleData(test)


write.csv(MAcrophageMarks, "Shut_etAl_MAcrophageMarks.csv")

getwd()

HeartMacro@misc

Keep = test@assays$integrated@scale.data

Keep = Keep[rownames(Keep) %in% top10$gene,]

 test@assays$integrated@scale.data = Keep

Idents(test) = 'Mouse_Human_Fact'



 test@assays$integrated@scale.data = scale(Keep)

DoHeatmap(test,features = top10$gene, raster = TRUE)

plot = DoHeatmap(test,features = top10$gene, raster = TRUE)
ggsave(plot =  plot, filename= "IntegratedMacrophageMarkers_Top100_Raster_Cheat.pdf", width = 20 , height = 20)

plot + scale_fill_viridis(option =  "magma")



table(Macro_Shu$Label)

# ============================================================
# Supplementary plotting for Shu et al. + integrated macrophage object
# Generates:
# 1. Shu et al. labelled CellType UMAP
# 2. Shu et al. collagen gene FeaturePlot UMAP panel
# 3. Shu et al. collagen DotPlot in VIC / TentitiveMacFib / M1 / M2
# 4. Integrated object labelled macrophage-state UMAP
# 5. Integrated object identity-marker FeaturePlot UMAP panel
# 6. Integrated object identity-marker DotPlot
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# ------------------------------------------------------------
# Output directory
# ------------------------------------------------------------
outdir <- publication_path("PAPER_OUTPUT_ROOT", "Paper_Tables/Shu_EtAl_QC/Supplementary_Figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------
present_features <- function(object, features, assay = NULL, label = "features") {
  if (!is.null(assay) && assay %in% names(object@assays)) {
    feature_universe <- rownames(object[[assay]])
  } else {
    feature_universe <- rownames(object)
  }

  keep <- features[features %in% feature_universe]
  missing <- setdiff(features, keep)

  message("\n", label, " retained: ", paste(keep, collapse = ", "))
  if (length(missing) > 0) {
    message(label, " missing and skipped: ", paste(missing, collapse = ", "))
  }

  return(keep)
}

feature_panel_height <- function(n_features, ncol = 4, row_height = 2.35, min_height = 4.5) {
  max(min_height, ceiling(n_features / ncol) * row_height)
}

plot_feature_umap_panel <- function(
    object,
    features,
    filename,
    title,
    cols,
    ncol = 4,
    reduction = "umap",
    pt.size = 0.08,
    width = 8.27
) {
  if (length(features) == 0) {
    warning("No features available for: ", title)
    return(NULL)
  }

  plot_list <- FeaturePlot(
    object = object,
    features = features,
    reduction = reduction,
    order = TRUE,
    min.cutoff = "q05",
    max.cutoff = "q95",
    cols = cols,
    pt.size = pt.size,
    combine = FALSE
  )

  plot_list <- lapply(plot_list, function(p) {
    p +
      theme_classic(base_size = 10) +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(
          size = 10,
          face = "bold",
          hjust = 0.5
        ),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)
      )
  })

  final_plot <- wrap_plots(plot_list, ncol = ncol) +
    plot_annotation(
      title = title,
      theme = theme(
        plot.title = element_text(
          size = 13,
          face = "bold",
          hjust = 0.5
        )
      )
    )

  ggsave(
    filename = file.path(outdir, filename),
    plot = final_plot,
    width = width,
    height = feature_panel_height(length(features), ncol = ncol),
    units = "in"
  )

  return(final_plot)
}

# ------------------------------------------------------------
# Use RNA assay for expression visualisation where available
# ------------------------------------------------------------
if ("RNA" %in% names(Macro_Shu@assays)) {
  DefaultAssay(Macro_Shu) <- "RNA"
}

if ("RNA" %in% names(Macro.integrated@assays)) {
  DefaultAssay(Macro.integrated) <- "RNA"
}

# ------------------------------------------------------------
# Palette-derived expression colours
# ------------------------------------------------------------
collagen_expr_cols <- c("grey92", "#42977E")  # Fibroblast-like green
marker_expr_cols   <- c("grey92", "#A53D99")  # Macrophage purple

# ============================================================
# 1. Shu et al. labelled CellType UMAP
# ============================================================
p_shu_celltype_umap <- DimPlot(
  Macro_Shu,
  reduction = "umap",
  group.by = "Label",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.08
) +
  theme_classic(base_size = 11) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
  ) +
  ggtitle("Shu et al. reference atlas cell-type annotation")

ggsave(
  filename = file.path(outdir, "Supp_ShU_CellType_UMAP.pdf"),
  plot = p_shu_celltype_umap,
  width = 8.27,
  height = 7.2,
  units = "in"
)

# ============================================================
# 2. Shu et al. collagen gene expression on the full UMAP
# ============================================================
collagen_genes_shu <- c(
  "COL1A1",
  "COL1A2",
  "COL3A1",
  "COL4A1",
  "COL4A2",
  "COL5A1",
  "COL5A2",
  "COL6A1",
  "COL6A2"
)

collagen_genes_shu_present <- present_features(
  object = Macro_Shu,
  features = collagen_genes_shu,
  assay = "RNA",
  label = "Shu collagen genes"
)

p_shu_collagen_umaps <- plot_feature_umap_panel(
  object = Macro_Shu,
  features = collagen_genes_shu_present,
  filename = "Supp_ShU_Collagen_FeatureUMAPs.pdf",
  title = "Collagen gene expression across the Shu et al. atlas",
  cols = collagen_expr_cols,
  ncol = 3,
  pt.size = 0.08
)

# ============================================================
# 3. Shu collagen DotPlot in the Fibroblast / MacFib / M1 / M2 comparison
# ============================================================
#shu_dot_celltypes <- c("VIC", "TentitiveMacFib", "M1", "M2")

Shu_Collagen_Dot <- subset(
  Macro_Shu,
)

Shu_Collagen_Dot$Label <- factor(
  Shu_Collagen_Dot$Label,
)

p_shu_collagen_dot <- DotPlot(
  Shu_Collagen_Dot,
  features = collagen_genes_shu_present,
  group.by = "Label",
  scale = FALSE
) +
  coord_flip() +
  scale_color_gradient(
    low = "grey92",
    high = "#42977E"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "italic"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "right"
  ) +
  ggtitle("Collagen gene expression in Shu et al. fibroblast-like and myeloid populations")

ggsave(
  filename = file.path(outdir, "Supp_ShU_Collagen_DotPlot_FibMac.pdf"),
  plot = p_shu_collagen_dot,
  width = 8.27,
  height = 5.8,
  units = "in"
)

# ============================================================
# 4. Define / order integrated macrophage identity variable
# ============================================================
if ("Mouse_Human" %in% colnames(Macro.integrated@meta.data)) {
  Macro.integrated$Mouse_Human_Fact <- factor(
    Macro.integrated$Mouse_Human,
    levels = c(
      "M2_Mouse",
      "M2_Human",
      "M1_Mouse",
      "M1_Human",
      "M_Ex_Mouse",
      "Remodeling_Human",
      "Monocyte_Human",
      "Dividing_Human"
    )
  )
}

integrated_group_var <- if ("Mouse_Human_Fact" %in% colnames(Macro.integrated@meta.data)) {
  "Mouse_Human_Fact"
} else if ("Mouse_Human" %in% colnames(Macro.integrated@meta.data)) {
  "Mouse_Human"
} else if ("Label" %in% colnames(Macro.integrated@meta.data)) {
  "Label"
} else {
  stop("Could not find Mouse_Human_Fact, Mouse_Human, or Label in Macro.integrated metadata.")
}

# ============================================================
# 5. Integrated object labelled UMAP
# ============================================================
p_integrated_umap <- DimPlot(
  Macro.integrated,
  reduction = "umap",
  group.by = integrated_group_var,
  label = TRUE,
  repel = TRUE,
  pt.size = 0.10
) +
  theme_classic(base_size = 11) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
  ) +
  ggtitle("Integrated mouse–human macrophage state annotation")

ggsave(
  filename = file.path(outdir, "Supp_Integrated_MacrophageState_UMAP.pdf"),
  plot = p_integrated_umap,
  width = 8.27,
  height = 7.2,
  units = "in"
)

# ============================================================
# 6. Integrated object identity-marker expression UMAPs
#    Uses the same annotation markers already present in the notebook
# ============================================================
identity_markers_integrated <- c(
  "CCL3L1",
  "CCL4L2",
  "IL1B",
  "RGS1",
  "PLTP",
  "LGMN",
  "CD163",
  "DCN",
  "EREG",
  "TOP2A"
)

identity_markers_integrated_present <- present_features(
  object = Macro.integrated,
  features = identity_markers_integrated,
  assay = "RNA",
  label = "Integrated identity markers"
)

p_integrated_marker_umaps <- plot_feature_umap_panel(
  object = Macro.integrated,
  features = identity_markers_integrated_present,
  filename = "Supp_Integrated_IdentityMarker_FeatureUMAPs.pdf",
  title = "Identity-marker expression in the integrated macrophage object",
  cols = marker_expr_cols,
  ncol = 3,
  pt.size = 0.10
)

# ============================================================
# 7. Integrated object identity-marker DotPlot
# ============================================================
p_integrated_marker_dot <- DotPlot(
  Macro.integrated,
  features = identity_markers_integrated_present,
  group.by = integrated_group_var,
  scale = FALSE
) +
  coord_flip() +
  scale_color_gradient(
    low = "grey92",
    high = "#A53D99"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "italic"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "right"
  ) +
  ggtitle("Identity-marker expression across integrated macrophage states")

ggsave(
  filename = file.path(outdir, "Supp_Integrated_IdentityMarker_DotPlot.pdf"),
  plot = p_integrated_marker_dot,
  width = 8.27,
  height = 6.2,
  units = "in"
)

# ============================================================
# 8. Print plots to notebook
# ============================================================
p_shu_celltype_umap
p_shu_collagen_umaps
p_shu_collagen_dot
p_integrated_umap
p_integrated_marker_umaps
p_integrated_marker_dot
