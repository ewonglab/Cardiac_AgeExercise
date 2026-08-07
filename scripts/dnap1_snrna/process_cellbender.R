# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Process CellBender-corrected AAV dnAP-1 snRNA-seq matrices.
# Steps: 1) Load corrected matrices; 2) perform QC and integration; 3) assign cell types and save the analysis object.
# =======================================
# Source reformatted from dnAP-1_snRNA/Process-Cellbender.ipynb.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))





library(Seurat)
library(tidyverse)

base_dir <- publication_path("DNAP1_SNRNA_ROOT", "data/CellBender_OCM")
pools <- c("Pool_MS", "Pool_SW", "Pool_VW")

# Prefer the filtered CellBender output for Seurat.
# Handles both naming styles:
#   old naming: Sample_2_cellbender_filtered_filtered.h5
#   cleaner naming: Sample_2_cellbender_filtered.h5
read_cellbender_10x_h5 <- function(filename, feature_type_keep = "Gene Expression") {
  h5 <- H5File$new(filename, mode = "r")
  on.exit(h5$close_all())
  
  # Find the actual 10x matrix group.
  # Standard CellRanger/CellBender output should have /matrix.
  if ("matrix" %in% names(h5)) {
    mat_group <- h5[["matrix"]]
  } else {
    stop(
      "Could not find /matrix group in H5 file.\n",
      "This file may not be a CellRanger-format CellBender output.\n",
      "Run inspect_h5(filename) and check the H5 structure."
    )
  }
  
  required <- c("data", "indices", "indptr", "shape", "barcodes", "features")
  missing_required <- setdiff(required, names(mat_group))
  
  if (length(missing_required) > 0) {
    stop(
      "The /matrix group is missing required fields: ",
      paste(missing_required, collapse = ", ")
    )
  }
  
  data <- mat_group[["data"]][]
  indices <- mat_group[["indices"]][]
  indptr <- mat_group[["indptr"]][]
  shape <- mat_group[["shape"]][]
  barcodes <- mat_group[["barcodes"]][]
  
  features <- mat_group[["features"]]
  
  feature_ids <- if ("id" %in% names(features)) {
    features[["id"]][]
  } else {
    paste0("feature_", seq_len(shape[1]))
  }
  
  feature_names <- if ("name" %in% names(features)) {
    features[["name"]][]
  } else {
    feature_ids
  }
  
  feature_types <- if ("feature_type" %in% names(features)) {
    features[["feature_type"]][]
  } else {
    rep("Gene Expression", length(feature_names))
  }
  
  # Build sparse matrix.
  # 10x H5 stores CSC matrix with:
  #   rows = features
  #   columns = barcodes
  counts <- sparseMatrix(
    i = indices + 1,
    p = indptr,
    x = data,
    dims = shape
  )
  
  rownames(counts) <- make.unique(feature_names)
  colnames(counts) <- barcodes
  
  # Subset to Gene Expression if present.
  if (!is.null(feature_type_keep)) {
    keep <- feature_types == feature_type_keep
    
    if (any(keep)) {
      counts <- counts[keep, , drop = FALSE]
    } else {
      warning(
        "No features found with feature_type == '", feature_type_keep, "'. ",
        "Keeping all features instead. Available feature types: ",
        paste(unique(feature_types), collapse = ", ")
      )
    }
  }
  
  return(counts)
}
find_cellbender_h5 <- function(sample_dir, sample) {
  candidates <- c(
    file.path(sample_dir, paste0(sample, "_cellbender_filtered_filtered.h5"))
  )
  
  existing <- candidates[file.exists(candidates)]
  
  if (length(existing) == 0) {
    warning("No CellBender H5 found for: ", sample_dir)
    return(NA_character_)
  }
  
  existing[[1]]
}

process_pool <- function(pool) {
  sample_outs_dir <- file.path(base_dir, pool)
  
  if (!dir.exists(sample_outs_dir)) {
    warning("Missing pool directory: ", sample_outs_dir)
    return(list())
  }
  
  samples <- list.dirs(sample_outs_dir, full.names = FALSE, recursive = FALSE)
  
  seurat_list <- lapply(samples, function(sample) {
    sample_dir <- file.path(sample_outs_dir, sample)
    h5_file <- find_cellbender_h5(sample_dir, sample)
    
    if (is.na(h5_file)) {
      return(NULL)
    }
    
    message("Reading: ", h5_file)
    
    counts <- read_cellbender_10x_h5(h5_file)
    message("here")
    
    
    sobj <- CreateSeuratObject(
      counts = counts,
      project = sample
    )
    
    sobj$Orig_Pool <- pool
    sobj$Sample_ID <- sample
    sobj$CellBender_H5 <- h5_file
    
    
    sobj
  })
  
  seurat_list <- seurat_list[!vapply(seurat_list, is.null, logical(1))]
  
  seurat_list
}

all_seurat_objects <- lapply(pools, process_pool)
all_seurat_objects <- unlist(all_seurat_objects, recursive = FALSE)

if (length(all_seurat_objects) == 0) {
  stop("No CellBender Seurat objects were created. Check base_dir and output file names.")
}

barcode_prefixes <- vapply(all_seurat_objects, function(x) {
  paste0(unique(x$Orig_Pool), "_", unique(x$Sample_ID))
}, character(1))

if (length(all_seurat_objects) == 1) {
  master_seurat <- all_seurat_objects[[1]]
  master_seurat <- RenameCells(master_seurat, add.cell.id = barcode_prefixes[[1]])
} else {
  master_seurat <- merge(
    x = all_seurat_objects[[1]],
    y = all_seurat_objects[2:length(all_seurat_objects)],
    add.cell.ids = barcode_prefixes,
    project = "AAV_SingleCell_CellBender"
  )
}

rm(all_seurat_objects)
gc()

head(master_seurat@meta.data)
table(master_seurat$Orig_Pool)
table(master_seurat$Sample_ID)

# Optional save
saveRDS(
  master_seurat,
  file = file.path(base_dir, "AAV_SingleCell_CellBender_master_seurat.rds")
)






master_seurat[["percent.mt"]] <- PercentageFeatureSet(master_seurat, pattern = "^mt-")
table(master_seurat$nFeature_RNA < 500)
table(master_seurat$nFeature_RNA > 5000)
table(master_seurat$percent.mt > 5)

master_seurat_filtered <- subset(
  master_seurat,
  subset = nFeature_RNA > 500 & 
           nFeature_RNA < 5000 & 
           percent.mt < 5
)

library(tidyverse)

# 1. Load the raw scores generated by your Python script
doublet_data <- read_csv(publication_path("DNAP1_SNRNA_ROOT", "Code/all_samples_scrublet_scores.csv"))

# 2. Create a lookup table specifying the threshold for EVERY unique Pool + Sample combination
# Note: Replace these dummy combinations with the actual 10 from your PDFs
custom_thresholds <- tribble(
    ~Pool_ID,   ~Sample_ID, ~Manual_Threshold,
    "Pool_MS",  "Sample_2", 0.092,
    "Pool_MS",  "Sample_3", 0.062,
    "Pool_MS",  "Sample_5", 0.107,
    "Pool_MS",  "Sample_6", 0.088,
    "Pool_SW",  "Sample_1", 0.167,
    "Pool_SW",  "Sample_4", 0.195,
    "Pool_VW",  "Sample_2", 0.094,
    "Pool_VW",  "Sample_3", 0.05,
    "Pool_VW",  "Sample_5", 0.098,
    "Pool_VW",  "Sample_6", 0.075
  # Add the remaining 4 samples here...
)

# 3. Join the thresholds to your data using BOTH columns to ensure a perfect match
adjusted_doublets <- doublet_data %>%
  dplyr::left_join(custom_thresholds, by = c("Pool_ID", "Sample_ID")) %>%
  dplyr::mutate(
    Manual_Doublet = if_else(Scrublet_Score >= Manual_Threshold, TRUE, FALSE)
  )

# 4. Verify the new, realistic doublet rates per sample
summary_stats <- adjusted_doublets %>%
  dplyr::group_by(Pool_ID, Sample_ID) %>%
  dplyr::summarise(
    Total_Cells = n(),
    Doublets = sum(Manual_Doublet),
    Percent_Doublet = (Doublets / Total_Cells) * 100,
    .groups = "drop"
  )

print(summary_stats)

library(scDblFinder)

library(Seurat)
library(tidyverse)
library(patchwork)

# 1. Format the doublet data to match Seurat's internal cell IDs
# Assuming 'adjusted_doublets' is the dataframe from the previous step
doublets_for_seurat <- adjusted_doublets %>%
  dplyr::mutate(Seurat_Barcode = paste0(Pool_ID, "_", Sample_ID, "_", Barcode)) %>%
  tibble::column_to_rownames("Seurat_Barcode") %>%
  # Only keep the numerical score and the True/False call for the metadata
  dplyr::select(Scrublet_Score, Manual_Doublet) 

# 2. Add the clean metadata to your master Seurat object
master_seurat <- AddMetaData(master_seurat, metadata = doublets_for_seurat)

# 3. Standard Preprocessing (Normalisation, Feature Selection, Scaling)
master_seurat <- NormalizeData(master_seurat)
master_seurat <- FindVariableFeatures(master_seurat, selection.method = "vst", nfeatures = 2000)
master_seurat <- ScaleData(master_seurat)

# 4. Dimensionality Reduction
master_seurat <- RunPCA(master_seurat, features = VariableFeatures(object = master_seurat))

# Run the UMAP (using the top 30 principal components)
# Note: If the 'RcppAnnoy' package is not installed, this will trigger the fallback warning
# but will still successfully generate the UMAP, just a bit slower.
master_seurat <- RunUMAP(master_seurat, dims = 1:30)

# 5. Visualise the Results
# Plot A: Distinct True/False calls, split by pool
plot_calls <- DimPlot(master_seurat, 
                      group.by = "Manual_Doublet", 
                      split.by = "Orig_Pool",
                      cols = c("FALSE" = "#E5E5E5", "TRUE" = "#D55E00"),
                      pt.size = 0.5) + 
              ggtitle("Manually Curated Doublets")

# Plot B: The continuous raw Scrublet scores 
plot_scores <- FeaturePlot(master_seurat, 
                           features = "Scrublet_Score", 
                           split.by = "Orig_Pool",
                           cols = c("#E5E5E5", "#990000"),
                           pt.size = 0.5)

# Display the plots stacked on top of each other
plot_calls / plot_scores

DimPlot(master_seurat, shuffle = TRUE,
                      group.by = "Manual_Doublet", 

                      pt.size = 0.5) 

library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
library(BiocParallel)
library(tidyverse)

# 1. Create a unique identifier for each 10x run
master_seurat$Pool_Sample <- paste0(master_seurat$Orig_Pool, "_", master_seurat$Sample_ID)

# 2. Collapse the fragmented Seurat v5 layers back into a single matrix
master_seurat <- JoinLayers(master_seurat)

# 3. Convert the Seurat object to a SingleCellExperiment object
sce <- as.SingleCellExperiment(master_seurat)

# 4. Set up parallel processing
# Adjust the workers to match your HPC job allocation
bp <- MulticoreParam(workers = 8)

# 5. Run scDblFinder (calculating independently per sample)
sce <- scDblFinder(sce, samples = "Pool_Sample", BPPARAM = bp)

# 6. Map the results directly back to the master Seurat object
master_seurat$scDblFinder_Class <- sce$scDblFinder.class
master_seurat$scDblFinder_Score <- sce$scDblFinder.score

# 7. Generate a summary table to check the new estimated doublet rates
scDblFinder_summary <- master_seurat@meta.data %>%
  dplyr::group_by(Orig_Pool, Sample_ID) %>%
  dplyr::summarise(
    Total_Cells = n(),
    Doublets_Called = sum(scDblFinder_Class == "doublet"),
    Percent_Doublet = (Doublets_Called / Total_Cells) * 100,
    .groups = "drop"
  )

print(scDblFinder_summary)

library(Seurat)
library(tidyverse)
library(patchwork)

# 1. Extract the metadata for comparison
meta <- master_seurat@meta.data

# 2. Standardise the labels for a clean cross-tabulation
# Scrublet is currently TRUE/FALSE, scDblFinder is "singlet"/"doublet"
meta <- meta %>%
  dplyr::mutate(
    Scrublet_Call = if_else(Manual_Doublet == TRUE, "Doublet", "Singlet"),
    # Convert "doublet" to "Doublet" for exact matching
    scDblFinder_Call = stringr::str_to_title(scDblFinder_Class) 
  )

# 3. Generate the Confusion Matrix
confusion_matrix <- table(Scrublet = meta$Scrublet_Call, scDblFinder = meta$scDblFinder_Call)

cat("\n--- Confusion Matrix ---\n")
print(confusion_matrix)

# 4. Calculate overall agreement percentages
total_cells <- nrow(meta)
agree_singlets <- confusion_matrix["Singlet", "Singlet"]
agree_doublets <- confusion_matrix["Doublet", "Doublet"]

overall_agreement <- ((agree_singlets + agree_doublets) / total_cells) * 100

cat(paste0("\nOverall Agreement: ", round(overall_agreement, 2), "%\n"))

# 5. Create a Consensus category in the Seurat object for plotting
master_seurat@meta.data <- meta %>%
  dplyr::mutate(
    Consensus_Call = dplyr::case_when(
      Scrublet_Call == "Doublet" & scDblFinder_Call == "Doublet" ~ "Doublet (Both)",
      Scrublet_Call == "Doublet" & scDblFinder_Call == "Singlet" ~ "Doublet (Scrublet Only)",
      Scrublet_Call == "Singlet" & scDblFinder_Call == "Doublet" ~ "Doublet (scDblFinder Only)",
      TRUE ~ "Singlet (Both)"
    )
  )

# 6. Visualise the comparison on the UMAP
plot_scrublet <- DimPlot(master_seurat, 
                         group.by = "Scrublet_Call", 
                         cols = c("Singlet" = "#E5E5E5", "Doublet" = "#D55E00"), 
                         pt.size = 0.5) + 
                 ggtitle("Scrublet (Manual)")

plot_scDbl <- DimPlot(master_seurat, 
                      group.by = "scDblFinder_Call", 
                      cols = c("Singlet" = "#E5E5E5", "Doublet" = "#0072B2"), 
                      pt.size = 0.5) + 
              ggtitle("scDblFinder")

# Plotting the consensus highlights exactly where the algorithms disagree
plot_consensus <- DimPlot(master_seurat, 
                          group.by = "Consensus_Call", 
                          cols = c("Singlet (Both)" = "#E5E5E5", 
                                   "Doublet (Both)" = "#009E73", 
                                   "Doublet (Scrublet Only)" = "#D55E00", 
                                   "Doublet (scDblFinder Only)" = "#0072B2"), 
                          pt.size = 0.5) + 
                  ggtitle("Consensus Map")

# Display the plots: Top row individual algorithms, bottom row consensus
(plot_scrublet | plot_scDbl) / plot_consensus

master_seurat_clean <- subset(master_seurat, subset = scDblFinder_Class == "singlet")

# Optional: Clean up the metadata by dropping the old Scrublet columns to keep the object lightweight
master_seurat_clean@meta.data <- master_seurat_clean@meta.data %>%
  dplyr::select(-Scrublet_Score, -Manual_Doublet, -Scrublet_Call, -Consensus_Call)

# Verify the final cell count
master_seurat_clean

library(Seurat)
library(patchwork)

# 1. Re-calculate normalisation and variable features on the clean data
master_seurat_clean <- NormalizeData(master_seurat_clean)
master_seurat_clean <- FindVariableFeatures(master_seurat_clean, selection.method = "vst", nfeatures = 2000)

# 2. Scale the data
# This step standardises the expression values for each gene across all cells
master_seurat_clean <- ScaleData(master_seurat_clean)

# 3. Linear Dimensionality Reduction (PCA)
master_seurat_clean <- RunPCA(master_seurat_clean, features = VariableFeatures(object = master_seurat_clean))

# 4. Clustering
# Using the top 30 principal components as established in the previous run.
# A resolution of 0.5 is a standard starting point for ~80,000 cells. 
# You can increase this (e.g., 0.8 or 1.0) later if you need finer sub-clusters.
master_seurat_clean <- FindNeighbors(master_seurat_clean, dims = 1:30)
master_seurat_clean <- FindClusters(master_seurat_clean, resolution = 0.5)

# 5. Non-linear Dimensionality Reduction (UMAP)
master_seurat_clean <- RunUMAP(master_seurat_clean, dims = 1:30)

# 6. Visualisation
# Plot A: The standard UMAP coloured and labelled by the new Seurat clusters
plot_clusters <- DimPlot(master_seurat_clean, 
                         reduction = "umap", 
                         label = TRUE, 
                         pt.size = 0.5) + 
                 ggtitle("Cleaned UMAP Clusters")

# Plot B: The UMAP coloured by the original pools to check for batch effects
plot_pools <- DimPlot(master_seurat_clean, 
                      reduction = "umap", 
                      group.by = "Orig_Pool", 
                      pt.size = 0.1, 
                      alpha = 0.5) + 
              ggtitle("Distribution by AAV Pool")

# Display the plots side-by-side for comparison
plot_clusters | plot_pools



plot_clusters

FeaturePlot(master_seurat_clean, features = c("Cd163"),
                         reduction = "umap", 
                         label = TRUE, 
                         pt.size = 0.5) 

FeaturePlot(master_seurat_clean, features = c("Col3a1"),
                         reduction = "umap", 
                         label = TRUE, 
                         pt.size = 0.5) 

library(Seurat)
library(patchwork)

# 1. Subset the master object based strictly on your manual Scrublet calls
# We use the 'Scrublet_Call' column generated during the consensus step
master_seurat_scrublet <- subset(master_seurat, subset = Scrublet_Call == "Singlet")

# Optional: Clean up the metadata so it doesn't carry the scDblFinder data into this specific object
master_seurat_scrublet@meta.data <- master_seurat_scrublet@meta.data %>%
  dplyr::select(-scDblFinder_Class, -scDblFinder_Score, -Consensus_Call)

# 2. Re-calculate normalisation and variable features on this specific subset
master_seurat_scrublet <- NormalizeData(master_seurat_scrublet)
master_seurat_scrublet <- FindVariableFeatures(master_seurat_scrublet, selection.method = "vst", nfeatures = 2000)

# 3. Scale the data
master_seurat_scrublet <- ScaleData(master_seurat_scrublet)

# 4. Dimensionality Reduction & Clustering
master_seurat_scrublet <- RunPCA(master_seurat_scrublet, features = VariableFeatures(object = master_seurat_scrublet))
master_seurat_scrublet <- FindNeighbors(master_seurat_scrublet, dims = 1:30)
master_seurat_scrublet <- FindClusters(master_seurat_scrublet, resolution = 0.5)

# 5. Non-linear Dimensionality Reduction (UMAP)
master_seurat_scrublet <- RunUMAP(master_seurat_scrublet, dims = 1:30)

# 6. Visualisation
# Plot A: The standard UMAP coloured and labelled by the new Seurat clusters
plot_scrublet_clusters <- DimPlot(master_seurat_scrublet, 
                                  reduction = "umap", 
                                  label = TRUE, 
                                  pt.size = 0.5) + 
                          ggtitle("Scrublet-Cleaned UMAP Clusters")

# Plot B: The UMAP coloured by the original pools to check for batch effects
plot_scrublet_pools <- DimPlot(master_seurat_scrublet, 
                               reduction = "umap", 
                               group.by = "Orig_Pool", 
                               pt.size = 0.1, 
                               alpha = 0.5) + 
                       ggtitle("Distribution by AAV Pool")

# Display the plots side-by-side for comparison
plot_scrublet_clusters | plot_scrublet_pools

plot_scrublet_clusters



library(ggplot2)
library(tidyr)
library(Seurat)

plot_custom_qc <- function(seurat_obj, group.by = NULL) {
  # 1. Fetch the data using Seurat's internal extractor (works regardless of layers)
  vars <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  if (!is.null(group.by)) {
    vars <- c(vars, group.by)
  }
  
  plot_data <- FetchData(seurat_obj, vars = vars)
  
  # 2. Reshape into long format for ggplot
  plot_data_long <- plot_data %>%
    pivot_longer(cols = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                 names_to = "Metric", 
                 values_to = "Value")
  
  # 3. Build the plot
  # We use facet_wrap with 'scales = "free"' because the Y-axis ranges 
  # are completely different for counts vs percentages.
  p <- ggplot(plot_data_long, aes(x = "Cell", y = Value, fill = "Cell")) +
    geom_violin(trim = TRUE, scale = "width", fill = "#69b3a2", alpha = 0.7) +
    geom_boxplot(width = 0.1, color = "black", alpha = 0.5) +
    facet_wrap(~Metric, scales = "free", ncol = 3) +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = "none"
    ) +
    labs(y = "Value", title = "Quality Control Metrics")
  
  return(p)
}

my_qc_plot <- plot_custom_qc(master_seurat_clean)

my_qc_plot





library(Seurat)
library(patchwork)

# 1. Subset the master object based strictly on your manual Scrublet calls
# We use the 'Scrublet_Call' column generated during the consensus step
master_seurat_scrublet <- subset(master_seurat, subset = Scrublet_Call == "Singlet")

# Optional: Clean up the metadata so it doesn't carry the scDblFinder data into this specific object
master_seurat_scrublet@meta.data <- master_seurat_scrublet@meta.data %>%
  dplyr::select(-scDblFinder_Class, -scDblFinder_Score, -Consensus_Call)

# 2. Re-calculate normalisation and variable features on this specific subset
master_seurat_scrublet <- NormalizeData(master_seurat_scrublet)
master_seurat_scrublet <- FindVariableFeatures(master_seurat_scrublet, selection.method = "vst", nfeatures = 2000)

# 3. Scale the data
master_seurat_scrublet <- ScaleData(master_seurat_scrublet)

# 4. Dimensionality Reduction & Clustering
master_seurat_scrublet <- RunPCA(master_seurat_scrublet, features = VariableFeatures(object = master_seurat_scrublet))
master_seurat_scrublet <- FindNeighbors(master_seurat_scrublet, dims = 1:30)
master_seurat_scrublet <- FindClusters(master_seurat_scrublet, resolution = 0.5)

# 5. Non-linear Dimensionality Reduction (UMAP)
master_seurat_scrublet <- RunUMAP(master_seurat_scrublet, dims = 1:30)

# 6. Visualisation
# Plot A: The standard UMAP coloured and labelled by the new Seurat clusters
plot_scrublet_clusters <- DimPlot(master_seurat_scrublet, 
                                  reduction = "umap", 
                                  label = TRUE, 
                                  pt.size = 0.5) + 
                          ggtitle("Scrublet-Cleaned UMAP Clusters")

# Plot B: The UMAP coloured by the original pools to check for batch effects
plot_scrublet_pools <- DimPlot(master_seurat_scrublet, 
                               reduction = "umap", 
                               group.by = "Orig_Pool", 
                               pt.size = 0.1, 
                               alpha = 0.5) + 
                       ggtitle("Distribution by AAV Pool")

# Display the plots side-by-side for comparison
plot_scrublet_clusters | plot_scrublet_pools

plot_scrublet_clusters

############################################################
# Marker identification + first-pass cell type annotation
# Object: master_seurat_scrublet
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(tibble)
  library(stringr)
  library(purrr)
})

############################################################
# 0. Settings
############################################################



outdir <- "Scrublet_Cleaned_Marker_Annotation"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Cluster column to use
# Change this only if your clustering is stored under a different metadata column


Idents(master_seurat_scrublet) <- "seurat_clusters"

############################################################
# 1. Choose assay for marker testing
############################################################



# Prefer SCT if present and currently active; otherwise use RNA
# This preserves your likely existing workflow

  
  message("Using RNA assay for marker testing.")
  DefaultAssay(master_seurat_scrublet) <- "RNA"
  marker_assay <- "RNA"

############################################################
# 2. Identify cluster markers
############################################################

markers_all <- FindAllMarkers(
  object = master_seurat_scrublet,
  assay = marker_assay,
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.10,
  logfc.threshold = 0.5,
  return.thresh = 0.05,
  verbose = TRUE
)

markers_all <- markers_all %>%
  arrange(cluster, desc(p_val_adj))

write_csv(
  markers_all,
  file.path("All_Cluster_Markers_Scrublet_Cleaned.csv")
)

head(markers_all)

setwd(publication_path("DNAP1_SNRNA_ROOT", "Process/"))

saveRDS(master_seurat_scrublet , "CellBedner_Doublet_removed.rds")

master_seurat_scrublet$seurat_clusters

############################################################
# 3. Top markers per cluster
############################################################
fc_col ='avg_log2FC'
cluster_col = 'seurat_clusters'
top10_markers <- markers_all %>%
  group_by(cluster) %>%
  slice_max(order_by = .data[[fc_col]], n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(as.numeric(as.character(cluster)), desc(.data[[fc_col]]))



write_csv(
  top10_markers,
  file.path( "Top10_Markers_Per_Cluster_Scrublet_Cleaned.csv")
)


print(top10_markers %>% select(cluster, gene, all_of(fc_col), pct.1, pct.2, p_val_adj), n = Inf)

############################################################
# 4. Canonical marker panels for heart annotation
############################################################

marker_sets <- list(
  Cardiomyocyte = c(
    "Tnnt2", "Tnni3", "Actc1", "Myh6", "Myh7", "Ryr2", "Pln"
  ),
  
  Fibroblast = c(
    "Col1a1", "Col1a2", "Col3a1", "Dcn", "Lum", "Pdgfra", "Cthrc1"
  ),
  
  Endothelial = c(
    "Pecam1", "Cdh5", "Kdr", "Emcn", "Klf2", "Klf4", "Esam"
  ),
  
  Lymphatic_EC = c(
    "Prox1", "Lyve1", "Pdpn", "Flt4", "Ccl21a"
  ),
  
  Endocardial_EC = c(
    "Npr3", "Nrg1", "Plvap", "Gja5"
  ),
  
  Coronary_EC = c(
    "Fabp4", "Aplnr", "Car4", "Cd36", "Kdr"
  ),
  
  Smooth_Muscle = c(
    "Acta2", "Tagln", "Myh11", "Cnn1", "Myl9"
  ),
  
  Pericytes = c(
    "Rgs5", "Pdgfrb", "Cspg4", "Kcnj8", "Abcc9", "Notch3"
  ),
  
  Macrophage = c(
    "Adgre1", "Cd68", "C1qa", "C1qb", "C1qc", "Lyz2", "Lgals3"
  ),
  
  Monocyte_Inflammatory = c(
    "Ly6c2", "S100a8", "S100a9", "Ccr2", "Plac8", "Il1b"
  ),
  
  Dendritic = c(
    "Flt3", "Clec9a", "Xcr1", "Itgax", "Cd74", "H2-Ab1"
  ),
  
  T_cell = c(
    "Cd3d", "Cd3e", "Trac", "Ltb", "Il7r"
  ),
  
  B_cell = c(
    "Ms4a1", "Cd79a", "Cd79b", "Cd74", "H2-Aa"
  ),
  
  NK = c(
    "Nkg7", "Klrb1c", "Gzmb", "Prf1", "Ccl5"
  ),
  
  Neuronal = c(
    "Tubb3", "Snap25", "Syt1", "Elavl3", "Rbfox3"
  ),
  
  Mesothelial = c(
    "Upk3b", "Wt1", "Msln", "Krt19", "Aldh1a2"
  ),
  
  Cycling = c(
    "Mki67", "Top2a", "Cenpf", "Birc5", "Ube2c"
  )
)

# Keep only genes present in object
marker_sets_present <- map(
  marker_sets,
  ~ intersect(.x, rownames(master_seurat_scrublet))
)

marker_sets_present <- marker_sets_present[
  lengths(marker_sets_present) > 0
]

marker_gene_table <- tibble(
  CellType = rep(names(marker_sets_present), lengths(marker_sets_present)),
  Gene = unlist(marker_sets_present, use.names = FALSE)
)

write_csv(
  marker_gene_table,
  file.path( "Marker_Panel_Genes_Present_In_Object.csv")
)

############################################################
# 5. Dot plot of major marker genes
############################################################

dotplot_features <- unique(unlist(marker_sets_present, use.names = FALSE))

p_dot_all <- DotPlot(
  master_seurat_scrublet,
  features = dotplot_features,
  assay = marker_assay,
  group.by = cluster_col
) +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold")
  ) +
  ggtitle("Canonical Heart and Immune Markers by Cluster")

ggsave(
  file.path( "DotPlot_All_Canonical_Markers_By_Cluster.pdf"),
  p_dot_all,
  width = 26,
  height = 10
)

############################################################
# 6. More readable grouped dot plots
############################################################

plot_marker_group <- function(feature_vec, title_text, filename, width = 10, height = 6) {
  
  genes_present <- intersect(feature_vec, rownames(master_seurat_scrublet))
  
  if (length(genes_present) == 0) {
    message("Skipping ", title_text, ": no genes present.")
    return(NULL)
  }
  
  p <- DotPlot(
    master_seurat_scrublet,
    features = genes_present,
    assay = marker_assay,
    group.by = cluster_col
  ) +
    RotatedAxis() +
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(face = "bold")
    ) +
    ggtitle(title_text)
  
  ggsave(
    file.path(filename),
    p,
    width = width,
    height = height
  )
  
  return(p)
}

p_core_lineages <- plot_marker_group(
  feature_vec = unique(c(
    marker_sets$Cardiomyocyte,
    marker_sets$Fibroblast,
    marker_sets$Endothelial,
    marker_sets$Smooth_Muscle,
    marker_sets$Pericytes
  )),
  title_text = "Core Cardiac Stromal and Vascular Markers",
  filename = "DotPlot_Core_Cardiac_Lineages.pdf",
  width = 18,
  height = 8
)

p_immune <- plot_marker_group(
  feature_vec = unique(c(
    marker_sets$Macrophage,
    marker_sets$Monocyte_Inflammatory,
    marker_sets$Dendritic,
    marker_sets$T_cell,
    marker_sets$B_cell,
    marker_sets$NK
  )),
  title_text = "Immune Cell Markers",
  filename = "DotPlot_Immune_Lineages.pdf",
  width = 20,
  height = 8
)

p_specialised <- plot_marker_group(
  feature_vec = unique(c(
    marker_sets$Lymphatic_EC,
    marker_sets$Endocardial_EC,
    marker_sets$Coronary_EC,
    marker_sets$Mesothelial,
    marker_sets$Neuronal,
    marker_sets$Cycling
  )),
  title_text = "Specialised and Rare Population Markers",
  filename = "DotPlot_Specialised_Rare_Lineages.pdf",
  width = 18,
  height = 8
)

############################################################
# 7. Heatmap of top cluster markers
############################################################

top_heatmap_genes <- top10_markers %>%
  pull(gene) %>%
  unique()

top_heatmap_genes <- intersect(top_heatmap_genes, rownames(master_seurat_scrublet))

p_heatmap <- DoHeatmap(
  master_seurat_scrublet,
  features = top_heatmap_genes,
  assay = marker_assay,
  group.by = cluster_col,
  size = 3
) +
  ggtitle("Top 10 Positive Markers per Cluster")

ggsave(
  file.path( "Heatmap_Top10_Cluster_Markers.pdf"),
  p_heatmap,
  width = 18,
  height = 24
)

############################################################
# 8. Module score-based first-pass annotation
############################################################
# AddModuleScore calculates cell-level expression program scores
# for each supplied gene set, then we summarise by cluster.

plot_scrublet_clusters
marker_sets <- list(
  Cardiomyocyte = c(
    "Tnnt2", "Tnni3", "Actc1", "Myh6", "Myh7", "Ryr2", "Pln"
  ),
  
  Fibroblast = c(
    "Col1a1", "Col1a2", "Col3a1", "Dcn", "Lum", "Pdgfra", "Cthrc1"
  ),
  
  Endothelial = c(
    "Pecam1", "Cdh5", "Kdr", "Emcn", "Klf2", "Klf4", "Esam"
  ),
  
  Lymphatic_EC = c(
    "Prox1", "Lyve1", "Pdpn", "Flt4", "Ccl21a"
  ),
  
  Endocardial_EC = c(
    "Npr3", "Nrg1", "Plvap", "Gja5"
  ),
  
  Coronary_EC = c(
    "Fabp4", "Aplnr", "Car4", "Cd36", "Kdr"
  ),
  
  Smooth_Muscle = c(
    "Acta2", "Tagln", "Myh11", "Cnn1", "Myl9"
  ),
  
  Pericytes = c(
    "Rgs5", "Pdgfrb", "Cspg4", "Kcnj8", "Abcc9", "Notch3"
  ),
  
  Macrophage = c(
    "Adgre1", "Cd68", "C1qa", "C1qb", "C1qc", "Lyz2", "Lgals3"
  ),
  
  Monocyte_Inflammatory = c(
    "Ly6c2", "S100a8", "S100a9", "Ccr2", "Plac8", "Il1b"
  ),
  
  Dendritic = c(
    "Flt3", "Clec9a", "Xcr1", "Itgax", "Cd74", "H2-Ab1"
  ),
  
  T_cell = c(
    "Cd3d", "Cd3e", "Trac", "Ltb", "Il7r"
  ),
  
  B_cell = c(
    "Ms4a1", "Cd79a", "Cd79b", "Cd74", "H2-Aa"
  ),
  
  NK = c(
    "Nkg7", "Klrb1c", "Gzmb", "Prf1", "Ccl5"
  ),
  
  Neuronal = c(
    "Tubb3", "Snap25", "Syt1", "Elavl3", "Rbfox3"
  ),
  
  Mesothelial = c(
    "Upk3b", "Wt1", "Msln", "Krt19", "Aldh1a2"
  ),
  
  Cycling = c(
    "Mki67", "Top2a", "Cenpf", "Birc5", "Ube2c"
  )
)

marker_sets <- list(
  Cardiomyocyte = c(
    0,8,13,18,19
  ),
  
  Fibroblast = c(
    1,3
  ),
  
  Endothelial = c(
   4,6,7,9,12
  ),
  Smooth_Muscle = c(
   14
  ),
  
  Pericytes = c(
   5
  ),
  
  Macrophage = c(
   2
  ),
  
  Neutrophil = c(
    15
  ),
  
  Dendritic = c(
   17
  ),
  
  T_cell = c(
    11
  ),
  
  B_cell = c(
   10
  ),
  
  Atrial_cardio = c(
    19
  ),
  
  adipocytes = c(
    21
  ),
  
  Mesothelial = c(
    16
  ),
  
  Cycling = c(
    20
  )
)

DT::datatable(markers_all, filter = 'top')

marker_sets_refined <- list(
  Ventricular_CM = c(0),
  Male_biased_CM = c(8),
  Atrial_conduction_like_CM = c(13),
  Stress_Remodelling_CM = c(18),
  Atrial_CM = c(19),

  Fibroblast = c(1, 3),

  Capillary_Coronary_EC = c(4),
  Venous_Endocardial_like_EC = c(6),
  Transitional_EC = c(7),
  Arterial_Coronary_EC = c(9),
  Lymphatic_EC = c(12),

  Pericytes = c(5),
  Smooth_Muscle = c(14),

  Macrophage = c(2),
  Neutrophil_Inflammatory_Myeloid = c(15),
  Dendritic_cDC1 = c(17),

  T_cell = c(11),
  B_cell = c(10),

  Mesothelial = c(16),
  Cycling_Myeloid = c(20),
  Adipocyte = c(21),
  Schwann_Peripheral_Glia = c(22)
)

# If you want broad labels only:
marker_sets_broad <- list(
  Cardiomyocyte = c(0, 8, 13, 18, 19),
  Fibroblast = c(1, 3),
  Endothelial = c(4, 6, 7, 9, 12),
  Pericytes = c(5),
  Smooth_Muscle = c(14),
  Macrophage = c(2),
  Neutrophil = c(15),
  Dendritic = c(17),
  T_cell = c(11),
  B_cell = c(10),
  Mesothelial = c(16),
  Cycling = c(20),
  Adipocyte = c(21),
  Schwann_Peripheral_Glia = c(22)
)

# Useful checks for the cardiomyocyte split.
# Replace master_seurat with your object name if different.

cm_features <- c(
  # broad CM
  "Tnnt2", "Actc1", "Ttn", "Myh6", "Myh7", "Ryr2",

  # ventricular / mature CM
  "Myl2", "Myl3", "Lrrc10", "Mypn", "Pdk4", "Plin5", "Ucp3", "Ppargc1a",

  # atrial working CM
  "Myl7", "Myl4", "Sln", "Nppa", "Bmp10", "Mybphl",

  # atrial/conduction-like markers
  "Rxfp1", "Fgf12", "Ryr3", "Tbx5", "Cacna2d2", "Cacna1h",
  "Kcnj3", "Hcn1", "Gja5", "Nfasc", "Adcy1",

  # stress/remodelling CM
  "Nppb", "Ankrd1", "Xirp1", "Gadd45g", "Cited4", "Rrad", "Hspb7",

  # sex-driven cluster check
  "Uty", "Ddx3y", "Kdm5d", "Eif2s3y"
)



DotPlot(
  master_seurat_scrublet,
  features = cm_features,
  group.by = "seurat_clusters"
)

# Add refined labels to the object.
###############################################################################
# Add refined cell type annotations safely
###############################################################################

cluster_to_celltype <- c(
  "0"  = "Ventricular_CM",
  "1"  = "Fibroblast",
  "2"  = "Macrophage",
  "3"  = "Fibroblast",
  "4"  = "Capillary_Coronary_EC",
  "5"  = "Pericytes",
  "6"  = "Venous_Endocardial_like_EC",
  "7"  = "Transitional_EC",
  "8"  = "Male_biased_CM",
  "9"  = "Arterial_Coronary_EC",
  "10" = "B_cell",
  "11" = "T_cell",
  "12" = "Lymphatic_EC",
  "13" = "Atrial_conduction_like_CM",
  "14" = "Smooth_Muscle",
  "15" = "Neutrophil_Inflammatory_Myeloid",
  "16" = "Mesothelial",
  "17" = "Dendritic_cDC1",
  "18" = "Stress_Remodelling_CM",
  "19" = "Atrial_CM",
  "20" = "Cycling_Myeloid",
  "21" = "Adipocyte",
  "22" = "Schwann_Peripheral_Glia"
)

# Check cluster column exists
stopifnot("seurat_clusters" %in% colnames(master_seurat_scrublet@meta.data))

# Convert cluster IDs to character
cluster_ids <- as.character(master_seurat_scrublet$seurat_clusters)

# Map cluster IDs to cell type labels
celltype_refined <- cluster_to_celltype[cluster_ids]

# Remove names inherited from cluster_to_celltype
celltype_refined <- unname(celltype_refined)

# Add cell barcodes as names so Seurat can match correctly
names(celltype_refined) <- colnames(master_seurat_scrublet)

# Check for missing annotations
missing_clusters <- sort(unique(cluster_ids[is.na(celltype_refined)]))

if (length(missing_clusters) > 0) {
  stop(
    "These clusters are missing from cluster_to_celltype: ",
    paste(missing_clusters, collapse = ", ")
  )
}

# Add to metadata
master_seurat <- AddMetaData(
  object = master_seurat_scrublet,
  metadata = celltype_refined,
  col.name = "celltype_refined"
)

# Confirm it worked
table(master_seurat_scrublet$seurat_clusters, master_seurat_scrublet$celltype_refined)

DimPlot(
  master_seurat_scrublet,
  group.by = "celltype_refined",
  label = TRUE,
  repel = TRUE
)


DimPlot(
  master_seurat_scrublet,
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
)

# Directly compare CM clusters 13 and 19 against other cardiomyocytes.
cm_clusters <- c("0", "8", "13", "18", "19")

cm_obj <- subset(
  master_seurat_scrublet,
  subset = seurat_clusters %in% cm_clusters
)

Idents(cm_obj) <- "seurat_clusters"

markers_13_vs_19 <- FindMarkers(
  cm_obj,
  ident.1 = "13",
  ident.2 = "19",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

markers_13_vs_other_cm <- FindMarkers(
  cm_obj,
  ident.1 = "13",
  ident.2 = setdiff(cm_clusters, "13"),
  logfc.threshold = 0.25,
  min.pct = 0.1
)

markers_19_vs_other_cm <- FindMarkers(
  cm_obj,
  ident.1 = "19",
  ident.2 = setdiff(cm_clusters, "19"),
  logfc.threshold = 0.25,
  min.pct = 0.1
)

head(markers_13_vs_19[order(markers_13_vs_19$avg_log2FC, decreasing = TRUE), ], 30)
head(markers_13_vs_19[order(markers_13_vs_19$avg_log2FC, decreasing = FALSE), ], 30)

write.csv(markers_13_vs_19, "CM_cluster13_vs_cluster19_markers.csv")
write.csv(markers_13_vs_other_cm, "CM_cluster13_vs_other_CM_markers.csv")
write.csv(markers_19_vs_other_cm, "CM_cluster19_vs_other_CM_markers.csv")

###############################################################################
# Add mouse-level metadata to Seurat object using Sample_number
###############################################################################

library(dplyr)
library(stringr)
library(tibble)

# Check existing metadata columns
colnames(master_seurat_scrublet@meta.data)

# This assumes your current Seurat metadata has Sample_ID like:
# Sample_1, Sample_2, Sample_3, etc.
# If your column is named differently, change Sample_ID below.

master_seurat_scrublet$Sample_number <- str_extract(
  as.character(master_seurat_scrublet$Sample_ID),
  "\\d+$"
) %>%
  as.integer()

# Mouse/sample metadata
sample_metadata <- tibble::tribble(
  ~Sample_number, ~Mouse_ID,                 ~Group,                ~Experimental_Sample_ID, ~Sex,      ~Tube_label,       ~Process_time,
  1,              "000330//94//30300194",    "AAV9-cTnT",           "EW#6_1",                "Male",    "EW#6_1_Heart",    "9:52am",
  2,              "000333//44//30300244",    "AAV9-AP-1_Dosage_2",  "EW#6_2",                "Female",  "EW#6_2_Heart",    "8:57am",
  3,              "000334//45//30300245",    "AAV9-AP-1_Dosage_2",  "EW#6_3",                "Female",  "EW#6_3_Heart",    "9:25am",
  4,              "000335//79//30300379",    "AAV9-AP-1_Dosage_2",  "EW#6_4",                "Female",  "EW#6_4_Heart",    "9:55am",
  5,              "000337//33//30300510",    "AAV9-cTnT",           "EW#6_5",                "Female",  "EW#6_5_Heart",    "8:59am",
  6,              "000338//1//30300515",     "AAV9-cTnT",           "EW#6_6",                "Female",  "EW#6_6_Heart",    "9:27am"
)

# Check sample numbers present in Seurat object
table(master_seurat_scrublet$Sample_ID, master_seurat_scrublet$Sample_number)

# Check for Seurat samples missing from metadata
missing_from_metadata <- setdiff(
  sort(unique(master_seurat_scrublet$Sample_number)),
  sample_metadata$Sample_number
)

if (length(missing_from_metadata) > 0) {
  stop(
    "These Sample_number values are present in the Seurat object but missing from sample_metadata: ",
    paste(missing_from_metadata, collapse = ", ")
  )
}

# Check for metadata samples not present in Seurat object
missing_from_seurat <- setdiff(
  sample_metadata$Sample_number,
  sort(unique(master_seurat$Sample_number))
)

if (length(missing_from_seurat) > 0) {
  warning(
    "These Sample_number values are in sample_metadata but not present in the Seurat object: ",
    paste(missing_from_seurat, collapse = ", ")
  )
}

# Create metadata dataframe with cell barcodes preserved
meta_new <- master_seurat_scrublet@meta.data %>%
  rownames_to_column("cell_barcode") %>%
  left_join(sample_metadata, by = "Sample_number") %>%
  column_to_rownames("cell_barcode")

# Safety check: cell order must match exactly
stopifnot(identical(rownames(meta_new), colnames(master_seurat_scrublet)))

# Add columns back to Seurat metadata
cols_to_add <- c(
  "Sample_number",
  "Mouse_ID",
  "Group",
  "Experimental_Sample_ID",
  "Sex",
  "Tube_label",
  "Process_time"
)

# master_seurat_scrublet <- AddMetaData(
#   object = master_seurat_scrublet,
#   metadata = meta_new[, cols_to_add, drop = FALSE]
# )
colnames(master_seurat_scrublet@meta.data)
# Check result
table(master_seurat_scrublet$Sample_number, master_seurat_scrublet$Orig_Pool)
table(master_seurat_scrublet$Sample_number, master_seurat_scrublet$Sex)
table(master_seurat_scrublet$Sample_ID, master_seurat_scrublet$Experimental_Sample_ID)

head(master_seurat_scrublet@meta.data[, c(
  "Orig_Pool",
  "Sample_ID",
  "Sample_number",
  "Mouse_ID",
  "Group",
  "Experimental_Sample_ID",
  "Sex",
  "Tube_label",
  "Process_time"
)])

###############################################################################
# QC and marker checks for CM clusters 13 vs 19
###############################################################################

library(Seurat)
library(dplyr)
library(ggplot2)

# Define CM clusters
cm_clusters <- c("0", "8", "13", "18", "19")

# Make sure cluster IDs are character
master_seurat_scrublet$cluster_chr <- as.character(master_seurat_scrublet$seurat_clusters)

# Subset cardiomyocyte-labelled clusters
cm_obj <- subset(
  master_seurat_scrublet,
  subset = cluster_chr %in% cm_clusters
)

Idents(master_seurat_scrublet) <- "cluster_chr"
Idents(cm_obj) <- "cluster_chr"

###############################################################################
# 1. QC by cluster
###############################################################################

qc_features <- intersect(
  c("nCount_RNA", "nFeature_RNA", "percent.mt", "Scrublet_score", "doublet_score"),
  colnames(master_seurat_scrublet@meta.data)
)



# QC table
qc_table <- master_seurat_scrublet@meta.data %>%
  mutate(cluster_chr = as.character(seurat_clusters)) %>%
  filter(cluster_chr %in% cm_clusters) %>%
  group_by(cluster_chr) %>%
  summarise(
    n_cells = n(),
    median_nCount_RNA = median(nCount_RNA, na.rm = TRUE),
    median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE),
    median_percent_mt = if ("percent.mt" %in% colnames(master_seurat_scrublet@meta.data)) median(percent.mt, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )

print(qc_table)

###############################################################################
# 2. Sample composition
###############################################################################

if ("Sample_ID" %in% colnames(master_seurat_scrublet@meta.data)) {
  print(table(master_seurat_scrublet$cluster_chr, master_seurat_scrublet$Sample_ID))
  
  sample_prop <- master_seurat_scrublet@meta.data %>%
    mutate(cluster_chr = as.character(seurat_clusters)) %>%
    filter(cluster_chr %in% cm_clusters) %>%
    dplyr::count(cluster_chr, Sample_ID) %>%
    group_by(cluster_chr) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()
  
  print(sample_prop)
  
  print(
    ggplot(sample_prop, aes(x = cluster_chr, y = prop, fill = Sample_ID)) +
      geom_col() +
      theme_classic() +
      labs(
        x = "Cluster",
        y = "Proportion within cluster",
        title = "Sample composition of cardiomyocyte clusters"
      )
  )
}

###############################################################################
# 3. DotPlot: CM, atrial/conduction, stress, contamination markers
###############################################################################

marker_check <- c(
  # broad CM
  "Tnnt2", "Ttn", "Actc1", "Myh6", "Myh7", "Ryr2",
  
  # ventricular
  "Myl2", "Myl3", "Lrrc10", "Mypn", "Ppargc1a",
  
  # atrial
  "Myl4", "Myl7", "Sln", "Nppa", "Bmp10", "Mybphl",
  
  # conduction / ion-channel-like
  "Rxfp1", "Fgf12", "Ryr3", "Tbx5", "Cacna2d2", "Cacna1h",
  "Kcnj3", "Hcn1", "Gja5", "Gja3", "Rnf207", "Wnk2", "Dsc2",
  
  # stress/remodelling
  "Nppb", "Ankrd1", "Xirp1", "Gadd45g", "Cited4", "Rrad", "Hspb7",
  
  # endothelial contamination
  "Pecam1", "Cdh5", "Kdr", "Ramp2", "Gpihbp1",
  
  # fibroblast / mural / adipocyte contamination
  "Dcn", "Bgn", "Mgp", "Vim", "Fabp4", "Col1a1", "Col1a2",
  
  # myeloid / RBC
  "Lyz2", "C1qa", "C1qb", "Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt"
)

marker_check <- marker_check[marker_check %in% rownames(master_seurat_scrublet)]

print(
  DotPlot(
    master_seurat_scrublet,
    features = marker_check,
    group.by = "cluster_chr",
    idents = cm_clusters
  ) +
    RotatedAxis() +
    theme_classic() +
    labs(title = "Cardiomyocyte cluster marker validation")
)

###############################################################################
# 4. FeaturePlot for cluster 13 vs 19 markers
###############################################################################

FeaturePlot(
  master_seurat_scrublet,
  features = c(
    "Fgf12", "Rxfp1", "Ryr3", "Tbx5", "Kcnj3", "Rnf207",
    "Myl7", "Sln", "Nppa", "Myl4",
    "Pecam1", "Cdh5", "Kdr", "Dcn", "Lyz2", "Hba-a1"
  )[c(
    "Fgf12", "Rxfp1", "Ryr3", "Tbx5", "Kcnj3", "Rnf207",
    "Myl7", "Sln", "Nppa", "Myl4",
    "Pecam1", "Cdh5", "Kdr", "Dcn", "Lyz2", "Hba-a1"
  ) %in% rownames(master_seurat_scrublet)],
  ncol = 4,
  order = TRUE
)

###############################################################################
# 5. Co-expression contamination check in cluster 19
###############################################################################

cluster19 <- subset(master_seurat_scrublet, subset = cluster_chr == "19")

contam_markers <- c("Pecam1", "Cdh5", "Kdr", "Dcn", "Bgn", "Mgp", "Lyz2", "Hba-a1", "Hbb-bs")
cm_markers <- c("Tnnt2", "Ttn", "Actc1", "Myh6", "Myl7", "Sln", "Nppa")

contam_markers <- contam_markers[contam_markers %in% rownames(cluster19)]
cm_markers <- cm_markers[cm_markers %in% rownames(cluster19)]

if (length(contam_markers) > 0 && length(cm_markers) > 0) {
  expr <- FetchData(cluster19, vars = c(cm_markers, contam_markers))
  
  contam_positive <- rowSums(expr[, contam_markers, drop = FALSE] > 0) > 0
  cm_positive <- rowSums(expr[, cm_markers, drop = FALSE] > 0) > 0
  
  cat("\nCluster 19 cells:", nrow(expr), "\n")
  cat("CM marker positive:", sum(cm_positive), "(", round(mean(cm_positive) * 100, 2), "% )\n")
  cat("Contamination marker positive:", sum(contam_positive), "(", round(mean(contam_positive) * 100, 2), "% )\n")
  cat("Both CM and contamination marker positive:", sum(cm_positive & contam_positive), "(", round(mean(cm_positive & contam_positive) * 100, 2), "% )\n")
  
  contam_summary <- data.frame(
    marker = contam_markers,
    pct_positive_cluster19 = colMeans(expr[, contam_markers, drop = FALSE] > 0) * 100
  ) %>%
    arrange(desc(pct_positive_cluster19))
  
  print(contam_summary)
}

###############################################################################
# 6. Temporary annotation
###############################################################################

cluster_to_celltype_refined <- c(
  "0"  = "Ventricular_CM",
  "8"  = "Male_biased_CM",
  "13" = "Atrial_conduction_like_CM",
  "18" = "Stress_Remodelling_CM",
  "19" = "Atrial_metabolic_CM_check_quality"
)

cm_obj$cm_subtype_check <- unname(cluster_to_celltype_refined[as.character(cm_obj$cluster_chr)])

print(
  DimPlot(
    cm_obj,
    group.by = "cm_subtype_check",
    label = TRUE,
    repel = TRUE
  )
)

###############################################################################
# Cluster 19 contamination severity check
###############################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)

cm_clusters <- c("0", "8", "13", "18", "19")

master_seurat_scrublet$cluster_chr <- as.character(master_seurat_scrublet$seurat_clusters)

cm_markers <- c("Tnnt2", "Ttn", "Actc1", "Myh6", "Myh7", "Ryr2", "Myl7", "Sln", "Nppa", "Myl4")
contam_markers <- c("Mgp", "Dcn", "Bgn", "Kdr", "Pecam1", "Cdh5", "Lyz2", "Hba-a1", "Hbb-bs")

cm_markers <- cm_markers[cm_markers %in% rownames(master_seurat_scrublet)]
contam_markers <- contam_markers[contam_markers %in% rownames(master_seurat_scrublet)]

check_obj <- subset(master_seurat_scrublet, subset = cluster_chr %in% cm_clusters)

# Add module scores
check_obj <- AddModuleScore(
  check_obj,
  features = list(cm_markers),
  name = "CM_score"
)

check_obj <- AddModuleScore(
  check_obj,
  features = list(contam_markers),
  name = "Contam_score"
)

# Seurat appends 1 to score names
check_obj$CM_score <- check_obj$CM_score1
check_obj$Contam_score <- check_obj$Contam_score1

###############################################################################
# 1. Compare contamination scores across CM clusters
###############################################################################



###############################################################################
# 2. DotPlot of CM vs contamination markers across CM clusters
###############################################################################

marker_check <- c(
  cm_markers,
  contam_markers
)

p2 <- DotPlot(
  check_obj,
  features = marker_check,
  group.by = "cluster_chr"
) +
  RotatedAxis() +
  theme_classic() +
  labs(title = "CM clusters: cardiomyocyte vs contaminant marker expression")

print(p2)

###############################################################################
# 3. Is cluster 19 only 1-UMI ambient noise?
###############################################################################

cluster19 <- subset(check_obj, subset = cluster_chr == "19")

expr_counts <- GetAssayData(cluster19, assay = DefaultAssay(cluster19), slot = "counts")

contam_present <- contam_markers[contam_markers %in% rownames(expr_counts)]

contam_summary <- lapply(contam_present, function(gene) {
  vals <- as.numeric(expr_counts[gene, ])
  data.frame(
    gene = gene,
    pct_gt0 = mean(vals > 0) * 100,
    pct_gt1 = mean(vals > 1) * 100,
    pct_gt2 = mean(vals > 2) * 100,
    mean_count_positive = ifelse(any(vals > 0), mean(vals[vals > 0]), NA_real_),
    max_count = max(vals)
  )
}) %>%
  bind_rows() %>%
  arrange(desc(pct_gt0))

print(contam_summary)

###############################################################################
# 4. Check whether contaminant-positive cluster 19 cells have high nCount/nFeature
###############################################################################

cluster19_meta <- cluster19@meta.data %>%
  mutate(
    cell_barcode = rownames(cluster19@meta.data)
  )

cluster19_expr <- FetchData(
  cluster19,
  vars = c(cm_markers, contam_markers)
)

cluster19_meta$CM_positive <- rowSums(cluster19_expr[, cm_markers, drop = FALSE] > 0) > 0
cluster19_meta$Contam_positive <- rowSums(cluster19_expr[, contam_markers, drop = FALSE] > 0) > 0

p3 <- ggplot(cluster19_meta, aes(x = Contam_positive, y = nCount_RNA)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  theme_classic() +
  labs(
    title = "Cluster 19 nCount_RNA by contaminant marker positivity",
    x = "Any contaminant marker detected",
    y = "nCount_RNA"
  )

print(p3)

p4 <- ggplot(cluster19_meta, aes(x = Contam_positive, y = nFeature_RNA)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  theme_classic() +
  labs(
    title = "Cluster 19 nFeature_RNA by contaminant marker positivity",
    x = "Any contaminant marker detected",
    y = "nFeature_RNA"
  )

print(p4)

###############################################################################
# 5. Optional: flag cluster 19 as low-confidence
###############################################################################

master_seurat_scrublet$celltype_refined_checked <- master_seurat_scrublet$celltype_refined

master_seurat_scrublet$celltype_refined_checked[
  master_seurat_scrublet$cluster_chr == "19"
] <- "Atrial_CM_mixed_low_confidence"

table(master_seurat_scrublet$cluster_chr, master_seurat_scrublet$celltype_refined_checked)

###############################################################################
# Remove low-confidence cluster 19, recluster, and reannotate from scratch
###############################################################################

library(Seurat)
library(dplyr)
library(ggplot2)

###############################################################################
# 1. Save original object before removing anything
###############################################################################

master_seurat_scrublet$cluster_chr_original <- as.character(master_seurat_scrublet$seurat_clusters)

saveRDS(
  master_seurat_scrublet,
  file = "master_seurat_scrublet_with_cluster19_flagged_before_reclustering.rds"
)

###############################################################################
# 2. Remove cluster 19
###############################################################################

master_clean <- subset(
  master_seurat_scrublet,
  subset = cluster_chr_original != "19"
)

cat("\nOriginal cluster counts:\n")
print(table(master_seurat_scrublet$cluster_chr_original))

cat("\nClean object cluster counts after removing original cluster 19:\n")
print(table(master_clean$cluster_chr_original))

###############################################################################
# 3. Clear old dimensional reductions, graphs, neighbours and cluster calls
###############################################################################

master_clean@reductions <- list()
master_clean@graphs <- list()
master_clean@neighbors <- list()

old_cluster_cols <- grep(
  pattern = "snn_res|seurat_clusters",
  x = colnames(master_clean@meta.data),
  value = TRUE
)

if (length(old_cluster_cols) > 0) {
  master_clean@meta.data[, old_cluster_cols] <- NULL
}

###############################################################################
# 4. Reprocess from raw RNA counts
###############################################################################

DefaultAssay(master_clean) <- "RNA"

master_clean <- NormalizeData(
  master_clean,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

master_clean <- FindVariableFeatures(
  master_clean,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = TRUE
)

master_clean <- ScaleData(
  master_clean,
  features = VariableFeatures(master_clean),
  verbose = TRUE
)

master_clean <- RunPCA(
  master_clean,
  features = VariableFeatures(master_clean),
  npcs = 50,
  verbose = TRUE
)

###############################################################################
# 5. Choose dims
###############################################################################

ElbowPlot(master_clean, ndims = 50)

# Adjust this after checking ElbowPlot.
dims_use <- 1:30

###############################################################################
# 6. Recluster and rerun UMAP
###############################################################################

master_clean <- FindNeighbors(
  master_clean,
  dims = dims_use,
  verbose = TRUE
)

master_clean <- FindClusters(
  master_clean,
  resolution = 0.5,
  verbose = TRUE
)

master_clean <- RunUMAP(
  master_clean,
  dims = dims_use,
  verbose = TRUE
)

###############################################################################
# 7. Check the new clustering
###############################################################################

DimPlot(
  master_clean,
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
)

if ("Sample_ID" %in% colnames(master_clean@meta.data)) {
  DimPlot(
    master_clean,
    group.by = "Sample_ID"
  )
}

if ("Group" %in% colnames(master_clean@meta.data)) {
  DimPlot(
    master_clean,
    group.by = "Group"
  )
}

###############################################################################
# 8. Recalculate markers for new clusters
###############################################################################

Idents(master_clean) <- "seurat_clusters"

markers_clean <- FindAllMarkers(
  master_clean,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25
)

markers_clean <- markers_clean %>%
  arrange(cluster, desc(avg_log2FC))

write.csv(
  markers_clean,
  file = "All_Cluster_Markers_cluster19_removed_reclustered.csv",
  row.names = FALSE
)

top20_markers_clean <- markers_clean %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20, with_ties = FALSE) %>%
  ungroup()

write.csv(
  top20_markers_clean,
  file = "Top20_Cluster_Markers_cluster19_removed_reclustered.csv",
  row.names = FALSE
)

###############################################################################
# 9. Marker validation DotPlot for major cell types
###############################################################################

marker_check <- c(
  # Cardiomyocyte
  "Tnnt2", "Ttn", "Actc1", "Myh6", "Myh7", "Ryr2",
  "Myl2", "Myl3", "Myl4", "Myl7", "Sln", "Nppa",
  "Rxfp1", "Fgf12", "Ryr3", "Tbx5", "Hcn1", "Kcnj3",
  "Nppb", "Ankrd1", "Xirp1",

  # Fibroblast
  "Dcn", "Lum", "Col1a1", "Col1a2", "Pdgfra",

  # Endothelial
  "Pecam1", "Cdh5", "Kdr", "Emcn", "Klf2",

  # Lymphatic
  "Prox1", "Pdpn", "Lyve1", "Flt4",

  # Pericyte / smooth muscle
  "Rgs5", "Pdgfrb", "Kcnj8", "Acta2", "Tagln", "Myh11",

  # Immune
  "Lyz2", "Adgre1", "C1qa", "C1qb", "Cd68",
  "S100a8", "S100a9", "Csf3r",
  "Xcr1", "Itgae", "Clec9a",
  "Cd3d", "Cd3e", "Nkg7",
  "Ms4a1", "Cd79a", "Cd74",

  # Mesothelial / adipocyte / Schwann
  "Upk3b", "Upk1b", "Muc16", "Wt1",
  "Plin1", "Cidec", "Pck1",
  "Sox10", "Plp1", "Mpz", "Pmp22",

  # Low-quality / contamination checks
  "Hba-a1", "Hbb-bs", "Mgp", "Bgn"
)

marker_check <- marker_check[marker_check %in% rownames(master_clean)]

DotPlot(
  master_clean,
  features = marker_check,
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  theme_classic() +
  labs(title = "Marker validation after removing cluster 19 and reclustering")

###############################################################################
# 10. Save clean reclustered object
###############################################################################

saveRDS(
  master_clean,
  file = "master_seurat_scrublet_cluster19_removed_reclustered.rds"
)

###############################################################################
# 11. Optional: recluster cardiomyocytes only after global reclustering
###############################################################################

cm_marker_genes <- c("Tnnt2", "Ttn", "Actc1", "Myh6", "Myh7", "Ryr2")
cm_marker_genes <- cm_marker_genes[cm_marker_genes %in% rownames(master_clean)]

# After inspecting markers_clean, replace this with the new CM cluster IDs.
# Example only:
# new_cm_clusters <- c("0", "7", "12", "15")
new_cm_clusters <- c()

if (length(new_cm_clusters) > 0) {
  master_clean$cluster_chr_new <- as.character(master_clean$seurat_clusters)

  cm_clean <- subset(
    master_clean,
    subset = cluster_chr_new %in% new_cm_clusters
  )

  cm_clean@reductions <- list()
  cm_clean@graphs <- list()
  cm_clean@neighbors <- list()

  DefaultAssay(cm_clean) <- "RNA"

  cm_clean <- NormalizeData(cm_clean, verbose = TRUE)
  cm_clean <- FindVariableFeatures(cm_clean, selection.method = "vst", nfeatures = 3000, verbose = TRUE)
  cm_clean <- ScaleData(cm_clean, features = VariableFeatures(cm_clean), verbose = TRUE)
  cm_clean <- RunPCA(cm_clean, features = VariableFeatures(cm_clean), npcs = 30, verbose = TRUE)

  cm_dims_use <- 1:20

  cm_clean <- FindNeighbors(cm_clean, dims = cm_dims_use, verbose = TRUE)
  cm_clean <- FindClusters(cm_clean, resolution = 0.3, verbose = TRUE)
  cm_clean <- RunUMAP(cm_clean, dims = cm_dims_use, verbose = TRUE)

  DimPlot(
    cm_clean,
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE
  )

  cm_markers_clean <- FindAllMarkers(
    cm_clean,
    only.pos = TRUE,
    min.pct = 0.1,
    logfc.threshold = 0.25
  )

  write.csv(
    cm_markers_clean,
    file = "CM_only_markers_cluster19_removed_reclustered.csv",
    row.names = FALSE
  )

  saveRDS(
    cm_clean,
    file = "cm_clean_cluster19_removed_reclustered.rds"
  )
}
