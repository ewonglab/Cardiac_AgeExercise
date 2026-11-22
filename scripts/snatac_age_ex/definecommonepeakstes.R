# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Define common ATAC peak sets and prepare downstream differential analyses.
# Steps: 1) Load per-sample/cell-type peaks; 2) merge and annotate common regions; 3) save peak sets and summaries.
# =======================================
# Source reformatted from snATAC_AgeEx/DefineCommonePeakStes.ipynb.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

# Redo analysis
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
  #legend.title = element_blank(),
  legend.text = element_text(size = 12),
  legend.position = "right",
  
  # Specify the panel/grid appearance, no grid line 
  panel.background = element_blank(),#element_rect(fill = "white", color = "black"),
  panel.grid.major = element_blank(),#element_line(color = "lightgray", size = 0.2),
  panel.grid.minor = element_blank(),#element_line(color = "lightgray", size = 0.2),
  
  # Specify the plot margins
  plot.margin = margin(0, 0, 0, 0, "cm"),
  
  #set the ratio of the plot 
  #aspect.ratio=1
)



setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2024-08-02.MACS_CallFinal/"))
peaks = readRDS("All_peaks.rds")

Macrophage = readRDS("M1M2_peaks.rds")







# Exploratory check of the fibroblast consensus peaks. Fibroblast_consensus.bed is written by the
# consensus-peak export further down ("5. Final Aggregation and Export"), so it does not exist on a
# first run. The check is therefore only performed when the file is already present.
refined_peaks_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_RefinedPeaks")
fib <- NULL
if (file.exists(file.path(refined_peaks_dir, "Fibroblast_consensus.bed"))) {
  setwd(refined_peaks_dir)
  list.files(pattern = ".bed")
  fib = read.table('Fibroblast_consensus.bed')
}

base_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Output/")
bindetect_dir <- file.path(base_dir, 'BINDetect')

# Output Directory
out_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Footprints/")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# TF List (AP-1 Family)
# Ensure these match the folder names inside BINDetect exactly
AP1_TFs <- c(
  "Fos_Homo_sapiens_DBD_1.00_M08069_2.00"
 
)
comp_dirs <- list.dirs(bindetect_dir, recursive = FALSE)
cell_types <- sub("_comp_all$", "", basename(comp_dirs[grep("_comp_all$", comp_dirs)]))

message("Found Cell Types: ", paste(cell_types, collapse = ", "))
comp_dirs

# CHeck the bed 

Fos_bind = read.table(header = TRUE,publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Output/BINDetect/Fibroblast_comp_all/Fos_Homo_sapiens_DBD_1.00_M08069_2.00/Fos_Homo_sapiens_DBD_1.00_M08069_2.00_overview.txt"))

Fos_bind$start = as.character(Fos_bind$TFBS_start)

DT::datatable(Fos_bind, filter = 'top')

if (!is.null(fib)) fib$start = as.character(fib$V2)
76693404
76693615
76693736

if (!is.null(fib)) DT::datatable(fib, filter = 'top')

76693404



library(GenomicRanges)
library(tidyverse)
library(stringr)
library(rtracklayer)

# 1. Load Data
setwd(publication_path("HEART_SNRNA_ROOT", "analysis/2024-08-02.MACS_CallFinal/"))
peaks_main <- readRDS("All_peaks.rds")
peaks_macro <- readRDS("M1M2_peaks.rds")

# Combine lists
peaks_raw <- c(peaks_main, peaks_macro)

# --- THE FIX: Clean the list ---
# Keep only elements that are actually GRanges objects. 
# This removes NULLs or corrupted entries.
peaks <- peaks_raw[sapply(peaks_raw, function(x) is(x, "GRanges"))]

message(paste("Original list length:", length(peaks_raw)))
message(paste("Cleaned list length:", length(peaks)))


# 2. Build Metadata (Now on the cleaned list)
meta_df <- tibble(
  list_index = seq_along(peaks),
  full_ident = map_chr(peaks, ~ as.character(unique(.x$ident))[1])
) %>%
  mutate(
    # Extract Cell Type
    cell_type = str_extract(full_ident, "^(.*?)_(?=MAE|MAS|YE|YS)"),
    cell_type = str_remove(cell_type, "_$"),
    
    # Extract Condition (MAE, MAS, YE, YS)
    condition = str_extract(full_ident, "(MAE|MAS|YE|YS)"),
    
    # Create Group ID
    group_id = paste(cell_type, condition, sep = "_")
  )

# Remove any rows where parsing failed (e.g. if an object didn't have the right name format)
meta_df <- meta_df %>% filter(!is.na(cell_type) & !is.na(condition))

# QC: Check counts
print(table(meta_df$group_id))


# 3. Consensus Function
get_consensus_peaks <- function(subset_indices, all_peaks_list, min_reps = 3) {
  
  # Safe subset
  gr_subset <- as(all_peaks_list[subset_indices], "GRangesList")
  all_ranges_flat <- unlist(gr_subset, use.names = FALSE)
  
  # Reduce to find universe
  union_peaks <- GenomicRanges::reduce(all_ranges_flat)
  
  # Count overlaps
  ov <- findOverlaps(union_peaks, all_ranges_flat)
  n_peaks_per_sample <- elementNROWS(gr_subset)
  sample_id_vector <- rep(seq_along(gr_subset), times = n_peaks_per_sample)
  
  counts <- tibble(
    peak_idx = queryHits(ov),
    sample_id = sample_id_vector[subjectHits(ov)]
  ) %>%
    distinct() %>% 
    group_by(peak_idx) %>%
    summarise(n_reps = n())
  
  # Filter
  keep_idx <- counts %>% 
    dplyr::filter(n_reps >= min_reps) %>% 
    pull(peak_idx)
  
  return(union_peaks[keep_idx])
}


# 4. Processing Loop (Consensus per Condition)
condition_consensus_list <- lapply(unique(meta_df$group_id), function(gid) {
  
  indices <- meta_df %>% 
    dplyr::filter(group_id == gid) %>% 
    pull(list_index)
  
  if(length(indices) < 3) {
    message(paste("Skipping", gid, "- fewer than 3 replicates."))
    return(NULL)
  }
  
  res <- get_consensus_peaks(indices, peaks, min_reps = 3)
  return(res)
})
names(condition_consensus_list) <- unique(meta_df$group_id)
condition_consensus_list <- condition_consensus_list[!sapply(condition_consensus_list, is.null)]


# 5. Final Aggregation and Export
out_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_RefinedPeaks")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

unique_cell_types <- unique(meta_df$cell_type)

for (ctype in unique_cell_types) {
  
  message(paste("Finalizing Cell Type:", ctype))
  
  # Identify relevant groups for this cell type
  relevant_groups <- grep(paste0("^", ctype, "_"), names(condition_consensus_list), value = TRUE)
  
  if (length(relevant_groups) == 0) {
    message("No valid consensus peaks found for ", ctype)
    next
  }
  
  # Merge the condition-specific consensus peaks
  combined_gr_list <- GRangesList(condition_consensus_list[relevant_groups])
  combined_gr_flat <- unlist(combined_gr_list, use.names = FALSE)
  
  # Final Reduce: Combine valid peaks from different conditions
  final_peaks <- GenomicRanges::reduce(combined_gr_flat)
  
  # Export
  safe_name <- gsub(" ", "_", ctype)
  filename <- file.path(out_dir, paste0(safe_name, "_consensus.bed"))
  export.bed(final_peaks, con = filename)
}

message("Processing Complete.")

# Load necessary libraries
library(tools)

# --- Configuration ---

# Input Directories (From previous steps)
base_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Output/")
bindetect_dir <- file.path(base_dir, 'BINDetect')

# Output Directory
out_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Footprints/")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# TF List (AP-1 Family)
# Ensure these match the folder names inside BINDetect exactly
AP1_TFs <- c(
  "Fos_Homo_sapiens_DBD_1.00_M08069_2.00",
  "Jun_M01816_2.00", 
  "Junb_M00760_2.00",
  "Jund_M01822_2.00", 
  "Fosl1_M01805_2.00",
  "Fosb_Homo_sapiens_DBD_1.00_M04272_2.00", 
  "Bach2_Homo_sapiens_DBD_0.95_M08065_2.00",
  "Nfe2l1_Homo_sapiens_DBD_1.00_M04008_2.00", 
  "Nfe2l2_M08077_2.00",
  "Nfe2_Homo_sapiens_DBD_0.95_M07813_2.00"
)

# Detect Cell Types from the BINDetect folder structure
# Looks for folders ending in "_comp_all"
comp_dirs <- list.dirs(bindetect_dir, recursive = FALSE)
cell_types <- sub("_comp_all$", "", basename(comp_dirs[grep("_comp_all$", comp_dirs)]))

message("Found Cell Types: ", paste(cell_types, collapse = ", "))

# --- Main Function ---

makeFootprints <- function(TFs, file_suffix, cell_types) {
  
  conditions <- c('MAS', 'MAE', 'YS', 'YE')
  
  for (cell_type in cell_types) {
    
    message(paste0("\n>>> Processing: ", cell_type, " <<<"))
    
    # Define paths
    ct_bindetect_dir <- file.path(bindetect_dir, paste0(cell_type, '_comp_all'))
    ct_out_dir <- file.path(out_dir, cell_type)
    if (!dir.exists(ct_out_dir)) dir.create(ct_out_dir)
    
    # 1. Aggregate Bound Sites per Condition (Optional but useful)
    # This creates a combined BED of all AP1 sites bound in MAS, then MAE, etc.
    for (con in conditions) {
      outfile <- file.path(ct_out_dir, paste0(file_suffix, "_", con, "_bound.bed"))
      
      # Initialize empty file
      system(paste(">", outfile))
      
      # Cat relevant TF beds
      for (TF in TFs) {
        tf_bed <- file.path(ct_bindetect_dir, TF, "beds", paste0(TF, "_", con, "_bound.bed"))
        
        if (file.exists(tf_bed)) {
          system(paste("cat", tf_bed, ">>", outfile))
        } else {
          warning(paste("Missing BED for", TF, con))
        }
      }
      
      # Sort
      if (file.info(outfile)$size > 0) {
        system(paste("sort -k1,1 -k2,2n", outfile, "-o", outfile))
      }
    }
    
    # 2. Create Master Region File for Plotting
    # Merges all sites from the list TFs into one file
    master_bed <- file.path(ct_out_dir, paste0(file_suffix, "_all_merged.bed"))
    system(paste(">", master_bed))
    
    for (TF in TFs) {
      # Use the 'all.bed' which contains all potential sites (bound and unbound)
      # Or use the bound union if you only want bound sites. 
      # Standard PlotAggregate usage often takes the union of bound sites across conditions.
      # Here we use 'all.bed' as per your original script logic (merging TF_all.bed)
      
      tf_all_bed <- file.path(ct_bindetect_dir, TF, "beds", paste0(TF, "_all.bed"))
      
      if (file.exists(tf_all_bed)) {
        system(paste("cat", tf_all_bed, ">>", master_bed))
      }
    }
    
    # Sort and Merge
    if (file.info(master_bed)$size > 0) {
      system(paste("sort -k1,1 -k2,2n", master_bed, "-o", master_bed))
      system(paste("bedtools merge -i", master_bed, ">", paste0(master_bed, ".tmp")))
      file.rename(paste0(master_bed, ".tmp"), master_bed)
    } else {
      warning(paste("No regions found for", cell_type))
      next
    }
    
    # 3. Define Signals (BigWigs)
    # We look for the ATACorrect output (corrected.bw)
    get_signal <- function(cond) {
      # Expected Path: .../Output/Coronary_EC_MAS/Coronary_EC_MAS_corrected.bw
      path <- file.path(base_dir, paste0(cell_type, "_", cond), paste0(cell_type, "_", cond, "_corrected.bw"))
      
      # Fallback: check for _merged_corrected.bw if using old naming
      if (!file.exists(path)) {
        path <- file.path(base_dir, paste0(cell_type, "_", cond), paste0(cell_type, "_", cond, "_merged_corrected.bw"))
      }
      return(path)
    }
    
    signals <- c(get_signal('MAS'), get_signal('MAE'), get_signal('YE'), get_signal('YS'))
    
    # Check if signals exist
    if (!all(file.exists(signals))) {
      warning("Missing signal files for ", cell_type, ". Skipping plot.")
      print(signals[!file.exists(signals)])
      next
    }
    
    # 4. Run TOBIAS PlotAggregate
    pdf_out <- file.path(ct_out_dir, paste0(file_suffix, "_footprint_comparisons.pdf"))
    txt_out <- file.path(ct_out_dir, paste0(file_suffix, "_footprint_comparisons.txt"))
    
    cmd <- paste(
      "TOBIAS PlotAggregate",
      "--TFBS", master_bed,
      "--signals", paste(signals, collapse = " "),
      "--output", pdf_out,
      "--output-txt", txt_out,
      "--share_y both",
      "--plot_boundaries",
      "--signal-on-x"
    )
    
    message("Running PlotAggregate...")
    system(cmd)
  }
}

# --- Execution ---

# Check if bedtools is installed/loaded
if (Sys.which("bedtools") == "") {
  stop("Error: 'bedtools' not found. Please run 'module load bedtools' before running R.")
}

# Run the function
makeFootprints(AP1_TFs, "AP1_Aggregated", cell_types)



library(ggplot2)
library(dplyr)
library(stringr)

# --- Configuration ---
# Update to the directory where we just generated the footprints
base_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Footprints")
setwd(base_dir)

# Define Colors
Group_Pallate <- c('#E41A1C', '#A020F0', '#00008B', '#377EB8')
Group_order <- c("MAS", "MAE", "YS", "YE")
group_colours <- setNames(Group_Pallate, Group_order)

# --- Functions ---

loadPlotAggregate <- function(file) {
  # Read the TOBIAS output text file
  Tef_comparison <- read.table(file, col.names = c('bigWig', 'bedfile', 'average.signal'))
  
  # Parse the "average.signal" column (comma-separated string -> numeric vector)
  # Each row becomes 120 rows (positions -59 to 60)
  TOBIAS_footprint_df <- data.frame(
    bedfile = rep(Tef_comparison$bedfile, each = 120),
    original_path = rep(Tef_comparison$bigWig, each = 120),
    average.signal = as.numeric(do.call('c', strsplit(Tef_comparison$average.signal, ','))),
    position = rep(c(-59:60), nrow(Tef_comparison))
  )
  
  # Extract Metadata from the file path/name
  # Logic: " .../Coronary_EC_MAS_corrected.bw" or similar
  
  # 1. Extract Condition (MAE, MAS, YE, YS)
  # We look for one of the specific strings in the path
  TOBIAS_footprint_df$condition <- str_extract(basename(TOBIAS_footprint_df$original_path), "(MAE|MAS|YE|YS)")
  
  # 2. Extract Cell Type
  # Assumes structure: CellType_Condition...
  # We grab everything before the condition
  TOBIAS_footprint_df$cell_type <- str_extract(basename(TOBIAS_footprint_df$original_path), "^(.*?)_(?=MAE|MAS|YE|YS)")
  TOBIAS_footprint_df$cell_type <- str_remove(TOBIAS_footprint_df$cell_type, "_$") # Remove trailing underscore
  
  return(TOBIAS_footprint_df)
}

FootprintPlot <- function(df, colour = 'condition') {
  ggplot(df, aes(x = position, y = average.signal, colour = .data[[colour]])) +
    scale_colour_manual(values = group_colours) +
    theme_classic() + 
    # Use the bedfile name as the title if multiple BEDs exist, otherwise facet by cell type later
    # facet_wrap(vars(bedfile)) + 
    theme(
      legend.direction = 'horizontal', 
      legend.position = 'bottom',
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    # Heavy smoothing as requested (polynomial degree 20)
    geom_smooth(linewidth = 1, se = FALSE, method = 'lm', formula = y ~ poly(x, 20))
}

# --- Execution ---

# 1. Find all comparison text files (Recursive search to look inside cell type folders)
# Pattern matches the suffix we defined in the previous step
footprint_files <- list.files(pattern = 'AP1_Aggregated_footprint_comparisons.txt$', recursive = TRUE, full.names = TRUE)

message(paste("Found", length(footprint_files), "footprint files."))

# 2. Load and Combine Data
if (length(footprint_files) > 0) {
  load_all_AP1 <- do.call('rbind', lapply(footprint_files, loadPlotAggregate))
  
  # 3. Filter Data
  # Remove YE condition
  load_all_AP1 <- load_all_AP1 %>% filter(condition != "YE")
  
  # Remove specific cell types
  load_all_AP1 <- load_all_AP1 %>% filter(!cell_type %in% c("M1", "M2", "LSC"))
  
  # Ensure Factor Order for Colors
  load_all_AP1$condition <- factor(load_all_AP1$condition, levels = Group_order)
  
  # 4. Plot
  p <- FootprintPlot(load_all_AP1, colour = 'condition') + 
       facet_wrap(vars(cell_type), scales = 'free_y', nrow = 2) +
       labs(y = "Tn5 Cut Signal (Normalized)", x = "Distance from Center (bp)")
  
  print(p)
  
  # 5. Save
  # Ensure output directory exists
  out_file <- publication_path("LEGACY_ATAC_RESULTS_ROOT", "FOS_lmsmoothed_allCellTypes_allsites.pdf")
  dir_create_path <- dirname(out_file)
  if (!dir.exists(dir_create_path)) dir.create(dir_create_path, recursive = TRUE)
  
  ggsave(out_file, plot = p, width = 297, height = 210, units = 'mm') # A4 Landscape roughly
  message(paste("Saved plot to:", out_file))
  
} else {
  warning("No footprint files found. Check directory and previous TOBIAS run.")
}

library(dplyr)
library(stringr)
library(readr)

# --- Configuration ---
# Set this to the parent BINDetect folder containing all your cell type subfolders
base_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Output/BINDetect")

# --- 1. Find all result files ---
# Recursive search for 'bindetect_results.txt' inside any subfolder
result_files <- list.files(base_dir, pattern = "bindetect_results.txt", recursive = TRUE, full.names = TRUE)

message(paste("Found", length(result_files), "result files."))

# --- 2. Function to Process Single File ---
process_tobias_file <- function(filepath) {
  
  # Read the text file (TOBIAS outputs are tab-separated)
  df <- read.table(filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Extract Cell Type from the parent folder name
  # Path looks like: .../BINDetect/Coronary_EC_comp_all/bindetect_results.txt
  parent_folder <- basename(dirname(filepath))
  cell_type_name <- str_remove(parent_folder, "_comp_all$")
  
  # Clean the Motif Name (Remove the ID in parentheses if present)
  # e.g., "Fos(MA0476.1)" -> "Fos"
    df$name <- df$name %>%
  str_remove("^\\(") %>%        # Step 1: Remove leading (
  str_remove("(_|\\().*") %>%   # Step 2: Cut off at metadata
  str_remove("\\)$")    
  # Add Metadata columns
  df$cell_type <- cell_type_name
  df$method <- "TOBIAS"
  df$ct_expressed <- NA # Placeholder if you need this column later
  df$X <- 1:nrow(df)    # Row index
  
  return(df)
}

# --- 3. Read and Combine All ---
# lapply runs the function on every file, bind_rows puts them into one big table
TOBIAS_combined <- bind_rows(lapply(result_files, process_tobias_file))

# --- 4. Final Post-Processing ---
TOBIAS_combined <- TOBIAS_combined %>%
  # Flip the sign for MAS_MAE change as requested
  mutate(MAS_MAE_change = -MAS_MAE_change) %>%
  # Filter out specific cell types
  filter(!cell_type %in% c("M1", "M2", "LSC"))

# Check the result
print(table(TOBIAS_combined$cell_type))
head(TOBIAS_combined)

# Optional: Save the master table
base_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Output/BINDetect")
write.csv(TOBIAS_combined, file.path(base_dir, "bindetect_results_all_merged.csv"), row.names = FALSE)

base_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Output/BINDetect")

# --- 1. Find all result files ---
# Recursive search for 'bindetect_results.txt' inside any subfolder
result_files <- list.files(base_dir, pattern = "bindetect_results.txt", recursive = TRUE, full.names = TRUE)

result_files

library(dplyr)
library(tidyr)
library(stringr)

# --- Configuration ---
base_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Output/BINDetect")
out_dir  <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Classification")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Get list of cell types automatically
comp_dirs <- list.dirs(base_dir, recursive = FALSE)
cells <- sub("_comp_all$", "", basename(comp_dirs[grep("_comp_all$", comp_dirs)]))

# --- Helper Function: Import Motif AND Peak Info ---
import_bed_safe <- function(path) {
  if (file.exists(path) && file.info(path)$size > 0) {
    tryCatch({
      # Read the file (No header, tab separated)
      # We rely on TOBIAS output structure:
      # Cols 1-3: Motif Coords
      # Cols 7-9: Peak Coords (from the user-provided peak file)
      df <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      
      # Ensure file has enough columns (at least 9 for Peak info)
      if (ncol(df) < 9) {
        warning(paste("File", path, "does not have peak columns (cols 7-9). Skipping."))
        return(NULL)
      }
      
      # Extract & Rename
      # Motif Info
      res <- df[, 1:3]
      colnames(res) <- c("seqnames", "start", "end")
      
      # Peak Info (Parent Peak)
      res$peak_chr   <- df[, 7]
      res$peak_start <- df[, 8]
      res$peak_end   <- df[, 9]
      
      # Create Identifiers
      # 1. Motif ID (Small footprint, e.g. 10bp) -> Used for "Is this motif bound in condition X?"
      res$motif_identifier <- paste(res$seqnames, res$start, res$end, sep = "-")
      
      # 2. Peak ID (Large consensus peak, e.g. 500bp) -> Used for "Who else is in this peak?"
      res$peak_identifier  <- paste(res$peak_chr, res$peak_start, res$peak_end, sep = "-")
      
      return(res)
      
    }, error = function(e) {
      warning(paste("Failed to read:", path, "-", e$message))
      return(NULL)
    })
  } else {
    return(NULL)
  }
}

# --- Main Loop ---
# Assuming TOBIAS_combined is loaded, or get unique motifs from folder names
# Motifs <- sort(unique(TOBIAS_combined$output_prefix)) 
# Fallback if TOBIAS_combined isn't loaded: get from one cell folder
example_cell <- cells[1]
example_path <- file.path(base_dir, paste0(example_cell, "_comp_all"))
Motifs <- list.dirs(example_path, full.names = FALSE, recursive = FALSE)
Motifs <- Motifs[Motifs != ""] # Remove empty

for (Motif in Motifs) {
  
  # Check output existence
  out_file <- file.path(out_dir, paste0(Motif, "_containingPeaks_Classified.rds"))
  if (file.exists(out_file)) {
     #message(paste("Skipping", Motif, "- Output exists.")) 
    next 
  }
  
  message(paste("Processing Motif:", Motif))
  save_list <- list()
  
  for (cell in cells) {
    
    # Define paths
    motif_dir <- file.path(base_dir, paste0(cell, "_comp_all"), Motif, "beds")
    f_mae <- file.path(motif_dir, paste0(Motif, "_MAE_bound.bed"))
    f_mas <- file.path(motif_dir, paste0(Motif, "_MAS_bound.bed"))
    f_ys  <- file.path(motif_dir, paste0(Motif, "_YS_bound.bed"))
    
    # Import Data (Now includes Peak Info)
    df_mae <- import_bed_safe(f_mae)
    df_mas <- import_bed_safe(f_mas)
    df_ys  <- import_bed_safe(f_ys)
    
    if (is.null(df_mae) && is.null(df_mas) && is.null(df_ys)) next
    
    # Combine Unique Motif Footprints
    # We define "uniqueness" by the MOTIF footprint, not the peak
    all_peaks <- bind_rows(df_mae, df_mas, df_ys) %>%
      distinct(motif_identifier, .keep_all = TRUE)
    
    # Annotate Conditions (Based on MOTIF presence)
    if (!is.null(df_mae)) all_peaks$in_MAE <- all_peaks$motif_identifier %in% df_mae$motif_identifier else all_peaks$in_MAE <- FALSE
    if (!is.null(df_mas)) all_peaks$in_MAS <- all_peaks$motif_identifier %in% df_mas$motif_identifier else all_peaks$in_MAS <- FALSE
    if (!is.null(df_ys))  all_peaks$in_YS  <- all_peaks$motif_identifier %in% df_ys$motif_identifier  else all_peaks$in_YS  <- FALSE
    
    # Classify
    all_peaks <- all_peaks %>%
      mutate(
        category = case_when(
          in_MAE & in_MAS & in_YS ~ "In All Conditions",
          in_MAE & in_YS & !in_MAS ~ "Lost AgeSed",      
          in_MAE & in_MAS & !in_YS ~ "Lost YoungSed",    
          in_MAS & in_YS & !in_MAE ~ "Lost AgeEx",       
          in_YS & !in_MAE & !in_MAS ~ "In Young Sedentary Only",
          in_MAE & !in_YS & !in_MAS ~ "In Aged Ex Only",
          in_MAS & !in_YS & !in_MAE ~ "In Aged Sed Only",
          TRUE ~ "Complex/Other"
        ),
        cell_type = cell,
        motif_id = Motif
      )
    
    save_list[[cell]] <- all_peaks
  }
  
  if (length(save_list) > 0) {
    final_df <- bind_rows(save_list)
    saveRDS(final_df, out_file)
    message(paste("Saved:", out_file))
  } else {
    message(paste("No binding sites found for", Motif))
  }
}

message("Processing Complete.")

library(dplyr)
library(tidyr)
library(stringr)

# --- Configuration ---
base_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Output/BINDetect")
out_dir  <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Classification")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Get list of cell types automatically
comp_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
cells <- sub("_comp_all$", "", basename(comp_dirs[grep("_comp_all$", comp_dirs)]))

# --- Helpers ---
import_bed_safe <- function(path, score_name = NULL) {
  if (!file.exists(path) || file.info(path)$size == 0) return(NULL)

  tryCatch({
    df <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")

    # Need motif coords + peak coords
    if (ncol(df) < 9) {
      warning(paste("File", path, "does not have peak columns (cols 7-9). Skipping."))
      return(NULL)
    }

    res <- df[, 1:3, drop = FALSE]
    colnames(res) <- c("seqnames", "start", "end")

    res$peak_chr   <- df[, 7]
    res$peak_start <- df[, 8]
    res$peak_end   <- df[, 9]

    res$motif_identifier <- paste(res$seqnames, res$start, res$end, sep = "-")
    res$peak_identifier  <- paste(res$peak_chr, res$peak_start, res$peak_end, sep = "-")

    # TOBIAS bed usually has "score" in col 5 (bed6: chrom start end name score strand)
    if (!is.null(score_name) && ncol(df) >= 5) {
      res[[score_name]] <- suppressWarnings(as.numeric(df[, 5]))
    } else if (!is.null(score_name)) {
      res[[score_name]] <- NA_real_
    }

    res
  }, error = function(e) {
    warning(paste("Failed to read:", path, "-", e$message))
    NULL
  })
}

import_bindetect_scores_safe <- function(motif_root) {
  # Look for a file in the motif folder that contains MAS_score/MAE_score/YS_score columns
  if (!dir.exists(motif_root)) return(NULL)

  all_paths <- list.files(motif_root, full.names = TRUE, recursive = FALSE)
  if (length(all_paths) == 0) return(NULL)

  files <- all_paths[!file.info(all_paths)$isdir]
  if (length(files) == 0) return(NULL)

  for (f in files) {
    hdr <- tryCatch(readLines(f, n = 1), error = function(e) "")
    if (!nzchar(hdr)) next

    if (grepl("MAS_score", hdr, fixed = TRUE) &&
        grepl("MAE_score", hdr, fixed = TRUE) &&
        grepl("YS_score",  hdr, fixed = TRUE)) {

      df <- tryCatch(
        read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", fill = TRUE),
        error = function(e) NULL
      )
      if (is.null(df) || nrow(df) == 0) next

      # Normalise names
      colnames(df) <- gsub("\\.", "_", colnames(df))

      # Identify coordinate columns
      nms_l <- tolower(colnames(df))
      pick_col <- function(candidates) {
        hit <- which(nms_l %in% candidates)
        if (length(hit) > 0) return(colnames(df)[hit[1]])
        NA_character_
      }

      chr_col <- pick_col(c("seqnames", "chrom", "chr"))
      st_col  <- pick_col(c("start"))
      en_col  <- pick_col(c("end"))

      # If not found, assume first 3 columns are coords
      if (is.na(chr_col) || is.na(st_col) || is.na(en_col)) {
        chr_col <- colnames(df)[1]
        st_col  <- colnames(df)[2]
        en_col  <- colnames(df)[3]
      }

      if (!all(c("MAS_score", "MAE_score", "YS_score") %in% colnames(df))) next

      res <- df %>%
        transmute(
          seqnames = .data[[chr_col]],
          start    = .data[[st_col]],
          end      = .data[[en_col]],
          MAS_score = suppressWarnings(as.numeric(MAS_score)),
          MAE_score = suppressWarnings(as.numeric(MAE_score)),
          YS_score  = suppressWarnings(as.numeric(YS_score))
        ) %>%
        mutate(
          motif_identifier = paste(seqnames, start, end, sep = "-")
        ) %>%
        distinct(motif_identifier, .keep_all = TRUE)

      return(res)
    }
  }

  NULL
}

# --- Motif discovery (from one example cell) ---
example_cell <- cells[1]
example_path <- file.path(base_dir, paste0(example_cell, "_comp_all"))
Motifs <- list.dirs(example_path, full.names = FALSE, recursive = FALSE)
Motifs <- Motifs[Motifs != ""]

# --- Main Loop ---
for (Motif in Motifs) {

  out_file <- file.path(out_dir, paste0(Motif, "_containingPeaks_Classified_binding.rds"))
  if (file.exists(out_file)) next

  message(paste("Processing Motif:", Motif))
  save_list <- list()

  for (cell in cells) {

    motif_root <- file.path(base_dir, paste0(cell, "_comp_all"), Motif)

    # BEDs (bound sites per condition)
    motif_beds_dir <- file.path(motif_root, "beds")
    f_mae <- file.path(motif_beds_dir, paste0(Motif, "_MAE_bound.bed"))
    f_mas <- file.path(motif_beds_dir, paste0(Motif, "_MAS_bound.bed"))
    f_ys  <- file.path(motif_beds_dir, paste0(Motif, "_YS_bound.bed"))

    df_mae <- import_bed_safe(f_mae, score_name = "MAE_score")
    df_mas <- import_bed_safe(f_mas, score_name = "MAS_score")
    df_ys  <- import_bed_safe(f_ys,  score_name = "YS_score")

    if (is.null(df_mae) && is.null(df_mas) && is.null(df_ys)) next

    # Union of motif footprints (unique by motif footprint coords)
    sites <- bind_rows(
      df_mae %>% select(seqnames, start, end, peak_chr, peak_start, peak_end, motif_identifier, peak_identifier),
      df_mas %>% select(seqnames, start, end, peak_chr, peak_start, peak_end, motif_identifier, peak_identifier),
      df_ys  %>% select(seqnames, start, end, peak_chr, peak_start, peak_end, motif_identifier, peak_identifier)
    ) %>%
      distinct(motif_identifier, .keep_all = TRUE)

    # Presence flags (by motif footprint)
    sites$in_MAE <- if (!is.null(df_mae)) sites$motif_identifier %in% df_mae$motif_identifier else FALSE
    sites$in_MAS <- if (!is.null(df_mas)) sites$motif_identifier %in% df_mas$motif_identifier else FALSE
    sites$in_YS  <- if (!is.null(df_ys))  sites$motif_identifier %in% df_ys$motif_identifier  else FALSE

    # Scores: prefer a BINDetect per-site score file if present; otherwise use score from each bound BED
    score_df <- import_bindetect_scores_safe(motif_root)

    if (!is.null(score_df)) {
      sites <- sites %>%
        left_join(score_df %>% select(motif_identifier, MAS_score, MAE_score, YS_score), by = "motif_identifier")
    } else {
      mae_scores <- if (!is.null(df_mae)) df_mae %>% distinct(motif_identifier, .keep_all = TRUE) %>% select(motif_identifier, MAE_score) else NULL
      mas_scores <- if (!is.null(df_mas)) df_mas %>% distinct(motif_identifier, .keep_all = TRUE) %>% select(motif_identifier, MAS_score) else NULL
      ys_scores  <- if (!is.null(df_ys))  df_ys  %>% distinct(motif_identifier, .keep_all = TRUE) %>% select(motif_identifier, YS_score)  else NULL

      if (!is.null(mae_scores)) sites <- left_join(sites, mae_scores, by = "motif_identifier")
      if (!is.null(mas_scores)) sites <- left_join(sites, mas_scores, by = "motif_identifier")
      if (!is.null(ys_scores))  sites <- left_join(sites, ys_scores,  by = "motif_identifier")
    }

    # Classify
    sites <- sites %>%
      mutate(
        category = case_when(
          in_MAE & in_MAS & in_YS ~ "In All Conditions",
          in_MAE & in_YS & !in_MAS ~ "Lost AgeSed",
          in_MAE & in_MAS & !in_YS ~ "Lost YoungSed",
          in_MAS & in_YS & !in_MAE ~ "Lost AgeEx",
          in_YS & !in_MAE & !in_MAS ~ "In Young Sedentary Only",
          in_MAE & !in_YS & !in_MAS ~ "In Aged Ex Only",
          in_MAS & !in_YS & !in_MAE ~ "In Aged Sed Only",
          TRUE ~ "Complex/Other"
        ),
        cell_type = cell,
        motif_id = Motif
      )

    save_list[[cell]] <- sites
  }

  if (length(save_list) > 0) {
    final_df <- bind_rows(save_list)
    saveRDS(final_df, out_file)
    message(paste("Saved:", out_file))
  } else {
    message(paste("No binding sites found for", Motif))
  }
}

message("Processing Complete.")

# Define Paths
base_dir   <- publication_path("HEART_SNRNA_ROOT", "analysis")
tobias_dir <- file.path(base_dir, "2026-01-19.TOBIAS_Classification") # NEW PATH
mast_dir   <- file.path(base_dir, "2023-11-06.MAST_Final_Fixed")
cicero_dir <- file.path(base_dir, "2023-10-31.Cicero")
gtf_path   <- publication_path("CELLRANGER_GEX_REFERENCE", "genes/genes.gtf")

# Output for intermediate saves
out_dir <- file.path(cicero_dir, "2026-01-20_Recreated_Links")
if(!dir.exists(out_dir)) dir.create(out_dir)

# ==============================================================================
# 2. Load NEW TOBIAS Data (Classified)
# ==============================================================================
message("Loading new TOBIAS classification files...")

tobias_files <- list.files(tobias_dir, pattern = "_containingPeaks_Classified_binding.rds$", full.names = TRUE)

# Read and combine all TOBIAS files
Tobias_List <- lapply(tobias_files, function(f) {
  # Extract motif name from filename
  # Filename format: "MotifName_..._containingPeaks_Classified.rds"
  motif_name <- basename(f) %>% str_remove("_containingPeaks_Classified_binding.rds")
  
  # Read RDS (This is a single DF or a list of DFs depending on previous steps)
  dat <- readRDS(f)
  
  # Ensure it is a data frame
  if (!is.data.frame(dat)) {
    # If it's a list (split by cell), bind it back
    dat <- bind_rows(dat)
  }
  
  if (nrow(dat) > 0) {
    dat$motif <- motif_name
    return(dat)
  } else {
    return(NULL)
  }
})

Tobias_df <- bind_rows(Tobias_List)

# Clean Cell Types
Tobias_df$cell_type <- ifelse(Tobias_df$cell_type == "M1M2", "AllMacrophage", Tobias_df$cell_type)
Tobias_df$cell_type <- ifelse(Tobias_df$cell_type == "MSC", "HSPC", Tobias_df$cell_type)
# Rename to match downstream logic
Tobias_df <- Tobias_df %>% dplyr::rename(cell = cell_type, peak_coord_Age = peak_identifier)

message(paste("Loaded TOBIAS data. Total rows:", nrow(Tobias_df)))
out_dir <- file.path(cicero_dir, "2026-01-20_Recreated_Links")
saveRDS(Tobias_df, file.path(out_dir, "Tobias_df_Classified_Full_scores.rds"))

head(out_dir)

library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(scales)

# --- 1. Load Data ---
# TOBIAS_combined <- read.csv(...) # Load if needed

# --- 2. Define TFs of Interest ---
Rev_TF <- c('Fos','Fosb','Fosl1','Fosl2',
            'Junb','Jund','Jun',
            'Atf2','Atf3','Atf4','Batf','Batf3',
            'Maf','Mafa','Mafb','Mafg')

# --- 3. Filter Data ---
plot_df <- TOBIAS_combined %>%
  filter(name %in% Rev_TF) %>%
  filter(!cell_type %in% c("M1", "M2", "LSC")) 

# --- 4. Cluster Rows ---
matrix_df <- plot_df %>%
  select(cell_type, name, MAS_YS_change) %>%
  pivot_wider(names_from = cell_type, values_from = MAS_YS_change, values_fill = 0) %>%
  as.data.frame()

rownames(matrix_df) <- matrix_df$name
matrix_df$name <- NULL

ht <- Heatmap(as.matrix(matrix_df), cluster_columns = FALSE, cluster_rows = TRUE)
ht_draw <- draw(ht)
ordered_tfs <- rownames(matrix_df)[row_order(ht_draw)]

plot_df$name <- factor(plot_df$name, levels = ordered_tfs)

# --- 5. Define Theme (Rotated X-Axis 90 degrees) ---
theme_dotplot <- theme_bw() + theme(
  legend.title = element_blank(),
  # Rotate 90 degrees, align vertical center to tick
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
  axis.text.y = element_text(color = "black"),
  panel.grid.major = element_line(color = "grey90"),
  panel.grid.minor = element_blank()
)

# --- 6. Plot 1: Aging Effect ---
Age_FC_AP1 <- ggplot(plot_df, aes(x = cell_type, y = name, 
                                  colour = MAS_YS_change, 
                                  size = -log10(MAS_YS_pvalue))) +
  geom_point() +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", 
    midpoint = 0, 
    limits = c(-0.3, 0.3), 
    name = "Differential Binding"
  ) +
  scale_size_continuous(limits = c(0, 170), name = "-log10(P)") +
  labs(title = "Aging Effect (Old vs Young)", x = "", y = "") +
  theme_dotplot

# --- 7. Plot 2: Exercise Effect ---
Ex_FC_AP1 <- ggplot(plot_df, aes(x = cell_type, y = name, 
                                 colour = MAS_MAE_change, 
                                 size = -log10(MAS_MAE_pvalue))) +
  geom_point() +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", 
    midpoint = 0, 
    limits = c(-0.3, 0.3), 
    name = "Differential Binding"
  ) +
  scale_size_continuous(limits = c(0, 170), name = "-log10(P)") +
  labs(title = "Exercise Effect (Ex vs Sed)", x = "", y = "") +
  theme_dotplot

# --- 8. Output ---
print(Age_FC_AP1)
print(Ex_FC_AP1)

# Save
out_path <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Footprints/")
ggsave(file.path(out_path, "DotPlot_AP1_Aging.pdf"), Age_FC_AP1, width = 7, height = 6)
ggsave(file.path(out_path, "DotPlot_AP1_Exercise.pdf"), Ex_FC_AP1, width = 7, height = 6)

base_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Output/BINDetect")
TOBIAS_combined = read.csv(file.path(base_dir, "bindetect_results_all_merged.csv"))
head(TOBIAS_combined)

# Ceck BACh1 


library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(scales)

# --- 1. Load Data ---
# TOBIAS_combined <- read.csv(...) # Load if needed

# --- 2. Define TFs of Interest ---
Rev_TF <- c('Fos','Fosb','Fosl1','Fosl2',
            'Junb','Jund','Jun',
            'Atf2','Atf3','Atf4','Batf','Batf3',
            'Maf','Mafa','Mafb','Mafg',"Bach1")

# --- 3. Filter Data ---
plot_df <- TOBIAS_combined %>%
  filter(name %in% Rev_TF) %>%
  filter(!cell_type %in% c("M1", "M2", "LSC")) 

# --- 4. Cluster Rows ---
matrix_df <- plot_df %>%
  select(cell_type, name, MAS_YS_change) %>%
  pivot_wider(names_from = cell_type, values_from = MAS_YS_change, values_fill = 0) %>%
  as.data.frame()

rownames(matrix_df) <- matrix_df$name
matrix_df$name <- NULL

ht <- Heatmap(as.matrix(matrix_df), cluster_columns = FALSE, cluster_rows = TRUE)
ht_draw <- draw(ht)
ordered_tfs <- rownames(matrix_df)[row_order(ht_draw)]

plot_df$name <- factor(plot_df$name, levels = ordered_tfs)

# --- 5. Define Theme (Rotated X-Axis 90 degrees) ---
theme_dotplot <- theme_bw() + theme(
  legend.title = element_blank(),
  # Rotate 90 degrees, align vertical center to tick
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
  axis.text.y = element_text(color = "black"),
  panel.grid.major = element_line(color = "grey90"),
  panel.grid.minor = element_blank()
)

# --- 6. Plot 1: Aging Effect ---
Age_FC_AP1 <- ggplot(plot_df, aes(x = cell_type, y = name, 
                                  colour = MAS_YS_change, 
                                  size = -log10(MAS_YS_pvalue))) +
  geom_point() +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", 
    midpoint = 0, 
    limits = c(-0.3, 0.3), 
    name = "Differential Binding"
  ) +
  scale_size_continuous(limits = c(0, 170), name = "-log10(P)") +
  labs(title = "Aging Effect (Old vs Young)", x = "", y = "") +
  theme_dotplot

# --- 7. Plot 2: Exercise Effect ---
Ex_FC_AP1 <- ggplot(plot_df, aes(x = cell_type, y = name, 
                                 colour = MAS_MAE_change, 
                                 size = -log10(MAS_MAE_pvalue))) +
  geom_point() +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", 
    midpoint = 0, 
    limits = c(-0.3, 0.3), 
    name = "Differential Binding"
  ) +
  scale_size_continuous(limits = c(0, 170), name = "-log10(P)") +
  labs(title = "Exercise Effect (Ex vs Sed)", x = "", y = "") +
  theme_dotplot

# --- 8. Output ---
print(Age_FC_AP1)
print(Ex_FC_AP1)

# Save
out_path <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Footprints/")
#ggsave(file.path(out_path, "DotPlot_AP1_Aging.pdf"), Age_FC_AP1, width = 7, height = 6)
#ggsave(file.path(out_path, "DotPlot_AP1_Exercise.pdf"), Ex_FC_AP1, width = 7, height = 6)



library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

# ==============================================================================
# 1. Setup & Filtering
# ==============================================================================

# Define AP-1 TFs
Rev_TF <- c('Fos','Fosb','Fosl1','Fosl2',
            'Junb','Jund','Jun',
            'Atf2','Atf3','Atf4','Batf','Batf3',
            'Maf','Mafa','Mafb','Mafg')

# 1. Define Candidate Pool (Stable Motifs Only)
# We look for motifs that are NOT AP-1 and have very low change scores
candidates_pool <- TOBIAS_combined %>%
  filter(!cell_type %in% c("M1", "M2", "LSC")) %>%
  filter(!name %in% Rev_TF) %>%
  group_by(name) %>%
  summarise(
    # Stability Score: Sum of absolute changes (lower = more stable)
    stability = mean(abs(MAS_YS_change) + abs(MAS_MAE_change), na.rm = TRUE),
    n_types = n_distinct(cell_type)
  ) %>%
  # Strict Stability Filter (Ensure they are true controls)
  filter(stability < 0.04) %>%
  pull(name)

cat("Found", length(candidates_pool), "stable candidate motifs.\n")

# ==============================================================================
# 2. Prepare Data Matrices for Matching
# ==============================================================================

# We need a matrix of Counts: Rows = Motifs, Cols = Cell Types
count_matrix <- TOBIAS_combined %>%
  filter(!cell_type %in% c("M1", "M2", "LSC")) %>%
  filter(name %in% c(Rev_TF, candidates_pool)) %>%
  select(name, cell_type, MAS_bound) %>%
  # Handle missing values (fill with 0 if motif not found in a cell type)
  pivot_wider(names_from = cell_type, values_from = MAS_bound, values_fill = 0) %>%
  column_to_rownames("name")

# Normalize counts (Log10) to prevent 'Super Binders' from dominating the distance metric
# (Differences in millions matter less than log-fold differences)
log_matrix <- log10(count_matrix + 1)

# Split into Target (AP-1) and Candidates
target_mat <- log_matrix[rownames(log_matrix) %in% Rev_TF, ]
candidate_mat <- log_matrix[rownames(log_matrix) %in% candidates_pool, ]

# ==============================================================================
# 3. The Matching Algorithm
# ==============================================================================
# For each AP-1, find the candidate with the smallest Euclidean distance

results_list <- list()
used_controls <- c() # Track used controls to prevent duplicates if desired

for (ap1 in rownames(target_mat)) {
  
  # Get the profile for this specific AP-1 motif
  target_vec <- target_mat[ap1, ]
  
  # Calculate distance to ALL candidates
  # Dist = Sum of squared differences across all cell types
  dists <- apply(candidate_mat, 1, function(x) sum((x - target_vec)^2))
  
  # Sort by closest match
  matches <- sort(dists)
  
  # Pick the best match that hasn't been used yet (Greedy matching)
  # If you don't mind duplicates, remove the '!names %in% used_controls' check
  available_matches <- matches[!names(matches) %in% used_controls]
  
  if(length(available_matches) > 0) {
    best_match <- names(available_matches)[1]
    score <- available_matches[1]
    
    # Mark as used
    used_controls <- c(used_controls, best_match)
    
    results_list[[ap1]] <- data.frame(
      AP1_Motif = ap1,
      Control_Motif = best_match,
      Distance_Score = score,
      stringsAsFactors = FALSE
    )
  } else {
    warning(paste("No unique control found for", ap1))
  }
}

# Compile Matching Table
matched_pairs <- do.call(rbind, results_list)
print(matched_pairs)

# Save the pair list
write.csv(matched_pairs, file.path(out_path, "AP1_Matched_Control_Pairs.csv"), row.names = FALSE)

# ==============================================================================
# 4. Visualization: Verify the Pairing
# ==============================================================================

# Prepare data for plotting
# We pull the counts for the specific pairs selected
plot_df <- TOBIAS_combined %>%
  filter(!cell_type %in% c("M1", "M2", "LSC")) %>%
  select(name, cell_type, MAS_bound)

# Join with our matching table
# This creates a "long" format with pairs side-by-side
pair_plot_data <- matched_pairs %>%
  inner_join(plot_df, by = c("AP1_Motif" = "name")) %>%
  rename(AP1_Count = MAS_bound) %>%
  inner_join(plot_df, by = c("Control_Motif" = "name", "cell_type")) %>%
  rename(Control_Count = MAS_bound)

# --- Plot A: Correlation Scatter ---
# Ideally, points should fall on the diagonal line (Perfect 1:1 match)
p_corr <- ggplot(pair_plot_data, aes(x = AP1_Count, y = Control_Count)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  geom_point(aes(color = cell_type), alpha = 0.6) +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(~AP1_Motif, scales = "free") + # One plot per AP-1 motif
  theme_bw() +
  labs(title = "Quality Control: AP-1 vs Selected Partner",
       subtitle = "Points should cluster along the diagonal line",
       x = "AP-1 Peak Count", y = "Control Peak Count")

print(p_corr)
ggsave(file.path(out_path, "QC_Matched_Pairs_Scatter.pdf"), p_corr, width = 12, height = 10)



library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

# --- 1. Prepare Data ---
# Combine AP-1 and the NEW Control TFs into one dataframe
# Assumes 'Control_TF' is available from your previous step
plot_data <- TOBIAS_combined %>%
  filter(name %in% c(matched_pairs$AP1_Motif, matched_pairs$Control_Motif)) %>%
  filter(!cell_type %in% c("M1", "M2", "LSC")) %>%
  mutate(Group = ifelse(name %in% Rev_TF, "AP-1", "Control"))

# ==============================================================================
# PLOT 1: Distribution of Prevalence (Peak Counts)
# ==============================================================================
# This answers: "Are the counts similar across cell types?"

p_dist <- ggplot(plot_data, aes(x = cell_type, y = MAS_bound, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  # Add points to show real data spread
  geom_point(position = position_jitterdodge(jitter.width = 0.2), 
             size = 0.5, alpha = 0.4, color = "black") +
  scale_fill_manual(values = c("AP-1" = "firebrick", "Control" = "grey70")) +
  scale_y_log10(labels = scales::comma) + 
  theme_bw() +
  labs(title = "Prevalence Distribution: AP-1 vs Matched Controls",
       subtitle = "Controls should now match the abundance of AP-1",
       y = "Number of Bound Peaks (Log10)",
       x = "Cell Type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

print(p_dist)


# ==============================================================================
# PLOT 2 & 3: Changes in Binding (Dot Plots for Controls)
# ==============================================================================
# This answers: "How do the controls change with Age/Exercise?" (Should be stable)

# Subset just the controls for the Dot Plots
ctrl_df <- plot_data %>% filter(Group == "Control")

# --- Cluster Rows for Visualization ---
# This orders the TFs so similar patterns sit next to each other
matrix_df <- ctrl_df %>%
  dplyr::select(cell_type, name, MAS_YS_change) %>%
  tidyr::pivot_wider(names_from = cell_type, values_from = MAS_YS_change, values_fill = 0) %>%
  as.data.frame()
rownames(matrix_df) <- matrix_df$name
matrix_df$name <- NULL

ht <- Heatmap(as.matrix(matrix_df), cluster_columns = FALSE, cluster_rows = TRUE)
ht_draw <- draw(ht)
ordered_tfs <- rownames(matrix_df)[row_order(ht_draw)]
ctrl_df$name <- factor(ctrl_df$name, levels = ordered_tfs)

# --- Define Theme ---
theme_dotplot <- theme_bw() + theme(
  legend.title = element_text(size=10),
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
  axis.text.y = element_text(color = "black"),
  panel.grid.major = element_line(color = "grey90"),
  panel.grid.minor = element_blank()
)

# --- Plot Aging Changes ---
p_age <- ggplot(ctrl_df, aes(x = cell_type, y = name, 
                             colour = MAS_YS_change, 
                             size = -log10(MAS_YS_pvalue))) +
  geom_point() +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", 
    midpoint = 0, 
    limits = c(-0.3, 0.3), # Fixed scale matching AP-1
    name = "Diff. Binding"
  ) +
  scale_size_continuous(limits = c(0, 170), name = "-log10(P)") +
  labs(title = "Control Stability: Aging (Old vs Young)", 
       subtitle = "Ideally mostly white (no change)", x = "", y = "") +
  theme_dotplot

# --- Plot Exercise Changes ---
p_ex <- ggplot(ctrl_df, aes(x = cell_type, y = name, 
                            colour = MAS_MAE_change, 
                            size = -log10(MAS_MAE_pvalue))) +
  geom_point() +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", 
    midpoint = 0, 
    limits = c(-0.3, 0.3), # Fixed scale matching AP-1
    name = "Diff. Binding"
  ) +
  scale_size_continuous(limits = c(0, 170), name = "-log10(P)") +
  labs(title = "Control Stability: Exercise (Ex vs Sed)", 
       subtitle = "Ideally mostly white (no change)", x = "", y = "") +
  theme_dotplot

print(p_age)
print(p_ex)

# ==============================================================================
# SAVE
# ==============================================================================
out_path <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Footprints/")
ggsave(file.path(out_path, "Boxplot_Prevalence_AP1_vs_Control.pdf"), p_dist, width = 8, height = 6)
ggsave(file.path(out_path, "DotPlot_Controls_Aging_Matched.pdf"), p_age, width = 7, height = 6)
ggsave(file.path(out_path, "DotPlot_Controls_Exercise_Matched.pdf"), p_ex, width = 7, height = 6)



library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

# ==============================================================================
# 1. Setup & Path Discovery
# ==============================================================================

base_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Output/BINDetect")
out_file <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Classification/AP1_Matched_Control_Peak_Coordinates.csv")

# 1. Define Short Names (Target List)
Rev_TF <- c('Fos','Fosb','Fosl1','Fosl2',
            'Junb','Jund','Jun',
            'Atf2','Atf3','Atf4','Batf','Batf3',
            'Maf','Mafa','Mafb','Mafg')

# Load Matched Pairs to get Control Short Names
pair_file <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Footprints/AP1_Matched_Control_Pairs.csv")
if (file.exists(pair_file)) {
  matched_pairs <- read.csv(pair_file)
  Control_TF <- unique(matched_pairs$Control_Motif)
} else {
  stop("Matched Pairs file not found. Please run the pairing step first.")
}

Target_Short_Names <- unique(c(Rev_TF, Control_TF))

# 2. Get Cell Types
comp_dirs <- list.dirs(base_dir, recursive = FALSE)
cells <- sub("_comp_all$", "", basename(comp_dirs[grep("_comp_all$", comp_dirs)]))
cells <- cells[!cells %in% c("M1", "M2", "LSC", "M1M2")]

# ==============================================================================
# 2. Build Name Map (Short Name -> Full Directory Name)
# ==============================================================================
example_cell_dir <- file.path(base_dir, paste0(cells[1], "_comp_all"))
all_motif_folders <- list.dirs(example_cell_dir, full.names = FALSE, recursive = FALSE)

get_full_name <- function(short, all_folders) {
  pattern <- paste0("^", short, "(_|$)")
  matches <- all_folders[grepl(pattern, all_folders)]
  if (length(matches) == 0) return(NA)
  if (length(matches) > 1) return(matches[which.min(nchar(matches))]) 
  return(matches)
}

name_map <- data.frame(Short_Name = Target_Short_Names, stringsAsFactors = FALSE)
name_map$Full_Name <- sapply(name_map$Short_Name, get_full_name, all_folders = all_motif_folders)
valid_map <- name_map %>% filter(!is.na(Full_Name))

cat("Mapped", nrow(valid_map), "motifs to full folder names.\n")

# ==============================================================================
# 3. Helper Function
# ==============================================================================
import_bed_safe <- function(path) {
  if (file.exists(path) && file.info(path)$size > 0) {
    tryCatch({
      df <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      if (ncol(df) < 9) return(NULL)
      
      res <- df[, 1:3]
      colnames(res) <- c("TFBS chr", "TFBS start", "TFBS end")
      res$`Peak chr`   <- df[, 7]
      res$`Peak start` <- df[, 8]
      res$`Peak End`   <- df[, 9]
      return(res)
    }, error = function(e) return(NULL))
  } else {
    return(NULL)
  }
}

# ==============================================================================
# 4. Main Loop
# ==============================================================================
results_list <- list()
counter <- 1

cat("Extracting coordinates...\n")

for (cell in cells) {
  cell_path <- file.path(base_dir, paste0(cell, "_comp_all"))
  
  for (i in 1:nrow(valid_map)) {
    short_id <- valid_map$Short_Name[i]
    full_id  <- valid_map$Full_Name[i]
    
    motif_path <- file.path(cell_path, full_id)
    bed_dir <- file.path(motif_path, "beds")
    
    # Check files using FULL ID
    f_mae <- file.path(bed_dir, paste0(full_id, "_MAE_bound.bed"))
    f_mas <- file.path(bed_dir, paste0(full_id, "_MAS_bound.bed"))
    f_ys  <- file.path(bed_dir, paste0(full_id, "_YS_bound.bed"))
    
    df_mae <- import_bed_safe(f_mae)
    df_mas <- import_bed_safe(f_mas)
    df_ys  <- import_bed_safe(f_ys)
    
    if (is.null(df_mae) && is.null(df_mas) && is.null(df_ys)) next
    
    combined <- bind_rows(df_mae, df_mas, df_ys) %>%
      distinct(`TFBS chr`, `TFBS start`, `TFBS end`, .keep_all = TRUE)
    
    combined$Motif <- short_id
    combined$Full_Motif_ID <- full_id
    combined$`cell type` <- cell
    combined$Group <- ifelse(short_id %in% Rev_TF, "AP-1", "Control")
    
    results_list[[counter]] <- combined
    counter <- counter + 1
  }
}

# ==============================================================================
# 5. Save Output (With Rename)
# ==============================================================================
if (length(results_list) > 0) {
  final_df <- bind_rows(results_list)
  
  # --- RENAME MSC to HSPC ---
  # Check if MSC exists and replace it
  final_df <- final_df %>%
    mutate(`cell type` = ifelse(`cell type` == "MSC", "HSPC", `cell type`))
  
  # Add Pairing ID
  pair_map <- matched_pairs %>% select(AP1_Motif, Control_Motif)
  
  ap1_part <- final_df %>% filter(Group == "AP-1") %>% mutate(Pair_ID = Motif)
  ctrl_part <- final_df %>% filter(Group == "Control") %>%
    left_join(pair_map, by = c("Motif" = "Control_Motif")) %>%
    mutate(Pair_ID = AP1_Motif) %>%
    select(-AP1_Motif)
  
  final_df <- bind_rows(ap1_part, ctrl_part) %>%
    select(Motif, Group, Pair_ID, `cell type`, `Peak chr`, `Peak start`, `Peak End`, `TFBS chr`, `TFBS start`, `TFBS end`)

  cat("Total Unique Sites:", nrow(final_df), "\n")
  cat("Renamed 'MSC' to 'HSPC'.\n")
  
  fwrite(final_df, out_file)
  cat("Saved coordinates to:", out_file, "\n")
  
} else {
  warning("No sites found.")
}

library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

# ==============================================================================
# 1. Setup & Path Discovery
# ==============================================================================

base_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Output/BINDetect")
out_file <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Classification/AP1_Matched_Control_Peak_Coordinates_Classified.csv")

# 1. Define Short Names (Target List)
Rev_TF <- c('Fos','Fosb','Fosl1','Fosl2',
            'Junb','Jund','Jun',
            'Atf2','Atf3','Atf4','Batf','Batf3',
            'Maf','Mafa','Mafb','Mafg')

# Load Matched Pairs to get Control Short Names
pair_file <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Footprints/AP1_Matched_Control_Pairs.csv")
if (file.exists(pair_file)) {
  matched_pairs <- read.csv(pair_file)
  Control_TF <- unique(matched_pairs$Control_Motif)
} else {
  stop("Matched Pairs file not found. Please run the pairing step first.")
}

Target_Short_Names <- unique(c(Rev_TF, Control_TF))

# 2. Get Cell Types
comp_dirs <- list.dirs(base_dir, recursive = FALSE)
cells <- sub("_comp_all$", "", basename(comp_dirs[grep("_comp_all$", comp_dirs)]))
cells <- cells[!cells %in% c("M1", "M2", "LSC", "M1M2")]

# ==============================================================================
# 2. Build Name Map (Short Name -> Full Directory Name)
# ==============================================================================
example_cell_dir <- file.path(base_dir, paste0(cells[1], "_comp_all"))
all_motif_folders <- list.dirs(example_cell_dir, full.names = FALSE, recursive = FALSE)

get_full_name <- function(short, all_folders) {
  pattern <- paste0("^", short, "(_|$)")
  matches <- all_folders[grepl(pattern, all_folders)]
  if (length(matches) == 0) return(NA)
  if (length(matches) > 1) return(matches[which.min(nchar(matches))]) 
  return(matches)
}

name_map <- data.frame(Short_Name = Target_Short_Names, stringsAsFactors = FALSE)
name_map$Full_Name <- sapply(name_map$Short_Name, get_full_name, all_folders = all_motif_folders)
valid_map <- name_map %>% filter(!is.na(Full_Name))

cat("Mapped", nrow(valid_map), "motifs to full folder names.\n")

# ==============================================================================
# 3. Helper Function
# ==============================================================================
import_bed_safe <- function(path) {
  if (file.exists(path) && file.info(path)$size > 0) {
    tryCatch({
      df <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      if (ncol(df) < 9) return(NULL)
      
      res <- df[, 1:3]
      colnames(res) <- c("TFBS chr", "TFBS start", "TFBS end")
      res$`Peak chr`   <- df[, 7]
      res$`Peak start` <- df[, 8]
      res$`Peak End`   <- df[, 9]
      
      # Create unique identifier for merging logic later
      res$site_id <- paste(res$`TFBS chr`, res$`TFBS start`, res$`TFBS end`, sep="-")
      
      return(res)
    }, error = function(e) return(NULL))
  } else {
    return(NULL)
  }
}

# ==============================================================================
# 4. Main Loop with Classification
# ==============================================================================
results_list <- list()
counter <- 1

cat("Extracting coordinates and classifying binding patterns...\n")

for (cell in cells) {
  cell_path <- file.path(base_dir, paste0(cell, "_comp_all"))
  
  for (i in 1:nrow(valid_map)) {
    short_id <- valid_map$Short_Name[i]
    full_id  <- valid_map$Full_Name[i]
    
    motif_path <- file.path(cell_path, full_id)
    bed_dir <- file.path(motif_path, "beds")
    
    # Check files using FULL ID
    f_mae <- file.path(bed_dir, paste0(full_id, "_MAE_bound.bed"))
    f_mas <- file.path(bed_dir, paste0(full_id, "_MAS_bound.bed"))
    f_ys  <- file.path(bed_dir, paste0(full_id, "_YS_bound.bed"))
    
    df_mae <- import_bed_safe(f_mae)
    df_mas <- import_bed_safe(f_mas)
    df_ys  <- import_bed_safe(f_ys)
    
    if (is.null(df_mae) && is.null(df_mas) && is.null(df_ys)) next
    
    # --- Classification Logic ---
    # 1. Combine all potential sites found in any condition
    combined <- bind_rows(df_mae, df_mas, df_ys) %>%
      distinct(`TFBS chr`, `TFBS start`, `TFBS end`, .keep_all = TRUE)
    
    # 2. Check presence in each specific condition
    if (!is.null(df_mae)) combined$in_MAE <- combined$site_id %in% df_mae$site_id else combined$in_MAE <- FALSE
    if (!is.null(df_mas)) combined$in_MAS <- combined$site_id %in% df_mas$site_id else combined$in_MAS <- FALSE
    if (!is.null(df_ys))  combined$in_YS  <- combined$site_id %in% df_ys$site_id  else combined$in_YS  <- FALSE
    
    # 3. Assign Category
    combined <- combined %>%
      mutate(category = case_when(
        in_MAE & in_MAS & in_YS ~ "In All Conditions",
        
        # Lost with Age (Present in YS, Lost in MAS/MAE)
        in_YS & !in_MAS & !in_MAE ~ "In Young Sedentary Only",
        
        # Gained with Age (Absent in YS, Present in MAS/MAE)
        !in_YS & in_MAS & in_MAE ~ "Gained in Age (Both)",
        
        # Exercise specific restoration/loss patterns
        in_MAE & in_YS & !in_MAS ~ "Lost AgeSed (Rescued by Ex)", 
        in_MAE & in_MAS & !in_YS ~ "Lost YoungSed (Gained in Age)", 
        in_MAS & in_YS & !in_MAE ~ "Lost AgeEx",       
        
        # Exclusive conditions
        in_MAE & !in_YS & !in_MAS ~ "In Aged Ex Only",
        in_MAS & !in_YS & !in_MAE ~ "In Aged Sed Only",
        
        TRUE ~ "Complex/Other"
      ))

    # Add Metadata
    combined$Motif <- short_id
    combined$Full_Motif_ID <- full_id
    combined$`cell type` <- cell
    combined$Group <- ifelse(short_id %in% Rev_TF, "AP-1", "Control")
    
    results_list[[counter]] <- combined
    counter <- counter + 1
  }
}

# ==============================================================================
# 5. Save Output (With Rename)
# ==============================================================================
if (length(results_list) > 0) {
  final_df <- bind_rows(results_list)
  
  # --- RENAME MSC to HSPC ---
  final_df <- final_df %>%
    mutate(`cell type` = ifelse(`cell type` == "MSC", "HSPC", `cell type`))
  
  # Add Pairing ID
  pair_map <- matched_pairs %>% select(AP1_Motif, Control_Motif)
  
  ap1_part <- final_df %>% filter(Group == "AP-1") %>% mutate(Pair_ID = Motif)
  ctrl_part <- final_df %>% filter(Group == "Control") %>%
    left_join(pair_map, by = c("Motif" = "Control_Motif")) %>%
    mutate(Pair_ID = AP1_Motif) %>%
    select(-AP1_Motif)
  
  # Final Selection (Now includes 'category')
  final_df <- bind_rows(ap1_part, ctrl_part) %>%
    select(Motif, Group, Pair_ID, category, `cell type`, 
           `Peak chr`, `Peak start`, `Peak End`, 
           `TFBS chr`, `TFBS start`, `TFBS end`)

  cat("Total Unique Sites:", nrow(final_df), "\n")
  cat("Renamed 'MSC' to 'HSPC'.\n")
  cat("Preview of Categories:\n")
  print(table(final_df$category))
  
  fwrite(final_df, out_file)
  cat("Saved coordinates to:", out_file, "\n")
  
} else {
  warning("No sites found.")
}

library(dplyr)
library(stringr)
library(readr)

# --- 1. Load and Merge All Classification Results ---
indir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Classification")
rds_files <- list.files(indir, pattern = "_containingPeaks_Classified.rds", full.names = TRUE)

message(paste("Loading", length(rds_files), "motif files..."))
master_binding_df <- bind_rows(lapply(rds_files, readRDS))

# --- 2. Clean Motif Names ---
master_binding_df <- master_binding_df %>%
  mutate(
    clean_motif = motif_id %>% 
      str_remove("^\\(") %>%       
      str_remove("(_|\\().*") %>%   
      str_remove("\\)$")            
  )

# --- 3. Updated Function: Group by PARENT PEAK ---
create_cobinding_report <- function(df, condition_col, out_filename) {
  
  message(paste("Generating report for:", condition_col))
  
  # Filter: Keep rows where the specific condition is TRUE
  subset_df <- df %>% 
    filter(.data[[condition_col]] == TRUE)
  
  # --- KEY CHANGE HERE ---
  # We group by the 'peak_identifier' (The broad consensus peak)
  # rather than the 'identifier' (The tiny TFBS footprint).
  
  co_binding_df <- subset_df %>%
    group_by(cell_type, peak_identifier, peak_chr, peak_start, peak_end) %>%
    summarize(
      motif_count = n_distinct(clean_motif),
      # This collapses all unique motifs found in this ONE peak into a list
      motifs = paste(sort(unique(clean_motif)), collapse = ","),
      .groups = "drop"
    )
  
  # Save
  out_path <- file.path(indir, out_filename)
  write_csv(co_binding_df, out_path)
  message(paste("Saved:", out_path))
}

# --- 4. Run for Each Condition ---

# Young Sedentary (YS)
create_cobinding_report(master_binding_df, "in_YS", "YS_Tobias_Cobinding.csv")

# Aged Sedentary (MAS)
create_cobinding_report(master_binding_df, "in_MAS", "MAS_Tobias_Cobinding.csv")

# Aged Exercise (MAE)
create_cobinding_report(master_binding_df, "in_MAE", "MAE_Tobias_Cobinding.csv")

library(dplyr)
library(stringr)
library(readr)
library(tidyr) # Required for pivot_wider

# ==============================================================================
# 1. Load and Merge All Classification Results
# ==============================================================================
indir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Classification")
rds_files <- list.files(indir, pattern = "_containingPeaks_Classified.rds", full.names = TRUE)

message(paste("Loading", length(rds_files), "motif files..."))
master_binding_df <- bind_rows(lapply(rds_files, readRDS))

# ==============================================================================
# 2. Clean Motif Names
# ==============================================================================
master_binding_df <- master_binding_df %>%
  mutate(
    clean_motif = motif_id %>% 
      str_remove("^\\(") %>%        
      str_remove("(_|\\().*") %>%    
      str_remove("\\)$")            
  )

# ==============================================================================
# 3. New Function: Create Peak x Motif Count Matrix
# ==============================================================================
create_motif_matrix <- function(df, condition_col, out_filename) {
  
  message(paste("Generating count matrix for:", condition_col))
  
  # 1. Filter for the specific condition (e.g., in_YS)
  subset_df <- df %>% 
    filter(.data[[condition_col]] == TRUE)
  
  # 2. Create the Matrix
  # - Unite peak and celltype to make a unique row ID
  # - Count occurrences of each motif per peak
  # - Pivot to wide format (Motifs become columns)
  motif_matrix <- subset_df %>%
    mutate(peak_celltype = paste(peak_identifier, cell_type, sep = "_")) %>%
    count(peak_celltype, clean_motif) %>%
    pivot_wider(
      names_from = clean_motif, 
      values_from = n, 
      values_fill = 0
    )
  
  # 3. Save
  out_path <- file.path(indir, out_filename)
  write_csv(motif_matrix, out_path)
  message(paste("Saved matrix with dimensions:", paste(dim(motif_matrix), collapse = "x")))
  message(paste("File saved to:", out_path, "\n"))
}

# ==============================================================================
# 4. Run for Each Condition
# ==============================================================================

# Young Sedentary (YS)
create_motif_matrix(master_binding_df, "in_YS", "YS_Motif_Count_Matrix.csv")

# Aged Sedentary (MAS)
create_motif_matrix(master_binding_df, "in_MAS", "MAS_Motif_Count_Matrix.csv")

# Aged Exercise (MAE)
create_motif_matrix(master_binding_df, "in_MAE", "MAE_Motif_Count_Matrix.csv")

# Check the results 
cob_MAS = read.csv(publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Classification/MAS_Motif_Count_Matrix.csv"))
cob_MAE = read.csv(publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Classification/MAE_Motif_Count_Matrix.csv"))
cob_YS = read.csv(publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Classification/YS_Motif_Count_Matrix.csv"))

table(cob_MAS > 1)
table(cob_MAE > 1)
table(cob_YS > 1)

# Check the results 
cob_MAS = read.csv(publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Classification/MAS_Tobias_Cobinding.csv"))
cob_MAE = read.csv(publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Classification/MAE_Tobias_Cobinding.csv"))
cob_YS = read.csv(publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Classification/YS_Tobias_Cobinding.csv"))

head(cob_MAS)
head(cob_MAE)

library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)

# --- Configuration ---
base_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Classification")
plot_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Classification")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# --- Cell Type Order (for sorting) ---
target_col_order <- c(
  "Cardiomyocyte", "Coronary_EC", "Endocardial_EC", 
  "Lymphatic_EC", "Fibroblast", "MSC", "M1M2", "Pericytes"
)

# --- Helper: Diversity Selection (Limits Sox family to 3) ---
select_diverse_top_hits <- function(sorted_genes, n_limit = 5, family_cap = 3, family_pattern = "^Sox") {
  selected <- c()
  family_count <- 0
  for (gene in names(sorted_genes)) {
    if (length(selected) >= n_limit) break
    is_family <- str_detect(gene, family_pattern)
    if (is_family) {
      if (family_count < family_cap) {
        selected <- c(selected, gene)
        family_count <- family_count + 1
      }
    } else {
      selected <- c(selected, gene)
    }
  }
  return(selected)
}

# --- PART 1: Analysis Function ---

get_partner_data <- function(csv_filename, condition_label, anchor_pattern = "Fos") {
  
  message(paste("Analyzing:", condition_label))
  infile <- file.path(base_dir, csv_filename)
  if (!file.exists(infile)) return(NULL)
  
  data <- read_csv(infile, show_col_types = FALSE)
  
  # Filter for Anchor
  anchor_peaks <- data %>%
    filter(str_detect(motifs, anchor_pattern)) %>%
    filter(!cell_type %in% c("M1", "M2", "LSC"))
  
  if (nrow(anchor_peaks) < 10) return(NULL)
  
  partner_stats <- anchor_peaks %>%
    separate_rows(motifs, sep = ",") %>%
    rename(partner = motifs) %>%
    filter(!str_detect(partner, anchor_pattern)) %>%
    filter(!partner %in% c("Jun", "Junb", "Jund")) %>% 
    group_by(cell_type, partner) %>%
    summarise(count = n(), .groups = "drop_last") %>%
    mutate(
      total_anchor_peaks = n_distinct(anchor_peaks$peak_identifier[anchor_peaks$cell_type == unique(cell_type)]),
      freq = count / total_anchor_peaks
    ) %>%
    ungroup()
  
  # Matrices
  mat_freq <- partner_stats %>%
    select(cell_type, partner, freq) %>%
    pivot_wider(names_from = cell_type, values_from = freq, values_fill = 0) %>%
    column_to_rownames("partner") %>%
    as.matrix()
  
  mat_count <- partner_stats %>%
    select(cell_type, partner, count) %>%
    pivot_wider(names_from = cell_type, values_from = count, values_fill = 0) %>%
    column_to_rownames("partner") %>%
    as.matrix()
  
  # Filter low signal
  keep_partners <- rownames(mat_freq)[apply(mat_freq, 1, max) > 0.05]
  mat_filtered <- mat_freq[keep_partners, , drop = FALSE]
  
  if (nrow(mat_filtered) == 0) return(NULL)
  
  mat_zscore <- t(scale(t(mat_filtered)))
  mat_zscore[is.na(mat_zscore)] <- 0
  
  # --- SELECTION & TRACKING (Needed for Sorting) ---
  hit_tracker <- data.frame(Gene = character(), Origin_Cell = character(), Z_Score_At_Selection = numeric(), stringsAsFactors = FALSE)
  
  # Per-Cell Selection
  for (cell in colnames(mat_zscore)) {
    cell_vals <- sort(mat_zscore[, cell], decreasing = TRUE)
    selected_tfs <- select_diverse_top_hits(cell_vals, n_limit = 5, family_cap = 3)
    
    if(length(selected_tfs) > 0) {
      hit_tracker <- rbind(hit_tracker, data.frame(
        Gene = selected_tfs, 
        Origin_Cell = cell, 
        Z_Score_At_Selection = mat_zscore[selected_tfs, cell]
      ))
    }
  }
  
  return(list(
    label = condition_label,
    mat_freq = mat_freq,      
    mat_count = mat_count,    
    hit_tracker = hit_tracker
  ))
}

# --- PART 2: Execute & Build Tracker ---

results <- list()
results[["YS"]]  <- get_partner_data("YS_Tobias_Cobinding.csv", "YS")
results[["MAS"]] <- get_partner_data("MAS_Tobias_Cobinding.csv", "MAS")
results[["MAE"]] <- get_partner_data("MAE_Tobias_Cobinding.csv", "MAE")

# Combine trackers
all_hits_long <- bind_rows(
  results[["YS"]]$hit_tracker %>% mutate(Cond = "YS", Cond_Order = 1),
  results[["MAS"]]$hit_tracker %>% mutate(Cond = "MAS", Cond_Order = 2),
  results[["MAE"]]$hit_tracker %>% mutate(Cond = "MAE", Cond_Order = 3)
)

# Get unique genes
all_top_hits <- unique(all_hits_long$Gene)
message(paste("Total Unique Partners:", length(all_top_hits)))

# --- PART 3: Strict Sorting Logic (Internal) ---
# We keep this logic to ensure the heatmap has a nice diagonal structure,
# even though we aren't printing the origin on the label anymore.

gene_sort_df <- all_hits_long %>%
  # Map Origin_Cell to Index
  mutate(Cell_Index = match(Origin_Cell, target_col_order)) %>%
  # Handle unknown cells (put at end)
  mutate(Cell_Index = ifelse(is.na(Cell_Index), 99, Cell_Index)) %>%
  
  # Group by Gene to summarize
  group_by(Gene) %>%
  summarize(
    # Min Index ensures if it's selected by CM and Fib, it stays with CM
    Primary_Sort = min(Cell_Index),
    # Min Condition ensures if it's in YS and MAS, it sorts with YS
    Secondary_Sort = min(Cond_Order),
    # Max Z-score for tie-breaking
    Tertiary_Sort = max(Z_Score_At_Selection),
    .groups = "drop"
  ) %>%
  # Apply Sorting
  arrange(Primary_Sort, Secondary_Sort, desc(Tertiary_Sort))

diagonal_row_order <- gene_sort_df$Gene

message("Strict diagonal order established.")

# --- PART 4: Plotting Function ---

plot_final_heatmap <- function(result_obj, master_genes, master_order, col_order) {
  
  if (is.null(result_obj)) return()
  
  message(paste("Plotting:", result_obj$label))
  
  # Align Data
  plot_mat <- matrix(0, nrow = length(master_genes), ncol = ncol(result_obj$mat_freq))
  rownames(plot_mat) <- master_genes
  colnames(plot_mat) <- colnames(result_obj$mat_freq)
  genes_present <- intersect(rownames(result_obj$mat_freq), master_genes)
  plot_mat[genes_present, ] <- result_obj$mat_freq[genes_present, ]
  
  mat_scaled <- t(scale(t(plot_mat)))
  mat_scaled[is.na(mat_scaled)] <- 0
  
  # Clean Labels: Gene (Count)
  new_names <- sapply(rownames(mat_scaled), function(tf) {
    count_str <- "(0)"
    if (tf %in% rownames(result_obj$mat_count)) {
      max_cell <- names(which.max(plot_mat[tf, ]))
      cnt <- result_obj$mat_count[tf, max_cell]
      count_str <- paste0("(", cnt, ")")
    }
    # No origin string here
    return(paste(tf, count_str))
  })
  
  # Reorder using STRICT master order
  valid_rows <- intersect(master_order, rownames(mat_scaled))
  mat_final <- mat_scaled[valid_rows, , drop = FALSE]
  labels_final <- new_names[valid_rows]
  
  # Reorder Columns
  valid_cols <- intersect(col_order, colnames(mat_final))
  if(length(valid_cols) > 0) mat_final <- mat_final[, valid_cols, drop = FALSE]

  # Plot
  out_pdf <- file.path(plot_dir, paste0("Fos_Partners_Unified_", result_obj$label, ".pdf"))
  
  max_abs <- max(abs(mat_final))
  limit <- min(max_abs, 2.5) 
  col_fun <- colorRamp2(c(-limit, 0, limit), c("blue", "white", "red"))
  
  pdf(out_pdf, width = 8, height = 12) # Adjusted width since labels are shorter
  ht <- Heatmap(mat_final, 
                name = "Z-Score", 
                col = col_fun,
                cluster_rows = FALSE,     
                row_labels = labels_final,
                cluster_columns = FALSE,  
                column_order = valid_cols, 
                column_title = paste0(result_obj$label, " (Strict Diagonal)"),
                row_names_gp = gpar(fontsize = 9),
                column_names_gp = gpar(fontsize = 10, fontface = "bold")
  )
  draw(ht)
  dev.off()
}

# --- PART 5: Generate All Plots ---
for (res in results) {
  plot_final_heatmap(res, all_top_hits, diagonal_row_order, target_col_order)
}

message("All unified plots generated.")

library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)

# --- Configuration ---
base_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Classification")
plot_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Classification")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# --- Cell Type Order ---
target_col_order <- c(
  "Cardiomyocyte", "Coronary_EC", "Endocardial_EC", 
  "Lymphatic_EC", "Fibroblast", "MSC", "M1M2", "Pericytes"
)

# --- PART 1: Data Prep (Recreating Exclusive File with Peak Logic) ---
rds_files <- list.files(base_dir, pattern = "_containingPeaks_Classified.rds", full.names = TRUE)
exclusive_outfile <- file.path(base_dir, "MAS_Exclusive_Cobinding.csv")

# We Force Re-creation to ensure Peak Logic is applied
if (TRUE) { 
  message("Generating Exclusive AgedSed Dataset (Using Peak Logic)...")
  
  # 1. Load All Data
  master_df <- bind_rows(lapply(rds_files, readRDS))
  
  # 2. Clean Motif Names
  master_df <- master_df %>%
    mutate(clean_motif = motif_id %>% str_remove("^\\(") %>% str_remove("(_|\\().*") %>% str_remove("\\)$"))
  
  # 3. Filter strictly for "In Aged Sed Only"
  exclusive_df <- master_df %>%
    filter(category == "In Aged Sed Only") %>%
    filter(!cell_type %in% c("M1", "M2", "LSC"))
  
  # 4. Aggregate by PARENT PEAK (peak_identifier)
  # This matches the logic of the Global file: identifying partners in the same peak
  exclusive_cobinding <- exclusive_df %>%
    group_by(cell_type, peak_identifier, peak_chr, peak_start, peak_end) %>%
    summarize(
      motifs = paste(sort(unique(clean_motif)), collapse = ","),
      .groups = "drop"
    )
  
  write_csv(exclusive_cobinding, exclusive_outfile)
  message("Saved Exclusive MAS data.")
}

# --- Helper: Diversity Selection (Limits Sox family to 3) ---
select_diverse_top_hits <- function(sorted_genes, n_limit = 5, family_cap = 3, family_pattern = "^Sox") {
  selected <- c()
  family_count <- 0
  for (gene in names(sorted_genes)) {
    if (length(selected) >= n_limit) break
    is_family <- str_detect(gene, family_pattern)
    if (is_family) {
      if (family_count < family_cap) {
        selected <- c(selected, gene)
        family_count <- family_count + 1
      }
    } else {
      selected <- c(selected, gene)
    }
  }
  return(selected)
}

# --- PART 2: Analysis Function ---

get_partner_data <- function(csv_filename, condition_label, anchor_pattern = "Fos") {
  
  message(paste("Analyzing:", condition_label))
  infile <- file.path(base_dir, csv_filename)
  if (!file.exists(infile)) return(NULL)
  
  data <- read_csv(infile, show_col_types = FALSE)
  
  # Filter for Anchor
  anchor_peaks <- data %>%
    filter(str_detect(motifs, anchor_pattern)) %>%
    filter(!cell_type %in% c("M1", "M2", "LSC"))
  
  if (nrow(anchor_peaks) < 10) return(NULL)
  
  # Calculate Co-binding
  partner_stats <- anchor_peaks %>%
    separate_rows(motifs, sep = ",") %>%
    rename(partner = motifs) %>%
    filter(!str_detect(partner, anchor_pattern)) %>%
    filter(!partner %in% c("Jun", "Junb", "Jund")) %>% 
    group_by(cell_type, partner) %>%
    summarise(count = n(), .groups = "drop_last") %>%
    mutate(
      total_anchor_peaks = n_distinct(anchor_peaks$peak_identifier[anchor_peaks$cell_type == unique(cell_type)]),
      freq = count / total_anchor_peaks
    ) %>%
    ungroup()
  
  # Matrices
  mat_freq <- partner_stats %>%
    select(cell_type, partner, freq) %>%
    pivot_wider(names_from = cell_type, values_from = freq, values_fill = 0) %>%
    column_to_rownames("partner") %>%
    as.matrix()
  
  mat_count <- partner_stats %>%
    select(cell_type, partner, count) %>%
    pivot_wider(names_from = cell_type, values_from = count, values_fill = 0) %>%
    column_to_rownames("partner") %>%
    as.matrix()
  
  # Filter low signal
  keep_partners <- rownames(mat_freq)[apply(mat_freq, 1, max) > 0.05]
  mat_filtered <- mat_freq[keep_partners, , drop = FALSE]
  
  if (nrow(mat_filtered) == 0) return(NULL)
  
  mat_zscore <- t(scale(t(mat_filtered)))
  mat_zscore[is.na(mat_zscore)] <- 0
  
  # --- SELECTION & TRACKING ---
  hit_tracker <- data.frame(Gene = character(), Origin_Cell = character(), Z_Score_At_Selection = numeric(), stringsAsFactors = FALSE)
  
  for (cell in colnames(mat_zscore)) {
    cell_vals <- sort(mat_zscore[, cell], decreasing = TRUE)
    selected_tfs <- select_diverse_top_hits(cell_vals, n_limit = 5, family_cap = 3)
    
    if(length(selected_tfs) > 0) {
      hit_tracker <- rbind(hit_tracker, data.frame(
        Gene = selected_tfs, 
        Origin_Cell = cell, 
        Z_Score_At_Selection = mat_zscore[selected_tfs, cell]
      ))
    }
  }
  
  return(list(
    label = condition_label,
    mat_freq = mat_freq,      
    mat_count = mat_count,    
    hit_tracker = hit_tracker
  ))
}

# --- PART 3: Execute Analysis ---

results <- list()
# Run on Global MAS
results[["Global"]] <- get_partner_data("MAS_Tobias_Cobinding.csv", "Global_MAS")
# Run on Exclusive MAS
results[["Exclusive"]] <- get_partner_data("MAS_Exclusive_Cobinding.csv", "Exclusive_MAS")

# Combine trackers
all_hits_long <- bind_rows(
  results[["Global"]]$hit_tracker %>% mutate(Cond = "Global", Cond_Order = 1),
  results[["Exclusive"]]$hit_tracker %>% mutate(Cond = "Exclusive", Cond_Order = 2)
)

# Build Origin Labels (e.g., "[Global]" or "[Global; Exclusive]")
origin_map_df <- all_hits_long %>%
  group_by(Gene) %>%
  summarize(
    Origin_Label = paste0("[", paste(unique(Cond), collapse = "; "), "]"),
    .groups = "drop"
  )
origin_map <- setNames(origin_map_df$Origin_Label, origin_map_df$Gene)
all_top_hits <- origin_map_df$Gene

message(paste("Total Unique Partners:", length(all_top_hits)))

# --- PART 4: Sorting (Cell Type -> Condition -> Score) ---

gene_sort_df <- all_hits_long %>%
  # Map Origin_Cell to Index
  mutate(Cell_Index = match(Origin_Cell, target_col_order)) %>%
  mutate(Cell_Index = ifelse(is.na(Cell_Index), 99, Cell_Index)) %>%
  
  group_by(Gene) %>%
  summarize(
    Primary_Sort = min(Cell_Index),
    Secondary_Sort = min(Cond_Order),
    Tertiary_Sort = max(Z_Score_At_Selection),
    .groups = "drop"
  ) %>%
  arrange(Primary_Sort, Secondary_Sort, desc(Tertiary_Sort))

diagonal_row_order <- gene_sort_df$Gene

message("Strict diagonal order established.")

# --- PART 5: Plotting Function ---

plot_final_heatmap <- function(result_obj, master_genes, master_order, col_order, origins) {
  
  if (is.null(result_obj)) return()
  
  message(paste("Plotting:", result_obj$label))
  
  # Align Data
  plot_mat <- matrix(0, nrow = length(master_genes), ncol = ncol(result_obj$mat_freq))
  rownames(plot_mat) <- master_genes
  colnames(plot_mat) <- colnames(result_obj$mat_freq)
  genes_present <- intersect(rownames(result_obj$mat_freq), master_genes)
  plot_mat[genes_present, ] <- result_obj$mat_freq[genes_present, ]
  
  mat_scaled <- t(scale(t(plot_mat)))
  mat_scaled[is.na(mat_scaled)] <- 0
  
  # Create Labels: Gene (Count) [Origin]
  new_names <- sapply(rownames(mat_scaled), function(tf) {
    count_str <- "(0)"
    if (tf %in% rownames(result_obj$mat_count)) {
      max_cell <- names(which.max(plot_mat[tf, ]))
      cnt <- result_obj$mat_count[tf, max_cell]
      count_str <- paste0("(", cnt, ")")
    }
    
    # Add Origin Annotation
    origin_str <- ifelse(tf %in% names(origins), origins[tf], "")
    
    return(paste(tf, count_str, origin_str))
  })
  
  # Reorder
  valid_rows <- intersect(master_order, rownames(mat_scaled))
  mat_final <- mat_scaled[valid_rows, , drop = FALSE]
  labels_final <- new_names[valid_rows]
  
  valid_cols <- intersect(col_order, colnames(mat_final))
  if(length(valid_cols) > 0) mat_final <- mat_final[, valid_cols, drop = FALSE]

  out_pdf <- file.path(plot_dir, paste0("Fos_Partners_", result_obj$label, "_Annotated.pdf"))
  
  max_abs <- max(abs(mat_final))
  limit <- min(max_abs, 2.5) 
  col_fun <- colorRamp2(c(-limit, 0, limit), c("blue", "white", "red"))
  
  pdf(out_pdf, width = 9, height = 12) 
  ht <- Heatmap(mat_final, 
                name = "Z-Score", 
                col = col_fun,
                cluster_rows = FALSE,     
                row_labels = labels_final,
                cluster_columns = FALSE,  
                column_order = valid_cols, 
                column_title = paste0(result_obj$label, " (Sorted by Cell & Origin)"),
                row_names_gp = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 10, fontface = "bold")
  )
  draw(ht)
  dev.off()
}

# --- PART 6: Generate Plots ---
for (res in results) {
  plot_final_heatmap(res, all_top_hits, diagonal_row_order, target_col_order, origin_map)
}

message("Annotated plots generated.")

library(dplyr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(knitr) # Optional, for nicer table printing

# --- 1. Load Data (If needed) ---
# results <- list()
# results[["Global"]] <- get_partner_data("MAS_Tobias_Cobinding.csv", "Global_MAS")
# results[["Exclusive"]] <- get_partner_data("MAS_Exclusive_Cobinding.csv", "Exclusive_MAS")

# --- 2. Statistical Comparison Function with Summary Printing ---

compare_partner_enrichment <- function(res_global, res_exclusive, cell_types) {
  
  message("Running Fisher's Exact Tests for Differential Co-binding...")
  
  stats_list <- list()
  summary_list <- list() # To store N and counts
  
  # Align Genes
  all_genes <- unique(c(rownames(res_global$mat_count), rownames(res_exclusive$mat_count)))
  
  for (cell in cell_types) {
    
    # --- 1. Recover Total Peak Counts (N) ---
    # Global
    g_counts <- res_global$mat_count[, cell]
    g_freqs  <- res_global$mat_freq[, cell]
    valid_idx_g <- which(g_counts > 0)[1] 
    
    if (is.na(valid_idx_g)) {
      n_global <- 0
    } else {
      n_global <- round(g_counts[valid_idx_g] / g_freqs[valid_idx_g])
    }
    
    # Exclusive
    e_counts <- res_exclusive$mat_count[, cell]
    e_freqs  <- res_exclusive$mat_freq[, cell]
    valid_idx_e <- which(e_counts > 0)[1]
    
    if (is.na(valid_idx_e)) {
      n_exclusive <- 0
    } else {
      n_exclusive <- round(e_counts[valid_idx_e] / e_freqs[valid_idx_e])
    }
    
    # Store N for this cell type
    current_summary <- data.frame(
      Cell_Type = cell,
      Total_N_Global = n_global,
      Total_N_Exclusive = n_exclusive,
      stringsAsFactors = FALSE
    )
    
    # Skip if data missing
    if (n_global < 10 || n_exclusive < 10) {
      current_summary$Sig_Hits <- 0
      summary_list[[cell]] <- current_summary
      next
    }
    
    # --- 2. Iterate through Genes ---
    for (gene in all_genes) {
      
      k_global <- ifelse(gene %in% names(g_counts), g_counts[gene], 0)
      k_excl   <- ifelse(gene %in% names(e_counts), e_counts[gene], 0)
      
      if (k_global < 5 && k_excl < 5) next
      
      mat <- matrix(c(k_excl, n_exclusive - k_excl, 
                      k_global, n_global - k_global), 
                    nrow = 2, byrow = TRUE)
      
      ft <- fisher.test(mat)
      
      # Log2FC (with pseudocount)
      f_excl <- (k_excl + 0.1) / (n_exclusive + 1)
      f_glob <- (k_global + 0.1) / (n_global + 1)
      l2fc <- log2(f_excl / f_glob)
      
      stats_list[[length(stats_list) + 1]] <- data.frame(
        Cell_Type = cell,
        Gene = gene,
        Count_Global = k_global,
        Count_Exclusive = k_excl,
        P_Value = ft$p.value,
        Log2FC = l2fc
      )
    }
    
    # Store temp summary (will update sig hits count later)
    summary_list[[cell]] <- current_summary
  }
  
  # --- 3. Process Stats & Finalize Summary ---
  df_stats <- bind_rows(stats_list) %>%
    mutate(FDR = p.adjust(P_Value, method = "BH"))
  
  # Count Significant Hits per Cell (FDR < 0.05 AND Log2FC > 1)
  sig_counts <- df_stats %>%
    filter(FDR < 0.05 & abs(Log2FC) > 1.0) %>% # Strict Filter
    group_by(Cell_Type) %>%
    summarise(Sig_Hits = n())
  
  # Merge Counts into Summary Table
  final_summary <- bind_rows(summary_list) %>%
    left_join(sig_counts, by = "Cell_Type") %>%
    mutate(Sig_Hits = ifelse(is.na(Sig_Hits), 0, Sig_Hits))
  
  # --- 4. PRINT THE TABLE ---
  message("\n========================================================")
  message("          Sample Sizes & Significant Differences         ")
  message("          (Sig: FDR < 0.05, |Log2FC| > 1.0)              ")
  message("========================================================")
  print(final_summary)
  message("========================================================\n")
  
  return(df_stats)
}

# --- 3. Run Comparison ---
target_cols <- c("Cardiomyocyte", "Coronary_EC", "Endocardial_EC", 
                 "Lymphatic_EC", "Fibroblast", "MSC", "M1M2", "Pericytes")

diff_stats <- compare_partner_enrichment(results[["Global"]], results[["Exclusive"]], target_cols)

# --- 4. Plotting (Rest of pipeline) ---

# Filter Top 50 Extreme Changes
sig_hits <- diff_stats %>%
  filter(FDR < 0.05 & abs(Log2FC) > 1.0)

if (nrow(sig_hits) > 50) {
  top_genes <- sig_hits %>%
    group_by(Gene) %>%
    summarize(Max_Change = max(abs(Log2FC))) %>%
    arrange(desc(Max_Change)) %>%
    head(50) %>%
    pull(Gene)
  sig_hits <- sig_hits %>% filter(Gene %in% top_genes)
}

# Prepare Matrices
plot_genes <- unique(sig_hits$Gene)

mat_l2fc <- diff_stats %>%
  filter(Gene %in% plot_genes) %>%
  dplyr::select(Gene, Cell_Type, Log2FC) %>%
  pivot_wider(names_from = Cell_Type, values_from = Log2FC, values_fill = 0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

mat_sig <- diff_stats %>%
  filter(Gene %in% plot_genes) %>%
  dplyr::select(Gene, Cell_Type, FDR) %>%
  mutate(Star = ifelse(FDR < 0.001, "***", ifelse(FDR < 0.01, "**", ifelse(FDR < 0.05, "*", "")))) %>%
  dplyr::select(Gene, Cell_Type, Star) %>%
  pivot_wider(names_from = Cell_Type, values_from = Star, values_fill = "") %>%
  column_to_rownames("Gene") %>%
  as.matrix()

# Sorting
get_max_col_index <- function(gene_vals, col_names, target_order) {
  max_col <- col_names[which.max(abs(gene_vals))]
  match(max_col, target_order)
}
sort_df <- data.frame(Gene = rownames(mat_l2fc))
sort_df$Max_Index <- apply(mat_l2fc, 1, function(x) get_max_col_index(x, colnames(mat_l2fc), target_cols))
sort_df$Max_Val <- apply(mat_l2fc, 1, max)
sort_df <- sort_df %>% arrange(Max_Index, desc(Max_Val))

mat_l2fc_sorted <- mat_l2fc[sort_df$Gene, , drop = FALSE]
mat_sig_sorted <- mat_sig[sort_df$Gene, , drop = FALSE]

valid_cols <- intersect(target_cols, colnames(mat_l2fc_sorted))
mat_l2fc_sorted <- mat_l2fc_sorted[, valid_cols, drop = FALSE]
mat_sig_sorted <- mat_sig_sorted[, valid_cols, drop = FALSE]

# Plot
out_pdf <- file.path(plot_dir, "Fos_Partners_Global_vs_Exclusive_Diff_Stats.pdf")
pdf(out_pdf, width = 9, height = 12)

col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
ht <- Heatmap(mat_l2fc_sorted,
              name = "Log2FC", 
              heatmap_legend_param = list(title = "Log2FC\n(Exclusive/Global)"),
              col = col_fun,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              column_order = valid_cols,
              column_title = "Differential Co-binding (Top 50 Extreme Changes)",
              row_names_gp = gpar(fontsize = 9),
              column_names_gp = gpar(fontsize = 10, fontface = "bold"),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(mat_sig_sorted[i, j], x, y, gp = gpar(fontsize = 8))
              }
)
draw(ht)
dev.off()

library(GenomicRanges)
library(ComplexHeatmap)
library(UpSetR)
library(stringr)
library(dplyr)

# --- Configuration ---
indir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Classification")
out_dir <- file.path(indir, "Plots_UpSet")
if (!dir.exists(out_dir)) dir.create(out_dir)

# --- 1. Define the UpSet Function (Any Condition / Union) ---
create_union_upset <- function(gr_list, motif_name) {
  
  # Remove empty cell types
  non_empty_indices <- sapply(gr_list, function(gr) length(gr) > 0)
  
  if (all(!non_empty_indices)) {
    message(paste("Skipping", motif_name, "- No peaks found."))
    return(NULL)
  }
  
  gr_list <- gr_list[non_empty_indices]
  
  if (length(gr_list) < 2) {
    message("Fewer than 2 cells have data. Skipping UpSet.")
    return(NULL)
  }
  
  message(paste("Generating Peak-Based Union UpSet for:", motif_name))
  
  # Step 1: Consensus Peaks (Using PEAK coordinates)
  # reduce() merges overlapping regions across cell types into a single consensus set
  all_ranges <- unlist(GRangesList(gr_list))
  consensus_peaks <- GenomicRanges::reduce(all_ranges)
  
  # Step 2: Build Binary Matrix
  mat <- matrix(0, nrow = length(consensus_peaks), ncol = length(gr_list))
  colnames(mat) <- names(gr_list)
  
  for (cell in names(gr_list)) {
    hits <- findOverlaps(consensus_peaks, gr_list[[cell]])
    mat[queryHits(hits), cell] <- 1
  }
  
  upset_df <- as.data.frame(mat)
  
  # Step 3: Plot (Black/Neutral color for "Any Condition")
  p <- upset(upset_df, 
             sets = colnames(upset_df), 
             order.by = "freq", 
             empty.intersections = "on",
             nsets = length(colnames(upset_df)),
             main.bar.color = "black",
             sets.bar.color = "black", 
             text.scale = 1.5
  )
  
  return(p)
}

# --- 2. Main Execution (Fos Only) ---

all_files <- list.files(indir, pattern = "_containingPeaks_Classified.rds$", full.names = TRUE)
fos_file <- all_files[grepl("/Fos_", all_files) | grepl("^Fos_", basename(all_files))]

if (length(fos_file) > 0) {
  fos_file <- fos_file[1]
  motif_name <- basename(fos_file) %>% str_remove("_containingPeaks_Classified.rds")
  
  df_all <- readRDS(fos_file)
  
  if (!is.null(df_all) && nrow(df_all) > 0) {
    
    df_list <- split(df_all, df_all$cell_type)
    
    # --- CRITICAL FIX APPLIED HERE ---
    gr_list <- lapply(df_list, function(df) {
      
      # FIX: Use dplyr::rename explicitly
      df_safe <- df %>%
        dplyr::rename(
          motif_start = start,
          motif_end = end,
          motif_seqnames = seqnames
       
        )
      
      # Create GRanges using PEAK coordinates
      makeGRangesFromDataFrame(df_safe, 
                               keep.extra.columns = TRUE, 
                               seqnames.field = "peak_chr",   
                               start.field    = "peak_start", 
                               end.field      = "peak_end"    
      )
    })
    
    # Cleanup unwanted cell types
    gr_list$M1 <- NULL
    gr_list$M2 <- NULL
    gr_list$LSC <- NULL
    
    # Run for "Any Condition" (Union)
    upset_plot <- create_union_upset(gr_list, motif_name)
    
    if (!is.null(upset_plot)) {
      pdf_file <- file.path(out_dir, paste0(motif_name, "_UpSet_AnyCondition_PEAKS.pdf"))
      
      pdf(pdf_file, width = 10, height = 6, onefile = FALSE)
      print(upset_plot)
      dev.off()
      message(paste("Saved:", pdf_file))
    }
  }
}

head(df)

library(GenomicRanges)
library(ComplexHeatmap)
library(UpSetR)
library(stringr)
library(dplyr)

# --- Configuration ---
indir <- publication_path("HEART_SNRNA_ROOT", "analysis/2026-01-19.TOBIAS_Classification")
out_dir <- file.path(indir, "Plots_UpSet")
if (!dir.exists(out_dir)) dir.create(out_dir)

# --- 1. Define the UpSet Function ---
create_condition_upset <- function(gr_list, target_category, motif_name) {
  
  message(paste("Filtering for category:", target_category))
  
  gr_list_filtered <- lapply(gr_list, function(gr) {
    if ("category" %in% colnames(mcols(gr))) {
      return(gr[which(mcols(gr)$category == target_category)])
    } else {
      return(GRanges())
    }
  })
  
  non_empty_indices <- sapply(gr_list_filtered, function(gr) length(gr) > 0)
  
  if (all(!non_empty_indices)) {
    message(paste("Skipping - No peaks found for category:", target_category))
    return(NULL)
  }
  
  gr_list_filtered <- gr_list_filtered[non_empty_indices]
  
  if (length(gr_list_filtered) < 2) {
    message("Fewer than 2 cells have this category. Skipping UpSet.")
    return(NULL)
  }
  
  # Step 2: Consensus Peaks (Using PEAK coordinates)
  all_ranges <- unlist(GRangesList(gr_list_filtered))
  consensus_peaks <- GenomicRanges::reduce(all_ranges)
  
  # Step 3: Binary Matrix
  mat <- matrix(0, nrow = length(consensus_peaks), ncol = length(gr_list_filtered))
  colnames(mat) <- names(gr_list_filtered)
  
  for (cell in names(gr_list_filtered)) {
    hits <- findOverlaps(consensus_peaks, gr_list_filtered[[cell]])
    mat[queryHits(hits), cell] <- 1
  }
  
  upset_df <- as.data.frame(mat)
  
  # Step 4: Plot
  p <- upset(upset_df, 
             sets = colnames(upset_df), 
             order.by = "freq", 
             empty.intersections = "on",
             nsets = length(colnames(upset_df)),
             main.bar.color = "black",
             sets.bar.color = "#E41A1C", # Red for Aged Sedentary
             text.scale = 1.5
  )
  
  return(p)
}

# --- 2. Main Execution (Fos Only) ---

all_files <- list.files(indir, pattern = "_containingPeaks_Classified.rds$", full.names = TRUE)
fos_file <- all_files[grepl("/Fos_", all_files) | grepl("^Fos_", basename(all_files))]

if (length(fos_file) > 0) {
  fos_file <- fos_file[1]
  motif_name <- basename(fos_file) %>% str_remove("_containingPeaks_Classified.rds")
  
  df_all <- readRDS(fos_file)
  
  if (!is.null(df_all) && nrow(df_all) > 0) {
    
    df_list <- split(df_all, df_all$cell_type)
    
    # --- CRITICAL FIX APPLIED HERE ---
    gr_list <- lapply(df_list, function(df) {
      
      # FIX: Use dplyr::rename explicitly to avoid conflict with plyr/reshape
      df_safe <- df %>%
        dplyr::rename(
          motif_start = start,
          motif_end = end,
          motif_seqnames = seqnames,
        )
      
      # Create GRanges using PEAK coordinates
      makeGRangesFromDataFrame(df_safe, 
                               keep.extra.columns = TRUE, 
                               seqnames.field = "peak_chr",   
                               start.field    = "peak_start", 
                               end.field      = "peak_end"    
      )
    })
    
    # Cleanup unwanted cell types
    gr_list$M1 <- NULL
    gr_list$M2 <- NULL
    gr_list$LSC <- NULL
    
    # Run for "In Aged Sed Only" (Flipped/Gained Peaks)
    upset_plot <- create_condition_upset(gr_list, target_category = "In Aged Sed Only", motif_name)
    
    if (!is.null(upset_plot)) {
      pdf_file <- file.path(out_dir, paste0(motif_name, "_UpSet_AgedSedOnly_PEAKS.pdf"))
      
      pdf(pdf_file, width = 10, height = 6, onefile = FALSE)
      print(upset_plot)
      dev.off()
      message(paste("Saved:", pdf_file))
    }
  }
}

# Guarded: the identical rename is already applied earlier in this script (right after the TOBIAS
# tables are bound), and repeating it fails because cell_type / peak_identifier no longer exist.
if ("cell_type" %in% colnames(Tobias_df)) {
  Tobias_df <- Tobias_df %>% dplyr::rename(cell = cell_type, peak_coord_Age = peak_identifier)
}

# ==============================================================================
# 1. Setup and Libraries
# ==============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(GenomicRanges)
  library(rtracklayer)
  library(openxlsx)
  library(EnhancedVolcano)
})

# Define Paths
base_dir   <- publication_path("HEART_SNRNA_ROOT", "analysis")
tobias_dir <- file.path(base_dir, "2026-01-19.TOBIAS_Classification") # NEW PATH
mast_dir   <- file.path(base_dir, "2023-11-06.MAST_Final_Fixed")
cicero_dir <- file.path(base_dir, "2023-10-31.Cicero")
gtf_path   <- publication_path("CELLRANGER_GEX_REFERENCE", "genes/genes.gtf")

# Output for intermediate saves
out_dir <- file.path(cicero_dir, "2026-01-20_Recreated_Links")
if(!dir.exists(out_dir)) dir.create(out_dir)

# ==============================================================================
# 2. Load NEW TOBIAS Data (Classified)
# ==============================================================================
message("Loading new TOBIAS classification files...")

tobias_files <- list.files(tobias_dir, pattern = "_containingPeaks_Classified.rds$", full.names = TRUE)

# Read and combine all TOBIAS files
Tobias_List <- lapply(tobias_files, function(f) {
  # Extract motif name from filename
  # Filename format: "MotifName_..._containingPeaks_Classified.rds"
  motif_name <- basename(f) %>% str_remove("_containingPeaks_Classified.rds")
  
  # Read RDS (This is a single DF or a list of DFs depending on previous steps)
  dat <- readRDS(f)
  
  # Ensure it is a data frame
  if (!is.data.frame(dat)) {
    # If it's a list (split by cell), bind it back
    dat <- bind_rows(dat)
  }
  
  if (nrow(dat) > 0) {
    dat$motif <- motif_name
    return(dat)
  } else {
    return(NULL)
  }
})

Tobias_df <- bind_rows(Tobias_List)

# Clean Cell Types
Tobias_df$cell_type <- ifelse(Tobias_df$cell_type == "M1M2", "AllMacrophage", Tobias_df$cell_type)
Tobias_df$cell_type <- ifelse(Tobias_df$cell_type == "MSC", "HSPC", Tobias_df$cell_type)
# Rename to match downstream logic
Tobias_df <- Tobias_df %>% dplyr::rename(cell = cell_type, peak_coord_Age = peak_identifier)

message(paste("Loaded TOBIAS data. Total rows:", nrow(Tobias_df)))

saveRDS(Tobias_df, file.path(out_dir, "Tobias_df_Classified_Full.rds"))

# ==============================================================================
# 3. Load and Process MAST (DEG) Data
# ==============================================================================
message("Loading MAST DEG results...")

mast_files <- list.files(mast_dir, pattern = 'MAST.*csv', recursive = TRUE, full.names = TRUE)
mast_files <- mast_files[!grepl("Global", mast_files)]

All <- bind_rows(lapply(mast_files, function(f) {
  df <- read.csv(f)
  clean_name <- str_remove_all(basename(f), "MAST|.csv")
  condition  <- str_extract(clean_name, "Aged.*|Young.*")
  cell       <- str_remove(clean_name, "_Aged.*|_Young.*")
  
  df$cell <- cell
  df$condition <- condition
  return(df)
}))

# Format Identifiers
All$state <- paste(All$primerid, All$cell)

# Split and Filter Intersection
Age <- All %>% filter(condition == 'Aged_Sed-Young_Sed')
Ex  <- All %>% filter(condition == 'Aged_Ex-Aged_Sed')
common <- intersect(Age$state, Ex$state)

Age <- Age %>% filter(state %in% common) %>% arrange(state)
Ex  <- Ex %>% filter(state %in% common) %>% arrange(state)

colnames(Ex) <- paste0(colnames(Ex), ".Ex")
colnames(Age) <- paste0(colnames(Age), ".Age")

Age_Ex <- cbind(Ex, Age)
Age_Ex <- Age_Ex[rowSums(is.na(Age_Ex)) != ncol(Age_Ex), ]

# Define Directions
FC <- 0.5
sig_age <- Age_Ex$fdr.Age < 0.05
sig_ex  <- Age_Ex$fdr.Ex < 0.05
fc_age  <- Age_Ex$coef.Age
fc_ex   <- Age_Ex$coef.Ex

Age_Ex$Direction <- case_when(
  sig_age & abs(fc_age) > FC & (!sig_ex | abs(fc_ex) < FC) ~ "Age DEG",
  sig_ex  & abs(fc_ex) > FC  & (!sig_age | abs(fc_age) < FC) ~ "Ex DEG",
  sig_age & sig_ex & ((fc_age > FC & fc_ex > FC) | (fc_age < -FC & fc_ex < -FC)) ~ "Exercise exacerbates",
  sig_age & sig_ex & ((fc_age < -FC & fc_ex > FC) | (fc_age > FC & fc_ex < -FC)) ~ "Exercise reverses",
  TRUE ~ "NSig"
)

message("DEG Processing Complete.")

# ==============================================================================
# 4. Process Cicero Connections
# ==============================================================================
message("Loading Cicero connections...")

conns <- readRDS(file.path(cicero_dir, "conns.rds"))
IntConns <- subset(conns, coaccess != 0)

# Add Self-Loops (Peak1 == Peak2)
all_peaks <- unique(IntConns$Peak1)
self_loops <- data.frame(Peak1 = all_peaks, Peak2 = all_peaks, coaccess = 1)
IntConns <- rbind(IntConns, self_loops)
IntConns$intcons_index <- 1:nrow(IntConns)

# Helper: Parse "chr-start-end"
parse_peaks <- function(peak_strings) {
  peak_strings <- as.character(peak_strings)
  # Robust regex for chr-start-end
  # Matches "chr1-100-200"
  parts <- str_match(peak_strings, "^(chr[^:-]+)[-:](\\d+)[-](\\d+)$")
  
  GRanges(seqnames = parts[,2],
          ranges = IRanges(start = as.numeric(parts[,3]), end = as.numeric(parts[,4])))
}

# Load GTF for Promoters
gtf <- import(gtf_path)
gtf_genes <- gtf[gtf$type == "gene"]
promoters_gr <- promoters(gtf_genes, upstream = 1000, downstream = 1000)
mcols(promoters_gr) <- mcols(promoters_gr)[, c("gene_name", "gene_id")]

message("Identifying Promoter-Anchored Interactions...")

# Convert Cicero Peaks to GRanges
p1_gr <- parse_peaks(IntConns$Peak1)
p2_gr <- parse_peaks(IntConns$Peak2)

# Overlap Peak 1 with Promoters
ov1 <- findOverlaps(p1_gr, promoters_gr)
df1 <- data.frame(
  intcons_index = IntConns$intcons_index[queryHits(ov1)],
  gene_id = promoters_gr$gene_name[subjectHits(ov1)],
  candidate_peak = IntConns$Peak2[queryHits(ov1)], # The OTHER end is the candidate enhancer
  stringsAsFactors = FALSE
)

# Overlap Peak 2 with Promoters
ov2 <- findOverlaps(p2_gr, promoters_gr)
df2 <- data.frame(
  intcons_index = IntConns$intcons_index[queryHits(ov2)],
  gene_id = promoters_gr$gene_name[subjectHits(ov2)],
  candidate_peak = IntConns$Peak1[queryHits(ov2)], # The OTHER end is the candidate enhancer
  stringsAsFactors = FALSE
)

# Combine -> This DF contains all interactions where at least one end is a promoter
promoter_intcons_df <- bind_rows(df1, df2) %>% distinct()

# Pre-parse candidate peaks (the potential enhancers)
candidate_peaks_gr <- parse_peaks(promoter_intcons_df$candidate_peak)
promoter_intcons_df$cand_chr <- as.character(seqnames(candidate_peaks_gr))
promoter_intcons_df$cand_start <- start(candidate_peaks_gr)
promoter_intcons_df$cand_end <- end(candidate_peaks_gr)

head(promoter_intcons_df)

# ==============================================================================
# 5. Map TOBIAS Peaks to Cicero Links (Per Cell Type)
# ==============================================================================
cells <- intersect(unique(Tobias_df$cell), unique(Age_Ex$cell.Ex))
result_list <- list()
message(cells)
for (ct in cells) {
  message(paste("Linking data for cell type:", ct))
  
  # 1. Filter Data for this Cell
  tobias_sub <- Tobias_df %>% filter(cell == ct)
  deg_sub    <- Age_Ex %>% filter(cell.Ex == ct)
  
  if(nrow(tobias_sub) == 0) next
  
  # 2. Convert TOBIAS peaks to GRanges
  # Tobias_df uses "peak_coord_Age" (which is the peak_identifier from the classification file)
  tob_gr <- parse_peaks(tobias_sub$peak_coord_Age)
  
  # 3. Create GRanges for the Candidate Enhancers in Cicero
  # We already parsed these into candidate_peaks_gr, but we need to subset corresponding to promoter_intcons_df
  # For efficiency, we just use the pre-calculated GRanges
  # But we need to overlap 'candidate_peaks_gr' with 'tob_gr'
  
  # findOverlaps: Query = Cicero Candidates, Subject = TOBIAS Peaks
  hits <- findOverlaps(candidate_peaks_gr, tob_gr)
  
  if (length(hits) == 0) {
    message("  No overlaps found.")
    next
  }
  
  # 4. Build the Linked Data Frame
  # Start with the Cicero Promoter Links
  linked_df <- promoter_intcons_df[queryHits(hits), ]
  
  # Add the TOBIAS info (Subject Hits)
  tob_matches <- tobias_sub[subjectHits(hits), ]
  
  # Combine columns
  # We use cbind because the indexing is perfectly aligned by the Hits object
  combined_df <- cbind(linked_df, tob_matches)
  
  # 5. Merge with DEG Data
  # Join on Gene ID
  final_df <- combined_df %>%
    inner_join(deg_sub, by = c("gene_id" = "primerid.Age")) %>%
    # Add Co-access Score
    left_join(IntConns %>% select(intcons_index, coaccess), by = "intcons_index")
  
  # Save intermediate
  if (nrow(final_df) > 0) {
    saveRDS(final_df, file.path(out_dir, paste0(ct, "_MappedTogether_Classified.rds")))
    result_list[[ct]] <- final_df
  }
}

# ==============================================================================
# 6. Combine and Save Final Result
# ==============================================================================
if (length(result_list) > 0) {
  All_final_linked_Tobias_df <- bind_rows(result_list)
  
  save_path <- file.path(cicero_dir, "All_TOBIAS_LINKS_Classified.rds")
  saveRDS(All_final_linked_Tobias_df, save_path)
  message(paste("Saved Final Recreated Links to:", save_path))
  
} else {
  message("No linked data found.")
}


# Define Paths
base_dir   <- publication_path("HEART_SNRNA_ROOT", "analysis")
tobias_dir <- file.path(base_dir, "2026-01-19.TOBIAS_Classification") # NEW PATH
mast_dir   <- file.path(base_dir, "2023-11-06.MAST_Final_Fixed")
cicero_dir <- file.path(base_dir, "2023-10-31.Cicero")
gtf_path   <- publication_path("CELLRANGER_GEX_REFERENCE", "genes/genes.gtf")

# Output for intermediate saves
out_dir <- file.path(cicero_dir, "2026-01-20_Recreated_Links")
if(!dir.exists(out_dir)) dir.create(out_dir)
out_dir





# ==============================================================================
# 1. Setup and Libraries
# ==============================================================================
suppressPackageStartupMessages({
  library(rtracklayer)
  library(openxlsx)
  library(stringr)
  library(EnhancedVolcano)
  library(tidyverse) # Includes dplyr, ggplot2, tidyr, readr
  library(GenomicRanges)
  library(rGREAT)
})

# Define a default theme for plots
theme_mytheme <- theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white"))

# Paths
base_dir   <- publication_path("HEART_SNRNA_ROOT", "analysis")
mast_dir   <- file.path(base_dir, "2023-11-06.MAST_Final_Fixed")
cicero_dir <- file.path(base_dir, "2023-10-31.Cicero")
cisbp_path <- file.path(base_dir, "2024-05-24.ChromBPMotif/moodsUnion_Cisbp/EnhancerOrder/genes_expressed_by_celltype.csv")
gtf_path   <- publication_path("CELLRANGER_GEX_REFERENCE", "genes/genes.gtf")

# ==============================================================================
# 2. Load Gene Annotations and Expression Data
# ==============================================================================
message("Loading gene annotations...")
expressed_genes <- read.csv(cisbp_path)

# Load GTF and define TSS
gene_coordinates <- import(gtf_path)
genes <- gene_coordinates[gene_coordinates$type == 'gene']

tss <- GRanges(
  seqnames = seqnames(genes),
  ranges = IRanges(
    start = ifelse(strand(genes) == "+", start(genes), end(genes)),
    end   = ifelse(strand(genes) == "+", start(genes), end(genes))
  ),
  strand = strand(genes),
  gene_id = mcols(genes)$gene_id,
  gene_name = mcols(genes)$gene_name
)

# ==============================================================================
# 3. Load and Process MAST Results (DEGs)
# ==============================================================================
message("Loading MAST results...")
mast_files <- list.files(mast_dir, pattern = 'MAST.*csv', recursive = TRUE, full.names = TRUE)
mast_files <- mast_files[!grepl("Global", mast_files)]

All <- bind_rows(lapply(mast_files, function(f) {
  df <- read.csv(f)
  clean_name <- str_remove_all(basename(f), "MAST|.csv")
  condition  <- str_extract(clean_name, "Aged.*|Young.*")
  cell       <- str_remove(clean_name, "_Aged.*|_Young.*")
  
  df$cell <- cell
  df$condition <- condition
  return(df)
}))

All$cell_cond <- paste(All$cell, All$condition)
All$state <- paste(All$primerid, All$cell)

Age <- All %>% filter(condition == 'Aged_Sed-Young_Sed')
Ex  <- All %>% filter(condition == 'Aged_Ex-Aged_Sed')
common <- intersect(Age$state, Ex$state)

Age <- Age %>% filter(state %in% common) %>% arrange(state)
Ex  <- Ex %>% filter(state %in% common) %>% arrange(state)

colnames(Ex) <- paste0(colnames(Ex), ".Ex")
colnames(Age) <- paste0(colnames(Age), ".Age")

Age_Ex <- cbind(Ex, Age)
Age_Ex <- Age_Ex[rowSums(is.na(Age_Ex)) != ncol(Age_Ex), ]

# ==============================================================================
# 4. Categorize DEGs
# ==============================================================================
message("Categorizing DEGs...")
FC <- 0.5
sig_age <- Age_Ex$fdr.Age < 0.05
sig_ex  <- Age_Ex$fdr.Ex < 0.05
fc_age  <- Age_Ex$coef.Age
fc_ex   <- Age_Ex$coef.Ex

Age_Ex$Direction <- case_when(
  sig_age & fc_age > FC  & (!sig_ex | abs(fc_ex) < FC) ~ "Age Up",
  sig_age & fc_age < -FC & (!sig_ex | abs(fc_ex) < FC) ~ "Age Down",
  sig_ex  & fc_ex > FC   & (!sig_age | abs(fc_age) < FC) ~ "Ex Up",
  sig_ex  & fc_ex < -FC  & (!sig_age | abs(fc_age) < FC) ~ "Ex Down",
  sig_age & sig_ex & fc_age > FC  & fc_ex > FC  ~ "Pro-Age AgeUp",
  sig_age & sig_ex & fc_age < -FC & fc_ex < -FC ~ "Pro-Age AgeDown",
  sig_age & sig_ex & fc_age < -FC & fc_ex > FC  ~ "Anti-Age AgeDown",
  sig_age & sig_ex & fc_age > FC  & fc_ex < -FC ~ "Anti-Age AgeUp",
  TRUE ~ "NSig"
)

# Sub-categories
Age_Ex$AgeDirection <- "NSig"
Age_Ex$AgeDirection[sig_age & fc_age > FC]  <- "Age DEG Up"
Age_Ex$AgeDirection[sig_age & fc_age < -FC] <- "Age DEG Down"

Age_Ex$ExDirection <- "NSig"
Age_Ex$ExDirection[sig_ex & fc_ex > FC]  <- "Ex DEG Up"
Age_Ex$ExDirection[sig_ex & fc_ex < -FC] <- "Ex DEG Down"

print(table(Age_Ex$Direction))

# ==============================================================================
# 5. Helper Function: Annotate Nearest Genes
# ==============================================================================
annotate_regions_both_sides_Mast <- function(celltype, gr, expressed_genes, gene_gr) {
  cell_lookup <- case_when(
    celltype == "MSC" ~ "HSPC",
    celltype == "LSC" ~ "M_Ex",
    celltype == "Macrophage" ~ "M1M2",
    TRUE ~ celltype
  )
  genes_of_interest <- expressed_genes %>% filter(CellType == cell_lookup) %>% pull(Gene)
  expressed_gene_gr <- gene_gr[gene_gr$gene_name %in% genes_of_interest]
  
  nearest_left <- precede(gr, expressed_gene_gr, ignore.strand = TRUE)
  valid_left <- !is.na(nearest_left)
  mcols(gr)$distance_to_left_gene <- NA
  if(any(valid_left)){
    mcols(gr)$distance_to_left_gene[valid_left] <- abs(start(gr[valid_left]) - end(expressed_gene_gr[nearest_left[valid_left]]))
  }
  
  nearest_right <- follow(gr, expressed_gene_gr, ignore.strand = TRUE)
  valid_right <- !is.na(nearest_right)
  mcols(gr)$distance_to_right_gene <- NA
  if(any(valid_right)){
    mcols(gr)$distance_to_right_gene[valid_right] <- abs(start(expressed_gene_gr[nearest_right[valid_right]]) - end(gr[valid_right]))
  }
  return(gr)
}

# ==============================================================================
# 6. Load NEW Recreated Links (Classified)
# ==============================================================================
message("Loading NEW Classified Links...")

# Path to the file we just created
links_file <- file.path(cicero_dir, "All_TOBIAS_LINKS_Classified.rds")

if (!file.exists(links_file)) {
  stop("New links file not found! Did you run the recreation script?")
}

AllMerged_df <- readRDS(links_file)

# Ensure we are only looking at Fos (should already be filtered, but double check)
    AllMerged_df$motif <- AllMerged_df$motif %>%
  str_remove("^\\(") %>%        # Step 1: Remove leading (
  str_remove("(_|\\().*") %>%   # Step 2: Cut off at metadata
  str_remove("\\)$")    
AllMerged_df <- AllMerged_df %>% filter(motif == "Fos")

# Add Promoter Annotations if missing (though they were added in creation)
AllMerged_df$Promoter_Annotation <- "Distal_link"
message(paste("Loaded Links with Categories:", paste(unique(AllMerged_df$category), collapse=", ")))



table(AllMerged_df$category)



# ==============================================================================
# 7. Enrichment Analysis Function
# ==============================================================================
run_enrichment_test <- function(links_df, deg_df, category_col, category_val) {
  
  # Set interest flag in DEG table
  deg_local <- deg_df
  deg_local$interest <- ifelse(deg_local[[category_col]] == category_val, "DEG", "Nsig")
  
  # Identify shared cell column
  cell_col <- if ("cell.Age" %in% names(links_df)) "cell.Age" else "cell"
  
  # Function to compute stats for one cell
  compute_cell_stats <- function(cell_check) {
    # 1. Genes linked to FOS in this cell
    linked_genes <- links_df %>%
      filter(.data[[cell_col]] == cell_check) %>%
      pull(gene_id) %>% unique()
    
    # 2. Prepare contingency table
    ax <- deg_local %>%
      filter(.data[[cell_col]] == cell_check) %>%
      transmute(
        gene_id = primerid.Age,
        DEG_Status = ifelse(interest == "DEG", "DEG", "Other"),
        Link_Status = ifelse(gene_id %in% linked_genes, "Linked", "Unlinked")
      )
    
    if(nrow(ax) == 0) return(NULL)
    
    # 3. Counts
    a <- sum(ax$DEG_Status == "DEG"   & ax$Link_Status == "Linked")
    b <- sum(ax$DEG_Status == "DEG"   & ax$Link_Status == "Unlinked")
    c <- sum(ax$DEG_Status == "Other" & ax$Link_Status == "Linked")
    d <- sum(ax$DEG_Status == "Other" & ax$Link_Status == "Unlinked")
    
    # 4. Correction & Stats
    add05 <- any(c(a,b,c,d) == 0)
    if (add05) { a <- a + 0.5; b <- b + 0.5; c <- c + 0.5; d <- d + 0.5 }
    
    or    <- (a*d)/(b*c)
    logor <- log(or)
    se    <- sqrt(1/a + 1/b + 1/c + 1/d)
    z     <- logor / se
    p_two <- 2 * pnorm(-abs(z))
    
    data.frame(
      cell = cell_check,
      OR = or,
      p_two_sided = p_two,
      a_DEG_Linked = a, b_DEG_Unlinked = b, c_Other_Linked = c, d_Other_Unlinked = d,
      stringsAsFactors = FALSE
    )
  }
  
  cells <- sort(unique(links_df[[cell_col]]))
  res_df <- bind_rows(lapply(cells, compute_cell_stats))
  
  if(nrow(res_df) > 0) {
    res_df$DEG_Class <- category_val
    res_df$cell <- ifelse(res_df$cell == "AllMacrophage", "Macrophage", res_df$cell)
  }
  return(res_df)
}

# ==============================================================================
# 8. Main Execution Loop
# ==============================================================================

# Iterate through NEW FOS categories (e.g., "In Aged Sed Only", "Global_MAS")
#for (Fos_type in unique(AllMerged_df$category)) {
 Fos_type = "In Aged Sed Only"
  message(paste("Processing FOS type:", Fos_type))
  
  # Filter Links for this category
  FOS_Links_Sub <- AllMerged_df %>%
    filter(abs(coaccess) > 0.2, Promoter_Annotation != "Promoter_promoter", category == Fos_type)
  
  if (nrow(FOS_Links_Sub) == 0) {
    message("  No links found for this category.")
    next
  }

  res_list <- list()
  
  # 1. Main Direction Categories
  cats_main <- unique(Age_Ex$Direction)
  res_list[["main"]] <- bind_rows(lapply(cats_main, function(cat) {
    run_enrichment_test(FOS_Links_Sub, Age_Ex, "Direction", cat)
  }))
  
  # 2. Age Specific Categories
  cats_age <- unique(Age_Ex$AgeDirection)
  res_list[["age"]] <- bind_rows(lapply(cats_age, function(cat) {
    run_enrichment_test(FOS_Links_Sub, Age_Ex, "AgeDirection", cat)
  }))
  cats_ex <- unique(Age_Ex$ExDirection)
  res_list[["ex"]] <- bind_rows(lapply(cats_ex, function(cat) {
    run_enrichment_test(FOS_Links_Sub, Age_Ex, "ExDirection", cat)
  }))
  # Combine
  All_Results <- bind_rows(res_list)
  if(nrow(All_Results) == 0) next
  
  # Prepare Plotting Data
  forest_df2 <- All_Results %>%
    mutate(
      p_bin = cut(p_two_sided, breaks = c(-Inf, 1e-3, 1e-2, 5e-2, Inf), labels = c("<0.001", "<0.01", "<0.05", ">0.05"), right = TRUE),
      logOR = log2(OR)
    )

  # Palette
  pal <- c("<0.001" = "#A81B1B", "<0.01" = "#F99349", "<0.05" = "#3A85A8", ">0.05" = "#808080")
  
  # Plot
  lim <- max(abs(forest_df2$logOR[is.finite(forest_df2$logOR)]), na.rm = TRUE)
  if(is.infinite(lim)) lim <- 4
  
  p <- ggplot(forest_df2, aes(y = cell, x = logOR)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_point(aes(color = p_bin), size = 3) + 
    facet_wrap(~DEG_Class, scales = "free_x", ncol = 2) + 
    xlim(-lim, lim) +
    scale_color_manual(values = pal) +
    labs(title = paste0("Enrichment: ", Fos_type), x = "Log2 Odds Ratio") + 
    theme_mytheme
    print(p)
  # Save
  ggsave(file.path(tobias_dir, paste0(Fos_type, "_DEG_enrichment_AP1.pdf")), plot = p, width = 10, height = 20)
  write.csv(forest_df2, file.path(tobias_dir, paste0(Fos_type, "_DEG_enrichment_AP1.csv")))
#}

library(ggplot2)
library(dplyr)

# 1. Filter for the specific classes of interest
# Note: We use the names we standardized in the previous script 
# (Age DEG Up -> "Age Up", Ex DEG Up -> "Exercise Up", etc.)
target_classes <- c("Age Up", "Age Down", "Ex Up", "Ex Down")

heatmap_data <- forest_df2 %>%
  filter(DEG_Class %in% target_classes) %>%
  mutate(DEG_Class = factor(DEG_Class, levels = target_classes)) # Enforce order

# 2. Plot
p_heat <- ggplot(heatmap_data, aes(x = DEG_Class, y = cell, fill = logOR)) +
  geom_tile(color = "white", lwd = 0.5) + # White borders
  
  # Add stars for significance
  geom_text(aes(label = ifelse(p_two_sided < 0.001, "***", 
                        ifelse(p_two_sided < 0.01, "**",
                        ifelse(p_two_sided < 0.05, "*", "")))), 
            vjust = 0.75, size = 5, color = "black") +
  
  # Color scale: Blue (Depleted) -> White -> Red (Enriched)
  scale_fill_gradient2(
    low = "#3A85A8", mid = "white", high = "#A81B1B", 
    midpoint = 0, name = "Log2(OR)"
  ) +
  
  # Clean layout
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title = element_blank()
  ) +
  labs(title = paste0("Enrichment of DEGs in ", Fos_type, " FOS Peaks"))

# 3. Save
print(p_heat)
ggsave(file.path(cicero_dir, paste0(Fos_type, "_Heatmap_Enrichment_MainCats.pdf")), 
       plot = p_heat, width = 6, height = 8)



# ==============================================================================
# 7. Enrichment Analysis Function
# ==============================================================================
run_enrichment_test <- function(links_df, deg_df, category_col, category_val) {
  
  # Set interest flag in DEG table
  deg_local <- deg_df
  deg_local$interest <- ifelse(deg_local[[category_col]] == category_val, "DEG", "Nsig")
  
  # Identify shared cell column
  cell_col <- if ("cell.Age" %in% names(links_df)) "cell.Age" else "cell"
  
  # Function to compute stats for one cell
  compute_cell_stats <- function(cell_check) {
    # 1. Genes linked to FOS in this cell
    linked_genes <- links_df %>%
      filter(.data[[cell_col]] == cell_check) %>%
      pull(gene_id) %>% unique()
    
    # 2. Prepare contingency table
    ax <- deg_local %>%
      filter(.data[[cell_col]] == cell_check) %>%
      transmute(
        gene_id = primerid.Age,
        DEG_Status = ifelse(interest == "DEG", "DEG", "Other"),
        Link_Status = ifelse(gene_id %in% linked_genes, "Linked", "Unlinked")
      )
    
    if(nrow(ax) == 0) return(NULL)
    
    # 3. Counts
    a <- sum(ax$DEG_Status == "DEG"   & ax$Link_Status == "Linked")
    b <- sum(ax$DEG_Status == "DEG"   & ax$Link_Status == "Unlinked")
    c <- sum(ax$DEG_Status == "Other" & ax$Link_Status == "Linked")
    d <- sum(ax$DEG_Status == "Other" & ax$Link_Status == "Unlinked")
    
    # 4. Correction & Stats
    add05 <- any(c(a,b,c,d) == 0)
    if (add05) { a <- a + 0.5; b <- b + 0.5; c <- c + 0.5; d <- d + 0.5 }
    
    or    <- (a*d)/(b*c)
    logor <- log(or)
    se    <- sqrt(1/a + 1/b + 1/c + 1/d)
    z     <- logor / se
    p_two <- 2 * pnorm(-abs(z))
    
    data.frame(
      cell = cell_check,
      OR = or,
      p_two_sided = p_two,
      a_DEG_Linked = a, b_DEG_Unlinked = b, c_Other_Linked = c, d_Other_Unlinked = d,
      stringsAsFactors = FALSE
    )
  }
  
  cells <- sort(unique(links_df[[cell_col]]))
  res_df <- bind_rows(lapply(cells, compute_cell_stats))
  
  if(nrow(res_df) > 0) {
    res_df$DEG_Class <- category_val
    res_df$cell <- ifelse(res_df$cell == "AllMacrophage", "Macrophage", res_df$cell)
  }
  return(res_df)
}

# ==============================================================================
# 8. Main Execution Loop
# ==============================================================================

# Iterate through NEW FOS categories (e.g., "In Aged Sed Only", "Global_MAS")
#for (Fos_type in unique(AllMerged_df$category)) {
 Fos_type = "All Fos Peaks"
  message(paste("Processing FOS type:", Fos_type))
  
  # Filter Links for this category
  FOS_Links_Sub <- AllMerged_df %>%
    filter(abs(coaccess) > 0.2, Promoter_Annotation != "Promoter_promoter")
  
  if (nrow(FOS_Links_Sub) == 0) {
    message("  No links found for this category.")
    next
  }

  res_list <- list()
  
  # 1. Main Direction Categories
  cats_main <- unique(Age_Ex$Direction)
  res_list[["main"]] <- bind_rows(lapply(cats_main, function(cat) {
    run_enrichment_test(FOS_Links_Sub, Age_Ex, "Direction", cat)
  }))
  
  # 2. Age Specific Categories
  cats_age <- unique(Age_Ex$AgeDirection)
  res_list[["age"]] <- bind_rows(lapply(cats_age, function(cat) {
    run_enrichment_test(FOS_Links_Sub, Age_Ex, "AgeDirection", cat)
  }))
  cats_ex <- unique(Age_Ex$ExDirection)
  res_list[["ex"]] <- bind_rows(lapply(cats_ex, function(cat) {
    run_enrichment_test(FOS_Links_Sub, Age_Ex, "ExDirection", cat)
  }))
  # Combine
  All_Results <- bind_rows(res_list)
  if(nrow(All_Results) == 0) next
  
  # Prepare Plotting Data
  forest_df2 <- All_Results %>%
    mutate(
      p_bin = cut(p_two_sided, breaks = c(-Inf, 1e-3, 1e-2, 5e-2, Inf), labels = c("<0.001", "<0.01", "<0.05", ">0.05"), right = TRUE),
      logOR = log2(OR)
    )

  # Palette
  pal <- c("<0.001" = "#A81B1B", "<0.01" = "#F99349", "<0.05" = "#3A85A8", ">0.05" = "#808080")
  
  # Plot
  lim <- max(abs(forest_df2$logOR[is.finite(forest_df2$logOR)]), na.rm = TRUE)
  if(is.infinite(lim)) lim <- 4
  
  p <- ggplot(forest_df2, aes(y = cell, x = logOR)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_point(aes(color = p_bin), size = 3) + 
    facet_wrap(~DEG_Class, scales = "free_x", ncol = 2) + 
    xlim(-lim, lim) +
    scale_color_manual(values = pal) +
    labs(title = paste0("Enrichment: ", Fos_type), x = "Log2 Odds Ratio") + 
    theme_mytheme
    print(p)
  # Save
  ggsave(file.path(tobias_dir, paste0(Fos_type, "_DEG_enrichment_AP1.pdf")), plot = p, width = 10, height = 20)
  write.csv(forest_df2, file.path(tobias_dir, paste0(Fos_type, "_DEG_enrichment_AP1.csv")))
#}

library(ggplot2)
library(dplyr)

# 1. Filter for the specific classes of interest
# Note: We use the names we standardized in the previous script 
# (Age DEG Up -> "Age Up", Ex DEG Up -> "Exercise Up", etc.)
target_classes <- c("Age Up", "Age Down", "Ex Up", "Ex Down")

heatmap_data <- forest_df2 %>%
  filter(DEG_Class %in% target_classes) %>%
  mutate(DEG_Class = factor(DEG_Class, levels = target_classes)) # Enforce order

# 2. Plot
p_heat <- ggplot(heatmap_data, aes(x = DEG_Class, y = cell, fill = logOR)) +
  geom_tile(color = "white", lwd = 0.5) + # White borders
  
  # Add stars for significance
  geom_text(aes(label = ifelse(p_two_sided < 0.001, "***", 
                        ifelse(p_two_sided < 0.01, "**",
                        ifelse(p_two_sided < 0.05, "*", "")))), 
            vjust = 0.75, size = 5, color = "black") +
  
  # Color scale: Blue (Depleted) -> White -> Red (Enriched)
  scale_fill_gradient2(
    low = "#3A85A8", mid = "white", high = "#A81B1B", 
    midpoint = 0, name = "Log2(OR)"
  ) +
  
  # Clean layout
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title = element_blank()
  ) +
  labs(title = paste0("Enrichment of DEGs in ", Fos_type, " FOS Peaks"))

# 3. Save
print(p_heat)
ggsave(file.path(cicero_dir, paste0(Fos_type, "_Heatmap_Enrichment_AllAP1.pdf")), 
       plot = p_heat, width = 6, height = 8)



# ==============================================================================
# 1. Setup and Libraries
# ==============================================================================
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(BSgenome.Mmusculus.UCSC.mm10)
})

# Define Paths
base_dir   <- publication_path("HEART_SNRNA_ROOT", "analysis")
tobias_dir <- file.path(base_dir, "2026-01-19.TOBIAS_Classification") 
hic_dir    <- file.path(base_dir, "2023-10-31.Cicero")
out_dir    <- file.path(tobias_dir, "HiC_Integration")

if(!dir.exists(out_dir)) dir.create(out_dir)

# ==============================================================================
# 2. Recreate Hi-C Windows (1Mb Window, 100kb Step)
# ==============================================================================
message("Loading Hi-C Data and Creating Windows...")

# Load Liftedover Hi-C Links
df_hic <- readRDS(file.path(hic_dir, "Liftedover.rds"))

# Create GRanges for Hi-C Anchors
Gr_Peak1 <- GRanges(seqnames = df_hic$peak1_seqnames,
                    ranges = IRanges(start = df_hic$peak1_start, end = df_hic$peak1_end))
Gr_Peak2 <- GRanges(seqnames = df_hic$peak2_seqnames,
                    ranges = IRanges(start = df_hic$peak2_start, end = df_hic$peak2_end))

# Assign Link IDs
Gr_Peak1$link_id <- 1:nrow(df_hic)
Gr_Peak2$link_id <- 1:nrow(df_hic)
HICs <- c(Gr_Peak1, Gr_Peak2)
seqlevelsStyle(HICs) <- "UCSC"

# Create Sliding Windows 
chrom_sizes <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)
chrom_sizes_df <- data.frame(seqnames = names(chrom_sizes), seqlengths = as.numeric(chrom_sizes))
standard_chromosomes <- paste0('chr', c(1:19, 'X', 'Y'))
chrom_sizes_df <- chrom_sizes_df[chrom_sizes_df$seqnames %in% standard_chromosomes, ]

window_size <- 1e6  # 1Mb
step_size   <- 1e5  # 100kb

create_windows <- function(chrom, length) {
  starts <- seq(1, length - window_size + 1, by = step_size)
  ends <- starts + window_size - 1
  GRanges(seqnames = chrom, ranges = IRanges(start = starts, end = ends))
}

gr_windows_list <- lapply(1:nrow(chrom_sizes_df), function(i) {
  create_windows(chrom_sizes_df$seqnames[i], chrom_sizes_df$seqlengths[i])
})
gr_windows <- do.call(c, gr_windows_list)

# Count Links per Window
message("Counting Hi-C links per window (this may take a moment)...")
overlaps <- findOverlaps(gr_windows, HICs)

overlaps_df <- data.frame(
  window_id = queryHits(overlaps),
  link_id = HICs$link_id[subjectHits(overlaps)]
)

counts_df <- overlaps_df %>%
  group_by(window_id) %>%
  summarize(hic_density = n_distinct(link_id))

gr_windows$hic_density <- 0
gr_windows$hic_density[counts_df$window_id] <- counts_df$hic_density

message("Hi-C Window Map Created.")

# ==============================================================================
# 3. Load NEW Classified TOBIAS Peaks
# ==============================================================================
message("Loading Classified TOBIAS Peaks...")

# Find Fos Files
tobias_files <- list.files(tobias_dir, pattern = "_containingPeaks_Classified.rds$", full.names = TRUE)
fos_files <- tobias_files[grepl("Fos", tobias_files)]

# Load and Combine
Tobias_List <- lapply(fos_files, function(f) {
  motif_name <- basename(f) %>% str_remove("_containingPeaks_Classified.rds")
  dat <- readRDS(f)
  
  if (!is.data.frame(dat)) {
    dat <- bind_rows(dat)
  }
  
  if (nrow(dat) > 0) {
    dat$motif <- motif_name
    return(dat)
  } else {
    return(NULL)
  }
})

Tobias_df <- bind_rows(Tobias_List)

# Filter for Fos only and valid categories
Tobias_df <- Tobias_df %>% 
  filter(grepl("Fos", motif)) %>%
  filter(!is.na(category))

# --- FIX: Rename reserved columns before converting to GRanges ---
# Columns like 'start', 'end', 'width', 'strand' cannot be metadata columns in GRanges
Tobias_df_Safe <- Tobias_df %>%
  dplyr::rename(
    motif_seqnames = seqnames,
    motif_start    = start,
    motif_end      = end
  )

# Convert TOBIAS Peaks to GRanges
# We use PEAK coordinates for the range, and keep the rest (now renamed) as metadata
tobias_gr <- makeGRangesFromDataFrame(Tobias_df_Safe, 
                                      keep.extra.columns = TRUE, 
                                      seqnames.field = "peak_chr", 
                                      start.field = "peak_start", 
                                      end.field = "peak_end")
seqlevelsStyle(tobias_gr) <- "UCSC"

# ==============================================================================
# 4. Overlap Peaks with Hi-C Windows
# ==============================================================================
message("Annotating Peaks with Hi-C Density...")

# Find which window each peak falls into
hits <- findOverlaps(tobias_gr, gr_windows)

# Extract density. If a peak overlaps multiple windows (due to sliding), take MEAN.
density_mapping <- data.frame(
  peak_idx = queryHits(hits),
  density  = gr_windows$hic_density[subjectHits(hits)]
) %>%
  group_by(peak_idx) %>%
  summarize(mean_hic_density = mean(density))

# Add to original dataframe
Tobias_df$HiC_Density <- NA
Tobias_df$HiC_Density[density_mapping$peak_idx] <- density_mapping$mean_hic_density

# Filter out peaks that didn't overlap any valid window
Tobias_df_Clean <- Tobias_df %>% filter(!is.na(HiC_Density))

# ==============================================================================
# 5. Analyze & Plot
# ==============================================================================
message("Generating Plots...")

# Define Categories of Interest
target_cats <- c("In Aged Sed Only", "In All Conditions", "Global_MAS", "Global_YS")
Plot_Data <- Tobias_df_Clean %>% 
  filter(category %in% target_cats)

# Plot 1: Boxplot of Link Density by Category
p1 <- ggplot(Plot_Data, aes(x = category, y = HiC_Density, fill = category)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.3) +
  theme_bw() +
  labs(title = "Hi-C Connectivity of Fos Binding Sites",
       subtitle = "Are dynamic peaks in regions with more looping?",
       y = "Hi-C Link Density (Links per 1Mb Window)",
       x = "Peak Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

print(p1)
ggsave(file.path(out_dir, "Fos_HiC_Density_Boxplot.pdf"), plot = p1, width = 6, height = 6)

# Plot 2: Density Curve
p2 <- ggplot(Plot_Data, aes(x = HiC_Density, color = category, fill = category)) +
  geom_density(alpha = 0.3) +
  theme_bw() +
  labs(title = "Distribution of Hi-C Density",
       x = "Hi-C Link Density",
       y = "Frequency")

print(p2)
ggsave(file.path(out_dir, "Fos_HiC_Density_Distribution.pdf"), plot = p2, width = 7, height = 5)

# Save Summary Stats
summary_stats <- Tobias_df_Clean %>%
  group_by(category) %>%
  summarize(
    Mean_Density = mean(HiC_Density),
    Median_Density = median(HiC_Density),
    Peak_Count = n()
  )

write.csv(summary_stats, file.path(out_dir, "HiC_Density_Summary.csv"))
print(summary_stats)

message("Analysis Complete.")

# ==============================================================================
# 1. Setup & Define Hotspots (Poisson Method)
# ==============================================================================
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(BSgenome.Mmusculus.UCSC.mm10)
})

# Define Paths
base_dir   <- publication_path("HEART_SNRNA_ROOT", "analysis")
tobias_dir <- file.path(base_dir, "2026-01-19.TOBIAS_Classification")
hic_dir    <- file.path(base_dir, "2023-10-31.Cicero")
out_dir    <- file.path(tobias_dir, "HiC_Integration")
if(!dir.exists(out_dir)) dir.create(out_dir)

# --- 1A. Recreate Windows (Needed for Poisson Test) ---
# We need the full distribution of densities to calculate the Poisson background
message("Recreating Hi-C windows for statistical testing...")

df_hic <- readRDS(file.path(hic_dir, "Liftedover.rds"))

# Create anchors
Gr_Peak1 <- GRanges(seqnames = df_hic$peak1_seqnames, ranges = IRanges(start = df_hic$peak1_start, end = df_hic$peak1_end))
Gr_Peak2 <- GRanges(seqnames = df_hic$peak2_seqnames, ranges = IRanges(start = df_hic$peak2_start, end = df_hic$peak2_end))
HICs <- c(Gr_Peak1, Gr_Peak2)
HICs$link_id <- rep(1:nrow(df_hic), 2)
seqlevelsStyle(HICs) <- "UCSC"

# Sliding Windows (1Mb, 100kb step)
chrom_sizes <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)
standard_chromosomes <- paste0('chr', c(1:19, 'X', 'Y'))
chrom_sizes <- chrom_sizes[names(chrom_sizes) %in% standard_chromosomes]

create_windows <- function(chrom, length) {
  starts <- seq(1, length - 1e6 + 1, by = 1e5)
  GRanges(seqnames = chrom, ranges = IRanges(start = starts, width = 1e6))
}
gr_windows <- do.call(c, lapply(names(chrom_sizes), function(n) create_windows(n, chrom_sizes[[n]])))

# Count Links
overlaps <- findOverlaps(gr_windows, HICs)
counts <- data.frame(win = queryHits(overlaps), link = HICs$link_id[subjectHits(overlaps)]) %>%
  group_by(win) %>% summarise(n = n_distinct(link))
gr_windows$hic_density <- 0
gr_windows$hic_density[counts$win] <- counts$n

# --- 1B. NEW: Poisson Test for Hotspots ---
message("Calculating Poisson Statistics for Hotspots...")

# 1. Calculate Lambda (Expected background rate)
# We exclude regions with 0 signal if they represent unmappable regions, 
# but for a conservative test, including them lowers lambda (making it easier to be significant).
# Usually, it is safer to use the mean of non-zero regions or the global mean.
# Here we use the global mean to define "Genome-wide Enriched".
lambda_global <- mean(gr_windows$hic_density)

# 2. Poisson Test
# P(X >= k) = P(X > k-1) = ppois(k-1, lambda, lower.tail=FALSE)
p_values <- ppois(gr_windows$hic_density - 1, lambda = lambda_global, lower.tail = FALSE)

# 3. FDR Correction
fdr_values <- p.adjust(p_values, method = "BH")

# 4. Define Hotspots (Significantly Enriched)
# FDR < 0.05 AND Count > Mean (to ensure enrichment, not just deviation)
hotspots_HIC <- gr_windows[fdr_values < 0.05 & gr_windows$hic_density > lambda_global]

message(paste0("Defined ", length(hotspots_HIC), " Hotspots using Poisson Test (FDR < 0.05)."))
message(paste0("Background Lambda: ", round(lambda_global, 2), " links/window."))
message(paste0("Min Density in Hotspots: ", min(hotspots_HIC$hic_density)))

# ==============================================================================
# 2. Prepare TOBIAS Data & Clean Cells
# ==============================================================================
# Load data (assuming Tobias_df_Clean exists from previous steps, or reloading here for safety)
# For this script to be standalone-ish, I'll rely on the classification files
tobias_files <- list.files(tobias_dir, pattern = "_containingPeaks_Classified.rds$", full.names = TRUE)
fos_files <- tobias_files[grepl("Fos", tobias_files)]

Tobias_List <- lapply(fos_files, function(f) {
  dat <- readRDS(f)
  if (!is.data.frame(dat)) dat <- bind_rows(dat)
  if (nrow(dat) > 0) { dat$motif <- basename(f) %>% str_remove("_containingPeaks_Classified.rds"); return(dat) }
})
Age_Ex_ATAC <- bind_rows(Tobias_List) %>% filter(grepl("Fos", motif), !is.na(category))

# Clean Cell Names
Age_Ex_ATAC$cell <- case_when(
  Age_Ex_ATAC$cell_type == "M1M2" ~ "AllMacrophage",
  Age_Ex_ATAC$cell_type == "Macrophage" ~ "AllMacrophage",
  Age_Ex_ATAC$cell_type == "MSC" ~ "HSPC",
  Age_Ex_ATAC$cell_type == "LSC" ~ "MacFib",
  TRUE ~ Age_Ex_ATAC$cell_type
)

# Fix Columns for GRanges
Age_Ex_ATAC <- Age_Ex_ATAC %>% dplyr::rename(peak_chr = peak_chr, peak_start = peak_start, peak_end = peak_end)

# ==============================================================================
# 3. Main Analysis Loop (Fisher Test)
# ==============================================================================
AP_Plot <- list()
store_results <- list()

categories <- unique(Age_Ex_ATAC$category)
categories <- categories[!is.na(categories)]

for (cat in categories) {
  
  message(paste("Analyzing Category:", cat))
  store <- list()
  
  for (cell_extract in unique(Age_Ex_ATAC$cell)) {
    
    # A. Test Set
    test_df <- Age_Ex_ATAC %>% filter(category == cat & cell == cell_extract)
    if (nrow(test_df) < 10) next
    
    Flip <- makeGRangesFromDataFrame(test_df, seqnames.field = "peak_chr", start.field = "peak_start", end.field = "peak_end")
    
    # B. Background Set
    if (cat == 'In All Conditions') {
      bg_df <- Age_Ex_ATAC %>% filter(category != cat & cell == cell_extract)
    } else {
      bg_df <- Age_Ex_ATAC %>% filter(category == 'In All Conditions' & cell == cell_extract)
    }
    if (nrow(bg_df) < 10) next
    
    Background_gr <- makeGRangesFromDataFrame(bg_df, seqnames.field = "peak_chr", start.field = "peak_start", end.field = "peak_end")
    
    # C. Overlaps
    num_overlaps_test <- length(unique(queryHits(findOverlaps(Flip, hotspots_HIC))))
    num_non_overlaps_test <- length(Flip) - num_overlaps_test
    
    num_overlaps_bg <- length(unique(queryHits(findOverlaps(Background_gr, hotspots_HIC))))
    num_non_overlaps_bg <- length(Background_gr) - num_overlaps_bg
    
    # D. Fisher Test
    contingency_table <- matrix(c(num_overlaps_test, num_non_overlaps_test, num_overlaps_bg, num_non_overlaps_bg), nrow = 2, byrow = TRUE)
    store[[cell_extract]] <- list(fisher.test(contingency_table), contingency_table)
  }
  
  # ==========================================================================
  # 4. Plotting
  # ==========================================================================
  if (length(store) == 0) next
  
  results_df <- do.call(rbind, lapply(names(store), function(cell_type) {
    ft <- store[[cell_type]][[1]]
    data.frame(cell = cell_type, OR = log2(ft$estimate), p.value = ft$p.value, stringsAsFactors = FALSE)
  }))
  
  # Cleanup & Limits
  results_df <- results_df %>% filter(!grepl("unk|Cardiac|Smooth|Dendri|NK|M1|M2", cell))
  
  # Handle Inf
  max_val <- max(abs(results_df$OR[is.finite(results_df$OR)]), na.rm=TRUE)
  results_df$OR[is.infinite(results_df$OR) & results_df$OR > 0] <- max_val + 0.5
  results_df$OR[is.infinite(results_df$OR) & results_df$OR < 0] <- -(max_val + 0.5)
  
  results_df <- results_df %>% mutate(p.value_category = case_when(
    p.value < 0.001 ~ "< 0.001", p.value < 0.01 ~ "< 0.01", p.value < 0.05 ~ "< 0.05", TRUE ~ "NS"
  ))
  
  # Theme
  theme_mytheme <- theme_bw() + theme(panel.grid.minor = element_blank(), axis.text = element_text(color="black"))
  
  # Plot
  p <- ggplot(results_df, aes(x = OR, y = cell)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_segment(aes(x = 0, xend = OR, y = cell, yend = cell), color = "grey") +
    geom_point(aes(color = p.value_category), size = 5) +
    scale_color_manual(values = c("< 0.001" = "red", "< 0.01" = "orange", "< 0.05" = "yellow", "NS" = "black")) +
    labs(title = paste0(cat, ": Hi-C Hotspot Enrichment (Poisson)"), x = "Log2 Odds Ratio", y = "Cell Type", color = "P-value") +
    theme_mytheme + coord_cartesian(xlim = c(-1.5, 1.5))
  
  safe_cat <- str_replace_all(cat, " ", "_")
  ggsave(file.path(out_dir, paste0("HiC_Enrichment_Poisson_", safe_cat, ".pdf")), plot = p, width = 7, height = 6)
  AP_Plot[[cat]] <- p
}

message("Poisson-defined Hotspot Analysis Complete.")

# ==============================================================================
# 1. Setup and Inputs
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(ggplot2)
  library(ggpubr)
  library(GenomicRanges)
  library(viridisLite)
  library(stringr)
})

# --- Paths ---
base_dir   <- publication_path("HEART_SNRNA_ROOT")
atac_path  <- file.path(base_dir, "data/2022-10-17.SharedObject/2024-07-16_Full_ATAC_MatchAnnotated.rds")
tobias_dir <- file.path(base_dir, "analysis/2026-01-19.TOBIAS_Classification")

# --- User Controls ---
motif_pattern <- "Fos"          
meta_condition_col <- "age_condition" 
assay_name <- "ATAC"
peak_slot <- "data"             

# ==============================================================================
# 2. Load Data & Clean Labels
# ==============================================================================
message("Loading ATAC-seq object...")
ATAC <- readRDS(atac_path)
DefaultAssay(ATAC) <- assay_name

message("Applying Cell Type updates (LSC->HSPC, M1/M2->Macrophage)...")
# Update ATAC Object
ATAC$CellType <- as.character(ATAC$CellType)
ATAC$CellType[ATAC$CellType %in% c("M1", "M2")] <- "Macrophage"
ATAC$CellType[ATAC$CellType == "LSC"] <- "HSPC"

message("Loading Classified TOBIAS Peaks...")
all_files <- list.files(tobias_dir, pattern = "_containingPeaks_Classified.rds$", full.names = TRUE)
fos_file <- all_files[grepl(paste0("^", motif_pattern), basename(all_files))][1]
if (is.na(fos_file)) stop("Could not find the Fos Classification file!")

tobias_dat <- readRDS(fos_file)
if (!is.data.frame(tobias_dat)) tobias_dat <- bind_rows(tobias_dat)

# Update TOBIAS Data labels
tobias_dat$cell_type[tobias_dat$cell_type %in% c("M1", "M2", "M1M2")] <- "Macrophage"
tobias_dat$cell_type[tobias_dat$cell_type == "LSC"] <- "HSPC"

# ==============================================================================
# 3. Accessibility Analysis Loop
# ==============================================================================
make_gr <- function(df) {
  if (nrow(df) == 0) return(GRanges())
  GRanges(seqnames = df$peak_chr, ranges = IRanges(start = df$peak_start, end = df$peak_end))
}

overlap_to_df <- function(obj, gr_set, label, cell_type) {
  if (length(gr_set) == 0L) return(NULL)
  signac_peaks <- obj@assays[[assay_name]]@ranges
  hits <- findOverlaps(signac_peaks, gr_set, ignore.strand = TRUE)
  if (length(hits) == 0L) return(NULL)
  
  peak_names <- rownames(obj)[queryHits(hits)]
  acc <- GetAssayData(obj, assay = assay_name, slot = peak_slot)
  acc_sub <- acc[intersect(peak_names, rownames(acc)), , drop = FALSE]
  if (nrow(acc_sub) == 0L) return(NULL)
  
  df <- as.data.frame(as.matrix(acc_sub))
  df$peak <- rownames(df)
  long <- pivot_longer(df, -peak, names_to = "cell", values_to = "accessibility")
  
  cond <- obj[[meta_condition_col]][,1]
  names(cond) <- colnames(obj)
  long$condition <- cond[long$cell]
  long$CellType <- cell_type
  
  out <- as.data.table(long)[, .(mean_accessibility = mean(accessibility, na.rm = TRUE)),
                             by = .(CellType, peak, condition)]
  out$group <- label
  return(out)
}

per_peak_condition_all <- list()
cell_types <- sort(unique(ATAC$CellType))
all_categories <- unique(tobias_dat$category)[!is.na(unique(tobias_dat$category))]

for (cell_check in cell_types) {
  if (!cell_check %in% tobias_dat$cell_type) next
  message("Processing: ", cell_check)
  
  obj_ct <- subset(ATAC, subset = CellType == cell_check)
  if (ncol(obj_ct) == 0) next
  
  cell_results <- list()
  for (cat in all_categories) {
    df_cat <- tobias_dat %>% filter(cell_type == cell_check, category == cat)
    gr_cat <- make_gr(df_cat)
    if (length(gr_cat) > 0) {
      acc_cat <- overlap_to_df(obj_ct, gr_cat, cat, cell_check)
      if (!is.null(acc_cat)) cell_results[[cat]] <- acc_cat
    }
  }
  if (length(cell_results) > 0) per_peak_condition_all[[cell_check]] <- bind_rows(cell_results)
}

acc_combined <- bind_rows(per_peak_condition_all)
if (nrow(acc_combined) == 0L) stop("No overlapping peaks found.")

# ==============================================================================
# 4. Statistical Testing (Compare against 'In All Conditions')
# ==============================================================================
message("Running Statistical Tests against 'In All Conditions'...")

# Reference Category
ref_cat <- "In All Conditions"

# Container for results
stats_list <- list()

# Get unique combos of CellType and Condition
combos <- acc_combined %>% select(CellType, condition) %>% distinct()

for (i in 1:nrow(combos)) {
  ct <- combos$CellType[i]
  cond <- combos$condition[i]
  
  # Get data for this specific square (Cell + Condition)
  subset_data <- acc_combined %>% filter(CellType == ct, condition == cond)
  
  # Check if Reference exists
  ref_data <- subset_data %>% filter(group == ref_cat)
  if (nrow(ref_data) < 3) next # Need at least 3 peaks to test
  
  # Test every other category against Reference
  other_cats <- setdiff(unique(subset_data$group), ref_cat)
  
  for (test_cat in other_cats) {
    test_data <- subset_data %>% filter(group == test_cat)
    if (nrow(test_data) < 3) next
    
    # Wilcoxon Rank Sum Test
    res <- wilcox.test(test_data$mean_accessibility, ref_data$mean_accessibility)
    
    stats_list[[length(stats_list) + 1]] <- data.frame(
      CellType = ct,
      condition = cond,
      group = test_cat,
      p.value = res$p.value
    )
  }
}

stats_df <- bind_rows(stats_list)

# Add stars
stats_df <- stats_df %>%
  mutate(stars = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE ~ ""
  ))

# ==============================================================================
# 5. Visualization with Stats
# ==============================================================================

# --- Plot 1: Boxplot ---
p_box <- ggplot(acc_combined, aes(x = group, y = mean_accessibility, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(condition ~ CellType, scales = "free_y", switch = "y") +
  labs(x = NULL, y = "Mean accessibility (TF-IDF)", title = paste("Accessibility at", motif_pattern, "sites")) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

if (exists("theme_mytheme")) p_box <- p_box + theme_mytheme

# --- Plot 2: Heatmap with Significance Stars ---
summ_heat <- acc_combined %>%
  group_by(group, CellType, condition) %>%
  summarise(mean_accessibility = mean(mean_accessibility, na.rm = TRUE), .groups = "drop")

# Merge Stats into Heatmap Data
summ_heat <- summ_heat %>%
  left_join(stats_df, by = c("CellType", "condition", "group")) %>%
  mutate(stars = ifelse(is.na(stars), "", stars)) # Fill NAs (e.g. Ref category) with empty

# Shared Limits
vmin <- min(summ_heat$mean_accessibility, na.rm = TRUE)
vmax <- max(summ_heat$mean_accessibility, na.rm = TRUE)

p_heat <- ggplot(summ_heat, aes(x = condition, y = CellType, fill = mean_accessibility)) +
  geom_tile(width = 0.95, height = 0.95, colour = NA) +
  
  # Overlay Stars
  geom_text(aes(label = stars), color = "white", size = 5, vjust = 0.75) +
  
  scale_fill_gradientn(
    colours = viridis(256),
    limits = c(vmin, vmax),
    name = "Mean\nAccess."
  ) +
  facet_wrap(~ group, nrow = 2) + 
  coord_fixed() +
  labs(x = "Condition", y = "Cell type", 
       title = paste("Accessibility Heatmap -", motif_pattern),
       subtitle = "Significance vs 'In All Conditions' (* P<0.05)") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold", size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

if (exists("theme_mytheme")) p_heat <- p_heat + theme_mytheme

# ==============================================================================
# 6. Output
# ==============================================================================
print(p_box)
print(p_heat)

ggsave(paste0(motif_pattern, "_AllCats_Accessibility_Boxplot.pdf"), plot = p_box, width = 16, height = 12)
ggsave(paste0(motif_pattern, "_AllCats_Accessibility_Heatmap_Stats.pdf"), plot = p_heat, width = 12, height = 10)

# Save Stats Table
write.csv(stats_df, paste0(motif_pattern, "_Accessibility_Stats.csv"), row.names = FALSE)

message("Analysis Complete. Heatmap now includes significance stars.")



# Read in all Fos Sites message(paste("Loaded TOBIAS data. Total rows:", nrow(Tobias_df)))
out_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2023-10-31.Cicero/2026-01-20_Recreated_Links")
TOBIAS = readRDS(file.path(out_dir, "/Tobias_df_Classified_Full.rds"))



library(dplyr)
# Read in all Fos Sites message(paste("Loaded TOBIAS data. Total rows:", nrow(Tobias_df)))
out_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2023-10-31.Cicero/2026-01-20_Recreated_Links")
TOBIAS = readRDS(file.path(out_dir, "/Tobias_df_Classified_Full.rds"))
FOS = TOBIAS %>% filter(motif == 'Fos_Homo_sapiens_DBD_1.00_M08069_2.00')

library(dplyr)
# Read in all Fos Sites message(paste("Loaded TOBIAS data. Total rows:", nrow(Tobias_df)))
out_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2023-10-31.Cicero/2026-01-20_Recreated_Links")
TOBIAS = readRDS(file.path(out_dir, "/Tobias_df_Classified_Full.rds"))
FOS = TOBIAS %>% filter(motif == 'Fos_Homo_sapiens_DBD_1.00_M08069_2.00')

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggalluvial)
  library(forcats)
  library(purrr)
})

# ---- inputs ----
# expects: TOBIAS with columns: celltype, in_YS, in_MAS, in_MAE (TRUE/FALSE)
df <- FOS
df$celltype = df$cell
stopifnot(all(c("celltype","in_YS","in_MAS","in_MAE") %in% colnames(df)))

# ---- your palette ----
RNA_pallate = c(
  '#A53D99','#F99349','#87CEEB','#FDDA62','#808080','#BEBADA',
  "#A81B1B","#101369","#606692","#3A85A8","#42977E","#591f0a",
  "#3b7818","#B6D66A","#1CE382","#BF862B","#CC6A6F","#E086B5",
  "#5F1296","#5F1296","#B784A7"
)

to_state <- function(x) {
  x <- as.logical(x)
  dplyr::case_when(
    is.na(x) ~ "NA",
    x ~ "TRUE",
    TRUE ~ "FALSE"
  )
}

df2 <- df %>%
  mutate(
    in_YS  = as.logical(in_YS),
    in_MAS = as.logical(in_MAS),
    in_MAE = as.logical(in_MAE)
  )

allu <- df2 %>%
  mutate(.id = row_number()) %>%
  transmute(
    celltype,
    .id,
    YS  = to_state(in_YS),
    MAS = to_state(in_MAS),
    MAE = to_state(in_MAE)
  ) %>%
  group_by(celltype, YS, MAS, MAE) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    traj = dplyr::case_when(
      YS=="TRUE"  & MAS=="TRUE"  & MAE=="TRUE"  ~ "TTT (TRUE all)",            # grey
      YS=="FALSE" & MAS=="TRUE"  & MAE=="FALSE" ~ "FTF (TRUE only MAS)",       # green (highlight)
      YS=="TRUE"  & MAS=="FALSE" & MAE=="TRUE"  ~ "TFT (TRUE YS+MAE)",         # purple (highlight)
      YS=="FALSE" & MAS=="TRUE"  & MAE=="TRUE"  ~ "FTT (TRUE MAS+MAE)",        # red
      YS=="TRUE"  & MAS=="TRUE"  & MAE=="FALSE" ~ "TTF (drops at MAE)",
      YS=="TRUE"  & MAS=="FALSE" & MAE=="FALSE" ~ "TFF (TRUE only YS)",
      YS=="FALSE" & MAS=="FALSE" & MAE=="TRUE"  ~ "FFT (TRUE only MAE)",
      YS=="FALSE" & MAS=="FALSE" & MAE=="FALSE" ~ "FFF (FALSE all)",
      TRUE ~ "has NA"
    )
  )

traj_levels <- c(
  "FTF (TRUE only MAS)",
  "TFT (TRUE YS+MAE)",
  "FTT (TRUE MAS+MAE)",
  "TTT (TRUE all)",
  "TTF (drops at MAE)",
  "TFF (TRUE only YS)",
  "FFT (TRUE only MAE)",
  "FFF (FALSE all)",
  "has NA"
)

# requested:
# - TTT grey
# - FTF green
# - FTT red
# plus sensible colours to make FTF and TFT pop
traj_cols <- c(
  "TTT (TRUE all)"       = "grey",      # grey (your palette)
  "FTF (TRUE only MAS)"  = 'green',# green highlight (#1CE382)
  "FTT (TRUE MAS+MAE)"   = 'red', # red (#A81B1B)
  "TFT (TRUE YS+MAE)"    = 'dark green',# purple (#5F1296)
  "TTF (drops at MAE)"   = 'blue', # orange
  "TFF (TRUE only YS)"   = 'dark red', # dark blue
  "FFT (TRUE only MAE)"  = 'dark blue', # yellow
  "FFF (FALSE all)"      = RNA_pallate[6], # neutral-ish
  "has NA"               = RNA_pallate[18] # pink
)

allu <- allu %>%
  mutate(
    traj = factor(traj, levels = traj_levels),
    YS  = factor(YS,  levels = c("FALSE","TRUE","NA")),
    MAS = factor(MAS, levels = c("FALSE","TRUE","NA")),
    MAE = factor(MAE, levels = c("FALSE","TRUE","NA"))
  )

plot_one <- function(ct) {
  d <- allu %>% filter(celltype == ct)
  ggplot(d, aes(axis1 = YS, axis2 = MAS, axis3 = MAE, y = n)) +
    geom_alluvium(aes(fill = traj), width = 0.22, alpha = 0.90) +
    geom_stratum(width = 0.22, colour = "black", linewidth = 0.3, fill = "white") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c("in_YS","in_MAS","in_MAE"), expand = c(0.05, 0.05)) +
    scale_fill_manual(values = traj_cols, drop = FALSE) +
    labs(title = paste0("Alluvial: ", ct), x = NULL, y = "Count", fill = "Trajectory") +
    theme_classic(base_size = 12) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

plots <- lapply(sort(unique(allu$celltype)), plot_one)
names(plots) <- sort(unique(allu$celltype))

for (ct in names(plots)) print(plots[[ct]])

# optional save:
# dir.create("alluvial_by_celltype", showWarnings = FALSE)
# for (ct in names(plots)) {
#   ggsave(file.path("alluvial_by_celltype", paste0("alluvial_", gsub("[^A-Za-z0-9_]+","_", ct), ".pdf")),
#          plots[[ct]], width = 6.5, height = 4.5)
# }


getwd()
tobias_dir <- file.path(base_dir, "analysis/2026-01-19.TOBIAS_Classification")  # base_dir is HEART_SNRNA_ROOT here (set above), so "analysis/" is required
setwd(tobias_dir)
out_pdf <- "TOBIAS_FOS_alluvial_by_celltype.pdf"

pdf(out_pdf, width = 7, height = 5, onefile = TRUE)
for (ct in names(plots)) {
  print(plots[[ct]])
}
dev.off()

message("Saved: ", normalizePath(out_pdf))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggalluvial)
})

# ---- inputs ----
# expects: df with columns: celltype, in_YS, in_MAS, in_MAE (TRUE/FALSE)
df <- FOS
df$celltype <- df$cell
stopifnot(all(c("celltype","in_YS","in_MAS","in_MAE") %in% colnames(df)))

# ---- your palette (kept) ----
RNA_pallate = c(
  '#A53D99','#F99349','#87CEEB','#FDDA62','#808080','#BEBADA',
  "#A81B1B","#101369","#606692","#3A85A8","#42977E","#591f0a",
  "#3b7818","#B6D66A","#1CE382","#BF862B","#CC6A6F","#E086B5",
  "#5F1296","#5F1296","#B784A7"
)

to_state <- function(x) {
  x <- as.logical(x)
  dplyr::case_when(
    is.na(x) ~ "NA",
    x ~ "TRUE",
    TRUE ~ "FALSE"
  )
}

df2 <- df %>%
  mutate(
    in_YS  = as.logical(in_YS),
    in_MAS = as.logical(in_MAS),
    in_MAE = as.logical(in_MAE)
  )

# ----------------------------
# Build pooled alluvial (ALL cell types together)
# ----------------------------
allu_all <- df2 %>%
  transmute(
    YS  = to_state(in_YS),
    MAS = to_state(in_MAS),
    MAE = to_state(in_MAE)
  ) %>%
  count(YS, MAS, MAE, name = "n") %>%
  mutate(
    traj = dplyr::case_when(
      YS=="TRUE"  & MAS=="TRUE"  & MAE=="TRUE"  ~ "TTT (TRUE all)",            # grey
      YS=="FALSE" & MAS=="TRUE"  & MAE=="FALSE" ~ "FTF (TRUE only MAS)",       # green
      YS=="TRUE"  & MAS=="FALSE" & MAE=="TRUE"  ~ "TFT (TRUE YS+MAE)",         # purple-ish
      YS=="FALSE" & MAS=="TRUE"  & MAE=="TRUE"  ~ "FTT (TRUE MAS+MAE)",        # red
      YS=="TRUE"  & MAS=="TRUE"  & MAE=="FALSE" ~ "TTF (drops at MAE)",
      YS=="TRUE"  & MAS=="FALSE" & MAE=="FALSE" ~ "TFF (TRUE only YS)",
      YS=="FALSE" & MAS=="FALSE" & MAE=="TRUE"  ~ "FFT (TRUE only MAE)",
      YS=="FALSE" & MAS=="FALSE" & MAE=="FALSE" ~ "FFF (FALSE all)",
      TRUE ~ "has NA"
    )
  )

traj_levels <- c(
  "FTF (TRUE only MAS)",
  "TFT (TRUE YS+MAE)",
  "FTT (TRUE MAS+MAE)",
  "TTT (TRUE all)",
  "TTF (drops at MAE)",
  "TFF (TRUE only YS)",
  "FFT (TRUE only MAE)",
  "FFF (FALSE all)",
  "has NA"
)

traj_cols <- c(
  "TTT (TRUE all)"       = "grey",
  "FTF (TRUE only MAS)"  = "green",
  "FTT (TRUE MAS+MAE)"   = "red",
  "TFT (TRUE YS+MAE)"    = "darkgreen",
  "TTF (drops at MAE)"   = "blue",
  "TFF (TRUE only YS)"   = "darkred",
  "FFT (TRUE only MAE)"  = "darkblue",
  "FFF (FALSE all)"      = RNA_pallate[6],
  "has NA"               = RNA_pallate[18]
)

allu_all <- allu_all %>%
  mutate(
    traj = factor(traj, levels = traj_levels),
    YS  = factor(YS,  levels = c("FALSE","TRUE","NA")),
    MAS = factor(MAS, levels = c("FALSE","TRUE","NA")),
    MAE = factor(MAE, levels = c("FALSE","TRUE","NA"))
  )

# Add legend labels with pooled n and %
leg_df <- allu_all %>%
  group_by(traj) %>%
  summarise(n_traj = sum(n), .groups = "drop") %>%
  mutate(
    pct_traj = 100 * n_traj / sum(n_traj),
    lab = paste0(as.character(traj), " (n=", n_traj, ", ", sprintf("%.1f", pct_traj), "%)")
  )
traj_labs <- setNames(leg_df$lab, as.character(leg_df$traj))

# ----------------------------
# ONE plot
# ----------------------------
p_all <- ggplot(allu_all, aes(axis1 = YS, axis2 = MAS, axis3 = MAE, y = n)) +
  geom_alluvium(aes(fill = traj), width = 0.22, alpha = 0.90) +
  geom_stratum(width = 0.22, colour = "black", linewidth = 0.3, fill = "white") +
  geom_text(
    stat = "stratum",
    aes(label = paste0(after_stat(stratum), "\n", sprintf("%.1f%%", after_stat(prop) * 100))),
    size = 3,
    lineheight = 0.95
  ) +
  scale_x_discrete(limits = c("in_YS","in_MAS","in_MAE"), expand = c(0.05, 0.05)) +
  scale_fill_manual(values = traj_cols, breaks = names(traj_labs), labels = unname(traj_labs), drop = FALSE) +
  labs(
    title = "Alluvial (ALL cell types pooled): in_YS → in_MAS → in_MAE",
    x = NULL, y = "Count", fill = "Trajectory"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.text = element_text(size = 9)
  )

print(p_all)

# optional save
out_pdf <- "alluvial_presence_ALL_celltypes.pdf"
pdf(out_pdf, width = 7.2, height = 4.8, onefile = TRUE)
print(p_all)
dev.off()
message("Saved: ", normalizePath(out_pdf))

out_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2023-10-31.Cicero/2026-01-20_Recreated_Links")
TOBIAS = readRDS(file.path(out_dir, "/Tobias_df_Classified_Full_scores.rds"))
FOS = TOBIAS %>% filter(motif == 'Fos_Homo_sapiens_DBD_1.00_M08069_2.00')

head(FOS)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggalluvial)
})

# ----------------------------
# Load + subset
# ----------------------------
out_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2023-10-31.Cicero/2026-01-20_Recreated_Links")
TOBIAS <- readRDS(file.path(out_dir, "/Tobias_df_Classified_Full_scores.rds"))
FOS <- TOBIAS %>% filter(motif == "Fos_Homo_sapiens_DBD_1.00_M08069_2.00")

df <- FOS
df$celltype <- df$cell

stopifnot(all(c("celltype","YS_score","MAS_score","MAE_score") %in% colnames(df)))

# ----------------------------
# Palette (yours)
# ----------------------------
RNA_pallate = c(
  '#A53D99','#F99349','#87CEEB','#FDDA62','#808080','#BEBADA',
  "#A81B1B","#101369","#606692","#3A85A8","#42977E","#591f0a",
  "#3b7818","#B6D66A","#1CE382","#BF862B","#CC6A6F","#E086B5",
  "#5F1296","#5F1296","#B784A7"
)

# ----------------------------
# Thresholded direction categories
# ----------------------------
thresh <- 0.1

dir_state <- function(from, to, thresh = 0.1) {
  dplyr::case_when(
    is.na(from) | is.na(to) ~ "NA",
    (to - from) >=  thresh  ~ "Increasing",
    (to - from) <= -thresh  ~ "Decreasing",
    TRUE ~ "NoChange"
  )
}

df2 <- df %>%
  mutate(
    YS_score  = as.numeric(YS_score),
    MAS_score = as.numeric(MAS_score),
    MAE_score = as.numeric(MAE_score),
    YS_to_MAS  = dir_state(YS_score,  MAS_score, thresh = thresh),
    MAS_to_MAE = dir_state(MAS_score, MAE_score, thresh = thresh)
  )

# ----------------------------
# Alluvial counts + within-celltype percentages
# ----------------------------
allu <- df2 %>%
  count(celltype, YS_to_MAS, MAS_to_MAE, name = "n") %>%
  group_by(celltype) %>%
  mutate(
    pct = 100 * n / sum(n),
    traj = paste0(YS_to_MAS, " -> ", MAS_to_MAE)
  ) %>%
  ungroup() %>%
  mutate(
    YS_to_MAS  = factor(YS_to_MAS,  levels = c("Decreasing","NoChange","Increasing","NA")),
    MAS_to_MAE = factor(MAS_to_MAE, levels = c("Decreasing","NoChange","Increasing","NA"))
  )

traj_levels <- c(
  "Increasing -> Increasing",
  "Decreasing -> Decreasing",
  "Increasing -> Decreasing",
  "Decreasing -> Increasing",
  "NoChange -> Increasing",
  "NoChange -> Decreasing",
  "Increasing -> NoChange",
  "Decreasing -> NoChange",
  "NoChange -> NoChange",
  "NA -> Increasing","NA -> Decreasing","NA -> NoChange",
  "Increasing -> NA","Decreasing -> NA","NoChange -> NA",
  "NA -> NA"
)

traj_cols <- c(
  "Increasing -> Increasing" = RNA_pallate[15],
  "Decreasing -> Decreasing" = RNA_pallate[7],
  "Increasing -> Decreasing" = RNA_pallate[19],
  "Decreasing -> Increasing" = RNA_pallate[2],
  "NoChange -> Increasing"   = RNA_pallate[3],
  "NoChange -> Decreasing"   = RNA_pallate[4],
  "Increasing -> NoChange"   = RNA_pallate[1],
  "Decreasing -> NoChange"   = RNA_pallate[10],
  "NoChange -> NoChange"     = RNA_pallate[5],
  "NA -> Increasing"         = RNA_pallate[6],
  "NA -> Decreasing"         = RNA_pallate[6],
  "NA -> NoChange"           = RNA_pallate[6],
  "Increasing -> NA"         = RNA_pallate[6],
  "Decreasing -> NA"         = RNA_pallate[6],
  "NoChange -> NA"           = RNA_pallate[6],
  "NA -> NA"                 = RNA_pallate[6]
)

allu <- allu %>%
  mutate(traj = factor(traj, levels = traj_levels))

# ----------------------------
# Legend labels WITH n and % (per celltype) computed inside plot_one
# ----------------------------
plot_one <- function(ct) {
  d <- allu %>% filter(celltype == ct)

  # build per-celltype legend labels: "Trajectory (n=..., x.x%)"
  leg_df <- d %>%
    group_by(traj) %>%
    summarise(n_traj = sum(n), .groups = "drop") %>%
    mutate(
      pct_traj = 100 * n_traj / sum(n_traj),
      lab = paste0(as.character(traj), " (n=", n_traj, ", ", sprintf("%.1f", pct_traj), "%)")
    )

  # named vector for labels (levels not present won't be shown)
  traj_labs <- setNames(leg_df$lab, as.character(leg_df$traj))

  ggplot(d, aes(axis1 = YS_to_MAS, axis2 = MAS_to_MAE, y = n)) +
    geom_alluvium(aes(fill = traj), width = 0.25, alpha = 0.90) +
    geom_stratum(width = 0.25, colour = "black", linewidth = 0.3, fill = "white") +
    # internal stratum labels: state + percent of total at that axis
    geom_text(
      stat = "stratum",
      aes(label = paste0(after_stat(stratum), "\n", sprintf("%.1f%%", after_stat(prop) * 100))),
      size = 2.8,
      lineheight = 0.95
    ) +
    scale_x_discrete(limits = c("YS→MAS", "MAS→MAE"), expand = c(0.05, 0.05)) +
    scale_fill_manual(
      values = traj_cols,
      breaks = names(traj_labs),
      labels = unname(traj_labs),
      drop = FALSE
    ) +
    labs(
      title = paste0("Alluvial (threshold ", thresh, "): ", ct),
      x = NULL, y = "Count of sites", fill = "Trajectory"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.text = element_text(size = 9)
    )
}

plots <- lapply(sort(unique(allu$celltype)), plot_one)
names(plots) <- sort(unique(allu$celltype))
for (ct in names(plots)) print(plots[[ct]])

# optional save all to one pdf
 out_pdf <- paste0("alluvial_direction_counts", thresh, ".pdf")
pdf(out_pdf, width = 7, height = 5, onefile = TRUE)
 for (ct in names(plots)) print(plots[[ct]])
 dev.off()
message("Saved: ", normalizePath(out_pdf))



# YES: use ONE per-celltype threshold computed from BOTH transitions' |delta|
# Then apply that same threshold to YS→MAS and MAS→MAE.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggalluvial)
})

# ----------------------------
# Load + subset
# ----------------------------
out_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2023-10-31.Cicero/2026-01-20_Recreated_Links")
TOBIAS <- readRDS(file.path(out_dir, "/Tobias_df_Classified_Full_scores.rds"))
FOS <- TOBIAS %>% filter(motif == "Fos_Homo_sapiens_DBD_1.00_M08069_2.00")

df <- FOS
df$celltype <- df$cell
stopifnot(all(c("celltype","YS_score","MAS_score","MAE_score") %in% colnames(df)))

# ----------------------------
# Palette (yours)
# ----------------------------
RNA_pallate = c(
  '#A53D99','#F99349','#87CEEB','#FDDA62','#808080','#BEBADA',
  "#A81B1B","#101369","#606692","#3A85A8","#42977E","#591f0a",
  "#3b7818","#B6D66A","#1CE382","#BF862B","#CC6A6F","#E086B5",
  "#5F1296","#5F1296","#B784A7"
)

# ----------------------------
# One threshold per celltype (based on pooled |delta|)
# ----------------------------
q_abs   <- 0.25
min_thr <- 0.0

dir_from_delta <- function(delta, thr) {
  dplyr::case_when(
    is.na(delta) | is.na(thr) ~ "NA",
    delta >=  thr ~ "Increasing",
    delta <= -thr ~ "Decreasing",
    TRUE ~ "NoChange"
  )
}

df2 <- df %>%
  mutate(
    YS_score  = as.numeric(YS_score),
    MAS_score = as.numeric(MAS_score),
    MAE_score = as.numeric(MAE_score),
    d_YS_MAS  = MAS_score - YS_score,
    d_MAS_MAE = MAE_score - MAS_score,
    abs_d_YS_MAS  = abs(d_YS_MAS),
    abs_d_MAS_MAE = abs(d_MAS_MAE)
  )

thr_tbl <- df2 %>%
  group_by(celltype) %>%
  summarise(
    thr = max(
      min_thr,
      as.numeric(quantile(c(abs_d_YS_MAS, abs_d_MAS_MAE), probs = q_abs, na.rm = TRUE))
    ),
    n_sites = n(),
    .groups = "drop"
  )

print(thr_tbl %>% arrange(celltype))

df3 <- df2 %>%
  left_join(thr_tbl, by = "celltype") %>%
  mutate(
    YS_to_MAS  = dir_from_delta(d_YS_MAS,  thr),
    MAS_to_MAE = dir_from_delta(d_MAS_MAE, thr)
  )

# sanity: percent called changing (should be similar-ish across transitions now)
sanity <- df3 %>%
  group_by(celltype) %>%
  summarise(
    pct_change_YS_MAS  = 100 * mean(YS_to_MAS %in% c("Increasing","Decreasing"), na.rm = TRUE),
    pct_change_MAS_MAE = 100 * mean(MAS_to_MAE %in% c("Increasing","Decreasing"), na.rm = TRUE),
    .groups = "drop"
  )
print(sanity %>% arrange(desc(pct_change_YS_MAS)))

# ----------------------------
# Alluvial counts + within-celltype % (y axis = counts)
# ----------------------------
allu <- df3 %>%
  count(celltype, YS_to_MAS, MAS_to_MAE, name = "n") %>%
  group_by(celltype) %>%
  mutate(
    pct = 100 * n / sum(n),
    traj = paste0(YS_to_MAS, " -> ", MAS_to_MAE)
  ) %>%
  ungroup() %>%
  mutate(
    YS_to_MAS  = factor(YS_to_MAS,  levels = c("Decreasing","NoChange","Increasing","NA")),
    MAS_to_MAE = factor(MAS_to_MAE, levels = c("Decreasing","NoChange","Increasing","NA"))
  )

traj_levels <- c(
  "Increasing -> Increasing",
  "Decreasing -> Decreasing",
  "Increasing -> Decreasing",
  "Decreasing -> Increasing",
  "NoChange -> Increasing",
  "NoChange -> Decreasing",
  "Increasing -> NoChange",
  "Decreasing -> NoChange",
  "NoChange -> NoChange",
  "NA -> Increasing","NA -> Decreasing","NA -> NoChange",
  "Increasing -> NA","Decreasing -> NA","NoChange -> NA",
  "NA -> NA"
)

traj_cols <- c(
  "Increasing -> Increasing" = 'red',
  "Decreasing -> Decreasing" = 'dark red',
  "Increasing -> Decreasing" = 'green',
  "Decreasing -> Increasing" = 'dark green',
  "NoChange -> Increasing"   = 'blue',
  "NoChange -> Decreasing"   = 'light blue',
  "Increasing -> NoChange"   = 'pink',
  "Decreasing -> NoChange"   = 'orange',
  "NoChange -> NoChange"     = 'grey'

)

allu <- allu %>%
  mutate(traj = factor(traj, levels = traj_levels))
# Add a RIGHT-SIDE "final" axis that shows trajectory % (and n), while keeping MAS→MAE % labels too.
# Drop this in AFTER you build `allu` and BEFORE `plot_one()`.

allu <- allu %>%
  mutate(
    traj_chr = as.character(traj),
    traj_abbr = dplyr::recode(
      traj_chr,
      "Increasing -> Increasing" = "Inc→Inc",
      "Decreasing -> Decreasing" = "Dec→Dec",
      "Increasing -> Decreasing" = "Inc→Dec",
      "Decreasing -> Increasing" = "Dec→Inc",
      "NoChange -> Increasing"   = "NC→Inc",
      "NoChange -> Decreasing"   = "NC→Dec",
      "Increasing -> NoChange"   = "Inc→NC",
      "Decreasing -> NoChange"   = "Dec→NC",
      "NoChange -> NoChange"     = "NC→NC",
      "NA -> Increasing"         = "NA→Inc",
      "NA -> Decreasing"         = "NA→Dec",
      "NA -> NoChange"           = "NA→NC",
      "Increasing -> NA"         = "Inc→NA",
      "Decreasing -> NA"         = "Dec→NA",
      "NoChange -> NA"           = "NC→NA",
      "NA -> NA"                 = "NA→NA",
      .default = traj_chr
    ),
    traj_abbr = factor(traj_abbr, levels = unique(dplyr::recode(
      traj_levels,
      "Increasing -> Increasing" = "Inc→Inc",
      "Decreasing -> Decreasing" = "Dec→Dec",
      "Increasing -> Decreasing" = "Inc→Dec",
      "Decreasing -> Increasing" = "Dec→Inc",
      "NoChange -> Increasing"   = "NC→Inc",
      "NoChange -> Decreasing"   = "NC→Dec",
      "Increasing -> NoChange"   = "Inc→NC",
      "Decreasing -> NoChange"   = "Dec→NC",
      "NoChange -> NoChange"     = "NC→NC",
      "NA -> Increasing"         = "NA→Inc",
      "NA -> Decreasing"         = "NA→Dec",
      "NA -> NoChange"           = "NA→NC",
      "Increasing -> NA"         = "Inc→NA",
      "Decreasing -> NA"         = "Dec→NA",
      "NoChange -> NA"           = "NC→NA",
      "NA -> NA"                 = "NA→NA",
      .default = traj_levels
    ))),
    traj_chr = NULL
  )

plot_one <- function(ct) {
  d <- allu %>% filter(celltype == ct)

  # legend labels with n and % (per celltype)
  leg_df <- d %>%
    group_by(traj) %>%
    summarise(n_traj = sum(n), .groups = "drop") %>%
    mutate(
      pct_traj = 100 * n_traj / sum(n_traj),
      lab = paste0(as.character(traj), " (n=", n_traj, ", ", sprintf("%.1f", pct_traj), "%)")
    )
  traj_labs <- setNames(leg_df$lab, as.character(leg_df$traj))

  # celltype-specific pooled threshold for title (from df3)
  thr_ct <- unique(df3$thr[df3$celltype == ct])
  thr_ct <- if (length(thr_ct) == 1) thr_ct else NA_real_

  ggplot(d, aes(axis1 = YS_to_MAS, axis2 = MAS_to_MAE, axis3 = traj_abbr, y = n)) +
    geom_alluvium(aes(fill = traj), width = 0.22, alpha = 0.90) +
    geom_stratum(width = 0.22, colour = "black", linewidth = 0.3, fill = "white") +
    geom_text(
      stat = "stratum",
      aes(label = ifelse(
        after_stat(x) == 3,
        paste0(sprintf("%.1f%%", after_stat(prop) * 100), "\n(n=", after_stat(count), ")"),
        paste0(after_stat(stratum), "\n", sprintf("%.1f%%", after_stat(prop) * 100))
      )),
      size = 2.7,
      lineheight = 0.95
    ) +
    scale_x_discrete(limits = c("YS->MAS", "MAS->MAE", "Final (trajectory)"), expand = c(0.05, 0.05)) +
    scale_fill_manual(
      values = traj_cols,
      breaks = names(traj_labs),
      labels = unname(traj_labs),
      drop = FALSE
    ) +
    labs(
      title = paste0("Alluvial ",ct," thr=", sprintf("%.3f", thr_ct), "): ", ct),
      x = NULL, y = "Count of sites", fill = "Trajectory"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.text = element_text(size = 9)
    )
}

plots <- lapply(sort(unique(allu$celltype)), plot_one)
names(plots) <- sort(unique(allu$celltype))
for (ct in names(plots)) print(plots[[ct]])

 out_pdf <- paste0("alluvial_direction_counts_split", thresh, ".pdf")
pdf(out_pdf, width = 7, height = 5, onefile = TRUE)
 for (ct in names(plots)) print(plots[[ct]])
 dev.off()
message("Saved: ", normalizePath(out_pdf))

# Plot ALL cell types together in ONE alluvial
# Key change: collapse counts over celltype (but keep celltype as a stratum fill if you want)
# This keeps your per-celltype thresholding upstream (df3), then aggregates trajectories across cell types.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggalluvial)
})

# ----------------------------
# Load + subset
# ----------------------------
out_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2023-10-31.Cicero/2026-01-20_Recreated_Links")
TOBIAS <- readRDS(file.path(out_dir, "/Tobias_df_Classified_Full_scores.rds"))
FOS <- TOBIAS %>% filter(motif == "Fos_Homo_sapiens_DBD_1.00_M08069_2.00")

df <- FOS
df$celltype <- df$cell
stopifnot(all(c("celltype","YS_score","MAS_score","MAE_score") %in% colnames(df)))

# ----------------------------
# One threshold per celltype (based on pooled |delta|)
# ----------------------------
q_abs   <- 0.25
min_thr <- 0.0

dir_from_delta <- function(delta, thr) {
  dplyr::case_when(
    is.na(delta) | is.na(thr) ~ "NA",
    delta >=  thr ~ "Increasing",
    delta <= -thr ~ "Decreasing",
    TRUE ~ "NoChange"
  )
}

df2 <- df %>%
  mutate(
    YS_score  = as.numeric(YS_score),
    MAS_score = as.numeric(MAS_score),
    MAE_score = as.numeric(MAE_score),
    d_YS_MAS  = MAS_score - YS_score,
    d_MAS_MAE = MAE_score - MAS_score,
    abs_d_YS_MAS  = abs(d_YS_MAS),
    abs_d_MAS_MAE = abs(d_MAS_MAE)
  )

thr_tbl <- df2 %>%
  group_by(celltype) %>%
  summarise(
    thr = max(
      min_thr,
      as.numeric(quantile(c(abs_d_YS_MAS, abs_d_MAS_MAE), probs = q_abs, na.rm = TRUE))
    ),
    n_sites = n(),
    .groups = "drop"
  )

df3 <- df2 %>%
  left_join(thr_tbl, by = "celltype") %>%
  mutate(
    YS_to_MAS  = dir_from_delta(d_YS_MAS,  thr),
    MAS_to_MAE = dir_from_delta(d_MAS_MAE, thr)
  )

# ----------------------------
# Alluvial counts pooled across ALL cell types
# ----------------------------
traj_levels <- c(
  "Increasing -> Increasing",
  "Decreasing -> Decreasing",
  "Increasing -> Decreasing",
  "Decreasing -> Increasing",
  "NoChange -> Increasing",
  "NoChange -> Decreasing",
  "Increasing -> NoChange",
  "Decreasing -> NoChange",
  "NoChange -> NoChange",
  "NA -> Increasing","NA -> Decreasing","NA -> NoChange",
  "Increasing -> NA","Decreasing -> NA","NoChange -> NA",
  "NA -> NA"
)

traj_cols <- c(
  "Increasing -> Increasing" = "red",
  "Decreasing -> Decreasing" = "darkred",
  "Increasing -> Decreasing" = "green",
  "Decreasing -> Increasing" = "darkgreen",
  "NoChange -> Increasing"   = "blue",
  "NoChange -> Decreasing"   = "lightblue",
  "Increasing -> NoChange"   = "pink",
  "Decreasing -> NoChange"   = "orange",
  "NoChange -> NoChange"     = "grey",
  "NA -> Increasing"         = "grey70",
  "NA -> Decreasing"         = "grey70",
  "NA -> NoChange"           = "grey70",
  "Increasing -> NA"         = "grey70",
  "Decreasing -> NA"         = "grey70",
  "NoChange -> NA"           = "grey70",
  "NA -> NA"                 = "grey85"
)

allu_all <- df3 %>%
  count(YS_to_MAS, MAS_to_MAE, name = "n") %>%
  mutate(
    YS_to_MAS  = factor(YS_to_MAS,  levels = c("Decreasing","NoChange","Increasing","NA")),
    MAS_to_MAE = factor(MAS_to_MAE, levels = c("Decreasing","NoChange","Increasing","NA")),
    traj = paste0(as.character(YS_to_MAS), " -> ", as.character(MAS_to_MAE)),
    traj = factor(traj, levels = traj_levels)
  ) %>%
  mutate(
    traj_abbr = dplyr::recode(
      as.character(traj),
      "Increasing -> Increasing" = "Inc→Inc",
      "Decreasing -> Decreasing" = "Dec→Dec",
      "Increasing -> Decreasing" = "Inc→Dec",
      "Decreasing -> Increasing" = "Dec→Inc",
      "NoChange -> Increasing"   = "NC→Inc",
      "NoChange -> Decreasing"   = "NC→Dec",
      "Increasing -> NoChange"   = "Inc→NC",
      "Decreasing -> NoChange"   = "Dec→NC",
      "NoChange -> NoChange"     = "NC→NC",
      "NA -> Increasing"         = "NA→Inc",
      "NA -> Decreasing"         = "NA→Dec",
      "NA -> NoChange"           = "NA→NC",
      "Increasing -> NA"         = "Inc→NA",
      "Decreasing -> NA"         = "Dec→NA",
      "NoChange -> NA"           = "NC→NA",
      "NA -> NA"                 = "NA→NA",
      .default = as.character(traj)
    )
  )

# Legend labels with n and % (pooled)
leg_df <- allu_all %>%
  group_by(traj) %>%
  summarise(n_traj = sum(n), .groups = "drop") %>%
  mutate(
    pct_traj = 100 * n_traj / sum(n_traj),
    lab = paste0(as.character(traj), " (n=", n_traj, ", ", sprintf("%.1f", pct_traj), "%)")
  )
traj_labs <- setNames(leg_df$lab, as.character(leg_df$traj))

# ----------------------------
# ONE plot
# ----------------------------
p_all <- ggplot(allu_all, aes(axis1 = YS_to_MAS, axis2 = MAS_to_MAE, axis3 = traj_abbr, y = n)) +
  geom_alluvium(aes(fill = traj), width = 0.22, alpha = 0.90) +
  geom_stratum(width = 0.22, colour = "black", linewidth = 0.3, fill = "white") +
  geom_text(
    stat = "stratum",
    aes(label = ifelse(
      after_stat(x) == 3,
      paste0(sprintf("%.1f%%", after_stat(prop) * 100), "\n(n=", after_stat(count), ")"),
      paste0(after_stat(stratum), "\n", sprintf("%.1f%%", after_stat(prop) * 100))
    )),
    size = 2.7,
    lineheight = 0.95
  ) +
  scale_x_discrete(limits = c("YS->MAS", "MAS->MAE", "Final (trajectory)"), expand = c(0.05, 0.05)) +
  scale_fill_manual(
    values = traj_cols,
    breaks = names(traj_labs),
    labels = unname(traj_labs),
    drop = FALSE
  ) +
  labs(
    title = paste0("Alluvial (ALL cell types pooled) | per-celltype thr = q", q_abs),
    x = NULL, y = "Count of sites", fill = "Trajectory"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.text = element_text(size = 9)
  )

print(p_all)

# ----------------------------
# Save single PDF
# ----------------------------
out_pdf <- paste0("alluvial_direction_counts_ALL_celltypes_q", q_abs, ".pdf")
pdf(out_pdf, width = 7.5, height = 5.2, onefile = TRUE)
print(p_all)
dev.off()
message("Saved: ", normalizePath(out_pdf))

# ----------------------------
# OPTIONAL: if you want each celltype "represented" but still in ONE plot:
# Use facets instead of separate plots:
#   + facet_wrap(~celltype)
# BUT you asked for not splitting, so the pooled plot above is the direct answer.
# ----------------------------

# Plot ALL cell types together in ONE alluvial
# Key change: collapse counts over celltype (but keep celltype as a stratum fill if you want)
# This keeps your per-celltype thresholding upstream (df3), then aggregates trajectories across cell types.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggalluvial)
})

# ----------------------------
# Load + subset
# ----------------------------
out_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2023-10-31.Cicero/2026-01-20_Recreated_Links")
TOBIAS <- readRDS(file.path(out_dir, "/Tobias_df_Classified_Full_scores.rds"))
FOS <- TOBIAS %>% filter(motif == "Nfic_Homo_sapiens_DBD_0.98_M08164_2.00")

df <- FOS
df$celltype <- df$cell
stopifnot(all(c("celltype","YS_score","MAS_score","MAE_score") %in% colnames(df)))

# ----------------------------
# One threshold per celltype (based on pooled |delta|)
# ----------------------------
q_abs   <- 0.25
min_thr <- 0.0

dir_from_delta <- function(delta, thr) {
  dplyr::case_when(
    is.na(delta) | is.na(thr) ~ "NA",
    delta >=  thr ~ "Increasing",
    delta <= -thr ~ "Decreasing",
    TRUE ~ "NoChange"
  )
}

df2 <- df %>%
  mutate(
    YS_score  = as.numeric(YS_score),
    MAS_score = as.numeric(MAS_score),
    MAE_score = as.numeric(MAE_score),
    d_YS_MAS  = MAS_score - YS_score,
    d_MAS_MAE = MAE_score - MAS_score,
    abs_d_YS_MAS  = abs(d_YS_MAS),
    abs_d_MAS_MAE = abs(d_MAS_MAE)
  )

thr_tbl <- df2 %>%
  group_by(celltype) %>%
  summarise(
    thr = max(
      min_thr,
      as.numeric(quantile(c(abs_d_YS_MAS, abs_d_MAS_MAE), probs = q_abs, na.rm = TRUE))
    ),
    n_sites = n(),
    .groups = "drop"
  )

df3 <- df2 %>%
  left_join(thr_tbl, by = "celltype") %>%
  mutate(
    YS_to_MAS  = dir_from_delta(d_YS_MAS,  thr),
    MAS_to_MAE = dir_from_delta(d_MAS_MAE, thr)
  )

# ----------------------------
# Alluvial counts pooled across ALL cell types
# ----------------------------
traj_levels <- c(
  "Increasing -> Increasing",
  "Decreasing -> Decreasing",
  "Increasing -> Decreasing",
  "Decreasing -> Increasing",
  "NoChange -> Increasing",
  "NoChange -> Decreasing",
  "Increasing -> NoChange",
  "Decreasing -> NoChange",
  "NoChange -> NoChange",
  "NA -> Increasing","NA -> Decreasing","NA -> NoChange",
  "Increasing -> NA","Decreasing -> NA","NoChange -> NA",
  "NA -> NA"
)

traj_cols <- c(
  "Increasing -> Increasing" = "red",
  "Decreasing -> Decreasing" = "darkred",
  "Increasing -> Decreasing" = "green",
  "Decreasing -> Increasing" = "darkgreen",
  "NoChange -> Increasing"   = "blue",
  "NoChange -> Decreasing"   = "lightblue",
  "Increasing -> NoChange"   = "pink",
  "Decreasing -> NoChange"   = "orange",
  "NoChange -> NoChange"     = "grey",
  "NA -> Increasing"         = "grey70",
  "NA -> Decreasing"         = "grey70",
  "NA -> NoChange"           = "grey70",
  "Increasing -> NA"         = "grey70",
  "Decreasing -> NA"         = "grey70",
  "NoChange -> NA"           = "grey70",
  "NA -> NA"                 = "grey85"
)

allu_all <- df3 %>%
  count(YS_to_MAS, MAS_to_MAE, name = "n") %>%
  mutate(
    YS_to_MAS  = factor(YS_to_MAS,  levels = c("Decreasing","NoChange","Increasing","NA")),
    MAS_to_MAE = factor(MAS_to_MAE, levels = c("Decreasing","NoChange","Increasing","NA")),
    traj = paste0(as.character(YS_to_MAS), " -> ", as.character(MAS_to_MAE)),
    traj = factor(traj, levels = traj_levels)
  ) %>%
  mutate(
    traj_abbr = dplyr::recode(
      as.character(traj),
      "Increasing -> Increasing" = "Inc→Inc",
      "Decreasing -> Decreasing" = "Dec→Dec",
      "Increasing -> Decreasing" = "Inc→Dec",
      "Decreasing -> Increasing" = "Dec→Inc",
      "NoChange -> Increasing"   = "NC→Inc",
      "NoChange -> Decreasing"   = "NC→Dec",
      "Increasing -> NoChange"   = "Inc→NC",
      "Decreasing -> NoChange"   = "Dec→NC",
      "NoChange -> NoChange"     = "NC→NC",
      "NA -> Increasing"         = "NA→Inc",
      "NA -> Decreasing"         = "NA→Dec",
      "NA -> NoChange"           = "NA→NC",
      "Increasing -> NA"         = "Inc→NA",
      "Decreasing -> NA"         = "Dec→NA",
      "NoChange -> NA"           = "NC→NA",
      "NA -> NA"                 = "NA→NA",
      .default = as.character(traj)
    )
  )

# Legend labels with n and % (pooled)
leg_df <- allu_all %>%
  group_by(traj) %>%
  summarise(n_traj = sum(n), .groups = "drop") %>%
  mutate(
    pct_traj = 100 * n_traj / sum(n_traj),
    lab = paste0(as.character(traj), " (n=", n_traj, ", ", sprintf("%.1f", pct_traj), "%)")
  )
traj_labs <- setNames(leg_df$lab, as.character(leg_df$traj))

# ----------------------------
# ONE plot
# ----------------------------
p_all <- ggplot(allu_all, aes(axis1 = YS_to_MAS, axis2 = MAS_to_MAE, axis3 = traj_abbr, y = n)) +
  geom_alluvium(aes(fill = traj), width = 0.22, alpha = 0.90) +
  geom_stratum(width = 0.22, colour = "black", linewidth = 0.3, fill = "white") +
  geom_text(
    stat = "stratum",
    aes(label = ifelse(
      after_stat(x) == 3,
      paste0(sprintf("%.1f%%", after_stat(prop) * 100), "\n(n=", after_stat(count), ")"),
      paste0(after_stat(stratum), "\n", sprintf("%.1f%%", after_stat(prop) * 100))
    )),
    size = 2.7,
    lineheight = 0.95
  ) +
  scale_x_discrete(limits = c("YS->MAS", "MAS->MAE", "Final (trajectory)"), expand = c(0.05, 0.05)) +
  scale_fill_manual(
    values = traj_cols,
    breaks = names(traj_labs),
    labels = unname(traj_labs),
    drop = FALSE
  ) +
  labs(
    title = paste0("Alluvial (ALL cell types pooled) | per-celltype thr = q", q_abs),
    x = NULL, y = "Count of sites", fill = "Trajectory"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.text = element_text(size = 9)
  )

print(p_all)

# ----------------------------
# Save single PDF
# ----------------------------
out_pdf <- paste0("NFIC_control_alluvial_direction_counts_ALL_celltypes_q", q_abs, ".pdf")
pdf(out_pdf, width = 7.5, height = 5.2, onefile = TRUE)
print(p_all)
dev.off()
message("Saved: ", normalizePath(out_pdf))

# ----------------------------
# OPTIONAL: if you want each celltype "represented" but still in ONE plot:
# Use facets instead of separate plots:
#   + facet_wrap(~celltype)
# BUT you asked for not splitting, so the pooled plot above is the direct answer.
# ----------------------------

# YES: use ONE per-celltype threshold computed from BOTH transitions' |delta|
# Then apply that same threshold to YS→MAS and MAS→MAE.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggalluvial)
})

# ----------------------------
# Load + subset
# ----------------------------
out_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2023-10-31.Cicero/2026-01-20_Recreated_Links")
TOBIAS <- readRDS(file.path(out_dir, "/Tobias_df_Classified_Full_scores.rds"))
FOS <- TOBIAS %>% filter(motif == "Nfic_Homo_sapiens_DBD_0.98_M08164_2.00")

df <- FOS
df$celltype <- df$cell
stopifnot(all(c("celltype","YS_score","MAS_score","MAE_score") %in% colnames(df)))

# ----------------------------
# Palette (yours)
# ----------------------------
RNA_pallate = c(
  '#A53D99','#F99349','#87CEEB','#FDDA62','#808080','#BEBADA',
  "#A81B1B","#101369","#606692","#3A85A8","#42977E","#591f0a",
  "#3b7818","#B6D66A","#1CE382","#BF862B","#CC6A6F","#E086B5",
  "#5F1296","#5F1296","#B784A7"
)

# ----------------------------
# One threshold per celltype (based on pooled |delta|)
# ----------------------------
q_abs   <- 0.25
min_thr <- 0.0

dir_from_delta <- function(delta, thr) {
  dplyr::case_when(
    is.na(delta) | is.na(thr) ~ "NA",
    delta >=  thr ~ "Increasing",
    delta <= -thr ~ "Decreasing",
    TRUE ~ "NoChange"
  )
}

df2 <- df %>%
  mutate(
    YS_score  = as.numeric(YS_score),
    MAS_score = as.numeric(MAS_score),
    MAE_score = as.numeric(MAE_score),
    d_YS_MAS  = MAS_score - YS_score,
    d_MAS_MAE = MAE_score - MAS_score,
    abs_d_YS_MAS  = abs(d_YS_MAS),
    abs_d_MAS_MAE = abs(d_MAS_MAE)
  )

thr_tbl <- df2 %>%
  group_by(celltype) %>%
  summarise(
    thr = max(
      min_thr,
      as.numeric(quantile(c(abs_d_YS_MAS, abs_d_MAS_MAE), probs = q_abs, na.rm = TRUE))
    ),
    n_sites = n(),
    .groups = "drop"
  )

print(thr_tbl %>% arrange(celltype))

df3 <- df2 %>%
  left_join(thr_tbl, by = "celltype") %>%
  mutate(
    YS_to_MAS  = dir_from_delta(d_YS_MAS,  thr),
    MAS_to_MAE = dir_from_delta(d_MAS_MAE, thr)
  )

# sanity: percent called changing (should be similar-ish across transitions now)
sanity <- df3 %>%
  group_by(celltype) %>%
  summarise(
    pct_change_YS_MAS  = 100 * mean(YS_to_MAS %in% c("Increasing","Decreasing"), na.rm = TRUE),
    pct_change_MAS_MAE = 100 * mean(MAS_to_MAE %in% c("Increasing","Decreasing"), na.rm = TRUE),
    .groups = "drop"
  )
print(sanity %>% arrange(desc(pct_change_YS_MAS)))

# ----------------------------
# Alluvial counts + within-celltype % (y axis = counts)
# ----------------------------
allu <- df3 %>%
  count(celltype, YS_to_MAS, MAS_to_MAE, name = "n") %>%
  group_by(celltype) %>%
  mutate(
    pct = 100 * n / sum(n),
    traj = paste0(YS_to_MAS, " -> ", MAS_to_MAE)
  ) %>%
  ungroup() %>%
  mutate(
    YS_to_MAS  = factor(YS_to_MAS,  levels = c("Decreasing","NoChange","Increasing","NA")),
    MAS_to_MAE = factor(MAS_to_MAE, levels = c("Decreasing","NoChange","Increasing","NA"))
  )

traj_levels <- c(
  "Increasing -> Increasing",
  "Decreasing -> Decreasing",
  "Increasing -> Decreasing",
  "Decreasing -> Increasing",
  "NoChange -> Increasing",
  "NoChange -> Decreasing",
  "Increasing -> NoChange",
  "Decreasing -> NoChange",
  "NoChange -> NoChange",
  "NA -> Increasing","NA -> Decreasing","NA -> NoChange",
  "Increasing -> NA","Decreasing -> NA","NoChange -> NA",
  "NA -> NA"
)

traj_cols <- c(
  "Increasing -> Increasing" = 'red',
  "Decreasing -> Decreasing" = 'dark red',
  "Increasing -> Decreasing" = 'green',
  "Decreasing -> Increasing" = 'dark green',
  "NoChange -> Increasing"   = 'blue',
  "NoChange -> Decreasing"   = 'light blue',
  "Increasing -> NoChange"   = 'pink',
  "Decreasing -> NoChange"   = 'orange',
  "NoChange -> NoChange"     = 'grey'

)

allu <- allu %>%
  mutate(traj = factor(traj, levels = traj_levels))
# Add a RIGHT-SIDE "final" axis that shows trajectory % (and n), while keeping MAS→MAE % labels too.
# Drop this in AFTER you build `allu` and BEFORE `plot_one()`.

allu <- allu %>%
  mutate(
    traj_chr = as.character(traj),
    traj_abbr = dplyr::recode(
      traj_chr,
      "Increasing -> Increasing" = "Inc→Inc",
      "Decreasing -> Decreasing" = "Dec→Dec",
      "Increasing -> Decreasing" = "Inc→Dec",
      "Decreasing -> Increasing" = "Dec→Inc",
      "NoChange -> Increasing"   = "NC→Inc",
      "NoChange -> Decreasing"   = "NC→Dec",
      "Increasing -> NoChange"   = "Inc→NC",
      "Decreasing -> NoChange"   = "Dec→NC",
      "NoChange -> NoChange"     = "NC→NC",
      "NA -> Increasing"         = "NA→Inc",
      "NA -> Decreasing"         = "NA→Dec",
      "NA -> NoChange"           = "NA→NC",
      "Increasing -> NA"         = "Inc→NA",
      "Decreasing -> NA"         = "Dec→NA",
      "NoChange -> NA"           = "NC→NA",
      "NA -> NA"                 = "NA→NA",
      .default = traj_chr
    ),
    traj_abbr = factor(traj_abbr, levels = unique(dplyr::recode(
      traj_levels,
      "Increasing -> Increasing" = "Inc→Inc",
      "Decreasing -> Decreasing" = "Dec→Dec",
      "Increasing -> Decreasing" = "Inc→Dec",
      "Decreasing -> Increasing" = "Dec→Inc",
      "NoChange -> Increasing"   = "NC→Inc",
      "NoChange -> Decreasing"   = "NC→Dec",
      "Increasing -> NoChange"   = "Inc→NC",
      "Decreasing -> NoChange"   = "Dec→NC",
      "NoChange -> NoChange"     = "NC→NC",
      "NA -> Increasing"         = "NA→Inc",
      "NA -> Decreasing"         = "NA→Dec",
      "NA -> NoChange"           = "NA→NC",
      "Increasing -> NA"         = "Inc→NA",
      "Decreasing -> NA"         = "Dec→NA",
      "NoChange -> NA"           = "NC→NA",
      "NA -> NA"                 = "NA→NA",
      .default = traj_levels
    ))),
    traj_chr = NULL
  )

plot_one <- function(ct) {
  d <- allu %>% filter(celltype == ct)

  # legend labels with n and % (per celltype)
  leg_df <- d %>%
    group_by(traj) %>%
    summarise(n_traj = sum(n), .groups = "drop") %>%
    mutate(
      pct_traj = 100 * n_traj / sum(n_traj),
      lab = paste0(as.character(traj), " (n=", n_traj, ", ", sprintf("%.1f", pct_traj), "%)")
    )
  traj_labs <- setNames(leg_df$lab, as.character(leg_df$traj))

  # celltype-specific pooled threshold for title (from df3)
  thr_ct <- unique(df3$thr[df3$celltype == ct])
  thr_ct <- if (length(thr_ct) == 1) thr_ct else NA_real_

  ggplot(d, aes(axis1 = YS_to_MAS, axis2 = MAS_to_MAE, axis3 = traj_abbr, y = n)) +
    geom_alluvium(aes(fill = traj), width = 0.22, alpha = 0.90) +
    geom_stratum(width = 0.22, colour = "black", linewidth = 0.3, fill = "white") +
    geom_text(
      stat = "stratum",
      aes(label = ifelse(
        after_stat(x) == 3,
        paste0(sprintf("%.1f%%", after_stat(prop) * 100), "\n(n=", after_stat(count), ")"),
        paste0(after_stat(stratum), "\n", sprintf("%.1f%%", after_stat(prop) * 100))
      )),
      size = 2.7,
      lineheight = 0.95
    ) +
    scale_x_discrete(limits = c("YS->MAS", "MAS->MAE", "Final (trajectory)"), expand = c(0.05, 0.05)) +
    scale_fill_manual(
      values = traj_cols,
      breaks = names(traj_labs),
      labels = unname(traj_labs),
      drop = FALSE
    ) +
    labs(
      title = paste0("Alluvial ",ct," thr=", sprintf("%.3f", thr_ct), "): ", ct),
      x = NULL, y = "Count of sites", fill = "Trajectory"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.text = element_text(size = 9)
    )
}

plots <- lapply(sort(unique(allu$celltype)), plot_one)
names(plots) <- sort(unique(allu$celltype))
for (ct in names(plots)) print(plots[[ct]])

 out_pdf <- paste0("NFic_control_alluvial_direction_counts_split", thresh, ".pdf")
pdf(out_pdf, width = 7, height = 5, onefile = TRUE)
 for (ct in names(plots)) print(plots[[ct]])
 dev.off()
message("Saved: ", normalizePath(out_pdf))

# Loop reversal logic across ALL motifs, compute reversal proportions, and test which motifs are most reversal-enriched
# Assumes TOBIAS has: motif, cell (or celltype), YS_score, MAS_score, MAE_score

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ----------------------------
# Load
# ----------------------------
out_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2023-10-31.Cicero/2026-01-20_Recreated_Links")
TOBIAS <- readRDS(file.path(out_dir, "/Tobias_df_Classified_Full_scores.rds"))

stopifnot(all(c("motif","cell","YS_score","MAS_score","MAE_score") %in% colnames(TOBIAS)))

# ----------------------------
# Params
# ----------------------------
q_abs <- 0.25      # per-celltype threshold quantile using pooled |delta| (within each motif)
min_thr <- 0.0     # hard floor on threshold
min_n_for_test <- 50  # minimum n (sites changing in BOTH steps) to include motif in test
top_n_plot <- 40

# ----------------------------
# Prep deltas
# ----------------------------
df_all <- TOBIAS %>%
  transmute(
    motif = as.character(motif),
    celltype = as.character(cell),
    YS_score  = suppressWarnings(as.numeric(YS_score)),
    MAS_score = suppressWarnings(as.numeric(MAS_score)),
    MAE_score = suppressWarnings(as.numeric(MAE_score))
  ) %>%
  filter(!is.na(motif), !is.na(celltype)) %>%
  mutate(
    d_YS_MAS  = MAS_score - YS_score,      # ageing step
    d_MAS_MAE = MAE_score - MAS_score,     # exercise step
    abs_d_YS_MAS  = abs(d_YS_MAS),
    abs_d_MAS_MAE = abs(d_MAS_MAE)
  )

# ----------------------------
# Per (motif x celltype) threshold from pooled |delta| across both steps
# ----------------------------
thr_tbl <- df_all %>%
  group_by(motif, celltype) %>%
  summarise(
    thr = max(
      min_thr,
      as.numeric(quantile(c(abs_d_YS_MAS, abs_d_MAS_MAE), probs = q_abs, na.rm = TRUE))
    ),
    n_sites = n(),
    .groups = "drop"
  )

df3 <- df_all %>%
  left_join(thr_tbl, by = c("motif","celltype")) %>%
  mutate(
    chg_age = is.finite(d_YS_MAS)  & is.finite(thr) & (abs(d_YS_MAS)  >= thr),
    chg_ex  = is.finite(d_MAS_MAE) & is.finite(thr) & (abs(d_MAS_MAE) >= thr),

    inc_age = is.finite(d_YS_MAS)  & is.finite(thr) & (d_YS_MAS  >=  thr),
    dec_age = is.finite(d_YS_MAS)  & is.finite(thr) & (d_YS_MAS  <= -thr),
    inc_ex  = is.finite(d_MAS_MAE) & is.finite(thr) & (d_MAS_MAE >=  thr),
    dec_ex  = is.finite(d_MAS_MAE) & is.finite(thr) & (d_MAS_MAE <= -thr),

    rev_strict = (inc_age & dec_ex) | (dec_age & inc_ex),
    same_strict = (inc_age & inc_ex) | (dec_age & dec_ex)
  )

# ----------------------------
# Summarise per motif
# Denominator for "reversal proportion" defaults to sites that change in BOTH steps (more meaningful)
# ----------------------------
motif_tbl <- df3 %>%
  group_by(motif) %>%
  summarise(
    n_total = n(),
    n_age_change  = sum(chg_age, na.rm = TRUE),
    n_ex_change   = sum(chg_ex, na.rm = TRUE),
    n_both_change = sum(chg_age & chg_ex, na.rm = TRUE),

    n_rev  = sum(rev_strict, na.rm = TRUE),
    n_same = sum(same_strict, na.rm = TRUE),

    prop_rev_total      = n_rev / n_total,
    prop_rev_given_both = ifelse(n_both_change > 0, n_rev / n_both_change, NA_real_),
    prop_same_given_both = ifelse(n_both_change > 0, n_same / n_both_change, NA_real_),
    .groups = "drop"
  )

# ----------------------------
# Enrichment test: motifs with more reversal than global background
# Uses prop.test (fast). Alternative = greater (reversal-enriched).
# Testing on n_both_change (sites dynamic in BOTH steps).
# ----------------------------
test_tbl <- motif_tbl %>%
  filter(is.finite(prop_rev_given_both), n_total >= min_n_for_test)

p0 <- with(test_tbl, sum(n_rev) / sum(n_total))  # global reversal rate among "both-change" sites

pvals <- mapply(
  function(k, n) stats::prop.test(x = k, n = n, p = p0, alternative = "greater", correct = FALSE)$p.value,
  test_tbl$n_rev, test_tbl$n_total
)

test_tbl <- test_tbl %>%
  mutate(
    p = pvals,
    fdr = p.adjust(p, method = "BH"),
    enrich = prop_rev_given_both / p0,
    log2_enrich = log2(enrich),
    delta = prop_rev_given_both - p0
  ) %>%
  arrange(fdr, desc(prop_rev_given_both))

print(head(test_tbl, 30))
message("Global reversal rate (among both-change sites): ", sprintf("%.4f", p0))

# ----------------------------
# Save table
# ----------------------------
out_tsv <- paste0("motif_reversal_summary_q", q_abs, "_minN", min_n_for_test, ".tsv")
readr::write_tsv(test_tbl, out_tsv)
message("Saved: ", normalizePath(out_tsv))

# ----------------------------
# Plot: top motifs by reversal proportion (among both-change sites)
# ----------------------------
plot_df <- test_tbl %>%
  slice_head(n = top_n_plot) %>%
  mutate(motif = factor(motif, levels = rev(motif)))

p_top <- ggplot(plot_df, aes(x = motif, y = prop_rev_given_both)) +
  geom_col() +
  geom_hline(yintercept = p0, linetype = "dashed") +
  coord_flip() +
  labs(
    title = paste0("Top motifs by reversal proportion (q=", q_abs, ", min both-change n=", min_n_for_test, ")"),
    subtitle = paste0("Dashed = global reversal rate (", sprintf("%.3f", p0), "). Tested with prop.test (greater), BH FDR."),
    x = NULL,
    y = "Reversal proportion (n_rev / n_total)"
  ) +
  theme_classic()

print(p_top)

out_pdf <- paste0("top_motif_reversal_barplot_q", q_abs, "_minN", min_n_for_test, ".pdf")
ggsave(out_pdf, p_top, width = 7.5, height = 9)
message("Saved: ", normalizePath(out_pdf))

# Test motifs enriched for AgeUp+ExDown relative to AgeDown+ExUp
# Age step: YS -> MAS  (d_YS_MAS)
# Ex  step: MAS -> MAE (d_MAS_MAE)
# Uses thresholding on |delta| to call Inc/Dec per step, then compares the TWO reversal directions:
#   AUED = AgeUp & ExDown
#   ADED = AgeDown & ExUp
#
# Output: per-motif counts, prop_AUED, Fisher OR (motif vs others) and BH FDR.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ----------------------------
# Load
# ----------------------------
out_dir <- publication_path("HEART_SNRNA_ROOT", "analysis/2023-10-31.Cicero/2026-01-20_Recreated_Links")
TOBIAS <- readRDS(file.path(out_dir, "Tobias_df_Classified_Full_scores.rds"))

stopifnot(all(c("motif","cell","YS_score","MAS_score","MAE_score") %in% colnames(TOBIAS)))

# ----------------------------
# Params
# ----------------------------
q_abs <- 0.25
min_thr <- 0.0
min_n_rev <- 50              # minimum reversal-direction sites (AUED+ADED) per motif to test
thr_mode <- "celltype"       # "celltype" or "motif_celltype"
out_prefix <- paste0("AUED_enrichment_q", q_abs, "_", thr_mode, "_minN", min_n_rev)

# ----------------------------
# Prep deltas
# ----------------------------
df_all <- TOBIAS %>%
  transmute(
    motif = as.character(motif),
    celltype = as.character(cell),
    YS_score  = suppressWarnings(as.numeric(YS_score)),
    MAS_score = suppressWarnings(as.numeric(MAS_score)),
    MAE_score = suppressWarnings(as.numeric(MAE_score))
  ) %>%
  mutate(
    d_YS_MAS  = MAS_score - YS_score,
    d_MAS_MAE = MAE_score - MAS_score,
    abs_d_YS_MAS  = abs(d_YS_MAS),
    abs_d_MAS_MAE = abs(d_MAS_MAE)
  ) %>%
  filter(!is.na(motif), !is.na(celltype))

# ----------------------------
# Threshold table
# ----------------------------
if (thr_mode == "celltype") {
  thr_tbl <- df_all %>%
    group_by(celltype) %>%
    summarise(
      thr = max(min_thr, as.numeric(quantile(c(abs_d_YS_MAS, abs_d_MAS_MAE), probs = q_abs, na.rm = TRUE))),
      .groups = "drop"
    )
  df_thr <- df_all %>% left_join(thr_tbl, by = "celltype")
} else if (thr_mode == "motif_celltype") {
  thr_tbl <- df_all %>%
    group_by(motif, celltype) %>%
    summarise(
      thr = max(min_thr, as.numeric(quantile(c(abs_d_YS_MAS, abs_d_MAS_MAE), probs = q_abs, na.rm = TRUE))),
      .groups = "drop"
    )
  df_thr <- df_all %>% left_join(thr_tbl, by = c("motif","celltype"))
} else {
  stop("thr_mode must be 'celltype' or 'motif_celltype'")
}

# ----------------------------
# Call direction per step and define the two reversal directions
# ----------------------------
df_calls <- df_thr %>%
  mutate(
    inc_age = is.finite(d_YS_MAS)  & is.finite(thr) & (d_YS_MAS  >=  thr),
    dec_age = is.finite(d_YS_MAS)  & is.finite(thr) & (d_YS_MAS  <= -thr),
    inc_ex  = is.finite(d_MAS_MAE) & is.finite(thr) & (d_MAS_MAE >=  thr),
    dec_ex  = is.finite(d_MAS_MAE) & is.finite(thr) & (d_MAS_MAE <= -thr),

    AUED = inc_age & dec_ex,   # Age Up, Ex Down
    ADED = dec_age & inc_ex    # Age Down, Ex Up
  )

# Universe restricted to the two reversal directions only (ignores same-direction and no-change)
rev_univ <- df_calls %>%
  filter(AUED | ADED) %>%
  transmute(motif, AUED = AUED, ADED = ADED)

# Per-motif counts
motif_cnt <- rev_univ %>%
  group_by(motif) %>%
  summarise(
    n_AUED = sum(AUED, na.rm = TRUE),
    n_ADED = sum(ADED, na.rm = TRUE),
    n_rev_dir = n_AUED + n_ADED,
    prop_AUED = ifelse(n_rev_dir > 0, n_AUED / n_rev_dir, NA_real_),
    .groups = "drop"
  ) %>%
  filter(n_rev_dir >= min_n_rev)

# Global AUED fraction within this universe (used for binomial baseline)
p0 <- with(motif_cnt, sum(n_AUED) / sum(n_rev_dir))
message("Global AUED fraction (within AUED+ADED universe): ", sprintf("%.4f", p0))

# Totals for "others" in Fisher tests
TOT_AUED <- sum(motif_cnt$n_AUED)
TOT_ADED <- sum(motif_cnt$n_ADED)

# ----------------------------
# Tests
# 1) Fisher: is AUED:ADED ratio higher in motif than in all other motifs
# 2) Binomial: is prop_AUED higher than global p0 (similar question, different framing)
# ----------------------------
test_tbl <- motif_cnt %>%
  rowwise() %>%
  mutate(
    other_AUED = TOT_AUED - n_AUED,
    other_ADED = TOT_ADED - n_ADED,

    fisher = list(fisher.test(matrix(c(n_AUED, n_ADED, other_AUED, other_ADED), nrow = 2), alternative = "greater")),
    OR = unname(fisher$estimate),
    p_fisher = fisher$p.value,

    p_binom = binom.test(x = n_AUED, n = n_rev_dir, p = p0, alternative = "greater")$p.value
  ) %>%
  ungroup() %>%
  mutate(
    fdr_fisher = p.adjust(p_fisher, method = "BH"),
    fdr_binom  = p.adjust(p_binom,  method = "BH"),
    log2OR = log2(OR),
    delta = prop_AUED - p0
  ) %>%
  arrange(fdr_fisher, desc(prop_AUED))

print(head(test_tbl, 30))

# ----------------------------
# Save table
# ----------------------------
out_tsv <- paste0(out_prefix, "_motif_tests.tsv")
write.table(test_tbl, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message("Saved: ", normalizePath(out_tsv))

# ----------------------------
# Plot: top motifs enriched for AUED (by Fisher FDR then effect size)
# ----------------------------
top_n <- 40
plot_df <- test_tbl %>%
  arrange(fdr_fisher, desc(delta), desc(log2OR)) %>%
  slice_head(n = top_n) %>%
  mutate(motif = factor(motif, levels = rev(motif)),
         sig = ifelse(fdr_fisher < 0.05, "FDR<0.05", "ns"))

p <- ggplot(plot_df, aes(x = motif, y = prop_AUED)) +
  geom_col() +
  geom_hline(yintercept = p0, linetype = "dashed") +
  coord_flip() +
  labs(
    title = "Motifs enriched for AgeUp + ExDown relative to AgeDown + ExUp",
    subtitle = paste0("Universe = AUED+ADED only. Dashed = global AUED fraction (", sprintf("%.3f", p0), "). ",
                      "Test = Fisher (motif vs others), BH FDR."),
    x = NULL,
    y = "prop_AUED = n_AUED / (n_AUED + n_ADED)"
  ) +
  theme_classic()

print(p)

out_pdf <- paste0(out_prefix, "_top", top_n, "_barplot.pdf")
ggsave(out_pdf, p, width = 8.5, height = 10)
message("Saved: ", normalizePath(out_pdf))
