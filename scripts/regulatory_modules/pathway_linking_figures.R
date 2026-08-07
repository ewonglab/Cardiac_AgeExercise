#!/usr/bin/env Rscript
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Link AP-1 footprints to ECM and IL-1 pathway modules.
# Steps: 1) Load TOBIAS footprints and reference TSSs; 2) assign nearby genes/modules; 3) export regression and exercise figures.
# =======================================
# Source reformatted from emily_plots/pathway_plots/22_fig5CE_replot.R.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))
# ==============================================================================
# 22_fig5CE_replot.R
# Panel C: Regression (ECM, IL-1 only) — all fibroblast AP-1 sites
# Panel E: Exercise effect dotplot — all 19 modules, all cell types, clustered
# ==============================================================================

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(GenomicRanges)

setwd(publication_path("PEAK2PATHWAY_ROOT"))
out_dir <- publication_path("PAPER_OUTPUT_ROOT", "regulatory_modules")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. Load and annotate
# ==============================================================================

cat("Loading TOBIAS data...\n")
tob <- as.data.table(readRDS("Tobias_df_Classified_Full.rds"))
ap1_members <- c("Fos", "Fosb", "Fosl1", "Fosl2", "Jun", "Junb", "Jund")
ap1 <- tob[name %in% ap1_members]
rm(tob); gc()

cat("Annotating nearest gene...\n")
tss <- readRDS("mm10_refGene_tss.rds")
site_gr <- GRanges(seqnames = ap1$seqnames,
                    ranges = IRanges(start = ap1$start, end = ap1$end))
ni <- nearest(site_gr, tss)
ap1$nearest_gene <- mcols(tss)$gene_symbol[ni]
ap1$dist_to_tss <- distance(site_gr, tss[ni])

nearby <- ap1[dist_to_tss <= 200000]
ct_exclude <- c("MacFib", "HSPC", "M1", "M2")
nearby <- nearby[!CellType %in% ct_exclude]

nearby[, aging_delta := MAS_score - YS_score]
nearby[, exercise_delta := MAE_score - MAS_score]

# ==============================================================================
# 2. Gene modules (all 19)
# ==============================================================================

gene_modules <- list(
  ECM_Remodeling = c(
    "Col1a1", "Col1a2", "Col3a1", "Col4a1", "Col6a1", "Col6a3", "Fn1",
    "Lama2", "Lama4", "Mmp2", "Mmp3", "Mmp9", "Mmp14", "Timp1", "Timp2",
    "Timp3", "Lox", "Loxl2", "Dcn", "Bgn", "Vcan", "Postn", "Thbs1",
    "Thbs2", "Sparc", "Tnc", "Tgfb1", "Tgfb2", "Tgfb3", "Tgfbr1",
    "Tgfbr2", "Adamts2"),
  IL1_Pathway = c(
    "Il1a", "Il1b", "Il18", "Il33", "Il36a", "Il36b", "Il36g", "Il36rn",
    "Il1r1", "Il1r2", "Il1rap", "Il1rl1", "Il1rl2", "Il18r1", "Il18rap", "Il1rn",
    "Myd88", "Irak1", "Irak2", "Irak4", "Traf6", "Map3k7", "Mapk14", "Mapk8",
    "Nfkb1", "Rela", "Fos", "Jun", "Junb",
    "Nlrp3", "Pycard", "Casp1", "Casp4",
    "Ptgs2", "Ccl2", "Cxcl1", "Cxcl2", "Mmp3", "Mmp9", "Icam1", "Vcam1", "Sele"),
  Inflammation = c(
    "Tnf", "Il6", "Il1b", "Il1a", "Il18", "Ccl2", "Ccl3", "Ccl4", "Ccl5",
    "Ccl7", "Cxcl1", "Cxcl2", "Cxcl10", "Cxcl12", "Ccr2", "Cxcr4",
    "Nfkb1", "Nfkb2", "Rela", "Ikbkb", "Nlrp3", "Pycard", "Casp1",
    "Ptgs2", "Nos2", "Mmp9", "Il10", "Tgfb1", "Tlr2", "Tlr4"),
  Mechanotransduction = c(
    "Itgb1", "Itgb3", "Itga5", "Itga7", "Itgav", "Vcl", "Tln1", "Pxn",
    "Ptk2", "Yap1", "Wwtr1", "Tead1", "Tead4", "Lats1", "Lats2", "Piezo1",
    "Piezo2", "Trpv4", "Rock1", "Rock2", "Rhoa", "Dmd", "Dag1", "Sgcd",
    "Ttn", "Cd44", "Sdc1", "Sdc4", "Actn1", "Flna"),
  Angiogenesis = c(
    "Vegfa", "Vegfb", "Vegfc", "Flt1", "Kdr", "Flt4", "Angpt1", "Angpt2",
    "Tek", "Pecam1", "Cdh5", "Nos3", "Eng", "Dll4", "Notch1", "Ednra",
    "Ednrb", "Edn1", "Ace", "Ace2", "Hif1a", "Epas1", "Nrp1", "Nrp2",
    "Pdgfb", "Pdgfrb", "Pgf"),
  Oxidative_Stress = c(
    "Sod1", "Sod2", "Sod3", "Cat", "Gpx1", "Gpx4", "Gsr", "Gclc", "Gclm",
    "Gstm1", "Txn1", "Txnrd1", "Prdx1", "Prdx2", "Nfe2l2", "Keap1",
    "Hmox1", "Nqo1", "Cybb", "Cyba", "Nox4", "Hif1a", "Epas1", "Egln1",
    "Egln3", "Bnip3", "Gadd45a", "Gadd45b", "Ldha", "Pdk1"),
  Mitochondrial = c(
    "Ppargc1a", "Ppargc1b", "Tfam", "Nrf1", "Ndufs1", "Ndufv1", "Sdha",
    "Sdhb", "Uqcrc1", "Cox5a", "Cox7a1", "Atp5f1a", "Dnm1l", "Fis1",
    "Mff", "Mfn1", "Mfn2", "Opa1", "Pink1", "Prkn", "Bnip3l", "Fundc1",
    "Hspd1", "Hspe1", "Lonp1", "Clpp", "Cpt1a", "Cpt1b", "Cpt2", "Acadl",
    "Acadm", "Cs", "Idh2"),
  Cardiac_Hypertrophy = c(
    "Nppa", "Nppb", "Myh7", "Acta1", "Gata4", "Nkx2-5", "Mef2a", "Mef2c",
    "Hand2", "Tbx5", "Myocd", "Srf", "Ankrd1", "Ctgf", "Tgfb1"),
  Senescence = c(
    "Cdkn1a", "Cdkn2a", "Tp53", "Rb1", "Serpine1", "Igfbp3", "Igfbp5",
    "Ccnd1", "Bmi1", "Hmga2", "Gadd45b", "Gadd45a", "Trp53bp1", "Pml"),
  Calcium_Handling = c(
    "Ryr2", "Atp2a2", "Pln", "Casq2", "Cacna1c", "Slc8a1", "Camk2d",
    "Calm1", "Calm2", "Calm3", "Trdn", "Jph2", "Hrc", "Mybpc3"),
  Fibrosis_Signaling = c(
    "Smad2", "Smad3", "Smad4", "Smad7", "Ctgf", "Pdgfa", "Pdgfb",
    "Pdgfra", "Pdgfrb", "Wnt3a", "Wnt5a", "Fzd1", "Fzd2", "Axin2",
    "Dkk1", "Wisp1"),
  MAPK_Signaling = c(
    "Map2k1", "Map2k2", "Mapk1", "Mapk3", "Mapk8", "Mapk9", "Mapk14",
    "Map3k1", "Map3k7", "Raf1", "Braf", "Sos1", "Grb2", "Kras", "Hras",
    "Dusp1", "Dusp4", "Dusp6", "Dusp10", "Spry1", "Spry2"),
  NFkB_Signaling = c(
    "Nfkb1", "Nfkb2", "Rela", "Relb", "Rel", "Ikbkb", "Ikbkg", "Ikbke",
    "Chuk", "Nfkbia", "Nfkbib", "Tnfaip3", "Traf2", "Traf6", "Ripk1",
    "Birc2", "Birc3", "Bcl3", "Card11"),
  TGFb_Signaling = c(
    "Tgfb1", "Tgfb2", "Tgfb3", "Tgfbr1", "Tgfbr2", "Smad2", "Smad3",
    "Smad4", "Smad6", "Smad7", "Acvrl1", "Acvr1", "Bmpr1a", "Bmpr2",
    "Bmp2", "Bmp4", "Bmp7", "Id1", "Id2", "Id3"),
  Apoptosis = c(
    "Bax", "Bak1", "Bcl2", "Bcl2l1", "Mcl1", "Bid", "Bad", "Bim",
    "Casp3", "Casp7", "Casp8", "Casp9", "Apaf1", "Cycs", "Diablo",
    "Fas", "Fasl", "Fadd", "Tradd", "Xiap"),
  Cell_Cycle = c(
    "Ccnd1", "Ccnd2", "Ccne1", "Ccna2", "Ccnb1", "Cdk1", "Cdk2", "Cdk4",
    "Cdk6", "Cdkn1a", "Cdkn1b", "Cdkn2a", "Rb1", "E2f1", "E2f3",
    "Myc", "Pcna", "Mki67", "Top2a"),
  EMT_EndMT = c(
    "Snai1", "Snai2", "Twist1", "Twist2", "Zeb1", "Zeb2", "Vim", "Fn1",
    "Cdh1", "Cdh2", "Cdh11", "Acta2", "Tagln", "S100a4", "Mmp2", "Mmp9",
    "Tgfb1", "Tgfbr1", "Wnt5a", "Notch1"),
  Autophagy = c(
    "Atg5", "Atg7", "Atg12", "Atg16l1", "Becn1", "Map1lc3b", "Map1lc3a",
    "Sqstm1", "Ulk1", "Ulk2", "Pik3c3", "Lamp1", "Lamp2", "Bnip3",
    "Bnip3l", "Tfeb", "Mtor", "Rptor", "Tsc1", "Tsc2"),
  Wnt_Signaling = c(
    "Wnt3a", "Wnt5a", "Wnt7a", "Wnt11", "Ctnnb1", "Lef1", "Tcf7",
    "Tcf7l2", "Axin2", "Apc", "Gsk3b", "Dvl1", "Dvl2", "Fzd1", "Fzd2",
    "Fzd7", "Sfrp1", "Sfrp2", "Dkk1", "Rspo1")
)

ct_labels <- c(
  Fibroblast = "Fibroblast", Cardiomyocyte = "CM",
  AllMacrophage = "Macrophage", Coronary_EC = "Coronary EC",
  Endocardial_EC = "Endocardial EC", Lymphatic_EC = "Lymphatic EC",
  Pericytes = "Pericyte"
)

mod_labels <- c(
  ECM_Remodeling = "ECM Remodeling", IL1_Pathway = "IL-1 Pathway",
  MAPK_Signaling = "MAPK", TGFb_Signaling = "TGFb",
  Inflammation = "Inflammation", Mechanotransduction = "Mechanotransduction",
  Senescence = "Senescence", Cardiac_Hypertrophy = "Cardiac Hypertrophy",
  Oxidative_Stress = "Oxidative Stress", Angiogenesis = "Angiogenesis",
  Mitochondrial = "Mitochondrial", Calcium_Handling = "Calcium Handling",
  Fibrosis_Signaling = "Fibrosis", NFkB_Signaling = "NF-kB",
  Apoptosis = "Apoptosis", Cell_Cycle = "Cell Cycle",
  EMT_EndMT = "EMT / EndMT", Autophagy = "Autophagy",
  Wnt_Signaling = "Wnt"
)

# ==============================================================================
# 3. Assign module per site (first match)
# ==============================================================================

cat("Assigning modules...\n")
mod_lookup <- stack(gene_modules)
setnames(setDT(mod_lookup), c("gene", "module"))
# For genes in multiple modules, keep first
mod_lookup <- mod_lookup[!duplicated(gene)]

nearby[, module := mod_lookup$module[match(nearest_gene, mod_lookup$gene)]]

cat(sprintf("  Sites with module: %d / %d\n",
            sum(!is.na(nearby$module)), nrow(nearby)))

# ==============================================================================
# 4. Panel C: Regression (ECM and IL-1 only, all fibroblast sites)
# ==============================================================================

cat("\n=== Panel C: Regression ===\n")

reg_data <- nearby[CellType == "Fibroblast" & !is.na(module) &
                   module %in% c("ECM_Remodeling", "IL1_Pathway")]
reg_data[, mod_disp := ifelse(module == "ECM_Remodeling", "ECM", "IL-1")]
reg_data$mod_disp <- factor(reg_data$mod_disp, levels = c("ECM", "IL-1"))

reg_stats <- reg_data[, .(
  n = .N,
  slope = coef(lm(exercise_delta ~ aging_delta))[2],
  p = summary(lm(exercise_delta ~ aging_delta))$coefficients[2, 4]
), by = mod_disp]
reg_stats[, label := sprintf("slope=%.2f, P=%.1e, n=%d", slope, p, n)]

cat("Regression stats:\n")
print(reg_stats)

mod_colours <- c("ECM" = "#D9A05B", "IL-1" = "#8073AC")

# Separate plots for ECM and IL-1
for (mod in c("ECM", "IL-1")) {
  mod_data <- reg_data[mod_disp == mod]
  mod_stats <- reg_stats[mod_disp == mod]
  mod_col <- mod_colours[mod]

  p <- ggplot(mod_data, aes(x = aging_delta, y = exercise_delta)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_point(size = 0.8, alpha = 0.5, colour = mod_col) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, alpha = 0.15,
                colour = mod_col, fill = mod_col) +
    geom_text(data = mod_stats, aes(label = label),
              x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
              size = 2.2, colour = mod_col) +
    labs(x = "Aging delta (MAS \u2013 YS)",
         y = "Exercise delta (MAE \u2013 MAS)",
         title = paste0(mod, ": exercise correction\nscales with aging effect")) +
    theme_minimal(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", size = 8),
      axis.text = element_text(size = 7),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(5, 5, 5, 5)
    )

  mod_slug <- gsub("-", "", mod)
  ggsave(file.path(out_dir, paste0("fig5C_regression_", mod_slug, ".pdf")),
         p, width = 2.8, height = 2.8, bg = "white")
  ggsave(file.path(out_dir, paste0("fig5C_regression_", mod_slug, ".png")),
         p, width = 2.8, height = 2.8, dpi = 300, bg = "white")
  cat(sprintf("Saved fig5C_regression_%s.pdf/png\n", mod_slug))
}

# ==============================================================================
# 5. Panel E: Exercise effect dotplot — all 19 modules, all cell types
# ==============================================================================

cat("\n=== Panel E: Exercise dotplot (19 modules) ===\n")

# Z-score exercise delta within each cell type to normalize for score scale differences
nearby[, exercise_delta_z := (exercise_delta - mean(exercise_delta)) / sd(exercise_delta),
       by = CellType]

celltypes <- sort(unique(nearby$CellType))
results_list <- list()

for (ct in celltypes) {
  ct_data <- nearby[CellType == ct]

  for (mod_name in names(gene_modules)) {
    in_mod <- ct_data$nearest_gene %in% gene_modules[[mod_name]]
    n_mod <- sum(in_mod)
    if (n_mod < 3) {
      results_list[[length(results_list) + 1]] <- data.frame(
        cell_type = ct, module = mod_name, n_sites = n_mod,
        mean_exercise = NA, mean_exercise_z = NA, p_exercise = NA)
      next
    }

    mod_vals <- ct_data$exercise_delta[in_mod]
    mod_vals_z <- ct_data$exercise_delta_z[in_mod]

    # Wilcoxon signed-rank test: is exercise delta != 0?
    wt <- wilcox.test(mod_vals, mu = 0)

    results_list[[length(results_list) + 1]] <- data.frame(
      cell_type = ct, module = mod_name, n_sites = n_mod,
      mean_exercise = mean(mod_vals), mean_exercise_z = mean(mod_vals_z),
      p_exercise = wt$p.value)
  }
}

ex_res <- do.call(rbind, results_list)
rownames(ex_res) <- NULL

# BH correction within each cell type
ex_res <- ex_res %>%
  group_by(cell_type) %>%
  mutate(p_adj = p.adjust(p_exercise, method = "BH")) %>%
  ungroup()

write.csv(ex_res, file.path(out_dir, "exercise_effect_19modules.csv"), row.names = FALSE)
cat("Saved exercise_effect_19modules.csv\n")

# Print significant
cat("\n=== Significant exercise effects (P < 0.05) ===\n")
sig <- ex_res %>% filter(p_exercise < 0.05) %>% arrange(p_exercise)
print(as.data.frame(sig), row.names = FALSE)

# --- Heatmap with significance stars ---

plot_data <- ex_res %>%
  filter(!is.na(p_exercise)) %>%
  mutate(
    ct_disp = ifelse(cell_type %in% names(ct_labels), ct_labels[cell_type], cell_type),
    mod_disp = ifelse(module %in% names(mod_labels), mod_labels[module], module),
    ex_cap = pmax(pmin(mean_exercise_z, 0.4), -0.4),
    sig_star = case_when(
      p_adj < 0.05  ~ "***",
      p_exercise < 0.01 ~ "**",
      p_exercise < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# Cluster on z-scored exercise delta values
clust_mat <- plot_data %>%
  select(mod_disp, ct_disp, mean_exercise_z) %>%
  mutate(mean_exercise_z = replace_na(mean_exercise_z, 0)) %>%
  pivot_wider(names_from = ct_disp, values_from = mean_exercise_z, values_fill = 0) %>%
  as.data.frame()
rownames(clust_mat) <- clust_mat$mod_disp
clust_mat$mod_disp <- NULL
mat <- as.matrix(clust_mat)

row_clust <- hclust(dist(mat), method = "ward.D2")
col_clust <- hclust(dist(t(mat)), method = "ward.D2")

# Remove modules with no significance in any cell type (+ Fibrosis)
mods_with_sig <- plot_data %>%
  filter(sig_star != "") %>%
  pull(mod_disp) %>%
  unique()
mods_with_sig <- setdiff(mods_with_sig, "Fibrosis")

plot_filt <- plot_data %>% filter(mod_disp %in% mods_with_sig)

# Re-cluster on filtered data
clust_mat2 <- plot_filt %>%
  select(mod_disp, ct_disp, mean_exercise_z) %>%
  mutate(mean_exercise_z = replace_na(mean_exercise_z, 0)) %>%
  pivot_wider(names_from = ct_disp, values_from = mean_exercise_z, values_fill = 0) %>%
  as.data.frame()
rownames(clust_mat2) <- clust_mat2$mod_disp
clust_mat2$mod_disp <- NULL
mat2 <- as.matrix(clust_mat2)

row_clust2 <- hclust(dist(mat2), method = "ward.D2")
col_clust2 <- hclust(dist(t(mat2)), method = "ward.D2")

plot_filt$mod_disp <- factor(plot_filt$mod_disp, levels = rev(rownames(mat2)[row_clust2$order]))
plot_filt$ct_disp  <- factor(plot_filt$ct_disp, levels = colnames(mat2)[col_clust2$order])

# --- Muted terracotta/green color scale (z-scored) ---
muted_scale <- scale_fill_gradient2(
  low = "#7A9E7E", mid = "#FAF8F5", high = "#C08B70", midpoint = 0,
  limits = c(-0.4, 0.4), name = "Exercise\ndelta (z)",
  guide = guide_colorbar(barwidth = 0.5, barheight = 3)
)

# --- Plot 1: Full heatmap (rows with any significance, no Fibrosis) ---

p_e1 <- ggplot(plot_filt, aes(x = ct_disp, y = mod_disp)) +
  geom_tile(aes(fill = ex_cap), colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sig_star), size = 2.5, colour = "black", vjust = 0.75) +
  muted_scale +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       title = "Exercise effect on AP-1 footprints",
       caption = "* P<0.05   ** P<0.01   *** FDR<0.05") +
  theme_minimal(base_size = 8) +
  theme(
    plot.title = element_text(face = "bold", size = 9),
    plot.caption = element_text(size = 5.5, colour = "grey40", hjust = 0),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 6.5),
    panel.grid = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.key.size = unit(0.35, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )

ggsave(file.path(out_dir, "fig5E_exercise_by_celltype.pdf"),
       p_e1, width = 3.3, height = 3.5, bg = "white")
ggsave(file.path(out_dir, "fig5E_exercise_by_celltype.png"),
       p_e1, width = 3.3, height = 3.5, dpi = 300, bg = "white")
cat("Saved fig5E_exercise_by_celltype.pdf/png\n")

# --- Plot 2: Significant-only (blank if ns) ---

plot_sig <- plot_filt %>%
  mutate(
    ex_cap_sig = ifelse(sig_star == "", NA_real_, ex_cap)
  )

p_e2 <- ggplot(plot_sig, aes(x = ct_disp, y = mod_disp)) +
  geom_tile(aes(fill = ex_cap_sig), colour = "grey85", linewidth = 0.3) +
  geom_text(aes(label = sig_star), size = 2.5, colour = "black", vjust = 0.75) +
  scale_fill_gradient2(
    low = "#7A9E7E", mid = "#FAF8F5", high = "#C08B70", midpoint = 0,
    limits = c(-0.4, 0.4), name = "Exercise\ndelta (z)",
    na.value = "white",
    guide = guide_colorbar(barwidth = 0.5, barheight = 3)
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       title = "Exercise effect on AP-1 footprints",
       subtitle = "Significant effects only",
       caption = "* P<0.05   ** P<0.01   *** FDR<0.05") +
  theme_minimal(base_size = 8) +
  theme(
    plot.title = element_text(face = "bold", size = 9),
    plot.subtitle = element_text(size = 7, colour = "grey40"),
    plot.caption = element_text(size = 5.5, colour = "grey40", hjust = 0),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 6.5),
    panel.grid = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.key.size = unit(0.35, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )

ggsave(file.path(out_dir, "fig5E_exercise_sig_only.pdf"),
       p_e2, width = 3.3, height = 3.5, bg = "white")
ggsave(file.path(out_dir, "fig5E_exercise_sig_only.png"),
       p_e2, width = 3.3, height = 3.5, dpi = 300, bg = "white")
cat("Saved fig5E_exercise_sig_only.pdf/png\n")

cat("Done.\n")
