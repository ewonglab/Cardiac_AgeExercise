#!/usr/bin/env Rscript
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Run log-fold-change-ranked Hallmark GSEA for the anti-IL-1β RNA contrast.
# Steps: 1) Load the anti-IL-1β edgeR object; 2) construct ranked genes and Hallmark sets; 3) export GSEA results, plots and a method summary.
# =======================================
# Source reformatted from anti-il1code/run_hallmark_gsea_logfc_rank.R.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(clusterProfiler)
  library(msigdbr)
})

project_dir <- publication_path("H3K27AC_PROJECT_ROOT")
input_rds <- publication_path("ANTI_IL1B_RNA_EDGER_OBJECT")
out_dir <- file.path(project_dir, "results", "rnaseq_hallmark_gsea_logfc_rank")
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260728)

get_msigdbr_table <- function(collection, subcollection = NULL) {
  tryCatch(
    msigdbr::msigdbr(
      species = "Mus musculus",
      db_species = "MM",
      collection = collection,
      subcollection = subcollection
    ),
    error = function(e) {
      msigdbr::msigdbr(
        species = "Mus musculus",
        category = collection,
        subcategory = subcollection
      )
    }
  )
}

clean_hallmark_label <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  x <- gsub("\\bTnfa\\b", "TNFA", x)
  x <- gsub("\\bNfkb\\b", "NFKB", x)
  x <- gsub("\\bIl6\\b", "IL6", x)
  x <- gsub("\\bJak\\b", "JAK", x)
  x <- gsub("\\bStat3\\b", "STAT3", x)
  x <- gsub("\\bIl2\\b", "IL2", x)
  x <- gsub("\\bStat5\\b", "STAT5", x)
  x <- gsub("\\bG2m\\b", "G2M", x)
  x <- gsub("\\bE2f\\b", "E2F", x)
  x <- gsub("\\bKras\\b", "KRAS", x)
  x <- gsub("\\bDna\\b", "DNA", x)
  x
}

compact_hallmark_label <- function(id, clean_term) {
  out <- clean_term
  out[id == "HALLMARK_OXIDATIVE_PHOSPHORYLATION"] <- "OXPHOS"
  out[id == "HALLMARK_FATTY_ACID_METABOLISM"] <- "Fatty acid metab."
  out[id == "HALLMARK_BILE_ACID_METABOLISM"] <- "Bile acid metab."
  out[id == "HALLMARK_INTERFERON_GAMMA_RESPONSE"] <- "IFN-gamma response"
  out[id == "HALLMARK_INTERFERON_ALPHA_RESPONSE"] <- "IFN-alpha response"
  out[id == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"] <- "EMT"
  out[id == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"] <- "TNFA/NFKB signaling"
  out[id == "HALLMARK_IL6_JAK_STAT3_SIGNALING"] <- "IL6/JAK/STAT3"
  out[id == "HALLMARK_G2M_CHECKPOINT"] <- "G2M checkpoint"
  out[id == "HALLMARK_E2F_TARGETS"] <- "E2F targets"
  out[id == "HALLMARK_KRAS_SIGNALING_UP"] <- "KRAS signaling up"
  out
}

message("Reading edgeR object")
obj <- readRDS(input_rds)
de <- as.data.table(obj$de)
required_cols <- c("gene_id", "gene_symbol", "logFC")
missing_cols <- setdiff(required_cols, names(de))
if (length(missing_cols) > 0) {
  stop("DE table is missing required column(s): ", paste(missing_cols, collapse = ", "))
}

rank_dt <- de[
  !is.na(gene_symbol) &
    gene_symbol != "" &
    gene_symbol != "NA" &
    is.finite(logFC),
  .(
    gene_id = as.character(gene_id),
    gene_symbol = as.character(gene_symbol),
    logFC = as.numeric(logFC)
  )
]
rank_dt[, abs_logFC := abs(logFC)]

# If a symbol appears more than once, keep the transcript/gene row with the
# largest absolute logFC. P values are deliberately not used in the rank.
setorder(rank_dt, gene_symbol, -abs_logFC, gene_id)
rank_dt <- rank_dt[!duplicated(gene_symbol)]
gene_list <- rank_dt$logFC
names(gene_list) <- rank_dt$gene_symbol
gene_list <- sort(gene_list, decreasing = TRUE)

fwrite(rank_dt[order(-logFC)], file.path(out_dir, "ranked_gene_list_symbol_logFC_table.csv"))
write.table(
  data.table(gene = names(gene_list), rank = as.numeric(gene_list)),
  file.path(out_dir, "ranked_gene_list_symbol_logFC.rnk"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

message("Preparing mouse Hallmark gene sets")
msig_h <- as.data.table(get_msigdbr_table("MH"))
term2gene <- unique(msig_h[!is.na(gene_symbol) & gene_symbol != "", .(
  term = as.character(gs_name),
  gene = as.character(gene_symbol)
)])
term2name <- unique(msig_h[, .(
  term = as.character(gs_name),
  name = as.character(gs_name)
)])

message("Running Hallmark GSEA with logFC rank")
gsea_obj <- clusterProfiler::GSEA(
  geneList = gene_list,
  TERM2GENE = term2gene,
  TERM2NAME = term2name,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE,
  seed = TRUE,
  by = "fgsea"
)

res <- as.data.table(as.data.frame(gsea_obj))
if (nrow(res) > 0) {
  res[, `:=`(
    collection = "MSigDB_Hallmark",
    rank_method = "logFC",
    rank_description = "edgeR log2 fold-change only",
    absNES = abs(NES),
    clean_term = clean_hallmark_label(Description),
    direction = fifelse(NES >= 0, "Higher in anti-IL-1b", "Lower in anti-IL-1b")
  )]
  res[, compact_term := compact_hallmark_label(ID, clean_term)]
  setorder(res, p.adjust, pvalue, -absNES)
}

saveRDS(gsea_obj, file.path(out_dir, "MSigDB_Hallmark_logFC_GSEA_result.rds"))
fwrite(res, file.path(out_dir, "MSigDB_Hallmark_logFC_GSEA_results.csv"))

message("Writing plots")
fill_values <- c(
  "Lower in anti-IL-1b" = "#3E6FA3",
  "Higher in anti-IL-1b" = "#C85050"
)

sig_counts <- res[, .N, by = .(direction, significant = p.adjust < 0.05)]
fwrite(sig_counts, file.path(out_dir, "MSigDB_Hallmark_logFC_GSEA_significance_counts.csv"))

down <- res[NES < 0][order(p.adjust, pvalue, NES)][seq_len(min(.N, 12))]
up_sig <- res[NES > 0 & p.adjust < 0.05][order(p.adjust, pvalue, -NES)]
plot_dt <- rbindlist(list(down, up_sig), fill = TRUE)
plot_dt[, compact_term := factor(compact_term, levels = plot_dt[order(NES), compact_term])]
plot_dt[, direction := factor(direction, levels = c("Lower in anti-IL-1b", "Higher in anti-IL-1b"))]

p_compact <- ggplot(plot_dt, aes(x = NES, y = compact_term, fill = direction)) +
  geom_vline(xintercept = 0, linewidth = 0.22, color = "grey45") +
  geom_col(width = 0.7, color = NA) +
  scale_fill_manual(values = fill_values, drop = FALSE) +
  scale_x_continuous(
    breaks = c(-2, -1, 0, 1),
    limits = c(-2.55, 1.85),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(x = "NES", y = NULL, fill = NULL) +
  theme_classic(base_size = 6.8) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 5.7, color = "black", margin = margin(r = 1.2)),
    axis.text.x = element_text(size = 5.9, color = "black"),
    axis.title.x = element_text(size = 6.5, color = "black", margin = margin(t = 1.5)),
    axis.line.x = element_line(linewidth = 0.22, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.22, color = "black"),
    legend.position = "top",
    legend.justification = "left",
    legend.key.size = unit(2.6, "mm"),
    legend.text = element_text(size = 5.7),
    legend.margin = margin(0, 0, 0.5, 0),
    legend.box.margin = margin(0, 0, 0, -2),
    plot.margin = margin(1.5, 2.2, 1.5, 1.5)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ggsave(
  file.path(plot_dir, "hallmark_logFC_rank_top_terms_ultracompact_manuscript.png"),
  p_compact,
  width = 2.65,
  height = 2.75,
  dpi = 600,
  bg = "white"
)
ggsave(
  file.path(plot_dir, "hallmark_logFC_rank_top_terms_ultracompact_manuscript.pdf"),
  p_compact,
  width = 2.65,
  height = 2.75,
  bg = "white",
  useDingbats = FALSE
)

plot_dt[, clean_term := factor(clean_term, levels = plot_dt[order(NES), clean_term])]
p_full_labels <- ggplot(plot_dt, aes(x = NES, y = clean_term, fill = direction)) +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "grey45") +
  geom_col(width = 0.68, color = NA) +
  scale_fill_manual(values = fill_values, drop = FALSE) +
  scale_x_continuous(
    breaks = c(-2, -1, 0, 1),
    limits = c(-2.55, 1.85),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(x = "NES", y = NULL, fill = NULL) +
  theme_classic(base_size = 7.2) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 6.1, color = "black", margin = margin(r = 1.5)),
    axis.text.x = element_text(size = 6.4, color = "black"),
    axis.title.x = element_text(size = 7, color = "black", margin = margin(t = 2)),
    axis.line.x = element_line(linewidth = 0.25, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.25, color = "black"),
    legend.position = "top",
    legend.justification = "left",
    legend.key.size = unit(3.2, "mm"),
    legend.text = element_text(size = 6.3),
    legend.margin = margin(0, 0, 1, 0),
    legend.box.margin = margin(0, 0, 0, -2),
    plot.margin = margin(2, 3, 2, 2)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ggsave(
  file.path(plot_dir, "hallmark_logFC_rank_top_terms_compact_full_labels.png"),
  p_full_labels,
  width = 3.55,
  height = 2.9,
  dpi = 600,
  bg = "white"
)
ggsave(
  file.path(plot_dir, "hallmark_logFC_rank_top_terms_compact_full_labels.pdf"),
  p_full_labels,
  width = 3.55,
  height = 2.9,
  bg = "white",
  useDingbats = FALSE
)

fwrite(plot_dt[order(NES), .(
  ID,
  term = clean_term,
  compact_term,
  NES,
  pvalue,
  FDR = p.adjust,
  setSize,
  direction,
  core_enrichment
)], file.path(out_dir, "hallmark_logFC_rank_top_terms_for_manuscript.csv"))

summary_lines <- c(
  "# Hallmark GSEA, log2FC Ranked",
  "",
  paste0("- Bulk object: `", input_rds, "`."),
  "- Ranking metric: edgeR `logFC` only. No p-value weighting is used in the rank.",
  "- Positive rank/NES means higher expression in anti-IL-1b relative to IgG.",
  paste0("- Ranked unique gene symbols: ", length(gene_list), "."),
  paste0("- Hallmark terms tested: ", nrow(res), "."),
  paste0("- Significant lower terms at FDR < 0.05: ", res[NES < 0 & p.adjust < 0.05, .N], "."),
  paste0("- Significant higher terms at FDR < 0.05: ", res[NES > 0 & p.adjust < 0.05, .N], "."),
  "",
  "Top lower terms:",
  paste0(
    "- ",
    paste(res[NES < 0][order(p.adjust, pvalue, NES)][1:8, clean_term], collapse = "; "),
    "."
  ),
  "",
  "Top higher terms:",
  paste0(
    "- ",
    paste(res[NES > 0][order(p.adjust, pvalue, -NES)][1:4, clean_term], collapse = "; "),
    "."
  ),
  "",
  "Main outputs:",
  "- `MSigDB_Hallmark_logFC_GSEA_results.csv`",
  "- `ranked_gene_list_symbol_logFC.rnk`",
  "- `plots/hallmark_logFC_rank_top_terms_ultracompact_manuscript.png`",
  "- `plots/hallmark_logFC_rank_top_terms_ultracompact_manuscript.pdf`",
  "- `plots/hallmark_logFC_rank_top_terms_compact_full_labels.png`",
  "- `plots/hallmark_logFC_rank_top_terms_compact_full_labels.pdf`"
)
writeLines(summary_lines, file.path(out_dir, "hallmark_logFC_rank_GSEA_summary.md"))

cat("Wrote logFC-ranked Hallmark GSEA outputs to ", out_dir, "\n", sep = "")
