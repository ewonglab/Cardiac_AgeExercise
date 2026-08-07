#!/usr/bin/env Rscript
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Plot cross-contrast anti-IL-1β, aging and exercise GSEA concordance summaries.
# Steps: 1) Load the combined GSEA table; 2) calculate family-wise q values and selection sets; 3) export pathway heatmaps and source tables.
# =======================================
# Source reformatted from anti-il1code/plot_edgeR_gsea_all_terms_q020_two_tests.R.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(cowplot)
})

project_dir <- publication_path("H3K27AC_PROJECT_ROOT")
in_file <- file.path(
  project_dir, "results/antiil1b_age_exercise_agent_explore",
  "edgeR_pseudobulk_hallmark_go_gsea_mimicry",
  "clusterProfiler_Hallmark_GO_GSEA_edgeR_pseudobulk_ranked_FC_long.csv"
)
out_dir <- file.path(
  project_dir, "results/antiil1b_age_exercise_agent_explore",
  "edgeR_GSEA_all_terms_q020_two_tests"
)
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

q_cutoff <- 0.20
go_rows_per_page <- 36L

clean_term <- function(x) {
  out <- gsub("^(HALLMARK|GOBP|GOCC|GOMF)_", "", x)
  out <- gsub("_", " ", out)
  tools::toTitleCase(tolower(out))
}

clean_hallmark <- function(x) {
  x <- clean_term(x)
  replacements <- c(
    "Tnfa" = "TNFalpha", "Nfkb" = "NFkB", "Il2" = "IL2", "Il6" = "IL6",
    "Jak" = "JAK", "Stat3" = "STAT3", "Pi3k" = "PI3K", "Mtorc1" = "mTORC1",
    "Mtor" = "mTOR", "G2m" = "G2M", "E2f" = "E2F", "Myc" = "MYC",
    "Kras" = "KRAS", "Uv" = "UV"
  )
  for (from in names(replacements)) x <- gsub(from, replacements[[from]], x, fixed = TRUE)
  x
}

wrap_label <- function(x, width = 45) {
  vapply(strwrap(x, width = width, simplify = FALSE), paste, character(1), collapse = "\n")
}

fmt_q <- function(q) {
  ifelse(
    is.na(q), "q=NA",
    ifelse(q < 0.001, "q<0.001",
           ifelse(q < 0.01, paste0("q=", formatC(q, format = "f", digits = 3)),
                  paste0("q=", formatC(q, format = "f", digits = 2))))
  )
}

theme_pub <- function(base_size = 6) {
  theme_classic(base_size = base_size, base_family = "ArialMT") +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1.0, margin = margin(b = 1)),
      plot.subtitle = element_text(size = base_size - 0.2, color = "grey25", margin = margin(b = 1)),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.2, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.3),
      legend.key.width = unit(3.2, "mm"),
      legend.key.height = unit(2.7, "mm"),
      panel.spacing = unit(1, "mm"),
      plot.margin = margin(2, 2, 2, 2)
    )
}

save_plot <- function(p, stub, width, height) {
  qpdf <- function(filename, width, height, ...) {
    quartz(type = "pdf", file = filename, width = width, height = height, family = "ArialMT")
  }
  ggsave(file.path(plot_dir, paste0(stub, ".pdf")), p, width = width, height = height,
         units = "in", device = qpdf, bg = "white")
  ggsave(file.path(plot_dir, paste0(stub, ".png")), p, width = width, height = height,
         units = "in", dpi = 600, bg = "white")
}

gsea <- fread(in_file)
gsea[, comparison := factor(
  rank_name,
  levels = c("anti_IL1b", "edgeR_age", "edgeR_exercise"),
  labels = c("anti-IL1b", "age", "exercise")
)]

# Hallmark q values are adjusted only across the Hallmark family.
gsea[collection_label == "Hallmark", q_family := p.adjust(pvalue, method = "BH"),
     by = rank_name]

# GO q values are adjusted across all GO BP/CC/MF terms together.
gsea[collection_label != "Hallmark", q_family := p.adjust(pvalue, method = "BH"),
     by = rank_name]

wide_q <- dcast(
  gsea,
  collection_label + term ~ rank_name,
  value.var = "q_family"
)
setnames(
  wide_q,
  c("anti_IL1b", "edgeR_age", "edgeR_exercise"),
  c("q_anti", "q_age", "q_exercise")
)
wide_q[, n_tests_q020 := rowSums(
  cbind(q_anti < q_cutoff, q_age < q_cutoff, q_exercise < q_cutoff),
  na.rm = TRUE
)]

selected_keys <- wide_q[n_tests_q020 >= 2, .(collection_label, term, n_tests_q020)]
selected <- merge(
  gsea,
  selected_keys,
  by = c("collection_label", "term"),
  all = FALSE
)
selected[, significance_border := fifelse(
  q_family < 0.05, "q < 0.05",
  fifelse(q_family < 0.10, "q < 0.10",
          fifelse(q_family < 0.20, "q < 0.20", NA_character_))
)]
selected[, tile_label := paste0(sprintf("%.2f", NES), "\n", fmt_q(q_family))]
selected[, NES_limited := pmax(pmin(NES, 2.6), -2.6)]

hallmark <- copy(selected[collection_label == "Hallmark"])
hallmark[, term_display := clean_hallmark(term)]

go <- copy(selected[collection_label != "Hallmark"])
go[, ontology := factor(
  collection_label,
  levels = c("GO BP", "GO CC", "GO MF"),
  labels = c("BP", "CC", "MF")
)]
go[, term_display := paste0(ontology, ": ", clean_term(term))]

cluster_order <- function(dt) {
  mat <- dcast(dt, term_display ~ comparison, value.var = "NES")
  row_names <- mat$term_display
  x <- as.matrix(mat[, -"term_display"])
  rownames(x) <- row_names
  if (nrow(x) <= 2) return(row_names)
  rownames(x)[hclust(dist(x), method = "complete")$order]
}

hallmark_order <- cluster_order(hallmark)
go_order <- cluster_order(go)

make_heatmap <- function(dt, term_order, title, subtitle, base_size = 6,
                         y_text_size = 5.4, label_size = 1.55) {
  x <- copy(dt)
  x[, term_display := factor(term_display, levels = rev(term_order))]
  x[, comparison := factor(comparison, levels = c("anti-IL1b", "age", "exercise"))]

  ggplot(x, aes(comparison, term_display)) +
    geom_tile(
      aes(fill = NES_limited),
      color = "grey88", linewidth = 0.16, width = 0.9, height = 0.9
    ) +
    geom_tile(
      data = x[!is.na(significance_border)],
      aes(color = significance_border),
      fill = NA, linewidth = 0.44, width = 0.9, height = 0.9
    ) +
    geom_text(
      aes(label = tile_label),
      family = "ArialMT", size = label_size, lineheight = 0.82
    ) +
    scale_fill_gradient2(
      low = "#315AA6", mid = "white", high = "#C93F42",
      midpoint = 0, limits = c(-2.6, 2.6),
      breaks = c(-2, -1, 0, 1, 2), name = "NES"
    ) +
    scale_color_manual(
      values = c("q < 0.05" = "black", "q < 0.10" = "grey35", "q < 0.20" = "grey65"),
      name = "FDR"
    ) +
    guides(
      fill = guide_colorbar(order = 1, barheight = unit(18, "mm")),
      color = guide_legend(order = 2, override.aes = list(fill = NA, linewidth = 0.55))
    ) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_pub(base_size) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      axis.text.y = element_text(size = y_text_size),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right"
    )
}

p_hallmark <- make_heatmap(
  hallmark,
  hallmark_order,
  "edgeR pseudobulk Hallmarks significant in at least two comparisons",
  "All Hallmarks with q < 0.20 in >=2 tests; no directional filtering; clustered by NES",
  base_size = 6,
  y_text_size = 5.3,
  label_size = 1.5
)
save_plot(
  p_hallmark,
  "Hallmark_all_q020_in_at_least_two_tests",
  4.75,
  max(3.0, 0.225 * length(hallmark_order) + 0.9)
)

hallmark_compact <- copy(hallmark)
hallmark_compact[, tile_label := paste0(
  sprintf("%.2f", NES),
  "\n",
  fifelse(
    q_family < 0.001, "q<.001",
    fifelse(
      q_family < 0.01,
      paste0("q=", sub("^0", "", formatC(q_family, format = "f", digits = 3))),
      paste0("q=", sub("^0", "", formatC(q_family, format = "f", digits = 2)))
    )
  )
)]

p_hallmark_compact <- make_heatmap(
  hallmark_compact,
  hallmark_order,
  NULL,
  NULL,
  base_size = 5.5,
  y_text_size = 5.0,
  label_size = 1.05
) +
  scale_x_discrete(expand = expansion(add = 0)) +
  scale_y_discrete(expand = expansion(add = 0)) +
  coord_fixed(ratio = 1, clip = "off") +
  theme(
    axis.text.x = element_text(
      angle = 42, hjust = 1, vjust = 1,
      size = 5.0, margin = margin(t = 0.5)
    ),
    axis.text.y = element_text(
      size = 5.0, margin = margin(r = 0.7)
    ),
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(0.5, "mm"),
    legend.spacing.y = unit(0.6, "mm"),
    legend.key.width = unit(2.5, "mm"),
    legend.key.height = unit(2.4, "mm"),
    plot.margin = margin(0.4, 0.4, 0.4, 0.4)
  ) +
  scale_color_manual(
    values = c("q < 0.05" = "black", "q < 0.10" = "grey35", "q < 0.20" = "grey65"),
    labels = c("q<.05", "q<.10", "q<.20"),
    name = "FDR"
  ) +
  guides(
    fill = guide_colorbar(
      order = 1,
      barheight = unit(13, "mm"),
      barwidth = unit(2.3, "mm"),
      title.position = "top"
    ),
    color = guide_legend(
      order = 2,
      override.aes = list(fill = NA, linewidth = 0.5),
      title.position = "top",
      keyheight = unit(2.4, "mm")
    )
  )
save_plot(
  p_hallmark_compact,
  "Hallmark_all_q020_in_at_least_two_tests_compact_square",
  2.75,
  4.05
)

# Main-figure candidate: apply a direction-independent q < 0.05 rule first,
# then display the selected terms in blocks according to their NES directions.
hallmark_q05_wide <- dcast(
  gsea[collection_label == "Hallmark"],
  term ~ rank_name,
  value.var = c("NES", "q_family")
)
hallmark_q05_wide[, n_tests_q005 := rowSums(cbind(
  q_family_anti_IL1b < 0.05,
  q_family_edgeR_age < 0.05,
  q_family_edgeR_exercise < 0.05
), na.rm = TRUE)]
hallmark_q05_wide <- hallmark_q05_wide[n_tests_q005 >= 2]
hallmark_q05_wide[, direction_group := fifelse(
  sign(NES_anti_IL1b) == sign(NES_edgeR_exercise) &
    sign(NES_anti_IL1b) == -sign(NES_edgeR_age),
  "Anti-IL1b and exercise concordant;\nage NES opposite",
  "Other directional patterns"
)]

hallmark_q05 <- merge(
  gsea[collection_label == "Hallmark"],
  hallmark_q05_wide[, .(term, n_tests_q005, direction_group)],
  by = "term",
  all = FALSE
)
hallmark_q05[, term_display := clean_hallmark(term)]
hallmark_q05[, significance_border := fifelse(
  q_family < 0.05, "q < 0.05",
  fifelse(q_family < 0.10, "q < 0.10",
          fifelse(q_family < 0.20, "q < 0.20", NA_character_))
)]
hallmark_q05[, tile_label := paste0(
  sprintf("%.2f", NES),
  "\n",
  fifelse(
    q_family < 0.001, "q<.001",
    fifelse(
      q_family < 0.01,
      paste0("q=", sub("^0", "", formatC(q_family, format = "f", digits = 3))),
      paste0("q=", sub("^0", "", formatC(q_family, format = "f", digits = 2)))
    )
  )
)]
hallmark_q05[, NES_limited := pmax(pmin(NES, 2.6), -2.6)]

expected_label <- "Anti-IL1b and exercise concordant;\nage NES opposite"
other_label <- "Other directional patterns"
hallmark_q05_expected <- hallmark_q05[direction_group == expected_label]
hallmark_q05_other <- hallmark_q05[direction_group == other_label]
expected_order <- cluster_order(hallmark_q05_expected)
other_order <- cluster_order(hallmark_q05_other)

make_compact_block <- function(dt, term_order, block_title, show_legend = FALSE) {
  x <- copy(dt)
  x[, term_display := factor(term_display, levels = rev(term_order))]
  x[, comparison := factor(
    rank_name,
    levels = c("anti_IL1b", "edgeR_age", "edgeR_exercise"),
    labels = c("anti-IL1b", "age", "exercise")
  )]

  ggplot(x, aes(comparison, term_display)) +
    geom_tile(
      aes(fill = NES_limited),
      color = "grey88", linewidth = 0.14, width = 0.91, height = 0.91
    ) +
    geom_tile(
      data = x[!is.na(significance_border)],
      aes(color = significance_border),
      fill = NA, linewidth = 0.42, width = 0.91, height = 0.91
    ) +
    geom_text(
      aes(label = tile_label),
      family = "ArialMT", size = 1.15, lineheight = 0.82
    ) +
    scale_fill_gradient2(
      low = "#315AA6", mid = "white", high = "#C93F42",
      midpoint = 0, limits = c(-2.6, 2.6),
      breaks = c(-2, -1, 0, 1, 2), name = "NES"
    ) +
    scale_color_manual(
      values = c(
        "q < 0.05" = "black",
        "q < 0.10" = "grey35",
        "q < 0.20" = "grey65"
      ),
      limits = c("q < 0.05", "q < 0.10", "q < 0.20"),
      breaks = c("q < 0.05", "q < 0.10", "q < 0.20"),
      labels = c("<0.05", "<0.10", "<0.20"),
      name = "FDR",
      drop = FALSE
    ) +
    scale_x_discrete(expand = expansion(add = 0)) +
    scale_y_discrete(expand = expansion(add = 0)) +
    coord_fixed(ratio = 1, clip = "off") +
    labs(title = block_title, x = NULL, y = NULL) +
    theme_pub(5.5) +
    theme(
      plot.title = element_text(
        size = 5.4, face = "bold", lineheight = 0.92,
        hjust = 0, margin = margin(b = 1.2)
      ),
      axis.text.x = element_text(
        angle = 40, hjust = 1, vjust = 1,
        size = 5.0, margin = margin(t = 0.4)
      ),
      axis.text.y = element_text(
        size = 5.0, margin = margin(r = 0.6)
      ),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      legend.position = if (show_legend) "right" else "none",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(0.5, "mm"),
      legend.spacing.y = unit(0.6, "mm"),
      legend.key.width = unit(2.5, "mm"),
      legend.key.height = unit(2.4, "mm"),
      plot.margin = margin(0.4, 0.4, 0.4, 0.4)
    ) +
    guides(
      fill = guide_colorbar(
        order = 1, barheight = unit(12, "mm"), barwidth = unit(2.3, "mm"),
        title.position = "top"
      ),
      color = guide_legend(
        order = 2, override.aes = list(fill = NA, linewidth = 0.5),
        title.position = "top", keyheight = unit(2.4, "mm")
      )
    )
}

p_q05_expected <- make_compact_block(
  hallmark_q05_expected,
  expected_order,
  "anti-IL1b = exercise;\nage NES opposite"
)
p_q05_other <- make_compact_block(
  hallmark_q05_other,
  other_order,
  "other directional patterns"
)
p_q05_legend <- make_compact_block(
  hallmark_q05_expected,
  expected_order,
  NULL,
  show_legend = TRUE
)
q05_legend <- get_legend(p_q05_legend)

p_hallmark_q05_blocks <- ggdraw() +
  draw_plot(p_q05_expected, x = 0.00, y = 0.00, width = 0.52, height = 1.00) +
  draw_plot(p_q05_other, x = 0.52, y = 0.31, width = 0.34, height = 0.69) +
  draw_grob(q05_legend, x = 0.86, y = 0.04, width = 0.14, height = 0.92)

save_plot(
  p_hallmark_q05_blocks,
  "Hallmark_q005_in_at_least_two_tests_direction_blocks",
  5.2,
  2.35
)

hallmark_q05_table <- dcast(
  hallmark_q05,
  term + term_display + direction_group + n_tests_q005 ~ rank_name,
  value.var = c("NES", "q_family")
)
hallmark_q05_table[, display_order := fifelse(
  direction_group == expected_label,
  match(term_display, expected_order),
  match(term_display, other_order)
)]
setorder(hallmark_q05_table, direction_group, display_order)
fwrite(
  hallmark_q05_table,
  file.path(out_dir, "Hallmark_q005_in_at_least_two_tests_direction_blocks.csv")
)

go_pages <- split(
  go_order,
  ceiling(seq_along(go_order) / go_rows_per_page)
)
go_pdf <- file.path(plot_dir, "GO_all_q020_in_at_least_two_tests_complete.pdf")
quartz(
  type = "pdf", file = go_pdf, width = 5.65, height = 8.2,
  family = "ArialMT"
)
for (i in seq_along(go_pages)) {
  page_terms <- go_pages[[i]]
  page_dt <- go[term_display %in% page_terms]
  p_page <- make_heatmap(
    page_dt,
    page_terms,
    "edgeR pseudobulk GO terms significant in at least two comparisons",
    paste0(
      "GO-only q < 0.20 in >=2 tests; no directional filtering; clustered by NES | ",
      "page ", i, " of ", length(go_pages)
    ),
    base_size = 5.6,
    y_text_size = 4.6,
    label_size = 1.35
  )
  print(p_page)
  ggsave(
    file.path(plot_dir, sprintf("GO_all_q020_in_at_least_two_tests_page_%02d.png", i)),
    p_page, width = 5.65, height = 8.2, units = "in", dpi = 450, bg = "white"
  )
}
dev.off()

hallmark_table <- dcast(
  hallmark,
  collection_label + term + term_display + n_tests_q020 ~ rank_name,
  value.var = c("NES", "q_family")
)
go_table <- dcast(
  go,
  collection_label + ontology + term + term_display + n_tests_q020 ~ rank_name,
  value.var = c("NES", "q_family")
)
hallmark_table[, cluster_order := match(term_display, hallmark_order)]
go_table[, cluster_order := match(term_display, go_order)]
setorder(hallmark_table, cluster_order)
setorder(go_table, cluster_order)

fwrite(
  hallmark_table,
  file.path(out_dir, "Hallmark_all_q020_in_at_least_two_tests.csv")
)
fwrite(
  go_table,
  file.path(out_dir, "GO_all_q020_in_at_least_two_tests.csv")
)

counts <- rbindlist(list(
  selected_keys[collection_label == "Hallmark", .(
    family = "Hallmark",
    n_terms = .N,
    n_q020_in_two = sum(n_tests_q020 == 2),
    n_q020_in_three = sum(n_tests_q020 == 3)
  )],
  selected_keys[collection_label != "Hallmark", .(
    family = "All GO",
    n_terms = .N,
    n_q020_in_two = sum(n_tests_q020 == 2),
    n_q020_in_three = sum(n_tests_q020 == 3)
  )],
  selected_keys[collection_label != "Hallmark", .(
    family = collection_label,
    n_terms = .N,
    n_q020_in_two = sum(n_tests_q020 == 2),
    n_q020_in_three = sum(n_tests_q020 == 3)
  ), by = collection_label][, collection_label := NULL]
), use.names = TRUE, fill = TRUE)
fwrite(counts, file.path(out_dir, "selected_term_counts.csv"))

summary_lines <- c(
  "# Uncurated edgeR pseudobulk GSEA heatmaps",
  "",
  "Selection rule: q < 0.20 in at least two of anti-IL1b, age, and exercise.",
  "No directional pattern was used for selection.",
  "Hallmark FDR was calculated across all Hallmark terms per comparison.",
  "GO FDR was calculated across all GO BP/CC/MF terms together per comparison.",
  "Rows were ordered by complete-linkage hierarchical clustering of the three NES values.",
  "",
  "## Counts",
  capture.output(print(counts))
)
writeLines(summary_lines, file.path(out_dir, "analysis_summary.md"))

message("Wrote uncurated GSEA heatmaps to: ", out_dir)
