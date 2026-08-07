#!/usr/bin/env Rscript
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Assess anti-IL-1β RNA responsiveness of genes linked to dynamic versus stable AP-1 sites.
# Steps: 1) Load RNA, age/exercise and AP-1 linkage inputs; 2) calculate gene-level contrasts and adjusted models; 3) export summaries and supporting plots.
# =======================================
# Source reformatted from anti-il1code/test_ap1_dynamic_vs_stable_gene_responsiveness.R.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

set.seed(20260730)

project_dir <- publication_path("H3K27AC_PROJECT_ROOT")
rna_rds <- publication_path("ANTI_IL1B_RNA_EDGER_OBJECT")
age_ex_rds <- publication_path("AGE_EX_DEG_OBJECT")
ap1_pairs_file <- file.path(
  project_dir,
  "results",
  "rnaseq_tobias_ap1_100kb_age_effect_overlap",
  "TSS_100kb_AP1_overlap_pairs_with_age_effect.csv"
)
nearest_assignments_file <- file.path(
  project_dir,
  "results",
  "rnaseq_tobias_ap1_antiil1b_refinements",
  "nearest_expressed_gene_site_assignments.csv"
)

out_dir <- file.path(project_dir, "results", "rnaseq_ap1_dynamic_vs_stable_gene_responsiveness")
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

major_cell_types <- c(
  "Cardiomyocyte", "Fibroblast", "Coronary_EC", "Endocardial_EC",
  "Lymphatic_EC", "Pericytes", "HSPC", "Smooth_Muscle", "AllMacrophage",
  "B_cell", "T_cell"
)

class_defs <- data.table(
  class_key = c(
    "age_changed",
    "exercise_changed",
    "age_gain",
    "exercise_loss",
    "age_gain_exercise_loss",
    "age_loss",
    "exercise_gain",
    "age_loss_exercise_gain",
    "any_dynamic"
  ),
  class_label = c(
    "Age changed AP-1",
    "Exercise changed AP-1",
    "Age-gain AP-1",
    "Exercise-loss AP-1",
    "Age-gain/exercise-loss AP-1",
    "Age-loss AP-1",
    "Exercise-gain AP-1",
    "Age-loss/exercise-gain AP-1",
    "Any dynamic AP-1"
  ),
  expected_sign = c(NA, NA, -1, -1, -1, 1, 1, 1, NA)
)
class_levels <- class_defs$class_label

distance_thresholds <- c(10000, 25000, 50000, 100000)
distance_labels <- c("10 kb", "25 kb", "50 kb", "100 kb")

clean_gene <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "NA"] <- NA_character_
  x
}

short_ct <- function(x) gsub("_", " ", x)

fmt_p <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(x < 0.001, formatC(x, format = "e", digits = 1), as.character(signif(x, 2)))
  )
}

safe_wilcox <- function(x, y, alternative = "greater") {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 5 || length(y) < 5 || length(unique(c(x, y))) < 2) {
    return(NA_real_)
  }
  suppressWarnings(wilcox.test(x, y, alternative = alternative, exact = FALSE)$p.value)
}

safe_wilcox_one <- function(x, alternative = "greater") {
  x <- x[is.finite(x)]
  if (length(x) < 5 || length(unique(x)) < 2) {
    return(NA_real_)
  }
  suppressWarnings(wilcox.test(x, mu = 0, alternative = alternative, exact = FALSE)$p.value)
}

safe_fisher_greater <- function(dynamic_success, dynamic_fail, stable_success, stable_fail) {
  if (any(!is.finite(c(dynamic_success, dynamic_fail, stable_success, stable_fail)))) {
    return(list(or = NA_real_, pvalue = NA_real_))
  }
  if ((dynamic_success + dynamic_fail) < 5 || (stable_success + stable_fail) < 5) {
    return(list(or = NA_real_, pvalue = NA_real_))
  }
  ft <- fisher.test(
    matrix(
      c(dynamic_success, dynamic_fail, stable_success, stable_fail),
      nrow = 2,
      byrow = TRUE
    ),
    alternative = "greater"
  )
  list(or = unname(ft$estimate), pvalue = ft$p.value)
}

safe_lm_group <- function(dt, y_col) {
  d <- copy(dt[is.finite(get(y_col)) & is.finite(logCPM)])
  if (nrow(d) < 20 || uniqueN(d$group) < 2 || uniqueN(d[[y_col]]) < 2) {
    return(list(beta = NA_real_, pvalue = NA_real_))
  }
  d[, group_dynamic := as.integer(group == "Dynamic")]
  d[, z_logCPM := as.numeric(scale(logCPM))]
  fit <- tryCatch(lm(get(y_col) ~ group_dynamic + z_logCPM, data = d), error = function(e) NULL)
  if (is.null(fit)) {
    return(list(beta = NA_real_, pvalue = NA_real_))
  }
  co <- summary(fit)$coefficients
  if (!"group_dynamic" %in% rownames(co)) {
    return(list(beta = NA_real_, pvalue = NA_real_))
  }
  list(beta = unname(co["group_dynamic", "Estimate"]), pvalue = unname(co["group_dynamic", "Pr(>|t|)"]))
}

safe_glm_group <- function(dt) {
  d <- copy(dt[is.finite(logCPM)])
  if (nrow(d) < 20 || uniqueN(d$group) < 2 || uniqueN(d$nominal_P0.05) < 2) {
    return(list(beta = NA_real_, pvalue = NA_real_))
  }
  d[, group_dynamic := as.integer(group == "Dynamic")]
  d[, z_logCPM := as.numeric(scale(logCPM))]
  fit <- tryCatch(
    glm(nominal_P0.05 ~ group_dynamic + z_logCPM, data = d, family = binomial()),
    error = function(e) NULL,
    warning = function(w) suppressWarnings(glm(nominal_P0.05 ~ group_dynamic + z_logCPM, data = d, family = binomial()))
  )
  if (is.null(fit)) {
    return(list(beta = NA_real_, pvalue = NA_real_))
  }
  co <- summary(fit)$coefficients
  if (!"group_dynamic" %in% rownames(co)) {
    return(list(beta = NA_real_, pvalue = NA_real_))
  }
  list(beta = unname(co["group_dynamic", "Estimate"]), pvalue = unname(co["group_dynamic", "Pr(>|z|)"]))
}

safe_lm_count_competition <- function(dt, y_col, dynamic_count_col) {
  d <- copy(dt[is.finite(get(y_col)) & is.finite(logCPM)])
  if (nrow(d) < 30 || uniqueN(d[[dynamic_count_col]]) < 2 || uniqueN(d$stable_site_count) < 2) {
    return(list(
      beta_dynamic = NA_real_, p_dynamic = NA_real_,
      beta_stable = NA_real_, p_stable = NA_real_,
      beta_dynamic_minus_stable = NA_real_, p_dynamic_minus_stable = NA_real_
    ))
  }
  d[, z_dynamic := as.numeric(scale(log1p(get(dynamic_count_col))))]
  d[, z_stable := as.numeric(scale(log1p(stable_site_count)))]
  d[, z_logCPM := as.numeric(scale(logCPM))]
  fit <- tryCatch(lm(get(y_col) ~ z_dynamic + z_stable + z_logCPM, data = d), error = function(e) NULL)
  if (is.null(fit)) {
    return(list(
      beta_dynamic = NA_real_, p_dynamic = NA_real_,
      beta_stable = NA_real_, p_stable = NA_real_,
      beta_dynamic_minus_stable = NA_real_, p_dynamic_minus_stable = NA_real_
    ))
  }
  co <- summary(fit)$coefficients
  vc <- vcov(fit)
  needed <- c("z_dynamic", "z_stable")
  if (!all(needed %in% rownames(co)) || !all(needed %in% rownames(vc))) {
    return(list(
      beta_dynamic = NA_real_, p_dynamic = NA_real_,
      beta_stable = NA_real_, p_stable = NA_real_,
      beta_dynamic_minus_stable = NA_real_, p_dynamic_minus_stable = NA_real_
    ))
  }
  diff_est <- unname(co["z_dynamic", "Estimate"] - co["z_stable", "Estimate"])
  diff_se <- sqrt(vc["z_dynamic", "z_dynamic"] + vc["z_stable", "z_stable"] - 2 * vc["z_dynamic", "z_stable"])
  diff_t <- diff_est / diff_se
  diff_p <- 2 * pt(abs(diff_t), df = fit$df.residual, lower.tail = FALSE)
  list(
    beta_dynamic = unname(co["z_dynamic", "Estimate"]),
    p_dynamic = unname(co["z_dynamic", "Pr(>|t|)"]),
    beta_stable = unname(co["z_stable", "Estimate"]),
    p_stable = unname(co["z_stable", "Pr(>|t|)"]),
    beta_dynamic_minus_stable = diff_est,
    p_dynamic_minus_stable = diff_p
  )
}

theme_pub <- function(base_size = 6) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = 7),
      plot.subtitle = element_text(size = 5.4, color = "grey30"),
      axis.title = element_text(size = 6),
      axis.text = element_text(size = 5.4, color = "black"),
      axis.line = element_line(linewidth = 0.18, color = "black"),
      axis.ticks = element_line(linewidth = 0.16, color = "black"),
      legend.title = element_text(size = 5.2),
      legend.text = element_text(size = 5.2),
      legend.key.height = unit(2.4, "mm"),
      legend.key.width = unit(3.0, "mm"),
      panel.spacing = unit(1.2, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 5.5),
      plot.margin = margin(2, 2, 2, 2)
    )
}

message("Reading anti-IL1b RNA")
obj <- readRDS(rna_rds)
anti <- as.data.table(obj$de)[, .(
  gene = clean_gene(gene_symbol),
  gene_id = as.character(gene_id),
  logFC = as.numeric(logFC),
  logCPM = as.numeric(logCPM),
  F = as.numeric(F),
  PValue = as.numeric(PValue),
  FDR = as.numeric(FDR)
)]
anti <- anti[!is.na(gene) & is.finite(logFC)]
anti[, abs_logFC := abs(logFC)]
setorder(anti, gene, PValue, -abs_logFC)
anti <- anti[!duplicated(gene)]
anti[, `:=`(
  abs_logFC = abs(logFC),
  neg_log10_p = -log10(pmax(PValue, 1e-300)),
  nominal_P0.05 = is.finite(PValue) & PValue < 0.05
)]

message("Reading cell-type expressed/tested universes")
age_ex <- as.data.table(readRDS(age_ex_rds))[, .(
  gene = clean_gene(gene),
  cell_type = as.character(CellType)
)]
age_ex <- unique(age_ex[!is.na(gene) & cell_type %in% major_cell_types & gene %in% anti$gene])
expressed_universe <- lapply(split(age_ex$gene, age_ex$cell_type), unique)

message("Reading AP-1 site-gene data")
raw_pairs <- fread(ap1_pairs_file)
raw_pairs[, `:=`(
  gene = clean_gene(gene_symbol),
  cell_type = as.character(cell_type),
  bound_in_young_sedentary = as.logical(bound_in_young_sedentary),
  bound_in_aged_sedentary = as.logical(bound_in_aged_sedentary),
  bound_in_aged_exercised = as.logical(bound_in_aged_exercised),
  distance_to_tss_bp = as.numeric(distance_to_tss_bp)
)]
raw_pairs <- raw_pairs[!is.na(gene) & cell_type %in% names(expressed_universe) & gene %in% anti$gene]
raw_pairs[, abs_distance_to_tss_bp := abs(distance_to_tss_bp)]

nearest_pairs <- fread(nearest_assignments_file)
nearest_pairs[, `:=`(
  gene = clean_gene(gene),
  cell_type = as.character(cell_type),
  abs_distance_to_tss_bp = as.numeric(abs_distance_to_tss_bp),
  bound_in_young_sedentary = as.logical(bound_in_young_sedentary),
  bound_in_aged_sedentary = as.logical(bound_in_aged_sedentary),
  bound_in_aged_exercised = as.logical(bound_in_aged_exercised)
)]
nearest_pairs <- nearest_pairs[!is.na(gene) & cell_type %in% names(expressed_universe) & gene %in% anti$gene]

add_classes <- function(dt) {
  dt[, `:=`(
    ap1_any = TRUE,
    ap1_stable_bound = bound_in_young_sedentary == TRUE &
      bound_in_aged_sedentary == TRUE &
      bound_in_aged_exercised == TRUE,
    ap1_age_gain = bound_in_aged_sedentary == TRUE & bound_in_young_sedentary != TRUE,
    ap1_age_loss = bound_in_young_sedentary == TRUE & bound_in_aged_sedentary != TRUE,
    ap1_exercise_gain = bound_in_aged_exercised == TRUE & bound_in_aged_sedentary != TRUE,
    ap1_exercise_loss = bound_in_aged_sedentary == TRUE & bound_in_aged_exercised != TRUE
  )]
  dt[, `:=`(
    ap1_age_changed = ap1_age_gain | ap1_age_loss,
    ap1_exercise_changed = ap1_exercise_gain | ap1_exercise_loss,
    ap1_age_gain_exercise_loss = ap1_age_gain & ap1_exercise_loss,
    ap1_age_loss_exercise_gain = ap1_age_loss & ap1_exercise_gain
  )]
  dt[, ap1_any_dynamic := ap1_any & !ap1_stable_bound]
  dt
}

raw_pairs <- add_classes(raw_pairs)
nearest_pairs <- add_classes(nearest_pairs)

assignment_inputs <- list(
  "Nearest expressed site" = nearest_pairs,
  "All TSS +/-100kb pairs" = raw_pairs
)

make_gene_level <- function(pairs_dt, assignment_label, distance_bp) {
  rows <- list()
  keys <- c(
    age_changed = "ap1_age_changed",
    exercise_changed = "ap1_exercise_changed",
    age_gain = "ap1_age_gain",
    exercise_loss = "ap1_exercise_loss",
    age_gain_exercise_loss = "ap1_age_gain_exercise_loss",
    age_loss = "ap1_age_loss",
    exercise_gain = "ap1_exercise_gain",
    age_loss_exercise_gain = "ap1_age_loss_exercise_gain",
    any_dynamic = "ap1_any_dynamic",
    stable = "ap1_stable_bound",
    any = "ap1_any"
  )
  for (ct in names(expressed_universe)) {
    universe <- expressed_universe[[ct]]
    base <- anti[gene %in% universe, .(gene, logFC, logCPM, PValue, FDR, abs_logFC, neg_log10_p, nominal_P0.05)]
    sub <- pairs_dt[cell_type == ct & abs_distance_to_tss_bp <= distance_bp]
    sub <- sub[gene %in% universe]
    for (key in names(keys)) {
      col <- keys[[key]]
      gm <- sub[get(col) == TRUE, .(
        site_count = uniqueN(site_key),
        min_abs_distance_to_tss_bp = min(abs_distance_to_tss_bp, na.rm = TRUE)
      ), by = gene]
      setnames(gm, c("site_count", "min_abs_distance_to_tss_bp"), c(paste0(key, "_site_count"), paste0(key, "_min_distance")))
      base <- merge(base, gm, by = "gene", all.x = TRUE, sort = FALSE)
    }
    count_cols <- grep("_site_count$", names(base), value = TRUE)
    distance_cols <- grep("_min_distance$", names(base), value = TRUE)
    for (col in count_cols) {
      base[is.na(get(col)), (col) := 0L]
    }
    for (col in distance_cols) {
      base[!is.finite(get(col)), (col) := NA_real_]
    }
    base[, `:=`(
      cell_type = ct,
      assignment = assignment_label,
      distance_threshold_bp = distance_bp,
      distance_label = paste0(distance_bp / 1000, " kb")
    )]
    rows[[ct]] <- base
  }
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

message("Building gene-level AP-1 dynamic/stable summaries")
gene_level <- rbindlist(lapply(names(assignment_inputs), function(mode) {
  rbindlist(lapply(distance_thresholds, function(distance_bp) {
    make_gene_level(assignment_inputs[[mode]], mode, distance_bp)
  }), use.names = TRUE, fill = TRUE)
}), use.names = TRUE, fill = TRUE)
fwrite(gene_level, file.path(out_dir, "gene_level_dynamic_stable_ap1_metrics.csv"))

message("Running dynamic-vs-stable contrast tests")
contrast_rows <- list()
for (mode in unique(gene_level$assignment)) {
  for (distance_bp in distance_thresholds) {
    for (ct in unique(gene_level$cell_type)) {
      base <- gene_level[assignment == mode & distance_threshold_bp == distance_bp & cell_type == ct]
      if (nrow(base) == 0) next
      stable_strict <- base[stable_site_count > 0 & any_dynamic_site_count == 0]
      stable_relaxed <- base[stable_site_count > 0]
      for (i in seq_len(nrow(class_defs))) {
        key <- class_defs$class_key[i]
        label <- class_defs$class_label[i]
        expected_sign <- class_defs$expected_sign[i]
        count_col <- paste0(key, "_site_count")
        if (!count_col %in% names(base)) next

        dyn <- base[get(count_col) > 0]
        stable_control <- stable_strict[get(count_col) == 0]
        control_label <- "Stable only"
        if (nrow(stable_control) < 25) {
          stable_control <- stable_relaxed[get(count_col) == 0]
          control_label <- "Stable-bound excluding class"
        }
        if (nrow(dyn) < 10 || nrow(stable_control) < 10) next

        test_dt <- rbind(
          dyn[, .(gene, logFC, logCPM, abs_logFC, neg_log10_p, nominal_P0.05, group = "Dynamic")],
          stable_control[, .(gene, logFC, logCPM, abs_logFC, neg_log10_p, nominal_P0.05, group = "Stable")]
        )

        ft <- safe_fisher_greater(
          sum(dyn$nominal_P0.05, na.rm = TRUE),
          sum(!dyn$nominal_P0.05, na.rm = TRUE),
          sum(stable_control$nominal_P0.05, na.rm = TRUE),
          sum(!stable_control$nominal_P0.05, na.rm = TRUE)
        )
        lm_abs <- safe_lm_group(test_dt, "abs_logFC")
        lm_rank <- safe_lm_group(test_dt, "neg_log10_p")
        glm_nominal <- safe_glm_group(test_dt)

        dyn_expected <- NA_real_
        stable_expected <- NA_real_
        p_expected_gt_stable <- NA_real_
        p_expected_gt_zero <- NA_real_
        frac_expected_dynamic <- NA_real_
        frac_expected_stable <- NA_real_
        lm_expected_beta <- NA_real_
        lm_expected_p <- NA_real_
        if (is.finite(expected_sign)) {
          test_dt[, expected_signed_logFC := expected_sign * logFC]
          dyn_expected <- median(test_dt[group == "Dynamic", expected_signed_logFC], na.rm = TRUE)
          stable_expected <- median(test_dt[group == "Stable", expected_signed_logFC], na.rm = TRUE)
          p_expected_gt_stable <- safe_wilcox(
            test_dt[group == "Dynamic", expected_signed_logFC],
            test_dt[group == "Stable", expected_signed_logFC],
            alternative = "greater"
          )
          p_expected_gt_zero <- safe_wilcox_one(
            test_dt[group == "Dynamic", expected_signed_logFC],
            alternative = "greater"
          )
          frac_expected_dynamic <- mean(test_dt[group == "Dynamic", expected_signed_logFC] > 0, na.rm = TRUE)
          frac_expected_stable <- mean(test_dt[group == "Stable", expected_signed_logFC] > 0, na.rm = TRUE)
          lm_expected <- safe_lm_group(test_dt, "expected_signed_logFC")
          lm_expected_beta <- lm_expected$beta
          lm_expected_p <- lm_expected$pvalue
        }

        contrast_rows[[paste(mode, distance_bp, ct, key, sep = "\r")]] <- data.table(
          assignment = mode,
          distance_threshold_bp = distance_bp,
          distance_label = paste0(distance_bp / 1000, " kb"),
          cell_type = ct,
          class_key = key,
          class_label = label,
          expected_sign = expected_sign,
          stable_control = control_label,
          dynamic_genes = nrow(dyn),
          stable_genes = nrow(stable_control),
          dynamic_sites = sum(dyn[[count_col]], na.rm = TRUE),
          stable_sites = sum(stable_control$stable_site_count, na.rm = TRUE),
          median_abs_logFC_dynamic = median(dyn$abs_logFC, na.rm = TRUE),
          median_abs_logFC_stable = median(stable_control$abs_logFC, na.rm = TRUE),
          delta_median_abs_logFC = median(dyn$abs_logFC, na.rm = TRUE) - median(stable_control$abs_logFC, na.rm = TRUE),
          abs_logFC_wilcox_gt_stable_p = safe_wilcox(dyn$abs_logFC, stable_control$abs_logFC, alternative = "greater"),
          lm_abs_logFC_beta_dynamic_vs_stable = lm_abs$beta,
          lm_abs_logFC_p_dynamic_vs_stable = lm_abs$pvalue,
          median_neg_log10p_dynamic = median(dyn$neg_log10_p, na.rm = TRUE),
          median_neg_log10p_stable = median(stable_control$neg_log10_p, na.rm = TRUE),
          delta_median_neg_log10p = median(dyn$neg_log10_p, na.rm = TRUE) - median(stable_control$neg_log10_p, na.rm = TRUE),
          neg_log10p_wilcox_gt_stable_p = safe_wilcox(dyn$neg_log10_p, stable_control$neg_log10_p, alternative = "greater"),
          lm_neg_log10p_beta_dynamic_vs_stable = lm_rank$beta,
          lm_neg_log10p_p_dynamic_vs_stable = lm_rank$pvalue,
          nominal_dynamic = sum(dyn$nominal_P0.05, na.rm = TRUE),
          nominal_stable = sum(stable_control$nominal_P0.05, na.rm = TRUE),
          nominal_fraction_dynamic = mean(dyn$nominal_P0.05, na.rm = TRUE),
          nominal_fraction_stable = mean(stable_control$nominal_P0.05, na.rm = TRUE),
          nominal_fisher_OR = ft$or,
          nominal_fisher_gt_stable_p = ft$pvalue,
          glm_nominal_beta_dynamic_vs_stable = glm_nominal$beta,
          glm_nominal_p_dynamic_vs_stable = glm_nominal$pvalue,
          median_expected_signed_dynamic = dyn_expected,
          median_expected_signed_stable = stable_expected,
          delta_median_expected_signed = dyn_expected - stable_expected,
          expected_signed_wilcox_gt_stable_p = p_expected_gt_stable,
          expected_signed_dynamic_gt0_p = p_expected_gt_zero,
          fraction_expected_dynamic = frac_expected_dynamic,
          fraction_expected_stable = frac_expected_stable,
          lm_expected_signed_beta_dynamic_vs_stable = lm_expected_beta,
          lm_expected_signed_p_dynamic_vs_stable = lm_expected_p
        )
      }
    }
  }
}
contrast_tests <- rbindlist(contrast_rows, use.names = TRUE, fill = TRUE)
if (nrow(contrast_tests) > 0) {
  contrast_tests[, `:=`(
    abs_logFC_wilcox_FDR = p.adjust(abs_logFC_wilcox_gt_stable_p, method = "BH"),
    lm_abs_logFC_FDR = p.adjust(lm_abs_logFC_p_dynamic_vs_stable, method = "BH"),
    neg_log10p_wilcox_FDR = p.adjust(neg_log10p_wilcox_gt_stable_p, method = "BH"),
    lm_neg_log10p_FDR = p.adjust(lm_neg_log10p_p_dynamic_vs_stable, method = "BH"),
    nominal_fisher_FDR = p.adjust(nominal_fisher_gt_stable_p, method = "BH"),
    glm_nominal_FDR = p.adjust(glm_nominal_p_dynamic_vs_stable, method = "BH"),
    expected_signed_wilcox_FDR = p.adjust(expected_signed_wilcox_gt_stable_p, method = "BH"),
    expected_signed_dynamic_gt0_FDR = p.adjust(expected_signed_dynamic_gt0_p, method = "BH"),
    lm_expected_signed_FDR = p.adjust(lm_expected_signed_p_dynamic_vs_stable, method = "BH")
  )]
}
fwrite(contrast_tests, file.path(out_dir, "dynamic_vs_stable_contrast_tests.csv"))

message("Running count-competition models")
model_rows <- list()
for (mode in unique(gene_level$assignment)) {
  for (distance_bp in distance_thresholds) {
    for (ct in unique(gene_level$cell_type)) {
      base <- gene_level[assignment == mode & distance_threshold_bp == distance_bp & cell_type == ct]
      if (nrow(base) == 0) next
      for (i in seq_len(nrow(class_defs))) {
        key <- class_defs$class_key[i]
        label <- class_defs$class_label[i]
        expected_sign <- class_defs$expected_sign[i]
        count_col <- paste0(key, "_site_count")
        if (!count_col %in% names(base)) next

        abs_model <- safe_lm_count_competition(base, "abs_logFC", count_col)
        rank_model <- safe_lm_count_competition(base, "neg_log10_p", count_col)

        exp_model <- list(
          beta_dynamic = NA_real_, p_dynamic = NA_real_,
          beta_stable = NA_real_, p_stable = NA_real_,
          beta_dynamic_minus_stable = NA_real_, p_dynamic_minus_stable = NA_real_
        )
        if (is.finite(expected_sign)) {
          base[, expected_signed_logFC_tmp := expected_sign * logFC]
          exp_model <- safe_lm_count_competition(base, "expected_signed_logFC_tmp", count_col)
          base[, expected_signed_logFC_tmp := NULL]
        }

        model_rows[[paste(mode, distance_bp, ct, key, sep = "\r")]] <- data.table(
          assignment = mode,
          distance_threshold_bp = distance_bp,
          distance_label = paste0(distance_bp / 1000, " kb"),
          cell_type = ct,
          class_key = key,
          class_label = label,
          expected_sign = expected_sign,
          tested_genes = nrow(base),
          dynamic_linked_genes = sum(base[[count_col]] > 0, na.rm = TRUE),
          stable_linked_genes = sum(base$stable_site_count > 0, na.rm = TRUE),
          abs_beta_dynamic_count = abs_model$beta_dynamic,
          abs_p_dynamic_count = abs_model$p_dynamic,
          abs_beta_stable_count = abs_model$beta_stable,
          abs_p_stable_count = abs_model$p_stable,
          abs_beta_dynamic_minus_stable = abs_model$beta_dynamic_minus_stable,
          abs_p_dynamic_minus_stable = abs_model$p_dynamic_minus_stable,
          rank_beta_dynamic_count = rank_model$beta_dynamic,
          rank_p_dynamic_count = rank_model$p_dynamic,
          rank_beta_stable_count = rank_model$beta_stable,
          rank_p_stable_count = rank_model$p_stable,
          rank_beta_dynamic_minus_stable = rank_model$beta_dynamic_minus_stable,
          rank_p_dynamic_minus_stable = rank_model$p_dynamic_minus_stable,
          expected_beta_dynamic_count = exp_model$beta_dynamic,
          expected_p_dynamic_count = exp_model$p_dynamic,
          expected_beta_stable_count = exp_model$beta_stable,
          expected_p_stable_count = exp_model$p_stable,
          expected_beta_dynamic_minus_stable = exp_model$beta_dynamic_minus_stable,
          expected_p_dynamic_minus_stable = exp_model$p_dynamic_minus_stable
        )
      }
    }
  }
}
count_models <- rbindlist(model_rows, use.names = TRUE, fill = TRUE)
if (nrow(count_models) > 0) {
  count_models[, `:=`(
    abs_FDR_dynamic_minus_stable = p.adjust(abs_p_dynamic_minus_stable, method = "BH"),
    rank_FDR_dynamic_minus_stable = p.adjust(rank_p_dynamic_minus_stable, method = "BH"),
    expected_FDR_dynamic_minus_stable = p.adjust(expected_p_dynamic_minus_stable, method = "BH")
  )]
}
fwrite(count_models, file.path(out_dir, "dynamic_vs_stable_count_competition_models.csv"))

message("Preparing plot-ready gene categories")
plot_categories <- rbindlist(lapply(c("AllMacrophage", "Fibroblast", "Cardiomyocyte"), function(ct) {
  d <- gene_level[
    assignment == "Nearest expressed site" &
      distance_threshold_bp == 100000 &
      cell_type == ct
  ]
  if (nrow(d) == 0) return(NULL)
  d[, category := fifelse(
    age_gain_exercise_loss_site_count > 0,
    "Age-gain/exercise-loss",
    fifelse(
      exercise_loss_site_count > 0,
      "Exercise-loss",
      fifelse(
        age_gain_site_count > 0,
        "Age-gain",
        fifelse(
          any_dynamic_site_count > 0,
          "Other dynamic",
          fifelse(stable_site_count > 0, "Stable only", "No AP-1")
        )
      )
    )
  )]
  d[, .(gene, cell_type, category, logFC, abs_logFC, neg_log10_p, nominal_P0.05)]
}), use.names = TRUE, fill = TRUE)
plot_categories <- plot_categories[
  category %in% c("Stable only", "Other dynamic", "Age-gain", "Exercise-loss", "Age-gain/exercise-loss")
]
plot_categories[, category := factor(
  category,
  levels = c("Stable only", "Other dynamic", "Age-gain", "Exercise-loss", "Age-gain/exercise-loss")
)]
plot_categories[, cell_label := factor(short_ct(cell_type), levels = c("AllMacrophage", "Fibroblast", "Cardiomyocyte"))]
fwrite(plot_categories, file.path(out_dir, "plot_gene_categories_nearest_100kb_key_celltypes.csv"))

message("Making plots")
heat_dt <- contrast_tests[
  assignment == "Nearest expressed site" &
    distance_threshold_bp == 100000 &
    class_key %in% c(
      "age_changed", "exercise_changed", "age_gain", "exercise_loss",
      "age_gain_exercise_loss", "age_loss", "exercise_gain",
      "age_loss_exercise_gain", "any_dynamic"
    )
]
heat_dt[, `:=`(
  cell_label = factor(short_ct(cell_type), levels = rev(short_ct(unique(cell_type)))),
  class_label = factor(class_label, levels = rev(class_levels)),
  abs_sig = abs_logFC_wilcox_FDR < 0.10,
  rank_sig = neg_log10p_wilcox_FDR < 0.10,
  expected_sig = expected_signed_wilcox_FDR < 0.10
)]

p_abs_heat <- ggplot(heat_dt, aes(x = class_label, y = cell_label, fill = delta_median_abs_logFC)) +
  geom_tile(color = "white", linewidth = 0.12) +
  geom_point(data = heat_dt[abs_logFC_wilcox_gt_stable_p < 0.05], shape = 16, size = 0.6, color = "black") +
  geom_point(data = heat_dt[abs_sig == TRUE], shape = 21, size = 1.3, stroke = 0.35, fill = NA, color = "black") +
  scale_fill_gradient2(
    low = "#3E6FA3",
    mid = "white",
    high = "#C85050",
    midpoint = 0,
    limits = c(-0.05, 0.05),
    oob = squish,
    name = "delta median\n|log2FC|"
  ) +
  labs(
    title = "Dynamic AP-1 genes versus stable-bound AP-1 genes",
    subtitle = "Nearest expressed gene assignment, TSS +/-100 kb; dots P < 0.05, rings FDR < 0.10",
    x = NULL,
    y = NULL
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5.0))
ggsave(file.path(plot_dir, "dynamic_vs_stable_abs_logFC_heatmap_nearest_100kb.pdf"), p_abs_heat, width = 6.6, height = 3.8, device = "pdf")
ggsave(file.path(plot_dir, "dynamic_vs_stable_abs_logFC_heatmap_nearest_100kb.png"), p_abs_heat, width = 6.6, height = 3.8, dpi = 450)

expected_heat <- heat_dt[is.finite(expected_sign)]
p_expected_heat <- ggplot(expected_heat, aes(x = class_label, y = cell_label, fill = delta_median_expected_signed)) +
  geom_tile(color = "white", linewidth = 0.12) +
  geom_point(data = expected_heat[expected_signed_wilcox_gt_stable_p < 0.05], shape = 16, size = 0.6, color = "black") +
  geom_point(data = expected_heat[expected_sig == TRUE], shape = 21, size = 1.3, stroke = 0.35, fill = NA, color = "black") +
  scale_fill_gradient2(
    low = "#3E6FA3",
    mid = "white",
    high = "#C85050",
    midpoint = 0,
    limits = c(-0.08, 0.08),
    oob = squish,
    name = "delta median\nexpected log2FC"
  ) +
  labs(
    title = "Expected-direction anti-IL1b response at dynamic AP-1 genes",
    subtitle = "Positive means stronger expected anti-IL1b movement than stable-only AP-1 genes",
    x = NULL,
    y = NULL
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5.0))
ggsave(file.path(plot_dir, "dynamic_vs_stable_expected_direction_heatmap_nearest_100kb.pdf"), p_expected_heat, width = 6.0, height = 3.8, device = "pdf")
ggsave(file.path(plot_dir, "dynamic_vs_stable_expected_direction_heatmap_nearest_100kb.png"), p_expected_heat, width = 6.0, height = 3.8, dpi = 450)

p_box_abs <- ggplot(plot_categories, aes(x = category, y = abs_logFC)) +
  geom_boxplot(width = 0.58, outlier.size = 0.12, linewidth = 0.18, fill = "grey92", color = "grey25") +
  stat_summary(fun = median, geom = "point", size = 0.75, color = "#8E2E2E") +
  facet_wrap(~cell_label, nrow = 1, scales = "free_y") +
  labs(
    title = "Anti-IL1b response magnitude by AP-1 site class",
    subtitle = "Nearest expressed AP-1-linked genes within TSS +/-100 kb",
    x = NULL,
    y = "|anti-IL1b log2FC|"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 5.0))
ggsave(file.path(plot_dir, "key_celltypes_dynamic_vs_stable_abs_logFC_boxplot.pdf"), p_box_abs, width = 6.8, height = 3.1, device = "pdf")
ggsave(file.path(plot_dir, "key_celltypes_dynamic_vs_stable_abs_logFC_boxplot.png"), p_box_abs, width = 6.8, height = 3.1, dpi = 450)

down_box <- plot_categories[category %in% c("Stable only", "Age-gain", "Exercise-loss", "Age-gain/exercise-loss")]
down_box[, expected_down_score := -logFC]
p_box_down <- ggplot(down_box, aes(x = category, y = expected_down_score)) +
  geom_hline(yintercept = 0, linewidth = 0.16, color = "grey55") +
  geom_boxplot(width = 0.58, outlier.size = 0.12, linewidth = 0.18, fill = "grey92", color = "grey25") +
  stat_summary(fun = median, geom = "point", size = 0.75, color = "#8E2E2E") +
  facet_wrap(~cell_label, nrow = 1, scales = "free_y") +
  labs(
    title = "Expected-lower AP-1 dynamic genes",
    subtitle = "Higher values mean anti-IL1b moves genes lower; compared with stable-only AP-1 genes",
    x = NULL,
    y = "-anti-IL1b log2FC"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 5.0))
ggsave(file.path(plot_dir, "key_celltypes_expected_lower_dynamic_vs_stable_boxplot.pdf"), p_box_down, width = 6.8, height = 3.1, device = "pdf")
ggsave(file.path(plot_dir, "key_celltypes_expected_lower_dynamic_vs_stable_boxplot.png"), p_box_down, width = 6.8, height = 3.1, dpi = 450)

model_plot_dt <- count_models[
  assignment == "Nearest expressed site" &
    distance_threshold_bp == 100000 &
    class_key %in% c("age_changed", "exercise_changed", "age_gain", "exercise_loss", "age_gain_exercise_loss", "any_dynamic")
]
model_plot_dt[, `:=`(
  cell_label = factor(short_ct(cell_type), levels = rev(short_ct(unique(cell_type)))),
  class_label = factor(class_label, levels = rev(class_levels)),
  sig = abs_FDR_dynamic_minus_stable < 0.10
)]
p_model <- ggplot(model_plot_dt, aes(x = abs_beta_dynamic_minus_stable, y = class_label)) +
  geom_vline(xintercept = 0, linewidth = 0.18, color = "grey55") +
  geom_point(aes(size = -log10(pmax(abs_p_dynamic_minus_stable, 1e-300)), fill = sig), shape = 21, stroke = 0.25) +
  facet_wrap(~cell_label, ncol = 4) +
  scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "#C85050"), name = "FDR < 0.10") +
  scale_size_continuous(range = c(0.6, 2.4), name = "-log10 P") +
  labs(
    title = "Dynamic site-count effect beyond stable site-count",
    subtitle = "Linear model: |anti-IL1b log2FC| ~ dynamic count + stable count + logCPM",
    x = "beta(dynamic count) - beta(stable count)",
    y = NULL
  ) +
  theme_pub()
ggsave(file.path(plot_dir, "dynamic_vs_stable_count_competition_abs_logFC_model.pdf"), p_model, width = 6.8, height = 4.4, device = "pdf")
ggsave(file.path(plot_dir, "dynamic_vs_stable_count_competition_abs_logFC_model.png"), p_model, width = 6.8, height = 4.4, dpi = 450)

message("Writing summary")
best_abs <- contrast_tests[
  assignment == "Nearest expressed site" &
    distance_threshold_bp == 100000 &
    delta_median_abs_logFC > 0
][order(abs_logFC_wilcox_gt_stable_p)][1:min(.N, 12)]

best_expected <- contrast_tests[
  assignment == "Nearest expressed site" &
    distance_threshold_bp == 100000 &
    is.finite(expected_sign) &
    delta_median_expected_signed > 0
][order(expected_signed_wilcox_gt_stable_p)][1:min(.N, 12)]

best_models <- count_models[
  assignment == "Nearest expressed site" &
    distance_threshold_bp == 100000 &
    abs_beta_dynamic_minus_stable > 0
][order(abs_p_dynamic_minus_stable)][1:min(.N, 12)]

summary_lines <- c(
  "# AP-1 dynamic versus stable gene responsiveness",
  "",
  "Question: are genes linked to age/exercise-dynamic AP-1 sites more anti-IL1b-responsive than genes linked only to stable-bound AP-1 sites?",
  "",
  "## Tests run",
  "- Dynamic-vs-stable Wilcoxon tests for |anti-IL1b log2FC| and -log10(P).",
  "- Fisher/logistic tests for enrichment of nominal anti-IL1b-responsive genes (P < 0.05).",
  "- Expected-direction tests for directional classes: age-gain/exercise-loss should move lower with anti-IL1b; age-loss/exercise-gain should move higher.",
  "- Count-competition models asking whether dynamic AP-1 site count predicts response beyond stable AP-1 site count and logCPM.",
  "- Repeated all tests for nearest-expressed-site and all TSS +/-100 kb assignments, at 10, 25, 50 and 100 kb.",
  "",
  "## Best absolute-response dynamic-vs-stable results at nearest expressed TSS +/-100 kb",
  paste(capture.output(print(best_abs[, .(
    cell_type, class_label, dynamic_genes, stable_genes,
    delta_median_abs_logFC, abs_logFC_wilcox_gt_stable_p, abs_logFC_wilcox_FDR,
    nominal_fraction_dynamic, nominal_fraction_stable, nominal_fisher_OR, nominal_fisher_gt_stable_p
  )], nrows = 12)), collapse = "\n"),
  "",
  "## Best expected-direction dynamic-vs-stable results at nearest expressed TSS +/-100 kb",
  paste(capture.output(print(best_expected[, .(
    cell_type, class_label, dynamic_genes, stable_genes,
    delta_median_expected_signed, expected_signed_wilcox_gt_stable_p, expected_signed_wilcox_FDR,
    fraction_expected_dynamic, fraction_expected_stable
  )], nrows = 12)), collapse = "\n"),
  "",
  "## Best dynamic-count versus stable-count models at nearest expressed TSS +/-100 kb",
  paste(capture.output(print(best_models[, .(
    cell_type, class_label, dynamic_linked_genes, stable_linked_genes,
    abs_beta_dynamic_minus_stable, abs_p_dynamic_minus_stable, abs_FDR_dynamic_minus_stable
  )], nrows = 12)), collapse = "\n"),
  "",
  "## Main output files",
  "- dynamic_vs_stable_contrast_tests.csv",
  "- dynamic_vs_stable_count_competition_models.csv",
  "- gene_level_dynamic_stable_ap1_metrics.csv",
  "- plot_gene_categories_nearest_100kb_key_celltypes.csv",
  "- plots/dynamic_vs_stable_abs_logFC_heatmap_nearest_100kb.pdf",
  "- plots/dynamic_vs_stable_expected_direction_heatmap_nearest_100kb.pdf",
  "- plots/key_celltypes_dynamic_vs_stable_abs_logFC_boxplot.pdf",
  "- plots/key_celltypes_expected_lower_dynamic_vs_stable_boxplot.pdf",
  "- plots/dynamic_vs_stable_count_competition_abs_logFC_model.pdf"
)
writeLines(summary_lines, file.path(out_dir, "dynamic_vs_stable_gene_responsiveness_summary.md"))

cat(out_dir, "\n", sep = "")
