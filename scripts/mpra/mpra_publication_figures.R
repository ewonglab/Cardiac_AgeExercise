# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Generate the selected MPRA publication panels.
# Steps: 1) Load the final MPRA result table; 2) calculate group and matched summaries; 3) export individual and combined figures.
# =======================================
# Source reformatted from emily_plots/MPRA/mpra_publication_figures.R.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

# MPRA Publication Figures
# Key message: Dynamic binding sites (age up, age down, ex up, and age-up
# sites reverted by exercise) share enhanced sequence-context-dependent
# regulatory activity. This resides in unannotated sequence elements,
# not core TF binding motifs.

library(tidyverse)
library(ggplot2)
library(ggpubr)

setwd(publication_path("MPRA_FIGURE_ROOT"))
df <- read.csv("MPRA_Results_with_TSS_distance.csv")
df$Group <- factor(df$Group, levels = c("Stable", "Age Up", "Age Down", "Ex Up", "Ex Down"))

# -------------------------------------------------------
# Define readable subcategory labels
# -------------------------------------------------------
df$subcat_label <- case_when(
  df$category == "In All Conditions" ~ "Stable\n(all conditions)",
  df$category == "In Aged Sed" & df$Group == "Age Up"   ~ "In Aged Sed\n(Age Up)",
  df$category == "In Aged Sed" & df$Group == "Ex Down"  ~ "In Aged Sed\n(Ex Down)",
  df$category == "Lost YoungSed"     ~ "Lost Young Sed\n(Age Up)",
  df$category == "In Young Sedentary" ~ "In Young Sed\n(Age Down)",
  df$category == "Lost AgeSed" & df$Group == "Age Down" ~ "Lost Aged Sed\n(Age Down)",
  df$category == "Lost AgeSed" & df$Group == "Ex Up"    ~ "Lost Aged Sed\n(Ex Up)",
  df$category == "In Aged Ex"        ~ "In Aged Ex\n(Ex Up)",
  df$category == "Lost AgeEx"        ~ "Lost Aged Ex\n(Ex Down)"
)

df$subcat_label <- factor(df$subcat_label, levels = c(
  "Stable\n(all conditions)",
  "In Aged Sed\n(Age Up)", "Lost Young Sed\n(Age Up)",
  "In Young Sed\n(Age Down)", "Lost Aged Sed\n(Age Down)",
  "In Aged Ex\n(Ex Up)", "Lost Aged Sed\n(Ex Up)",
  "In Aged Sed\n(Ex Down)", "Lost Aged Ex\n(Ex Down)"
))

# Colour by biological behaviour
df$behaviour <- case_when(
  df$category == "In All Conditions" ~ "Stable",
  df$category == "Lost AgeEx"        ~ "Context-dependent\n(no sequence signal)",
  TRUE                               ~ "Sequence-dependent\n(dynamic)"
)
df$behaviour <- factor(df$behaviour, levels = c(
  "Stable", "Sequence-dependent\n(dynamic)", "Context-dependent\n(no sequence signal)"
))

behaviour_colors <- c(
  "Stable" = "#808080",
  "Sequence-dependent\n(dynamic)" = "#E64B35",
  "Context-dependent\n(no sequence signal)" = "#3C5488"
)

group_colors <- c(
  "Stable" = "#808080", "Age Up" = "#E64B35", "Age Down" = "#4DBBD5",
  "Ex Up" = "#00A087", "Ex Down" = "#3C5488"
)

# ================================================================
# PANEL A: Subcategory-level delta vs full shuffle
# Shows all subcategories — "In Aged Sed" appears twice (Age Up &
# Ex Down) with identical values, demonstrating these are the same
# sequences. "Lost AgeEx" is the clear outlier.
# ================================================================

# Significance vs Stable
stable_delta <- df$delta_vs_FullShuff[df$Group == "Stable"]
subcat_pvals <- df %>%
  filter(subcat_label != "Stable\n(all conditions)") %>%
  group_by(subcat_label) %>%
  summarise(
    mean_delta = mean(delta_vs_FullShuff),
    p = wilcox.test(delta_vs_FullShuff, stable_delta)$p.value,
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(p_label = case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  ))

cat("=== Subcategory p-values vs Stable ===\n")
print(subcat_pvals, width = Inf)

panelA <- ggplot(df, aes(x = subcat_label, y = delta_vs_FullShuff, fill = behaviour)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_boxplot(width = 0.6, outlier.size = 0.5, alpha = 0.85) +
  scale_fill_manual(values = behaviour_colors) +
  # Add n labels at bottom
  stat_summary(fun.data = function(x) data.frame(y = -3.2, label = paste0("n=", length(x))),
               geom = "text", size = 2.8, color = "grey30") +
  labs(
    x = "", y = "Delta log2(alpha)\n(Endogenous - Full Shuffle)",
    fill = ""
  ) +
  coord_cartesian(ylim = c(-3.5, 4)) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "top",
    legend.text = element_text(size = 9),
    plot.margin = margin(5, 10, 5, 15)
  )

ggsave("Fig_panelA_subcategory_delta.pdf", panelA, width = 10, height = 5.5)

# ================================================================
# PANEL B: Context vs Motif — paired group-level comparison
# Shows that full shuffle disrupts activity but motif scramble does not
# ================================================================

df_comp <- df %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    mean_full  = mean(delta_vs_FullShuff),
    se_full    = sd(delta_vs_FullShuff) / sqrt(n()),
    mean_motif = mean(delta_vs_MotifShuff),
    se_motif   = sd(delta_vs_MotifShuff) / sqrt(n()),
    .groups = "drop"
  )

df_comp_long <- bind_rows(
  df_comp %>% transmute(Group, comparison = "Entire sequence\nshuffled",
                         mean_delta = mean_full, se = se_full),
  df_comp %>% transmute(Group, comparison = "Motifs only\nscrambled",
                         mean_delta = mean_motif, se = se_motif)
) %>%
  mutate(comparison = factor(comparison, levels = c("Entire sequence\nshuffled", "Motifs only\nscrambled")))

panelB <- ggplot(df_comp_long, aes(x = Group, y = mean_delta, fill = comparison)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = mean_delta - se, ymax = mean_delta + se),
                position = position_dodge(0.7), width = 0.15, linewidth = 0.4) +
  scale_fill_manual(values = c("Entire sequence\nshuffled" = "#E64B35",
                                "Motifs only\nscrambled" = "#91D1C2")) +
  labs(
    x = "", y = "Mean delta log2(alpha)\n(Endogenous - Control)",
    fill = ""
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "top",
    legend.text = element_text(size = 9)
  )

ggsave("Fig_panelB_context_vs_motif.pdf", panelB, width = 7, height = 5)

# ================================================================
# PANEL C: Matched analysis summary — forest plot style
# Per-group matched results (GC + TSS distance matched)
# ================================================================

library(MatchIt)

# Run per-group matching and collect results
# Use "Age Up / Ex Down" (In Aged Sed only) instead of full Ex Down
matched_results <- list()

for (grp in c("Age Up", "Age Down", "Ex Up")) {
  df_pair <- df %>% filter(Group %in% c("Stable", grp))
  df_pair$is_dynamic <- ifelse(df_pair$Group == grp, 1, 0)

  set.seed(42)
  m_pair <- matchit(is_dynamic ~ GC_content + log10_dist_TSS, data = df_pair,
                    method = "nearest", ratio = 1, caliper = 0.2)
  df_m <- match.data(m_pair)

  s <- df_m$delta_vs_FullShuff[df_m$Group == "Stable"]
  d <- df_m$delta_vs_FullShuff[df_m$Group == grp]

  matched_results[[grp]] <- tibble(
    Group = grp,
    n_matched = length(d),
    mean_dynamic = mean(d),
    mean_stable  = mean(s),
    diff = mean(d) - mean(s),
    se_diff = sqrt(var(d)/length(d) + var(s)/length(s)),
    p = wilcox.test(d, s)$p.value
  )
}

# Age Up / Ex Down: match only "In Aged Sed" sequences (age-gained, exercise-reverted)
df_aged_sed_exdown <- df %>% filter(category == "In Aged Sed", Group == "Ex Down")
df_pair_aued <- bind_rows(df %>% filter(Group == "Stable"), df_aged_sed_exdown)
df_pair_aued$is_dynamic <- ifelse(df_pair_aued$Group == "Ex Down", 1, 0)

set.seed(42)
m_aued <- matchit(is_dynamic ~ GC_content + log10_dist_TSS, data = df_pair_aued,
                  method = "nearest", ratio = 1, caliper = 0.2)
df_m_aued <- match.data(m_aued)

s_aued <- df_m_aued$delta_vs_FullShuff[df_m_aued$Group == "Stable"]
d_aued <- df_m_aued$delta_vs_FullShuff[df_m_aued$Group == "Ex Down"]

matched_results[["Age Up / Ex Down"]] <- tibble(
  Group = "Age Up / Ex Down",
  n_matched = length(d_aued),
  mean_dynamic = mean(d_aued),
  mean_stable  = mean(s_aued),
  diff = mean(d_aued) - mean(s_aued),
  se_diff = sqrt(var(d_aued)/length(d_aued) + var(s_aued)/length(s_aued)),
  p = wilcox.test(d_aued, s_aued)$p.value
)

mr <- bind_rows(matched_results)
mr$Group <- factor(mr$Group, levels = rev(c("Age Up", "Age Down", "Ex Up", "Age Up / Ex Down")))
mr$sig <- ifelse(mr$p < 0.05, "p < 0.05", "ns")

cat("\n=== Matched analysis summary ===\n")
print(mr, width = Inf)

panelC <- ggplot(mr, aes(x = diff, y = Group, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = diff - 1.96 * se_diff, xmax = diff + 1.96 * se_diff),
                 height = 0.2, linewidth = 0.6) +
  geom_point(size = 3.5) +
  scale_color_manual(values = c("p < 0.05" = "#E64B35", "ns" = "#808080")) +
  geom_text(aes(label = sprintf("p = %.3f", p)), hjust = -0.15, vjust = -0.8, size = 3.2,
            color = "grey30") +
  geom_text(aes(label = sprintf("n = %d pairs", n_matched)), x = max(mr$diff + 1.96 * mr$se_diff) + 0.02,
            hjust = 0, size = 3, color = "grey50") +
  labs(
    x = "Difference in delta log2(alpha)\n(Dynamic - Matched Stable)",
    y = "",
    color = ""
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 11)
  )

ggsave("Fig_panelC_matched_forest.pdf", panelC, width = 7, height = 3.5)

# ================================================================
# PANEL D: "In Aged Sed" — same sequences, two contexts
# Directly shows the 246 shared sequences behave identically
# whether classified as Age Up or Ex Down
# ================================================================

df_aged_sed <- df %>%
  filter(category == "In Aged Sed") %>%
  select(Coords, Group, delta_vs_FullShuff, delta_vs_MotifShuff, Endo_residual)

# These are paired — same Coords in both groups
df_paired <- df_aged_sed %>%
  pivot_wider(id_cols = Coords, names_from = Group, values_from = delta_vs_FullShuff)

panelD <- ggplot(df_paired, aes(x = `Age Up`, y = `Ex Down`)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.5, size = 1.5, color = "#E64B35") +
  stat_cor(method = "pearson", size = 3.5, label.x.npc = 0.05, label.y.npc = 0.95) +
  labs(
    x = "Delta log2(alpha) when classified as Age Up",
    y = "Delta log2(alpha) when classified as Ex Down",
    title = "\"In Aged Sed\" sequences (n = 246)",
    subtitle = "Same sequences — age-gained sites reverted by exercise"
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "grey30"))

ggsave("Fig_panelD_shared_sequences.pdf", panelD, width = 5, height = 5)

# ================================================================
# COMBINED FIGURE: A + B top row, C + D bottom row
# ================================================================
combined <- ggarrange(
  panelA, panelB,
  panelC, panelD,
  labels = c("A", "B", "C", "D"),
  ncol = 2, nrow = 2,
  heights = c(1.1, 0.9),
  widths = c(1.2, 1)
)

ggsave("Fig_MPRA_combined.pdf", combined, width = 14, height = 11)

# ================================================================
# Print summary statistics for the text
# ================================================================
cat("\n\n=== STATISTICS FOR TEXT ===\n\n")

# Group-level delta vs full shuffle
cat("Delta vs Full Shuffle by Group:\n")
df %>%
  group_by(Group) %>%
  summarise(mean = mean(delta_vs_FullShuff), se = sd(delta_vs_FullShuff)/sqrt(n()),
            n = n(), .groups = "drop") %>%
  print(width = Inf)

# "In Aged Sed" in both contexts
cat("\n'In Aged Sed' sequences (n=246):\n")
cat(sprintf("  As Age Up:  mean delta = %.3f\n", mean(df$delta_vs_FullShuff[df$category == "In Aged Sed" & df$Group == "Age Up"])))
cat(sprintf("  As Ex Down: mean delta = %.3f\n", mean(df$delta_vs_FullShuff[df$category == "In Aged Sed" & df$Group == "Ex Down"])))
cat(sprintf("  Correlation: r = %.3f\n", cor(df_paired$`Age Up`, df_paired$`Ex Down`)))

# Lost AgeEx
cat("\n'Lost AgeEx' sequences (n=101):\n")
lost_ageex <- df %>% filter(category == "Lost AgeEx")
cat(sprintf("  Mean delta = %.3f (vs Stable mean = %.3f)\n",
            mean(lost_ageex$delta_vs_FullShuff), mean(df$delta_vs_FullShuff[df$Group == "Stable"])))
test <- wilcox.test(lost_ageex$delta_vs_FullShuff, df$delta_vs_FullShuff[df$Group == "Stable"])
cat(sprintf("  vs Stable: p = %.2e\n", test$p.value))

# Dynamic (excluding Lost AgeEx) vs Stable
cat("\nAll dynamic sites excluding Lost AgeEx:\n")
dynamic_no_lostageex <- df %>% filter(Group != "Stable", category != "Lost AgeEx")
test2 <- wilcox.test(dynamic_no_lostageex$delta_vs_FullShuff, stable_delta)
cat(sprintf("  n = %d, mean delta = %.3f vs Stable = %.3f, p = %.2e\n",
            nrow(dynamic_no_lostageex), mean(dynamic_no_lostageex$delta_vs_FullShuff),
            mean(stable_delta), test2$p.value))

cat("\nFigures saved:\n")
cat("  Fig_panelA_subcategory_delta.pdf\n")
cat("  Fig_panelB_context_vs_motif.pdf\n")
cat("  Fig_panelC_matched_forest.pdf\n")
cat("  Fig_panelD_shared_sequences.pdf\n")
cat("  Fig_MPRA_combined.pdf\n")
