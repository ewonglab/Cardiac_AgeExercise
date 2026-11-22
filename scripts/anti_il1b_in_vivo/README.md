# In-vivo anti-IL-1β workflow

These scripts cover the anti-IL-1β bulk-RNA and H3K27ac follow-up analyses.

## Configuration

Set these variables in `config/.env` before running the workflow:

- `H3K27AC_PROJECT_ROOT`: root containing the declared `results/` hierarchy.
- `ANTI_IL1B_RNA_EDGER_OBJECT`: anti-IL-1β bulk-RNA edgeR RDS object.
- `AGE_EX_DEG_OBJECT`: age/exercise differential-expression RDS object used by the supporting AP-1 gene-response analysis.

## Execution order

1. `overlap_h3k27ac_and_antiil1b_rna_promoters.R` creates promoter and nearby H3K27ac annotations from anti-IL-1β RNA/H3K27ac differential-result tables.
2. `plot_promoter_h3k27ac_rna_quintiles.R` consumes the promoter annotation and produces the supplementary RNA-quintile figure plus source tables.
3. `run_hallmark_gsea_logfc_rank.R` writes the `MSigDB_Hallmark_logFC_GSEA_results.csv` table consumed by Panel B of `make_final_antiil1b_manuscript_panels.R`. The assembler maps the standard `p.adjust` field to its local `FDR` display field without changing values; `test_allcells_top250_alternative_gene_set_statistics.R` generates the directional top-250 gene-set results used in the manuscript.
4. `plot_edgeR_gsea_all_terms_q020_two_tests.R` draws the cross-contrast GSEA summaries from the combined GSEA input table.
5. `plot_pooled_ap1_directional_excluding_mac_subsets_hspc.R` calculates the dynamic-versus-stable AP-1 H3K27ac table. `make_final_antiil1b_manuscript_panels.R` reads that table and the other final tables to assemble the manuscript panels.
6. `test_ap1_dynamic_vs_stable_gene_responsiveness.R` is a supporting AP-1 gene-response analysis; it is not required for the supplementary promoter plot.
