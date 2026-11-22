# 10-week c-FOS CUT&Tag workflow

These scripts cover the 10-week c-FOS reversal analysis, known-motif enrichment, and selected publication figures. Set the `CUTTAG_10_WEEK_*`, `PEAK2PATHWAY_ROOT`, `PAPER_OUTPUT_ROOT` and `HOMER_BIN` variables in `config/.env`, export them, then run from the repository root.

Suggested order:

1. Submit `homer_standard_background_10wk.sh` to a PBS-compatible environment after staging the nf-core/CUT&Run output below `CUTTAG_10_WEEK_ROOT`.
2. Run `cuttag_reversal_10wk.R` to make consensus sets, counts and differential/reversal summaries.
3. Run the `plot_cfos_*`, `plot_chr10_peak_reversal.R` and `cuttag_tobias_pathway_overlap.R` scripts after their required upstream tables are staged below `CUTTAG_10_WEEK_ANALYSIS_ROOT` and `PEAK2PATHWAY_ROOT`.
