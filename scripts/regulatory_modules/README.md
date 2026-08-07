# Regulatory-module and AP-1 linking analyses

Configure and export `SIGNALCRE_ROOT`, `PEAK2PATHWAY_ROOT` and `PAPER_OUTPUT_ROOT` from `config/.env` before running these scripts from the `publication_code` directory.

- `module_gc_content.R` calculates GC percentage for existing NMF-module peak assignments. It expects `outputs_pooled/chromatin_modules.csv` below `SIGNALCRE_ROOT` and writes `outputs_tobias/module_gc_content.csv` there.
- `tre_vs_cre_footprinting.R` compares TRE (FOS/JUN) and CRE (CREB/ATF) footprint trajectories using `Tobias_df_Classified_Full.rds` staged below `PEAK2PATHWAY_ROOT`.
- `pathway_linking_figures.R` links AP-1 sites to ECM and IL-1 gene modules using the same TOBIAS object and `mm10_refGene_tss.rds`.

The NMF module-construction workflow and animal-level GC-redistribution inference are not present in the supplied sources; the first script begins from module assignments generated upstream.
