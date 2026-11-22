# AlphaGenome MPRA sequence-context workflow

This directory contains the AlphaGenome MPRA sequence-context workflow: 21 Python source files and two small input fixtures.

## Required configuration and inputs

Copy `config/.env.example` to `config/.env`, populate the variables below, and export it from the repository root before running a script.

| Variable | Required contents |
| --- | --- |
| `MPRA_DATA_ROOT` | `corrected_data/MPRA_Results_with_TSS_distance.csv` and `deep_learning/results/sequences_cache.fa` |
| `MPRA_KOO_PROJECT` | `bundle/` (frozen encoder checkpoint, MPRA heads and MPRA inputs) and `koorun/` (the original model helper, `out/embeddings.npy` and `out/jaspar_panel.meme`) |
| `MPRA_SEQLET_TASK_ROOT` | A writable analysis root. The workflow reads and writes its intermediate tables and figures in `outputs/` below this location. |

The original activity/ISM code also requires the `alphagenome_encoder_ft` package. Install the packages listed in `environment/requirements_mpra_seqlet_grammar.txt` in an environment compatible with the frozen AlphaGenome checkpoint.

## Model training

The AlphaGenome MPRA model-training notebook is available as a [Google Colab notebook](https://colab.research.google.com/drive/1BQK3NK4TWkoi-9fDRWXXr9wqAbTiLzXW#scrollTo=6TDg6mlHVGHb).

## Run order

Run from the repository root after exporting configuration:

```bash
set -a; source config/.env; set +a

python scripts/mpra_seqlet_grammar/original_activity_ism/validate.py
python scripts/mpra_seqlet_grammar/original_activity_ism/ism_counterfactual.py
python scripts/mpra_seqlet_grammar/original_activity_ism/fig_f2cA_fisher.py
python scripts/mpra_seqlet_grammar/original_activity_ism/ism_positions_full.py
python scripts/mpra_seqlet_grammar/original_activity_ism/motif_top20_countmatch.py

python scripts/mpra_seqlet_grammar/contrast_seqlets/current_mpra_contrast_benchmark.py
python scripts/mpra_seqlet_grammar/contrast_seqlets/train_nested_contrast_probes.py
python scripts/mpra_seqlet_grammar/contrast_seqlets/mpra_contrast_followup.py
python scripts/mpra_seqlet_grammar/contrast_seqlets/context_probe_attribution.py
python scripts/mpra_seqlet_grammar/contrast_seqlets/context_module_disruption.py
python scripts/mpra_seqlet_grammar/contrast_seqlets/context_module_epistasis_rescue.py
python scripts/mpra_seqlet_grammar/contrast_seqlets/context_module_consolidate.py
python scripts/mpra_seqlet_grammar/contrast_seqlets/analyze_context_module_pairs.py
python scripts/mpra_seqlet_grammar/contrast_seqlets/redo_uncorrected_motif_analysis.py
python scripts/mpra_seqlet_grammar/contrast_seqlets/spacing_epistasis_scan.py

python scripts/mpra_seqlet_grammar/figures/make_context_manuscript_figure.py
python scripts/mpra_seqlet_grammar/figures/make_manuscript_evidence_figures.py
```

`spacing_epistasis_scan.py` imports the `spacing_grammar_sweep.py` helper from its own directory.


