# AlphaGenome MPRA sequence-context workflow

This directory is the publication copy of the supplied AlphaGenome MPRA seqlet-grammar bundle. It contains 21 Python source files and two small input fixtures only. The original bundle remains unchanged in `RealCode/MPRA_seqlet_grammar_code_bundle/`; its source checksums are retained in `provenance/mpra_seqlet_grammar_source_SHA256SUMS.txt`.

## Required configuration and inputs

Copy `config/.env.example` to `config/.env`, populate the variables below, and export it from the `publication_code` directory before running a script.

| Variable | Required contents |
| --- | --- |
| `MPRA_DATA_ROOT` | `corrected_data/MPRA_Results_with_TSS_distance.csv` and `deep_learning/results/sequences_cache.fa` |
| `MPRA_KOO_PROJECT` | `bundle/` (frozen encoder checkpoint, MPRA heads and MPRA inputs) and `koorun/` (the original model helper, `out/embeddings.npy` and `out/jaspar_panel.meme`) |
| `MPRA_SEQLET_TASK_ROOT` | A writable analysis root. The workflow reads and writes its intermediate tables and figures in `outputs/` below this location. |

The original activity/ISM code also requires the locally supplied `alphagenome_encoder_ft` package. Install the packages listed in `environment/requirements_mpra_seqlet_grammar.txt` in an environment compatible with the frozen AlphaGenome checkpoint. The checkpoint, fold-specific heads, encoder package, embeddings, motif panel, raw MPRA data and generated intermediate tables are intentionally not bundled here.

## Model-training reference

The original AlphaGenome MPRA model-training notebook is retained externally as a read-only reference: [Google Colab notebook](https://colab.research.google.com/drive/1BQK3NK4TWkoi-9fDRWXXr9wqAbTiLzXW#scrollTo=6TDg6mlHVGHb). It is linked only; no notebook content was downloaded, copied, or modified.

## Run order

Run from the `publication_code` directory after exporting configuration:

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

`spacing_epistasis_scan.py` imports the bundled `spacing_grammar_sweep.py` helper from its own directory; no duplicate `work/` directory is required.

## Data-free checks

The supplied fixtures allow two checks without the model or raw data:

```bash
python scripts/mpra_seqlet_grammar/audit_followup_panel.py
python scripts/mpra_seqlet_grammar/figures/plot_top_interacting_seqlet_pairs.py \
  --output-prefix "$MPRA_SEQLET_TASK_ROOT/outputs/Fig_top_interacting_seqlet_pairs"
```

The first reproduces the 144 panel entries and 142 unique coordinates. The second regenerates the top-pair plot from the included pair table. Pre-rendered figure files are excluded from this repository because they are generated outputs.

## Analysis terminology

The 142 unique coordinates are the *mechanistic follow-up panel*, not 142 high-confidence loci: the 144 entries comprise 96 mechanism-reliable entries, 24 deliberate model-failure controls and 24 AP-1-responsive controls, with two overlapping coordinates. Attribution and single-seqlet perturbation use all 142 coordinates; epistasis and rescue use the 91 mechanism-reliable loci with complete high- and low-attribution seqlet sets.
