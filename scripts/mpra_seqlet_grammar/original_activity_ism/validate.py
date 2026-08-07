# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Validate the MPRA activity head and ensemble on the held-out fold.
# Steps: 1) Load held-out sequences; 2) score single and ensemble models; 3) report correlations and ISM timing estimates.
# =======================================
"""Sanity check: reproduce the notebook's held-out (fold 10) performance on CPU,
and time inference so we can gauge ISM feasibility."""
import csv, time
import numpy as np
from koo_model import (BUNDLE, load_single, load_ensemble,
                       predict_single, predict_ensemble, pearson, spearman, device)

print("device:", device)
TSV = BUNDLE / "log2_alpha_Endo_train_data.tsv"
TEST_FOLD = "10"

rows = [r for r in csv.DictReader(open(TSV), delimiter="\t") if r["fold"] == TEST_FOLD]
cores = [r["seq"].strip().upper() for r in rows]
y = np.array([float(r["mean_value"]) for r in rows], np.float32)
print(f"{len(cores)} test sequences, {len(cores[0])} bp")

print("loading single model...")
t0 = time.time()
model = load_single()
print(f"  loaded in {time.time()-t0:.1f}s")

t0 = time.time()
p_single = predict_single(model, cores)
dt = time.time() - t0
print(f"single: Pearson {pearson(y,p_single):.3f}  Spearman {spearman(y,p_single):.3f}"
      f"   ({dt:.1f}s for {len(cores)} seq = {1000*dt/len(cores):.1f} ms/seq)")

print("loading ensemble heads...")
base, heads = load_ensemble(model)
t0 = time.time()
p_ens = predict_ensemble(base, heads, cores)
dt = time.time() - t0
print(f"{len(heads)}-fold ensemble: Pearson {pearson(y,p_ens):.3f}  Spearman {spearman(y,p_ens):.3f}"
      f"   ({dt:.1f}s = {1000*dt/len(cores):.1f} ms/seq)")

# ISM cost projection: 170 positions x 3 alt bases = 510 mutant forwards / sequence
print(f"\nISM cost: ~{510*1000*dt/len(cores)/1000:.1f}s per full-170bp sequence (ensemble, CPU)")
