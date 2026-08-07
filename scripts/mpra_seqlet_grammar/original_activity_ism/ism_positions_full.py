# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Run full in-silico saturation mutagenesis across all MPRA enhancers.
# Steps: 1) Load the ensemble and sequences; 2) score all single-base substitutions; 3) save per-position importance arrays.
# =======================================
"""Full in-silico saturation mutagenesis on ALL 1653 enhancers (not the 310 subset), saving the
per-position importance matrix (imp, 1653x170) so the motif analyses (FigF4/F5/F7) can be re-run at
full per-subgroup power. imp[p] = wild-type predicted activity - most-deactivating single-base
substitution at p. Runs on MPS (9-fold AlphaGenome ensemble)."""
import os, sys, time
import numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from koo_model import load_ensemble, predict_ensemble, device
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from config.paths import publication_path

MPRA = Path(publication_path("MPRA_DATA_ROOT"))
OUT = Path(publication_path("MPRA_SEQLET_TASK_ROOT", "outputs", "original_activity_ism")); BASES = "ACGT"
df = __import__("pandas").read_csv(MPRA / "corrected_data" / "MPRA_Results_with_TSS_distance.csv").drop_duplicates("Coords").reset_index(drop=True)
seqs = {}; nm = None
for l in open(MPRA / "deep_learning" / "results" / "sequences_cache.fa"):
    if l.startswith(">"): nm = l.strip()[1:]
    else: seqs[nm] = l.strip().upper()
sequences = [seqs[c] for c in df["Coords"]]; grp = df["Group"].to_numpy(); coords = df["Coords"].to_numpy()
L = len(sequences[0]); N = len(sequences)
base, heads = load_ensemble()
def predict(x): return np.asarray(predict_ensemble(base, heads, x))
def mutate(s, p, b): a = list(s); a[p] = b; return "".join(a)
print(f"device={device}; full ISM on N={N} enhancers", flush=True)

IMP = np.zeros((N, L), dtype=np.float32); t0 = time.time()
for n in range(N):
    s = sequences[n]; wt = float(predict([s])[0])
    muts, meta = [], []
    for p in range(L):
        for b in BASES:
            if b != s[p]: muts.append(mutate(s, p, b)); meta.append(p)
    preds = predict(muts); by = {}
    for p, pr in zip(meta, preds): by.setdefault(p, []).append(float(pr))
    IMP[n] = np.array([wt - min(by[p]) for p in range(L)])
    if (n+1) % 50 == 0:
        el = time.time()-t0; print(f"  {n+1}/{N} ({el:.0f}s, ETA {el/(n+1)*(N-n-1)/60:.0f}min)", flush=True)
        np.savez(OUT / "ism_imp_full.npz", imp=IMP[:n+1], groups=grp[:n+1], coords=coords[:n+1], n_done=n+1)
np.savez(OUT / "ism_imp_full.npz", imp=IMP, groups=grp, coords=coords, n_done=N)
print(f"Saved out/ism_imp_full.npz ({time.time()-t0:.0f}s)", flush=True)
