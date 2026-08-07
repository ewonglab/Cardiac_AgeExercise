# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Run balanced full ISM and minimal-edit counterfactual analysis.
# Steps: 1) Load the activity ensemble and sequences; 2) score substitutions and edits; 3) annotate AP-1 reliance and export results.
# =======================================
"""FULL single-base ISM + minimal-edit counterfactual, dynamic vs stable.

For a balanced sample of enhancers (Koo 9-fold ensemble, MPS):
  1. WT predicted activity; a composition-preserving shuffle -> compositional floor.
     arrangement_activity = WT - shuffle = the activity that depends on ARRANGEMENT
     (not bulk composition). This is what "grammar" would live in.
  2. FULL ISM: every position x 3 alt bases. Per-position deactivation importance =
     WT - min_alt(pred).
  3. Minimal-edit counterfactual: greedily apply the most-deactivating substitutions
     (ranked by ISM), cumulatively, and predict at each step. Measure how many edits it
     takes to remove 50% / 90% of arrangement_activity. FEW edits -> concentrated grammar;
     MANY diffuse edits -> distributed compositional code.
  4. Annotate the top-important positions: AP-1 core / homopolymer / other.
Compare dynamic vs stable on all of the above.
"""
import os, sys, json, time, re
import numpy as np
from pathlib import Path
from scipy import stats
sys.path.insert(0, str(Path(__file__).resolve().parent))
from koo_model import load_ensemble, predict_ensemble, device
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from config.paths import publication_path

MPRA = Path(publication_path("MPRA_DATA_ROOT"))
DATA_FILE = MPRA / "corrected_data" / "MPRA_Results_with_TSS_distance.csv"
SEQ_CACHE = MPRA / "deep_learning" / "results" / "sequences_cache.fa"
OUT = Path(publication_path("MPRA_SEQLET_TASK_ROOT", "outputs", "original_activity_ism")); (OUT / "pub").mkdir(parents=True, exist_ok=True)

SEED = 42
STAB_N = int(os.environ.get("STAB_N", "150"))
SUB_N = int(os.environ.get("SUB_N", "40"))
N_SHUF = int(os.environ.get("N_SHUF", "5"))
KGRID = [0, 1, 2, 3, 5, 8, 12, 17, 25, 35, 50, 70]
BASES = "ACGT"

import pandas as pd
df = pd.read_csv(DATA_FILE).drop_duplicates("Coords").reset_index(drop=True)
seqs = {}; name = None
for line in open(SEQ_CACHE):
    if line.startswith(">"): name = line.strip()[1:]
    else: seqs[name] = line.strip().upper()
sequences = [seqs[c] for c in df["Coords"]]
grp = df["Group"].to_numpy()
L = len(sequences[0])

rng = np.random.default_rng(SEED)
sample = []
for g, n in [("Stable", STAB_N), ("Age Up", SUB_N), ("Age Down", SUB_N), ("Ex Up", SUB_N), ("Ex Down", SUB_N)]:
    idx = np.where(grp == g)[0]
    pick = rng.choice(idx, min(n, len(idx)), replace=False)
    sample.extend((int(i), g) for i in pick)
print(f"sample: {len(sample)} seqs ({STAB_N} stable + subgroups); device={device}", flush=True)

base, heads = load_ensemble()
def predict(x): return np.asarray(predict_ensemble(base, heads, x))

AP1_RE = re.compile("TGA[CG]TCA")
def ap1_positions(s):
    pos = set()
    for mm in AP1_RE.finditer(s):
        pos.update(range(mm.start(), mm.end()))
    rc = s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    for mm in AP1_RE.finditer(rc):
        pos.update(range(L - mm.end(), L - mm.start()))   # map revcomp hit back to fwd coords
    return pos
def homopolymer_positions(s):
    pos = set()
    for mm in re.finditer(r"(A{4,}|T{4,}|G{4,}|C{4,})", s):
        pos.update(range(mm.start(), mm.end()))
    return pos


def mutate(s, edits):  # edits: list of (pos, base)
    a = list(s)
    for p, b in edits: a[p] = b
    return "".join(a)


rows = []
t0 = time.time()
for n, (i, g) in enumerate(sample):
    s = sequences[i]
    wt = float(predict([s])[0])
    # compositional floor: average over N_SHUF mononucleotide shuffles
    shuf_seqs = []
    for k in range(N_SHUF):
        r = np.random.default_rng(SEED * 100 + n * 10 + k)
        shuf_seqs.append("".join(r.permutation(list(s))))
    shuf = float(predict(shuf_seqs).mean())
    arr = wt - shuf                                  # arrangement-dependent activity

    # FULL ISM: 170 x 3 alts
    muts, meta = [], []
    for p in range(L):
        for b in BASES:
            if b != s[p]:
                muts.append(mutate(s, [(p, b)])); meta.append((p, b))
    preds = predict(muts)
    imp = np.zeros(L); best_alt = [None] * L
    by_pos = {}
    for (p, b), pr in zip(meta, preds):
        by_pos.setdefault(p, []).append((float(pr), b))
    for p in range(L):
        mn, mb = min(by_pos[p])                       # most-deactivating alt
        imp[p] = wt - mn                              # >0 means mutating drops activity
        best_alt[p] = mb

    # greedy cumulative knockout by importance
    order = np.argsort(-imp)
    edits = [(int(p), best_alt[p]) for p in order]
    kmax = min(max(KGRID), L)
    cum_seqs = [mutate(s, edits[:k]) for k in KGRID if k <= kmax]
    cum_pred = predict(cum_seqs)
    curve = {k: float(cum_pred[j]) for j, k in enumerate([k for k in KGRID if k <= kmax])}

    def edits_to_frac(frac):                          # edits to remove frac of arrangement_activity
        if arr <= 1e-6: return np.nan
        target = wt - frac * arr
        ks = sorted(curve)
        for k in ks:
            if curve[k] <= target: return k
        return np.nan

    ap1 = ap1_positions(s); homo = homopolymer_positions(s)
    top = order[:5]
    ann = {"ap1": int(sum(p in ap1 for p in top)), "homo": int(sum(p in homo for p in top)),
           "other": int(sum((p not in ap1 and p not in homo) for p in top))}
    # gini of importance (positive part)
    iv = np.clip(imp, 0, None); iv.sort(); c = np.cumsum(iv)
    gini = float((L + 1 - 2 * np.sum(c) / c[-1]) / L) if c[-1] > 0 else 0.0

    rows.append(dict(group=g, wt=wt, shuf=shuf, arrangement=arr,
                     edits_50=edits_to_frac(0.5), edits_90=edits_to_frac(0.9),
                     gini_imp=gini, top_ann=ann,
                     frac_arr_top5=float(np.clip(imp, 0, None)[top].sum() / arr) if arr > 1e-6 else np.nan))
    if (n + 1) % 25 == 0:
        print(f"  {n+1}/{len(sample)} ({time.time()-t0:.0f}s)", flush=True)

# aggregate dynamic vs stable
isdyn = np.array([r["group"] != "Stable" for r in rows])
def arr_of(key): return np.array([r[key] for r in rows], float)
def nmean(x):
    x = x[np.isfinite(x)]
    return float(x.mean()) if x.size else float("nan")
print("\n=== ISM counterfactual: dynamic vs stable ===", flush=True)
# fraction of sequences whose native arrangement is activating (WT > shuffle)
arr = arr_of("arrangement")
fpos = {"stable": float((arr[~isdyn] > 0).mean()), "dynamic": float((arr[isdyn] > 0).mean())}
print(f"  frac_arr_positive   stable={fpos['stable']:.2f}  dynamic={fpos['dynamic']:.2f}  "
      f"(activating native arrangement)", flush=True)
agg = {"frac_arr_positive": fpos}
for key in ["wt", "shuf", "arrangement", "edits_50", "edits_90", "gini_imp", "frac_arr_top5"]:
    v = arr_of(key); m = np.isfinite(v)
    ps = stats.mannwhitneyu(v[m & isdyn], v[m & ~isdyn]).pvalue if (m & isdyn).sum() and (m & ~isdyn).sum() else np.nan
    agg[key] = {"stable": nmean(v[~isdyn]), "dynamic": nmean(v[isdyn]), "mwu_p": float(ps),
                "n_stable": int((m & ~isdyn).sum()), "n_dynamic": int((m & isdyn).sum())}
    print(f"  {key:14s} stable={agg[key]['stable']:+.3f}(n={agg[key]['n_stable']})  "
          f"dynamic={agg[key]['dynamic']:+.3f}(n={agg[key]['n_dynamic']})  p={ps:.3f}", flush=True)
# annotation of top positions
for lab in ["ap1", "homo", "other"]:
    s_ = np.mean([r["top_ann"][lab] for r in rows if r["group"] == "Stable"])
    d_ = np.mean([r["top_ann"][lab] for r in rows if r["group"] != "Stable"])
    print(f"  top5 {lab:5s}/seq  stable={s_:.2f}  dynamic={d_:.2f}")

json.dump({"n": len(rows), "params": dict(STAB_N=STAB_N, SUB_N=SUB_N, N_SHUF=N_SHUF, KGRID=KGRID),
           "aggregate": agg,
           "top_annotation": {lab: {"stable": float(np.mean([r["top_ann"][lab] for r in rows if r["group"] == "Stable"])),
                                    "dynamic": float(np.mean([r["top_ann"][lab] for r in rows if r["group"] != "Stable"]))}
                              for lab in ["ap1", "homo", "other"]},
           "rows": rows}, open(OUT / "ism_counterfactual.json", "w"), indent=2, default=float)

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
                     "font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 300})
fig, ax = plt.subplots(1, 3, figsize=(9, 3))
for k, key in enumerate(["arrangement", "edits_50", "gini_imp"]):
    sv = arr_of(key)[~isdyn]; dv = arr_of(key)[isdyn]
    sv = sv[np.isfinite(sv)]; dv = dv[np.isfinite(dv)]
    data = [sv if sv.size else np.array([np.nan]), dv if dv.size else np.array([np.nan])]
    bp = ax[k].boxplot(data, positions=[0, 1], widths=0.6, showfliers=False, patch_artist=True,
                       medianprops=dict(color="#222", lw=1))
    for patch, cc in zip(bp["boxes"], ["#5A6B8C", "#C44E52"]): patch.set_facecolor(cc); patch.set_alpha(.85)
    ax[k].set_xticks([0, 1]); ax[k].set_xticklabels([f"Stable\n(n={sv.size})", f"Dynamic\n(n={dv.size})"])
    ax[k].set_title(key)
fig.suptitle("Full-ISM minimal-edit counterfactual (dynamic vs stable)", fontsize=9)
fig.tight_layout(); fig.savefig(OUT / "pub" / "FigE_ism_counterfactual.png", bbox_inches="tight")
print(f"\nSaved: {OUT/'ism_counterfactual.json'} and pub/FigE_ism_counterfactual.png  (total {time.time()-t0:.0f}s)", flush=True)
