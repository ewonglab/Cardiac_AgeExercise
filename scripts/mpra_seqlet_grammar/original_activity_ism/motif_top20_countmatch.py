# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Assess motif reliance in the top 20 ISM positions with site-count matching.
# Steps: 1) Load full ISM arrays and sequences; 2) scan motifs; 3) standardise dynamic-versus-stable effects and export statistics.
# =======================================
"""(b) Count-matched top-20 motif comparison, pooled Dynamic vs Stable. For each TF we ask whether
dynamic enhancers RELY on that motif more/less among their top-20 ISM positions, controlling for how
many of that TF's sites the enhancer actually contains (so we measure reliance-given-presence, not
mere motif abundance). Metric per enhancer = # of its top-20 ISM positions inside that TF's
footprint. Adjusted difference = direct standardization over the TF's own site-count strata (0,1,2,3+);
bootstrap 95% CI (5000). Raw difference tested with two-sided Mann-Whitney U. AP-1 also reported as
binary >=1-in-top-20 (extends FigF2 from top-5 to top-20)."""
import re, json
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
from Bio.Seq import Seq
from pyjaspar import jaspardb
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from config.paths import publication_path

MPRA = Path(publication_path("MPRA_DATA_ROOT"))
OUT = Path(publication_path("MPRA_SEQLET_TASK_ROOT", "outputs", "original_activity_ism")); AP1_RE = re.compile("TGA[CG]TCA"); L = 170; K = 20; NBOOT = 5000
PANEL = ["RUNX1", "RUNX2", "ETS1", "ELK1", "FLI1", "NFATC1", "NFATC2", "TEAD1", "GATA4",
         "CEBPA", "CEBPB", "CEBPD", "SMAD3", "RELA", "EGR1", "TCF21", "SRF", "MEF2C"]; THR = 0.80

d = np.load(OUT / "ism_imp_full.npz", allow_pickle=True)   # full 1653 ISM (was 310 subset)
IMP = d["imp"]; groups = d["groups"].astype(str); coords = d["coords"].astype(str)
seqmap = {}; nm = None
for l in open(MPRA / "deep_learning" / "results" / "sequences_cache.fa"):
    if l.startswith(">"): nm = l.strip()[1:]
    else: seqmap[nm] = l.strip().upper()
S = [seqmap[c] for c in coords]; isdyn = groups != "Stable"
from collections import Counter
c = Counter("".join(S)); tot = sum(c[b] for b in "ACGT"); bg = {b: c[b]/tot for b in "ACGT"}
jdb = jaspardb(release="JASPAR2024"); pssms = {}
for k in PANEL:
    try: ms = jdb.fetch_motifs_by_name(k)
    except Exception: ms = None
    if ms: pssms[k] = ms[0].counts.normalize(pseudocounts=0.5).log_odds(bg)
TFS = ["AP1"] + list(pssms)

def ap1_cover(s):
    cov = np.zeros(L, bool)
    for m in AP1_RE.finditer(s): cov[m.start():m.end()] = True
    rc = s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    for m in AP1_RE.finditer(rc): cov[L-m.end():L-m.start()] = True
    return cov
def pwm_cover(ps, s):
    Lm = ps.length; smin, smax = ps.min, ps.max; rng = smax - smin; sq = Seq(s)
    both = np.nanmax(np.vstack([np.asarray(ps.calculate(sq), float),
                                np.asarray(ps.reverse_complement().calculate(sq), float)]), 0)
    rel = (both - smin) / rng; cov = np.zeros(L, bool)
    for i in range(len(rel)):
        if np.isfinite(rel[i]) and rel[i] >= THR: cov[i:i+Lm] = True
    return cov
def nseg(cov): return int((cov[0]) + np.sum(cov[1:] & ~cov[:-1]))   # # of binding sites (coverage runs)

count = {tf: np.zeros(len(S), int) for tf in TFS}      # # top-20 positions in TF
sites = {tf: np.zeros(len(S), int) for tf in TFS}      # # of TF sites in enhancer
for j, s in enumerate(S):
    cov = {"AP1": ap1_cover(s)}
    for tf in pssms: cov[tf] = pwm_cover(pssms[tf], s)
    top = np.argsort(-IMP[j])[:K]
    for tf in TFS: count[tf][j] = int(cov[tf][top].sum()); sites[tf][j] = nseg(cov[tf])
    if (j+1) % 100 == 0: print(f"  scanned {j+1}/{len(S)}", flush=True)

def std_diff(metric, dyn, st, use):
    num, den = 0.0, 0
    for lv in use:
        m = st == lv
        if not m.any(): continue
        ds, dd = m & ~dyn, m & dyn
        if not ds.any() or not dd.any(): continue
        num += m.sum() * (metric[dd].mean() - metric[ds].mean()); den += m.sum()
    return num/den if den else np.nan

boot = np.random.default_rng(42); n = len(S)
res = {"K": K, "n": n, "n_stable": int((~isdyn).sum()), "n_dynamic": int(isdyn.sum()), "TFs": {}}
print(f"\nn={n} (stable {(~isdyn).sum()}, dynamic {isdyn.sum()}); count-matched top-{K} reliance (Dynamic - Stable)")
print(f"{'TF':8s} {'rawS':>5s} {'rawD':>5s} {'rawDiff':>7s} {'MWUp':>6s} | {'adjDiff':>7s} {'95% CI':>16s} {'p~':>5s}")
for tf in TFS:
    metric = count[tf].astype(float); st = np.minimum(sites[tf], 3)
    rs, rd = metric[~isdyn].mean(), metric[isdyn].mean()
    mwu = stats.mannwhitneyu(metric[isdyn], metric[~isdyn]).pvalue
    strata = [lv for lv in [0, 1, 2, 3] if ((st == lv) & ~isdyn).any() and ((st == lv) & isdyn).any()]
    adj = std_diff(metric, isdyn, st, strata)
    bs = np.empty(NBOOT)
    for b in range(NBOOT):
        ix = boot.integers(0, n, n); bs[b] = std_diff(metric[ix], isdyn[ix], st[ix], strata)
    bs = bs[np.isfinite(bs)]; lo, hi = np.percentile(bs, [2.5, 97.5]); fgt = float(np.mean(bs > 0))
    p2 = 2 * min(fgt, 1 - fgt)
    EPS = 1e-9
    adj_sig = bool(bs.std() > EPS and (lo > EPS or hi < -EPS))  # CI strictly excludes 0, non-degenerate
    res["TFs"][tf] = dict(raw_stable=float(rs), raw_dynamic=float(rd), raw_diff=float(rd - rs), mwu_p=float(mwu),
                          adj_diff=float(adj), adj_ci=[float(lo), float(hi)], adj_p_approx=float(p2),
                          adj_sig=adj_sig, strata_used=strata)
    flag = " <-" if adj_sig else ""
    print(f"{tf:8s} {rs:5.2f} {rd:5.2f} {rd-rs:+7.2f} {mwu:6.3g} | {adj:+7.2f} [{lo:+.2f},{hi:+.2f}] {p2:5.2g}{flag}")

# AP-1 binary >=1 in top-20, count-matched (extends FigF2 to top-20)
apb = (count["AP1"] >= 1).astype(float); stA = np.minimum(sites["AP1"], 3)
strataA = [lv for lv in [1, 2, 3] if ((stA == lv) & ~isdyn).any() and ((stA == lv) & isdyn).any()]
adjA = std_diff(apb, isdyn, stA, strataA); bsA = np.empty(NBOOT)
for b in range(NBOOT):
    ix = boot.integers(0, n, n); bsA[b] = std_diff(apb[ix], isdyn[ix], stA[ix], strataA)
bsA = bsA[np.isfinite(bsA)]; ciA = np.percentile(bsA, [2.5, 97.5])
res["AP1_binary_top20"] = dict(raw_stable=float(apb[~isdyn].mean()), raw_dynamic=float(apb[isdyn].mean()),
                               adj_diff=float(adjA), adj_ci=[float(ciA[0]), float(ciA[1])],
                               frac_gt0=float(np.mean(bsA > 0)))
print(f"\nAP-1 binary >=1 in top-20 (count-matched): stable {apb[~isdyn].mean():.0%} dyn {apb[isdyn].mean():.0%}  "
      f"adj {adjA:+.3f} [{ciA[0]:+.3f},{ciA[1]:+.3f}]")
json.dump(res, open(OUT / "motif_top20_countmatch.json", "w"), indent=2)

# forest plot of adjusted differences (count metric), sorted
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
                     "font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 300})
order = sorted(TFS, key=lambda t: res["TFs"][t]["adj_diff"])
y = np.arange(len(order))
fig, ax = plt.subplots(figsize=(5.6, 0.34*len(order)+1))
ax.axvline(0, color="#333", ls="--", lw=.9)
for yi, tf in zip(y, order):
    r = res["TFs"][tf]; col = "#B0302E" if r["adj_sig"] else "#8a8a8a"
    ax.plot(r["adj_ci"], [yi, yi], color=col, lw=1.5); ax.plot(r["adj_diff"], yi, "o", color=col, ms=5)
ax.set_yticks(y); ax.set_yticklabels(order, fontsize=7.5)
ax.set_xlabel("Dynamic − Stable, count-matched\n(mean # of top-20 ISM positions in motif)")
ax.set_title(f"Count-matched top-{K} motif reliance (bootstrap 95% CI)\nred = 95% CI excludes 0",
             loc="left", fontweight="bold", fontsize=8.5)
fig.tight_layout()
for e in ("pdf", "png"): fig.savefig(OUT / "pub" / f"FigF5_top20_countmatched.{e}", bbox_inches="tight")
print("Saved motif_top20_countmatch.json and pub/FigF5_top20_countmatched.png")
