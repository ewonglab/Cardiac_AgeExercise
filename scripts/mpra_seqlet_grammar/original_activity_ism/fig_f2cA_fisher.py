# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Test AP-1 reliance after stratifying ISM results by AP-1 site count.
# Steps: 1) Load balanced ISM results; 2) run Fisher and CMH tests with bootstrapping; 3) write the figure panel and statistics.
# =======================================
"""FigF2cA (the easy-to-read panel): fraction with AP-1 in top-5 ISM within each AP-1-site-count
stratum, Stable vs Dynamic, bootstrap 95% CI error bars, with a per-stratum Fisher exact p on top
(exact 2x2 test, class x AP-1-in-top5). Pooled AP-1-count-adjusted effect summarised by the
Cochran-Mantel-Haenszel stratified test (closed-form) alongside the bootstrap standardized diff."""
import re, json
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from config.paths import publication_path

MPRA = Path(publication_path("MPRA_DATA_ROOT"))
OUT = Path(publication_path("MPRA_SEQLET_TASK_ROOT", "outputs", "original_activity_ism")); AP1_RE = re.compile("TGA[CG]TCA"); L = 170
SEED = 42; STAB_N, SUB_N = 150, 40; NBOOT = 5000

df = pd.read_csv(MPRA / "corrected_data" / "MPRA_Results_with_TSS_distance.csv").drop_duplicates("Coords").reset_index(drop=True)
seqs = {}; nm = None
for l in open(MPRA / "deep_learning" / "results" / "sequences_cache.fa"):
    if l.startswith(">"): nm = l.strip()[1:]
    else: seqs[nm] = l.strip().upper()
sequences = [seqs[c] for c in df["Coords"]]; grp = df["Group"].to_numpy()
rng = np.random.default_rng(SEED); sample = []
for g, n in [("Stable", STAB_N), ("Age Up", SUB_N), ("Age Down", SUB_N), ("Ex Up", SUB_N), ("Ex Down", SUB_N)]:
    idx = np.where(grp == g)[0]; sample.extend((int(i), g) for i in rng.choice(idx, min(n, len(idx)), replace=False))
rows = json.load(open(OUT / "ism_counterfactual.json"))["rows"]

def n_ap1(s):
    spans = [(m.start(), m.end()) for m in AP1_RE.finditer(s)]
    rc = s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    spans += [(L - m.end(), L - m.start()) for m in AP1_RE.finditer(rc)]
    spans.sort(); cnt = 0; end = -1
    for a, b in spans:
        if a >= end: cnt += 1; end = b
        elif b > end: end = b
    return cnt

isdyn = np.array([g != "Stable" for _, g in sample])
nap1 = np.array([n_ap1(sequences[i]) for i, _ in sample]); st = np.minimum(nap1, 3)
ap1top = np.array([rows[k]["top_ann"]["ap1"] >= 1 for k in range(len(rows))])
levels = [1, 2]   # informative strata (n_ap1=0 -> all-zero, degenerate)

# per-stratum 2x2 Fisher + bootstrap CI on each fraction
def frac(m): return ap1top[m].mean() if m.any() else np.nan
info = {}
for lv in levels:
    ms, md = (st == lv) & ~isdyn, (st == lv) & isdyn
    sp, sn = int(ap1top[ms].sum()), int((~ap1top[ms]).sum())
    dp, dn = int(ap1top[md].sum()), int((~ap1top[md]).sum())
    OR, p = stats.fisher_exact([[dp, dn], [sp, sn]])   # dynamic vs stable
    info[lv] = dict(f_s=frac(ms), f_d=frac(md), n_s=int(ms.sum()), n_d=int(md.sum()),
                    table=[[dp, dn], [sp, sn]], fisher_or=float(OR), fisher_p=float(p))
    print(f"n_ap1={lv}: stable {sp}/{sp+sn}={info[lv]['f_s']:.0%}  dynamic {dp}/{dp+dn}={info[lv]['f_d']:.0%}  Fisher OR={OR:.2f} p={p:.3g}")

# Cochran-Mantel-Haenszel across informative strata (adjusted for n_ap1)
def cmh(tables):  # tables: list of [[a,b],[c,d]] (exposed=dynamic top row)
    num = 0.0; Evar = 0.0; mh_num = 0.0; mh_den = 0.0
    for (a, b), (c, d) in tables:
        Ni = a + b + c + d
        if Ni == 0: continue
        num += a - (a + b) * (a + c) / Ni
        Evar += (a + b) * (c + d) * (a + c) * (b + d) / (Ni**2 * (Ni - 1)) if Ni > 1 else 0
        mh_num += a * d / Ni; mh_den += b * c / Ni
    chi = (abs(num) - 0.5)**2 / Evar if Evar > 0 else np.nan
    p = stats.chi2.sf(chi, 1); OR = mh_num / mh_den if mh_den > 0 else np.nan
    return chi, p, OR
chi, p_cmh, or_cmh = cmh([info[lv]["table"] for lv in levels])
print(f"\nCochran-Mantel-Haenszel (adjust for AP-1 count): chi2={chi:.2f}, p={p_cmh:.3g}, MH OR={or_cmh:.2f}")

# bootstrap CI on each fraction + standardized diff (for the bars / cross-check)
def std_diff(at, dyn, s, use):
    num, den = 0.0, 0
    for lv in use:
        m = s == lv
        if not m.any(): continue
        ds, dd = m & ~dyn, m & dyn
        if not ds.any() or not dd.any(): continue
        num += m.sum() * (at[dd].mean() - at[ds].mean()); den += m.sum()
    return num / den if den else np.nan
boot = np.random.default_rng(SEED); n = len(rows)
bS = {lv: [] for lv in levels}; bD = {lv: [] for lv in levels}; b_pos = np.empty(NBOOT)
for j in range(NBOOT):
    ix = boot.integers(0, n, n); at, dy, sn = ap1top[ix], isdyn[ix], st[ix]
    b_pos[j] = std_diff(at, dy, sn, [1, 2, 3])
    for lv in levels:
        bS[lv].append(at[(sn == lv) & ~dy].mean() if ((sn == lv) & ~dy).any() else np.nan)
        bD[lv].append(at[(sn == lv) & dy].mean() if ((sn == lv) & dy).any() else np.nan)
def ci(v): v = np.asarray(v); v = v[np.isfinite(v)]; return np.percentile(v, [2.5, 97.5])
ciS = {lv: ci(bS[lv]) for lv in levels}; ciD = {lv: ci(bD[lv]) for lv in levels}; ci_pos = ci(b_pos)
std_pos = std_diff(ap1top, isdyn, st, [1, 2, 3])
json.dump(dict(strata={str(lv): info[lv] for lv in levels},
               cmh=dict(chi2=float(chi), p=float(p_cmh), mh_or=float(or_cmh)),
               bootstrap_std_diff=dict(est=float(std_pos), ci=[float(ci_pos[0]), float(ci_pos[1])])),
          open(OUT / "fig_f2cA_fisher.json", "w"), indent=2)

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
                     "font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 300})
C_S, C_D = "#5A6B8C", "#C0504D"
fig, ax = plt.subplots(figsize=(4.6, 3.8))
xc = np.arange(len(levels)); tags = [f"{lv} site" + ("s" if lv > 1 else "") for lv in levels]
for offx, fr, cis, col, lab in [(-0.2, [info[lv]["f_s"] for lv in levels], [ciS[lv] for lv in levels], C_S, "Stable"),
                                (0.2, [info[lv]["f_d"] for lv in levels], [ciD[lv] for lv in levels], C_D, "Dynamic")]:
    los = [fr[i] - cis[i][0] for i in range(len(levels))]; his = [cis[i][1] - fr[i] for i in range(len(levels))]
    ax.bar(xc + offx, fr, 0.38, color=col, alpha=.85, label=lab, yerr=[los, his], error_kw=dict(lw=1, capsize=3))
ymax = 1.12
for i, lv in enumerate(levels):
    y = max(info[lv]["f_s"] + (ciS[lv][1] - info[lv]["f_s"]), info[lv]["f_d"] + (ciD[lv][1] - info[lv]["f_d"])) + 0.04
    ax.plot([xc[i]-0.2, xc[i]-0.2, xc[i]+0.2, xc[i]+0.2], [y, y+.02, y+.02, y], lw=.8, color="#333")
    ax.text(xc[i], y+.02, f"Fisher p={info[lv]['fisher_p']:.3g}\n(n={info[lv]['n_s']}/{info[lv]['n_d']})",
            ha="center", va="bottom", fontsize=6.6)
ax.set_xticks(xc); ax.set_xticklabels(tags); ax.set_ylim(0, ymax)
ax.set_xlabel("AP-1 sites in enhancer (matched)"); ax.set_ylabel("Fraction with AP-1 in top-5 ISM")
ax.legend(frameon=False, fontsize=7.5, loc="upper left")
ax.set_title(f"AP-1 reliance, count-matched\nper-stratum Fisher; pooled CMH p={p_cmh:.3g} (OR={or_cmh:.2f})",
             loc="left", fontweight="bold", fontsize=8.5)
fig.tight_layout()
for e in ("pdf", "png"): fig.savefig(OUT / "pub" / f"FigF2cA_fisher.{e}", bbox_inches="tight")
print(f"\nbootstrap standardized diff (cross-check) {std_pos:+.3f} CI[{ci_pos[0]:+.3f},{ci_pos[1]:+.3f}]")
print("Saved fig_f2cA_fisher.json and pub/FigF2cA_fisher.png")
