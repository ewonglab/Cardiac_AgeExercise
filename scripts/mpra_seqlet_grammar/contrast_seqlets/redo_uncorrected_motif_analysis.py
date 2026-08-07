#!/usr/bin/env python3
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Analyse motif occupancy while retaining native sequence composition.
# Steps: 1) Load sequences and module outputs; 2) scan JASPAR motifs; 3) test stable-versus-dynamic and seqlet enrichment.
# =======================================
"""Primary motif analysis treating sequence composition as biological signal."""

from pathlib import Path
import json
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from config.paths import publication_path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(publication_path("MPRA_SEQLET_TASK_ROOT"))
OUT = ROOT / "outputs"
PROJECT = Path(publication_path("MPRA_KOO_PROJECT"))
MEME = PROJECT / "koorun/out/jaspar_panel.meme"
RELATIVE_THRESHOLD = 0.80
BASES = np.array(list("ACGT"))
BASE_INDEX = {b: i for i, b in enumerate("ACGT")}
CLASS_ORDER = ["Stable", "Age-only", "Exercise-only", "Reversal"]
CLASS_COLORS = ["#6b7280", "#2563eb", "#059669", "#7c3aed"]


def read_fasta(path):
    records, key, buf = {}, None, []
    for raw in path.open():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if key is not None:
                records[key] = "".join(buf)
            kind, coord = line[1:].split("|", 1)
            key, buf = (kind, coord), []
        else:
            buf.append(line)
    if key is not None:
        records[key] = "".join(buf)
    return records


def core(seq):
    return seq[15:-15].upper()


def parse_meme(path):
    lines = path.read_text().splitlines()
    motifs = {}
    i = 0
    while i < len(lines):
        if lines[i].startswith("MOTIF "):
            name = lines[i].split(maxsplit=1)[1]
            i += 1
            while i < len(lines) and "letter-probability matrix" not in lines[i]:
                i += 1
            width = int(re.search(r"\bw=\s*(\d+)", lines[i]).group(1))
            motifs[name] = np.asarray([[float(x) for x in row.split()[:4]]
                                       for row in lines[i + 1:i + 1 + width]], float)
            i += width
        i += 1
    return motifs


def motif_gene(name):
    return re.sub(r"_MA\d+\.\d+$", "", name)


def bh(pvalues):
    p = np.asarray(pvalues, float)
    order = np.argsort(p)
    ranked = p[order]
    qrank = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(qrank)
    out[order] = np.minimum(qrank, 1)
    return out


def scan_matrix(encoded, pwm):
    """Return best score, relative score, hit count, strand and offset per sequence."""
    logodds = np.log2((pwm + 1e-5) / 0.25)
    n, length = encoded.shape
    width = len(pwm)
    if width > length:
        return None
    minimum = float(np.min(logodds, axis=1).sum())
    maximum = float(np.max(logodds, axis=1).sum())
    threshold = minimum + RELATIVE_THRESHOLD * (maximum - minimum)
    best = np.full(n, -np.inf)
    best_strand = np.full(n, "+", dtype="U1")
    best_offset = np.zeros(n, int)
    hit_count = np.zeros(n, int)
    rc = logodds[::-1, ::-1]
    arange = np.arange(width)[None, :]
    for start in range(length - width + 1):
        plus = logodds[arange, encoded[:, start:start + width]].sum(axis=1)
        minus = rc[arange, encoded[:, start:start + width]].sum(axis=1)
        at_position = np.maximum(plus, minus)
        hit_count += at_position >= threshold
        take = at_position > best
        best[take] = at_position[take]
        best_offset[take] = start
        best_strand[take] = np.where(plus[take] >= minus[take], "+", "-")
    relative = (best - minimum) / max(1e-9, maximum - minimum)
    return best, relative, hit_count, best_strand, best_offset


def encode(sequences):
    return np.asarray([[BASE_INDEX[b] for b in s] for s in sequences], dtype=np.int8)


def resolved_class(groups):
    if groups == {"Stable"}:
        return "Stable"
    if ("Age Up" in groups and "Ex Down" in groups) or ("Age Down" in groups and "Ex Up" in groups):
        return "Reversal"
    if groups & {"Age Up", "Age Down"}:
        return "Age-only"
    if groups & {"Ex Up", "Ex Down"}:
        return "Exercise-only"
    return "Other"


def main():
    motifs = parse_meme(MEME)
    fa = read_fasta(PROJECT / "bundle/MPRA_AP1_fibroblast.fasta")
    source_all = pd.read_csv(PROJECT / "bundle/MPRA_Results_with_TSS_distance.csv")
    classes = source_all.groupby("Coords")["Group"].agg(set).map(resolved_class)
    source = source_all.drop_duplicates("Coords").set_index("Coords")
    source["resolved_class"] = classes
    source["is_dynamic"] = source["resolved_class"] != "Stable"
    source["measured_context"] = source["log2_alpha_MotifShuff"] - source["log2_alpha_FullShuff"]
    source["measured_total"] = source["log2_alpha_Endo"] - source["log2_alpha_FullShuff"]

    # Treat GC as a biological phenotype rather than a nuisance covariate.
    stable = source[source.resolved_class == "Stable"]
    dynamic = source[source.resolved_class != "Stable"]
    gc_stats = {
        "stable_n": int(len(stable)),
        "dynamic_n": int(len(dynamic)),
        "stable_mean_GC": float(stable.GC_content.mean()),
        "dynamic_mean_GC": float(dynamic.GC_content.mean()),
        "dynamic_minus_stable_GC": float(dynamic.GC_content.mean() - stable.GC_content.mean()),
        "mannwhitney_p": float(stats.mannwhitneyu(dynamic.GC_content, stable.GC_content).pvalue),
        "class_kruskal_p": float(stats.kruskal(*[source.loc[source.resolved_class == c, "GC_content"]
                                                  for c in CLASS_ORDER]).pvalue),
        "spearman_GC_vs_context": float(stats.spearmanr(source.GC_content, source.measured_context).statistic),
        "spearman_GC_vs_context_p": float(stats.spearmanr(source.GC_content, source.measured_context).pvalue),
        "spearman_GC_vs_total": float(stats.spearmanr(source.GC_content, source.measured_total).statistic),
        "spearman_GC_vs_total_p": float(stats.spearmanr(source.GC_content, source.measured_total).pvalue),
    }

    # A. Uncorrected motif assignment in the exact 12-bp high modules.
    pair = pd.read_csv(OUT / "context_module_pair_results.csv")
    module_parts = []
    for side in [1, 2]:
        z = pair[["Coords", f"module_{side}", f"module_{side}_start", f"module_{side}_end",
                  f"module_{side}_sequence", f"module_{side}_direct_effect",
                  f"module_{side}_attribution_rank"]].copy()
        z.columns = ["Coords", "letter", "start_0based", "end_0based_exclusive", "sequence",
                     "direct_absolute_effect", "original_attribution_rank"]
        module_parts.append(z)
    modules = pd.concat(module_parts, ignore_index=True).drop_duplicates(["Coords", "letter"]).reset_index(drop=True)
    module_encoded = encode(modules.sequence.astype(str).tolist())
    per_motif = []
    for jaspar, pwm in motifs.items():
        scanned = scan_matrix(module_encoded, pwm)
        if scanned is None:
            continue
        raw, relative, count, strand, offset = scanned
        for i in range(len(modules)):
            per_motif.append({
                "Coords": modules.loc[i, "Coords"], "letter": modules.loc[i, "letter"],
                "sequence": modules.loc[i, "sequence"], "motif": motif_gene(jaspar),
                "jaspar_motif": jaspar, "raw_log_odds_score": float(raw[i]),
                "relative_PWM_score": float(relative[i]), "hit_count_at_80pct": int(count[i]),
                "strand": strand[i], "motif_offset_in_module": int(offset[i]),
                "passes_80pct_relative_threshold": bool(relative[i] >= RELATIVE_THRESHOLD),
            })
    module_matches = pd.DataFrame(per_motif)
    module_matches.to_csv(OUT / "context_module_individual_motif_matches_uncorrected.csv", index=False)
    assignments = []
    for (coord, letter), z in module_matches.groupby(["Coords", "letter"]):
        z = z.sort_values("relative_PWM_score", ascending=False)
        passing = z[z.passes_80pct_relative_threshold]
        top = z.iloc[0]
        names = list(dict.fromkeys(passing.motif.astype(str)))
        base = modules[(modules.Coords == coord) & (modules.letter == letter)].iloc[0]
        assignments.append({
            "Coords": coord, "letter": letter, "start_0based": int(base.start_0based),
            "end_0based_exclusive": int(base.end_0based_exclusive), "sequence": base.sequence,
            "direct_absolute_effect": float(base.direct_absolute_effect),
            "top_motif": top.motif, "top_jaspar_motif": top.jaspar_motif,
            "top_relative_PWM_score": float(top.relative_PWM_score), "top_strand": top.strand,
            "top_motif_offset": int(top.motif_offset_in_module),
            "n_motifs_passing_80pct": int(len(passing)),
            "motifs_passing_80pct": ";".join(names),
            "motif_assignment": str(top.motif) if top.relative_PWM_score >= RELATIVE_THRESHOLD
                                else f"unresolved (top {top.motif}, relative={top.relative_PWM_score:.2f})",
        })
    assignments = pd.DataFrame(assignments)
    assignments.to_csv(OUT / "context_module_motif_assignments_uncorrected.csv", index=False)

    annotated_pairs = pair.copy()
    for side in [1, 2]:
        lookup = assignments.add_prefix(f"module_{side}_motif_").rename(columns={
            f"module_{side}_motif_Coords": "Coords", f"module_{side}_motif_letter": f"module_{side}"})
        annotated_pairs = annotated_pairs.merge(lookup, on=["Coords", f"module_{side}"], how="left")
    annotated_pairs["motif_pair_name"] = (annotated_pairs.module_1_motif_motif_assignment.astype(str)
                                           + " + " + annotated_pairs.module_2_motif_motif_assignment.astype(str))
    annotated_pairs.to_csv(OUT / "context_module_pairs_with_motifs_uncorrected.csv", index=False)
    annotated_pairs.head(30).to_csv(OUT / "context_module_pair_top_hits_with_motifs_uncorrected.csv", index=False)

    # B. Raw individual-PWM enrichment in high modules versus low controls.
    candidate = pd.read_csv(OUT / "context_probe_candidate_modules.csv")
    local_sequences = []
    for r in candidate.itertuples(index=False):
        seq = core(fa[("MOTIF_SHUF", str(r.Coords))])
        center = (int(r.start_0based) + int(r.end_0based_exclusive)) // 2
        lo, hi = max(0, center - 11), min(170, center + 11)
        local_sequences.append(seq[lo:hi])
    lengths = pd.Series([len(x) for x in local_sequences])
    enrichment_rows = []
    for jaspar, pwm in motifs.items():
        raw = np.full(len(candidate), np.nan)
        for length in sorted(lengths.unique()):
            idx = np.flatnonzero(lengths.to_numpy() == length)
            scanned = scan_matrix(encode([local_sequences[i] for i in idx]), pwm)
            if scanned is not None:
                raw[idx] = scanned[0]
        tmp = candidate[["Coords", "module_type"]].copy()
        tmp["score"] = raw
        by = tmp.groupby(["Coords", "module_type"])["score"].mean().unstack().dropna()
        test = stats.wilcoxon(by.high, by.low_control, alternative="greater")
        enrichment_rows.append({
            "motif": motif_gene(jaspar), "jaspar_motif": jaspar, "n_loci": int(len(by)),
            "high_mean_raw_PWM": float(by.high.mean()), "low_mean_raw_PWM": float(by.low_control.mean()),
            "high_minus_low": float(by.high.mean() - by.low_control.mean()),
            "paired_wilcoxon_greater_p": float(test.pvalue),
        })
    enrichment = pd.DataFrame(enrichment_rows)
    enrichment["BH_q"] = bh(enrichment.paired_wilcoxon_greater_p)
    enrichment = enrichment.sort_values(["BH_q", "paired_wilcoxon_greater_p"])
    enrichment.to_csv(OUT / "motif_enrichment_functional_modules_uncorrected.csv", index=False)
    # The prespecified 15-family panel is also reported. Family-level testing is
    # less redundant than testing all 41 closely related individual matrices.
    family_enrichment = pd.read_csv(OUT / "context_module_motif_enrichment.csv")
    family_enrichment = (family_enrichment[family_enrichment.analysis == "raw_PWM"]
                         .dropna(subset=["family"])[["family", "n_loci", "high_mean", "low_mean",
                                                    "paired_wilcoxon_greater_p", "BH_q"]]
                         .sort_values(["BH_q", "paired_wilcoxon_greater_p"]))
    family_enrichment.to_csv(OUT / "motif_family_enrichment_functional_modules_uncorrected.csv", index=False)

    # C. Raw motif occupancy in complete native inserts: dynamic versus stable.
    coords = source.index.astype(str).tolist()
    native_encoded = encode([core(fa[("ORIG", c)]) for c in coords])
    dynamic_mask = source.loc[coords, "is_dynamic"].to_numpy(bool)
    occupancy_rows = []
    occupancy_long = []
    for jaspar, pwm in motifs.items():
        raw, relative, count, _, _ = scan_matrix(native_encoded, pwm)
        d, s = count[dynamic_mask], count[~dynamic_mask]
        test = stats.mannwhitneyu(d, s, alternative="two-sided")
        occupancy_rows.append({
            "motif": motif_gene(jaspar), "jaspar_motif": jaspar,
            "dynamic_mean_hits": float(d.mean()), "stable_mean_hits": float(s.mean()),
            "dynamic_minus_stable_hits": float(d.mean() - s.mean()),
            "dynamic_presence_fraction": float(np.mean(d > 0)),
            "stable_presence_fraction": float(np.mean(s > 0)),
            "mannwhitney_two_sided_p": float(test.pvalue),
        })
    occupancy = pd.DataFrame(occupancy_rows)
    occupancy["BH_q"] = bh(occupancy.mannwhitney_two_sided_p)
    occupancy["enriched_in"] = np.where(occupancy.dynamic_minus_stable_hits > 0, "dynamic", "stable")
    occupancy = occupancy.sort_values(["BH_q", "mannwhitney_two_sided_p"])
    occupancy.to_csv(OUT / "motif_enrichment_dynamic_vs_stable_uncorrected.csv", index=False)

    summary = {
        "primary_analysis_policy": "No sequence-composition adjustment; GC and oligonucleotide composition are treated as biological signal.",
        "relative_PWM_presence_threshold": RELATIVE_THRESHOLD,
        "gc_biology": gc_stats,
        "functional_module_families_BH_q_lt_0_05": family_enrichment[family_enrichment.BH_q < .05].to_dict("records"),
        "functional_module_motifs_BH_q_lt_0_05": enrichment[enrichment.BH_q < .05].to_dict("records"),
        "dynamic_vs_stable_native_motifs_BH_q_lt_0_05": occupancy[occupancy.BH_q < .05].to_dict("records"),
    }
    with open(OUT / "motif_analysis_uncorrected_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    # Figure: GC biology, raw module enrichment, native dynamic/stable occupancy,
    # and motif names on the strongest affected pairs.
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42})
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.8), constrained_layout=True)
    axa, axb, axc, axd = axes.ravel()
    gc_data = [source.loc[source.resolved_class == c, "GC_content"] for c in CLASS_ORDER]
    bp = axa.boxplot(gc_data, labels=[x.replace("-only", "") for x in CLASS_ORDER], patch_artist=True,
                     showfliers=False)
    for box, color in zip(bp["boxes"], CLASS_COLORS):
        box.set_facecolor(color); box.set_alpha(.75)
    axa.set_ylabel("Native insert GC fraction")
    axa.set_title("Lower GC is part of the dynamic-enhancer phenotype")
    axa.text(.98, .03, f"dynamic − stable = {gc_stats['dynamic_minus_stable_GC']:.3f}\nP = {gc_stats['mannwhitney_p']:.1e}",
             transform=axa.transAxes, ha="right", va="bottom")

    top_enrich = family_enrichment.copy()
    top_enrich["high_minus_low"] = top_enrich["high_mean"] - top_enrich["low_mean"]
    top_enrich = top_enrich.sort_values("high_minus_low").tail(12)
    colors = ["#d95f02" if q < .05 else "#9ca3af" for q in top_enrich.BH_q]
    axb.barh(np.arange(len(top_enrich)), top_enrich.high_minus_low, color=colors)
    axb.set_yticks(np.arange(len(top_enrich)), top_enrich.family, fontsize=7.5)
    axb.axvline(0, color="#6b7280", lw=.8)
    axb.set_xlabel("Raw PWM score: high module − low control")
    axb.set_title("Uncorrected motif-family enrichment in functional modules")
    nkx_row = top_enrich[top_enrich.family == "NKX2-5"]
    if not nkx_row.empty:
        idx = int(np.flatnonzero(top_enrich.family.to_numpy() == "NKX2-5")[0])
        axb.text(float(nkx_row.high_minus_low.iloc[0]) + .02, idx,
                 f"q={float(nkx_row.BH_q.iloc[0]):.3f}", va="center", fontsize=7.5)

    sig_occ = occupancy[occupancy.BH_q < .05].copy()
    if len(sig_occ) < 12:
        sig_occ = occupancy.head(12).copy()
    else:
        sig_occ = sig_occ.head(16)
    sig_occ = sig_occ.sort_values("dynamic_minus_stable_hits")
    colors = ["#7c3aed" if x > 0 else "#6b7280" for x in sig_occ.dynamic_minus_stable_hits]
    axc.barh(np.arange(len(sig_occ)), sig_occ.dynamic_minus_stable_hits, color=colors)
    axc.set_yticks(np.arange(len(sig_occ)), sig_occ.motif, fontsize=7.5)
    axc.axvline(0, color="#374151", lw=.8)
    axc.set_xlabel("Mean motif hits: dynamic − stable")
    axc.set_title("Native motif occupancy differs with enhancer class")

    top_pairs = annotated_pairs.head(12).sort_values("mean_absolute_epsilon")
    labels = [f"{r.Coords} {r.interaction}\n{r.motif_pair_name}" for _, r in top_pairs.iterrows()]
    axd.barh(np.arange(len(top_pairs)), top_pairs.mean_absolute_epsilon, color="#d95f02", alpha=.8)
    axd.set_yticks(np.arange(len(top_pairs)), labels, fontsize=6.2)
    axd.set_xlabel("Mean |epistasis ε|")
    axd.set_title("Uncorrected candidate motif names on top pairs")
    for label, ax in zip("ABCD", axes.ravel()):
        ax.text(-.14, 1.07, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")
    fig.suptitle("Enhancer composition and motif content as biological signal", fontsize=14, fontweight="bold")
    fig.savefig(OUT / "Fig_motif_analysis_uncorrected.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "Fig_motif_analysis_uncorrected.pdf", bbox_inches="tight")

    # Revised motif-pair report without composition correction.
    lines = [
        "# Uncorrected motif names for impacted enhancer-module pairs", "",
        "## Analysis policy", "",
        "Sequence composition is treated as part of enhancer biology. No GC, mononucleotide or dinucleotide correction is applied. Motif presence is called using the raw JASPAR PWM and an 80% relative-score threshold. Multiple-testing correction across motif enrichment tests is retained because it addresses the number of hypotheses, not sequence composition.", "",
        f"Dynamic loci are {abs(gc_stats['dynamic_minus_stable_GC'])*100:.2f} percentage points lower in GC than stable loci ({gc_stats['dynamic_mean_GC']:.3f} versus {gc_stats['stable_mean_GC']:.3f}; P = {gc_stats['mannwhitney_p']:.2g}). Lower GC is also associated with larger context dependence (Spearman ρ = {gc_stats['spearman_GC_vs_context']:.3f}, P = {gc_stats['spearman_GC_vs_context_p']:.2g}).", "",
        "## Strongest affected pairs with uncorrected motif names", "",
        "| Rank | Locus | Class | Pair | Motif pair | Gap | Mean |ε| | Rescue |",
        "|---:|---|---|---|---|---:|---:|---:|",
    ]
    for _, r in annotated_pairs.head(25).iterrows():
        lines.append(f"| {int(r.epistasis_rank)} | {r.Coords} | {r.resolved_class} | {r.interaction} | {r.motif_pair_name} | {int(r.edge_to_edge_gap_bp)} bp | {r.mean_absolute_epsilon:.3f} | {r.rescue_high_gap_closed:.1%} |")
    lines += ["", "## Raw motif enrichment in functional modules", "",
              f"In the prespecified 15-family analysis, NKX2-5 is enriched in high-impact modules (raw paired P = {family_enrichment.iloc[0].paired_wilcoxon_greater_p:.3g}, BH q = {family_enrichment.iloc[0].BH_q:.3g}). In the more redundant 41-individual-matrix panel the same NKX2-5 signal is P = {enrichment.iloc[0].paired_wilcoxon_greater_p:.3g}, q = {enrichment.iloc[0].BH_q:.3g}.", "",
              "| Motif | High mean PWM | Low mean PWM | Difference | BH q |",
              "|---|---:|---:|---:|---:|"]
    for _, r in enrichment.head(15).iterrows():
        lines.append(f"| {r.motif} | {r.high_mean_raw_PWM:.3f} | {r.low_mean_raw_PWM:.3f} | {r.high_minus_low:.3f} | {r.BH_q:.3g} |")
    lines += ["", "## Interpretation", "",
              "The reported motif and GC differences are total biological associations. They may reflect TF recognition, DNA shape, nucleosome preference and broader enhancer composition together. Direct motif substitutions in the next MPRA remain necessary to distinguish these mechanisms, but the motif signals are not removed from the primary result.", ""]
    report_text = "\n".join(lines)
    (OUT / "context_module_pair_motif_report_uncorrected.md").write_text(report_text)
    # Replace the earlier motif-pair report so the primary report follows the requested policy.
    (OUT / "context_module_pair_motif_report.md").write_text(report_text)
    print(f"uncorrected module assignments={len(assignments)}; functional family q<.05={(family_enrichment.BH_q < .05).sum()}; native class motif q<.05={(occupancy.BH_q < .05).sum()}")


if __name__ == "__main__":
    main()
