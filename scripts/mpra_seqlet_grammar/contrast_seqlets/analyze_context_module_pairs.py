#!/usr/bin/env python3
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Summarise pair-level context-seqlet epistasis and rescue effects.
# Steps: 1) Load consolidated module results; 2) calculate pair statistics; 3) export class-comparison figures and tables.
# =======================================
"""Resolve pair-level context-module epistasis and rescue."""

from pathlib import Path
import json
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
PAIRS = ["AB", "AC", "BC"]
CLASSES = ["Stable", "Age-only", "Exercise-only", "Reversal"]
HIGH = "#d95f02"
LOW = "#4c78a8"


def bh_adjust(pvalues):
    p = np.asarray(pvalues, float)
    order = np.argsort(p)
    ranked = p[order]
    q = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1)
    return out


def safe_wilcoxon(x, y=None, alternative="two-sided"):
    try:
        return float(stats.wilcoxon(x, y, alternative=alternative).pvalue)
    except ValueError:
        return 1.0


def selected_module_map(summary, detail):
    rows = []
    sequence_lookup = (detail.sort_values("shuffle_replicate")
                       .drop_duplicates(["Coords", "module_type", "start_0based", "end_0based_exclusive"])
                       .set_index(["Coords", "module_type", "start_0based", "end_0based_exclusive"])["original_window"])
    for (coord, module_type), part in summary.groupby(["Coords", "module_type"]):
        chosen = part.nlargest(3, "mean_absolute_effect").reset_index(drop=True)
        if len(chosen) != 3:
            continue
        for i, letter in enumerate("ABC"):
            r = chosen.iloc[i]
            key = (coord, module_type, int(r.start_0based), int(r.end_0based_exclusive))
            rows.append({
                "Coords": coord,
                "module_type": module_type,
                "letter": letter,
                "original_attribution_rank": int(r.module_rank),
                "start_0based": int(r.start_0based),
                "end_0based_exclusive": int(r.end_0based_exclusive),
                "sequence": str(sequence_lookup.loc[key]),
                "direct_absolute_effect": float(r.mean_absolute_effect),
                "attribution_abs_sum": float(r.attribution_abs_sum),
            })
    return pd.DataFrame(rows)


def main():
    epi_raw = pd.read_csv(OUT / "context_module_epistasis.csv")
    rescue = pd.read_csv(OUT / "context_module_rescue.csv")
    disruption = pd.read_csv(OUT / "context_module_disruption_summary.csv")
    detail = pd.read_csv(OUT / "context_module_disruption_all_shuffles.csv")
    modules = selected_module_map(disruption, detail)
    module_idx = modules.set_index(["Coords", "module_type", "letter"])

    # One observation per locus and pair, after averaging the three sequence-shuffle realizations.
    pair_agg = (epi_raw[epi_raw["interaction"].isin(PAIRS)]
                .groupby(["Coords", "resolved_class", "module_type", "interaction"])
                .agg(mean_epsilon=("epsilon_score", "mean"),
                     mean_loss_interaction=("loss_interaction", "mean"),
                     mean_absolute_epsilon=("absolute_epsilon", "mean"),
                     mean_fractional_nonadditivity=("fractional_nonadditivity", "mean"),
                     epistasis_sd=("epsilon_score", "std"),
                     sign_consistency=("loss_interaction", lambda x: max(np.mean(x > 0), np.mean(x < 0))))
                .reset_index())

    high = pair_agg[pair_agg["module_type"] == "high"].copy()
    low = (pair_agg[pair_agg["module_type"] == "low_control"]
           .drop(columns=["resolved_class", "module_type"])
           .rename(columns={c: "low_" + c for c in pair_agg.columns
                            if c not in {"Coords", "resolved_class", "module_type", "interaction"}}))
    pairs = high.merge(low, on=["Coords", "interaction"], how="inner")

    # Map A/B/C to the exact selected 12-bp modules. A is the strongest directly
    # validated module, B the second strongest and C the third strongest.
    annotations = []
    for _, row in pairs.iterrows():
        a, b = row["interaction"]
        ma = module_idx.loc[(row["Coords"], "high", a)]
        mb = module_idx.loc[(row["Coords"], "high", b)]
        gap = max(0, max(ma.start_0based, mb.start_0based) - min(ma.end_0based_exclusive, mb.end_0based_exclusive))
        annotations.append({
            "module_1": a,
            "module_1_start": int(ma.start_0based),
            "module_1_end": int(ma.end_0based_exclusive),
            "module_1_sequence": ma.sequence,
            "module_1_direct_effect": float(ma.direct_absolute_effect),
            "module_1_attribution_rank": int(ma.original_attribution_rank),
            "module_2": b,
            "module_2_start": int(mb.start_0based),
            "module_2_end": int(mb.end_0based_exclusive),
            "module_2_sequence": mb.sequence,
            "module_2_direct_effect": float(mb.direct_absolute_effect),
            "module_2_attribution_rank": int(mb.original_attribution_rank),
            "edge_to_edge_gap_bp": int(gap),
            "center_distance_bp": float(abs((ma.start_0based + ma.end_0based_exclusive) / 2 -
                                            (mb.start_0based + mb.end_0based_exclusive) / 2)),
        })
    pairs = pd.concat([pairs.reset_index(drop=True), pd.DataFrame(annotations)], axis=1)

    # Add pair and component rescue measurements.
    rescue_wide = rescue.pivot_table(index=["Coords", "module_type"], columns="combo",
                                     values=["absolute_gap_closed", "signed_recovery_fraction"])
    rescue_rows = []
    for _, row in pairs.iterrows():
        coord, pair = row.Coords, row.interaction
        a, b = pair
        h = rescue_wide.loc[(coord, "high")]
        l = rescue_wide.loc[(coord, "low_control")]
        rescue_rows.append({
            "rescue_high_gap_closed": float(h["absolute_gap_closed"][pair]),
            "rescue_low_gap_closed": float(l["absolute_gap_closed"][pair]),
            "rescue_high_minus_low": float(h["absolute_gap_closed"][pair] - l["absolute_gap_closed"][pair]),
            "rescue_high_signed": float(h["signed_recovery_fraction"][pair]),
            "rescue_high_increment_over_best_single": float(
                h["absolute_gap_closed"][pair] - max(h["absolute_gap_closed"][a], h["absolute_gap_closed"][b])
            ),
            "rescue_high_nonadditivity": float(
                h["signed_recovery_fraction"][pair]
                - h["signed_recovery_fraction"][a]
                - h["signed_recovery_fraction"][b]
            ),
        })
    pairs = pd.concat([pairs, pd.DataFrame(rescue_rows)], axis=1)
    pairs["interaction_mode"] = np.where(pairs["mean_loss_interaction"] > 0, "synergy", "redundancy")
    cutoff = pairs["mean_absolute_epsilon"].quantile(0.90)
    pairs["robust_top_decile_pair"] = (
        (pairs["mean_absolute_epsilon"] >= cutoff)
        & (pairs["mean_absolute_epsilon"] > pairs["low_mean_absolute_epsilon"])
        & (pairs["sign_consistency"] >= 2 / 3)
    )
    pairs["strong_rescue_pair"] = (
        (pairs["rescue_high_gap_closed"] >= 0.50)
        & (pairs["rescue_high_gap_closed"] > pairs["rescue_low_gap_closed"])
    )
    pairs["epistasis_rank"] = pairs["mean_absolute_epsilon"].rank(method="min", ascending=False).astype(int)
    pairs = pairs.sort_values("epistasis_rank")
    pairs.to_csv(OUT / "context_module_pair_results.csv", index=False)
    pairs.head(30).to_csv(OUT / "context_module_pair_top_hits.csv", index=False)

    # Pair-type statistics.
    pair_stats = {}
    epi_ps = []
    for pair in PAIRS:
        z = pairs[pairs["interaction"] == pair]
        p = safe_wilcoxon(z["mean_absolute_epsilon"], z["low_mean_absolute_epsilon"], "greater")
        epi_ps.append(p)
        class_arrays = [z.loc[z["resolved_class"] == c, "mean_absolute_epsilon"] for c in CLASSES]
        pair_stats[pair] = {
            "n_loci": int(len(z)),
            "high_mean_abs_epsilon": float(z["mean_absolute_epsilon"].mean()),
            "low_mean_abs_epsilon": float(z["low_mean_absolute_epsilon"].mean()),
            "high_low_ratio": float(z["mean_absolute_epsilon"].mean() / z["low_mean_absolute_epsilon"].mean()),
            "high_vs_low_p": p,
            "synergy_fraction": float(np.mean(z["mean_loss_interaction"] > 0)),
            "signed_interaction_mean": float(z["mean_loss_interaction"].mean()),
            "signed_interaction_two_sided_p": safe_wilcoxon(z["mean_loss_interaction"]),
            "class_kruskal_p": float(stats.kruskal(*class_arrays).pvalue),
        }
    for pair, q in zip(PAIRS, bh_adjust(epi_ps)):
        pair_stats[pair]["high_vs_low_BH_q"] = float(q)

    # Repeated-measures comparison of AB, AC and BC among high modules.
    epi_wide = pairs.pivot(index="Coords", columns="interaction", values="mean_absolute_epsilon")[PAIRS]
    friedman = stats.friedmanchisquare(*(epi_wide[p] for p in PAIRS))
    posthoc = []
    for a, b in [("AB", "AC"), ("AB", "BC"), ("AC", "BC")]:
        posthoc.append({"comparison": f"{a}>{b}", "p": safe_wilcoxon(epi_wide[a], epi_wide[b], "greater")})
    posthoc_q = bh_adjust([x["p"] for x in posthoc])
    for x, q in zip(posthoc, posthoc_q):
        x["BH_q"] = float(q)

    # Distance sensitivity. The nearest/farthest comparison is paired within locus.
    distance_stats = {}
    for module_type, prefix in [("high", ""), ("low_control", "low_")]:
        if module_type == "high":
            z = pairs[["Coords", "interaction", "edge_to_edge_gap_bp", "mean_absolute_epsilon"]].copy()
            z.columns = ["Coords", "interaction", "gap", "effect"]
        else:
            # Low windows have different locations, so calculate their distances separately.
            low_rows = []
            low_pairs = pair_agg[pair_agg["module_type"] == "low_control"]
            for _, row in low_pairs.iterrows():
                a, b = row.interaction
                ma = module_idx.loc[(row.Coords, "low_control", a)]
                mb = module_idx.loc[(row.Coords, "low_control", b)]
                gap = max(0, max(ma.start_0based, mb.start_0based) - min(ma.end_0based_exclusive, mb.end_0based_exclusive))
                low_rows.append({"Coords": row.Coords, "interaction": row.interaction,
                                 "gap": gap, "effect": row.mean_absolute_epsilon})
            z = pd.DataFrame(low_rows)
        near = z.sort_values(["Coords", "gap"]).groupby("Coords").first()
        far = z.sort_values(["Coords", "gap"]).groupby("Coords").last()
        distance_stats[module_type] = {
            "nearest_mean_gap_bp": float(near["gap"].mean()),
            "farthest_mean_gap_bp": float(far["gap"].mean()),
            "nearest_mean_abs_epsilon": float(near["effect"].mean()),
            "farthest_mean_abs_epsilon": float(far["effect"].mean()),
            "nearest_gt_farthest_p": safe_wilcoxon(near["effect"], far["effect"], "greater"),
            "fraction_loci_nearest_gt_farthest": float(np.mean(near["effect"] > far["effect"])),
        }

    # Within-locus rank correlation between high-module separation and epistasis.
    gap_wide = pairs.pivot(index="Coords", columns="interaction", values="edge_to_edge_gap_bp")[PAIRS]
    rank_correlations = []
    for coord in epi_wide.index:
        rg = stats.rankdata(gap_wide.loc[coord].to_numpy())
        re = stats.rankdata(epi_wide.loc[coord].to_numpy())
        rank_correlations.append(float(np.corrcoef(rg, re)[0, 1]))
    distance_stats["high_within_locus_rank"] = {
        "mean_spearman": float(np.nanmean(rank_correlations)),
        "median_spearman": float(np.nanmedian(rank_correlations)),
        "fraction_negative": float(np.mean(np.asarray(rank_correlations) < 0)),
    }

    # Distance bins for visualization.
    bins = [-1, 12, 36, 80, 1000]
    labels = ["0–12", "13–36", "37–80", "81+"]
    high_bins = pairs.assign(distance_bin=pd.cut(pairs.edge_to_edge_gap_bp, bins, labels=labels))
    distance_bins = (high_bins.groupby("distance_bin", observed=True)["mean_absolute_epsilon"]
                     .agg(["count", "mean", "sem"]).reset_index())

    # Rescue statistics by pair.
    rescue_stats = {}
    rescue_ps = []
    rescue_nonadd_ps = []
    for pair in PAIRS:
        z = pairs[pairs["interaction"] == pair]
        p = safe_wilcoxon(z["rescue_high_gap_closed"], z["rescue_low_gap_closed"], "greater")
        p_nonadd = safe_wilcoxon(z["rescue_high_nonadditivity"])
        rescue_ps.append(p)
        rescue_nonadd_ps.append(p_nonadd)
        rescue_stats[pair] = {
            "high_mean_gap_closed": float(z["rescue_high_gap_closed"].mean()),
            "high_median_gap_closed": float(z["rescue_high_gap_closed"].median()),
            "low_mean_gap_closed": float(z["rescue_low_gap_closed"].mean()),
            "high_vs_low_p": p,
            "fraction_high_gt_low": float(np.mean(z["rescue_high_gap_closed"] > z["rescue_low_gap_closed"])),
            "mean_increment_over_best_single": float(z["rescue_high_increment_over_best_single"].mean()),
            "mean_signed_rescue_nonadditivity": float(z["rescue_high_nonadditivity"].mean()),
            "signed_rescue_nonadditivity_two_sided_p": p_nonadd,
        }
    for pair, q, qn in zip(PAIRS, bh_adjust(rescue_ps), bh_adjust(rescue_nonadd_ps)):
        rescue_stats[pair]["high_vs_low_BH_q"] = float(q)
        rescue_stats[pair]["signed_rescue_nonadditivity_BH_q"] = float(qn)

    summary = {
        "pair_definition": "A, B and C are ranked by directly validated module effect within each locus; A is strongest.",
        "n_loci": int(pairs.Coords.nunique()),
        "n_pair_observations": int(len(pairs)),
        "epistasis_by_pair": pair_stats,
        "epistasis_pair_type_friedman": {"statistic": float(friedman.statistic), "p": float(friedman.pvalue)},
        "epistasis_pair_type_posthoc": posthoc,
        "distance": distance_stats,
        "rescue_by_pair": rescue_stats,
        "n_robust_top_decile_pairs": int(pairs.robust_top_decile_pair.sum()),
        "n_strong_rescue_pairs": int(pairs.strong_rescue_pair.sum()),
    }
    with open(OUT / "context_module_pair_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    # Figure.
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42})
    fig, axes = plt.subplots(2, 2, figsize=(10.6, 8.2), constrained_layout=True)
    axa, axb, axc, axd = axes.ravel()
    xpos = np.arange(3)
    hmean = [pair_stats[p]["high_mean_abs_epsilon"] for p in PAIRS]
    lmean = [pair_stats[p]["low_mean_abs_epsilon"] for p in PAIRS]
    axa.bar(xpos - .18, hmean, .36, color=HIGH, label="High-attribution modules")
    axa.bar(xpos + .18, lmean, .36, color=LOW, label="Low controls")
    axa.set_xticks(xpos, ["A+B", "A+C", "B+C"])
    axa.set_ylabel("Mean |epistasis ε|")
    axa.set_title("A+B, the two strongest modules, interact most")
    axa.legend(frameon=False, fontsize=8)
    for i, p in enumerate(PAIRS):
        axa.text(i, max(hmean[i], lmean[i]) + .004, f"q={pair_stats[p]['high_vs_low_BH_q']:.2g}", ha="center", fontsize=7.5)

    # high and low distance-bin curves, calculating low bins here.
    axb.errorbar(np.arange(4), distance_bins["mean"], yerr=distance_bins["sem"],
                 marker="o", linewidth=2, color=HIGH, capsize=3, label="High-attribution modules")
    low_for_plot = []
    low_pairs = pair_agg[pair_agg.module_type == "low_control"]
    for _, row in low_pairs.iterrows():
        a, b = row.interaction
        ma, mb = module_idx.loc[(row.Coords, "low_control", a)], module_idx.loc[(row.Coords, "low_control", b)]
        gap = max(0, max(ma.start_0based, mb.start_0based) - min(ma.end_0based_exclusive, mb.end_0based_exclusive))
        low_for_plot.append({"gap": gap, "effect": row.mean_absolute_epsilon})
    low_for_plot = pd.DataFrame(low_for_plot)
    low_for_plot["distance_bin"] = pd.cut(low_for_plot.gap, bins, labels=labels)
    low_bins = low_for_plot.groupby("distance_bin", observed=True).effect.agg(["mean", "sem"]).reindex(labels)
    axb.errorbar(np.arange(4), low_bins["mean"], yerr=low_bins["sem"], marker="o", linewidth=2,
                 color=LOW, capsize=3, label="Low controls")
    axb.set_xticks(np.arange(4), labels)
    axb.set_xlabel("Edge-to-edge module separation (bp)")
    axb.set_ylabel("Mean |epistasis ε|")
    axb.set_title("Non-additivity decays with module separation")
    axb.legend(frameon=False, fontsize=8)

    rh = [rescue_stats[p]["high_mean_gap_closed"] for p in PAIRS]
    rl = [rescue_stats[p]["low_mean_gap_closed"] for p in PAIRS]
    axc.bar(xpos - .18, rh, .36, color=HIGH, label="High-attribution modules")
    axc.bar(xpos + .18, rl, .36, color=LOW, label="Low controls")
    axc.set_xticks(xpos, ["A+B", "A+C", "B+C"])
    axc.set_ylabel("Mean fraction of context gap closed")
    axc.set_title("Rescue requires the strongest module A")
    for i, p in enumerate(PAIRS):
        axc.text(i, max(rh[i], rl[i]) + .025, f"q={rescue_stats[p]['high_vs_low_BH_q']:.2g}", ha="center", fontsize=7.5)

    top = pairs.head(12).sort_values("mean_absolute_epsilon")
    labels_top = [f"{r.Coords} {r.interaction}" for _, r in top.iterrows()]
    colors = ["#b91c1c" if x == "synergy" else "#2563eb" for x in top.interaction_mode]
    axd.barh(np.arange(len(top)), top.mean_absolute_epsilon, color=colors, alpha=.8)
    axd.set_yticks(np.arange(len(top)), labels_top, fontsize=6.8)
    axd.set_xlabel("Mean |epistasis ε|")
    axd.set_title("Strongest individual module pairs")
    axd.text(.98, .03, "red: synergistic loss\nblue: redundancy", transform=axd.transAxes,
             ha="right", va="bottom", fontsize=7.5)
    for label, ax in zip("ABCD", axes.ravel()):
        ax.text(-.13, 1.07, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")
    fig.suptitle("Pair-level combinatorial grammar and context rescue", fontsize=14, fontweight="bold")
    fig.savefig(OUT / "Fig_context_pair_grammar.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "Fig_context_pair_grammar.pdf", bbox_inches="tight")

    # Human-readable pair report with the top ten exact pairs.
    top10 = pairs.head(10)
    lines = [
        "# Pair-level combinatorial grammar and rescue", "",
        "## Answer", "",
        "A, B and C are the three directly validated 12-bp context modules at each locus, ranked by direct perturbation effect (A strongest). All three pair types show more absolute epistasis than matched low-attribution controls, but A+B is the strongest interaction class. Rescue is significant for A+B and A+C, but not B+C. This means rescue is largely anchored by the strongest module rather than requiring a unique B+C grammar.", "",
        "The sign of disruption epistasis is not consistent: approximately half of pairs show synergistic loss and half show redundancy. The robust feature is the magnitude of non-additivity, not a universal AND-gate.", "",
        "## Pair-type tests", "",
        "| Pair | Interpretation | High mean | Low-control mean | Ratio | BH q | Synergistic fraction |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    meanings = {"AB": "two strongest modules", "AC": "strongest + third", "BC": "second + third"}
    for p in PAIRS:
        s = pair_stats[p]
        lines.append(f"| {p} | {meanings[p]} | {s['high_mean_abs_epsilon']:.4f} | {s['low_mean_abs_epsilon']:.4f} | {s['high_low_ratio']:.2f}× | {s['high_vs_low_BH_q']:.2g} | {s['synergy_fraction']:.1%} |")
    lines += ["", f"AB, AC and BC differ overall by a repeated-measures Friedman test (P = {friedman.pvalue:.3g}); AB exceeds AC and BC after correction, while AC and BC are indistinguishable.", "",
              "## Spacing component", "",
              f"Within a locus, the nearest high-module pair had mean |ε| = {distance_stats['high']['nearest_mean_abs_epsilon']:.4f}, versus {distance_stats['high']['farthest_mean_abs_epsilon']:.4f} for the farthest pair (paired P = {distance_stats['high']['nearest_gt_farthest_p']:.2g}). The nearer pair was stronger in {distance_stats['high']['fraction_loci_nearest_gt_farthest']:.1%} of loci. Low-attribution controls also show distance decay, so proximity is a general locality effect that is amplified, rather than unique, in causal modules.", "",
              "## Pairwise rescue", "",
              "| Pair restored | High mean gap closed | Low-control mean | BH q | Interpretation |",
              "|---|---:|---:|---:|---|",
    ]
    for p in PAIRS:
        s = rescue_stats[p]
        interp = "supported" if s["high_vs_low_BH_q"] < .05 else "not supported"
        lines.append(f"| {p} | {s['high_mean_gap_closed']:.1%} | {s['low_mean_gap_closed']:.1%} | {s['high_vs_low_BH_q']:.2g} | {interp} |")
    lines += ["", "Pair rescue is generally subadditive relative to the sum of the two single-module rescues, consistent with saturation or redundancy. A pair also does not reliably outperform its best single member. Therefore the rescue result supports modular sufficiency, but not obligate pairwise synergy.", "",
              "## Strongest individual affected pairs", "",
              "These are ranked hypotheses, not per-locus significance calls; the three shuffle realizations are technical sequence variants, not biological replicates.", "",
              "| Rank | Locus | Class | Pair | Module 1 (offset; sequence) | Module 2 (offset; sequence) | Gap | Mean |ε| | Mode | Pair rescue |",
              "|---:|---|---|---|---|---|---:|---:|---|---:|",
    ]
    for _, r in top10.iterrows():
        lines.append(
            f"| {int(r.epistasis_rank)} | {r.Coords} | {r.resolved_class} | {r.interaction} | "
            f"{int(r.module_1_start)}–{int(r.module_1_end)}; `{r.module_1_sequence}` | "
            f"{int(r.module_2_start)}–{int(r.module_2_end)}; `{r.module_2_sequence}` | "
            f"{int(r.edge_to_edge_gap_bp)} bp | {r.mean_absolute_epsilon:.3f} | {r.interaction_mode} | {r.rescue_high_gap_closed:.1%} |"
        )
    lines += ["", "## Motif-pair interpretation", "",
              "These should currently be described as sequence-module pairs, not named TF–TF pairs. The earlier composition-controlled scan of 15 canonical motif families found no significant family enrichment. In particular, the data do not justify labelling the pairs AP-1–NKX2-5 or another specific TF grammar without additional motif perturbations or binding data.", ""]
    (OUT / "context_module_pair_report.md").write_text("\n".join(lines))
    print(f"Wrote {len(pairs)} locus-pair rows; top-decile robust={pairs.robust_top_decile_pair.sum()}")


if __name__ == "__main__":
    main()
