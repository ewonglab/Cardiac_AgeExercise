#!/usr/bin/env python3
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Assemble the context-seqlet manuscript figure.
# Steps: 1) Load attribution, disruption, epistasis and rescue summaries; 2) render the six-panel figure; 3) save publication outputs.
# =======================================
"""Create a consolidated manuscript-style figure for the context-module analyses."""

from pathlib import Path
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from config.paths import publication_path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
from scipy.stats import pearsonr


ROOT = Path(publication_path("MPRA_SEQLET_TASK_ROOT"))
OUT = ROOT / "outputs"
SRC = Path(publication_path("MPRA_KOO_PROJECT", "bundle"))

CLASS_ORDER = ["Stable", "Age-only", "Exercise-only", "Reversal"]
COLORS = {
    "Stable": "#6b7280",
    "Age-only": "#2563eb",
    "Exercise-only": "#059669",
    "Reversal": "#7c3aed",
}
HIGH = "#d95f02"
LOW = "#4c78a8"


def resolved_classes():
    dat = pd.read_csv(SRC / "MPRA_Results_with_TSS_distance.csv")
    sets = dat.groupby("Coords")["Group"].agg(set)

    def classify(groups):
        if groups == {"Stable"}:
            return "Stable"
        if ("Age Up" in groups and "Ex Down" in groups) or ("Age Down" in groups and "Ex Up" in groups):
            return "Reversal"
        if groups & {"Age Up", "Age Down"}:
            return "Age-only"
        if groups & {"Ex Up", "Ex Down"}:
            return "Exercise-only"
        return "Other"

    return sets.map(classify).rename("resolved_class")


def style_boxplot(bp, color_map=None):
    for i, box in enumerate(bp["boxes"]):
        if color_map:
            box.set_facecolor(color_map[i])
        box.set_alpha(0.75)
        box.set_edgecolor("#374151")
    for key in ["whiskers", "caps", "medians"]:
        for artist in bp[key]:
            artist.set_color("#374151")


def panel_label(ax, label):
    ax.text(-0.14, 1.08, label, transform=ax.transAxes, fontsize=13, fontweight="bold", va="top")


def main():
    pred = pd.read_csv(OUT / "mpra_nested_contrast_predictions.csv").rename(columns={
        "Coords": "coords",
        "measured__context_Motif_minus_Full": "measured",
        "nested_chr_oof__context_Motif_minus_Full": "predicted",
    })
    pred = pred.join(resolved_classes(), on="coords")
    attr = pd.read_csv(OUT / "context_probe_attribution_summary.csv")
    direct = pd.read_csv(OUT / "context_module_disruption_summary.csv")
    epi = pd.read_csv(OUT / "context_module_epistasis.csv")
    rescue = pd.read_csv(OUT / "context_module_rescue.csv")
    locus = pd.read_csv(OUT / "context_module_locus_level_results.csv")
    with open(OUT / "context_module_consolidated_summary.json") as fh:
        consolidated = json.load(fh)["class_tests_independent_locus_unit"]
    with open(OUT / "context_module_disruption_summary.json") as fh:
        disruption_stats = json.load(fh)
    with open(OUT / "context_module_epistasis_rescue_summary.json") as fh:
        epistasis_stats = json.load(fh)

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
    })
    fig, axes = plt.subplots(3, 2, figsize=(11.2, 12.2), constrained_layout=True)
    axa, axb, axc, axd, axe, axf = axes.ravel()

    # A: held-out context prediction.
    for cls in CLASS_ORDER:
        z = pred[pred["resolved_class"] == cls]
        axa.scatter(z["measured"], z["predicted"], s=10, alpha=0.45, linewidth=0,
                    color=COLORS[cls], label=f"{cls} (n={len(z)})")
    lims = [min(pred["measured"].min(), pred["predicted"].min()),
            max(pred["measured"].max(), pred["predicted"].max())]
    axa.plot(lims, lims, linestyle="--", linewidth=1, color="#6b7280")
    axa.set(xlabel="Measured context effect (motif shuffle − full shuffle)",
            ylabel="Chromosome-held-out prediction",
            title="Contrast-aware model predicts context dependence")
    axa.legend(frameon=False, fontsize=7, loc="upper left")
    r = pearsonr(pred["measured"], pred["predicted"]).statistic
    r2 = 1 - np.sum((pred["measured"] - pred["predicted"]) ** 2) / np.sum((pred["measured"] - pred["measured"].mean()) ** 2)
    axa.text(0.98, 0.03, f"r = {r:.3f}\nR² = {r2:.3f}",
             transform=axa.transAxes, ha="right", va="bottom")
    panel_label(axa, "A")

    # B: attribution concentration.
    vals = [attr.loc[attr["resolved_class"] == cls, "attribution_gini_motif"].dropna() for cls in CLASS_ORDER]
    bp = axb.boxplot(vals, labels=[x.replace("-only", "") for x in CLASS_ORDER], patch_artist=True,
                     widths=0.65, showfliers=False)
    style_boxplot(bp, [COLORS[c] for c in CLASS_ORDER])
    rng = np.random.default_rng(3)
    for i, (cls, y) in enumerate(zip(CLASS_ORDER, vals), 1):
        axb.scatter(rng.normal(i, 0.055, len(y)), y, s=8, alpha=0.4, color=COLORS[cls], linewidth=0)
    axb.set(ylabel="Attribution Gini coefficient", title="Context dependence is similarly distributed")
    axb.text(0.98, 0.97, "class P = 0.80", transform=axb.transAxes, ha="right", va="top")
    panel_label(axb, "B")

    # C: direct local perturbation, paired by locus.
    wide = (locus.set_index("Coords")[["direct_low_control", "direct_high"]]
            .rename(columns={"direct_low_control": "low", "direct_high": "high"}).dropna())
    for _, row in wide.iterrows():
        axc.plot([1, 2], [row["low"], row["high"]], color="#9ca3af", alpha=0.22, linewidth=0.65)
    axc.boxplot([wide["low"], wide["high"]], positions=[1, 2], widths=0.48, showfliers=False,
                patch_artist=True, boxprops=dict(facecolor="white", edgecolor="#374151"),
                medianprops=dict(color="#111827"), whiskerprops=dict(color="#374151"), capprops=dict(color="#374151"))
    axc.scatter([1, 2], [wide["low"].mean(), wide["high"].mean()], s=34, color=[LOW, HIGH], zorder=4)
    axc.set_xticks([1, 2], ["Low-attribution\ncontrol", "High-attribution\nmodule"])
    axc.set(ylabel="|Change in predicted context effect|", title="Local shuffles validate nominated modules")
    axc.text(0.03, 0.97,
             f"{disruption_stats['attribution_validation']['high_to_low_ratio']:.2f}× larger\npaired P = {disruption_stats['attribution_validation']['paired_wilcoxon_greater_p']:.1e}",
             transform=axc.transAxes, va="top")
    panel_label(axc, "C")

    # D: epistasis.
    ewide = (epi.groupby(["Coords", "module_type"])["absolute_epsilon"].mean().unstack()
             .rename(columns={"low_control": "low"}).dropna())
    for _, row in ewide.iterrows():
        axd.plot([1, 2], [row["low"], row["high"]], color="#9ca3af", alpha=0.26, linewidth=0.7)
    axd.boxplot([ewide["low"], ewide["high"]], positions=[1, 2], widths=0.48, showfliers=False,
                patch_artist=True, boxprops=dict(facecolor="white", edgecolor="#374151"),
                medianprops=dict(color="#111827"), whiskerprops=dict(color="#374151"), capprops=dict(color="#374151"))
    axd.scatter([1, 2], [ewide["low"].mean(), ewide["high"].mean()], s=34, color=[LOW, HIGH], zorder=4)
    axd.set_xticks([1, 2], ["Low-attribution\ncontrol", "High-attribution\nmodule"])
    axd.set(ylabel="Mean |non-additive interaction (ε)|", title="Nominated modules interact non-additively")
    axd.text(0.03, 0.97,
             f"{epistasis_stats['epistasis']['high_to_low_ratio']:.2f}× larger\npaired P = {epistasis_stats['epistasis']['paired_wilcoxon_greater_p']:.1e}",
             transform=axd.transAxes, va="top")
    panel_label(axd, "D")

    # E: rescue curve.
    summary = (rescue.groupby(["module_type", "n_modules"])["absolute_gap_closed"]
               .agg(["mean", "sem"]).reset_index())
    for typ, color, label in [("low_control", LOW, "Low-attribution control"), ("high", HIGH, "High-attribution module")]:
        z = summary[summary["module_type"] == typ]
        axe.errorbar(z["n_modules"], z["mean"], yerr=z["sem"], marker="o", markersize=5,
                     capsize=3, linewidth=2, color=color, label=label)
    axe.axhline(0, color="#9ca3af", linewidth=0.8)
    axe.set_xticks([1, 2, 3])
    axe.set(xlabel="Number of native modules restored into full shuffle",
            ylabel="Fraction of predicted WT–full-shuffle gap closed",
            title="Three short modules restore much of the context effect")
    axe.legend(frameon=False, fontsize=8)
    axe.text(0.98, 0.05,
             f"3 modules: {epistasis_stats['rescue']['high_triple_gap_closed_median']:.0%} median rescue\npaired P = {epistasis_stats['rescue']['paired_wilcoxon_greater_p']:.1e}",
             transform=axe.transAxes, ha="right", va="bottom")
    panel_label(axe, "E")

    # F: class means, normalized to the overall mean for each assay.
    outcomes = [
        ("Direct module\ndisruption", "direct_high", consolidated["direct_high_effect"]),
        ("Non-additive\ninteraction", "absolute_epsilon_high", consolidated["high_abs_epistasis"]),
        ("Three-module\nrescue", "rescue_high", consolidated["high_triple_rescue"]),
    ]
    mat = np.full((len(outcomes), len(CLASS_ORDER)), np.nan)
    for i, (_, col, _) in enumerate(outcomes):
        overall = locus[col].dropna().mean()
        for j, cls in enumerate(CLASS_ORDER):
            mat[i, j] = locus.loc[locus["resolved_class"] == cls, col].dropna().mean() / overall
    im = axf.imshow(mat, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=0.70, vcenter=1, vmax=1.30), aspect="auto")
    axf.set_xticks(np.arange(len(CLASS_ORDER)), [x.replace("-only", "") for x in CLASS_ORDER])
    axf.set_yticks(np.arange(len(outcomes)), [x[0] for x in outcomes])
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            axf.text(j, i, f"{mat[i, j]:.2f}×", ha="center", va="center", fontsize=8,
                     color="white" if abs(mat[i, j] - 1) > 0.18 else "#111827")
        stat = outcomes[i][2]
        axf.text(3.65, i + 0.30,
                 f"class P={stat['kruskal_p_locus_level']:.2f}\nadj. P={stat['class_partial_F_controlling_abs_effect_GC_Endo']['p']:.2f}",
                 ha="left", va="center", fontsize=7.5, clip_on=False)
    axf.set_title("Module architecture does not separate biological classes")
    axf.set_xlabel("Class mean / overall mean")
    axf.spines[["left", "bottom"]].set_visible(False)
    cbar = fig.colorbar(im, ax=axf, fraction=0.035, pad=0.16)
    cbar.set_label("Relative class mean")
    panel_label(axf, "F")

    fig.suptitle("Perturbation-aware AlphaGenome reveals a shared, distributed context code",
                 fontsize=15, fontweight="bold")
    fig.savefig(OUT / "Fig_context_modules_manuscript.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "Fig_context_modules_manuscript.pdf", bbox_inches="tight")
    print(OUT / "Fig_context_modules_manuscript.png")


if __name__ == "__main__":
    main()
