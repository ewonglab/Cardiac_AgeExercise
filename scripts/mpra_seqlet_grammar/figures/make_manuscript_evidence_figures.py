#!/usr/bin/env python3
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Assemble MPRA sequence-context evidence figures.
# Steps: 1) Load final analysis tables; 2) draw evidence panels; 3) save manuscript-ready figure files.
# =======================================
"""Create manuscript-ready evidence figures for the MPRA grammar study."""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from config.paths import publication_path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(publication_path("MPRA_SEQLET_TASK_ROOT"))
OUT = ROOT / "outputs"
PROJECT = Path(publication_path("MPRA_KOO_PROJECT"))
RAW_RESULTS = PROJECT / "bundle/MPRA_Results_with_TSS_distance.csv"

# Restrained manuscript palette. Blue is used for dynamic/high/supported evidence;
# rust is used only to distinguish the second factor/readout.
INK = "#30343B"
MID = "#7C828A"
LIGHT = "#D7D9DC"
PALE = "#EEF0F1"
BLUE = "#58788D"
BLUE_LIGHT = "#B7C7D0"
RUST = "#9A6B57"
RUST_LIGHT = "#D2B8AC"
OLIVE = "#7B8265"
WHITE = "#FFFFFF"

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "axes.titlesize": 9.5,
        "axes.labelsize": 8.5,
        "axes.edgecolor": MID,
        "axes.linewidth": 0.75,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.color": INK,
        "ytick.color": INK,
        "text.color": INK,
        "axes.labelcolor": INK,
        "figure.facecolor": WHITE,
        "savefig.facecolor": WHITE,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor=WHITE)
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", facecolor=WHITE)
    plt.close(fig)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.14,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
        ha="left",
        clip_on=False,
    )


def resolved_class(groups: set[str]) -> str:
    if groups == {"Stable"}:
        return "Stable"
    if ("Age Up" in groups and "Ex Down" in groups) or (
        "Age Down" in groups and "Ex Up" in groups
    ):
        return "Reversal"
    if groups & {"Age Up", "Age Down"}:
        return "Age-only"
    if groups & {"Ex Up", "Ex Down"}:
        return "Exercise-only"
    return "Other"


def format_p(value: float) -> str:
    if value < 0.001:
        exponent = int(math.floor(math.log10(value)))
        coeff = value / (10**exponent)
        return rf"$P={coeff:.1f}\times10^{{{exponent}}}$"
    return rf"$P={value:.3f}$"


def style_axis(ax: plt.Axes, grid: str | None = "y") -> None:
    if grid:
        ax.grid(axis=grid, color=LIGHT, linewidth=0.55, alpha=0.75)
        ax.set_axisbelow(True)


def sequence_bar(
    ax: plt.Axes,
    y: float,
    mode: str,
    label: str,
    ap1_start: float = 0.40,
    ap1_width: float = 0.17,
) -> None:
    x0, width, height = 0.17, 0.73, 0.30
    if mode == "full":
        rect = patches.Rectangle(
            (x0, y),
            width,
            height,
            facecolor=PALE,
            edgecolor=MID,
            linewidth=0.9,
            hatch="////",
        )
        ax.add_patch(rect)
    else:
        ax.add_patch(
            patches.Rectangle(
                (x0, y),
                width,
                height,
                facecolor="#E3E4E6",
                edgecolor=MID,
                linewidth=0.9,
            )
        )
        for x in np.linspace(x0 + 0.03, x0 + width - 0.03, 16):
            ax.plot([x, x + 0.012], [y + 0.08, y + 0.22], color="#B7BABE", lw=0.7)
        ap1_x = x0 + width * ap1_start
        if mode == "native":
            ax.add_patch(
                patches.Rectangle(
                    (ap1_x, y),
                    width * ap1_width,
                    height,
                    facecolor=BLUE,
                    edgecolor=INK,
                    linewidth=0.7,
                )
            )
            ax.text(
                ap1_x + width * ap1_width / 2,
                y + height / 2,
                "AP-1",
                color=WHITE,
                ha="center",
                va="center",
                fontsize=7,
                fontweight="bold",
            )
        else:
            ax.add_patch(
                patches.Rectangle(
                    (ap1_x, y),
                    width * ap1_width,
                    height,
                    facecolor=PALE,
                    edgecolor=INK,
                    linewidth=0.7,
                    hatch="xxxx",
                )
            )
    ax.text(0.02, y + height / 2, label, ha="left", va="center", fontsize=8.5)


def figure_1() -> None:
    pred = pd.read_csv(OUT / "mpra_nested_contrast_predictions.csv")
    pred["state"] = np.where(pred["Group"].eq("Stable"), "Stable", "Dynamic")
    contrasts = [
        (
            "total",
            "measured__total_Endo_minus_Full",
            "nested_chr_oof__total_Endo_minus_Full",
            "Native − full",
            0.772,
            0.596,
        ),
        (
            "ap1",
            "measured__AP1_Endo_minus_Motif",
            "nested_chr_oof__AP1_Endo_minus_Motif",
            "Native − AP-1 shuffle",
            0.431,
            0.183,
        ),
        (
            "context",
            "measured__context_Motif_minus_Full",
            "nested_chr_oof__context_Motif_minus_Full",
            "AP-1 shuffle − full",
            0.754,
            0.568,
        ),
    ]
    raw = pd.read_csv(RAW_RESULTS)
    classes = raw.groupby("Coords")["Group"].agg(set).map(resolved_class)
    gc = raw.drop_duplicates("Coords").set_index("Coords")[["GC_content"]].copy()
    gc["class"] = classes
    gc["state"] = np.where(gc["class"].eq("Stable"), "Stable", "Dynamic")

    fig = plt.figure(figsize=(11.2, 7.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, width_ratios=[1.28, 1, 1], wspace=0.25, hspace=0.30)

    ax = fig.add_subplot(gs[0, 0])
    panel_label(ax, "A")
    sequence_bar(ax, 2.48, "native", "Native")
    sequence_bar(ax, 1.68, "motif", "AP-1 shuffle")
    sequence_bar(ax, 0.88, "full", "Full shuffle")
    ax.annotate(
        "",
        xy=(0.94, 2.63),
        xytext=(0.94, 1.83),
        arrowprops=dict(arrowstyle="<->", color=BLUE, lw=1.2),
    )
    ax.text(0.97, 2.23, "AP-1\ncomponent", va="center", fontsize=7.7, color=BLUE)
    ax.annotate(
        "",
        xy=(0.94, 1.83),
        xytext=(0.94, 1.03),
        arrowprops=dict(arrowstyle="<->", color=RUST, lw=1.2),
    )
    ax.text(0.97, 1.43, "remaining\ncontext", va="center", fontsize=7.7, color=RUST)
    ax.text(
        0.17,
        0.34,
        "Subtraction partitions total native sequence dependence\n"
        "without regressing away genomic composition.",
        fontsize=8,
        va="top",
    )
    ax.set_xlim(0, 1.18)
    ax.set_ylim(0.12, 3.12)
    ax.axis("off")
    ax.set_title("Experimental contrasts isolate two sequence components", loc="left", pad=4)

    ax = fig.add_subplot(gs[1, 0])
    panel_label(ax, "B")
    positions = np.arange(3)
    pvals = [0.026, 0.425, 0.043]
    for i, (_, measured, _, label, _, _) in enumerate(contrasts):
        stable = pred.loc[pred.state.eq("Stable"), measured].dropna().to_numpy()
        dynamic = pred.loc[pred.state.eq("Dynamic"), measured].dropna().to_numpy()
        bp = ax.boxplot(
            [stable, dynamic],
            positions=[i - 0.17, i + 0.17],
            widths=0.25,
            patch_artist=True,
            showfliers=False,
            medianprops=dict(color=INK, lw=1.0),
            whiskerprops=dict(color=MID, lw=0.8),
            capprops=dict(color=MID, lw=0.8),
        )
        bp["boxes"][0].set(facecolor=LIGHT, edgecolor=MID)
        bp["boxes"][1].set(facecolor=BLUE_LIGHT, edgecolor=BLUE)
        ax.scatter(i - 0.17, stable.mean(), marker="D", s=22, color=INK, zorder=3)
        ax.scatter(i + 0.17, dynamic.mean(), marker="D", s=22, color=BLUE, zorder=3)
        ymax = max(np.quantile(stable, 0.95), np.quantile(dynamic, 0.95))
        ax.text(i, ymax + 0.10, format_p(pvals[i]), ha="center", va="bottom", fontsize=7.2)
    ax.axhline(0, color=INK, lw=0.7)
    ax.set_xticks(positions, ["Native −\nfull", "Native −\nAP-1 shuf.", "AP-1 shuf. −\nfull"])
    ax.set_ylabel("Measured MPRA contrast (log₂)")
    ax.set_title("Dynamic enhancers depend more on non-AP-1 context", loc="left")
    style_axis(ax)
    ax.legend(
        [
            patches.Patch(facecolor=LIGHT, edgecolor=MID),
            patches.Patch(facecolor=BLUE_LIGHT, edgecolor=BLUE),
        ],
        ["Stable", "Dynamic"],
        frameon=False,
        ncol=2,
        loc="lower right",
        fontsize=7.5,
    )

    for j, (_, measured, predicted, label, reported_r, reported_r2) in enumerate(contrasts):
        row, col = (0, 1) if j == 0 else ((0, 2) if j == 1 else (1, 1))
        ax = fig.add_subplot(gs[row, col])
        panel_label(ax, chr(ord("C") + j))
        z = pred[[measured, predicted]].dropna()
        ax.scatter(
            z[measured],
            z[predicted],
            s=7,
            alpha=0.26,
            color=BLUE if j != 1 else MID,
            linewidths=0,
            rasterized=True,
        )
        lo = float(min(z.min()))
        hi = float(max(z.max()))
        pad = (hi - lo) * 0.04
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color=LIGHT, lw=1.0)
        slope, intercept = np.polyfit(z[measured], z[predicted], 1)
        xx = np.array([lo, hi])
        ax.plot(xx, intercept + slope * xx, color=INK, lw=1.1)
        ax.text(
            0.04,
            0.96,
            rf"$r={reported_r:.3f}$" + "\n" + rf"$R^2={reported_r2:.3f}$",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.25", fc=WHITE, ec=LIGHT, alpha=0.92),
        )
        ax.set_xlabel("Measured contrast")
        ax.set_ylabel("Chromosome-held-out prediction")
        ax.set_title(label, loc="left")
        style_axis(ax)

    ax = fig.add_subplot(gs[1, 2])
    panel_label(ax, "F")
    stable_gc = gc.loc[gc.state.eq("Stable"), "GC_content"].dropna().to_numpy()
    dynamic_gc = gc.loc[gc.state.eq("Dynamic"), "GC_content"].dropna().to_numpy()
    parts = ax.violinplot(
        [stable_gc, dynamic_gc],
        positions=[0, 1],
        widths=0.75,
        showextrema=False,
    )
    for body, color in zip(parts["bodies"], [LIGHT, BLUE_LIGHT]):
        body.set_facecolor(color)
        body.set_edgecolor(MID if color == LIGHT else BLUE)
        body.set_alpha(0.95)
    bp = ax.boxplot(
        [stable_gc, dynamic_gc],
        positions=[0, 1],
        widths=0.20,
        showfliers=False,
        patch_artist=True,
        boxprops=dict(facecolor=WHITE, edgecolor=INK, lw=0.8),
        medianprops=dict(color=INK, lw=1.1),
        whiskerprops=dict(color=INK, lw=0.75),
        capprops=dict(color=INK, lw=0.75),
    )
    ax.scatter([0, 1], [stable_gc.mean(), dynamic_gc.mean()], marker="D", s=25, color=[INK, BLUE])
    ax.text(
        0.5,
        max(np.quantile(stable_gc, 0.99), np.quantile(dynamic_gc, 0.99)) + 0.012,
        r"$P=1.31\times10^{-6}$",
        ha="center",
        fontsize=7.5,
    )
    ax.set_xticks([0, 1], ["Stable\n0.517", "Dynamic\n0.497"])
    ax.set_ylabel("Native GC fraction")
    ax.set_title("Lower GC is part of the dynamic state", loc="left")
    style_axis(ax)

    fig.suptitle(
        "Sequence context, not AP-1 loss alone, distinguishes dynamic enhancers",
        x=0.01,
        y=1.01,
        ha="left",
        fontsize=13,
        fontweight="bold",
    )
    save(fig, "Fig1_sequence_context_and_prediction")


def comparison_scatter(
    ax: plt.Axes,
    frame: pd.DataFrame,
    xcol: str,
    ycol: str,
    xlabel: str,
    ylabel: str,
    title: str,
    annotation: str,
) -> None:
    z = frame[[xcol, ycol]].dropna()
    ax.scatter(
        z[xcol],
        z[ycol],
        s=18,
        color=BLUE,
        alpha=0.55,
        edgecolor=WHITE,
        linewidth=0.25,
        rasterized=True,
    )
    low = min(z[xcol].quantile(0.01), z[ycol].quantile(0.01))
    high = max(z[xcol].quantile(0.99), z[ycol].quantile(0.99))
    pad = max((high - low) * 0.08, 0.02)
    ax.plot([low - pad, high + pad], [low - pad, high + pad], color=MID, ls="--", lw=0.9)
    ax.set_xlim(low - pad, high + pad)
    ax.set_ylim(low - pad, high + pad)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left")
    ax.text(
        0.04,
        0.96,
        annotation,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.6,
        bbox=dict(boxstyle="round,pad=0.25", fc=WHITE, ec=LIGHT, alpha=0.94),
    )
    style_axis(ax)


def figure_2() -> None:
    locus = pd.read_csv(OUT / "context_module_locus_level_results.csv")
    rescue = pd.read_csv(OUT / "context_module_rescue.csv")
    rescue_locus = (
        rescue.groupby(["Coords", "module_type", "n_modules"], as_index=False)[
            "absolute_gap_closed"
        ]
        .mean()
        .replace([np.inf, -np.inf], np.nan)
        .dropna()
    )
    summary = (
        rescue_locus.groupby(["module_type", "n_modules"])["absolute_gap_closed"]
        .agg(["mean", "sem"])
        .reset_index()
    )

    fig = plt.figure(figsize=(12.0, 7.6), constrained_layout=False)
    gs = fig.add_gridspec(2, 4, height_ratios=[0.82, 1.18], hspace=0.29, wspace=0.26)
    fig.subplots_adjust(left=0.07, right=0.985, top=0.91, bottom=0.10, wspace=0.48, hspace=0.50)

    ax = fig.add_subplot(gs[0, :])
    panel_label(ax, "A")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Attribution-guided perturbation asks whether context is modular and rescuable", loc="left")
    x0, y0, width, height = 0.05, 0.57, 0.58, 0.16
    ax.add_patch(patches.Rectangle((x0, y0), width, height, fc=PALE, ec=MID, lw=1))
    for x in np.linspace(x0 + 0.02, x0 + width - 0.02, 24):
        ax.plot([x, x + 0.007], [y0 + 0.04, y0 + 0.12], color="#B4B7BB", lw=0.7)
    ap1 = patches.Rectangle((x0 + 0.24, y0), 0.08, height, fc=INK, ec=INK)
    ax.add_patch(ap1)
    ax.text(x0 + 0.28, y0 + height / 2, "AP-1", color=WHITE, ha="center", va="center", fontsize=7)
    highs = [(0.11, "A"), (0.39, "B"), (0.50, "C")]
    lows = [(0.18, "L1"), (0.34, "L2"), (0.54, "L3")]
    for relx, label in highs:
        ax.add_patch(
            patches.Rectangle(
                (x0 + relx * width, y0 - 0.03),
                0.043,
                height + 0.06,
                fc=BLUE_LIGHT,
                ec=BLUE,
                lw=1.2,
            )
        )
        ax.text(
            x0 + relx * width + 0.0215,
            y0 + height + 0.055,
            label,
            color=BLUE,
            ha="center",
            va="bottom",
            fontsize=7.5,
            fontweight="bold",
        )
    for relx, label in lows:
        ax.add_patch(
            patches.Rectangle(
                (x0 + relx * width, y0 + 0.01),
                0.034,
                height - 0.02,
                fc=WHITE,
                ec=MID,
                lw=0.9,
                hatch="..",
            )
        )
        ax.text(
            x0 + relx * width + 0.017,
            y0 - 0.05,
            label,
            color=MID,
            ha="center",
            va="top",
            fontsize=6.8,
        )
    ax.text(x0, y0 + height + 0.19, "170-bp enhancer", fontsize=8.5, fontweight="bold")
    ax.text(x0, y0 - 0.16, "High-attribution 12-bp modules", color=BLUE, fontsize=8)
    ax.text(x0 + 0.28, y0 - 0.16, "Matched low-attribution controls", color=MID, fontsize=8)
    ax.annotate(
        "",
        xy=(0.70, 0.66),
        xytext=(0.64, 0.66),
        arrowprops=dict(arrowstyle="->", color=INK, lw=1.1),
    )
    boxes = [
        (0.72, 0.68, "single-module\nloss", "direct effect"),
        (0.82, 0.68, "pair loss\nA+B", "epistasis"),
        (0.92, 0.68, "restore\nA→A+B→ABC", "rescue"),
    ]
    for x, y, headline, sub in boxes:
        ax.add_patch(
            patches.FancyBboxPatch(
                (x - 0.045, y - 0.12),
                0.09,
                0.24,
                boxstyle="round,pad=0.012",
                fc=WHITE,
                ec=BLUE if "restore" in headline else MID,
                lw=1.0,
            )
        )
        ax.text(x, y + 0.025, headline, ha="center", va="center", fontsize=7)
        ax.text(x, y - 0.075, sub, ha="center", va="center", fontsize=6.6, color=MID)
    ax.text(
        0.05,
        0.16,
        "Interpretation: a distributed enhancer grammar predicts that high-attribution modules "
        "show larger loss, stronger non-additivity, and dose-dependent restoration.",
        fontsize=8.4,
        va="center",
    )

    ax = fig.add_subplot(gs[1, 0])
    panel_label(ax, "B")
    comparison_scatter(
        ax,
        locus,
        "direct_low_control",
        "direct_high",
        "Low-module effect",
        "High-module effect",
        "Direct perturbation",
        "2.58× larger\n" + r"$P=2.97\times10^{-21}$",
    )

    ax = fig.add_subplot(gs[1, 1])
    panel_label(ax, "C")
    comparison_scatter(
        ax,
        locus,
        "absolute_epsilon_low_control",
        "absolute_epsilon_high",
        "Low-pair |ε|",
        "High-pair |ε|",
        "Pairwise epistasis",
        "2.10× larger\n" + r"$P=4.26\times10^{-10}$",
    )

    ax = fig.add_subplot(gs[1, 2])
    panel_label(ax, "D")
    for module_type, color, label in [
        ("low_control", MID, "Low controls"),
        ("high", BLUE, "High modules"),
    ]:
        z = summary[summary.module_type.eq(module_type)].sort_values("n_modules")
        ax.errorbar(
            z.n_modules,
            100 * z["mean"],
            yerr=100 * 1.96 * z["sem"],
            marker="o",
            ms=5,
            lw=1.5,
            color=color,
            label=label,
            capsize=2,
        )
    high_vals = [28.4, 41.2, 54.2]
    for x, y in zip([1, 2, 3], high_vals):
        ax.text(x, y + 5.5, f"{y:.1f}%", color=BLUE, ha="center", fontsize=7)
    ax.set_xticks([1, 2, 3])
    ax.set_xlabel("Modules restored")
    ax.set_ylabel("Context gap closed (%)")
    ax.set_title("Dose-dependent rescue", loc="left")
    ax.legend(frameon=False, fontsize=7.2, loc="lower right")
    style_axis(ax)

    ax = fig.add_subplot(gs[1, 3])
    panel_label(ax, "E")
    order = ["Stable", "Age-only", "Exercise-only", "Reversal"]
    values = [locus.loc[locus.resolved_class.eq(c), "rescue_high"].dropna().to_numpy() for c in order]
    bp = ax.boxplot(
        values,
        positions=np.arange(len(order)),
        widths=0.58,
        patch_artist=True,
        showfliers=False,
        boxprops=dict(facecolor=BLUE_LIGHT, edgecolor=BLUE, lw=0.9),
        medianprops=dict(color=INK, lw=1.1),
        whiskerprops=dict(color=MID, lw=0.8),
        capprops=dict(color=MID, lw=0.8),
    )
    rng = np.random.default_rng(12)
    for i, vals in enumerate(values):
        ax.scatter(
            rng.normal(i, 0.055, size=len(vals)),
            vals,
            s=7,
            alpha=0.25,
            color=BLUE,
            linewidths=0,
            rasterized=True,
        )
    ax.set_xticks(np.arange(len(order)), ["Stable", "Age", "Exercise", "Reversal"], rotation=28, ha="right")
    ax.set_ylabel("Three-module rescue")
    ax.set_title("Shared across response classes", loc="left")
    ax.text(
        0.04,
        0.96,
        "Class test\n" + r"$P=0.735$",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.7,
        bbox=dict(boxstyle="round,pad=0.25", fc=WHITE, ec=LIGHT),
    )
    style_axis(ax)

    fig.suptitle(
        "Dynamic enhancer function is distributed across interacting sequence modules",
        x=0.01,
        y=0.98,
        ha="left",
        fontsize=13,
        fontweight="bold",
    )
    save(fig, "Fig2_distributed_modules_and_rescue")


def spacing_panel(
    ax: plt.Axes,
    spacing: pd.DataFrame,
    coord: str,
    title: str,
    peak_note: str,
    show_legend: bool = False,
    show_native: bool = True,
) -> None:
    z = spacing[(spacing.Coords.eq(coord)) & (spacing.design_mode.eq("motif_core"))].copy()
    z = z.sort_values("target_gap_bp")
    ax.plot(
        z.target_gap_bp,
        z.loss_synergy_context_Motif_minus_Full,
        marker="o",
        ms=3.2,
        lw=1.25,
        color=BLUE,
        label="Context model",
    )
    ax.plot(
        z.target_gap_bp,
        z.loss_synergy_mpra_activity,
        marker="s",
        ms=3.0,
        lw=1.1,
        color=RUST,
        label="MPRA activity head",
    )
    ax.axhline(0, color=LIGHT, lw=0.8)
    native = int(z.native_motif_gap_bp.iloc[0])
    if show_native and native <= z.target_gap_bp.max():
        ax.axvline(native, color=INK, lw=0.85, ls=":")
        ax.text(
            native,
            ax.get_ylim()[1] if ax.get_ylim()[1] else 0,
            f" native {native}",
            rotation=90,
            va="top",
            ha="left",
            fontsize=6.5,
            color=INK,
        )
    ax.text(
        0.98,
        0.96,
        peak_note,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=7.2,
        bbox=dict(boxstyle="round,pad=0.22", fc=WHITE, ec=LIGHT, alpha=0.92),
    )
    ax.set_xlabel("Edge-to-edge motif gap (bp)")
    ax.set_ylabel("Loss synergy")
    ax.set_title(title, loc="left")
    if show_legend:
        ax.legend(frameon=False, fontsize=7, loc="lower right")
    style_axis(ax)


def figure_3() -> None:
    motif = pd.read_csv(OUT / "motif_enrichment_dynamic_vs_stable_uncorrected.csv")
    pair = pd.read_csv(OUT / "context_module_pair_results.csv")
    spacing = pd.read_csv(OUT / "spacing_epistasis_by_gap.csv")

    fig = plt.figure(figsize=(11.3, 8.1), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, height_ratios=[0.94, 1.15], hspace=0.31, wspace=0.28)

    ax = fig.add_subplot(gs[0, 0])
    panel_label(ax, "A")
    names = ["TCF7L2", "GATA4", "MEF2C", "NFATC3", "ETV4"]
    z = motif.set_index("motif").loc[names].reset_index()
    y = np.arange(len(z))[::-1]
    colors = [BLUE if value > 0 else RUST for value in z.dynamic_minus_stable_hits]
    ax.barh(y, z.dynamic_minus_stable_hits, color=colors, edgecolor=WHITE, height=0.62)
    ax.axvline(0, color=INK, lw=0.8)
    ax.set_yticks(y, z.motif)
    ax.set_xlabel("Dynamic − stable motif hits / enhancer")
    ax.set_title("Uncorrected native-sequence enrichment", loc="left")
    q_map = dict(zip(z.motif, z.BH_q))
    for yi, row in zip(y, z.itertuples()):
        small_bar = abs(row.dynamic_minus_stable_hits) < 0.30
        text_x = (
            row.dynamic_minus_stable_hits - 0.025
            if small_bar and row.dynamic_minus_stable_hits < 0
            else row.dynamic_minus_stable_hits / 2
        )
        ax.text(
            text_x,
            yi,
            f"q={q_map[row.motif]:.3g}",
            va="center",
            ha="right" if small_bar and row.dynamic_minus_stable_hits < 0 else "center",
            fontsize=6.8,
            color=INK if small_bar else WHITE,
        )
    style_axis(ax, "x")

    ax = fig.add_subplot(gs[0, 1])
    panel_label(ax, "B")
    pair_order = ["AB", "AC", "BC"]
    x = np.arange(3)
    high = [pair.loc[pair.interaction.eq(k), "mean_absolute_epsilon"].mean() for k in pair_order]
    low = [pair.loc[pair.interaction.eq(k), "low_mean_absolute_epsilon"].mean() for k in pair_order]
    ax.bar(x - 0.17, low, width=0.32, color=LIGHT, edgecolor=MID, label="Low controls")
    ax.bar(x + 0.17, high, width=0.32, color=BLUE_LIGHT, edgecolor=BLUE, label="High modules")
    ratios = np.array(high) / np.array(low)
    for i, ratio in enumerate(ratios):
        ax.text(i, max(high[i], low[i]) + 0.004, f"{ratio:.2f}×", ha="center", fontsize=7)
    ax.set_xticks(x, pair_order)
    ax.set_xlabel("Ranked high-module pair")
    ax.set_ylabel("Mean |ε|")
    ax.set_title("All three pairs are non-additive", loc="left")
    ax.legend(frameon=False, fontsize=7)
    style_axis(ax)

    ax = fig.add_subplot(gs[0, 2])
    panel_label(ax, "C")
    high_r = [pair.loc[pair.interaction.eq(k), "rescue_high_gap_closed"].mean() for k in pair_order]
    low_r = [pair.loc[pair.interaction.eq(k), "rescue_low_gap_closed"].mean() for k in pair_order]
    ax.bar(x - 0.17, 100 * np.array(low_r), width=0.32, color=LIGHT, edgecolor=MID)
    ax.bar(x + 0.17, 100 * np.array(high_r), width=0.32, color=BLUE_LIGHT, edgecolor=BLUE)
    qvals = [3.24e-6, 0.00174, 0.38]
    for i, q in enumerate(qvals):
        ax.text(
            i,
            max(100 * high_r[i], 100 * low_r[i]) + 3,
            "q<0.01" if q < 0.01 else "n.s.",
            ha="center",
            fontsize=7,
        )
    ax.set_xticks(x, pair_order)
    ax.set_xlabel("Ranked high-module pair")
    ax.set_ylabel("Pair rescue (% gap closed)")
    ax.set_title("AB and AC give the clearest rescue", loc="left")
    style_axis(ax)

    ax = fig.add_subplot(gs[1, 0])
    panel_label(ax, "D")
    spacing_panel(
        ax,
        spacing,
        "chr12:81788643-81788812",
        "MEF2C–TBX5",
        "both readouts peak at 9 bp\ncurve agreement r=0.962",
        show_legend=True,
    )

    ax = fig.add_subplot(gs[1, 1])
    panel_label(ax, "E")
    spacing_panel(
        ax,
        spacing,
        "chr3:103612737-103612906",
        "TCF7L2–RUNX3, locus 1",
        "both readouts peak at 4 bp\ncurve agreement r=0.913",
    )

    ax = fig.add_subplot(gs[1, 2])
    panel_label(ax, "F")
    spacing_panel(
        ax,
        spacing,
        "chr12:15566211-15566380",
        "TCF7L2–RUNX3, locus 2",
        "context peak at 4 bp;\nactivity is less concordant",
        show_native=False,
    )
    ax.set_xlim(-1, 30)
    ax.text(
        0.98,
        0.04,
        "native gap = 87 bp (off-axis)",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.8,
        color=MID,
    )

    fig.suptitle(
        "Candidate motif pairs exhibit position-dependent combinatorial grammar",
        x=0.01,
        y=1.01,
        ha="left",
        fontsize=13,
        fontweight="bold",
    )
    save(fig, "Fig3_motifs_pairs_and_spacing")


def load_atom_table(cif_path: Path) -> pd.DataFrame:
    spec = importlib.util.spec_from_file_location(
        "af3_geometry", ROOT / "work/analyze_af3_geometry.py"
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module.atom_table(cif_path)


def project_structure(atoms: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    selection = atoms[
        ((atoms.chain.isin(["A", "B", "C"])) & atoms.atom.eq("CA"))
        | ((atoms.chain.isin(["D", "E"])) & atoms.atom.eq("P"))
    ].copy()
    coords = selection[["x", "y", "z"]].to_numpy(float)
    center = coords.mean(axis=0)
    _, _, vt = np.linalg.svd(coords - center, full_matrices=False)
    basis = vt[:2].T
    selection[["u", "v"]] = (coords - center) @ basis
    return selection, center, basis


def figure_4() -> None:
    aggregate = pd.read_csv(
        OUT
        / "alphafold3_spacing_results/analysis/alphafold3_combined_job_metrics.csv"
    )
    cif = (
        OUT
        / "alphafold3_spacing_results/jobs/"
        "tcf7l2_runx3_locus1_gap04_intact/"
        "fold_tcf7l2_runx3_locus1_gap04_intact_model_0.cif"
    )
    atoms = load_atom_table(cif)
    projection, center, basis = project_structure(atoms)

    fig = plt.figure(figsize=(12.0, 8.3), constrained_layout=False)
    gs = fig.add_gridspec(
        3,
        3,
        height_ratios=[1.06, 1.0, 0.35],
        width_ratios=[1.0, 1.15, 1.0],
        hspace=0.26,
        wspace=0.27,
    )
    fig.subplots_adjust(left=0.065, right=0.985, top=0.91, bottom=0.08, wspace=0.40, hspace=0.52)

    ax = fig.add_subplot(gs[0, 0])
    panel_label(ax, "A")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("AlphaFold 3 was used as an orthogonal plausibility filter", loc="left")
    steps = [
        (0.12, 0.71, "In-silico\nspacing peak", "4 bp"),
        (0.50, 0.71, "AF3 jobs", "19 designs"),
        (0.88, 0.71, "Replicate\nstructures", "95 models"),
    ]
    for i, (x, y, label, value) in enumerate(steps):
        ax.add_patch(
            patches.FancyBboxPatch(
                (x - 0.105, y - 0.14),
                0.21,
                0.28,
                boxstyle="round,pad=0.014",
                fc=PALE if i != 1 else BLUE_LIGHT,
                ec=BLUE if i != 2 else MID,
                lw=1.1,
            )
        )
        ax.text(x, y + 0.035, label, ha="center", va="center", fontsize=7.7)
        ax.text(x, y - 0.078, value, ha="center", va="center", fontsize=9, fontweight="bold")
        if i < len(steps) - 1:
            ax.annotate(
                "",
                xy=(steps[i + 1][0] - 0.12, y),
                xytext=(x + 0.12, y),
                arrowprops=dict(arrowstyle="->", color=INK, lw=1.0),
            )
    ax.text(
        0.04,
        0.30,
        "Fixed protein domains + designed dsDNA\n"
        "• five models per job\n"
        "• pair confidence, PAE and geometry\n"
        "• motif engagement checked explicitly",
        va="top",
        fontsize=8,
        linespacing=1.45,
    )

    ax = fig.add_subplot(gs[0, 1:])
    panel_label(ax, "B")
    chain_style = {
        "A": (BLUE, 2.0, "TCF7L2"),
        "B": (RUST, 2.0, "RUNX3"),
        "C": (OLIVE, 1.5, "CBFβ"),
        "D": (MID, 1.25, "DNA"),
        "E": (MID, 1.25, None),
    }
    for chain, (color, lw, label) in chain_style.items():
        z = projection[projection.chain.eq(chain)].sort_values("resid")
        ax.plot(z.u, z.v, color=color, lw=lw, alpha=0.95, label=label)
        if chain in ["D", "E"]:
            ax.scatter(z.u, z.v, s=4, color=color, alpha=0.7)
    contacts = [
        ("A", 32, "B", 40, "K352–T89"),
        ("A", 32, "B", 89, "K352–D138"),
    ]
    for c1, r1, c2, r2, label in contacts:
        a = atoms[(atoms.chain.eq(c1)) & (atoms.resid.eq(r1)) & (atoms.atom.eq("CA"))]
        b = atoms[(atoms.chain.eq(c2)) & (atoms.resid.eq(r2)) & (atoms.atom.eq("CA"))]
        if not a.empty and not b.empty:
            pa = (a[["x", "y", "z"]].to_numpy()[0] - center) @ basis
            pb = (b[["x", "y", "z"]].to_numpy()[0] - center) @ basis
            ax.plot([pa[0], pb[0]], [pa[1], pb[1]], color=INK, ls="--", lw=0.9)
            mid = (pa + pb) / 2
            ax.text(mid[0], mid[1], label, fontsize=6.5, ha="center", va="bottom")
            ax.scatter([pa[0], pb[0]], [pa[1], pb[1]], s=18, color=INK, zorder=5)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Representative chr3 4-bp model: both TFs engage DNA and contact", loc="left")
    ax.legend(
        frameon=False,
        ncol=4,
        fontsize=7.3,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.03),
    )

    locus1 = aggregate[
        aggregate.job_name.str.match(r"tcf7l2_runx3_locus1_gap\d+_intact")
    ].sort_values("gap_bp")
    ax = fig.add_subplot(gs[1, 0])
    panel_label(ax, "C")
    x = np.arange(len(locus1))
    med = locus1.median_primary_tf_pair_iptm.to_numpy()
    lower = med - locus1.q25_primary_tf_pair_iptm.to_numpy()
    upper = locus1.q75_primary_tf_pair_iptm.to_numpy() - med
    colors = [BLUE if gap == 4 else LIGHT for gap in locus1.gap_bp]
    ax.bar(x, med, color=colors, edgecolor=[BLUE if gap == 4 else MID for gap in locus1.gap_bp])
    ax.errorbar(x, med, yerr=np.vstack([lower, upper]), fmt="none", ecolor=INK, lw=0.8, capsize=2)
    ax.set_xticks(x, locus1.gap_bp.astype(int))
    ax.set_xlabel("Designed gap (bp)")
    ax.set_ylabel("Median pair ipTM")
    ax.set_title("Pair confidence peaks at 4 bp", loc="left")
    ax.text(
        0.97,
        0.95,
        "4 bp: 0.64\npeak − off-peak: 0.585",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=7.2,
        bbox=dict(boxstyle="round,pad=0.22", fc=WHITE, ec=LIGHT),
    )
    style_axis(ax)

    ax = fig.add_subplot(gs[1, 1])
    panel_label(ax, "D")
    contact = locus1.median_primary_contact_max.to_numpy()
    bars = ax.bar(
        x,
        contact,
        color=[BLUE if gap == 4 else LIGHT for gap in locus1.gap_bp],
        edgecolor=[BLUE if gap == 4 else MID for gap in locus1.gap_bp],
    )
    for bar, n in zip(bars, locus1.median_primary_tf_pair_residue_pair_contacts):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.025,
            f"{int(n)} pairs",
            ha="center",
            fontsize=6.4,
        )
    ax.set_xticks(x, locus1.gap_bp.astype(int))
    ax.set_xlabel("Designed gap (bp)")
    ax.set_ylabel("Maximum inter-TF contact probability")
    ax.set_ylim(0, 0.88)
    ax.set_title("Physical contacts collapse off-peak", loc="left")
    style_axis(ax)

    ax = fig.add_subplot(gs[1, 2])
    panel_label(ax, "E")
    factorial_names = [
        "tcf7l2_runx3_locus1_gap04_intact",
        "tcf7l2_runx3_locus1_gap04_anchor_only",
        "tcf7l2_runx3_locus1_gap04_moved_only",
        "tcf7l2_runx3_locus1_gap04_neither",
    ]
    factor = aggregate.set_index("job_name").loc[factorial_names].reset_index()
    labels = ["Both motifs", "RUNX only", "TCF only", "Neither"]
    vals = factor.median_primary_contact_max.to_numpy()
    colors = [BLUE, RUST_LIGHT, BLUE_LIGHT, LIGHT]
    ax.bar(np.arange(4), vals, color=colors, edgecolor=[BLUE, RUST, BLUE, MID])
    ax.set_xticks(np.arange(4), labels, rotation=24, ha="right")
    ax.set_ylabel("Maximum inter-TF contact probability")
    ax.set_ylim(0, 0.84)
    ax.set_title("Both motifs are required for contact rescue", loc="left")
    for i, value in enumerate(vals):
        ax.text(i, value + 0.025, f"{value:.2f}", ha="center", fontsize=7.2)
    ax.text(
        0.97,
        0.95,
        "factorial rescue = 0.54\n7 physical residue pairs",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=7.1,
        bbox=dict(boxstyle="round,pad=0.22", fc=WHITE, ec=LIGHT),
    )
    style_axis(ax)

    ax = fig.add_subplot(gs[2, :])
    panel_label(ax, "F")
    ax.axis("off")
    statements = [
        (
            0.02,
            "SUPPORTED",
            "chr3 TCF7L2–RUNX3",
            "4-bp activity peak + motif-anchored TF–TF contact",
            BLUE,
        ),
        (
            0.36,
            "NOT GENERALIZED",
            "chr12 TCF7L2–RUNX3",
            "apparent contacts occur away from intended motif sites",
            MID,
        ),
        (
            0.69,
            "POSITIVE-CONTROL FAILURE",
            "MEF2C–TBX5",
            "AF3 did not recover the 9-bp functional optimum",
            RUST,
        ),
    ]
    for x0, status, pair_name, statement, color in statements:
        ax.add_patch(
            patches.FancyBboxPatch(
                (x0, 0.10),
                0.29,
                0.76,
                boxstyle="round,pad=0.012",
                fc=WHITE,
                ec=color,
                lw=1.0,
            )
        )
        ax.text(x0 + 0.015, 0.70, status, color=color, fontsize=7.1, fontweight="bold")
        ax.text(x0 + 0.015, 0.48, pair_name, fontsize=8.2, fontweight="bold")
        ax.text(x0 + 0.015, 0.23, statement, fontsize=7.3)

    fig.suptitle(
        "AlphaFold supports one spacing-dependent assembly but does not validate spacing by itself",
        x=0.01,
        y=0.98,
        ha="left",
        fontsize=13,
        fontweight="bold",
    )
    save(fig, "Fig4_alphafold_spacing_validation")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    figure_1()
    figure_2()
    figure_3()
    figure_4()
    print("Created Figures 1–4 as PNG and PDF.")


if __name__ == "__main__":
    main()
