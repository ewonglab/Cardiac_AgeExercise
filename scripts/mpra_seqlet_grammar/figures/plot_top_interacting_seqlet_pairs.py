#!/usr/bin/env python3
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Plot the strongest interacting seqlet pairs from a versioned table.
# Steps: 1) Read pair statistics; 2) rank the top pairs; 3) write PNG and PDF horizontal-bar figures.
# =======================================
"""Reproduce the horizontal bar chart of the top interacting seqlet pairs."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pairs",
        type=Path,
        default=Path(__file__).resolve().parents[1]
        / "example_inputs"
        / "context_module_pairs_with_motifs_uncorrected.csv",
        help="Pair-level table produced by redo_uncorrected_motif_analysis.py.",
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=Path(__file__).resolve().parent / "Fig_top_interacting_seqlet_pairs",
    )
    parser.add_argument("--top-n", type=int, default=9)
    args = parser.parse_args()

    pairs = pd.read_csv(args.pairs)
    required = {"epistasis_rank", "motif_pair_name", "mean_absolute_epsilon"}
    missing = required.difference(pairs.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    top = (
        pairs.sort_values(["epistasis_rank", "mean_absolute_epsilon"], ascending=[True, False])
        .head(args.top_n)
        .sort_values("mean_absolute_epsilon", ascending=True)
    )

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 12,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, ax = plt.subplots(figsize=(10.6, 5.35))
    y = range(len(top))
    ax.barh(
        y,
        top["mean_absolute_epsilon"],
        color="#9A6B57",
        edgecolor="none",
        height=0.72,
    )
    ax.set_yticks(list(y), top["motif_pair_name"])
    ax.set_xlabel(r"Mean $|\mathrm{epistasis}\ \varepsilon|$", labelpad=16)
    ax.set_ylabel("")
    ax.set_xlim(left=0)
    ax.tick_params(axis="y", labelsize=11.5)
    ax.tick_params(axis="x", labelsize=10.5)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)
    fig.subplots_adjust(left=0.54, right=0.98, bottom=0.22, top=0.98)

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output_prefix.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(args.output_prefix.with_suffix(".pdf"), facecolor="white")
    print(args.output_prefix.with_suffix(".png"))


if __name__ == "__main__":
    main()
