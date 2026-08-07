#!/usr/bin/env python3
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Test direct disruption of attribution-selected context modules.
# Steps: 1) Load held-out probe and selected loci; 2) create deterministic local shuffles; 3) score effects and save validation outputs.
# =======================================
"""Directly perturb attribution-selected context modules with local mono shuffles."""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from config.paths import publication_path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from scipy import stats
from sklearn.model_selection import GroupKFold


PROJECT = Path(publication_path("MPRA_KOO_PROJECT"))
TASK = Path(publication_path("MPRA_SEQLET_TASK_ROOT"))
OUT = TASK / "outputs"
sys.path.insert(0, str(PROJECT / "koorun"))
from koo_model import CONSTRUCT_MODE, device, load_single  # noqa: E402
from alphagenome_pytorch.utils.sequence import sequence_to_onehot_tensor  # noqa: E402


N_SHUFFLES = 5
BATCH = int(os.environ.get("BATCH", "64"))


def read_fasta(path: Path) -> dict[tuple[str, str], str]:
    records = {}
    key = None
    buf = []
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


def core(seq: str) -> str:
    return seq[15:-15].upper()


def seed_for(*parts: object) -> int:
    text = "|".join(map(str, parts)).encode()
    return int.from_bytes(hashlib.sha256(text).digest()[:8], "little")


def local_shuffles(seq: str, start: int, end: int, n: int) -> list[str]:
    window = seq[start:end]
    rng = np.random.default_rng(seed_for(seq, start, end))
    variants = []
    seen = {window}
    chars = np.array(list(window))
    for _ in range(1000):
        shuffled = "".join(rng.permutation(chars))
        if shuffled in seen:
            continue
        seen.add(shuffled)
        assert Counter(shuffled) == Counter(window)
        variants.append(seq[:start] + shuffled + seq[end:])
        if len(variants) == n:
            break
    return variants


@torch.no_grad()
def score_sequences(
    base,
    sequences: list[str],
    fold: np.ndarray,
    weights: np.ndarray,
    intercepts: np.ndarray,
) -> np.ndarray:
    out = []
    t0 = time.time()
    for start in range(0, len(sequences), BATCH):
        chunk = sequences[start : start + BATCH]
        assembled = base.construct_spec.assemble_sequences(chunk, mode=CONSTRUCT_MODE)
        onehot = torch.stack(
            [
                sequence_to_onehot_tensor(s, dtype=torch.float32, device=device)
                for s in assembled
            ]
        )
        embedding = base.encode(onehot).reshape(len(chunk), -1)
        f = fold[start : start + len(chunk)]
        w = torch.tensor(weights[f], dtype=torch.float32, device=device)
        b = torch.tensor(intercepts[f], dtype=torch.float32, device=device)
        out.append(((embedding * w).sum(dim=1) + b).float().cpu().numpy())
        print(
            f"  scored {min(start+BATCH, len(sequences))}/{len(sequences)} "
            f"({time.time()-t0:.0f}s)",
            flush=True,
        )
    return np.concatenate(out)


def main() -> None:
    modules = pd.read_csv(OUT / "context_probe_candidate_modules.csv")
    attr_summary = pd.read_csv(OUT / "context_probe_attribution_summary.csv").set_index("Coords")
    fa = read_fasta(PROJECT / "bundle/MPRA_AP1_fibroblast.fasta")

    models = np.load(OUT / "mpra_nested_contrast_models.npz", allow_pickle=True)
    target_names = list(models["target_names"].astype(str))
    ti = target_names.index("context_Motif_minus_Full")
    weights = models["outer_coef"][ti]
    intercepts = models["outer_intercept"][ti]
    model_coords = models["coords"].astype(str)
    source = (
        pd.read_csv(PROJECT / "bundle/MPRA_Results_with_TSS_distance.csv")
        .drop_duplicates("Coords")
        .set_index("Coords")
    )
    model_chrom = source.loc[model_coords, "chr"].astype(str).to_numpy()
    assignment = np.full(len(model_coords), -1, int)
    for f, (_, test) in enumerate(
        GroupKFold(5).split(np.zeros(len(model_coords)), groups=model_chrom)
    ):
        assignment[test] = f
    fold_map = {c: int(assignment[i]) for i, c in enumerate(model_coords)}

    records = []
    sequences = []
    folds = []
    for row in modules.itertuples(index=False):
        coord = str(row.Coords)
        seq = core(fa[("MOTIF_SHUF", coord)])
        variants = local_shuffles(seq, int(row.start_0based), int(row.end_0based_exclusive), N_SHUFFLES)
        for replicate, mutant in enumerate(variants, 1):
            sequences.append(mutant)
            folds.append(fold_map[coord])
            records.append(
                {
                    "Coords": coord,
                    "resolved_class": row.resolved_class,
                    "Group": row.Group,
                    "panel_membership": row.panel_membership,
                    "module_type": row.module_type,
                    "module_rank": int(row.module_rank),
                    "start_0based": int(row.start_0based),
                    "end_0based_exclusive": int(row.end_0based_exclusive),
                    "attribution_abs_sum": float(row.attribution_abs_sum),
                    "shuffle_replicate": replicate,
                    "original_window": seq[int(row.start_0based) : int(row.end_0based_exclusive)],
                    "mutant_window": mutant[int(row.start_0based) : int(row.end_0based_exclusive)],
                    "outer_fold": fold_map[coord] + 1,
                }
            )
    print(
        f"device={device}; modules={len(modules)}; mutable variants={len(sequences)}; batch={BATCH}",
        flush=True,
    )
    base = load_single()
    predicted = score_sequences(
        base, sequences, np.asarray(folds, int), weights, intercepts
    )
    del base

    detail = pd.DataFrame(records)
    detail["mutant_context_score"] = predicted
    detail["baseline_context_score"] = detail["Coords"].map(
        attr_summary["context_score_motif"]
    )
    detail["signed_loss"] = detail["baseline_context_score"] - detail["mutant_context_score"]
    detail["absolute_effect"] = detail["signed_loss"].abs()
    detail.to_csv(OUT / "context_module_disruption_all_shuffles.csv", index=False)

    keys = [
        "Coords",
        "resolved_class",
        "Group",
        "panel_membership",
        "module_type",
        "module_rank",
        "start_0based",
        "end_0based_exclusive",
        "attribution_abs_sum",
        "outer_fold",
    ]
    agg = (
        detail.groupby(keys, dropna=False)
        .agg(
            n_shuffles=("absolute_effect", "size"),
            mean_signed_loss=("signed_loss", "mean"),
            sd_signed_loss=("signed_loss", "std"),
            mean_absolute_effect=("absolute_effect", "mean"),
            max_absolute_effect=("absolute_effect", "max"),
        )
        .reset_index()
    )
    agg.to_csv(OUT / "context_module_disruption_summary.csv", index=False)

    # Direct test: within-locus high-attribution modules versus low-attribution controls.
    by_locus = (
        agg.groupby(["Coords", "resolved_class", "module_type"])["mean_absolute_effect"]
        .mean()
        .unstack("module_type")
        .dropna()
        .reset_index()
    )
    wilcoxon = stats.wilcoxon(
        by_locus["high"], by_locus["low_control"], alternative="greater"
    )
    rho = stats.spearmanr(agg["attribution_abs_sum"], agg["mean_absolute_effect"])

    class_stats = {}
    high = agg[agg["module_type"] == "high"]
    arrays = []
    for cls in ["Stable", "Age-only", "Exercise-only", "Reversal"]:
        values = high.loc[high["resolved_class"] == cls, "mean_absolute_effect"].to_numpy()
        arrays.append(values)
        class_stats[cls] = {
            "n_modules": int(len(values)),
            "mean": float(values.mean()),
            "median": float(np.median(values)),
        }
    class_stats["kruskal_p"] = float(stats.kruskal(*arrays).pvalue)

    summary = {
        "n_loci": int(by_locus["Coords"].nunique()),
        "n_modules": int(len(agg)),
        "n_mutant_sequences": int(len(detail)),
        "local_shuffle": {
            "window_bp": int((agg["end_0based_exclusive"] - agg["start_0based"]).median()),
            "replicates_requested": N_SHUFFLES,
            "preserves": "exact window mononucleotide counts",
        },
        "attribution_validation": {
            "high_mean": float(by_locus["high"].mean()),
            "low_control_mean": float(by_locus["low_control"].mean()),
            "high_to_low_ratio": float(
                by_locus["high"].mean() / (by_locus["low_control"].mean() + 1e-12)
            ),
            "paired_wilcoxon_greater_p": float(wilcoxon.pvalue),
            "spearman_attribution_vs_effect": float(rho.statistic),
            "spearman_p": float(rho.pvalue),
        },
        "class_high_module_effect": class_stats,
    }
    with (OUT / "context_module_disruption_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.6))
    for _, row in by_locus.iterrows():
        axes[0].plot(
            [0, 1],
            [row["low_control"], row["high"]],
            color="#999",
            alpha=0.25,
            lw=0.6,
        )
    axes[0].scatter(
        np.zeros(len(by_locus)), by_locus["low_control"], s=8, alpha=0.45, color="#4C72B0"
    )
    axes[0].scatter(
        np.ones(len(by_locus)), by_locus["high"], s=8, alpha=0.45, color="#C44E52"
    )
    axes[0].set_xticks([0, 1], ["low-attribution\ncontrols", "high-attribution\nmodules"])
    axes[0].set_ylabel("mean |context-score change|")
    axes[0].set_title(f"Direct local disruption\npaired P={wilcoxon.pvalue:.1e}")

    classes = ["Stable", "Age-only", "Exercise-only", "Reversal"]
    data = [
        high.loc[high["resolved_class"] == cls, "mean_absolute_effect"] for cls in classes
    ]
    axes[1].boxplot(data, tick_labels=classes, showfliers=False)
    axes[1].tick_params(axis="x", rotation=25)
    axes[1].set_ylabel("high-module |context-score change|")
    axes[1].set_title(f"High modules by class\nKruskal P={class_stats['kruskal_p']:.2g}")
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "context_module_disruption.png", dpi=240, bbox_inches="tight")
    fig.savefig(OUT / "context_module_disruption.pdf", bbox_inches="tight")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
