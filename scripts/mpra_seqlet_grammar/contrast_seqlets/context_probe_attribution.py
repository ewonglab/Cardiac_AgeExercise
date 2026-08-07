#!/usr/bin/env python3
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Compute held-out context-probe attributions and nominate seqlets.
# Steps: 1) Load contrast-probe inputs; 2) calculate base-level gradients; 3) select high and low attribution windows for testing.
# =======================================
"""Gradient attribution and candidate-module mapping for the context contrast probe."""

from __future__ import annotations

import json
import os
import sys
import time
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


BATCH = int(os.environ.get("BATCH", "8"))
WINDOW = 12
N_MODULES = 3


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


def gini(x: np.ndarray) -> float:
    x = np.abs(np.asarray(x, float))
    if not np.any(x):
        return 0.0
    x = np.sort(x)
    n = len(x)
    return float((2 * np.sum(np.arange(1, n + 1) * x) / (n * x.sum())) - (n + 1) / n)


def pick_windows(
    importance: np.ndarray,
    excluded: np.ndarray,
    strongest: bool,
    n: int = N_MODULES,
) -> list[tuple[int, int, float]]:
    candidates = []
    for start in range(0, len(importance) - WINDOW + 1):
        end = start + WINDOW
        if excluded[start:end].any():
            continue
        candidates.append((start, end, float(importance[start:end].sum())))
    candidates.sort(key=lambda x: x[2], reverse=strongest)
    selected = []
    for item in candidates:
        start, end, _ = item
        if all(end <= a or start >= b for a, b, _ in selected):
            selected.append(item)
        if len(selected) == n:
            break
    return selected


def main() -> None:
    panel = pd.read_csv(OUT / "mpra_in_silico_followup_loci.csv")
    membership = panel.groupby("Coords")["panel_type"].agg(lambda x: ";".join(sorted(set(x))))
    panel = panel.drop_duplicates("Coords").set_index("Coords")
    panel["panel_membership"] = membership
    panel = panel.reset_index()
    coords = panel["Coords"].astype(str).to_numpy()

    fa = read_fasta(PROJECT / "bundle/MPRA_AP1_fibroblast.fasta")
    seq = {
        kind: [core(fa[(kind, c)]) for c in coords]
        for kind in ["ORIG", "MOTIF_SHUF"]
    }
    models = np.load(OUT / "mpra_nested_contrast_models.npz", allow_pickle=True)
    names = list(models["target_names"].astype(str))
    ti = names.index("context_Motif_minus_Full")
    model_coords = models["coords"].astype(str)
    source = (
        pd.read_csv(PROJECT / "bundle/MPRA_Results_with_TSS_distance.csv")
        .drop_duplicates("Coords")
        .set_index("Coords")
    )
    model_chrom = source.loc[model_coords, "chr"].astype(str).to_numpy()
    outer_assignment = np.full(len(model_coords), -1, dtype=int)
    for fold, (_, test) in enumerate(
        GroupKFold(5).split(np.zeros(len(model_coords)), groups=model_chrom)
    ):
        outer_assignment[test] = fold
    fold_map = {c: int(outer_assignment[i]) for i, c in enumerate(model_coords)}
    locus_folds = np.array([fold_map[c] for c in coords], dtype=int)
    outer_weight = models["outer_coef"][ti]
    outer_intercept = models["outer_intercept"][ti]

    base = load_single()
    left = len(base.construct_spec.left_adapter or "")
    core_slice = slice(left, left + 170)
    print(f"device={device}; loci={len(coords)}; batch={BATCH}; core={core_slice}", flush=True)

    attrs = {}
    scores = {}
    t0 = time.time()
    torch.set_grad_enabled(True)
    for kind in ["ORIG", "MOTIF_SHUF"]:
        onehot_all = np.zeros((len(coords), 170, 4), np.float32)
        hyp_all = np.zeros_like(onehot_all)
        score_all = np.zeros(len(coords), np.float32)
        for start in range(0, len(coords), BATCH):
            chunk = seq[kind][start : start + BATCH]
            assembled = base.construct_spec.assemble_sequences(chunk, mode=CONSTRUCT_MODE)
            oh = torch.stack(
                [
                    sequence_to_onehot_tensor(s, dtype=torch.float32, device=device)
                    for s in assembled
                ]
            )
            oh.requires_grad_(True)
            embedding = base.encode(oh).reshape(len(chunk), -1)
            fold_chunk = locus_folds[start : start + len(chunk)]
            weight = torch.tensor(
                outer_weight[fold_chunk], dtype=torch.float32, device=device
            )
            intercept = torch.tensor(
                outer_intercept[fold_chunk], dtype=torch.float32, device=device
            )
            out = (embedding * weight).sum(dim=1) + intercept
            base.zero_grad(set_to_none=True)
            if oh.grad is not None:
                oh.grad = None
            out.sum().backward()
            grad = oh.grad.detach().cpu().numpy()[:, core_slice, :]
            onehot = oh.detach().cpu().numpy()[:, core_slice, :]
            stop = start + len(chunk)
            onehot_all[start:stop] = onehot
            hyp_all[start:stop] = grad
            score_all[start:stop] = out.detach().cpu().numpy()
            print(
                f"  {kind}: {stop}/{len(coords)} ({time.time()-t0:.0f}s)",
                flush=True,
            )
        hyp_all -= hyp_all.mean(axis=2, keepdims=True)
        attrs[kind] = {
            "onehot": onehot_all,
            "hypothetical": hyp_all,
            "contrib": hyp_all * onehot_all,
        }
        scores[kind] = score_all
    torch.set_grad_enabled(False)
    del base

    np.savez_compressed(
        OUT / "context_probe_panel_attributions.npz",
        coords=coords,
        orig_onehot=attrs["ORIG"]["onehot"],
        orig_hypothetical=attrs["ORIG"]["hypothetical"],
        orig_contrib=attrs["ORIG"]["contrib"],
        motif_onehot=attrs["MOTIF_SHUF"]["onehot"],
        motif_hypothetical=attrs["MOTIF_SHUF"]["hypothetical"],
        motif_contrib=attrs["MOTIF_SHUF"]["contrib"],
        orig_score=scores["ORIG"],
        motif_score=scores["MOTIF_SHUF"],
        outer_fold=locus_folds,
    )

    modules = []
    summary_rows = []
    motif_actual = attrs["MOTIF_SHUF"]["contrib"].sum(axis=2)
    orig_actual = attrs["ORIG"]["contrib"].sum(axis=2)
    for i, coord in enumerate(coords):
        edit_mask = np.array(list(seq["ORIG"][i])) != np.array(list(seq["MOTIF_SHUF"][i]))
        importance = np.abs(motif_actual[i])
        high = pick_windows(importance, edit_mask, strongest=True)
        low = pick_windows(importance, edit_mask, strongest=False)
        for module_type, chosen in [("high", high), ("low_control", low)]:
            for rank, (start, end, value) in enumerate(chosen, 1):
                modules.append(
                    {
                        "Coords": coord,
                        "resolved_class": panel.loc[i, "resolved_class"],
                        "Group": panel.loc[i, "Group"],
                        "panel_membership": panel.loc[i, "panel_membership"],
                        "module_type": module_type,
                        "module_rank": rank,
                        "start_0based": start,
                        "end_0based_exclusive": end,
                        "sequence_motif_shuf": seq["MOTIF_SHUF"][i][start:end],
                        "attribution_abs_sum": value,
                        "overlaps_AP1_shuffle": bool(edit_mask[start:end].any()),
                    }
                )
        order = np.argsort(-importance)
        frac = np.cumsum(importance[order]) / (importance.sum() + 1e-12)
        summary_rows.append(
            {
                "Coords": coord,
                "resolved_class": panel.loc[i, "resolved_class"],
                "panel_membership": panel.loc[i, "panel_membership"],
                "context_score_orig": float(scores["ORIG"][i]),
                "context_score_motif": float(scores["MOTIF_SHUF"][i]),
                "outer_fold": int(locus_folds[i] + 1),
                "attribution_gini_motif": gini(importance),
                "positions_for_50pct": int(np.searchsorted(frac, 0.5) + 1),
                "n_AP1_shuffled_bases": int(edit_mask.sum()),
            }
        )
    modules = pd.DataFrame(modules)
    summaries = pd.DataFrame(summary_rows)
    modules.to_csv(OUT / "context_probe_candidate_modules.csv", index=False)
    summaries.to_csv(OUT / "context_probe_attribution_summary.csv", index=False)

    # Class comparison is descriptive; perturbation validation follows next.
    class_stats = {}
    for metric in ["attribution_gini_motif", "positions_for_50pct"]:
        class_stats[metric] = {}
        arrays = []
        for cls in ["Stable", "Age-only", "Exercise-only", "Reversal"]:
            values = summaries.loc[summaries["resolved_class"] == cls, metric].to_numpy(float)
            arrays.append(values)
            class_stats[metric][cls] = {
                "n": int(len(values)),
                "mean": float(values.mean()),
                "median": float(np.median(values)),
            }
        class_stats[metric]["kruskal_p"] = float(stats.kruskal(*arrays).pvalue)

    with (OUT / "context_probe_attribution_summary.json").open("w") as handle:
        json.dump(
            {
                "n_loci": int(len(coords)),
                "window": WINDOW,
                "modules_per_locus": N_MODULES,
                "class_stats": class_stats,
            },
            handle,
            indent=2,
        )

    # Compact concentration figure and representative attribution maps.
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
    classes = ["Stable", "Age-only", "Exercise-only", "Reversal"]
    data = [
        summaries.loc[summaries["resolved_class"] == cls, "attribution_gini_motif"]
        for cls in classes
    ]
    axes[0].boxplot(data, tick_labels=classes, showfliers=False)
    axes[0].set_ylabel("context attribution Gini")
    axes[0].tick_params(axis="x", rotation=25)
    axes[0].set_title("Attribution concentration")

    # One high-confidence representative per class, chosen by high-window score.
    hm = []
    labels = []
    high_score = modules[modules["module_type"] == "high"].groupby("Coords")[
        "attribution_abs_sum"
    ].sum()
    for cls in classes:
        candidates = summaries.loc[summaries["resolved_class"] == cls, "Coords"]
        chosen = high_score.loc[high_score.index.intersection(candidates)].idxmax()
        idx = int(np.where(coords == chosen)[0][0])
        v = motif_actual[idx]
        hm.append(v / (np.max(np.abs(v)) + 1e-12))
        labels.append(cls)
    im = axes[1].imshow(np.vstack(hm), aspect="auto", cmap="coolwarm", vmin=-1, vmax=1)
    axes[1].set_yticks(range(4), labels)
    axes[1].set_xlabel("position in 170-bp element")
    axes[1].set_title("Representative signed context attributions")
    fig.colorbar(im, ax=axes[1], fraction=0.04, pad=0.02)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "context_probe_attribution.png", dpi=240, bbox_inches="tight")
    fig.savefig(OUT / "context_probe_attribution.pdf", bbox_inches="tight")
    print(json.dumps(class_stats, indent=2), flush=True)


if __name__ == "__main__":
    main()
