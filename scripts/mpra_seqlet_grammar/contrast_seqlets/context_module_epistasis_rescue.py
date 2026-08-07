#!/usr/bin/env python3
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Test combinatorial seqlet disruption and composition-preserving rescue.
# Steps: 1) Generate matched single/double/triple perturbations; 2) score held-out probes; 3) calculate epistasis and rescue summaries.
# =======================================
"""Combinatorial context-module disruption and composition-preserving rescue."""

from __future__ import annotations

import json
import os
import sys
import time
from collections import Counter
from itertools import combinations
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


BATCH = int(os.environ.get("BATCH", "64"))
N_DISRUPTION_REPLICATES = 3
LETTERS = ["A", "B", "C"]
COMBOS = ["A", "B", "C", "AB", "AC", "BC", "ABC"]


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


def replace_windows(seq: str, modules: dict[str, dict], combo: str, field: str) -> str:
    chars = list(seq)
    for letter in combo:
        module = modules[letter]
        start, end = module["start"], module["end"]
        replacement = module[field]
        assert len(replacement) == end - start
        chars[start:end] = replacement
    return "".join(chars)


def balanced_rescue(full: str, target: str, modules: dict[str, dict], combo: str) -> str:
    """Set chosen positions to target bases using swaps, preserving total mono counts."""
    assert Counter(full) == Counter(target)
    chars = list(full)
    chosen = set()
    for letter in combo:
        chosen.update(range(modules[letter]["start"], modules[letter]["end"]))
    for pos in sorted(chosen):
        want = target[pos]
        if chars[pos] == want:
            continue
        candidates = [q for q in range(len(chars)) if q not in chosen and chars[q] == want]
        if not candidates:
            raise RuntimeError(f"No compensating base for rescue at {pos}")
        q = candidates[0]
        chars[pos], chars[q] = chars[q], chars[pos]
    rescued = "".join(chars)
    assert Counter(rescued) == Counter(full)
    assert all(rescued[p] == target[p] for p in chosen)
    return rescued


@torch.no_grad()
def score_sequences(base, sequences, fold, weights, intercepts):
    out = []
    t0 = time.time()
    for start in range(0, len(sequences), BATCH):
        chunk = sequences[start : start + BATCH]
        assembled = base.construct_spec.assemble_sequences(chunk, mode=CONSTRUCT_MODE)
        onehot = torch.stack(
            [sequence_to_onehot_tensor(s, dtype=torch.float32, device=device) for s in assembled]
        )
        embedding = base.encode(onehot).reshape(len(chunk), -1)
        f = fold[start : start + len(chunk)]
        w = torch.tensor(weights[f], dtype=torch.float32, device=device)
        b = torch.tensor(intercepts[f], dtype=torch.float32, device=device)
        out.append(((embedding * w).sum(1) + b).float().cpu().numpy())
        print(
            f"  scored {min(start+BATCH,len(sequences))}/{len(sequences)} ({time.time()-t0:.0f}s)",
            flush=True,
        )
    return np.concatenate(out)


def main() -> None:
    detail = pd.read_csv(OUT / "context_module_disruption_all_shuffles.csv")
    agg = pd.read_csv(OUT / "context_module_disruption_summary.csv")
    attr = pd.read_csv(OUT / "context_probe_attribution_summary.csv").set_index("Coords")
    mechanism_coords = set(
        pd.read_csv(OUT / "mpra_in_silico_followup_loci.csv")
        .query("panel_type == 'mechanism-reliable'")["Coords"]
        .astype(str)
    )
    detail["Coords"] = detail["Coords"].astype(str)
    agg["Coords"] = agg["Coords"].astype(str)
    attr.index = attr.index.astype(str)

    # Keep loci with three directly tested high modules and three low controls.
    selected_modules = {}
    for coord in sorted(mechanism_coords):
        selected_modules[coord] = {}
        ok = True
        for module_type in ["high", "low_control"]:
            part = agg[(agg["Coords"] == coord) & (agg["module_type"] == module_type)]
            part = part.nlargest(3, "mean_absolute_effect")
            if len(part) < 3:
                ok = False
                break
            selected_modules[coord][module_type] = part
        if not ok:
            selected_modules.pop(coord, None)
    coords = sorted(selected_modules)

    fa = read_fasta(PROJECT / "bundle/MPRA_AP1_fibroblast.fasta")
    models = np.load(OUT / "mpra_nested_contrast_models.npz", allow_pickle=True)
    names = list(models["target_names"].astype(str))
    ti = names.index("context_Motif_minus_Full")
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

    sequences = []
    folds = []
    records = []
    baseline = {}
    module_meta = {}

    for coord in coords:
        motif = core(fa[("MOTIF_SHUF", coord)])
        full = core(fa[("SHUF_FULL", coord)])
        baseline[coord] = {
            "motif": float(attr.loc[coord, "context_score_motif"]),
            "full": None,
        }
        module_meta[coord] = {}
        for module_type in ["high", "low_control"]:
            part = selected_modules[coord][module_type].reset_index(drop=True)
            mods = {}
            for j, row in part.iterrows():
                letter = LETTERS[j]
                start, end = int(row.start_0based), int(row.end_0based_exclusive)
                mods[letter] = {"start": start, "end": end, "target": motif[start:end]}
            module_meta[coord][module_type] = mods

            # Three matched combinatorial disruption replicates.
            for replicate in range(1, N_DISRUPTION_REPLICATES + 1):
                repmods = {}
                for letter in LETTERS:
                    m = mods[letter]
                    row = detail[
                        (detail["Coords"] == coord)
                        & (detail["module_type"] == module_type)
                        & (detail["start_0based"] == m["start"])
                        & (detail["end_0based_exclusive"] == m["end"])
                        & (detail["shuffle_replicate"] == replicate)
                    ]
                    if row.empty:
                        raise RuntimeError(f"Missing disruption replicate {coord} {module_type} {letter}")
                    repmods[letter] = dict(m, mutant=str(row.iloc[0]["mutant_window"]))
                for combo in COMBOS:
                    sequences.append(replace_windows(motif, repmods, combo, "mutant"))
                    folds.append(fold_map[coord])
                    records.append(
                        {
                            "Coords": coord,
                            "resolved_class": str(attr.loc[coord, "resolved_class"]),
                            "design": "disruption",
                            "module_type": module_type,
                            "replicate": replicate,
                            "combo": combo,
                            "n_modules": len(combo),
                            "outer_fold": fold_map[coord] + 1,
                        }
                    )

            # Composition-preserving rescue from the measured full shuffle.
            for combo in COMBOS:
                sequences.append(balanced_rescue(full, motif, mods, combo))
                folds.append(fold_map[coord])
                records.append(
                    {
                        "Coords": coord,
                        "resolved_class": str(attr.loc[coord, "resolved_class"]),
                        "design": "rescue",
                        "module_type": module_type,
                        "replicate": 1,
                        "combo": combo,
                        "n_modules": len(combo),
                        "outer_fold": fold_map[coord] + 1,
                    }
                )
        # Score the full baseline once per locus.
        sequences.append(full)
        folds.append(fold_map[coord])
        records.append(
            {
                "Coords": coord,
                "resolved_class": str(attr.loc[coord, "resolved_class"]),
                "design": "full_baseline",
                "module_type": "none",
                "replicate": 1,
                "combo": "none",
                "n_modules": 0,
                "outer_fold": fold_map[coord] + 1,
            }
        )

    print(f"device={device}; loci={len(coords)}; sequences={len(sequences)}", flush=True)
    base = load_single()
    scores = score_sequences(
        base, sequences, np.asarray(folds, int), weights, intercepts
    )
    del base
    scored = pd.DataFrame(records)
    scored["context_score"] = scores
    full_scores = scored[scored["design"] == "full_baseline"].set_index("Coords")[
        "context_score"
    ]
    scored["motif_baseline_score"] = scored["Coords"].map(
        {c: baseline[c]["motif"] for c in coords}
    )
    scored["full_baseline_score"] = scored["Coords"].map(full_scores)
    scored.to_csv(OUT / "context_module_epistasis_rescue_scores.csv", index=False)

    # Pairwise and triple epistasis on disruption variants.
    epi_rows = []
    disruption = scored[scored["design"] == "disruption"]
    for (coord, module_type, rep), part in disruption.groupby(
        ["Coords", "module_type", "replicate"]
    ):
        s = part.set_index("combo")["context_score"]
        wt = float(part["motif_baseline_score"].iloc[0])
        for a, b in combinations(LETTERS, 2):
            pair = a + b
            epsilon = float(s[pair] - s[a] - s[b] + wt)
            loss_interaction = -epsilon
            denom = abs(wt - s[a]) + abs(wt - s[b]) + 1e-9
            epi_rows.append(
                {
                    "Coords": coord,
                    "resolved_class": part["resolved_class"].iloc[0],
                    "module_type": module_type,
                    "replicate": int(rep),
                    "interaction": pair,
                    "epsilon_score": epsilon,
                    "loss_interaction": loss_interaction,
                    "absolute_epsilon": abs(epsilon),
                    "fractional_nonadditivity": abs(epsilon) / denom,
                }
            )
        triple_loss_interaction = float(
            (s["A"] + s["B"] + s["C"] - s["ABC"] - 2 * wt)
        )
        epi_rows.append(
            {
                "Coords": coord,
                "resolved_class": part["resolved_class"].iloc[0],
                "module_type": module_type,
                "replicate": int(rep),
                "interaction": "ABC",
                "epsilon_score": -triple_loss_interaction,
                "loss_interaction": triple_loss_interaction,
                "absolute_epsilon": abs(triple_loss_interaction),
                "fractional_nonadditivity": abs(triple_loss_interaction)
                / (
                    abs(wt - s["A"])
                    + abs(wt - s["B"])
                    + abs(wt - s["C"])
                    + 1e-9
                ),
            }
        )
    epi = pd.DataFrame(epi_rows)
    epi.to_csv(OUT / "context_module_epistasis.csv", index=False)

    # Rescue recovery, including exact mononucleotide-preserving low-module controls.
    rescue = scored[scored["design"] == "rescue"].copy()
    gap = rescue["motif_baseline_score"] - rescue["full_baseline_score"]
    rescue["signed_recovery_fraction"] = (
        rescue["context_score"] - rescue["full_baseline_score"]
    ) / gap.replace(0, np.nan)
    rescue["absolute_gap_closed"] = 1 - (
        rescue["motif_baseline_score"] - rescue["context_score"]
    ).abs() / gap.abs().replace(0, np.nan)
    rescue.to_csv(OUT / "context_module_rescue.csv", index=False)

    # Locus-level comparisons high versus low.
    epi_locus = (
        epi.groupby(["Coords", "resolved_class", "module_type"])["absolute_epsilon"]
        .mean()
        .unstack("module_type")
        .dropna()
        .reset_index()
    )
    epi_test = stats.wilcoxon(
        epi_locus["high"], epi_locus["low_control"], alternative="greater"
    )
    triple = rescue[rescue["combo"] == "ABC"].pivot_table(
        index=["Coords", "resolved_class"],
        columns="module_type",
        values="absolute_gap_closed",
    ).dropna().reset_index()
    rescue_test = stats.wilcoxon(
        triple["high"], triple["low_control"], alternative="greater"
    )

    high_epi = epi[epi["module_type"] == "high"]
    class_arrays = []
    class_stats = {}
    for cls in ["Stable", "Age-only", "Exercise-only", "Reversal"]:
        values = high_epi.loc[
            high_epi["resolved_class"] == cls, "absolute_epsilon"
        ].to_numpy()
        class_arrays.append(values)
        class_stats[cls] = {
            "n": int(len(values)),
            "mean_abs_epsilon": float(values.mean()),
            "median_abs_epsilon": float(np.median(values)),
        }
    class_stats["kruskal_p"] = float(stats.kruskal(*class_arrays).pvalue)

    summary = {
        "n_loci": int(len(coords)),
        "n_scored_sequences": int(len(scored)),
        "epistasis": {
            "high_mean_abs_epsilon": float(epi_locus["high"].mean()),
            "low_control_mean_abs_epsilon": float(epi_locus["low_control"].mean()),
            "high_to_low_ratio": float(
                epi_locus["high"].mean() / (epi_locus["low_control"].mean() + 1e-12)
            ),
            "paired_wilcoxon_greater_p": float(epi_test.pvalue),
            "high_loss_synergy_fraction": float(
                np.mean(high_epi["loss_interaction"] > 0)
            ),
            "high_loss_redundancy_fraction": float(
                np.mean(high_epi["loss_interaction"] < 0)
            ),
            "class_comparison": class_stats,
        },
        "rescue": {
            "high_triple_gap_closed_mean": float(triple["high"].mean()),
            "high_triple_gap_closed_median": float(triple["high"].median()),
            "low_triple_gap_closed_mean": float(triple["low_control"].mean()),
            "low_triple_gap_closed_median": float(triple["low_control"].median()),
            "paired_wilcoxon_greater_p": float(rescue_test.pvalue),
            "preserves": "exact full-sequence mononucleotide counts by compensating swaps",
        },
    }
    with (OUT / "context_module_epistasis_rescue_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)

    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.7))
    for _, row in epi_locus.iterrows():
        axes[0].plot(
            [0, 1], [row["low_control"], row["high"]], color="#999", alpha=0.25, lw=0.6
        )
    axes[0].scatter(np.zeros(len(epi_locus)), epi_locus["low_control"], s=8, alpha=.5)
    axes[0].scatter(np.ones(len(epi_locus)), epi_locus["high"], s=8, alpha=.5, color="#C44E52")
    axes[0].set_xticks([0, 1], ["low controls", "causal modules"])
    axes[0].set_ylabel("mean |epistasis ε|")
    axes[0].set_title(f"Module non-additivity\npaired P={epi_test.pvalue:.1e}")

    rescue_plot = rescue.groupby(["module_type", "n_modules"])["absolute_gap_closed"].agg(
        ["mean", "sem"]
    )
    for module_type, color, label in [
        ("high", "#C44E52", "causal modules"),
        ("low_control", "#4C72B0", "low controls"),
    ]:
        d = rescue_plot.loc[module_type]
        axes[1].errorbar(d.index, d["mean"], yerr=d["sem"], marker="o", color=color, label=label)
    axes[1].axhline(0, color="#777", lw=.6)
    axes[1].set_xticks([1, 2, 3])
    axes[1].set_xlabel("modules restored")
    axes[1].set_ylabel("fraction of motif↔full gap closed")
    axes[1].set_title(f"Composition-preserving rescue\ntriple paired P={rescue_test.pvalue:.1e}")
    axes[1].legend(frameon=False, fontsize=7)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "context_module_epistasis_rescue.png", dpi=240, bbox_inches="tight")
    fig.savefig(OUT / "context_module_epistasis_rescue.pdf", bbox_inches="tight")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
