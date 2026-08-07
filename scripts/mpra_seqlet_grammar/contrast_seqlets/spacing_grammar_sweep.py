#!/usr/bin/env python3
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Run broad AP-1-centred motif-spacing sensitivity analyses.
# Steps: 1) Construct composition-preserving spacing variants; 2) score model responses; 3) compare curves across enhancer classes and folds.
# =======================================
"""Motif-preserving spacing sweeps for MPRA-derived enhancer grammar.

The primary perturbation moves an intact PWM footprint by block transposition.
The sensitivity analysis moves the complete 12-bp attribution module containing
that footprint.  Both operations preserve the exact 170-bp mononucleotide
composition because no base is added, removed, or substituted.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations
from pathlib import Path
import json
import os
import re
import sys
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.model_selection import GroupKFold
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from config.paths import publication_path


PROJECT = Path(publication_path("MPRA_KOO_PROJECT"))
TASK = Path(publication_path("MPRA_SEQLET_TASK_ROOT"))
OUT = TASK / "outputs"
MEME = PROJECT / "koorun/out/jaspar_panel.meme"
sys.path.insert(0, str(PROJECT / "koorun"))
from koo_model import CONSTRUCT_MODE, device, load_single  # noqa: E402
from alphagenome_pytorch.utils.sequence import sequence_to_onehot_tensor  # noqa: E402


BATCH = int(os.environ.get("BATCH", "64"))
DRY_RUN = bool(int(os.environ.get("DRY_RUN", "0")))
ANALYZE_ONLY = bool(int(os.environ.get("ANALYZE_ONLY", "0")))
DISCOVERY_RANK_MAX = int(os.environ.get("DISCOVERY_RANK_MAX", "60"))
RELATIVE_THRESHOLD = 0.80
FINE_GAPS = list(range(0, 31)) + [32, 34, 36, 40, 44, 48, 56, 64]
COARSE_GAPS = [0, 2, 4, 6, 8, 10, 12, 16, 20, 24, 28, 32, 40, 48, 56, 64, 80]
TARGET_RANKS = {1, 2, 4, 6, 15}
BASES = "ACGT"
BASE_INDEX = {b: i for i, b in enumerate(BASES)}


def read_fasta(path: Path) -> dict[tuple[str, str], str]:
    records: dict[tuple[str, str], str] = {}
    key = None
    buf: list[str] = []
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


def parse_meme(path: Path) -> dict[str, np.ndarray]:
    lines = path.read_text().splitlines()
    motifs: dict[str, np.ndarray] = {}
    i = 0
    while i < len(lines):
        if lines[i].startswith("MOTIF "):
            name = lines[i].split(maxsplit=1)[1]
            i += 1
            while i < len(lines) and "letter-probability matrix" not in lines[i]:
                i += 1
            width = int(re.search(r"\bw=\s*(\d+)", lines[i]).group(1))
            motifs[name] = np.asarray(
                [[float(x) for x in row.split()[:4]] for row in lines[i + 1 : i + 1 + width]],
                float,
            )
            i += width
        i += 1
    return motifs


def motif_gene(name: str) -> str:
    return re.sub(r"_MA\d+\.\d+$", "", str(name))


def transpose_block(seq: str, start: int, end: int, target_start: int) -> str:
    """Move seq[start:end] to target_start in final coordinates."""
    width = end - start
    if not (0 <= start < end <= len(seq)):
        raise ValueError("invalid source block")
    if not (0 <= target_start <= len(seq) - width):
        raise ValueError("invalid target block")
    if target_start == start:
        return seq
    # Overlapping source/target placements are ambiguous block rotations and are
    # deliberately excluded from the causal spacing panel.
    if max(start, target_start) < min(end, target_start + width):
        raise ValueError("overlapping source and target block")
    block = seq[start:end]
    remainder = seq[:start] + seq[end:]
    # ``target_start`` is defined in final-sequence coordinates.  Inserting the
    # block at that index in the shortened remainder places the first block base
    # at exactly the requested final coordinate, irrespective of move direction.
    insertion = target_start
    result = remainder[:insertion] + block + remainder[insertion:]
    assert len(result) == len(seq)
    assert Counter(result) == Counter(seq)
    assert result[target_start : target_start + width] == block
    return result


def scan_pwm(seq: str, pwm: np.ndarray) -> tuple[float, int]:
    encoded = np.asarray([BASE_INDEX[b] for b in seq], int)
    logodds = np.log2((pwm + 1e-5) / 0.25)
    rc = logodds[::-1, ::-1]
    minimum = float(np.min(logodds, axis=1).sum())
    maximum = float(np.max(logodds, axis=1).sum())
    threshold = minimum + RELATIVE_THRESHOLD * (maximum - minimum)
    best = -np.inf
    hits = 0
    width = len(pwm)
    for pos in range(len(seq) - width + 1):
        plus = float(logodds[np.arange(width), encoded[pos : pos + width]].sum())
        minus = float(rc[np.arange(width), encoded[pos : pos + width]].sum())
        value = max(plus, minus)
        best = max(best, value)
        hits += int(value >= threshold)
    relative = (best - minimum) / max(1e-9, maximum - minimum)
    return float(relative), int(hits)


def build_pair_record(row: pd.Series, motifs: dict[str, np.ndarray]) -> dict:
    members = []
    for side in [1, 2]:
        jaspar = str(row[f"module_{side}_motif_top_jaspar_motif"])
        if jaspar not in motifs:
            raise KeyError(jaspar)
        module_start = int(row[f"module_{side}_start"])
        module_end = int(row[f"module_{side}_end"])
        offset = int(row[f"module_{side}_motif_top_motif_offset"])
        motif_start = module_start + offset
        motif_end = motif_start + len(motifs[jaspar])
        if motif_end > module_end:
            raise ValueError(f"motif outside module: {row.Coords} side {side}")
        members.append(
            {
                "side": side,
                "letter": str(row[f"module_{side}"]),
                "motif": motif_gene(jaspar),
                "jaspar": jaspar,
                "strand": str(row[f"module_{side}_motif_top_strand"]),
                "relative": float(row[f"module_{side}_motif_top_relative_PWM_score"]),
                "motif_start": motif_start,
                "motif_end": motif_end,
                "module_start": module_start,
                "module_end": module_end,
                "motif_offset": offset,
            }
        )
    members.sort(key=lambda x: (x["motif_start"], x["motif_end"]))
    upstream, downstream = members
    native_gap = downstream["motif_start"] - upstream["motif_end"]
    if native_gap < 0:
        raise ValueError(f"overlapping motif footprints: {row.Coords}")
    canonical = " + ".join(sorted([upstream["motif"], downstream["motif"]]))
    return {
        "Coords": str(row.Coords),
        "resolved_class": str(row.resolved_class),
        "epistasis_rank": int(row.epistasis_rank),
        "mean_absolute_epsilon": float(row.mean_absolute_epsilon),
        "rescue_high_gap_closed": float(row.rescue_high_gap_closed),
        "motif_pair_name": str(row.motif_pair_name),
        "canonical_pair": canonical,
        "upstream": upstream,
        "downstream": downstream,
        "native_motif_gap_bp": int(native_gap),
    }


def variant_for_gap(
    seq: str,
    pair: dict,
    gap: int,
    anchor: str,
    mode: str,
) -> tuple[str, int, int]:
    """Return sequence and final moved motif coordinates."""
    up, down = pair["upstream"], pair["downstream"]
    if anchor == "upstream":
        fixed, moved = up, down
        desired_motif_start = fixed["motif_end"] + gap
    elif anchor == "downstream":
        fixed, moved = down, up
        desired_motif_start = fixed["motif_start"] - gap - (moved["motif_end"] - moved["motif_start"])
    else:
        raise ValueError(anchor)

    if mode == "motif_core":
        source_start, source_end = moved["motif_start"], moved["motif_end"]
        target_start = desired_motif_start
        final_motif_start = target_start
    elif mode == "module_12bp":
        source_start, source_end = moved["module_start"], moved["module_end"]
        target_start = desired_motif_start - moved["motif_offset"]
        final_motif_start = target_start + moved["motif_offset"]
        # Do not place one attribution module across the other module.
        fixed_start, fixed_end = fixed["module_start"], fixed["module_end"]
        if max(target_start, fixed_start) < min(target_start + 12, fixed_end):
            raise ValueError("module overlap")
    else:
        raise ValueError(mode)

    result = transpose_block(seq, source_start, source_end, target_start)
    width = moved["motif_end"] - moved["motif_start"]
    final_motif_end = final_motif_start + width
    expected = seq[moved["motif_start"] : moved["motif_end"]]
    assert result[final_motif_start:final_motif_end] == expected

    # The fixed motif remains at its original absolute coordinates because the
    # moved block stays on the same side of it.
    fixed_seq = seq[fixed["motif_start"] : fixed["motif_end"]]
    assert result[fixed["motif_start"] : fixed["motif_end"]] == fixed_seq
    return result, int(final_motif_start), int(final_motif_end)


def fold_map_for_coords(coords: list[str]) -> dict[str, int]:
    models = np.load(OUT / "mpra_nested_contrast_models.npz", allow_pickle=True)
    model_coords = models["coords"].astype(str)
    source = (
        pd.read_csv(PROJECT / "bundle/MPRA_Results_with_TSS_distance.csv")
        .drop_duplicates("Coords")
        .set_index("Coords")
    )
    chrom = source.loc[model_coords, "chr"].astype(str).to_numpy()
    assignment = np.full(len(model_coords), -1, int)
    for fold, (_, test) in enumerate(GroupKFold(5).split(np.zeros(len(model_coords)), groups=chrom)):
        assignment[test] = fold
    lookup = {coord: int(assignment[i]) for i, coord in enumerate(model_coords)}
    return {coord: lookup[coord] for coord in coords}


@torch.no_grad()
def score_unique_sequences(
    sequences: list[str], folds: np.ndarray
) -> dict[str, np.ndarray]:
    models = np.load(OUT / "mpra_nested_contrast_models.npz", allow_pickle=True)
    names = list(models["target_names"].astype(str))
    weights = models["outer_coef"]
    intercepts = models["outer_intercept"]
    base = load_single()
    output = {name: [] for name in names}
    output["mpra_activity"] = []
    t0 = time.time()
    for start in range(0, len(sequences), BATCH):
        chunk = sequences[start : start + BATCH]
        chunk_folds = folds[start : start + len(chunk)]
        assembled = base.construct_spec.assemble_sequences(chunk, mode=CONSTRUCT_MODE)
        onehot = torch.stack(
            [sequence_to_onehot_tensor(s, dtype=torch.float32, device=device) for s in assembled]
        )
        encoder_output = base.encode(onehot)
        output["mpra_activity"].append(
            base.head(encoder_output).float().cpu().numpy().reshape(-1)
        )
        embedding = encoder_output.reshape(len(chunk), -1)
        for target_index, name in enumerate(names):
            w = torch.tensor(weights[target_index, chunk_folds], dtype=torch.float32, device=device)
            b = torch.tensor(intercepts[target_index, chunk_folds], dtype=torch.float32, device=device)
            output[name].append(((embedding * w).sum(1) + b).float().cpu().numpy())
        print(
            f"  scored {min(start+BATCH, len(sequences))}/{len(sequences)} "
            f"({time.time()-t0:.0f}s)",
            flush=True,
        )
    del base
    return {key: np.concatenate(value) for key, value in output.items()}


def centered_curve(frame: pd.DataFrame, value: str) -> pd.Series:
    curve = frame.groupby("target_gap_bp")[value].mean().sort_index()
    return curve - curve.mean()


def pairwise_curve_agreement(curves: list[pd.Series]) -> float:
    values = []
    for a, b in combinations(curves, 2):
        common = a.index.intersection(b.index)
        if len(common) < 6:
            continue
        x, y = a.loc[common].to_numpy(), b.loc[common].to_numpy()
        if np.std(x) == 0 or np.std(y) == 0:
            continue
        values.append(float(np.corrcoef(x, y)[0, 1]))
    return float(np.mean(values)) if values else np.nan


def helical_amplitude(gaps: np.ndarray, values: np.ndarray, period: float = 10.5) -> float:
    x = np.column_stack(
        [np.ones(len(gaps)), gaps, np.sin(2 * np.pi * gaps / period), np.cos(2 * np.pi * gaps / period)]
    )
    coef, *_ = np.linalg.lstsq(x, values, rcond=None)
    return float(np.sqrt(coef[2] ** 2 + coef[3] ** 2))


def summarize(raw: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    score = "context_Motif_minus_Full"
    curve_rows = []
    keys = ["Coords", "epistasis_rank", "canonical_pair", "grid", "design_mode", "anchor"]
    for key, part in raw.groupby(keys):
        curve = part.groupby("target_gap_bp", as_index=False).agg(
            context_score=(score, "mean"), mpra_activity=("mpra_activity", "mean")
        )
        if len(curve) < 6:
            continue
        native_gap = int(part.native_motif_gap_bp.iloc[0])
        baseline_context = float(part.native_context_Motif_minus_Full.iloc[0])
        baseline_activity = float(part.native_mpra_activity.iloc[0])
        best = curve.loc[curve.context_score.idxmax()]
        activity_best = curve.loc[curve.mpra_activity.idxmax()]
        curve_rows.append(
            {
                **dict(zip(keys, key)),
                "resolved_class": str(part.resolved_class.iloc[0]),
                "motif_pair_name": str(part.motif_pair_name.iloc[0]),
                "n_gaps": int(len(curve)),
                "native_motif_gap_bp": native_gap,
                "mean_absolute_epsilon": float(part.mean_absolute_epsilon.iloc[0]),
                "rescue_high_gap_closed": float(part.rescue_high_gap_closed.iloc[0]),
                "best_context_gap_bp": int(best.target_gap_bp),
                "best_activity_gap_bp": int(activity_best.target_gap_bp),
                "context_amplitude": float(curve.context_score.max() - curve.context_score.min()),
                "activity_amplitude": float(curve.mpra_activity.max() - curve.mpra_activity.min()),
                "best_context_minus_native": float(best.context_score - baseline_context),
                "best_activity_minus_native": float(activity_best.mpra_activity - baseline_activity),
                "native_context_percentile": float(np.mean(curve.context_score <= baseline_context)),
                "native_activity_percentile": float(np.mean(curve.mpra_activity <= baseline_activity)),
                "context_activity_curve_r": float(
                    stats.pearsonr(curve.context_score, curve.mpra_activity).statistic
                    if np.std(curve.context_score) and np.std(curve.mpra_activity)
                    else np.nan
                ),
                "context_helical_10_5bp_amplitude": helical_amplitude(
                    curve.target_gap_bp.to_numpy(float), curve.context_score.to_numpy(float)
                ),
                "activity_helical_10_5bp_amplitude": helical_amplitude(
                    curve.target_gap_bp.to_numpy(float), curve.mpra_activity.to_numpy(float)
                ),
            }
        )
    curve_summary = pd.DataFrame(curve_rows)

    pair_rows = []
    pair_keys = ["Coords", "epistasis_rank", "canonical_pair", "grid"]
    for key, part in raw.groupby(pair_keys):
        curves_context = [centered_curve(z, score) for _, z in part.groupby(["design_mode", "anchor"])]
        curves_activity = [centered_curve(z, "mpra_activity") for _, z in part.groupby(["design_mode", "anchor"])]
        meta_context = pd.concat(curves_context, axis=1).mean(axis=1).dropna()
        meta_activity = pd.concat(curves_activity, axis=1).mean(axis=1).dropna()
        common = meta_context.index.intersection(meta_activity.index)
        pair_rows.append(
            {
                **dict(zip(pair_keys, key)),
                "resolved_class": str(part.resolved_class.iloc[0]),
                "motif_pair_name": str(part.motif_pair_name.iloc[0]),
                "native_motif_gap_bp": int(part.native_motif_gap_bp.iloc[0]),
                "mean_absolute_epsilon": float(part.mean_absolute_epsilon.iloc[0]),
                "rescue_high_gap_closed": float(part.rescue_high_gap_closed.iloc[0]),
                "n_curves": int(len(curves_context)),
                "context_curve_agreement": pairwise_curve_agreement(curves_context),
                "activity_curve_agreement": pairwise_curve_agreement(curves_activity),
                "meta_best_context_gap_bp": int(meta_context.idxmax()),
                "meta_best_activity_gap_bp": int(meta_activity.idxmax()),
                "meta_context_amplitude": float(meta_context.max() - meta_context.min()),
                "meta_activity_amplitude": float(meta_activity.max() - meta_activity.min()),
                "meta_context_activity_r": float(
                    np.corrcoef(meta_context.loc[common], meta_activity.loc[common])[0, 1]
                    if len(common) >= 6 else np.nan
                ),
            }
        )
    pair_summary = pd.DataFrame(pair_rows)

    # Motif-family meta-curves across loci in the discovery grid.  Each curve is
    # centered and scaled before averaging so a single high-amplitude locus does
    # not dominate the inferred shared spacing preference.
    family_rows = []
    discovery = raw[raw.grid == "coarse_discovery"]
    for canonical, part in discovery.groupby("canonical_pair"):
        if part.Coords.nunique() < 2:
            continue
        context_curves = []
        activity_curves = []
        locus_context = []
        for (coord, mode, anchor), z in part.groupby(["Coords", "design_mode", "anchor"]):
            c = centered_curve(z, score)
            a = centered_curve(z, "mpra_activity")
            if c.std() > 0:
                c = c / c.std()
            if a.std() > 0:
                a = a / a.std()
            context_curves.append(c.rename((coord, mode, anchor)))
            activity_curves.append(a.rename((coord, mode, anchor)))
        context_matrix = pd.concat(context_curves, axis=1)
        activity_matrix = pd.concat(activity_curves, axis=1)
        meta_context = context_matrix.mean(axis=1).dropna()
        meta_activity = activity_matrix.mean(axis=1).dropna()
        for coord in part.Coords.unique():
            cols = [c for c in context_matrix.columns if c[0] == coord]
            locus_context.append(context_matrix[cols].mean(axis=1).rename(coord))
        locus_matrix = pd.concat(locus_context, axis=1)
        loo_best = []
        for coord in locus_matrix.columns:
            remaining = locus_matrix.drop(columns=coord)
            if remaining.shape[1]:
                loo_best.append(int(remaining.mean(axis=1).idxmax()))
        common = meta_context.index.intersection(meta_activity.index)
        family_rows.append(
            {
                "canonical_pair": canonical,
                "n_loci": int(part.Coords.nunique()),
                "n_curves": int(len(context_curves)),
                "best_context_gap_bp": int(meta_context.idxmax()),
                "best_activity_gap_bp": int(meta_activity.idxmax()),
                "leave_one_locus_out_best_gaps": ";".join(map(str, loo_best)),
                "locus_context_curve_agreement": pairwise_curve_agreement(locus_context),
                "all_context_curve_agreement": pairwise_curve_agreement(context_curves),
                "all_activity_curve_agreement": pairwise_curve_agreement(activity_curves),
                "meta_context_activity_r": float(
                    np.corrcoef(meta_context.loc[common], meta_activity.loc[common])[0, 1]
                    if len(common) >= 6 else np.nan
                ),
                "mean_pair_epistasis": float(part.groupby("Coords").mean_absolute_epsilon.first().mean()),
                "max_pair_epistasis": float(part.groupby("Coords").mean_absolute_epsilon.first().max()),
                "mean_pair_rescue": float(part.groupby("Coords").rescue_high_gap_closed.first().mean()),
            }
        )
    family_summary = pd.DataFrame(family_rows)
    if len(family_summary):
        family_summary = family_summary.sort_values(
            ["n_loci", "locus_context_curve_agreement", "max_pair_epistasis"],
            ascending=[False, False, False],
        )
    return curve_summary, pair_summary, family_summary


def make_figure(raw: pd.DataFrame, family_summary: pd.DataFrame) -> None:
    score = "context_Motif_minus_Full"
    targets = raw[(raw.grid == "fine_targeted") & raw.epistasis_rank.isin(sorted(TARGET_RANKS))]
    ranks = sorted(targets.epistasis_rank.unique())
    fig, axes = plt.subplots(len(ranks), 2, figsize=(10.5, 2.4 * len(ranks)), squeeze=False)
    colors = {"motif_core": "#2563eb", "module_12bp": "#d97706"}
    for ri, rank in enumerate(ranks):
        part = targets[targets.epistasis_rank == rank]
        title = f"rank {rank}: {part.motif_pair_name.iloc[0]} ({part.resolved_class.iloc[0]})"
        for ci, value in enumerate([score, "mpra_activity"]):
            ax = axes[ri, ci]
            for (mode, anchor), z in part.groupby(["design_mode", "anchor"]):
                curve = z.groupby("target_gap_bp")[value].mean().sort_index()
                ax.plot(
                    curve.index,
                    curve - curve.mean(),
                    marker="o",
                    ms=2.5,
                    lw=1,
                    color=colors[mode],
                    alpha=0.9 if anchor == "upstream" else 0.45,
                    label=f"{mode}, anchor {anchor}" if ri == 0 else None,
                )
            ax.axvline(part.native_motif_gap_bp.iloc[0], color="#111", ls="--", lw=0.8)
            ax.set_title(title if ci == 0 else "fine-tuned MPRA activity")
            ax.set_xlabel("edge-to-edge motif gap (bp)")
            ax.set_ylabel("centered model score")
            ax.spines[["top", "right"]].set_visible(False)
    if ranks:
        axes[0, 0].legend(frameon=False, fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(OUT / "Fig_spacing_grammar_targeted.png", dpi=240, bbox_inches="tight")
    fig.savefig(OUT / "Fig_spacing_grammar_targeted.pdf", bbox_inches="tight")
    plt.close(fig)

    if len(family_summary):
        top_families = family_summary.head(8).canonical_pair.tolist()
        fig, axes = plt.subplots(2, 4, figsize=(12, 6), squeeze=False)
        for ax, canonical in zip(axes.ravel(), top_families):
            part = raw[(raw.grid == "coarse_discovery") & (raw.canonical_pair == canonical)]
            curves = []
            for _, z in part.groupby(["Coords", "design_mode", "anchor"]):
                c = centered_curve(z, score)
                if c.std() > 0:
                    c = c / c.std()
                curves.append(c)
                ax.plot(c.index, c, color="#94a3b8", lw=0.6, alpha=0.35)
            meta = pd.concat(curves, axis=1).mean(axis=1)
            ax.plot(meta.index, meta, color="#7c3aed", marker="o", ms=3, lw=1.8)
            ax.set_title(f"{canonical}\n{part.Coords.nunique()} loci")
            ax.set_xlabel("motif gap (bp)")
            ax.set_ylabel("standardized context score")
            ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        fig.savefig(OUT / "Fig_spacing_grammar_discovery.png", dpi=240, bbox_inches="tight")
        fig.savefig(OUT / "Fig_spacing_grammar_discovery.pdf", bbox_inches="tight")
        plt.close(fig)


def main() -> None:
    motifs = parse_meme(MEME)
    fasta = read_fasta(PROJECT / "bundle/MPRA_AP1_fibroblast.fasta")
    pairs_df = pd.read_csv(OUT / "context_module_pairs_with_motifs_uncorrected.csv")

    selected_rows: list[tuple[str, pd.Series, list[int]]] = []
    for _, row in pairs_df.iterrows():
        if int(row.epistasis_rank) in TARGET_RANKS:
            selected_rows.append(("fine_targeted", row, FINE_GAPS))
        good_motifs = (
            float(row.module_1_motif_top_relative_PWM_score) >= RELATIVE_THRESHOLD
            and float(row.module_2_motif_top_relative_PWM_score) >= RELATIVE_THRESHOLD
        )
        # Discovery is restricted to strong pair interactions.  Rescue alone can
        # be high for weak pairs because it also contains additive restoration.
        discovery_priority = (
            int(row.epistasis_rank) <= DISCOVERY_RANK_MAX
            or bool(row.robust_top_decile_pair)
        )
        if good_motifs and discovery_priority:
            selected_rows.append(("coarse_discovery", row, COARSE_GAPS))

    records = []
    baseline_records = {}
    seen_designs = set()
    for grid, row, gaps in selected_rows:
        try:
            pair = build_pair_record(row, motifs)
        except ValueError as exc:
            print(f"skip {row.Coords} rank {row.epistasis_rank}: {exc}", flush=True)
            continue
        seq = core(fasta[("MOTIF_SHUF", pair["Coords"])])
        baseline_records[pair["Coords"]] = seq
        for anchor in ["upstream", "downstream"]:
            for mode in ["motif_core", "module_12bp"]:
                for gap in gaps:
                    design_key = (grid, pair["Coords"], pair["epistasis_rank"], anchor, mode, gap)
                    if design_key in seen_designs:
                        continue
                    seen_designs.add(design_key)
                    try:
                        variant, moved_start, moved_end = variant_for_gap(seq, pair, gap, anchor, mode)
                    except ValueError:
                        continue
                    up = pair["upstream"]
                    down = pair["downstream"]
                    rel1, hits1 = scan_pwm(variant, motifs[up["jaspar"]])
                    rel2, hits2 = scan_pwm(variant, motifs[down["jaspar"]])
                    base_rel1, base_hits1 = scan_pwm(seq, motifs[up["jaspar"]])
                    base_rel2, base_hits2 = scan_pwm(seq, motifs[down["jaspar"]])
                    records.append(
                        {
                            **{k: v for k, v in pair.items() if k not in {"upstream", "downstream"}},
                            "grid": grid,
                            "design_mode": mode,
                            "anchor": anchor,
                            "target_gap_bp": int(gap),
                            "moved_motif_start_0based": moved_start,
                            "moved_motif_end_0based_exclusive": moved_end,
                            "upstream_motif": up["motif"],
                            "downstream_motif": down["motif"],
                            "upstream_native_strand": up["strand"],
                            "downstream_native_strand": down["strand"],
                            "sequence_170bp": variant,
                            "exact_composition_preserved": Counter(variant) == Counter(seq),
                            "upstream_best_relative_PWM": rel1,
                            "downstream_best_relative_PWM": rel2,
                            "upstream_PWM_hit_count": hits1,
                            "downstream_PWM_hit_count": hits2,
                            "upstream_PWM_hit_count_changed": hits1 != base_hits1,
                            "downstream_PWM_hit_count_changed": hits2 != base_hits2,
                            "upstream_baseline_best_relative_PWM": base_rel1,
                            "downstream_baseline_best_relative_PWM": base_rel2,
                        }
                    )
    raw = pd.DataFrame(records)
    if raw.empty:
        raise RuntimeError("No valid spacing variants")
    print(
        f"device={device}; variants={len(raw)}; unique variants={raw.sequence_170bp.nunique()}; "
        f"loci={raw.Coords.nunique()}; pairs={raw[['Coords','epistasis_rank']].drop_duplicates().shape[0]}",
        flush=True,
    )
    raw.to_csv(OUT / "spacing_grammar_variants_unscored.csv", index=False)
    if DRY_RUN:
        return
    if ANALYZE_ONLY:
        scored_path = OUT / "spacing_grammar_model_scores.csv"
        if not scored_path.exists():
            raise FileNotFoundError(scored_path)
        raw = pd.read_csv(scored_path)
        curve_summary, pair_summary, family_summary = summarize(raw)
        curve_summary.to_csv(OUT / "spacing_grammar_curve_summary.csv", index=False)
        pair_summary.to_csv(OUT / "spacing_grammar_pair_summary.csv", index=False)
        family_summary.to_csv(OUT / "spacing_grammar_family_discovery_summary.csv", index=False)
        make_figure(raw, family_summary)
        summary = {
            "analysis_policy": (
                "No sequence-composition correction. Spacing variants preserve the exact 170-bp "
                "mononucleotide composition by block transposition."
            ),
            "n_variant_rows": int(len(raw)),
            "n_unique_sequences": int(raw.sequence_170bp.nunique()),
            "n_loci": int(raw.Coords.nunique()),
            "n_locus_pairs": int(raw[["Coords", "epistasis_rank"]].drop_duplicates().shape[0]),
            "targeted_pair_summaries": pair_summary[pair_summary.grid == "fine_targeted"].to_dict("records"),
            "recurring_family_summaries": family_summary.to_dict("records"),
        }
        with (OUT / "spacing_grammar_summary.json").open("w") as handle:
            json.dump(summary, handle, indent=2)
        print(json.dumps(summary, indent=2), flush=True)
        return

    fold_lookup = fold_map_for_coords(sorted(raw.Coords.unique()))
    unique = raw[["Coords", "sequence_170bp"]].drop_duplicates().reset_index(drop=True)
    unique["outer_fold"] = unique.Coords.map(fold_lookup)
    scores = score_unique_sequences(
        unique.sequence_170bp.tolist(), unique.outer_fold.to_numpy(int)
    )
    for name, values in scores.items():
        unique[name] = values
    raw = raw.merge(unique, on=["Coords", "sequence_170bp"], how="left", validate="many_to_one")

    baseline = pd.DataFrame(
        {"Coords": list(baseline_records), "sequence_170bp": list(baseline_records.values())}
    ).drop_duplicates("Coords")
    baseline["outer_fold"] = baseline.Coords.map(fold_lookup)
    missing = baseline.merge(
        unique, on=["Coords", "sequence_170bp"], how="left", suffixes=("", "_scored")
    )
    score_names = ["mpra_activity", "total_Endo_minus_Full", "AP1_Endo_minus_Motif", "context_Motif_minus_Full"]
    need = missing[missing[score_names[0]].isna()][["Coords", "sequence_170bp", "outer_fold"]]
    if len(need):
        extra = score_unique_sequences(need.sequence_170bp.tolist(), need.outer_fold.to_numpy(int))
        for name, values in extra.items():
            need[name] = values
        baseline_scored = pd.concat(
            [missing[missing[score_names[0]].notna()][["Coords"] + score_names], need[["Coords"] + score_names]],
            ignore_index=True,
        )
    else:
        baseline_scored = missing[["Coords"] + score_names]
    baseline_scored = baseline_scored.drop_duplicates("Coords")
    baseline_scored = baseline_scored.rename(
        columns={name: f"native_{name}" for name in score_names}
    )
    raw = raw.merge(baseline_scored, on="Coords", how="left", validate="many_to_one")
    raw.to_csv(OUT / "spacing_grammar_model_scores.csv", index=False)

    curve_summary, pair_summary, family_summary = summarize(raw)
    curve_summary.to_csv(OUT / "spacing_grammar_curve_summary.csv", index=False)
    pair_summary.to_csv(OUT / "spacing_grammar_pair_summary.csv", index=False)
    family_summary.to_csv(OUT / "spacing_grammar_family_discovery_summary.csv", index=False)
    make_figure(raw, family_summary)

    summary = {
        "analysis_policy": (
            "No sequence-composition correction. Spacing variants preserve the exact 170-bp "
            "mononucleotide composition by block transposition."
        ),
        "n_variant_rows": int(len(raw)),
        "n_unique_sequences": int(raw.sequence_170bp.nunique()),
        "n_loci": int(raw.Coords.nunique()),
        "n_locus_pairs": int(raw[["Coords", "epistasis_rank"]].drop_duplicates().shape[0]),
        "targeted_pair_summaries": pair_summary[pair_summary.grid == "fine_targeted"].to_dict("records"),
        "recurring_family_summaries": family_summary.to_dict("records"),
    }
    with (OUT / "spacing_grammar_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
