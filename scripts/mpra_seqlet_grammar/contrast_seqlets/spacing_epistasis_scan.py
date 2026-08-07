#!/usr/bin/env python3
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Run targeted factorial spacing scans for nominated motif pairs.
# Steps: 1) Construct moved/disrupted motif states; 2) score context and activity readouts; 3) calculate gap-specific epistasis.
# =======================================
"""Four-state spacing scan that separates pair grammar from position effects.

At every gap this script scores both motifs intact, each single-site disruption,
and the double disruption.  Motif disruptions are within-footprint nucleotide
permutations, so every construct retains the exact 170-bp base composition.
The factorial interaction is

    epsilon(gap) = both - anchor_only - moved_only + neither.

Positive epsilon indicates synergistic activity of the intact pair at that gap.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
import hashlib
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from config.paths import publication_path


TASK = Path(publication_path("MPRA_SEQLET_TASK_ROOT"))
OUT = TASK / "outputs"
sys.path.insert(0, str(Path(__file__).resolve().parent))
import spacing_grammar_sweep as sg  # noqa: E402


BATCH = int(os.environ.get("BATCH", "64"))
ANALYZE_ONLY = bool(int(os.environ.get("ANALYZE_ONLY", "0")))
REUSE_EXISTING = bool(int(os.environ.get("REUSE_EXISTING", "0")))
TARGET_RANKS = {1, 2, 4, 6, 15}
DISCOVERY_RANKS = {3, 5, 11, 12, 39, 42}
SELECTED_RANKS = TARGET_RANKS | DISCOVERY_RANKS
FINE_GAPS = list(range(0, 31)) + [32, 34, 36, 40, 44, 48, 56, 64]
COARSE_GAPS = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 28, 32, 36, 40, 48, 56, 64]
STATES = ["both_intact", "anchor_only", "moved_only", "neither"]


def weakest_permutation(seq: str, pwm: np.ndarray, jaspar: str) -> tuple[str, float]:
    """Find a low-scoring exact-composition permutation of a motif footprint."""
    digest = hashlib.sha256(f"{jaspar}|{seq}".encode()).digest()
    rng = np.random.default_rng(int.from_bytes(digest[:8], "little"))
    best_seq = seq
    best_score = sg.scan_pwm(seq, pwm)[0]
    seen = {seq}
    # Widths are 4–11 bp. Five thousand draws are ample for short repeated motifs
    # and give a reproducible low-affinity permutation for longer footprints.
    for _ in range(5000):
        candidate = "".join(rng.permutation(list(seq)))
        if candidate in seen:
            continue
        seen.add(candidate)
        score = sg.scan_pwm(candidate, pwm)[0]
        if score < best_score:
            best_seq, best_score = candidate, score
    if best_seq == seq:
        raise RuntimeError(f"Could not disrupt {jaspar} footprint {seq}")
    assert Counter(best_seq) == Counter(seq)
    return best_seq, float(best_score)


def replace(seq: str, start: int, end: int, value: str) -> str:
    assert len(value) == end - start
    result = seq[:start] + value + seq[end:]
    assert len(result) == len(seq)
    return result


def build_designs() -> pd.DataFrame:
    motifs = sg.parse_meme(sg.MEME)
    fasta = sg.read_fasta(sg.PROJECT / "bundle/MPRA_AP1_fibroblast.fasta")
    pairs = pd.read_csv(OUT / "context_module_pairs_with_motifs_uncorrected.csv")
    pairs = pairs[pairs.epistasis_rank.isin(SELECTED_RANKS)].copy()
    records = []
    for _, row in pairs.iterrows():
        pair = sg.build_pair_record(row, motifs)
        native = sg.core(fasta[("MOTIF_SHUF", pair["Coords"])])
        up, down = pair["upstream"], pair["downstream"]
        up_seq = native[up["motif_start"] : up["motif_end"]]
        down_seq = native[down["motif_start"] : down["motif_end"]]
        up_disrupted, up_disrupted_rel = weakest_permutation(up_seq, motifs[up["jaspar"]], up["jaspar"])
        down_disrupted, down_disrupted_rel = weakest_permutation(
            down_seq, motifs[down["jaspar"]], down["jaspar"]
        )
        grid = FINE_GAPS if int(row.epistasis_rank) in TARGET_RANKS else COARSE_GAPS
        # Always include the exact native motif-footprint gap, including long
        # native arrangements outside the standard discovery grid.
        gaps = sorted(set(grid + [int(pair["native_motif_gap_bp"])]))
        for mode in ["motif_core", "module_12bp"]:
            for gap in gaps:
                try:
                    intact, moved_start, moved_end = sg.variant_for_gap(
                        native, pair, gap, "upstream", mode
                    )
                except ValueError:
                    continue
                anchor_start, anchor_end = up["motif_start"], up["motif_end"]
                moved_original = intact[moved_start:moved_end]
                assert moved_original == down_seq
                variants = {
                    "both_intact": intact,
                    "anchor_only": replace(intact, moved_start, moved_end, down_disrupted),
                    "moved_only": replace(intact, anchor_start, anchor_end, up_disrupted),
                    "neither": replace(
                        replace(intact, anchor_start, anchor_end, up_disrupted),
                        moved_start,
                        moved_end,
                        down_disrupted,
                    ),
                }
                for state, sequence in variants.items():
                    assert Counter(sequence) == Counter(native)
                    records.append(
                        {
                            "Coords": pair["Coords"],
                            "resolved_class": pair["resolved_class"],
                            "epistasis_rank": pair["epistasis_rank"],
                            "motif_pair_name": pair["motif_pair_name"],
                            "canonical_pair": pair["canonical_pair"],
                            "mean_absolute_epsilon_original": pair["mean_absolute_epsilon"],
                            "rescue_high_gap_closed": pair["rescue_high_gap_closed"],
                            "native_motif_gap_bp": pair["native_motif_gap_bp"],
                            "upstream_motif": up["motif"],
                            "downstream_motif": down["motif"],
                            "upstream_native_strand": up["strand"],
                            "downstream_native_strand": down["strand"],
                            "upstream_native_footprint": up_seq,
                            "downstream_native_footprint": down_seq,
                            "upstream_disrupted_footprint": up_disrupted,
                            "downstream_disrupted_footprint": down_disrupted,
                            "upstream_disrupted_relative_PWM": up_disrupted_rel,
                            "downstream_disrupted_relative_PWM": down_disrupted_rel,
                            "design_mode": mode,
                            "anchor": "upstream",
                            "target_gap_bp": int(gap),
                            "state": state,
                            "sequence_170bp": sequence,
                            "exact_composition_preserved": True,
                        }
                    )
    return pd.DataFrame(records)


def calculate_epistasis(raw: pd.DataFrame) -> pd.DataFrame:
    values = [
        "context_Motif_minus_Full",
        "mpra_activity",
        "total_Endo_minus_Full",
        "AP1_Endo_minus_Motif",
    ]
    index = [
        "Coords",
        "resolved_class",
        "epistasis_rank",
        "motif_pair_name",
        "canonical_pair",
        "mean_absolute_epsilon_original",
        "rescue_high_gap_closed",
        "native_motif_gap_bp",
        "design_mode",
        "target_gap_bp",
    ]
    output = None
    for value in values:
        wide = raw.pivot_table(index=index, columns="state", values=value).reset_index()
        wide[f"epsilon_{value}"] = (
            wide["both_intact"] - wide["anchor_only"] - wide["moved_only"] + wide["neither"]
        )
        # This matches the established pair-disruption convention used in the
        # project: negative factorial epsilon means the double hit loses more
        # activity than additive expectation, so loss-synergy is -epsilon.
        wide[f"loss_synergy_{value}"] = -wide[f"epsilon_{value}"]
        wide[f"anchor_effect_{value}"] = wide["both_intact"] - wide["moved_only"]
        wide[f"moved_effect_{value}"] = wide["both_intact"] - wide["anchor_only"]
        keep = index + [
            f"epsilon_{value}",
            f"loss_synergy_{value}",
            f"anchor_effect_{value}",
            f"moved_effect_{value}",
        ]
        output = wide[keep] if output is None else output.merge(wide[keep], on=index, how="inner")
    return output


def correlation(a: pd.Series, b: pd.Series) -> float:
    common = a.dropna().index.intersection(b.dropna().index)
    if len(common) < 6 or a.loc[common].std() == 0 or b.loc[common].std() == 0:
        return np.nan
    return float(np.corrcoef(a.loc[common], b.loc[common])[0, 1])


def summarize(epi: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    ec = "epsilon_context_Motif_minus_Full"
    ea = "epsilon_mpra_activity"
    pair_rows = []
    meta_curves = {}
    for (coord, rank), part in epi.groupby(["Coords", "epistasis_rank"]):
        context_curves = {}
        activity_curves = {}
        for mode, z in part.groupby("design_mode"):
            context_curves[mode] = z.set_index("target_gap_bp")[ec].sort_index()
            activity_curves[mode] = z.set_index("target_gap_bp")[ea].sort_index()
        common_context = pd.concat(context_curves, axis=1).mean(axis=1)
        common_activity = pd.concat(activity_curves, axis=1).mean(axis=1)
        native_gap = int(part.native_motif_gap_bp.iloc[0])
        native_context = (
            float(common_context.loc[native_gap]) if native_gap in common_context.index else np.nan
        )
        native_activity = (
            float(common_activity.loc[native_gap]) if native_gap in common_activity.index else np.nan
        )
        pair_rows.append(
            {
                "Coords": coord,
                "epistasis_rank": int(rank),
                "resolved_class": str(part.resolved_class.iloc[0]),
                "motif_pair_name": str(part.motif_pair_name.iloc[0]),
                "canonical_pair": str(part.canonical_pair.iloc[0]),
                "native_motif_gap_bp": native_gap,
                "original_pair_absolute_epsilon": float(part.mean_absolute_epsilon_original.iloc[0]),
                "original_pair_rescue": float(part.rescue_high_gap_closed.iloc[0]),
                "n_gaps": int(len(common_context)),
                "core_module_context_curve_r": correlation(
                    context_curves.get("motif_core", pd.Series(dtype=float)),
                    context_curves.get("module_12bp", pd.Series(dtype=float)),
                ),
                "core_module_activity_curve_r": correlation(
                    activity_curves.get("motif_core", pd.Series(dtype=float)),
                    activity_curves.get("module_12bp", pd.Series(dtype=float)),
                ),
                "context_activity_epistasis_curve_r": correlation(common_context, common_activity),
                "best_loss_synergy_context_gap_bp": int(common_context.idxmin()),
                "best_loss_synergy_activity_gap_bp": int(common_activity.idxmin()),
                "best_redundancy_context_gap_bp": int(common_context.idxmax()),
                "best_redundancy_activity_gap_bp": int(common_activity.idxmax()),
                "strongest_context_interaction_gap_bp": int(common_context.abs().idxmax()),
                "strongest_activity_interaction_gap_bp": int(common_activity.abs().idxmax()),
                "max_context_loss_synergy": float(-common_context.min()),
                "max_activity_loss_synergy": float(-common_activity.min()),
                "context_epsilon_range": float(common_context.max() - common_context.min()),
                "activity_epsilon_range": float(common_activity.max() - common_activity.min()),
                "native_context_epsilon": native_context,
                "native_activity_epsilon": native_activity,
                "native_context_loss_synergy_percentile": (
                    float(np.mean(common_context >= native_context)) if np.isfinite(native_context) else np.nan
                ),
                "native_activity_loss_synergy_percentile": (
                    float(np.mean(common_activity >= native_activity)) if np.isfinite(native_activity) else np.nan
                ),
            }
        )
        meta_curves[(coord, int(rank), "context")] = common_context
        meta_curves[(coord, int(rank), "activity")] = common_activity
    pair_summary = pd.DataFrame(pair_rows).sort_values("epistasis_rank")

    family_rows = []
    for canonical, members in pair_summary.groupby("canonical_pair"):
        if len(members) < 2:
            continue
        context_curves = []
        activity_curves = []
        for r in members.itertuples(index=False):
            c = meta_curves[(r.Coords, r.epistasis_rank, "context")]
            a = meta_curves[(r.Coords, r.epistasis_rank, "activity")]
            if c.std() > 0:
                c = c / c.std()
            if a.std() > 0:
                a = a / a.std()
            context_curves.append(c.rename(r.Coords))
            activity_curves.append(a.rename(r.Coords))
        context_matrix = pd.concat(context_curves, axis=1)
        activity_matrix = pd.concat(activity_curves, axis=1)
        meta_context = context_matrix.mean(axis=1)
        meta_activity = activity_matrix.mean(axis=1)
        locus_context_r = []
        locus_activity_r = []
        for a, b in __import__("itertools").combinations(range(len(context_curves)), 2):
            locus_context_r.append(correlation(context_curves[a], context_curves[b]))
            locus_activity_r.append(correlation(activity_curves[a], activity_curves[b]))
        loo_context_best = []
        loo_activity_best = []
        for coord in context_matrix.columns:
            other_c = context_matrix.drop(columns=coord)
            other_a = activity_matrix.drop(columns=coord)
            if other_c.shape[1]:
                loo_context_best.append(int(other_c.mean(axis=1).idxmin()))
                loo_activity_best.append(int(other_a.mean(axis=1).idxmin()))
        family_rows.append(
            {
                "canonical_pair": canonical,
                "n_loci": int(len(members)),
                "locus_context_curve_r": float(np.nanmean(locus_context_r)),
                "locus_activity_curve_r": float(np.nanmean(locus_activity_r)),
                "meta_context_activity_r": correlation(meta_context, meta_activity),
                "best_loss_synergy_context_gap_bp": int(meta_context.idxmin()),
                "best_loss_synergy_activity_gap_bp": int(meta_activity.idxmin()),
                "leave_one_locus_out_context_best_gaps": ";".join(map(str, loo_context_best)),
                "leave_one_locus_out_activity_best_gaps": ";".join(map(str, loo_activity_best)),
                "mean_original_pair_absolute_epsilon": float(members.original_pair_absolute_epsilon.mean()),
                "mean_original_pair_rescue": float(members.original_pair_rescue.mean()),
            }
        )
    family_summary = pd.DataFrame(family_rows)
    return pair_summary, family_summary


def make_figure(epi: pd.DataFrame, pair_summary: pd.DataFrame) -> None:
    ec = "epsilon_context_Motif_minus_Full"
    ea = "epsilon_mpra_activity"
    order = pair_summary.epistasis_rank.tolist()
    fig, axes = plt.subplots(len(order), 2, figsize=(10.5, 2.15 * len(order)), squeeze=False)
    colors = {"motif_core": "#2563eb", "module_12bp": "#d97706"}
    for ri, rank in enumerate(order):
        part = epi[epi.epistasis_rank == rank]
        native_gap = int(part.native_motif_gap_bp.iloc[0])
        for ci, value in enumerate([ec, ea]):
            ax = axes[ri, ci]
            for mode, z in part.groupby("design_mode"):
                z = z.sort_values("target_gap_bp")
                ax.plot(
                    z.target_gap_bp,
                    z[value],
                    marker="o",
                    ms=2.5,
                    lw=1,
                    color=colors[mode],
                    label=mode if ri == 0 else None,
                )
            ax.axhline(0, color="#64748b", lw=0.7)
            ax.axvline(native_gap, color="#111827", ls="--", lw=0.8)
            ax.set_xlabel("edge-to-edge motif gap (bp)")
            ax.set_ylabel("factorial ε (negative = loss synergy)")
            if ci == 0:
                ax.set_title(
                    f"rank {rank}: {part.motif_pair_name.iloc[0]} ({part.resolved_class.iloc[0]})"
                )
            else:
                ax.set_title("fine-tuned MPRA activity")
            ax.spines[["top", "right"]].set_visible(False)
    axes[0, 0].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "Fig_spacing_epistasis_scan.png", dpi=240, bbox_inches="tight")
    fig.savefig(OUT / "Fig_spacing_epistasis_scan.pdf", bbox_inches="tight")
    plt.close(fig)


def finish(raw: pd.DataFrame) -> None:
    epi = calculate_epistasis(raw)
    pair_summary, family_summary = summarize(epi)
    epi.to_csv(OUT / "spacing_epistasis_by_gap.csv", index=False)
    pair_summary.to_csv(OUT / "spacing_epistasis_pair_summary.csv", index=False)
    family_summary.to_csv(OUT / "spacing_epistasis_family_summary.csv", index=False)
    make_figure(epi, pair_summary)
    summary = {
        "analysis_policy": (
            "No sequence-composition correction. All spacing and motif-disruption constructs "
            "preserve exact 170-bp mononucleotide composition."
        ),
        "factorial_interaction": "both - anchor_only - moved_only + neither",
        "n_construct_rows": int(len(raw)),
        "n_unique_sequences": int(raw.sequence_170bp.nunique()),
        "pair_summaries": pair_summary.to_dict("records"),
        "recurring_pair_summaries": family_summary.to_dict("records"),
    }
    with (OUT / "spacing_epistasis_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
    print(json.dumps(summary, indent=2), flush=True)


def main() -> None:
    if ANALYZE_ONLY:
        finish(pd.read_csv(OUT / "spacing_epistasis_model_scores.csv"))
        return
    design = build_designs()
    design.to_csv(OUT / "spacing_epistasis_constructs.csv", index=False)
    unique = design[["Coords", "sequence_170bp"]].drop_duplicates().reset_index(drop=True)
    fold_lookup = sg.fold_map_for_coords(sorted(unique.Coords.unique()))
    unique["outer_fold"] = unique.Coords.map(fold_lookup)
    score_names = [
        "total_Endo_minus_Full",
        "AP1_Endo_minus_Motif",
        "context_Motif_minus_Full",
        "mpra_activity",
    ]
    if REUSE_EXISTING and (OUT / "spacing_epistasis_model_scores.csv").exists():
        old = pd.read_csv(OUT / "spacing_epistasis_model_scores.csv")
        old = old[["Coords", "sequence_170bp"] + score_names].drop_duplicates(
            ["Coords", "sequence_170bp"]
        )
        unique = unique.merge(old, on=["Coords", "sequence_170bp"], how="left")
        missing = unique[unique[score_names[0]].isna()].copy()
        print(f"reusing {len(unique)-len(missing)} scores; scoring {len(missing)} new sequences", flush=True)
        if len(missing):
            scores = sg.score_unique_sequences(
                missing.sequence_170bp.tolist(), missing.outer_fold.to_numpy(int)
            )
            for name, values in scores.items():
                missing[name] = values
            unique = pd.concat(
                [unique[unique[score_names[0]].notna()], missing], ignore_index=True
            )
    else:
        scores = sg.score_unique_sequences(
            unique.sequence_170bp.tolist(), unique.outer_fold.to_numpy(int)
        )
        for name, values in scores.items():
            unique[name] = values
    raw = design.merge(unique, on=["Coords", "sequence_170bp"], how="left", validate="many_to_one")
    raw.to_csv(OUT / "spacing_epistasis_model_scores.csv", index=False)
    finish(raw)


if __name__ == "__main__":
    main()
