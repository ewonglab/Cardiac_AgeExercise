#!/usr/bin/env python3
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Consolidate validated context modules and annotate their motifs.
# Steps: 1) Load perturbation outputs and sequence inputs; 2) calculate locus-level summaries; 3) write motif and module tables.
# =======================================
"""Locus-level statistics and motif annotation for validated context modules."""

from __future__ import annotations

import json
import re
from collections import Counter
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from config.paths import publication_path

import numpy as np
import pandas as pd
from scipy import stats


PROJECT = Path(publication_path("MPRA_KOO_PROJECT"))
TASK = Path(publication_path("MPRA_SEQLET_TASK_ROOT"))
OUT = TASK / "outputs"


def read_fasta(path: Path) -> dict[tuple[str, str], str]:
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


def core(seq: str) -> str:
    return seq[15:-15].upper()


def parse_meme(path: Path) -> dict[str, np.ndarray]:
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
            matrix = []
            for row in lines[i + 1 : i + 1 + width]:
                matrix.append([float(x) for x in row.split()[:4]])
            motifs[name] = np.asarray(matrix, float)
            i += width
        i += 1
    return motifs


def pwm_best(seq: str, pwm: np.ndarray) -> float:
    bases = {b: i for i, b in enumerate("ACGT")}
    logodds = np.log2((pwm + 1e-5) / 0.25)
    rc = logodds[::-1, ::-1]
    best = -np.inf
    for mat in [logodds, rc]:
        w = len(mat)
        for start in range(len(seq) - w + 1):
            score = sum(mat[j, bases[seq[start + j]]] for j in range(w))
            best = max(best, score)
    return float(best)


def bh(pvalues: list[float]) -> list[float]:
    p = np.asarray(pvalues, float)
    order = np.argsort(p)
    q = np.empty_like(p)
    running = 1.0
    for rank in range(len(p) - 1, -1, -1):
        idx = order[rank]
        running = min(running, p[idx] * len(p) / (rank + 1))
        q[idx] = running
    return q.tolist()


def partial_class_f(y: np.ndarray, covariates: np.ndarray, classes: np.ndarray) -> dict[str, float]:
    # Rank-transform the outcome; compare covariates-only with covariates + class dummies.
    yr = stats.rankdata(y) / (len(y) + 1)
    levels = ["Age-only", "Exercise-only", "Reversal"]
    dummy = np.column_stack([(classes == level).astype(float) for level in levels])
    reduced = np.column_stack([np.ones(len(y)), covariates])
    full = np.column_stack([reduced, dummy])
    br, *_ = np.linalg.lstsq(reduced, yr, rcond=None)
    bf, *_ = np.linalg.lstsq(full, yr, rcond=None)
    sse_r = np.sum((yr - reduced @ br) ** 2)
    sse_f = np.sum((yr - full @ bf) ** 2)
    q = full.shape[1] - reduced.shape[1]
    df = len(y) - full.shape[1]
    f = ((sse_r - sse_f) / q) / (sse_f / df)
    return {"F": float(f), "df1": int(q), "df2": int(df), "p": float(stats.f.sf(f, q, df))}


def main() -> None:
    source = (
        pd.read_csv(PROJECT / "bundle/MPRA_Results_with_TSS_distance.csv")
        .drop_duplicates("Coords")
        .set_index("Coords")
    )
    nested = pd.read_csv(OUT / "mpra_nested_contrast_predictions.csv").set_index("Coords")
    attr = pd.read_csv(OUT / "context_probe_attribution_summary.csv").set_index("Coords")
    disruption = pd.read_csv(OUT / "context_module_disruption_summary.csv")
    epi = pd.read_csv(OUT / "context_module_epistasis.csv")
    rescue = pd.read_csv(OUT / "context_module_rescue.csv")
    for frame in [source, nested, attr]:
        frame.index = frame.index.astype(str)
    for frame in [disruption, epi, rescue]:
        frame["Coords"] = frame["Coords"].astype(str)

    # One independent row per enhancer.
    direct = disruption.groupby(["Coords", "resolved_class", "module_type"])[
        "mean_absolute_effect"
    ].mean().unstack("module_type")
    epiloc = epi.groupby(["Coords", "resolved_class", "module_type"])[
        ["absolute_epsilon", "loss_interaction"]
    ].mean().unstack("module_type")
    epiloc.columns = [f"{a}_{b}" for a, b in epiloc.columns]
    triple = rescue[rescue["combo"] == "ABC"].pivot_table(
        index=["Coords", "resolved_class"], columns="module_type", values="absolute_gap_closed"
    )
    direct.columns = [f"direct_{c}" for c in direct.columns]
    triple.columns = [f"rescue_{c}" for c in triple.columns]
    locus = direct.join(epiloc).join(triple).reset_index()
    locus["measured_context"] = locus["Coords"].map(
        nested["measured__context_Motif_minus_Full"]
    )
    locus["predicted_context_oof"] = locus["Coords"].map(
        nested["nested_chr_oof__context_Motif_minus_Full"]
    )
    locus["attribution_gini"] = locus["Coords"].map(attr["attribution_gini_motif"])
    locus["GC_content"] = locus["Coords"].map(source["GC_content"])
    locus["Endo_activity"] = locus["Coords"].map(source["log2_alpha_Endo"])
    locus.to_csv(OUT / "context_module_locus_level_results.csv", index=False)

    outcomes = {
        "direct_high_effect": "direct_high",
        "high_abs_epistasis": "absolute_epsilon_high",
        "high_triple_rescue": "rescue_high",
    }
    class_report = {}
    cov = np.column_stack(
        [
            np.abs(locus["measured_context"].to_numpy(float)),
            locus["GC_content"].to_numpy(float),
            locus["Endo_activity"].to_numpy(float),
        ]
    )
    for label, column in outcomes.items():
        ok = locus[column].notna().to_numpy()
        arrays, rows = [], {}
        for cls in ["Stable", "Age-only", "Exercise-only", "Reversal"]:
            values = locus.loc[ok & (locus["resolved_class"] == cls), column].to_numpy(float)
            arrays.append(values)
            rows[cls] = {
                "n_loci": int(len(values)),
                "mean": float(values.mean()),
                "median": float(np.median(values)),
            }
        rows["kruskal_p_locus_level"] = float(stats.kruskal(*arrays).pvalue)
        rows["class_partial_F_controlling_abs_effect_GC_Endo"] = partial_class_f(
            locus.loc[ok, column].to_numpy(float),
            cov[ok],
            locus.loc[ok, "resolved_class"].to_numpy(),
        )
        class_report[label] = rows

    # PWM-family annotation of high versus low modules, using 22-bp local contexts.
    modules = pd.read_csv(OUT / "context_probe_candidate_modules.csv")
    modules["Coords"] = modules["Coords"].astype(str)
    fa = read_fasta(PROJECT / "bundle/MPRA_AP1_fibroblast.fasta")
    motifs = parse_meme(PROJECT / "koorun/out/jaspar_panel.meme")
    family_members = {
        "AP1": ["FOS::JUN", "JUN_", "JUNB", "FOSL1", "FOSL2", "BATF", "JDP2"],
        "RUNX": ["RUNX"],
        "ETS": ["ETS1", "ELK", "GABPA", "FLI1", "ETV"],
        "CEBP": ["CEBP"],
        "NFAT": ["NFAT"],
        "GATA": ["GATA"],
        "MEF2": ["MEF2"],
        "TEAD": ["TEAD"],
        "TBX": ["TBX"],
        "NKX2-5": ["NKX2-5"],
        "SRF": ["SRF"],
        "ESRRB": ["ESRRB"],
        "SP/KLF/EGR": ["SP1", "KLF", "EGR1"],
        "TCF7L2": ["TCF7L2"],
        "CTCF": ["CTCF"],
    }
    pwm_rows = []
    for row in modules.itertuples(index=False):
        seq = core(fa[("MOTIF_SHUF", str(row.Coords))])
        center = (int(row.start_0based) + int(row.end_0based_exclusive)) // 2
        lo, hi = max(0, center - 11), min(170, center + 11)
        local = seq[lo:hi]
        scores = {name: pwm_best(local, pwm) for name, pwm in motifs.items() if len(pwm) <= len(local)}
        record = {
            "Coords": str(row.Coords),
            "module_type": row.module_type,
            "module_rank": int(row.module_rank),
            "start_0based": int(row.start_0based),
            "end_0based_exclusive": int(row.end_0based_exclusive),
            "local_sequence": local,
        }
        for family, prefixes in family_members.items():
            member = [v for name, v in scores.items() if any(name.startswith(p) for p in prefixes)]
            record[family] = max(member) if member else np.nan
        pwm_rows.append(record)
    pwm = pd.DataFrame(pwm_rows)
    pwm.to_csv(OUT / "context_module_pwm_scores.csv", index=False)

    motif_results = []
    for family in family_members:
        by = pwm.groupby(["Coords", "module_type"])[family].mean().unstack("module_type").dropna()
        test = stats.wilcoxon(by["high"], by["low_control"], alternative="greater")
        motif_results.append(
            {
                "family": family,
                "n_loci": int(len(by)),
                "high_mean": float(by["high"].mean()),
                "low_mean": float(by["low_control"].mean()),
                "paired_wilcoxon_greater_p": float(test.pvalue),
            }
        )
    q = bh([r["paired_wilcoxon_greater_p"] for r in motif_results])
    for row, value in zip(motif_results, q):
        row["BH_q"] = value
    motif_results.sort(key=lambda x: x["BH_q"])

    # Composition control: residualize every PWM-family score on local mono- and
    # dinucleotide frequencies, then repeat the within-locus high-vs-low test.
    feature_rows = []
    for seq in pwm["local_sequence"]:
        n = len(seq)
        feat = [n]
        feat += [seq.count(b) / n for b in "ACGT"]
        feat += [seq.count(a + b) / max(1, n - 1) for a in "ACGT" for b in "ACGT"]
        feature_rows.append(feat)
    xcomp = np.column_stack([np.ones(len(pwm)), np.asarray(feature_rows, float)])
    controlled = []
    for family in family_members:
        y = pwm[family].to_numpy(float)
        coef, *_ = np.linalg.lstsq(xcomp, y, rcond=None)
        residual = y - xcomp @ coef
        pwm[f"{family}_composition_residual"] = residual
        by = (
            pwm.assign(_residual=residual)
            .groupby(["Coords", "module_type"])["_residual"]
            .mean()
            .unstack("module_type")
            .dropna()
        )
        test = stats.wilcoxon(by["high"], by["low_control"], alternative="greater")
        controlled.append(
            {
                "family": family,
                "n_loci": int(len(by)),
                "high_residual_mean": float(by["high"].mean()),
                "low_residual_mean": float(by["low_control"].mean()),
                "paired_wilcoxon_greater_p": float(test.pvalue),
            }
        )
    q_controlled = bh([r["paired_wilcoxon_greater_p"] for r in controlled])
    for row, value in zip(controlled, q_controlled):
        row["BH_q"] = value
    controlled.sort(key=lambda x: x["BH_q"])
    pwm.to_csv(OUT / "context_module_pwm_scores.csv", index=False)
    raw_df = pd.DataFrame(motif_results)
    raw_df["analysis"] = "raw_PWM"
    ctrl_df = pd.DataFrame(controlled)
    ctrl_df["analysis"] = "mono_dinucleotide_residual"
    pd.concat([raw_df, ctrl_df], ignore_index=True, sort=False).to_csv(
        OUT / "context_module_motif_enrichment.csv", index=False
    )

    report = {
        "n_loci": int(len(locus)),
        "class_tests_independent_locus_unit": class_report,
        "motif_family_enrichment_high_vs_low_raw": motif_results,
        "motif_family_enrichment_high_vs_low_composition_controlled": controlled,
    }
    with (OUT / "context_module_consolidated_summary.json").open("w") as handle:
        json.dump(report, handle, indent=2)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
