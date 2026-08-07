#!/usr/bin/env python3
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Perform robustness diagnostics and select contrast-model follow-up loci.
# Steps: 1) Load benchmark outputs; 2) calculate cross-validation and bootstrap metrics; 3) write selected-locus summaries.
# =======================================
"""Robustness, subclass diagnostics, and locus selection for MPRA contrast models."""

from __future__ import annotations

import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from config.paths import publication_path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


PROJECT = Path(publication_path("MPRA_KOO_PROJECT"))
TASK = Path(publication_path("MPRA_SEQLET_TASK_ROOT"))
OUT = TASK / "outputs"
ALPHAS = [300.0, 1000.0, 3000.0, 10000.0, 30000.0]


def metrics(y: np.ndarray, p: np.ndarray) -> dict[str, float | int]:
    ok = np.isfinite(y) & np.isfinite(p)
    y, p = y[ok], p[ok]
    return {
        "n": int(len(y)),
        "pearson": float(stats.pearsonr(y, p).statistic),
        "spearman": float(stats.spearmanr(y, p).statistic),
        "r2": float(r2_score(y, p)),
        "rmse": float(mean_squared_error(y, p) ** 0.5),
        "sign_accuracy": float(np.mean(np.sign(y) == np.sign(p))),
    }


def cv_ridge(x: np.ndarray, y: np.ndarray, chrom: np.ndarray, alpha: float) -> np.ndarray:
    pred = np.full(len(y), np.nan)
    for tr, te in GroupKFold(5).split(x, y, groups=chrom):
        model = make_pipeline(StandardScaler(), Ridge(alpha=alpha))
        model.fit(x[tr], y[tr])
        pred[te] = model.predict(x[te])
    return pred


def block_bootstrap_ci(
    y: np.ndarray, p: np.ndarray, chrom: np.ndarray, nboot: int = 3000
) -> dict[str, list[float]]:
    rng = np.random.default_rng(20260722)
    levels = np.unique(chrom)
    rs, r2s = [], []
    by = {c: np.flatnonzero(chrom == c) for c in levels}
    for _ in range(nboot):
        chosen = rng.choice(levels, len(levels), replace=True)
        idx = np.concatenate([by[c] for c in chosen])
        rs.append(stats.pearsonr(y[idx], p[idx]).statistic)
        r2s.append(r2_score(y[idx], p[idx]))
    return {
        "pearson_95ci": [float(x) for x in np.percentile(rs, [2.5, 97.5])],
        "r2_95ci": [float(x) for x in np.percentile(r2s, [2.5, 97.5])],
    }


def zscore(x: pd.Series) -> pd.Series:
    sd = x.std(ddof=0)
    return (x - x.mean()) / (sd if sd > 0 else 1.0)


def main() -> None:
    pred = pd.read_csv(OUT / "mpra_contrast_benchmark_predictions.csv")
    pred = pred.pivot_table(
        index=["Coords", "Group", "chr", "training_fold"],
        columns="contrast",
        values=["measured", "existing_probe_oof", "retrained_delta_probe_chr_oof"],
    )
    pred.columns = [f"{a}__{b}" for a, b in pred.columns]
    pred = pred.reset_index()

    raw = pd.read_csv(PROJECT / "bundle/MPRA_Results_with_TSS_distance.csv")
    unique = raw.drop_duplicates("Coords").set_index("Coords")
    pred["replicated"] = pred["Coords"].map(unique["replicated"]).astype(bool)

    # Resolve overlapping age/exercise calls rather than assigning the first duplicate row.
    coord_groups = raw.groupby("Coords")["Group"].agg(lambda x: set(x))

    def resolved(gs: set[str]) -> str:
        if "Stable" in gs:
            return "Stable"
        age = any(g.startswith("Age") for g in gs)
        ex = any(g.startswith("Ex") for g in gs)
        if age and ex:
            return "Reversal"
        return "Age-only" if age else "Exercise-only"

    pred["resolved_class"] = pred["Coords"].map(coord_groups.map(resolved))
    expected = {"Stable": 457, "Age-only": 523, "Exercise-only": 314, "Reversal": 359}
    observed = pred["resolved_class"].value_counts().to_dict()
    if observed != expected:
        raise RuntimeError(f"Resolved class counts differ: {observed}")

    cache = np.load(OUT / "mpra_triplet_embeddings_motif_full.npz")
    cache_coords = cache["coords"].astype(str)
    source = (
        pd.read_csv(PROJECT / "bundle/MPRA_Results_with_TSS_distance.csv")
        .drop_duplicates("Coords")
        .reset_index(drop=True)
    )
    eorig_all = np.load(PROJECT / "koorun/out/embeddings.npy")
    eorig_map = {c: eorig_all[i] for i, c in enumerate(source["Coords"])}
    emap = {c: i for i, c in enumerate(cache_coords)}
    eorig = np.vstack([eorig_map[c] for c in pred["Coords"]]).astype(np.float32)
    emotif = np.vstack([cache["motif"][emap[c]] for c in pred["Coords"]]).astype(np.float32)
    efull = np.vstack([cache["full"][emap[c]] for c in pred["Coords"]]).astype(np.float32)
    emb = {"ORIG": eorig, "MOTIF_SHUF": emotif, "SHUF_FULL": efull}

    definitions = {
        "total_Endo_minus_Full": ("ORIG", "SHUF_FULL"),
        "AP1_Endo_minus_Motif": ("ORIG", "MOTIF_SHUF"),
        "context_Motif_minus_Full": ("MOTIF_SHUF", "SHUF_FULL"),
    }
    chrom = pred["chr"].to_numpy()
    report: dict[str, object] = {
        "resolved_class_counts": observed,
        "alpha_sensitivity": {},
        "WT_embedding_baseline": {},
        "subclass_performance": {},
        "replicated_only": {},
        "blocked_bootstrap": {},
    }

    # Alpha sensitivity and a key ablation: WT embedding versus paired embedding change.
    for name, (a, b) in definitions.items():
        y = pred[f"measured__{name}"].to_numpy(float)
        xdelta = emb[a] - emb[b]
        report["alpha_sensitivity"][name] = {}
        for alpha in ALPHAS:
            p = cv_ridge(xdelta, y, chrom, alpha)
            report["alpha_sensitivity"][name][str(int(alpha))] = metrics(y, p)
        p_wt = cv_ridge(eorig, y, chrom, 3000.0)
        report["WT_embedding_baseline"][name] = metrics(y, p_wt)

        pe = pred[f"existing_probe_oof__{name}"].to_numpy(float)
        pr = pred[f"retrained_delta_probe_chr_oof__{name}"].to_numpy(float)
        report["blocked_bootstrap"][name] = {
            "existing": block_bootstrap_ci(y, pe, chrom),
            "retrained": block_bootstrap_ci(y, pr, chrom),
        }
        report["subclass_performance"][name] = {}
        for cls in ["Stable", "Age-only", "Exercise-only", "Reversal"]:
            m = pred["resolved_class"].to_numpy() == cls
            report["subclass_performance"][name][cls] = {
                "existing": metrics(y[m], pe[m]),
                "retrained": metrics(y[m], pr[m]),
            }
        m = pred["replicated"].to_numpy(bool)
        report["replicated_only"][name] = {
            "n": int(m.sum()),
            "existing": metrics(y[m], pe[m]),
            "retrained": metrics(y[m], pr[m]),
        }

    # Select mechanism-reliable and deliberate model-failure loci from the context contrast.
    name = "context_Motif_minus_Full"
    ycol = f"measured__{name}"
    pcol = f"retrained_delta_probe_chr_oof__{name}"
    work = pred.copy()
    work["residual_abs"] = (work[ycol] - work[pcol]).abs()
    work["sign_concordant"] = np.sign(work[ycol]) == np.sign(work[pcol])
    work["effect_direction"] = np.where(work[ycol] >= 0, "activating", "repressive")
    work["reliability_score"] = (
        zscore(work[ycol].abs())
        + zscore(work[pcol].abs())
        - zscore(work["residual_abs"])
    )

    selected = []
    for cls in ["Stable", "Age-only", "Exercise-only", "Reversal"]:
        for direction in ["activating", "repressive"]:
            part = work[
                (work["resolved_class"] == cls)
                & (work["effect_direction"] == direction)
                & work["replicated"]
                & work["sign_concordant"]
            ].nlargest(12, "reliability_score")
            part = part.assign(panel_type="mechanism-reliable")
            selected.append(part)
        failure = work[
            (work["resolved_class"] == cls) & work["replicated"]
        ].nlargest(6, "residual_abs")
        failure = failure.assign(panel_type="model-failure-control")
        selected.append(failure)

    # AP-1-responsive controls are selected independently, since AP-1 does not separate classes.
    ap1 = "AP1_Endo_minus_Motif"
    ap = pred.copy()
    ap["ap1_abs"] = ap[f"measured__{ap1}"].abs()
    ap["ap1_residual_abs"] = (
        ap[f"measured__{ap1}"] - ap[f"retrained_delta_probe_chr_oof__{ap1}"]
    ).abs()
    ap["ap1_concordant"] = (
        np.sign(ap[f"measured__{ap1}"])
        == np.sign(ap[f"retrained_delta_probe_chr_oof__{ap1}"])
    )
    ap["reliability_score"] = zscore(ap["ap1_abs"]) - zscore(ap["ap1_residual_abs"])
    ap = ap[ap["replicated"] & ap["ap1_concordant"]].nlargest(24, "reliability_score")
    ap = ap.assign(panel_type="AP1-responsive-control")
    selected.append(ap)

    panel = pd.concat(selected, ignore_index=True).drop_duplicates(
        ["Coords", "panel_type"]
    )
    keep = [
        "panel_type",
        "Coords",
        "resolved_class",
        "Group",
        "chr",
        "replicated",
        "effect_direction",
        "reliability_score",
        "sign_concordant",
        "residual_abs",
    ]
    for name in definitions:
        keep += [
            f"measured__{name}",
            f"existing_probe_oof__{name}",
            f"retrained_delta_probe_chr_oof__{name}",
        ]
    keep = [c for c in keep if c in panel.columns]
    panel[keep].to_csv(OUT / "mpra_in_silico_followup_loci.csv", index=False)
    report["candidate_panel_counts"] = panel["panel_type"].value_counts().to_dict()

    with (OUT / "mpra_contrast_followup.json").open("w") as handle:
        json.dump(report, handle, indent=2)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
