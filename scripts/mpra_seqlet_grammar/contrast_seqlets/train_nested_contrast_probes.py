#!/usr/bin/env python3
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Train chromosome-held-out ridge probes for MPRA perturbation contrasts.
# Steps: 1) Load frozen representations and measured contrasts; 2) fit nested cross-validation models; 3) save coefficients and predictions.
# =======================================
"""Train manuscript-ready perturbation-aware probes with nested outer evaluation.

Alpha is selected by generalized cross-validation using training data only. The
outer test folds hold out complete chromosomes. Effective coefficients are saved
on the raw flattened encoder scale for direct scoring and attribution.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from config.paths import publication_path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import RidgeCV
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler


PROJECT = Path(publication_path("MPRA_KOO_PROJECT"))
TASK = Path(publication_path("MPRA_SEQLET_TASK_ROOT"))
OUT = TASK / "outputs"
ALPHAS = np.array([100.0, 300.0, 1000.0, 3000.0, 10000.0, 30000.0, 100000.0])
TARGETS = {
    "total_Endo_minus_Full": ("ORIG", "SHUF_FULL"),
    "AP1_Endo_minus_Motif": ("ORIG", "MOTIF_SHUF"),
    "context_Motif_minus_Full": ("MOTIF_SHUF", "SHUF_FULL"),
}


def metrics(y: np.ndarray, p: np.ndarray) -> dict[str, float | int]:
    slope, intercept = np.polyfit(p, y, 1)
    return {
        "n": int(len(y)),
        "pearson": float(stats.pearsonr(y, p).statistic),
        "spearman": float(stats.spearmanr(y, p).statistic),
        "r2": float(r2_score(y, p)),
        "rmse": float(mean_squared_error(y, p) ** 0.5),
        "sign_accuracy": float(np.mean(np.sign(y) == np.sign(p))),
        "calibration_slope": float(slope),
        "calibration_intercept": float(intercept),
    }


def raw_coefficients(scaler: StandardScaler, ridge: RidgeCV) -> tuple[np.ndarray, float]:
    coef = np.asarray(ridge.coef_, dtype=float).reshape(-1)
    effective = coef / scaler.scale_
    intercept = float(ridge.intercept_ - np.dot(coef, scaler.mean_ / scaler.scale_))
    return effective.astype(np.float32), intercept


def main() -> None:
    long = pd.read_csv(OUT / "mpra_contrast_benchmark_predictions.csv")
    base = long.drop_duplicates("Coords")[["Coords", "Group", "chr", "training_fold"]]
    base = base.reset_index(drop=True)
    coords = base["Coords"].astype(str).to_numpy()
    chrom = base["chr"].astype(str).to_numpy()
    ymap = {
        name: long[long["contrast"] == name]
        .set_index("Coords")
        .loc[coords, "measured"]
        .to_numpy(float)
        for name in TARGETS
    }

    cache = np.load(OUT / "mpra_triplet_embeddings_motif_full.npz")
    cache_index = {str(c): i for i, c in enumerate(cache["coords"].astype(str))}
    source = (
        pd.read_csv(PROJECT / "bundle/MPRA_Results_with_TSS_distance.csv")
        .drop_duplicates("Coords")
        .reset_index(drop=True)
    )
    eorig_all = np.load(PROJECT / "koorun/out/embeddings.npy")
    orig_index = {str(c): i for i, c in enumerate(source["Coords"].astype(str))}
    emb = {
        "ORIG": np.vstack([eorig_all[orig_index[c]] for c in coords]).astype(np.float32),
        "MOTIF_SHUF": np.vstack([cache["motif"][cache_index[c]] for c in coords]).astype(np.float32),
        "SHUF_FULL": np.vstack([cache["full"][cache_index[c]] for c in coords]).astype(np.float32),
    }

    names = list(TARGETS)
    outer_coef = np.zeros((len(names), 5, 3072), dtype=np.float32)
    outer_intercept = np.zeros((len(names), 5), dtype=np.float32)
    outer_alpha = np.zeros((len(names), 5), dtype=np.float32)
    final_coef = np.zeros((len(names), 3072), dtype=np.float32)
    final_intercept = np.zeros(len(names), dtype=np.float32)
    final_alpha = np.zeros(len(names), dtype=np.float32)
    nested_predictions = base.copy()
    report: dict[str, object] = {
        "design": {
            "outer_cv": "5-fold GroupKFold by chromosome",
            "alpha_selection": "RidgeCV generalized cross-validation on outer-training data only",
            "alpha_grid": ALPHAS.tolist(),
            "unit": "parent enhancer; paired embedding difference",
        },
        "targets": {},
    }

    outer = list(GroupKFold(5).split(np.zeros(len(coords)), groups=chrom))
    for t, name in enumerate(names):
        a, b = TARGETS[name]
        x = (emb[a] - emb[b]).astype(np.float64)
        y = ymap[name]
        pred = np.full(len(y), np.nan)
        folds = []
        print(f"training {name}", flush=True)
        for k, (train, test) in enumerate(outer):
            scaler = StandardScaler().fit(x[train])
            xt = scaler.transform(x[train])
            xv = scaler.transform(x[test])
            model = RidgeCV(alphas=ALPHAS, cv=None, gcv_mode="auto").fit(xt, y[train])
            pred[test] = model.predict(xv)
            coef, intercept = raw_coefficients(scaler, model)
            outer_coef[t, k] = coef
            outer_intercept[t, k] = intercept
            outer_alpha[t, k] = float(model.alpha_)
            folds.append(
                {
                    "fold": k + 1,
                    "alpha": float(model.alpha_),
                    "test_chromosomes": sorted(set(chrom[test])),
                    "metrics": metrics(y[test], pred[test]),
                }
            )
            print(f"  fold {k+1}: alpha={model.alpha_:.0f} R2={r2_score(y[test], pred[test]):.3f}", flush=True)

        scaler = StandardScaler().fit(x)
        model = RidgeCV(alphas=ALPHAS, cv=None, gcv_mode="auto").fit(scaler.transform(x), y)
        coef, intercept = raw_coefficients(scaler, model)
        final_coef[t] = coef
        final_intercept[t] = intercept
        final_alpha[t] = float(model.alpha_)
        nested_predictions[f"measured__{name}"] = y
        nested_predictions[f"nested_chr_oof__{name}"] = pred
        report["targets"][name] = {
            "nested_metrics": metrics(y, pred),
            "outer_folds": folds,
            "final_alpha": float(model.alpha_),
            "outer_alpha_counts": {
                str(int(a)): int(np.sum(outer_alpha[t] == a)) for a in np.unique(outer_alpha[t])
            },
            "outer_weight_cosine_mean": float(
                np.mean(
                    [
                        np.dot(outer_coef[t, i], outer_coef[t, j])
                        / (
                            np.linalg.norm(outer_coef[t, i])
                            * np.linalg.norm(outer_coef[t, j])
                            + 1e-12
                        )
                        for i in range(5)
                        for j in range(i + 1, 5)
                    ]
                )
            ),
        }

    np.savez_compressed(
        OUT / "mpra_nested_contrast_models.npz",
        target_names=np.array(names),
        coords=coords,
        outer_coef=outer_coef,
        outer_intercept=outer_intercept,
        outer_alpha=outer_alpha,
        final_coef=final_coef,
        final_intercept=final_intercept,
        final_alpha=final_alpha,
    )
    nested_predictions.to_csv(OUT / "mpra_nested_contrast_predictions.csv", index=False)
    with (OUT / "mpra_nested_contrast_models.json").open("w") as handle:
        json.dump(report, handle, indent=2)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
