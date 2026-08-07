#!/usr/bin/env python3
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Benchmark existing and retrained MPRA contrast probes.
# Steps: 1) Build paired construct representations; 2) perform chromosome-held-out evaluation; 3) export performance and prediction tables.
# =======================================
"""Benchmark existing and retrained AlphaGenome probes on paired MPRA perturbations.

Read-only with respect to the source project. Outputs are written to this task's
outputs directory. The key unit of validation is the parent enhancer: no parent
appears in both train and test. Chromosome-held-out Ridge benchmarks are used for
new contrast heads.
"""

from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
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
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold, KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


PROJECT = Path(publication_path("MPRA_KOO_PROJECT"))
TASK = Path(publication_path("MPRA_SEQLET_TASK_ROOT"))
OUT = TASK / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
KOORUN = PROJECT / "koorun"
sys.path.insert(0, str(KOORUN))

from koo_model import CONSTRUCT_MODE, device, load_single  # noqa: E402
from alphagenome_pytorch.utils.sequence import sequence_to_onehot_tensor  # noqa: E402


KINDS = ["ORIG", "MOTIF_SHUF", "SHUF_FULL"]
MEASURED = {
    "ORIG": "log2_alpha_Endo",
    "MOTIF_SHUF": "log2_alpha_MotifShuff",
    "SHUF_FULL": "log2_alpha_FullShuff",
}
CONTRASTS = {
    "total_Endo_minus_Full": ("ORIG", "SHUF_FULL"),
    "AP1_Endo_minus_Motif": ("ORIG", "MOTIF_SHUF"),
    "context_Motif_minus_Full": ("MOTIF_SHUF", "SHUF_FULL"),
}
ALPHA = 3000.0


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
            key = (kind, coord)
            buf = []
        else:
            buf.append(line)
    if key is not None:
        records[key] = "".join(buf)
    return records


def core(seq: str) -> str:
    return seq[15:-15].upper()


def safe_corr(x: np.ndarray, y: np.ndarray, kind: str = "pearson") -> float:
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return float("nan")
    if kind == "spearman":
        return float(stats.spearmanr(x[ok], y[ok]).statistic)
    return float(stats.pearsonr(x[ok], y[ok]).statistic)


def metric_block(y: np.ndarray, p: np.ndarray) -> dict[str, float | int]:
    ok = np.isfinite(y) & np.isfinite(p)
    yy, pp = y[ok], p[ok]
    if len(yy) < 3:
        return {"n": int(len(yy))}
    slope, intercept = np.polyfit(pp, yy, 1)
    return {
        "n": int(len(yy)),
        "pearson": safe_corr(yy, pp),
        "spearman": safe_corr(yy, pp, "spearman"),
        "r2": float(r2_score(yy, pp)),
        "rmse": float(mean_squared_error(yy, pp) ** 0.5),
        "calibration_slope_measured_on_predicted": float(slope),
        "calibration_intercept": float(intercept),
        "sign_accuracy": float(np.mean(np.sign(yy) == np.sign(pp))),
    }


def load_head_states() -> list[dict[str, np.ndarray]]:
    out = []
    for p in sorted(
        (PROJECT / "bundle/ensemble").glob("fold_*_head.pt"),
        key=lambda x: int(x.stem.split("_")[1]),
    ):
        state = torch.load(p, map_location="cpu", weights_only=True)
        out.append({k: v.detach().cpu().numpy() for k, v in state.items()})
    return out


def head_predict(flat_embedding: np.ndarray, state: dict[str, np.ndarray]) -> np.ndarray:
    # The encoder stores (B, 1536, 2); MPRAHead transposes this to (B, 2, 1536)
    # before its position-wise LayerNorm and flatten operation.
    x = (
        flat_embedding.reshape(len(flat_embedding), 1536, 2)
        .transpose(0, 2, 1)
        .astype(np.float64)
    )
    mu = x.mean(axis=-1, keepdims=True)
    var = x.var(axis=-1, keepdims=True)
    x = (x - mu) / np.sqrt(var + 1e-5)
    x = x * state["norm.weight"][None, None, :] + state["norm.bias"][None, None, :]
    x = x.reshape(len(x), -1)
    return x @ state["output_layer.weight"].reshape(-1) + float(
        state["output_layer.bias"].reshape(-1)[0]
    )


@torch.no_grad()
def encode_sequences(model, seqs: list[str], batch_size: int = 64) -> np.ndarray:
    blocks = []
    for i in range(0, len(seqs), batch_size):
        chunk = seqs[i : i + batch_size]
        assembled = model.construct_spec.assemble_sequences(chunk, mode=CONSTRUCT_MODE)
        onehot = torch.stack(
            [
                sequence_to_onehot_tensor(s, dtype=torch.float32, device=device)
                for s in assembled
            ]
        )
        emb = model.encode(onehot).float().cpu().numpy().reshape(len(chunk), -1)
        blocks.append(emb)
        print(f"encoded {min(i + batch_size, len(seqs))}/{len(seqs)}", flush=True)
    return np.concatenate(blocks).astype(np.float32)


def cv_ridge(
    x: np.ndarray, y: np.ndarray, chrom: np.ndarray, alpha: float = ALPHA
) -> tuple[np.ndarray, dict[str, float | int]]:
    pred = np.full(len(y), np.nan, dtype=float)
    cv = GroupKFold(n_splits=5)
    for train, test in cv.split(x, y, groups=chrom):
        model = make_pipeline(StandardScaler(), Ridge(alpha=alpha))
        model.fit(x[train], y[train])
        pred[test] = model.predict(x[test])
    return pred, metric_block(y, pred)


def kmer_matrix(seqs: list[str], k: int) -> np.ndarray:
    alphabet = "ACGT"
    kmers = [""]
    for _ in range(k):
        kmers = [p + b for p in kmers for b in alphabet]
    lookup = {v: i for i, v in enumerate(kmers)}
    x = np.zeros((len(seqs), len(kmers)), dtype=np.float32)
    for r, s in enumerate(seqs):
        for j in range(len(s) - k + 1):
            x[r, lookup[s[j : j + k]]] += 1.0
        x[r] /= max(1, len(s) - k + 1)
    return x


def permutation_p_group_difference(
    values: np.ndarray, dynamic: np.ndarray, rng: np.random.Generator, nperm: int = 10000
) -> float:
    obs = float(values[dynamic].mean() - values[~dynamic].mean())
    hits = 0
    for _ in range(nperm):
        lab = rng.permutation(dynamic)
        diff = float(values[lab].mean() - values[~lab].mean())
        hits += abs(diff) >= abs(obs)
    return float((hits + 1) / (nperm + 1))


def main() -> None:
    print(f"device={device}", flush=True)
    df = (
        pd.read_csv(PROJECT / "bundle/MPRA_Results_with_TSS_distance.csv")
        .drop_duplicates("Coords")
        .reset_index(drop=True)
    )
    fasta = read_fasta(PROJECT / "bundle/MPRA_AP1_fibroblast.fasta")
    complete = [
        c for c in df["Coords"] if all((kind, c) in fasta for kind in KINDS)
    ]
    df = df.set_index("Coords").loc[complete].reset_index()
    coords = df["Coords"].to_numpy()
    cores = {kind: [core(fasta[(kind, c)]) for c in coords] for kind in KINDS}

    train = pd.read_csv(PROJECT / "bundle/log2_alpha_Endo_train_data.tsv", sep="\t")
    seq_to_fold: dict[str, list[int]] = defaultdict(list)
    for seq, fold in zip(train["seq"], train["fold"]):
        seq_to_fold[str(seq).upper()].append(int(fold))
    folds = np.array(
        [seq_to_fold[s][0] if len(set(seq_to_fold[s])) == 1 else -1 for s in cores["ORIG"]]
    )
    if np.any(folds < 1):
        raise RuntimeError(f"Could not uniquely map {(folds < 1).sum()} sequences to folds")
    chrom = df["chr"].to_numpy()
    dynamic = df["Group"].to_numpy() != "Stable"
    measured = {kind: df[col].to_numpy(float) for kind, col in MEASURED.items()}

    # Endogenous embeddings already exist and are in the same drop-duplicate coordinate order.
    existing_df = (
        pd.read_csv(PROJECT / "bundle/MPRA_Results_with_TSS_distance.csv")
        .drop_duplicates("Coords")
        .reset_index(drop=True)
    )
    existing_e = np.load(PROJECT / "koorun/out/embeddings.npy")
    existing_map = {c: existing_e[i] for i, c in enumerate(existing_df["Coords"])}
    embeddings: dict[str, np.ndarray] = {
        "ORIG": np.vstack([existing_map[c] for c in coords]).astype(np.float32)
    }
    cache = OUT / "mpra_triplet_embeddings_motif_full.npz"
    if cache.exists():
        z = np.load(cache)
        if list(z["coords"].astype(str)) == list(coords.astype(str)):
            embeddings["MOTIF_SHUF"] = z["motif"].astype(np.float32)
            embeddings["SHUF_FULL"] = z["full"].astype(np.float32)
            print(f"loaded cached shuffled embeddings {cache}", flush=True)
    if "MOTIF_SHUF" not in embeddings:
        base = load_single()
        for kind in ["MOTIF_SHUF", "SHUF_FULL"]:
            print(f"extracting {kind} embeddings", flush=True)
            embeddings[kind] = encode_sequences(base, cores[kind])
        np.savez_compressed(
            cache,
            coords=coords.astype(str),
            motif=embeddings["MOTIF_SHUF"],
            full=embeddings["SHUF_FULL"],
        )
        del base

    # Existing heads: ensemble and parent-fold OOF predictions.
    states = load_head_states()
    pred_allheads: dict[str, np.ndarray] = {}
    pred_oof: dict[str, np.ndarray] = {}
    for kind in KINDS:
        per_head = np.vstack([head_predict(embeddings[kind], s) for s in states])
        pred_allheads[kind] = per_head.mean(axis=0)
        oof = np.empty(len(df), dtype=float)
        for i, fold in enumerate(folds):
            # Heads 1..9 each omit their named fold; every head omits fold 10.
            oof[i] = per_head[fold - 1, i] if fold <= 9 else per_head[:, i].mean()
        pred_oof[kind] = oof

    result: dict[str, object] = {
        "n": int(len(df)),
        "n_by_group": df["Group"].value_counts().to_dict(),
        "device": str(device),
        "existing_probe": {"activity": {}, "contrasts": {}, "fold10": {}},
        "retrained_contrast_heads": {},
        "kmer_delta_baselines": {},
        "group_effects": {},
        "notes": [
            "Existing head OOF uses fold_i_head for parents in folds 1-9 and the nine-head ensemble for fold 10.",
            "New Ridge probes use chromosome-held-out five-fold validation and alpha=3000 fixed a priori from the existing project analyses.",
            "The frozen encoder is not updated in this benchmark.",
        ],
    }

    for kind in KINDS:
        result["existing_probe"]["activity"][kind] = metric_block(
            measured[kind], pred_oof[kind]
        )
        m10 = folds == 10
        result["existing_probe"]["fold10"][kind] = metric_block(
            measured[kind][m10], pred_oof[kind][m10]
        )

    contrast_rows = []
    retrained_preds: dict[str, np.ndarray] = {}
    for name, (a, b) in CONTRASTS.items():
        y = measured[a] - measured[b]
        p_existing = pred_oof[a] - pred_oof[b]
        result["existing_probe"]["contrasts"][name] = metric_block(y, p_existing)

        x_delta = embeddings[a] - embeddings[b]
        p_delta, met_delta = cv_ridge(x_delta, y, chrom)
        result["retrained_contrast_heads"][name] = met_delta
        retrained_preds[name] = p_delta

        # K-mer change baselines test whether the embedding adds beyond local word counts.
        result["kmer_delta_baselines"][name] = {}
        for k in [2, 3, 4]:
            xk = kmer_matrix(cores[a], k) - kmer_matrix(cores[b], k)
            _, mk = cv_ridge(xk, y, chrom, alpha=10.0)
            result["kmer_delta_baselines"][name][str(k)] = mk

        measured_gap = float(y[dynamic].mean() - y[~dynamic].mean())
        existing_gap = float(p_existing[dynamic].mean() - p_existing[~dynamic].mean())
        retrained_gap = float(p_delta[dynamic].mean() - p_delta[~dynamic].mean())
        result["group_effects"][name] = {
            "measured_stable_mean": float(y[~dynamic].mean()),
            "measured_dynamic_mean": float(y[dynamic].mean()),
            "measured_dynamic_minus_stable": measured_gap,
            "measured_permutation_p": permutation_p_group_difference(
                y, dynamic, np.random.default_rng(20260722)
            ),
            "existing_probe_stable_mean": float(p_existing[~dynamic].mean()),
            "existing_probe_dynamic_mean": float(p_existing[dynamic].mean()),
            "existing_probe_dynamic_minus_stable": existing_gap,
            "retrained_probe_stable_mean": float(p_delta[~dynamic].mean()),
            "retrained_probe_dynamic_mean": float(p_delta[dynamic].mean()),
            "retrained_probe_dynamic_minus_stable": retrained_gap,
        }
        for i in range(len(df)):
            contrast_rows.append(
                {
                    "Coords": coords[i],
                    "Group": df.loc[i, "Group"],
                    "chr": chrom[i],
                    "training_fold": int(folds[i]),
                    "contrast": name,
                    "measured": float(y[i]),
                    "existing_probe_oof": float(p_existing[i]),
                    "retrained_delta_probe_chr_oof": float(p_delta[i]),
                }
            )

    # Fold-head identity diagnostic: named head should be the held-out head for folds 1..9.
    fold_diagnostic = {}
    e = embeddings["ORIG"]
    per_head = np.vstack([head_predict(e, s) for s in states])
    for f in range(1, 10):
        m = folds == f
        named = safe_corr(measured["ORIG"][m], per_head[f - 1, m])
        others = [safe_corr(measured["ORIG"][m], per_head[h, m]) for h in range(9) if h != f - 1]
        fold_diagnostic[str(f)] = {
            "n": int(m.sum()),
            "named_head_r": named,
            "other_heads_mean_r": float(np.nanmean(others)),
        }
    result["fold_head_diagnostic"] = fold_diagnostic

    # Output table with all raw activity predictions as well.
    activity_table = df[["Coords", "Group", "chr", "replicated"]].copy()
    activity_table["training_fold"] = folds
    for kind in KINDS:
        activity_table[f"measured_{kind}"] = measured[kind]
        activity_table[f"pred_existing_oof_{kind}"] = pred_oof[kind]
        activity_table[f"pred_existing_allheads_{kind}"] = pred_allheads[kind]
    activity_table.to_csv(OUT / "mpra_triplet_model_predictions.csv", index=False)
    pd.DataFrame(contrast_rows).to_csv(
        OUT / "mpra_contrast_benchmark_predictions.csv", index=False
    )
    with (OUT / "mpra_contrast_benchmark.json").open("w") as handle:
        json.dump(result, handle, indent=2)

    # Compact figure: measured vs predicted contrasts, existing vs retrained.
    fig, axes = plt.subplots(3, 2, figsize=(8.2, 10.2))
    palette = np.where(dynamic, "#C44E52", "#4C566A")
    for row, (name, (a, b)) in enumerate(CONTRASTS.items()):
        y = measured[a] - measured[b]
        panels = [
            (pred_oof[a] - pred_oof[b], "existing native-activity probe"),
            (retrained_preds[name], "retrained delta probe"),
        ]
        for col, (p, title) in enumerate(panels):
            ax = axes[row, col]
            ax.scatter(y, p, s=7, alpha=0.32, c=palette, linewidths=0)
            lo = float(min(np.nanpercentile(y, 1), np.nanpercentile(p, 1)))
            hi = float(max(np.nanpercentile(y, 99), np.nanpercentile(p, 99)))
            ax.plot([lo, hi], [lo, hi], ls="--", lw=0.7, color="#555")
            ax.axhline(0, lw=0.45, color="#999")
            ax.axvline(0, lw=0.45, color="#999")
            ax.set_xlabel("measured paired effect (log2)")
            ax.set_ylabel("predicted paired effect (log2)")
            met = metric_block(y, p)
            ax.set_title(
                f"{name}\n{title}\nr={met.get('pearson', np.nan):.2f}; "
                f"rho={met.get('spearman', np.nan):.2f}; R2={met.get('r2', np.nan):.2f}",
                fontsize=8,
            )
            ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(
        "Can the current AlphaGenome model predict measured MPRA perturbation effects?",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(OUT / "mpra_contrast_benchmark.png", dpi=240, bbox_inches="tight")
    fig.savefig(OUT / "mpra_contrast_benchmark.pdf", bbox_inches="tight")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
    main()
