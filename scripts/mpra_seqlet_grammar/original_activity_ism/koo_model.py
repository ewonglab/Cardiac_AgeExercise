# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Provide shared loading and prediction helpers for the frozen AlphaGenome MPRA model.
# Steps: 1) Load the encoder and fold heads; 2) encode MPRA inserts; 3) return single-head or ensemble activity predictions.
# =======================================
"""Shared loader + predict helpers for Koo's AlphaGenome encoder+probe MPRA model.

Wraps the single-head checkpoint and the 9-fold ensemble (shared frozen encoder,
one linear head per training split) behind a plain `predict(cores) -> np.ndarray`
interface so our existing MPRA analyses (ISM flanking, spacing grammar) can call it.
"""
from pathlib import Path
import os
import copy
import numpy as np
import torch

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from config.paths import publication_path

from alphagenome_encoder_ft import EncoderMPRAModel
from alphagenome_pytorch.utils.sequence import sequence_to_onehot_tensor

BUNDLE = Path(publication_path("MPRA_KOO_PROJECT", "bundle"))
SINGLE_CKPT = BUNDLE / "single" / "best.pt"
ENSEMBLE_DIR = BUNDLE / "ensemble"
CONSTRUCT_MODE = "adapters"   # variable core flanked by the MPRA adapters, as in the notebook
BASES = "ACGT"

# Device: KOO_DEVICE env overrides (e.g. "cuda", "mps", "cpu"); else CUDA > CPU.
# MPS is opt-in (KOO_DEVICE=mps) because some encoder ops may fall back to CPU.
_dev = os.environ.get("KOO_DEVICE")
if _dev:
    device = torch.device(_dev)
elif torch.cuda.is_available():
    device = torch.device("cuda")
else:
    device = torch.device("cpu")
torch.set_grad_enabled(False)


def load_single():
    m = EncoderMPRAModel.from_checkpoint(SINGLE_CKPT, device=device)
    m.eval()
    return m


def load_ensemble(base=None):
    """Return (base_model, [heads]) — one head per fold on the single shared encoder."""
    base = base or load_single()
    head_files = sorted(ENSEMBLE_DIR.glob("fold_*_head.pt"),
                        key=lambda p: int(p.stem.split("_")[1]))
    heads = []
    for hf in head_files:
        h = copy.deepcopy(base.head)
        h.load_state_dict(torch.load(hf, map_location="cpu", weights_only=True))
        h.to(device).eval()
        heads.append(h)
    return base, heads


@torch.no_grad()
def predict_single(model, cores, batch_size=64):
    out = []
    for i in range(0, len(cores), batch_size):
        p = model.predict_sequences([c.strip().upper() for c in cores[i:i + batch_size]],
                                    construct_mode=CONSTRUCT_MODE)
        out.append(p.float().cpu().numpy().reshape(-1))
    return np.concatenate(out)


@torch.no_grad()
def predict_ensemble(base, heads, cores, batch_size=64):
    """Run the shared frozen encoder once per sequence, average all fold heads."""
    preds = []
    for i in range(0, len(cores), batch_size):
        chunk = [s.strip().upper() for s in cores[i:i + batch_size]]
        asm = base.construct_spec.assemble_sequences(chunk, mode=CONSTRUCT_MODE)
        onehot = torch.stack([sequence_to_onehot_tensor(s, dtype=torch.float32, device=device)
                              for s in asm])
        emb = base.encode(onehot)                        # frozen encoder, run once
        fold = np.stack([h(emb).float().cpu().numpy().reshape(-1) for h in heads])
        preds.append(fold.mean(axis=0))
    return np.concatenate(preds)


# --- metrics ---
def pearson(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.corrcoef(a, b)[0, 1])


def spearman(a, b):
    from scipy.stats import rankdata
    return pearson(rankdata(a), rankdata(b))
