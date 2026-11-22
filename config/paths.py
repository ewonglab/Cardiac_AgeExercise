# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Resolve configured Python workflow paths from named environment variables.
# Steps: 1) Read a configured root. 2) Fail clearly when it is unset. 3) Join optional path components reproducibly.
# =======================================
"""Resolve publication workflow paths from environment variables."""
from __future__ import annotations

import os
from pathlib import Path


def publication_path(variable: str, *parts: str) -> str:
    root = os.environ.get(variable)
    if not root:
        raise RuntimeError(f"Set {variable} in config/.env before running this workflow.")
    return str(Path(root, *parts))
