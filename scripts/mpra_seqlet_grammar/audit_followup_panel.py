#!/usr/bin/env python3
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Audit the composition of the mechanistic MPRA follow-up panel.
# Steps: 1) Load the versioned panel table; 2) count membership and unique loci; 3) print class and direction summaries.
# =======================================
"""Report the composition of the 142-coordinate mechanistic follow-up panel."""

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
PANEL = ROOT / "example_inputs" / "mpra_in_silico_followup_loci.csv"


def main() -> None:
    panel = pd.read_csv(PANEL)
    memberships = panel.groupby("Coords")["panel_type"].agg(
        lambda values: ";".join(sorted(set(values)))
    )
    mechanism = panel[panel["panel_type"].eq("mechanism-reliable")]

    print(f"Panel entries: {len(panel)}")
    print(f"Unique coordinates: {panel['Coords'].nunique()}")
    print("\nEntries by panel type:")
    print(panel["panel_type"].value_counts().to_string())
    print("\nUnique coordinates by combined membership:")
    print(memberships.value_counts().to_string())
    print("\nMechanism-reliable entries by enhancer class and effect direction:")
    print(
        pd.crosstab(
            mechanism["resolved_class"], mechanism["effect_direction"]
        ).to_string()
    )


if __name__ == "__main__":
    main()
