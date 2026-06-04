#!/usr/bin/env python3
"""CLI wrapper -- build thermo-source comparison figures from on-disk CSVs.

Reads the canonical CSVs from ``results/thermo_sources/`` and writes PNGs
to ``reports/figures/thermo_sources/`` by delegating to the
``thermo_source_figures`` library.  All matplotlib code lives in that
library; this script just loads inputs and dispatches.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

# Make the project scripts importable when the file is invoked directly.
SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import thermo_source_figures as tsf  # noqa: E402

ROOT = Path("/scratch/ctaylor/core_models_analysis")
IN_DIR = ROOT / "results" / "thermo_sources"
OUT_DIR = ROOT / "reports" / "figures" / "thermo_sources"


def load_inputs():
    """Load the long, coverage, overrides, and direction DataFrames."""
    long_df = pd.read_csv(IN_DIR / "panel_fba_long.csv")
    coverage = {s: pd.read_csv(IN_DIR / f"coverage_{s}.csv") for s in tsf.SOURCES}
    overrides = {s: pd.read_csv(IN_DIR / f"overrides_{s}.csv") for s in tsf.SOURCES}
    rxn_dirs = {
        s: pd.read_csv(IN_DIR / f"rxn_directions_{s}.csv") for s in tsf.SOURCES
    }
    return long_df, coverage, overrides, rxn_dirs


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    long_df, coverage, overrides, rxn_dirs = load_inputs()

    produced = [
        tsf.make_fig_grower_counts(long_df, OUT_DIR / "fig_grower_counts.png"),
        tsf.make_fig_mean_flux(long_df, OUT_DIR / "fig_mean_flux.png"),
        tsf.make_fig_flux_violin(long_df, OUT_DIR / "fig_flux_violin.png"),
        tsf.make_fig_per_model_heatmap(long_df, OUT_DIR / "fig_per_model_heatmap.png"),
        tsf.make_fig_coverage_per_source(coverage, OUT_DIR / "fig_coverage_per_source.png"),
        tsf.make_fig_override_transitions(overrides, OUT_DIR / "fig_override_transitions.png"),
        tsf.make_fig_flux_vs_baseline_scatter(long_df, OUT_DIR / "fig_flux_vs_baseline_scatter.png"),
        tsf.make_fig_dg_distribution_per_source(rxn_dirs, OUT_DIR / "fig_dg_distribution_per_source.png"),
    ]

    bad = [str(p) for p in produced if not p.exists() or p.stat().st_size == 0]
    if bad:
        print("MISSING OR EMPTY:", bad, file=sys.stderr)
        return 1
    for p in produced:
        print(f"OK {p} ({p.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
