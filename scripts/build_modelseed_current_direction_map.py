#!/usr/bin/env python3
"""Build the reaction-direction map for "what is currently in ModelSEED":
eQuilibrator's own reversibility call, falling back to the reaction's existing
(Group-Contribution-backed) canonical reversibility when eQuilibrator has no
data for that reaction.

This is NOT a reimplementation -- it calls ``estimate_one(rxn_entry, 'EQ')``
directly from the live ModelSEEDDatabase checkout's
``Estimate_Reaction_Reversibility.py``, which is exactly what running that
script's ``EQ`` CLI mode against the live database does:

  * if the reaction has eQuilibrator data, run the unmodified DEFAULT_HEURISTICS
    cascade against the canonical top-level DeltaG (gated on eQuilibrator
    eligibility) -- ``top_level_energy('EQ')``;
  * else, if the reaction has Group Contribution data, fall back to the
    reaction's current stored ``reversibility`` field (``_incomplete_decision``);
  * else, undefined ("?") -- excluded from the output map, same convention as
    the other four thermo-source maps in this pipeline.

Output: ``results/thermo_source_fba_all_models/rxn_directions_modelseed_current.json``
(``{rxn_id: operator}``) and a coverage CSV alongside it.
"""
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR", "/scratch/ctaylor/core_models_analysis"))
MSDB_ROOT = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase"))
OUT_DIR = ANALYSIS_DIR / "results" / "thermo_source_fba_all_models"

sys.path.insert(0, str(MSDB_ROOT / "Libs" / "Python"))
sys.path.insert(0, str(MSDB_ROOT / "Scripts" / "Thermodynamics"))

from BiochemPy import Reactions  # noqa: E402
from Estimate_Reaction_Reversibility import estimate_one  # noqa: E402


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("loading reactions from live ModelSEEDDatabase checkout...")
    reactions = Reactions().loadReactions()
    print(f"  {len(reactions)} reactions loaded")

    direction_map = {}
    rows = []
    n_via_eq = n_via_gc_fallback = n_undefined = 0
    for rxn_id in sorted(reactions):
        rxn_entry = reactions[rxn_id]
        if rxn_entry.get("status") == "EMPTY":
            continue
        status, operator, source_label = estimate_one(rxn_entry, "EQ")
        has_data = operator in (">", "<", "=")
        if has_data:
            direction_map[rxn_id] = operator
            if source_label == "eQuilibrator":
                n_via_eq += 1
            else:
                n_via_gc_fallback += 1
        else:
            n_undefined += 1
        rows.append({
            "rxn_id": rxn_id,
            "has_modelseed_current": has_data,
            "op_modelseed_current": operator if has_data else "",
            "via": ("eQuilibrator" if source_label == "eQuilibrator"
                     else "GC_fallback" if has_data else "undefined"),
        })

    out_path = OUT_DIR / "rxn_directions_modelseed_current.json"
    with open(out_path, "w") as fh:
        json.dump(direction_map, fh, indent=1, sort_keys=True)
    print(f"wrote {out_path} ({len(direction_map)} reactions: "
          f"{n_via_eq} via eQuilibrator, {n_via_gc_fallback} via GC fallback, "
          f"{n_undefined} undefined/excluded)")

    coverage_csv = OUT_DIR / "rxn_source_coverage_modelseed_current.csv"
    with open(coverage_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {coverage_csv} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
