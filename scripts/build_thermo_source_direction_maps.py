#!/usr/bin/env python3
"""Build four fresh, mutually-consistent reaction-direction maps, one per
thermodynamic source, using the *unmodified* ModelSEED reversibility cascade.

Unlike the older per-source artifacts elsewhere in this repo (some of which
rerun the cascade, some of which just overlay a pre-computed operator from a
CSV generated from a different git ref), this script:

  * reads reactions/compounds directly from the *live, currently checked-out*
    ModelSEEDDatabase working tree via ``BiochemPy`` -- the same loader
    ``Estimate_Reaction_Reversibility.py`` itself uses -- so results are
    guaranteed current;
  * runs the exact same ``DEFAULT_HEURISTICS`` cascade
    (``reversibility_heuristics.py``) for all four sources, varying only the
    energy input;
  * uses ``per_source_energy(label)`` for Group Contribution / eQuilibrator /
    dGPredictor -- i.e. that source's *own* dG -- not ``top_level_energy``,
    which would use the canonical top-level dG (merely gated by source
    eligibility) and, for eQuilibrator, silently fall back to Group
    Contribution's reversibility when eQuilibrator itself has no data. Both
    behaviors would break "only source X" purity.

Outputs (under ``results/thermo_source_fba_all_models/``):
  * ``rxn_directions_{modelseed,group_contribution,equilibrator,dgpredictor}.json``
    -- ``{rxn_id: operator}`` for reactions where that source has usable data.
  * ``rxn_source_coverage.csv`` -- one row per non-EMPTY reaction: which
    sources have data, and what operator each computed.
  * ``cpd_source_coverage.csv`` -- one row per compound: which sources have a
    defined formation energy (dGPredictor never does -- it's a direct
    reaction-DeltaG predictor with no compound-level analog; its column is
    always False, documented rather than omitted).
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

from BiochemPy import Reactions, Compounds  # noqa: E402
from reversibility_heuristics import (  # noqa: E402
    DEFAULT_HEURISTICS,
    run_reversibility,
    top_level_energy,
    per_source_energy,
    SENTINEL_DG,
)

# source key -> (energy source callable, thermodynamics-sublist label or None)
REACTION_SOURCES = {
    "modelseed": (top_level_energy(""), None),
    "group_contribution": (per_source_energy("Group contribution"), "Group contribution"),
    "equilibrator": (per_source_energy("eQuilibrator"), "eQuilibrator"),
    "dgpredictor": (per_source_energy("dGPredictor"), "dGPredictor"),
}


def _compound_has_energy(cpd_entry: dict, label: str | None) -> bool:
    """True iff the compound carries a non-sentinel dG for this source.

    ``label is None`` means the canonical top-level ``deltag``/``deltagerr``
    (the "original ModelSEED" source)."""
    if label is None:
        dg = cpd_entry.get("deltag")
    else:
        thermo = cpd_entry.get("thermodynamics") or {}
        pair = thermo.get(label)
        dg = pair[0] if pair else None
    if dg is None:
        return False
    try:
        dg = float(dg)
    except (TypeError, ValueError):
        return False
    return dg != SENTINEL_DG


def build_reaction_maps(reactions: dict, dgp_mask: set | None = None) -> tuple[dict, list]:
    """Returns (``{source: {rxn_id: operator}}``, coverage_rows).

    ``dgp_mask`` (build_dgpredictor_kegg_mask.py) lists reactions whose stored
    dGPredictor value was predicted from a KEGG reaction that is not theirs.
    They are withheld from the ``dgpredictor`` map only: the reaction then has
    no dGPredictor-defined direction, exactly as if dGPredictor had never
    covered it, and downstream FBA falls back to that source's normal
    no-coverage behaviour. Group Contribution and eQuilibrator are untouched.
    """
    dgp_mask = dgp_mask or set()
    maps = {src: {} for src in REACTION_SOURCES}
    coverage_rows = []
    n_masked = 0
    for rxn_id in sorted(reactions):
        rxn_entry = reactions[rxn_id]
        if rxn_entry.get("status") == "EMPTY":
            continue
        row = {"rxn_id": rxn_id}
        for src, (energy_source, _label) in REACTION_SOURCES.items():
            if src == "dgpredictor" and rxn_id in dgp_mask:
                row["has_dgpredictor"] = False
                row["op_dgpredictor"] = ""
                n_masked += 1
                continue
            status, operator, _src_label = run_reversibility(
                rxn_entry, energy_source, DEFAULT_HEURISTICS)
            has_data = operator is not None
            row[f"has_{src}"] = has_data
            row[f"op_{src}"] = operator or ""
            if has_data:
                maps[src][rxn_id] = operator
        coverage_rows.append(row)
    if n_masked:
        print(f"  dgpredictor: withheld {n_masked} reaction(s) with an unvouched "
              f"KEGG mapping (dgpredictor_kegg_mask.json)")
    return maps, coverage_rows


def build_compound_coverage(compounds: dict) -> list:
    rows = []
    for cpd_id in sorted(compounds):
        cpd_entry = compounds[cpd_id]
        row = {
            "cpd_id": cpd_id,
            "has_modelseed": _compound_has_energy(cpd_entry, None),
            "has_group_contribution": _compound_has_energy(cpd_entry, "Group contribution"),
            "has_equilibrator": _compound_has_energy(cpd_entry, "eQuilibrator"),
            "has_dgpredictor": False,  # dGPredictor has no compound-level formation energies
        }
        rows.append(row)
    return rows


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("loading reactions from live ModelSEEDDatabase checkout...")
    reactions = Reactions().loadReactions()
    print(f"  {len(reactions)} reactions loaded")

    print("loading compounds...")
    compounds = Compounds().loadCompounds()
    print(f"  {len(compounds)} compounds loaded")

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from build_dgpredictor_kegg_mask import load_mask
    maps, rxn_coverage_rows = build_reaction_maps(reactions, dgp_mask=load_mask())
    for src, m in maps.items():
        out_path = OUT_DIR / f"rxn_directions_{src}.json"
        with open(out_path, "w") as fh:
            json.dump(m, fh, indent=1, sort_keys=True)
        print(f"  {src}: {len(m)} reactions with a usable direction -> {out_path}")

    coverage_csv = OUT_DIR / "rxn_source_coverage.csv"
    fieldnames = ["rxn_id"] + [f"{p}_{src}" for src in REACTION_SOURCES for p in ("has", "op")]
    with open(coverage_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rxn_coverage_rows)
    print(f"wrote {coverage_csv} ({len(rxn_coverage_rows)} rows)")

    cpd_rows = build_compound_coverage(compounds)
    cpd_csv = OUT_DIR / "cpd_source_coverage.csv"
    with open(cpd_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(cpd_rows[0].keys()))
        w.writeheader()
        w.writerows(cpd_rows)
    print(f"wrote {cpd_csv} ({len(cpd_rows)} rows)")

    n_nonempty = len(rxn_coverage_rows)
    n_compounds = len(cpd_rows)
    summary = {
        "n_reactions_total": len(reactions),
        "n_reactions_nonempty": n_nonempty,
        "n_compounds_total": n_compounds,
        "reaction_coverage": {src: len(m) for src, m in maps.items()},
        "compound_coverage": {
            "modelseed": sum(1 for r in cpd_rows if r["has_modelseed"]),
            "group_contribution": sum(1 for r in cpd_rows if r["has_group_contribution"]),
            "equilibrator": sum(1 for r in cpd_rows if r["has_equilibrator"]),
            "dgpredictor": 0,
        },
    }
    with open(OUT_DIR / "direction_maps_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
