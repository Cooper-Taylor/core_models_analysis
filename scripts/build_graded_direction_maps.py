#!/usr/bin/env python3
"""Reaction-direction maps for the graded thermo-source comparison.

Rebuilds the per-source direction maps of ``build_thermo_source_direction_maps.py``
against a NEWER snapshot and adds the graded/recommended variant:

  gc              Group Contribution's own dG only          (Convention A rebuild)
  eq              eQuilibrator's own dG only
  dgpms           dGPredictor-ModelSEED's own dG only        (the retrain, not
                  the KEGG-keyed legacy source -- structurally immune to the
                  17,271-reaction mis-mapping, so no mask is needed)
  graded          per reaction, the dG of the best-GRADED source, any grade
  graded_trusted  same, but reactions whose best grade is BRONZE are dropped
                  entirely -- the model keeps its native bound there

WHY A NEW SCRIPT RATHER THAN AN EDIT. The original reads the live
``/scratch/ctaylor/ModelSEEDDatabase`` working tree, which has neither the
``dGPredictor-ModelSEED`` label nor the Convention A Group Contribution
rebuild, and its outputs back the 2026-08-03 report. This one reads a dev
snapshot so that the three sources match the ones the grades were fitted on.
The earlier artifacts are left untouched.

WHAT IS HELD CONSTANT. The cascade. ``DEFAULT_HEURISTICS`` is imported from the
local ModelSEEDDatabase checkout -- byte-identical in order to the snapshot's
own copy (verified: atp_synthase, abc_transporter, stored_bounds, mmdeltag_band,
low_energy, default) -- so the only thing varying between the five maps is which
dG feeds in.

Outputs: ``results/thermo_grades_fba/rxn_directions_<variant>.json`` +
``rxn_source_coverage.csv`` + ``direction_maps_summary.json``.
"""
from __future__ import annotations

import csv
import glob
import json
import os
import sys
from pathlib import Path

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
# Data: a dev snapshot carrying all three source labels.
MSDB_DATA = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/tmp/devsnap2"))
# Code: the local checkout (the snapshot ships no Libs/Python).
MSDB_CODE = Path(os.environ.get("MSDB_CODE", "/scratch/ctaylor/ModelSEEDDatabase"))
OUT_DIR = ANALYSIS_DIR / "results" / "thermo_grades_fba"

sys.path.insert(0, str(MSDB_CODE / "Scripts" / "Thermodynamics"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from reversibility_heuristics import (  # noqa: E402
    DEFAULT_HEURISTICS, SENTINEL_DG, per_source_energy, run_reversibility,
)
from grade_thermo_sources import recommended_energy_map  # noqa: E402

PER_SOURCE = {
    "gc": "Group contribution",
    "eq": "eQuilibrator",
    "dgpms": "dGPredictor-ModelSEED",
}
# variant -> (grade floor, use the held-out grades?)
GRADED = {
    "graded": ("BRONZE", False),
    "graded_trusted": ("SILVER", False),
    # No TECRDB source and no measurement-derived grade. The only variant that
    # can be scored against TECRDB without circularity -- ``graded`` reproduces
    # the experiment by construction wherever the experiment exists.
    "graded_heldout": ("BRONZE", True),
}
VARIANTS = list(PER_SOURCE) + list(GRADED)


def load_reactions() -> dict:
    """``{rxn_id: entry}`` from the snapshot, normalised the way BiochemPy does.

    ``BiochemPy.Reactions.loadReactions`` only replaces ``None`` with the string
    ``"null"`` (the cascade relies on that for ``notes``), so this reproduces it
    without needing the snapshot to ship ``Libs/Python``.
    """
    out = {}
    for path in sorted(glob.glob(str(MSDB_DATA / "Biochemistry" / "reaction_*.json"))):
        for entry in json.load(open(path)):
            for key, val in list(entry.items()):
                if val is None:
                    entry[key] = "null"
                elif isinstance(val, list):
                    entry[key] = ["null" if v is None else v for v in val]
                elif isinstance(val, dict):
                    entry[key] = {k: ("null" if v is None else v) for k, v in val.items()}
            out[entry["id"]] = entry
    return out


def graded_energy(energy_map: dict):
    """Energy source over ``{rxn_id: (dg, dge, source_label)}`` from the grades."""
    def resolve(rxn_entry):
        hit = energy_map.get(rxn_entry["id"])
        if hit is None:
            return None, None, "graded"
        dg, dge, label = hit
        if dg == SENTINEL_DG:
            return None, None, label
        return dg, dge, label
    return resolve


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    reactions = load_reactions()
    n_nonempty = sum(1 for r in reactions.values() if r.get("status") != "EMPTY")
    print(f"{len(reactions)} reactions from {MSDB_DATA} ({n_nonempty} non-EMPTY)")

    sources = {k: per_source_energy(lbl) for k, lbl in PER_SOURCE.items()}
    chosen = {}
    for variant, (floor, heldout) in GRADED.items():
        emap = recommended_energy_map(min_grade=floor, heldout=heldout)
        chosen[variant] = emap
        sources[variant] = graded_energy(emap)
        print(f"  {variant:15s} best-graded source available for {len(emap)} reactions "
              f"(floor {floor}{', held out' if heldout else ''})")

    maps = {v: {} for v in VARIANTS}
    picked = {v: {} for v in GRADED}
    rows = []
    for rxn_id in sorted(reactions):
        entry = reactions[rxn_id]
        if entry.get("status") == "EMPTY":
            continue
        row = {"rxn_id": rxn_id}
        for v in VARIANTS:
            status, operator, _ = run_reversibility(entry, sources[v], DEFAULT_HEURISTICS)
            row[f"has_{v}"] = operator is not None
            row[f"op_{v}"] = operator or ""
            row[f"status_{v}"] = status or ""
            if operator is not None:
                maps[v][rxn_id] = operator
        for v in GRADED:
            hit = chosen[v].get(rxn_id)
            row[f"src_{v}"] = hit[2] if hit else ""
            picked[v][rxn_id] = hit[2] if hit else ""
        rows.append(row)

    for v in VARIANTS:
        path = OUT_DIR / f"rxn_directions_{v}.json"
        json.dump(maps[v], open(path, "w"), indent=1, sort_keys=True)
        ops = {}
        for op in maps[v].values():
            ops[op] = ops.get(op, 0) + 1
        print(f"  {v:15s} {len(maps[v]):6d} directions   " +
              "  ".join(f"{o}:{n}" for o, n in sorted(ops.items())))

    cov = OUT_DIR / "rxn_source_coverage.csv"
    fields = ["rxn_id"] + [f"{p}_{v}" for v in VARIANTS for p in ("has", "op", "status")] \
        + [f"src_{v}" for v in GRADED]
    with open(cov, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {cov} ({len(rows)} rows)")

    summary = {
        "msdb_data": str(MSDB_DATA), "msdb_code": str(MSDB_CODE),
        "n_reactions_nonempty": len(rows),
        "coverage": {v: len(maps[v]) for v in VARIANTS},
        "operator_mix": {v: {op: sum(1 for o in maps[v].values() if o == op)
                             for op in ("=", ">", "<", "?")} for v in VARIANTS},
        "graded_source_mix": {
            v: {s: sum(1 for x in picked[v].values() if x == s)
                for s in sorted(set(picked[v].values()) - {""})} for v in GRADED},
    }
    json.dump(summary, open(OUT_DIR / "direction_maps_summary.json", "w"), indent=2)
    print(json.dumps(summary["graded_source_mix"], indent=1))


if __name__ == "__main__":
    main()
