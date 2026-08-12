#!/usr/bin/env python3
"""Characterise the reaction directions the core models ship with.

The ``implicit`` FBA variant applies no override -- it runs each Kegg2 core
model on the bounds baked in when the model was built. Everywhere else in this
analysis that variant is only a growth baseline; this script asks what those
bounds actually SAY, and whether they are right.

Three questions:

  1. What direction does each model assign natively, and do models agree?
     (They do -- unanimously, on all 239 core reactions -- which means the
     implicit bounds are effectively one global map and can be scored like any
     other direction source.)
  2. How accurate is that map against the experimental reference, i.e. the
     cascade run on measured TECRDB energies?
  3. Where do the thermodynamic variants move it, and in which direction?

A reaction's native operator is read from its bounds:
    lb < 0 and ub > 0  -> '='      lb >= 0 and ub > 0 -> '>'
    lb < 0 and ub <= 0 -> '<'      otherwise          -> 'blocked'

Outputs (results/thermo_grades_fba/):
    implicit_directions.tsv       one row per core reaction
    implicit_summary.json         the aggregate tables
"""
from __future__ import annotations

import collections
import json
import os
import sys
from pathlib import Path

import pandas as pd

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
MSDB_CODE = Path(os.environ.get("MSDB_CODE", "/scratch/ctaylor/ModelSEEDDatabase"))
MODELS_DIR = ANALYSIS_DIR / "data" / "core_models_kegg2"
DATA = ANALYSIS_DIR / "results" / "thermo_grades_fba"
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(MSDB_CODE / "Scripts" / "Thermodynamics"))

VARIANTS = ["gc", "eq", "dgpms", "graded", "graded_trusted", "graded_heldout"]


def native_operator(lb: float, ub: float) -> str:
    if lb < 0 and ub > 0:
        return "="
    if lb >= 0 and ub > 0:
        return ">"
    if lb < 0 and ub <= 0:
        return "<"
    return "blocked"


def scan_models() -> tuple:
    """{rxn: Counter(operator -> n_models)} plus the global (model, rxn) mix."""
    from seed_annotation import seed_id
    per_rxn = collections.defaultdict(collections.Counter)
    mix = collections.Counter()
    n_models = 0
    for path in sorted(MODELS_DIR.glob("*.json")):
        d = json.load(open(path))
        n_models += 1
        seen = {}
        for r in d.get("reactions", []):
            s = seed_id(r)
            if s:
                seen[s] = native_operator(float(r.get("lower_bound", 0)),
                                          float(r.get("upper_bound", 0)))
        for s, o in seen.items():
            per_rxn[s][o] += 1
            mix[o] += 1
    return per_rxn, mix, n_models


def tecrdb_reference() -> dict:
    from analyze_graded_fba import tecrdb_reference as _ref
    return _ref()


def main() -> None:
    per_rxn, mix, n_models = scan_models()
    print(f"{n_models} models, {len(per_rxn)} distinct core reactions")

    total = sum(mix.values())
    print(f"\nnative bound mix over all {total} (model, reaction) pairs:")
    for o, c in mix.most_common():
        print(f"  {o:8s} {c:7d}  {c/total:6.1%}")

    consistency = collections.Counter()
    modal = {}
    for s, c in per_rxn.items():
        top, n = c.most_common(1)[0]
        modal[s] = top
        f = n / sum(c.values())
        consistency["unanimous" if f == 1 else
                    (">=90%" if f >= 0.9 else (">=50%" if f >= 0.5 else "split"))] += 1
    print("\nper-reaction consistency of the native direction across models:")
    for k in ("unanimous", ">=90%", ">=50%", "split"):
        print(f"  {k:10s} {consistency[k]:4d} of {len(per_rxn)}")

    modal_mix = collections.Counter(modal.values())
    print("\nmodal native operator per core reaction:")
    for o, c in modal_mix.most_common():
        print(f"  {o:8s} {c:4d}  {c/len(modal):6.1%}")

    ref = tecrdb_reference()
    cov = pd.read_csv(DATA / "rxn_source_coverage.csv", low_memory=False).set_index("rxn_id")
    scoreable = [r for r in modal if r in ref and ref[r][1] == "stereo_exact"]
    acc = {}
    n = ok = 0
    mismatch = collections.Counter()
    for r in scoreable:
        n += 1
        if modal[r] == ref[r][0]:
            ok += 1
        else:
            mismatch[f"{modal[r]} -> {ref[r][0]}"] += 1
    acc["implicit"] = {"n": n, "correct": ok, "accuracy": ok / n if n else None}
    print(f"\ndirection accuracy vs the experimental reference (n={n} core reactions):")
    print(f"  {'implicit (native bounds)':26s} {ok:3d}/{n:3d}  {ok/n:6.1%}")
    for v in VARIANTS:
        n2 = ok2 = 0
        for r in scoreable:
            if r in cov.index and cov.at[r, f"has_{v}"]:
                n2 += 1
                ok2 += (cov.at[r, f"op_{v}"] == ref[r][0])
        acc[v] = {"n": n2, "correct": ok2, "accuracy": ok2 / n2 if n2 else None}
        print(f"  {v:26s} {ok2:3d}/{n2:3d}  {ok2/n2:6.1%}")
    print(f"  implicit mismatches (native -> experiment): {dict(mismatch)}")

    agree = {}
    for v in VARIANTS:
        n3 = a = 0
        for r in modal:
            if r in cov.index and cov.at[r, f"has_{v}"] and modal[r] != "blocked":
                n3 += 1
                a += (cov.at[r, f"op_{v}"] == modal[r])
        agree[v] = {"n": n3, "agree": a, "frac": a / n3 if n3 else None}
    print("\nagreement with the native bound (core reactions the variant covers):")
    for v, d in agree.items():
        print(f"  {v:26s} {d['agree']:3d}/{d['n']:3d}  {d['frac']:6.1%}")

    transitions = collections.Counter()
    for r in modal:
        if r in cov.index and cov.at[r, "has_graded"] and modal[r] != "blocked":
            transitions[f"{modal[r]} -> {cov.at[r, 'op_graded']}"] += 1
    print("\nnative -> graded transitions (core reactions):")
    for k, c in sorted(transitions.items(), key=lambda x: -x[1]):
        print(f"  {k:12s} {c:4d}")

    rows = []
    for r in sorted(modal):
        row = {"rxn": r, "n_core_models": sum(per_rxn[r].values()),
               "native_op": modal[r],
               "native_unanimous": len(per_rxn[r]) == 1,
               "tecrdb_direction": ref[r][0] if r in ref else None,
               "tecrdb_match_tier": ref[r][1] if r in ref else None}
        for v in VARIANTS:
            row[f"op_{v}"] = cov.at[r, f"op_{v}"] if (
                r in cov.index and cov.at[r, f"has_{v}"]) else ""
        rows.append(row)
    pd.DataFrame(rows).to_csv(DATA / "implicit_directions.tsv", sep="\t", index=False)
    json.dump({"n_models": n_models, "n_core_reactions": len(per_rxn),
               "pairwise_mix": dict(mix), "consistency": dict(consistency),
               "modal_mix": dict(modal_mix), "direction_accuracy": acc,
               "implicit_mismatches": dict(mismatch),
               "agreement_with_native": agree,
               "native_to_graded": dict(transitions)},
              open(DATA / "implicit_summary.json", "w"), indent=1)
    print(f"\nwrote {DATA}/implicit_directions.tsv, implicit_summary.json")


if __name__ == "__main__":
    main()
