#!/usr/bin/env python3
"""Investigate WHY variant 3.1 (eQuilibrator reversibility index) breaks so many models.

Three questions:
  Q1. Which reactions are responsible? (single-reaction knock-in across all models)
  Q2. Does the panel's ~73% grower break-rate generalize to all ~3.5K growers?
  Q3. Does 3.1 flip any non-growers into growth?

Method for Q1: for each 3.1-changed reaction that appears in >=1 core model,
build a reversibility map = baseline with ONLY that one reaction set to its
3.1 direction, run FBA across the models that contain it, and record which
baseline-growers lose growth. The reactions whose single change breaks the
most growers are the culprits. A greedy union over these single-reaction
break-sets shows how much of 3.1's total breakage a handful of reactions
explain.

Outputs JSON to results/variant_31_breakage.json and prints a summary.
Reads (no recompute) where possible:
  - site/data/baseline.json                  (baseline rev map)
  - site/data/variants/3.1.json              (3.1 diffs)
  - site/data/all_models_baseline_fba.json   (baseline grower status)
  - site/data/all_models_variant_fba__3.1.json (full-3.1 grower status)
  - site/data/all_models_rxnsets.json        (which models contain each rxn)
"""
from __future__ import annotations
import json
import os
import sys
import time
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
ROOT = SCRIPTS.parent
SITE = ROOT / "site" / "data"
MSDB = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase"))
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(MSDB / "Libs" / "Python"))

import growth_heuristics as gh  # noqa: E402


def load_json(p):
    return json.loads(Path(p).read_text())


def main():
    t0 = time.time()
    baseline = load_json(SITE / "baseline.json")["map"]
    diffs = {d["rxn"]: (d["base"], d["new"]) for d in load_json(SITE / "variants/3.1.json")["diffs"]}

    base_fba = {r["model_id"]: r for r in load_json(SITE / "all_models_baseline_fba.json")}
    v31_fba = {r["model_id"]: r for r in load_json(SITE / "all_models_variant_fba__3.1.json")}
    rxnsets = load_json(SITE / "all_models_rxnsets.json")  # {model_id: [rxn,...]}

    base_growers = {m for m, r in base_fba.items() if r.get("grows")}
    # Full-3.1 broken set: baseline grower that 3.1 turns into a non-grower.
    # v31_fba only contains models with >=1 changed rxn; models absent kept baseline status.
    def grows_under_31(m):
        if m in v31_fba:
            return bool(v31_fba[m].get("grows"))
        return bool(base_fba.get(m, {}).get("grows"))
    full_broken = {m for m in base_growers if not grows_under_31(m)}
    full_rescued = {m for m in base_fba if (not base_fba[m].get("grows")) and grows_under_31(m)}

    print(f"baseline growers: {len(base_growers)}")
    print(f"3.1 full breakage (grower->non-grower): {len(full_broken)}")
    print(f"3.1 full rescue   (non-grower->grower): {len(full_rescued)}")

    # Candidate reactions: changed by 3.1 AND present in >=1 model.
    rxn_to_models = {}
    for m, rxns in rxnsets.items():
        for r in rxns:
            if r in diffs:
                rxn_to_models.setdefault(r, set()).add(m)
    candidates = sorted(rxn_to_models, key=lambda r: -len(rxn_to_models[r]))
    print(f"candidate reactions (changed & present): {len(candidates)}")

    # Single-reaction knock-in: baseline + only this reaction's 3.1 dir.
    per_rxn = {}
    for r in candidates:
        b, n = diffs[r]
        models_with = sorted(rxn_to_models[r] & base_growers)  # only growers can break
        if not models_with:
            per_rxn[r] = {"base": b, "new": n, "n_models_with_rxn": len(rxn_to_models[r]),
                          "n_growers_with_rxn": 0, "break_set": [], "n_break": 0}
            continue
        eff = dict(baseline)
        eff[r] = n
        res = gh.run_panel(models_with, reversibility_map=eff, baseline_map=None, n_workers=64)
        broke = sorted(m for m in (rr["model_id"] for rr in res if not rr.get("grows")))
        # confirm these were growers (they are, models_with is grower-filtered)
        per_rxn[r] = {
            "base": b, "new": n,
            "n_models_with_rxn": len(rxn_to_models[r]),
            "n_growers_with_rxn": len(models_with),
            "n_break": len(broke),
            "break_set": broke,
        }
        print(f"  {r:12s} {b}->{n}  growers_with={len(models_with):5d}  "
              f"single-rxn breaks={len(broke):5d}  ({time.time()-t0:.0f}s)", flush=True)

    # Greedy cover of full_broken using single-reaction break-sets.
    remaining = set(full_broken)
    greedy = []
    sets = {r: set(per_rxn[r]["break_set"]) & full_broken for r in candidates}
    while remaining:
        best = max(candidates, key=lambda r: len(sets[r] & remaining))
        gain = len(sets[best] & remaining)
        if gain == 0:
            break
        greedy.append({"rxn": best, "marginal": gain,
                       "cumulative": len(full_broken) - len(remaining) + gain})
        remaining -= sets[best]
    explained = len(full_broken) - len(remaining)

    print(f"\nGreedy single-reaction cover of the {len(full_broken)} broken growers:")
    for g in greedy[:10]:
        b, n = diffs[g["rxn"]]
        print(f"  +{g['marginal']:5d}  cum {g['cumulative']:5d}  {g['rxn']} ({b}->{n})")
    print(f"explained by single-reaction knock-ins: {explained}/{len(full_broken)} "
          f"= {100*explained/max(1,len(full_broken)):.1f}%")
    print(f"unexplained (combination-only breaks): {len(remaining)}")

    out = {
        "baseline_growers": len(base_growers),
        "full_broken": len(full_broken),
        "full_rescued": len(full_rescued),
        "rescued_models": sorted(full_rescued),
        "candidates": [
            {"rxn": r, **{k: v for k, v in per_rxn[r].items() if k != "break_set"}}
            for r in candidates
        ],
        "greedy_cover": greedy,
        "explained_by_singles": explained,
        "unexplained": len(remaining),
        "elapsed_s": round(time.time() - t0, 1),
    }
    (ROOT / "results" / "variant_31_breakage.json").write_text(json.dumps(out, indent=2))
    print(f"\nwrote results/variant_31_breakage.json  ({out['elapsed_s']}s)")


if __name__ == "__main__":
    main()
