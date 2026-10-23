#!/usr/bin/env python3
"""Database-wide: how often does the compound-summed drG (the only thing
FullThermoPkg can build) disagree with the reaction-level drG the DB ships?

Also reports the overlap with the reaction pool actually used by the
core_models_kegg2 panel.
"""
import glob, json, statistics
from collections import Counter
from pathlib import Path

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
import _paths as _p

BIOCHEM = str(_p.BIOCHEM)
MODELS = _p.MODELS
SENTINEL = 10000000
KCAL_TO_KJ = 4.184

cpd_dg = {}
for f in sorted(glob.glob(f"{BIOCHEM}/compound_*.json")):
    for c in json.load(open(f)):
        v = c.get("deltag")
        cpd_dg[c["id"]] = None if v in (None, SENTINEL) else v

results = {}
for f in sorted(glob.glob(f"{BIOCHEM}/reaction_*.json")):
    for r in json.load(open(f)):
        stoich = r.get("stoichiometry")
        if not stoich:
            continue
        total, ok = 0.0, True
        for term in stoich:
            dg = cpd_dg.get(term.get("compound"))
            if dg is None:
                ok = False; break
            total += float(term["coefficient"]) * dg
        if not ok:
            continue
        rl = r.get("deltag")
        if rl in (None, SENTINEL):
            continue
        results[r["id"]] = (total * KCAL_TO_KJ, rl * KCAL_TO_KJ)

diffs = {k: abs(a - b) for k, (a, b) in results.items()}
n = len(diffs)
print(f"MSDB reactions with both a compound-summable drG and a reaction-level drG: {n}\n")
bands = [(5, "< 5"), (50, "5 - 50"), (200, "50 - 200"), (1000, "200 - 1000"), (1e9, "> 1000")]
prev = 0
for hi, label in bands:
    c = sum(1 for d in diffs.values() if prev <= d < hi)
    print(f"  |compound-sum - reaction-level|  {label:>12s} kJ/mol : {c:6d}  ({100*c/n:5.1f}%)")
    prev = hi
bad = {k for k, d in diffs.items() if d > 50}
print(f"\ncorrupted (> 50 kJ/mol disagreement): {len(bad)} ({100*len(bad)/n:.1f}%)")

# which of those are actually used by the core model panel?
pool = Counter()
files = sorted(MODELS.glob("GCF_*.json"))
sample = files[:400]
for mf in sample:
    try:
        m = json.load(open(mf))
    except Exception:
        continue
    for rxn in m.get("reactions", []):
        rid = rxn["id"].split("_")[0]
        if rid.startswith("rxn"):
            pool[rid] += 1
print(f"\ncore-model panel sample: {len(sample)} models, {len(pool)} distinct rxn ids")
hit = {r: pool[r] for r in pool if r in bad}
print(f"  distinct corrupted reactions present in the panel: {len(hit)}")
if hit:
    frac = [100 * v / len(sample) for v in hit.values()]
    print(f"  they appear in a median of {statistics.median(frac):.0f}% of models "
          f"(max {max(frac):.0f}%)")
    models_hit = set()
    for mf in sample:
        try:
            m = json.load(open(mf))
        except Exception:
            continue
        ids = {rxn["id"].split("_")[0] for rxn in m.get("reactions", [])}
        if ids & bad:
            models_hit.add(mf.name)
    print(f"  models containing >=1 corrupted reaction: {len(models_hit)}/{len(sample)} "
          f"({100*len(models_hit)/len(sample):.0f}%)")
    top = sorted(hit.items(), key=lambda kv: -kv[1])[:10]
    print(f"\n  {'reaction':12s} {'models':>7s} {'cmp-sum kJ':>12s} {'rxn-level':>10s}")
    for rid, cnt in top:
        a, b = results[rid]
        print(f"  {rid:12s} {cnt:7d} {a:12.1f} {b:10.1f}")

json.dump({"corrupted": sorted(bad)},
          open(str(_p.RESULTS / "dbwide_corrupted.json"), "w"))
