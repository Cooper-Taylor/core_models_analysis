#!/usr/bin/env python3
"""Cross-check the drG that FullThermoPkg would build from compound deltag
against the reaction-level drG values the database already ships.

FullThermoPkg never reads a reaction drG; it sums compound potentials. If that
sum disagrees with the DB's own reaction-level numbers, the compound scalars
are not usable as a tFBA input.
"""
import glob, json, statistics, sys
from pathlib import Path

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
import _paths as _p

BIOCHEM = str(_p.BIOCHEM)
MODEL = sys.argv[1] if len(sys.argv) > 1 else str(_p.DEFAULT_MODEL)
SENTINEL = 10000000
KCAL_TO_KJ = 4.184

cpds = {}
for f in sorted(glob.glob(f"{BIOCHEM}/compound_*.json")):
    for c in json.load(open(f)):
        cpds[c["id"]] = c
rxns = {}
for f in sorted(glob.glob(f"{BIOCHEM}/reaction_*.json")):
    for r in json.load(open(f)):
        rxns[r["id"]] = r

model = json.load(open(MODEL))
recs = []
for rxn in model["reactions"]:
    rid = rxn["id"].split("_")[0]
    if not rid.startswith("rxn"):
        continue
    dbr = rxns.get(rid)
    if dbr is None:
        continue
    # drG summed from compound scalars -- exactly what FullThermoPkg does
    total, ok = 0.0, True
    for met_id, s in rxn["metabolites"].items():
        c = cpds.get(met_id.split("_")[0])
        if c is None or c.get("deltag") in (None, SENTINEL):
            ok = False
            break
        total += s * c["deltag"]
    if not ok:
        continue
    th = dbr.get("thermodynamics") or {}
    rec = {"id": rid, "from_compounds": total * KCAL_TO_KJ}
    rl = dbr.get("deltag")
    rec["reaction_level"] = rl * KCAL_TO_KJ if rl not in (None, SENTINEL) else None
    for src in ("Group contribution", "eQuilibrator"):
        v = th.get(src)
        rec[src] = v[0] * KCAL_TO_KJ if v and v[0] not in (None, SENTINEL) else None
    recs.append(rec)

print(f"reactions compared: {len(recs)}\n")


def summarize(key):
    pairs = [(r["from_compounds"], r[key]) for r in recs if r.get(key) is not None]
    if not pairs:
        print(f"{key}: no overlap")
        return
    diffs = [abs(a - b) for a, b in pairs]
    within = sum(1 for d in diffs if d < 5)
    print(f"vs {key:22s} n={len(pairs):4d}  median|diff|={statistics.median(diffs):8.1f} kJ/mol  "
          f"max={max(diffs):9.1f}  within 5 kJ/mol: {within} ({100*within/len(pairs):.0f}%)")


for k in ("reaction_level", "Group contribution", "eQuilibrator"):
    summarize(k)

mag = [abs(r["from_compounds"]) for r in recs]
print(f"\n|drG| summed from compound scalars: median={statistics.median(mag):.1f} "
      f"max={max(mag):.1f} kJ/mol")
over = sum(1 for m in mag if m > 200)
print(f"reactions with |drG| > 200 kJ/mol (physically implausible): {over}/{len(mag)} "
      f"({100*over/len(mag):.0f}%)")

eq = [abs(r["eQuilibrator"]) for r in recs if r.get("eQuilibrator") is not None]
if eq:
    print(f"\nfor comparison, |drG| from the DB's reaction-level eQuilibrator values: "
          f"median={statistics.median(eq):.1f} max={max(eq):.1f} kJ/mol")
    print(f"  of those, > 200 kJ/mol: {sum(1 for m in eq if m>200)}/{len(eq)}")

worst = sorted(recs, key=lambda r: -abs(r["from_compounds"]))[:10]
print(f"\n{'reaction':12s} {'from compounds':>15s} {'rxn-level':>11s} {'GC':>10s} {'eQuilibrator':>13s}")
for r in worst:
    f = lambda v: f"{v:10.1f}" if v is not None else "         -"
    print(f"{r['id']:12s} {r['from_compounds']:15.1f} {f(r['reaction_level'])} "
          f"{f(r['Group contribution'])} {f(r['eQuilibrator'])}")
