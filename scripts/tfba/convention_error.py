#!/usr/bin/env python3
"""Per-reaction error introduced by the mixed-convention `deltag` scalar.

FullThermoPkg models the proton explicitly (cpd00067 at 1e-7 M, and the DB
gives it deltag = 0), so it needs every compound in Convention B.

For a Convention-A compound the stored value is low by n_H * 9.539 kcal/mol,
so the error contributed to a reaction's drG is

    mix_error = -9.539 * sum_{i in A} s_i * n_H_i      [kcal/mol]

Only confidently-A-classified compounds are counted, which makes this a LOWER
BOUND on the true convention error: drifted compounds add more.
"""
import glob, json, re, sys
from collections import Counter
from pathlib import Path

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
import _paths as _p

BIOCHEM = str(_p.BIOCHEM)
MODEL = sys.argv[1] if len(sys.argv) > 1 else str(_p.DEFAULT_MODEL)
SENTINEL = 10000000
SHIFT = 9.539
KCAL_TO_KJ = 4.184


def n_hydrogen(formula):
    if not formula or re.search(r"[RXZ\*]", str(formula)):
        return None
    m = re.search(r"H(\d*)(?![a-z])", str(formula))
    if not m:
        return 0
    return int(m.group(1)) if m.group(1) else 1


cpds = {}
for f in sorted(glob.glob(f"{BIOCHEM}/compound_*.json")):
    for c in json.load(open(f)):
        cpds[c["id"]] = c


def classify(c):
    dg = c.get("deltag")
    gc = (c.get("thermodynamics") or {}).get("Group contribution")
    nh = n_hydrogen(c.get("formula"))
    if dg in (None, SENTINEL) or not gc or gc[0] in (None, SENTINEL) or nh is None:
        return "unknown", nh
    if nh == 0:
        return "A=B", nh
    if abs(dg - gc[0]) < 0.5:
        return "A", nh
    if abs(dg - (gc[0] + nh * SHIFT)) < 1.5:
        return "B", nh
    return "drifted", nh


model = json.load(open(MODEL))
conv_tally = Counter()
rows = []
for rxn in model["reactions"]:
    if rxn["id"][:3] in ("EX_", "SK_", "DM_") or rxn["id"].startswith("bio"):
        continue
    a_h_flow = 0.0          # sum of s_i * n_H_i over Convention-A compounds
    drifted = 0
    unknown = 0
    dg_mixed = 0.0
    complete = True
    for met_id, s in rxn["metabolites"].items():
        msid = met_id.split("_")[0]
        c = cpds.get(msid)
        if c is None or c.get("deltag") in (None, SENTINEL):
            complete = False
            continue
        dg_mixed += s * c["deltag"]
        conv, nh = classify(c)
        conv_tally[conv] += 1
        if conv == "A" and nh:
            a_h_flow += s * nh
        elif conv == "drifted":
            drifted += 1
        elif conv == "unknown":
            unknown += 1
    err = -SHIFT * a_h_flow
    rows.append({
        "id": rxn["id"], "name": rxn.get("name", "")[:40],
        "dg_mixed_kJ": dg_mixed * KCAL_TO_KJ,
        "mix_error_kJ": err * KCAL_TO_KJ,
        "drifted_cpds": drifted, "unknown_cpds": unknown,
        "complete": complete,
    })

usable = [r for r in rows if r["complete"]]
print(f"model: {Path(MODEL).stem}")
print(f"internal reactions: {len(rows)}   with all deltag present: {len(usable)}\n")

print("compound-convention occurrences across those reactions:")
for k, v in conv_tally.most_common():
    print(f"   {k:10s} {v:6d}")

big = [r for r in usable if abs(r["mix_error_kJ"]) > 1.0]
print(f"\nreactions with a convention-mix error > 1 kJ/mol: {len(big)} / {len(usable)}"
      f"  ({100*len(big)/max(len(usable),1):.0f}%)")
flipped = [r for r in big
           if r["dg_mixed_kJ"] * (r["dg_mixed_kJ"] - r["mix_error_kJ"]) < 0]
print(f"reactions where removing the error FLIPS the sign of drG: {len(flipped)}")

big.sort(key=lambda r: -abs(r["mix_error_kJ"]))
print(f"\n{'reaction':14s} {'name':40s} {'drG(mixed)':>11s} {'error':>9s} {'drG(corr)':>10s}")
for r in big[:15]:
    print(f"{r['id']:14s} {r['name']:40s} {r['dg_mixed_kJ']:11.1f} "
          f"{r['mix_error_kJ']:9.1f} {r['dg_mixed_kJ']-r['mix_error_kJ']:10.1f}")

out = _p.RESULTS / "convention_error.json"
json.dump(rows, open(out, "w"), indent=2)
print(f"\nwrote {out}")
