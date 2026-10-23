"""Classify the convention of each compound's top-level `deltag` scalar.

Convention A (MFAToolkit / Chris Henry 2010): H contribution left IN,  dfG(H+) = -9.5 kcal/mol
Convention B (Alberty / eQuilibrator):        H contribution removed,  dfG(H+) =  0.0 kcal/mol
Relationship:  dG_B = dG_A + n_H * 9.539

`thermodynamics["Group contribution"]` is guaranteed Convention A, so it is the
reference against which the top-level scalar is classified.
"""
import glob, json, re, sys
from collections import Counter

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
import _paths as _p

SENTINEL = 10000000
SHIFT = 9.539
BIOCHEM = str(_p.BIOCHEM)


def n_hydrogen(formula):
    if not formula or formula in ("null", "None"):
        return None
    # reject formulas with R-groups / wildcards
    if re.search(r"[RXZ\*]", formula):
        return None
    m = re.search(r"H(\d*)(?![a-z])", formula)
    if not m:
        return 0
    return int(m.group(1)) if m.group(1) else 1


rows = []
for f in sorted(glob.glob(f"{BIOCHEM}/compound_*.json")):
    rows.extend(json.load(open(f)))

tally = Counter()
examples = {"A": [], "B": [], "neither": []}
for c in rows:
    dg = c.get("deltag")
    th = c.get("thermodynamics") or {}
    gc = th.get("Group contribution")
    if dg in (None, SENTINEL) or not gc or gc[0] in (None, SENTINEL):
        tally["unclassifiable (missing deltag or GC)"] += 1
        continue
    nh = n_hydrogen(c.get("formula"))
    if nh is None:
        tally["unclassifiable (no usable formula)"] += 1
        continue
    d_A = abs(dg - gc[0])
    d_B = abs(dg - (gc[0] + nh * SHIFT))
    if d_A < 0.5 and d_B < 0.5:
        tally["ambiguous (n_H=0, conventions coincide)"] += 1
    elif d_A < 0.5:
        tally["Convention A"] += 1
        if len(examples["A"]) < 5:
            examples["A"].append((c["id"], c["name"][:24], c["formula"], dg, gc[0]))
    elif d_B < 0.5:
        tally["Convention B"] += 1
        if len(examples["B"]) < 5:
            examples["B"].append((c["id"], c["name"][:24], c["formula"], dg, gc[0]))
    else:
        tally["neither (drifted)"] += 1
        if len(examples["neither"]) < 5:
            examples["neither"].append((c["id"], c["name"][:24], c["formula"], dg, gc[0]))

print(f"total compounds: {len(rows)}")
for k, v in tally.most_common():
    print(f"  {k:46s} {v:7d}")

classified = tally["Convention A"] + tally["Convention B"]
if classified:
    print(f"\nof the {classified} unambiguously classifiable (n_H>0, both values present):")
    print(f"  Convention A : {tally['Convention A']:6d}  ({100*tally['Convention A']/classified:.1f}%)")
    print(f"  Convention B : {tally['Convention B']:6d}  ({100*tally['Convention B']/classified:.1f}%)")

for k in ("A", "B", "neither"):
    if examples[k]:
        print(f"\nexamples — {k}:")
        for e in examples[k]:
            print(f"   {e[0]} {e[1]:24s} {e[2]:16s} deltag={e[3]:10.3f} GC={e[4]:10.3f}")
