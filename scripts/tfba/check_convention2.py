"""Classify convention by the per-hydrogen shift s = (deltag - GC) / n_H.

s ~ 0.0   -> Convention A (GC/MFAToolkit, dfG(H+) = -9.5 kcal/mol)
s ~ 9.5   -> Convention B (Alberty/eQuilibrator, dfG(H+) = 0)
anything else -> drifted / indeterminate
"""
import glob, json, re
from collections import Counter

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
import _paths as _p

SENTINEL = 10000000
BIOCHEM = str(_p.BIOCHEM)


def n_hydrogen(formula):
    if not formula or formula in ("null", "None") or re.search(r"[RXZ\*]", formula):
        return None
    m = re.search(r"H(\d*)(?![a-z])", formula)
    if not m:
        return 0
    return int(m.group(1)) if m.group(1) else 1


rows = []
for f in sorted(glob.glob(f"{BIOCHEM}/compound_*.json")):
    rows.extend(json.load(open(f)))

bins = Counter()
shifts = []
for c in rows:
    dg, th = c.get("deltag"), (c.get("thermodynamics") or {})
    gc = th.get("Group contribution")
    if dg in (None, SENTINEL) or not gc or gc[0] in (None, SENTINEL):
        continue
    nh = n_hydrogen(c.get("formula"))
    if not nh:          # None or 0
        continue
    s = (dg - gc[0]) / nh
    shifts.append((s, c["id"], c["name"], c["formula"], dg, gc[0]))
    if abs(s) < 0.2:
        bins["A  (s~0.0, Convention A)"] += 1
    elif 9.2 <= s <= 9.8:
        bins["B  (s~9.5, Convention B)"] += 1
    elif 8.5 <= s < 9.2:
        bins["B' (s~8.5-9.2, partial transform)"] += 1
    else:
        bins["?  (other shift)"] += 1

tot = len(shifts)
print(f"compounds with both deltag and GC and n_H>0: {tot}\n")
for k, v in sorted(bins.items()):
    print(f"  {k:40s} {v:6d}  ({100*v/tot:5.1f}%)")

print("\n--- key central-metabolism compounds ---")
want = ["cpd00002 ATP", "cpd00008 ADP", "cpd00038 GTP", "cpd00031 GDP",
        "cpd00003 NAD", "cpd00004 NADH", "cpd00001 H2O", "cpd00067 H+",
        "cpd00009 Pi", "cpd00023 Glu", "cpd00027 Glc", "cpd00020 Pyr"]
idx = {r[1]: r for r in shifts}
raw = {c["id"]: c for c in rows}
print(f"{'id':10s} {'name':22s} {'formula':16s} {'nH':>3s} {'deltag':>9s} {'GC':>9s} {'shift/H':>8s}  conv")
for w in want:
    cid = w.split()[0]
    c = raw.get(cid)
    if not c:
        continue
    gc = (c.get("thermodynamics") or {}).get("Group contribution")
    nh = n_hydrogen(c.get("formula"))
    dg = c.get("deltag")
    if gc is None or gc[0] in (None, SENTINEL) or dg in (None, SENTINEL):
        print(f"{cid:10s} {c['name'][:22]:22s} {str(c.get('formula'))[:16]:16s} -- no GC/deltag --")
        continue
    s = (dg - gc[0]) / nh if nh else float("nan")
    conv = "A" if abs(s) < 0.2 else ("B" if 8.5 <= s <= 9.8 else "?")
    if nh == 0:
        conv = "A=B"
    print(f"{cid:10s} {c['name'][:22]:22s} {str(c.get('formula'))[:16]:16s} {nh:3d} "
          f"{dg:9.2f} {gc[0]:9.2f} {s:8.2f}  {conv}")
