#!/usr/bin/env python3
"""Which reactions does the TMFA constraint hard-lock, and in which direction?

In FullThermoPkg a reaction's drG is
    drG = sum_i s_i * (dGf_i + R*T*ln c_i + z_i*F*psi_compartment + err_i)
and the revbin coupling forces sign(drG) to oppose the flux direction.

The only freedom available to satisfy that sign is the log-concentration range
(default 1e-6 .. 0.02 M, i.e. R*T*ln(0.02/1e-6) = 24.5 kJ/mol per metabolite)
plus the per-metabolite dG error slack (default +/-5 kJ/mol). If the standard
drG exceeds that total swing, the direction is locked and no flux solution can
override it.
"""
import glob, json, math, sys
from pathlib import Path

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
import _paths as _p

BIOCHEM = str(_p.BIOCHEM)
MODEL = sys.argv[1] if len(sys.argv) > 1 else str(_p.DEFAULT_MODEL)
SENTINEL = 10000000
KCAL_TO_KJ = 4.184
R_KJ = 8.314462618 / 1000.0
T = 298.0
RT = R_KJ * T
CONC_SWING = RT * (math.log(0.02) - math.log(1e-6))   # 24.5 kJ/mol per metabolite
ERR = 5.0                                             # kJ/mol per metabolite

cpds, rxn_dg = {}, {}
for f in sorted(glob.glob(f"{BIOCHEM}/compound_*.json")):
    for c in json.load(open(f)):
        cpds[c["id"]] = c
for f in sorted(glob.glob(f"{BIOCHEM}/reaction_*.json")):
    for r in json.load(open(f)):
        v = r.get("deltag")
        try:
            v = float(v)
        except (TypeError, ValueError):
            v = None
        rxn_dg[r["id"]] = None if v is None or v == SENTINEL else v * KCAL_TO_KJ

model = json.load(open(MODEL))
locked_asis, locked_fixed, wrong_asis, wrong_fixed = [], [], [], []
n = 0
for rxn in model["reactions"]:
    if rxn["id"][:3] in ("EX_", "SK_", "DM_") or rxn["id"].startswith("bio"):
        continue
    dg_kcal, swing, ok = 0.0, 0.0, True
    for met_id, s in rxn["metabolites"].items():
        msid = met_id.split("_")[0]
        c = cpds.get(msid)
        if c is None or c.get("deltag") in (None, SENTINEL):
            ok = False
            break
        dg_kcal += s * c["deltag"]
        if msid not in ("cpd00001", "cpd00067"):   # water fixed, H+ pinned at pH 7
            swing += abs(s) * CONC_SWING
        swing += abs(s) * ERR
    if not ok:
        continue
    n += 1
    dg_fixed = dg_kcal * KCAL_TO_KJ          # correct kcal -> kJ
    dg_asis = dg_kcal / KCAL_TO_KJ           # what the shipped code computes
    ref = rxn_dg.get(rxn["id"].split("_")[0])
    for dg, locked, wrong in ((dg_asis, locked_asis, wrong_asis),
                              (dg_fixed, locked_fixed, wrong_fixed)):
        if abs(dg) > swing:
            rec = {"id": rxn["id"], "dg": dg, "swing": swing, "ref": ref,
                   "name": rxn.get("name", "")[:38]}
            locked.append(rec)
            if ref is not None and dg * ref < 0:
                wrong.append(rec)

print(f"model {Path(MODEL).stem}: {n} internal reactions with complete compound deltag\n")
for label, locked, wrong in (("as shipped  (deltag / 4.184)", locked_asis, wrong_asis),
                             ("unit-fixed  (deltag * 4.184)", locked_fixed, wrong_fixed)):
    print(f"{label}")
    print(f"   direction hard-locked by thermodynamics : {len(locked):3d} / {n} "
          f"({100*len(locked)/n:.0f}%)")
    print(f"   ...of which locked OPPOSITE to the DB's own reaction-level drG : {len(wrong)}")
    if wrong:
        for r in sorted(wrong, key=lambda x: -abs(x["dg"]))[:8]:
            print(f"      {r['id']:14s} {r['name']:38s} tFBA drG={r['dg']:9.1f}  "
                  f"swing=+/-{r['swing']:6.1f}  DB drG={r['ref']:7.1f}")
    print()
