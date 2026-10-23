#!/usr/bin/env python3
"""FullThermoPkg end-to-end with NO API shims -- stock modelseedpy 0.4.2
against a legacy-layout ModelSEEDDatabase (the pre-2022-05-20 layout).

The only patch applied is the speed one (closed-form max_energy_magnitude in
place of 2 MILP solves per reaction); it changes no numbers.
"""
import sys, time
from pathlib import Path
import cobra
from cobra.io import load_json_model
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager

sys.path.insert(0, str(Path(__file__).parent))
import fast_thermo
fast_thermo.install()

LEGACY = str(_p.LEGACY_MSDB)
MEDIA = _p.MEDIA
MODEL = str(_p.DEFAULT_MODEL)

model = load_json_model(MODEL)
model.solver = "glpk"
MINIMAL = {"cpd00027": -10.0, "cpd00007": -20.0}
OPEN = {"cpd00001","cpd00067","cpd00009","cpd00013","cpd00048","cpd00205","cpd00971",
        "cpd00099","cpd00030","cpd00034","cpd00058","cpd00149","cpd10515","cpd10516",
        "cpd00254","cpd00063","cpd00011"}
for r in model.reactions:
    if r.id.startswith("EX_"):
        m = list(r.metabolites)
        if len(m) == 1:
            c = m[0].id.split("_")[0]
            r.upper_bound = 1000.0
            r.lower_bound = MINIMAL.get(c, -1000.0 if c in OPEN else 0.0)

base = model.optimize()
print(f"plain FBA growth : {base.objective_value:.4f}  ({base.status})")

saved = model.objective
t0 = time.time()
pkgmgr = MSPackageManager.get_pkg_mgr(model)
ft = pkgmgr.getpkg("FullThermoPkg")
ft.build_package({"modelseed_path": LEGACY}, verbose=False)
print(f"build_package OK : {time.time()-t0:.1f}s   <-- STOCK CODE, NO SHIM")

sp = pkgmgr.getpkg("SimpleThermoPkg")
print(f"  potc constraints (metabolites priced) : {len(ft.constraints.get('potc', {}))}")
print(f"  thermo constraints (reactions coupled): {len(sp.constraints.get('thermo', {}))}")
print(f"  logconc vars: {len(ft.variables.get('logconc', {}))}   "
      f"dgerr vars: {len(ft.variables.get('dgerr', {}))}   "
      f"revbin vars: {len(pkgmgr.getpkg('RevBinPkg').variables.get('revbin', {}))}")
print(f"  objective after build_package == original? "
      f"{str(model.objective.expression) == str(saved.expression)}")

# ---- make the implied potential bounds explicit, then solve -----------------
from scipy.constants import R, kilo

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
import _paths as _p
RT = R / kilo * ft.parameters["temperature"]
for met_id, potc in ft.constraints.get("potc", {}).items():
    pot = sp.variables["potential"].get(met_id)
    if pot is None:
        continue
    lo = hi = potc.lb
    lc = ft.variables["logconc"].get(met_id)
    if lc is not None:
        lo += RT * lc.lb; hi += RT * lc.ub
    de = ft.variables["dgerr"].get(met_id)
    if de is not None:
        lo += de.lb; hi += de.ub
    pot.lb, pot.ub = lo, hi

try:
    model.solver.configuration.timeout = 600
except Exception as e:
    print("  (no solver timeout support:", e, ")")
t0 = time.time()
sol = model.optimize()
print(f"\ntFBA solve: status={sol.status}  growth={sol.objective_value}  "
      f"({time.time()-t0:.0f}s)")
if sol.status == "optimal":
    changed = [(r.id, base.fluxes[r.id], sol.fluxes[r.id]) for r in model.reactions
               if abs(base.fluxes[r.id] - sol.fluxes[r.id]) > 1e-6]
    print(f"  growth vs FBA : {100*sol.objective_value/base.objective_value:.1f}%")
    print(f"  reactions with changed flux : {len(changed)}")
    import math
    print(f"\n  example metabolite concentrations the solver chose (M):")
    for mid in ["cpd00002_c0", "cpd00027_c0", "cpd00020_c0"]:
        v = ft.variables["logconc"].get(mid)
        if v is not None and v.primal is not None:
            print(f"    {mid:14s} {math.exp(v.primal):.3e}")
    print(f"\n  example reaction drG (kJ/mol) from the solved potentials:")
    for rid in ["rxn00777_c0", "rxn00216_c0", "rxn00459_c0"]:
        if rid not in [r.id for r in model.reactions]:
            continue
        rx = model.reactions.get_by_id(rid)
        g = sum(s * sp.variables["potential"][m.id].primal
                for m, s in rx.metabolites.items()
                if m.id in sp.variables["potential"])
        print(f"    {rid:14s} {g:10.2f}   flux={sol.fluxes[rid]:8.3f}")
