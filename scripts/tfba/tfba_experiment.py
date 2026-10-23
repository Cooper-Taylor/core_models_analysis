#!/usr/bin/env python3
"""Run FBA vs thermodynamic FBA (TMFA) on a core model.

modelseedpy 0.4.2's FullThermoPkg cannot run against its own biochem loader,
so two compatibility shims are installed first:

  1. `ModelSEEDDatabase` (what `from_local` returns) has no `get_seed_compound`;
     that method lives on the older `ModelSEEDBiochem` class.
  2. `ModelSEEDCompound2` exposes `delta_g` / `delta_g_error`, but the package
     reads `.deltag` / `.deltagerr`.

A third switch controls the kcal->kJ conversion. ModelSEED delta_g is kcal/mol.
The package computes `deltag / calorie`; converting kcal to kJ is `* calorie`.
Setting the module-level `calorie` to 1/4.184 makes the package's own
expression produce the correct kJ/mol value, so both behaviours can be run
through identical code.
"""
import argparse, json, math, sys, time
from pathlib import Path

import cobra
from cobra.io import load_json_model
from modelseedpy.biochem import from_local
from modelseedpy.fbapkg import fullthermopkg
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
import _paths as _p

FAST = False
TIGHTEN = False
MSDB = str(_p.MSDB)
MEDIA = _p.MEDIA
SENTINEL = 10000000
KCAL_TO_KJ = 4.184


class _CpdAdapter:
    """Presents ModelSEEDCompound2 under the attribute names FullThermoPkg wants."""
    __slots__ = ("deltag", "deltagerr", "id")

    def __init__(self, cpd):
        self.id = cpd.id
        dg = cpd.delta_g
        self.deltag = SENTINEL if dg is None else dg
        err = getattr(cpd, "delta_g_error", None)
        self.deltagerr = SENTINEL if err is None else err


def install_db_shim():
    """Patch the class, not an instance: build_package calls from_local() itself
    via FBAHelper.get_modelseed_db_api, so an instance-level shim is discarded."""
    from modelseedpy.biochem.modelseed_biochem import ModelSEEDDatabase

    def get_seed_compound(self, seed_id):
        try:
            return _CpdAdapter(self.compounds.get_by_id(seed_id))
        except KeyError:
            return None

    ModelSEEDDatabase.get_seed_compound = get_seed_compound


def load_model(path, medium, solver="highs"):
    model = load_json_model(path)
    model.solver = solver
    media = {l.strip() for l in open(MEDIA) if l.strip()}
    for rxn in model.reactions:
        if not rxn.id.startswith("EX_"):
            continue
        mets = list(rxn.metabolites)
        if len(mets) != 1:
            continue
        cpd = mets[0].id.split("_")[0]
        rxn.upper_bound = 1000.0
        if medium == "complete":
            rxn.lower_bound = -1000.0 if cpd in media else 0.0
        else:  # glucose minimal, aerobic
            MINIMAL = {
                "cpd00027": -10.0,   # D-glucose  (the limiting carbon source)
                "cpd00007": -20.0,   # O2
                "cpd00001": -1000.0, "cpd00067": -1000.0,   # water, proton
                "cpd00009": -1000.0, "cpd00013": -1000.0,   # phosphate, ammonia
                "cpd00048": -1000.0, "cpd00205": -1000.0,   # sulfate, potassium
                "cpd00971": -1000.0, "cpd00099": -1000.0,   # sodium, chloride
                "cpd00030": -1000.0, "cpd00034": -1000.0,   # Mn2+, Zn2+
                "cpd00058": -1000.0, "cpd00149": -1000.0,   # Cu2+, Co2+
                "cpd10515": -1000.0, "cpd10516": -1000.0,   # Fe2+, Fe3+
                "cpd00254": -1000.0, "cpd00063": -1000.0,   # Mg, Ca
                "cpd00011": -1000.0,                        # CO2
            }
            rxn.lower_bound = MINIMAL.get(cpd, 0.0)
    return model


def thermo_coverage(model, db):
    """How many metabolites carry a usable delta_g, and how many reactions are
    therefore actually constrained vs silently left free."""
    have, missing, no_msid = set(), set(), set()
    for met in model.metabolites:
        msid = met.id.split("_")[0]
        if not msid.startswith("cpd"):
            no_msid.add(met.id); continue
        try:
            cpd = db.compounds.get_by_id(msid)
        except KeyError:
            missing.add(met.id); continue
        if cpd.delta_g is None or cpd.delta_g == SENTINEL:
            missing.add(met.id)
        else:
            have.add(met.id)
    internal = [r for r in model.reactions if r.id[:3] not in ("EX_", "SK_", "DM_")]
    fully, partly = [], []
    for r in internal:
        bad = [m.id for m in r.metabolites if m.id in missing or m.id in no_msid]
        (partly if bad else fully).append(r.id)
    return {
        "metabolites_total": len(model.metabolites),
        "metabolites_with_dg": len(have),
        "metabolites_without_dg": len(missing) + len(no_msid),
        "internal_reactions": len(internal),
        "reactions_all_mets_have_dg": len(fully),
        "reactions_with_a_free_potential": len(partly),
        "unconstrained_reaction_ids": sorted(partly),
        "missing_metabolite_ids": sorted(missing | no_msid),
    }


def tighten_potentials(ft, pkgmgr):
    """Make the potential bounds that potc already implies explicit.

    potc pins   potential = constant + (R*T)*logconc + dgerr
    with logconc and dgerr both box-bounded, so each potential variable has an
    exact implied range. The package leaves them at the +/-100000 default, which
    gives the branch-and-bound a hopeless LP relaxation. Writing the implied
    range onto the variable removes no feasible point -- it only tells the
    solver what the constraints already say.
    """
    from scipy.constants import R, kilo
    sp = pkgmgr.getpkg("SimpleThermoPkg")
    RT = R / kilo * ft.parameters["temperature"]
    n = 0
    for met_id, potc in ft.constraints.get("potc", {}).items():
        pot = sp.variables["potential"].get(met_id)
        if pot is None:
            continue
        const = potc.lb                      # constraint is equality: lb == ub
        lo = hi = const
        lc = ft.variables["logconc"].get(met_id)
        if lc is not None:
            lo += RT * lc.lb
            hi += RT * lc.ub
        de = ft.variables["dgerr"].get(met_id)
        if de is not None:
            lo += de.lb
            hi += de.ub
        pot.lb, pot.ub = lo, hi
        n += 1
    return n


def run(path, medium, mode, verbose=False, solver="highs"):
    """mode: 'fba' | 'asis' (deltag/4.184) | 'fixed' (deltag*4.184)"""
    model = load_model(path, medium, solver)
    base = model.optimize()
    result = {"model": Path(path).stem, "medium": medium, "mode": mode,
              "fba_status": base.status, "fba_growth": base.objective_value}
    if mode == "fba":
        result["fluxes"] = base.fluxes.to_dict()
        return result, model

    install_db_shim()
    if FAST:
        sys.path.insert(0, str(Path(__file__).parent))
        import fast_thermo
        fast_thermo.install()
    # the package computes `deltag / calorie`; kcal->kJ needs `* 4.184`
    fullthermopkg.calorie = KCAL_TO_KJ if mode == "asis" else 1.0 / KCAL_TO_KJ

    saved_objective = model.objective          # SimpleThermoPkg clobbers this
    t0 = time.time()
    pkgmgr = MSPackageManager.get_pkg_mgr(model)
    ft = pkgmgr.getpkg("FullThermoPkg")
    ft.build_package({"modelseed_path": MSDB}, verbose=verbose)
    result["build_seconds"] = round(time.time() - t0, 1)
    if TIGHTEN:
        result["potentials_tightened"] = tighten_potentials(ft, pkgmgr)
    model.objective = saved_objective

    sol = model.optimize()
    result["tfba_status"] = sol.status
    result["tfba_growth"] = sol.objective_value
    result["n_potc_constraints"] = len(ft.constraints.get("potc", {}))
    result["n_thermo_constraints"] = len(
        pkgmgr.getpkg("SimpleThermoPkg").constraints.get("thermo", {}))
    if sol.status == "optimal":
        result["fluxes"] = sol.fluxes.to_dict()
    return result, model


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default=str(_p.DEFAULT_MODEL))
    ap.add_argument("--medium", default="minimal", choices=["complete", "minimal"])
    ap.add_argument("--mode", default="asis", choices=["fba", "asis", "fixed"])
    ap.add_argument("--coverage", action="store_true")
    ap.add_argument("--solver", default="glpk")
    ap.add_argument("--tighten", action="store_true",
                    help="write potc-implied bounds onto the potential variables")
    ap.add_argument("--fast", action="store_true",
                    help="closed-form max_energy_magnitude instead of 2 MILPs per reaction")
    ap.add_argument("--out")
    a = ap.parse_args()

    FAST = a.fast
    TIGHTEN = a.tighten
    if a.coverage:
        m = load_model(a.model, a.medium, a.solver)
        cov = thermo_coverage(m, from_local(MSDB))
        print(json.dumps({k: v for k, v in cov.items() if not k.endswith("_ids")}, indent=2))
        if a.out:
            json.dump(cov, open(a.out, "w"), indent=2)
        sys.exit(0)

    res, model = run(a.model, a.medium, a.mode, solver=a.solver)
    flux = res.pop("fluxes", None)
    print(json.dumps(res, indent=2, default=str))
    if a.out:
        res["fluxes"] = flux
        json.dump(res, open(a.out, "w"), indent=2, default=str)
