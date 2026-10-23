#!/usr/bin/env python3
"""Flux effect of the thermodynamic direction locks, swept over many core models.

This is the TMFA constraint reduced to what it actually does to a reaction that
cannot escape its sign: fix the direction. It ignores the cross-reaction
coupling that the full MILP captures (metabolite potentials are shared), so it
is a LOWER bound on the effect -- but it is fast enough to run over the panel
and it isolates the contribution of the delta-G data itself.

Three direction sources are compared:
  compound-sum  : sum_i s_i * dGf_i from the compound `deltag` scalars, which is
                  the only thing FullThermoPkg can build
  reaction-level: the DB's own reaction `deltag`
  eQuilibrator  : the DB's per-source eQuilibrator reaction value
"""
import glob, json, math, sys
from pathlib import Path

import cobra
from cobra.io import load_json_model

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
import _paths as _p

BIOCHEM = str(_p.BIOCHEM)
MEDIA = _p.MEDIA
MODELS = _p.MODELS
SENTINEL = 10000000
KCAL_TO_KJ = 4.184
RT = 8.314462618 / 1000.0 * 298.0
CONC_SWING = RT * (math.log(0.02) - math.log(1e-6))
ERR = 5.0

MINIMAL = {"cpd00027": -10.0, "cpd00007": -20.0}
OPEN = {"cpd00001", "cpd00067", "cpd00009", "cpd00013", "cpd00048", "cpd00205",
        "cpd00971", "cpd00099", "cpd00030", "cpd00034", "cpd00058", "cpd00149",
        "cpd10515", "cpd10516", "cpd00254", "cpd00063", "cpd00011"}

cpd_dg, rxn_dg, rxn_eq = {}, {}, {}
for f in sorted(glob.glob(f"{BIOCHEM}/compound_*.json")):
    for c in json.load(open(f)):
        v = c.get("deltag")
        cpd_dg[c["id"]] = None if v in (None, SENTINEL) else float(v)
for f in sorted(glob.glob(f"{BIOCHEM}/reaction_*.json")):
    for r in json.load(open(f)):
        try:
            v = float(r.get("deltag"))
        except (TypeError, ValueError):
            v = None
        rxn_dg[r["id"]] = None if v is None or v == SENTINEL else v * KCAL_TO_KJ
        eq = (r.get("thermodynamics") or {}).get("eQuilibrator")
        rxn_eq[r["id"]] = eq[0] * KCAL_TO_KJ if eq and eq[0] not in (None, SENTINEL) else None


def apply_medium(model):
    for rxn in model.reactions:
        if not rxn.id.startswith("EX_"):
            continue
        mets = list(rxn.metabolites)
        if len(mets) != 1:
            continue
        cpd = mets[0].id.split("_")[0]
        rxn.upper_bound = 1000.0
        rxn.lower_bound = MINIMAL.get(cpd, -1000.0 if cpd in OPEN else 0.0)


def locks(model, source, scale=KCAL_TO_KJ):
    """Return {reaction_id: 'fwd'|'rev'} for reactions the thermo data pins."""
    out = {}
    for rxn in model.reactions:
        if rxn.id[:3] in ("EX_", "SK_", "DM_") or rxn.id.startswith("bio"):
            continue
        rid = rxn.id.split("_")[0]
        if source == "compound-sum":
            dg, swing, ok = 0.0, 0.0, True
            for met, s in rxn.metabolites.items():
                v = cpd_dg.get(met.id.split("_")[0])
                if v is None:
                    ok = False
                    break
                dg += s * v
                if met.id.split("_")[0] not in ("cpd00001", "cpd00067"):
                    swing += abs(s) * CONC_SWING
                swing += abs(s) * ERR
            if not ok:
                continue
            dg *= scale
        else:
            dg = (rxn_dg if source == "reaction-level" else rxn_eq).get(rid)
            if dg is None:
                continue
            swing = sum(abs(s) * (CONC_SWING + ERR) for s in rxn.metabolites.values())
        if abs(dg) > swing:
            out[rxn.id] = "fwd" if dg < 0 else "rev"
    return out


def growth_with(path, source=None, scale=KCAL_TO_KJ):
    model = load_json_model(path)
    model.solver = "glpk"
    apply_medium(model)
    n = 0
    if source:
        for rid, d in locks(model, source, scale).items():
            r = model.reactions.get_by_id(rid)
            if d == "fwd":
                r.lower_bound = max(r.lower_bound, 0.0)
            else:
                r.upper_bound = min(r.upper_bound, 0.0)
            n += 1
    s = model.optimize()
    return (s.objective_value if s.status == "optimal" else 0.0), s.status, n


if __name__ == "__main__":
    n_models = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    files = sorted(MODELS.glob("GCF_*.json"))[:n_models]
    rows = []
    for mf in files:
        try:
            base, st, _ = growth_with(mf)
            rec = {"model": mf.stem, "base": base, "status": st}
            for src, scale in (("compound-sum", KCAL_TO_KJ),
                               ("compound-sum-asis", 1.0 / KCAL_TO_KJ),
                               ("reaction-level", KCAL_TO_KJ),
                               ("eQuilibrator", KCAL_TO_KJ)):
                s = "compound-sum" if src.startswith("compound-sum") else src
                g, stt, nl = growth_with(mf, s, scale)
                rec[src] = g
                rec[src + "_nlocks"] = nl
            rows.append(rec)
        except Exception as e:
            rows.append({"model": mf.stem, "error": str(e)[:80]})
    good = [r for r in rows if "error" not in r and r["base"] > 1e-6]
    print(f"models run: {len(rows)}  with baseline growth: {len(good)}\n")
    print(f"{'direction source':22s} {'mean locks':>10s} {'lost growth':>12s} {'no growth':>10s} {'mean % of base':>15s}")
    for src in ("compound-sum-asis", "compound-sum", "reaction-level", "eQuilibrator"):
        nl = sum(r[src + "_nlocks"] for r in good) / len(good)
        dead = sum(1 for r in good if r[src] < 1e-6)
        red = sum(1 for r in good if r[src] < r["base"] - 1e-6)
        pct = sum(100 * r[src] / r["base"] for r in good) / len(good)
        print(f"{src:22s} {nl:10.1f} {red:12d} {dead:10d} {pct:14.1f}%")
    json.dump(rows, open(str(_p.RESULTS / "lock_sweep.json"), "w"), indent=2)
