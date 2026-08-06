#!/usr/bin/env python3
"""Per-reaction, per-direction influence dataset for ALL core models (~5,683).

For every model, starting from its DEFAULT bounds, each internal reaction is set
ONE AT A TIME to each of the four ModelSEED direction options and the model is
re-solved, isolating that single reaction's influence:

    "<"  reverse     -> (-1000, 0)
    ">"  forward     -> (0, 1000)
    "="  reversible  -> (-1000, 1000)
    "?"  unknown/off -> (0, 0)          (knocked out)

Captured per (model, reaction, direction):
  * growth      -- max biomass flux (FBA objective value)
  * status      -- solver status
  * n_active    -- reactions carrying flux (|v|>tol) at the optimum (structural)
  * egc_atp/redox/mass -- on a CLOSED copy of the model with the reaction forced
      to this direction, the max rate a representative ATP / redox / mass drain
      can run (>0 => this direction lets the closed model make energy/mass from
      nothing = an energy-generating cycle becomes possible)
  * shadow_price duals -- full metabolite dual values (sparse: |dual|>tol), the
      marginal value/scarcity of each metabolite at the growth optimum

Also, per model: the full ``find_flux_loops`` (12-probe) EGC enumeration on the
baseline direction map (actual loop reactions), via growth_heuristics.flux_loops_one.

Output (partitioned Parquet under ``results/reaction_effects_all/`` — one shard
file per model so the run is parallel + restartable):
  effects/<model_id>.parquet        (model_id, rxn, base, default_dir, dir,
                                      growth, status, n_active, egc_atp,
                                      egc_redox, egc_mass)   ~2.8M rows total
  shadow_prices/<model_id>.parquet  (model_id, rxn, dir, metabolite,
                                      shadow_price)          ~180M rows total (~4 GB)
  model_flux_loops.jsonl            one JSON line per model: baseline EGC loops
  manifest.json                     run metadata + schema

Prereqs: cobra + modelseedpy (for the model-level find_flux_loops). Pure-cobra
for the per-reaction sweeps. Run:  python3 scripts/build_reaction_effects_all.py --workers 32
"""

from __future__ import annotations

import argparse
import datetime
import json
import sys
import time
import multiprocessing as mp
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPTS = Path(__file__).resolve().parent
ANALYSIS = SCRIPTS.parent
MODELS_DIR = ANALYSIS / "data" / "core_models_kegg2"
SITE_DATA = ANALYSIS / "site" / "data"
OUT_DIR = ANALYSIS / "results" / "reaction_effects_all"
IDS_SUMMARY = SITE_DATA / "all_models_rxnsets_summary.json"  # keys = all 5,683 ids

sys.path.insert(0, str(SCRIPTS))
from direction_change_template_eval import direction_from_bounds  # noqa: E402

OPTIONS = ["<", ">", "=", "?"]
OPTION_BOUNDS = {"<": (-1000.0, 0.0), ">": (0.0, 1000.0),
                 "=": (-1000.0, 1000.0), "?": (0.0, 0.0)}
TOL = 1e-6

# One representative drain per EGC probe group (ModelSEED cpd ids, product-side
# drain). The model-level find_flux_loops pass uses the full 12-probe battery;
# the per-(reaction,direction) sweep uses these three for speed.
PROBES = {
    "atp":   {"cpd00002": -1, "cpd00001": -1, "cpd00008": 1, "cpd00009": 1, "cpd00067": 1},
    "redox": {"cpd00004": -1, "cpd00003": 1, "cpd11749": 1},   # H2 absent -> usually 0
    "mass":  {"cpd00011": -1},                                  # CO2 sink
}

_CFG: dict = {}


def _init(baseline_map, params):
    import cobra
    import logging
    cobra.Configuration().solver = "glpk"
    logging.getLogger("cobra").setLevel(logging.ERROR)
    logging.getLogger("modelseedpy").setLevel(logging.ERROR)
    _CFG["baseline_map"] = baseline_map
    _CFG["params"] = params


def _find_biomass(model):
    for rid in ("bio1", "bio2", "biomass", "Biomass"):
        if rid in model.reactions:
            return model.reactions.get_by_id(rid)
    for r in model.reactions:
        if r.id.lower().startswith("bio") and not r.id.startswith("SK_"):
            return r
    return None


def _is_boundary(rid):
    return (rid.startswith(("EX_", "DM_", "SK_"))
            or rid.startswith("bio") or "biomass" in rid.lower())


def _internal_rxns(model):
    return [r for r in model.reactions if not _is_boundary(r.id)]


def _seed_base(r):
    from seed_annotation import seed_id
    return seed_id(r) or ""


def _build_closed_with_probes(model_id, comp="c0"):
    """Fresh model, fully closed (all EX_/DM_/SK_/bio at 0,0), plus the three
    representative drain probes wired to present metabolites."""
    import cobra
    m = cobra.io.load_json_model(str(MODELS_DIR / f"{model_id}.json"))
    for r in m.reactions:
        if _is_boundary(r.id):
            r.bounds = (0.0, 0.0)
    mets = {x.id for x in m.metabolites}
    probe_ids = {}
    for grp, stoich in PROBES.items():
        need = [f"{c}_{comp}" for c in stoich]
        if not all(mm in mets for mm in need):
            probe_ids[grp] = None
            continue
        rid = f"PROBE_{grp}"
        rxn = cobra.Reaction(rid)
        rxn.bounds = (0.0, 1000.0)
        rxn.add_metabolites({m.metabolites.get_by_id(f"{c}_{comp}"): coef
                             for c, coef in stoich.items()})
        m.add_reactions([rxn])
        probe_ids[grp] = rid
    return m, probe_ids


def eval_model(model_id):
    try:
        import cobra
        p = _CFG["params"]

        # ---- growth + shadow-price sweep (OPEN / default-bounds model) ----
        m = cobra.io.load_json_model(str(MODELS_DIR / f"{model_id}.json"))
        bio = _find_biomass(m)
        if bio is None:
            return {"model_id": model_id, "error": "no_biomass"}
        m.objective = bio
        m.objective_direction = "max"
        base_sol = m.optimize()
        base_flux = float(base_sol.objective_value) if base_sol.status == "optimal" else 0.0

        internal = _internal_rxns(m)
        eff_rows = []
        dual_rows = []
        for r in internal:
            base_dir = direction_from_bounds(r.lower_bound, r.upper_bound)
            seed = _seed_base(r)
            for opt in OPTIONS:
                lb, ub = OPTION_BOUNDS[opt]
                with m:
                    r.bounds = (lb, ub)
                    sol = m.optimize()
                    if sol.status == "optimal":
                        g = float(sol.objective_value or 0.0)
                        fl = sol.fluxes.values
                        n_active = int(np.count_nonzero(np.abs(fl) > TOL))
                        if p["duals"]:
                            sp = sol.shadow_prices
                            spv = sp.values
                            for j in np.where(np.abs(spv) > TOL)[0]:
                                dual_rows.append((model_id, r.id, opt,
                                                  sp.index[j], float(spv[j])))
                    else:
                        g = float("nan")
                        n_active = -1
                eff_rows.append([model_id, r.id, seed, base_dir, opt, g,
                                 sol.status, n_active])

        # ---- EGC-rate sweep (CLOSED model + 3 representative drains) ----
        egc = defaultdict(lambda: {"atp": 0.0, "redox": 0.0, "mass": 0.0})
        cm, probe_ids = _build_closed_with_probes(model_id)
        cm_internal = [r for r in cm.reactions if not _is_boundary(r.id)
                       and not r.id.startswith("PROBE_")]
        for grp, pid in probe_ids.items():
            if pid is None:
                continue
            cm.objective = pid
            cm.objective_direction = "max"
            for r in cm_internal:
                for opt in OPTIONS:
                    lb, ub = OPTION_BOUNDS[opt]
                    with cm:
                        r.bounds = (lb, ub)
                        v = cm.slim_optimize()
                    egc[(r.id, opt)][grp] = float(v) if v == v else 0.0  # nan->0

        # merge EGC rates onto the effect rows
        for row in eff_rows:
            rr = egc.get((row[1], row[4]), {})
            row.extend([round(rr.get("atp", 0.0), 4),
                        round(rr.get("redox", 0.0), 4),
                        round(rr.get("mass", 0.0), 4)])

        # ---- model-level full loop enumeration (baseline map, 12 probes) ----
        loops = None
        if p["model_loops"]:
            import growth_heuristics as gh
            loops = gh.flux_loops_one(model_id, _CFG["baseline_map"],
                                      objective="all", max_loops_per_probe=5)

        _write_shards(model_id, eff_rows, dual_rows if p["duals"] else [])
        return {"model_id": model_id, "base_flux": round(base_flux, 6),
                "n_rxn": len(internal), "n_dual_rows": len(dual_rows),
                "loops": loops, "error": None}
    except Exception as e:
        import traceback
        return {"model_id": model_id, "error": f"{type(e).__name__}: {e}",
                "trace": traceback.format_exc()[-400:]}


EFFECT_COLS = ["model_id", "rxn", "base", "default_dir", "dir", "growth",
               "status", "n_active", "egc_atp", "egc_redox", "egc_mass"]
DUAL_COLS = ["model_id", "rxn", "dir", "metabolite", "shadow_price"]


def _write_shards(model_id, eff_rows, dual_rows):
    (OUT_DIR / "effects").mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(eff_rows, columns=EFFECT_COLS)
    for c in ("growth", "egc_atp", "egc_redox", "egc_mass"):
        df[c] = df[c].astype("float32")
    df["n_active"] = df["n_active"].astype("int32")
    df.to_parquet(OUT_DIR / "effects" / f"{model_id}.parquet",
                  compression="zstd", index=False)
    if dual_rows:
        (OUT_DIR / "shadow_prices").mkdir(parents=True, exist_ok=True)
        dd = pd.DataFrame(dual_rows, columns=DUAL_COLS)
        dd["shadow_price"] = dd["shadow_price"].astype("float32")
        dd.to_parquet(OUT_DIR / "shadow_prices" / f"{model_id}.parquet",
                      compression="zstd", index=False)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--workers", type=int, default=max(1, min(mp.cpu_count() - 1, 32)))
    ap.add_argument("--limit", type=int, default=0, help="limit number of models (smoke runs)")
    ap.add_argument("--no-duals", action="store_true", help="skip shadow-price duals")
    ap.add_argument("--no-loops", action="store_true", help="skip model-level find_flux_loops")
    args = ap.parse_args(argv)

    ids = sorted(json.loads(IDS_SUMMARY.read_text()))
    if args.limit:
        ids = ids[: args.limit]
    baseline_map = json.loads((SITE_DATA / "baseline.json").read_text())["map"]
    params = {"duals": not args.no_duals, "model_loops": not args.no_loops}
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[reff] {len(ids)} models, workers={args.workers}, duals={params['duals']}, "
          f"model_loops={params['model_loops']}", flush=True)
    t0 = time.time()
    loops_path = OUT_DIR / "model_flux_loops.jsonl"
    errors = []
    n_ok = n_dual_rows = 0
    nw = max(1, min(args.workers, len(ids)))
    with open(loops_path, "w") as lf:
        if nw == 1:
            _init(baseline_map, params)
            it = (eval_model(m) for m in ids)
        else:
            ctx = mp.get_context("fork")
            pool = ctx.Pool(nw, initializer=_init, initargs=(baseline_map, params))
            it = pool.imap_unordered(eval_model, ids, chunksize=1)
        for i, rec in enumerate(it, 1):
            if rec.get("error"):
                errors.append((rec["model_id"], rec["error"]))
            else:
                n_ok += 1
                n_dual_rows += rec.get("n_dual_rows", 0)
                if rec.get("loops") is not None:
                    lf.write(json.dumps({"model_id": rec["model_id"],
                                         "base_flux": rec["base_flux"],
                                         "loops": rec["loops"]}) + "\n")
            if i % 100 == 0 or i == len(ids):
                print(f"  {i}/{len(ids)}  ok={n_ok} err={len(errors)}  "
                      f"dual_rows={n_dual_rows:,}  {i/max(time.time()-t0,1e-9):.1f} models/s",
                      flush=True)
        if nw != 1:
            pool.close(); pool.join()

    manifest = {
        "generated_at": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "n_models": n_ok,
        "n_errors": len(errors),
        "options": OPTIONS,
        "option_bounds": {k: list(v) for k, v in OPTION_BOUNDS.items()},
        "option_note": "'?' = unknown, tested as off/(0,0)",
        "egc_probes_per_reaction": {g: f"drain {list(s)}" for g, s in PROBES.items()},
        "egc_note": "egc_* = max drain rate on the CLOSED model with the reaction forced "
                    "to that direction (>0 => energy-generating cycle possible). redox "
                    "usually 0 (no H2/cpd11749 in core models). Model-level model_flux_loops.jsonl "
                    "has the full 12-probe find_flux_loops enumeration for the baseline map.",
        "effects_schema": EFFECT_COLS,
        "shadow_prices_schema": DUAL_COLS,
        "n_dual_rows": n_dual_rows,
        "duals": params["duals"],
        "model_loops": params["model_loops"],
    }
    (OUT_DIR / "manifest.json").write_text(json.dumps(manifest, indent=2))
    dt = time.time() - t0
    print(f"[reff] DONE {n_ok} ok, {len(errors)} errors in {dt/60:.1f} min; "
          f"{n_dual_rows:,} dual rows -> {OUT_DIR}", flush=True)
    if errors:
        print("  first errors:", errors[:5], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
