#!/usr/bin/env python3
"""Regenerate ``results/results.csv`` with every model's reaction bounds
re-bound from a ModelSEED reversibility map, instead of the on-disk
(template-time) bounds that ``analyze_growth.py`` reads verbatim.

Motivation: the reversibility cascade was changed (Heuristics-Review fixes
H2 + H3 adopted as the new baseline in
``ModelSEEDDatabase/Scripts/Thermodynamics/Estimate_Reaction_Reversibility.py``).
To make the descriptive 100-model panel reflect the adopted cascade, the
grower set must be recomputed under it before re-running ``select_diverse.py``.

The reversibility map is read straight from the (freshly regenerated) MSDB
reaction JSONs' canonical ``reversibility`` field, so this growth panel is
consistent with the committed database directions. core_models_kegg2/*.json
are never written -- the rebind is in-memory only.

Schema matches ``analyze_growth.py`` exactly so all downstream consumers
(select_diverse.py, site build) keep working unchanged.
"""
from __future__ import annotations

import csv
import logging
import multiprocessing as mp
import os
import sys
from pathlib import Path

import cobra
from cobra.io import load_json_model

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR", "/scratch/ctaylor/core_models_analysis"))
MSDB_ROOT = os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase")
MODELS_DIR = ANALYSIS_DIR / "data" / "core_models_kegg2"
MEDIA_FILE = Path(MSDB_ROOT + "/Media/KBaseMedia.cpd")
RESULTS_CSV = ANALYSIS_DIR / "results" / "results.csv"

sys.path.insert(0, str(ANALYSIS_DIR / "scripts"))
from seed_annotation import seed_id  # noqa: E402  (normalizes _c-suffix)

GROWTH_THRESHOLD = 1e-6
FIELDS = [
    "model_id", "n_metabolites", "n_reactions", "n_genes",
    "biomass_rxn", "n_exchanges_total", "n_exchanges_open",
    "status", "growth_flux", "grows", "error",
]


def load_media_compounds(path: Path) -> set:
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}


def load_reversibility_map() -> dict:
    """Canonical ModelSEED reversibility, straight from the reaction JSONs."""
    sys.path.insert(0, MSDB_ROOT + "/Libs/Python")
    from BiochemPy import Reactions
    rxns = Reactions().loadReactions()
    return {rid: rxns[rid].get("reversibility") for rid in rxns}


def _bounds_for_rev(rev: str, default_bound: float = 1000.0):
    if rev == ">":
        return 0.0, default_bound
    if rev == "<":
        return -default_bound, 0.0
    return -default_bound, default_bound   # '=' and '?' -> reversible


# --- worker state (pickled once per child) ---------------------------------
_STATE: dict = {}


def _init_worker(rev_map, media_cpds):
    _STATE["rev_map"] = rev_map
    _STATE["media"] = media_cpds
    logging.getLogger("cobra").setLevel(logging.ERROR)


def _apply_media(model, media_cpds) -> int:
    open_count = 0
    for rxn in model.reactions:
        if not rxn.id.startswith("EX_"):
            continue
        mets = list(rxn.metabolites.keys())
        if len(mets) != 1:
            continue
        cpd_id = mets[0].id.split("_")[0]
        if cpd_id in media_cpds:
            rxn.lower_bound = -1000.0
            open_count += 1
        else:
            rxn.lower_bound = 0.0
        if rxn.upper_bound < 1000.0:
            rxn.upper_bound = 1000.0
    return open_count


def _find_biomass(model):
    if "bio1" in model.reactions:
        return model.reactions.get_by_id("bio1")
    for rid in ("bio2", "biomass", "Biomass"):
        if rid in model.reactions:
            return model.reactions.get_by_id(rid)
    for r in model.reactions:
        if r.id.lower().startswith("bio") and not r.id.startswith("SK_"):
            return r
    return None


def _rebind(model, rev_map) -> None:
    """Full rebind: every non-exchange reaction whose seed id is in the map
    gets bounds from its cascade direction. Exchanges/sinks/biomass untouched
    so the media step stays valid."""
    for rxn in model.reactions:
        seed = seed_id(rxn)
        if not seed:
            continue
        rev = rev_map.get(seed)
        if rev is None:
            continue
        rxn.lower_bound, rxn.upper_bound = _bounds_for_rev(rev)


def analyze_one(path_str: str) -> dict:
    path = Path(path_str)
    res = {k: 0 for k in FIELDS}
    res.update({"model_id": path.stem, "biomass_rxn": "", "status": "",
                "growth_flux": 0.0, "grows": False, "error": ""})
    try:
        model = load_json_model(str(path))
        res["n_metabolites"] = len(model.metabolites)
        res["n_reactions"] = len(model.reactions)
        res["n_genes"] = len(model.genes)
        res["n_exchanges_total"] = sum(1 for r in model.reactions if r.id.startswith("EX_"))

        _rebind(model, _STATE["rev_map"])
        res["n_exchanges_open"] = _apply_media(model, _STATE["media"])

        bio = _find_biomass(model)
        if bio is None:
            res["status"] = "no_biomass"
            return res
        res["biomass_rxn"] = bio.id
        model.objective = bio
        sol = model.optimize()
        res["status"] = sol.status
        flux = float(sol.objective_value) if sol.objective_value is not None else 0.0
        res["growth_flux"] = flux
        res["grows"] = (sol.status == "optimal") and (flux > GROWTH_THRESHOLD)
    except Exception as e:  # noqa: BLE001
        res["error"] = f"{type(e).__name__}: {e}"
    return res


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=max(1, min(mp.cpu_count() - 1, 64)))
    ap.add_argument("--out", default=str(RESULTS_CSV))
    args = ap.parse_args()

    print("loading reversibility map from MSDB reaction JSONs ...", flush=True)
    rev_map = load_reversibility_map()
    print(f"  {len(rev_map)} reactions in map", flush=True)
    media = load_media_compounds(MEDIA_FILE)
    print(f"  media has {len(media)} compounds", flush=True)

    model_paths = sorted(str(p) for p in MODELS_DIR.glob("*.json"))
    print(f"found {len(model_paths)} models; running with {args.workers} workers", flush=True)

    rows = []
    grew = 0
    ctx = mp.get_context("fork")
    with ctx.Pool(args.workers, initializer=_init_worker, initargs=(rev_map, media)) as pool:
        for i, res in enumerate(pool.imap_unordered(analyze_one, model_paths, chunksize=4), 1):
            rows.append(res)
            grew += bool(res["grows"])
            if i % 500 == 0:
                print(f"  {i}/{len(model_paths)}  growing={grew}", flush=True)

    rows.sort(key=lambda r: r["model_id"])
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)

    print(f"\nDone. {len(rows)} models; {grew} grow (>{GROWTH_THRESHOLD}).", flush=True)
    print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
