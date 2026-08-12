#!/usr/bin/env python3
"""FBA over all Kegg2 core models under each thermodynamic source in isolation
and under the graded/recommended source, using the maps from
``build_graded_direction_maps.py``.

Variants (all run on the same model, same media, same cascade):

    implicit        no override at all -- the bounds baked into the model file
    gc              Group Contribution only
    eq              eQuilibrator only
    dgpms           dGPredictor-ModelSEED only
    graded          per reaction, the best-graded source (any grade)
    graded_trusted  same, but BRONZE-best reactions get no call at all

Mirrors ``run_thermo_source_fba_all_models.py`` -- same per-model loading, same
``growth_heuristics`` overlay/media/biomass path, same "reactions the source has
no opinion on keep their native bounds" policy -- so the numbers are directly
comparable to the 2026-08-03 sweep. Two things are added: per-variant counts of
bounds actually CHANGED (as opposed to merely touched), and the per-reaction
operator table needed to attribute growth differences.

Outputs under ``results/thermo_grades_fba/``:
    model_results.csv    one row per model
    summary_stats.json   combined totals + growth counts per variant
    manifest.json        run metadata
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
MODELS_DIR = ANALYSIS_DIR / "data" / "core_models_kegg2"
DATA_DIR = ANALYSIS_DIR / "results" / "thermo_grades_fba"
SCRIPTS = ANALYSIS_DIR / "scripts"
sys.path.insert(0, str(SCRIPTS))

MAP_SOURCES = ("gc", "eq", "dgpms", "graded", "graded_trusted", "graded_heldout")
FBA_SOURCES = ("implicit",) + MAP_SOURCES


def _load_direction_maps() -> dict:
    return {src: json.load(open(DATA_DIR / f"rxn_directions_{src}.json"))
            for src in MAP_SOURCES}


def _load_coverage_sets() -> dict:
    cov = {src: set() for src in MAP_SOURCES}
    with open(DATA_DIR / "rxn_source_coverage.csv") as fh:
        for row in csv.DictReader(fh):
            for src in MAP_SOURCES:
                if row[f"has_{src}"] == "True":
                    cov[src].add(row["rxn_id"])
    return cov


_state: dict = {}


def _init_worker(maps, cov):
    _state["maps"] = maps
    _state["cov"] = cov
    logging.getLogger("cobra").setLevel(logging.ERROR)


def eval_model(model_id: str) -> dict:
    import growth_heuristics as gh
    from cobra.io import load_json_model
    from seed_annotation import seed_id

    maps, cov = _state["maps"], _state["cov"]
    res = {"model_id": model_id, "error": ""}
    try:
        base = load_json_model(str(MODELS_DIR / f"{model_id}.json"))
        rxn_ids = {s for s in (seed_id(r) for r in base.reactions) if s}
        cpd_ids = {m.id.rsplit("_", 1)[0] for m in base.metabolites}
        res["n_unique_reactions"] = len(rxn_ids)
        res["n_unique_compounds"] = len(cpd_ids)
        for src in MAP_SOURCES:
            res[f"n_reactions_with_direction_{src}"] = len(rxn_ids & cov[src])

        # native bounds, keyed by SEED id, to count real changes per variant
        native = {}
        for r in base.reactions:
            s = seed_id(r)
            if s:
                native[s] = (r.lower_bound, r.upper_bound)

        for src in FBA_SOURCES:
            model = base.copy()
            if src == "implicit":
                stats = {"touched": 0}
            else:
                stats = gh.override_bounds(model, maps[src])
            changed = 0
            for r in model.reactions:
                s = seed_id(r)
                if s and s in native and (r.lower_bound, r.upper_bound) != native[s]:
                    changed += 1
            gh.apply_media(model)
            bio = gh.find_biomass_reaction(model)
            res[f"fba_n_overrides_{src}"] = stats["touched"]
            res[f"fba_n_bounds_changed_{src}"] = changed
            if bio is None:
                res[f"fba_status_{src}"] = "no_biomass"
                res[f"fba_grows_{src}"] = False
                res[f"fba_growth_flux_{src}"] = 0.0
                continue
            model.objective = bio
            sol = model.optimize()
            flux = float(sol.objective_value) if sol.objective_value is not None else 0.0
            res[f"fba_status_{src}"] = sol.status
            res[f"fba_growth_flux_{src}"] = flux
            res[f"fba_grows_{src}"] = (sol.status == "optimal") and (flux > gh.GROWTH_THRESHOLD)

        res["_rxn_ids"] = sorted(rxn_ids)
        res["_cpd_ids"] = sorted(cpd_ids)
    except Exception as e:  # noqa: BLE001 -- one bad model must not kill the sweep
        res["error"] = f"{type(e).__name__}: {e}"
    return res


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--workers", type=int, default=max(1, min(mp.cpu_count() - 1, 32)))
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args(argv)

    t0 = time.time()
    ids = sorted(p.stem for p in MODELS_DIR.glob("*.json"))
    if args.limit:
        ids = ids[: args.limit]
    maps = _load_direction_maps()
    cov = _load_coverage_sets()
    print(f"{len(ids)} models")
    for src in MAP_SOURCES:
        print(f"  {src:15s} {len(maps[src]):6d} directions")
    print("  implicit        native on-disk bounds, no override")

    logging.getLogger("cobra").setLevel(logging.ERROR)
    rows, errors = [], []
    union_rxn, union_cpd = set(), set()
    union_by_src = {src: set() for src in MAP_SOURCES}
    ctx = mp.get_context("fork")
    with ctx.Pool(args.workers, initializer=_init_worker, initargs=(maps, cov)) as pool:
        for i, rec in enumerate(pool.imap_unordered(eval_model, ids, chunksize=4), 1):
            if rec.get("error"):
                errors.append(rec)
                print(f"  ERROR {rec['model_id']}: {rec['error']}")
            else:
                r_ids = set(rec.pop("_rxn_ids"))
                c_ids = set(rec.pop("_cpd_ids"))
                union_rxn |= r_ids
                union_cpd |= c_ids
                for src in MAP_SOURCES:
                    union_by_src[src] |= (r_ids & cov[src])
                rows.append(rec)
            if i % 1000 == 0 or i == len(ids):
                print(f"  [{i}/{len(ids)}]")

    rows.sort(key=lambda r: r["model_id"])
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with open(DATA_DIR / "model_results.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote model_results.csv ({len(rows)} rows)")

    summary = {
        "n_models": len(ids), "n_ok": len(rows), "n_errors": len(errors),
        "errors": [{"model_id": e["model_id"], "error": e["error"]} for e in errors][:20],
        "combined_unique_reactions_all_models": len(union_rxn),
        "combined_unique_compounds_all_models": len(union_cpd),
        "combined_reactions_with_direction_by_source": {
            src: len(union_by_src[src]) for src in MAP_SOURCES},
        "fba_growth_totals_by_source": {
            src: sum(1 for r in rows if r.get(f"fba_grows_{src}")) for src in FBA_SOURCES},
        "median_bounds_changed_by_source": {
            src: sorted(r[f"fba_n_bounds_changed_{src}"] for r in rows)[len(rows) // 2]
            for src in FBA_SOURCES},
    }
    json.dump(summary, open(DATA_DIR / "summary_stats.json", "w"), indent=2)
    print(json.dumps(summary, indent=2))
    json.dump({"n_models": len(ids), "n_errors": len(errors), "workers": args.workers,
               "elapsed_s": round(time.time() - t0, 2),
               "fba_sources": list(FBA_SOURCES)},
              open(DATA_DIR / "manifest.json", "w"), indent=2)


if __name__ == "__main__":
    main()
