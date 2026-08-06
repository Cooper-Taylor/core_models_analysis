#!/usr/bin/env python3
"""Full FBA pipeline over all 5,683 Kegg2 core models, one run per
thermodynamic/direction source, plus the compound/reaction inventory +
coverage stats needed alongside it.

Sources with a genuine, global {rxn_id: operator} direction map
(``MAP_SOURCES``), built by ``build_thermo_source_direction_maps.py`` /
``build_modelseed_current_direction_map.py``:
  * ``modelseed``            -- canonical top-level DeltaG (unfiltered)
  * ``group_contribution``   -- Group Contribution's own DeltaG
  * ``equilibrator``         -- eQuilibrator (2.0)'s own DeltaG
  * ``dgpredictor``          -- dGPredictor (ModelSEED-retrained)'s own DeltaG
  * ``modelseed_current``    -- "what's currently in ModelSEED": eQuilibrator,
    falling back to the reaction's existing (Group-Contribution-backed)
    stored reversibility when eQuilibrator has no data (``estimate_one(rxn, 'EQ')``)

Plus one source with NO override at all (``implicit``): the reaction bounds
already baked into each Kegg2 core model file, i.e. whatever direction was
implicitly assigned when that model was built/gap-filled. This is
necessarily per-model rather than a global map, so it is handled as a
FBA-only special case (``FBA_SOURCES`` but not ``MAP_SOURCES``).

Per model (loaded once):
  * unique reactions/compounds = distinct base ModelSEED IDs
    (``seed_annotation.seed_id`` for reactions; compartment-suffix-stripped
    metabolite id for compounds) -- see reports/thermoComparison/THERMO_SOURCE_FBA_PIPELINE.md
    for why boundary/biomass pseudo-reactions are excluded.
  * for each MAP_SOURCES source, intersect the model's unique-reaction set
    against that source's coverage to get "# reactions with a defined
    direction under source X"; same for compounds against compound-energy
    coverage where applicable (dGPredictor and modelseed_current have no
    compound-level energies -- CPD_COVERAGE_SOURCES excludes them).
  * FBA is run once per FBA_SOURCES source: apply that source's direction map
    on top of the model's native bounds (``growth_heuristics.override_bounds``
    leaves reactions the source has no opinion on untouched; ``implicit``
    applies no override at all), restrict uptake to the standard KBase media
    (``apply_media``), optimize biomass.

Outputs (under ``results/thermo_source_fba_all_models/``):
  * ``model_results.csv`` -- one row per model.
  * ``summary_stats.json`` -- combined-across-all-models totals plus
    per-source FBA growth totals.
  * ``manifest.json`` -- run metadata.
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

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR", "/scratch/ctaylor/core_models_analysis"))
MODELS_DIR = ANALYSIS_DIR / "data" / "core_models_kegg2"
DATA_DIR = ANALYSIS_DIR / "results" / "thermo_source_fba_all_models"
SCRIPTS = ANALYSIS_DIR / "scripts"
sys.path.insert(0, str(SCRIPTS))

MAP_SOURCES = ("modelseed", "group_contribution", "equilibrator", "dgpredictor", "modelseed_current")
CPD_COVERAGE_SOURCES = ("modelseed", "group_contribution", "equilibrator", "dgpredictor")
FBA_SOURCES = MAP_SOURCES + ("implicit",)


def _load_direction_maps() -> dict:
    maps = {}
    for src in MAP_SOURCES:
        with open(DATA_DIR / f"rxn_directions_{src}.json") as fh:
            maps[src] = json.load(fh)
    return maps


def _load_coverage_sets() -> tuple[dict, dict]:
    """Returns (rxn_coverage_by_source -> set of rxn_id, cpd_coverage_by_source -> set of cpd_id)."""
    rxn_cov = {src: set() for src in MAP_SOURCES}
    with open(DATA_DIR / "rxn_source_coverage.csv") as fh:
        for row in csv.DictReader(fh):
            for src in ("modelseed", "group_contribution", "equilibrator", "dgpredictor"):
                if row[f"has_{src}"] == "True":
                    rxn_cov[src].add(row["rxn_id"])
    with open(DATA_DIR / "rxn_source_coverage_modelseed_current.csv") as fh:
        for row in csv.DictReader(fh):
            if row["has_modelseed_current"] == "True":
                rxn_cov["modelseed_current"].add(row["rxn_id"])

    cpd_cov = {src: set() for src in CPD_COVERAGE_SOURCES}
    with open(DATA_DIR / "cpd_source_coverage.csv") as fh:
        for row in csv.DictReader(fh):
            for src in CPD_COVERAGE_SOURCES:
                key = "has_modelseed" if src == "modelseed" else f"has_{src}"
                if row[key] == "True":
                    cpd_cov[src].add(row["cpd_id"])
    return rxn_cov, cpd_cov


_worker_state: dict = {}


def _init_worker(direction_maps, rxn_coverage, cpd_coverage):
    _worker_state["maps"] = direction_maps
    _worker_state["rxn_cov"] = rxn_coverage
    _worker_state["cpd_cov"] = cpd_coverage
    logging.getLogger("cobra").setLevel(logging.ERROR)


def eval_model(model_id: str) -> dict:
    import growth_heuristics as gh
    from cobra.io import load_json_model
    from seed_annotation import seed_id

    direction_maps = _worker_state["maps"]
    rxn_cov = _worker_state["rxn_cov"]
    cpd_cov = _worker_state["cpd_cov"]

    res = {"model_id": model_id, "error": ""}
    try:
        base_model = load_json_model(str(MODELS_DIR / f"{model_id}.json"))

        rxn_ids = set()
        for rxn in base_model.reactions:
            s = seed_id(rxn)
            if s:
                rxn_ids.add(s)
        cpd_ids = {met.id.rsplit("_", 1)[0] for met in base_model.metabolites}

        res["n_unique_reactions"] = len(rxn_ids)
        res["n_unique_compounds"] = len(cpd_ids)

        for src in MAP_SOURCES:
            res[f"n_reactions_with_direction_{src}"] = len(rxn_ids & rxn_cov[src])
        for src in CPD_COVERAGE_SOURCES:
            res[f"n_compounds_with_energy_{src}"] = len(cpd_ids & cpd_cov[src])

        for src in FBA_SOURCES:
            model = base_model.copy()
            if src == "implicit":
                stats = {"touched": 0, "unchanged": 0, "no_anno": 0}
            else:
                stats = gh.override_bounds(model, direction_maps[src])
            gh.apply_media(model)
            bio_rxn = gh.find_biomass_reaction(model)
            if bio_rxn is None:
                res[f"fba_status_{src}"] = "no_biomass"
                res[f"fba_grows_{src}"] = False
                res[f"fba_growth_flux_{src}"] = 0.0
                res[f"fba_n_overrides_{src}"] = stats["touched"]
                continue
            model.objective = bio_rxn
            sol = model.optimize()
            flux = float(sol.objective_value) if sol.objective_value is not None else 0.0
            res[f"fba_status_{src}"] = sol.status
            res[f"fba_growth_flux_{src}"] = flux
            res[f"fba_grows_{src}"] = (sol.status == "optimal") and (flux > gh.GROWTH_THRESHOLD)
            res[f"fba_n_overrides_{src}"] = stats["touched"]

        res["_rxn_ids"] = sorted(rxn_ids)
        res["_cpd_ids"] = sorted(cpd_ids)
    except Exception as e:
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
    print(f"{len(ids)} models to process")

    direction_maps = _load_direction_maps()
    rxn_cov, cpd_cov = _load_coverage_sets()
    for src in MAP_SOURCES:
        cpd_n = len(cpd_cov[src]) if src in cpd_cov else "n/a"
        print(f"  {src}: {len(direction_maps[src])} rxn directions, "
              f"{len(rxn_cov[src])} rxn coverage, {cpd_n} cpd coverage")
    print("  implicit: no direction map (uses each model's native on-disk bounds)")

    logging.getLogger("cobra").setLevel(logging.ERROR)
    ctx = mp.get_context("fork")
    rows = []
    errors = []
    global_rxn_union = set()
    global_cpd_union = set()
    global_rxn_union_by_src = {src: set() for src in MAP_SOURCES}
    global_cpd_union_by_src = {src: set() for src in CPD_COVERAGE_SOURCES}

    with ctx.Pool(args.workers, initializer=_init_worker,
                  initargs=(direction_maps, rxn_cov, cpd_cov)) as pool:
        for i, rec in enumerate(pool.imap_unordered(eval_model, ids, chunksize=4), 1):
            if rec.get("error"):
                errors.append(rec)
                print(f"  [{i}/{len(ids)}] ERROR {rec['model_id']}: {rec['error']}")
            else:
                rxn_ids = set(rec.pop("_rxn_ids"))
                cpd_ids = set(rec.pop("_cpd_ids"))
                global_rxn_union |= rxn_ids
                global_cpd_union |= cpd_ids
                for src in MAP_SOURCES:
                    global_rxn_union_by_src[src] |= (rxn_ids & rxn_cov[src])
                for src in CPD_COVERAGE_SOURCES:
                    global_cpd_union_by_src[src] |= (cpd_ids & cpd_cov[src])
                rows.append(rec)
            if i % 500 == 0 or i == len(ids):
                print(f"  [{i}/{len(ids)}] processed")

    rows.sort(key=lambda r: r["model_id"])
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = DATA_DIR / "model_results.csv"
    if rows:
        fieldnames = list(rows[0].keys())
        with open(csv_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(rows)
    print(f"wrote {csv_path} ({len(rows)} rows)")

    summary = {
        "n_models": len(ids),
        "n_ok": len(rows),
        "n_errors": len(errors),
        "errors": [{"model_id": e["model_id"], "error": e["error"]} for e in errors],
        "combined_unique_reactions_all_models": len(global_rxn_union),
        "combined_unique_compounds_all_models": len(global_cpd_union),
        "combined_reactions_with_direction_by_source": {
            src: len(global_rxn_union_by_src[src]) for src in MAP_SOURCES
        },
        "combined_compounds_with_energy_by_source": {
            src: len(global_cpd_union_by_src[src]) for src in CPD_COVERAGE_SOURCES
        },
        "fba_growth_totals_by_source": {
            src: sum(1 for r in rows if r.get(f"fba_grows_{src}")) for src in FBA_SOURCES
        },
    }
    with open(DATA_DIR / "summary_stats.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))

    manifest = {
        "n_models": len(ids),
        "n_errors": len(errors),
        "workers": args.workers,
        "elapsed_s": round(time.time() - t0, 2),
        "map_sources": list(MAP_SOURCES),
        "fba_sources": list(FBA_SOURCES),
    }
    with open(DATA_DIR / "manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
