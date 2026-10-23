#!/usr/bin/env python3
"""How many of a set of ModelSEED reactions are present in, and active in, the
FBA solutions of the core and genome-scale model panels.

Uses the repo's own cma model-set API (model sets, media, biomass policy and
the seed-key mapping from a model reaction id to its MSDB rxn id). Nothing from
ModelSEEDpy's fbapkg is involved.

"Active" is measured on the **pFBA** flux vector, matching the convention in
scripts/build_flux_distributions.py: plain FBA picks an arbitrary vertex among
alternate optima, so whether a given reaction carries flux in it is not
reproducible. The plain-FBA count is recorded alongside for comparison.

Usage::

    python3 presence_and_activity.py                       # all three panels
    python3 presence_and_activity.py --sets core_kegg2 --limit 50
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))          # scripts/, so `cma` imports
from cma import models as M  # noqa: E402
from cma import paths  # noqa: E402

#: The reaction list under test, shipped next to this script.
TARGETS_FILE = HERE / "target_reactions.tsv"
OUT = paths.results("direction_change_targets")

PANELS = ("core_kegg2", "ms2_gsm", "ms2_gsm_auxo")
ACTIVE_TOL = 1e-9        # build_flux_distributions.py
GROWTH_THRESHOLD = 1e-6  # repo convention

_SET = None
_TARGETS: set[str] = set()


def load_targets() -> dict[str, str]:
    out = {}
    with open(TARGETS_FILE) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["reaction_id"].strip()] = row["change"].strip()
    return out


def _init_worker(set_name: str, targets: list[str]) -> None:
    global _SET, _TARGETS
    _SET = M.get(set_name)
    _TARGETS = set(targets)
    os.environ.setdefault("OMP_NUM_THREADS", "1")


def analyse_one(model_id: str) -> dict:
    """Presence and pFBA activity of the target reactions in one model."""
    from cobra.flux_analysis import pfba

    rec = {"model_id": model_id, "error": "", "grows": False,
           "present": [], "active_pfba": [], "active_fba": []}
    try:
        ms = _SET
        model = M.load_model(ms, model_id)

        # model reaction id -> MSDB rxn id, via the set's seed-key policy
        seed_of = {}
        for r in model.reactions:
            sid = M.seed_id(r, ms.seed_key)
            if sid in _TARGETS:
                seed_of[r.id] = sid
        rec["present"] = sorted(set(seed_of.values()))
        if not seed_of:
            return rec                      # nothing to measure in this model

        M.apply_media(model, M.media(ms.media))
        bio = M.find_biomass(model, ms.biomass)
        if bio is None:
            rec["error"] = "no_biomass"
            return rec
        model.objective = bio

        sol = model.optimize()
        if sol.status != "optimal" or sol.objective_value is None:
            rec["error"] = f"status:{sol.status}"
            return rec
        if sol.objective_value <= GROWTH_THRESHOLD:
            return rec                      # non-grower has no flux solution
        rec["grows"] = True

        rec["active_fba"] = sorted({sid for rid, sid in seed_of.items()
                                    if abs(sol.fluxes[rid]) > ACTIVE_TOL})
        fluxes = pfba(model).fluxes
        rec["active_pfba"] = sorted({sid for rid, sid in seed_of.items()
                                     if abs(fluxes[rid]) > ACTIVE_TOL})
    except Exception as exc:                # noqa: BLE001
        rec["error"] = f"{type(exc).__name__}: {exc}"[:120]
    return rec


def run_panel(set_name: str, targets: dict[str, str], limit: int | None,
              jobs: int, keep: list | None = None) -> dict:
    ms = M.get(set_name)
    ids = ms.model_ids()
    if limit:
        ids = ids[:limit]
    t0 = time.time()
    recs = []
    with ProcessPoolExecutor(max_workers=jobs, initializer=_init_worker,
                             initargs=(set_name, list(targets))) as ex:
        for i, rec in enumerate(ex.map(analyse_one, ids, chunksize=4), 1):
            recs.append(rec)
            if i % 250 == 0:
                print(f"  [{set_name}] {i}/{len(ids)}  "
                      f"{time.time()-t0:.0f}s", flush=True)

    n_present = {r: 0 for r in targets}
    #: of the models containing the reaction, how many have non-zero biomass.
    #: A non-grower has no flux solution at all, so this is the denominator
    #: that `active` should really be read against.
    n_present_growing = {r: 0 for r in targets}
    n_active_pfba = {r: 0 for r in targets}
    n_active_fba = {r: 0 for r in targets}
    growers = 0
    errors = 0
    models_with_any = 0
    models_with_active = 0
    for rec in recs:
        if rec["error"]:
            errors += 1
        if rec["grows"]:
            growers += 1
        if rec["present"]:
            models_with_any += 1
        if rec["active_pfba"]:
            models_with_active += 1
        for r in rec["present"]:
            n_present[r] += 1
            if rec["grows"]:
                n_present_growing[r] += 1
        for r in rec["active_pfba"]:
            n_active_pfba[r] += 1
        for r in rec["active_fba"]:
            n_active_fba[r] += 1

    if keep is not None:
        # only models that contain at least one target are worth keeping
        keep.extend({"model_id": r["model_id"], "grows": r["grows"],
                     "present": r["present"], "active_pfba": r["active_pfba"],
                     "active_fba": r["active_fba"]}
                    for r in recs if r["present"])

    return {
        "set": set_name, "media": ms.media, "n_models": len(ids),
        "n_growers": growers, "n_errors": errors,
        "models_containing_any_target": models_with_any,
        "models_with_any_target_active": models_with_active,
        "targets_present_in_ge1_model": sum(1 for v in n_present.values() if v),
        "targets_active_in_ge1_model": sum(1 for v in n_active_pfba.values() if v),
        "targets_active_fba_ge1_model": sum(1 for v in n_active_fba.values() if v),
        "per_reaction": {r: {"present": n_present[r],
                             "present_growing": n_present_growing[r],
                             "active_pfba": n_active_pfba[r],
                             "active_fba": n_active_fba[r]} for r in targets},
        "seconds": round(time.time() - t0, 1),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sets", nargs="*", default=list(PANELS))
    ap.add_argument("--limit", type=int)
    ap.add_argument("-j", "--jobs", type=int, default=min(32, os.cpu_count() or 8))
    ap.add_argument("--tag", default="")
    ap.add_argument("--dump-records", action="store_true",
                    help="also write per_model_records.json (per-model detail)")
    a = ap.parse_args()

    targets = load_targets()
    print(f"{len(targets)} target reactions\n")
    OUT.mkdir(parents=True, exist_ok=True)

    out = {}
    records: dict[str, list] = {}
    for s in a.sets:
        print(f"=== {s}")
        keep = [] if a.dump_records else None
        out[s] = run_panel(s, targets, a.limit, a.jobs, keep)
        if keep is not None:
            records[s] = keep
        r = out[s]
        print(f"  models {r['n_models']}  growers {r['n_growers']}  errors {r['n_errors']}")
        print(f"  targets present in >=1 model : {r['targets_present_in_ge1_model']}/{len(targets)}")
        print(f"  targets active  in >=1 model : {r['targets_active_in_ge1_model']}/{len(targets)} (pFBA)")
        print(f"  ({r['seconds']}s)\n")

    path = OUT / f"presence_activity{a.tag}.json"
    json.dump({"targets": targets, "panels": out}, open(path, "w"), indent=2)
    print(f"wrote {path}")
    if records:
        rp = OUT / f"per_model_records{a.tag}.json"
        json.dump(records, open(rp, "w"))
        print(f"wrote {rp}")


if __name__ == "__main__":
    main()
