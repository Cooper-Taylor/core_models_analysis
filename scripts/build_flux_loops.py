#!/usr/bin/env python3
"""Energy-generating-cycle (EGC) detection per reversibility variant, per panel model.

For the baseline direction map and each reversibility variant, re-bounds every
panel model to that map, CLOSES the model, and runs ``find_flux_loops`` (via
``growth_heuristics.flux_loops_one``) for the full probe battery (atp+redox+mass,
12 probes).  Reports the energy-generating cycles each heuristic induces vs. the
baseline -- a good reversibility heuristic should not manufacture free-ATP/redox
loops.

Auto-scoped: a variant is evaluated only on models that contain at least one
reaction the variant changes (a variant that changes no in-model reaction cannot
change that model's EGCs -> delta 0, omitted).

Output (under ``site/data/``):
  - ``panel_flux_loops.json``
      {"meta": {...}, "models": {mid: {"baseline": {...},
                                       "variants": {tag: {...}}}},
       "global": {"by_variant": {tag: {...}}}}

Prerequisite: ``modelseedpy`` must be installed (find_flux_loops needs
MSModelUtil + its ReactionUse MILP package).  Inputs: site/data/baseline.json,
site/data/variants/*.json, site/data/panel_rxnsets.json, results/selected_ids.txt.
Runs the closed-model MILP battery per (model, variant); parallelize with --workers.
"""

from __future__ import annotations

import argparse
import datetime
import json
import logging
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path
from typing import Optional

SCRIPTS = Path(__file__).resolve().parent
ANALYSIS_ROOT = SCRIPTS.parent
SITE_DATA = ANALYSIS_ROOT / "site" / "data"
MSDB_ROOT = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase"))
PANEL_FILE = ANALYSIS_ROOT / "results" / "selected_ids.txt"
OUT_FILE = SITE_DATA / "panel_flux_loops.json"

sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(MSDB_ROOT / "Libs" / "Python"))

import growth_heuristics as gh

_CFG: dict = {}


def _init(baseline_map, variant_diffs, params):
    _CFG["baseline_map"] = baseline_map
    _CFG["variant_diffs"] = variant_diffs  # {tag: [{rxn,new}, ...]}
    _CFG["params"] = params
    logging.getLogger("cobra").setLevel(logging.ERROR)
    logging.getLogger("modelseedpy").setLevel(logging.ERROR)


def _work(task):
    model_id, tag = task
    eff_map = dict(_CFG["baseline_map"])
    if tag != "baseline":
        for d in _CFG["variant_diffs"][tag]:
            eff_map[d["rxn"]] = d["new"]
    p = _CFG["params"]
    res = gh.flux_loops_one(
        model_id, eff_map,
        objective=p["objective"],
        max_loops_per_probe=p["max_loops"],
    )
    return (model_id, tag, res)


def _probe_signature(loop_list):
    """Order-independent signature of a probe's loop set: {frozenset((id,dir))}."""
    return frozenset(
        frozenset((r["id"], r["dir"]) for r in loop["reactions"])
        for loop in loop_list
    )


def main(argv: Optional[list] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--workers", type=int, default=max(1, min(mp.cpu_count() - 1, 16)))
    ap.add_argument("--limit", type=int, default=None, help="limit number of panel models (smoke runs)")
    ap.add_argument("--variants", type=str, default=None, help="comma-separated variant subset")
    ap.add_argument("--max-loops", type=int, default=5)
    ap.add_argument("--objective", type=str, default="all", choices=["all", "atp", "redox", "mass"])
    args = ap.parse_args(argv)

    baseline_map = json.loads((SITE_DATA / "baseline.json").read_text())["map"]
    panel_ids = PANEL_FILE.read_text().split()
    if args.limit:
        panel_ids = panel_ids[: args.limit]
    panel_set = set(panel_ids)
    panel_rxnsets = json.loads((SITE_DATA / "panel_rxnsets.json").read_text())

    # Load variant diffs
    variant_diffs = {}
    for vfile in sorted((SITE_DATA / "variants").glob("*.json")):
        tag = vfile.stem
        if args.variants and tag not in args.variants.split(","):
            continue
        variant_diffs[tag] = json.loads(vfile.read_text())["diffs"]
    tags = sorted(variant_diffs)

    # Build task list: baseline for all panel models; each variant only on models
    # containing >=1 of its changed reactions (auto-scope).  Track n_changed_here.
    tasks = [(mid, "baseline") for mid in panel_ids]
    n_changed_here = {}  # (mid, tag) -> int
    for tag in tags:
        changed = {d["rxn"] for d in variant_diffs[tag]}
        for mid in panel_ids:
            present = changed & set(panel_rxnsets.get(mid, []))
            if present:
                tasks.append((mid, tag))
                n_changed_here[(mid, tag)] = len(present)

    params = {"objective": args.objective, "max_loops": args.max_loops}
    print(f"[egc] {len(panel_ids)} models x {len(tags)} variants -> {len(tasks)} "
          f"(model,variant) tasks; objective={args.objective} max_loops={args.max_loops}", flush=True)

    t0 = time.time()
    nw = max(1, min(args.workers, len(tasks)))
    if nw == 1:
        _init(baseline_map, variant_diffs, params)
        results = [_work(t) for t in tasks]
    else:
        ctx = mp.get_context("fork")
        with ctx.Pool(nw, initializer=_init, initargs=(baseline_map, variant_diffs, params)) as pool:
            results = list(pool.imap_unordered(_work, tasks, chunksize=1))
    dt = time.time() - t0

    # Organize: models[mid] = {baseline, variants:{tag}}
    per_model = {mid: {"baseline": None, "variants": {}} for mid in panel_ids}
    n_err = 0
    for mid, tag, res in results:
        if res.get("status") != "ok":
            n_err += 1
        if tag == "baseline":
            per_model[mid]["baseline"] = res
        else:
            per_model[mid]["variants"][tag] = res

    # Assemble output: variant deltas vs baseline; store variant loops only where
    # the per-probe loop set differs from baseline (keeps the file small).
    models_out = {}
    glob = {t: {"n_models_scoped": 0, "n_models_new_egc": 0,
                "new_atp": 0, "new_redox": 0, "new_mass": 0, "net_total_delta": 0}
            for t in tags}
    for mid in panel_ids:
        base = per_model[mid]["baseline"]
        if base is None:
            continue
        base_block = {"n_loops_total": base["n_loops_total"], "by_group": base["by_group"],
                      "by_probe": base["by_probe"], "loops": base["loops"]}
        base_sig = {p: _probe_signature(ll) for p, ll in base["loops"].items()}
        variants_out = {}
        for tag, r in per_model[mid]["variants"].items():
            delta = {g: r["by_group"][g] - base["by_group"][g] for g in ("atp", "redox", "mass")}
            delta["total"] = r["n_loops_total"] - base["n_loops_total"]
            all_probes = set(r["by_probe"]) | set(base["by_probe"])
            new_probes = sorted(p for p in all_probes
                                if r["by_probe"].get(p, 0) > base["by_probe"].get(p, 0))
            resolved_probes = sorted(p for p in all_probes
                                     if base["by_probe"].get(p, 0) > 0 and r["by_probe"].get(p, 0) == 0)
            # Keep only probes whose loop set differs from baseline.
            diff_loops = {p: ll for p, ll in r["loops"].items()
                          if _probe_signature(ll) != base_sig.get(p, frozenset())}
            variants_out[tag] = {
                "in_scope": True,
                "n_changed_here": n_changed_here.get((mid, tag), 0),
                "n_loops_total": r["n_loops_total"],
                "by_group": r["by_group"],
                "by_probe": r["by_probe"],
                "delta_vs_baseline": delta,
                "new_probes": new_probes,
                "resolved_probes": resolved_probes,
                "loops": diff_loops,
            }
            gv = glob[tag]
            gv["n_models_scoped"] += 1
            if delta["total"] > 0:
                gv["n_models_new_egc"] += 1
            for g in ("atp", "redox", "mass"):
                gv[f"new_{g}"] += max(0, delta[g])
            gv["net_total_delta"] += delta["total"]
        models_out[mid] = {"baseline": base_block, "variants": variants_out}

    from kbutillib.ms_fba_utils import EGC_PROBE_CATALOG
    probe_groups = {g: [p["name"] for p in EGC_PROBE_CATALOG[g]] for g in ("atp", "redox", "mass")}
    meta = {
        "objective": args.objective,
        "max_loops_per_probe": args.max_loops,
        "fraction_of_optimum": 0.999,
        "probe_groups": probe_groups,
        "n_models": len(models_out),
        "n_variants": len(tags),
        "generated_at": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "redox_note": "H2 (cpd11749) is absent in core panel models, so redox probes "
                      "structurally cannot fire (report 0) -- expected, not a bug.",
        "closed_model_note": "All EX_/DM_/SK_ exchanges + real biomass are closed; the "
                             "small ATP-maintenance 'bio2' reaction is left open (it is "
                             "reused as the ATP probe). +delta loops = heuristic manufactured "
                             "new energy-generating cycles (bad).",
    }
    OUT_FILE.write_text(json.dumps({"meta": meta, "models": models_out, "global": {"by_variant": glob}},
                                   separators=(",", ":")))
    print(f"[egc] done in {dt:.1f}s ({len(models_out)} models, {n_err} task errors); "
          f"wrote {OUT_FILE.name} ({OUT_FILE.stat().st_size/1024:.0f} KB)", flush=True)


if __name__ == "__main__":
    main()
