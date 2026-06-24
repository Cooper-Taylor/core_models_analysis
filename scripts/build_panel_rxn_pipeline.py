#!/usr/bin/env python3
"""Precompute the per-reaction perturbation pipeline for the Panel Models tab.

For every (panel model, non-baseline variant) pair, decompose the variant's
whole-model ATP-flux effect into per-reaction contributions:

  * **single** -- flip each variant-changed reaction (that exists in the model)
    one at a time vs the heuristic baseline, recording its Δ ATP flux;
  * **cumulative** -- sort those reactions by |Δ ATP flux| descending and
    re-apply them in that order, recording the running Δ ATP flux.

The heavy lifting (load model once, perturb in memory) lives in
``growth_heuristics.rxn_pipeline_one``; this script just builds the task list,
parallelizes, and writes the result.

Output (under ``site/data/``):

  - ``panel_model_rxn_pipeline.json``
      {model_id: {variant_tag: {base_flux, singles:[{rxn,delta}...],
                                cumulative:[{rxn,delta}...]}}}
      singles + cumulative share the same |delta|-descending reaction order;
      cumulative[-1].delta == that variant's whole-model delta (cross-check).

Inputs (all produced by build_site_data.py): site/data/baseline.json (the
heuristic-baseline reversibility map), site/data/variants/*.json (per-variant
diffs), site/data/panel_rxnsets.json (model -> seed reaction set), and the
panel id list results/selected_ids.txt.

Runtime: ~1-2 min on the panel (each (model,variant) pair loads its model once
and solves ~2*n_changed small LPs in memory).
"""

from __future__ import annotations

import argparse
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
MODELS_DIR = ANALYSIS_ROOT / "data" / "core_models_kegg2"
MSDB_ROOT = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase"))
PANEL_FILE = ANALYSIS_ROOT / "results" / "selected_ids.txt"
OUT_FILE = SITE_DATA / "panel_model_rxn_pipeline.json"

sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(MSDB_ROOT / "Libs" / "Python"))

import growth_heuristics as gh
from seed_annotation import seed_id  # normalizes the _c-suffix bug


def model_rxnset(model_id: str) -> set:
    """SEED reaction-id set for a model, parsed straight from its JSON.

    Sourced from the model JSON (not site/data/panel_rxnsets.json, whose
    kbcache origin is missing 17 of the 100 panel models) so every panel
    model is covered. Matches the seed_id canonicalization that
    growth_heuristics.override_bounds uses on the cobra model."""
    p = MODELS_DIR / f"{model_id}.json"
    if not p.exists():
        return set()
    d = json.loads(p.read_text())
    return {s for r in d.get("reactions", []) if (s := seed_id(r))}


def load_variant_diffs() -> dict:
    """Return {tag: {seed_rxn: variant_direction}} for every non-baseline variant."""
    out = {}
    for vfile in sorted((SITE_DATA / "variants").glob("*.json")):
        tag = vfile.stem
        if tag == "baseline":
            continue
        payload = json.loads(vfile.read_text())
        out[tag] = {d["rxn"]: d["new"] for d in payload.get("diffs", [])}
    return out


# ---- worker (baseline map shared once via initializer) --------------------
_BASELINE: dict = {}


def _init(baseline_map):
    _BASELINE["map"] = baseline_map
    logging.getLogger("cobra").setLevel(logging.ERROR)


def _work(task):
    model_id, tag, changed_dirs = task
    res = gh.rxn_pipeline_one(model_id, _BASELINE["map"], changed_dirs)
    return model_id, tag, res


def _round(x, n=4):
    return round(float(x), n)


def main(argv: Optional[list] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--workers", type=int,
                    default=max(1, min(mp.cpu_count() - 1, 16)))
    ap.add_argument("--variants", nargs="*", default=None,
                    help="Only process these variant tags (default: all)")
    args = ap.parse_args(argv)

    baseline_map = json.loads((SITE_DATA / "baseline.json").read_text())["map"]
    panel_ids = PANEL_FILE.read_text().split()
    rxnsets = {m: model_rxnset(m) for m in panel_ids}
    n_empty = sum(1 for m in panel_ids if not rxnsets[m])
    print(f"[pipeline] rxnsets from model JSON ({n_empty} panel models with no reactions)",
          flush=True)
    variant_diffs = load_variant_diffs()
    if args.variants:
        variant_diffs = {t: d for t, d in variant_diffs.items() if t in args.variants}
    print(f"[pipeline] {len(panel_ids)} panel models x {len(variant_diffs)} variants",
          flush=True)

    # one task per (model, variant) pair that actually has changed reactions
    tasks = []
    for mid in panel_ids:
        model_rxns = rxnsets.get(mid, set())
        for tag, diff in variant_diffs.items():
            changed = {r: d for r, d in diff.items() if r in model_rxns}
            if changed:
                tasks.append((mid, tag, changed))
    print(f"[pipeline] {len(tasks)} (model,variant) tasks with >=1 changed rxn",
          flush=True)

    t0 = time.time()
    n_workers = max(1, min(args.workers, len(tasks) or 1))
    if n_workers == 1:
        _init(baseline_map)
        results = [_work(t) for t in tasks]
    else:
        ctx = mp.get_context("fork")
        with ctx.Pool(n_workers, initializer=_init, initargs=(baseline_map,)) as pool:
            results = list(pool.imap_unordered(_work, tasks, chunksize=4))
    dt = time.time() - t0

    out: dict = {}
    n_err = 0
    for mid, tag, res in results:
        if res.get("error"):
            n_err += 1
            continue
        out.setdefault(mid, {})[tag] = {
            "base_flux": _round(res["base_flux"]),
            "singles": [{"rxn": s["rxn"], "delta": _round(s["delta"])}
                        for s in res["singles"]],
            "cumulative": [{"rxn": c["rxn"], "delta": _round(c["delta"])}
                           for c in res["cumulative"]],
        }

    OUT_FILE.write_text(json.dumps(out, separators=(",", ":")))
    size_kb = OUT_FILE.stat().st_size / 1024
    print(f"[pipeline] done in {dt:.1f}s ({len(tasks)} tasks, {n_err} errors); "
          f"wrote {OUT_FILE.name} ({len(out)} models, {size_kb:.0f} KB)", flush=True)


if __name__ == "__main__":
    main()
