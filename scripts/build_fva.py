#!/usr/bin/env python3
"""Flux variability analysis (FVA) of each panel model at near-optimal growth.

Classifies every reaction as blocked / flux-forced (obligate for optimal growth)
/ flexible (Mahadevan & Schilling 2003). Per-model sweep lives in
``growth_heuristics.fva_one``.

Output (under ``site/data/``):
  - ``panel_fva.json``
      {"models": {model_id: {base_flux, n_blocked, n_forced, n_flexible,
                            reactions:[{rxn,min,max,span,kind}]}},
       "global": [{rxn, n_forced_models, n_blocked_models, max_span}...]}

Inputs: site/data/baseline.json + results/selected_ids.txt. FVA is the slowest
sweep (~1-2 min for the 100-model panel); run in the background.
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
MSDB_ROOT = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase"))
PANEL_FILE = ANALYSIS_ROOT / "results" / "selected_ids.txt"
OUT_FILE = SITE_DATA / "panel_fva.json"

sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(MSDB_ROOT / "Libs" / "Python"))

import growth_heuristics as gh

_CFG: dict = {}


def _init(baseline_map, fraction, top_n):
    _CFG["map"] = baseline_map
    _CFG["fraction"] = fraction
    _CFG["top_n"] = top_n
    logging.getLogger("cobra").setLevel(logging.ERROR)


def _work(model_id):
    return gh.fva_one(model_id, _CFG["map"], fraction=_CFG["fraction"], top_n=_CFG["top_n"])


def _r(x, n=4):
    return round(float(x), n)


def main(argv: Optional[list] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--workers", type=int, default=max(1, min(mp.cpu_count() - 1, 16)))
    ap.add_argument("--fraction", type=float, default=0.99)
    ap.add_argument("--top-n", type=int, default=60)
    ap.add_argument("--global-n", type=int, default=60)
    args = ap.parse_args(argv)

    baseline_map = json.loads((SITE_DATA / "baseline.json").read_text())["map"]
    panel_ids = PANEL_FILE.read_text().split()
    print(f"[FVA] {len(panel_ids)} panel models @ fraction_of_optimum={args.fraction}", flush=True)

    t0 = time.time()
    nw = max(1, min(args.workers, len(panel_ids)))
    if nw == 1:
        _init(baseline_map, args.fraction, args.top_n)
        results = [_work(m) for m in panel_ids]
    else:
        ctx = mp.get_context("fork")
        with ctx.Pool(nw, initializer=_init, initargs=(baseline_map, args.fraction, args.top_n)) as pool:
            results = list(pool.imap_unordered(_work, panel_ids, chunksize=1))
    dt = time.time() - t0

    models, tally = {}, {}
    n_err = 0
    for res in results:
        if res.get("error"):
            n_err += 1
            continue
        mid = res["model_id"]
        rxns = [{"rxn": d["rxn"], "min": _r(d["min"]), "max": _r(d["max"]),
                 "span": _r(d["span"]), "kind": d["kind"]} for d in res["reactions"]]
        models[mid] = {"base_flux": _r(res["base_flux"]), "n_blocked": res["n_blocked"],
                       "n_forced": res["n_forced"], "n_flexible": res["n_flexible"],
                       "reactions": rxns}
        for d in res["reactions"]:
            t = tally.setdefault(d["rxn"], {"f": 0, "b": 0, "span": 0.0})
            if d["kind"] == "flux_forced":
                t["f"] += 1
            elif d["kind"] == "blocked":
                t["b"] += 1
            t["span"] = max(t["span"], abs(d["span"]))

    glob = [{"rxn": r, "n_forced_models": t["f"], "n_blocked_models": t["b"],
             "max_span": _r(t["span"])} for r, t in tally.items() if t["f"] or t["b"]]
    glob.sort(key=lambda d: (-d["n_forced_models"], -d["n_blocked_models"], d["rxn"]))
    glob = glob[:args.global_n]

    OUT_FILE.write_text(json.dumps({"models": models, "global": glob}, separators=(",", ":")))
    print(f"[FVA] done in {dt:.1f}s ({len(models)} models, {n_err} errors); "
          f"wrote {OUT_FILE.name} ({OUT_FILE.stat().st_size/1024:.0f} KB)", flush=True)


if __name__ == "__main__":
    main()
