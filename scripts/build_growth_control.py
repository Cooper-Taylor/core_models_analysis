#!/usr/bin/env python3
"""Find each panel model's growth-controlling reactions (knockout essentiality).

The classic FBA single-reaction-deletion analysis (Edwards & Palsson 2000; Orth,
Thiele & Palsson 2010): for every panel model, block each reaction in turn and
re-solve biomass, recording the growth change. Reactions whose removal collapses
growth are **essential** (key to keeping growth HIGH); reactions whose removal
*raises* growth are **growth-limiting** (key to keeping growth LOW). Also records
the flux each reaction carries at the growth optimum and its LP reduced cost. The
per-reaction sweep lives in ``growth_heuristics.growth_control_one`` (model loaded
once, ~one extra LP per reaction, in memory).

Output (under ``site/data/``):

  - ``panel_growth_control.json``
      {"models": {model_id: {base_flux, n_tested, n_essential, n_limiting,
                             reactions:[{rxn, ko_delta, flux_opt, reduced_cost,
                                         kind}...]}},
       "global": [{rxn, n_essential_models, n_limiting_models, max_abs_ko,
                   mean_ko}...]}

Inputs: site/data/baseline.json + results/selected_ids.txt. Run after
build_site_data.py.
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
OUT_FILE = SITE_DATA / "panel_growth_control.json"

sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(MSDB_ROOT / "Libs" / "Python"))

import growth_heuristics as gh

_BASELINE: dict = {}


def _init(baseline_map, top_n):
    _BASELINE["map"] = baseline_map
    _BASELINE["top_n"] = top_n
    logging.getLogger("cobra").setLevel(logging.ERROR)


def _work(model_id):
    return gh.growth_control_one(model_id, _BASELINE["map"], top_n=_BASELINE["top_n"])


def _r(x, n=4):
    return round(float(x), n)


def main(argv: Optional[list] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--workers", type=int, default=max(1, min(mp.cpu_count() - 1, 16)))
    ap.add_argument("--top-n", type=int, default=90)
    ap.add_argument("--global-n", type=int, default=60)
    args = ap.parse_args(argv)

    baseline_map = json.loads((SITE_DATA / "baseline.json").read_text())["map"]
    panel_ids = PANEL_FILE.read_text().split()
    print(f"[ctrl] {len(panel_ids)} panel models; single-reaction knockout sweep", flush=True)

    t0 = time.time()
    nw = max(1, min(args.workers, len(panel_ids)))
    if nw == 1:
        _init(baseline_map, args.top_n)
        results = [_work(m) for m in panel_ids]
    else:
        ctx = mp.get_context("fork")
        with ctx.Pool(nw, initializer=_init, initargs=(baseline_map, args.top_n)) as pool:
            results = list(pool.imap_unordered(_work, panel_ids, chunksize=2))
    dt = time.time() - t0

    models, tally, met_tally = {}, {}, {}
    n_err = n_solves = 0
    for res in results:
        if res.get("error"):
            n_err += 1
            continue
        mid = res["model_id"]
        n_solves += res.get("n_tested", 0)
        rxns = [{"rxn": d["rxn"], "ko_delta": _r(d["ko_delta"]), "flux_opt": _r(d["flux_opt"]),
                 "reduced_cost": _r(d["reduced_cost"]), "kind": d["kind"]} for d in res["reactions"]]
        mets = [{"met": m["met"], "name": m["name"], "shadow_price": _r(m["shadow_price"])}
                for m in res.get("metabolites", [])]
        models[mid] = {"base_flux": _r(res["base_flux"]), "n_tested": res["n_tested"],
                       "n_essential": res["n_essential"], "n_limiting": res["n_limiting"],
                       "reactions": rxns, "metabolites": mets}
        for d in res["reactions"]:
            t = tally.setdefault(d["rxn"], {"ess": 0, "lim": 0, "kos": []})
            if d["kind"] == "essential":
                t["ess"] += 1
            elif d["kind"] == "limiting":
                t["lim"] += 1
            t["kos"].append(d["ko_delta"])
        for m in res.get("metabolites", []):
            mt = met_tally.setdefault(m["met"], {"name": m["name"], "n": 0, "sps": []})
            mt["n"] += 1
            mt["sps"].append(abs(m["shadow_price"]))

    glob = [{"rxn": r, "n_essential_models": t["ess"], "n_limiting_models": t["lim"],
             "max_abs_ko": _r(max(abs(x) for x in t["kos"])),
             "mean_ko": _r(sum(t["kos"]) / len(t["kos"]))}
            for r, t in tally.items() if t["ess"] or t["lim"]]
    glob.sort(key=lambda d: (-d["n_essential_models"], -d["max_abs_ko"], d["rxn"]))
    glob = glob[:args.global_n]

    met_glob = [{"met": m, "name": v["name"], "n_models": v["n"],
                 "max_abs_sp": _r(max(v["sps"])), "mean_abs_sp": _r(sum(v["sps"]) / len(v["sps"]))}
                for m, v in met_tally.items()]
    met_glob.sort(key=lambda d: (-d["n_models"], -d["max_abs_sp"], d["met"]))
    met_glob = met_glob[:args.global_n]

    OUT_FILE.write_text(json.dumps({"models": models, "global": glob, "metabolites_global": met_glob},
                                   separators=(",", ":")))
    print(f"[ctrl] done in {dt:.1f}s ({len(models)} models, {n_err} errors, "
          f"~{n_solves} knockouts); wrote {OUT_FILE.name} ({OUT_FILE.stat().st_size/1024:.0f} KB)",
          flush=True)
    if glob:
        print(f"[ctrl] most broadly essential: {glob[0]['rxn']} "
              f"(essential in {glob[0]['n_essential_models']} models)")


if __name__ == "__main__":
    main()
