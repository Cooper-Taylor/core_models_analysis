#!/usr/bin/env python3
"""Find synthetic-lethal / synthetic-sick reaction pairs in each panel model.

Extends single-reaction essentiality to *pairs*: reactions that are dispensable
alone but jointly essential (Suthers et al. 2009; Fast-SL, Pratapa et al. 2015).
The per-model pairwise double-knockout sweep (restricted to the top flux-carrying
individually-non-essential candidates) lives in
``growth_heuristics.synthetic_lethal_one``.

Output (under ``site/data/``):
  - ``panel_synthetic_lethal.json``
      {"models": {model_id: {base_flux, n_candidates, n_pairs,
                            pairs:[{a,b,joint_delta,epistasis,single_a,single_b}]}},
       "global": [{pair:"a+b", a, b, n_models, max_abs_joint}...]}

Inputs: site/data/baseline.json + results/selected_ids.txt. Run after build_site_data.py.
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
OUT_FILE = SITE_DATA / "panel_synthetic_lethal.json"

sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(MSDB_ROOT / "Libs" / "Python"))

import growth_heuristics as gh

_BASELINE: dict = {}


def _init(baseline_map, n_cand, top_n):
    _BASELINE["map"] = baseline_map
    _BASELINE["n_cand"] = n_cand
    _BASELINE["top_n"] = top_n
    logging.getLogger("cobra").setLevel(logging.ERROR)


def _work(model_id):
    return gh.synthetic_lethal_one(model_id, _BASELINE["map"],
                                   n_cand=_BASELINE["n_cand"], top_n=_BASELINE["top_n"])


def _r(x, n=4):
    return round(float(x), n)


def main(argv: Optional[list] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--workers", type=int, default=max(1, min(mp.cpu_count() - 1, 16)))
    ap.add_argument("--n-cand", type=int, default=35)
    ap.add_argument("--top-n", type=int, default=40)
    ap.add_argument("--global-n", type=int, default=50)
    args = ap.parse_args(argv)

    baseline_map = json.loads((SITE_DATA / "baseline.json").read_text())["map"]
    panel_ids = PANEL_FILE.read_text().split()
    print(f"[SL] {len(panel_ids)} panel models; pairwise double-knockout over top-{args.n_cand} carriers",
          flush=True)

    t0 = time.time()
    nw = max(1, min(args.workers, len(panel_ids)))
    if nw == 1:
        _init(baseline_map, args.n_cand, args.top_n)
        results = [_work(m) for m in panel_ids]
    else:
        ctx = mp.get_context("fork")
        with ctx.Pool(nw, initializer=_init, initargs=(baseline_map, args.n_cand, args.top_n)) as pool:
            results = list(pool.imap_unordered(_work, panel_ids, chunksize=2))
    dt = time.time() - t0

    models, tally = {}, {}
    n_err = 0
    for res in results:
        if res.get("error"):
            n_err += 1
            continue
        mid = res["model_id"]
        pairs = [{"a": p["a"], "b": p["b"], "joint_delta": _r(p["joint_delta"]),
                  "epistasis": _r(p["epistasis"]), "single_a": _r(p["single_a"]),
                  "single_b": _r(p["single_b"])} for p in res["pairs"]]
        models[mid] = {"base_flux": _r(res["base_flux"]), "n_candidates": res["n_candidates"],
                       "n_pairs": res["n_pairs"], "pairs": pairs}
        for p in res["pairs"]:
            key = "+".join(sorted((p["a"], p["b"])))
            t = tally.setdefault(key, {"n": 0, "joints": []})
            t["n"] += 1
            t["joints"].append(p["joint_delta"])

    glob = [{"pair": k, "a": k.split("+")[0], "b": k.split("+")[1], "n_models": t["n"],
             "max_abs_joint": _r(max(abs(x) for x in t["joints"]))}
            for k, t in tally.items()]
    glob.sort(key=lambda d: (-d["n_models"], -d["max_abs_joint"], d["pair"]))
    glob = glob[:args.global_n]

    OUT_FILE.write_text(json.dumps({"models": models, "global": glob}, separators=(",", ":")))
    n_with = sum(1 for m in models.values() if m["n_pairs"])
    print(f"[SL] done in {dt:.1f}s ({len(models)} models, {n_err} errors, "
          f"{n_with} with >=1 SL pair); wrote {OUT_FILE.name} ({OUT_FILE.stat().st_size/1024:.0f} KB)",
          flush=True)


if __name__ == "__main__":
    main()
