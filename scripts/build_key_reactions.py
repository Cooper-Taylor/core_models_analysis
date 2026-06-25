#!/usr/bin/env python3
"""Find each panel model's "key" reactions by direction-sensitivity.

For every panel model, probe every SEED reaction: flip it to each non-baseline
direction ('>', '<', '=') one at a time vs the cascade baseline and re-solve
growth. Reactions whose direction change causes a large |Δ growth| are the
model's "key" reactions -- a single (or small subset of) reaction(s) that
control whether the network grows. The per-reaction sweep is the heavy lifting
in ``growth_heuristics.key_reactions_one`` (model loaded once, ~hundreds of LP
solves per model, in memory).

Output (under ``site/data/``):

  - ``panel_key_reactions.json``
      {"models": {model_id: {base_flux, n_tested,
                             reactions:[{rxn, base_dir, best_dir, best_delta,
                                         by_dir:{dir:Δ}, severity}...]}},
       "global": [{rxn, n_models, max_severity, mean_severity}...]}
      Per-model reactions are sorted by severity (max |Δ growth|) descending and
      capped; "global" tallies reactions that are key across many panel models.

Inputs: site/data/baseline.json (heuristic-baseline reversibility map) and the
panel id list results/selected_ids.txt. Run any time after build_site_data.py.
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
OUT_FILE = SITE_DATA / "panel_key_reactions.json"

sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(MSDB_ROOT / "Libs" / "Python"))

import growth_heuristics as gh

_BASELINE: dict = {}


def _init(baseline_map, top_n):
    _BASELINE["map"] = baseline_map
    _BASELINE["top_n"] = top_n
    logging.getLogger("cobra").setLevel(logging.ERROR)


def _work(model_id):
    return gh.key_reactions_one(model_id, _BASELINE["map"], top_n=_BASELINE["top_n"])


def _round(x, n=4):
    return round(float(x), n)


def main(argv: Optional[list] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--workers", type=int, default=max(1, min(mp.cpu_count() - 1, 16)))
    ap.add_argument("--top-n", type=int, default=75,
                    help="max key reactions kept per model")
    ap.add_argument("--global-n", type=int, default=60,
                    help="max reactions in the cross-panel 'global' tally")
    args = ap.parse_args(argv)

    baseline_map = json.loads((SITE_DATA / "baseline.json").read_text())["map"]
    panel_ids = PANEL_FILE.read_text().split()
    print(f"[key] {len(panel_ids)} panel models; probing every reaction in all directions",
          flush=True)

    t0 = time.time()
    n_workers = max(1, min(args.workers, len(panel_ids)))
    if n_workers == 1:
        _init(baseline_map, args.top_n)
        results = [_work(m) for m in panel_ids]
    else:
        ctx = mp.get_context("fork")
        with ctx.Pool(n_workers, initializer=_init, initargs=(baseline_map, args.top_n)) as pool:
            results = list(pool.imap_unordered(_work, panel_ids, chunksize=2))
    dt = time.time() - t0

    models: dict = {}
    n_err = 0
    n_solves = 0
    tally: dict = {}  # seed -> {"n": int, "sevs": [..]}
    for res in results:
        if res.get("error"):
            n_err += 1
            continue
        mid = res["model_id"]
        n_solves += res.get("n_tested", 0)
        rxns = []
        for r in res["reactions"]:
            rxns.append({
                "rxn": r["rxn"],
                "base_dir": r["base_dir"],
                "best_dir": r["best_dir"],
                "best_delta": _round(r["best_delta"]),
                "by_dir": {k: _round(v) for k, v in r["by_dir"].items()},
                "severity": _round(r["severity"]),
            })
            t = tally.setdefault(r["rxn"], {"n": 0, "sevs": []})
            t["n"] += 1
            t["sevs"].append(r["severity"])
        models[mid] = {
            "base_flux": _round(res["base_flux"]),
            "n_tested": res.get("n_tested", 0),
            "reactions": rxns,
        }

    glob = [
        {"rxn": rxn,
         "n_models": t["n"],
         "max_severity": _round(max(t["sevs"])),
         "mean_severity": _round(sum(t["sevs"]) / len(t["sevs"]))}
        for rxn, t in tally.items()
    ]
    glob.sort(key=lambda d: (-d["n_models"], -d["max_severity"], d["rxn"]))
    glob = glob[:args.global_n]

    OUT_FILE.write_text(json.dumps({"models": models, "global": glob}, separators=(",", ":")))
    size_kb = OUT_FILE.stat().st_size / 1024
    print(f"[key] done in {dt:.1f}s ({len(models)} models, {n_err} errors, "
          f"~{n_solves} reaction-direction solves); wrote {OUT_FILE.name} ({size_kb:.0f} KB)",
          flush=True)
    if glob:
        top = glob[0]
        print(f"[key] most-frequently-key reaction: {top['rxn']} in {top['n_models']} models "
              f"(max |Δ growth| {top['max_severity']})")


if __name__ == "__main__":
    main()
