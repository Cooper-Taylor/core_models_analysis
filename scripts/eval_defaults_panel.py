#!/usr/bin/env python3
"""Evaluate every ReversibilityConfig variant against the *new* fixed baseline
(H2+H3 adopted) on the *new* 100-model descriptive panel, to inform the
choice of default cascade configuration.

For each variant it reports:
  - n_changed   : reactions whose EQ direction differs from baseline (all 56k)
  - transitions : the dominant baseline->variant direction shifts
  - panel grow-flips: grew->not / not->grew on the 100 panel models
  - panel mean biomass flux vs baseline

Self-contained: computes cascade maps via reversibility_lib and panel FBA via
growth_heuristics (full rebind). Does not touch the notebook kbcache or any
on-disk model JSON.
"""
from __future__ import annotations

import json
import os
import sys
import time
from collections import Counter
from pathlib import Path

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR", "/scratch/ctaylor/core_models_analysis"))
MSDB_ROOT = os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase")
sys.path.insert(0, str(ANALYSIS_DIR / "scripts"))
sys.path.insert(0, MSDB_ROOT + "/Libs/Python")

import reversibility_lib as lib       # noqa: E402
import variant_catalog as vc          # noqa: E402
import growth_heuristics as gh        # noqa: E402
from BiochemPy import Reactions       # noqa: E402

N_WORKERS = int(os.environ.get("EVAL_WORKERS", "32"))


def cascade_map(rxns, cfg) -> dict:
    """EQ-level cascade with the GC-first pre-roll, returning {rxn_id: rev}."""
    out = lib.run_cascade(rxns, db_level="EQ", cfg=cfg, gc_first=True)
    return {r: rev for r, (status, rev) in out.items()}


def main():
    print("loading reactions ...", flush=True)
    rxns = Reactions().loadReactions()
    panel_ids = (ANALYSIS_DIR / "results" / "selected_ids.txt").read_text().split()
    print(f"panel size: {len(panel_ids)}", flush=True)

    print("baseline cascade ...", flush=True)
    base_map = cascade_map(rxns, lib.ReversibilityConfig())

    print("baseline panel FBA ...", flush=True)
    base_fba = gh.run_panel(panel_ids, reversibility_map=base_map,
                            baseline_map=None, n_workers=N_WORKERS)
    base_by_id = {r["model_id"]: r for r in base_fba}
    base_growers = sum(1 for r in base_fba if r["grows"])
    base_mean = sum(r["growth_flux"] for r in base_fba) / len(base_fba)
    print(f"  baseline growers on panel: {base_growers}/{len(panel_ids)}  mean_flux={base_mean:.4g}", flush=True)

    report = {
        "panel_size": len(panel_ids),
        "baseline": {"panel_growers": base_growers, "panel_mean_flux": base_mean},
        "variants": {},
    }

    for v in vc.VARIANTS:
        tag = v["tag"]
        if tag == "baseline":
            continue
        t0 = time.time()
        cfg = v["cfg"]()
        vmap = cascade_map(rxns, cfg)
        # all-DB direction changes vs baseline
        nchg = 0
        trans = Counter()
        for r in base_map:
            b, n = base_map[r], vmap.get(r)
            if b != n:
                nchg += 1
                trans[(b, n)] += 1
        # panel FBA
        vfba = gh.run_panel(panel_ids, reversibility_map=vmap,
                            baseline_map=None, n_workers=N_WORKERS)
        g2n = n2g = 0
        vmean = sum(r["growth_flux"] for r in vfba) / len(vfba)
        vgrow = 0
        for r in vfba:
            b = base_by_id[r["model_id"]]
            vgrow += bool(r["grows"])
            if b["grows"] and not r["grows"]:
                g2n += 1
            if r["grows"] and not b["grows"]:
                n2g += 1
        top = ", ".join(f"{k[0]}->{k[1]}:{c}" for k, c in trans.most_common(4))
        report["variants"][tag] = {
            "n_changed_rxns": nchg,
            "top_transitions": top,
            "panel_growers": vgrow,
            "panel_grew_to_not": g2n,
            "panel_not_to_grew": n2g,
            "panel_net_grow_delta": vgrow - base_growers,
            "panel_mean_flux": vmean,
        }
        print(f"[{tag}] nchg={nchg} grow={vgrow}({vgrow-base_growers:+d}) "
              f"g2n={g2n} n2g={n2g} mean_flux={vmean:.4g} "
              f"({time.time()-t0:.0f}s)  {top}", flush=True)

    out = ANALYSIS_DIR / "results" / "defaults_panel_eval.json"
    out.write_text(json.dumps(report, indent=2))
    print(f"\nwrote {out}", flush=True)


if __name__ == "__main__":
    main()
