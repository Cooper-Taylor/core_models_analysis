#!/usr/bin/env python3
"""Summarise presence_activity.json into the headline counts."""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from cma import paths  # noqa: E402

d = json.load(open(paths.results("direction_change_targets") / "presence_activity.json"))
targets: dict[str, str] = d["targets"]
panels: dict[str, dict] = d["panels"]

LABEL = {"core_kegg2": "core (core_kegg2)",
         "ms2_gsm": "full / GSM (glucose minimal)",
         "ms2_gsm_auxo": "full / GSM (auxotrophy media)"}

explicit = {r for r, c in targets.items() if c != "no change"}
uncertain = {r for r, c in targets.items() if c == "no change"}
N = len(targets)

print(f"{N} target reactions "
      f"({len(explicit)} explicit direction change, "
      f"{len(uncertain)} labelled 'no change' = uncertainty-driven)\n")

print(f"{'panel':34s} {'models':>7s} {'growers':>8s} {'present':>9s} {'active':>8s}")
for k, p in panels.items():
    print(f"{LABEL.get(k, k):34s} {p['n_models']:7d} {p['n_growers']:8d} "
          f"{p['targets_present_in_ge1_model']:6d}/{N} "
          f"{p['targets_active_in_ge1_model']:5d}/{N}")

# union across panels
present_any, active_any, active_fba_any = set(), set(), set()
for p in panels.values():
    for r, v in p["per_reaction"].items():
        if v["present"]:
            present_any.add(r)
        if v["active_pfba"]:
            active_any.add(r)
        if v["active_fba"]:
            active_fba_any.add(r)

print(f"\n{'union over all panels':34s} {'':7s} {'':8s} "
      f"{len(present_any):6d}/{N} {len(active_any):5d}/{N}")
print(f"  never present in any model : {N - len(present_any)}")
print(f"  present but never carrying flux : {len(present_any - active_any)}")

print("\n--- split by category (union over panels) ---")
for name, group in (("explicit direction change", explicit),
                    ("labelled 'no change'", uncertain)):
    print(f"  {name:28s} present {len(group & present_any):2d}/{len(group):2d}   "
          f"active {len(group & active_any):2d}/{len(group):2d}")

print(f"\n--- plain FBA vs pFBA (union) ---")
print(f"  active under pFBA      : {len(active_any)}")
print(f"  active under plain FBA : {len(active_fba_any)}")
extra = active_fba_any - active_any
if extra:
    print(f"  active only in plain FBA (alternate-optimum / loop artefacts): "
          f"{len(extra)} -> {sorted(extra)}")

print("\n--- per-reaction, reactions active somewhere ---")
print(f"{'reaction':10s} {'change':10s} " +
      " ".join(f"{LABEL.get(k,k).split('(')[0].strip()[:12]:>26s}" for k in panels))
hdr = " " * 21 + " ".join(f"{'present/active':>26s}" for _ in panels)
print(hdr)
rows = sorted(active_any, key=lambda r: -sum(panels[k]["per_reaction"][r]["active_pfba"]
                                             for k in panels))
for r in rows:
    cells = []
    for k in panels:
        v = panels[k]["per_reaction"][r]
        cells.append(f"{v['present']:11d} /{v['active_pfba']:12d}")
    print(f"{r:10s} {targets[r][:10]:10s} " + " ".join(cells))

print("\n--- present but never active anywhere ---")
quiet = sorted(present_any - active_any)
for r in quiet:
    tot = sum(panels[k]["per_reaction"][r]["present"] for k in panels)
    print(f"  {r:10s} {targets[r][:10]:10s} present in {tot} model-instances, never carries flux")

print("\n--- how much of each panel is touched ---")
for k, p in panels.items():
    print(f"  {LABEL.get(k,k):34s} models containing >=1 target: "
          f"{p['models_containing_any_target']:5d}/{p['n_models']}   "
          f"with >=1 target active: {p['models_with_any_target_active']:5d}")
