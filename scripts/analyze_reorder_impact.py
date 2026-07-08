#!/usr/bin/env python3
"""Analyze what the reversibility heuristic REORDER changed, before vs after.

The reorder (ATP-synthase + ABC-transporter heuristics moved ahead of the
general MdeltaG stored-bounds rule) is applied in both MSDB
(Estimate_Reaction_Reversibility DEFAULT_HEURISTICS, branch claude-changes) and
in this repo's port (reversibility_lib.estimate_one). This script isolates the
reorder's effect by running the port cascade with the PRE-reorder code (HEAD of
reversibility_lib.py) vs the current (post-reorder) code, on the same MSDB data:

  1. Direction transitions: every reaction whose baseline cascade direction
     changed, with its name, old->new direction, and the heuristic that now
     fires (ATPS / ABCT).
  2. Model footprint: how many of the 5,683 core models (and the 100-model
     descriptive panel) contain each changed reaction.
  3. Panel FBA impact: 100-model growth/flux with the OLD baseline map vs the
     NEW baseline map (full rebind, baseline_map=None), so only the changed
     reactions can move flux.

Writes results/reorder_impact.json and prints a summary. Read-only on MSDB and
core_models_kegg2.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from collections import Counter

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(SCRIPTS)
MSDB = os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, os.path.join(MSDB, "Libs", "Python"))

import reversibility_lib as new_lib            # post-reorder (current working tree)
import growth_heuristics as gh
from seed_annotation import normalize_seed_id
from BiochemPy import Reactions

OUT = os.path.join(ROOT, "results", "reorder_impact.json")
PANEL_IDS = os.path.join(ROOT, "results", "selected_ids.txt")
HEAD_TMP = os.path.join(SCRIPTS, "_revlib_head.py")


def _load_old_lib():
    """Import the pre-reorder reversibility_lib from git HEAD as a sibling module."""
    src = subprocess.run(
        ["git", "show", "HEAD:scripts/reversibility_lib.py"],
        cwd=ROOT, capture_output=True, text=True, check=True).stdout
    with open(HEAD_TMP, "w") as fh:
        fh.write(src)
    import _revlib_head as old_lib  # noqa: E402
    return old_lib


def _cascade_map(lib, rxns):
    """{rxn_id: (status, direction)} from a fresh load (run_cascade mutates the
    reversibility field via its GC pre-pass, so callers pass a fresh dict)."""
    return lib.run_cascade(rxns, db_level="EQ", cfg=lib.ReversibilityConfig(),
                           gc_first=True)


def main():
    old_lib = _load_old_lib()
    try:
        new_casc = _cascade_map(new_lib, Reactions().loadReactions())
        rxns = Reactions().loadReactions()      # fresh (un-mutated) for names + old pass
        names = {rid: r.get("name", "") for rid, r in rxns.items()}
        is_tx = {rid: r.get("is_transport") in (1, "1", True) for rid, r in rxns.items()}
        old_casc = _cascade_map(old_lib, rxns)
    finally:
        if os.path.exists(HEAD_TMP):
            os.remove(HEAD_TMP)

    old_map = {r: rev for r, (_s, rev) in old_casc.items()}
    new_map = {r: rev for r, (_s, rev) in new_casc.items()}
    new_status = {r: s for r, (s, _rev) in new_casc.items()}

    changed = {r: (old_map[r], new_map[r]) for r in new_map if old_map[r] != new_map[r]}
    print(f"[reorder] reactions: {len(new_map)}  changed: {len(changed)}")
    print("[reorder] transitions:",
          dict(Counter(f"{a}->{b}" for a, b in changed.values())))

    def heuristic_of(status):
        if status.startswith("ATPS"):
            return "ATP synthase"
        if status.startswith("ABCT"):
            return "ABC transporter"
        return status.split(":")[0]

    by_heur = Counter(heuristic_of(new_status[r]) for r in changed)
    print("[reorder] deciding heuristic (new):", dict(by_heur))

    # --- model footprint over all 5,683 core models + the 100-model panel ---
    panel = set(open(PANEL_IDS).read().split())
    changed_set = set(changed)
    all_count = Counter()
    panel_count = Counter()
    model_files = sorted(gh.MODELS_DIR.glob("*.json"))
    for f in model_files:
        mid = f.stem
        seeds = set()
        for rxn in json.load(open(f)).get("reactions", []):
            anno = rxn.get("annotation") or {}
            s = anno.get("seed.reaction")
            if isinstance(s, list):
                s = s[0] if s else None
            s = normalize_seed_id(s)
            if s in changed_set:
                seeds.add(s)
        for s in seeds:
            all_count[s] += 1
            if mid in panel:
                panel_count[s] += 1
    # models containing >=1 changed reaction (all + panel)
    models_with_change_all = set()
    models_with_change_panel = set()
    for f in model_files:
        mid = f.stem
        for rxn in json.load(open(f)).get("reactions", []):
            anno = rxn.get("annotation") or {}
            s = anno.get("seed.reaction")
            if isinstance(s, list):
                s = s[0] if s else None
            s = normalize_seed_id(s)
            if s in changed_set:
                models_with_change_all.add(mid)
                if mid in panel:
                    models_with_change_panel.add(mid)
                break

    reactions = []
    for r in sorted(changed):
        o, n = changed[r]
        reactions.append({
            "rxn_id": r, "name": names.get(r, ""),
            "old": o, "new": n,
            "heuristic": heuristic_of(new_status[r]),
            "status": new_status[r],
            "is_transport": bool(is_tx.get(r)),
            "n_models_all": all_count.get(r, 0),
            "n_models_panel": panel_count.get(r, 0),
        })

    # --- panel FBA: OLD baseline map vs NEW baseline map (full rebind) ---
    panel_ids = sorted(panel)
    print(f"[reorder] panel FBA on {len(panel_ids)} models (old vs new baseline)...")
    res_old = {x["model_id"]: x for x in gh.run_panel(panel_ids, reversibility_map=old_map,
                                                      baseline_map=None, n_workers=16)}
    res_new = {x["model_id"]: x for x in gh.run_panel(panel_ids, reversibility_map=new_map,
                                                      baseline_map=None, n_workers=16)}
    fba_rows, n_grow_flip, n_flux_change = [], 0, 0
    EPS = 1e-6
    for mid in panel_ids:
        b, v = res_old.get(mid, {}), res_new.get(mid, {})
        bf, vf = float(b.get("growth_flux", 0) or 0), float(v.get("growth_flux", 0) or 0)
        bg, vg = bool(b.get("grows")), bool(v.get("grows"))
        d = vf - bf
        if bg != vg:
            n_grow_flip += 1
        if abs(d) > EPS:
            n_flux_change += 1
        if bg != vg or abs(d) > EPS:
            fba_rows.append({"model_id": mid, "old_flux": bf, "new_flux": vf,
                             "delta_flux": d, "old_grows": bg, "new_grows": vg})
    fba_rows.sort(key=lambda x: -abs(x["delta_flux"]))

    summary = {
        "n_reactions_total": len(new_map),
        "n_changed": len(changed),
        "transitions": dict(Counter(f"{a}->{b}" for a, b in changed.values())),
        "deciding_heuristic": dict(by_heur),
        "n_models_all_with_change": len(models_with_change_all),
        "n_models_panel_with_change": len(models_with_change_panel),
        "panel_fba": {
            "n_models": len(panel_ids),
            "n_grow_flips": n_grow_flip,
            "n_flux_changes": n_flux_change,
            "changed_models": fba_rows,
        },
        "reactions": reactions,
    }
    with open(OUT, "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"[reorder] models containing >=1 changed rxn: "
          f"{len(models_with_change_all)}/{len(model_files)} all, "
          f"{len(models_with_change_panel)}/{len(panel)} panel")
    print(f"[reorder] panel FBA: {n_grow_flip} grow-flips, {n_flux_change} flux changes "
          f"(of {len(panel_ids)})")
    print(f"[reorder] wrote {OUT}")


if __name__ == "__main__":
    main()
