#!/usr/bin/env python3
"""Per-reaction direction-sensitivity of growth for the 100-model panel, with the
four heuristic direction calls overlaid.

For every panel model, starting from the model's DEFAULT bounds, each reaction is
changed ONE AT A TIME to each of the four ModelSEED direction options and the
model's biomass growth (max flux) is re-solved:

    "<"  reverse      -> bounds (-1000, 0)
    ">"  forward      -> bounds (0, 1000)
    "="  reversible   -> bounds (-1000, 1000)
    "?"  unknown/off  -> bounds (0, 0)      <-- knocked out (user-chosen semantics)

All other reactions stay at their default bounds, so each number is the isolated
effect of that single reaction's direction on growth.  For each reaction we also
record where the four heuristics "send" it:

    default    -> the reaction's own in-place direction (from its bounds)
    jankowski  -> Jankowski_2008 column (group contribution)   \\
    flamholz   -> Flamholz_2012 column (eQuilibrator)           } literature TSV
    opus       -> LLM_Opus_4.8 column (Claude Opus 4.8)         /

Output: ``site/data/reaction_direction_effects_panel.json`` consumed by the
"Reaction-direction heuristics" chart under the site's Panel Models tab.

Pure cobra (no modelseedpy needed).  Run:
    PY=/mnt/homes/ctaylor/conda/miniforge3/envs/core_models_analysis/bin/python
    $PY scripts/build_reaction_direction_effects.py [--workers 32]
"""

from __future__ import annotations

import argparse
import datetime
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from direction_change_template_eval import build_base_to_model_index, direction_from_bounds
from template_quality_heuristics import load_literature_maps

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR", "/scratch/ctaylor/core_models_analysis"))
RESULTS_DIR = ANALYSIS_DIR / "results"
SITE_DATA = ANALYSIS_DIR / "site" / "data"
MODELS_DIR = Path("/scratch/ctaylor/core_models_kegg2")
PANEL_IDS = RESULTS_DIR / "selected_ids.txt"

# The 4 options the user asked to test, with '?' = off/knockout.
OPTION_BOUNDS = {
    "<": (-1000.0, 0.0),
    ">": (0.0, 1000.0),
    "=": (-1000.0, 1000.0),
    "?": (0.0, 0.0),
}
OPTIONS = ["<", ">", "=", "?"]
SCHEME_COLS = {"jankowski": "Jankowski_2008", "flamholz": "Flamholz_2012", "opus": "LLM_Opus_4.8"}

_LIT = None


def _find_biomass(model):
    for rid in ("bio1", "bio2", "biomass", "Biomass"):
        if rid in model.reactions:
            return model.reactions.get_by_id(rid)
    for r in model.reactions:
        if r.id.lower().startswith("bio") and not r.id.startswith("SK_"):
            return r
    return None


def _num(v):
    return round(float(v), 6) if (v is not None and v == v) else None  # nan -> None


def _init(lit):
    global _LIT
    import cobra
    cobra.Configuration().solver = "glpk"
    _LIT = lit


def eval_model(model_id):
    try:
        import cobra
        model = cobra.io.load_json_model(str(MODELS_DIR / f"{model_id}.json"))
        bio = _find_biomass(model)
        if bio is None:
            return {"model_id": model_id, "error": "no biomass reaction"}
        model.objective = bio
        model.objective_direction = "max"
        base_flux = _num(model.slim_optimize()) or 0.0

        base2mdl = build_base_to_model_index(model)
        jank = _LIT[SCHEME_COLS["jankowski"]]
        flam = _LIT[SCHEME_COLS["flamholz"]]
        opus = _LIT[SCHEME_COLS["opus"]]

        reactions = []
        for base, mids in base2mdl.items():
            for rid in mids:
                r = model.reactions.get_by_id(rid)
                default_dir = direction_from_bounds(r.lower_bound, r.upper_bound)
                g = {}
                for opt in OPTIONS:
                    lb, ub = OPTION_BOUNDS[opt]
                    try:
                        with model:
                            r.bounds = (lb, ub)
                            g[opt] = _num(model.slim_optimize())
                    except Exception:
                        g[opt] = None
                reactions.append({
                    "rxn": rid,
                    "base": base,
                    "dirs": {
                        "default": default_dir,
                        "jankowski": jank.get(base, "NA"),
                        "flamholz": flam.get(base, "NA"),
                        "opus": opus.get(base, "NA"),
                    },
                    "g": g,
                })
        # order by direction-sensitivity (spread of growth across options), desc
        def spread(rec):
            vals = [v for v in rec["g"].values() if v is not None]
            return (max(vals) - min(vals)) if vals else 0.0
        reactions.sort(key=spread, reverse=True)
        return {
            "model_id": model_id,
            "base_flux": round(base_flux, 6),
            "biomass_id": bio.id,
            "n_reactions": len(reactions),
            "reactions": reactions,
        }
    except Exception as exc:
        import traceback
        return {"model_id": model_id, "error": f"{type(exc).__name__}: {exc}",
                "trace": traceback.format_exc()[-500:]}


def build_global(models):
    """Cross-panel tally: reactions where the heuristics most often disagree, and
    reactions most often essential-when-off."""
    from collections import defaultdict
    disagree = defaultdict(int)
    off_ess = defaultdict(int)
    seen = defaultdict(int)
    for m in models.values():
        for rec in m["reactions"]:
            base = rec["base"]
            seen[base] += 1
            calls = {v for k, v in rec["dirs"].items() if v not in ("NA",)}
            if len({c for c in calls if c != "?"}) > 1:
                disagree[base] += 1
            goff = rec["g"].get("?")
            if m["base_flux"] > 1e-6 and goff is not None and goff < 1e-6:
                off_ess[base] += 1
    top_dis = sorted(disagree.items(), key=lambda kv: -kv[1])[:20]
    top_off = sorted(off_ess.items(), key=lambda kv: -kv[1])[:20]
    return {
        "most_disagreed": [{"base": b, "n_models": n} for b, n in top_dis],
        "most_off_essential": [{"base": b, "n_models": n} for b, n in top_off],
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--workers", type=int, default=32)
    p.add_argument("--limit", type=int, default=0)
    args = p.parse_args(argv)

    ids = PANEL_IDS.read_text().split()
    if args.limit:
        ids = ids[: args.limit]
    SITE_DATA.mkdir(parents=True, exist_ok=True)
    out_path = SITE_DATA / "reaction_direction_effects_panel.json"

    lit = load_literature_maps()
    print(f"panel models={len(ids)} workers={args.workers} "
          f"lit={{{', '.join(f'{k}:{len(v)}' for k, v in lit.items())}}}", flush=True)

    import multiprocessing as mp
    t0 = time.time()
    models = {}
    errors = []
    with mp.Pool(args.workers, initializer=_init, initargs=(lit,)) as pool:
        for i, rec in enumerate(pool.imap_unordered(eval_model, ids, chunksize=1), 1):
            if "error" in rec:
                errors.append((rec["model_id"], rec["error"]))
            else:
                models[rec["model_id"]] = {k: rec[k] for k in
                                           ("base_flux", "biomass_id", "n_reactions", "reactions")}
            if i % 20 == 0 or i == len(ids):
                print(f"  {i}/{len(ids)}  {i/max(time.time()-t0,1e-9):.1f} models/s", flush=True)

    out = {
        "generated": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "n_models": len(models),
        "n_errors": len(errors),
        "options": OPTIONS,
        "option_bounds": {k: list(v) for k, v in OPTION_BOUNDS.items()},
        "option_note": "'?' = unknown, tested as blocked/off (0,0)",
        "schemes": ["default", "jankowski", "flamholz", "opus"],
        "scheme_labels": {
            "default": "Default (model bounds)",
            "jankowski": "Jankowski (group contribution)",
            "flamholz": "Flamholz 2012 (eQuilibrator)",
            "opus": "Claude Opus 4.8",
        },
        "global": build_global(models),
        "models": models,
    }
    out_path.write_text(json.dumps(out))
    print(f"DONE: {len(models)} models, {len(errors)} errors -> {out_path} "
          f"({out_path.stat().st_size/1e6:.1f} MB)", flush=True)
    if errors:
        print("  errors:", errors[:5], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
