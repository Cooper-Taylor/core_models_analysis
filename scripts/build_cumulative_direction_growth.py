#!/usr/bin/env python3
"""Cumulative growth trajectory of a heuristic's reaction-direction changes.

For each panel model and each heuristic (Jankowski / Flamholz / Opus 4.8):

  1. Starting from the model's DEFAULT bounds, take every reaction whose heuristic
     direction differs from the model's own direction, and measure the *individual*
     effect on biomass growth of making only that one change (single-reaction diff
     vs the default baseline).
  2. Sort those reactions by that individual growth effect (largest first).
  3. Apply the changes ONE ON ANOTHER in that sorted order and record the model's
     growth after each cumulative prefix:
        [baseline, r1, r1+r2, r1+r2+r3, ...]

Both the ranking (step 1) and the cumulative application (step 3) are driven by the
KBUtilLib feature we have been using -- ``MSTemplateUtils.diff_template_evaluation``:
  * step 1 uses ``mode="independent"`` (each change vs the shared default baseline),
  * step 3 uses ``mode="cumulative"`` (each change vs the previous state).
We subclass the offline runner so ``_evaluate_model_quality`` reports biomass flux
(a single ``slim_optimize``), and read the growth of the baseline + each step out of
the sequence of evaluation calls diff_template_evaluation performs.

Output: ``site/data/cumulative_direction_growth_panel.json`` -> the "Cumulative
direction changes" line chart under the site's Panel Models tab.

    PY=/mnt/homes/ctaylor/conda/miniforge3/envs/core_models_analysis/bin/python
    $PY scripts/build_cumulative_direction_growth.py [--workers 48]
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
from direction_change_template_eval import (
    bounds_for_direction, build_base_to_model_index, direction_from_bounds, _find_biomass,
)
from template_quality_heuristics import load_literature_maps

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR", "/scratch/ctaylor/core_models_analysis"))
RESULTS_DIR = ANALYSIS_DIR / "results"
SITE_DATA = ANALYSIS_DIR / "site" / "data"
MODELS_DIR = Path("/scratch/ctaylor/core_models_kegg2")
PANEL_IDS = RESULTS_DIR / "selected_ids.txt"

SCHEMES = [
    {"key": "jankowski", "col": "Jankowski_2008", "label": "Jankowski (group contribution)"},
    {"key": "flamholz", "col": "Flamholz_2012", "label": "Flamholz 2012 (eQuilibrator)"},
    {"key": "opus", "col": "LLM_Opus_4.8", "label": "Claude Opus 4.8"},
]

_RUNNER = None
_MU = None
_LIT = None


def _init(lit):
    global _RUNNER, _MU, _LIT
    import cobra
    cobra.Configuration().solver = "glpk"
    from direction_change_template_eval import make_runner_class
    cls, MU = make_runner_class()

    class GrowthReportingEval(cls):
        # Report ONLY biomass growth so diff_template_evaluation is fast, and record
        # each evaluation's growth into self._trace (baseline first, then each step).
        def _evaluate_model_quality(self, mdlutl, rich_media=None, minimal_media=None):
            cm = mdlutl.model
            bio = _find_biomass(cm)
            g = 0.0
            if bio is not None:
                cm.objective = bio
                cm.objective_direction = "max"
                v = cm.slim_optimize()
                g = float(v) if (v is not None and v == v) else 0.0
            self._trace.append(round(g, 6))
            return {"growth": round(g, 6), "template_metadata": {"id": getattr(cm, "id", "m")}}

    r = GrowthReportingEval.build()
    r._trace = []
    _RUNNER, _MU, _LIT = r, MU, lit


def eval_model(model_id):
    try:
        import cobra
        model = cobra.io.load_json_model(str(MODELS_DIR / f"{model_id}.json"))
        base2mdl = build_base_to_model_index(model)
        rxn_by_id = {r.id: r for r in model.reactions}

        # explicit baseline growth (all reactions at default)
        mc = model.copy()
        bio = _find_biomass(mc)
        base_growth = 0.0
        if bio is not None:
            mc.objective = bio
            mc.objective_direction = "max"
            v = mc.slim_optimize()
            base_growth = round(float(v), 6) if (v is not None and v == v) else 0.0

        out = {"model_id": model_id, "baseline": base_growth, "heuristics": {}}

        for s in SCHEMES:
            dmap = _LIT[s["col"]]
            perts, meta = [], []
            for base, mids in base2mdl.items():
                hd = dmap.get(base)
                if not hd or hd in ("NA", "?"):
                    continue
                b = bounds_for_direction(hd)
                if b is None:
                    continue
                for mid in mids:
                    r = rxn_by_id[mid]
                    dd = direction_from_bounds(r.lower_bound, r.upper_bound)
                    if hd == dd:
                        continue  # heuristic agrees with the model's own direction
                    perts.append({"op": "modify", "reaction_id": mid,
                                  "lower_bound": b[0], "upper_bound": b[1]})
                    meta.append({"rxn": mid, "base": base, "to": hd, "from": dd})

            if not perts:
                out["heuristics"][s["key"]] = {"n_changed": 0, "reactions": [],
                                               "cumulative": [base_growth]}
                continue

            # (1) individual effects vs default baseline  -- diff_template_evaluation independent
            _RUNNER._trace = []
            _RUNNER.diff_template_evaluation(model, perts, mode="independent", baseline_report=None)
            tr = _RUNNER._trace                       # [baseline, ind_1, ind_2, ...]
            deltas = [round(tr[1 + i] - tr[0], 6) for i in range(len(perts))]

            # (2) sort by individual growth effect, largest first
            order = sorted(range(len(perts)), key=lambda i: deltas[i], reverse=True)
            sorted_perts = [perts[i] for i in order]
            sorted_meta = [{**meta[i], "delta": deltas[i]} for i in order]

            # (3) apply one-on-another in that order  -- diff_template_evaluation cumulative
            _RUNNER._trace = []
            _RUNNER.diff_template_evaluation(model, sorted_perts, mode="cumulative", baseline_report=None)
            cum = _RUNNER._trace                      # [baseline, cum_1, cum_12, ...]

            out["heuristics"][s["key"]] = {
                "n_changed": len(perts),
                "reactions": sorted_meta,
                "cumulative": [round(x, 6) for x in cum],
                "final": round(cum[-1], 6),
                "peak": round(max(cum), 6),
                "peak_at": int(max(range(len(cum)), key=lambda i: cum[i])),
            }
        return out
    except Exception as exc:
        import traceback
        return {"model_id": model_id, "error": f"{type(exc).__name__}: {exc}",
                "trace": traceback.format_exc()[-500:]}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--workers", type=int, default=48)
    p.add_argument("--limit", type=int, default=0)
    args = p.parse_args(argv)

    ids = PANEL_IDS.read_text().split()
    if args.limit:
        ids = ids[: args.limit]
    SITE_DATA.mkdir(parents=True, exist_ok=True)
    out_path = SITE_DATA / "cumulative_direction_growth_panel.json"

    lit = load_literature_maps()
    print(f"panel models={len(ids)} workers={args.workers}", flush=True)

    import multiprocessing as mp
    t0 = time.time()
    models, errors = {}, []
    with mp.Pool(args.workers, initializer=_init, initargs=(lit,)) as pool:
        for i, rec in enumerate(pool.imap_unordered(eval_model, ids, chunksize=1), 1):
            if "error" in rec:
                errors.append((rec["model_id"], rec["error"]))
            else:
                models[rec["model_id"]] = {"baseline": rec["baseline"], "heuristics": rec["heuristics"]}
            if i % 20 == 0 or i == len(ids):
                print(f"  {i}/{len(ids)}  {i/max(time.time()-t0,1e-9):.1f} models/s", flush=True)

    out = {
        "generated": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "n_models": len(models),
        "n_errors": len(errors),
        "schemes": SCHEMES,
        "note": "cumulative[i] = biomass growth after applying the i highest-individual-effect "
                "direction changes (sorted desc), one on another, via "
                "MSTemplateUtils.diff_template_evaluation(mode='cumulative').",
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
