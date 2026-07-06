#!/usr/bin/env python3
"""Evaluate every core model's *template quality* under four reaction-direction
heuristics, using KBUtilLib's new ``MSTemplateUtils`` battery, and emit a JSON
snapshot for the core_models_analysis website.

For each core model we build the FBA/FVA quality report (via
``MSTemplateUtils._evaluate_model_quality``: reaction classification, closed-mode
loops, Biolog growth, producible/consumable-metabolite sweeps) under four
direction-assignment schemes, then measure how the qualities change:

  (1) default      -- the model's own in-place bounds (no override)
  (2) jankowski    -- Jankowski (group contribution, Henry-2007 feasibility rule);
                      column ``Jankowski_2008`` of the literature TSV
  (3) flamholz     -- Flamholz 2012 (eQuilibrator reversibility index);
                      column ``Flamholz_2012``
  (4) opus         -- Claude Opus 4.8 directionality; column ``LLM_Opus_4.8``

Sources (2)-(4) come from ``results/reaction_directions_literature_vs_llm.tsv``
(produced by ``estimate_directions_literature.py``); each method uses its OWN
directionality rule, not the MSDB cascade.  A direction is applied by rewriting
the matching in-model reaction's bounds ( ">"->(0,1000), "<"->(-1000,0),
"="->(-1000,1000) ); "NA"/"?" (no call / uncertain) leave the model's bound alone.

The heavy lifting (offline construction, objective shim, per-stage isolation that
avoids the GLPK process-abort) is reused from ``direction_change_template_eval``'s
``OfflineTemplateEval`` -- see that module and
``reports/TEMPLATE_DIRECTION_EVAL.md`` / the KBUtilLib audit for why the shims
exist.

Output: ``site/data/template_quality_<scope>.json`` (+ a raw ``.jsonl`` for
resumability).  Run:

    PY=/mnt/homes/ctaylor/conda/miniforge3/envs/core_models_analysis/bin/python
    $PY scripts/template_quality_heuristics.py --panel            # 100-model panel
    $PY scripts/template_quality_heuristics.py --all --workers 96  # all core models
"""

from __future__ import annotations

import argparse
import csv
import datetime
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR", "/scratch/ctaylor/core_models_analysis"))
RESULTS_DIR = ANALYSIS_DIR / "results"
SITE_DATA = ANALYSIS_DIR / "site" / "data"
MODELS_DIR = Path("/scratch/ctaylor/core_models_kegg2")
LIT_TSV = RESULTS_DIR / "reaction_directions_literature_vs_llm.tsv"
PANEL_IDS = RESULTS_DIR / "selected_ids.txt"

SCHEMES = [
    {"key": "default", "label": "Default (model bounds)", "col": None},
    {"key": "jankowski", "label": "Jankowski (group contribution)", "col": "Jankowski_2008"},
    {"key": "flamholz", "label": "Flamholz 2012 (eQuilibrator)", "col": "Flamholz_2012"},
    {"key": "opus", "label": "Claude Opus 4.8", "col": "LLM_Opus_4.8"},
]
METRICS = [
    ("dead", "dead reactions"),
    ("forward_only", "forward-only"),
    ("reverse_only", "reverse-only"),
    ("reversible", "reversible"),
    ("essential", "essential (union)"),
    ("closed_mode", "closed-mode loops"),
    ("producible_complete", "producible (complete)"),
    ("producible_glucose", "producible (glucose-min)"),
    ("consumable", "consumable (complete)"),
    ("growth", "biomass flux"),
]


# ---------------------------------------------------------------------------
# Literature direction maps
# ---------------------------------------------------------------------------
def load_literature_maps():
    """Return {col: {base_rxn_id: direction}} for the source columns (NA dropped)."""
    cols = [s["col"] for s in SCHEMES if s["col"]]
    maps = {c: {} for c in cols}
    with LIT_TSV.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rid = row["rxn_id"]
            for c in cols:
                v = (row.get(c) or "").strip()
                if v and v != "NA":
                    maps[c][rid] = v
    return maps


# ---------------------------------------------------------------------------
# Per-model evaluation (runs in a worker)
# ---------------------------------------------------------------------------
_RUNNER = None
_MU = None
_LIT = None
_BOUNDS = None
_BASEIDX = None
_GROWTH = None


def _init_worker(lit_maps):
    """Build the (expensive) MSTemplateUtils runner once per worker process."""
    global _RUNNER, _MU, _LIT, _BOUNDS, _BASEIDX, _GROWTH
    import cobra
    cobra.Configuration().solver = "glpk"  # safe given per-stage isolation
    from direction_change_template_eval import (
        make_runner_class, bounds_for_direction, build_base_to_model_index, growth_of,
    )
    cls, MU = make_runner_class()
    _RUNNER = cls.build()
    _MU = MU
    _LIT = lit_maps
    _BOUNDS = bounds_for_direction
    _BASEIDX = build_base_to_model_index
    _GROWTH = growth_of


def _count(report, *path):
    node = report
    for k in path:
        node = node.get(k, {})
    return node.get("count", 0) if isinstance(node, dict) else 0


def _summarize(report, growth):
    rich = report.get("reaction_classes", {}).get("rich", {})
    return {
        "dead": _count(rich, "dead"),
        "forward_only": _count(rich, "forward_only"),
        "reverse_only": _count(rich, "reverse_only"),
        "reversible": _count(rich, "reversible"),
        "essential": _count(rich, "essential", "union"),
        "closed_mode": _count(report, "closed_mode_reactions"),
        "producible_complete": _count(report, "producible_metabolites", "complete"),
        "producible_glucose": _count(report, "producible_metabolites", "glucose_minimal"),
        "consumable": _count(report, "consumable_metabolites", "complete"),
        "growth": round(growth, 6),
    }


def _diff_counts(before, after):
    """Count changed categories / added / removed between two reports."""
    from kbutillib.ms_template_utils import _compute_diff
    delta = _compute_diff(before, after)
    n_cat = n_add = n_rem = 0
    for k, v in delta.items():
        if isinstance(v, dict) and ("added" in v or "removed" in v):
            a, r = len(v.get("added", [])), len(v.get("removed", []))
            if a or r:
                n_cat += 1
                n_add += a
                n_rem += r
    return {"categories_changed": n_cat, "added": n_add, "removed": n_rem}


def eval_one_model(model_id):
    """Evaluate one model under all 4 schemes.  Returns a record dict (or error)."""
    try:
        import cobra
        path = MODELS_DIR / f"{model_id}.json"
        model = cobra.io.load_json_model(str(path))
        base2mdl = _BASEIDX(model)
        n_rxn = sum(1 for r in model.reactions)

        out = {"model_id": model_id, "n_reactions": n_rxn, "schemes": {}}

        # (1) default
        default_report = _RUNNER._evaluate_model_quality(_MU.get(model.copy()))
        default_growth = _GROWTH(model.copy())
        out["schemes"]["default"] = {
            **_summarize(default_report, default_growth),
            "n_overridden": 0,
        }

        # (2)-(4) literature/AI schemes
        for scheme in SCHEMES[1:]:
            dmap = _LIT[scheme["col"]]
            working = model.copy()
            n_over = 0
            for base, mids in base2mdl.items():
                d = dmap.get(base)
                if not d or d == "?":
                    continue
                b = _BOUNDS(d)
                if b is None:
                    continue
                for mid in mids:
                    r = working.reactions.get_by_id(mid)
                    if (r.lower_bound, r.upper_bound) != b:
                        r.bounds = b  # atomic
                        n_over += 1
            report = _RUNNER._evaluate_model_quality(_MU.get(working))
            growth = _GROWTH(working.copy())
            rec = {**_summarize(report, growth), "n_overridden": n_over}
            rec["diff_vs_default"] = _diff_counts(default_report, report)
            rec["growth_delta"] = round(growth - default_growth, 6)
            out["schemes"][scheme["key"]] = rec
        return out
    except Exception as exc:  # keep the batch alive; record the failure
        import traceback
        return {"model_id": model_id, "error": f"{type(exc).__name__}: {exc}",
                "trace": traceback.format_exc()[-800:]}


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def list_models(scope):
    if scope == "panel":
        ids = PANEL_IDS.read_text().split()
    else:
        ids = sorted(p.stem for p in MODELS_DIR.glob("*.json"))
    return ids


def aggregate(records):
    """Build the site JSON: per-model metrics + cross-model aggregate."""
    ok = [r for r in records if "error" not in r]
    scheme_keys = [s["key"] for s in SCHEMES]
    metric_keys = [m[0] for m in METRICS]

    models = {}
    for r in ok:
        models[r["model_id"]] = {
            "n_reactions": r.get("n_reactions"),
            "schemes": r["schemes"],
        }

    # cross-model aggregate: mean of each metric per scheme; growth-change counts
    agg = {}
    for sk in scheme_keys:
        rows = [r["schemes"][sk] for r in ok if sk in r["schemes"]]
        if not rows:
            continue
        means = {mk: round(sum(x.get(mk, 0) for x in rows) / len(rows), 3) for mk in metric_keys}
        entry = {"mean": means, "n_models": len(rows)}
        if sk != "default":
            entry["mean_overridden"] = round(sum(x.get("n_overridden", 0) for x in rows) / len(rows), 2)
            entry["models_growth_gt"] = sum(1 for x in rows if x.get("growth_delta", 0) > 1e-6)
            entry["models_growth_lt"] = sum(1 for x in rows if x.get("growth_delta", 0) < -1e-6)
            entry["mean_categories_changed"] = round(
                sum(x.get("diff_vs_default", {}).get("categories_changed", 0) for x in rows) / len(rows), 2)
        agg[sk] = entry

    return {
        "generated": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "n_models": len(ok),
        "n_errors": len(records) - len(ok),
        "schemes": SCHEMES,
        "metrics": [{"key": k, "label": lbl} for k, lbl in METRICS],
        "aggregate": agg,
        "models": models,
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = p.add_mutually_exclusive_group()
    g.add_argument("--panel", action="store_true", help="100-model panel (default)")
    g.add_argument("--all", action="store_true", help="every core model (~5,683)")
    p.add_argument("--workers", type=int, default=64)
    p.add_argument("--limit", type=int, default=0, help="cap #models (debug)")
    p.add_argument("--resume", action="store_true", help="skip models already in the .jsonl")
    args = p.parse_args(argv)

    scope = "all" if args.all else "panel"
    ids = list_models(scope)
    if args.limit:
        ids = ids[: args.limit]

    SITE_DATA.mkdir(parents=True, exist_ok=True)
    jsonl_path = SITE_DATA / f"template_quality_{scope}.jsonl"
    json_path = SITE_DATA / f"template_quality_{scope}.json"

    done = set()
    if args.resume and jsonl_path.exists():
        for line in jsonl_path.read_text().splitlines():
            try:
                done.add(json.loads(line)["model_id"])
            except Exception:
                pass
    todo = [m for m in ids if m not in done]
    print(f"scope={scope} models={len(ids)} todo={len(todo)} workers={args.workers}", flush=True)

    lit = load_literature_maps()
    print("literature maps:", {c: len(m) for c, m in lit.items()}, flush=True)

    import multiprocessing as mp
    t0 = time.time()
    n_done = 0
    mode = "a" if (args.resume and jsonl_path.exists()) else "w"
    with jsonl_path.open(mode) as out, \
            mp.Pool(args.workers, initializer=_init_worker, initargs=(lit,),
                    maxtasksperchild=25) as pool:
        for rec in pool.imap_unordered(eval_one_model, todo, chunksize=1):
            out.write(json.dumps(rec) + "\n")
            out.flush()
            n_done += 1
            if n_done % 50 == 0 or n_done == len(todo):
                rate = n_done / max(time.time() - t0, 1e-9)
                eta = (len(todo) - n_done) / max(rate, 1e-9)
                print(f"  {n_done}/{len(todo)}  {rate:.1f} models/s  ETA {eta/60:.1f} min",
                      flush=True)

    # aggregate everything in the jsonl (covers resumed runs too)
    records = [json.loads(l) for l in jsonl_path.read_text().splitlines() if l.strip()]
    site = aggregate(records)
    site["scope"] = scope
    try:
        site["panel_ids"] = PANEL_IDS.read_text().split()
    except Exception:
        site["panel_ids"] = []
    json_path.write_text(json.dumps(site))
    errs = [r["model_id"] for r in records if "error" in r]
    print(f"DONE: {site['n_models']} models, {site['n_errors']} errors -> {json_path}", flush=True)
    if errs:
        print(f"  first errors: {errs[:5]}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
