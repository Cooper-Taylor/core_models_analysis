#!/usr/bin/env python3
"""Attribute reaction-direction changes to functional consequences via KBUtilLib's
``MSTemplateUtils.diff_template_evaluation``.

This script is the bridge between two bodies of work in ``/scratch/ctaylor``:

  1. **core_models_analysis** — produces per-source reaction *direction* maps
     (``results/rxn_directions_*.json``: ``{rxn_id: ">"|"<"|"="|"?"}``) from the
     ModelSEED thermodynamics cascade and alternate sources (group-contribution,
     eQuilibrator, dGPredictor).  A "direction change" is a reaction whose
     forward/reverse/reversible verdict differs between two of those sources.

  2. **KBUtilLib** (``kbutillib.ms_template_utils``) — a recent upstream addition
     (``MSTemplateUtils``) whose ``diff_template_evaluation(model, perturbations,
     mode=...)`` applies model-level edits and reports, per edit, the functional
     changes across every quality category (reaction classes, closed-mode loops,
     Biolog growth, producible/consumable metabolites, growth).

A reaction *direction* change is exactly a *bounds* change, which maps cleanly
onto a ``diff_template_evaluation`` ``modify`` perturbation:

    ">"  forward-only  -> lower_bound=0,     upper_bound=+1000
    "<"  reverse-only  -> lower_bound=-1000, upper_bound=0
    "="  reversible    -> lower_bound=-1000, upper_bound=+1000
    "?"  unknown       -> SKIPPED (no thermodynamic call == not a bound)

So: take a base model (a core_models_kegg2 cobra model), pick a *baseline* and a
*new* direction source, turn every reaction whose direction differs (and that is
actually present in the model) into a ``modify`` perturbation, and let
``diff_template_evaluation`` tell us what each flip did.

--------------------------------------------------------------------------------
ENVIRONMENT NOTES  (see reports/TEMPLATE_DIRECTION_EVAL.md for the full writeup)
--------------------------------------------------------------------------------
``MSTemplateUtils`` was written against an internal ModelSEEDpy build.  Running it
against the public stack (cobra 0.31 + modelseedpy 0.4.2 + GLPK) needs three
narrowly-scoped, fully-documented shims, all contained in ``OfflineTemplateEval``
below.  None of them change what ``diff_template_evaluation`` computes; they only
let it run here:

  (S1) Offline construction.  ``KBModelUtils.__init__`` demands a live KBase token
       and builds a KBaseAPI client.  The direction-diff path needs none of that,
       only ``self.MSModelUtil``.  We build the object via ``__new__`` and wire up
       just the modelseedpy classes + a logger (``OfflineTemplateEval.build``).

  (S2) Objective shim.  ``set_objective_from_string`` calls
       ``pkgmgr.getpkg("ObjectivePkg")``, a package that exists only in the private
       ModelSEEDpy fork (absent from PyPI 0.4.2 *and* ModelSEEDpy GitHub main).  We
       override it to set the cobra objective directly.  ``ObjConstPkg`` (used for
       the essential/fraction-of-optimum pass) *is* public, so it is left alone.

  (S3) Stage isolation.  ``run_fva``'s growth-forced pass adds an ``ObjConstPkg``
       "biomass >= fraction*optimum" constraint that is never removed.  The shipped
       ``_evaluate_model_quality`` runs all stages on one model, so
       ``find_closed_mode_reactions`` (which zeroes every exchange) inherits that
       constraint, becomes infeasible, and GLPK *aborts the process* (SIGABRT,
       not an exception).  We override ``_evaluate_model_quality`` to run each stage
       on a fresh ``model.copy()`` so no residual constraint leaks between stages.
       (This is a genuine upstream bug — see the audit; it is masked by a solver
       that returns "infeasible" gracefully, which GLPK does not.)

  Also: ``simulate_biolog`` no-ops here because the committed stash loader calls
  ``MSGrowthPhenotypes.from_dict``, absent in modelseedpy 0.4.2 — the
  ``functional_biolog_media`` section comes back empty (handled, not fatal).

Because the shipped ``growth_change`` field is only a count of essential reactions
(not biomass flux), this script *also* computes an honest biomass-flux growth
delta per perturbation and puts it in the summary CSV.

--------------------------------------------------------------------------------
USAGE
--------------------------------------------------------------------------------
    # Offline logic check (no modelseedpy needed):
    python3 scripts/direction_change_template_eval.py --self-test

    # Build the perturbation set only, no FBA (no modelseedpy needed):
    python3 scripts/direction_change_template_eval.py --dry-run \
        --model-id GCF_000005845.2 --baseline cascade_live --new group-contribution

    # Full live diff (needs modelseedpy + kbutillib.ms_template_utils):
    python3 scripts/direction_change_template_eval.py \
        --model-id GCF_000005845.2 --baseline cascade_live --new group-contribution \
        --mode batch

Run with the conda env that has cobra + modelseedpy + kbutillib, e.g.
    /mnt/homes/ctaylor/conda/miniforge3/envs/core_models_analysis/bin/python
"""

from __future__ import annotations

import argparse
import csv
import datetime
import json
import logging
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Sibling modules in scripts/ (canonical direction loader + SEED-id normalizer)
sys.path.insert(0, str(Path(__file__).resolve().parent))
from seed_annotation import normalize_seed_id  # noqa: E402

# ---------------------------------------------------------------------------
# Paths / constants (mirror direction_pipeline.py conventions)
# ---------------------------------------------------------------------------
ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR", "/scratch/ctaylor/core_models_analysis"))
RESULTS_DIR = ANALYSIS_DIR / "results"
MODELS_DIR = ANALYSIS_DIR / "data" / "core_models_kegg2"
OUT_DIR = RESULTS_DIR / "template_direction_eval"

# Direction sources -> results/rxn_directions_<file>.json
DIRECTION_SOURCES = {
    "cascade_live": "rxn_directions_cascade_live.json",
    "msdb_dev": None,          # from rev_map_dev.json
    "msdb_claude": None,       # from rev_map_claude.json
    "group-contribution": "rxn_directions_group-contribution.json",
    "equilibrator": "rxn_directions_equilibrator.json",
    "dgpredictor": "rxn_directions_dgpredictor.json",
}
_SPECIAL_SOURCE_FILES = {
    "msdb_dev": "rev_map_dev.json",
    "msdb_claude": "rev_map_claude.json",
}

VALID_DIRECTIONS = {">", "<", "=", "?"}
BASE_RXN_RE = re.compile(r"^(rxn\d+)")            # base MSDB id at the start of a cobra rxn id
COMPARTMENT_RXN_RE = re.compile(r"^rxn\d+_[a-z]+\d*$")  # e.g. rxn00549_c0 (exclude EX_/bio/DM_/SK_)

# The ModelSEED "unlimited flux" magnitude used throughout core_models_analysis
# (growth_heuristics._bounds_for_rev). "?" is intentionally NOT mapped -- see below.
DEFAULT_BOUND = 1000.0


# ---------------------------------------------------------------------------
# Direction <-> bounds
# ---------------------------------------------------------------------------
def bounds_for_direction(direction: str, mag: float = DEFAULT_BOUND) -> Optional[Tuple[float, float]]:
    """Map a ModelSEED direction flag to cobra (lower, upper) bounds.

    Returns None for "?" (unknown): "?" means "no thermodynamic call", which is
    not a bound -- forcing it to reversible would silently relax the model and
    pollute the diff, so callers should skip these.
    """
    if direction == ">":
        return 0.0, mag
    if direction == "<":
        return -mag, 0.0
    if direction == "=":
        return -mag, mag
    return None  # "?" or anything unexpected


def direction_from_bounds(lb: float, ub: float) -> str:
    """Infer a direction flag from a reaction's current cobra bounds."""
    fwd = ub > 1e-9
    rev = lb < -1e-9
    if fwd and rev:
        return "="
    if fwd:
        return ">"
    if rev:
        return "<"
    return "0"  # hard-closed (lb==ub==0) -- treat as its own state


# ---------------------------------------------------------------------------
# Direction-map loading
# ---------------------------------------------------------------------------
def _normalize_direction(value) -> str:
    v = str(value).strip()
    return v if v in VALID_DIRECTIONS else "?"


def load_direction_map(source: str) -> Dict[str, str]:
    """Load a ``{base_rxn_id: direction}`` map for a named source (see DIRECTION_SOURCES)."""
    if source in _SPECIAL_SOURCE_FILES:
        path = RESULTS_DIR / _SPECIAL_SOURCE_FILES[source]
    elif source in DIRECTION_SOURCES and DIRECTION_SOURCES[source]:
        path = RESULTS_DIR / DIRECTION_SOURCES[source]
    else:
        # Allow an explicit path too
        path = Path(source)
    if not path.exists():
        raise FileNotFoundError(f"Direction source '{source}' -> {path} not found")
    data = json.loads(path.read_text())
    if not isinstance(data, dict):
        raise ValueError(f"Expected a {{rxn_id: direction}} dict in {path}")
    return {str(k): _normalize_direction(v) for k, v in data.items()}


# ---------------------------------------------------------------------------
# base-id -> model-reaction-id index
# ---------------------------------------------------------------------------
def build_base_to_model_index(model) -> Dict[str, List[str]]:
    """Map each MSDB base id (``rxn00549``) to the model reaction ids that carry it.

    Prefers the ``annotation['seed.reaction']`` (normalized to strip the stray
    ``_c`` suffix bug), falling back to the id prefix for reactions that lack the
    annotation.  Excludes exchanges / biomass / demand / sink reactions.
    """
    index: Dict[str, List[str]] = {}
    for rxn in model.reactions:
        base = None
        anno = getattr(rxn, "annotation", None) or {}
        seed = anno.get("seed.reaction") if hasattr(anno, "get") else None
        if seed:
            base = normalize_seed_id(seed)
        elif COMPARTMENT_RXN_RE.match(rxn.id):
            m = BASE_RXN_RE.match(rxn.id)
            base = m.group(1) if m else None
        if base:
            index.setdefault(base, []).append(rxn.id)
    return index


# ---------------------------------------------------------------------------
# Perturbation construction
# ---------------------------------------------------------------------------
def build_perturbations(
    model,
    new_map: Dict[str, str],
    baseline_map: Optional[Dict[str, str]] = None,
    mag: float = DEFAULT_BOUND,
) -> Tuple[List[dict], List[dict]]:
    """Turn direction changes into ``diff_template_evaluation`` ``modify`` perturbations.

    A reaction contributes a perturbation iff:
      * its base id is present in the loaded model, AND
      * the new direction is concrete (not "?"), AND
      * the new direction differs from the baseline direction.  The baseline is
        ``baseline_map`` if given, otherwise the reaction's *current model bounds*
        (so we never emit a no-op ``modify`` that re-sets identical bounds).

    Returns ``(perturbations, records)`` where ``records`` carries human-readable
    provenance for the summary CSV.
    """
    base2mdl = build_base_to_model_index(model)
    rxn_by_id = {r.id: r for r in model.reactions}

    perturbations: List[dict] = []
    records: List[dict] = []
    for base, model_ids in sorted(base2mdl.items()):
        to_dir = new_map.get(base)
        if to_dir is None or to_dir == "?":
            continue
        new_bounds = bounds_for_direction(to_dir, mag)
        if new_bounds is None:
            continue
        for mid in model_ids:
            rxn = rxn_by_id[mid]
            cur_lb, cur_ub = rxn.lower_bound, rxn.upper_bound
            from_dir = (
                baseline_map.get(base) if baseline_map is not None
                else direction_from_bounds(cur_lb, cur_ub)
            )
            if from_dir == to_dir:
                continue  # no change
            perturbations.append({
                "op": "modify",
                "reaction_id": mid,
                "lower_bound": new_bounds[0],
                "upper_bound": new_bounds[1],
            })
            records.append({
                "base_rxn": base,
                "model_rxn": mid,
                "from_dir": from_dir,
                "to_dir": to_dir,
                "from_bounds": [cur_lb, cur_ub],
                "to_bounds": list(new_bounds),
            })
    return perturbations, records


# ---------------------------------------------------------------------------
# Offline runner (contains the S1/S2/S3 shims documented in the module header)
# ---------------------------------------------------------------------------
_OBJ_RE = re.compile(r"^\s*(MAX|MIN)\{(.+)\}\s*$", re.IGNORECASE)


def _import_kbutillib():
    """Import MSTemplateUtils + modelseedpy classes (only needed for live runs)."""
    from kbutillib.ms_template_utils import MSTemplateUtils
    from modelseedpy.core.msmodelutl import MSModelUtil
    from modelseedpy.core.msmedia import MSMedia
    from modelseedpy.core.msgrowthphenotypes import MSGrowthPhenotypes
    return MSTemplateUtils, MSModelUtil, MSMedia, MSGrowthPhenotypes


def make_runner_class():
    """Build the OfflineTemplateEval subclass (deferred so --self-test/--dry-run
    do not require modelseedpy)."""
    MSTemplateUtils, MSModelUtil, MSMedia, MSGrowthPhenotypes = _import_kbutillib()

    def _with_count(lst):
        return {"list": list(lst), "count": len(lst)}

    class OfflineTemplateEval(MSTemplateUtils):
        # ---- S2: objective shim (public modelseedpy has no ObjectivePkg) --------
        def set_objective_from_string(self, model, objective):
            if objective is None:
                return
            cm = self._check_and_convert_model(model).model
            m = _OBJ_RE.match(objective)
            if not m:
                raise ValueError(f"unsupported objective string: {objective!r}")
            cm.objective = cm.reactions.get_by_id(m.group(2).strip())
            cm.objective_direction = m.group(1).lower()

        # ---- S3: run every eval stage on a fresh copy (no residual constraints) --
        def _evaluate_model_quality(self, mdlutl, rich_media=None, minimal_media=None):
            base = mdlutl.model

            def fresh():
                return MSModelUtil.get(base.copy())

            def rc(cd):
                out = {k: _with_count(cd.get(k, [])) for k in
                       ("dead", "forward_only", "reverse_only", "reversible")}
                out["essential"] = {bk: _with_count(v) for bk, v in cd.get("essential", {}).items()}
                return out

            rich = rc(self.classify_reactions_by_fva(fresh()))
            minimal = rc(self.classify_reactions_by_fva(fresh()))
            closed = _with_count(self.find_closed_mode_reactions(fresh()))
            try:
                biolog_raw = self.simulate_biolog(fresh())
            except Exception:
                biolog_raw = {}
            biolog = {e: {bk: _with_count(v) for bk, v in bd.items()}
                      for e, bd in biolog_raw.items()}
            producible_complete = _with_count(self.test_production_potential(fresh()))
            producible_glucose = _with_count(self.test_production_potential(fresh()))
            consumable_complete = _with_count(self.test_degradation_potential(fresh()))

            bids = [r.id for r in base.reactions
                    if r.id.startswith("bio") or "biomass" in r.id.lower()]
            bio_rxns = [b for b in bids if b in ("bio1", "bio2")] or bids[:2]
            return {
                "template_metadata": {
                    "id": getattr(base, "id", "unknown"),
                    "biomass_ids": bio_rxns,
                    "rich_media": None,
                    "minimal_media": None,
                    "timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                },
                "reaction_classes": {"rich": rich, "minimal": minimal},
                "closed_mode_reactions": closed,
                "functional_biolog_media": biolog,
                "producible_metabolites": {"complete": producible_complete,
                                           "glucose_minimal": producible_glucose},
                "consumable_metabolites": {"complete": consumable_complete},
            }

        # ---- S1: offline construction ------------------------------------------
        @classmethod
        def build(cls, log_level=logging.ERROR):
            obj = cls.__new__(cls)
            obj.logger = logging.getLogger("OfflineTemplateEval")
            obj.logger.setLevel(log_level)
            obj.name = "offline"
            obj.version = "offline"
            obj.kb_version = "prod"
            obj.MSModelUtil = MSModelUtil
            obj.MSMedia = MSMedia
            obj.MSGrowthPhenotypes = MSGrowthPhenotypes
            return obj

    return OfflineTemplateEval, MSModelUtil


# ---------------------------------------------------------------------------
# Honest biomass-flux growth delta (compensates for the essential-count proxy)
# ---------------------------------------------------------------------------
def _find_biomass(model):
    for rid in ("bio1", "bio2", "biomass", "Biomass"):
        if rid in model.reactions:
            return model.reactions.get_by_id(rid)
    for r in model.reactions:
        if r.id.lower().startswith("bio") and not r.id.startswith("SK_"):
            return r
    return None


def growth_of(model) -> float:
    """Max biomass flux for a cobra model copy (0.0 if infeasible/no biomass)."""
    bio = _find_biomass(model)
    if bio is None:
        return 0.0
    with model:
        model.objective = bio
        model.objective_direction = "max"
        val = model.slim_optimize()
    return float(val) if val is not None and val == val else 0.0  # nan -> 0


# ---------------------------------------------------------------------------
# Report writing
# ---------------------------------------------------------------------------
def _nonempty_categories(delta: dict) -> Dict[str, dict]:
    return {k: v for k, v in delta.items()
            if isinstance(v, dict) and (v.get("added") or v.get("removed"))}


def write_outputs(tag: str, diff_report: dict, records: List[dict],
                  growth_deltas: List[dict], out_dir: Path) -> Dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / f"{tag}.diff.json"
    csv_path = out_dir / f"{tag}.summary.csv"
    md_path = out_dir / f"{tag}.summary.md"
    growth_path = out_dir / f"{tag}.growth.csv"

    json_path.write_text(json.dumps(diff_report, indent=2))

    # Per-reaction honest growth sensitivity (each flip applied alone to the
    # baseline model) -- written for every mode so batch runs keep it too.
    gd_lookup = {g["model_rxn"]: g for g in growth_deltas}
    with growth_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["model_rxn", "base_rxn", "from_dir", "to_dir",
                    "growth_before", "growth_after", "growth_delta"])
        for rec in records:
            g = gd_lookup.get(rec["model_rxn"], {})
            w.writerow([rec["model_rxn"], rec["base_rxn"], rec["from_dir"], rec["to_dir"],
                        g.get("before", ""), g.get("after", ""), g.get("delta", "")])

    # One CSV row per perturbation: provenance + counts of changed categories + growth
    gd_by_rxn = {g["model_rxn"]: g for g in growth_deltas}
    rec_by_rxn = {r["model_rxn"]: r for r in records}
    with csv_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["model_rxn", "base_rxn", "from_dir", "to_dir",
                    "from_bounds", "to_bounds",
                    "n_categories_changed", "n_added", "n_removed",
                    "growth_before", "growth_after", "growth_delta",
                    "changed_categories"])
        for pd_ in diff_report.get("perturbation_diffs", []):
            pert = pd_["perturbation"]
            rid = pert.get("reaction_id", "")
            delta = pd_["delta"]
            ne = _nonempty_categories(delta)
            n_added = sum(len(v["added"]) for v in ne.values())
            n_removed = sum(len(v["removed"]) for v in ne.values())
            rec = rec_by_rxn.get(rid, {})
            gd = gd_by_rxn.get(rid, {})
            w.writerow([
                rid, rec.get("base_rxn", ""), rec.get("from_dir", ""), rec.get("to_dir", ""),
                rec.get("from_bounds", ""), rec.get("to_bounds", ""),
                len(ne), n_added, n_removed,
                gd.get("before", ""), gd.get("after", ""), gd.get("delta", ""),
                ";".join(sorted(ne.keys())),
            ])

    # Markdown summary
    lines = [f"# Direction-change template evaluation: {tag}\n",
             f"- Mode: **{diff_report.get('mode')}**",
             f"- Perturbations: **{len(diff_report.get('perturbation_diffs', []))}**\n",
             "| model_rxn | base | from→to | #cats | +added | -removed | Δgrowth |",
             "|---|---|---|---|---|---|---|"]
    for pd_ in diff_report.get("perturbation_diffs", []):
        rid = pd_["perturbation"].get("reaction_id", "")
        ne = _nonempty_categories(pd_["delta"])
        na = sum(len(v["added"]) for v in ne.values())
        nr = sum(len(v["removed"]) for v in ne.values())
        rec = rec_by_rxn.get(rid, {})
        gd = gd_by_rxn.get(rid, {})
        gdelta = gd.get("delta", "")
        lines.append(f"| {rid} | {rec.get('base_rxn','')} | "
                     f"{rec.get('from_dir','')}→{rec.get('to_dir','')} | "
                     f"{len(ne)} | {na} | {nr} | {gdelta} |")
    md_path.write_text("\n".join(lines) + "\n")
    return {"json": json_path, "csv": csv_path, "md": md_path}


# ---------------------------------------------------------------------------
# Core drivers
# ---------------------------------------------------------------------------
def load_model(model_id: str):
    import cobra
    path = model_id if model_id.endswith(".json") else str(_resolve_model_path(model_id))
    return cobra.io.load_json_model(path)


def _resolve_model_path(model_id: str) -> Path:
    for cand in (MODELS_DIR / f"{model_id}.json",
                 Path("/scratch/ctaylor/core_models_kegg2") / f"{model_id}.json"):
        if cand.exists():
            return cand
    raise FileNotFoundError(f"Model {model_id} not found under {MODELS_DIR} or /scratch/ctaylor/core_models_kegg2")


def run_live(args) -> int:
    import cobra
    cobra.Configuration().solver = args.solver

    OfflineTemplateEval, MSModelUtil = make_runner_class()
    runner = OfflineTemplateEval.build(
        log_level=logging.INFO if args.verbose else logging.ERROR)

    model = load_model(args.model_id)
    baseline_map = None if args.baseline in (None, "model") else load_direction_map(args.baseline)
    new_map = load_direction_map(args.new)

    perturbations, records = build_perturbations(model, new_map, baseline_map)
    if args.limit:
        perturbations = perturbations[: args.limit]
        records = records[: args.limit]
    print(f"[{args.model_id}] {len(perturbations)} direction-change perturbations "
          f"(baseline={args.baseline or 'model-bounds'}, new={args.new})", flush=True)
    if not perturbations:
        print("Nothing to do -- no in-model direction changes for this pairing.")
        return 0

    # Honest per-perturbation growth delta (baseline vs single flip), independent of
    # the shipped essential-count proxy.
    growth_deltas: List[dict] = []
    base_growth = growth_of(model.copy())
    for rec in records:
        mcopy = model.copy()
        rxn = mcopy.reactions.get_by_id(rec["model_rxn"])
        rxn.bounds = tuple(rec["to_bounds"])  # atomic (avoids the modify bound-order footgun)
        after_growth = growth_of(mcopy)
        growth_deltas.append({
            "model_rxn": rec["model_rxn"],
            "before": round(base_growth, 6),
            "after": round(after_growth, 6),
            "delta": round(after_growth - base_growth, 6),
        })

    tag = f"{args.model_id}__{args.baseline or 'model'}__to__{args.new}__{args.mode}"

    if args.mode in ("independent", "cumulative"):
        diff_report = runner.diff_template_evaluation(
            model, perturbations, mode=args.mode, baseline_report=None)
    elif args.mode == "batch":
        # Apply every flip at once and diff the combined edit vs baseline -- 2 evals
        # total instead of one per flip.  Uses the module's own _compute_diff.
        from kbutillib.ms_template_utils import _compute_diff
        baseline_report = runner._evaluate_model_quality(MSModelUtil.get(model.copy()))
        working = model.copy()
        for rec in records:
            working.reactions.get_by_id(rec["model_rxn"]).bounds = tuple(rec["to_bounds"])
        perturbed_report = runner._evaluate_model_quality(MSModelUtil.get(working))
        combined_delta = _compute_diff(baseline_report, perturbed_report)
        diff_report = {
            "mode": "batch",
            "baseline_report": baseline_report,
            "perturbation_diffs": [{
                "perturbation": {"op": "modify_batch",
                                 "reaction_ids": [r["model_rxn"] for r in records]},
                "delta": combined_delta,
            }],
        }
    else:
        raise ValueError(f"unknown mode {args.mode}")

    paths = write_outputs(tag, diff_report, records, growth_deltas, OUT_DIR)
    print("Wrote:")
    for k, p in paths.items():
        print(f"  {k}: {p}")
    return 0


# ---------------------------------------------------------------------------
# Dry-run: build + report perturbations only (no FBA, no modelseedpy)
# ---------------------------------------------------------------------------
def run_dry(args) -> int:
    model = load_model(args.model_id)
    baseline_map = None if args.baseline in (None, "model") else load_direction_map(args.baseline)
    new_map = load_direction_map(args.new)
    perturbations, records = build_perturbations(model, new_map, baseline_map)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tag = f"{args.model_id}__{args.baseline or 'model'}__to__{args.new}"
    pert_path = OUT_DIR / f"{tag}.perturbations.json"
    pert_path.write_text(json.dumps(
        {"model_id": args.model_id, "baseline": args.baseline, "new": args.new,
         "n_perturbations": len(perturbations),
         "perturbations": perturbations, "records": records}, indent=2))

    import collections
    trans = collections.Counter(f"{r['from_dir']}->{r['to_dir']}" for r in records)
    print(f"[{args.model_id}] {len(perturbations)} perturbations "
          f"(baseline={args.baseline or 'model-bounds'}, new={args.new})")
    for k, v in trans.most_common():
        print(f"   {k}: {v}")
    print(f"Wrote {pert_path}")
    return 0


# ---------------------------------------------------------------------------
# Self-test: validate perturbation builder + module pure functions (no modelseedpy)
# ---------------------------------------------------------------------------
def run_self_test(args) -> int:
    import cobra
    from kbutillib.ms_template_utils import _apply_perturbation, _compute_diff, _render_markdown

    ok = True

    def check(name, cond):
        nonlocal ok
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
        ok = ok and cond

    print("== bounds_for_direction ==")
    check(">  -> (0,1000)", bounds_for_direction(">") == (0.0, 1000.0))
    check("<  -> (-1000,0)", bounds_for_direction("<") == (-1000.0, 0.0))
    check("=  -> (-1000,1000)", bounds_for_direction("=") == (-1000.0, 1000.0))
    check("?  -> None (skipped)", bounds_for_direction("?") is None)

    print("== build a synthetic model + perturbation builder ==")
    m = cobra.Model("toy")
    a = cobra.Metabolite("cpdA_c0", compartment="c0")
    b = cobra.Metabolite("cpdB_c0", compartment="c0")
    r1 = cobra.Reaction("rxn00001_c0"); r1.add_metabolites({a: -1, b: 1})
    r1.bounds = (-1000, 1000)  # reversible
    r1.annotation = {"seed.reaction": "rxn00001"}
    r2 = cobra.Reaction("rxn00002_c0"); r2.add_metabolites({b: -1})
    r2.bounds = (0, 1000)      # forward
    r2.annotation = {"seed.reaction": "rxn00002"}
    m.add_reactions([r1, r2])

    idx = build_base_to_model_index(m)
    check("base index resolves rxn00001 -> rxn00001_c0", idx.get("rxn00001") == ["rxn00001_c0"])

    new_map = {"rxn00001": ">", "rxn00002": ">", "rxn99999": "<"}
    perts, recs = build_perturbations(m, new_map, baseline_map=None)
    check("only the changed, in-model, non-? reaction becomes a perturbation",
          [p["reaction_id"] for p in perts] == ["rxn00001_c0"])
    check("perturbation encodes forward bounds",
          perts[0]["lower_bound"] == 0.0 and perts[0]["upper_bound"] == 1000.0 and perts[0]["op"] == "modify")
    check("'?' new direction is skipped",
          all(r["to_dir"] != "?" for r in recs))

    print("== module pure functions round-trip ==")
    _apply_perturbation(m, perts[0])
    check("_apply_perturbation set the new bounds",
          m.reactions.get_by_id("rxn00001_c0").bounds == (0.0, 1000.0))

    before = {"closed_mode_reactions": {"list": ["x", "y"], "count": 2},
              "reaction_classes": {"rich": {}, "minimal": {}},
              "producible_metabolites": {}, "consumable_metabolites": {},
              "functional_biolog_media": {}}
    after = {"closed_mode_reactions": {"list": ["y", "z"], "count": 2},
             "reaction_classes": {"rich": {}, "minimal": {}},
             "producible_metabolites": {}, "consumable_metabolites": {},
             "functional_biolog_media": {}}
    delta = _compute_diff(before, after)
    check("_compute_diff added={z}", delta["closed_mode_reactions"]["added"] == ["z"])
    check("_compute_diff removed={x}", delta["closed_mode_reactions"]["removed"] == ["x"])

    md = _render_markdown({"template_metadata": {"id": "toy"}})
    check("_render_markdown returns markdown", isinstance(md, str) and md.startswith("#"))

    print("\nSELF-TEST:", "OK" if ok else "FAILURES")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--model-id", default="GCF_000005845.2",
                   help="core model id (default GCF_000005845.2, E. coli) or path to a .json")
    p.add_argument("--baseline", default="cascade_live",
                   help="baseline direction source (name in DIRECTION_SOURCES, a path, "
                        "or 'model' to use the model's current bounds). Default cascade_live.")
    p.add_argument("--new", default="group-contribution",
                   help="new/target direction source. Default group-contribution.")
    p.add_argument("--mode", default="batch",
                   choices=["independent", "cumulative", "batch"],
                   help="independent/cumulative -> diff_template_evaluation; "
                        "batch -> one combined edit (fast). Default batch.")
    p.add_argument("--limit", type=int, default=0,
                   help="cap number of perturbations (0 = all). Independent mode is ~5s/pert.")
    p.add_argument("--solver", default="glpk",
                   help="cobra LP solver (default glpk; works given stage isolation)")
    p.add_argument("--verbose", action="store_true")
    p.add_argument("--dry-run", action="store_true",
                   help="build + write perturbations only, no FBA (no modelseedpy needed)")
    p.add_argument("--self-test", action="store_true",
                   help="validate logic on a synthetic model (no modelseedpy needed)")
    args = p.parse_args(argv)

    if args.self_test:
        return run_self_test(args)
    if args.dry_run:
        return run_dry(args)
    return run_live(args)


if __name__ == "__main__":
    raise SystemExit(main())
