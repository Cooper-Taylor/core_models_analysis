"""
Growth pipeline that re-bounds each panel model's reactions using a chosen
reversibility map, then runs FBA maximizing ATP production.

This complements ``analyze_growth.py``.  That script reads the on-disk
``lower_bound`` / ``upper_bound`` of every reaction as-is.  This module lets
the notebook overlay a fresh ``{rxn_id: reversibility}`` map (e.g. the output
of ``reversibility_lib.run_cascade(cfg=...)``) so we can quantify how each
heuristic change moves model ATP production.

We never write back to ``core_models_kegg2/*.json`` -- all overrides live in
memory.
"""

from __future__ import annotations
import os

import logging
import multiprocessing as mp
from pathlib import Path
from typing import Optional

import cobra
from cobra.io import load_json_model

from seed_annotation import seed_id  # normalizes the _c-suffix bug

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR", "/scratch/ctaylor/core_models_analysis"))
MODELS_DIR = ANALYSIS_DIR / "data" / "core_models_kegg2"
MEDIA_FILE = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase") + "/Media/KBaseMedia.cpd")
GROWTH_THRESHOLD = 1e-6
# Uniform uptake cap (mmol/gDW/h) on every open media exchange. The objective is
# ATP production (see ``ensure_atp_objective``); without a finite uptake cap the
# ATP-hydrolysis demand saturates at its own upper bound for every model, so the
# metric is reported per unit substrate uptake ("ATP per unit total uptake").
UPTAKE_BOUND = 1.0


def load_media_compounds(path: Path = MEDIA_FILE) -> set:
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}


# Lazy module-level singleton -- the multiprocessing workers re-import this
# module per child process and load the media list once.
_MEDIA_CPDS: Optional[set] = None


def _media_cpds() -> set:
    global _MEDIA_CPDS
    if _MEDIA_CPDS is None:
        _MEDIA_CPDS = load_media_compounds()
    return _MEDIA_CPDS


def apply_media(model: cobra.Model, media_cpds: Optional[set] = None) -> int:
    """Restrict uptake to the media compounds, each capped at ``-UPTAKE_BOUND``.

    Open media exchanges get a finite uptake lower bound (``-UPTAKE_BOUND``)
    rather than the usual ``-1000`` so the ATP-production objective yields a
    finite, comparable number per unit substrate uptake (see ``UPTAKE_BOUND``);
    secretion (upper bound) stays open at ``1000``."""
    if media_cpds is None:
        media_cpds = _media_cpds()
    open_count = 0
    for rxn in model.reactions:
        if not rxn.id.startswith("EX_"):
            continue
        mets = list(rxn.metabolites.keys())
        if len(mets) != 1:
            continue
        cpd_id = mets[0].id.split("_")[0]
        if cpd_id in media_cpds:
            rxn.lower_bound = -UPTAKE_BOUND
            open_count += 1
        else:
            rxn.lower_bound = 0.0
        if rxn.upper_bound < 1000.0:
            rxn.upper_bound = 1000.0
    return open_count


def find_biomass_reaction(model: cobra.Model):
    """Same precedence as ``analyze_growth.find_biomass_reaction``."""
    if "bio1" in model.reactions:
        return model.reactions.get_by_id("bio1")
    for rid in ("bio2", "biomass", "Biomass"):
        if rid in model.reactions:
            return model.reactions.get_by_id(rid)
    for r in model.reactions:
        if r.id.lower().startswith("bio") and not r.id.startswith("SK_"):
            return r
    return None


# ATP-production objective. These core models are scored on their capacity to
# PRODUCE ATP (not biomass): we add an ATP-hydrolysis demand reaction and
# maximize it (the standard MSATPCorrection-style "max ATP production" test):
#   ATP[c] + H2O[c] -> ADP[c] + Pi[c] + H+[c]
# The max flux is the network's ATP-regeneration rate on the applied media.
ATP_DEMAND_RXN = "DM_atp_c0"
_ATP_STOICH = {
    "cpd00002_c0": -1.0,  # ATP
    "cpd00001_c0": -1.0,  # H2O
    "cpd00008_c0": 1.0,   # ADP
    "cpd00009_c0": 1.0,   # Pi
    "cpd00067_c0": 1.0,   # H+
}


def ensure_atp_objective(model: cobra.Model):
    """Add (once) the ATP-hydrolysis demand reaction ``DM_atp_c0`` and return it.

    The ``DM_`` prefix keeps it out of the ``ignore_bounds`` rebinder and
    ``override_bounds`` (no seed id), so reversibility variants still apply only
    to internal reactions. Returns ``None`` if the model lacks any of the five
    ATP-test metabolites (then it can't be scored for ATP production)."""
    if ATP_DEMAND_RXN in model.reactions:
        return model.reactions.get_by_id(ATP_DEMAND_RXN)
    mets = {}
    for cid, coef in _ATP_STOICH.items():
        if cid not in model.metabolites:
            return None
        mets[model.metabolites.get_by_id(cid)] = coef
    rxn = cobra.Reaction(ATP_DEMAND_RXN, name="ATP production (hydrolysis demand)",
                         lower_bound=0.0, upper_bound=1000.0)
    rxn.add_metabolites(mets)
    model.add_reactions([rxn])
    return rxn


def _bounds_for_rev(rev: str, default_bound: float = 1000.0):
    """Map a ModelSEED reversibility flag to cobra bounds.

    ``?`` is treated as reversible (the conservative choice -- matches the
    way ModelSEED handles unknown direction when building a model).
    """
    if rev == ">":
        return 0.0, default_bound
    if rev == "<":
        return -default_bound, 0.0
    return -default_bound, default_bound


def override_bounds(model: cobra.Model, reversibility_map: dict,
                    only_changed_vs_msdb: Optional[dict] = None) -> dict:
    """Rewrite every model reaction's bounds from ``reversibility_map``.

    Only reactions whose ``annotation['seed.reaction']`` key is in the map
    get touched -- exchange reactions, sinks, and the biomass reaction are
    left alone so the media-application step stays valid.

    SEED ids are normalized via ``seed_annotation.seed_id`` so that the
    17 cobra reactions whose annotation carries a stray ``_c`` suffix
    (e.g. ``rxn11322_c``) are correctly matched against the cascade's
    bare ``rxn11322`` key. Without the normalization, those overrides
    were silently skipped.

    If ``only_changed_vs_msdb`` is supplied (the *baseline* reversibility map
    that the on-disk model bounds already reflect), only reactions whose new
    flag differs from the baseline are rewritten.  This minimizes the FBA
    perturbation so the diff isolates the heuristic change.

    Returns a dict of stats: ``{'touched': N, 'unchanged': M, 'no_anno': K}``.
    """
    touched = 0
    unchanged = 0
    no_anno = 0
    for rxn in model.reactions:
        seed = seed_id(rxn)
        if not seed:
            no_anno += 1
            continue
        new_rev = reversibility_map.get(seed)
        if new_rev is None:
            continue
        if only_changed_vs_msdb is not None:
            base_rev = only_changed_vs_msdb.get(seed)
            if base_rev == new_rev:
                unchanged += 1
                continue
        lb, ub = _bounds_for_rev(new_rev)
        rxn.lower_bound = lb
        rxn.upper_bound = ub
        touched += 1
    return {"touched": touched, "unchanged": unchanged, "no_anno": no_anno}


def fba_one(model_id: str, reversibility_map: Optional[dict] = None,
            baseline_map: Optional[dict] = None,
            ignore_bounds: bool = False) -> dict:
    """Apply media, optionally rebound, run FBA maximizing ATP production.

    - ``reversibility_map`` is ``None``  -- keep the on-disk bounds.
    - ``reversibility_map`` is a dict   -- rewrite bounds before solving.
    - ``baseline_map`` lets you rewrite only where the new map *differs* from
      the baseline (= what the on-disk model bounds already encode).
    - ``ignore_bounds`` strips every non-exchange reaction back to
      ``(-1000, 1000)`` first.  Used by the "all-reversible" control.
    """
    model_path = MODELS_DIR / f"{model_id}.json"
    res = {
        "model_id": model_id,
        "status": "",
        "growth_flux": 0.0,
        "grows": False,
        "biomass_rxn": "",
        "n_overrides": 0,
        "n_unchanged_vs_baseline": 0,
        "error": "",
    }
    try:
        model = load_json_model(str(model_path))
        if ignore_bounds:
            for r in model.reactions:
                if not r.id.startswith(("EX_", "SK_", "DM_", "bio")):
                    r.lower_bound = -1000.0
                    r.upper_bound = 1000.0
        if reversibility_map is not None:
            stats = override_bounds(model, reversibility_map,
                                    only_changed_vs_msdb=baseline_map)
            res["n_overrides"] = stats["touched"]
            res["n_unchanged_vs_baseline"] = stats["unchanged"]

        apply_media(model)
        atp_rxn = ensure_atp_objective(model)
        if atp_rxn is None:
            res["status"] = "no_atp_metabolites"
            return res
        res["biomass_rxn"] = atp_rxn.id  # objective reaction (key kept for compatibility)
        model.objective = atp_rxn
        sol = model.optimize()
        res["status"] = sol.status
        flux = float(sol.objective_value) if sol.objective_value is not None else 0.0
        res["growth_flux"] = flux  # now ATP-production flux
        res["grows"] = (sol.status == "optimal") and (flux > GROWTH_THRESHOLD)  # produces ATP
    except Exception as e:
        res["error"] = f"{type(e).__name__}: {e}"
    return res


def rxn_pipeline_one(model_id: str, baseline_map: dict, changed_dirs: dict) -> dict:
    """Decompose a variant's ATP-flux effect on one model into per-reaction parts.

    Rebinds the model to ``baseline_map`` (the heuristic-baseline cascade), then,
    for each reaction in ``changed_dirs`` (``{seed_rxn_id: variant_direction}``)
    that is present in the model:

      * **single** -- flip only that reaction to its variant direction, solve, and
        record ``delta = ATP_flux - base_flux`` (then restore the baseline bounds);
      * **cumulative** -- after sorting the single-reaction deltas by ``|delta|``
        descending, re-apply them in that order *without* reverting, solving after
        each, to show how the direction changes compound.

    The model JSON is loaded once; all perturbations happen in memory. Returns
    ``{model_id, base_flux, singles:[{rxn,delta}...], cumulative:[{rxn,delta}...]}``
    with ``singles`` and ``cumulative`` in the same (|delta|-descending) reaction
    order, or ``{model_id, error}`` on failure. ``cumulative[-1]`` (all changes
    applied) equals the variant's whole-model delta -- a built-in cross-check."""
    try:
        model = load_json_model(str(MODELS_DIR / f"{model_id}.json"))
        override_bounds(model, baseline_map)
        apply_media(model)
        atp_rxn = ensure_atp_objective(model)
        if atp_rxn is None:
            return {"model_id": model_id, "error": "no_atp_metabolites"}
        model.objective = atp_rxn

        def _flux():
            sol = model.optimize()
            if sol.status == "optimal" and sol.objective_value is not None:
                return float(sol.objective_value)
            return 0.0

        base_flux = _flux()

        # seed reaction id -> the model reaction(s) carrying that SEED annotation
        seed_rxns: dict = {}
        for rxn in model.reactions:
            s = seed_id(rxn)
            if s in changed_dirs:
                seed_rxns.setdefault(s, []).append(rxn)

        singles = []
        for seed, new_dir in changed_dirs.items():
            rs = seed_rxns.get(seed)
            if not rs:
                continue  # reaction not actually present in this model
            backup = [(r, r.lower_bound, r.upper_bound) for r in rs]
            lb, ub = _bounds_for_rev(new_dir)
            for r in rs:
                r.lower_bound, r.upper_bound = lb, ub
            singles.append({"rxn": seed, "delta": _flux() - base_flux})
            for r, lo, hi in backup:
                r.lower_bound, r.upper_bound = lo, hi

        singles.sort(key=lambda d: (-abs(d["delta"]), d["rxn"]))

        cumulative = []
        for d in singles:
            seed = d["rxn"]
            lb, ub = _bounds_for_rev(changed_dirs[seed])
            for r in seed_rxns.get(seed, []):
                r.lower_bound, r.upper_bound = lb, ub
            cumulative.append({"rxn": seed, "delta": _flux() - base_flux})

        return {"model_id": model_id, "base_flux": base_flux,
                "singles": singles, "cumulative": cumulative}
    except Exception as e:
        return {"model_id": model_id, "error": f"{type(e).__name__}: {e}"}


_worker_kwargs_cache: dict = {}


def _init_worker(reversibility_map, baseline_map, ignore_bounds):
    _worker_kwargs_cache["rmap"] = reversibility_map
    _worker_kwargs_cache["bmap"] = baseline_map
    _worker_kwargs_cache["ignore"] = ignore_bounds
    logging.getLogger("cobra").setLevel(logging.ERROR)


def _worker_run(model_id):
    return fba_one(
        model_id,
        reversibility_map=_worker_kwargs_cache.get("rmap"),
        baseline_map=_worker_kwargs_cache.get("bmap"),
        ignore_bounds=_worker_kwargs_cache.get("ignore", False),
    )


def run_panel(model_ids, reversibility_map: Optional[dict] = None,
              baseline_map: Optional[dict] = None,
              ignore_bounds: bool = False,
              n_workers: Optional[int] = None) -> list:
    """Run :func:`fba_one` across a list of model IDs in parallel.

    Pickling the reversibility maps once via ``initializer`` keeps each
    cross-process send small even when the map has 56K entries.
    """
    logging.getLogger("cobra").setLevel(logging.ERROR)
    n_workers = n_workers or max(1, min(mp.cpu_count() - 1, 16))
    if n_workers == 1:
        _init_worker(reversibility_map, baseline_map, ignore_bounds)
        return [_worker_run(mid) for mid in model_ids]
    ctx = mp.get_context("spawn") if mp.get_start_method(allow_none=True) == "spawn" \
        else mp.get_context("fork")
    with ctx.Pool(
        n_workers,
        initializer=_init_worker,
        initargs=(reversibility_map, baseline_map, ignore_bounds),
    ) as pool:
        return list(pool.imap_unordered(_worker_run, model_ids, chunksize=4))


def diff_panel(baseline_results: list, variant_results: list,
               growth_eps: float = 1e-6) -> dict:
    """Compare two ``run_panel`` outputs.

    Returns counts + the changed rows for use in the notebook tables.
    """
    base_by_id = {r["model_id"]: r for r in baseline_results}
    rows = []
    grow_change = 0
    flux_changes = 0
    for v in variant_results:
        b = base_by_id.get(v["model_id"])
        if b is None:
            continue
        d = {
            "model_id": v["model_id"],
            "baseline_grows": b["grows"],
            "variant_grows": v["grows"],
            "baseline_flux": b["growth_flux"],
            "variant_flux": v["growth_flux"],
            "delta_flux": v["growth_flux"] - b["growth_flux"],
            "n_overrides": v.get("n_overrides", 0),
        }
        rows.append(d)
        if b["grows"] != v["grows"]:
            grow_change += 1
        if abs(d["delta_flux"]) > growth_eps:
            flux_changes += 1
    return {
        "rows": rows,
        "n_models": len(rows),
        "n_grow_change": grow_change,
        "n_flux_change": flux_changes,
    }


def reversibility_diff(baseline_map: dict, variant_map: dict) -> dict:
    """Return per-reaction reversibility deltas between two maps."""
    diffs = []
    seen = set(baseline_map) | set(variant_map)
    for rxn in sorted(seen):
        b = baseline_map.get(rxn)
        v = variant_map.get(rxn)
        if b != v:
            diffs.append({"rxn": rxn, "baseline": b, "variant": v})
    counts = {}
    for d in diffs:
        key = (d["baseline"], d["variant"])
        counts[key] = counts.get(key, 0) + 1
    return {"n_changed": len(diffs), "by_transition": counts, "diffs": diffs}
