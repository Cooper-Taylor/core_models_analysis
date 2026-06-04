"""
Growth pipeline that re-bounds each panel model's reactions using a chosen
reversibility map, then runs FBA on the biomass reaction.

This complements ``analyze_growth.py``.  That script reads the on-disk
``lower_bound`` / ``upper_bound`` of every reaction as-is.  This module lets
the notebook overlay a fresh ``{rxn_id: reversibility}`` map (e.g. the output
of ``reversibility_lib.run_cascade(cfg=...)``) so we can quantify how each
heuristic change moves model growth.

We never write back to ``core_models_kegg2/*.json`` -- all overrides live in
memory.
"""

from __future__ import annotations

import logging
import multiprocessing as mp
from pathlib import Path
from typing import Optional

import cobra
from cobra.io import load_json_model

ANALYSIS_DIR = Path("/scratch/ctaylor/core_models_analysis")
MODELS_DIR = ANALYSIS_DIR / "data" / "core_models_kegg2"
MEDIA_FILE = Path("/scratch/ctaylor/ModelSEEDDatabase/Media/KBaseMedia.cpd")
GROWTH_THRESHOLD = 1e-6


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
    """Mirror of ``analyze_growth.apply_media`` -- restrict uptake to media."""
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
            rxn.lower_bound = -1000.0
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
        seed = rxn.annotation.get("seed.reaction") if rxn.annotation else None
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
    """Apply media, optionally rebound, run FBA on biomass.

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
        bio_rxn = find_biomass_reaction(model)
        if bio_rxn is None:
            res["status"] = "no_biomass"
            return res
        res["biomass_rxn"] = bio_rxn.id
        model.objective = bio_rxn
        sol = model.optimize()
        res["status"] = sol.status
        flux = float(sol.objective_value) if sol.objective_value is not None else 0.0
        res["growth_flux"] = flux
        res["grows"] = (sol.status == "optimal") and (flux > GROWTH_THRESHOLD)
    except Exception as e:
        res["error"] = f"{type(e).__name__}: {e}"
    return res


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
