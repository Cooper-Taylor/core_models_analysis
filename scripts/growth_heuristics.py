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


def rxn_pipeline_one(model_id: str, baseline_map: dict, changed_dirs: dict) -> dict:
    """Decompose a variant's growth-flux effect on one model into per-reaction parts.

    Rebinds the model to ``baseline_map`` (the heuristic-baseline cascade), then,
    for each reaction in ``changed_dirs`` (``{seed_rxn_id: variant_direction}``)
    that is present in the model:

      * **single** -- flip only that reaction to its variant direction, solve, and
        record ``delta = growth_flux - base_flux`` (then restore the baseline bounds);
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
        bio_rxn = find_biomass_reaction(model)
        if bio_rxn is None:
            return {"model_id": model_id, "error": "no_biomass"}
        model.objective = bio_rxn

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


_ALL_DIRS = (">", "<", "=")


def key_reactions_one(model_id: str, baseline_map: dict, top_n: int = 75) -> dict:
    """Rank a model's reactions by how much *changing their direction* disrupts growth.

    Variant-agnostic single-reaction sensitivity probe. Rebinds the model to the
    cascade baseline (``baseline_map``) and records baseline growth, then for each
    SEED reaction in the model flips it to **every other** direction in
    ``('>', '<', '=')`` one at a time (keeping all others at baseline), re-solving
    growth each time. A reaction is "key" when some direction change causes a large
    |Δ growth|. The model JSON is loaded once; all flips happen in memory, so each
    reaction costs ~2 LP solves (typically a few hundred solves per model).

    Returns ``{model_id, base_flux, n_tested, reactions:[{rxn, base_dir, best_dir,
    best_delta, by_dir:{dir:Δ}, severity}...]}`` sorted by ``severity`` (max |Δ|)
    descending and capped at ``top_n``; or ``{model_id, error}`` on failure."""
    try:
        model = load_json_model(str(MODELS_DIR / f"{model_id}.json"))
        override_bounds(model, baseline_map)
        apply_media(model)
        bio_rxn = find_biomass_reaction(model)
        if bio_rxn is None:
            return {"model_id": model_id, "error": "no_biomass"}
        model.objective = bio_rxn

        def _flux():
            sol = model.optimize()
            if sol.status == "optimal" and sol.objective_value is not None:
                return float(sol.objective_value)
            return 0.0

        base_flux = _flux()

        seed_rxns: dict = {}
        for rxn in model.reactions:
            s = seed_id(rxn)
            if s:
                seed_rxns.setdefault(s, []).append(rxn)

        results = []
        for seed, rs in seed_rxns.items():
            base_dir = baseline_map.get(seed)
            if base_dir is None:
                continue  # reaction outside the cascade map: no defined baseline dir
            base_bounds = _bounds_for_rev(base_dir)
            backup = [(r, r.lower_bound, r.upper_bound) for r in rs]
            by_dir = {}
            for d in _ALL_DIRS:
                if _bounds_for_rev(d) == base_bounds:
                    continue  # equivalent to baseline (e.g. '=' vs '?') -> no change
                lb, ub = _bounds_for_rev(d)
                for r in rs:
                    r.lower_bound, r.upper_bound = lb, ub
                by_dir[d] = _flux() - base_flux
                for r, lo, hi in backup:
                    r.lower_bound, r.upper_bound = lo, hi
            if not by_dir:
                continue
            best_dir = max(by_dir, key=lambda k: abs(by_dir[k]))
            results.append({
                "rxn": seed, "base_dir": base_dir, "best_dir": best_dir,
                "best_delta": by_dir[best_dir], "by_dir": by_dir,
                "severity": abs(by_dir[best_dir]),
            })

        n_tested = len(results)
        results.sort(key=lambda x: (-x["severity"], x["rxn"]))
        keep = [r for r in results if r["severity"] > GROWTH_THRESHOLD][:top_n]
        return {"model_id": model_id, "base_flux": base_flux,
                "n_tested": n_tested, "reactions": keep}
    except Exception as e:
        return {"model_id": model_id, "error": f"{type(e).__name__}: {e}"}


def growth_control_one(model_id: str, baseline_map: dict, top_n: int = 90) -> dict:
    """Per-reaction *growth control*: single-reaction knockout Δgrowth + flux
    carried at the growth optimum + LP reduced cost.

    The classic FBA reaction-deletion analysis (Edwards & Palsson 2000): rebind
    the model to the cascade baseline + media, solve biomass once (capturing the
    optimal flux distribution and reduced costs), then BLOCK each SEED reaction
    (lb=ub=0) in turn and re-solve. A reaction is **essential** when blocking
    collapses growth (ko_delta < 0) and **growth-limiting** when blocking raises
    growth (ko_delta > 0, the reaction was siphoning flux away from biomass);
    ``flux_opt`` is the net flux it carries at the optimum and ``reduced_cost``
    its LP marginal. Model loaded once; ~one extra LP per reaction.

    Returns ``{model_id, base_flux, n_tested, n_essential, n_limiting,
    reactions:[{rxn, ko_delta, flux_opt, reduced_cost, kind}...]}`` sorted by
    ``|ko_delta|`` then ``|flux_opt|`` (capped at ``top_n``), or ``{model_id, error}``."""
    try:
        model = load_json_model(str(MODELS_DIR / f"{model_id}.json"))
        override_bounds(model, baseline_map)
        apply_media(model)
        bio_rxn = find_biomass_reaction(model)
        if bio_rxn is None:
            return {"model_id": model_id, "error": "no_biomass"}
        model.objective = bio_rxn

        sol = model.optimize()
        ok = sol.status == "optimal" and sol.objective_value is not None
        base_flux = float(sol.objective_value) if ok else 0.0
        fluxes = sol.fluxes if ok else None
        rcosts = getattr(sol, "reduced_costs", None) if ok else None

        def _flux():
            so = model.optimize()
            if so.status == "optimal" and so.objective_value is not None:
                return float(so.objective_value)
            return 0.0

        seed_rxns: dict = {}
        for rxn in model.reactions:
            s = seed_id(rxn)
            if s:
                seed_rxns.setdefault(s, []).append(rxn)

        results = []
        for seed, rs in seed_rxns.items():
            backup = [(r, r.lower_bound, r.upper_bound) for r in rs]
            for r in rs:
                r.lower_bound, r.upper_bound = 0.0, 0.0
            ko_delta = _flux() - base_flux
            for r, lo, hi in backup:
                r.lower_bound, r.upper_bound = lo, hi
            flux_opt, rc = 0.0, 0.0
            if fluxes is not None:
                for r in rs:
                    try:
                        flux_opt += float(fluxes[r.id])
                    except Exception:
                        pass
                    if rcosts is not None:
                        try:
                            v = float(rcosts[r.id])
                            if abs(v) > abs(rc):
                                rc = v
                        except Exception:
                            pass
            kind = ("essential" if ko_delta < -GROWTH_THRESHOLD
                    else "limiting" if ko_delta > GROWTH_THRESHOLD else "neutral")
            results.append({"rxn": seed, "ko_delta": ko_delta, "flux_opt": flux_opt,
                            "reduced_cost": rc, "kind": kind})

        n_tested = len(results)
        n_ess = sum(1 for d in results if d["kind"] == "essential")
        n_lim = sum(1 for d in results if d["kind"] == "limiting")
        results.sort(key=lambda d: (-abs(d["ko_delta"]), -abs(d["flux_opt"]), d["rxn"]))
        keep = [d for d in results
                if abs(d["ko_delta"]) > GROWTH_THRESHOLD or abs(d["flux_opt"]) > 1e-6][:top_n]

        # Limiting metabolites: LP shadow prices (duals) at the growth optimum --
        # the marginal change in growth per unit change in a metabolite's mass
        # balance. Large |shadow price| = the metabolite pool constrains growth.
        mets = []
        sp = getattr(sol, "shadow_prices", None) if ok else None
        if sp is not None:
            name_of = {m.id: (m.name or m.id) for m in model.metabolites}
            for mid_, val in sp.items():
                try:
                    v = float(val)
                except Exception:
                    continue
                if abs(v) > 1e-9:
                    mets.append({"met": mid_, "name": name_of.get(mid_, mid_), "shadow_price": v})
            mets.sort(key=lambda d: -abs(d["shadow_price"]))
            mets = mets[:25]

        return {"model_id": model_id, "base_flux": base_flux, "n_tested": n_tested,
                "n_essential": n_ess, "n_limiting": n_lim, "reactions": keep,
                "metabolites": mets}
    except Exception as e:
        return {"model_id": model_id, "error": f"{type(e).__name__}: {e}"}


def synthetic_lethal_one(model_id: str, baseline_map: dict,
                         n_cand: int = 35, top_n: int = 40) -> dict:
    """Find synthetic-lethal / synthetic-sick reaction *pairs* in one model.

    Single-reaction essentiality misses reactions that are dispensable alone but
    jointly essential. Rebind to the cascade baseline + media, solve biomass, and
    among the reactions that carry flux at the optimum keep those that are
    **individually non-essential** (single-knockout barely moves growth) as
    candidates (top ``n_cand`` by |flux|). Then knock out every candidate *pair*
    and record ``joint_delta`` (Δgrowth of the double knockout) and
    ``epistasis = joint_delta - (single_a + single_b)`` (synergy beyond the two
    singles). Pairs with a strongly negative joint_delta are synthetic-lethal/sick.

    Model loaded once; ~one LP per single + one per pair. Returns
    ``{model_id, base_flux, n_candidates, n_pairs, pairs:[{a,b,joint_delta,
    epistasis,single_a,single_b}...]}`` sorted by most-negative joint_delta
    (capped at ``top_n``), or ``{model_id, error}``."""
    import itertools
    try:
        model = load_json_model(str(MODELS_DIR / f"{model_id}.json"))
        override_bounds(model, baseline_map)
        apply_media(model)
        bio_rxn = find_biomass_reaction(model)
        if bio_rxn is None:
            return {"model_id": model_id, "error": "no_biomass"}
        model.objective = bio_rxn
        sol = model.optimize()
        ok = sol.status == "optimal" and sol.objective_value is not None
        base_flux = float(sol.objective_value) if ok else 0.0
        if not ok or base_flux <= GROWTH_THRESHOLD:
            return {"model_id": model_id, "base_flux": base_flux,
                    "n_candidates": 0, "n_pairs": 0, "pairs": []}
        fluxes = sol.fluxes

        def _flux():
            so = model.optimize()
            if so.status == "optimal" and so.objective_value is not None:
                return float(so.objective_value)
            return 0.0

        seed_rxns: dict = {}
        for rxn in model.reactions:
            s = seed_id(rxn)
            if s:
                seed_rxns.setdefault(s, []).append(rxn)

        # flux-carrying seeds, with single-KO delta
        carriers = []
        for seed, rs in seed_rxns.items():
            fx = sum(abs(float(fluxes[r.id])) for r in rs if r.id in fluxes.index) \
                if hasattr(fluxes, "index") else 0.0
            if fx <= 1e-6:
                continue
            carriers.append((seed, fx, rs))
        carriers.sort(key=lambda x: -x[1])
        single = {}
        cand = []
        for seed, fx, rs in carriers:
            backup = [(r, r.lower_bound, r.upper_bound) for r in rs]
            for r in rs:
                r.lower_bound, r.upper_bound = 0.0, 0.0
            d = _flux() - base_flux
            for r, lo, hi in backup:
                r.lower_bound, r.upper_bound = lo, hi
            single[seed] = d
            if d > -GROWTH_THRESHOLD:  # individually non-essential
                cand.append((seed, rs))
            if len(cand) >= n_cand:
                break

        pairs = []
        for (sa, ra), (sb, rb) in itertools.combinations(cand, 2):
            rxs = ra + rb
            backup = [(r, r.lower_bound, r.upper_bound) for r in rxs]
            for r in rxs:
                r.lower_bound, r.upper_bound = 0.0, 0.0
            joint = _flux() - base_flux
            for r, lo, hi in backup:
                r.lower_bound, r.upper_bound = lo, hi
            if joint < -GROWTH_THRESHOLD:  # the pair reduces growth (singles ~did not)
                pairs.append({"a": sa, "b": sb, "joint_delta": joint,
                              "epistasis": joint - (single[sa] + single[sb]),
                              "single_a": single[sa], "single_b": single[sb]})
        pairs.sort(key=lambda d: (d["joint_delta"], d["epistasis"]))
        return {"model_id": model_id, "base_flux": base_flux,
                "n_candidates": len(cand), "n_pairs": len(pairs),
                "pairs": pairs[:top_n]}
    except Exception as e:
        return {"model_id": model_id, "error": f"{type(e).__name__}: {e}"}


def fva_one(model_id: str, baseline_map: dict, fraction: float = 0.99,
            top_n: int = 60) -> dict:
    """Flux variability analysis of one model at near-optimal growth.

    Rebind to the cascade baseline + media, then compute the min/max flux each
    reaction can carry while keeping growth at >= ``fraction`` of optimum
    (Mahadevan & Schilling 2003). Per SEED reaction classify:
      * **blocked**     -- |min|,|max| < 1e-6 (cannot carry flux)
      * **flux_forced** -- min*max > 0 (range excludes 0: obligate for growth)
      * **flexible**    -- range spans 0
    Returns ``{model_id, base_flux, n_blocked, n_forced, n_flexible,
    reactions:[{rxn,min,max,span,kind}...]}`` (flux-forced + widest-range first,
    capped at ``top_n``), or ``{model_id, error}``."""
    from cobra.flux_analysis import flux_variability_analysis
    try:
        model = load_json_model(str(MODELS_DIR / f"{model_id}.json"))
        override_bounds(model, baseline_map)
        apply_media(model)
        bio_rxn = find_biomass_reaction(model)
        if bio_rxn is None:
            return {"model_id": model_id, "error": "no_biomass"}
        model.objective = bio_rxn
        sol = model.optimize()
        base_flux = float(sol.objective_value) if (sol.status == "optimal" and sol.objective_value is not None) else 0.0
        # processes=1: cobra FVA would otherwise spawn its own pool, which fails
        # inside our (daemonic) multiprocessing workers.
        fr = flux_variability_analysis(model, fraction_of_optimum=fraction, processes=1)

        # aggregate cobra reaction ranges onto SEED ids (net min/max across copies)
        per_seed = {}
        for rxn in model.reactions:
            s = seed_id(rxn)
            if not s or rxn.id not in fr.index:
                continue
            lo = float(fr.loc[rxn.id, "minimum"])
            hi = float(fr.loc[rxn.id, "maximum"])
            cur = per_seed.get(s)
            if cur is None:
                per_seed[s] = [lo, hi]
            else:
                cur[0] = min(cur[0], lo)
                cur[1] = max(cur[1], hi)
        rows = []
        n_b = n_f = n_x = 0
        for s, (lo, hi) in per_seed.items():
            if abs(lo) < 1e-6 and abs(hi) < 1e-6:
                kind = "blocked"; n_b += 1
            elif lo * hi > 1e-12:
                kind = "flux_forced"; n_f += 1
            else:
                kind = "flexible"; n_x += 1
            rows.append({"rxn": s, "min": lo, "max": hi, "span": hi - lo, "kind": kind})
        # flux-forced first (obligate), then widest range
        rows.sort(key=lambda d: (0 if d["kind"] == "flux_forced" else 1, -d["span"], d["rxn"]))
        return {"model_id": model_id, "base_flux": base_flux,
                "n_blocked": n_b, "n_forced": n_f, "n_flexible": n_x,
                "reactions": rows[:top_n]}
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
