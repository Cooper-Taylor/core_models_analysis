"""The built-in pipeline stages.

Each is one unit of "info" ``beginPipeline`` can produce. List them with
``beginPipeline --list``; select with ``--only a,b``; the default is all.

Adding an output is one ``@stage`` function here -- no new driver script, no
argparse, no multiprocessing scaffold.
"""

from __future__ import annotations

import logging
import multiprocessing as mp
from collections import Counter

from . import models as models_reg
from .pipeline import Ctx, stage

GROWTH_THRESHOLD = 1e-6


# ---------------------------------------------------------------------------
# helpers shared by the FBA stages
# ---------------------------------------------------------------------------
_W: dict = {}


def _init_worker(model_set, media, dmaps, only_changed_vs):
    logging.getLogger("cobra").setLevel(logging.ERROR)
    _W["ms"] = model_set
    _W["media"] = media
    _W["dmaps"] = dmaps
    _W["base"] = only_changed_vs


def _solve(model_id: str) -> dict:
    """Load once, then solve under each direction map. The shared prologue."""
    ms, media, dmaps, base = _W["ms"], _W["media"], _W["dmaps"], _W["base"]
    row = {"model_id": model_id, "runs": {}}
    try:
        template = models_reg.load_model(ms, model_id)
    except Exception as exc:  # noqa: BLE001
        row["error"] = f"{type(exc).__name__}: {exc}"
        return row
    for label, dmap in dmaps.items():
        m = template.copy()
        touched = 0
        if dmap:
            for rxn in m.reactions:
                seed = models_reg.seed_id(rxn, ms.seed_key)
                if not seed:
                    continue
                op = dmap.get(seed)
                if op is None:
                    continue
                if base is not None and base.get(seed) == op:
                    continue
                # '?' (no opinion) becomes fully reversible, matching
                # growth_heuristics._bounds_for_rev and ModelSEED's own handling
                # of unknown direction. It is the permissive choice: a map full
                # of '?' REMOVES constraints rather than adding them.
                if op == ">":
                    rxn.lower_bound, rxn.upper_bound = 0.0, 1000.0
                elif op == "<":
                    rxn.lower_bound, rxn.upper_bound = -1000.0, 0.0
                else:
                    rxn.lower_bound, rxn.upper_bound = -1000.0, 1000.0
                touched += 1
        models_reg.apply_media(m, media)
        bio = models_reg.find_biomass(m, ms.biomass)
        if bio is None:
            row["runs"][label] = {"status": "no_biomass", "growth": 0.0,
                                  "grows": False, "n_overrides": touched}
            continue
        m.objective = bio
        sol = m.optimize()
        flux = float(sol.objective_value) if sol.objective_value is not None else 0.0
        row["runs"][label] = {
            "status": sol.status, "growth": flux,
            "grows": sol.status == "optimal" and flux > GROWTH_THRESHOLD,
            "n_overrides": touched, "biomass_rxn": bio.id,
        }
    return row


def _run_panel(ctx: Ctx, dmaps: dict) -> list:
    ids = ctx.model_ids[: ctx.limit] if ctx.limit else ctx.model_ids
    n_workers = max(1, min(ctx.workers, len(ids)))
    args = (ctx.model_set, ctx.media, dmaps, None)
    if n_workers == 1:
        _init_worker(*args)
        return [_solve(i) for i in ids]
    mpctx = mp.get_context("fork")
    with mpctx.Pool(n_workers, initializer=_init_worker, initargs=args) as pool:
        return list(pool.imap_unordered(_solve, ids, chunksize=4))


def _all_dmaps(ctx: Ctx) -> dict:
    """Every direction set to solve under: on-disk bounds, plus each input."""
    d = {"on_disk": {}}
    d.update(ctx.direction_maps)
    d.update(ctx.results.get("heuristics", {}).get("maps", {}))
    return d


# ---------------------------------------------------------------------------
# inventory stages (cheap, no FBA)
# ---------------------------------------------------------------------------
@stage("inputs", "Echo every resolved input: model set, media, maps, heuristics, sources",
       tags=("inventory",))
def _inputs(ctx: Ctx) -> dict:
    ms = ctx.model_set
    return {
        "model_set": {
            "name": ms.name, "dir": str(ms.models_dir()), "glob": ms.id_glob,
            "loader": ms.loader, "n_models_available": len(ms.model_ids()),
            "n_models_selected": len(ctx.model_ids[: ctx.limit] if ctx.limit else ctx.model_ids),
            "media": ms.media, "results_namespace": ms.results_namespace,
        },
        "media": {"name": ctx.media.name, "n_compounds": len(ctx.media.loader())},
        "direction_maps": {k: {"n_reactions": len(v),
                               "operators": dict(Counter(v.values()))}
                           for k, v in ctx.direction_maps.items()},
        "heuristics": sorted(ctx.heuristics),
        "thermo_sources": {k: {"n_reactions": len(v)} for k, v in ctx.thermo.items()},
        "workers": ctx.workers, "limit": ctx.limit,
    }


@stage("model_summary", "Size distribution, compartments and SEED-id coverage of the models",
       tags=("inventory",))
def _model_summary(ctx: Ctx) -> dict:
    ms = ctx.model_set
    ids = ctx.model_ids[: ctx.limit] if ctx.limit else ctx.model_ids
    sample = ids[: min(len(ids), ctx.options.get("summary_sample", 200))]
    n_rxn, n_met, seeds, cpts, bios = [], [], set(), Counter(), Counter()
    no_seed = 0
    unreadable = []
    for mid in sample:
        try:
            m = models_reg.load_model(ms, mid)
        except Exception as exc:  # noqa: BLE001
            # growth tolerates an unreadable model; this stage used to abort the
            # whole run on the first one.
            unreadable.append({"model_id": mid, "error": f"{type(exc).__name__}: {exc}"})
            continue
        n_rxn.append(len(m.reactions))
        n_met.append(len(m.metabolites))
        for r in m.reactions:
            s = models_reg.seed_id(r, ms.seed_key)
            if s:
                seeds.add(s)
            else:
                no_seed += 1
        for mm in m.metabolites:
            cpts[mm.compartment] += 1
        b = models_reg.find_biomass(m, ms.biomass)
        bios[b.id if b else "NONE"] += 1
    n_rxn.sort()
    return {
        "_n": len(sample),
        "n_models_scanned": len(sample),
        "n_models_selected": len(ids),
        "sampled": len(sample) < len(ids),
        "sampling_note": (f"statistics cover {len(sample)} of {len(ids)} models; "
                          "raise with --option summary_sample=N")
        if len(sample) < len(ids) else None,
        "reactions_per_model": {"min": n_rxn[0], "median": n_rxn[len(n_rxn) // 2],
                                "max": n_rxn[-1]} if n_rxn else {},
        "metabolites_per_model": {"min": min(n_met), "max": max(n_met)} if n_met else {},
        "distinct_seed_reactions": len(seeds),
        "reactions_without_seed_id": no_seed,
        "compartments": dict(cpts),
        "biomass_reactions": dict(bios),
        "n_unreadable": len(unreadable),
        "unreadable": unreadable[:20],
    }


@stage("direction_summary", "Operator composition and pairwise agreement of the direction maps",
       needs=("heuristics",), requires_any=("direction_maps", "heuristics"),
       tags=("inventory", "directions"))
def _direction_summary(ctx: Ctx) -> dict:
    maps = dict(ctx.direction_maps)
    maps.update(ctx.results.get("heuristics", {}).get("maps", {}))
    per = {k: {"n_reactions": len(v), "operators": dict(Counter(v.values()))}
           for k, v in maps.items()}
    labels = sorted(maps)
    agree = {}
    for i, a in enumerate(labels):
        for b in labels[i + 1:]:
            shared = set(maps[a]) & set(maps[b])
            same = sum(1 for k in shared if maps[a][k] == maps[b][k])
            agree[f"{a} vs {b}"] = {
                "shared": len(shared), "agree": same,
                "disagree": len(shared) - same,
                "pct_agree": round(100.0 * same / len(shared), 2) if shared else None,
            }
    return {"_n": len(per), "per_map": per, "pairwise_agreement": agree}


# ---------------------------------------------------------------------------
# heuristics -> direction maps
# ---------------------------------------------------------------------------
@stage("heuristics", "Run each heuristic's cascade over MSDB to produce a direction map",
       requires=("heuristics",), tags=("directions",))
def _heuristics(ctx: Ctx) -> dict:
    import reversibility_lib as lib

    from . import paths

    paths.add_msdb_to_syspath()
    import export_thermo_variants as ex

    rxns = ex.load_msdb_reactions(use_cache=True)
    out, maps = {}, {}
    for label, cfg in ctx.heuristics.items():
        res = lib.run_cascade(rxns, db_level=ctx.options.get("db_level", "EQ"), cfg=cfg)
        ops = {r: v[1] if isinstance(v, (tuple, list)) else v for r, v in res.items()}
        maps[f"heuristic:{label}"] = ops
        out[label] = {"n_reactions": len(ops), "operators": dict(Counter(ops.values()))}
    return {"_n": len(out), "per_heuristic": out, "maps": maps}


@stage("thermo_summary", "Coverage and ΔG′° distribution of each thermodynamic source",
       requires=("thermo",), tags=("inventory", "thermo"))
def _thermo_summary(ctx: Ctx) -> dict:
    out = {}
    for label, table in ctx.thermo.items():
        dgs = sorted(v[0] for v in table.values())
        errs = sorted(v[1] for v in table.values())
        n = len(dgs)
        # eQuilibrator carries a few thousand records whose uncertainty is a
        # placeholder several orders of magnitude above any real value. A mean
        # over those is meaningless and reverses the ranking of the sources, so
        # report the median and say how many outliers there are rather than
        # hiding them inside an average.
        huge = [e for e in errs if e > 1e3]
        out[label] = {
            "n_reactions": n,
            "dg_kcal_per_mol": {"min": dgs[0], "p25": dgs[n // 4], "median": dgs[n // 2],
                                "p75": dgs[3 * n // 4], "max": dgs[-1]} if n else {},
            "uncertainty_kcal_per_mol": {
                "median": round(errs[n // 2], 3), "p75": round(errs[3 * n // 4], 3),
                "max": round(errs[-1], 3),
            } if n else {},
            "n_uncertainty_over_1000": len(huge),
            "note": (f"{len(huge)} record(s) carry a placeholder uncertainty > 1000 "
                     "kcal/mol; the median is reported because a mean is dominated "
                     "by them") if huge else None,
        }
    return {"_n": len(out), "per_source": out}


# ---------------------------------------------------------------------------
# FBA stages
# ---------------------------------------------------------------------------
@stage("growth", "FBA biomass flux for every model under every direction set",
       needs=("heuristics",), heavy=True, tags=("fba",))
def _growth(ctx: Ctx) -> dict:
    dmaps = _all_dmaps(ctx)
    rows = _run_panel(ctx, dmaps)
    n_failed = sum(1 for r in rows if r.get("error"))
    per_variant = {}
    for label in dmaps:
        vals = [r["runs"][label] for r in rows if label in r.get("runs", {})]
        grows = sum(1 for v in vals if v["grows"])
        fluxes = [v["growth"] for v in vals if v["grows"]]
        per_variant[label] = {
            "n_models": len(vals), "n_models_attempted": len(rows),
            "n_failed_to_load": n_failed, "n_growers": grows,
            "frac_growers": round(grows / len(vals), 4) if vals else None,
            "mean_growth_of_growers": round(sum(fluxes) / len(fluxes), 4) if fluxes else 0.0,
            "mean_overrides": round(sum(v["n_overrides"] for v in vals) / len(vals), 1)
            if vals else 0,
        }
    warnings = []
    for label, v in per_variant.items():
        if label != "on_disk" and v["mean_overrides"] == 0:
            warnings.append(
                f"{label!r} overrode 0 reactions in every model -- its reaction ids "
                "probably do not match this model set's seed.reaction annotations, "
                "so its growth numbers are identical to on_disk by accident")
    for w in warnings:
        print(f"    WARNING: {w}", flush=True)
    return {"_n": len(rows), "per_variant": per_variant, "warnings": warnings,
            "per_model": sorted(rows, key=lambda r: r["model_id"]),
            "errors": [r for r in rows if r.get("error")]}


@stage("growth_diff", "Which models change grow/no-grow between direction sets",
       needs=("growth",), tags=("fba", "compare"))
def _growth_diff(ctx: Ctx) -> dict:
    g = ctx.results.get("growth") or {}
    rows = g.get("per_model", [])
    labels = sorted(g.get("per_variant", {}))
    if ctx.baseline and ctx.baseline not in labels:
        raise ValueError(
            f"--baseline {ctx.baseline!r} is not one of the direction sets in this run "
            f"({', '.join(labels)}). Silently diffing against something else would "
            "misreport every number in this table."
        )
    base = ctx.baseline or ("on_disk" if "on_disk" in labels else
                            (labels[0] if labels else None))
    if base is None:
        return {"_n": 0, "note": "no growth results to diff"}
    out = {}
    for label in labels:
        if label == base:
            continue
        gained, lost, changed = [], [], []
        for r in rows:
            a, b = r["runs"].get(base), r["runs"].get(label)
            if not a or not b:
                continue
            if b["grows"] and not a["grows"]:
                gained.append(r["model_id"])
            elif a["grows"] and not b["grows"]:
                lost.append(r["model_id"])
            elif abs(a["growth"] - b["growth"]) > 1e-6:
                changed.append({"model_id": r["model_id"],
                                "from": round(a["growth"], 6), "to": round(b["growth"], 6)})
        cap = ctx.options.get("list_cap", 200)
        out[label] = {"baseline": base, "n_gained_growth": len(gained),
                      "n_lost_growth": len(lost), "n_flux_changed": len(changed),
                      "gained": gained[:cap], "lost": lost[:cap], "changed": changed[:cap],
                      "lists_truncated_to": cap if max(len(gained), len(lost),
                                                       len(changed)) > cap else None}
    return {"_n": len(out), "baseline": base, "per_variant": out}


@stage("direction_changes",
       "Reactions whose direction differs between the supplied direction sets",
       requires_any=("direction_maps", "heuristics"), tags=("directions", "compare"))
def _direction_changes(ctx: Ctx) -> dict:
    """Which reactions the direction sets disagree about, and how.

    This is a comparison of the MAPS, not a measurement of growth. It was
    originally named ``reaction_effects`` and described as "per-reaction growth
    impact", which it never computed -- and it declared a dependency on
    ``growth``, so asking for it forced the whole FBA sweep to produce a table
    that never used the result. Measuring a genuine per-reaction growth impact
    needs one FBA per reaction per model; that is a separate stage, not yet
    written.
    """
    maps = _all_dmaps(ctx)
    maps = {k: v for k, v in maps.items() if v}
    labels = sorted(maps)
    if len(labels) < 2:
        return {"_n": 0,
                "note": ("needs at least two non-empty direction sets to compare; "
                         f"got {labels or 'none'}"),
                "sets_compared": labels}

    ref = ctx.baseline if ctx.baseline in labels else labels[0]
    others = [lab for lab in labels if lab != ref]
    changed: dict = {}
    for label in others:
        for rxn, op in maps[label].items():
            base_op = maps[ref].get(rxn)
            if base_op is not None and base_op != op:
                changed.setdefault(rxn, {})[label] = [base_op, op]

    transitions = Counter()
    for d in changed.values():
        for _label, (a, b) in d.items():
            transitions[f"{a} -> {b}"] += 1
    top = sorted(changed.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    n_top = ctx.options.get("top_n", 200)
    return {
        "_n": len(changed),
        "reference_set": ref,
        "sets_compared": others,
        "n_reactions_differing": len(changed),
        "transition_counts": dict(transitions.most_common()),
        "top_reactions": [{"rxn_id": r, "transitions": d} for r, d in top[:n_top]],
        "truncated_to": n_top if len(top) > n_top else None,
    }


@stage("report", "A human-readable markdown summary of every stage that ran",
       tags=("report",))
def _report(ctx: Ctx) -> dict:
    lines = ["# Pipeline run", ""]
    inp = ctx.results.get("inputs", {})
    if inp:
        ms = inp.get("model_set", {})
        lines += [
            "## Inputs", "",
            f"- **Model set**: `{ms.get('name')}` — {ms.get('n_models_selected')} of "
            f"{ms.get('n_models_available')} models from `{ms.get('dir')}`",
            f"- **Media**: {inp.get('media', {}).get('name')} "
            f"({inp.get('media', {}).get('n_compounds')} compounds)",
            f"- **Direction maps**: {', '.join(inp.get('direction_maps', {})) or 'none'}",
            f"- **Heuristics**: {', '.join(inp.get('heuristics', [])) or 'none'}",
            f"- **Thermo sources**: {', '.join(inp.get('thermo_sources', {})) or 'none'}",
            "",
        ]
    ms_sum = ctx.results.get("model_summary")
    if ms_sum:
        r = ms_sum.get("reactions_per_model", {})
        lines += ["## Models", "",
                  f"- Reactions per model: min {r.get('min')}, median {r.get('median')}, "
                  f"max {r.get('max')}",
                  f"- Distinct SEED reactions: {ms_sum.get('distinct_seed_reactions')}",
                  f"- Reactions with no SEED id: {ms_sum.get('reactions_without_seed_id')}", ""]
    g = ctx.results.get("growth")
    if g:
        lines += ["## Growth", "", "| direction set | growers | fraction | mean flux | overrides |",
                  "|---|---:|---:|---:|---:|"]
        for label, v in g.get("per_variant", {}).items():
            lines.append(f"| {label} | {v['n_growers']}/{v['n_models']} | "
                         f"{v['frac_growers']} | {v['mean_growth_of_growers']} | "
                         f"{v['mean_overrides']} |")
        lines.append("")
    d = ctx.results.get("growth_diff")
    if d and d.get("per_variant"):
        lines += [f"## Change vs `{d.get('baseline')}`", "",
                  "| direction set | gained | lost | flux changed |", "|---|---:|---:|---:|"]
        for label, v in d["per_variant"].items():
            lines.append(f"| {label} | {v['n_gained_growth']} | {v['n_lost_growth']} | "
                         f"{v['n_flux_changed']} |")
        lines.append("")
    ds = ctx.results.get("direction_summary")
    if ds and ds.get("pairwise_agreement"):
        lines += ["## Direction-map agreement", "", "| pair | shared | agree | % |",
                  "|---|---:|---:|---:|"]
        for pair, v in ds["pairwise_agreement"].items():
            lines.append(f"| {pair} | {v['shared']} | {v['agree']} | {v['pct_agree']} |")
        lines.append("")
    text = "\n".join(lines)
    (ctx.out_dir / "REPORT.md").write_text(text)
    return {"_n": len(lines), "path": str(ctx.out_dir / "REPORT.md"), "markdown": text}
