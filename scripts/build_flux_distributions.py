#!/usr/bin/env python3
"""Internal-flux distributions for the ModelSEED v2 model sets, loops discounted.

The question this answers: a genome-scale model's biomass flux is ~1.2
mmol gDW-1 h-1, which *sounds* small, but biomass is one reaction with ~60
precursors -- so how active are the internal reactions that feed it? This
script measures the distribution of internal reaction fluxes across every
model in each panel, and reports median and IQR.

THE LOOP PROBLEM
A flux-balance solution may contain thermodynamically infeasible cycles: sets
of reactions carrying flux around a closed loop with no net exchange, driven
only by the arbitrary +-1000 bound rather than by any thermodynamic force.
Their magnitudes are an artefact of that bound, so including them in a median
measures the bound, not the biology. Two filters are computed per model:

``A`` -- *the +-1000 test.* FVA on the model with its medium applied; a
    reaction is flagged when it can reach ``>= +LOOP_BOUND`` AND
    ``<= -LOOP_BOUND``. Cheap and intuitive: nothing bounded by real substrate
    supply can run a thousand units in both directions.

``B`` -- *the closed-exchange test* (default). Every boundary reaction is shut,
    then FVA; anything that can still carry flux must be cycling internally,
    because no matter can enter or leave. This is the standard type-III /
    infeasible-loop definition (Schellenberger, Lewis & Palsson 2011).

**A and B are not interchangeable, and which one is right depends on the
medium.** Measured here, not assumed:

  * On the genome-scale panels (uptake -10) A is a strict *subset* of B --
    sound but incomplete, catching 79 of 140 loop reactions on the first
    model. Any reaction that can run 1000 units when nothing is supplied
    faster than 10 is certainly looping; but a loop member can also be capped
    below 1000 in one direction and so escape A.
  * On the core panel the complete medium supplies every compound at **-1000**,
    so a reaction can legitimately carry 1000 units without any loop at all.
    There A is no longer even a subset of B (23 flagged against B's 14): it
    over-flags, and the quantity it measures is the uptake bound.

B is therefore the default, because it is the only one of the two that does
not change meaning with the medium. A is computed and reported alongside it
because it is the criterion originally asked for and the comparison is the
point. ``--filter`` selects which one the headline numbers use.

THE FLUX VECTOR IS pFBA, NOT PLAIN FBA
Plain FBA picks an arbitrary vertex among the alternate optima, so its
individual flux values are not reproducible and a median over them means
little. pFBA (Lewis et al. 2010) minimises total absolute flux subject to
optimal biomass, which is unique in practice and is the conservative choice
here: if even the minimal-flux solution shows internal reactions running well
above the biomass flux, the point holds a fortiori.

Only reactions carrying flux are counted. Including the zeros would make every
median 0 -- most of a genome-scale model is inactive on any one medium -- which
measures model size rather than metabolic activity. ``ACTIVE_TOL`` is the cut.

Outputs land in ``results/flux_distributions/<panel>/``:

  per_model.tsv     one row per model: biomass flux, reaction/active/loop
                    counts, and the flux quantiles under each filter
  fluxes.npz        every active internal |flux| in the panel, with the model
                    index and a per-value loop flag, so any filter can be
                    re-applied at plot time without re-solving
  summary.json      pooled quantiles per filter, and the loop accounting

Usage::

    python3 scripts/build_flux_distributions.py                  # all panels
    python3 scripts/build_flux_distributions.py --sets ms2_gsm -j 64
    python3 scripts/build_flux_distributions.py --limit 20       # smoke test
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
from cma import models as M  # noqa: E402
from cma import paths  # noqa: E402

#: Model sets to sweep, in figure order. These are the same three panels the
#: baseline biomass figures use (see build_baseline_biomass_taxonomy.py).
PANELS = ("core_kegg2", "ms2_gsm", "ms2_gsm_auxo")

#: A reaction counts as carrying flux above this. Below it, pFBA's value is
#: solver noise rather than metabolism.
ACTIVE_TOL = 1e-9
#: Biomass flux above which a model is treated as growing (repo convention,
#: matches growth_heuristics.GROWTH_THRESHOLD).
GROWTH_THRESHOLD = 1e-6
#: The arbitrary bound a looping reaction runs away to. Reactions reaching
#: +this and -this in FVA are flagged by filter A.
LOOP_BOUND = 1000.0
#: Tolerance on hitting LOOP_BOUND exactly.
LOOP_EPS = 1e-2
#: Under filter B, flux this far from zero with every boundary shut is a loop.
CLOSED_TOL = 1e-6

#: Bit flags stored per flux value in fluxes.npz.
FLAG_LOOP_A = 1
FLAG_LOOP_B = 2

QUANTILES = (0.05, 0.25, 0.50, 0.75, 0.95)


# ---------------------------------------------------------------------------
# per-model work
# ---------------------------------------------------------------------------

_SET = None  # per-worker cache, set by _init_worker


def _init_worker(set_name: str) -> None:
    global _SET
    _SET = M.get(set_name)
    # One LP per reaction per direction; cobra would otherwise try to spawn its
    # own pool inside this already-forked worker.
    os.environ.setdefault("OMP_NUM_THREADS", "1")


def analyse_one(model_id: str) -> dict:
    """Flux distribution and loop membership for one model.

    Returns a record with the per-model summary plus three parallel arrays
    (``flux``, ``flag``) describing every active internal reaction.
    """
    from cobra.flux_analysis import flux_variability_analysis, pfba

    ms = _SET
    rec: dict = {"model_id": model_id, "error": "", "grows": False,
                 "biomass_flux": 0.0, "n_reactions": 0, "n_internal": 0,
                 "n_active": 0, "n_loop_a": 0, "n_loop_b": 0,
                 "flux": np.empty(0, dtype=np.float32),
                 "flag": np.empty(0, dtype=np.uint8)}
    try:
        model = M.load_model(ms, model_id)
        M.apply_media(model, M.media(ms.media))
        bio = M.find_biomass(model, ms.biomass)
        if bio is None:
            rec["error"] = "no_biomass"
            return rec
        model.objective = bio
        rec["n_reactions"] = len(model.reactions)

        sol = model.optimize()
        if sol.status != "optimal" or sol.objective_value is None:
            rec["error"] = f"status:{sol.status}"
            return rec
        rec["biomass_flux"] = float(sol.objective_value)
        if rec["biomass_flux"] <= GROWTH_THRESHOLD:
            return rec           # a non-grower has no flux distribution
        rec["grows"] = True

        # Internal = not a boundary reaction (cobra classifies exchange/demand/
        # sink itself) and not a biomass reaction.
        boundary = {r.id for r in model.boundary}
        internal = [r.id for r in model.reactions
                    if r.id not in boundary and not r.id.startswith("bio")]
        rec["n_internal"] = len(internal)

        # --- the flux vector: parsimonious, so it is reproducible -----------
        fluxes = pfba(model).fluxes

        # --- filter A: can it run LOOP_BOUND in BOTH directions? ------------
        fva = flux_variability_analysis(model, fraction_of_optimum=0.0,
                                        processes=1)
        loop_a = set(fva.index[(fva["maximum"] >= LOOP_BOUND - LOOP_EPS)
                               & (fva["minimum"] <= -LOOP_BOUND + LOOP_EPS)])

        # --- filter B: with nothing entering or leaving, what still moves? --
        # The context manager reverts the bounds, which is far cheaper than
        # re-reading the model from disk.
        with model:
            for rid in boundary:
                model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
            fva_c = flux_variability_analysis(model, fraction_of_optimum=0.0,
                                              processes=1)
        loop_b = set(fva_c.index[(fva_c["maximum"].abs() > CLOSED_TOL)
                                 | (fva_c["minimum"].abs() > CLOSED_TOL)])

        rec["n_loop_a"] = len(loop_a & set(internal))
        rec["n_loop_b"] = len(loop_b & set(internal))

        vals, flags = [], []
        for rid in internal:
            v = abs(float(fluxes[rid]))
            if v <= ACTIVE_TOL:
                continue
            vals.append(v)
            flags.append((FLAG_LOOP_A if rid in loop_a else 0)
                         | (FLAG_LOOP_B if rid in loop_b else 0))
        rec["flux"] = np.asarray(vals, dtype=np.float32)
        rec["flag"] = np.asarray(flags, dtype=np.uint8)
        rec["n_active"] = len(vals)
    except Exception as exc:                                  # noqa: BLE001
        rec["error"] = f"{type(exc).__name__}: {exc}"
    return rec


# ---------------------------------------------------------------------------
# aggregation
# ---------------------------------------------------------------------------

def quantiles(arr: np.ndarray) -> dict:
    if arr.size == 0:
        return {f"q{int(q * 100):02d}": None for q in QUANTILES} | {"n": 0}
    out = {f"q{int(q * 100):02d}": float(np.quantile(arr, q)) for q in QUANTILES}
    out["n"] = int(arr.size)
    out["mean"] = float(arr.mean())
    out["max"] = float(arr.max())
    return out


def sweep(set_name: str, jobs: int, limit: int | None, out_root: Path) -> dict:
    ms = M.get(set_name)
    ids = ms.model_ids()
    if limit:
        ids = ids[:limit]
    out_dir = out_root / set_name
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{set_name}] {len(ids):,} models, {jobs} workers", flush=True)
    t0 = time.time()
    rows, all_flux, all_flag, all_idx, all_bio = [], [], [], [], []
    done = 0
    with ProcessPoolExecutor(max_workers=jobs, initializer=_init_worker,
                             initargs=(set_name,)) as pool:
        for rec in pool.map(analyse_one, ids, chunksize=1):
            done += 1
            if done % 250 == 0 or done == len(ids):
                el = time.time() - t0
                print(f"[{set_name}] {done:,}/{len(ids):,}  {el:,.0f}s  "
                      f"({done / el:.1f}/s)", flush=True)
            rows.append({k: rec[k] for k in
                         ("model_id", "error", "grows", "biomass_flux",
                          "n_reactions", "n_internal", "n_active",
                          "n_loop_a", "n_loop_b")})
            if rec["flux"].size:
                mi = len(all_bio)
                all_flux.append(rec["flux"])
                all_flag.append(rec["flag"])
                all_idx.append(np.full(rec["flux"].size, mi, dtype=np.int32))
                all_bio.append(rec["biomass_flux"])

    flux = np.concatenate(all_flux) if all_flux else np.empty(0, np.float32)
    flag = np.concatenate(all_flag) if all_flag else np.empty(0, np.uint8)
    idx = np.concatenate(all_idx) if all_idx else np.empty(0, np.int32)
    bio = np.asarray(all_bio, dtype=np.float64)

    np.savez_compressed(out_dir / "fluxes.npz", flux=flux, flag=flag,
                        model_index=idx, model_biomass=bio)

    with open(out_dir / "per_model.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    keep_a = (flag & FLAG_LOOP_A) == 0
    keep_b = (flag & FLAG_LOOP_B) == 0
    growers = [r for r in rows if r["grows"]]
    summary = {
        "model_set": set_name,
        "media": ms.media,
        "n_models": len(ids),
        "n_growers": len(growers),
        "n_errors": sum(1 for r in rows if r["error"]),
        "n_flux_values": int(flux.size),
        "loops": {
            "a_plus_minus_1000": int((~keep_a).sum()),
            "b_closed_exchange": int((~keep_b).sum()),
            "a_and_b": int(((~keep_a) & (~keep_b)).sum()),
            "a_only": int(((~keep_a) & keep_b).sum()),
            "b_only": int((keep_a & (~keep_b)).sum()),
            "a_is_subset_of_b": bool(((~keep_a) & keep_b).sum() == 0),
        },
        "pooled": {
            "unfiltered": quantiles(flux),
            "filter_a": quantiles(flux[keep_a]),
            "filter_b": quantiles(flux[keep_b]),
        },
        # The headline framing: internal flux measured in units of that
        # model's own biomass flux.
        "pooled_per_biomass": {
            "unfiltered": quantiles(flux / bio[idx]) if flux.size else quantiles(flux),
            "filter_a": quantiles(flux[keep_a] / bio[idx[keep_a]]) if flux.size else quantiles(flux),
            "filter_b": quantiles(flux[keep_b] / bio[idx[keep_b]]) if flux.size else quantiles(flux),
        },
        "biomass_flux": quantiles(bio),
        "elapsed_s": round(time.time() - t0, 1),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"[{set_name}] done in {summary['elapsed_s']:,.0f}s -> {out_dir}",
          flush=True)
    return summary


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sets", nargs="*", default=list(PANELS),
                    help=f"model sets to sweep (default: {' '.join(PANELS)})")
    ap.add_argument("-j", "--jobs", type=int, default=64)
    ap.add_argument("--limit", type=int, default=None,
                    help="only the first N models of each set (smoke test)")
    ap.add_argument("--out", default=None,
                    help="output root (default results/flux_distributions)")
    args = ap.parse_args()

    out_root = (Path(args.out) if args.out
                else paths.analysis_dir() / "results" / "flux_distributions")
    out_root.mkdir(parents=True, exist_ok=True)

    index = {}
    for name in args.sets:
        index[name] = sweep(name, args.jobs, args.limit, out_root)
    (out_root / "index.json").write_text(json.dumps(index, indent=2))
    print(f"wrote {out_root}/index.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
