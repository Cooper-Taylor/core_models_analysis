#!/usr/bin/env python3
"""Choose the largest set of reactions on which eQuilibrator and
dGPredictor-ModelSEED can be trusted to agree, by explicit optimisation.

THE PROBLEM
-----------
Let R be the reconciled key subset, x_i = dG_eq(i), y_i = dG_dgp(i). For S subset R:

    maximise    |S|                              coverage
    subject to  CCC(S)          >= c*            agreement about the line y = x
                RMSE(S)         <= E*            error magnitude, kcal/mol
                |slope(S) - 1|  <= delta         no systematic scaling
                |S|             >= N_min         coverage floor
                per-|dG|-decile retention >= rho stratum floor

Coverage is the maximand and quality the constraint: coverage is what you want
to be greedy about, quality is what you need guaranteed.

WHY NOT MAXIMISE PEARSON r
--------------------------
Because it is degenerate here, measurably. On the same pool:

    |dG_eq| >  50 :  n=1451,  r=0.774,  median |delta| = 11.91
    |dG_eq| <= 10 :  n=7835,  r=0.366,  median |delta| =  2.52

The set that agrees 4.7x better in kcal/mol scores less than half the
correlation, because r rewards spread rather than agreement. An optimiser told
to maximise r would throw away the near-zero bulk where most of metabolism sits.

So agreement is measured by Lin's concordance correlation coefficient, which
penalises scatter AND departure from the 1:1 line, backed by an explicit slope
guard. Plain r is still reported, as a diagnostic only.

THE STRATUM FLOOR
-----------------
Without it, any quality bar can be met by deleting a whole energy regime -- the
same degeneracy re-entering through the back door. Requiring a minimum retention
in every |dG| decile makes that inadmissible.

TWO SOLVERS
-----------
oracle  Selects directly on |x-y|. That is circular -- it is the outcome we are
        trying to guarantee -- so it is a BOUND, never a deliverable. It is
        exactly solvable: for fixed k the min-RMSE subset is the k smallest
        residuals and prefix-RMSE is monotone in k, so sort ascending and take
        the longest feasible prefix. O(n log n).

rule    A threshold vector over features knowable WITHOUT seeing the
        disagreement, so it applies to reactions outside the comparison set.
        Fitted by coordinate ascent from a coarse grid, cross-checked with
        differential evolution, and 5-fold cross-validated -- without CV the
        thresholds are fitted to noise and the reported quality is optimistic.

The oracle/rule gap is itself a result: it quantifies how much of the
disagreement is simply not predictable from selection-time information.

OUTPUTS (results/eq_vs_dgpms/)
    consensus_frontier.tsv   the trade-off surface, oracle and rule
    consensus_rule.json      fitted thresholds + CV metrics; read by load_selector
    consensus_selected.tsv   the reaction list at the chosen operating point

CONSUMERS import ``load_selector()``, which mirrors the contract of
``build_dgpredictor_kegg_mask.load_mask()``.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
DATA = Path(os.environ.get("EQDGP_OUT", str(ANALYSIS_DIR / "results" / "eq_vs_dgpms")))
KEY = DATA / "key_subset_classified.tsv"
RULE_JSON = DATA / "consensus_rule.json"

RNG = np.random.default_rng(20260806)
N_FOLDS = 5

# Operating point actually shipped. RMSE <= 3 kcal/mol is ~2x the +/-2.0 kcal/mol
# reversible band the ModelSEED cascade uses on mMdeltaG
# (reversibility_heuristics.py:327), so a set meeting it will rarely flip a
# direction call on disagreement alone. CCC 0.95 and the slope band are the
# anti-degeneracy guards rather than independently tuned targets.
TARGET = {"rmse": 3.0, "ccc": 0.95, "slope_tol": 0.15, "stratum_floor": 0.0}

# The stratum floor defaults OFF, and that is a considered choice, not a
# convenience. It was insurance against an optimiser satisfying a quality bar by
# deleting an energy regime -- a real hazard when MAXIMISING r, which rewards
# spread. This objective maximises COVERAGE, which pushes the opposite way: it
# keeps everything it can. Measured, a floor of 0.15 makes every bar below
# RMSE <= 5 infeasible, and the decile that fails is |dG| 88.9-451.8 at 6.7%
# retention -- the O2 / quinone regime, where dGPredictor is both uncertain and
# demonstrably wrong (report section 2). Discarding it is the correct answer, not
# gaming. So retention is REPORTED per decile always, and the floor is swept to
# quantify its cost, rather than silently forcing a worse solution.
STRATUM_FLOOR_SWEEP = 0.15

# Continuous features: keep a reaction when feature <= threshold. Bounds are the
# search range; None upper bound means "use the observed max" (i.e. inactive).
CONT_FEATURES = {
    "dgp_uncertainty": (0.5, 200.0),
    "eq_uncertainty": (0.05, 100.0),
    "abs_dg_eq": (1.0, 2000.0),
    "n_participants": (2.0, 40.0),
    "max_carbon": (1.0, 120.0),
    "total_arom_rings": (0.0, 40.0),
    "abs_net_proton": (0.0, 20.0),
}


# ----------------------------------------------------------------- metrics
def ccc(x: np.ndarray, y: np.ndarray) -> float:
    """Lin's concordance correlation coefficient: agreement about y = x.

    Unlike Pearson r this is maximal only when the points lie on the identity
    line, so it cannot be inflated by choosing a widely-spread subset.
    """
    if len(x) < 3:
        return float("nan")
    vx, vy = x.var(), y.var()
    denom = vx + vy + (x.mean() - y.mean()) ** 2
    if denom <= 0:
        return float("nan")
    return float(2 * np.cov(x, y, ddof=0)[0, 1] / denom)


def metrics(x: np.ndarray, y: np.ndarray) -> dict:
    if len(x) < 3:
        return {"n": len(x), "ccc": np.nan, "r": np.nan, "slope": np.nan,
                "rmse": np.nan, "median_absdiff": np.nan}
    d = x - y
    slope = float(np.polyfit(x, y, 1)[0]) if x.std() > 0 else np.nan
    return {
        "n": int(len(x)),
        "ccc": ccc(x, y),
        "r": float(np.corrcoef(x, y)[0, 1]) if x.std() > 0 and y.std() > 0 else np.nan,
        "slope": slope,
        "rmse": float(np.sqrt(np.mean(d ** 2))),
        "median_absdiff": float(np.median(np.abs(d))),
    }


def feasible(m: dict, tgt: dict, n_min: int = 0) -> bool:
    return (m["n"] >= n_min
            and np.isfinite(m["ccc"]) and m["ccc"] >= tgt["ccc"]
            and m["rmse"] <= tgt["rmse"]
            and np.isfinite(m["slope"]) and abs(m["slope"] - 1.0) <= tgt["slope_tol"])


def stratum_retention(mask: np.ndarray, strata: np.ndarray) -> np.ndarray:
    """Retained fraction per |dG| decile. The optimiser must not be allowed to
    satisfy a quality bar by deleting an entire energy regime."""
    return np.array([mask[strata == b].mean() if (strata == b).any() else 1.0
                     for b in np.unique(strata)])


# ------------------------------------------------------------------ oracle
def oracle_frontier(x: np.ndarray, y: np.ndarray, bars: list[float]) -> pd.DataFrame:
    """Exact max-coverage-at-RMSE bound. Circular by construction: a benchmark."""
    res = np.abs(x - y)
    order = np.argsort(res, kind="stable")
    cum = np.cumsum(res[order] ** 2)
    prefix_rmse = np.sqrt(cum / np.arange(1, len(res) + 1))
    assert np.all(np.diff(prefix_rmse) >= -1e-12), "prefix RMSE must be monotone"

    rows = []
    for E in bars:
        feas = np.flatnonzero(prefix_rmse <= E)
        if len(feas) == 0:
            continue
        k = int(feas[-1]) + 1
        idx = order[:k]
        rows.append({"solver": "oracle", "rmse_bar": E,
                     **metrics(x[idx], y[idx]),
                     "coverage": k / len(x)})
    return pd.DataFrame(rows)


# -------------------------------------------------------------------- rule
def _apply(df: pd.DataFrame, theta: dict, drop_classes: set[str]) -> np.ndarray:
    """Keep a reaction unless a feature we KNOW exceeds its threshold.

    A missing value is treated as "no information on this axis", not as a
    violation: eq_uncertainty is absent for 620 reactions, and excluding them
    on ignorance would silently cost coverage the rule never chose to give up.
    """
    m = np.ones(len(df), bool)
    for f, t in theta.items():
        col = df[f].to_numpy(float)
        m &= ~(np.isfinite(col) & (col > t))
    if drop_classes:
        m &= ~df["chem_class"].isin(drop_classes).to_numpy()
    return m


def fit_rule(df: pd.DataFrame, tgt: dict, x: np.ndarray, y: np.ndarray,
             strata: np.ndarray, n_rounds: int = 6, verbose: bool = True,
             de_refine: bool = True, de_maxiter: int = 60, de_popsize: int = 15,
             de_seed: int = 7):
    """Coordinate ascent maximising coverage subject to the quality constraints.

    Returns (theta, dropped_classes) or None when the bar is unreachable.

    Direction matters. Starting from all thresholds inactive is infeasible at
    tight bars, and no SINGLE coordinate move can restore feasibility from
    there -- coordinate ascent then sits at the infeasible start forever and
    silently returns it. So we first SEED a feasible point by tightening sigma
    (the strongest single predictor) as little as possible, and only then loosen
    coordinate-wise to buy back coverage.
    """
    classes = sorted(df["chem_class"].dropna().unique())
    inactive = {f: float(np.nanmax(df[f].to_numpy(float))) for f in CONT_FEATURES}
    drop: set[str] = set()

    def score(th, dr):
        m = _apply(df, th, dr)
        if m.sum() < 50:
            return -1.0, None, m
        mm = metrics(x[m], y[m])
        ok = feasible(mm, tgt)
        if ok and tgt.get("stratum_floor", 0.0) > 0 and \
                stratum_retention(m, strata).min() < tgt["stratum_floor"]:
            ok = False
        return (m.sum() / len(df) if ok else -1.0), mm, m

    # ---- seed: loosest sigma that is feasible, optionally after class drops
    sig = np.unique(np.nanquantile(df["dgp_uncertainty"].to_numpy(float),
                                   np.linspace(0.01, 1.0, 200)))
    theta = None
    for dr_try in ({}, {"Redox: quinone / quinol"}):
        for t in sig[::-1]:                       # loosest first
            trial = dict(inactive); trial["dgp_uncertainty"] = float(t)
            if score(trial, set(dr_try))[0] > 0:
                theta, drop = trial, set(dr_try)
                break
        if theta is not None:
            break
    if theta is None:
        return None                                # bar genuinely unreachable

    best = score(theta, drop)[0]
    if verbose:
        print(f"    seed: sigma<={theta['dgp_uncertainty']:.3g}"
              f"{' + drop quinones' if drop else ''} -> coverage {best:.4f}")

    for rnd in range(n_rounds):
        improved = False
        # tighten continuous thresholds
        for f, (lo, hi) in CONT_FEATURES.items():
            vals = df[f].to_numpy(float)
            grid = np.unique(np.nanquantile(vals[np.isfinite(vals)],
                                            np.linspace(0.02, 1.0, 60)))
            grid = np.clip(grid, lo, hi)
            cur = theta[f]
            for t in grid:
                trial = dict(theta); trial[f] = float(t)
                s, _, _ = score(trial, drop)
                if s > best + 1e-9:
                    best, theta, improved = s, trial, True
            if theta[f] != cur and verbose:
                print(f"    round {rnd}: {f} <= {theta[f]:.3g}  -> coverage {best:.4f}")
        # class exclusions
        for c in classes:
            if c in drop:
                continue
            trial = drop | {c}
            s, _, _ = score(theta, trial)
            if s > best + 1e-9:
                best, drop, improved = s, trial, True
                if verbose:
                    print(f"    round {rnd}: drop class '{c}'         -> coverage {best:.4f}")
        if not improved:
            break

    # ---- joint refinement.
    # Coordinate ascent can only move one threshold at a time, so from a feasible
    # point it can never LOOSEN sigma while TIGHTENING something else -- and that
    # is exactly the trade that pays here. Measured at RMSE <= 3, this stage lifts
    # coverage 23.9% -> 29.3% (+22%) by relaxing sigma from 7.9 to 14.4 and paying
    # for it with eq_uncertainty, |dG|, participant-count, aromatic-ring and
    # proton-balance caps. Differential evolution searches all thresholds jointly.
    if de_refine:
        from scipy.optimize import differential_evolution
        feats = list(CONT_FEATURES)
        bounds = [(CONT_FEATURES[f][0], inactive[f]) for f in feats]

        def neg(v):
            trial = {f: float(val) for f, val in zip(feats, v)}
            return -max(score(trial, drop)[0], 0.0)

        res = differential_evolution(neg, bounds, seed=de_seed, maxiter=de_maxiter,
                                     popsize=de_popsize, tol=1e-8, polish=False,
                                     workers=1, updating="immediate")
        cand = {f: float(v) for f, v in zip(feats, res.x)}
        if score(cand, drop)[0] > best + 1e-9:
            if verbose:
                print(f"    DE refine: coverage {best:.4f} -> {score(cand, drop)[0]:.4f}")
            theta = cand
    return theta, drop


def cross_validate(df: pd.DataFrame, tgt: dict, x: np.ndarray, y: np.ndarray,
                   strata: np.ndarray) -> pd.DataFrame:
    """Fit thresholds on 4/5 of the data, score on the held-out fifth."""
    idx = RNG.permutation(len(df))
    folds = np.array_split(idx, N_FOLDS)
    rows = []
    for i, test in enumerate(folds):
        train = np.concatenate([f for j, f in enumerate(folds) if j != i])
        got = fit_rule(df.iloc[train], tgt, x[train], y[train],
                       strata[train], verbose=False)
        if got is None:
            rows.append({"fold": i, "train_coverage": np.nan, "test_coverage": np.nan})
            continue
        th, dr = got
        mte = _apply(df.iloc[test], th, dr)
        mtr = _apply(df.iloc[train], th, dr)
        rows.append({"fold": i,
                     "train_coverage": mtr.mean(), **{f"train_{k}": v for k, v in
                                                      metrics(x[train][mtr], y[train][mtr]).items()},
                     "test_coverage": mte.mean(), **{f"test_{k}": v for k, v in
                                                     metrics(x[test][mte], y[test][mte]).items()}})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ public
def load_selector(path: Path | None = None) -> Callable[[pd.DataFrame], np.ndarray]:
    """Return a predicate mapping a feature frame -> boolean keep-mask.

    Mirrors build_dgpredictor_kegg_mask.load_mask(): if the fitted rule is
    missing the caller still runs, but unfiltered, and the warning is loud.
    """
    path = path or RULE_JSON
    if not path.exists():
        print(f"  WARNING: consensus rule not found at {path}; keeping ALL reactions. "
              f"Run optimize_thermo_consensus.py first.")
        return lambda df: np.ones(len(df), bool)
    obj = json.loads(path.read_text())
    theta = {k: float(v) for k, v in obj["thresholds"].items()}
    drop = set(obj.get("dropped_classes", []))

    def predicate(df: pd.DataFrame) -> np.ndarray:
        d = df.copy()
        if "abs_dg_eq" not in d and "dg_eq" in d:
            d["abs_dg_eq"] = d["dg_eq"].abs()
        if "abs_net_proton" not in d and "net_proton" in d:
            d["abs_net_proton"] = d["net_proton"].abs()
        missing = [f for f in theta if f not in d.columns]
        if missing:
            raise KeyError(f"selector needs columns {missing}")
        return _apply(d, theta, drop)

    return predicate


# -------------------------------------------------------------------- main
def main() -> None:
    df = pd.read_csv(KEY, sep="\t", low_memory=False)
    df["abs_dg_eq"] = df["dg_eq"].abs()
    df["abs_net_proton"] = df["net_proton"].abs()
    x = df["dg_eq"].to_numpy(float)
    y = df["dg_dgp"].to_numpy(float)
    strata = pd.qcut(df["abs_dg_eq"], 10, labels=False, duplicates="drop").to_numpy()
    print(f"pool: {len(df)} reactions;  unfiltered {metrics(x, y)}")

    # ---- oracle bound
    bars = [1, 1.5, 2, 2.5, 3, 4, 5, 7.5, 10, 15, 20]
    orc = oracle_frontier(x, y, bars)
    print("\n=== ORACLE (exact bound; selects on the outcome, so not a deliverable) ===")
    print(orc[["rmse_bar", "n", "coverage", "ccc", "r", "slope", "median_absdiff"]]
          .to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    # ---- rule frontier
    print("\n=== RULE (generalisable) — fitting across the frontier ===")
    rule_rows, fitted = [], {}
    for E in bars:
        tgt = dict(TARGET); tgt["rmse"] = float(E)
        got = fit_rule(df, tgt, x, y, strata, verbose=False)
        if got is None:
            print(f"  RMSE<={E:<5g} NO FEASIBLE RULE at CCC>={tgt['ccc']}, "
                  f"|slope-1|<={tgt['slope_tol']}")
            rule_rows.append({"solver": "rule", "rmse_bar": E, "n": 0,
                              "coverage": 0.0, "feasible": False})
            continue
        th, dr = got
        m = _apply(df, th, dr)
        mm = metrics(x[m], y[m])
        ret = stratum_retention(m, strata)
        ok = feasible(mm, tgt) and (tgt["stratum_floor"] <= 0
                                    or ret.min() >= tgt["stratum_floor"])
        assert ok, f"fit_rule returned an infeasible solution at RMSE<={E}"
        rule_rows.append({"solver": "rule", "rmse_bar": E, **mm,
                          "coverage": m.mean(), "feasible": True,
                          "min_decile_retention": float(ret.min()),
                          "top_decile_retention": float(ret[-1])})
        fitted[E] = (th, dr)
        orc_n = orc.loc[orc.rmse_bar == E, "n"]
        rec = (mm["n"] / orc_n.iat[0]) if len(orc_n) and orc_n.iat[0] else np.nan
        print(f"  RMSE<={E:<5g} n={mm['n']:6d} ({m.mean():5.1%})  CCC={mm['ccc']:.3f}  "
              f"slope={mm['slope']:.2f}  med|d|={mm['median_absdiff']:5.2f}  "
              f"recovers {rec:5.1%} of oracle")

    # What does insisting on the stratum floor cost? Sweep it separately so the
    # trade-off is visible instead of buried in an infeasibility.
    print(f"\n=== cost of forcing a {STRATUM_FLOOR_SWEEP:.0%} floor in EVERY |dG| decile ===")
    floor_rows = []
    for E in bars:
        tgt = dict(TARGET); tgt["rmse"] = float(E)
        tgt["stratum_floor"] = STRATUM_FLOOR_SWEEP
        got = fit_rule(df, tgt, x, y, strata, verbose=False, de_refine=False)
        if got is None:
            print(f"  RMSE<={E:<5g} infeasible with the floor")
            floor_rows.append({"solver": "rule_floored", "rmse_bar": E,
                               "n": 0, "coverage": 0.0, "feasible": False})
            continue
        m = _apply(df, *got); mm = metrics(x[m], y[m])
        base = next((r["n"] for r in rule_rows if r["rmse_bar"] == E), None)
        floor_rows.append({"solver": "rule_floored", "rmse_bar": E, **mm,
                           "coverage": m.mean(), "feasible": True})
        print(f"  RMSE<={E:<5g} n={mm['n']:6d} ({m.mean():5.1%})"
              + (f"   vs {base} unfloored" if base else ""))

    frontier = pd.concat([orc, pd.DataFrame(rule_rows), pd.DataFrame(floor_rows)],
                         ignore_index=True)
    frontier.to_csv(DATA / "consensus_frontier.tsv", sep="\t", index=False,
                    float_format="%.5f")

    # ---- ship the target operating point
    if TARGET["rmse"] not in fitted:
        raise SystemExit(f"no feasible rule at the target RMSE <= {TARGET['rmse']}; "
                         f"loosen TARGET rather than shipping a rule that misses it")
    th, dr = fitted[TARGET["rmse"]]
    keep = _apply(df, th, dr)
    mm = metrics(x[keep], y[keep])
    ret = stratum_retention(keep, strata)
    print(f"\n=== OPERATING POINT (RMSE <= {TARGET['rmse']}) ===")
    print(f"  thresholds: " + ", ".join(f"{k}<={v:.3g}" for k, v in th.items()
                                        if v < np.nanmax(df[k].to_numpy(float))))
    print(f"  dropped classes: {sorted(dr) if dr else 'none'}")
    print(f"  n={mm['n']} ({keep.mean():.1%})  CCC={mm['ccc']:.4f}  r={mm['r']:.4f}  "
          f"slope={mm['slope']:.3f}  RMSE={mm['rmse']:.3f}  median|d|={mm['median_absdiff']:.3f}")
    print("  per-|dG|-decile retention (the solution's shape):")
    for b, r in enumerate(ret):
        lo = df.loc[strata == b, "abs_dg_eq"].min()
        hi = df.loc[strata == b, "abs_dg_eq"].max()
        bar = "#" * int(round(r * 40))
        print(f"    decile {b}  |dG| {lo:7.2f}-{hi:8.2f}  {r:6.1%} {bar}")

    print("\n=== 5-fold cross-validation (thresholds refitted per fold) ===")
    cv = cross_validate(df, TARGET, x, y, strata)
    print(cv[["fold", "train_coverage", "train_ccc", "train_rmse",
              "test_coverage", "test_ccc", "test_rmse"]]
          .to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print(f"  mean held-out: coverage {cv.test_coverage.mean():.1%}  "
          f"CCC {cv.test_ccc.mean():.4f}  RMSE {cv.test_rmse.mean():.3f}")

    RULE_JSON.write_text(json.dumps({
        "target": TARGET,
        "thresholds": {k: float(v) for k, v in th.items()},
        "dropped_classes": sorted(dr),
        "in_sample": mm,
        "cv_mean_test": {"coverage": float(cv.test_coverage.mean()),
                         "ccc": float(cv.test_ccc.mean()),
                         "rmse": float(cv.test_rmse.mean())},
        "note": "Selects on agreement between two estimators, not on correctness.",
    }, indent=1))
    out = df.loc[keep, ["rxn", "name", "ec", "chem_class", "dg_eq", "dg_dgp",
                        "dgp_uncertainty", "eq_uncertainty"]].copy()
    out["absdiff"] = (out.dg_eq - out.dg_dgp).abs()
    out.sort_values("absdiff").to_csv(DATA / "consensus_selected.tsv", sep="\t",
                                      index=False, float_format="%.4f")
    print(f"\nwrote {RULE_JSON.name}, consensus_frontier.tsv, consensus_selected.tsv")


if __name__ == "__main__":
    main()
