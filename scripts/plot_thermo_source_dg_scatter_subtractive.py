#!/usr/bin/env python3
"""The error-filtered ΔG′° scatter plots, run SUBTRACTIVELY: at each threshold T
the CONFIDENT reactions are removed and what is left is drawn.

plot_thermo_source_dg_scatter_filtered.py keeps the reactions both sources are
confident about (both σ ≤ T, or both ê ≤ T) and asks "how well do the sources
agree once we only trust the confident ones?". This script asks the mirror-image
question -- "what is left over, and what does it look like?" -- by deleting those
same reactions at T ∈ {0.5, 1, 2, 5, 10, 20} kcal/mol and plotting the residue.
Same three pairs, same two criteria (σ and ê), same pair of ellipses (95%
concentration + minimum-volume enclosing), same reversibility-transition
palette. Everything except the mask is inherited from that script by import, so
the two families of figures cannot drift apart.

TWO REMOVAL RULES
-----------------
any   (default; the exact complement of the additive plots)
      Remove reaction i iff max(crit_a, crit_b) ≤ T -- exactly the set the
      additive panel at T draws. Keep iff AT LEAST ONE source admits more than
      T. n_keep here + n_keep in the additive panel = the pair's base set, so the
      two figures partition the data and the counts add up.

both  (the strict residue)
      Remove reaction i iff min(crit_a, crit_b) ≤ T -- i.e. one confident source
      is enough to retire the reaction. Keep iff BOTH sources admit more than T:
      the reactions where nobody claims to know the answer. This empties out much
      faster than `any` and is the interesting one at large T.

NaN CRITERION = "MORE THAN T"
-----------------------------
ê is NaN where the assignment model refuses the source outright (eQuilibrator
sentinels are already dropped at load; dGPredictor-ModelSEED on the
quinone/quinol couple). In the additive script NaN fails the filter and the
reaction is excluded. Here NaN is treated as exceeding every T -- a refused
source is the opposite of a confident one -- so a NaN reaction is never removed
under `any`, and under `both` it goes only when its PARTNER source is confident.
That is what makes `any` an exact complement. Each panel prints how many of its points
are riding on a NaN criterion, and slice_counts_subtractive.tsv carries the count
as n_nan_crit, because at large T under the `ehat` criterion the survivors are
increasingly just the refusals.

OUTPUTS  reports/thermoComparison/figures/thermo_source_dg_scatter_subtractive/
    grid_unfiltered.png                       the 3 pairs, nothing removed (reference)
    grid_<rule>_<crit>_<mode>.png             3 pairs x 6 thresholds, 8 combinations
    <rule>/<crit>/<mode>/dg_scatter_<a>_vs_<b>_<crit>_gt<T>.png   144 single panels
    slice_counts_subtractive.tsv              n kept / n removed / r / median |Δ| per slice
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_thermo_source_dg_scatter_filtered import (  # noqa: E402
    ANALYSIS_DIR, CRIT_COLUMN, CRIT_SHORT, CRITERIA, EXTREME_CUTOFF, INK_MUTED,
    INK_PRIMARY, INK_SECONDARY, MODES, MSDB_DATA, PAIRS, SOURCES, STAT_KEYS,
    SURFACE, THRESHOLDS, attach_ehat, base_set, draw_panel, ellipse_note,
    load_table,
)

OUT_DIR = (ANALYSIS_DIR / "reports" / "thermoComparison" / "figures"
           / "thermo_source_dg_scatter_subtractive")
# Re-parsing every reaction JSON and re-running the cascade three times over
# costs ~2.5 minutes, which is most of this script's runtime and all of the cost
# of any follow-up question about the slices. --cache parks the assembled table
# so those are seconds. Opt-in, and NOT invalidated automatically: delete it (or
# run without --cache) after MSDB_DATA or the cascade changes.
CACHE_PATH = ANALYSIS_DIR / "results" / "thermo_scatter_cache" / "source_table.pkl"

RULES = ["any", "both"]
RULE_LABEL = {
    "any": ("removing every reaction both sources are confident about "
            "(the exact complement of the additive panels)"),
    "both": ("removing every reaction EITHER source is confident about — "
             "only reactions no source vouches for survive"),
}
RULE_SHORT = {
    "any":  "max({a}, {b}) > T",
    "both": "min({a}, {b}) > T",
}
CRIT_LABEL_SUB = {
    "sigma": "source-reported σ (predicted error)",
    "ehat":  "calibrated ê (expected error vs TECRDB)",
}


def surviving_mask(sub: pd.DataFrame, a: str, b: str, crit: str, tol: float,
                   rule: str) -> np.ndarray:
    """True for the reactions that SURVIVE the subtraction at tolerance ``tol``.

    NaN propagates as "exceeds T": ``NaN <= tol`` is False, so a refused source
    never counts as confident. Under `any` that makes the reaction unremovable;
    under `both` its confident partner can still retire it.
    """
    col = CRIT_COLUMN[crit]
    ea = sub[f"{col}_{a}"].to_numpy(float)
    eb = sub[f"{col}_{b}"].to_numpy(float)
    with np.errstate(invalid="ignore"):
        confident_a, confident_b = ea <= tol, eb <= tol
    removed = (confident_a & confident_b) if rule == "any" else (confident_a | confident_b)
    return ~removed


def nan_count(sub: pd.DataFrame, a: str, b: str, crit: str, keep: np.ndarray) -> int:
    col = CRIT_COLUMN[crit]
    bad = (~np.isfinite(sub[f"{col}_{a}"].to_numpy(float))
           | ~np.isfinite(sub[f"{col}_{b}"].to_numpy(float)))
    return int((keep & bad).sum())


def pair_ceiling(sub: pd.DataFrame, a: str, b: str, crit: str, rule: str) -> float:
    """Smallest T that empties the panel: the largest value the removal rule has
    to clear. Per reaction that is max(crit_a, crit_b) under `any` and
    min(crit_a, crit_b) under `both`; the ceiling is the largest of those over
    the pair. NaN is +inf here for the same reason it is elsewhere -- a refused
    source is never confident -- so under `any` one NaN makes the reaction
    unremovable (ceiling inf), while under `both` its finite partner still
    decides."""
    col = CRIT_COLUMN[crit]
    ea = sub[f"{col}_{a}"].to_numpy(float)
    eb = sub[f"{col}_{b}"].to_numpy(float)
    ea = np.where(np.isnan(ea), np.inf, ea)
    eb = np.where(np.isnan(eb), np.inf, eb)
    worst = np.maximum(ea, eb) if rule == "any" else np.minimum(ea, eb)
    return float(worst.max()) if worst.size else float("nan")


def load_or_cache(cache: Path | None) -> pd.DataFrame:
    """The dg / σ / ê / cascade-operator table, optionally read from or written
    to a pickle. See CACHE_PATH on invalidation: there is none."""
    if cache is not None and cache.exists():
        print(f"reusing cached source table {cache}")
        return pd.read_pickle(cache)
    print(f"loading reactions from {MSDB_DATA} ...")
    df = attach_ehat(load_table())
    if cache is not None:
        cache.parent.mkdir(parents=True, exist_ok=True)
        df.to_pickle(cache)
        print(f"cached source table to {cache}")
    return df


def sub_title(a: str, b: str, crit: str, tol: float | None, rule: str, st: dict,
              n_base: int, compact: bool = False) -> str:
    head = f"{SOURCES[a][1]} vs {SOURCES[b][1]}"
    if tol is None:
        return (f"{head}\nnothing removed · n = {st['n_keep']:,} · "
                f"r = {st['r']:.2f} · median |Δ| = {st['median_absdiff']:.2f} kcal/mol")
    frac = f"n = {st['n_keep']:,} of {n_base:,} ({st['n_keep']/max(n_base,1):.1%} kept)"
    stats = f"r = {st['r']:.2f} · median |Δ| = {st['median_absdiff']:.2f} kcal/mol"
    cut = f"{CRIT_SHORT[crit]} ≤ {tol:g} removed"
    if compact:
        return f"{cut}\n{frac}\n{stats}"
    return f"{head} · {cut}\n{frac} · {stats}"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--rules", nargs="+", choices=RULES, default=RULES)
    ap.add_argument("--no-singles", action="store_true",
                    help="write only the grid figures, skip the 144 single panels")
    ap.add_argument("--cache", nargs="?", const=str(CACHE_PATH), default=None,
                    metavar="PATH",
                    help="read the source table from this pickle if it exists, "
                         "write it there otherwise (no staleness check — delete "
                         "it after the data or the cascade changes)")
    args = ap.parse_args()
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    df = load_or_cache(Path(args.cache) if args.cache else None)

    bases, n_extreme = {}, {}
    for a, b in PAIRS:
        bases[(a, b)], n_extreme[(a, b)] = base_set(df, a, b)
        print(f"  {SOURCES[a][1]} ∩ {SOURCES[b][1]}: {len(bases[(a,b)]):,} reactions "
              f"({n_extreme[(a,b)]} dropped as |ΔG′°| > {EXTREME_CUTOFF:g})")

    records = []

    # ---- reference: nothing removed -----------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(21, 7.0), dpi=130)
    fig.patch.set_facecolor(SURFACE)
    for ax, (a, b) in zip(axes, PAIRS):
        sub = bases[(a, b)]
        st = draw_panel(ax, sub, a, b, np.ones(len(sub), dtype=bool), "filtered")
        ax.set_title(sub_title(a, b, "none", None, "any", st, len(sub)),
                     color=INK_PRIMARY, fontsize=11.5, pad=10)
        ax.set_xlabel(SOURCES[a][2], color=INK_SECONDARY, fontsize=10)
        ax.set_ylabel(SOURCES[b][2], color=INK_SECONDARY, fontsize=10)
        records.append({"pair": f"{a}_vs_{b}", "rule": "none", "criterion": "none",
                        "tolerance": np.nan, "n_base": len(sub), "n_nan_crit": 0,
                        **{k: st[k] for k in STAT_KEYS}})
    fig.suptitle("ΔG′° agreement between thermodynamic sources — nothing removed\n"
                 "filled ellipse = 95% concentration (empirical Mahalanobis quantile; "
                 "major axis = orthogonal fit) · dashed = minimum-volume enclosing "
                 "ellipse (Löwner–John hull proxy)",
                 color=INK_PRIMARY, fontsize=12.5, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    p = out_dir / "grid_unfiltered.png"
    fig.savefig(p, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {p}")

    # ---- the threshold grids -------------------------------------------------
    for rule in args.rules:
        for crit in CRITERIA:
            for mode in MODES:
                excl = f"Removed ({CRIT_SHORT[crit]} ≤ T)"
                empty = "every reaction was\nremoved at this threshold"
                fig, axes = plt.subplots(len(PAIRS), len(THRESHOLDS),
                                         figsize=(4.45 * len(THRESHOLDS),
                                                  4.95 * len(PAIRS)),
                                         dpi=120)
                fig.patch.set_facecolor(SURFACE)
                for ri, (a, b) in enumerate(PAIRS):
                    sub = bases[(a, b)]
                    row_lim = None
                    if mode == "context":
                        pts = np.concatenate([sub[f"dg_{a}"].to_numpy(float),
                                              sub[f"dg_{b}"].to_numpy(float)])
                        lo, hi = float(pts.min()), float(pts.max())
                        pad = 0.09 * (hi - lo)
                        row_lim = (lo - pad, hi + pad)
                    for ci, tol in enumerate(THRESHOLDS):
                        ax = axes[ri][ci]
                        keep = surviving_mask(sub, a, b, crit, tol, rule)
                        st = draw_panel(ax, sub, a, b, keep, mode,
                                        xlim=row_lim, ylim=row_lim,
                                        show_legend=True, compact=True,
                                        excluded_label=excl, empty_note=empty)
                        ax.set_title(sub_title(a, b, crit, tol, rule, st, len(sub),
                                               compact=True),
                                     color=INK_PRIMARY, fontsize=9, pad=6)
                        if ci == 0:
                            ax.text(-0.30, 0.5,
                                    f"{SOURCES[a][1]}\nvs\n{SOURCES[b][1]}",
                                    transform=ax.transAxes, ha="center", va="center",
                                    rotation=90, fontsize=11, color=INK_PRIMARY)
                        if ri == len(PAIRS) - 1:
                            ax.set_xlabel(SOURCES[a][2], color=INK_SECONDARY, fontsize=9)
                        if ci == 0:
                            ax.set_ylabel(SOURCES[b][2], color=INK_SECONDARY, fontsize=9)
                        if mode == "filtered":
                            records.append({
                                "pair": f"{a}_vs_{b}", "rule": rule, "criterion": crit,
                                "tolerance": tol, "n_base": len(sub),
                                "n_nan_crit": nan_count(sub, a, b, crit, keep),
                                **{k: st[k] for k in STAT_KEYS}})
                sub_note = ("only the surviving reactions are drawn; axes rescale to them"
                            if mode == "filtered" else
                            "every reaction is drawn — removed reactions grayed out, "
                            "survivors in color and inside the ellipses; axes shared per row")
                fig.suptitle(
                    f"ΔG′° source agreement, subtractive on {CRIT_LABEL_SUB[crit]} — "
                    f"{RULE_LABEL[rule]}\n{sub_note}",
                    color=INK_PRIMARY, fontsize=15, y=0.995)
                fig.tight_layout(rect=(0.022, 0, 1, 0.95))
                p = out_dir / f"grid_{rule}_{crit}_{mode}.png"
                fig.savefig(p, facecolor=SURFACE)
                plt.close(fig)
                print(f"wrote {p}")

                if args.no_singles:
                    continue
                single_dir = out_dir / rule / crit / mode
                single_dir.mkdir(parents=True, exist_ok=True)
                for a, b in PAIRS:
                    sub = bases[(a, b)]
                    row_lim = None
                    if mode == "context":
                        pts = np.concatenate([sub[f"dg_{a}"].to_numpy(float),
                                              sub[f"dg_{b}"].to_numpy(float)])
                        lo, hi = float(pts.min()), float(pts.max())
                        pad = 0.09 * (hi - lo)
                        row_lim = (lo - pad, hi + pad)
                    for tol in THRESHOLDS:
                        fig, ax = plt.subplots(figsize=(7.2, 6.6), dpi=150)
                        fig.patch.set_facecolor(SURFACE)
                        keep = surviving_mask(sub, a, b, crit, tol, rule)
                        st = draw_panel(ax, sub, a, b, keep, mode,
                                        xlim=row_lim, ylim=row_lim,
                                        excluded_label=excl, empty_note=empty)
                        ax.set_title(sub_title(a, b, crit, tol, rule, st, len(sub)),
                                     color=INK_PRIMARY, fontsize=11.5, pad=10)
                        ax.set_xlabel(SOURCES[a][2], color=INK_SECONDARY, fontsize=10)
                        ax.set_ylabel(SOURCES[b][2], color=INK_SECONDARY, fontsize=10)
                        keeps = RULE_SHORT[rule].format(a=CRIT_SHORT[crit] + "_A",
                                                        b=CRIT_SHORT[crit] + "_B")
                        where = ("BOTH sources" if rule == "any" else "EITHER source")
                        note = (f"subtraction: every reaction with "
                                f"{CRIT_SHORT[crit]} ≤ {tol:g} kcal/mol in {where} is "
                                f"removed — survivors satisfy {keeps}\n"
                                + ellipse_note(st))
                        n_nan = nan_count(sub, a, b, crit, keep)
                        if n_nan:
                            note += (f"\n{n_nan:,} survivor(s) have no "
                                     f"{CRIT_SHORT[crit]} at all (source refused); "
                                     f"they are counted as exceeding every threshold")
                        if n_extreme[(a, b)]:
                            note += (f"\n{n_extreme[(a,b)]} reaction(s) with |ΔG′°| > "
                                     f"{EXTREME_CUTOFF:g} kcal/mol excluded from the "
                                     f"base set")
                        if st["n_keep"] == 0:
                            ceil = pair_ceiling(sub, a, b, crit, rule)
                            note += (f"\nempty because the threshold clears every "
                                     f"reaction; the last to go is at "
                                     f"{ceil:.2f} kcal/mol")
                        ax.text(0.5, -0.125, note,
                                transform=ax.transAxes, ha="center", va="top",
                                fontsize=7.5, color=INK_MUTED, linespacing=1.6)
                        fig.tight_layout()
                        fig.subplots_adjust(bottom=0.135 + 0.028 * (note.count("\n") + 1))
                        fig.savefig(
                            single_dir / f"dg_scatter_{a}_vs_{b}_{crit}_gt{tol:g}.png",
                            facecolor=SURFACE)
                        plt.close(fig)
                print(f"  wrote {len(PAIRS) * len(THRESHOLDS)} single panels to {single_dir}")

    counts = pd.DataFrame(records)
    counts_path = out_dir / "slice_counts_subtractive.tsv"
    counts.to_csv(counts_path, sep="\t", index=False, float_format="%.4f")
    print(f"wrote {counts_path}")
    print(counts.to_string(index=False))


if __name__ == "__main__":
    main()
