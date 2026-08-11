#!/usr/bin/env python3
"""Pairwise DeltaG scatter plots across the three thermodynamic sources, sliced by
how much error each source ADMITS to on the reaction, with a covering oval.

This is the error-filtered companion to plot_thermo_source_dg_scatter.py. Same
three pairs -- Group Contribution, eQuilibrator, dGPredictor-ModelSEED (the
retrained, ModelSEED-finetuned variant on dev, NOT the original KEGG-mapped
dGPredictor) -- but each pair is redrawn six times, keeping only the reactions
where both sources are confident to within a tolerance T in
{0.5, 1, 2, 5, 10, 20} kcal/mol.

TWO NOTIONS OF "ERROR", PLOTTED SEPARATELY
------------------------------------------
sigma  ("predicted error", criterion 1)
    The uncertainty each source writes into its OWN record --
    ``thermodynamics[label] = [dg, dge, operator]``, second field. This is the
    source's self-reported error bar, on its own scale, never calibrated
    against anything. Keep reaction i iff max(sigma_a(i), sigma_b(i)) <= T.

ehat   ("tolerance", criterion 2)
    The CALIBRATED expected |error vs TECRDB| from
    optimize_thermo_source_assignment.py: each source's sigma pushed through an
    isotonic regression fitted on 802 stereo-exact TECRDB matches, so the number
    is in real kcal/mol-of-actual-error rather than in that source's private
    units. Read from results/eq_vs_dgpms/source_assignment.tsv. Keep reaction i
    iff max(ehat_a(i), ehat_b(i)) <= T.
    ehat is NaN -- and the reaction therefore always fails -- where the
    assignment model refuses the source outright: eQuilibrator sentinel
    uncertainties (sigma > 100 kcal/mol, the source disclaiming an estimate) and
    dGPredictor-ModelSEED on the quinone/quinol couple (52.8% sign-wrong).

The two are NOT interchangeable and the same T means different things to each:
sigma <= 2 is "the source says 2"; ehat <= 2 is "the source has been measured to
be off by about 2".

TWO RENDERINGS PER SLICE
------------------------
filtered  only the reactions that pass are drawn, axes rescaled to them.
context   every reaction in the pair's base set is drawn; the ones that fail the
          threshold are grayed out and pushed behind, the ones that pass keep
          their categorical color and get the oval. Axes are shared down a
          column-group so the shrinking of the passing set is visible as an
          actual shrinking, not hidden by autoscale.

THE OVAL
--------
Minimum-volume enclosing ellipse (Khachiyan's algorithm) over the points that
pass the threshold -- the tightest oval that genuinely covers all of them, not a
confidence ellipse. Its area is the honest picture of how much of the ΔG'° plane
a "confident" subset still occupies.

Point color is the reversibility transition between the two sources, identical
in definition and palette to plot_thermo_source_dg_scatter.py: each source's own
ΔG'° run through the unmodified cascade (DEFAULT_HEURISTICS via
per_source_energy), collapsed to reversible ("=") vs irreversible (">"/"<").

OUTPUTS  reports/thermoComparison/figures/thermo_source_dg_scatter_filtered/
    grid_unfiltered.png                   the 3 pairs, no threshold, with oval
    grid_<crit>_<mode>.png                3 pairs x 6 thresholds, 4 combinations
    <crit>/<mode>/dg_scatter_<a>_vs_<b>_<crit>_le<T>.png      72 single panels
    slice_counts.tsv                      n passing / n excluded / r per slice
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Ellipse

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
# Data comes from the dev snapshot: the only checkout carrying the retrained
# dGPredictor-ModelSEED energies. Code (the cascade) comes from the working
# ModelSEEDDatabase checkout.
MSDB_DATA = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/tmp/devsnap"))
MSDB_CODE = Path(os.environ.get("MSDB_CODE_ROOT", "/scratch/ctaylor/ModelSEEDDatabase"))
ASSIGN_TSV = Path(os.environ.get(
    "EQDGP_ASSIGNMENT",
    str(ANALYSIS_DIR / "results" / "eq_vs_dgpms" / "source_assignment.tsv")))
OUT_DIR = (ANALYSIS_DIR / "reports" / "thermoComparison" / "figures"
           / "thermo_source_dg_scatter_filtered")

sys.path.insert(0, str(MSDB_CODE / "Scripts" / "Thermodynamics"))
from reversibility_heuristics import (  # noqa: E402
    DEFAULT_HEURISTICS, run_reversibility, per_source_energy,
)

# key -> (thermodynamics subkey, short label, axis title)
SOURCES = {
    "GC": ("Group contribution", "Group contribution",
           "Group Contribution ΔG′° (kcal/mol)"),
    "EQ": ("eQuilibrator", "eQuilibrator",
           "eQuilibrator (2.0) ΔG′° (kcal/mol)"),
    "DGPMS": ("dGPredictor-ModelSEED", "dGPredictor-ModelSEED",
              "dGPredictor-ModelSEED ΔG′° (kcal/mol)"),
}
PAIRS = [("GC", "EQ"), ("GC", "DGPMS"), ("EQ", "DGPMS")]
THRESHOLDS = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
CRITERIA = ["sigma", "ehat"]
MODES = ["filtered", "context"]

CRIT_LABEL = {
    "sigma": "source-reported σ (predicted error)",
    "ehat":  "calibrated ê (expected error vs TECRDB)",
}
CRIT_SHORT = {"sigma": "σ", "ehat": "ê"}
CRIT_COLUMN = {"sigma": "sig", "ehat": "ehat"}    # column prefix in the table

EQ_SENTINEL = 100.0     # kcal/mol; eQuilibrator's "I have no estimate" marker
# A handful of aggregate/polymer reactions carry genuine but chemically
# implausible magnitudes (rxn05017, Group Contribution ~15,900 kcal/mol). They
# are dropped from the base set here rather than merely clipped, because a
# single such point would set the oval on its own. Count is reported per panel.
EXTREME_CUTOFF = 1500.0

# Palette carried over unchanged from plot_thermo_source_dg_scatter.py: three
# validated categorical hues plus a neutral gray for "No change". Threshold-
# excluded points use a much lighter neutral, drawn smaller and behind, so the
# two grays are separated by lightness and size rather than hue.
CATEGORY_ORDER = [
    "No change",
    "Reversible → Irreversible",
    "Irreversible → Reversible",
    "Irreversible → Irreversible",
]
CATEGORY_COLOR = {
    "No change": "#9c9a94",
    "Reversible → Irreversible": "#2a78d6",
    "Irreversible → Reversible": "#eb6834",
    "Irreversible → Irreversible": "#1baf7a",
}
EXCLUDED_COLOR = "#dedbd1"
OVAL_COLOR = "#5a3fb0"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
BASELINE = "#c3c2b7"
SURFACE = "#fcfcfb"


# ------------------------------------------------------------------ geometry
def mvee(points: np.ndarray, tol: float = 1e-4, max_iter: int = 5000):
    """Minimum-volume enclosing ellipse (Khachiyan). Returns (center, width,
    height, angle_deg) for matplotlib, or None if degenerate.

    Points are isotropically rescaled before the solve -- isotropic scaling maps
    ellipses to ellipses and preserves which one is minimal, so this only buys
    conditioning.
    """
    P = np.asarray(points, dtype=float)
    P = P[np.isfinite(P).all(axis=1)]
    if len(P) < 3:
        return None
    scale = float(np.abs(P).max()) or 1.0
    P = P / scale
    n, d = P.shape
    Q = np.column_stack([P, np.ones(n)]).T          # (d+1, n)
    u = np.full(n, 1.0 / n)
    for _ in range(max_iter):
        X = (Q * u) @ Q.T
        try:
            Xinv = np.linalg.inv(X)
        except np.linalg.LinAlgError:
            return None
        M = np.einsum("ij,ij->j", Q, Xinv @ Q)
        j = int(np.argmax(M))
        maxM = M[j]
        step = (maxM - d - 1.0) / ((d + 1.0) * (maxM - 1.0))
        if not np.isfinite(step) or step <= 0:
            break
        new_u = (1.0 - step) * u
        new_u[j] += step
        if np.linalg.norm(new_u - u) < tol:
            u = new_u
            break
        u = new_u
    c = P.T @ u
    try:
        A = np.linalg.inv((P.T * u) @ P - np.outer(c, c)) / d
    except np.linalg.LinAlgError:
        return None
    evals, evecs = np.linalg.eigh(A)
    if np.any(evals <= 0) or not np.all(np.isfinite(evals)):
        return None
    axes = 1.0 / np.sqrt(evals)                      # semi-axis lengths
    order = np.argsort(axes)[::-1]
    axes, evecs = axes[order], evecs[:, order]
    angle = float(np.degrees(np.arctan2(evecs[1, 0], evecs[0, 0])))
    return (c * scale, 2.0 * axes[0] * scale, 2.0 * axes[1] * scale, angle)


def draw_oval(ax, points: np.ndarray, label: str | None = None):
    ell = mvee(points)
    if ell is None:
        return None
    c, w, h, ang = ell
    patch = Ellipse(tuple(c), w, h, angle=ang, facecolor=OVAL_COLOR, alpha=0.055,
                    edgecolor=OVAL_COLOR, linewidth=1.4, linestyle="-",
                    zorder=4, label=label)
    ax.add_patch(patch)
    return ell


# ------------------------------------------------------------------- loading
def classify(op_a: str, op_b: str) -> str:
    rev_a, rev_b = op_a == "=", op_b == "="
    if rev_a and rev_b:
        return "No change"
    if rev_a and not rev_b:
        return "Reversible → Irreversible"
    if not rev_a and rev_b:
        return "Irreversible → Reversible"
    return "Irreversible → Irreversible"


def load_table() -> pd.DataFrame:
    """One row per non-EMPTY reaction: dg / sigma / cascade operator per source."""
    rows = []
    n_sentinel = 0
    for path in sorted(glob.glob(str(MSDB_DATA / "Biochemistry" / "reaction_*.json"))):
        for entry in json.load(open(path)):
            if entry.get("status") == "EMPTY":
                continue
            thermo = entry.get("thermodynamics") or {}
            row = {"rxn": entry["id"], "name": entry.get("name", "")}
            for key, (subkey, _, _) in SOURCES.items():
                triple = thermo.get(subkey)
                if not triple or len(triple) < 3 or triple[2] in (None, "?"):
                    continue
                try:
                    dg, sig = float(triple[0]), abs(float(triple[1]))
                except (TypeError, ValueError):
                    continue
                if key == "EQ" and sig > EQ_SENTINEL:
                    # the source itself declares no estimate -- not a data point
                    n_sentinel += 1
                    continue
                _, op, _ = run_reversibility(entry, per_source_energy(subkey),
                                             DEFAULT_HEURISTICS)
                if op is None:
                    continue
                row[f"dg_{key}"] = dg
                row[f"sig_{key}"] = sig
                row[f"op_{key}"] = op
            rows.append(row)
    df = pd.DataFrame(rows).set_index("rxn")
    print(f"  {len(df):,} non-EMPTY reactions; dropped {n_sentinel:,} eQuilibrator "
          f"sentinel value(s) (σ > {EQ_SENTINEL:g})")
    for key, (_, label, _) in SOURCES.items():
        print(f"    {label:24s} {df.get(f'dg_{key}', pd.Series(dtype=float)).notna().sum():,} usable ΔG′°")
    return df


def attach_ehat(df: pd.DataFrame) -> pd.DataFrame:
    a = pd.read_csv(ASSIGN_TSV, sep="\t").set_index("rxn")
    for key in SOURCES:
        df[f"ehat_{key}"] = a[f"ehat_{key}"].reindex(df.index)
    n = {k: int(df[f"ehat_{k}"].notna().sum()) for k in SOURCES}
    print(f"  ê attached from {ASSIGN_TSV.name}: " +
          ", ".join(f"{SOURCES[k][1]} {v:,}" for k, v in n.items()))
    return df


def base_set(df: pd.DataFrame, a: str, b: str) -> pd.DataFrame:
    """Reactions both sources cover, with implausible magnitudes removed."""
    m = df[f"dg_{a}"].notna() & df[f"dg_{b}"].notna()
    sub = df[m]
    keep = (sub[f"dg_{a}"].abs() <= EXTREME_CUTOFF) & (sub[f"dg_{b}"].abs() <= EXTREME_CUTOFF)
    return sub[keep], int((~keep).sum())


def passing_mask(sub: pd.DataFrame, a: str, b: str, crit: str, tol: float) -> np.ndarray:
    col = CRIT_COLUMN[crit]
    ea = sub[f"{col}_{a}"].to_numpy(float)
    eb = sub[f"{col}_{b}"].to_numpy(float)
    with np.errstate(invalid="ignore"):
        return (ea <= tol) & (eb <= tol)     # NaN ê (refused source) fails


# ------------------------------------------------------------------ plotting
def draw_panel(ax, sub: pd.DataFrame, a: str, b: str, keep: np.ndarray,
               mode: str, xlim=None, ylim=None, show_legend: bool = True,
               compact: bool = False) -> dict:
    xs = sub[f"dg_{a}"].to_numpy(float)
    ys = sub[f"dg_{b}"].to_numpy(float)
    cats = np.array([classify(oa, ob)
                     for oa, ob in zip(sub[f"op_{a}"], sub[f"op_{b}"])], dtype=object)

    n_keep = int(keep.sum())
    n_drop = int((~keep).sum())
    xk, yk = xs[keep], ys[keep]
    r = float(np.corrcoef(xk, yk)[0, 1]) if n_keep > 1 and xk.std() and yk.std() else float("nan")
    med = float(np.median(np.abs(xk - yk))) if n_keep else float("nan")
    # Reactions both sources put at exactly 0 kcal/mol -- transport, some
    # isomerisations, anything whose net group change cancels. They are real
    # entries, but they are also free agreement: they inflate r and pull the
    # median |delta| to 0, and low-sigma slices are heavily enriched for them.
    n_zero = int(((xk == 0.0) & (yk == 0.0)).sum()) if n_keep else 0
    frac_zero = n_zero / n_keep if n_keep else float("nan")

    ax.set_facecolor(SURFACE)
    if xlim is None or ylim is None:
        pts = np.concatenate([xs, ys]) if mode == "context" else np.concatenate([xk, yk])
        if pts.size == 0:
            pts = np.array([-1.0, 1.0])
        lo, hi = float(pts.min()), float(pts.max())
        if hi - lo < 1e-9:
            lo, hi = lo - 1.0, hi + 1.0
        pad = 0.09 * (hi - lo)
        xlim = ylim = (lo - pad, hi + pad)

    ax.plot(xlim, xlim, linestyle="--", linewidth=1.1, color=BASELINE, zorder=1)

    counts = Counter(cats[keep])
    if mode == "context" and n_drop:
        ax.scatter(xs[~keep], ys[~keep], s=5, linewidths=0, color=EXCLUDED_COLOR,
                   alpha=0.55, zorder=2,
                   label=f"Excluded by threshold ({n_drop:,})")
    for cat in CATEGORY_ORDER:
        idx = keep & (cats == cat)
        if not idx.any():
            continue
        ax.scatter(xs[idx], ys[idx], s=13 if not compact else 9, linewidths=0,
                   color=CATEGORY_COLOR[cat], alpha=0.68, zorder=3,
                   label=f"{cat} ({counts.get(cat, 0):,})")
    if n_keep >= 3:
        draw_oval(ax, np.column_stack([xk, yk]),
                  label="Covering oval (min-volume ellipse)")

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.grid(True, color=GRIDLINE, linewidth=0.7, zorder=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(BASELINE)
    ax.tick_params(colors=INK_MUTED, labelsize=8 if compact else 9)
    if n_keep == 0:
        ax.text(0.5, 0.5, "no reaction qualifies\nat this threshold",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=11 if not compact else 9.5, color=INK_MUTED)
    elif n_zero:
        ax.text(0.985, 0.03,
                f"{frac_zero:.0%} of the plotted points are ΔG′° = 0 in both sources",
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=7 if compact else 8, color=INK_MUTED)
    if show_legend:
        leg = ax.legend(loc="upper left", frameon=False, fontsize=8.5,
                        labelcolor=INK_PRIMARY, markerscale=1.7,
                        title="Reversibility transition", title_fontsize=8.5)
        leg.get_title().set_color(INK_SECONDARY)
    return {"n_keep": n_keep, "n_drop": n_drop, "r": r, "median_absdiff": med,
            "n_both_zero": n_zero, "frac_both_zero": frac_zero,
            "xlim": xlim, "ylim": ylim}


def slice_title(a: str, b: str, crit: str, tol: float | None, st: dict,
                n_base: int, compact: bool = False) -> str:
    head = f"{SOURCES[a][1]} vs {SOURCES[b][1]}"
    if tol is None:
        return (f"{head}\nno error filter · n = {st['n_keep']:,} · "
                f"r = {st['r']:.2f} · median |Δ| = {st['median_absdiff']:.2f} kcal/mol")
    frac = f"n = {st['n_keep']:,} of {n_base:,} ({st['n_keep']/max(n_base,1):.1%})"
    stats = f"r = {st['r']:.2f} · median |Δ| = {st['median_absdiff']:.2f} kcal/mol"
    if compact:
        # grid panels are ~4.4 in wide; three short lines beat two long ones,
        # which otherwise run under the neighbouring panel's title
        return f"{CRIT_SHORT[crit]} ≤ {tol:g} kcal/mol\n{frac}\n{stats}"
    return f"{head} · {CRIT_SHORT[crit]} ≤ {tol:g} kcal/mol\n{frac} · {stats}"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--no-singles", action="store_true",
                    help="write only the grid figures, skip the 72 single panels")
    args = ap.parse_args()
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"loading reactions from {MSDB_DATA} ...")
    df = attach_ehat(load_table())

    # A threshold below a source's floor makes every panel in that column empty,
    # which is a property of the source, not of the reactions.
    print("  floors (lowest value any reaction attains):")
    for key, (_, label, _) in SOURCES.items():
        print(f"    {label:24s} σ ≥ {df[f'sig_{key}'].min():.2f}   "
              f"ê ≥ {df[f'ehat_{key}'].min():.2f} kcal/mol")

    bases, n_extreme = {}, {}
    for a, b in PAIRS:
        bases[(a, b)], n_extreme[(a, b)] = base_set(df, a, b)
        print(f"  {SOURCES[a][1]} ∩ {SOURCES[b][1]}: {len(bases[(a,b)]):,} reactions "
              f"({n_extreme[(a,b)]} dropped as |ΔG′°| > {EXTREME_CUTOFF:g})")

    records = []

    # ---- reference: the same three pairs with no error filter at all --------
    fig, axes = plt.subplots(1, 3, figsize=(21, 7.0), dpi=130)
    fig.patch.set_facecolor(SURFACE)
    for ax, (a, b) in zip(axes, PAIRS):
        sub = bases[(a, b)]
        st = draw_panel(ax, sub, a, b, np.ones(len(sub), dtype=bool), "filtered")
        ax.set_title(slice_title(a, b, "none", None, st, len(sub)),
                     color=INK_PRIMARY, fontsize=11.5, pad=10)
        ax.set_xlabel(SOURCES[a][2], color=INK_SECONDARY, fontsize=10)
        ax.set_ylabel(SOURCES[b][2], color=INK_SECONDARY, fontsize=10)
        records.append({"pair": f"{a}_vs_{b}", "criterion": "none", "tolerance": np.nan,
                        "n_base": len(sub), **{k: st[k] for k in
                                               ("n_keep", "n_drop", "r", "median_absdiff",
                                                "n_both_zero", "frac_both_zero")}})
    fig.suptitle("ΔG′° agreement between thermodynamic sources — no error filter, "
                 "with minimum-volume covering oval",
                 color=INK_PRIMARY, fontsize=14, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    p = out_dir / "grid_unfiltered.png"
    fig.savefig(p, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {p}")

    # ---- the 4 threshold grids ---------------------------------------------
    for crit in CRITERIA:
        for mode in MODES:
            fig, axes = plt.subplots(len(PAIRS), len(THRESHOLDS),
                                     figsize=(4.45 * len(THRESHOLDS), 4.95 * len(PAIRS)),
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
                    keep = passing_mask(sub, a, b, crit, tol)
                    st = draw_panel(ax, sub, a, b, keep, mode,
                                    xlim=row_lim, ylim=row_lim,
                                    show_legend=(ri == 0 and ci == 0), compact=True)
                    ax.set_title(slice_title(a, b, crit, tol, st, len(sub),
                                             compact=True),
                                 color=INK_PRIMARY, fontsize=9, pad=6)
                    if ci == 0:
                        ax.text(-0.30, 0.5, f"{SOURCES[a][1]}\nvs\n{SOURCES[b][1]}",
                                transform=ax.transAxes, ha="center", va="center",
                                rotation=90, fontsize=11, color=INK_PRIMARY)
                    if ri == len(PAIRS) - 1:
                        ax.set_xlabel(SOURCES[a][2], color=INK_SECONDARY, fontsize=9)
                    if ci == 0:
                        ax.set_ylabel(SOURCES[b][2], color=INK_SECONDARY, fontsize=9)
                    if mode == "filtered":
                        records.append({"pair": f"{a}_vs_{b}", "criterion": crit,
                                        "tolerance": tol, "n_base": len(sub),
                                        **{k: st[k] for k in ("n_keep", "n_drop",
                                                              "r", "median_absdiff",
                                                              "n_both_zero",
                                                              "frac_both_zero")}})
            sub_note = ("only the passing reactions are drawn; axes rescale to them"
                        if mode == "filtered" else
                        "every reaction is drawn — failing reactions grayed out, "
                        "passing reactions in color and inside the oval; axes shared per row")
            fig.suptitle(
                f"ΔG′° source agreement filtered by {CRIT_LABEL[crit]} — "
                f"both sources ≤ T\n{sub_note}",
                color=INK_PRIMARY, fontsize=15, y=0.995)
            fig.tight_layout(rect=(0.022, 0, 1, 0.955))
            p = out_dir / f"grid_{crit}_{mode}.png"
            fig.savefig(p, facecolor=SURFACE)
            plt.close(fig)
            print(f"wrote {p}")

            if args.no_singles:
                continue
            single_dir = out_dir / crit / mode
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
                    keep = passing_mask(sub, a, b, crit, tol)
                    st = draw_panel(ax, sub, a, b, keep, mode,
                                    xlim=row_lim, ylim=row_lim)
                    ax.set_title(slice_title(a, b, crit, tol, st, len(sub)),
                                 color=INK_PRIMARY, fontsize=11.5, pad=10)
                    ax.set_xlabel(SOURCES[a][2], color=INK_SECONDARY, fontsize=10)
                    ax.set_ylabel(SOURCES[b][2], color=INK_SECONDARY, fontsize=10)
                    note = (f"filter: {CRIT_LABEL[crit]} ≤ {tol:g} kcal/mol "
                            f"for BOTH sources")
                    if n_extreme[(a, b)]:
                        note += (f"\n{n_extreme[(a,b)]} reaction(s) with |ΔG′°| > "
                                 f"{EXTREME_CUTOFF:g} kcal/mol excluded from the base set")
                    if st["n_keep"] == 0:
                        floor = min(sub[f"{CRIT_COLUMN[crit]}_{k}"].min()
                                    for k in (a, b))
                        note += (f"\nempty because the lower of the two sources' "
                                 f"{CRIT_SHORT[crit]} floors is {floor:.2f} kcal/mol")
                    ax.text(0.5, -0.125, note,
                            transform=ax.transAxes, ha="center", va="top",
                            fontsize=7.5, color=INK_MUTED, linespacing=1.6)
                    fig.tight_layout()
                    fig.subplots_adjust(bottom=0.135 + 0.028 * (note.count("\n") + 1))
                    fig.savefig(single_dir / f"dg_scatter_{a}_vs_{b}_{crit}_le{tol:g}.png",
                                facecolor=SURFACE)
                    plt.close(fig)
            print(f"  wrote {len(PAIRS) * len(THRESHOLDS)} single panels to {single_dir}")

    counts = pd.DataFrame(records)
    counts_path = out_dir / "slice_counts.tsv"
    counts.to_csv(counts_path, sep="\t", index=False, float_format="%.4f")
    print(f"wrote {counts_path}")
    print(counts.to_string(index=False))


if __name__ == "__main__":
    main()
