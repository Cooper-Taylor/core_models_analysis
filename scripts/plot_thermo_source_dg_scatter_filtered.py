#!/usr/bin/env python3
"""Pairwise DeltaG scatter plots across the three thermodynamic sources, sliced by
how much error each source ADMITS to on the reaction, with two ellipses.

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
    units. Read from results/eq_vs_dgpms_gcA/source_assignment.tsv -- refitted
    against the Convention A GC rebuild, since ehat_GC derived from the previous
    GC values no longer describes these ones. Keep reaction i iff
    max(ehat_a(i), ehat_b(i)) <= T.
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
          their categorical color and get the ellipses. Axes are shared down a
          column-group so the shrinking of the passing set is visible as an
          actual shrinking, not hidden by autoscale.

THE TWO ELLIPSES
----------------
Every panel carries two, and they answer different questions.

filled   95% CONCENTRATION ELLIPSE. Covariance shape, radius the EMPIRICAL 95th
         percentile of the Mahalanobis distance d^2 = (z-c)' S^-1 (z-c), so the
         coverage claim is distribution-free and exact rather than normal-theory:
         the legend prints the realised fraction inside. Its major axis is the
         principal axis of S, i.e. the total-least-squares / orthogonal
         regression fit -- the right fit here because BOTH coordinates carry
         error -- and its tilt is therefore readable against 45 degrees, which is
         perfect agreement. The semi-minor axis is the 95% orthogonal-residual
         scale. This is the ellipse to quote about the cloud.

outline  MINIMUM-VOLUME ENCLOSING ELLIPSE (Loewner-John, by Khachiyan's
         algorithm): the unique minimiser of -log det A subject to
         (z_i - c)' A (z_i - c) <= 1 for every point. Uniqueness is John's 1948
         theorem; in the plane it is fixed by at most d(d+3)/2 = 5 of the points
         (the legend prints how many are actually active), and shrinking it by a
         factor of d = 2 about its centre lands it inside the convex hull:
         (1/2)E is a subset of conv(P), itself a subset of E. So it is a smooth
         two-parameter stand-in for the convex hull -- an EXTREMES statistic with
         a breakdown point of zero, not a description of the cloud. On the
         unfiltered panels it runs 99-166x the area of the 95% ellipse, and
         deleting its 5 active points shrinks it by 37-51%. Read it as "how far
         into the plane does this subset reach", nothing more.

Both areas, tilts and the realised coverage go into slice_counts.tsv.

Point color is the reversibility transition between the two sources, identical
in definition and palette to plot_thermo_source_dg_scatter.py: each source's own
ΔG'° run through the unmodified cascade (DEFAULT_HEURISTICS via
per_source_energy), collapsed to reversible ("=") vs irreversible (">"/"<").

DATA PROVENANCE
---------------
ModelSEED dev @ 49563c6f (2026-08-07 GC rebuild landed). Group Contribution is
Convention A there; eQuilibrator and dGPredictor-ModelSEED are still Convention
B. That mismatch is a real concern for a cross-source scatter, and it was
tested: the GC-vs-eQuilibrator residual regressed on net H+ has slope -2.67
kcal/mol per proton, nowhere near the +/-9.539 a systematic A-vs-B gap would
produce, and the slope SHRANK from -3.39 under the old mixed-convention GC.
Consistent A vs consistent B cancels for mass-balanced reactions; what the
rebuild fixed was convention MIXING inside single reactions. So the panels are
drawn on raw stored values with no convention correction.

OUTPUTS  reports/thermoComparison/figures/thermo_source_dg_scatter_filtered/
    grid_unfiltered.png                   the 3 pairs, no threshold, both ellipses
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
# dGPredictor-ModelSEED energies. devsnap2 is dev @ 49563c6f, i.e. AFTER
# ad34d6ab "Rebuild GC energies under Convention A" -- Group Contribution values
# there differ from the earlier devsnap (34992d39) on 53% of reactions and cover
# 1,501 more. Code (the cascade) comes from the working ModelSEEDDatabase
# checkout.
MSDB_DATA = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/tmp/devsnap2"))
MSDB_CODE = Path(os.environ.get("MSDB_CODE_ROOT", "/scratch/ctaylor/ModelSEEDDatabase"))
ASSIGN_TSV = Path(os.environ.get(
    "EQDGP_ASSIGNMENT",
    str(ANALYSIS_DIR / "results" / "eq_vs_dgpms_gcA" / "source_assignment.tsv")))
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

# Everything draw_panel() reports that belongs in slice_counts.tsv, in order.
STAT_KEYS = ("n_keep", "n_drop", "r", "median_absdiff", "n_both_zero",
             "frac_both_zero", "conc_area", "conc_tilt", "conc_cover",
             "conc_semimajor", "conc_semiminor", "mvee_area", "mvee_tilt",
             "mvee_support")

EQ_SENTINEL = 100.0     # kcal/mol; eQuilibrator's "I have no estimate" marker
# A handful of aggregate/polymer reactions carry genuine but chemically
# implausible magnitudes (rxn05017, Group Contribution ~15,900 kcal/mol). They
# are dropped from the base set here rather than merely clipped, because a
# single such point would set the enclosing ellipse on its own -- it is fixed by
# at most 5 points. Count is reported per panel.
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
OVAL_COLOR = "#5a3fb0"        # the 95% concentration ellipse: filled, solid edge
MVEE_COLOR = "#9b8fd0"        # the enclosing ellipse: unfilled, dashed, thinner
MVEE_DASH = (0, (5, 3))
CONC_FRAC = 0.95
# Below this many points the 95th percentile of the Mahalanobis distance is
# essentially the maximum -- the "95% ellipse" would just be a covering ellipse
# with a misleading label -- so it is not drawn. The MVEE still is: it is exactly
# a covering ellipse and says so.
CONC_MIN_N = 20
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


def concentration_ellipse(points: np.ndarray, frac: float = CONC_FRAC) -> dict:
    """Empirical concentration ellipse: the shape of the covariance S, scaled to
    the ``frac`` quantile of the Mahalanobis distance of the points themselves.

    Returns {"ellipse": (center, width, height, angle_deg) or None, "cover":
    realised fraction inside, "area", "tilt"}. Because the radius is an EMPIRICAL
    quantile rather than sqrt(chi2_2,0.95) = 2.4477, the coverage statement holds
    without assuming bivariate normality -- ``cover`` is the number to print, and
    it exceeds ``frac`` when the distances tie (the ΔG′° = 0 clusters do exactly
    this, and can collapse the ellipse to a point: ellipse is then None).

    The major axis is the leading principal axis of S, which is the
    total-least-squares / orthogonal-regression fit line. That is the right fit
    for these panels -- ordinary least squares would assume the x source is
    error-free -- and it is why the tilt is worth reading against 45 degrees.
    """
    P = np.asarray(points, dtype=float)
    P = P[np.isfinite(P).all(axis=1)]
    out = {"ellipse": None, "cover": float("nan"), "area": float("nan"),
           "tilt": float("nan"), "semi_major": float("nan"),
           "semi_minor": float("nan")}
    if len(P) < CONC_MIN_N:
        return out
    c = P.mean(axis=0)
    S = np.cov(P.T)
    if not np.all(np.isfinite(S)):
        return out
    evals, evecs = np.linalg.eigh(S)
    d2 = np.einsum("ij,jk,ik->i", P - c, np.linalg.pinv(S), P - c)
    k2 = float(np.quantile(d2, frac))
    out["cover"] = float((d2 <= k2).mean())
    if k2 <= 0 or evals.min() <= 0:      # singular or degenerate: nothing to draw
        out["area"] = 0.0
        return out
    a, b = np.sqrt(k2 * evals[1]), np.sqrt(k2 * evals[0])   # eigh sorts ascending
    angle = float(np.degrees(np.arctan2(evecs[1, 1], evecs[0, 1])))
    out.update(ellipse=(c, 2.0 * a, 2.0 * b, angle), area=float(np.pi * a * b),
               tilt=angle % 180.0, semi_major=float(a), semi_minor=float(b))
    return out


def ellipse_support(points: np.ndarray, ell, tol: float = 1e-3) -> int:
    """How many DISTINCT points sit on the boundary of ``ell``.

    The exact MVEE has an active set of at most d(d+3)/2 = 5 in the plane, but
    this is not that number and must not be labelled as if it were, for two
    reasons. Khachiyan is run to tol=1e-4, so the returned ellipse is
    approximate and several points can land within 0.1% of its boundary without
    being active in the exact program (observed: 6 distinct points at q >= 0.999
    on eQ vs dGPMS at sigma <= 20). And reactions genuinely share ΔG′° values, so
    a single contact location can be several rows -- hence the de-duplication,
    which is what keeps this number near the theoretical bound instead of 9.
    """
    if ell is None:
        return 0
    c, w, h, ang = ell
    t = np.radians(ang)
    R = np.array([[np.cos(t), np.sin(t)], [-np.sin(t), np.cos(t)]])
    P = np.asarray(points, dtype=float)
    U = (P - c) @ R.T
    q = (U[:, 0] / (w / 2.0)) ** 2 + (U[:, 1] / (h / 2.0)) ** 2
    on = q >= 1.0 - tol
    if not on.any():
        return 0
    return int(len(np.unique(np.round(P[on], 6), axis=0)))


def proxy_handle(kind: str, label: str):
    """Legend entry that does not depend on the host panel having drawn anything."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    if kind == "conc":
        return Patch(facecolor=OVAL_COLOR, alpha=0.25, edgecolor=OVAL_COLOR,
                     linewidth=1.4, label=label)
    if kind == "mvee":
        return Patch(facecolor="none", edgecolor=MVEE_COLOR, linewidth=1.1,
                     linestyle=MVEE_DASH, label=label)
    color = EXCLUDED_COLOR if kind == "excluded" else CATEGORY_COLOR[kind]
    return Line2D([], [], marker="o", linestyle="none", markersize=5.5,
                  markerfacecolor=color, markeredgecolor="none", label=label)


def pair_floor(sub: pd.DataFrame, a: str, b: str, crit: str) -> float:
    """Smallest T at which ANY reaction can pass, i.e. min over reactions of
    max(crit_a, crit_b). Not min(crit_a.min(), crit_b.min()) -- the filter needs
    both sources under T on the SAME reaction, so a source's own floor says
    nothing on its own."""
    col = CRIT_COLUMN[crit]
    worse = np.maximum(sub[f"{col}_{a}"].to_numpy(float),
                       sub[f"{col}_{b}"].to_numpy(float))
    worse = worse[np.isfinite(worse)]
    return float(worse.min()) if worse.size else float("nan")


def draw_ellipses(ax, points: np.ndarray) -> dict:
    """Draw the 95% concentration ellipse (filled) and the minimum-volume
    enclosing ellipse (dashed outline), and return both sets of numbers."""
    conc = concentration_ellipse(points)
    if conc["ellipse"] is not None:
        c, w, h, ang = conc["ellipse"]
        # Where the sources agree this ellipse is a NEEDLE -- 95% of the reactions
        # sit within a few kcal/mol of the fit line while spanning hundreds along
        # it -- so it renders as a thick diagonal stroke buried in the point
        # cloud. That is the honest geometry (and the whole point next to the
        # balloon-shaped MVEE), so rather than widen it, give it a white halo and
        # put the half-width in the legend.
        import matplotlib.patheffects as pe
        ax.add_patch(Ellipse(tuple(c), w, h, angle=ang, facecolor=OVAL_COLOR,
                             alpha=0.10, edgecolor=OVAL_COLOR, linewidth=1.7,
                             linestyle="-", zorder=6,
                             path_effects=[pe.withStroke(linewidth=3.4,
                                                         foreground="white")]))
    ell = mvee(points)
    if ell is not None:
        c, w, h, ang = ell
        ax.add_patch(Ellipse(tuple(c), w, h, angle=ang, facecolor="none",
                             edgecolor=MVEE_COLOR, linewidth=1.0,
                             linestyle=MVEE_DASH, zorder=4))
    return {
        "conc_area": conc["area"], "conc_tilt": conc["tilt"],
        "conc_cover": conc["cover"], "conc_semimajor": conc["semi_major"],
        "conc_semiminor": conc["semi_minor"],
        "mvee_area": float(np.pi * ell[1] * ell[2] / 4.0) if ell else float("nan"),
        "mvee_tilt": (ell[3] % 180.0) if ell else float("nan"),
        "mvee_support": ellipse_support(points, ell),
    }


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
               compact: bool = False, excluded_label: str = "Excluded by threshold",
               empty_note: str = "no reaction qualifies\nat this threshold") -> dict:
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
                   label=f"{excluded_label} ({n_drop:,})")
    for cat in CATEGORY_ORDER:
        idx = keep & (cats == cat)
        if not idx.any():
            continue
        ax.scatter(xs[idx], ys[idx], s=13 if not compact else 9, linewidths=0,
                   color=CATEGORY_COLOR[cat], alpha=0.68, zorder=3,
                   label=f"{cat} ({counts.get(cat, 0):,})")
    geom = {"conc_area": float("nan"), "conc_tilt": float("nan"),
            "conc_cover": float("nan"), "conc_semimajor": float("nan"),
            "conc_semiminor": float("nan"), "mvee_area": float("nan"),
            "mvee_tilt": float("nan"), "mvee_support": 0}
    if n_keep >= 3:
        geom = draw_ellipses(ax, np.column_stack([xk, yk]))

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.grid(True, color=GRIDLINE, linewidth=0.7, zorder=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(BASELINE)
    ax.tick_params(colors=INK_MUTED, labelsize=8 if compact else 9)
    if n_keep == 0:
        ax.text(0.5, 0.5, empty_note,
                transform=ax.transAxes, ha="center", va="center",
                fontsize=11 if not compact else 9.5, color=INK_MUTED)
    elif n_zero:
        ax.text(0.985, 0.03,
                f"{frac_zero:.0%} of the plotted points are ΔG′° = 0 in both sources",
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=7 if compact else 8, color=INK_MUTED)
    if show_legend:
        # Built from proxy handles rather than from whatever this panel happens
        # to have drawn: the legend host panel can be empty (a threshold below
        # the pair's floor), and ax.legend() would then silently emit nothing.
        handles = ([proxy_handle("excluded", f"{excluded_label} ({n_drop:,})")]
                   if mode == "context" else [])
        handles += [proxy_handle(cat, f"{cat} ({counts.get(cat, 0):,})")
                    for cat in CATEGORY_ORDER]
        if np.isfinite(geom["conc_cover"]):
            if np.isfinite(geom["conc_semiminor"]):
                # the half-width is the readable quantity: 95% of the reactions
                # lie within it of the orthogonal fit line
                shape = (f"±{geom['conc_semiminor']:.1f} @ {geom['conc_tilt']:.0f}°"
                         if compact else
                         f"±{geom['conc_semiminor']:.1f} kcal/mol about the "
                         f"{geom['conc_tilt']:.1f}° orthogonal fit")
            else:
                shape = "degenerate"
            head = "95% conc. ellipse" if compact else "95% concentration ellipse"
            handles.append(proxy_handle(
                "conc", f"{head} ({geom['conc_cover']:.1%} inside, {shape})"))
        if np.isfinite(geom["mvee_area"]):
            handles.append(proxy_handle(
                "mvee", ("Min-volume enclosing ellipse "
                         f"({geom['mvee_support']} boundary pts)" if compact else
                         "Min-volume enclosing ellipse (hull proxy, "
                         f"{geom['mvee_support']} distinct points on its boundary)")))
        fs = 6.6 if compact else 8.5
        leg = ax.legend(handles=handles, loc="upper left", frameon=False,
                        fontsize=fs, labelcolor=INK_PRIMARY,
                        title="Reversibility transition", title_fontsize=fs,
                        labelspacing=0.32 if compact else 0.5,
                        handletextpad=0.4 if compact else 0.8,
                        borderpad=0.2 if compact else 0.4)
        leg.get_title().set_color(INK_SECONDARY)
    return {"n_keep": n_keep, "n_drop": n_drop, "r": r, "median_absdiff": med,
            "n_both_zero": n_zero, "frac_both_zero": frac_zero,
            "xlim": xlim, "ylim": ylim, **geom}


def ellipse_note(st: dict) -> str:
    """The mathematics of the two ellipses, with this panel's numbers in it, for
    the footnote of a full-size single panel."""
    import textwrap
    wrap = lambda s: "\n".join(textwrap.fill(ln, 112) for ln in s.split("\n"))
    if np.isfinite(st["conc_tilt"]):
        line = (f"filled = {CONC_FRAC:.0%} concentration ellipse — covariance shape, "
                f"radius the empirical {CONC_FRAC:.0%} quantile of Mahalanobis d², so "
                f"{st['conc_cover']:.1%} of the points lie inside; major axis is the "
                f"orthogonal (total-least-squares) fit at {st['conc_tilt']:.1f}° "
                f"(45° = agreement), semi-axes {st['conc_semimajor']:,.0f} × "
                f"{st['conc_semiminor']:,.1f} kcal/mol — i.e. {CONC_FRAC:.0%} of the "
                f"reactions sit within {st['conc_semiminor']:,.1f} kcal/mol of that "
                f"line — area {st['conc_area']:,.0f} (kcal/mol)²")
    elif np.isfinite(st["conc_cover"]):
        line = (f"filled = {CONC_FRAC:.0%} concentration ellipse — degenerate here "
                f"({st['conc_cover']:.1%} of the points coincide), nothing drawn")
    else:
        line = (f"no {CONC_FRAC:.0%} concentration ellipse: fewer than {CONC_MIN_N} "
                f"points, at which size its radius is just the maximum")
    if np.isfinite(st["mvee_area"]):
        ca = st.get("conc_area", float("nan"))
        ratio = (f", {st['mvee_area'] / ca:,.0f}× its area"
                 if np.isfinite(ca) and ca > 0 else "")
        line += (f"\ndashed = minimum-volume enclosing ellipse — argmin −log det A "
                 f"s.t. (zᵢ−c)ᵀA(zᵢ−c) ≤ 1 ∀i; unique (John 1948), ½E ⊆ conv(P) ⊆ E, "
                 f"and at most 5 points can be active in the plane — "
                 f"{st['mvee_support']} distinct points of {st['n_keep']:,} sit on this "
                 f"boundary (Khachiyan solved to 1e-4, so near-active points are not "
                 f"resolved from active ones); area {st['mvee_area']:,.0f}{ratio}")
    return wrap(line)


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
                        "n_base": len(sub), **{k: st[k] for k in STAT_KEYS}})
    fig.suptitle("ΔG′° agreement between thermodynamic sources — no error filter\n"
                 "filled ellipse = 95% concentration (empirical Mahalanobis quantile; "
                 "major axis = orthogonal fit) · dashed = minimum-volume enclosing "
                 "ellipse (Löwner–John hull proxy)",
                 color=INK_PRIMARY, fontsize=12.5, y=0.995)
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
                                    show_legend=True, compact=True)
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
                                        **{k: st[k] for k in STAT_KEYS}})
            sub_note = ("only the passing reactions are drawn; axes rescale to them"
                        if mode == "filtered" else
                        "every reaction is drawn — failing reactions grayed out, "
                        "passing reactions in color and inside the ellipses; axes shared per row")
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
                            f"for BOTH sources\n" + ellipse_note(st))
                    if n_extreme[(a, b)]:
                        note += (f"\n{n_extreme[(a,b)]} reaction(s) with |ΔG′°| > "
                                 f"{EXTREME_CUTOFF:g} kcal/mol excluded from the base set")
                    if st["n_keep"] == 0:
                        floor = pair_floor(sub, a, b, crit)
                        note += (f"\nempty because no reaction has both sources under "
                                 f"{CRIT_SHORT[crit]} = {tol:g}; the best any reaction "
                                 f"achieves is {floor:.2f} kcal/mol")
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
