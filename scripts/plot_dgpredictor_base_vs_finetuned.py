#!/usr/bin/env python3
"""ΔG′° scatter: the ORIGINAL (KEGG-trained) dGPredictor against the
ModelSEED-fine-tuned retrain, in the house style of
plot_thermo_source_dg_scatter.py.

Both series come from the same reaction records on ModelSEED ``dev``
(``/scratch/ctaylor/tmp/devsnap2``, dev @ 49563c6f), where they sit side by side
as two entries in the same ``thermodynamics`` dict:

    "dGPredictor"            [dg, sigma, operator]   original, KEGG-mapped
    "dGPredictor-ModelSEED"  [dg, sigma, operator]   freiburgermsu retrain

so the comparison needs no unit conversion, no re-prediction, and no ID
crosswalk: it is the same reaction, the same kcal/mol, the same stored-operator
convention, differing only in which dGPredictor produced the number.

Point colour is the reversibility TRANSITION between the two, defined exactly as
in plot_thermo_source_dg_scatter.py: each source's own ΔG′° is pushed through
the unmodified ModelSEED cascade (``reversibility_heuristics.DEFAULT_HEURISTICS``
via ``per_source_energy``) and the resulting call is collapsed to reversible
("=") vs irreversible (">"/"<"). "No change" is neutral gray rather than a 4th
categorical hue -- same reasoning, same palette.

THREE FIGURES
-------------
fig1_base_vs_finetuned.png
    All co-covered reactions, one panel.

fig2_split_by_kegg_mask.png
    The same points split into two panels by whether the ORIGINAL dGPredictor's
    value is trustworthy at all. ``results/thermo_agreement/dgpredictor_kegg_mask.json``
    lists the reactions that were staged a KEGG reaction id ModelSEED does not
    list as an alias -- the mis-mapping defect that drives most of the original
    predictor's disagreement with everything else. Splitting on it separates
    "the retrain changed the chemistry" from "the retrain removed a data defect",
    which the pooled panel cannot distinguish.

fig3_sigma.png
    Reported uncertainty, original vs retrain, on log axes. The retrain's sigma
    is the one that is actually calibrated (rho = +0.672 vs |eQ - dGP|), so the
    shift visible here is the substantive difference between the two models,
    not a cosmetic one.

Each scatter is written twice: once on axes that hold every point, and once
(``*_zoom.png``) clipped to +/- ZOOM_LIMIT kcal/mol, which is where essentially
all of the mass sits. The zoom is a rendering choice only -- every statistic
printed on a panel is computed over the panel's full point set, and the number
of off-scale points is stated on the figure.

Outputs to reports/thermoComparison/figures/dgpredictor_base_vs_finetuned/
alongside pair_stats.tsv.
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

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
# Data: the dev snapshot is the only checkout carrying BOTH dGPredictor variants.
MSDB_DATA = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/tmp/devsnap2"))
# Code (the cascade) comes from the working checkout.
MSDB_CODE = Path(os.environ.get("MSDB_CODE_ROOT", "/scratch/ctaylor/ModelSEEDDatabase"))
OUT_DIR = (ANALYSIS_DIR / "reports" / "thermoComparison" / "figures"
           / "dgpredictor_base_vs_finetuned")

sys.path.insert(0, str(MSDB_CODE / "Scripts" / "Thermodynamics"))
from reversibility_heuristics import (  # noqa: E402
    DEFAULT_HEURISTICS, run_reversibility, per_source_energy,
)
sys.path.insert(0, str(ANALYSIS_DIR / "scripts"))
from build_dgpredictor_kegg_mask import load_mask  # noqa: E402

BASE_KEY = "dGPredictor"
FT_KEY = "dGPredictor-ModelSEED"
BASE_AXIS = "dGPredictor (original, KEGG-trained) ΔG′° (kcal/mol)"
FT_AXIS = "dGPredictor-ModelSEED (fine-tuned) ΔG′° (kcal/mol)"

# Palette carried over unchanged from plot_thermo_source_dg_scatter.py.
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
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
BASELINE = "#c3c2b7"
SURFACE = "#fcfcfb"

# Same fixed, chemically-motivated cutoff used by the sibling scatter scripts:
# a handful of aggregate/polymer reactions carry genuine but implausible
# magnitudes and would otherwise crush the distribution onto a few pixels.
EXTREME_CUTOFF = 1500.0
# Axis clip for the *_zoom.png renderings. 250 kcal/mol holds >99% of both
# series; beyond it the panels are empty white space that squeezes the core.
ZOOM_LIMIT = 250.0


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
    """One row per non-EMPTY reaction: dg / sigma / cascade operator per variant."""
    rows = []
    for path in sorted(glob.glob(str(MSDB_DATA / "Biochemistry" / "reaction_*.json"))):
        for entry in json.load(open(path)):
            if entry.get("status") == "EMPTY":
                continue
            thermo = entry.get("thermodynamics") or {}
            row = {"rxn": entry["id"], "name": entry.get("name", "")}
            for key, subkey in (("BASE", BASE_KEY), ("FT", FT_KEY)):
                triple = thermo.get(subkey)
                if not triple or len(triple) < 3 or triple[2] in (None, "?"):
                    continue
                try:
                    dg, sig = float(triple[0]), abs(float(triple[1]))
                except (TypeError, ValueError):
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
    print(f"  {len(df):,} non-EMPTY reactions in {MSDB_DATA}")
    for key, subkey in (("BASE", BASE_KEY), ("FT", FT_KEY)):
        n = df.get(f"dg_{key}", pd.Series(dtype=float)).notna().sum()
        print(f"    {subkey:24s} {n:,} usable ΔG′°")
    return df


def _frame(ax) -> None:
    ax.grid(True, color=GRIDLINE, linewidth=0.8, zorder=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(BASELINE)
    ax.tick_params(colors=INK_MUTED, labelsize=9)


def draw_panel(ax, sub: pd.DataFrame, title: str, *, legend: bool = True,
               lims: tuple[float, float] | None = None) -> dict:
    xs = sub["dg_BASE"].to_numpy(float)
    ys = sub["dg_FT"].to_numpy(float)
    cats = [classify(a, b) for a, b in zip(sub["op_BASE"], sub["op_FT"])]
    counts = Counter(cats)

    r = np.corrcoef(xs, ys)[0, 1] if len(xs) > 1 else float("nan")
    rho = pd.Series(xs).corr(pd.Series(ys), method="spearman") if len(xs) > 1 else float("nan")
    d = ys - xs
    med = float(np.median(np.abs(d)))
    sign_flip = float(np.mean((np.sign(xs) * np.sign(ys)) < 0))

    if lims is None:
        lo = float(min(xs.min(), ys.min()))
        hi = float(max(xs.max(), ys.max()))
        pad = 0.06 * (hi - lo)
        lims = (lo - pad, hi + pad)

    ax.plot(lims, lims, linestyle="--", linewidth=1.2, color=BASELINE, zorder=1)
    for cat in CATEGORY_ORDER:
        idx = [i for i, c in enumerate(cats) if c == cat]
        if not idx:
            continue
        ax.scatter(xs[idx], ys[idx], s=14, linewidths=0, color=CATEGORY_COLOR[cat],
                   alpha=0.65, zorder=2, label=f"{cat} ({counts[cat]:,})")
    ax.set_xlim(*lims)
    ax.set_ylim(*lims)
    ax.set_xlabel(BASE_AXIS, color=INK_SECONDARY, fontsize=10)
    ax.set_ylabel(FT_AXIS, color=INK_SECONDARY, fontsize=10)
    ax.set_title(title, color=INK_PRIMARY, fontsize=11.5, pad=10)
    _frame(ax)
    if legend:
        lg = ax.legend(loc="upper left", frameon=False, fontsize=8.5,
                       labelcolor=INK_PRIMARY, markerscale=1.6,
                       title="Reversibility transition", title_fontsize=8.5)
        lg.get_title().set_color(INK_SECONDARY)
    n_offscale = int(np.sum((xs < lims[0]) | (xs > lims[1])
                            | (ys < lims[0]) | (ys > lims[1])))
    offscale_note = (f"\n{n_offscale} point(s) off-scale" if n_offscale else "")
    ax.text(0.98, 0.03,
            f"r = {r:.2f}   ρ = {rho:.2f}\nmedian |Δ| = {med:.2f} kcal/mol\n"
            f"sign flips = {100 * sign_flip:.1f}%{offscale_note}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5,
            color=INK_MUTED, linespacing=1.5)
    return {"n": int(len(sub)), "pearson_r": float(r), "spearman_rho": float(rho),
            "median_abs_delta": med, "frac_sign_flip": sign_flip,
            **{f"n_{c}": int(counts[c]) for c in CATEGORY_ORDER}}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"loading reaction records from {MSDB_DATA} ...")
    df = load_table()

    both = df[df["dg_BASE"].notna() & df["dg_FT"].notna()].copy()
    keep = ((both["dg_BASE"].abs() <= EXTREME_CUTOFF)
            & (both["dg_FT"].abs() <= EXTREME_CUTOFF))
    n_off = int((~keep).sum())
    both = both[keep]
    print(f"  co-covered: {len(both):,} reactions "
          f"({n_off} dropped with |ΔG′°| > {EXTREME_CUTOFF:g} kcal/mol)")
    print(f"  base-only : {int((df['dg_BASE'].notna() & df['dg_FT'].isna()).sum()):,}")
    print(f"  retrain-only: {int((df['dg_BASE'].isna() & df['dg_FT'].notna()).sum()):,}")

    stats = []
    mask = load_mask()
    is_bad = both.index.to_series().isin(mask).to_numpy()

    lo = float(min(both["dg_BASE"].min(), both["dg_FT"].min()))
    hi = float(max(both["dg_BASE"].max(), both["dg_FT"].max()))
    pad = 0.06 * (hi - lo)
    full_lims = (lo - pad, hi + pad)
    zoom_lims = (-ZOOM_LIMIT, ZOOM_LIMIT)

    for suffix, lims in (("", full_lims), ("_zoom", zoom_lims)):
        record = suffix == ""     # only tabulate stats once; they are lims-independent

        # ---- fig 1: pooled -------------------------------------------------
        fig, ax = plt.subplots(figsize=(7.2, 6.6), dpi=150)
        fig.patch.set_facecolor(SURFACE)
        ax.set_facecolor(SURFACE)
        s = draw_panel(ax, both,
                       f"Original dGPredictor vs ModelSEED fine-tuned retrain\n"
                       f"n = {len(both):,} reactions carrying both values "
                       f"(ModelSEED dev)", lims=lims)
        if record:
            stats.append({**s, "panel": "all"})
        ax.text(0.5, -0.115,
                f"{n_off} reaction(s) with |ΔG′°| > {EXTREME_CUTOFF:g} kcal/mol excluded "
                f"from the data; both values read from the same reaction record.",
                transform=ax.transAxes, ha="center", va="top", fontsize=7.5,
                color=INK_MUTED)
        fig.tight_layout()
        fig.subplots_adjust(bottom=0.13)
        p = out_dir / f"fig1_base_vs_finetuned{suffix}.png"
        fig.savefig(p, facecolor=SURFACE)
        plt.close(fig)
        print(f"wrote {p}")

        # ---- fig 2: split by the original predictor's KEGG mis-mapping -----
        fig, axes = plt.subplots(1, 2, figsize=(13.4, 6.6), dpi=150)
        fig.patch.set_facecolor(SURFACE)
        for ax, sel, title, tag in (
            (axes[0], ~is_bad,
             "KEGG id vouched by a ModelSEED alias\n(the original predictor's value is "
             "about the right reaction)", "vouched"),
            (axes[1], is_bad,
             "KEGG id inferred, not a ModelSEED alias\n(the original predictor's value "
             "belongs to a different reaction)", "mismapped"),
        ):
            ax.set_facecolor(SURFACE)
            sub = both[sel]
            s = draw_panel(ax, sub, f"{title}\nn = {len(sub):,}",
                           legend=(tag == "vouched"), lims=lims)
            if record:
                stats.append({**s, "panel": tag})
        fig.suptitle("Where the retrain differs from the original dGPredictor: chemistry "
                     "vs the KEGG mis-mapping defect",
                     color=INK_PRIMARY, fontsize=13, y=0.99)
        fig.tight_layout()
        p = out_dir / f"fig2_split_by_kegg_mask{suffix}.png"
        fig.savefig(p, facecolor=SURFACE)
        plt.close(fig)
        print(f"wrote {p}")

    # ---- fig 3: reported uncertainty ---------------------------------------
    fig, ax = plt.subplots(figsize=(6.8, 6.2), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    sb = both["sig_BASE"].to_numpy(float)
    sf = both["sig_FT"].to_numpy(float)
    ok = (sb > 0) & (sf > 0)
    slo = float(min(sb[ok].min(), sf[ok].min())) / 2.0
    shi = float(max(sb[ok].max(), sf[ok].max())) * 2.0
    ax.plot([slo, shi], [slo, shi], "--", linewidth=1.2, color=BASELINE, zorder=1)
    ax.scatter(sb[ok], sf[ok], s=10, linewidths=0, color="#5a3fb0", alpha=0.35, zorder=2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(slo, shi)
    ax.set_ylim(slo, shi)
    ax.set_xlabel("dGPredictor (original) reported σ (kcal/mol)",
                  color=INK_SECONDARY, fontsize=10)
    ax.set_ylabel("dGPredictor-ModelSEED reported σ (kcal/mol)",
                  color=INK_SECONDARY, fontsize=10)
    ax.set_title("Reported uncertainty: the retrain admits far more of it\n"
                 f"median σ  {np.median(sb[ok]):.2f} → {np.median(sf[ok]):.2f} kcal/mol"
                 f"   (n = {int(ok.sum()):,})",
                 color=INK_PRIMARY, fontsize=11.5, pad=10)
    _frame(ax)
    fig.tight_layout()
    p = out_dir / "fig3_sigma.png"
    fig.savefig(p, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {p}")

    tsv = out_dir / "pair_stats.tsv"
    pd.DataFrame(stats).set_index("panel").to_csv(tsv, sep="\t")
    print(f"wrote {tsv}")
    print(pd.DataFrame(stats).set_index("panel").to_string())


if __name__ == "__main__":
    main()
