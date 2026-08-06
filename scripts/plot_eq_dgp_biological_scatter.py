#!/usr/bin/env python3
"""eQuilibrator vs dGPredictor-ModelSEED DeltaG scatter, with the biological
question overlaid: are the reactions that actually matter to the core models in
the confident half of the comparison, or the doubtful half?

Panels
  A  all key-subset reactions, coloured by dGPredictor's own posterior sigma
     tier -- the structure that predicts disagreement.
  B  the same axes with the core-model reactions picked out: those present in
     the 5,683 Kegg2 core models, sized by how often flipping their direction
     changes predicted growth (results/reaction_effects_all/).
  C  agreement against model prevalence, which is the summary answer.

"Biologically significant" is operationalised two ways, both from this repo
rather than asserted: presence in the core models
(site/data/reaction_model_counts.json, "all" = how many of 5,683 contain it)
and direction-sensitivity of growth (fraction of models where swapping the
reaction's bound direction moves the FBA objective).

Palette: the project's existing validated triple plus neutral gray.
"""
from __future__ import annotations

import glob
import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
DATA = ANALYSIS_DIR / "results" / "eq_vs_dgpms"
OUT_DIR = ANALYSIS_DIR / "reports" / "thermoComparison" / "figures" / "eq_vs_dgpms"
EFFECTS = ANALYSIS_DIR / "results" / "reaction_effects_all" / "effects"
COUNTS = ANALYSIS_DIR / "site" / "data" / "reaction_model_counts.json"
CACHE = ANALYSIS_DIR / "results" / "eq_vs_dgpms" / "rxn_growth_sensitivity.tsv"

BLUE, ORANGE, AQUA, NEUTRAL = "#2a78d6", "#eb6834", "#1baf7a", "#9c9a94"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#898781"
GRID, BASE, SURF = "#e1e0d9", "#c3c2b7", "#fcfcfb"
W = 150.0


def style(ax):
    ax.set_facecolor(SURF)
    ax.grid(True, color=GRID, linewidth=0.8, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(BASE)
    ax.tick_params(colors=INK3, labelsize=9)


def _one(f):
    d = pd.read_parquet(f, columns=["base", "dir", "growth", "status"])
    d = d[d.status == "optimal"]
    if d.empty:
        return None
    g = d.groupby("base")["growth"]
    return (g.max() - g.min()).rename("span").reset_index()


def growth_sensitivity() -> pd.DataFrame:
    """Per reaction: in what fraction of core models does its DIRECTION change
    the FBA objective. Cached, because it reads 5,683 parquet files."""
    if CACHE.exists():
        return pd.read_csv(CACHE, sep="\t")
    files = sorted(glob.glob(str(EFFECTS / "*.parquet")))
    if not files:
        return pd.DataFrame(columns=["rxn", "frac_models_growth_sensitive"])
    with ProcessPoolExecutor(24) as ex:
        parts = [p for p in ex.map(_one, files, chunksize=40) if p is not None]
    allp = pd.concat(parts, ignore_index=True)
    agg = (allp.groupby("base")
           .agg(n_models_tested=("span", "size"),
                frac_models_growth_sensitive=("span", lambda s: float((s > 1e-6).mean())),
                median_growth_span=("span", "median"))
           .reset_index().rename(columns={"base": "rxn"}))
    agg.to_csv(CACHE, sep="\t", index=False)
    return agg


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    k = pd.read_csv(DATA / "key_subset_classified.tsv", sep="\t", low_memory=False)
    k["absdiff"] = (k.dg_eq - k.dg_dgp).abs()
    counts = json.load(open(COUNTS))
    k["prev"] = [counts.get(r, {}).get("all", 0) for r in k.rxn]
    k = k.merge(growth_sensitivity(), on="rxn", how="left")
    k["tier"] = np.where(k.dgp_uncertainty <= 3, "high",
                         np.where(k.dgp_uncertainty <= 30, "medium", "low"))

    fig, axes = plt.subplots(1, 3, figsize=(16.6, 5.6), dpi=150,
                             gridspec_kw={"width_ratios": [1, 1, 1.1]})
    fig.patch.set_facecolor(SURF)

    # ---- A: everything, by confidence tier
    ax = axes[0]; style(ax)
    ax.plot([-W, W], [-W, W], "--", lw=1.2, color=BASE, zorder=1)
    for lab, col, z in (("medium", NEUTRAL, 2), ("low", ORANGE, 3), ("high", AQUA, 4)):
        s = k[k.tier == lab]
        x, y = s.dg_eq.to_numpy(float), s.dg_dgp.to_numpy(float)
        m = (np.abs(x) <= W) & (np.abs(y) <= W)
        ax.scatter(x[m], y[m], s=7, lw=0, color=col, alpha=0.45, zorder=z,
                   label=f"σ {lab}  n={len(s):,}  med|Δ|={s.absdiff.median():.2f}")
    ax.set_xlim(-W, W); ax.set_ylim(-W, W)
    ax.set_xlabel("eQuilibrator ΔG′° (kcal/mol)", color=INK2, fontsize=10)
    ax.set_ylabel("dGPredictor-ModelSEED ΔG′° (kcal/mol)", color=INK2, fontsize=10)
    ax.legend(loc="upper left", frameon=False, fontsize=7.8, labelcolor=INK, markerscale=2.2)
    ax.set_title(f"A · all {len(k):,} reconciled reactions,\nby the model's own confidence",
                 color=INK, fontsize=10.5, pad=8)

    # ---- B: the biologically real ones
    ax = axes[1]; style(ax)
    core = k[k.prev > 0].copy()
    rest = k[k.prev == 0]
    ax.plot([-W, W], [-W, W], "--", lw=1.2, color=BASE, zorder=1)
    x, y = rest.dg_eq.to_numpy(float), rest.dg_dgp.to_numpy(float)
    m = (np.abs(x) <= W) & (np.abs(y) <= W)
    ax.scatter(x[m], y[m], s=5, lw=0, color=NEUTRAL, alpha=0.22, zorder=2,
               label=f"in no core model  n={len(rest):,}")
    x, y = core.dg_eq.to_numpy(float), core.dg_dgp.to_numpy(float)
    fr = core.frac_models_growth_sensitive.fillna(0).to_numpy()
    m = (np.abs(x) <= W) & (np.abs(y) <= W)
    ax.scatter(x[m], y[m], s=18 + 190 * fr[m], lw=0.4, edgecolor="white",
               color=BLUE, alpha=0.85, zorder=4,
               label=f"in ≥1 core model  n={len(core):,}")
    ax.set_xlim(-W, W); ax.set_ylim(-W, W)
    ax.set_xlabel("eQuilibrator ΔG′° (kcal/mol)", color=INK2, fontsize=10)
    ax.legend(loc="upper left", frameon=False, fontsize=7.8, labelcolor=INK, markerscale=1.4)
    ax.text(0.97, 0.05, "marker area ∝ fraction of models\nwhere its direction moves growth",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.2, color=INK3)
    ax.set_title(f"B · core-model reactions sit on the diagonal\n"
                 f"median |Δ| {core.absdiff.median():.2f} vs {rest.absdiff.median():.2f} kcal/mol",
                 color=INK, fontsize=10.5, pad=8)

    # ---- C: the summary answer
    ax = axes[2]; style(ax)
    bins = [(-1, 0, "none"), (0, 499, "1–499"), (499, 1999, "500–2k"), (1999, 10**9, "≥2,000")]
    labs, meds, hi, ns = [], [], [], []
    for lo, h, lab in bins:
        s = k[(k.prev > lo) & (k.prev <= h)]
        if len(s) < 5:
            continue
        labs.append(lab); ns.append(len(s))
        meds.append(s.absdiff.median())
        hi.append((s.tier == "high").mean() * 100)
    xs = np.arange(len(labs))
    ax.bar(xs - 0.2, meds, width=0.38, color=ORANGE, zorder=3, label="median |Δ| (kcal/mol)")
    ax2 = ax.twinx()
    ax2.bar(xs + 0.2, hi, width=0.38, color=AQUA, zorder=3, label="% in high-confidence tier")
    ax2.set_ylabel("% of reactions with σ ≤ 3", color=AQUA, fontsize=9.5)
    ax2.tick_params(axis="y", colors=AQUA, labelsize=9)
    ax2.spines["top"].set_visible(False)
    for sp in ("left", "bottom", "right"):
        ax2.spines[sp].set_color(BASE)
    ax2.set_ylim(0, 100)
    ax.set_ylabel("median |Δ| (kcal/mol)", color=ORANGE, fontsize=9.5)
    ax.tick_params(axis="y", colors=ORANGE)
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{l}\nn={n:,}" for l, n in zip(labs, ns)], fontsize=8.5, color=INK2)
    ax.set_xlabel("core models containing the reaction (of 5,683)", color=INK2, fontsize=10)
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper center", frameon=False, fontsize=8, labelcolor=INK)
    ax.set_title("C · the more biologically real a reaction,\nthe better the two agree",
                 color=INK, fontsize=10.5, pad=8)

    fig.suptitle("Do the two methods disagree on reactions that matter? Largely no — "
                 "the disagreement is in the periphery",
                 color=INK, fontsize=13, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    out = OUT_DIR / "fig6_biological_significance.png"
    fig.savefig(out, facecolor=SURF); plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
