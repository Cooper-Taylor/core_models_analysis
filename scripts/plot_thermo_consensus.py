#!/usr/bin/env python3
"""Figure for the consensus-selection optimisation.

  A  the trade-off surface: coverage against the RMSE bar, for the exact oracle
     bound and the generalisable rule, with the hand-picked sigma <= 20 baseline
     marked so the gain is visible.
  B  why the objective is CCC and not Pearson r -- the two subsets where r ranks
     them backwards.
  C  the selected set on the dG axes, against what was dropped.
  D  per-|dG|-decile retention: the shape of the solution, including the top
     decile the rule deliberately gives up.

Palette: the project's validated triple plus neutral gray. No new hues.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
DATA = Path(os.environ.get("EQDGP_OUT", str(ANALYSIS_DIR / "results" / "eq_vs_dgpms")))
OUT_DIR = Path(os.environ.get(
    "EQDGP_FIGS",
    str(ANALYSIS_DIR / "reports" / "thermoComparison" / "figures" / "eq_vs_dgpms")))

BLUE, ORANGE, AQUA, NEUTRAL = "#2a78d6", "#eb6834", "#1baf7a", "#9c9a94"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#898781"
GRID, BASE, SURF = "#e1e0d9", "#c3c2b7", "#fcfcfb"


def style(ax):
    ax.set_facecolor(SURF)
    ax.grid(True, color=GRID, linewidth=0.8, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(BASE)
    ax.tick_params(colors=INK3, labelsize=9)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from optimize_thermo_consensus import load_selector, metrics, stratum_retention

    fr = pd.read_csv(DATA / "consensus_frontier.tsv", sep="\t")
    rule_cfg = json.loads((DATA / "consensus_rule.json").read_text())
    d = pd.read_csv(DATA / "key_subset_classified.tsv", sep="\t", low_memory=False)
    d["abs_dg_eq"] = d.dg_eq.abs()
    d["abs_net_proton"] = d.net_proton.abs()
    x, y = d.dg_eq.to_numpy(float), d.dg_dgp.to_numpy(float)
    keep = load_selector()(d)
    strata = pd.qcut(d.abs_dg_eq, 10, labels=False, duplicates="drop").to_numpy()

    fig, axes = plt.subplots(1, 4, figsize=(19.2, 5.0), dpi=150,
                             gridspec_kw={"width_ratios": [1.15, 0.85, 1, 1]})
    fig.patch.set_facecolor(SURF)

    # ---- A frontier
    ax = axes[0]; style(ax)
    o = fr[fr.solver == "oracle"].sort_values("rmse_bar")
    r_ = fr[(fr.solver == "rule") & (fr.feasible == True)].sort_values("rmse_bar")
    f_ = fr[(fr.solver == "rule_floored") & (fr.feasible == True)].sort_values("rmse_bar")
    ax.plot(o.rmse_bar, o.coverage * 100, "-o", color=NEUTRAL, lw=2, ms=5, zorder=3,
            label="oracle bound (selects on the answer)")
    ax.plot(r_.rmse_bar, r_.coverage * 100, "-o", color=ORANGE, lw=2.4, ms=6, zorder=4,
            label="fitted rule (generalises)")
    if len(f_):
        ax.plot(f_.rmse_bar, f_.coverage * 100, "--s", color=AQUA, lw=1.6, ms=4,
                zorder=3, label="rule + 15% floor in every decile")
    bm = metrics(x[(d.dgp_uncertainty <= 20).to_numpy()],
                 y[(d.dgp_uncertainty <= 20).to_numpy()])
    ax.scatter([bm["rmse"]], [(d.dgp_uncertainty <= 20).mean() * 100], s=130,
               marker="*", color=BLUE, zorder=6, label="hand-picked σ ≤ 20")
    op = rule_cfg["in_sample"]
    ax.scatter([op["rmse"]], [op["n"] / len(d) * 100], s=95, marker="D",
               facecolor="none", edgecolor=INK, lw=1.6, zorder=7,
               label="shipped operating point")
    ax.set_xscale("log")
    ax.set_xlabel("RMSE bar  (kcal/mol, log)", color=INK2, fontsize=10)
    ax.set_ylabel("% of reactions retained", color=INK2, fontsize=10)
    ax.legend(loc="lower right", frameon=False, fontsize=7.6, labelcolor=INK)
    ax.set_title("A · coverage vs the error you accept", color=INK, fontsize=10.5, pad=8)

    # ---- B why not r
    ax = axes[1]; style(ax)
    hi = np.abs(x) > 50
    lo = np.abs(x) <= 10
    labs = ["|ΔG| > 50", "|ΔG| ≤ 10"]
    rs = [np.corrcoef(x[hi], y[hi])[0, 1], np.corrcoef(x[lo], y[lo])[0, 1]]
    ms = [np.median(np.abs(x[hi] - y[hi])), np.median(np.abs(x[lo] - y[lo]))]
    xs = np.arange(2)
    ax.bar(xs - 0.2, rs, width=0.38, color=BLUE, zorder=3, label="Pearson r")
    ax2 = ax.twinx()
    ax2.bar(xs + 0.2, ms, width=0.38, color=ORANGE, zorder=3,
            label="median |Δ| (kcal/mol)")
    ax2.set_ylabel("median |Δ| (kcal/mol)", color=ORANGE, fontsize=9)
    ax2.tick_params(axis="y", colors=ORANGE, labelsize=9)
    ax2.spines["top"].set_visible(False)
    for sp in ("left", "bottom", "right"):
        ax2.spines[sp].set_color(BASE)
    ax.set_ylabel("Pearson r", color=BLUE, fontsize=9)
    ax.tick_params(axis="y", colors=BLUE)
    ax.set_xticks(xs); ax.set_xticklabels(labs, fontsize=9, color=INK2)
    ax.set_ylim(0, 1.34)
    ax2.set_ylim(0, max(ms) * 1.34)
    for xi, (rv, mv) in enumerate(zip(rs, ms)):
        ax.text(xi - 0.2, rv + 0.03, f"{rv:.2f}", ha="center", fontsize=8, color=BLUE)
        ax2.text(xi + 0.2, mv + max(ms) * 0.03, f"{mv:.1f}", ha="center",
                 fontsize=8, color=ORANGE)
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper center", frameon=False, fontsize=7.4,
              labelcolor=INK, ncol=1, borderpad=0.1)
    ax.set_title("B · why not maximise r:\nit ranks these backwards",
                 color=INK, fontsize=10.5, pad=8)

    # ---- C selected set
    ax = axes[2]; style(ax)
    W = 150
    ax.plot([-W, W], [-W, W], "--", lw=1.2, color=BASE, zorder=1)
    m1 = (np.abs(x) <= W) & (np.abs(y) <= W) & ~keep
    m2 = (np.abs(x) <= W) & (np.abs(y) <= W) & keep
    ax.scatter(x[m1], y[m1], s=6, lw=0, color=NEUTRAL, alpha=0.25, zorder=2,
               label=f"dropped  n={int((~keep).sum()):,}")
    ax.scatter(x[m2], y[m2], s=8, lw=0, color=AQUA, alpha=0.6, zorder=3,
               label=f"selected  n={int(keep.sum()):,}")
    ax.set_xlim(-W, W); ax.set_ylim(-W, W)
    ax.set_xlabel("eQuilibrator ΔG′° (kcal/mol)", color=INK2, fontsize=10)
    ax.set_ylabel("dGPredictor-ModelSEED ΔG′° (kcal/mol)", color=INK2, fontsize=10)
    ax.legend(loc="upper left", frameon=False, fontsize=8, labelcolor=INK, markerscale=2)
    ax.set_title(f"C · selected set: CCC {op['ccc']:.3f}, RMSE {op['rmse']:.2f},\n"
                 f"median |Δ| {op['median_absdiff']:.2f} kcal/mol",
                 color=INK, fontsize=10.5, pad=8)

    # ---- D stratum shape
    ax = axes[3]; style(ax)
    ret = stratum_retention(keep, strata) * 100
    ys = np.arange(len(ret))
    ax.barh(ys, ret, color=[ORANGE if v < 15 else AQUA for v in ret], height=0.68, zorder=3)
    ax.axvline(15, color=BASE, ls="--", lw=1.3, zorder=2)
    labs = []
    for b in range(len(ret)):
        s = d.loc[strata == b, "abs_dg_eq"]
        labs.append(f"{s.min():.0f}–{s.max():.0f}")
    ax.set_yticks(ys); ax.set_yticklabels(labs, fontsize=8, color=INK2)
    ax.set_xlabel("% retained", color=INK2, fontsize=10)
    ax.set_ylabel("|ΔG′°| decile (kcal/mol)", color=INK2, fontsize=9.5)
    ax.text(16, len(ret) - 0.6, "15% floor", fontsize=7.5, color=INK3)
    ax.set_title("D · what the rule gives up:\nthe top energy decile, deliberately",
                 color=INK, fontsize=10.5, pad=8)

    fig.suptitle("Optimising which reactions to trust: maximise coverage subject to "
                 "CCC, RMSE and slope guarantees",
                 color=INK, fontsize=13, y=0.99)
    fig.text(0.5, 0.005,
             "The oracle selects directly on |eQ − dGP|, which is the outcome being "
             "guaranteed — it is a bound, not a usable rule. The fitted rule uses only "
             "features knowable before comparing\n(dGPredictor σ, eQuilibrator σ, |ΔG|, "
             "participant count, aromatic rings, proton balance) and is 5-fold "
             "cross-validated.",
             ha="center", va="bottom", fontsize=7.8, color=INK3, linespacing=1.5)
    fig.tight_layout(rect=(0, 0.05, 1, 0.92))
    out = OUT_DIR / "fig7_consensus_optimization.png"
    fig.savefig(out, facecolor=SURF); plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
