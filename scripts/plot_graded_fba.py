#!/usr/bin/env python3
"""Figures for the graded thermo-source core-model comparison.

Reads the tables written by ``analyze_graded_fba.py`` and emits three figures to
``reports/thermoComparison/figures/graded_fba/``:

  fig1_growth.png              % of the 5,683 core models that grow, per variant,
                               beside the permissiveness of each variant's calls.
                               Two panels, NOT a dual axis -- growth and
                               permissiveness are different measures and a shared
                               y-scale would imply a comparison that is not there.
  fig2_direction_accuracy.png  agreement with the experimental (TECRDB) reference
                               direction, on all matched reactions and on the
                               subset where the experiment says "directional".
                               ``graded``/``graded_trusted`` are omitted: they
                               USE TECRDB, so they score 100% by construction.
  fig3_core_grades.png         grade composition of the 239 core reactions, per
                               source, plus which source the graded map picked.

Palette and ink tokens are the ones already used by
``plot_thermo_source_growth_bars.py`` so the report's figures read as one set.
Categorical slots validated for the adjacent pairlist; the grade ramp is the
ordinal blue ramp (validated with --ordinal).
"""
from __future__ import annotations

import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
DATA = ANALYSIS_DIR / "results" / "thermo_grades_fba"
OUT_DIR = ANALYSIS_DIR / "reports" / "thermoComparison" / "figures" / "graded_fba"

INK_PRIMARY, INK_SECONDARY, INK_MUTED = "#0b0b0b", "#52514e", "#898781"
GRIDLINE, BASELINE, SURFACE = "#e1e0d9", "#c3c2b7", "#fcfcfb"
# categorical slots 1,2,3,4,7 -- validated light, worst adjacent CVD dE 9.1
COLOR = {"gc": "#2a78d6", "eq": "#eb6834", "dgpms": "#1baf7a",
         "implicit": "#eda100", "graded": "#4a3aa7",
         "graded_trusted": "#4a3aa7", "graded_heldout": "#4a3aa7",
         "TECRDB": "#e87ba4"}
# ordinal blue ramp, dark = best; neutral for "no source"
GRADE_COLOR = {"GOLD": "#184f95", "SILVER": "#3987e5", "BRONZE": "#86b6ef",
               "none": "#e1e0d9"}
SHORT = {"implicit": "model's own\nbounds", "gc": "Group\nContrib.",
         "eq": "eQuili-\nbrator", "dgpms": "dGPredictor\n-ModelSEED",
         "graded": "graded\n(recomm.)", "graded_trusted": "graded\n(SILVER)",
         "graded_heldout": "graded\n(held out)"}


def _style(ax, ylabel="", title=""):
    ax.set_facecolor(SURFACE)
    ax.grid(True, axis="y", color=GRIDLINE, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(BASELINE)
    ax.tick_params(colors=INK_MUTED, labelsize=9)
    if ylabel:
        ax.set_ylabel(ylabel, color=INK_SECONDARY, fontsize=10.5)
    if title:
        ax.set_title(title, color=INK_PRIMARY, fontsize=11.5, pad=12, loc="left")


def _bars(ax, labels, values, colors, fmt="{:.1%}", rot=0):
    x = np.arange(len(labels))
    # 4px-equivalent rounded data-ends are not available on mpl bars; a small
    # inter-bar gap (width 0.68) plays the surface-gap role instead.
    ax.bar(x, values, width=0.68, color=colors, edgecolor="none", zorder=2)
    for xi, v in zip(x, values):
        ax.annotate(fmt.format(v), (xi, v), textcoords="offset points",
                    xytext=(0, 4), ha="center", fontsize=9, color=INK_SECONDARY)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8.5, color=INK_SECONDARY, rotation=rot)


def fig_growth(growth: pd.DataFrame) -> None:
    order = ["implicit", "gc", "eq", "dgpms", "graded", "graded_trusted"]
    g = growth.set_index("variant").loc[order]
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 4.9))
    fig.patch.set_facecolor(SURFACE)

    _bars(axes[0], [SHORT[v] for v in order], g.pct_grows.to_numpy(),
          [COLOR[v] for v in order])
    axes[0].set_ylim(0, 0.80)
    _style(axes[0], "share of models with biomass flux > 0",
           "Core models that grow (n = 5,683)")

    # "model's own bounds" makes no direction calls at all, so it has no
    # permissiveness to plot -- omit the bar rather than draw a zero.
    order2 = [v for v in order if v != "implicit"]
    g2 = g.loc[order2]
    _bars(axes[1], [SHORT[v] for v in order2], g2.frac_reversible_core.to_numpy(),
          [COLOR[v] for v in order2])
    axes[1].set_ylim(0, 1.02)
    _style(axes[1], "share of direction calls that are reversible '='",
           "How permissive each variant is (core reactions)")

    for ax in axes:
        ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(xmax=1, decimals=0))
    fig.suptitle("Growth is not a quality metric on its own — a more permissive "
                 "map grows more models", color=INK_PRIMARY, fontsize=12.5,
                 x=0.008, ha="left", y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUT_DIR / "fig1_growth.png", dpi=200, facecolor=SURFACE)
    plt.close(fig)


def fig_accuracy(acc: pd.DataFrame) -> None:
    order = ["gc", "eq", "dgpms", "graded_heldout"]
    subsets = [("all_stereo_exact", "All matched reactions (n = 802)"),
               ("reference_directional", "Experiment says DIRECTIONAL (n = 155)")]
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 4.9), sharey=True)
    fig.patch.set_facecolor(SURFACE)
    for ax, (key, title) in zip(axes, subsets):
        sub = acc[(acc.subset == key)].set_index("variant").loc[order]
        _bars(ax, [f"{SHORT[v]}\nn={n}" for v, n in zip(order, sub.n_scored)],
              sub.accuracy.to_numpy(), [COLOR[v] for v in order])
        ax.set_ylim(0, 1.06)
        ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(xmax=1, decimals=0))
        _style(ax, "agreement with experimental direction" if ax is axes[0] else "", title)
    fig.suptitle("Direction accuracy against the TECRDB experimental reference\n"
                 "graded (recommended) and graded (SILVER floor) are omitted — "
                 "they use TECRDB, so they score 100% by construction",
                 color=INK_PRIMARY, fontsize=12, x=0.008, ha="left", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    fig.savefig(OUT_DIR / "fig2_direction_accuracy.png", dpi=200, facecolor=SURFACE)
    plt.close(fig)


def fig_core_grades(core: pd.DataFrame) -> None:
    srcs = [("grade_TECRDB", "TECRDB"), ("grade_EQ", "eQuilibrator"),
            ("grade_DGPMS", "dGPredictor\n-ModelSEED"), ("grade_GC", "Group\nContribution")]
    tiers = ["GOLD", "SILVER", "BRONZE", "none"]
    n_core = len(core)
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.2),
                             gridspec_kw={"width_ratios": [1.5, 1]})
    fig.patch.set_facecolor(SURFACE)

    ax = axes[0]
    x = np.arange(len(srcs))
    bottom = np.zeros(len(srcs))
    for tier in tiers:
        vals = np.array([
            int((core[col] == tier).sum()) if tier != "none"
            else int(core[col].isna().sum()) for col, _ in srcs], dtype=float)
        ax.bar(x, vals, width=0.66, bottom=bottom, color=GRADE_COLOR[tier],
               edgecolor=SURFACE, linewidth=1.6, zorder=2, label=tier)
        for xi, (v, b) in enumerate(zip(vals, bottom)):
            if v >= 14:
                ax.annotate(f"{int(v)}", (xi, b + v / 2), ha="center", va="center",
                            fontsize=9,
                            color="#ffffff" if tier in ("GOLD", "SILVER") else INK_SECONDARY)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels([lab for _, lab in srcs], fontsize=8.5, color=INK_SECONDARY)
    ax.set_ylim(0, n_core * 1.02)
    _style(ax, f"core reactions (of {n_core})", "Grade of each source on the core set")
    ax.legend(frameon=False, fontsize=8.5, labelcolor=INK_SECONDARY, ncol=4,
              loc="upper center", bbox_to_anchor=(0.5, -0.13))

    ax = axes[1]
    pick = core.loc[core.graded_pick.notna() & (core.graded_pick != ""), "graded_pick"]
    counts = pick.value_counts()
    key = {"TECRDB": "TECRDB", "eQuilibrator": "eq", "dGPredictor-ModelSEED": "dgpms",
           "Group contribution": "gc"}
    labels = list(counts.index)
    cols = [COLOR[key[k]] for k in labels]
    _bars(ax, [k.replace("dGPredictor-ModelSEED", "dGPredictor\n-ModelSEED")
               .replace("Group contribution", "Group\nContribution") for k in labels],
          counts.to_numpy().astype(float), cols, fmt="{:.0f}")
    ax.set_ylim(0, counts.max() * 1.22)
    _style(ax, "core reactions", "Which source the graded map used")

    fig.suptitle("The 239 reactions that appear in at least one core model",
                 color=INK_PRIMARY, fontsize=12.5, x=0.008, ha="left", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUT_DIR / "fig3_core_grades.png", dpi=200, facecolor=SURFACE)
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    growth = pd.read_csv(DATA / "variant_growth.tsv", sep="\t")
    acc = pd.read_csv(DATA / "direction_accuracy.tsv", sep="\t")
    core = pd.read_csv(DATA / "core_reaction_grades.tsv", sep="\t", low_memory=False)
    fig_growth(growth)
    fig_accuracy(acc)
    fig_core_grades(core)
    print(f"wrote 3 figures to {OUT_DIR}")


if __name__ == "__main__":
    main()
