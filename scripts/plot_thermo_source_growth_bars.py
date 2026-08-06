#!/usr/bin/env python3
"""Bar charts of biomass growth-flux statistics across all 5,683 Kegg2 core
models, one bar per direction/thermodynamic source:

  * Group Contribution, eQuilibrator (2.0), dGPredictor -- reaction direction
    from that source's own DeltaG, run through the unmodified ModelSEED
    heuristic cascade (see reports/thermoComparison/THERMO_SOURCE_FBA_PIPELINE.md).
  * Implicit -- no override at all; each model's own on-disk reaction bounds
    (whatever direction was implicitly baked in when the model was built).
  * ModelSEED (current) -- eQuilibrator, falling back to the reaction's
    existing Group-Contribution-backed reversibility when eQuilibrator has no
    data; i.e. what ``Estimate_Reaction_Reversibility.py EQ`` actually
    computes against the live database today.

Two figures:
  1. Median +/- standard deviation of growth flux over ALL 5,683 models
     (non-growing models contribute a flux of 0).
  2. Median +/- standard deviation of growth flux over only the models that
     grow under that source (n differs per bar -- annotated on each bar).

Each bar also overlays open circles for individual models whose growth flux
falls beyond median +/- OUTLIER_SIGMA standard deviations (jittered
horizontally to reduce overlap). A literal +/-1 std -- the whisker actually
drawn -- flags 30-60% of models per source on this skewed/bimodal data, so
OUTLIER_SIGMA defaults to 2, which gives a sparse, plottable outlier set.

Source: results/thermo_source_fba_all_models/model_results.csv
(built by run_thermo_source_fba_all_models.py).
Outputs: reports/thermoComparison/figures/thermo_source_dg_scatter/ (kept alongside the
pairwise DeltaG scatter plots from the same analysis).
"""
from __future__ import annotations

import csv
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR", "/scratch/ctaylor/core_models_analysis"))
RESULTS_CSV = ANALYSIS_DIR / "results" / "thermo_source_fba_all_models" / "model_results.csv"
OUT_DIR = ANALYSIS_DIR / "reports" / "thermoComparison" / "figures" / "thermo_source_growth_bars"

# Bar order + labels + categorical color (dataviz skill: bar/adjacent-pair
# charts clear CVD-safety for all 8 default slots, so 5 distinct hues here
# are safe -- unlike the 3-slot cap that applies to the all-pairs scatter
# plots elsewhere in this analysis).
SOURCES = [
    ("group_contribution", "Group\nContribution", "#2a78d6"),   # slot 1 blue
    ("equilibrator", "eQuilibrator\n(2.0)", "#eb6834"),          # slot 2 orange
    ("dgpredictor", "dGPredictor", "#1baf7a"),                    # slot 3 aqua
    ("implicit", "Implicit\n(on-disk)", "#eda100"),               # slot 4 yellow
    ("modelseed_current", "ModelSEED\n(current)", "#e87ba4"),    # slot 5 magenta
]

INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
BASELINE = "#c3c2b7"
SURFACE = "#fcfcfb"


def load_rows() -> list:
    with open(RESULTS_CSV) as fh:
        return list(csv.DictReader(fh))


# Outlier threshold for the open-circle markers: beyond median +/- OUTLIER_SIGMA
# standard deviations. A literal +/-1 std (the whisker actually drawn) flags
# 30-60% of models per source here -- these growth-flux distributions are
# skewed/bimodal (many exact-zero non-growers, a wide grower spread), so that
# threshold is not a small "outlier" set, it's most of the data. +/-2 std gives
# a sparse, plottable set (0-7% of models per source).
OUTLIER_SIGMA = 2.0
_rng = np.random.default_rng(42)  # fixed seed: deterministic jitter across reruns


def make_bar_chart(rows: list, only_growers: bool, out_path: Path, title: str) -> None:
    labels, medians, stds, ns, colors, outlier_vals_list = [], [], [], [], [], []
    for src, label, color in SOURCES:
        if only_growers:
            vals = np.array([float(r[f"fba_growth_flux_{src}"]) for r in rows
                              if r[f"fba_grows_{src}"] == "True"])
        else:
            vals = np.array([float(r[f"fba_growth_flux_{src}"]) for r in rows])
        med = np.median(vals) if len(vals) else 0.0
        std = np.std(vals, ddof=1) if len(vals) > 1 else 0.0
        labels.append(label)
        medians.append(med)
        stds.append(std)
        ns.append(len(vals))
        colors.append(color)
        outlier_vals_list.append(vals[(vals < med - OUTLIER_SIGMA * std) | (vals > med + OUTLIER_SIGMA * std)])

    medians = np.array(medians)
    stds = np.array(stds)
    # Growth flux cannot be negative -- clip the lower error whisker at 0
    # rather than letting it imply negative growth; document the asymmetry.
    lower_err = np.minimum(stds, medians)
    upper_err = stds

    fig, ax = plt.subplots(figsize=(7.5, 6), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    x = np.arange(len(labels))
    bars = ax.bar(x, medians, yerr=[lower_err, upper_err], capsize=5, width=0.6,
                    color=colors, edgecolor="none", zorder=2,
                    error_kw={"ecolor": INK_SECONDARY, "elinewidth": 1.3, "capthick": 1.3, "zorder": 3})

    # Outlier markers use a single neutral dark outline (not each bar's own
    # hue) -- most of these outliers land BELOW the bar's median (inside the
    # solid fill area), where a same-hue outline would be invisible against
    # the same-color fill behind it.
    y_top = float(max(medians + upper_err))
    for xi, outliers in zip(x, outlier_vals_list):
        if len(outliers) == 0:
            continue
        jitter = _rng.uniform(-0.22, 0.22, size=len(outliers))
        ax.scatter(xi + jitter, outliers, s=24, facecolors="none",
                    edgecolors=INK_PRIMARY, linewidths=1.1, alpha=0.85, zorder=5)
        y_top = max(y_top, float(outliers.max()))

    label_y = y_top * 1.06
    for xi, med, n, err, outliers in zip(x, medians, ns, upper_err, outlier_vals_list):
        ax.text(xi, label_y, f"n={n:,}\n{len(outliers)} outlier(s)",
                 ha="center", va="bottom", fontsize=8, color=INK_MUTED, linespacing=1.4)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9.5, color=INK_SECONDARY)
    ax.set_ylabel("Growth flux (biomass objective) -- median ± std dev",
                    color=INK_SECONDARY, fontsize=11)
    ax.set_title(title, color=INK_PRIMARY, fontsize=12.5, pad=14)
    ax.set_ylim(0, label_y * 1.16)
    ax.grid(True, axis="y", color=GRIDLINE, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(BASELINE)
    ax.tick_params(colors=INK_MUTED, labelsize=9)

    handle = plt.Line2D([0], [0], marker="o", linestyle="none", markersize=6,
                          markerfacecolor="none", markeredgecolor=INK_SECONDARY)
    ax.legend([handle], [f"model beyond median ± {OUTLIER_SIGMA:.0f}σ"],
               loc="upper left", frameon=False, fontsize=8.5, labelcolor=INK_SECONDARY)

    fig.tight_layout()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out_path}")
    for label, med, std, n, outliers in zip(labels, medians, stds, ns, outlier_vals_list):
        print(f"  {label.replace(chr(10), ' ')}: n={n:,} median={med:.2f} std={std:.2f} "
              f"outliers(beyond {OUTLIER_SIGMA:.0f}σ)={len(outliers)}")


def main() -> None:
    rows = load_rows()
    print(f"loaded {len(rows)} model rows from {RESULTS_CSV}")

    make_bar_chart(
        rows, only_growers=False,
        out_path=OUT_DIR / "growth_flux_median_std_all_models.png",
        title="Growth flux by direction source -- all 5,683 core models\n"
              "(non-growing models contribute a flux of 0)",
    )
    make_bar_chart(
        rows, only_growers=True,
        out_path=OUT_DIR / "growth_flux_median_std_growing_only.png",
        title="Growth flux by direction source -- growing core models only\n"
              "(n differs per source; see labels)",
    )


if __name__ == "__main__":
    main()
