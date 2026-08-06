#!/usr/bin/env python3
"""Bar charts of thermodynamic-data coverage across the 5,683 Kegg2 core
models, one bar per thermo source (Group Contribution, eQuilibrator (2.0),
dGPredictor, ModelSEED): what fraction / how many of a model's unique
reactions have a defined direction, and of its unique compounds have a
defined formation energy, under that source.

"ModelSEED" here is the canonical top-level source (`modelseed` column in
model_results.csv) rather than the eQuilibrator-with-GC-fallback
`modelseed_current` used for the growth-flux bar charts: the two are
numerically identical for reaction coverage in this checkout (same 202
reactions), but only `modelseed` has a compound-level coverage column
(there is no compound-level "EQ, fallback to GC" concept in this database --
dGPredictor and eQuilibrator/GC don't share a compound energy table the way
reactions do), so `modelseed` is the only one that covers both metabolites
and reactions consistently.

Four "per-model sample" charts (n = 5,683 models per bar, median +/- std dev,
open circles for models beyond median +/- 2 std dev, same convention as
plot_thermo_source_growth_bars.py):
  1. coverage_pct_reactions_per_model.png   -- % of a model's reactions covered
  2. coverage_abs_reactions_per_model.png   -- count of a model's reactions covered
  3. coverage_pct_compounds_per_model.png   -- % of a model's compounds covered
  4. coverage_abs_compounds_per_model.png   -- count of a model's compounds covered

Plus one single-value chart (no per-model sampling, no error bars):
  5. coverage_pct_combined_all_models.png -- % of the 239 combined-unique
     reactions / 182 combined-unique compounds (across all models) covered,
     grouped bars (reactions vs compounds) per source.

Source: results/thermo_source_fba_all_models/model_results.csv and
summary_stats.json (built by run_thermo_source_fba_all_models.py).
Outputs: reports/thermoComparison/figures/thermo_source_coverage_bars/.
"""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR", "/scratch/ctaylor/core_models_analysis"))
DATA_DIR = ANALYSIS_DIR / "results" / "thermo_source_fba_all_models"
RESULTS_CSV = DATA_DIR / "model_results.csv"
SUMMARY_JSON = DATA_DIR / "summary_stats.json"
OUT_DIR = ANALYSIS_DIR / "reports" / "thermoComparison" / "figures" / "thermo_source_coverage_bars"

# Same source -> color mapping as plot_thermo_source_growth_bars.py ("color
# follows the entity"): ModelSEED keeps the magenta slot used there for
# modelseed_current, since the two sources represent the same conceptual
# entity and are numerically identical for reaction coverage in this data.
SOURCES = [
    ("group_contribution", "Group\nContribution", "#2a78d6"),
    ("equilibrator", "eQuilibrator\n(2.0)", "#eb6834"),
    ("dgpredictor", "dGPredictor", "#1baf7a"),
    ("modelseed", "ModelSEED", "#e87ba4"),
]

INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
BASELINE = "#c3c2b7"
SURFACE = "#fcfcfb"

OUTLIER_SIGMA = 2.0
_rng = np.random.default_rng(42)


def load_rows() -> list:
    with open(RESULTS_CSV) as fh:
        return list(csv.DictReader(fh))


def make_sample_bar_chart(values_by_source: dict, out_path: Path, title: str, ylabel: str) -> None:
    """One bar per source: median +/- std dev over per-model samples, with
    open-circle outliers beyond median +/- OUTLIER_SIGMA std dev (identical
    convention to plot_thermo_source_growth_bars.py)."""
    labels, medians, stds, ns, colors, outlier_list = [], [], [], [], [], []
    for src, label, color in SOURCES:
        vals = values_by_source[src]
        med = np.median(vals)
        std = np.std(vals, ddof=1)
        labels.append(label)
        medians.append(med)
        stds.append(std)
        ns.append(len(vals))
        colors.append(color)
        outlier_list.append(vals[(vals < med - OUTLIER_SIGMA * std) | (vals > med + OUTLIER_SIGMA * std)])

    medians = np.array(medians)
    stds = np.array(stds)
    lower_err = np.minimum(stds, medians)  # cannot go below 0 (count/percentage)
    upper_err = stds

    fig, ax = plt.subplots(figsize=(7.5, 6), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    x = np.arange(len(labels))
    ax.bar(x, medians, yerr=[lower_err, upper_err], capsize=5, width=0.6,
            color=colors, edgecolor="none", zorder=2,
            error_kw={"ecolor": INK_SECONDARY, "elinewidth": 1.3, "capthick": 1.3, "zorder": 3})

    y_top = float(max(medians + upper_err))
    for xi, outliers in zip(x, outlier_list):
        if len(outliers) == 0:
            continue
        jitter = _rng.uniform(-0.22, 0.22, size=len(outliers))
        ax.scatter(xi + jitter, outliers, s=24, facecolors="none",
                    edgecolors=INK_PRIMARY, linewidths=1.1, alpha=0.85, zorder=5)
        y_top = max(y_top, float(outliers.max()))

    label_y = y_top * 1.06
    for xi, n, outliers in zip(x, ns, outlier_list):
        ax.text(xi, label_y, f"n={n:,}\n{len(outliers)} outlier(s)",
                 ha="center", va="bottom", fontsize=8, color=INK_MUTED, linespacing=1.4)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9.5, color=INK_SECONDARY)
    ax.set_ylabel(ylabel, color=INK_SECONDARY, fontsize=11)
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
    for label, med, std, n, outliers in zip(labels, medians, stds, ns, outlier_list):
        print(f"  {label.replace(chr(10), ' ')}: n={n:,} median={med:.2f} std={std:.2f} "
              f"outliers(beyond {OUTLIER_SIGMA:.0f}σ)={len(outliers)}")


def make_combined_pct_chart(summary: dict, out_path: Path) -> None:
    """Single-value grouped bars: % of the combined-unique reactions/compounds
    (across all 5,683 models) covered by each source -- no per-model sampling,
    no error bars."""
    n_rxn_total = summary["combined_unique_reactions_all_models"]
    n_cpd_total = summary["combined_unique_compounds_all_models"]
    rxn_by_src = summary["combined_reactions_with_direction_by_source"]
    cpd_by_src = summary["combined_compounds_with_energy_by_source"]

    labels = [lbl for _, lbl, _ in SOURCES]
    colors = [c for _, _, c in SOURCES]
    rxn_pct = [100 * rxn_by_src[src] / n_rxn_total for src, _, _ in SOURCES]
    cpd_pct = [100 * cpd_by_src.get(src, 0) / n_cpd_total for src, _, _ in SOURCES]

    fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    x = np.arange(len(labels))
    width = 0.32
    bars_rxn = ax.bar(x - width / 2, rxn_pct, width=width, color=colors,
                        edgecolor="none", zorder=2, label="Reactions")
    bars_cpd = ax.bar(x + width / 2, cpd_pct, width=width, color=colors,
                        edgecolor=INK_PRIMARY, linewidth=1.1, alpha=0.55, zorder=2,
                        hatch="///", label="Metabolites")

    for xi, v in zip(x - width / 2, rxn_pct):
        ax.text(xi, v + 1.5, f"{v:.0f}%", ha="center", va="bottom", fontsize=8.5, color=INK_MUTED)
    for xi, v in zip(x + width / 2, cpd_pct):
        ax.text(xi, v + 1.5, f"{v:.0f}%", ha="center", va="bottom", fontsize=8.5, color=INK_MUTED)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9.5, color=INK_SECONDARY)
    ax.set_ylabel("% of combined-unique reactions/compounds with defined energy",
                    color=INK_SECONDARY, fontsize=10.5)
    ax.set_title(f"Coverage across all combined core models\n"
                 f"({n_rxn_total} unique reactions, {n_cpd_total} unique compounds total)",
                 color=INK_PRIMARY, fontsize=12.5, pad=14)
    ax.set_ylim(0, 100)
    ax.grid(True, axis="y", color=GRIDLINE, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(BASELINE)
    ax.tick_params(colors=INK_MUTED, labelsize=9)

    rxn_handle = plt.Rectangle((0, 0), 1, 1, facecolor=INK_MUTED, edgecolor="none")
    cpd_handle = plt.Rectangle((0, 0), 1, 1, facecolor=INK_MUTED, edgecolor=INK_PRIMARY,
                                 linewidth=1.1, alpha=0.55, hatch="///")
    ax.legend([rxn_handle, cpd_handle], ["Reactions", "Metabolites"],
               loc="upper right", frameon=False, fontsize=9, labelcolor=INK_PRIMARY)

    fig.tight_layout()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out_path}")
    for label, rp, cp in zip(labels, rxn_pct, cpd_pct):
        print(f"  {label.replace(chr(10), ' ')}: reactions={rp:.1f}% compounds={cp:.1f}%")


def main() -> None:
    rows = load_rows()
    print(f"loaded {len(rows)} model rows from {RESULTS_CSV}")

    n_rxn = np.array([int(r["n_unique_reactions"]) for r in rows])
    n_cpd = np.array([int(r["n_unique_compounds"]) for r in rows])

    rxn_abs = {src: np.array([int(r[f"n_reactions_with_direction_{src}"]) for r in rows])
               for src, _, _ in SOURCES}
    cpd_abs = {src: np.array([int(r[f"n_compounds_with_energy_{src}"]) for r in rows])
               for src, _, _ in SOURCES}
    rxn_pct = {src: 100 * vals / n_rxn for src, vals in rxn_abs.items()}
    cpd_pct = {src: 100 * vals / n_cpd for src, vals in cpd_abs.items()}

    make_sample_bar_chart(
        rxn_pct, OUT_DIR / "coverage_pct_reactions_per_model.png",
        title="% of a model's reactions with a defined direction, per source\n(5,683 core models)",
        ylabel="% of model's unique reactions covered -- median ± std dev",
    )
    make_sample_bar_chart(
        rxn_abs, OUT_DIR / "coverage_abs_reactions_per_model.png",
        title="# of a model's reactions with a defined direction, per source\n(5,683 core models)",
        ylabel="# of model's unique reactions covered -- median ± std dev",
    )
    make_sample_bar_chart(
        cpd_pct, OUT_DIR / "coverage_pct_compounds_per_model.png",
        title="% of a model's compounds with a defined energy, per source\n(5,683 core models)",
        ylabel="% of model's unique compounds covered -- median ± std dev",
    )
    make_sample_bar_chart(
        cpd_abs, OUT_DIR / "coverage_abs_compounds_per_model.png",
        title="# of a model's compounds with a defined energy, per source\n(5,683 core models)",
        ylabel="# of model's unique compounds covered -- median ± std dev",
    )

    with open(SUMMARY_JSON) as fh:
        summary = json.load(fh)
    make_combined_pct_chart(summary, OUT_DIR / "coverage_pct_combined_all_models.png")


if __name__ == "__main__":
    main()
