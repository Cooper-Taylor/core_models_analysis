#!/usr/bin/env python3
"""Baseline biomass distributions and genome taxonomy for the ModelSEED v2 model sets.

Draws three PNGs under
``reports/baselineBiomassTaxonomy/figures/baseline_biomass_taxonomy/``:

  biomass_distribution.png   per-panel histogram of biomass flux among growers,
                             plus a median-normalised ECDF that compares the
                             three distributions' SHAPE on one axis.
  taxonomy_composition.png   genomes per phylum, core set vs genome-scale set,
                             beside the baseline grower fraction per phylum.
  biomass_by_phylum.png      per-panel median biomass flux by phylum, with the
                             interquartile range drawn as a bar.

All numbers come from ``results/baseline_biomass_taxonomy/`` -- written by
scripts/build_baseline_biomass_taxonomy.py, which is where the provenance of
every input is documented. Nothing is recomputed here and no LP is solved; this
script only draws what that build step already wrote.

WHY THE THREE PANELS NEVER SHARE AN X AXIS
Each panel is a different model set on a different media, so their absolute
fluxes live on different scales (core medians ~52, genome-scale ~1.2). Plotting
them against a common axis would compress two of the three into a spike and
imply a comparison the data does not support. Each histogram therefore gets its
own x axis, and the only cross-panel comparison drawn is the ECDF, where each
panel's flux is divided by that panel's own median so the curves describe
dispersion rather than magnitude.

  core        5,683 core models, ModelSEED complete media (347 compounds)
  full (GMM)  5,420 genome-scale models, glucose minimal media
  full (auxo) 5,420 genome-scale models, auxotrophy media

"Baseline" is the models exactly as they sit on disk: their own stored bounds,
their own biomass reaction, the media they were gap-filled against. No
direction map, thermodynamic source, heuristic or gap-fill is applied.

Colour: dataviz categorical slots 1-3 (blue/orange/aqua), assigned in fixed
order to core / full-GMM / full-auxo and never cycled. validate_palette.js
reports PASS on this triple under ``--pairs all``, light mode, with one
contrast WARN on aqua (2.74:1) -- the relief rule is met by the direct value
labels on every bar and by ``figure_stats.tsv`` beside the images.

Usage::

    python3 scripts/regen_figures.py baseline_biomass_taxonomy
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
DATA_DIR = ANALYSIS_DIR / "results" / "baseline_biomass_taxonomy"
OUT_DIR = (ANALYSIS_DIR / "reports" / "baselineBiomassTaxonomy" / "figures"
           / "baseline_biomass_taxonomy")

#: Panels in fixed drawing order -> (short label, media line, categorical slot).
PANELS = [
    ("core", "Core models", "complete media, 347 cpd", "#2a78d6"),
    ("full_gmm", "Genome-scale", "glucose minimal media", "#eb6834"),
    ("full_auxo", "Genome-scale", "auxotrophy media", "#1baf7a"),
]
PANEL_TITLE = {"core": "Core models",
               "full_gmm": "Genome-scale · glucose minimal",
               "full_auxo": "Genome-scale · auxotrophy"}

#: How many phyla to draw before folding the rest into "Other". Chosen so every
#: drawn phylum has >= 20 genomes in the union of the two sets.
TOP_PHYLA = 14
N_BINS = 60

INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
SURFACE = "#fcfcfb"
#: Neutral used where a mark carries no categorical identity ("Other", rules).
NEUTRAL = "#9c9a94"


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------


def load() -> tuple[dict, dict, list, dict]:
    with open(DATA_DIR / "summary.json") as fh:
        summary = json.load(fh)
    flux: dict[str, list[float]] = defaultdict(list)
    phylum_of: dict[str, str] = {}
    per_panel_phylum: dict[tuple, list] = defaultdict(list)
    with open(DATA_DIR / "per_model.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["grows"] != "True":
                continue
            val = float(row["growth"])
            flux[row["panel"]].append(val)
            phylum_of[row["genome_id"]] = row["phylum"]
            per_panel_phylum[(row["panel"], row["phylum"])].append(val)
    with open(DATA_DIR / "taxonomy_phylum.tsv") as fh:
        phyla = list(csv.DictReader(fh, delimiter="\t"))
    return summary, flux, phyla, per_panel_phylum


def style(ax) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRIDLINE)
    ax.tick_params(colors=INK_SECONDARY, labelsize=8.5, length=3)
    ax.grid(True, color=GRIDLINE, linewidth=0.6, alpha=0.9)
    ax.set_axisbelow(True)


# ---------------------------------------------------------------------------
# figure 1 - biomass distributions
# ---------------------------------------------------------------------------


def fig_distribution(summary: dict, flux: dict, out: Path) -> None:
    fig = plt.figure(figsize=(13.2, 7.4), facecolor=SURFACE)
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.92], hspace=0.50,
                          wspace=0.24, left=.062, right=.985, top=.845, bottom=.085)

    for col, (key, _lbl, media, colour) in enumerate(PANELS):
        ax = fig.add_subplot(gs[0, col])
        style(ax)
        vals = np.array(flux[key])
        rec = summary["panels"][key]
        ax.hist(vals, bins=N_BINS, color=colour, alpha=0.85,
                edgecolor=SURFACE, linewidth=0.4)
        for q, ls, lw in ((rec["median"], "-", 1.8), (rec["q1"], ":", 1.2),
                          (rec["q3"], ":", 1.2)):
            ax.axvline(q, color=INK_PRIMARY, linestyle=ls, linewidth=lw, alpha=.75)
        ax.set_title(PANEL_TITLE[key], color=INK_PRIMARY, fontsize=11.5,
                     fontweight="semibold", pad=15, loc="left")
        ax.text(0, 1.015,
                f"{rec['n_growers']:,} of {rec['n_models']:,} grow "
                f"({100 * rec['frac_growers']:.1f}%)  ·  {media}",
                transform=ax.transAxes, color=INK_SECONDARY, fontsize=8.6)
        # Direct labels satisfy the relief rule for the low-contrast slot.
        ax.text(.975, .95,
                f"median {rec['median']:.3g}\nIQR {rec['q1']:.3g}–{rec['q3']:.3g}\n"
                f"max {rec['max']:.3g}",
                transform=ax.transAxes, ha="right", va="top", fontsize=8.8,
                color=INK_PRIMARY, linespacing=1.5,
                bbox=dict(boxstyle="round,pad=0.42", facecolor=SURFACE,
                          edgecolor=GRIDLINE, linewidth=.8))
        ax.set_xlabel("biomass flux (mmol gDW⁻¹ h⁻¹)", color=INK_SECONDARY,
                      fontsize=9)
        if col == 0:
            ax.set_ylabel("models", color=INK_SECONDARY, fontsize=9)

    # ---- median-normalised ECDF: the only legitimate cross-panel view ----
    ax = fig.add_subplot(gs[1, :])
    style(ax)
    for key, _lbl, media, colour in PANELS:
        vals = np.sort(np.array(flux[key]))
        med = summary["panels"][key]["median"]
        x = vals / med
        y = np.arange(1, len(x) + 1) / len(x)
        ax.plot(x, y, color=colour, linewidth=2.0,
                label=f"{PANEL_TITLE[key]}  (n={len(x):,})")
        # Direct-label each curve at its own right-hand end. The three curves
        # are pinned together at the median by construction and only separate
        # in the upper tail, so the tail is the one place labels do not collide.
        ax.annotate(PANEL_TITLE[key], (x[-1], 1.0),
                    textcoords="offset points", xytext=(-4, 7),
                    ha="right", color=colour, fontsize=8.8,
                    fontweight="semibold")
    ax.axvline(1.0, color=NEUTRAL, linewidth=1.2, linestyle="--")
    ax.text(1.0, .035, " panel median", transform=ax.get_xaxis_transform(),
            ha="left", va="bottom", fontsize=8.2, color=INK_MUTED)
    ax.set_xlim(0, 2.6)
    ax.set_ylim(0, 1.06)
    ax.set_xlabel("biomass flux ÷ that panel's own median  "
                  "(dimensionless — removes the media scale difference)",
                  color=INK_SECONDARY, fontsize=9.5)
    ax.set_ylabel("cumulative fraction\nof growing models", color=INK_SECONDARY,
                  fontsize=9)
    ax.set_title("Shape comparison — how tightly each panel's growers cluster "
                 "around their own median",
                 color=INK_PRIMARY, fontsize=11, fontweight="semibold",
                 pad=10, loc="left")
    leg = ax.legend(loc="lower right", frameon=True, fontsize=8.8,
                    facecolor=SURFACE, edgecolor=GRIDLINE)
    for text in leg.get_texts():
        text.set_color(INK_SECONDARY)

    fig.suptitle("Baseline biomass flux — ModelSEED v2 core and genome-scale "
                 "models, as they sit on disk",
                 color=INK_PRIMARY, fontsize=14, fontweight="semibold",
                 x=.062, ha="left", y=.975)
    fig.text(.062, .935,
             "No direction map, thermodynamic source, heuristic or gap-fill "
             "applied — each model's own stored bounds, its own biomass "
             "reaction, and the media it was gap-filled against.\n"
             "The three panels sit on different media, so their flux scales "
             "are not comparable: only the bottom axis compares across them.",
             color=INK_MUTED, fontsize=9, va="top", linespacing=1.45)
    fig.savefig(out, dpi=170, facecolor=SURFACE)
    plt.close(fig)


# ---------------------------------------------------------------------------
# figure 2 - taxonomy composition and grower rate
# ---------------------------------------------------------------------------


def top_phyla(phyla: list) -> list:
    ranked = [p for p in phyla if p["phylum"] != "Unassigned"]
    ranked.sort(key=lambda p: -int(p["n_genomes_union"]))
    return ranked[:TOP_PHYLA]


def fig_taxonomy(summary: dict, phyla: list, out: Path) -> None:
    rows = top_phyla(phyla)
    rest = [p for p in phyla if p not in rows]
    other = {
        "phylum": f"Other ({len(rest)} phyla)",
        "n_genomes_union": sum(int(p["n_genomes_union"]) for p in rest),
        "n_core": sum(int(p["n_core"]) for p in rest),
        "n_full": sum(int(p["n_full"]) for p in rest),
        "core_growers": sum(int(p["core_growers"]) for p in rest),
        "full_gmm_growers": sum(int(p["full_gmm_growers"]) for p in rest),
        "full_auxo_growers": sum(int(p["full_auxo_growers"]) for p in rest),
    }
    rows = rows + [other]
    labels = [r["phylum"] for r in rows]
    y = np.arange(len(rows))[::-1]

    fig, (axL, axR) = plt.subplots(
        1, 2, figsize=(14.6, 7.9), facecolor=SURFACE,
        gridspec_kw=dict(width_ratios=[1.0, 1.0], wspace=.30,
                         left=.155, right=.985, top=.825, bottom=.125))
    style(axL)
    style(axR)

    # ---- left: genome counts, core vs genome-scale ----
    h = .38
    core_n = [int(r["n_core"]) for r in rows]
    full_n = [int(r["n_full"]) for r in rows]
    axL.barh(y + h / 2, core_n, height=h, color="#2a78d6", label="Core set",
             edgecolor=SURFACE, linewidth=1.0)
    axL.barh(y - h / 2, full_n, height=h, color=NEUTRAL,
             label="Genome-scale set", edgecolor=SURFACE, linewidth=1.0)
    for yy, nc, nf in zip(y, core_n, full_n):
        axL.text(nc + 22, yy + h / 2, f"{nc:,}", va="center", fontsize=8.2,
                 color=INK_SECONDARY)
        axL.text(nf + 22, yy - h / 2, f"{nf:,}", va="center", fontsize=8.2,
                 color=INK_SECONDARY)
    axL.set_yticks(y, labels, fontsize=9.4, color=INK_PRIMARY)
    axL.set_xlabel("genomes", color=INK_SECONDARY, fontsize=9.5)
    axL.set_xlim(0, max(core_n) * 1.14)
    axL.set_title("Genomes per phylum", color=INK_PRIMARY, fontsize=11.5,
                  fontweight="semibold", pad=24, loc="left")
    axL.text(0, 1.018,
             f"{summary['genome_sets']['core']:,} core · "
             f"{summary['genome_sets']['full']:,} genome-scale · "
             f"{summary['genome_sets']['shared']:,} shared",
             transform=axL.transAxes, color=INK_SECONDARY, fontsize=8.8)
    leg = axL.legend(loc="lower right", frameon=True, fontsize=9,
                     facecolor=SURFACE, edgecolor=GRIDLINE)
    for t in leg.get_texts():
        t.set_color(INK_SECONDARY)

    # ---- right: grower fraction, three panels ----
    series = [("core_growers", "n_core", "Core · complete media", "#2a78d6"),
              ("full_gmm_growers", "n_full", "Genome-scale · glucose min.", "#eb6834"),
              ("full_auxo_growers", "n_full", "Genome-scale · auxotrophy", "#1baf7a")]
    hh = .26
    for i, (gkey, nkey, label, colour) in enumerate(series):
        offs = (1 - i) * hh
        fracs, present = [], []
        for r in rows:
            denom = int(r[nkey])
            fracs.append(100 * int(r[gkey]) / denom if denom else np.nan)
            present.append(denom > 0)
        axR.barh(y + offs, [f if p else 0 for f, p in zip(fracs, present)],
                 height=hh * .92, color=colour, label=label,
                 edgecolor=SURFACE, linewidth=.8)
        for yy, f, p in zip(y, fracs, present):
            if not p:
                # Say it once per phylum, not once per absent series.
                if i == 1:
                    axR.text(1.5, yy, "— phylum absent from the genome-scale set",
                             va="center", fontsize=7.6, color=INK_MUTED,
                             style="italic")
            else:
                axR.text(f + 1.2, yy + offs, f"{f:.0f}%", va="center",
                         fontsize=7.8, color=INK_SECONDARY)
    axR.set_yticks(y, ["" for _ in labels])
    axR.set_xlim(0, 118)
    axR.set_xticks([0, 25, 50, 75, 100])
    axR.set_xlabel("% of that phylum's models with non-zero biomass flux",
                   color=INK_SECONDARY, fontsize=9.5)
    axR.set_title("Baseline grower rate per phylum", color=INK_PRIMARY,
                  fontsize=11.5, fontweight="semibold", pad=24, loc="left")
    axR.text(0, 1.018,
             "same genomes, three on-disk baselines",
             transform=axR.transAxes, color=INK_SECONDARY, fontsize=8.8)
    # Below the axes: every row carries bars out to ~100%, so any in-axes
    # placement would sit on top of data.
    leg = axR.legend(loc="upper right", bbox_to_anchor=(1.0, -0.075), ncol=3,
                     frameon=False, fontsize=8.8)
    for t in leg.get_texts():
        t.set_color(INK_SECONDARY)

    fig.suptitle("Taxonomic distribution of the ModelSEED v2 genomes, and how "
                 "baseline growth tracks it",
                 color=INK_PRIMARY, fontsize=14, fontweight="semibold",
                 x=.012, ha="left", y=.972)
    sk = {r["superkingdom"]: r for r in _superkingdom_rows()}
    bact, arch = sk.get("Bacteria", {}), sk.get("Archaea", {})
    fig.text(.012, .935,
             f"NCBI lineages for all {summary['genome_sets']['union']:,} distinct "
             f"genomes: {int(bact.get('n_genomes_union', 0)):,} Bacteria, "
             f"{int(arch.get('n_genomes_union', 0)):,} Archaea, across "
             f"{summary['taxonomy']['phylum']['n_distinct']} phyla, "
             f"{summary['taxonomy']['genus']['n_distinct']:,} genera and "
             f"{summary['taxonomy']['species']['n_distinct']:,} species.\n"
             f"Not one of the {int(arch.get('n_core', 0))} archaeal CORE models "
             f"carries flux to biomass on complete media; the genome-scale "
             f"versions of those same genomes nearly all do. Cyanobacteria were "
             f"dropped from the genome-scale set by the manuscript.",
             color=INK_MUTED, fontsize=9, va="top", linespacing=1.45)
    fig.savefig(out, dpi=170, facecolor=SURFACE)
    plt.close(fig)


def _superkingdom_rows() -> list:
    with open(DATA_DIR / "taxonomy_superkingdom.tsv") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------------------------------------------------------------------
# figure 3 - biomass by phylum
# ---------------------------------------------------------------------------


def fig_by_phylum(summary: dict, phyla: list, per_panel_phylum: dict,
                  out: Path) -> None:
    rows = top_phyla(phyla)
    labels = [r["phylum"] for r in rows]
    y = np.arange(len(rows))[::-1]

    fig, axes = plt.subplots(
        1, 3, figsize=(14.6, 7.2), facecolor=SURFACE,
        gridspec_kw=dict(wspace=.14, left=.155, right=.99, top=.835, bottom=.095))

    for ax, (key, _lbl, media, colour) in zip(axes, PANELS):
        style(ax)
        meds, los, his, ns = [], [], [], []
        for r in rows:
            vals = np.array(per_panel_phylum.get((key, r["phylum"]), []))
            ns.append(len(vals))
            if len(vals) == 0:
                meds.append(np.nan); los.append(np.nan); his.append(np.nan)
                continue
            meds.append(float(np.median(vals)))
            los.append(float(np.percentile(vals, 25)))
            his.append(float(np.percentile(vals, 75)))
        for yy, lo, hi, md, n in zip(y, los, his, meds, ns):
            if n == 0:
                ax.text(0.02, yy, "no growing model", va="center", fontsize=7.6,
                        color=INK_MUTED, style="italic",
                        transform=ax.get_yaxis_transform())
                continue
            ax.barh(yy, hi - lo, left=lo, height=.5, color=colour, alpha=.42,
                    edgecolor=SURFACE, linewidth=.8)
            ax.plot([md], [yy], "o", color=colour, markersize=7,
                    markeredgecolor=SURFACE, markeredgewidth=1.4, zorder=3)
            ax.text(hi, yy + .42, f"{md:.3g}", va="bottom", ha="right",
                    fontsize=7.4, color=INK_SECONDARY)
        ax.set_yticks(y, labels if ax is axes[0] else ["" for _ in labels],
                      fontsize=9.4, color=INK_PRIMARY)
        rec = summary["panels"][key]
        ax.axvline(rec["median"], color=INK_PRIMARY, linestyle="--",
                   linewidth=1.1, alpha=.55)
        ax.set_title(PANEL_TITLE[key], color=INK_PRIMARY, fontsize=11,
                     fontweight="semibold", pad=17, loc="left")
        ax.text(0, 1.012, f"{media} · panel median {rec['median']:.3g} (dashed)",
                transform=ax.transAxes, color=INK_SECONDARY, fontsize=8.3)
        ax.set_xlabel("biomass flux (mmol gDW⁻¹ h⁻¹)", color=INK_SECONDARY,
                      fontsize=9)
        ax.set_xlim(left=0)

    fig.suptitle("Baseline biomass flux by phylum — median (dot) and "
                 "interquartile range (bar), growers only",
                 color=INK_PRIMARY, fontsize=14, fontweight="semibold",
                 x=.012, ha="left", y=.955)
    fig.text(.012, .895,
             "Each panel keeps its own x axis: the three are different model "
             "sets on different media and their flux scales are not "
             "comparable. Compare the spread of phyla WITHIN a panel, never a "
             "value ACROSS panels.",
             color=INK_MUTED, fontsize=9)
    fig.savefig(out, dpi=170, facecolor=SURFACE)
    plt.close(fig)


# ---------------------------------------------------------------------------


def write_stats(summary: dict, phyla: list, out_dir: Path) -> None:
    """The numbers the figures assert, in a form prose can be rechecked against."""
    with open(out_dir / "figure_stats.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["section", "key", "value"])
        for key, rec in summary["panels"].items():
            for field in ("model_set", "media", "n_models", "n_growers",
                          "frac_growers", "min", "q1", "median", "mean", "q3",
                          "max", "sd", "cv"):
                w.writerow([f"panel:{key}", field, rec[field]])
        for field, val in summary["genome_sets"].items():
            w.writerow(["genomes", field, val])
        for rank, rec in summary["taxonomy"].items():
            w.writerow([f"taxonomy:{rank}", "n_distinct", rec["n_distinct"]])
            w.writerow([f"taxonomy:{rank}", "n_unassigned_genomes",
                        rec["n_unassigned_genomes"]])
        for r in top_phyla(phyla):
            for field in ("n_genomes_union", "n_core", "n_full",
                          "core_grower_frac", "full_gmm_grower_frac",
                          "full_auxo_grower_frac"):
                w.writerow([f"phylum:{r['phylum']}", field, r[field]])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-subdir", default=None,
                    help="write PNGs to reports/figures/<NAME>/ instead of the default")
    args = ap.parse_args()
    out_dir = OUT_DIR if args.out_subdir is None else (
        ANALYSIS_DIR / "reports" / "figures" / args.out_subdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    summary, flux, phyla, per_panel_phylum = load()
    fig_distribution(summary, flux, out_dir / "biomass_distribution.png")
    fig_taxonomy(summary, phyla, out_dir / "taxonomy_composition.png")
    fig_by_phylum(summary, phyla, per_panel_phylum,
                  out_dir / "biomass_by_phylum.png")
    write_stats(summary, phyla, out_dir)
    print(f"wrote {out_dir}/biomass_distribution.png")
    print(f"wrote {out_dir}/taxonomy_composition.png")
    print(f"wrote {out_dir}/biomass_by_phylum.png")
    print(f"wrote {out_dir}/figure_stats.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
