#!/usr/bin/env python3
"""Internal-flux distributions across the ModelSEED v2 model sets, loops discounted.

Three PNGs under ``reports/fluxDistributions/figures/flux_distributions/``:

  flux_distribution.png   the ask: median and IQR of internal reaction flux,
                          one panel per model set, with the loop-flagged flux
                          drawn behind as a ghost so the discard is visible.
  flux_vs_biomass.png     the same flux divided by that model's own biomass
                          flux. Dimensionless, so all three panels share one
                          axis, and 1.0 means "running as fast as biomass".
  loop_accounting.png     what the two loop filters remove and how far they
                          move the median -- the evidence for preferring one.

Everything is read from ``results/flux_distributions/<set>/`` (summary.json +
fluxes.npz), written by scripts/build_flux_distributions.py, which is where the
method and every threshold is documented. No LP is solved here.

WHY THE X AXIS IS LOG, AND WHY IT IS NOT CLIPPED
Internal flux spans many orders of magnitude in a genome-scale model: trace
biosynthetic steps run at 1e-5 while central carbon runs at 30. On a linear
axis the whole distribution collapses into the first pixel column. The decade
range is taken from the data itself (``log_range``), shared across the three
panels so they stay comparable, and nothing is clipped -- a tail that looks
empty is empty, not cropped. Median and IQR are computed on the raw values,
never on the logs.

WHY FILTER B IS THE DEFAULT
``A`` flags a reaction that FVA can push to +-1000 in both directions; ``B``
closes every boundary reaction and flags whatever can still carry flux. On the
genome-scale panels A is a strict subset of B. On the core panel it is not:
that medium supplies every compound at -1000, so a reaction can reach 1000
without looping, and A over-flags by ~4x. B is medium-independent, so it is
the default and A is drawn beside it. See the builder's docstring.

Colour: the repo's validated categorical triple (blue/orange/aqua), assigned
in fixed order to core / GMM / auxotrophy and never cycled -- the same
assignment plot_baseline_biomass_taxonomy.py uses, so a panel keeps its colour
across figures. ``validate_palette.js --pairs all --mode light`` reports ALL
CHECKS PASS with one contrast WARN on aqua (2.74:1); the relief that WARN
obligates is the direct median/IQR label carried on every panel plus
``figure_stats.tsv`` beside the images.

Usage::

    python3 scripts/regen_figures.py flux_distributions
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
DATA_DIR = Path(os.environ.get("FLUX_DIST_DIR",
                               ANALYSIS_DIR / "results" / "flux_distributions"))
OUT_DIR = (ANALYSIS_DIR / "reports" / "fluxDistributions" / "figures"
           / "flux_distributions")

#: Panels in fixed drawing order -> (set name, title, short title, media,
#: colour). ``media`` is the subtitle prefix and is empty where the title
#: already names the medium -- repeating it there pushed the subtitle under
#: the median/IQR callout.
PANELS = [
    ("core_kegg2", "Core models", "Core",
     "complete media, 347 cpd  ·  ", "#2a78d6"),
    ("ms2_gsm", "Genome-scale · glucose minimal", "GSM · glucose min.",
     "", "#eb6834"),
    ("ms2_gsm_auxo", "Genome-scale · auxotrophy", "GSM · auxotrophy",
     "", "#1baf7a"),
]

#: Which loop filter the headline numbers use. 'b' = closed-exchange test.
DEFAULT_FILTER = "b"
FLAG_LOOP_A, FLAG_LOOP_B = 1, 2

N_BINS = 90
#: Marks the biomass reference line on both flux figures.
BIOMASS_RED = "#b5342b"
#: A flux this close to the +-1000 bound is pinned against it rather than set
#: by metabolism. On the core panel's complete medium (uptake -1000) a large
#: share of reactions saturate, which caps that panel's whole distribution --
#: the figure has to say so rather than let the median be read as biology.
BOUND = 1000.0
BOUND_EPS = 1e-2
#: Above this share pinned at the bound, the panel gets an explicit warning.
BOUND_WARN_FRAC = 0.01
#: Filter B's colour wherever the two filters are drawn against each other.
FILTER_B_COL = "#6b4fa8"

INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
SURFACE = "#fcfcfb"
NEUTRAL = "#9c9a94"
GHOST = "#c9c7c0"


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------

def load() -> dict:
    """Per panel: the summary, the flux array, the flags, the per-value biomass."""
    out = {}
    for key, _t, _s, _m, _c in PANELS:
        d = DATA_DIR / key
        if not (d / "summary.json").exists():
            raise SystemExit(
                f"missing {d}/summary.json — run first:\n"
                f"  python3 scripts/build_flux_distributions.py --sets {key}")
        npz = np.load(d / "fluxes.npz")
        flux = npz["flux"].astype(np.float64)
        flag = npz["flag"]
        out[key] = {
            "summary": json.loads((d / "summary.json").read_text()),
            "flux": flux,
            "keep_a": (flag & FLAG_LOOP_A) == 0,
            "keep_b": (flag & FLAG_LOOP_B) == 0,
            "rel": flux / npz["model_biomass"][npz["model_index"]],
        }
    return out


def keep_mask(rec: dict, which: str) -> np.ndarray:
    return rec["keep_b"] if which == "b" else rec["keep_a"]


def log_range(arrays: list) -> tuple:
    """Whole-decade span covering every value in ``arrays``.

    Shared across panels so they stay comparable, and derived from the data so
    no tail is ever cropped: if an axis looks empty out there, it is empty.
    """
    lo = min(float(a.min()) for a in arrays if a.size)
    hi = max(float(a.max()) for a in arrays if a.size)
    return math.floor(math.log10(lo)), math.ceil(math.log10(hi))


def style(ax) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRIDLINE)
    ax.tick_params(colors=INK_SECONDARY, labelsize=8.5, length=3)
    ax.grid(True, color=GRIDLINE, linewidth=0.6, alpha=0.9)
    ax.set_axisbelow(True)


def mark_quartiles(ax, vals: np.ndarray, colour: str) -> tuple:
    """Median rule + shaded IQR, drawn in log space. Returns (q1, med, q3)."""
    q1, med, q3 = (float(np.quantile(vals, q)) for q in (0.25, 0.5, 0.75))
    ax.axvspan(np.log10(q1), np.log10(q3), color=colour, alpha=0.16, zorder=1)
    ax.axvline(np.log10(med), color=INK_PRIMARY, linewidth=1.9, alpha=0.85,
               zorder=5)
    for q in (q1, q3):
        ax.axvline(np.log10(q), color=INK_PRIMARY, linewidth=1.0, alpha=0.55,
                   linestyle=":", zorder=5)
    return q1, med, q3


def decade_axis(ax, lo: int, hi: int) -> None:
    span = hi - lo
    step = 1 if span <= 8 else 2
    ticks = [t for t in range(lo, hi + 1) if t % step == 0]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"$10^{{{t}}}$" for t in ticks])
    ax.set_xlim(lo, hi)


def callout(ax, text: str, colour: str, fontsize: float = 10.2) -> None:
    """The median/IQR label, on the title row where no data can reach it.

    This is also the relief the aqua slot's contrast WARN obligates, so it is
    carried on every panel rather than only where it happens to fit.
    """
    ax.text(1.0, 1.035, text, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=fontsize, color=INK_PRIMARY, fontweight="semibold",
            bbox=dict(boxstyle="round,pad=0.38", facecolor=SURFACE,
                      edgecolor=colour, linewidth=1.3))


def below_axes_legend(ax, ncol: int = 2) -> None:
    leg = ax.legend(loc="upper right", bbox_to_anchor=(1.0, -0.135), ncol=ncol,
                    frameon=False, fontsize=8.3)
    for t in leg.get_texts():
        t.set_color(INK_SECONDARY)


# ---------------------------------------------------------------------------
# figure 1 - the distributions
# ---------------------------------------------------------------------------

def fig_distribution(D: dict, which: str, out: Path) -> list:
    fig, axes = plt.subplots(1, 3, figsize=(15.8, 7.4), facecolor=SURFACE,
                             gridspec_kw=dict(wspace=.17, left=.058,
                                              right=.988, top=.705,
                                              bottom=.235))
    lo, hi = log_range([D[k]["flux"] for k, *_ in PANELS])
    bins = np.linspace(lo, hi, N_BINS + 1)
    stats = []
    for ax, (key, title, _short, media, colour) in zip(axes, PANELS):
        style(ax)
        rec = D[key]
        keep = keep_mask(rec, which)
        kept, dropped = rec["flux"][keep], rec["flux"][~keep]

        # Ghost first: what the filter removed, so the discard is visible
        # rather than merely asserted.
        if dropped.size:
            ax.hist(np.log10(dropped), bins=bins, color=GHOST, alpha=0.95,
                    zorder=1,
                    label=f"loop-flagged, discarded (n={dropped.size:,})")
        ax.hist(np.log10(kept), bins=bins, color=colour, alpha=0.88,
                edgecolor=SURFACE, linewidth=0.25, zorder=2,
                label=f"counted (n={kept.size:,})")
        q1, med, q3 = mark_quartiles(ax, kept, colour)

        # The comparison the figure exists to make.
        bio_med = rec["summary"]["biomass_flux"]["q50"]
        ax.axvline(np.log10(bio_med), color=BIOMASS_RED, linewidth=1.8,
                   linestyle="--", zorder=6)
        ax.text(np.log10(bio_med), 0.55,
                f" median biomass flux {bio_med:.3g} ",
                transform=ax.get_xaxis_transform(), rotation=90,
                ha="right", va="center", fontsize=8.2, color=BIOMASS_RED,
                fontweight="semibold", zorder=7,
                bbox=dict(boxstyle="round,pad=0.2", facecolor=SURFACE,
                          edgecolor="none", alpha=0.88))

        decade_axis(ax, lo, hi)
        ax.set_title(title, color=INK_PRIMARY, fontsize=11.5,
                     fontweight="semibold", pad=40, loc="left")
        # One line only: a second line would run under the callout, and the
        # counted-n is already stated in the legend.
        ax.text(0, 1.035,
                f"{media}{rec['summary']['n_growers']:,} growing models",
                transform=ax.transAxes, color=INK_SECONDARY, fontsize=8.5,
                va="bottom")
        callout(ax, f"median {med:.3g}\nIQR {q1:.3g} – {q3:.3g}", colour)
        ax.set_xlabel("internal reaction flux |v|  (mmol gDW⁻¹ h⁻¹, log scale)",
                      color=INK_SECONDARY, fontsize=9)
        below_axes_legend(ax)

        # Saturation check. A panel whose fluxes pile up on the +-1000 bound is
        # reporting its bound, not its metabolism, and must say so on its face.
        at_bound = float((kept >= BOUND - BOUND_EPS).mean())
        if at_bound > BOUND_WARN_FRAC:
            ax.text(0, -.285,
                    f"⚠  {100 * at_bound:.0f}% of these fluxes sit exactly on "
                    f"the ±1000 bound.\nUptake is −1000 on this medium, so the "
                    f"panel is saturated:\nread its median as the bound, not "
                    f"as biology.",
                    transform=ax.transAxes, ha="left", va="top",
                    fontsize=8.4, color=BIOMASS_RED, fontweight="semibold",
                    linespacing=1.45)
        stats.append((key, q1, med, q3, kept.size, dropped.size, bio_med,
                      at_bound))
    axes[0].set_ylabel("reaction–model pairs", color=INK_SECONDARY, fontsize=9)

    fig.suptitle("How hard do the internal reactions run? — median and IQR of "
                 "flux across every model, thermodynamic loops discounted",
                 color=INK_PRIMARY, fontsize=14, fontweight="semibold",
                 x=.058, ha="left", y=.975)
    fig.text(.058, .936,
             "Parsimonious FBA (pFBA) at maximum biomass, every growing model "
             "in the set pooled. Only reactions actually carrying flux are "
             "counted; exchange, demand, sink and biomass reactions are not "
             "internal and are excluded.\n"
             "Grey = fluxes discarded because the reaction still carries flux "
             "with every boundary closed, i.e. it sits on a thermodynamically "
             "infeasible cycle whose magnitude is set by the ±1000 bound "
             "rather than by biology.",
             color=INK_MUTED, fontsize=9, va="top", linespacing=1.45)
    fig.savefig(out, dpi=170, facecolor=SURFACE)
    plt.close(fig)
    return stats


# ---------------------------------------------------------------------------
# figure 2 - flux in units of the model's own biomass flux
# ---------------------------------------------------------------------------

def fig_vs_biomass(D: dict, which: str, out: Path) -> list:
    fig = plt.figure(figsize=(15.8, 8.6), facecolor=SURFACE)
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.95], hspace=.62,
                          wspace=.17, left=.058, right=.988, top=.78,
                          bottom=.075)
    rels = {k: D[k]["rel"][keep_mask(D[k], which)] for k, *_ in PANELS}
    lo, hi = log_range(list(rels.values()))
    bins = np.linspace(lo, hi, N_BINS + 1)
    stats = []
    for col, (key, title, _short, _media, colour) in enumerate(PANELS):
        ax = fig.add_subplot(gs[0, col])
        style(ax)
        rel = rels[key]
        ax.hist(np.log10(rel), bins=bins, color=colour, alpha=0.88,
                edgecolor=SURFACE, linewidth=0.25, zorder=2)
        q1, med, q3 = mark_quartiles(ax, rel, colour)
        ax.axvline(0.0, color=BIOMASS_RED, linewidth=1.8, linestyle="--",
                   zorder=6)
        frac = float((rel >= 1.0).mean())
        decade_axis(ax, lo, hi)
        ax.set_title(title, color=INK_PRIMARY, fontsize=11, pad=30,
                     fontweight="semibold", loc="left")
        ax.text(0, 1.035, f"{100 * frac:.1f}% run at or above biomass",
                transform=ax.transAxes, color=BIOMASS_RED, fontsize=8.8,
                va="bottom", fontweight="semibold")
        callout(ax, f"median {med:.3g}×\nIQR {q1:.3g} – {q3:.3g}×", colour,
                fontsize=9.6)
        ax.set_xlabel("internal flux ÷ that model's biomass flux",
                      color=INK_SECONDARY, fontsize=9)
        if col == 0:
            ax.set_ylabel("reaction–model pairs", color=INK_SECONDARY,
                          fontsize=9)
        stats.append((key, q1, med, q3, frac))

    # ---- pooled ECDF: dimensionless, so all three legitimately share it ----
    ax = fig.add_subplot(gs[1, :])
    style(ax)
    for key, title, _short, _media, colour in PANELS:
        rel = np.sort(rels[key])
        y = np.arange(1, rel.size + 1) / rel.size
        step = max(1, rel.size // 4000)        # thin for file size, not shape
        ax.plot(np.log10(rel[::step]), y[::step], color=colour, linewidth=2.2,
                label=f"{title}  (n={rel.size:,})")
    ax.axvline(0.0, color=BIOMASS_RED, linewidth=1.8, linestyle="--", zorder=6)
    ax.text(0.04, .035, " internal flux = biomass flux",
            transform=ax.get_xaxis_transform(), ha="left", va="bottom",
            fontsize=8.6, color=BIOMASS_RED, fontweight="semibold")
    decade_axis(ax, lo, hi)
    ax.set_ylim(0, 1.04)
    ax.set_xlabel("internal flux ÷ that model's own biomass flux  "
                  "(dimensionless — the one view all three panels share)",
                  color=INK_SECONDARY, fontsize=9.5)
    ax.set_ylabel("cumulative fraction\nof reaction fluxes",
                  color=INK_SECONDARY, fontsize=9)
    ax.set_title("Pooled — what fraction of internal reactions run slower, or "
                 "faster, than the biomass reaction they feed",
                 color=INK_PRIMARY, fontsize=11, fontweight="semibold",
                 pad=10, loc="left")
    leg = ax.legend(loc="upper left", frameon=True, fontsize=9,
                    facecolor=SURFACE, edgecolor=GRIDLINE)
    for t, (_k, _ti, _s, _m, colour) in zip(leg.get_texts(), PANELS):
        t.set_color(colour)

    fig.suptitle("Is a biomass flux of ~1 \"low\"? — internal flux measured in "
                 "units of each model's own biomass flux",
                 color=INK_PRIMARY, fontsize=14, fontweight="semibold",
                 x=.058, ha="left", y=.975)
    fig.text(.058, .936,
             "Dividing by the model's own biomass flux removes the media "
             "scale, so the three panels can be read against one another and "
             "against the dashed line at 1.0.\n"
             "A value of 10 means that reaction carries ten times the flux of "
             "the biomass reaction it supplies. Loop-flagged fluxes are "
             "already discarded.",
             color=INK_MUTED, fontsize=9, va="top", linespacing=1.45)
    fig.savefig(out, dpi=170, facecolor=SURFACE)
    plt.close(fig)
    return stats


# ---------------------------------------------------------------------------
# figure 3 - what the two loop filters do
# ---------------------------------------------------------------------------

def fig_loops(D: dict, out: Path) -> list:
    fig, (axL, axR) = plt.subplots(
        1, 2, figsize=(14.8, 6.8), facecolor=SURFACE,
        gridspec_kw=dict(wspace=.26, left=.115, right=.985, top=.70,
                         bottom=.145))
    style(axL)
    style(axR)
    y = np.arange(len(PANELS))[::-1]
    h = .34
    rows = []

    # ---- left: share of flux values each filter discards ----
    a_share, b_share = [], []
    for key, *_ in PANELS:
        rec = D[key]
        n = rec["flux"].size
        a_share.append(100 * float((~rec["keep_a"]).sum()) / n)
        b_share.append(100 * float((~rec["keep_b"]).sum()) / n)
    axL.barh(y + h / 2, a_share, height=h, color=NEUTRAL, edgecolor=SURFACE,
             linewidth=1.0, label="A   ±1000 in both directions")
    axL.barh(y - h / 2, b_share, height=h, color=FILTER_B_COL,
             edgecolor=SURFACE, linewidth=1.0,
             label="B   any flux with boundaries closed")
    for yy, a, b in zip(y, a_share, b_share):
        axL.text(a + .4, yy + h / 2, f"{a:.1f}%", va="center", fontsize=8.6,
                 color=INK_SECONDARY)
        axL.text(b + .4, yy - h / 2, f"{b:.1f}%", va="center", fontsize=8.6,
                 color=INK_SECONDARY)
    axL.set_yticks(y, [s for _k, _t, s, _m, _c in PANELS], fontsize=9.6,
                   color=INK_PRIMARY)
    axL.set_xlim(0, max(a_share + b_share) * 1.2)
    axL.set_xlabel("% of internal reaction fluxes flagged as loop",
                   color=INK_SECONDARY, fontsize=9.5)
    axL.set_title("What each filter discards", color=INK_PRIMARY,
                  fontsize=11.5, fontweight="semibold", pad=14, loc="left")
    leg = axL.legend(loc="lower right", frameon=True, fontsize=8.8,
                     facecolor=SURFACE, edgecolor=GRIDLINE)
    for t in leg.get_texts():
        t.set_color(INK_SECONDARY)

    # ---- right: how far the median moves, each panel on its own row ----
    for i, (key, _title, _short, _media, colour) in enumerate(PANELS):
        rec = D[key]
        meds = [float(np.median(rec["flux"])),
                float(np.median(rec["flux"][rec["keep_a"]])),
                float(np.median(rec["flux"][rec["keep_b"]]))]
        ratios = [m / meds[0] for m in meds]
        yy = y[i]
        # A connecting rule makes the three marks read as one panel's journey.
        axR.plot(ratios, [yy] * 3, "-", color=colour, linewidth=1.3, alpha=.5,
                 zorder=2)
        # Three distinct vertical offsets, not two: on the core row B lands
        # within 1% of the unfiltered median, so sharing an offset with
        # "no filter" would overprint it.
        for j, (lbl, mk, dy) in enumerate((("no filter", "o", 13),
                                           ("A", "s", -21),
                                           ("B", "D", 30))):
            axR.plot([ratios[j]], [yy], mk, color=colour, markersize=9,
                     markeredgecolor=SURFACE, markeredgewidth=1.5, zorder=3)
            axR.annotate(f"{lbl}  {meds[j]:.3g}", (ratios[j], yy),
                         textcoords="offset points", xytext=(0, dy),
                         ha="center", fontsize=7.9, color=INK_SECONDARY)
        subset = rec["summary"]["loops"]["a_is_subset_of_b"]
        rows.append((key, meds[0], meds[1], meds[2], subset))
        # Inside the axes, bottom-left: outside it would run off the figure.
        axR.text(.015, yy - .42,
                 "A ⊆ B — sound here, just incomplete" if subset else
                 "A ⊄ B — flags non-loops: this medium supplies 1000 units",
                 transform=axR.get_yaxis_transform(), ha="left", va="center",
                 fontsize=8.1, style="italic",
                 color=INK_MUTED if subset else BIOMASS_RED)
    axR.axvline(1.0, color=NEUTRAL, linewidth=1.3, linestyle="--", zorder=1)
    axR.set_yticks(y, ["" for _ in PANELS])
    axR.set_ylim(-0.75, len(PANELS) - 0.35)
    # Room for the centred labels on the outermost marks.
    lo_r = min(r for _k, m0, ma, mb, _s in rows for r in (ma / m0, mb / m0))
    axR.set_xlim(lo_r - 0.08, 1.09)
    axR.set_xlabel("median flux ÷ the unfiltered median",
                   color=INK_SECONDARY, fontsize=9.5)
    axR.set_title("How far the median moves", color=INK_PRIMARY, fontsize=11.5,
                  fontweight="semibold", pad=14, loc="left")

    fig.suptitle("The two loop filters are not interchangeable — which one is "
                 "correct depends on the medium",
                 color=INK_PRIMARY, fontsize=14, fontweight="semibold",
                 x=.012, ha="left", y=.975)
    fig.text(.012, .925,
             "A flags a reaction FVA can drive to +1000 AND −1000. B closes "
             "every exchange, demand and sink, then flags anything that can "
             "still carry flux — the standard infeasible-cycle test.\n"
             "On the genome-scale panels (uptake −10) A is a strict subset of "
             "B: sound, but it misses loop members capped below 1000. On the "
             "core panel the complete medium supplies every compound at "
             "−1000, so a reaction\n"
             "can reach 1000 without looping and A over-flags. B is "
             "medium-independent, and is the filter the other two figures use.",
             color=INK_MUTED, fontsize=9, va="top", linespacing=1.45)
    fig.savefig(out, dpi=170, facecolor=SURFACE)
    plt.close(fig)
    return rows


# ---------------------------------------------------------------------------

def write_stats(D: dict, dist, rel, loops, which: str, out_dir: Path) -> None:
    """Every number the three figures assert, so prose can be rechecked."""
    with open(out_dir / "figure_stats.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["section", "panel", "key", "value"])
        w.writerow(["config", "-", "headline_filter", which])
        for key, q1, med, q3, nk, nd, bio, at_bound in dist:
            for k, v in (("flux_q25", q1), ("flux_median", med),
                         ("flux_q75", q3), ("n_counted", nk),
                         ("n_discarded", nd), ("median_biomass_flux", bio),
                         ("frac_pinned_at_1000_bound", at_bound)):
                w.writerow(["distribution", key, k, f"{v:.6g}"])
        for key, q1, med, q3, frac in rel:
            for k, v in (("rel_q25", q1), ("rel_median", med),
                         ("rel_q75", q3), ("frac_at_or_above_biomass", frac)):
                w.writerow(["vs_biomass", key, k, f"{v:.6g}"])
        for key, m0, ma, mb, subset in loops:
            for k, v in (("median_unfiltered", m0), ("median_filter_a", ma),
                         ("median_filter_b", mb)):
                w.writerow(["loops", key, k, f"{v:.6g}"])
            w.writerow(["loops", key, "a_is_subset_of_b", subset])
            for k, v in D[key]["summary"]["loops"].items():
                w.writerow(["loops", key, k, v])
        for key, *_ in PANELS:
            s = D[key]["summary"]
            for k in ("n_models", "n_growers", "n_errors", "n_flux_values"):
                w.writerow(["panel", key, k, s[k]])


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--filter", choices=("a", "b"), default=DEFAULT_FILTER,
                    help="loop filter the headline numbers use (default: b)")
    ap.add_argument("--out-subdir", default=None)
    args = ap.parse_args()

    out_dir = OUT_DIR if args.out_subdir is None else (
        ANALYSIS_DIR / "reports" / "figures" / args.out_subdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    D = load()
    dist = fig_distribution(D, args.filter, out_dir / "flux_distribution.png")
    rel = fig_vs_biomass(D, args.filter, out_dir / "flux_vs_biomass.png")
    loops = fig_loops(D, out_dir / "loop_accounting.png")
    write_stats(D, dist, rel, loops, args.filter, out_dir)
    for name in ("flux_distribution.png", "flux_vs_biomass.png",
                 "loop_accounting.png", "figure_stats.tsv"):
        print(f"wrote {out_dir}/{name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
