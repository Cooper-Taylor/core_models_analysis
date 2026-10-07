#!/usr/bin/env python3
"""Pareto coverage curve: why classifying the most-cited enzymes first
front-loads reaction coverage.

Reads results/enzyme_functionality/ec_frequency_coverage.tsv (ec_rank 1..3883,
cumulative_unique_reactions_covered, is_classified 0/1), written by
build_enzyme_functionality_outputs.py. This is the THEORETICAL curve over all
3,883 unique complete-EC enzymes in frequency order -- it exists regardless of
how many have actually been classified yet. The classified prefix (is_classified
== 1) is shaded in the status-good green; the remaining queue is the neutral
sequential blue used elsewhere in this repo for "not yet resolved" magnitude.

Data: ModelSEED dev @ 078a395f (/scratch/ctaylor/tmp/devsnap_078a395f).
"""
from __future__ import annotations

import csv
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
DATA = ANALYSIS_DIR / "results" / "enzyme_functionality" / "ec_frequency_coverage.tsv"
OUT_DIR = ANALYSIS_DIR / "reports" / "enzymeFunctionality" / "figures"

INK, INK2, INK3 = "#0b0b0b", "#52514e", "#898781"
SURFACE, GRIDLINE = "#fcfcfb", "#e1e0d9"
STATUS_GOOD = "#0ca30c"
REMAINING_BLUE = "#9ec5f4"
LINE_BLUE = "#2a78d6"


def _tint(hex_colour: str, amount: float) -> tuple:
    r, g, b = (int(hex_colour[i:i + 2], 16) / 255 for i in (1, 3, 5))
    return tuple(c + (1.0 - c) * amount for c in (r, g, b))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ranks, cum, classified = [], [], []
    with open(DATA) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            ranks.append(int(row["ec_rank"]))
            cum.append(int(row["cumulative_unique_reactions_covered"]))
            classified.append(int(row["is_classified"]))

    n_ec_total = ranks[-1]
    total_reactions = cum[-1]
    n_classified = sum(classified)

    fig, ax = plt.subplots(figsize=(9.0, 5.8))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    ax.fill_between(ranks, cum, color=_tint(REMAINING_BLUE, 0.35), zorder=1)
    ax.plot(ranks, cum, color=LINE_BLUE, linewidth=2.0, zorder=3)

    if n_classified > 0:
        ax.fill_between(ranks[:n_classified], cum[:n_classified],
                         color=_tint(STATUS_GOOD, 0.78), zorder=2)
        ax.axvline(n_classified, color=STATUS_GOOD, linewidth=1.3, linestyle="--", zorder=2)
        ax.text(n_classified, total_reactions * 0.04,
                f" {n_classified:,} ECs classified so far", color=STATUS_GOOD,
                fontsize=9.0, ha="left", va="bottom", fontweight="semibold")

    # mark rank 500 (one literature-search batch) as a reference point
    if n_ec_total >= 500:
        y500 = cum[499]
        ax.plot([500], [y500], marker="o", color=INK, markersize=5.5, zorder=4)
        ax.annotate(f"top 500 ECs ({500 / n_ec_total:.0%} of {n_ec_total:,})\n"
                    f"→ {y500:,} reactions ({y500 / total_reactions:.0%} of {total_reactions:,})",
                    xy=(500, y500), xytext=(500 + n_ec_total * 0.04, y500 * 0.70),
                    color=INK2, fontsize=9.0, ha="left", va="top", linespacing=1.3,
                    arrowprops=dict(arrowstyle="-", color=INK3, linewidth=0.9))

    ax.grid(True, axis="y", color=GRIDLINE, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color(GRIDLINE)
    ax.spines["bottom"].set_color(GRIDLINE)
    ax.set_xlim(1, n_ec_total)
    ax.set_ylim(0, total_reactions * 1.05)
    ax.tick_params(colors=INK2, labelsize=9.6)
    ax.set_xlabel("unique complete-EC enzymes, ranked by citation count (most-cited first)",
                  color=INK2, fontsize=10.2)
    ax.set_ylabel("cumulative unique reactions covered", color=INK2, fontsize=10.2)

    fig.suptitle("Classifying the most-cited enzymes first front-loads reaction coverage",
                 x=0.09, y=0.975, ha="left", color=INK, fontsize=13.6, fontweight="semibold")
    ax.set_title(f"{n_ec_total:,} unique complete ECs cover {total_reactions:,} reactions "
                 f"(of 20,490 total; the rest have no complete EC)",
                 loc="left", color=INK2, fontsize=9.6, fontweight="normal", pad=12)

    fig.subplots_adjust(left=0.09, right=0.97, top=0.82, bottom=0.12)
    out_png = OUT_DIR / "enzyme_functionality_ec_coverage.png"
    fig.savefig(out_png, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
