#!/usr/bin/env python3
"""Explanatory 3x3 confusion matrix: recommended reversibility vs LLM ensemble.

Reads results/llm_vs_recommended/summary_counts.json (written by
build_llm_vs_recommended_mismatches.py) and draws every (recommended, LLM)
pair among the 20,490 reactions where both calls are directional. The diagonal
is agreement; the six off-diagonal cells are the six mismatch groups the
reaction-level detail is saved under in mismatches_by_group.json.

DESIGN NOTES
  * Diagonal cells use the fixed status-good green (dataviz palette) with an
    "AGREE" label -- a state, not a series, so it is exempt from the
    categorical-hue rule.
  * Off-diagonal cells use the single-hue sequential blue ramp (same palette),
    mapped on a log scale because the six counts span 62-6,098 (~100x): a
    linear map would make five of six cells look identical. Text colour flips
    to white above step 450 for contrast.
  * Every cell is direct-labelled with its own count, so there is no separate
    legend to cross-reference -- this is a 9-cell matrix, not a multi-series
    chart.
  * Numbers are read from the JSON, never retyped -- see
    llm_vs_recommended_confusion_values.tsv beside the image for the same
    figures in a diffable form.

Data: ModelSEED dev @ 078a395f (/scratch/ctaylor/tmp/devsnap_078a395f),
fetched 2026-10-07.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
SUMMARY = Path(os.environ.get("LLM_VS_RECOMMENDED_SUMMARY",
                              str(ANALYSIS_DIR / "results" / "llm_vs_recommended" / "summary_counts.json")))
OUT_DIR = Path(os.environ.get("LLM_VS_RECOMMENDED_OUT",
                              str(ANALYSIS_DIR / "reports" / "llmVsRecommended" / "figures")))

# --- shared chart tokens (same as plot_method_flowcharts.py / plot_graded_fba.py)
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#898781"
SURFACE, RULE = "#fcfcfb", "#c3c2b7"
STATUS_GOOD = "#0ca30c"
# sequential blue ramp, dataviz reference palette (steps 100..700)
SEQ_BLUE = {100: "#cde2fb", 150: "#b7d3f6", 200: "#9ec5f4", 250: "#86b6ef",
            300: "#6da7ec", 350: "#5598e7", 400: "#3987e5", 450: "#2a78d6",
            500: "#256abf", 550: "#1c5cab", 600: "#184f95", 650: "#104281", 700: "#0d366b"}
SEQ_STEPS = sorted(SEQ_BLUE)

DIRS = (">", "=", "<")
LABEL = {">": "forward", "=": "reversible", "<": "reverse"}
GROUP_KEY = {
    ("=", ">"): "reversible_vs_forward", ("=", "<"): "reversible_vs_reverse",
    (">", "="): "forward_vs_reversible", (">", "<"): "forward_vs_reverse",
    ("<", "="): "reverse_vs_reversible", ("<", ">"): "reverse_vs_forward",
}


def _tint(hex_colour: str, amount: float) -> tuple:
    r, g, b = (int(hex_colour[i:i + 2], 16) / 255 for i in (1, 3, 5))
    return tuple(c + (1.0 - c) * amount for c in (r, g, b))


def seq_colour(count: int, lo: int, hi: int) -> tuple:
    """Log-scaled position on the sequential ramp, steps 150..650."""
    t = (math.log1p(count) - math.log1p(lo)) / (math.log1p(hi) - math.log1p(lo))
    t = min(max(t, 0.0), 1.0)
    usable = [s for s in SEQ_STEPS if 150 <= s <= 650]
    idx = round(t * (len(usable) - 1))
    return SEQ_BLUE[usable[idx]], usable[idx]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary = json.load(open(SUMMARY))
    matrix = {}
    for key, count in summary["matrix"].items():
        rec, llm = key.split("|")
        matrix[(rec, llm)] = count

    off_diag_counts = [matrix[(r, l)] for r in DIRS for l in DIRS if r != l]
    lo, hi = min(off_diag_counts), max(off_diag_counts)

    fig, ax = plt.subplots(figsize=(10.5, 9.6))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    n_comp = summary["n_comparable"]
    ax.set_title("Recommended reversibility vs. LLM ensemble direction call",
                 color=INK, fontsize=15.5, loc="left", pad=20, fontweight="semibold")
    ax.text(0, 99.5,
            f"{n_comp:,} reactions where both calls are directional ('>' / '=' / '<')  ·  "
            f"ModelSEED dev @ 078a395f  ·  fetched 2026-10-07",
            color=INK2, fontsize=10.2, va="top", ha="left")

    # grid geometry: 3 columns (LLM), 3 rows (recommended), top-to-bottom = forward,reversible,reverse
    x0, y0, cell, gap = 22, 13, 20.0, 2.6
    col_x = {d: x0 + i * (cell + gap) for i, d in enumerate(DIRS)}
    row_y = {d: (y0 + 2 * (cell + gap)) - i * (cell + gap) for i, d in enumerate(DIRS)}

    # column headers (LLM ensemble)
    ax.text(x0 + 1.5 * cell + gap, row_y[">"] + cell + 10.5,
            "LLM ensemble direction call", color=INK, fontsize=12.0,
            ha="center", va="bottom", fontweight="semibold")
    for d in DIRS:
        ax.text(col_x[d] + cell / 2, row_y[">"] + cell + 2.0,
                f"{LABEL[d]}  ({d})", color=INK2, fontsize=10.6, ha="center", va="bottom")

    # row headers (recommended)
    ax.text(5.0, y0 + (cell + gap), "Recommended\nreversibility", color=INK, fontsize=12.0,
            ha="center", va="center", fontweight="semibold", rotation=90, linespacing=1.3)
    for d in DIRS:
        ax.text(x0 - 3.2, row_y[d] + cell / 2, f"{LABEL[d]} ({d})", color=INK2, fontsize=10.6,
                ha="right", va="center")

    for r in DIRS:
        for l in DIRS:
            x, y = col_x[l], row_y[r]
            count = matrix.get((r, l), 0)
            agree = (r == l)
            if agree:
                face = _tint(STATUS_GOOD, 0.86)
                edge = STATUS_GOOD
                text_colour = INK
                tag = "AGREE"
            else:
                face, step = seq_colour(count, lo, hi)
                edge = SEQ_BLUE[650]
                text_colour = "#ffffff" if step >= 450 else INK
                tag = GROUP_KEY[(r, l)]
            ax.add_patch(FancyBboxPatch(
                (x, y), cell, cell, boxstyle="round,pad=0,rounding_size=1.4",
                linewidth=1.4, edgecolor=edge, facecolor=face, zorder=3))
            ax.text(x + cell / 2, y + cell * 0.62, f"{count:,}", color=text_colour,
                    fontsize=17.5, ha="center", va="center", fontweight="bold", zorder=4)
            ax.text(x + cell / 2, y + cell * 0.30, tag, color=text_colour,
                    fontsize=7.6, ha="center", va="center", zorder=4, linespacing=1.3,
                    wrap=True)

    # diagonal-agreement swatch, left-aligned under the row headers
    ax.add_patch(FancyBboxPatch((0, 3.0), 3.2, 3.2, boxstyle="round,pad=0,rounding_size=0.3",
                                linewidth=1.0, edgecolor=STATUS_GOOD,
                                facecolor=_tint(STATUS_GOOD, 0.86), zorder=3))
    ax.text(4.2, 4.6, "diagonal =\nagreement", color=INK2, fontsize=8.4,
            ha="left", va="center", linespacing=1.25)

    # sequential legend (mismatch magnitude), spanning the grid width
    lx, ly, lw, lh = x0, 3.0, 3 * cell + 2 * gap, 3.2
    n_sw = 7
    usable = [s for s in SEQ_STEPS if 150 <= s <= 650]
    for i in range(n_sw):
        step = usable[round(i / (n_sw - 1) * (len(usable) - 1))]
        ax.add_patch(FancyBboxPatch(
            (lx + i * lw / n_sw, ly), lw / n_sw * 0.92, lh, boxstyle="round,pad=0,rounding_size=0.3",
            linewidth=0, facecolor=SEQ_BLUE[step], zorder=3))
    ax.text(lx, ly - 1.3, f"fewer mismatches ({lo:,})", color=INK3, fontsize=8.6, ha="left", va="top")
    ax.text(lx + lw, ly - 1.3, f"more mismatches ({hi:,})", color=INK3, fontsize=8.6, ha="right", va="top")

    excl = (f"excluded from this matrix: recommended '?' = {summary['n_recommended_unknown_excluded']:,}  ·  "
            f"no LLM entry = {summary['n_llm_missing_excluded']:,}  ·  LLM abstained '?' = "
            f"{summary['n_llm_abstained_excluded']:,}")
    ax.text(0, -3.6, excl, color=INK3, fontsize=8.6, ha="left", va="top")
    ax.text(0, -6.8,
            f"total mismatch: {summary['n_mismatch']:,} ({summary['n_mismatch'] / n_comp:.1%})   "
            f"total agreement: {summary['n_agree']:,} ({summary['n_agree'] / n_comp:.1%})",
            color=INK2, fontsize=9.2, ha="left", va="top", fontweight="semibold")

    fig.subplots_adjust(left=0.02, right=0.98, top=0.93, bottom=0.10)
    out_png = OUT_DIR / "llm_vs_recommended_confusion.png"
    fig.savefig(out_png, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out_png}")

    tsv_path = OUT_DIR / "llm_vs_recommended_confusion_values.tsv"
    with open(tsv_path, "w") as fh:
        fh.write("recommended\tllm\tcount\tgroup\n")
        for r in DIRS:
            for l in DIRS:
                group = "agree" if r == l else GROUP_KEY[(r, l)]
                fh.write(f"{r}\t{l}\t{matrix.get((r, l), 0)}\t{group}\n")
    print(f"wrote {tsv_path}")


if __name__ == "__main__":
    main()
