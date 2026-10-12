#!/usr/bin/env python3
"""Of the 2,824 charge-imbalance-only reactions, how many have a uniquely
verified, zero-blast-radius fix (rescaling an existing same-formula redox-
cofactor pair) vs. need individual literature review.

Reads results/ci_only_fixes/outcome_counts.tsv (from
scripts/fix_ci_only_reactions.py).

Data: ModelSEED dev @ 078a395f (/scratch/ctaylor/tmp/devsnap_078a395f),
confirmed == origin/dev HEAD, 2026-10-12.
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
DATA = ANALYSIS_DIR / "results" / "ci_only_fixes" / "outcome_counts.tsv"
OUT_DIR = ANALYSIS_DIR / "reports" / "chargeImbalanceFixes" / "figures"

INK, INK2, INK3 = "#0b0b0b", "#52514e", "#898781"
SURFACE, GRIDLINE = "#fcfcfb", "#e1e0d9"
SOLVED_COLOR = "#1baf7a"
AMBIGUOUS_COLOR = "#e0a526"
UNRESOLVED_COLOR = "#c3c2b7"
COLOR = {"solved": SOLVED_COLOR, "ambiguous": AMBIGUOUS_COLOR, "unresolved": UNRESOLVED_COLOR}
LABEL = {
    "solved": "uniquely fixed & verified",
    "ambiguous": "ambiguous (2+ valid pairs)",
    "unresolved": "needs individual review",
}


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    counts = {}
    with open(DATA) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            counts[r["outcome"]] = int(r["n"])
    total = sum(counts.values())
    order = ["solved", "ambiguous", "unresolved"]

    fig, ax = plt.subplots(figsize=(9.5, 2.6))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    left = 0.0
    for key in order:
        n = counts[key]
        ax.barh([0], [n], left=left, color=COLOR[key], zorder=3,
                edgecolor=SURFACE, linewidth=2)
        if n / total > 0.04:
            ax.text(left + n / 2, 0, f"{n:,}\n{n/total:.1%}", va="center", ha="center",
                    color="white" if key != "ambiguous" else INK, fontsize=10.2, fontweight="bold")
        left += n

    ax.set_xlim(0, total)
    ax.set_ylim(-0.6, 0.6)
    ax.set_yticks([])
    ax.set_xticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    handles = [plt.Rectangle((0, 0), 1, 1, color=COLOR[k]) for k in order]
    ax.legend(handles, [f"{LABEL[k]} ({counts[k]:,})" for k in order],
              loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=3,
              frameon=False, fontsize=9.4, labelcolor=INK2, handlelength=1.1)

    fig.suptitle("2,824 charge-imbalance-only reactions: how many does arithmetic alone fix?",
                 x=0.055, y=0.97, ha="left", color=INK, fontsize=13.0, fontweight="semibold")
    ax.set_title(
        "fixed = rescaling an existing redox-cofactor pair (NAD/NADH, ferredoxin, ...) already "
        "in the equation,\nverified to zero charge without disturbing the already-balanced mass",
        loc="left", color=INK2, fontsize=9.0, fontweight="normal", pad=10)

    fig.subplots_adjust(left=0.03, right=0.98, top=0.72, bottom=0.30)
    out_png = OUT_DIR / "ci_only_fix_outcomes.png"
    fig.savefig(out_png, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
