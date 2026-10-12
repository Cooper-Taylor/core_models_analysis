#!/usr/bin/env python3
"""Of the 389 uniquely-fixed charge-imbalance-only reactions, which redox-
cofactor-pair family needed its coefficient corrected.

Reads results/ci_only_fixes/pair_family_counts.tsv (from
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
DATA = ANALYSIS_DIR / "results" / "ci_only_fixes" / "pair_family_counts.tsv"
OUT_DIR = ANALYSIS_DIR / "reports" / "chargeImbalanceFixes" / "figures"

INK, INK2, INK3 = "#0b0b0b", "#52514e", "#898781"
SURFACE, GRIDLINE = "#fcfcfb", "#e1e0d9"
BAR_COLOR = "#1baf7a"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    with open(DATA) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append((r["family"], int(r["n"])))
    rows.sort(key=lambda r: r[1])
    total = sum(n for _, n in rows)

    fig, ax = plt.subplots(figsize=(9.0, 3.6))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    y = range(len(rows))
    ax.barh(list(y), [n for _, n in rows], color=BAR_COLOR, zorder=3)
    ax.set_yticks(list(y))
    ax.set_yticklabels([fam for fam, _ in rows], color=INK2, fontsize=10.4)
    ax.tick_params(axis="x", colors=INK2, labelsize=9.4)
    ax.grid(True, axis="x", color=GRIDLINE, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color(GRIDLINE)
    ax.set_xlabel("reactions fixed", color=INK2, fontsize=10.0)

    xmax = max(n for _, n in rows)
    ax.set_xlim(0, xmax * 1.18)
    for yi, (fam, n) in zip(y, rows):
        ax.text(n + xmax * 0.02, yi, f"{n:,} ({n/total:.0%})", va="center", ha="left",
                color=INK3, fontsize=9.2)

    fig.suptitle("The 389 verified fixes, by which cofactor pair was mis-scaled",
                 x=0.07, y=0.96, ha="left", color=INK, fontsize=13.0, fontweight="semibold")
    ax.set_title("coefficient of an oxidized/reduced couple corrected to its unique charge-balancing value",
                 loc="left", color=INK2, fontsize=9.2, fontweight="normal", pad=10)

    fig.subplots_adjust(left=0.22, right=0.95, top=0.80, bottom=0.13)
    out_png = OUT_DIR / "ci_only_fix_families.png"
    fig.savefig(out_png, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
