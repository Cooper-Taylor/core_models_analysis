#!/usr/bin/env python3
"""How much of each imbalance group reduces to one identified, tractable
pattern vs. genuinely unexplained residual (the "ease of fixing" chart).

Reads results/reaction_balance_flags/group_tractability.tsv. Each bar is
split into the identified-pattern share (status-good green -- NOT "fixed",
just "has a concrete, named explanation that narrows the fix") and the
residual share (neutral grey -- needs individual review). CPDFORMERROR's
"identified" share is explicitly NOT the same kind of win as the other three
-- see the caption and the report text -- it's a concentration finding, not
a fixability one, since these are generic protein cofactors by design.

Data: ModelSEED dev @ 078a395f (/scratch/ctaylor/tmp/devsnap_078a395f),
confirmed == origin/dev HEAD, 2026-10-09.
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
DATA = ANALYSIS_DIR / "results" / "reaction_balance_flags" / "group_tractability.tsv"
OUT_DIR = ANALYSIS_DIR / "reports" / "reactionBalanceFlags" / "figures"

INK, INK2, INK3 = "#0b0b0b", "#52514e", "#898781"
SURFACE, GRIDLINE = "#fcfcfb", "#e1e0d9"
IDENTIFIED_COLOR = "#1baf7a"
RESIDUAL_COLOR = "#c3c2b7"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    with open(DATA) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append({"group": r["group"], "n": int(r["n_total"]),
                         "identified": int(r["n_identified_pattern"]),
                         "residual": int(r["n_residual"]), "label": r["pattern_label"]})
    rows.sort(key=lambda r: r["n"])

    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    y = range(len(rows))
    ax.barh(list(y), [r["identified"] for r in rows], color=IDENTIFIED_COLOR, zorder=3,
            label="identified pattern")
    ax.barh(list(y), [r["residual"] for r in rows], left=[r["identified"] for r in rows],
            color=RESIDUAL_COLOR, zorder=3, label="unexplained residual")

    ax.set_yticks(list(y))
    ax.set_yticklabels([r["group"] for r in rows], color=INK2, fontsize=10.6)
    ax.tick_params(axis="x", colors=INK2, labelsize=9.6)
    ax.grid(True, axis="x", color=GRIDLINE, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color(GRIDLINE)
    ax.set_xlabel("reactions", color=INK2, fontsize=10.2)

    xmax = max(r["n"] for r in rows)
    ax.set_xlim(0, xmax * 1.05)
    for yi, r in zip(y, rows):
        pct = r["identified"] / r["n"]
        ax.text(r["identified"] / 2, yi, f"{r['identified']:,}\n({pct:.0%})", va="center", ha="center",
                color="white", fontsize=9.0, fontweight="bold")
        ax.text(r["n"] + xmax * 0.012, yi, r["label"], va="center", ha="left", color=INK3, fontsize=8.6)

    fig.suptitle("How much of each imbalance group has a named, concrete cause",
                 x=0.07, y=0.975, ha="left", color=INK, fontsize=13.4, fontweight="semibold")
    ax.set_title("green = reactions matching the pattern described below; grey = still needs individual review",
                 loc="left", color=INK2, fontsize=9.4, fontweight="normal", pad=12)
    ax.legend(loc="lower right", frameon=False, fontsize=9.2, labelcolor=INK2)

    fig.subplots_adjust(left=0.14, right=0.60, top=0.82, bottom=0.12)
    out_png = OUT_DIR / "reaction_balance_group_tractability.png"
    fig.savefig(out_png, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
