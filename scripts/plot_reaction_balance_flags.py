#!/usr/bin/env python3
"""Flag-combination breakdown across all 48,384 non-obsolete ModelSEED
reactions (the `status` field, see analyze_reaction_balance_flags.py).

Reads results/reaction_balance_flags/flag_counts.tsv. One bar per observed
status-token combination. `OK` gets the status-good green (it's a state, not
a defect); the two curator-override combinations (`CI+CK`, `CK+OK`) get a
distinguishing outline since they're the only reactions anyone has actually
reviewed; everything else is on the neutral->sequential-blue scale used
elsewhere in this project for "not yet resolved" magnitude.

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
DATA = ANALYSIS_DIR / "results" / "reaction_balance_flags" / "flag_counts.tsv"
OUT_DIR = ANALYSIS_DIR / "reports" / "reactionBalanceFlags" / "figures"

INK, INK2, INK3 = "#0b0b0b", "#52514e", "#898781"
SURFACE, GRIDLINE = "#fcfcfb", "#e1e0d9"
STATUS_GOOD = "#0ca30c"
SEQ_BLUE = {"CPDFORMERROR": "#184f95", "CI+MI": "#2a78d6", "CI": "#5598e7",
            "MI": "#86b6ef", "CI+CK": "#cde2fb", "CK+OK": "#cde2fb", "EMPTY": "#e1e0d9"}

PRETTY = {"OK": "OK (balanced)", "CPDFORMERROR": "CPDFORMERROR\n(uncheckable)",
          "CI+MI": "CI + MI\n(both)", "CI": "CI only\n(charge)", "MI": "MI only\n(mass)",
          "CI+CK": "CI + CK\n(charge, overridden)", "CK+OK": "OK + CK\n(confirmed)", "EMPTY": "EMPTY"}


def _tint(hex_colour: str, amount: float) -> tuple:
    r, g, b = (int(hex_colour[i:i + 2], 16) / 255 for i in (1, 3, 5))
    return tuple(c + (1.0 - c) * amount for c in (r, g, b))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    with open(DATA) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append({"combo": r["combination"], "n": int(r["n"]), "share": float(r["share"])})
    total = sum(r["n"] for r in rows)
    rows.sort(key=lambda r: r["n"])  # ascending -> largest at top of barh

    fig, ax = plt.subplots(figsize=(10.0, 6.0))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    y = range(len(rows))
    colors = [STATUS_GOOD if r["combo"] == "OK" else SEQ_BLUE.get(r["combo"], "#9ec5f4") for r in rows]
    edge_colors = [INK if "CK" in r["combo"] else "none" for r in rows]
    edge_widths = [1.8 if "CK" in r["combo"] else 0 for r in rows]
    ax.barh(list(y), [r["n"] for r in rows], color=colors, edgecolor=edge_colors,
            linewidth=edge_widths, zorder=3)

    ax.set_yticks(list(y))
    ax.set_yticklabels([PRETTY.get(r["combo"], r["combo"]) for r in rows], color=INK2, fontsize=10.0)
    ax.tick_params(axis="x", colors=INK2, labelsize=9.6)
    ax.set_xscale("symlog", linthresh=10)
    ax.grid(True, axis="x", color=GRIDLINE, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color(GRIDLINE)
    ax.set_xlabel("reactions (log scale) -- of 48,384 non-obsolete", color=INK2, fontsize=10.2)

    for yi, r in zip(y, rows):
        ax.text(r["n"] * 1.15 if r["n"] > 0 else 0.5, yi, f"{r['n']:,}  ({r['share']:.1%})",
                va="center", ha="left", color=INK, fontsize=9.8,
                fontweight="bold" if r["combo"] in ("OK", "CPDFORMERROR") else "normal")

    ax.set_xlim(0.5, total * 3.0)

    fig.suptitle("Every reaction status flag combination, database-wide",
                 x=0.09, y=0.975, ha="left", color=INK, fontsize=13.6, fontweight="semibold")
    ax.set_title("48,384 non-obsolete reactions  ·  bold outline = curator-reviewed (CK)",
                 loc="left", color=INK2, fontsize=9.6, fontweight="normal", pad=12)

    fig.subplots_adjust(left=0.24, right=0.90, top=0.84, bottom=0.10)
    out_png = OUT_DIR / "reaction_balance_flag_counts.png"
    fig.savefig(out_png, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
