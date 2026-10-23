#!/usr/bin/env python3
"""Per-reaction: in how many models does each direction-change reaction exist,
and in how many does it actually carry flux, for each of the three panels.

One row per reaction, three facets (core_kegg2, ms2_gsm glucose-minimal,
ms2_gsm_auxo). In each facet the pale bar is the number of models CONTAINING
the reaction and the saturated bar drawn on top is the subset of those whose
FBA solution carries flux through it. Active is a subset of present by
construction, so the two are nested rather than grouped.

Reads results/direction_change_targets/presence_activity.json, written by
scripts/direction_change/presence_and_activity.py.

Non-obvious choices
-------------------
* Only reactions present in at least one model are drawn (36 of 76); the other
  40 would be three empty rows each and are counted in the subtitle instead.
* Rows are sorted by total active models across panels, descending, so the
  reactions that actually matter are at the top.
* Counts, not fractions: the three panels are 5,683 / 5,420 / 5,420 models, so
  raw counts stay comparable and avoid implying a precision fractions do not
  have on reactions present in only a handful of models.
* "Active" is the **pFBA** flux vector (ACTIVE_TOL 1e-9), per
  build_flux_distributions.py; plain FBA picks an arbitrary vertex.
* A shared x-axis across facets, so bar lengths are comparable between panels.

Data: core_models_analysis @ cc44689. Sweep run 2026-10-23 over 16,523 models,
0 solver errors; a second run reproduced every count.
"""
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
from cma import paths  # noqa: E402

DATA = paths.results("direction_change_targets") / "presence_activity.json"
OUT_DIR = Path(os.environ.get(
    "DIRECTION_CHANGE_FIG_DIR",
    paths.results("direction_change_targets") / "figures"))

# dataviz reference palette, categorical slots 1-3. Validated light mode under
# --pairs all: CVD dE 9.2, normal-vision dE 24.0, all checks pass. Aqua carries
# a contrast WARN (2.74:1), relieved by the axis labels and the committed TSV.
PANELS = [
    ("core_kegg2", "core  ·  5,683 models", "#2a78d6"),
    ("ms2_gsm", "genome-scale  ·  glucose minimal  ·  5,420", "#eb6834"),
    ("ms2_gsm_auxo", "genome-scale  ·  auxotrophy  ·  5,420", "#1baf7a"),
]
PRESENT_ALPHA = 0.26
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#898781"
SURFACE, GRIDLINE = "#fcfcfb", "#e1e0d9"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    d = json.load(open(DATA))
    targets: dict[str, str] = d["targets"]
    panels: dict[str, dict] = d["panels"]

    present_any = [r for r in targets
                   if any(panels[k]["per_reaction"][r]["present"]
                          for k, _, _ in PANELS)]
    order = sorted(
        present_any,
        key=lambda r: (sum(panels[k]["per_reaction"][r]["active_pfba"]
                           for k, _, _ in PANELS),
                       sum(panels[k]["per_reaction"][r]["present"]
                           for k, _, _ in PANELS)))
    n = len(order)
    ypos = list(range(n))

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 0.30 * n + 2.6),
                             sharey=True, sharex=True)
    fig.patch.set_facecolor(SURFACE)

    xmax = max(panels[k]["per_reaction"][r]["present"]
               for k, _, _ in PANELS for r in order)

    for ax, (key, title, color) in zip(axes, PANELS):
        ax.set_facecolor(SURFACE)
        pres = [panels[key]["per_reaction"][r]["present"] for r in order]
        act = [panels[key]["per_reaction"][r]["active_pfba"] for r in order]
        ax.barh(ypos, pres, color=color, alpha=PRESENT_ALPHA, height=0.74,
                zorder=3, linewidth=0)
        # 2px surface ring where the saturated bar overlaps the pale one
        ax.barh(ypos, act, color=color, height=0.52, zorder=4,
                edgecolor=SURFACE, linewidth=1.4)
        ax.set_title(title, fontsize=10.2, color=INK, fontweight="bold",
                     loc="left", pad=9)
        ax.grid(axis="x", color=GRIDLINE, linewidth=0.7, zorder=0)
        ax.set_axisbelow(True)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.spines["bottom"].set_color(GRIDLINE)
        ax.tick_params(axis="both", colors=INK2, length=0, labelsize=9)
        ax.set_xlabel("models", fontsize=10, color=INK2)

    # The pale/solid distinction is the only 2-series encoding in the figure,
    # so it gets a real legend as well as the subtitle. Demonstrated in slot-1
    # blue; the nesting is identical in every facet. Placed low-left in the
    # core facet, where the bottom rows are empty.
    c = PANELS[0][2]
    axes[0].legend(handles=[
        plt.Rectangle((0, 0), 1, 1, color=c, alpha=PRESENT_ALPHA,
                      label="present in model"),
        plt.Rectangle((0, 0), 1, 1, color=c, label="carries flux (pFBA)"),
    ], loc="lower right", frameon=False, fontsize=9, labelcolor=INK2,
        handlelength=1.5, bbox_to_anchor=(1.0, 0.012))

    axes[0].set_ylim(-0.8, n - 0.2)
    axes[0].set_xlim(0, xmax * 1.04)
    axes[0].set_yticks(ypos)
    axes[0].set_yticklabels(
        [f"{r}  {targets[r] if targets[r] != 'no change' else '(uncert.)'}"
         for r in order], fontsize=8.4, color=INK, fontfamily="monospace")

    fig.suptitle("Direction-change reactions: how many models contain each, "
                 "and how many run it",
                 fontsize=13.5, color=INK, fontweight="bold", x=0.012,
                 ha="left", y=0.995)
    fig.text(0.012, 0.962,
             f"pale = models containing the reaction   ·   solid = models whose "
             f"pFBA solution carries flux through it   ·   "
             f"{n} of 76 reactions appear in any model; the other "
             f"{76 - n} are shown in no row",
             fontsize=9.3, color=INK3, ha="left")

    fig.tight_layout(rect=(0, 0, 1, 0.952))
    # PNG only: no figure PDFs are tracked in this repo.
    fig.savefig(OUT_DIR / "direction_change_by_reaction.png", dpi=200, facecolor=SURFACE)
    plt.close(fig)

    stats = OUT_DIR / "direction_change_by_reaction.tsv"
    with open(stats, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        head = ["reaction_id", "change"]
        for key, _, _ in PANELS:
            head += [f"present_{key}", f"active_{key}"]
        w.writerow(head)
        for r in reversed(order):
            row = [r, targets[r]]
            for key, _, _ in PANELS:
                v = panels[key]["per_reaction"][r]
                row += [v["present"], v["active_pfba"]]
            w.writerow(row)
    print(f"wrote {OUT_DIR}/direction_change_by_reaction.png (+tsv), {n} rows")


if __name__ == "__main__":
    main()
