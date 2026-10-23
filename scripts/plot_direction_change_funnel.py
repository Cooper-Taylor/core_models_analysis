#!/usr/bin/env python3
"""How many of the 76 direction-change reactions survive each stage: all 76 ->
present in at least one model -> carrying flux in at least one FBA solution.

Each stage is split into the two categories of the source list: 56 reactions
with an explicit direction change (``= -> >`` / ``= -> <``) and 20 labelled
"no change" that nonetheless move because of thermodynamic uncertainty.

Reads results/direction_change_targets/presence_activity.json, written by
scripts/direction_change/presence_and_activity.py.

Non-obvious choices
-------------------
* "present"/"active" are UNION over the three panels (core_kegg2, ms2_gsm,
  ms2_gsm_auxo), so a reaction counts once if any panel has it. The per-panel
  split is the other figure, ``direction_change_by_reaction``.
* "active" is measured on the **pFBA** flux vector (ACTIVE_TOL 1e-9), matching
  build_flux_distributions.py -- plain FBA picks an arbitrary vertex among
  alternate optima, so reaction-level activity in it is not reproducible.
* Stacking is used rather than grouping because the two categories partition
  each stage exactly.

Data: core_models_analysis @ cc44689; ModelSEED dev snapshot as pinned by
cma.paths.msdb("live"). Sweep run 2026-10-23 over 16,523 models, 0 errors.
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
    paths.analysis_dir() / "reports" / "directionChangeTargets" / "figures"))

# dataviz reference palette, categorical slots 1 and 2. Validated light mode:
# CVD dE 24.7, normal-vision dE 33.6, contrast >=3:1 -- all checks pass.
EXPLICIT_COLOR = "#2a78d6"   # slot 1, blue
UNCERTAIN_COLOR = "#eb6834"  # slot 2, orange
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#898781"
SURFACE, GRIDLINE = "#fcfcfb", "#e1e0d9"

EXPLICIT_LABEL = "explicit direction change  (= → > or = → <)"
UNCERTAIN_LABEL = "labelled “no change”  (moves on uncertainty)"

STAGES = [
    ("all targets", "the reaction list under test"),
    ("present in ≥1 model", "exists in core or genome-scale models"),
    ("active in ≥1 FBA solution", "carries flux under pFBA"),
]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    d = json.load(open(DATA))
    targets: dict[str, str] = d["targets"]
    panels: dict[str, dict] = d["panels"]

    explicit = {r for r, c in targets.items() if c != "no change"}
    uncertain = {r for r, c in targets.items() if c == "no change"}

    present, active = set(), set()
    for p in panels.values():
        for r, v in p["per_reaction"].items():
            if v["present"]:
                present.add(r)
            if v["active_pfba"]:
                active.add(r)

    rows = [
        (STAGES[0][0], STAGES[0][1], explicit, uncertain),
        (STAGES[1][0], STAGES[1][1], explicit & present, uncertain & present),
        (STAGES[2][0], STAGES[2][1], explicit & active, uncertain & active),
    ]

    fig, ax = plt.subplots(figsize=(10.5, 4.7))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    ypos = list(range(len(rows)))[::-1]
    for y, (_, _, ex, un) in zip(ypos, rows):
        ax.barh([y], [len(ex)], color=EXPLICIT_COLOR, zorder=3,
                edgecolor=SURFACE, linewidth=2, height=0.56)
        ax.barh([y], [len(un)], left=len(ex), color=UNCERTAIN_COLOR, zorder=3,
                edgecolor=SURFACE, linewidth=2, height=0.56)
        # direct labels inside each segment (also the relief the contrast
        # check asks for), plus the stage total at the end of the bar
        if len(ex):
            ax.text(len(ex) / 2, y, str(len(ex)), va="center", ha="center",
                    color="white", fontsize=11, fontweight="bold", zorder=4)
        if len(un):
            ax.text(len(ex) + len(un) / 2, y, str(len(un)), va="center",
                    ha="center", color="white", fontsize=11, fontweight="bold",
                    zorder=4)
        ax.text(len(ex) + len(un) + 1.2, y, f"{len(ex) + len(un)}",
                va="center", ha="left", color=INK, fontsize=12,
                fontweight="bold", zorder=4)

    ax.set_yticks(ypos)
    ax.set_yticklabels([f"{n}\n" for n, _, _, _ in rows], fontsize=11.5,
                       color=INK)
    for y, (_, sub, _, _) in zip(ypos, rows):
        ax.text(-1.4, y - 0.22, sub, va="center", ha="right", fontsize=9,
                color=INK3, style="italic")

    ax.set_xlim(0, 86)
    ax.set_xlabel("reactions", fontsize=10.5, color=INK2)
    ax.grid(axis="x", color=GRIDLINE, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(GRIDLINE)
    ax.tick_params(axis="both", colors=INK2, length=0)

    n_never = 76 - len(present)
    fig.tight_layout(rect=(0, 0.055, 1, 0.855))

    # Title and legend are figure-level, not axes-level: the axes starts a
    # quarter of the way in (wide stage labels), which clipped a two-column
    # legend anchored to it.
    fig.suptitle("Of 76 reactions whose direction changes, how many reach an "
                 "FBA solution?",
                 fontsize=13, color=INK, fontweight="bold", x=0.012, ha="left",
                 y=0.975)
    fig.legend(handles=[
        plt.Rectangle((0, 0), 1, 1, color=EXPLICIT_COLOR, label=EXPLICIT_LABEL),
        plt.Rectangle((0, 0), 1, 1, color=UNCERTAIN_COLOR, label=UNCERTAIN_LABEL),
    ], loc="upper left", bbox_to_anchor=(0.012, 0.925), ncol=2, frameon=False,
        fontsize=9.5, labelcolor=INK2, handlelength=1.5, columnspacing=2.4)
    fig.text(0.012, 0.022,
             f"{n_never} of the 76 are in no model at all; "
             f"{len(present - active)} more are present but never carry flux.",
             fontsize=9.5, color=INK3, ha="left")
    # PNG only: no figure PDFs are tracked in this repo.
    fig.savefig(OUT_DIR / "direction_change_funnel.png", dpi=200, facecolor=SURFACE)
    plt.close(fig)

    stats = OUT_DIR / "direction_change_funnel.tsv"
    with open(stats, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["stage", "explicit_change", "labelled_no_change", "total"])
        for name, _, ex, un in rows:
            w.writerow([name, len(ex), len(un), len(ex) + len(un)])
    print(f"wrote {OUT_DIR}/direction_change_funnel.png (+tsv)")


if __name__ == "__main__":
    main()
