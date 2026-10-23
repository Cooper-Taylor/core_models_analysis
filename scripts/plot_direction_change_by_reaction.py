#!/usr/bin/env python3
"""Per-reaction: in how many models does each direction-change reaction exist,
in how many of those does the model grow, and in how many does it carry flux.

Three nested counts per reaction per panel, drawn as nested bars of one hue:

    present  >=  growing  >=  active

* **present** - models containing the reaction at all (palest).
* **growing** - of those, the ones with non-zero biomass (mid). A non-grower
  has no flux solution, so this is the denominator ``active`` should be read
  against, not ``present``.
* **active**  - of those, the ones whose pFBA solution carries flux (solid).

Reads results/direction_change_targets/presence_activity.json, written by
scripts/direction_change/presence_and_activity.py.

Variants (see scripts/figures.tsv)
----------------------------------
``--orient vertical``   reactions down the y-axis, bars running across. Tall;
                        fits all 36 reactions comfortably. The default.
``--orient horizontal`` reactions along the x-axis, bars running up. Wide;
                        facets become stacked rows sharing the x-axis.
``--panel <set>``       restrict to one model set and drop the other facets,
                        showing only the reactions present in that set.

Non-obvious choices
-------------------
* Only reactions present in at least one of the panels being drawn are shown;
  the rest would be empty rows and are counted in the subtitle instead.
* Sorted by total active models, descending, so what matters is read first.
* Counts, not fractions: the panels are 5,683 / 5,420 / 5,420 models, so raw
  counts stay comparable without implying precision fractions would not have
  on reactions present in a handful of models.
* "Active" is the **pFBA** flux vector (ACTIVE_TOL 1e-9), per
  build_flux_distributions.py; plain FBA picks an arbitrary vertex.
* A shared count axis across facets, so bar lengths compare between panels.

Data: core_models_analysis @ 0a4f619. Sweep run 2026-10-23 over 16,523 models,
0 solver errors; repeated runs reproduced every count.
"""
from __future__ import annotations

import argparse
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
# a contrast WARN (2.74:1), relieved by the axis labels and the committed TSVs.
PANELS = [
    ("core_kegg2", "core", "core  ·  5,683 models", "#2a78d6"),
    ("ms2_gsm", "gsm", "genome-scale  ·  glucose minimal  ·  5,420", "#eb6834"),
    ("ms2_gsm_auxo", "gsm_auxo", "genome-scale  ·  auxotrophy  ·  5,420", "#1baf7a"),
]
#: Nested tiers, palest to solid. One hue per panel, three lightness steps.
TIERS = [
    ("present", 0.22, "present in model"),
    ("present_growing", 0.52, "…and the model grows"),
    ("active_pfba", 1.0, "…and it carries flux (pFBA)"),
]
BAR_SPAN = [0.80, 0.80, 0.52]   # thickness per tier; active drawn narrower
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#898781"
SURFACE, GRIDLINE = "#fcfcfb", "#e1e0d9"


def load(panel_keys: list[str]):
    d = json.load(open(DATA))
    targets: dict[str, str] = d["targets"]
    panels: dict[str, dict] = d["panels"]

    def tier(key, r, name):
        return panels[key]["per_reaction"][r].get(name, 0)

    present_any = [r for r in targets
                   if any(tier(k, r, "present") for k in panel_keys)]
    order = sorted(present_any,
                   key=lambda r: (sum(tier(k, r, "active_pfba") for k in panel_keys),
                                  sum(tier(k, r, "present") for k in panel_keys)))
    return targets, panels, order, tier


def row_label(r, targets):
    c = targets[r] if targets[r] != "no change" else "(uncert.)"
    return f"{r}  {c}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--orient", choices=("vertical", "horizontal"),
                    default="vertical")
    ap.add_argument("--panel", default="all",
                    choices=["all"] + [k for k, _, _, _ in PANELS])
    ap.add_argument("--out-name", default=None)
    a = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    chosen = [p for p in PANELS if a.panel in ("all", p[0])]
    keys = [k for k, _, _, _ in chosen]
    targets, panels, order, tier = load(keys)
    n = len(order)
    short_of = {k: s for k, s, _, _ in PANELS}
    name = a.out_name or ("direction_change_by_reaction"
                          + ("" if a.orient == "vertical" else "_horizontal")
                          + ("" if a.panel == "all" else "_" + short_of[a.panel]))

    cmax = max(tier(k, r, "present") for k in keys for r in order)
    pos = list(range(n))

    if a.orient == "vertical":
        fig, axes = plt.subplots(1, len(chosen),
                                 figsize=(4.6 * len(chosen) + 0.4, 0.30 * n + 2.6),
                                 sharey=True, sharex=True, squeeze=False)
        axes = axes[0]
    else:
        fig, axes = plt.subplots(len(chosen), 1,
                                 figsize=(max(0.40 * n + 3.4, 9.6),
                                          3.3 * len(chosen) + 1.6),
                                 sharex=True, sharey=True, squeeze=False)
        axes = axes[:, 0]
    fig.patch.set_facecolor(SURFACE)

    for ax, (key, _, title, color) in zip(axes, chosen):
        ax.set_facecolor(SURFACE)
        for (field, alpha, _), span in zip(TIERS, BAR_SPAN):
            vals = [tier(key, r, field) for r in order]
            kw = dict(color=color, alpha=alpha, zorder=3 + alpha,
                      edgecolor=SURFACE, linewidth=1.4 if alpha == 1.0 else 0)
            if a.orient == "vertical":
                ax.barh(pos, vals, height=span, **kw)
            else:
                ax.bar(pos, vals, width=span, **kw)
        ax.set_title(title, fontsize=10.2, color=INK, fontweight="bold",
                     loc="left", pad=9)
        ax.grid(axis="x" if a.orient == "vertical" else "y",
                color=GRIDLINE, linewidth=0.7, zorder=0)
        ax.set_axisbelow(True)
        hide = ("top", "right", "left") if a.orient == "vertical" else ("top", "right")
        for s in hide:
            ax.spines[s].set_visible(False)
        ax.spines["bottom"].set_color(GRIDLINE)
        if a.orient == "horizontal":
            ax.spines["left"].set_color(GRIDLINE)
        ax.tick_params(axis="both", colors=INK2, length=0, labelsize=9)

    labels = [row_label(r, targets) for r in order]
    if a.orient == "vertical":
        axes[0].set_ylim(-0.8, n - 0.2)
        axes[0].set_xlim(0, cmax * 1.04)
        axes[0].set_yticks(pos)
        axes[0].set_yticklabels(labels, fontsize=8.4, color=INK,
                                fontfamily="monospace")
        for ax in axes:
            ax.set_xlabel("models", fontsize=10, color=INK2)
    else:
        axes[-1].set_xlim(-0.8, n - 0.2)
        axes[0].set_ylim(0, cmax * 1.04)
        axes[-1].set_xticks(pos)
        axes[-1].set_xticklabels(labels, fontsize=8.0, color=INK,
                                 fontfamily="monospace", rotation=90)
        for ax in axes:
            ax.set_ylabel("models", fontsize=10, color=INK2)

    c0 = chosen[0][3]
    handles = [plt.Rectangle((0, 0), 1, 1, color=c0, alpha=al, label=lab)
               for _, al, lab in TIERS]

    title_of = {k: t for k, _, t, _ in PANELS}
    scope = ("across all three model sets" if a.panel == "all"
             else f"in the {title_of[a.panel].split('  ·')[0]} set")
    # Header offsets in INCHES converted to figure fraction: the two layouts
    # differ in height by 2x, so a fixed fraction collided the subtitle with
    # the title on the shorter one.
    h = fig.get_figheight()
    # In the horizontal layout the bars are sorted ascending and fill the
    # plot, so there is no free corner for a legend: it goes in the header
    # band instead. The vertical layout has empty space at lower right.
    if a.orient == "horizontal":
        fig.legend(handles=handles, loc="upper left",
                   bbox_to_anchor=(0.008, 1 - 0.62 / h), ncol=3, frameon=False,
                   fontsize=9, labelcolor=INK2, handlelength=1.5,
                   columnspacing=2.2)
    else:
        axes[0].legend(handles=handles, loc="lower right", frameon=False,
                       fontsize=9, labelcolor=INK2, handlelength=1.5)
    fig.suptitle("Direction-change reactions: models containing each, growing, "
                 "and running it",
                 fontsize=13.5, color=INK, fontweight="bold", x=0.008,
                 ha="left", va="top", y=1 - 0.10 / h)
    fig.text(0.008, 1 - 0.40 / h,
             f"nested bars: present ≥ growing ≥ active   ·   {n} of 76 "
             f"reactions appear {scope}; the other {76 - n} are not drawn",
             fontsize=9.3, color=INK3, ha="left", va="top")

    fig.tight_layout(rect=(0, 0, 1,
                           1 - (1.02 if a.orient == "horizontal" else 0.72) / h))
    fig.savefig(OUT_DIR / f"{name}.png", dpi=200, facecolor=SURFACE)
    plt.close(fig)

    with open(OUT_DIR / f"{name}.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        head = ["reaction_id", "change"]
        for key, short, _, _ in chosen:
            head += [f"present_{short}", f"growing_{short}", f"active_{short}"]
        w.writerow(head)
        for r in reversed(order):
            row = [r, targets[r]]
            for key, _, _, _ in chosen:
                row += [tier(key, r, "present"), tier(key, r, "present_growing"),
                        tier(key, r, "active_pfba")]
            w.writerow(row)
    print(f"wrote {OUT_DIR}/{name}.png (+tsv), {n} reactions, "
          f"{len(chosen)} panel(s), {a.orient}")


if __name__ == "__main__":
    main()
