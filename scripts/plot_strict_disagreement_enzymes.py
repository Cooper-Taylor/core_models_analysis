#!/usr/bin/env python3
"""Two figures for the strict-disagreement enzyme-functionality report.

Reads results/strict_disagreement_enzymes/summary.json (written by
build_strict_disagreement_outputs.py) and draws:

  strict_disagreement_functionality.png
      The question as asked: of the 1,370 reactions where ModelSEED's
      recommended direction and the LLM council's call point in strictly
      opposite directions, how many are run by a bifunctional or
      polyfunctional enzyme? Answer: very few -- which is itself the finding,
      so the chart has to make the smallness legible rather than hide it.

  strict_disagreement_direction.png
      Why those reactions actually disagree. Left panel: which side the
      enzymology supports, split by mismatch group, because the two groups
      behave completely differently and a pooled bar would mislead. Right
      panel: how often the ModelSEED equation itself is defective, with the
      database's own mass/charge-imbalance flag as an independent check.

DESIGN NOTES
  * Categorical hues are this repo's validated set (#2a78d6 / #eb6834 /
    #1baf7a); "unknown" and "not researched" are gray because absence of an
    answer is a state, not a series.
  * The functionality chart is a horizontal bar on a log x-axis: the counts
    span ~20 to ~500, and on a linear axis the bifunctional and polyfunctional
    bars -- the ones the reader came for -- would be invisible slivers. Each
    bar is direct-labelled with its count and percentage so the log axis never
    has to be mentally un-warped.
  * The verdict panel is 100% stacked per group, with group n in the tick
    label, because the comparison that matters is proportion-within-group
    (62 vs 1,308 reactions) not absolute height.
  * Every number is read from the JSON, never retyped; the same values are
    written beside the images as *_values.tsv for a diffable record.

Data: ModelSEED dev @ 078a395f (/scratch/ctaylor/tmp/devsnap_078a395f).
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
SUMMARY = ANALYSIS_DIR / "results" / "strict_disagreement_enzymes" / "summary.json"
OUT_DIR = ANALYSIS_DIR / "reports" / "strictDisagreementEnzymes" / "figures"

INK, INK2 = "#0b0b0b", "#52514e"
SURFACE, GRIDLINE = "#fcfcfb", "#e1e0d9"
GRAY, GRAY2 = "#c3c2b7", "#e1e0d9"
BLUE, ORANGE, GREEN = "#2a78d6", "#eb6834", "#1baf7a"

FUNC_ORDER = ["monofunctional", "bifunctional", "polyfunctional", "unknown", "not_researched"]
FUNC_COLOR = {"monofunctional": BLUE, "bifunctional": ORANGE,
              "polyfunctional": GREEN, "unknown": GRAY, "not_researched": GRAY2}
FUNC_LABEL = {"monofunctional": "monofunctional", "bifunctional": "bifunctional",
              "polyfunctional": "polyfunctional", "unknown": "unknown\n(no enzyme identifiable)",
              "not_researched": "not researched"}

VERDICT_ORDER = ["llm", "modelseed", "neither", "unclear"]
VERDICT_COLOR = {"llm": BLUE, "modelseed": ORANGE, "neither": GREEN, "unclear": GRAY}
VERDICT_LABEL = {"llm": "LLM council right", "modelseed": "ModelSEED right",
                 "neither": "genuinely reversible /\npool-determined", "unclear": "unclear"}

DEFECT_ORDER = ["none", "cofactors_stripped", "mass_unbalanced",
                "lumped_multienzyme", "source_disowned"]
DEFECT_LABEL = {"none": "equation sound", "cofactors_stripped": "cofactors stripped",
                "mass_unbalanced": "mass unbalanced",
                "lumped_multienzyme": "lumped multi-enzyme",
                "source_disowned": "source disowns entry"}
DEFECT_COLOR = {"none": GREEN, "cofactors_stripped": BLUE, "mass_unbalanced": ORANGE,
                "lumped_multienzyme": "#8c5ad1", "source_disowned": GRAY}


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(GRIDLINE)
    ax.tick_params(colors=INK2, labelsize=9)


def fig_functionality(s, out_tsv):
    counts = s["by_functionality"]
    total = s["total_reactions"]
    order = [c for c in FUNC_ORDER if counts.get(c)]
    vals = [counts[c] for c in order]

    fig, ax = plt.subplots(figsize=(9.6, 4.2))
    fig.patch.set_facecolor(SURFACE)
    style(ax)
    y = range(len(order))
    ax.barh(list(y), vals, height=0.6, color=[FUNC_COLOR[c] for c in order],
            edgecolor="none", zorder=2)
    ax.set_xscale("log")
    ax.grid(True, axis="x", color=GRIDLINE, linewidth=0.8, zorder=0)
    ax.set_yticks(list(y))
    ax.set_yticklabels([FUNC_LABEL[c] for c in order], fontsize=9.5, color=INK)
    ax.invert_yaxis()
    for i, (c, v) in enumerate(zip(order, vals)):
        ax.text(v * 1.08, i, f"{v:,}  ({100*v/total:.1f}%)", va="center",
                fontsize=9.5, color=INK, fontweight="bold" if c in
                ("bifunctional", "polyfunctional") else "normal")
    ax.set_xlim(right=max(vals) * 3.2)
    ax.set_xlabel("reactions  (log scale)", fontsize=9, color=INK2)
    multi = counts.get("bifunctional", 0) + counts.get("polyfunctional", 0)
    researched = s["researched"]
    ax.set_title(
        f"Enzyme functionality behind the {total:,} strict direction disagreements\n"
        f"Only {multi:,} ({100*multi/researched:.1f}%) are run by a "
        f"bi- or polyfunctional enzyme",
        fontsize=11.5, color=INK, loc="left", pad=12)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "strict_disagreement_functionality.png", dpi=200,
                facecolor=SURFACE)
    plt.close(fig)

    with open(out_tsv, "w") as fh:
        fh.write("category\treactions\tpct_of_total\n")
        for c in FUNC_ORDER:
            v = counts.get(c, 0)
            fh.write(f"{c}\t{v}\t{100*v/total:.2f}\n")


def fig_direction(s, out_tsv):
    by_group = s["verdict_by_group"]
    groups = [g for g in ("reverse_vs_forward", "forward_vs_reverse") if by_group.get(g)]
    defects = s["by_equation_defect"]
    enr = s["modelseed_imbalance_enrichment"]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.4, 4.6),
                                   gridspec_kw={"width_ratios": [1.15, 1.0]})
    fig.patch.set_facecolor(SURFACE)

    # --- left: verdict composition per mismatch group -------------------
    style(axL)
    for i, g in enumerate(groups):
        v = by_group[g]
        n = sum(v.values())
        left = 0.0
        for k in VERDICT_ORDER:
            if not v.get(k):
                continue
            w = 100.0 * v[k] / n
            axL.barh(i, w, left=left, height=0.55, color=VERDICT_COLOR[k],
                     edgecolor="none", zorder=2)
            if w >= 7:
                axL.text(left + w / 2, i, f"{w:.0f}%", ha="center", va="center",
                         fontsize=9.5, color="white", fontweight="bold")
            left += w
    axL.set_yticks(range(len(groups)))
    axL.set_yticklabels([f"{g.replace('_', ' ')}\n(n={sum(by_group[g].values()):,})"
                         for g in groups], fontsize=9.5, color=INK)
    axL.invert_yaxis()
    axL.set_xlim(0, 100)
    axL.set_xlabel("share of researched reactions in the group (%)", fontsize=9, color=INK2)
    axL.grid(True, axis="x", color=GRIDLINE, linewidth=0.8, zorder=0)
    axL.set_title("Which side does the enzymology support?", fontsize=11,
                  color=INK, loc="left", pad=10)
    handles = [plt.Rectangle((0, 0), 1, 1, color=VERDICT_COLOR[k]) for k in VERDICT_ORDER]
    axL.legend(handles, [VERDICT_LABEL[k] for k in VERDICT_ORDER],
               loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2,
               frameon=False, fontsize=8.5, labelcolor=INK2)

    # --- right: equation defect mix -------------------------------------
    style(axR)
    order = [d for d in DEFECT_ORDER if defects.get(d)]
    vals = [defects[d] for d in order]
    tot = sum(vals)
    axR.bar(range(len(order)), vals, width=0.6,
            color=[DEFECT_COLOR[d] for d in order], edgecolor="none", zorder=2)
    axR.set_xticks(range(len(order)))
    axR.set_xticklabels([DEFECT_LABEL[d] for d in order], fontsize=8.5,
                        color=INK, rotation=20, ha="right")
    axR.grid(True, axis="y", color=GRIDLINE, linewidth=0.8, zorder=0)
    for i, v in enumerate(vals):
        axR.text(i, v, f"{v}\n{100*v/tot:.0f}%", ha="center", va="bottom",
                 fontsize=8.5, color=INK)
    axR.set_ylim(top=max(vals) * 1.28)
    axR.set_ylabel("reactions", fontsize=9, color=INK2)
    axR.set_title(
        "Is the ModelSEED equation itself sound?\n"
        f"Its own imbalance flag fires on {enr['strict_disagreements']['pct']}% "
        f"here, {enr['all_non_obsolete']['pct']}% database-wide",
        fontsize=10.5, color=INK, loc="left", pad=10)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "strict_disagreement_direction.png", dpi=200,
                facecolor=SURFACE)
    plt.close(fig)

    with open(out_tsv, "w") as fh:
        fh.write("panel\tkey\tgroup\tvalue\n")
        for g in groups:
            for k, v in by_group[g].items():
                fh.write(f"verdict\t{k}\t{g}\t{v}\n")
        for d in DEFECT_ORDER:
            fh.write(f"defect\t{d}\t-\t{defects.get(d, 0)}\n")
        for key, blk in enr.items():
            if isinstance(blk, dict) and "pct" in blk:
                fh.write(f"imbalance\t{key}\t-\t{blk['flagged']}/{blk['n']} ({blk['pct']}%)\n")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(SUMMARY) as fh:
        s = json.load(fh)
    fig_functionality(s, OUT_DIR / "strict_disagreement_functionality_values.tsv")
    fig_direction(s, OUT_DIR / "strict_disagreement_direction_values.tsv")
    print(f"wrote 2 figures + 2 value TSVs to {OUT_DIR}")


if __name__ == "__main__":
    main()
