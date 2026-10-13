#!/usr/bin/env python3
"""Which side the enzymology supports, among the strict disagreements --
restricted to the OK/CI-only population (see
analyze_strict_disagreement_balanced_correctness.py for the filter and why).

Same 100%-stacked-bar design and verdict categories/colours as the left panel
of plot_strict_disagreement_enzymes.py's strict_disagreement_direction.png;
this is a standalone figure (not a second panel) because the right-hand
equation-defect panel in that figure doesn't carry over cleanly once the
mass-imbalanced reactions (its single largest defect category) have already
been removed from the population.

Reads results/strict_disagreement_enzymes/correctness_balanced.json (written
by analyze_strict_disagreement_balanced_correctness.py). Each group's bar
shows both populations stacked as grouped rows (all 1,370 vs. the 879
OK/CI-only) so the "did filtering change the verdict?" comparison is visible
in one chart rather than requiring two images side by side.

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
DATA = ANALYSIS_DIR / "results" / "strict_disagreement_enzymes" / "correctness_balanced.json"
OUT_DIR = ANALYSIS_DIR / "reports" / "strictDisagreementEnzymes" / "figures"

INK, INK2 = "#0b0b0b", "#52514e"
SURFACE, GRIDLINE = "#fcfcfb", "#e1e0d9"
GRAY = "#c3c2b7"
BLUE, ORANGE, GREEN = "#2a78d6", "#eb6834", "#1baf7a"

VERDICT_ORDER = ["llm", "modelseed", "neither", "unclear"]
VERDICT_COLOR = {"llm": BLUE, "modelseed": ORANGE, "neither": GREEN, "unclear": GRAY}
VERDICT_LABEL = {"llm": "LLM council right", "modelseed": "ModelSEED right",
                 "neither": "genuinely reversible /\npool-determined", "unclear": "unclear"}

GROUPS = ("reverse_vs_forward", "forward_vs_reverse")
POPULATIONS = [("all_1370", "all 1,370"), ("OK_CI_only", "OK/CI-only")]


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(GRIDLINE)
    ax.tick_params(colors=INK2, labelsize=9)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    data = json.load(open(DATA))

    rows = []  # (label, group, pop_key)
    for g in GROUPS:
        for pop_key, pop_label in POPULATIONS:
            rows.append((g, pop_key, pop_label))

    fig, ax = plt.subplots(figsize=(10.2, 5.6))
    fig.patch.set_facecolor(SURFACE)
    style(ax)

    for i, (g, pop_key, pop_label) in enumerate(rows):
        block = data["by_group"][pop_key][g]
        n = block["n"]
        left = 0.0
        for v in VERDICT_ORDER:
            w = block["pct"][v]
            if w <= 0:
                continue
            ax.barh(i, w, left=left, height=0.62, color=VERDICT_COLOR[v],
                    edgecolor="none", zorder=2)
            if w >= 6:
                ax.text(left + w / 2, i, f"{w:.0f}%", ha="center", va="center",
                        fontsize=9.3, color="white", fontweight="bold")
            left += w

    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([
        f"{g.replace('_', ' ')}  —  {pop_label}  (n={data['by_group'][pop_key][g]['n']:,})"
        for g, pop_key, pop_label in rows
    ], fontsize=9.3, color=INK)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("share of researched reactions in the group (%)", fontsize=9, color=INK2)
    ax.grid(True, axis="x", color=GRIDLINE, linewidth=0.8, zorder=0)

    # faint separator between the two mismatch groups
    ax.axhline(1.5, color=GRIDLINE, linewidth=1.0, zorder=1)

    n_kept, n_excl = data["n_kept"], data["n_excluded"]
    ax.set_title(
        "Which side does the enzymology support?\n"
        f"Excluding {n_excl:,} mass-imbalanced reactions (status MI / CI+MI)\n"
        f"barely moves the verdict — {n_kept:,} of 1,370 remain",
        fontsize=11.5, color=INK, loc="left", pad=12)
    handles = [plt.Rectangle((0, 0), 1, 1, color=VERDICT_COLOR[v]) for v in VERDICT_ORDER]
    ax.legend(handles, [VERDICT_LABEL[v] for v in VERDICT_ORDER],
              loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=4,
              frameon=False, fontsize=8.8, labelcolor=INK2)

    fig.tight_layout()
    out_png = OUT_DIR / "strict_disagreement_direction_balanced.png"
    fig.savefig(out_png, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out_png}")

    tsv_path = OUT_DIR / "strict_disagreement_direction_balanced_values.tsv"
    with open(tsv_path, "w") as fh:
        fh.write("population\tgroup\tn\t" + "\t".join(VERDICT_ORDER) + "\t" +
                  "\t".join(f"{v}_pct" for v in VERDICT_ORDER) + "\n")
        for pop_key, pop_label in POPULATIONS:
            for g in GROUPS:
                b = data["by_group"][pop_key][g]
                fh.write(f"{pop_key}\t{g}\t{b['n']}\t" +
                         "\t".join(str(b['counts'][v]) for v in VERDICT_ORDER) + "\t" +
                         "\t".join(str(b['pct'][v]) for v in VERDICT_ORDER) + "\n")
    print(f"wrote {tsv_path}")


if __name__ == "__main__":
    main()
