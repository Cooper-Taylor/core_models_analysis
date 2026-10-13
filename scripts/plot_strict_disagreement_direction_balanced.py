#!/usr/bin/env python3
"""Literature-determined verdict counts for the strict direction disagreements
-- all 1,370 vs. the OK/CI-only (mass-balanced) subset of 879 -- split into
the two mismatch groups, which must never be pooled (see
reports/strictDisagreementEnzymes/STRICT_DISAGREEMENT_ENZYMES.md section 1).

Reads results/strict_disagreement_enzymes/correctness_balanced.json (written
by analyze_strict_disagreement_balanced_correctness.py; no new research --
every verdict is `direction_conflict.enzymology_supports`, already collected
by the literature subagents behind strict_disagreement_enzymes.json).

DESIGN NOTES
  * Raw reaction COUNTS, not percentages (the earlier version of this figure
    used 100%-stacked bars; a user request specifically wanted the absolute
    numbers back, since a percentage hides how much smaller forward_vs_reverse
    is than reverse_vs_forward).
  * Two panels with INDEPENDENT y-axis scaling, not one shared-scale chart:
    forward_vs_reverse tops out at 31, reverse_vs_forward at 1,221 (~40x) --
    on a shared axis the smaller group's bars would be unreadable slivers.
  * Grouped bars (all 1,370 vs OK/CI-only) share each verdict category's hue;
    the filtered subset is drawn as a lighter tint + hatch of the same colour
    (it is literally a subset of "all", not an independent category), rather
    than introducing a second categorical hue dimension for "population".
  * The group name convention (`forward_vs_reverse` = ModelSEED says forward,
    LLM says reverse) is easy to misread backwards, so each panel's subtitle
    spells out "ModelSEED recommends X -- LLM council says Y" explicitly
    instead of relying on the reader to parse the underscore-joined key.

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
VERDICT_LABEL = {"llm": "LLM council\nright", "modelseed": "ModelSEED\nright",
                 "neither": "genuinely\nreversible", "unclear": "unclear"}

# (group_key, recommended_symbol, recommended_word, llm_symbol, llm_word)
PANELS = [
    ("forward_vs_reverse", ">", "forward", "<", "reverse"),
    ("reverse_vs_forward", "<", "reverse", ">", "forward"),
]
POPULATIONS = [("all_1370", "all 1,370", 1.0, False), ("OK_CI_only", "OK/CI-only (mass-balanced)", 0.55, True)]

# The group key's first word is always ModelSEED's recommended direction, the
# second is the LLM ensemble's -- easy to misread backwards, so pin it down.
for _group_key, _rec_sym, _rec_word, _llm_sym, _llm_word in PANELS:
    assert _group_key == f"{_rec_word}_vs_{_llm_word}", _group_key
    assert (_rec_sym, _llm_sym) in {(">", "<"), ("<", ">")}, _group_key


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(GRIDLINE)
    ax.tick_params(colors=INK2, labelsize=9.5)


def _tint(hex_colour: str, amount: float) -> tuple:
    r, g, b = (int(hex_colour[i:i + 2], 16) / 255 for i in (1, 3, 5))
    return tuple(c + (1.0 - c) * amount for c in (r, g, b))


def draw_panel(ax, data, group_key, rec_sym, rec_word, llm_sym, llm_word):
    style(ax)
    n_cat = len(VERDICT_ORDER)
    width = 0.34
    x = list(range(n_cat))

    max_val = 0
    for j, (pop_key, pop_label, tint_amount, hatched) in enumerate(POPULATIONS):
        block = data["by_group"][pop_key][group_key]
        offset = (j - 0.5) * width
        for i, v in enumerate(VERDICT_ORDER):
            count = block["counts"][v]
            pct = block["pct"][v]
            max_val = max(max_val, count)
            face = _tint(VERDICT_COLOR[v], tint_amount) if hatched else VERDICT_COLOR[v]
            ax.bar(x[i] + offset, count, width=width * 0.92, color=face,
                   edgecolor=VERDICT_COLOR[v], linewidth=1.1,
                   hatch="///" if hatched else None, zorder=2)
            ax.text(x[i] + offset, count + max(1, 0.015 * (max_val or 1)),
                    f"{count:,}", ha="center", va="bottom", fontsize=8.8,
                    color=INK, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([VERDICT_LABEL[v] for v in VERDICT_ORDER], fontsize=9.3, color=INK)
    ax.grid(True, axis="y", color=GRIDLINE, linewidth=0.8, zorder=0)
    ax.set_ylim(top=max_val * 1.22)
    ax.set_ylabel("reactions", fontsize=9.5, color=INK2)

    n_all = data["by_group"]["all_1370"][group_key]["n"]
    n_kept = data["by_group"]["OK_CI_only"][group_key]["n"]
    ax.set_title(
        f"ModelSEED recommends {rec_word} ({rec_sym})  —  "
        f"LLM council says {llm_word} ({llm_sym})\n"
        f"n = {n_all:,} (all)  /  {n_kept:,} (OK/CI-only)",
        fontsize=10.8, color=INK, loc="left", pad=10)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    data = json.load(open(DATA))

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.8, 6.2))
    fig.patch.set_facecolor(SURFACE)

    draw_panel(axL, data, *PANELS[0])
    draw_panel(axR, data, *PANELS[1])

    fig.suptitle(
        "Strict-disagreement literature verdicts: reaction counts, all 1,370 vs. mass-balanced only\n"
        "Group names read ModelSEED direction → LLM direction, in that order",
        fontsize=12.3, color=INK, x=0.03, ha="left", y=0.995, fontweight="semibold")

    legend_handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=GRAY, edgecolor=INK, linewidth=1.1),
        plt.Rectangle((0, 0), 1, 1, facecolor=_tint(GRAY, 0.55), edgecolor=INK,
                      linewidth=1.1, hatch="///"),
    ]
    fig.legend(legend_handles, ["all 1,370 (unfiltered)", "OK/CI-only, n=879 (mass-balanced subset)"],
               loc="lower center", ncol=2, frameon=False, fontsize=9.3,
               labelcolor=INK2, bbox_to_anchor=(0.5, -0.01))

    fig.tight_layout(rect=(0, 0.045, 1, 0.88))
    out_png = OUT_DIR / "strict_disagreement_direction_balanced.png"
    fig.savefig(out_png, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out_png}")

    tsv_path = OUT_DIR / "strict_disagreement_direction_balanced_values.tsv"
    with open(tsv_path, "w") as fh:
        fh.write("population\tgroup\tn\t" + "\t".join(VERDICT_ORDER) + "\t" +
                  "\t".join(f"{v}_pct" for v in VERDICT_ORDER) + "\n")
        for pop_key, pop_label, _, _ in POPULATIONS:
            for group_key, *_ in PANELS:
                b = data["by_group"][pop_key][group_key]
                fh.write(f"{pop_key}\t{group_key}\t{b['n']}\t" +
                         "\t".join(str(b['counts'][v]) for v in VERDICT_ORDER) + "\t" +
                         "\t".join(str(b['pct'][v]) for v in VERDICT_ORDER) + "\n")
    print(f"wrote {tsv_path}")


if __name__ == "__main__":
    main()
