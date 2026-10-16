#!/usr/bin/env python3
"""Funnel chart for find_strict_disagreement_redox_fixes.py's coverage cascade:
of the 219 CI-only strict disagreements, how many get a mechanically-verified
"fixed" equation, how many of those map to a real eQuilibrator compound, how
many give a physically plausible |dG| (<=200 kJ/mol -- the generic, templated
"aldehyde dehydrogenase [NAD(P)+]" combinatorial-expansion family gives
280-690 kJ/mol otherwise, a mapping artifact, not chemistry), and of the
plausible ones, how many give a confident directional call (vs near-
equilibrium) and which side (LLM / ModelSEED / neither) that call matches.

This is a coverage funnel, not a distribution -- every stage is a strict
subset of the one above it, so bars shrink monotonically left to right (drawn
top to bottom here). The terminal result is small by construction (7 of 219
reactions could be mechanically fixed AND independently priced by a live
thermodynamic calculation) -- the figure's job is to make that attrition
honest, not to hide it inside a single percentage.

Data: results/strict_disagreement_enzymes/redox_fix_summary.json, written by
find_strict_disagreement_redox_fixes.py (run in the `eq3` conda env). Same
ModelSEED dev @ 078a395f snapshot as the rest of this pipeline.
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
DATA = ANALYSIS_DIR / "results" / "strict_disagreement_enzymes" / "redox_fix_summary.json"
OUT_DIR = ANALYSIS_DIR / "reports" / "strictDisagreementEnzymes" / "figures"

INK, INK2, INK3 = "#0b0b0b", "#52514e", "#898781"
SURFACE, GRIDLINE = "#fcfcfb", "#e1e0d9"
BLUE = "#2a78d6"   # agrees with LLM
ORANGE = "#eb6834"  # agrees with ModelSEED (none, here)
GREY = "#c3c2b7"    # dropped at this stage / near-equilibrium


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    s = json.load(open(DATA))

    stages = [
        ("CI-only strict disagreements", s["n_ci_only"], GREY),
        ("mechanically-verified\ncofactor fix found", s["n_solved"], GREY),
        ("fix maps to a real\neQuilibrator compound", s["n_eq3_computed_for_fixed_equation"], GREY),
        ("...with a physically\nplausible |ΔG′°|", s["n_eq3_computed_for_fixed_equation"] - s["n_implausible_magnitude_excluded"], GREY),
    ]
    n_plausible = s["n_eq3_computed_for_fixed_equation"] - s["n_implausible_magnitude_excluded"]
    terminal = [
        ("near-equilibrium\n(no confident call)", s["n_near_equilibrium"], GREY),
        ("confident call →\nagrees with LLM", s["n_fixed_agrees_with_llm"], BLUE),
        ("confident call →\nagrees with ModelSEED", s["n_fixed_agrees_with_modelseed"], ORANGE),
    ]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.6, 5.0), gridspec_kw={"width_ratios": [1.3, 1]})
    fig.patch.set_facecolor(SURFACE)

    for ax, rows, xlabel, title in (
        (ax1, stages, "reactions (of 219)", "Coverage funnel"),
        (ax2, terminal, f"reactions (of {n_plausible} plausible)", "Of the physically plausible 7"),
    ):
        ax.set_facecolor(SURFACE)
        y = range(len(rows))
        xmax = max(n for _, n, _ in rows) or 1
        ax.barh(list(y), [n for _, n, _ in rows], color=[c for _, _, c in rows], zorder=3)
        ax.set_yticks(list(y))
        ax.set_yticklabels([lbl for lbl, _, _ in rows], color=INK2, fontsize=9.6)
        ax.invert_yaxis()
        for yi, (_, n, _) in zip(y, rows):
            ax.text(n + xmax * 0.03, yi, str(n), va="center", ha="left", color=INK, fontsize=10.5, fontweight="bold")
        ax.set_xlim(0, xmax * 1.22)
        ax.tick_params(axis="x", colors=INK2, labelsize=9.0)
        ax.grid(True, axis="x", color=GRIDLINE, linewidth=0.8, zorder=0)
        ax.set_axisbelow(True)
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
        ax.spines["bottom"].set_color(GRIDLINE)
        ax.set_xlabel(xlabel, color=INK2, fontsize=9.8)
        ax.set_title(title, loc="left", color=INK, fontsize=11.6, fontweight="semibold", pad=10)

    fig.suptitle("Can a mechanically-fixed reaction's real thermodynamics be computed,\n"
                 "and does it side with the LLM?",
                 x=0.045, y=0.995, ha="left", color=INK, fontsize=14.0, fontweight="semibold")
    fig.text(0.045, 0.015,
              "219 CI-only strict disagreements -> 103 get a verified, mass+charge-exact cofactor-bookkeeping fix (no new chemistry invented) ->\n"
              "29 map to a compound eQuilibrator recognises -> 7 give a plausible |ΔG′°| (the other 22 are one templated, mismapped reaction family) ->\n"
              "of 5 with a confident directional call, 5/5 (100%) agree with the LLM, 0 with ModelSEED's original (broken) recommendation.",
              color=INK3, fontsize=8.6, ha="left")

    fig.subplots_adjust(left=0.22, right=0.97, top=0.80, bottom=0.17, wspace=0.45)
    out_png = OUT_DIR / "strict_disagreement_redox_fix_funnel.png"
    fig.savefig(out_png, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
