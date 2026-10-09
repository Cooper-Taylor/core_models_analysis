#!/usr/bin/env python3
"""Biological pathway/network breakdown of the 248 "residual" strict
disagreements -- the reactions explained by NEITHER of the two systematic
thermo-source defect classes (see plot_strict_disagreement_sound_families.py).

Reads results/strict_disagreement_enzymes/residual_pathways.tsv (written by
analyze_strict_disagreement_residual_pathways.py). Bars are coloured by
provenance mix: how much of each category's count came from a curated
KEGG/MetaCyc `pathways` annotation vs. a weaker text-mining fallback (enzyme
name/rationale/organism context), shown as a stacked bar so the evidence
strength is visible at a glance. "No pathway identifiable" (grey) is its own
bar, not a 3rd provenance colour -- it's a leftover, not evidence.

Data: ModelSEED dev @ 078a395f, strict-disagreement literature pass
2026-10-08 (results/strict_disagreement_enzymes/).
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
DATA = ANALYSIS_DIR / "results" / "strict_disagreement_enzymes" / "residual_pathways.tsv"
OUT_DIR = ANALYSIS_DIR / "reports" / "strictDisagreementEnzymes" / "figures"

INK, INK2, INK3 = "#0b0b0b", "#52514e", "#898781"
SURFACE, GRIDLINE = "#fcfcfb", "#e1e0d9"
# categorical slot 3 (green) for curated pathway annotation, a lighter tint of
# the same hue for the weaker text-inferred evidence -- one series, one hue,
# light->dark by confidence, not two competing categorical colours. "No
# pathway identifiable" gets the neutral "no data" grey, same convention as
# the other figures in this report.
ANNOTATED_COLOR = "#1baf7a"
INFERRED_COLOR = "#a8e0c9"
NO_PATHWAY_COLOR = "#c3c2b7"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    with open(DATA) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append({"category": r["category"], "n": int(r["n"]),
                         "annotated": int(r["n_pathway_annotated"]), "inferred": int(r["n_text_inferred"])})
    total = sum(r["n"] for r in rows)
    rows.sort(key=lambda r: r["n"])  # ascending -> largest at top of barh

    fig, ax = plt.subplots(figsize=(11.5, 6.2))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    y = range(len(rows))
    no_pathway = [r["n"] if r["category"] == "No pathway identifiable" else 0 for r in rows]
    annotated = [r["annotated"] for r in rows]
    inferred = [r["inferred"] for r in rows]

    ax.barh(list(y), no_pathway, color=NO_PATHWAY_COLOR, zorder=3, label="no pathway identifiable")
    ax.barh(list(y), annotated, color=ANNOTATED_COLOR, zorder=3, label="KEGG/MetaCyc pathway annotation")
    ax.barh(list(y), inferred, left=annotated, color=INFERRED_COLOR, zorder=3,
            label="text-inferred (no curated annotation)")

    ax.set_yticks(list(y))
    ax.set_yticklabels([r["category"] for r in rows], color=INK2, fontsize=10.2)
    ax.tick_params(axis="x", colors=INK2, labelsize=9.6)
    ax.grid(True, axis="x", color=GRIDLINE, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color(GRIDLINE)
    ax.set_xlabel("reactions (of the 248 residual)", color=INK2, fontsize=10.2)

    xmax = max(r["n"] for r in rows)
    ax.set_xlim(0, xmax * 1.2)
    for yi, r in zip(y, rows):
        ax.text(r["n"] + xmax * 0.012, yi, f"{r['n']}  ({100*r['n']/total:.0f}%)",
                va="center", ha="left", color=INK, fontsize=9.8,
                fontweight="bold" if r["n"] >= 40 else "normal")

    fig.suptitle("The 248 residual reactions cluster into two biological networks",
                 x=0.09, y=0.975, ha="left", color=INK, fontsize=13.6, fontweight="semibold")
    ax.set_title("categorised by KEGG/MetaCyc pathway annotation, with a text-mining fallback",
                 loc="left", color=INK2, fontsize=9.6, fontweight="normal", pad=12)

    ax.legend(loc="lower right", frameon=False, fontsize=9.2, labelcolor=INK2)

    fig.subplots_adjust(left=0.37, right=0.97, top=0.84, bottom=0.10)
    out_png = OUT_DIR / "strict_disagreement_residual_pathways.png"
    fig.savefig(out_png, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
