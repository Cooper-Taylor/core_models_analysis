#!/usr/bin/env python3
"""The 5 systematic, non-overlapping clusters (+ residual) behind the 465
"LLM-right, equation-sound" strict disagreements.

Reads results/strict_disagreement_enzymes/sound_subset_clusters.tsv (written
by analyze_strict_disagreement_sound_subset.py). Each of the first 5 rows is a
reaction cluster claimed in priority order (no double-counting); the 6th row,
"Residual", is everything left over. Bars are grouped and coloured by DEFECT
CLASS, not by enzyme family -- two of the five clusters (ATP/phosphate-coupled
and quinone/quinol-coupled) were found by searching the reaction EQUATION for
a cofactor/motif signature, not the enzyme name, and cut across many
differently-named enzyme families, so an enzyme-family bar chart (the
previous version of this figure) couldn't show them.

  dGPredictor cluster (orange): SAM methyltransferase, PAPS sulfotransferase,
    ATP/phosphate-coupled group transfer -- mechanistically the same category
    (an activated, high-group-transfer-potential cofactor), so one bad group
    contribution for that bond type plausibly explains all three.
  eQuilibrator cluster (blue): glycoside hydrolase (broad -- glucosidase,
    glucuronidase, rhamnosidase, ... all the same GH-family mechanism),
    quinone/quinol-coupled redox (independently flagged elsewhere in this
    project as eQuilibrator's weak spot).
  Residual (grey): case-specific enzymology with no further single-mechanism
    pattern found -- its thermo grade/source/confidence profile is
    statistically indistinguishable from the full 504, and a tested
    secondary-metabolism/molecular-novelty hypothesis does NOT explain it
    (that signal concentrates in the clusters, not here -- see
    sound_subset_analysis.json -> systematic_clusters_residual).

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
DATA = ANALYSIS_DIR / "results" / "strict_disagreement_enzymes" / "sound_subset_clusters.tsv"
OUT_DIR = ANALYSIS_DIR / "reports" / "strictDisagreementEnzymes" / "figures"

INK, INK2, INK3 = "#0b0b0b", "#52514e", "#898781"
SURFACE, GRIDLINE = "#fcfcfb", "#e1e0d9"
# categorical slots: dGPredictor = slot 2 (orange), eQuilibrator = slot 1 (blue) --
# same pair used across the LLM_VS_RECOMMENDED / strict-disagreement figures.
# Residual gets the neutral "no data" grey, not a 3rd categorical hue -- it's
# a leftover bucket, not a defect class competing for identity.
CLASS_COLOR = {
    "dGPredictor: activated-cofactor group transfer": "#eb6834",
    "eQuilibrator: structural motif": "#2a78d6",
    "none": "#c3c2b7",
}
CLASS_LABEL = {
    "dGPredictor: activated-cofactor group transfer": "dGPredictor cluster",
    "eQuilibrator: structural motif": "eQuilibrator cluster",
    "none": "residual (case-specific)",
}
SYSTEMATIC_PURITY_THRESHOLD = 0.95

PRETTY = {
    "SAM-dependent methyltransferase": "SAM methyltransferase",
    "PAPS-dependent sulfotransferase": "PAPS sulfotransferase",
    "Glycoside hydrolase (broad)": "Glycoside hydrolase\n(glucosidase/glucuronidase/...)",
    "Quinone/quinol-coupled redox": "Quinone/quinol-coupled redox",
    "ATP/ADP/phosphate-coupled group transfer": "ATP/phosphate-coupled\ngroup transfer",
    "Residual (case-specific enzymology)": "Residual\n(case-specific enzymology)",
}


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    with open(DATA) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            purity = float(r["dominant_source_purity"]) if r["dominant_source_purity"] else None
            rows.append({"cluster": r["cluster"], "n": int(r["n"]), "defect_class": r["defect_class"],
                         "source": r["dominant_source"], "purity": purity})
    total = sum(r["n"] for r in rows)
    # Group bars by defect class (not just colour) so the two-defect-class
    # story reads spatially: eQuilibrator cluster at the bottom, dGPredictor
    # cluster above it, residual on top. Smallest-to-largest within each
    # group, since this list is consumed bottom-to-top by barh.
    residual = [r for r in rows if r["defect_class"] == "none"]
    eq_cluster = sorted([r for r in rows if r["defect_class"].startswith("eQuilibrator")], key=lambda r: r["n"])
    dgp_cluster = sorted([r for r in rows if r["defect_class"].startswith("dGPredictor")], key=lambda r: r["n"])
    ordered = eq_cluster + dgp_cluster + residual

    fig, ax = plt.subplots(figsize=(9.8, 6.2))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    y = range(len(ordered))
    colors = [CLASS_COLOR[r["defect_class"]] for r in ordered]
    edge_colors = [INK if (r["purity"] or 0) >= SYSTEMATIC_PURITY_THRESHOLD else "none" for r in ordered]
    edge_widths = [1.8 if (r["purity"] or 0) >= SYSTEMATIC_PURITY_THRESHOLD else 0 for r in ordered]
    ax.barh(list(y), [r["n"] for r in ordered], color=colors,
            edgecolor=edge_colors, linewidth=edge_widths, zorder=3)

    ax.set_yticks(list(y))
    ax.set_yticklabels([PRETTY.get(r["cluster"], r["cluster"]) for r in ordered], color=INK2, fontsize=10.0)
    ax.tick_params(axis="x", colors=INK2, labelsize=9.6)
    ax.grid(True, axis="x", color=GRIDLINE, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color(GRIDLINE)
    ax.set_xlabel("reactions (LLM-correct, equation-sound)", color=INK2, fontsize=10.2)

    xmax = max(r["n"] for r in ordered)
    ax.set_xlim(0, xmax * 1.30)
    for yi, r in zip(y, ordered):
        if r["purity"] is None:
            label = f"{r['n']}  ({100*r['n']/total:.0f}% of 465)"
        else:
            label = f"{r['n']}  ({r['purity']:.0%} {r['source']})"
        ax.text(r["n"] + xmax * 0.012, yi, label, va="center", ha="left", color=INK, fontsize=9.6,
                fontweight="bold" if (r["purity"] or 0) >= SYSTEMATIC_PURITY_THRESHOLD else "normal")

    fig.suptitle("Two thermo-source defect classes, and what's left over",
                 x=0.09, y=0.975, ha="left", color=INK, fontsize=13.6, fontweight="semibold")
    n_dgp = sum(r["n"] for r in rows if r["defect_class"].startswith("dGPredictor"))
    n_eq = sum(r["n"] for r in rows if r["defect_class"].startswith("eQuilibrator"))
    ax.set_title(f"465 LLM-correct, equation-sound reactions  ·  dGPredictor cluster "
                 f"{n_dgp} ({100*n_dgp/total:.0f}%)  ·  eQuilibrator cluster {n_eq} ({100*n_eq/total:.0f}%)",
                 loc="left", color=INK2, fontsize=9.6, fontweight="normal", pad=12)

    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in CLASS_COLOR.values()]
    ax.legend(handles, [CLASS_LABEL[k] for k in CLASS_COLOR],
              loc="lower right", frameon=False, fontsize=9.4, labelcolor=INK2)

    fig.subplots_adjust(left=0.24, right=0.97, top=0.83, bottom=0.10)
    out_png = OUT_DIR / "strict_disagreement_sound_families.png"
    fig.savefig(out_png, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
