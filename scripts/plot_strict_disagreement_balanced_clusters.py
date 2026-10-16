#!/usr/bin/env python3
"""Do the two systematic defect clusters survive a stricter population
definition? Three small-multiple panels, same bars/colours/order as
strict_disagreement_sound_families.png, one per population:

  all_1370    the original 465 LLM-correct/equation-sound reactions
              (equation_defect == "none" across all 1,370 strict disagreements)
  ok_ci_only  448 LLM-correct/equation-sound reactions within the 879 whose
              ModelSEED `status` is OK or CI-only (charge-imbalanced but mass
              already balanced) -- i.e. mass-imbalanced reactions (MI, CI+MI)
              excluded by the DATABASE'S OWN flag, not by literature judgment
  ok_only     438 LLM-correct/equation-sound reactions within the 660 whose
              status carries no flag at all (OK, nothing else)

These three "sound" populations are NOT nested subsets of each other (the
equation_defect label is a literature judgment, modelseed_status is an
automated check, and they disagree on 18 reactions -- see
analyze_strict_disagreement_balanced_populations.py's docstring), so this
figure is the actual test of whether the two-cluster story is an artifact of
one particular sound/unsound split. It isn't: SAM (81), PAPS (35), and
glycosidase (54) are IDENTICAL across all three panels, quinone drops by
exactly 1 (one quinone reaction is CI-flagged and so absent from ok_only),
and ATP (24) and the residual (465/248 -> 448/242 -> 438/233) move only by
the small number of reactions that carry an MI/CI+MI status despite a
literature-sound equation, or vice versa.

Data: results/strict_disagreement_enzymes/balanced_population_clusters.tsv,
written by analyze_strict_disagreement_balanced_populations.py. Same
ModelSEED dev @ 078a395f strict-disagreement literature pass (2026-10-08) as
strict_disagreement_sound_families.png -- no new research, just re-tabulated
over different populations.
"""
from __future__ import annotations

import csv
import os
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
DATA = ANALYSIS_DIR / "results" / "strict_disagreement_enzymes" / "balanced_population_clusters.tsv"
OUT_DIR = ANALYSIS_DIR / "reports" / "strictDisagreementEnzymes" / "figures"

INK, INK2, INK3 = "#0b0b0b", "#52514e", "#898781"
SURFACE, GRIDLINE = "#fcfcfb", "#e1e0d9"
# same categorical pair as strict_disagreement_sound_families.png: dGPredictor
# = orange, eQuilibrator = blue, residual = neutral grey (a leftover bucket,
# not a 3rd defect class competing for identity).
DEFECT_OF_CLUSTER = {
    "SAM-dependent methyltransferase": "dGPredictor: activated-cofactor group transfer",
    "PAPS-dependent sulfotransferase": "dGPredictor: activated-cofactor group transfer",
    "ATP/ADP/phosphate-coupled group transfer": "dGPredictor: activated-cofactor group transfer",
    "Glycoside hydrolase (broad)": "eQuilibrator: structural motif",
    "Quinone/quinol-coupled redox": "eQuilibrator: structural motif",
    "Residual (case-specific enzymology)": "none",
}
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
PRETTY = {
    "SAM-dependent methyltransferase": "SAM methyltransferase",
    "PAPS-dependent sulfotransferase": "PAPS sulfotransferase",
    "Glycoside hydrolase (broad)": "Glycoside hydrolase",
    "Quinone/quinol-coupled redox": "Quinone/quinol redox",
    "ATP/ADP/phosphate-coupled group transfer": "ATP/phosphate-coupled\ngroup transfer",
    "Residual (case-specific enzymology)": "Residual\n(case-specific)",
}
PANEL_ORDER = ["all_1370", "ok_ci_only", "ok_only"]
PANEL_TITLE = {
    "all_1370": "all 1,370: equation_defect == none\n(465 LLM-correct)",
    "ok_ci_only": "status OK or CI-only, n=879\n(448 LLM-correct, equation-sound)",
    "ok_only": "status OK only (no flags), n=660\n(438 LLM-correct, equation-sound)",
}
ROW_ORDER = [
    "Glycoside hydrolase (broad)", "Quinone/quinol-coupled redox",
    "SAM-dependent methyltransferase", "PAPS-dependent sulfotransferase",
    "ATP/ADP/phosphate-coupled group transfer", "Residual (case-specific enzymology)",
]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    panels: dict[str, dict[str, int]] = defaultdict(dict)
    totals: dict[str, int] = {}
    with open(DATA) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            panels[r["population"]][r["cluster"]] = int(r["n"])
            totals[r["population"]] = int(r["n_llm_sound_total"])

    fig, axes = plt.subplots(1, 3, figsize=(13.6, 5.4), sharex=True)
    fig.patch.set_facecolor(SURFACE)

    xmax = max(v for pop in panels.values() for v in pop.values())
    y = range(len(ROW_ORDER))
    for ax, pop in zip(axes, PANEL_ORDER):
        ax.set_facecolor(SURFACE)
        counts = [panels[pop].get(c, 0) for c in ROW_ORDER]
        colors = [CLASS_COLOR[DEFECT_OF_CLUSTER[c]] for c in ROW_ORDER]
        ax.barh(list(y), counts, color=colors, zorder=3)
        total = totals[pop]
        for yi, n in zip(y, counts):
            ax.text(n + xmax * 0.02, yi, f"{n} ({100*n/total:.0f}%)", va="center", ha="left",
                    color=INK, fontsize=8.8)
        ax.set_yticks(list(y))
        ax.set_yticklabels([PRETTY[c] for c in ROW_ORDER], color=INK2, fontsize=9.4)
        ax.tick_params(axis="x", colors=INK2, labelsize=9.0)
        ax.grid(True, axis="x", color=GRIDLINE, linewidth=0.8, zorder=0)
        ax.set_axisbelow(True)
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
        ax.spines["bottom"].set_color(GRIDLINE)
        ax.set_xlim(0, xmax * 1.38)
        ax.set_title(PANEL_TITLE[pop], color=INK2, fontsize=9.8, pad=10)
        ax.invert_yaxis()

    for ax in axes[1:]:
        ax.tick_params(axis="y", labelleft=False)

    fig.suptitle("The same two thermo-source defect clusters survive three different\n"
                  "“sound reaction” definitions",
                  x=0.065, y=0.99, ha="left", color=INK, fontsize=14.0, fontweight="semibold")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in CLASS_COLOR.values()]
    fig.legend(handles, [CLASS_LABEL[k] for k in CLASS_COLOR], loc="lower center",
               ncol=3, frameon=False, fontsize=9.6, labelcolor=INK2, bbox_to_anchor=(0.5, 0.0))

    fig.subplots_adjust(left=0.20, right=0.98, top=0.80, bottom=0.14, wspace=0.12)
    out_png = OUT_DIR / "strict_disagreement_balanced_clusters.png"
    fig.savefig(out_png, dpi=200, facecolor=SURFACE)
    plt.close(fig)

    with open(OUT_DIR / "strict_disagreement_balanced_clusters_values.tsv", "w") as fh:
        fh.write("population\tcluster\tn\tn_llm_sound_total\tshare\n")
        for pop in PANEL_ORDER:
            for c in ROW_ORDER:
                n = panels[pop].get(c, 0)
                fh.write(f"{pop}\t{c}\t{n}\t{totals[pop]}\t{round(n/totals[pop], 4)}\n")

    print(f"wrote {out_png}")
    print(f"wrote {OUT_DIR / 'strict_disagreement_balanced_clusters_values.tsv'}")


if __name__ == "__main__":
    main()
