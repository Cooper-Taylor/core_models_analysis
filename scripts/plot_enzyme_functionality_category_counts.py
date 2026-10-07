#!/usr/bin/env python3
"""Bar chart: ModelSEED reactions (of the 20,490-reaction LLM-vs-recommended
comparable set) per enzyme-functionality category.

Reads results/enzyme_functionality/category_counts.tsv (written by
build_enzyme_functionality_outputs.py). Re-run after every literature-search
batch -- monofunctional/bifunctional/polyfunctional grow and incomplete_ec
stays fixed at 11,853 (reactions with no complete EC annotation, never
searched).

DESIGN NOTES
  * Categorical slots 1/2/3 (blue/orange/green) are this repo's validated
    palette (see plot_graded_fba.py) for the three literature-derived
    outcomes. incomplete_ec and pending are "not applicable yet" statuses,
    not classification outcomes, so they get neutral grays instead of
    competing with A/B/C for categorical identity.
  * category_counts.tsv carries a "pending" row (reactions with >=1 complete
    EC where none of them are classified yet) specifically so this figure's
    total/percentage is always the true 20,490-reaction denominator, not
    just whatever has been classified so far.
  * Direct value labels on every bar; no legend needed for 5 labelled bars.

Data: ModelSEED dev @ 078a395f (/scratch/ctaylor/tmp/devsnap_078a395f).
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
DATA = ANALYSIS_DIR / "results" / "enzyme_functionality" / "category_counts.tsv"
OUT_DIR = ANALYSIS_DIR / "reports" / "enzymeFunctionality" / "figures"

INK, INK2 = "#0b0b0b", "#52514e"
SURFACE, GRIDLINE = "#fcfcfb", "#e1e0d9"
NO_DATA_GRAY = "#e1e0d9"
PENDING_GRAY = "#c3c2b7"
# categorical slots 1,2,3 (validated light, this repo's standard set)
CAT_COLOR = {"monofunctional": "#2a78d6", "bifunctional": "#eb6834",
             "polyfunctional": "#1baf7a", "incomplete_ec": NO_DATA_GRAY,
             "pending": PENDING_GRAY}
ORDER = ["monofunctional", "bifunctional", "polyfunctional", "incomplete_ec", "pending"]
PRETTY = {"monofunctional": "monofunctional\n(A)", "bifunctional": "bifunctional\n(B)",
          "polyfunctional": "polyfunctional\n(C)", "incomplete_ec": "incomplete EC\n(not searched)",
          "pending": "pending\n(later batches)"}


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    counts = {}
    with open(DATA) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            counts[row["category"]] = int(row["n_reactions"])

    total = sum(counts.values())
    values = [counts.get(c, 0) for c in ORDER]
    colors = [CAT_COLOR[c] for c in ORDER]

    fig, ax = plt.subplots(figsize=(8.0, 5.6))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    x = range(len(ORDER))
    ax.bar(x, values, width=0.62, color=colors, edgecolor="none", zorder=2)
    ax.grid(True, axis="y", color=GRIDLINE, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color(GRIDLINE)

    ax.set_xticks(list(x))
    ax.set_xticklabels([PRETTY[c] for c in ORDER], color=INK2, fontsize=10.2)
    ax.tick_params(axis="y", colors=INK2, labelsize=9.6)
    ax.set_ylabel("ModelSEED reactions", color=INK2, fontsize=10.6)

    ymax = max(values) if max(values) > 0 else 1
    for xi, v in zip(x, values):
        ax.text(xi, v + ymax * 0.015, f"{v:,}", ha="center", va="bottom",
                color=INK, fontsize=10.8, fontweight="semibold")

    n_resolved = total - counts.get("pending", 0) - counts.get("incomplete_ec", 0)
    n_searchable = total - counts.get("incomplete_ec", 0)
    pct_done = 100.0 * n_resolved / n_searchable if n_searchable else 0.0
    fig.suptitle("Enzyme functionality of the LLM-vs-recommended reaction set",
                 x=0.11, y=0.975, ha="left", color=INK, fontsize=13.6, fontweight="semibold")
    ax.set_title(f"{total:,} reactions  ·  {n_searchable:,} literature-searchable (>=1 complete EC)  ·  "
                 f"{pct_done:.1f}% resolved so far",
                 loc="left", color=INK2, fontsize=9.6, fontweight="normal", pad=12)

    fig.subplots_adjust(left=0.11, right=0.97, top=0.82, bottom=0.14)
    out_png = OUT_DIR / "enzyme_functionality_category_counts.png"
    fig.savefig(out_png, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
