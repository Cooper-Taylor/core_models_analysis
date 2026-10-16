#!/usr/bin/env python3
"""Inside the "Residual (case-specific enzymology)" bar of
strict_disagreement_balanced_clusters.png: ten smaller, almost perfectly
thermo-source-pure sub-clusters, covering 78 of the 248 residual reactions
(31.5%).

Why this figure is ONE panel and not three. The sibling figure
(strict_disagreement_balanced_clusters.png) draws the same bars once per
population -- all_1370 / ok_ci_only / ok_only -- and reads the stability of
the numbers as corroboration. Those three residuals are in fact STRICTLY
NESTED (233 subset-of 242 subset-of 248, proven in
analyze_strict_disagreement_residual_subclusters.py: every population is a
row-subset of the one before it, and every filter downstream is a pure
per-reaction predicate), so three panels would be three near-copies of the
same reactions. Drawing it once and marking the difference is the honest
encoding:

    solid segment  = reactions present in ALL THREE populations (the ok_only
                     count) -- i.e. ModelSEED's own status flag is clean
    pale segment   = the few reactions the stricter balance filters remove
                     (present in all_1370 but not in ok_only)

Only two of the ten sub-clusters lose anything at all, one reaction each
(cis-dihydrodiol dehydrogenase 7->6, dehalogenase 5->4); the "still
unexplained" remainder loses 13. That is what "the residual is the same set
under every balance filter" looks like drawn out.

Layout: the still-unexplained remainder (170) is ~9x the largest sub-cluster
(19), so putting it on the same axis compresses the ten bars the figure is
about into the left tenth of the panel. It gets its own full-width
composition strip at the top instead -- which is also the more honest place
for it, since it is a leftover, not an eleventh cluster -- and the main axis
is scaled to the sub-clusters alone. Both scales start at zero; neither is
broken.

Colour follows the sibling figures' categorical pair exactly -- dGPredictor
orange #eb6834, eQuilibrator blue #2a78d6 -- so a reader moving between the
two figures keeps the same source->hue mapping. Validated with the dataviz
skill's validate_palette.js: the two categorical hues PASS all six checks
(worst adjacent CVD dE 24.7 protan, 32.7 tritan; normal-vision dE 33.6). The
third colour, grey #c3c2b7, deliberately FAILS the categorical chroma/
lightness checks because it is NOT a competing category -- it marks the
leftover bucket, the same way the sibling figure greys its residual. Its
sub-3:1 contrast is relieved as the validator requires, by a direct count
label on every bar plus the values TSV written beside the PNG.

Sub-cluster definitions, membership and the nesting proof:
scripts/analyze_strict_disagreement_residual_subclusters.py ->
results/strict_disagreement_enzymes/residual_subclusters.tsv (this figure's
input) and residual_subclusters.json. Same ModelSEED dev @ 078a395f
strict-disagreement literature pass (2026-10-08) as every other figure in
this report -- no new research, no ModelSEED biochemistry read, just a finer
cut of the already-collected deliverable.
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
DATA = ANALYSIS_DIR / "results" / "strict_disagreement_enzymes" / "residual_subclusters.tsv"
OUT_DIR = ANALYSIS_DIR / "reports" / "strictDisagreementEnzymes" / "figures"

INK, INK2, INK3 = "#0b0b0b", "#52514e", "#898781"
# Thin marks: leave surface between rows so the two-line y labels never touch.
BAR_H = 0.58
SURFACE, GRIDLINE = "#fcfcfb", "#e1e0d9"

# Same categorical pair as strict_disagreement_sound_families.png and
# strict_disagreement_balanced_clusters.png. Grey = leftover, not a category.
CLASS_COLOR = {
    "dGPredictor": "#eb6834",
    "eQuilibrator": "#2a78d6",
    "none": "#c3c2b7",
}
# Pale tint of each, for the reactions a stricter balance filter removes.
CLASS_TINT = {
    "dGPredictor": "#f8c7b2",
    "eQuilibrator": "#b5d0f0",
    "none": "#e3e2dc",
}
CLASS_LABEL = {
    "dGPredictor": "dGPredictor-pure sub-cluster",
    "eQuilibrator": "eQuilibrator-pure sub-cluster",
    "none": "still unexplained",
}

# Grouped by defect class, largest first inside each -- so the two
# consecutive steps of the chlorocatechol pathway (cycloisomerase, dGP;
# dienelactone hydrolase, eQ) stay visible as the split pair they are.
ROW_ORDER = [
    "Ring aminohydrolase / hydrolytic deaminase",
    "Triple-bond hydratase (alkyne / nitrile)",
    "Chloromuconate cycloisomerase",
    "cis-Dihydrodiol dehydrogenase",
    "Epoxide hydrolase",
    "Hydrolytic / reductive dehalogenase",
    "Flavonol glycoside acyl / glucosyl transfer",
    "Dienelactone hydrolase",
    "HCCA isomerase (naphthalene, NahD)",
    "Acyl-CoA dehydrogenase (beta-oxidation)",
]
LEFTOVER_ROW = "Still unexplained"
PRETTY = {
    "Ring aminohydrolase / hydrolytic deaminase": "Ring aminohydrolase /\nhydrolytic deaminase",
    "Triple-bond hydratase (alkyne / nitrile)": "Triple-bond hydratase\n(alkyne / nitrile)",
    "Chloromuconate cycloisomerase": "Chloromuconate\ncycloisomerase",
    "cis-Dihydrodiol dehydrogenase": "cis-Dihydrodiol\ndehydrogenase",
    "Epoxide hydrolase": "Epoxide hydrolase",
    "Hydrolytic / reductive dehalogenase": "Hydrolytic / reductive\ndehalogenase",
    "Flavonol glycoside acyl / glucosyl transfer": "Flavonol glycoside\nacyl / glucosyl transfer",
    "Dienelactone hydrolase": "Dienelactone hydrolase",
    "HCCA isomerase (naphthalene, NahD)": "HCCA isomerase\n(naphthalene, NahD)",
    "Acyl-CoA dehydrogenase (beta-oxidation)": "Acyl-CoA dehydrogenase\n(β-oxidation)",
}


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows: dict[str, dict[str, dict]] = {}
    totals: dict[str, int] = {}
    with open(DATA) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.setdefault(r["population"], {})[r["cluster"]] = r
            totals[r["population"]] = int(r["n_residual_total"])

    all_p, ok_p = rows["all_1370"], rows["ok_only"]

    fig, (ax_top, ax) = plt.subplots(
        2, 1, figsize=(11.4, 6.9),
        gridspec_kw={"height_ratios": [1, 8.2], "hspace": 0.30})
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    ax_top.set_facecolor(SURFACE)

    # ---- composition strip: what the 248 residual reactions are made of ----
    n_resid = totals["all_1370"]
    n_dgp = sum(int(all_p[c]["n"]) for c in ROW_ORDER if all_p[c]["defect_class"] == "dGPredictor")
    n_eq = sum(int(all_p[c]["n"]) for c in ROW_ORDER if all_p[c]["defect_class"] == "eQuilibrator")
    n_left = int(all_p[LEFTOVER_ROW]["n"])
    left = 0.0
    for n, klass in ((n_dgp, "dGPredictor"), (n_eq, "eQuilibrator"), (n_left, "none")):
        ax_top.barh(0, n, left=left, color=CLASS_COLOR[klass], zorder=3,
                    edgecolor=SURFACE, linewidth=2.0)
        ax_top.text(left + n / 2, 0, str(n), va="center", ha="center",
                    color="#ffffff" if klass != "none" else INK, fontsize=9.4,
                    fontweight="semibold", zorder=4)
        left += n
    ax_top.set_xlim(0, n_resid)
    ax_top.set_ylim(-0.6, 0.6)
    ax_top.axis("off")
    ax_top.set_title(f"all {n_resid} residual reactions  —  "
                     f"{n_dgp + n_eq} sub-clustered ({(n_dgp + n_eq) / n_resid:.0%}), "
                     f"{n_left} with no further pattern",
                     color=INK2, fontsize=9.6, loc="left", pad=6)

    # ---- main panel: the ten sub-clusters, on their own scale -------------
    xmax = max(int(all_p[c]["n"]) for c in ROW_ORDER)
    y = list(range(len(ROW_ORDER)))

    for yi, cluster in zip(y, ROW_ORDER):
        n_all = int(all_p[cluster]["n"])
        n_ok = int(ok_p[cluster]["n"])
        klass = all_p[cluster]["defect_class"]
        # solid = present in all three populations; pale = shaved by a
        # stricter balance filter. 2px surface gap between the two fills.
        ax.barh(yi, n_ok, height=BAR_H, color=CLASS_COLOR[klass], zorder=3)
        if n_all > n_ok:
            ax.barh(yi, n_all - n_ok, left=n_ok, height=BAR_H, color=CLASS_TINT[klass],
                    zorder=3, edgecolor=SURFACE, linewidth=2.0)

        pur = f"{100*float(all_p[cluster]['dominant_source_purity']):.0f}% " \
              f"{all_p[cluster]['dominant_source']}"
        label = f"{n_all}  ({pur})" if n_all == n_ok else f"{n_all}→{n_ok}  ({pur})"
        ax.text(n_all + xmax * 0.04, yi, label, va="center", ha="left",
                color=INK, fontsize=9.2)

    ax.set_yticks(y)
    ax.set_yticklabels([PRETTY[c] for c in ROW_ORDER], color=INK2, fontsize=9.2)
    ax.tick_params(axis="x", colors=INK2, labelsize=9.0)
    ax.set_xlabel("reactions", color=INK2, fontsize=9.6)
    ax.grid(True, axis="x", color=GRIDLINE, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color(GRIDLINE)
    ax.set_xlim(0, xmax * 1.44)
    ax.invert_yaxis()

    fig.suptitle("The “case-specific” residual is not structureless: ten small,\n"
                 "source-pure sub-clusters cover 78 of its 248 reactions",
                 x=0.012, y=0.985, ha="left", color=INK, fontsize=14.2,
                 fontweight="semibold")
    fig.text(0.012, 0.862,
             "Residual of the 465 LLM-correct, equation-sound strict disagreements, sub-clustered on enzyme name by thermo-source purity.\n"
             "The three populations of strict_disagreement_balanced_clusters.png are strictly nested (233 ⊆ 242 ⊆ 248), so the bars are drawn\n"
             "once: solid = present in all three, pale = removed by a stricter balance filter.",
             ha="left", va="top", color=INK3, fontsize=9.0)

    handles = [plt.Rectangle((0, 0), 1, 1, color=CLASS_COLOR[k]) for k in CLASS_COLOR]
    handles.append(plt.Rectangle((0, 0), 1, 1, color=CLASS_TINT["none"]))
    labels = [CLASS_LABEL[k] for k in CLASS_COLOR] + ["removed by a stricter balance filter"]
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False,
               fontsize=9.3, labelcolor=INK2, bbox_to_anchor=(0.54, 0.0))

    fig.subplots_adjust(left=0.215, right=0.985, top=0.745, bottom=0.125)
    out_png = OUT_DIR / "strict_disagreement_residual_subclusters.png"
    fig.savefig(out_png, dpi=200, facecolor=SURFACE)
    plt.close(fig)

    out_tsv = OUT_DIR / "strict_disagreement_residual_subclusters_values.tsv"
    with open(out_tsv, "w") as fh:
        fh.write("cluster\tdefect_class\tn_all_1370\tn_ok_ci_only\tn_ok_only\t"
                 "dominant_source\tdominant_source_purity\tn_residual_all_1370\n")
        for c in ROW_ORDER + [LEFTOVER_ROW]:
            fh.write(f"{c}\t{all_p[c]['defect_class']}\t{all_p[c]['n']}\t"
                     f"{rows['ok_ci_only'][c]['n']}\t{ok_p[c]['n']}\t"
                     f"{all_p[c]['dominant_source']}\t{all_p[c]['dominant_source_purity']}\t"
                     f"{totals['all_1370']}\n")

    print(f"wrote {out_png}")
    print(f"wrote {out_tsv}")


if __name__ == "__main__":
    main()
