#!/usr/bin/env python3
"""Figures for the direction-influence analysis (reports/DIRECTION_INFLUENCE_V201.md).

    python3 scripts/plot_direction_influence.py
    python3 scripts/plot_direction_influence.py --out reports/figures/direction_influence

Four figures, each answering one question the report asks:

  fig1_growth_by_panel   what each source does to growth, on all three panels
  fig2_grading_vs_thermo how the evidence tiers relate to the three estimators
  fig3_coverage_effect   whether how much a source changes predicts what it does
  fig4_mechanism         which reactions are responsible for the collapse

Colour follows the validated reference palette: slots 1-3 (blue / orange / aqua)
for the three model panels, a single blue ramp where the job is magnitude, and
grey for the baseline. Every mark is direct-labelled, which is also the relief
the aqua slot needs at its contrast ratio.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
from cma import paths  # noqa: E402

# --- reference palette -----------------------------------------------------
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8a85"
GRID, SURFACE = "#e6e6e2", "#fcfcfb"
BLUE_RAMP = LinearSegmentedColormap.from_list("blues", ["#f2f7fd", "#2a78d6", "#17407a"])

PANELS = [("core_all", "Core models\ncomplete media", BLUE),
          ("gsm_gmm", "Genome-scale\nglucose minimal", ORANGE),
          ("gsm_auxo", "Genome-scale\nauxotrophy", AQUA)]

ORDER = ["v201_equilibrator", "v201_dgpredictor", "v201_group_contribution",
         "v201_llm_council", "claude_opus48",
         "v201_gold_only", "v201_silver_only", "v201_bronze_only",
         "v201_gold_silver", "v201_gold_silver_bronze"]
LABEL = {"v201_group_contribution": "Group contribution", "v201_dgpredictor": "dGPredictor",
         "v201_equilibrator": "eQuilibrator", "v201_llm_council": "LLM council",
         "claude_opus48": "Claude Opus 4.8", "v201_gold_only": "gold",
         "v201_silver_only": "silver", "v201_bronze_only": "bronze",
         "v201_gold_silver": "gold + silver",
         "v201_gold_silver_bronze": "gold + silver + bronze"}
FAMILY = {k: ("estimator" if k in ("v201_group_contribution", "v201_dgpredictor",
                                   "v201_equilibrator") else
              "LLM" if k in ("v201_llm_council", "claude_opus48") else "grade")
          for k in ORDER}


def load_results() -> dict:
    """Assemble every input the figures need, straight from results/.

    Depending on a hand-made temp file would make the registry entry
    unreproducible for anyone but its author.
    """
    inf = paths.results("influence_v201")
    sets = paths.results("direction_sets_v201")
    need = {
        "core_all": inf / "core_all", "gsm_gmm": inf / "gsm_gmm",
        "gsm_auxo": inf / "gsm_auxo",
    }
    D = {}
    for tag, d in need.items():
        g, gd = d / "growth.json", d / "growth_diff.json"
        if not g.exists():
            raise SystemExit(
                f"missing {g}.\nRun the sweep first, e.g.\n"
                f"  python3 beginPipeline --models core_kegg2 --directions v201_gold ... "
                f"--only growth,growth_diff --out {d}")
        D[tag] = {"growth": json.loads(g.read_text())["per_variant"],
                  "diff": json.loads(gd.read_text())["per_variant"]}
    for key, path in (("coverage", inf / "coverage.json"),
                      ("manifest", sets / "manifest.json"),
                      ("agreement", sets / "comparison" / "direction_summary.json")):
        if not path.exists():
            raise SystemExit(f"missing {path}")
        D[key] = json.loads(path.read_text())
    return D


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8, length=0)
    ax.grid(axis="x", color=GRID, lw=0.8)
    ax.set_axisbelow(True)


def fig1(D, out: Path):
    """Growth under each source, one panel per model set, against its baseline."""
    fig, axes = plt.subplots(1, 3, figsize=(13, 6.2), sharey=True)
    fig.patch.set_facecolor(SURFACE)
    y = range(len(ORDER))
    for ax, (tag, title, colour) in zip(axes, PANELS, strict=True):
        g = D[tag]["growth"]
        base = g["on_disk"]["frac_growers"]
        vals = [g[k]["frac_growers"] for k in ORDER]
        ax.barh(list(y), vals, height=0.62, color=colour, zorder=3)
        ax.axvline(base, color=MUTED, lw=1.6, ls="--", zorder=4)
        ax.text(base, -0.95, f"baseline {base:.0%}", color=INK2, fontsize=8,
                va="center", ha="center")
        for i, v in enumerate(vals):
            # nudge past the baseline rule when a bar ends on top of it
            off = 0.055 if abs(v - base) < 0.035 else 0.018
            ax.text(v + off, i, f"{v:.0%}", va="center", fontsize=8, color=INK)
        ax.set_title(title, fontsize=10, color=INK, pad=10)
        ax.set_xlim(0, 1.12)
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xticklabels(["0", "25%", "50%", "75%", "100%"])
        style(ax)
    axes[0].set_yticks(list(y))
    axes[0].set_yticklabels([LABEL[k] for k in ORDER], fontsize=9, color=INK)
    axes[0].invert_yaxis()
    fig.suptitle("Share of models that grow under each direction source",
                 fontsize=13, color=INK, x=0.5, y=0.985)
    fig.text(0.5, 0.945, "baseline is the bounds already written in the model files",
             ha="center", fontsize=9, color=INK2)
    fig.tight_layout(rect=(0, 0.02, 1, 0.93))
    fig.savefig(out / "fig1_growth_by_panel.png", dpi=200, facecolor=SURFACE)
    fig.savefig(out / "fig1_growth_by_panel.pdf", facecolor=SURFACE)
    plt.close(fig)


def fig2(D, out: Path):
    """Two facts about the tiers: who they agree with, and how committal they are."""
    ag = D["agreement"]["pairwise_agreement"]
    rows = ["v201_gold_only", "v201_silver_only", "v201_bronze_only",
            "v201_gold_silver", "v201_gold_silver_bronze"]
    cols = ["v201_equilibrator", "v201_dgpredictor", "v201_group_contribution"]

    def cell(a, b):
        e = ag.get(f"{a} vs {b}") or ag.get(f"{b} vs {a}")
        return (e["pct_agree"], e["shared"]) if e and e["pct_agree"] is not None else (None, 0)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5.4),
                                   gridspec_kw={"width_ratios": [1.0, 1.15]})
    fig.patch.set_facecolor(SURFACE)

    M = [[cell(r, c)[0] for c in cols] for r in rows]
    axL.imshow(M, cmap=BLUE_RAMP, vmin=50, vmax=100, aspect="auto")
    axL.set_xticks(range(len(cols)), [LABEL[c] for c in cols], fontsize=9, color=INK)
    axL.set_yticks(range(len(rows)), [LABEL[r] for r in rows], fontsize=9, color=INK)
    for i, r in enumerate(rows):
        for j, c in enumerate(cols):
            pct, shared = cell(r, c)
            # One threshold for both lines, so the sub-label never lands in a
            # mid-tone where neither ink nor white reads.
            dark_cell = pct > 74
            fg = "#ffffff" if dark_cell else INK
            axL.text(j, i - 0.08, f"{pct:.1f}%", ha="center", va="center",
                     fontsize=11, color=fg)
            axL.text(j, i + 0.24, f"n = {shared:,}", ha="center", va="center",
                     fontsize=7.5, color=fg, alpha=0.75)
    axL.set_title("Agreement with each estimator\n"
                  "where both have an opinion on the same reaction",
                  fontsize=10.5, color=INK, pad=12, loc="left")
    for sp in axL.spines.values():
        sp.set_visible(False)
    axL.tick_params(length=0)

    # Right: operator composition -- the tiers are markedly less committal.
    stat = {x["name"]: x for x in D["manifest"]["sets"]}
    bars = ["v201_equilibrator", "v201_dgpredictor", "v201_group_contribution",
            "v201_gold_only", "v201_silver_only", "v201_bronze_only"]
    ypos = list(range(len(bars)))
    fwd = [stat[b]["operators"].get(">", 0) / stat[b]["n_callable"] * 100 for b in bars]
    rev = [stat[b]["operators"].get("<", 0) / stat[b]["n_callable"] * 100 for b in bars]
    rvs = [stat[b]["operators"].get("=", 0) / stat[b]["n_callable"] * 100 for b in bars]
    axR.barh(ypos, fwd, height=0.6, color=BLUE, zorder=3, label="forward")
    axR.barh(ypos, rev, height=0.6, left=[f + 0.6 for f in fwd], color=ORANGE,
             zorder=3, label="reverse")
    axR.barh(ypos, rvs, height=0.6, left=[f + r + 1.2 for f, r in zip(fwd, rev, strict=True)],
             color=AQUA, zorder=3, label="reversible")
    for i, (f, r, v) in enumerate(zip(fwd, rev, rvs, strict=True)):
        axR.text(f / 2, i, f"{f:.0f}", ha="center", va="center", fontsize=8.5, color="#ffffff")
        axR.text(f + r + 1.2 + v / 2, i, f"{v:.0f}", ha="center", va="center",
                 fontsize=8.5, color="#ffffff")
    axR.set_yticks(ypos, [LABEL[b] for b in bars], fontsize=9, color=INK)
    axR.set_xlim(0, 103)
    axR.set_xticks([0, 25, 50, 75, 100], ["0", "25%", "50%", "75%", "100%"])
    axR.set_title("How committal each is\n"
                  "share of its calls that fix a direction rather than allow both",
                  fontsize=10.5, color=INK, pad=12, loc="left")
    axR.legend(frameon=False, fontsize=8.5, ncol=3, loc="lower right",
               bbox_to_anchor=(1.0, -0.19), labelcolor=INK2)
    style(axR)
    # imshow already draws row 0 at the top; only the bar chart needs flipping
    axR.invert_yaxis()
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(out / "fig2_grading_vs_thermo.png", dpi=200, facecolor=SURFACE)
    fig.savefig(out / "fig2_grading_vs_thermo.pdf", facecolor=SURFACE)
    plt.close(fig)


def fig3(D, out: Path):
    """Does how much a source changes predict what it does? Mostly not."""
    cov = D["coverage"]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2))
    fig.patch.set_facecolor(SURFACE)
    for ax, (tag, title, colour), covkey in zip(
            axes, [PANELS[0], PANELS[1]], ["core_kegg2", "ms2_gsm"], strict=True):
        g = D[tag]["growth"]
        pts = []
        for k in ORDER:
            c = cov[covkey]["per_source"].get(k)
            if c:
                pts.append((c["pct_instances_differing"], g[k]["frac_growers"] * 100, LABEL[k]))
        ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=75, color=colour,
                   zorder=3, edgecolor=SURFACE, linewidth=1.5)

        # Label only what can be read. Points piled on top of each other get one
        # shared annotation rather than six overprinted ones, and labels
        # alternate above/below so near-neighbours do not collide.
        flat = [q for q in pts if q[1] < 0.5]
        tall = sorted([q for q in pts if q[1] >= 0.5], key=lambda q: q[0])
        for idx, (x, y_, n) in enumerate(tall):
            dy = 9 if idx % 2 == 0 else -15
            ax.annotate(n, (x, y_), textcoords="offset points", xytext=(0, dy),
                        ha="center", fontsize=8, color=INK2)
        if flat:
            xs = [q[0] for q in flat]
            ax.annotate(f"{len(flat)} sources at 0%   ({min(xs):.0f}-{max(xs):.0f}% changed)",
                        (sum(xs) / len(xs), 0), textcoords="offset points", xytext=(0, -22),
                        ha="center", fontsize=8, color=INK2)

        base = g["on_disk"]["frac_growers"] * 100
        ax.axhline(base, color=MUTED, lw=1.4, ls="--", zorder=2)

        # Each panel gets its own y-range. A shared 0-100 axis squeezes the core
        # panel, whose sources all sit between 51% and 70%, into an unreadable band.
        ally = [q[1] for q in pts] + [base]
        lo, hi = min(ally), max(ally)
        pad = max(6.0, (hi - lo) * 0.28)
        # never show negative growers -- it is not a meaningful value
        ax.set_ylim(max(-0.05 * hi, lo - pad) if lo >= 0 else lo - pad, hi + pad)
        ax.text(ax.get_xlim()[1], base + pad * 0.12, "baseline ", ha="right", va="bottom",
                fontsize=8, color=INK2)
        ax.set_title(title.replace("\n", ", "), fontsize=10.5, color=INK, pad=10, loc="left")
        ax.set_xlabel("% of reaction occurrences the source changes", fontsize=9, color=INK2)
        ax.set_ylabel("% of models that grow", fontsize=9, color=INK2)
        style(ax)
        ax.grid(axis="y", color=GRID, lw=0.8)
    fig.suptitle("How much a source changes does not reliably predict what it does",
                 fontsize=12.5, color=INK, y=0.98)
    fig.text(0.5, 0.925,
             "Spearman rho = -0.03 on the core panel, -0.65 on genome-scale; "
             "eleven nested sets, so a tendency rather than a rule",
             ha="center", fontsize=8.5, color=INK2)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(out / "fig3_coverage_effect.png", dpi=200, facecolor=SURFACE)
    fig.savefig(out / "fig3_coverage_effect.pdf", facecolor=SURFACE)
    plt.close(fig)


LETHAL = [("rxn03108", 98.6, "Thiamine phosphomethylpyrimidine kinase"),
          ("rxn00359", 24.7, "3'-phosphoadenylyl-sulfate sulfohydrolase"),
          ("rxn00763", 16.4, "Glycerol:NAD+ oxidoreductase"),
          ("rxn00379", 12.3, "ATP sulfurylase"),
          ("rxn08801", 12.3, "Lysophospholipase L1"),
          ("rxn05119", 5.5, "L-aspartate:NAD+ oxidoreductase")]


def fig4(_D, out: Path):
    """One reaction carries the genome-scale collapse."""
    fig, ax = plt.subplots(figsize=(11.5, 4.6))
    fig.patch.set_facecolor(SURFACE)
    ys = list(range(len(LETHAL)))
    vals = [v for _, v, _ in LETHAL]
    cols = [ORANGE] + [BLUE] * (len(LETHAL) - 1)
    ax.barh(ys, vals, height=0.6, color=cols, zorder=3)
    for i, v in enumerate(vals):
        ax.text(v + 1.5, i, f"{v:.0f}%", va="center", fontsize=9.5, color=INK)
    ax.set_yticks(ys, [f"{rid}\n{name}" for rid, _, name in LETHAL],
                  fontsize=8.5, color=INK)
    ax.set_xlim(0, 112)
    ax.set_xticks([0, 25, 50, 75, 100], ["0", "25%", "50%", "75%", "100%"])
    ax.set_xlabel("share of tested genome-scale models in which constraining this one "
                  "reaction alone abolishes growth", fontsize=9, color=INK2)
    ax.set_title("A single reaction carries most of the genome-scale collapse",
                 fontsize=12.5, color=INK, pad=12, loc="left")
    ax.invert_yaxis()
    style(ax)
    fig.tight_layout()
    fig.savefig(out / "fig4_mechanism.png", dpi=200, facecolor=SURFACE)
    fig.savefig(out / "fig4_mechanism.pdf", facecolor=SURFACE)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path,
                    default=paths.reports("figures", "direction_influence"))
    ap.add_argument("--data", type=Path, default=None,
                    help="pre-assembled JSON; by default the result files are read directly")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    D = json.loads(args.data.read_text()) if args.data else load_results()
    for fn in (fig1, fig2, fig3, fig4):
        fn(D, args.out)
        print(f"  {fn.__name__} -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
