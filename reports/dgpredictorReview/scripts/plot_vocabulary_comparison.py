#!/usr/bin/env python3
"""Three panels on what the fine-tuned dGPredictor's fingerprint basis actually is,
next to the two hand-curated group vocabularies it is usually compared with.

Panel A  declared vocabulary vs the part of it the data can constrain, per method.
         Drawn as a dumbbell on a log position axis (NOT bars: bar length is only
         honest from zero, and these span four orders of magnitude).
Panel B  every one of the 1,415 learned fingerprints: how many training compounds
         support it, against the magnitude of the energy it was assigned.
Panel C  for all predictable ModelSEED reactions, the fraction of the reaction's
         group-change vector that lies outside the span of the training set.

Palette: the two categorical hues already validated for this repo's thermo
figures (validate_palette.js, light mode, surface #fcfcfb, --pairs all: PASS),
plus the neutral gray used for "not a category" endpoints.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

BLUE = "#2a78d6"
ORANGE = "#eb6834"
NEUTRAL = "#9c9a94"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
BASELINE = "#c3c2b7"
SURFACE = "#fcfcfb"


def frame(ax, *, grid_axis: str = "both") -> None:
    ax.grid(True, axis=grid_axis, color=GRIDLINE, linewidth=0.8, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(BASELINE)
    ax.tick_params(colors=INK_MUTED, labelsize=9)


def main() -> None:
    voc = json.loads((RES / "vocabulary_comparison.json").read_text())
    ident = json.loads((RES / "identifiability_and_confidence.json").read_text())
    rowsp = json.loads((RES / "rowspace_coverage.json").read_text())
    base = json.loads((RES / "base_fit_summary.json").read_text())
    shares = np.load(RES / "out_of_span_share.npy")

    ft = voc["finetuned_dgpredictor"]
    eqc = voc["equilibrator_current"]
    eq2 = voc["equilibrator_in_dgpredictor_repo"]
    gc = voc["modelseed_group_contribution"]

    # (label, declared vocabulary, data-constrained dimensions, note)
    ROWS = [
        ("dGPredictor\nModelSEED fine-tuned",
         ft["declared_vocabulary"]["total"],
         ident["identifiability"]["finetuned_dgpredictor"]["rank_of_reaction_feature_matrix"],
         f"{ft['learned_vocabulary']['total']:,} non-zero coefficients"),
        ("dGPredictor\noriginal (KEGG)",
         base["n1"] + base["n2"], base["rank_of_reaction_feature_matrix"],
         f"{base['used_total']:,} non-zero coefficients"),
        ("eQuilibrator\ncomponent-contribution 0.7",
         eqc["n_real_groups"] + eqc["n_placeholder_columns"],
         eqc["reaction_feature_matrix_rank"],
         f"{eqc['n_columns_used_in_training']:,} non-zero coefficients"),
        ("eQuilibrator\nCC 2.x (in dGPredictor repo)",
         eq2["n_columns"], eq2["reaction_feature_matrix_rank"],
         f"{eq2['n_columns_used_in_training']:,} non-zero coefficients"),
        ("ModelSEED\nGroup Contribution",
         gc["vocabulary_size"], None,
         "hand-curated named groups; training set not distributed"),
    ]

    fig = plt.figure(figsize=(17.0, 5.8), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.5, 1.0, 1.0], wspace=0.34,
                          left=0.115, right=0.985, top=0.80, bottom=0.21)

    # ---------------- Panel A ------------------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    ax.set_facecolor(SURFACE)
    ys = np.arange(len(ROWS))[::-1]
    ax.set_xscale("log")
    ax.set_xlim(45, 6e5)
    for y, (label, declared, rank, note) in zip(ys, ROWS):
        if rank is not None:
            ax.plot([rank, declared], [y, y], color=BASELINE, linewidth=2.0, zorder=1)
            ax.scatter([rank], [y], s=90, color=BLUE, zorder=3, linewidths=0)
            # value labels: rank BELOW the dot, declared ABOVE it, so the two
            # never collide even when the dumbbell is short (eQuilibrator rows)
            ax.annotate(f"{rank:,}", (rank, y), xytext=(0, -16),
                        textcoords="offset points", ha="center", fontsize=8.5,
                        color=BLUE)
        ax.scatter([declared], [y], s=90, color=NEUTRAL, zorder=3, linewidths=0)
        ax.annotate(f"{declared:,}", (declared, y), xytext=(0, 11),
                    textcoords="offset points", ha="center", fontsize=8.5,
                    color=INK_SECONDARY)
        if note:
            ax.annotate(note, (declared, y), xytext=(10, -3),
                        textcoords="offset points", ha="left", va="center",
                        fontsize=7.5, color=INK_MUTED)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in ROWS], fontsize=9, color=INK_SECONDARY)
    ax.set_ylim(-0.55, len(ROWS) - 0.15)
    ax.set_xlabel("number of groups / independent directions (log)",
                  color=INK_SECONDARY, fontsize=10)
    ax.set_title("A · Alphabet size vs what the data can constrain",
                 color=INK_PRIMARY, fontsize=11.5, loc="left", pad=22)
    frame(ax, grid_axis="x")
    ax.scatter([], [], s=90, color=NEUTRAL, label="declared vocabulary")
    ax.scatter([], [], s=90, color=BLUE,
               label="rank of the training matrix (directions the data fixes)")
    lg = ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.005), ncol=2,
                   frameon=False, fontsize=8.5, labelcolor=INK_PRIMARY,
                   borderpad=0.0, columnspacing=1.4, handletextpad=0.4)
    for t in lg.get_texts():
        t.set_linespacing(1.3)

    # ---------------- Panel B ------------------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    ax.set_facecolor(SURFACE)
    rows = [ln.split("\t") for ln in
            (RES / "learned_fingerprints.tsv").read_text().strip().split("\n")[1:]]
    supp = np.array([int(r[4]) for r in rows], float)
    mag = np.abs(np.array([float(r[2]) for r in rows]))
    ax.scatter(np.clip(supp, 0.7, None), np.clip(mag, 1e-3, None),
               s=12, linewidths=0, color=ORANGE, alpha=0.45, zorder=2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    n_singleton = int((supp <= 1).sum())
    ax.axvline(1.5, color=BASELINE, linewidth=1.2, linestyle="--", zorder=1)
    ax.annotate(f"{n_singleton} fingerprints seen in\nexactly one training compound",
                (1.6, ax.get_ylim()[1]), xytext=(4, -6), textcoords="offset points",
                fontsize=8, color=INK_MUTED, va="top", linespacing=1.4)
    ax.set_xlabel("training compounds containing the fingerprint (log)",
                  color=INK_SECONDARY, fontsize=10)
    ax.set_ylabel("|assigned energy| (kJ/mol, log)", color=INK_SECONDARY, fontsize=10)
    ax.set_title(f"B · The {len(rows):,} fingerprints with a non-zero energy",
                 color=INK_PRIMARY, fontsize=11.5, loc="left", pad=22)
    frame(ax)

    # ---------------- Panel C ------------------------------------------------
    ax = fig.add_subplot(gs[0, 2])
    ax.set_facecolor(SURFACE)
    ax.hist(shares, bins=40, color=ORANGE, alpha=0.85, zorder=2, edgecolor=SURFACE,
            linewidth=0.5)
    med = float(np.median(shares))
    ax.axvline(med, color=BLUE, linewidth=2.0, zorder=3)
    ax.annotate(f"median {med:.2f}", (med, ax.get_ylim()[1]), xytext=(6, -8),
                textcoords="offset points", fontsize=9, color=BLUE, va="top")
    ax.set_xlabel("fraction of the reaction's group-change vector\n"
                  "outside the training span", color=INK_SECONDARY, fontsize=10)
    ax.set_ylabel("reactions", color=INK_SECONDARY, fontsize=10)
    ax.set_title(f"C · Extrapolation ({len(shares):,} predicted reactions)",
                 color=INK_PRIMARY, fontsize=11.5, loc="left", pad=22)
    frame(ax)
    ax.annotate(f"{100 * float((shares > 0.5).mean()):.0f}% of reactions are\n"
                f"more than half outside",
                (0.97, 0.72), xycoords="axes fraction", ha="right", fontsize=8.5,
                color=INK_MUTED, linespacing=1.4)

    fig.suptitle("What the ModelSEED-fine-tuned dGPredictor actually learned, "
                 "and how its alphabet compares with the curated group vocabularies",
                 color=INK_PRIMARY, fontsize=13, y=0.975, x=0.008, ha="left")
    out = FIG / "vocabulary_comparison.png"
    fig.savefig(out, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")
    print(f"  rowspace: {json.dumps(rowsp['out_of_span_norm_fraction'])}")


if __name__ == "__main__":
    main()
