#!/usr/bin/env python3
"""Figures for the three-source DeltaG agreement analysis.

  fig_dg_agreement_by_kegg_trust.png
      THE EVIDENCE FOR THE MASK, and the only figure here drawn from the
      UNFILTERED table (reaction_features_unmasked.tsv): Group Contribution vs
      dGPredictor, split by whether the KEGG reaction id dGPredictor was
      actually run on is one ModelSEED itself lists as an alias of that
      reaction. The near-zero whole-database correlation is carried almost
      entirely by the inferred-mapping half, which is what
      build_dgpredictor_kegg_mask.py then withholds.

  fig_dg_agreement_by_reaction_class.png
      Post-filter. Same axes over the reactions that survive the mask, colored
      by reaction class -- redox chemistry vs group transfer -- which is where
      the remaining, genuinely chemical, disagreement lives.

  fig_agreement_vs_snr.png
      Post-filter. Pairwise correlation by propagated signal-to-noise decile,
      for all three source pairs.

Palette: the 3-hue categorical triple already in use in
plot_thermo_source_dg_scatter.py; re-validated for this figure set with
validate_palette.js (light, --pairs all): all checks PASS, with a contrast WARN
on the aqua slot that is discharged by the always-present legend and the
in-panel text labels.
"""
from __future__ import annotations

import glob
import json
import os
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

MSDB_ROOT = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase"))
ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
DATA_DIR = ANALYSIS_DIR / "results" / "thermo_agreement"
OUT_DIR = ANALYSIS_DIR / "reports" / "thermoComparison" / "figures" / "thermo_agreement"

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
NEUTRAL = "#9c9a94"
INK_PRIMARY, INK_SECONDARY, INK_MUTED = "#0b0b0b", "#52514e", "#898781"
GRIDLINE, BASELINE, SURFACE = "#e1e0d9", "#c3c2b7", "#fcfcfb"

# Chemically ordinary single-step reaction DeltaG'deg sits well inside this
# window; a few aggregate/polymer reactions do not, and are reported in the
# caption rather than silently dropped (same convention as the existing
# thermo_source_dg_scatter figures).
WINDOW = 250.0


def style(ax) -> None:
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRIDLINE, linewidth=0.8, zorder=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(BASELINE)
    ax.tick_params(colors=INK_MUTED, labelsize=9)


def identity_line(ax, lo: float, hi: float) -> None:
    ax.plot([lo, hi], [lo, hi], linestyle="--", linewidth=1.2,
            color=BASELINE, zorder=1)


def fig_kegg_trust(t: pd.DataFrame) -> Path:
    panels = [
        ("KEGG id is a ModelSEED alias\nof this reaction", t[t["kegg_vouched"] == 1], BLUE),
        ("KEGG id inferred\n(reaction lists no KEGG alias)", t[t["kegg_vouched"] == 0], ORANGE),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.6), dpi=150, sharex=True, sharey=True)
    fig.patch.set_facecolor(SURFACE)
    for ax, (label, sub, color) in zip(axes, panels):
        x, y = sub["dg_gc"].to_numpy(float), sub["dg_dgp"].to_numpy(float)
        r_all = np.corrcoef(x, y)[0, 1]
        rho = stats.spearmanr(x, y).statistic
        keep = (np.abs(x) <= WINDOW) & (np.abs(y) <= WINDOW)
        style(ax)
        identity_line(ax, -WINDOW, WINDOW)
        ax.scatter(x[keep], y[keep], s=10, linewidths=0, color=color,
                   alpha=0.45, zorder=2)
        ax.set_xlim(-WINDOW, WINDOW)
        ax.set_ylim(-WINDOW, WINDOW)
        ax.set_title(label, color=INK_PRIMARY, fontsize=11, pad=10)
        ax.set_xlabel("Group Contribution ΔG′° (kcal/mol)",
                      color=INK_SECONDARY, fontsize=10)
        ax.text(0.04, 0.96,
                f"n = {len(sub):,}\nPearson r = {r_all:+.3f}\nSpearman ρ = {rho:+.3f}\n"
                f"median |Δ| = {np.median(np.abs(x - y)):.1f} kcal/mol",
                transform=ax.transAxes, va="top", ha="left", fontsize=9,
                color=INK_SECONDARY, linespacing=1.6)
        ax.text(0.97, 0.04, f"{(~keep).sum()} point(s) outside ±{WINDOW:.0f}",
                transform=ax.transAxes, va="bottom", ha="right",
                fontsize=7.5, color=INK_MUTED)
    axes[0].set_ylabel("dGPredictor ΔG′° (kcal/mol)", color=INK_SECONDARY, fontsize=10)
    fig.suptitle("dGPredictor tracks Group Contribution only where the KEGG reaction "
                 "mapping is vouched for",
                 color=INK_PRIMARY, fontsize=13, y=0.99)
    fig.text(0.5, 0.005,
             "dGPredictor predicts from a KEGG reaction, so the stored ModelSEED value is only about "
             "the ModelSEED reaction\nwhen that KEGG mapping is correct. Same axes, same method, both panels.",
             ha="center", va="bottom", fontsize=8, color=INK_MUTED, linespacing=1.5)
    fig.tight_layout(rect=(0, 0.05, 1, 0.94))
    out = OUT_DIR / "fig_dg_agreement_by_kegg_trust.png"
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)
    return out


def fig_reaction_class(t: pd.DataFrame) -> Path:
    v = t.copy()
    redox = (v[["cof_nad", "cof_nadh", "cof_nadp", "cof_nadph"]].max(axis=1) == 1) | \
            v["ec_class"].fillna("none").str.contains("1")
    transfer = (v[["cof_atp", "cof_adp", "cof_amp", "cof_sam", "cof_sah"]].max(axis=1) == 1) | \
               (v["d_phosphoanhydride"] != 0) | (v["d_thioester"] != 0)
    v["klass"] = "Other"
    v.loc[transfer & ~redox, "klass"] = "Group transfer (phosphoryl / methyl / acyl)"
    v.loc[redox & ~transfer, "klass"] = "Redox (oxidoreductase, NAD(P)(H))"

    order = ["Other",
             "Redox (oxidoreductase, NAD(P)(H))",
             "Group transfer (phosphoryl / methyl / acyl)"]
    colors = {order[0]: NEUTRAL, order[1]: AQUA, order[2]: ORANGE}

    fig, ax = plt.subplots(figsize=(7.6, 6.8), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    style(ax)
    identity_line(ax, -WINDOW, WINDOW)
    lines = []
    for name in order:
        sub = v[v["klass"] == name]
        x, y = sub["dg_gc"].to_numpy(float), sub["dg_dgp"].to_numpy(float)
        keep = (np.abs(x) <= WINDOW) & (np.abs(y) <= WINDOW)
        ax.scatter(x[keep], y[keep], s=11, linewidths=0, color=colors[name],
                   alpha=0.55 if name != "Other" else 0.3,
                   label=f"{name}  (n = {len(sub):,})", zorder=2 if name == "Other" else 3)
        lines.append(f"{name.split(' (')[0]}:  r = {np.corrcoef(x, y)[0, 1]:+.2f}, "
                     f"median |Δ| = {np.median(np.abs(x - y)):.1f} kcal/mol")
    ax.set_xlim(-WINDOW, WINDOW)
    ax.set_ylim(-WINDOW, WINDOW)
    ax.set_xlabel("Group Contribution ΔG′° (kcal/mol)", color=INK_SECONDARY, fontsize=10)
    ax.set_ylabel("dGPredictor ΔG′° (kcal/mol)", color=INK_SECONDARY, fontsize=10)
    ax.set_title("After removing the mis-mapped values, redox reactions agree\n"
                 "and group-transfer reactions do not",
                 color=INK_PRIMARY, fontsize=12, pad=12)
    legend = ax.legend(loc="upper left", frameon=False, fontsize=8.5,
                       labelcolor=INK_PRIMARY, markerscale=1.8)
    legend.set_zorder(5)
    ax.text(0.98, 0.03, "\n".join(lines), transform=ax.transAxes,
            va="bottom", ha="right", fontsize=8, color=INK_SECONDARY, linespacing=1.7)
    fig.tight_layout()
    out = OUT_DIR / "fig_dg_agreement_by_reaction_class.png"
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)
    return out


def fig_snr(t: pd.DataFrame) -> Path:
    cpd_dg, cpd_err = {}, {}
    for path in sorted(glob.glob(str(MSDB_ROOT / "Biochemistry" / "compound_*.json"))):
        for entry in json.load(open(path)):
            val = entry.get("deltag")
            if isinstance(val, (int, float)) and abs(val) < 1e6:
                cpd_dg[entry["id"]] = float(val)
                err = entry.get("deltagerr")
                cpd_err[entry["id"]] = float(err) if isinstance(err, (int, float)) and abs(err) < 1e6 else 0.0
    stoich = {}
    for path in sorted(glob.glob(str(MSDB_ROOT / "Biochemistry" / "reaction_*.json"))):
        for entry in json.load(open(path)):
            vec = defaultdict(float)
            for item in entry.get("stoichiometry") or []:
                c = float(item.get("coefficient", 0) or 0)
                if c:
                    vec[item["compound"]] += c
            stoich[entry["id"]] = {k: val for k, val in vec.items() if val != 0}

    v = t.copy()
    snr = []
    for rid in v["rxn"]:
        vec = stoich.get(rid, {})
        if not vec or any(c not in cpd_dg for c in vec):
            snr.append(np.nan); continue
        net = sum(c * cpd_dg[k] for k, c in vec.items())
        sig = float(np.sqrt(sum((c * cpd_err[k]) ** 2 for k, c in vec.items())))
        snr.append(abs(net) / sig if sig > 0 else np.nan)
    v["snr"] = snr
    v = v[np.isfinite(v["snr"])].copy()
    v["bin"] = pd.qcut(v["snr"], 10, labels=False, duplicates="drop")

    pairs = [("gc", "eq", "Group Contribution vs eQuilibrator", BLUE),
             ("gc", "dgp", "Group Contribution vs dGPredictor", ORANGE),
             ("eq", "dgp", "eQuilibrator vs dGPredictor", AQUA)]
    fig, ax = plt.subplots(figsize=(8.2, 5.4), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    style(ax)
    xs = np.arange(1, v["bin"].nunique() + 1)
    for a, b, label, color in pairs:
        ys = [np.corrcoef(g[f"dg_{a}"], g[f"dg_{b}"])[0, 1]
              for _, g in v.groupby("bin")]
        # Three series that converge and cross at both ends: end-of-line direct
        # labels would overlap each other, so identity is carried by the legend
        # alone (still not color-alone -- the legend is always present).
        ax.plot(xs, ys, linewidth=2, color=color, marker="o", markersize=6,
                markeredgecolor=SURFACE, markeredgewidth=1.5, label=label, zorder=3)
    ax.set_xticks(xs)
    ax.set_xlabel("propagated signal-to-noise decile   |ΔG′°| / σ from compound-energy errors "
                  "(1 = lowest)", color=INK_SECONDARY, fontsize=10)
    ax.set_ylabel("Pearson r between the two sources", color=INK_SECONDARY, fontsize=10)
    ax.set_title("All three sources converge on reactions whose ΔG is large relative to\n"
                 "the error already in their compound energies",
                 color=INK_PRIMARY, fontsize=12, pad=12)
    ax.legend(loc="lower right", frameon=False, fontsize=8.5, labelcolor=INK_PRIMARY)
    ax.set_xlim(0.6, len(xs) + 0.4)
    fig.text(0.5, 0.01,
             f"Reactions surviving the dGPredictor KEGG mask (n = {len(v):,}). Decile 10 falls back "
             "because it is dominated\nby a small number of very large-|ΔG| aggregate reactions.",
             ha="center", va="bottom", fontsize=8, color=INK_MUTED, linespacing=1.5)
    fig.tight_layout(rect=(0, 0.09, 1, 1))
    out = OUT_DIR / "fig_agreement_vs_snr.png"
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)
    return out


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(DATA_DIR / "reaction_features.tsv", sep="\t", low_memory=False)
    t = df[df["n_sources"] == 3].copy()
    print(f"{len(t)} three-source reactions (post-mask)")

    # The mask-evidence panel must show what the mask removed, so it reads the
    # unfiltered table. Everything else reads the filtered one.
    unmasked_path = DATA_DIR / "reaction_features_unmasked.tsv"
    if unmasked_path.exists():
        u = pd.read_csv(unmasked_path, sep="\t", low_memory=False)
        u = u[u["n_sources"] == 3].copy()
        print(f"{len(u)} three-source reactions (pre-mask, for the evidence panel)")
        print("wrote", fig_kegg_trust(u))
    else:
        print(f"  SKIP fig_dg_agreement_by_kegg_trust: {unmasked_path} not built "
              f"(run build_thermo_agreement_features.py --no-dgp-mask)")
    for fn in (fig_reaction_class, fig_snr):
        print("wrote", fn(t))


if __name__ == "__main__":
    main()
