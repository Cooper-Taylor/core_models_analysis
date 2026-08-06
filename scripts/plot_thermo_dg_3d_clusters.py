#!/usr/bin/env python3
"""3D views of the three-source DeltaG space with the HDBSCAN clustering from
cluster_thermo_dg_3d.py drawn on top.

Both figures are drawn AFTER the dGPredictor KEGG mask is applied, so every
reaction shown has a dGPredictor value backed by a ModelSEED-vouched KEGG
mapping. (Before the mask, the dominant clusters were mis-mapping artifacts --
blocks of reactions sharing one broadcast KEGG prediction. Those are gone, and
what is left is chemistry.)

  fig_dg_3d_clusters.png
      Every reaction carrying a DeltaG'deg from Group Contribution,
      eQuilibrator and dGPredictor, plotted on those three axes, colored by the
      dominant chemistry of the HDBSCAN cluster it landed in -- redox,
      group transfer, other, or HDBSCAN noise.

  fig_dg_3d_cluster_panels.png
      Four clusters shown on their own, chosen to span the concordant-to-
      discordant range that survives the mask.

Axes are drawn on a symmetric-log scale, sgn(x)*log10(1+|x|), with ticks
relabelled in kcal/mol: the distribution runs from 0 to ~10^4 kcal/mol and a
linear axis would collapse all of metabolism onto one pixel.

Palette: 3 categorical hues (validated all-pairs, light mode) plus neutral gray
for the noise category, matching the convention already used in
plot_thermo_source_dg_scatter.py.
"""
from __future__ import annotations

import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
DATA_DIR = ANALYSIS_DIR / "results" / "thermo_agreement"
OUT_DIR = ANALYSIS_DIR / "reports" / "thermoComparison" / "figures" / "thermo_agreement"

BLUE, ORANGE, AQUA, NEUTRAL = "#2a78d6", "#eb6834", "#1baf7a", "#9c9a94"
INK_PRIMARY, INK_SECONDARY, INK_MUTED = "#0b0b0b", "#52514e", "#898781"
SURFACE, PANE = "#fcfcfb", "#f4f3ee"

SPACE = "raw"          # clustering whose labels are drawn
TICKS = [-10000, -1000, -100, -10, 0, 10, 100, 1000, 10000]
# Panel figure zoom. 99.3% of the three-source set sits inside +/-1000 kcal/mol
# on all three axes; the panels use that window so the clusters of interest are
# not reduced to a few pixels by a handful of 10^4 kcal/mol reactions.
PANEL_LIM = symlog_lim = np.log10(1 + 1000.0)


def symlog(a):
    return np.sign(a) * np.log10(1.0 + np.abs(a))


def style3d(ax) -> None:
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_pane_color((0.957, 0.953, 0.933, 1.0))
        axis._axinfo["grid"]["color"] = (0.88, 0.877, 0.855, 1.0)
        axis._axinfo["grid"]["linewidth"] = 0.7
    ax.tick_params(colors=INK_MUTED, labelsize=7.5, pad=0)
    t = symlog(np.array(TICKS, float))
    labels = [("0" if v == 0 else f"{v:,.0f}") for v in TICKS]
    for setter, lsetter in ((ax.set_xticks, ax.set_xticklabels),
                            (ax.set_yticks, ax.set_yticklabels),
                            (ax.set_zticks, ax.set_zticklabels)):
        setter(t)
        lsetter(labels)


def label_axes(ax, size: int = 8, short: bool = False) -> None:
    """Full source names on the single-axes figure; three-letter tags on the
    small multiples, where full names collide with the tick labels."""
    names = (("GC", "eQ", "dGP") if short else
             ("Group Contribution\nΔG′° (kcal/mol)",
              "eQuilibrator\nΔG′° (kcal/mol)",
              "dGPredictor\nΔG′° (kcal/mol)"))
    ax.set_xlabel(names[0], color=INK_SECONDARY, fontsize=size,
                  labelpad=6 if short else 2, linespacing=1.3)
    ax.set_ylabel(names[1], color=INK_SECONDARY, fontsize=size,
                  labelpad=6 if short else 2, linespacing=1.3)
    ax.set_zlabel(names[2], color=INK_SECONDARY, fontsize=size,
                  labelpad=6 if short else 2, linespacing=1.3)


def classify_clusters(clu: pd.DataFrame, feat: pd.DataFrame) -> pd.Series:
    """Per-reaction category: noise, or the dominant chemistry of its cluster.

    A cluster is called redox or group-transfer when a majority of its members
    carry that chemistry; ties and everything else fall to "other". Three
    categorical hues plus neutral gray for noise is the all-pairs-safe budget
    for a scatter (see module docstring).
    """
    cols = ["rxn", "cof_nad", "cof_nadh", "cof_nadp", "cof_nadph", "ec_class",
            "cof_atp", "cof_adp", "cof_amp", "cof_sam", "cof_sah",
            "d_phosphoanhydride", "d_thioester"]
    merged = clu.merge(feat[cols], on="rxn", how="left", suffixes=("", "_f"))
    redox = (merged[["cof_nad", "cof_nadh", "cof_nadp", "cof_nadph"]].max(axis=1) == 1) | \
            merged["ec_class"].fillna("none").astype(str).str.contains("1")
    transfer = (merged[["cof_atp", "cof_adp", "cof_amp", "cof_sam",
                        "cof_sah"]].max(axis=1) == 1) | \
               (merged["d_phosphoanhydride"] != 0) | (merged["d_thioester"] != 0)
    merged["_redox"] = (redox & ~transfer).astype(int)
    merged["_transfer"] = (transfer & ~redox).astype(int)
    rate_r = merged.groupby("cluster")["_redox"].mean()
    rate_t = merged.groupby("cluster")["_transfer"].mean()
    out = []
    for lab in merged["cluster"]:
        if lab == -1:
            out.append("Not clustered (HDBSCAN noise)")
        elif rate_r[lab] >= 0.5 and rate_r[lab] > rate_t[lab]:
            out.append("Cluster: mostly redox")
        elif rate_t[lab] >= 0.5 and rate_t[lab] > rate_r[lab]:
            out.append("Cluster: mostly group transfer")
        else:
            out.append("Cluster: mixed / other chemistry")
    return pd.Series(out, index=merged.index)


def fig_overview(clu: pd.DataFrame, feat: pd.DataFrame) -> Path:
    clu = clu.copy()
    clu["kind"] = classify_clusters(clu, feat).to_numpy()
    order = ["Not clustered (HDBSCAN noise)",
             "Cluster: mostly redox",
             "Cluster: mixed / other chemistry",
             "Cluster: mostly group transfer"]
    color = {order[0]: NEUTRAL, order[1]: AQUA, order[2]: BLUE, order[3]: ORANGE}

    fig = plt.figure(figsize=(9.6, 8.4), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor(SURFACE)
    for name in order:
        sub = clu[clu["kind"] == name]
        ax.scatter(symlog(sub["dg_gc"]), symlog(sub["dg_eq"]),
                   symlog(sub["dg_dgp"]),
                   s=4 if name == order[0] else 9, linewidths=0,
                   color=color[name], alpha=0.16 if name == order[0] else 0.70,
                   label=f"{name}  (n = {len(sub):,})", depthshade=False)
    style3d(ax)
    label_axes(ax, size=8.5)
    ax.view_init(elev=20, azim=-58)

    ax.legend(loc="upper left", frameon=False, fontsize=8.8,
              labelcolor=INK_PRIMARY, markerscale=2.6,
              bbox_to_anchor=(-0.06, 0.99))
    fig.suptitle("Three-source ΔG′° space after the KEGG mask:\n"
                 "what remains are chemistry clusters on the diagonal",
                 color=INK_PRIMARY, fontsize=13, y=0.975)
    fig.text(0.5, 0.02,
             f"All {len(clu):,} reactions carrying all three sources after the dGPredictor "
             "KEGG mask.\nAxes are symmetric-log, tick labels in kcal/mol. Per-cluster table: "
             "cluster_profiles_raw.tsv",
             ha="center", va="bottom", fontsize=8, color=INK_MUTED, linespacing=1.6)
    fig.subplots_adjust(left=0.02, right=0.97, top=0.92, bottom=0.13)
    out = OUT_DIR / "fig_dg_3d_clusters.png"
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)
    return out


def fig_panels(clu: pd.DataFrame, prof: pd.DataFrame, feat: pd.DataFrame) -> Path:
    """Four named clusters, each on its own axes (single series per panel)."""
    # Data-driven so the figure survives a re-clustering: the two largest
    # clusters the sources agree on (median |GC-dGP| < 2 kcal/mol) and the two
    # largest they disagree on (> 10), labelled from their own enriched
    # features rather than hardcoded cluster numbers.
    real = prof[prof["cluster"] >= 0]
    conc = real[real["mad_gc_dgp"] < 2].nlargest(2, "n")
    disc = real[real["mad_gc_dgp"] > 10].nlargest(2, "n")
    picks = []
    for frame, col, tag in ((conc, AQUA, "sources agree"), (disc, ORANGE, "sources disagree")):
        for _, row in frame.iterrows():
            feats = str(row.get("enriched", "") or "-").split(";")[:2]
            desc = ", ".join(f.split(" ")[0].replace("cof_", "").replace("d_", "")
                                .replace("_nz", "") for f in feats if f.strip() and f.strip() != "-")
            picks.append((int(row["cluster"]),
                          f"C{int(row['cluster'])} · {desc or 'mixed'}\n{tag}",
                          col))
    if not picks:
        print("  no clusters matched the panel criteria; skipping panel figure")
        return OUT_DIR / "fig_dg_3d_cluster_panels.png"

    pinfo = prof.set_index("cluster")
    fig = plt.figure(figsize=(12.0, 10.2), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    for i, (lab, title, col) in enumerate(picks[:4]):
        ax = fig.add_subplot(2, 2, i + 1, projection="3d")
        ax.set_facecolor(SURFACE)
        rest = clu[clu["cluster"] != lab]
        sub = clu[clu["cluster"] == lab]
        ax.scatter(symlog(rest["dg_gc"]), symlog(rest["dg_eq"]),
                   symlog(rest["dg_dgp"]), s=2, linewidths=0,
                   color=NEUTRAL, alpha=0.09, depthshade=False)
        ax.scatter(symlog(sub["dg_gc"]), symlog(sub["dg_eq"]),
                   symlog(sub["dg_dgp"]), s=20, linewidths=0,
                   color=col, alpha=0.85, depthshade=False)
        style3d(ax)
        label_axes(ax, size=8, short=True)
        for setter in (ax.set_xlim, ax.set_ylim, ax.set_zlim):
            setter(-PANEL_LIM, PANEL_LIM)
        ax.view_init(elev=20, azim=-58)
        row = pinfo.loc[lab]
        ax.set_title(
            f"{title}\n"
            f"n = {int(row['n'])} · KEGG vouched {row['frac_kegg_vouched']:.0%} · "
            f"median |GC − dGP| = {row['mad_gc_dgp']:.1f} kcal/mol\n"
            f"median ΔG′°:  GC {row['gc_median']:.1f}   eQ {row['eq_median']:.1f}   "
            f"dGP {row['dgp_median']:.1f}",
            color=INK_PRIMARY, fontsize=8.8, pad=-14, linespacing=1.5)
    fig.suptitle("Post-mask clusters: what the three sources agree on, and what they do not",
                 color=INK_PRIMARY, fontsize=12.5, y=0.992)
    fig.text(0.5, 0.012,
             "Axes GC / eQ / dGP are the three sources' ΔG′° in kcal/mol, symmetric-log, "
             "windowed to ±1,000 kcal/mol.\nGray is every other three-source reaction.",
             ha="center", va="bottom", fontsize=8, color=INK_MUTED, linespacing=1.5)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.90, bottom=0.055,
                        wspace=0.02, hspace=0.30)
    out = OUT_DIR / "fig_dg_3d_cluster_panels.png"
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)
    return out


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    clu = pd.read_csv(DATA_DIR / f"clusters_3d_{SPACE}.tsv", sep="\t")
    prof = pd.read_csv(DATA_DIR / f"cluster_profiles_{SPACE}.tsv", sep="\t")
    feat = pd.read_csv(DATA_DIR / "reaction_features.tsv", sep="\t", low_memory=False)
    print(f"{len(clu)} reactions, "
          f"{clu['cluster'].nunique() - 1} clusters, "
          f"{(clu['cluster'] == -1).mean():.1%} noise")
    print("wrote", fig_overview(clu, feat))
    print("wrote", fig_panels(clu, prof, feat))


if __name__ == "__main__":
    main()
