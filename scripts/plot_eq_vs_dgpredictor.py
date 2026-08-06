#!/usr/bin/env python3
"""Figures for the eQuilibrator vs dGPredictor comparison.

  fig_eq_dgp_reconciliation.png   the funnel from "both sources present" to the
      key subset, and what each filter removes.
  fig_eq_dgp_scatter.png          the key-subset scatter, plus the two distinct
      failure modes side by side (absolute disagreement vs sign disagreement).
  fig_eq_dgp_by_class.png         median |eQ - dGP| by reaction class.
  fig_eq_dgp_mechanisms.png       the three method-level effects that survive:
      O2 dose response, stereo blindness, common-metabolite anchoring.
  fig_eq_dgp_metabolites.png      per-metabolite offsets, best and worst.

Palette: the same 3-hue categorical triple + neutral gray already validated for
this project's figures (plot_thermo_source_dg_scatter.py, validate_palette.js
light --pairs all). No new hues are introduced.
"""
from __future__ import annotations

import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
DATA = Path(os.environ.get("EQDGP_OUT", str(ANALYSIS_DIR / "results" / "eq_vs_dgp")))
OUT_DIR = Path(os.environ.get("EQDGP_FIGS", str(ANALYSIS_DIR / "reports" / "thermoComparison" / "figures" / "eq_vs_dgp")))

BLUE, ORANGE, AQUA, NEUTRAL = "#2a78d6", "#eb6834", "#1baf7a", "#9c9a94"
INK_PRIMARY, INK_SECONDARY, INK_MUTED = "#0b0b0b", "#52514e", "#898781"
GRIDLINE, BASELINE, SURFACE = "#e1e0d9", "#c3c2b7", "#fcfcfb"
WINDOW = 250.0


def style(ax) -> None:
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRIDLINE, linewidth=0.8, zorder=0)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(BASELINE)
    ax.tick_params(colors=INK_MUTED, labelsize=9)


def fig_reconciliation(rec: pd.DataFrame, key: pd.DataFrame) -> Path:
    steps = [
        ("Both sources present", len(rec)),
        ("eQuilibrator compound match\nexact (full InChIKey)", int((rec["eq_all_tier1"] == 1).sum())),
        ("eQuilibrator has a real estimate\n(not a sentinel uncertainty)",
         int(((rec["eq_all_tier1"] == 1) & (rec["eq_sentinel"] == 0)).sum())),
        ("+ balanced, non-transport,\nstructures (key subset)", len(key)),
    ]
    fig, ax = plt.subplots(figsize=(9.2, 4.6), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    style(ax)
    ys = np.arange(len(steps))[::-1]
    vals = [s[1] for s in steps]
    ax.barh(ys, vals, color=[NEUTRAL, BLUE, AQUA, ORANGE], height=0.6, zorder=3)
    for y, v in zip(ys, vals):
        ax.text(v + max(vals) * 0.012, y, f"{v:,}", va="center", ha="left",
                fontsize=10, color=INK_PRIMARY)
    ax.set_yticks(ys)
    ax.set_yticklabels([s[0] for s in steps], fontsize=9, color=INK_SECONDARY)
    ax.set_xlim(0, max(vals) * 1.16)
    ax.set_xlabel("ModelSEED reactions", color=INK_SECONDARY, fontsize=10)
    ax.set_title("Reconciling eQuilibrator and dGPredictor on ModelSEED reaction identity",
                 color=INK_PRIMARY, fontsize=12, pad=12)
    fig.tight_layout()
    out = OUT_DIR / "fig_eq_dgp_reconciliation.png"
    fig.savefig(out, facecolor=SURFACE); plt.close(fig)
    return out


def fig_scatter(key: pd.DataFrame) -> Path:
    x, y = key["dg_eq"].to_numpy(float), key["dg_dgp"].to_numpy(float)
    o2 = key["cof_o2"].to_numpy() == 1
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 6.0), dpi=150)
    fig.patch.set_facecolor(SURFACE)

    ax = axes[0]
    style(ax)
    ax.plot([-WINDOW, WINDOW], [-WINDOW, WINDOW], "--", lw=1.2, color=BASELINE, zorder=1)
    keep = (np.abs(x) <= WINDOW) & (np.abs(y) <= WINDOW)
    ax.scatter(x[keep & ~o2], y[keep & ~o2], s=9, lw=0, color=BLUE, alpha=0.35,
               zorder=2, label=f"no O₂  (n = {int((~o2).sum()):,})")
    ax.scatter(x[keep & o2], y[keep & o2], s=11, lw=0, color=ORANGE, alpha=0.55,
               zorder=3, label=f"O₂-consuming  (n = {int(o2.sum()):,})")
    ax.set_xlim(-WINDOW, WINDOW); ax.set_ylim(-WINDOW, WINDOW)
    ax.set_xlabel("eQuilibrator ΔG′° (kcal/mol)", color=INK_SECONDARY, fontsize=10)
    ax.set_ylabel("dGPredictor ΔG′° (kcal/mol)", color=INK_SECONDARY, fontsize=10)
    ax.set_title(f"Key subset: n = {len(key):,},  r = {np.corrcoef(x, y)[0, 1]:.3f},\n"
                 f"median |Δ| = {np.median(np.abs(x - y)):.2f} kcal/mol",
                 color=INK_PRIMARY, fontsize=11, pad=10)
    ax.legend(loc="upper left", frameon=False, fontsize=8.5, labelcolor=INK_PRIMARY,
              markerscale=1.8)
    ax.text(0.97, 0.04, f"{int((~keep).sum())} point(s) outside ±{WINDOW:.0f}",
            transform=ax.transAxes, va="bottom", ha="right", fontsize=7.5, color=INK_MUTED)

    # The two failure modes are different questions and separate cleanly by |dG|.
    ax = axes[1]
    style(ax)
    bins = [0, 1, 2, 5, 10, 25, 50, 100, 1e9]
    labels = ["0–1", "1–2", "2–5", "5–10", "10–25", "25–50", "50–100", ">100"]
    absdg = np.abs(x)
    med_abs, sign_dis, ns = [], [], []
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = (absdg >= lo) & (absdg < hi)
        ns.append(int(m.sum()))
        med_abs.append(np.median(np.abs(x[m] - y[m])) if m.sum() else np.nan)
        sign_dis.append((np.sign(x[m]) != np.sign(y[m])).mean() * 100 if m.sum() else np.nan)
    xs = np.arange(len(labels))
    ax.plot(xs, med_abs, "-o", color=ORANGE, lw=2, ms=6, zorder=3,
            label="median |eQ − dGP|  (kcal/mol)")
    ax.set_ylabel("median |eQ − dGP|  (kcal/mol)", color=ORANGE, fontsize=10)
    ax.tick_params(axis="y", colors=ORANGE)
    ax2 = ax.twinx()
    ax2.plot(xs, sign_dis, "-s", color=BLUE, lw=2, ms=5, zorder=3,
             label="% disagreeing on sign")
    ax2.set_ylabel("% of reactions where the two disagree on SIGN",
                   color=BLUE, fontsize=10)
    ax2.tick_params(axis="y", colors=BLUE, labelsize=9)
    ax2.spines["top"].set_visible(False)
    for sp in ("left", "bottom", "right"):
        ax2.spines[sp].set_color(BASELINE)
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{l}\nn={n:,}" for l, n in zip(labels, ns)], fontsize=8)
    ax.set_xlabel("|eQuilibrator ΔG′°|  (kcal/mol)", color=INK_SECONDARY, fontsize=10)
    ax.set_title("Two different failure modes:\nabsolute error grows with |ΔG|, "
                 "direction error collapses onto ΔG ≈ 0",
                 color=INK_PRIMARY, fontsize=11, pad=10)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper center", frameon=False, fontsize=8.5,
              labelcolor=INK_PRIMARY)

    fig.tight_layout()
    out = OUT_DIR / "fig_eq_dgp_scatter.png"
    fig.savefig(out, facecolor=SURFACE); plt.close(fig)
    return out


def fig_by_class(cls: pd.DataFrame) -> Path:
    c = cls.sort_values("median_absdiff")
    fig, ax = plt.subplots(figsize=(9.0, 6.2), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    style(ax)
    ys = np.arange(len(c))[::-1]
    ax.hlines(ys, 0, c["median_absdiff"], color=BASELINE, lw=1.2, zorder=2)
    ax.scatter(c["median_absdiff"], ys, s=64, zorder=3,
               color=[ORANGE if v > 3 else AQUA for v in c["median_absdiff"]])
    for y, (_, row) in zip(ys, c.iterrows()):
        ax.text(row["median_absdiff"] + 0.12, y,
                f"n={int(row['n']):,}   r={row['r']:.2f}", va="center",
                fontsize=7.8, color=INK_MUTED)
    ax.set_yticks(ys)
    ax.set_yticklabels(c["class"], fontsize=9, color=INK_SECONDARY)
    ax.set_xlabel("median |eQuilibrator − dGPredictor|  (kcal/mol)",
                  color=INK_SECONDARY, fontsize=10)
    ax.set_xlim(0, c["median_absdiff"].max() * 1.45)
    ax.set_title("Where the two methods agree, by reaction class\n"
                 "(key subset, both sources on firm ground)",
                 color=INK_PRIMARY, fontsize=12, pad=12)
    fig.tight_layout()
    out = OUT_DIR / "fig_eq_dgp_by_class.png"
    fig.savefig(out, facecolor=SURFACE); plt.close(fig)
    return out


def fig_mechanisms(key: pd.DataFrame, stoich: dict) -> Path:
    fig, axes = plt.subplots(1, 3, figsize=(14.4, 4.9), dpi=150)
    fig.patch.set_facecolor(SURFACE)

    # (a) O2 dose response -- a fixed per-O2 offset, not a blow-up
    ax = axes[0]; style(ax)
    key = key.copy()
    key["n_o2"] = [-stoich.get(r, {}).get("cpd00007", 0.0) for r in key["rxn"]]
    groups = [0.0, 1.0, 2.0]
    meds, errs, ns = [], [], []
    for g in groups:
        s = key.loc[key["n_o2"] == g, "diff"]
        meds.append(s.median()); ns.append(len(s))
        errs.append([s.median() - s.quantile(0.25), s.quantile(0.75) - s.median()])
    ax.errorbar(groups, meds, yerr=np.array(errs).T, fmt="o-", color=ORANGE,
                lw=2, ms=8, capsize=4, zorder=3)
    ax.axhline(0, color=BASELINE, lw=1.2, ls="--", zorder=1)
    for g, m, n in zip(groups, meds, ns):
        ax.text(g, m - 1.9, f"n={n:,}", ha="center", fontsize=8, color=INK_MUTED)
    ax.set_xticks(groups)
    ax.set_xlabel("O₂ molecules consumed", color=INK_SECONDARY, fontsize=10)
    ax.set_ylabel("median (eQ − dGP)  (kcal/mol)", color=INK_SECONDARY, fontsize=10)
    ax.set_title("A fixed offset per O₂,\nnot a breakdown", color=INK_PRIMARY,
                 fontsize=11, pad=10)

    # (b) stereo blindness -- radius-1 fragments, no-stereo variant
    ax = axes[1]; style(ax)
    iso = key["ec_class"].fillna("none").astype(str).str.contains("5")
    frac = [(key.loc[iso, "dg_dgp"].abs() < 0.5).mean() * 100,
            (key.loc[~iso, "dg_dgp"].abs() < 0.5).mean() * 100]
    ax.bar([0, 1], frac, color=[ORANGE, NEUTRAL], width=0.55, zorder=3)
    for i, (f, n) in enumerate(zip(frac, [int(iso.sum()), int((~iso).sum())])):
        ax.text(i, f + 0.7, f"{f:.1f}%\nn={n:,}", ha="center", fontsize=9,
                color=INK_PRIMARY)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["EC 5\nisomerases", "everything\nelse"], fontsize=9,
                       color=INK_SECONDARY)
    ax.set_ylabel("% with |dGPredictor ΔG′°| < 0.5 kcal/mol",
                  color=INK_SECONDARY, fontsize=10)
    ax.set_ylim(0, max(frac) * 1.35)
    ax.set_title("Radius-1 fragments are stereo-blind:\nisomerases collapse to ΔG ≈ 0",
                 color=INK_PRIMARY, fontsize=11, pad=10)

    # (c) common-metabolite anchoring
    ax = axes[2]; style(ax)
    from collections import Counter
    counts = Counter(c for r in key["rxn"] for c in stoich.get(r, {}))
    common = {c for c, n in counts.items() if n >= 50}
    allc = np.array([all(c in common for c in stoich.get(r, {})) for r in key["rxn"]])
    data = [key.loc[allc, "absdiff"].to_numpy(), key.loc[~allc, "absdiff"].to_numpy()]
    bp = ax.boxplot(data, positions=[0, 1], widths=0.5, showfliers=False,
                    patch_artist=True, zorder=3)
    for patch, col in zip(bp["boxes"], [AQUA, NEUTRAL]):
        patch.set_facecolor(col); patch.set_edgecolor(INK_MUTED); patch.set_alpha(0.85)
    for el in ("medians", "whiskers", "caps"):
        for it in bp[el]:
            it.set_color(INK_SECONDARY)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"all reagents common\n(n={int(allc.sum()):,})",
                        f"otherwise\n(n={int((~allc).sum()):,})"],
                       fontsize=9, color=INK_SECONDARY)
    ax.set_ylabel("|eQ − dGP|  (kcal/mol)", color=INK_SECONDARY, fontsize=10)
    ax.set_ylim(0, 14)
    ax.text(0, data[0].mean() + 6.5, f"median\n{np.median(data[0]):.2f}", ha="center",
            fontsize=9, color=INK_PRIMARY)
    ax.text(1, np.median(data[1]) + 6.0, f"median\n{np.median(data[1]):.2f}", ha="center",
            fontsize=9, color=INK_PRIMARY)
    ax.set_title("Agreement is best where both\nmethods have measured chemistry",
                 color=INK_PRIMARY, fontsize=11, pad=10)

    fig.suptitle("What the two methods' construction predicts, tested",
                 color=INK_PRIMARY, fontsize=13, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = OUT_DIR / "fig_eq_dgp_mechanisms.png"
    fig.savefig(out, facecolor=SURFACE); plt.close(fig)
    return out


def fig_metabolites(prof: pd.DataFrame) -> Path:
    p = prof[prof["n_reactions"] >= 15].copy()
    worst = p.nlargest(14, "abs_offset").sort_values("offset")
    best = p[p["n_reactions"] >= 50].nsmallest(10, "abs_offset").sort_values("offset")
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.8), dpi=150,
                             gridspec_kw={"width_ratios": [1.15, 1]})
    fig.patch.set_facecolor(SURFACE)
    for ax, sub, title, col in (
        (axes[0], worst, "Metabolites the two methods value most differently", ORANGE),
        (axes[1], best, "Metabolites they agree on (≥50 reactions)", AQUA)):
        style(ax)
        ys = np.arange(len(sub))
        ax.barh(ys, sub["offset"], color=col, height=0.62, zorder=3, alpha=0.9)
        ax.axvline(0, color=BASELINE, lw=1.2, zorder=2)
        ax.set_yticks(ys)
        ax.set_yticklabels([f"{n}  (n={int(k)})" for n, k in
                            zip(sub["name"], sub["n_reactions"])],
                           fontsize=8.5, color=INK_SECONDARY)
        ax.set_xlabel("fitted per-metabolite offset, eQ − dGPredictor  (kcal/mol)",
                      color=INK_SECONDARY, fontsize=9.5)
        ax.set_title(title, color=INK_PRIMARY, fontsize=11, pad=10)
    fig.text(0.5, 0.005,
             "Least-squares attribution of the reaction-level disagreement to individual metabolites. "
             "Only 21.5% of the residual\ndisagreement is explained this way once eQuilibrator's "
             "no-estimate reactions are removed — most of the rest is reaction-specific.",
             ha="center", va="bottom", fontsize=8, color=INK_MUTED, linespacing=1.5)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    out = OUT_DIR / "fig_eq_dgp_metabolites.png"
    fig.savefig(out, facecolor=SURFACE); plt.close(fig)
    return out


def main() -> None:
    import glob, json
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rec = pd.read_csv(DATA / "reconciliation.tsv", sep="\t", low_memory=False)
    key = pd.read_csv(DATA / "key_subset.tsv", sep="\t", low_memory=False)
    cls = pd.read_csv(DATA / "class_breakdown.tsv", sep="\t")
    prof = pd.read_csv(DATA / "metabolite_profile.tsv", sep="\t")

    msdb = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase"))
    stoich = {}
    for path in sorted(glob.glob(str(msdb / "Biochemistry" / "reaction_*.json"))):
        for e in json.load(open(path)):
            stoich[e["id"]] = {i["compound"]: float(i.get("coefficient", 0) or 0)
                               for i in (e.get("stoichiometry") or [])}

    for fn, args in ((fig_reconciliation, (rec, key)), (fig_scatter, (key,)),
                     (fig_by_class, (cls,)), (fig_mechanisms, (key, stoich)),
                     (fig_metabolites, (prof,))):
        print("wrote", fn(*args))


if __name__ == "__main__":
    main()
