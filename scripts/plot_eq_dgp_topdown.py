#!/usr/bin/env python3
"""Top-down figure set for eQuilibrator vs dGPredictor-ModelSEED.

Ordered the way the report reads: chemistry first, then the dominant failure,
then the two quantities that need explaining before any metabolite claim is
believable.

  fig1_reaction_class.png   Which ORGANIC transformation classes carry the
      disagreement. Absolute and relative error side by side, because the
      ranking flips between them.
  fig2_quinone.png          The dominant class in detail: the quinone/quinol
      couple, its sign error, and the EC 1.x.5.x enzymes it lives in.
  fig3_sigma.png            What sigma IS and whether it can be trusted: the
      model's own posterior SD, its calibration curve, and the resulting tiers.
  fig4_gauge.png            What a "fitted offset" is and why a large one can
      mean nothing -- the null-space demonstration and the fitted-vs-observed
      plane that separates real disagreement from bookkeeping.
  fig5_metabolites.png      Metabolites, ranked by OBSERVED disagreement, with
      the fitted offset shown alongside so the two can be compared.

Palette: the project's existing validated 3-hue categorical triple plus neutral
gray. No new hues.
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
DATA = Path(os.environ.get("EQDGP_OUT", str(ANALYSIS_DIR / "results" / "eq_vs_dgpms")))
OUT_DIR = Path(os.environ.get("EQDGP_FIGS",
                              str(ANALYSIS_DIR / "reports" / "thermoComparison" / "figures" / "eq_vs_dgpms")))

BLUE, ORANGE, AQUA, NEUTRAL = "#2a78d6", "#eb6834", "#1baf7a", "#9c9a94"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#898781"
GRID, BASE, SURF = "#e1e0d9", "#c3c2b7", "#fcfcfb"


def style(ax):
    ax.set_facecolor(SURF)
    ax.grid(True, color=GRID, linewidth=0.8, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(BASE)
    ax.tick_params(colors=INK3, labelsize=9)


def short(s, n=38):
    return s if len(s) <= n else s[:n - 1] + "…"


# --------------------------------------------------------------------- fig 1
def fig1(cls: pd.DataFrame, baseline: float) -> Path:
    c = cls.sort_values("median_absdiff")
    fig, axes = plt.subplots(1, 2, figsize=(14.0, 7.0), dpi=150,
                             gridspec_kw={"width_ratios": [1.25, 1]})
    fig.patch.set_facecolor(SURF)
    ys = np.arange(len(c))

    ax = axes[0]; style(ax)
    ax.axvline(baseline, color=BASE, ls="--", lw=1.3, zorder=1)
    ax.text(baseline, len(c) - 0.3, f" subset median\n {baseline:.2f}", fontsize=8,
            color=INK3, va="top")
    cols = [ORANGE if v >= 2 * baseline else (AQUA if v < baseline else NEUTRAL)
            for v in c["median_absdiff"]]
    ax.hlines(ys, 0, c["median_absdiff"], color=BASE, lw=1.2, zorder=2)
    ax.scatter(c["median_absdiff"], ys, s=[28 + 3.2 * np.sqrt(n) for n in c["n"]],
               color=cols, zorder=3)
    for y, (_, r) in zip(ys, c.iterrows()):
        ax.text(r["median_absdiff"] * 1.06 + 0.4, y,
                f"n={int(r['n']):,}  {r['vs_baseline']:.1f}× baseline",
                va="center", fontsize=7.6, color=INK3)
    ax.set_yticks(ys); ax.set_yticklabels([short(x) for x in c["class"]],
                                          fontsize=8.8, color=INK2)
    ax.set_xscale("log")
    ax.set_xlim(1, c["median_absdiff"].max() * 3.2)
    ax.set_xlabel("median |eQuilibrator − dGPredictor-ModelSEED|  (kcal/mol, log)",
                  color=INK2, fontsize=10)
    ax.set_title("ABSOLUTE disagreement", color=INK, fontsize=11.5, pad=10)

    # Relative error reorders the list -- oxygenation drops out, because those
    # reactions are simply large.
    ax = axes[1]; style(ax)
    c2 = c.sort_values("rel_error")
    ys2 = np.arange(len(c2))
    cols2 = [ORANGE if v >= 1.0 else (AQUA if v < 0.6 else NEUTRAL) for v in c2["rel_error"]]
    ax.hlines(ys2, 0, c2["rel_error"] * 100, color=BASE, lw=1.2, zorder=2)
    ax.scatter(c2["rel_error"] * 100, ys2,
               s=[28 + 3.2 * np.sqrt(n) for n in c2["n"]], color=cols2, zorder=3)
    for y, (_, r) in zip(ys2, c2.iterrows()):
        ax.text(r["rel_error"] * 100 + 6, y, f"|ΔG| ≈ {r['median_abs_dg_eq']:.0f}",
                va="center", fontsize=7.4, color=INK3)
    ax.set_yticks(ys2); ax.set_yticklabels([short(x, 30) for x in c2["class"]],
                                           fontsize=8.2, color=INK2)
    ax.set_xlabel("median |Δ| as % of the reaction's own |ΔG′°|", color=INK2, fontsize=10)
    ax.set_xlim(0, min(c2["rel_error"].max() * 100 * 1.25, 420))
    ax.set_title("RELATIVE disagreement", color=INK, fontsize=11.5, pad=10)

    fig.suptitle("Where the two methods disagree, by organic transformation class",
                 color=INK, fontsize=13.5, y=0.985)
    fig.text(0.5, 0.008,
             "Each reaction gets exactly one class from a priority cascade over what bonds change "
             "(organic_reaction_types.py), not from its EC number.\nDot area ∝ number of reactions. "
             "Isomerisation's high relative error is an artefact of its near-zero ΔG′° (median 0.5 kcal/mol), "
             "not a real failure.",
             ha="center", va="bottom", fontsize=8, color=INK3, linespacing=1.5)
    fig.tight_layout(rect=(0, 0.05, 1, 0.94))
    out = OUT_DIR / "fig1_reaction_class.png"
    fig.savefig(out, facecolor=SURF); plt.close(fig)
    return out


# --------------------------------------------------------------------- fig 2
def fig2(key: pd.DataFrame, enz: pd.DataFrame) -> Path:
    q = key["chem_class"] == "Redox: quinone / quinol"
    fig, axes = plt.subplots(1, 3, figsize=(15.4, 5.2), dpi=150,
                             gridspec_kw={"width_ratios": [1.05, 1, 1.15]})
    fig.patch.set_facecolor(SURF)

    ax = axes[0]; style(ax)
    W = 200
    x, y = key["dg_eq"].to_numpy(float), key["dg_dgp"].to_numpy(float)
    k = (np.abs(x) <= W) & (np.abs(y) <= W)
    ax.plot([-W, W], [-W, W], "--", lw=1.2, color=BASE, zorder=1)
    ax.axhline(0, color=BASE, lw=0.8, zorder=1)
    ax.scatter(x[k & ~q], y[k & ~q], s=7, lw=0, color=NEUTRAL, alpha=0.3, zorder=2,
               label=f"all other classes (n={int((~q).sum()):,})")
    ax.scatter(x[k & q], y[k & q], s=16, lw=0, color=ORANGE, alpha=0.75, zorder=3,
               label=f"quinone / quinol (n={int(q.sum()):,})")
    ax.set_xlim(-W, W); ax.set_ylim(-W, W)
    ax.set_xlabel("eQuilibrator ΔG′° (kcal/mol)", color=INK2, fontsize=10)
    ax.set_ylabel("dGPredictor-ModelSEED ΔG′° (kcal/mol)", color=INK2, fontsize=10)
    ax.legend(loc="upper left", frameon=False, fontsize=8, labelcolor=INK, markerscale=1.7)
    ax.set_title("The quinone class sits off the diagonal,\noften on the wrong side of zero",
                 color=INK, fontsize=10.5, pad=8)

    ax = axes[1]; style(ax)
    sub = key[q]
    wrong = (np.sign(sub["dg_eq"]) != np.sign(sub["dg_dgp"])).mean() * 100
    other = (np.sign(key.loc[~q, "dg_eq"]) != np.sign(key.loc[~q, "dg_dgp"])).mean() * 100
    ax.bar([0, 1], [wrong, other], color=[ORANGE, NEUTRAL], width=0.55, zorder=3)
    for i, v in enumerate([wrong, other]):
        ax.text(i, v + 1.4, f"{v:.1f}%", ha="center", fontsize=11, color=INK)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["quinone /\nquinol", "everything\nelse"],
                                              fontsize=9, color=INK2)
    ax.set_ylabel("% where the two disagree on the SIGN of ΔG′°", color=INK2, fontsize=9.5)
    ax.set_ylim(0, max(wrong, other) * 1.32)
    ax.set_title("It is a direction error,\nnot just a magnitude error",
                 color=INK, fontsize=10.5, pad=8)

    ax = axes[2]; style(ax)
    e = enz.nlargest(7, "median_absdiff").sort_values("median_absdiff")
    ys = np.arange(len(e))
    ax.barh(ys, e["median_absdiff"], color=[ORANGE if "quinone" in c else NEUTRAL
                                            for c in e["dominant_class"]],
            height=0.62, zorder=3)
    ax.set_yticks(ys)
    ax.set_yticklabels([f"EC {r.ec3}  (n={int(r.n)})" for r in e.itertuples()],
                       fontsize=8.4, color=INK2)
    for y, r in zip(ys, e.itertuples()):
        ax.text(1.5, y, short(r.example_enzyme, 34), va="center", fontsize=7.2,
                color=SURF if r.median_absdiff > 25 else INK3)
    ax.set_xlabel("median |Δ| (kcal/mol)", color=INK2, fontsize=9.5)
    ax.set_title("EC 1.x.5.x — the nomenclature's own\n\"quinone as acceptor\" families",
                 color=INK, fontsize=10.5, pad=8)

    fig.suptitle("The dominant failure: the quinone / hydroquinone two-electron couple",
                 color=INK, fontsize=13, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    out = OUT_DIR / "fig2_quinone.png"
    fig.savefig(out, facecolor=SURF); plt.close(fig)
    return out


# --------------------------------------------------------------------- fig 3
def fig3(key: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.3), dpi=150)
    fig.patch.set_facecolor(SURF)

    ax = axes[0]; style(ax)
    d = key.dropna(subset=["dgp_uncertainty"]).copy()
    d["bin"] = pd.qcut(d["dgp_uncertainty"], 8, duplicates="drop")
    g = d.groupby("bin", observed=True).agg(rep=("dgp_uncertainty", "median"),
                                            obs=("absdiff", "median"),
                                            n=("absdiff", "size")).reset_index()
    lim = max(g["rep"].max(), g["obs"].max()) * 1.15
    ax.plot([0, lim], [0, lim], "--", lw=1.2, color=BASE, zorder=1,
            label="perfect calibration")
    ax.plot(g["rep"], g["obs"], "-o", color=ORANGE, lw=2, ms=7, zorder=3)
    for _, r in g.iterrows():
        ax.annotate(f"n={int(r['n']):,}", (r["rep"], r["obs"]), fontsize=7,
                    color=INK3, xytext=(4, -9), textcoords="offset points")
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_xlabel("σ the model REPORTS  (kcal/mol)", color=INK2, fontsize=10)
    ax.set_ylabel("|eQ − dGP| actually observed  (kcal/mol)", color=INK2, fontsize=10)
    ax.legend(loc="upper left", frameon=False, fontsize=8.5, labelcolor=INK)
    ax.set_title("σ is informative (ρ = 0.67) and conservative:\n"
                 "it over-states its error, but monotonically",
                 color=INK, fontsize=10.5, pad=8)

    ax = axes[1]; style(ax)
    tiers = [("σ ≤ 3", key["dgp_uncertainty"] <= 3, AQUA),
             ("3 < σ ≤ 30", key["dgp_uncertainty"].between(3, 30, "right"), NEUTRAL),
             ("σ > 30", key["dgp_uncertainty"] > 30, ORANGE)]
    W = 120
    for i, (lab, m, col) in enumerate(tiers):
        s = key[m.to_numpy()]
        xx, yy = s["dg_eq"].to_numpy(float), s["dg_dgp"].to_numpy(float)
        kk = (np.abs(xx) <= W) & (np.abs(yy) <= W)
        ax.scatter(xx[kk], yy[kk] + 0, s=6, lw=0, color=col, alpha=0.45, zorder=2 + i,
                   label=f"{lab}   n={len(s):,},  r={np.corrcoef(xx, yy)[0, 1]:.3f},  "
                         f"median |Δ|={np.median(np.abs(xx - yy)):.2f}")
    ax.plot([-W, W], [-W, W], "--", lw=1.2, color=BASE, zorder=1)
    ax.set_xlim(-W, W); ax.set_ylim(-W, W)
    ax.set_xlabel("eQuilibrator ΔG′° (kcal/mol)", color=INK2, fontsize=10)
    ax.set_ylabel("dGPredictor-ModelSEED ΔG′° (kcal/mol)", color=INK2, fontsize=10)
    ax.legend(loc="upper left", frameon=False, fontsize=7.8, labelcolor=INK, markerscale=2.2)
    ax.set_title("So it can be tiered on its own output,\nwith no external evidence",
                 color=INK, fontsize=10.5, pad=8)

    fig.suptitle("σ = dGPredictor-ModelSEED's own BayesianRidge posterior SD, "
                 "stored per reaction in ModelSEED",
                 color=INK, fontsize=12.5, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    out = OUT_DIR / "fig3_sigma.png"
    fig.savefig(out, facecolor=SURF); plt.close(fig)
    return out


# --------------------------------------------------------------------- fig 4
def fig4(met: pd.DataFrame, gauge: pd.DataFrame, baseline: float) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 6.0), dpi=150,
                             gridspec_kw={"width_ratios": [1, 1.15]})
    fig.patch.set_facecolor(SURF)

    ax = axes[0]; style(ax)
    ax.axis("off")
    ax.set_facecolor(SURF)
    txt = (
        "A “fitted offset” is NOT a regression of one method's ΔG on the other.\n\n"
        "For every reaction r take the disagreement\n"
        "        d$_r$ = ΔG$_{eq}$(r) − ΔG$_{dgp}$(r)\n"
        "and solve, by least squares over the stoichiometric matrix S,\n"
        "        d$_r$  ≈  Σ$_i$  ν$_{ir}$ · x$_i$\n"
        "where ν$_{ir}$ is compound i's coefficient in reaction r.\n"
        "x$_i$ is that compound's offset. Legitimate because both methods are\n"
        "additive over compound formation energies, so their difference is too.\n\n"
        "THE CATCH — the solution is not unique. S has a null space: any z with\n"
        "S·z = 0 can be added to x and every predicted reaction value is\n"
        "unchanged. Element conservation supplies such z for free. Adding a\n"
        "fixed amount per atom to every compound changes nothing:\n"
    )
    ax.text(0, 1, txt, transform=ax.transAxes, va="top", ha="left", fontsize=9.2,
            color=INK2, linespacing=1.65)
    tbl = "     element     max change over all 11,097 reactions\n" + "\n".join(
        f"         {r.element}            {r.max_abs_change_kcal:.1f} kcal/mol"
        f"      ({r.frac_reactions_unchanged:.0%} unchanged)"
        for r in gauge.itertuples())
    ax.text(0, 0.30, tbl, transform=ax.transAxes, va="top", ha="left",
            fontsize=8.6, color=INK, family="monospace", linespacing=1.6)
    ax.text(0, 0.06, "→ an individual offset can move by tens of kcal/mol with no\n"
                     "   observable consequence. A large fitted value proves nothing\n"
                     "   on its own; it has to be checked against the reactions.",
            transform=ax.transAxes, va="top", ha="left", fontsize=9.2,
            color=ORANGE, linespacing=1.65)
    ax.set_title("What a fitted offset is, and why it can lie",
                 color=INK, fontsize=11.5, pad=10, loc="left")

    ax = axes[1]; style(ax)
    m = met.dropna(subset=["fitted_offset", "observed_median_absdiff"]).copy()
    real = m["ratio_vs_baseline"] >= 2
    ax.axhline(baseline, color=BASE, ls="--", lw=1.2, zorder=1)
    ax.text(1.2, baseline * 1.06, f"subset baseline {baseline:.2f}", fontsize=8, color=INK3)
    ax.scatter(m.loc[~real, "fitted_offset"].abs(), m.loc[~real, "observed_median_absdiff"],
               s=26, lw=0, color=NEUTRAL, alpha=0.6, zorder=2,
               label="not distinguishable from baseline")
    ax.scatter(m.loc[real, "fitted_offset"].abs(), m.loc[real, "observed_median_absdiff"],
               s=44, lw=0, color=ORANGE, alpha=0.85, zorder=3,
               label="real disagreement (≥ 2× baseline)")
    # The prenyl-quinone variants share coordinates exactly (same reaction set),
    # so label the family once instead of stacking six strings on one point.
    seen_xy: set = set()
    for _, r in m.nlargest(14, "observed_median_absdiff").iterrows():
        xk, yk = round(abs(r["fitted_offset"]), 1), round(r["observed_median_absdiff"], 1)
        if (xk, yk) in seen_xy:
            continue
        seen_xy.add((xk, yk))
        lab = ("ubiquinone / quinol ×6" if "biquino" in str(r["name"]).lower()
               or "idecarenone" in str(r["name"]).lower() else short(str(r["name"]), 20))
        ax.annotate(lab, (abs(r["fitted_offset"]), r["observed_median_absdiff"]),
                    fontsize=7.4, color=INK, xytext=(-6, 7), textcoords="offset points",
                    ha="right")
        if len(seen_xy) >= 6:
            break
    trap = m[(m["fitted_offset"].abs() > 30) & (m["ratio_vs_baseline"] < 1.5)]
    for _, r in trap.nlargest(4, "fitted_offset").iterrows():
        ax.annotate(short(str(r["name"]), 18), (abs(r["fitted_offset"]),
                                                r["observed_median_absdiff"]),
                    fontsize=7.4, color=BLUE, xytext=(5, -10), textcoords="offset points")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(right=m["fitted_offset"].abs().max() * 4.5)
    ax.set_xlabel("|fitted offset|  (kcal/mol, log)", color=INK2, fontsize=10)
    ax.set_ylabel("OBSERVED median |eQ − dGP| over reactions\ncontaining that compound "
                  "(kcal/mol, log)", color=INK2, fontsize=9.5)
    ax.legend(loc="upper left", frameon=False, fontsize=8.2, labelcolor=INK, markerscale=1.4)
    ax.set_title("Bottom-right = large fitted offset, no real disagreement.\n"
                 "Those (blue labels) are gauge, not chemistry.",
                 color=INK, fontsize=10.5, pad=8)

    fig.tight_layout()
    out = OUT_DIR / "fig4_gauge.png"
    fig.savefig(out, facecolor=SURF); plt.close(fig)
    return out


# --------------------------------------------------------------------- fig 5
def fig5(met: pd.DataFrame, baseline: float) -> Path:
    top = met.nlargest(14, "observed_median_absdiff").sort_values("observed_median_absdiff")
    trap = met[(met["fitted_offset"].abs() > 25) & (met["ratio_vs_baseline"] < 1.5)] \
        .reindex(met["fitted_offset"].abs().sort_values(ascending=False).index) \
        .dropna(how="all").head(8).sort_values("fitted_offset")
    fig, axes = plt.subplots(1, 2, figsize=(14.2, 6.2), dpi=150)
    fig.patch.set_facecolor(SURF)

    for ax, sub, title, col in (
            (axes[0], top, "Metabolites with REAL disagreement\n(ranked by observed, not fitted)", ORANGE),
            (axes[1], trap, "The gauge trap: big fitted offset,\nno actual disagreement", BLUE)):
        style(ax)
        ys = np.arange(len(sub))
        ax.barh(ys - 0.19, sub["observed_median_absdiff"], height=0.36, color=col,
                zorder=3, label="observed median |eQ − dGP|")
        ax.barh(ys + 0.19, sub["fitted_offset"].abs(), height=0.36, color=NEUTRAL,
                zorder=3, alpha=0.85, label="|fitted offset|")
        ax.axvline(baseline, color=BASE, ls="--", lw=1.2, zorder=2)
        ax.set_yticks(ys)
        ax.set_yticklabels([f"{short(str(n), 26)}  (n={int(k)})"
                            for n, k in zip(sub["name"], sub["n_reactions"])],
                           fontsize=8.2, color=INK2)
        ax.set_xlabel("kcal/mol", color=INK2, fontsize=10)
        ax.legend(loc="lower right", frameon=False, fontsize=8.2, labelcolor=INK)
        ax.set_title(title, color=INK, fontsize=11, pad=10)

    fig.text(0.5, 0.005,
             "Left: fitted and observed agree, so the offset is chemistry — and it is one class, the "
             "quinone/hydroquinone couple.\nRight: fitted offset is large while the reactions containing "
             "the compound sit at or below the dashed baseline. Those offsets are bookkeeping.",
             ha="center", va="bottom", fontsize=8, color=INK3, linespacing=1.5)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    out = OUT_DIR / "fig5_metabolites.png"
    fig.savefig(out, facecolor=SURF); plt.close(fig)
    return out


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    key = pd.read_csv(DATA / "key_subset_classified.tsv", sep="\t", low_memory=False)
    key["absdiff"] = (key["dg_eq"] - key["dg_dgp"]).abs()
    cls = pd.read_csv(DATA / "reaction_class_breakdown.tsv", sep="\t")
    enz = pd.read_csv(DATA / "enzyme_breakdown.tsv", sep="\t")
    met = pd.read_csv(DATA / "metabolite_validated.tsv", sep="\t")
    gauge = pd.read_csv(DATA / "gauge_demo.tsv", sep="\t")
    baseline = key["absdiff"].median()
    print("wrote", fig1(cls, baseline))
    print("wrote", fig2(key, enz))
    print("wrote", fig3(key))
    print("wrote", fig4(met, gauge, baseline))
    print("wrote", fig5(met, baseline))


if __name__ == "__main__":
    main()
