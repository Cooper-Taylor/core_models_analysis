#!/usr/bin/env python3
"""Figure for the per-reaction source assignment.

  A  calibration: each source's expected error against its own reported sigma,
     with the TECRDB-covered (gold) range marked -- the two-tier fit exists
     because gold data stops well short of the range the model must work over.
  B  held-out experimental validation against the fixed-source baselines and
     the incumbent dev priority.
  C  coverage vs error tolerance, and which source gets chosen where.
"""
import json, os, sys
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
DATA = Path(os.environ.get("EQDGP_OUT", str(ANALYSIS_DIR / "results" / "eq_vs_dgpms")))
OUT_DIR = Path(os.environ.get("EQDGP_FIGS", str(
    ANALYSIS_DIR / "reports" / "thermoComparison" / "figures" / "eq_vs_dgpms")))
BLUE, ORANGE, AQUA, NEUTRAL = "#2a78d6", "#eb6834", "#1baf7a", "#9c9a94"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#898781"
GRID, BASE, SURF = "#e1e0d9", "#c3c2b7", "#fcfcfb"
COL = {"EQ": BLUE, "DGPMS": ORANGE, "GC": AQUA}
NAME = {"EQ": "eQuilibrator", "DGPMS": "dGPredictor-MS", "GC": "Group contribution"}

def style(ax):
    ax.set_facecolor(SURF); ax.grid(True, color=GRID, lw=0.8, zorder=0)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"): ax.spines[s].set_color(BASE)
    ax.tick_params(colors=INK3, labelsize=9)

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    mdl = json.loads((DATA / "source_assignment_models.json").read_text())
    fr = pd.read_csv(DATA / "source_assignment_frontier.tsv", sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(15.6, 5.0), dpi=150,
                             gridspec_kw={"width_ratios": [1.05, 1, 1.15]})
    fig.patch.set_facecolor(SURF)

    ax = axes[0]; style(ax)
    for k, m in mdl["models"].items():
        if m["kind"] != "isotonic": continue
        ax.plot(m["x"], m["y"], "-", lw=2.2, color=COL[k], zorder=3,
                label=f"{NAME[k]}  (gold n={m['n_gold']}, silver n={m['n_silver']:,})")
        if np.isfinite(m.get("gold_sigma_p90", np.nan)):
            ax.axvline(m["gold_sigma_p90"], color=COL[k], ls=":", lw=1.2, alpha=0.7)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("source's own reported σ  (kcal/mol, log)", color=INK2, fontsize=10)
    ax.set_ylabel("expected |error|  ê  (kcal/mol, log)", color=INK2, fontsize=10)
    ax.legend(loc="upper left", frameon=False, fontsize=7.4, labelcolor=INK)
    ax.set_title("A · calibration, isotonic in σ\n(dotted = where TECRDB coverage ends)",
                 color=INK, fontsize=10.5, pad=8)

    ax = axes[1]; style(ax)
    v = pd.DataFrame(mdl["validation"]).sort_values("mean_abs_err")
    ys = np.arange(len(v))
    cols = [ORANGE if "assignment" in s else NEUTRAL for s in v.strategy]
    ax.barh(ys - 0.19, v.median_abs_err, height=0.36, color=cols, zorder=3, alpha=0.55,
            label="median")
    ax.barh(ys + 0.19, v.mean_abs_err, height=0.36, color=cols, zorder=3, label="mean")
    ax.set_yticks(ys)
    ax.set_yticklabels([s.replace(" (this script)", " ←") for s in v.strategy],
                       fontsize=8, color=INK2)
    ax.set_xlabel("|chosen source − experiment|  (kcal/mol)", color=INK2, fontsize=10)
    ax.legend(loc="lower right", frameon=False, fontsize=8, labelcolor=INK)
    ax.set_title(f"B · held-out TECRDB (n={int(v.n_scored.max())})\nassignment wins on mean",
                 color=INK, fontsize=10.5, pad=8)

    ax = axes[2]; style(ax)
    f = fr[fr.tolerance < 1e8].sort_values("tolerance")
    bot = np.zeros(len(f))
    for k in ("EQ", "DGPMS", "GC"):
        ax.bar(np.arange(len(f)), f[f"n_{k}"], bottom=bot, color=COL[k], zorder=3,
               width=0.72, label=NAME[k])
        bot += f[f"n_{k}"].to_numpy()
    ax.axhline(3246, color=INK, ls="--", lw=1.3, zorder=4)
    ax.text(0.1, 3900, "consensus-subset approach (3,246)", fontsize=7.6, color=INK)
    ax.set_xticks(np.arange(len(f)))
    ax.set_xticklabels([f"{t:g}" for t in f.tolerance], fontsize=8.5)
    ax.set_xlabel("expected-error tolerance ê  (kcal/mol)", color=INK2, fontsize=10)
    ax.set_ylabel("reactions assigned", color=INK2, fontsize=10)
    ax.legend(loc="upper left", frameon=False, fontsize=8, labelcolor=INK)
    ax.set_title("C · coverage and which source wins\nceiling 32,466 = every reaction with any source",
                 color=INK, fontsize=10.5, pad=8)

    fig.suptitle("Per-reaction source assignment: use the source expected to be most "
                 "accurate, keep it if that is good enough",
                 color=INK, fontsize=13, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    out = OUT_DIR / "fig8_source_assignment.png"
    fig.savefig(out, facecolor=SURF); plt.close(fig)
    print("wrote", out)

if __name__ == "__main__":
    main()
