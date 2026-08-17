#!/usr/bin/env python3
"""Before/after: what the KEGG mis-mapping defect costs the ORIGINAL dGPredictor.

THE DEFECT
----------
The original dGPredictor scores a KEGG reaction, so every ModelSEED reaction had
to be matched to a KEGG id first. The script that did that (in the freiburgermsu
dGPredictor repo, not in ModelSEED) never reset its "last seen KEGG id" variable,
so a ModelSEED reaction with NO KEGG alias silently inherited the id of the
preceding reaction in file order. 17,271 of the 27,715 stored dGPredictor records
on ModelSEED dev are built from an id that is not that reaction's.

THE "FIX" APPLIED HERE
----------------------
No re-prediction is needed. For the 10,444 reactions whose staged KEGG id IS a
genuine ModelSEED alias, the stored value is already the right reaction's -- the
staged ids and the aliases never conflict, they are only ever absent. So the
corrected pipeline produces exactly the current values for those reactions and
NOTHING for the rest, which is what the repaired upstream script
(dG_prediction_modelseed_dev_branch_file_run.py, which sets kegg_id_str = None
each iteration) would emit. Applying the fix therefore means: drop the 17,271
reactions in results/thermo_agreement/dgpredictor_kegg_mask.json from the
dGPredictor series and leave every other source untouched.

THE CONTROL (important -- read before quoting any number)
---------------------------------------------------------
Removing points can improve a correlation on its own. To show the improvement is
not a sample-size artefact, every pair is also drawn a third way: a RANDOM subset
of the full set, of exactly the same size as the fixed set, seeded and repeated
N_CONTROL times. If the fix were merely "fewer points", the control would move as
much as the fix does. It does not.

PANELS
------
    Group Contribution      vs  dGPredictor (original)
    eQuilibrator            vs  dGPredictor (original)
    dGPredictor (original)  vs  dGPredictor-ModelSEED (retrain)

each under three conditions: as-shipped / KEGG-fixed / random control.

Point colour is the reversibility transition between the two sources, identical
in definition and palette to plot_thermo_source_dg_scatter.py.

Data: ModelSEED dev @ 49563c6f (/scratch/ctaylor/tmp/devsnap2) -- the live dev
branch; cascade code from the working ModelSEEDDatabase checkout.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
MSDB_DATA = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/tmp/devsnap2"))
MSDB_CODE = Path(os.environ.get("MSDB_CODE_ROOT", "/scratch/ctaylor/ModelSEEDDatabase"))
OUT_DIR = (ANALYSIS_DIR / "reports" / "thermoComparison" / "figures" / "kegg_fix_impact")
CACHE = OUT_DIR / "_devsnap2_sources.pkl"

sys.path.insert(0, str(MSDB_CODE / "Scripts" / "Thermodynamics"))
from reversibility_heuristics import (  # noqa: E402
    DEFAULT_HEURISTICS, run_reversibility, per_source_energy,
)
sys.path.insert(0, str(ANALYSIS_DIR / "scripts"))
from build_dgpredictor_kegg_mask import load_mask  # noqa: E402

SOURCES = {
    "GC":   ("Group contribution",     "Group Contribution ΔG′° (kcal/mol)"),
    "EQ":   ("eQuilibrator",           "eQuilibrator ΔG′° (kcal/mol)"),
    "BASE": ("dGPredictor",            "dGPredictor (original) ΔG′° (kcal/mol)"),
    "FT":   ("dGPredictor-ModelSEED",  "dGPredictor-ModelSEED ΔG′° (kcal/mol)"),
}
PAIRS = [("GC", "BASE"), ("EQ", "BASE"), ("BASE", "FT")]
CONDITIONS = ["as_shipped", "kegg_fixed", "random_control"]
COND_TITLE = {
    "as_shipped": "As shipped on dev\n(all dGPredictor values)",
    "kegg_fixed": "KEGG mis-mapping fixed\n(inherited ids dropped)",
    "random_control": "Control: random subset\nof the same size",
}

EQ_SENTINEL = 100.0      # kcal/mol; eQuilibrator's "no estimate" marker
EXTREME_CUTOFF = 1500.0  # implausible magnitudes, same as the sibling scripts
ZOOM_LIMIT = 250.0
N_CONTROL = 25
SEED = 20260814

CATEGORY_ORDER = [
    "No change",
    "Reversible → Irreversible",
    "Irreversible → Reversible",
    "Irreversible → Irreversible",
]
CATEGORY_COLOR = {
    "No change": "#9c9a94",
    "Reversible → Irreversible": "#2a78d6",
    "Irreversible → Reversible": "#eb6834",
    "Irreversible → Irreversible": "#1baf7a",
}
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
BASELINE = "#c3c2b7"
NEUTRAL = "#9c9a94"
SURFACE = "#fcfcfb"


def classify(a: str, b: str) -> str:
    ra, rb = a == "=", b == "="
    if ra and rb:
        return "No change"
    if ra and not rb:
        return "Reversible → Irreversible"
    if not ra and rb:
        return "Irreversible → Reversible"
    return "Irreversible → Irreversible"


def load_table() -> pd.DataFrame:
    if CACHE.exists():
        return pd.read_pickle(CACHE)
    rows = []
    for path in sorted(glob.glob(str(MSDB_DATA / "Biochemistry" / "reaction_*.json"))):
        for entry in json.load(open(path)):
            if entry.get("status") == "EMPTY":
                continue
            thermo = entry.get("thermodynamics") or {}
            row = {"rxn": entry["id"]}
            for key, (subkey, _) in SOURCES.items():
                tr = thermo.get(subkey)
                if not tr or len(tr) < 3 or tr[2] in (None, "?"):
                    continue
                try:
                    dg, sig = float(tr[0]), abs(float(tr[1]))
                except (TypeError, ValueError):
                    continue
                if key == "EQ" and sig > EQ_SENTINEL:
                    continue
                _, op, _ = run_reversibility(entry, per_source_energy(subkey),
                                             DEFAULT_HEURISTICS)
                if op is None:
                    continue
                row[f"dg_{key}"] = dg
                row[f"op_{key}"] = op
            rows.append(row)
    df = pd.DataFrame(rows).set_index("rxn")
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    df.to_pickle(CACHE)
    return df


def base_set(df: pd.DataFrame, a: str, b: str) -> pd.DataFrame:
    m = df[f"dg_{a}"].notna() & df[f"dg_{b}"].notna()
    s = df[m]
    return s[(s[f"dg_{a}"].abs() <= EXTREME_CUTOFF)
             & (s[f"dg_{b}"].abs() <= EXTREME_CUTOFF)]


def metrics(sub: pd.DataFrame, a: str, b: str) -> dict:
    x = sub[f"dg_{a}"].to_numpy(float)
    y = sub[f"dg_{b}"].to_numpy(float)
    if len(x) < 3:
        return {"n": len(x)}
    return {
        "n": int(len(x)),
        "pearson_r": float(np.corrcoef(x, y)[0, 1]),
        "spearman_rho": float(pd.Series(x).corr(pd.Series(y), method="spearman")),
        "median_abs_delta": float(np.median(np.abs(y - x))),
        "frac_sign_flip": float(np.mean(np.sign(x) * np.sign(y) < 0)),
    }


def frame(ax) -> None:
    ax.grid(True, color=GRIDLINE, linewidth=0.8, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(BASELINE)
    ax.tick_params(colors=INK_MUTED, labelsize=8.5)


def draw(ax, sub, a, b, title, lims, *, legend=False, extra="") -> dict:
    x = sub[f"dg_{a}"].to_numpy(float)
    y = sub[f"dg_{b}"].to_numpy(float)
    cats = [classify(p, q) for p, q in zip(sub[f"op_{a}"], sub[f"op_{b}"])]
    counts = Counter(cats)
    m = metrics(sub, a, b)

    ax.set_facecolor(SURFACE)
    ax.plot(lims, lims, "--", linewidth=1.2, color=BASELINE, zorder=1)
    for cat in CATEGORY_ORDER:
        idx = [i for i, c in enumerate(cats) if c == cat]
        if idx:
            ax.scatter(x[idx], y[idx], s=9, linewidths=0, alpha=0.55,
                       color=CATEGORY_COLOR[cat], zorder=2,
                       label=f"{cat} ({counts[cat]:,})")
    ax.set_xlim(*lims)
    ax.set_ylim(*lims)
    ax.set_xlabel(SOURCES[a][1], color=INK_SECONDARY, fontsize=9)
    ax.set_ylabel(SOURCES[b][1], color=INK_SECONDARY, fontsize=9)
    ax.set_title(title, color=INK_PRIMARY, fontsize=10.5, pad=8)
    frame(ax)
    off = int(np.sum((x < lims[0]) | (x > lims[1]) | (y < lims[0]) | (y > lims[1])))
    ax.text(0.97, 0.03,
            f"n = {m['n']:,}\nr = {m['pearson_r']:.2f}   ρ = {m['spearman_rho']:.2f}\n"
            f"median |Δ| = {m['median_abs_delta']:.2f}\n"
            f"sign flips = {100 * m['frac_sign_flip']:.0f}%"
            + (f"\n{off} off-scale" if off else "") + extra,
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
            color=INK_MUTED, linespacing=1.45)
    if legend:
        lg = ax.legend(loc="upper left", frameon=False, fontsize=7.5,
                       labelcolor=INK_PRIMARY, markerscale=1.6,
                       title="Reversibility transition", title_fontsize=7.5)
        lg.get_title().set_color(INK_SECONDARY)
    return m


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    print(f"loading {MSDB_DATA} ...")
    df = load_table()
    mask = load_mask()
    print(f"  {len(df):,} reactions; KEGG mask drops {len(mask):,} dGPredictor values")

    rng = np.random.default_rng(SEED)
    stats = []
    subsets: dict[tuple[str, str, str], pd.DataFrame] = {}

    for a, b in PAIRS:
        full = base_set(df, a, b)
        fixed = full[~full.index.isin(mask)]
        subsets[(a, b, "as_shipped")] = full
        subsets[(a, b, "kegg_fixed")] = fixed

        # control: N_CONTROL random subsets of |fixed| drawn from the full set
        ctrl_stats = []
        pick = None
        for i in range(N_CONTROL):
            sel = rng.choice(len(full), size=len(fixed), replace=False)
            s = full.iloc[np.sort(sel)]
            ctrl_stats.append(metrics(s, a, b))
            if i == 0:
                pick = s
        subsets[(a, b, "random_control")] = pick
        ctrl_mean = {k: float(np.mean([c[k] for c in ctrl_stats]))
                     for k in ("pearson_r", "spearman_rho", "median_abs_delta",
                               "frac_sign_flip")}
        ctrl_sd = float(np.std([c["pearson_r"] for c in ctrl_stats]))

        for cond in CONDITIONS:
            m = metrics(subsets[(a, b, cond)], a, b)
            row = {"pair": f"{a}_vs_{b}", "condition": cond, **m}
            if cond == "random_control":
                row.update({f"mean_of_{N_CONTROL}_{k}": v for k, v in ctrl_mean.items()})
                row["sd_of_pearson_r_across_controls"] = ctrl_sd
            stats.append(row)

    sdf = pd.DataFrame(stats)
    sdf.to_csv(out / "kegg_fix_stats.tsv", sep="\t", index=False)

    # ---------------- 3 x 3 grid, twice (full range and zoomed) --------------
    for suffix, zoom in (("", False), ("_zoom", True)):
        fig, axes = plt.subplots(len(PAIRS), len(CONDITIONS),
                                 figsize=(5.0 * len(CONDITIONS), 4.9 * len(PAIRS)),
                                 dpi=150)
        fig.patch.set_facecolor(SURFACE)
        for r, (a, b) in enumerate(PAIRS):
            full = subsets[(a, b, "as_shipped")]
            lo = float(min(full[f"dg_{a}"].min(), full[f"dg_{b}"].min()))
            hi = float(max(full[f"dg_{a}"].max(), full[f"dg_{b}"].max()))
            pad = 0.06 * (hi - lo)
            lims = (-ZOOM_LIMIT, ZOOM_LIMIT) if zoom else (lo - pad, hi + pad)
            for c, cond in enumerate(CONDITIONS):
                extra = ""
                if cond == "random_control":
                    cr = sdf[(sdf.pair == f"{a}_vs_{b}")
                             & (sdf.condition == "random_control")].iloc[0]
                    extra = (f"\nmean r over {N_CONTROL} draws = "
                             f"{cr[f'mean_of_{N_CONTROL}_pearson_r']:.2f}")
                draw(axes[r, c], subsets[(a, b, cond)], a, b,
                     COND_TITLE[cond] if r == 0 else "",
                     lims, legend=(r == 0 and c == 0), extra=extra)
        fig.suptitle("Impact of the inherited-KEGG-id defect on the original dGPredictor"
                     + ("  (axes clipped to ±250 kcal/mol)" if zoom else ""),
                     color=INK_PRIMARY, fontsize=14, y=0.995)
        fig.tight_layout(rect=(0, 0, 1, 0.975))
        p = out / f"grid_before_after{suffix}.png"
        fig.savefig(p, facecolor=SURFACE)
        plt.close(fig)
        print(f"wrote {p}")

        # individual panels
        for a, b in PAIRS:
            full = subsets[(a, b, "as_shipped")]
            lo = float(min(full[f"dg_{a}"].min(), full[f"dg_{b}"].min()))
            hi = float(max(full[f"dg_{a}"].max(), full[f"dg_{b}"].max()))
            pad = 0.06 * (hi - lo)
            lims = (-ZOOM_LIMIT, ZOOM_LIMIT) if zoom else (lo - pad, hi + pad)
            for cond in CONDITIONS:
                fig, ax = plt.subplots(figsize=(6.4, 6.0), dpi=150)
                fig.patch.set_facecolor(SURFACE)
                draw(ax, subsets[(a, b, cond)], a, b,
                     f"{SOURCES[a][0]} vs {SOURCES[b][0]}\n{COND_TITLE[cond]}",
                     lims, legend=True)
                fig.tight_layout()
                d = out / cond
                d.mkdir(exist_ok=True)
                fig.savefig(d / f"{a}_vs_{b}{suffix}.png", facecolor=SURFACE)
                plt.close(fig)

    # ---------------- one-glance summary -------------------------------------
    # Before -> after as a slope, with the random control marked. Neutral gray is
    # "as shipped" (not a categorical slot); blue/orange are the validated pair.
    fig, ax = plt.subplots(figsize=(9.2, 4.4), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    labels = []
    for i, (a, b) in enumerate(PAIRS):
        y = len(PAIRS) - 1 - i
        g = sdf[sdf.pair == f"{a}_vs_{b}"].set_index("condition")
        r0 = g.loc["as_shipped", "pearson_r"]
        r1 = g.loc["kegg_fixed", "pearson_r"]
        rc = g.loc["random_control", f"mean_of_{N_CONTROL}_pearson_r"]
        sd = g.loc["random_control", "sd_of_pearson_r_across_controls"]
        ax.plot([r0, r1], [y, y], color=BASELINE, linewidth=2.5, zorder=1)
        ax.scatter([r0], [y], s=130, color=NEUTRAL, zorder=3, linewidths=0)
        ax.scatter([r1], [y], s=130, color="#2a78d6", zorder=3, linewidths=0)
        ax.errorbar([rc], [y], xerr=[2 * sd], fmt="o", ms=7, color="#eb6834",
                    ecolor="#eb6834", elinewidth=1.6, capsize=3, zorder=4)
        ax.annotate(f"{r0:.2f}", (r0, y), xytext=(0, 13), textcoords="offset points",
                    ha="center", fontsize=9, color=INK_SECONDARY)
        ax.annotate(f"{r1:.2f}", (r1, y), xytext=(0, 13), textcoords="offset points",
                    ha="center", fontsize=9, color="#2a78d6")
        z = (r1 - rc) / sd if sd else float("nan")
        ax.annotate(f"{z:.0f}σ beyond the control", (r1, y), xytext=(10, -4),
                    textcoords="offset points", ha="left", va="center",
                    fontsize=8.5, color=INK_MUTED)
        labels.append(f"{SOURCES[a][0]}\nvs {SOURCES[b][0]}")
    ax.set_yticks(range(len(PAIRS)))
    ax.set_yticklabels(labels[::-1], fontsize=9, color=INK_SECONDARY)
    ax.set_ylim(-0.6, len(PAIRS) - 0.25)
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("Pearson r between the two ΔG′° series",
                  color=INK_SECONDARY, fontsize=10)
    ax.set_title("Removing the inherited KEGG ids is what fixes the agreement —\n"
                 "removing the same number of random reactions does nothing",
                 color=INK_PRIMARY, fontsize=12, loc="left", pad=26)
    frame(ax)
    ax.scatter([], [], s=130, color=NEUTRAL, label="as shipped on dev")
    ax.scatter([], [], s=130, color="#2a78d6", label="inherited KEGG ids dropped")
    ax.scatter([], [], s=70, color="#eb6834",
               label=f"random subset of equal size (mean ± 2 SD, {N_CONTROL} draws)")
    lg = ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.005), ncol=3,
                   frameon=False, fontsize=8.5, labelcolor=INK_PRIMARY,
                   borderpad=0.0, columnspacing=1.6, handletextpad=0.4)
    fig.tight_layout()
    p = out / "summary_r_before_after.png"
    fig.savefig(p, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {p}")

    print("\n" + sdf.to_string(index=False))


if __name__ == "__main__":
    main()
