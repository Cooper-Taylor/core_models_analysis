#!/usr/bin/env python
"""Build matplotlib figures from variant FBA results."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

RESULTS_DIR = Path("/scratch/ctaylor/core_models_analysis/results")
OUT_DIR = Path("/scratch/ctaylor/core_models_analysis/reports/figures")
OUT_DIR.mkdir(parents=True, exist_ok=True)

PANEL_CSV = RESULTS_DIR / "variant_panel_fba.csv"
SUMMARY_CSV = RESULTS_DIR / "variant_panel_summary.csv"
DIFF_CSV = RESULTS_DIR / "rev_diff_dev_vs_claude.csv"
DIFF_SUMMARY_JSON = RESULTS_DIR / "rev_diff_summary.json"

# Canonical variant order
VARIANT_ORDER = [
    "on_disk",
    "msdb_dev",
    "cascade_live",
    "msdb_claude",
    "branch_diff_only",
    "all_reversible",
    "all_forward",
    "flip_eq_to_gt",
]


def variant_colors(variants):
    """on_disk -> gray, others -> distinct tab10 colors."""
    palette = plt.get_cmap("tab10").colors
    colors = {}
    j = 0
    for v in variants:
        if v == "on_disk":
            colors[v] = "#888888"
        else:
            colors[v] = palette[j % len(palette)]
            j += 1
    return colors


def order_variants(variants):
    present = [v for v in VARIANT_ORDER if v in variants]
    extras = [v for v in variants if v not in present]
    return present + extras


def fig1_grower_counts(panel: pd.DataFrame, summary: pd.DataFrame, path: Path):
    summary = summary.set_index("variant")
    variants = order_variants(summary.index.tolist())
    counts = [int(summary.loc[v, "n_growers"]) for v in variants]
    colors = variant_colors(variants)
    bar_colors = [colors[v] for v in variants]

    fig, ax = plt.subplots(figsize=(8, 5))
    y_pos = np.arange(len(variants))
    bars = ax.barh(y_pos, counts, color=bar_colors, edgecolor="black", linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(variants)
    ax.invert_yaxis()
    ax.set_xlabel("Number of growers (out of 100)")
    ax.set_title("Growers in the 100-model panel by reaction-direction source")
    ax.set_xlim(0, 105)
    for bar, c in zip(bars, counts):
        ax.text(
            bar.get_width() + 1,
            bar.get_y() + bar.get_height() / 2,
            str(c),
            va="center",
            fontsize=9,
        )
    ax.grid(axis="x", linestyle=":", alpha=0.5)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig2_mean_flux(panel: pd.DataFrame, summary: pd.DataFrame, path: Path):
    # Mean flux over growers only (computed from panel to be safe)
    growers = panel[panel["grows"] == True]  # noqa: E712
    grouped = (
        growers.groupby("variant")["growth_flux"]
        .mean()
        .reindex(summary["variant"].tolist())
        .fillna(0.0)
    )
    variants = order_variants(grouped.index.tolist())
    values = [grouped[v] for v in variants]
    colors = variant_colors(variants)
    bar_colors = [colors[v] for v in variants]

    fig, ax = plt.subplots(figsize=(9, 5))
    x_pos = np.arange(len(variants))
    bars = ax.bar(x_pos, values, color=bar_colors, edgecolor="black", linewidth=0.5)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(variants, rotation=30, ha="right")
    ax.set_ylabel("Mean growth flux (growers only)")
    ax.set_title("Mean growth flux per variant (averaged over growers)")
    for bar, v in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(values) * 0.01,
            f"{v:.1f}",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    ax.grid(axis="y", linestyle=":", alpha=0.5)
    ax.set_ylim(0, max(values) * 1.15 if max(values) > 0 else 1)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig3_flux_distribution(panel: pd.DataFrame, path: Path):
    variants = order_variants(sorted(panel["variant"].unique()))
    offset = 1e-6
    data = []
    for v in variants:
        vals = panel.loc[panel["variant"] == v, "growth_flux"].to_numpy()
        vals = np.where(vals < 0, 0.0, vals) + offset
        data.append(vals)
    colors = variant_colors(variants)

    fig, ax = plt.subplots(figsize=(11, 6))
    positions = np.arange(1, len(variants) + 1)

    # Violins (handle zero-variance variants gracefully)
    for pos, vals, v in zip(positions, data, variants):
        if np.ptp(vals) == 0:
            # all equal — just plot a horizontal tick
            ax.hlines(vals[0], pos - 0.25, pos + 0.25, color=colors[v], lw=2)
        else:
            parts = ax.violinplot(
                [vals], positions=[pos], widths=0.7, showmeans=False, showmedians=True
            )
            for body in parts["bodies"]:
                body.set_facecolor(colors[v])
                body.set_edgecolor("black")
                body.set_alpha(0.6)
            for key in ("cmaxes", "cmins", "cbars", "cmedians"):
                if key in parts:
                    parts[key].set_color("black")
                    parts[key].set_linewidth(0.8)
        # strip overlay
        jitter = (np.random.RandomState(0).rand(len(vals)) - 0.5) * 0.25
        ax.scatter(
            np.full_like(vals, pos) + jitter,
            vals,
            s=8,
            color=colors[v],
            alpha=0.5,
            edgecolor="black",
            linewidth=0.2,
        )

    ax.set_yscale("log")
    ax.set_xticks(positions)
    ax.set_xticklabels(variants, rotation=30, ha="right")
    ax.set_ylabel("Growth flux (log scale, +1e-6 offset)")
    ax.set_title("Growth-flux distribution across the 100-model panel")
    ax.grid(axis="y", linestyle=":", alpha=0.4, which="both")
    caption = (
        "Caption: 'all_reversible' is the upper bound: every internal reaction "
        "is opened to (-1000, +1000), so flux is unconstrained by direction. "
        "Variants with all values at the 1e-6 floor are non-growers."
    )
    fig.text(0.5, -0.02, caption, ha="center", va="top", fontsize=9, wrap=True)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig4_grow_change(panel: pd.DataFrame, path: Path):
    # 2x2 confusion: rows on_disk grows, cols msdb_claude grows
    pivot = panel.pivot(index="model_id", columns="variant", values="grows")
    if "on_disk" not in pivot.columns or "msdb_claude" not in pivot.columns:
        raise SystemExit("Required variants missing for fig4")
    a = pivot["on_disk"].astype(bool)
    b = pivot["msdb_claude"].astype(bool)
    mat = np.zeros((2, 2), dtype=int)
    for av, bv in zip(a, b):
        mat[0 if av else 1, 0 if bv else 1] += 1
    # rows ordered True, False (consistent with labels below)
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(mat, cmap="Blues")
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(["msdb_claude grows", "msdb_claude no-grow"])
    ax.set_yticklabels(["on_disk grows", "on_disk no-grow"])
    ax.set_title("Growth change: on_disk -> msdb_claude (n=100)")
    for i in range(2):
        for j in range(2):
            color = "white" if mat[i, j] > mat.max() / 2 else "black"
            ax.text(
                j,
                i,
                str(mat[i, j]),
                ha="center",
                va="center",
                fontsize=20,
                color=color,
                fontweight="bold",
            )
    fig.colorbar(im, ax=ax, fraction=0.045, pad=0.04, label="model count")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig5_branch_diff_transitions(path: Path):
    # Read diff summary JSON for the by_transition counts (CSV may be empty)
    with open(DIFF_SUMMARY_JSON) as fh:
        summary = json.load(fh)
    transitions = summary.get("by_transition", {}) or {}

    fig, ax = plt.subplots(figsize=(8, 5))
    if not transitions:
        ax.text(
            0.5,
            0.5,
            "No reversibility transitions found:\n"
            f"dev and claude branches agree on all "
            f"{summary.get('n_total_reactions_dev', 0):,} reactions "
            f"(n_changed = {summary.get('n_changed', 0)}).",
            ha="center",
            va="center",
            fontsize=12,
            transform=ax.transAxes,
            bbox=dict(boxstyle="round", facecolor="#f5f5f5", edgecolor="#aaa"),
        )
        ax.set_axis_off()
        ax.set_title("Reversibility transitions: dev -> claude-changes")
    else:
        items = sorted(transitions.items(), key=lambda kv: kv[1], reverse=True)[:8]
        labels = [k for k, _ in items]
        counts = [v for _, v in items]
        y_pos = np.arange(len(labels))
        ax.barh(y_pos, counts, color="#4C78A8", edgecolor="black", linewidth=0.5)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels)
        ax.invert_yaxis()
        ax.set_xlabel("Number of reactions")
        ax.set_title("Reversibility transitions: dev -> claude-changes (top 8)")
        for i, c in enumerate(counts):
            ax.text(c + max(counts) * 0.01, i, str(c), va="center", fontsize=9)
        ax.grid(axis="x", linestyle=":", alpha=0.5)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig6_per_model_heatmap(panel: pd.DataFrame, path: Path):
    pivot = panel.pivot(index="model_id", columns="variant", values="growth_flux")
    variants = order_variants(pivot.columns.tolist())
    pivot = pivot[variants]
    # sort rows by on_disk flux
    if "on_disk" in pivot.columns:
        pivot = pivot.sort_values("on_disk", ascending=False)
    data = pivot.to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(8, 14))
    im = ax.imshow(data, aspect="auto", cmap="viridis")
    ax.set_xticks(np.arange(len(variants)))
    ax.set_xticklabels(variants, rotation=30, ha="right")
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=5)
    ax.set_title("Per-model growth flux across variants\n(rows sorted by on_disk flux)")
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Growth flux")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig7_cascade_live_vs_msdb_dev_scatter(panel: pd.DataFrame, path: Path):
    """Per-model scatter: x = msdb_dev growth_flux, y = cascade_live growth_flux."""
    pivot = panel.pivot(index="model_id", columns="variant", values="growth_flux")
    missing = [v for v in ("msdb_dev", "cascade_live") if v not in pivot.columns]
    if missing:
        raise SystemExit(f"Required variants missing for fig7 scatter: {missing}")
    sub = pivot[["msdb_dev", "cascade_live"]].dropna()
    x = sub["msdb_dev"].to_numpy(dtype=float)
    y = sub["cascade_live"].to_numpy(dtype=float)
    n = len(sub)

    eps = 1e-9
    n_above = int(np.sum(y - x > eps))
    n_below = int(np.sum(x - y > eps))
    n_on = int(np.sum(np.abs(y - x) <= eps))
    any_diff = (n_above + n_below) > 0

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(x, y, s=24, color="#4C78A8", edgecolor="black",
               linewidth=0.3, alpha=0.75)

    lo = float(min(x.min(), y.min(), 0.0))
    hi = float(max(x.max(), y.max(), 1.0))
    pad = 0.02 * (hi - lo) if hi > lo else 1.0
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad],
            color="#888", linestyle="--", linewidth=1.0, label="1:1")
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(lo - pad, hi + pad)

    ax.set_xlabel("msdb_dev growth flux")
    ax.set_ylabel("cascade_live growth flux")
    title_suffix = "differs" if any_diff else "matches"
    ax.set_title(
        f"cascade_live vs msdb_dev growth flux (n={n})\n"
        f"cascade_live {title_suffix} msdb_dev for at least one panel model"
    )

    txt = (
        f"above diag (cascade > dev): {n_above}\n"
        f"below diag (cascade < dev): {n_below}\n"
        f"on  diag (cascade = dev): {n_on}"
    )
    ax.text(
        0.03, 0.97, txt,
        transform=ax.transAxes, ha="left", va="top", fontsize=10,
        bbox=dict(boxstyle="round", facecolor="white", edgecolor="#888", alpha=0.85),
    )
    ax.grid(linestyle=":", alpha=0.5)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def verify(path: Path) -> None:
    if not path.exists():
        raise SystemExit(f"Missing output: {path}")
    if path.stat().st_size == 0:
        raise SystemExit(f"Empty output: {path}")


def main():
    panel = pd.read_csv(PANEL_CSV)
    summary = pd.read_csv(SUMMARY_CSV)

    outputs = [
        ("fig_grower_counts_by_variant.png", lambda p: fig1_grower_counts(panel, summary, p)),
        ("fig_mean_flux_by_variant.png", lambda p: fig2_mean_flux(panel, summary, p)),
        ("fig_flux_distribution_violin.png", lambda p: fig3_flux_distribution(panel, p)),
        ("fig_grow_change_msdb_dev_vs_claude.png", lambda p: fig4_grow_change(panel, p)),
        ("fig_branch_diff_transitions.png", lambda p: fig5_branch_diff_transitions(p)),
        ("fig_per_model_growth_heatmap.png", lambda p: fig6_per_model_heatmap(panel, p)),
        ("fig_cascade_live_vs_msdb_dev_scatter.png",
         lambda p: fig7_cascade_live_vs_msdb_dev_scatter(panel, p)),
    ]
    written = []
    for name, fn in outputs:
        out_path = OUT_DIR / name
        fn(out_path)
        verify(out_path)
        written.append(str(out_path))
        print(f"wrote {out_path} ({out_path.stat().st_size} bytes)")
    return written


if __name__ == "__main__":
    main()
