"""Library of figure-builders for the thermo-source comparison pipeline.

Each ``make_fig_*`` function takes already-loaded DataFrame(s) plus an
``out_path``, writes a 150-dpi PNG via the matplotlib ``Agg`` backend
(``tight_layout``), and returns the output ``Path``.

This module is the canonical home for the matplotlib code.  The CLI
``build_thermo_source_figures.py`` and the per-figure cells in
``notebooks/10_ThermoSourceComparison.ipynb`` both import from here so
that plot semantics (colors, titles, ordering, axis scales) stay in
exactly one place.

Inputs expected on disk:

  - ``panel_fba_long.csv``     -- one row per (variant, model_id)
  - ``panel_fba_summary.csv``  -- one row per variant
  - ``coverage_<slug>.csv``    -- per-model coverage stats per source
  - ``overrides_<slug>.csv``   -- per-model bound-class transitions (wide schema)
  - ``rxn_directions_<slug>.csv`` -- per-source operator + dG columns

All ``out_path`` arguments are written under
``PROJECT_ROOT / 'reports' / 'figures' / 'thermo_sources/'`` by the
callers; this module does no path validation -- the callers own that.
"""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

# ---------------------------------------------------------------------------
# Shared style
# ---------------------------------------------------------------------------
VARIANTS = ["kbase_baseline", "gc", "eq", "dgp"]
SOURCES = ["gc", "eq", "dgp"]

VARIANT_COLORS = {
    "kbase_baseline": "#7f7f7f",  # gray
    "gc": "#1f77b4",              # blue
    "eq": "#2ca02c",              # green
    "dgp": "#d62728",             # red
}

OP_COLORS = {
    ">": "#1f77b4",
    "=": "#7f7f7f",
    "<": "#d62728",
}

OVERRIDE_KEYS = [
    "n_fwd_to_rev",
    "n_rev_to_fwd",
    "n_fwd_to_reversible",
    "n_reversible_to_fwd",
    "n_rev_to_reversible",
    "n_reversible_to_rev",
]

OVERRIDE_COLORS = {
    "n_fwd_to_rev":        "#1f77b4",
    "n_rev_to_fwd":        "#ff7f0e",
    "n_fwd_to_reversible": "#2ca02c",
    "n_reversible_to_fwd": "#d62728",
    "n_rev_to_reversible": "#9467bd",
    "n_reversible_to_rev": "#8c564b",
}

DPI = 150
FLUX_OFFSET = 1e-6  # for log scaling of zero/near-zero fluxes


def _save(fig: plt.Figure, out_path: Path) -> Path:
    """Tight-layout, savefig at the shared DPI, close, return the path."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=DPI)
    plt.close(fig)
    return out_path


# ---------------------------------------------------------------------------
# Figure 1: grower counts per variant
# ---------------------------------------------------------------------------
def make_fig_grower_counts(summary_df: pd.DataFrame, out_path: Path) -> Path:
    """Bar chart: # growers per variant (from panel_fba_long via groupby).

    Accepts either the long-form DataFrame (one row per (variant, model_id))
    or anything that has ``variant`` and ``grows`` columns -- the function
    just groups on ``variant`` and sums ``grows``.
    """
    counts = (
        summary_df.groupby("variant")["grows"].sum().reindex(VARIANTS).astype(int)
    )
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    bars = ax.bar(
        [v for v in VARIANTS],
        counts.values,
        color=[VARIANT_COLORS[v] for v in VARIANTS],
        edgecolor="black",
        linewidth=0.5,
    )
    for b, v in zip(bars, counts.values):
        ax.text(
            b.get_x() + b.get_width() / 2,
            b.get_height() + 0.5,
            f"{v}",
            ha="center",
            va="bottom",
            fontsize=10,
        )
    ax.set_ylim(0, 105)
    ax.set_ylabel("Growers (out of 100 models)")
    ax.set_xlabel("Variant")
    ax.set_title("Growers per variant (panel n=100)")
    ax.axhline(100, color="black", linestyle=":", linewidth=0.7, alpha=0.6)
    return _save(fig, out_path)


# ---------------------------------------------------------------------------
# Figure 2: mean flux per variant (over growers only)
# ---------------------------------------------------------------------------
def make_fig_mean_flux(summary_df: pd.DataFrame, out_path: Path) -> Path:
    """Bar chart: mean growth flux per variant, restricted to growers."""
    growers = summary_df[summary_df["grows"]]
    means = (
        growers.groupby("variant")["growth_flux"].mean().reindex(VARIANTS)
    )
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    bars = ax.bar(
        VARIANTS,
        means.values,
        color=[VARIANT_COLORS[v] for v in VARIANTS],
        edgecolor="black",
        linewidth=0.5,
    )
    for b, v in zip(bars, means.values):
        ax.text(
            b.get_x() + b.get_width() / 2,
            b.get_height() + 0.5,
            f"{v:.2f}",
            ha="center",
            va="bottom",
            fontsize=10,
        )
    ax.set_ylabel("Mean growth flux (growers only)")
    ax.set_xlabel("Variant")
    ax.set_title("Mean growth flux per variant")
    return _save(fig, out_path)


# ---------------------------------------------------------------------------
# Figure 3: per-variant violin + strip plot, log y
# ---------------------------------------------------------------------------
def make_fig_flux_violin(long_df: pd.DataFrame, out_path: Path) -> Path:
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    # Pre-sort by (variant in VARIANTS, model_id ascending) so the per-variant
    # value sequence is deterministic regardless of which producer wrote
    # panel_fba_long.csv. The jitter RNG (seed=0) pairs each point with a
    # specific x-offset, so input row order is load-bearing for the rendered
    # PNG.
    order_map = {v: i for i, v in enumerate(VARIANTS)}
    sorted_df = long_df.copy()
    sorted_df["__variant_order"] = sorted_df["variant"].map(
        lambda v: order_map.get(v, len(VARIANTS))
    )
    sorted_df = sorted_df.sort_values(
        ["__variant_order", "variant", "model_id"],
        kind="mergesort",
    )
    data = []
    for v in VARIANTS:
        vals = (
            sorted_df.loc[sorted_df["variant"] == v, "growth_flux"]
            .fillna(0.0)
            .values
        )
        data.append(vals + FLUX_OFFSET)
    positions = np.arange(len(VARIANTS))
    parts = ax.violinplot(
        data,
        positions=positions,
        widths=0.8,
        showmeans=False,
        showmedians=True,
        showextrema=False,
    )
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor(VARIANT_COLORS[VARIANTS[i]])
        body.set_edgecolor("black")
        body.set_alpha(0.5)
    if "cmedians" in parts:
        parts["cmedians"].set_color("black")

    rng = np.random.default_rng(0)
    for i, vals in enumerate(data):
        jitter = rng.uniform(-0.12, 0.12, size=len(vals))
        ax.scatter(
            positions[i] + jitter,
            vals,
            s=8,
            color=VARIANT_COLORS[VARIANTS[i]],
            edgecolor="black",
            linewidth=0.2,
            alpha=0.7,
        )
    ax.set_yscale("log")
    ax.set_xticks(positions)
    ax.set_xticklabels(VARIANTS)
    ax.set_ylabel(f"growth_flux + {FLUX_OFFSET:g} (log scale)")
    ax.set_xlabel("Variant")
    ax.set_title("Per-model growth flux distribution")
    ax.grid(True, axis="y", which="both", linestyle=":", alpha=0.4)
    return _save(fig, out_path)


# ---------------------------------------------------------------------------
# Figure 4: per-model heatmap (100 rows x 4 variants), viridis log-scaled
# ---------------------------------------------------------------------------
def make_fig_per_model_heatmap(long_df: pd.DataFrame, out_path: Path) -> Path:
    pivot = long_df.pivot(index="model_id", columns="variant", values="growth_flux")
    pivot = pivot.reindex(columns=VARIANTS)
    pivot = pivot.sort_values("kbase_baseline", ascending=False)

    arr = pivot.fillna(0.0).values + FLUX_OFFSET
    vmin = max(arr.min(), FLUX_OFFSET)
    vmax = arr.max()

    fig, ax = plt.subplots(figsize=(5.5, 12.0))
    im = ax.imshow(
        arr,
        aspect="auto",
        cmap="viridis",
        norm=LogNorm(vmin=vmin, vmax=vmax),
        interpolation="nearest",
    )
    ax.set_xticks(np.arange(len(VARIANTS)))
    ax.set_xticklabels(VARIANTS, rotation=30, ha="right")
    ax.set_yticks([])
    ax.set_ylabel("Models (sorted by baseline growth flux, descending)")
    ax.set_title("Per-model growth flux across variants")
    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label(f"growth_flux + {FLUX_OFFSET:g} (log)")
    return _save(fig, out_path)


# ---------------------------------------------------------------------------
# Figure 5: coverage histograms (1x3 panel)
# ---------------------------------------------------------------------------
def make_fig_coverage_per_source(
    coverage_dfs: Mapping[str, pd.DataFrame],
    out_path: Path,
) -> Path:
    """1x3 histogram of per-model ``frac_covered`` for each source."""
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2), sharey=True)
    bins = np.linspace(0.0, 1.0, 21)
    for ax, src in zip(axes, SOURCES):
        df = coverage_dfs[src]
        vals = df["frac_covered"].dropna().values
        ax.hist(
            vals,
            bins=bins,
            color=VARIANT_COLORS[src],
            edgecolor="black",
            alpha=0.85,
        )
        mean = vals.mean() if len(vals) else float("nan")
        ax.axvline(mean, color="black", linestyle="--", linewidth=1.0)
        ax.set_title(
            f"Per-model fraction of reactions covered by source {src}\n"
            f"mean={mean:.3f}"
        )
        ax.set_xlabel("frac_covered")
        ax.set_xlim(0.0, 1.0)
        ax.grid(True, axis="y", linestyle=":", alpha=0.4)
    axes[0].set_ylabel("# models")
    return _save(fig, out_path)


# ---------------------------------------------------------------------------
# Figure 6: aggregate override transitions (1x3 stacked bars)
# ---------------------------------------------------------------------------
def make_fig_override_transitions(
    overrides_dfs: Mapping[str, pd.DataFrame],
    out_path: Path,
) -> Path:
    """1x3 stacked bars showing total transitions of each kind per source.

    Expects the wide schema produced by direction_pipeline.override_transitions
    (and the standalone CLI): one row per model with the six
    ``n_<before>_to_<after>`` columns.
    """
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.6), sharey=True)
    for ax, src in zip(axes, SOURCES):
        df = overrides_dfs[src]
        totals = df[OVERRIDE_KEYS].sum()
        bottom = 0.0
        for key in OVERRIDE_KEYS:
            val = float(totals[key])
            ax.bar(
                [src],
                [val],
                bottom=[bottom],
                color=OVERRIDE_COLORS[key],
                edgecolor="black",
                linewidth=0.4,
                label=key if ax is axes[0] else None,
            )
            if val > 0:
                ax.text(
                    0,
                    bottom + val / 2,
                    f"{int(val)}",
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="white" if val > totals.sum() * 0.04 else "black",
                )
            bottom += val
        ax.set_title(f"source = {src}\ntotal transitions = {int(totals.sum())}")
        if ax is axes[0]:
            ax.set_ylabel("# overrides (summed across panel models)")
    axes[0].legend(loc="upper left", fontsize=8, bbox_to_anchor=(1.02, 1.0))
    fig.suptitle("Aggregate override transitions per source (n=100 models)")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI)
    plt.close(fig)
    return out_path


# ---------------------------------------------------------------------------
# Figure 7: scatter source flux vs baseline flux per model (1x3)
# ---------------------------------------------------------------------------
def make_fig_flux_vs_baseline_scatter(
    long_df: pd.DataFrame, out_path: Path
) -> Path:
    pivot = long_df.pivot(index="model_id", columns="variant", values="growth_flux")
    pivot = pivot.reindex(columns=VARIANTS).fillna(0.0)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), sharex=True, sharey=True)
    baseline = pivot["kbase_baseline"].values
    eps = 1e-9
    for ax, src in zip(axes, SOURCES):
        y = pivot[src].values
        delta = y - baseline
        above = int((delta > eps).sum())
        below = int((delta < -eps).sum())
        equal = int(np.isclose(delta, 0.0, atol=eps).sum())
        mean_abs_delta = float(np.mean(np.abs(delta)))

        colors = np.where(
            delta > eps, "#1a9850",
            np.where(delta < -eps, "#d73027", "#7f7f7f"),
        )
        ax.scatter(
            baseline,
            y,
            c=colors,
            edgecolor="black",
            linewidth=0.3,
            s=28,
            alpha=0.85,
        )
        lo = 0.0
        hi = max(baseline.max(), y.max()) * 1.05
        ax.plot([lo, hi], [lo, hi], "k--", linewidth=1.0, alpha=0.7)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_xlabel("baseline growth_flux")
        ax.set_ylabel("source growth_flux")
        ax.set_title(
            f"{src}\nabove={above}, below={below}, equal={equal}\n"
            f"mean|delta|={mean_abs_delta:.3f}"
        )
        ax.grid(True, linestyle=":", alpha=0.4)
    fig.suptitle("Per-model growth flux: source vs baseline")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI)
    plt.close(fig)
    return out_path


# ---------------------------------------------------------------------------
# Figure 8: dG distribution per source, colored by operator
# ---------------------------------------------------------------------------
def make_fig_dg_distribution_per_source(
    direction_dfs: Mapping[str, pd.DataFrame],
    out_path: Path,
) -> Path:
    """1x3 stacked-histogram of dG by operator per source.

    Expects DataFrames with ``operator`` and ``dg`` columns (produced by
    direction_pipeline.snapshot_msdb_per_source).
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), sharey=True)
    # shared bins across panels for visual comparability
    all_dg = np.concatenate(
        [
            direction_dfs[s]["dg"].dropna().values
            for s in SOURCES
            if "dg" in direction_dfs[s].columns
        ]
    )
    lo = np.percentile(all_dg, 0.5)
    hi = np.percentile(all_dg, 99.5)
    bins = np.linspace(lo, hi, 61)

    for ax, src in zip(axes, SOURCES):
        df = direction_dfs[src].dropna(subset=["dg"])
        operators = ["<", "=", ">"]
        stacked = [df.loc[df["operator"] == op, "dg"].values for op in operators]
        colors = [OP_COLORS[op] for op in operators]
        ax.hist(
            stacked,
            bins=bins,
            stacked=True,
            color=colors,
            label=operators,
            edgecolor="black",
            linewidth=0.2,
        )
        ax.set_yscale("log")
        ax.set_xlabel("dG (kJ/mol)")
        ax.set_title(
            f"source={src}  n_rxns={len(df)}\n"
            f">: {(df['operator']=='>').sum()}, =: {(df['operator']=='=').sum()}, "
            f"<: {(df['operator']=='<').sum()}"
        )
        ax.axvline(0.0, color="black", linewidth=0.7, linestyle=":")
        ax.grid(True, axis="y", which="both", linestyle=":", alpha=0.4)
        ax.legend(title="operator", fontsize=8, loc="upper right")
    axes[0].set_ylabel("# reactions (log)")
    fig.suptitle("Per-source dG distribution by direction operator")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI)
    plt.close(fig)
    return out_path
