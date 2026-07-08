#!/usr/bin/env python
"""Build presentation-grade interactive figures (Plotly HTML + static PNG).

Six figures that tell the story "reaction-reversibility / thermodynamics ->
core-model ATP (growth) production":

  1. Impact tornado          - net change in panel growth flux per variant
  2. Leverage scatter        - reactions changed vs models that flip grow/no-grow
  3. Variant x model heatmap - per-model delta growth flux across variants
  4. Monte-Carlo band        - DG'-uncertainty propagated to growth flux + P(grows)
  5. P(direction) confidence - analytic P(forward) vs the cascade's deterministic call
  6. Cascade Sankey          - how 56k reactions flow heuristic -> direction

The pipeline only *reads* existing artifacts produced by earlier stages
(variant diff JSONs, the statistical panel, the live cascade) - it never
recomputes FBA, so it is fast and deterministic. It does not touch
ModelSEEDDatabase or core_models_kegg2.

Outputs
  reports/presentation/index.html         scrollable dashboard (all 6)
  reports/presentation/figN_*.html        one interactive figure each
  reports/presentation/png/figN_*.png     static export (needs kaleido+chrome)
  site/data/figures/figN_*.html           copies for embedding in the site
  site/data/figures/manifest.json         figure index for the site
  site/figures.html                       standalone viewer over the manifest

Usage
  python3 scripts/build_presentation_figures.py [--out DIR] [--figures 1,3,6|all]
                                                [--no-png] [--no-site] [--open]
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

BASE = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR", "/scratch/ctaylor/core_models_analysis"))
RESULTS = BASE / "results"
STATS = RESULTS / "statistical_panel"
VARIANTS_DIR = BASE / "site" / "data" / "variants"
CASCADE_CSV = RESULTS / "rxn_directions_cascade_live.csv"

DEFAULT_OUT = BASE / "reports" / "presentation"
SITE_FIG_DIR = BASE / "site" / "data" / "figures"
SITE_PAGE = BASE / "site" / "figures.html"

# Logical left-to-right / top-to-bottom variant order for the §3.x panel.
VARIANT_ORDER = [
    "3.1", "3.3", "3.3_wide", "3.5", "3.5_wide",
    "3.6", "3.7", "3.10_tight", "3.10_loose", "H4", "ai_opus48",
]
# Variants with a Monte-Carlo statistical panel (file stems).
STAT_VARIANTS = ["baseline", "3.5", "H4", "pforward_50", "pforward_95"]

EPS = 1e-6

PALETTE = [
    "#4C78A8", "#F58518", "#54A24B", "#E45756", "#72B7B2",
    "#EECA3B", "#B279A2", "#FF9DA6", "#9D755D", "#BAB0AC", "#8C6BB1",
]


# --------------------------------------------------------------------------- #
# Loaders
# --------------------------------------------------------------------------- #
def order_key(tag: str) -> int:
    return VARIANT_ORDER.index(tag) if tag in VARIANT_ORDER else len(VARIANT_ORDER)


def load_variant_jsons() -> list[dict]:
    """Load every variant diff JSON, ordered by VARIANT_ORDER."""
    items = []
    for f in glob.glob(str(VARIANTS_DIR / "*.json")):
        with open(f) as fh:
            items.append(json.load(fh))
    items.sort(key=lambda d: order_key(d.get("tag", "")))
    return items


def variant_panel_long(variants: list[dict]) -> pd.DataFrame:
    """One row per (variant, model) from each variant's embedded panel_fba."""
    rows = []
    for d in variants:
        tag = d.get("tag")
        for r in d.get("panel_fba", []) or []:
            rows.append(
                dict(
                    variant=tag,
                    model_id=r.get("model_id"),
                    baseline_flux=r.get("baseline_flux"),
                    variant_flux=r.get("variant_flux"),
                    delta_flux=r.get("delta_flux"),
                    baseline_grows=r.get("baseline_grows"),
                    variant_grows=r.get("variant_grows"),
                )
            )
    return pd.DataFrame(rows)


def variant_summary(variants: list[dict]) -> pd.DataFrame:
    """Per-variant aggregate metrics for the tornado + leverage figures."""
    rows = []
    for d in variants:
        dl = np.array([r.get("delta_flux", 0.0) for r in d.get("panel_fba", []) or []], float)
        n = len(dl)
        rows.append(
            dict(
                tag=d.get("tag"),
                title=d.get("title", ""),
                section=d.get("section", ""),
                n_changed=int(d.get("n_changed", 0)),
                n_models_flip=int(d.get("n_models_flip", 0)),
                n_models_flux_change=int(d.get("n_models_flux_change", 0)),
                panel_n=n,
                mean_delta=float(dl.mean()) if n else 0.0,
                median_delta=float(np.median(dl)) if n else 0.0,
                mean_abs_delta=float(np.abs(dl).mean()) if n else 0.0,
                n_up=int((dl > EPS).sum()),
                n_down=int((dl < -EPS).sum()),
            )
        )
    return pd.DataFrame(rows)


def load_stat_panel() -> tuple[dict[str, pd.DataFrame], pd.DataFrame | None]:
    """Per-variant MC distribution frames + the summary table (mean_p_grows)."""
    frames = {}
    for v in STAT_VARIANTS:
        f = STATS / f"panel_distribution__{v}__N50.csv"
        if f.exists():
            frames[v] = pd.read_csv(f)
    summ = STATS / "summary.csv"
    summary = pd.read_csv(summ) if summ.exists() else None
    return frames, summary


def load_pdirection() -> pd.DataFrame | None:
    """Analytic P(direction) joined to the cascade's deterministic call."""
    pf = STATS / "p_direction__baseline.csv"
    if not (pf.exists() and CASCADE_CSV.exists()):
        return None
    pdir = pd.read_csv(pf)
    casc = pd.read_csv(CASCADE_CSV)[["rxn_id", "reversibility"]]
    return pdir.merge(casc, on="rxn_id", how="left")


def load_cascade_flows() -> pd.DataFrame | None:
    """Cascade rows tagged with the heuristic stage that resolved them."""
    if not CASCADE_CSV.exists():
        return None
    casc = pd.read_csv(CASCADE_CSV)

    def stage(s: str) -> str:
        m = re.match(r"([A-Za-z]+)", str(s).strip())
        return m.group(1) if m else "other"

    casc["stage"] = casc["status"].map(stage)
    return casc


# Pretty labels.
STAGE_LABELS = {
    "MdeltaG": "ΔG bounds (H1)",
    "ATPS": "ATP synthase (H2)",
    "mMdeltaG": "mM ΔG band (H4)",
    "lowE": "low-energy cpds (H5)",
    "default": "default → reversible",
    "Incomplete": "incomplete thermo",
    "Empty": "empty reaction",
}
STAGE_ORDER = ["MdeltaG", "ATPS", "mMdeltaG", "lowE", "default", "Incomplete", "Empty"]
CLASS_LABELS = {">": "forward (>)", "<": "reverse (<)", "=": "reversible (=)", "?": "unknown (?)"}
CLASS_COLORS = {">": "#54A24B", "<": "#E45756", "=": "#4C78A8", "?": "#BAB0AC"}


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #
def make_fig1_impact_tornado(summary: pd.DataFrame) -> go.Figure:
    s = summary.sort_values("mean_delta")
    colors = ["#E45756" if v < 0 else "#4C78A8" for v in s["mean_delta"]]
    customdata = np.stack(
        [s["title"], s["section"], s["median_delta"], s["mean_abs_delta"],
         s["n_changed"], s["n_models_flip"], s["n_models_flux_change"],
         s["n_up"], s["n_down"], s["panel_n"]],
        axis=-1,
    )
    fig = go.Figure(
        go.Bar(
            x=s["mean_delta"], y=s["tag"], orientation="h",
            marker=dict(color=colors, line=dict(color="#333", width=0.5)),
            customdata=customdata,
            hovertemplate=(
                "<b>%{y}</b> — %{customdata[0]}<br>"
                "section %{customdata[1]}<br>"
                "mean Δflux: %{x:.2f} | median: %{customdata[2]:.2f} | "
                "mean|Δ|: %{customdata[3]:.2f}<br>"
                "reactions changed: %{customdata[4]:,}<br>"
                "models flux↑: %{customdata[7]} | flux↓: %{customdata[8]} "
                "| grow-flips: %{customdata[5]} | panel: %{customdata[9]}"
                "<extra></extra>"
            ),
        )
    )
    fig.add_vline(x=0, line_color="#444", line_width=1)
    # Fold n_changed into the y tick labels so it travels with each variant
    # and can never collide with the bars/axis (tiny negative bars otherwise do).
    ticktext = [f"{t}  ({n:,} rxns)" for t, n in zip(s["tag"], s["n_changed"])]
    fig.update_yaxes(tickmode="array", tickvals=list(s["tag"]), ticktext=ticktext)
    fig.update_layout(
        title="Fig 1 — Net change in panel growth flux per thermodynamic variant",
        xaxis_title="Mean Δ growth flux vs baseline (mmol·gDW⁻¹·h⁻¹), 100-model panel",
        yaxis_title="variant  (reactions changed)",
        template="plotly_white", height=520, margin=dict(l=170, r=50, t=70, b=60),
    )
    return fig


def make_fig2_leverage_scatter(summary: pd.DataFrame) -> go.Figure:
    s = summary.copy()
    sizeref = 2.0 * max(s["n_models_flux_change"].max(), 1) / (38.0 ** 2)
    fig = go.Figure()
    for i, (_, r) in enumerate(s.iterrows()):
        fig.add_trace(
            go.Scatter(
                x=[max(r["n_changed"], 1)], y=[r["n_models_flip"]],
                mode="markers+text",
                marker=dict(
                    size=max(r["n_models_flux_change"], 4), sizemode="area",
                    sizeref=sizeref, sizemin=5,
                    color=PALETTE[i % len(PALETTE)],
                    line=dict(color="#333", width=0.6), opacity=0.85,
                ),
                text=[r["tag"]], textposition="top center", textfont=dict(size=10),
                name=r["tag"],
                customdata=[[r["title"], r["section"], r["n_changed"],
                            r["n_models_flip"], r["n_models_flux_change"]]],
                hovertemplate=(
                    "<b>%{customdata[0]}</b> (%{customdata[1]})<br>"
                    "reactions changed: %{customdata[2]:,}<br>"
                    "models grow-flipped: %{customdata[3]}<br>"
                    "models flux-changed: %{customdata[4]}<extra></extra>"
                ),
            )
        )
    fig.update_layout(
        title="Fig 2 — Leverage: reactions changed vs biological impact",
        xaxis=dict(title="reactions changed vs baseline (log scale)", type="log"),
        yaxis=dict(title="models that flip grow ↔ no-grow"),
        template="plotly_white", height=560, showlegend=False,
        margin=dict(l=70, r=40, t=70, b=60),
        annotations=[dict(
            x=0.02, y=0.98, xref="paper", yref="paper", align="left", showarrow=False,
            text="bubble size = # models with any flux change<br>"
                 "top-left = high leverage (few edits, big effect)",
            font=dict(size=11, color="#555"),
            bgcolor="rgba(255,255,255,0.7)", bordercolor="#ccc", borderwidth=1,
        )],
    )
    return fig


def make_fig3_variant_model_heatmap(panel: pd.DataFrame) -> go.Figure:
    delta = panel.pivot_table(index="model_id", columns="variant", values="delta_flux")
    cols = [v for v in VARIANT_ORDER if v in delta.columns]
    delta = delta[cols]
    # order rows by total absolute movement (most active models on top)
    delta = delta.reindex(delta.abs().sum(axis=1).sort_values(ascending=False).index)

    raw = delta.to_numpy(dtype=float)
    # signed-log compresses the bimodal 0..~1000 range for color while keeping sign
    z = np.sign(raw) * np.log1p(np.abs(raw))
    zmax = float(np.nanmax(np.abs(z))) or 1.0
    ticks = [-1000, -100, -10, 0, 10, 100, 1000]
    tickz = [np.sign(t) * np.log1p(abs(t)) for t in ticks]

    fig = go.Figure(
        go.Heatmap(
            z=z, x=delta.columns.tolist(), y=delta.index.tolist(),
            customdata=raw,
            colorscale="RdBu", reversescale=True, zmid=0, zmin=-zmax, zmax=zmax,
            colorbar=dict(title="Δ flux", tickvals=tickz,
                          ticktext=[f"{t:+d}" for t in ticks]),
            hovertemplate=("model %{y}<br>variant %{x}<br>"
                           "Δ growth flux: %{customdata:.3f}<extra></extra>"),
        )
    )
    fig.update_layout(
        title="Fig 3 — Per-model Δ growth flux across variants "
              "(rows = panel models, signed-log color)",
        xaxis_title="variant", yaxis_title="panel model",
        template="plotly_white", height=900,
        margin=dict(l=140, r=40, t=70, b=60),
    )
    fig.update_yaxes(showticklabels=False)
    return fig


def make_fig4_mc_uncertainty(frames: dict[str, pd.DataFrame],
                             summary: pd.DataFrame | None) -> go.Figure:
    variants = [v for v in STAT_VARIANTS if v in frames]
    pg_by_variant = {}
    if summary is not None:
        for _, r in summary.iterrows():
            pg_by_variant[str(r["variant"])] = r.get("mean_p_grows")

    titles = []
    for v in variants:
        mp = pg_by_variant.get(v)
        titles.append(f"{v}   (mean P(grows) = {mp:.2f})" if mp is not None else v)

    fig = make_subplots(rows=len(variants), cols=1, subplot_titles=titles,
                        vertical_spacing=0.06)
    for i, v in enumerate(variants, start=1):
        df = frames[v].sort_values("q50").reset_index(drop=True)
        x = np.arange(1, len(df) + 1)
        fig.add_trace(
            go.Scatter(
                x=x, y=df["q50"], mode="markers",
                marker=dict(size=5, color=df["p_grows"], coloraxis="coloraxis"),
                error_y=dict(
                    type="data", symmetric=False,
                    array=(df["q95"] - df["q50"]).clip(lower=0),
                    arrayminus=(df["q50"] - df["q05"]).clip(lower=0),
                    thickness=0.8, width=0, color="rgba(120,120,120,0.5)",
                ),
                customdata=np.stack([df["model_id"], df["q05"], df["q95"], df["p_grows"]], -1),
                hovertemplate=("%{customdata[0]}<br>median flux %{y:.2f} "
                               "[%{customdata[1]:.2f}, %{customdata[2]:.2f}]<br>"
                               "P(grows) %{customdata[3]:.2f}<extra></extra>"),
                showlegend=False,
            ),
            row=i, col=1,
        )
        fig.update_yaxes(title_text="growth flux", row=i, col=1)
    fig.update_xaxes(title_text="panel model (sorted by median flux)", row=len(variants), col=1)
    fig.update_layout(
        title="Fig 4 — Growth flux under ΔG′° uncertainty (Monte-Carlo, N=50): "
              "90% CI band, colored by P(grows)",
        template="plotly_white", height=260 * len(variants),
        coloraxis=dict(colorscale="Viridis", cmin=0, cmax=1,
                       colorbar=dict(title="P(grows)")),
        margin=dict(l=70, r=40, t=80, b=60),
    )
    return fig


def make_fig5_pdirection(pdir: pd.DataFrame) -> tuple[go.Figure, str]:
    df = pdir.dropna(subset=["reversibility"]).copy()
    classes = [c for c in [">", "=", "<"] if c in set(df["reversibility"])]

    # disagreement: cascade hedged "=" but stats are confident either way
    hedged = df[(df["reversibility"] == "=") &
                ((df["p_forward"] >= 0.95) | (df["p_forward"] <= 0.05))]
    # cascade committed > / < but stats are not confident
    overcommit = df[(df["reversibility"].isin([">", "<"])) &
                    (df["p_forward"].between(0.05, 0.95))]
    note = (f"{len(df)} reactions have an analytic P(direction). "
            f"The cascade marks {(df['reversibility'] == '=').sum()} of them reversible (=), "
            f"yet {len(hedged)} of those are statistically confident "
            f"(P(forward) ≥ 0.95 or ≤ 0.05). "
            f"{len(overcommit)} reactions are called > or < despite 0.05 < P(forward) < 0.95.")

    fig = make_subplots(
        rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.12,
        row_heights=[0.45, 0.55],
        subplot_titles=("distribution of P(forward) by cascade call",
                        "per-reaction P(forward), grouped by cascade call"),
    )
    for c in classes:
        sub = df[df["reversibility"] == c]
        fig.add_trace(
            go.Histogram(x=sub["p_forward"], name=CLASS_LABELS[c],
                         marker_color=CLASS_COLORS[c], opacity=0.75,
                         xbins=dict(start=0, end=1, size=0.05), legendgroup=c),
            row=1, col=1,
        )
    rng = np.random.RandomState(0)
    for j, c in enumerate(classes):
        sub = df[df["reversibility"] == c]
        y = j + (rng.rand(len(sub)) - 0.5) * 0.6
        fig.add_trace(
            go.Scatter(
                x=sub["p_forward"], y=y, mode="markers", name=CLASS_LABELS[c],
                marker=dict(color=CLASS_COLORS[c], size=6, opacity=0.7,
                            line=dict(color="#333", width=0.3)),
                customdata=np.stack([sub["rxn_id"], sub["p_reverse"]], -1),
                hovertemplate=("%{customdata[0]}<br>P(forward) %{x:.3f} | "
                               "P(reverse) %{customdata[1]:.3f}<extra></extra>"),
                legendgroup=c, showlegend=False,
            ),
            row=2, col=1,
        )
    fig.update_layout(barmode="overlay")
    for xline in (0.05, 0.5, 0.95):
        fig.add_vline(x=xline, line_dash="dot", line_color="#999", line_width=1, row=2, col=1)
    fig.update_yaxes(tickvals=list(range(len(classes))),
                     ticktext=[CLASS_LABELS[c] for c in classes], row=2, col=1)
    fig.update_xaxes(title_text="P(forward)  (P(reverse) = 1 − P(forward); P(reversible) ≈ 0)",
                     range=[-0.02, 1.02], row=2, col=1)
    fig.update_yaxes(title_text="reactions", row=1, col=1)
    fig.update_layout(
        title="Fig 5 — Analytic P(direction) vs the cascade's deterministic call",
        template="plotly_white", height=620, margin=dict(l=80, r=40, t=80, b=70),
    )
    return fig, note


def make_fig6_cascade_sankey(casc: pd.DataFrame) -> go.Figure:
    total = len(casc)
    stages = [s for s in STAGE_ORDER if s in set(casc["stage"])]
    classes = [c for c in [">", "<", "=", "?"] if c in set(casc["reversibility"])]

    src_node = "all reactions"
    labels = [f"{src_node} ({total:,})"]
    node_color = ["#666"]
    stage_idx = {}
    for s in stages:
        stage_idx[s] = len(labels)
        labels.append(f"{STAGE_LABELS.get(s, s)} ({int((casc['stage'] == s).sum()):,})")
        node_color.append("#9ecae1")
    class_idx = {}
    for c in classes:
        class_idx[c] = len(labels)
        labels.append(CLASS_LABELS[c])
        node_color.append(CLASS_COLORS[c])

    src, tgt, val, lcol = [], [], [], []
    # all -> stage
    for s in stages:
        src.append(0); tgt.append(stage_idx[s])
        val.append(int((casc["stage"] == s).sum())); lcol.append("rgba(150,150,150,0.35)")
    # stage -> class
    grp = casc.groupby(["stage", "reversibility"]).size()
    for (s, c), n in grp.items():
        if s in stage_idx and c in class_idx:
            src.append(stage_idx[s]); tgt.append(class_idx[c])
            val.append(int(n))
            rc = CLASS_COLORS.get(c, "#ccc")
            lcol.append("rgba" + str(tuple(list(_hex_rgb(rc)) + [0.4])))

    fig = go.Figure(go.Sankey(
        arrangement="snap",
        node=dict(label=labels, color=node_color, pad=16, thickness=16,
                  line=dict(color="#333", width=0.5),
                  hovertemplate="%{label}<br>%{value:,} reactions<extra></extra>"),
        link=dict(source=src, target=tgt, value=val, color=lcol,
                  hovertemplate="%{source.label} → %{target.label}<br>"
                                "%{value:,} reactions<extra></extra>"),
    ))
    fig.update_layout(
        title=f"Fig 6 — How {total:,} reactions flow through the reversibility "
              "cascade to a direction",
        template="plotly_white", height=620, margin=dict(l=20, r=20, t=70, b=20),
        font=dict(size=12),
    )
    return fig


def _hex_rgb(h: str) -> tuple[int, int, int]:
    h = h.lstrip("#")
    return int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)


# --------------------------------------------------------------------------- #
# Output assembly
# --------------------------------------------------------------------------- #
CAPTIONS = {
    "fig1": "Net effect of each thermodynamic-heuristic variant on panel growth flux "
            "(a proxy for ATP-generating capacity). Permissive variants that relax "
            "directionality (3.3_wide, 3.6, 3.10_loose, 3.7, ai_opus48) inflate flux toward "
            "the solver bound (~1000) — biologically implausible — while the principled "
            "composite (H4) and tighter bands (3.5, 3.10_tight) leave it essentially "
            "unchanged or slightly reduced.",
    "fig2": "Each variant placed by how many reactions it changes (x) against how many "
            "panel models flip between growing and non-growing (y). Variants in the "
            "upper-left change few reactions yet move biology the most — the high-leverage "
            "reactions worth curating by hand.",
    "fig3": "Per-model change in growth flux for every variant. Colour is signed-log so "
            "both small shifts and order-1000 blow-ups are visible. Horizontal bands of "
            "colour are models that are thermodynamically sensitive across many variants; "
            "grey/white rows are inert.",
    "fig4": "Growth flux when ΔG′° is resampled from its uncertainty and pushed through the "
            "whole cascade + FBA (50 Monte-Carlo draws). Each point is a model's median "
            "flux with its 90% credible interval; colour is the probability the model grows "
            "at all. Composite/probabilistic variants (H4, pforward_*) show many models "
            "dropping to P(grows) ≈ 0 once uncertainty is honoured.",
    "fig5": "Where the cascade's hard direction call sits relative to the analytic "
            "P(forward). The marginal posterior almost never supports 'reversible', so the "
            "many reactions the cascade marks '=' that pile up near P(forward)=0 or 1 are "
            "directions the heuristic left on the table.",
    "fig6": "Every reaction flows from the full set, through the heuristic stage that "
            "resolved it, to its assigned direction. Two thirds fall to 'incomplete thermo' "
            "and end up unknown (?), showing how much directionality is still gated on "
            "missing thermodynamic data rather than on the heuristics themselves.",
}

FIG_TITLES = {
    "fig1": "Impact tornado — ATP/growth-flux change per variant",
    "fig2": "Leverage — reactions changed vs models flipped",
    "fig3": "Variant × model Δ-flux heatmap",
    "fig4": "Monte-Carlo growth-flux under ΔG′° uncertainty",
    "fig5": "P(direction) vs the cascade's deterministic call",
    "fig6": "Reversibility cascade Sankey",
}


def write_figure(fig: go.Figure, key: str, out_dir: Path, do_png: bool,
                 do_site: bool, png_ok: list[bool]) -> None:
    """Write a single figure as standalone HTML (+ PNG, + site copy)."""
    html_path = out_dir / f"{key}_{_slug(key)}.html"
    fig.write_html(str(html_path), include_plotlyjs="cdn", full_html=True)
    print(f"wrote {html_path}")

    if do_site:
        SITE_FIG_DIR.mkdir(parents=True, exist_ok=True)
        site_path = SITE_FIG_DIR / f"{key}.html"
        fig.write_html(str(site_path), include_plotlyjs="cdn", full_html=True)

    if do_png and png_ok[0]:
        png_dir = out_dir / "png"
        png_dir.mkdir(parents=True, exist_ok=True)
        png_path = png_dir / f"{key}.png"
        try:
            fig.write_image(str(png_path), width=1300, height=fig.layout.height or 700,
                            scale=2)
            print(f"wrote {png_path}")
        except Exception as exc:  # kaleido/chrome unavailable
            png_ok[0] = False
            print(f"WARNING: PNG export failed ({type(exc).__name__}: {exc}). "
                  "Skipping remaining PNGs; HTML still produced.", file=sys.stderr)


def _slug(key: str) -> str:
    return {
        "fig1": "impact_tornado", "fig2": "leverage_scatter",
        "fig3": "variant_model_heatmap", "fig4": "mc_uncertainty",
        "fig5": "pdirection", "fig6": "cascade_sankey",
    }[key]


def build_index_html(entries: list[dict], figs: dict[str, go.Figure], out_dir: Path) -> None:
    """One scrollable dashboard embedding every built figure."""
    blocks = []
    first = True
    for e in entries:
        key = e["key"]
        div = figs[key].to_html(
            full_html=False,
            include_plotlyjs="cdn" if first else False,
            default_width="100%",
        )
        first = False
        blocks.append(
            f'<section id="{key}">'
            f'<h2>{e["title"]}</h2>'
            f'<p class="caption">{e["caption"]}</p>'
            f'<div class="figwrap">{div}</div>'
            f'</section>'
        )
    nav = " · ".join(f'<a href="#{e["key"]}">{e["key"].upper()}</a>' for e in entries)
    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Thermodynamics → core-model ATP production</title>
<style>
 body {{ font-family: -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif;
        margin: 0; color: #1a1a1a; background: #fafafa; }}
 header {{ background: #1f3a5f; color: #fff; padding: 28px 40px; }}
 header h1 {{ margin: 0 0 6px; font-size: 23px; }}
 header p {{ margin: 0; opacity: 0.85; font-size: 14px; }}
 nav {{ position: sticky; top: 0; background: #fff; border-bottom: 1px solid #ddd;
        padding: 10px 40px; font-size: 13px; z-index: 10; }}
 nav a {{ color: #1f3a5f; text-decoration: none; margin-right: 4px; }}
 section {{ background: #fff; margin: 26px 40px; padding: 22px 26px;
           border: 1px solid #e6e6e6; border-radius: 8px; }}
 h2 {{ font-size: 18px; margin: 0 0 8px; color: #1f3a5f; }}
 .caption {{ color: #444; font-size: 13.5px; line-height: 1.5; max-width: 1000px;
            margin: 0 0 14px; }}
 footer {{ color: #777; font-size: 12px; padding: 24px 40px 48px; }}
 code {{ background: #f0f0f0; padding: 1px 4px; border-radius: 3px; }}
</style></head>
<body>
<header>
 <h1>Reaction thermodynamics → core-model ATP production</h1>
 <p>How reversibility heuristics reshape growth/ATP flux across a 100-model ModelSEED panel</p>
</header>
<nav>{nav}</nav>
{''.join(blocks)}
<footer>
 Generated by <code>scripts/build_presentation_figures.py</code> from existing analysis
 artifacts (variant diff JSONs, the N=50 statistical panel, and the live reversibility
 cascade). No FBA was recomputed. Growth flux is the biomass objective
 (mmol·gDW⁻¹·h⁻¹); ATP is the dominant biomass cofactor, so panel growth flux tracks
 ATP-generating capacity.
</footer>
</body></html>"""
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "index.html").write_text(html)
    print(f"wrote {out_dir / 'index.html'}")


def build_site_manifest(entries: list[dict]) -> None:
    SITE_FIG_DIR.mkdir(parents=True, exist_ok=True)
    manifest = [
        {"key": e["key"], "file": f'{e["key"]}.html', "title": e["title"],
         "caption": e["caption"]}
        for e in entries
    ]
    (SITE_FIG_DIR / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"wrote {SITE_FIG_DIR / 'manifest.json'}")

    page = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Presentation figures</title>
<style>
 body { font-family: -apple-system, Segoe UI, Roboto, Arial, sans-serif; margin: 0;
        background: #fafafa; color: #1a1a1a; }
 header { background: #1f3a5f; color: #fff; padding: 22px 32px; }
 header h1 { margin: 0; font-size: 20px; }
 section { background: #fff; margin: 22px 32px; padding: 18px 22px;
           border: 1px solid #e6e6e6; border-radius: 8px; }
 h2 { color: #1f3a5f; font-size: 17px; margin: 0 0 6px; }
 p.caption { color: #444; font-size: 13px; max-width: 1000px; }
 iframe { width: 100%; height: 760px; border: 1px solid #e0e0e0; border-radius: 6px; }
</style></head>
<body>
<header><h1>Presentation figures — thermodynamics → ATP production</h1></header>
<div id="figs"></div>
<script>
fetch('data/figures/manifest.json').then(r => r.json()).then(items => {
  const root = document.getElementById('figs');
  for (const it of items) {
    const s = document.createElement('section');
    s.innerHTML = `<h2>${it.title}</h2><p class="caption">${it.caption}</p>` +
                  `<iframe src="data/figures/${it.file}" loading="lazy"></iframe>`;
    root.appendChild(s);
  }
});
</script>
</body></html>"""
    SITE_PAGE.write_text(page)
    print(f"wrote {SITE_PAGE}")


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def parse_selection(sel: str) -> set[int]:
    if not sel or sel.lower() == "all":
        return set(range(1, 7))
    out = set()
    for part in sel.split(","):
        part = part.strip()
        if "-" in part:
            a, b = part.split("-")
            out.update(range(int(a), int(b) + 1))
        elif part:
            out.add(int(part))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(DEFAULT_OUT), help="output directory")
    ap.add_argument("--figures", default="all", help="e.g. 'all', '1,3,6', '1-4'")
    ap.add_argument("--no-png", action="store_true", help="skip static PNG export")
    ap.add_argument("--no-site", action="store_true", help="skip site/data/figures export")
    ap.add_argument("--open", action="store_true", help="print a file:// link at the end")
    args = ap.parse_args(argv)

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    want = parse_selection(args.figures)
    do_png = not args.no_png
    do_site = not args.no_site
    png_ok = [do_png]

    variants = load_variant_jsons()
    summary = variant_summary(variants)
    panel = variant_panel_long(variants)

    builders = {}
    if 1 in want:
        builders["fig1"] = lambda: make_fig1_impact_tornado(summary)
    if 2 in want:
        builders["fig2"] = lambda: make_fig2_leverage_scatter(summary)
    if 3 in want:
        builders["fig3"] = lambda: make_fig3_variant_model_heatmap(panel)
    if 4 in want:
        frames, stat_summary = load_stat_panel()
        if frames:
            builders["fig4"] = lambda: make_fig4_mc_uncertainty(frames, stat_summary)
        else:
            print("WARNING: no statistical-panel files; skipping Fig 4", file=sys.stderr)
    notes = {}
    if 5 in want:
        pdir = load_pdirection()
        if pdir is not None:
            def _f5():
                fig, note = make_fig5_pdirection(pdir)
                notes["fig5"] = note
                return fig
            builders["fig5"] = _f5
        else:
            print("WARNING: P(direction)/cascade inputs missing; skipping Fig 5", file=sys.stderr)
    if 6 in want:
        casc = load_cascade_flows()
        if casc is not None:
            builders["fig6"] = lambda: make_fig6_cascade_sankey(casc)
        else:
            print("WARNING: cascade CSV missing; skipping Fig 6", file=sys.stderr)

    figs: dict[str, go.Figure] = {}
    entries: list[dict] = []
    for key in ["fig1", "fig2", "fig3", "fig4", "fig5", "fig6"]:
        if key not in builders:
            continue
        fig = builders[key]()
        figs[key] = fig
        write_figure(fig, key, out_dir, do_png, do_site, png_ok)
        caption = CAPTIONS[key] + (f" {notes[key]}" if key in notes else "")
        entries.append({"key": key, "title": FIG_TITLES[key], "caption": caption})

    if not figs:
        print("No figures built.", file=sys.stderr)
        return 1

    build_index_html(entries, figs, out_dir)
    if do_site:
        build_site_manifest(entries)

    print(f"\nDone: {len(figs)} figure(s).")
    if args.open:
        print(f"open: file://{(out_dir / 'index.html').resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
