#!/usr/bin/env python3
"""Pairwise DeltaG scatter plots across all ModelSEED reactions, one per pair of
thermodynamic sources: Group Contribution vs eQuilibrator, Group Contribution vs
dGPredictor, eQuilibrator vs dGPredictor.

For each pair, only reactions where BOTH sources have a non-sentinel DeltaG are
plotted (the intersection of coverage). Points are colored by *reversibility
transition* -- how the two sources' own reversibility calls (each run through
the unmodified ModelSEED heuristic cascade, reversibility_heuristics.DEFAULT_HEURISTICS,
fed that source's own DeltaG via per_source_energy) compare, collapsing each
call are compared:

  * No change                    -- the two sources make the IDENTICAL call:
    both reversible ("="), or both irreversible in the SAME direction (both ">"
    or both "<"). Nothing is in dispute.
  * Reversible -> Irreversible   -- source A reversible, source B irreversible
  * Irreversible -> Reversible   -- source A irreversible, source B reversible
  * Irreversible -> Irreversible -- both irreversible, in OPPOSITE directions
    (">" vs "<"). RESERVED for this case: it is the only genuine direction
    conflict, and it is what actually changes a flux model.

That last category is stricter than it was before 2026-08, when it held every
both-irreversible pair regardless of direction -- which coloured perfect
agreement identically to a reversal and made it the largest category in every
panel. "No change" is correspondingly larger now, so it is drawn smaller and
more transparent with the three real transitions on top.

The legend shows the reaction count for each category in parentheses.
"No change" is rendered in neutral gray (no hue) rather than a 4th categorical
hue: per the dataviz skill, all-pairs comparisons such as scatter plots only
keep 3 categorical hues CVD-safe at once (validate_palette.js confirms PASS on
the 3-hue triple used here); a neutral, chroma-free gray for one category adds
a 4th legend entry without reintroducing that risk.

Outputs three PNGs under reports/thermoComparison/figures/thermo_source_dg_scatter/.

With ``--subset PATH`` (a JSON list of reaction ids) the same three plots are
restricted to that set of reactions -- used to redraw the comparison over only
the 239 reactions that actually appear across the combined core models, rather
than all ~19k ModelSEED reactions. ``--out-subdir NAME`` redirects the output
directory so a subset run does not overwrite the all-reactions figures.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR", "/scratch/ctaylor/core_models_analysis"))
MSDB_ROOT = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase"))
OUT_DIR = ANALYSIS_DIR / "reports" / "thermoComparison" / "figures" / "thermo_source_dg_scatter"

sys.path.insert(0, str(MSDB_ROOT / "Libs" / "Python"))
sys.path.insert(0, str(MSDB_ROOT / "Scripts" / "Thermodynamics"))

from BiochemPy import Reactions  # noqa: E402
from reversibility_heuristics import (  # noqa: E402
    DEFAULT_HEURISTICS, run_reversibility, per_source_energy,
)

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_dgpredictor_kegg_mask import load_mask  # noqa: E402

SOURCE_LABEL = {
    "group_contribution": "Group contribution",
    "equilibrator": "eQuilibrator",
    "dgpredictor": "dGPredictor",
}
AXIS_TITLE = {
    "group_contribution": "Group Contribution ΔG′° (kcal/mol)",
    "equilibrator": "eQuilibrator (2.0) ΔG′° (kcal/mol)",
    "dgpredictor": "dGPredictor ΔG′° (kcal/mol)",
}
PAIRS = [
    ("group_contribution", "equilibrator"),
    ("group_contribution", "dgpredictor"),
    ("equilibrator", "dgpredictor"),
]

# dataviz skill: categorical hues assigned in fixed order, validated all-pairs
# (scatter) safe for the first 3 slots only -- validate_palette.js confirms
# PASS on this exact triple, light mode, --pairs all. "No change" uses neutral
# gray (no hue) instead of a 4th categorical slot -- see module docstring.
CATEGORY_ORDER = [
    "No change",
    "Reversible → Irreversible",
    "Irreversible → Reversible",
    "Irreversible → Irreversible",
]
CATEGORY_COLOR = {
    "No change": "#9c9a94",                    # neutral gray, not a categorical slot
    "Reversible → Irreversible": "#2a78d6",     # slot 1, blue
    "Irreversible → Reversible": "#eb6834",     # slot 2, orange
    "Irreversible → Irreversible": "#1baf7a",   # slot 3, aqua
}
# Displayed legend text. The category NAME is unchanged; the parenthetical
# states the definition, which is not guessable from the name alone.
CATEGORY_LEGEND = {
    "No change": "No change (same call)",
    "Reversible → Irreversible": "Reversible → Irreversible",
    "Irreversible → Reversible": "Irreversible → Reversible",
    "Irreversible → Irreversible": "Irreversible → Irreversible (opposite direction)",
}
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
BASELINE = "#c3c2b7"
SURFACE = "#fcfcfb"


def classify(op_a: str, op_b: str) -> str:
    """Reversibility transition from source A to source B.

    "Irreversible -> Irreversible" is RESERVED for the case that matters: both
    sources call the reaction irreversible but in OPPOSITE directions ('>' vs
    '<'). Two sources that agree on the same irreversible direction have not
    disagreed about anything, so they land in "No change" alongside the pairs
    that both call it reversible.
    """
    if op_a == op_b:
        return "No change"                    # both '=', or the same '>' / '<'
    if op_a == "=":
        return "Reversible → Irreversible"
    if op_b == "=":
        return "Irreversible → Reversible"
    return "Irreversible → Irreversible"      # '>' vs '<' — direction reversed


def load_source_data(reactions: dict, dgp_mask: set[str] | None = None) -> dict:
    """Returns {source: {rxn_id: (dg, operator)}} for GC/EQ/DGP.

    ``dgp_mask`` is the set of reactions whose stored dGPredictor value is
    attributable to a KEGG reaction that is not theirs (see
    build_dgpredictor_kegg_mask.py). Those reactions are dropped from the
    dGPredictor series only -- their Group Contribution and eQuilibrator values
    are unaffected by the KEGG mapping and are kept, so the GC-vs-eQuilibrator
    panel is identical with and without the mask.

    A reaction is kept for a source only when that source's OWN stored
    thermodynamics triple (``thermodynamics[label] = [dg, dge, operator]``, as
    written directly into the ModelSEEDDatabase reaction record) has BOTH a
    non-sentinel dG AND a defined (non-"?") stored operator -- i.e. the source
    itself considers the reaction's direction defined. In this checkout the two
    conditions are redundant (sentinel dG <=> stored operator "?"), but the
    "?" check is kept explicit so a reaction is never plotted as though a
    source had an opinion on it when that source's own record says otherwise.
    """
    dgp_mask = dgp_mask or set()
    out = {src: {} for src in SOURCE_LABEL}
    n_excluded_q = {src: 0 for src in SOURCE_LABEL}
    n_masked = 0
    for rxn_id, rxn_entry in reactions.items():
        if rxn_entry.get("status") == "EMPTY":
            continue
        for src, label in SOURCE_LABEL.items():
            if src == "dgpredictor" and rxn_id in dgp_mask:
                n_masked += 1
                continue
            status, operator, _ = run_reversibility(
                rxn_entry, per_source_energy(label), DEFAULT_HEURISTICS)
            if operator is None:
                continue
            thermo = rxn_entry.get("thermodynamics") or {}
            pair = thermo.get(label)
            stored_op = pair[2] if pair and len(pair) > 2 else None
            if stored_op == "?" or stored_op is None:
                n_excluded_q[src] += 1
                continue
            out[src][rxn_id] = (float(pair[0]), operator)
    for src, n in n_excluded_q.items():
        print(f"  {src}: excluded {n} reaction(s) with an undefined ('?') stored operator")
    if n_masked:
        print(f"  dgpredictor: excluded {n_masked} reaction(s) whose KEGG mapping is "
              f"not vouched by a ModelSEED alias (dgpredictor_kegg_mask.json)")
    return out


def plot_pair(src_a: str, src_b: str, data: dict, out_dir: Path | None = None,
              subset_note: str = "") -> Path:
    out_dir = OUT_DIR if out_dir is None else out_dir
    common = sorted(set(data[src_a]) & set(data[src_b]))
    xs = np.array([data[src_a][r][0] for r in common])
    ys = np.array([data[src_b][r][0] for r in common])
    cats = [classify(data[src_a][r][1], data[src_b][r][1]) for r in common]
    from collections import Counter
    cat_counts = Counter(cats)  # over ALL n common reactions, not just the in-range subset plotted below

    r_all = np.corrcoef(xs, ys)[0, 1] if len(xs) > 1 else float("nan")

    # A handful of reactions carry genuine (non-sentinel) but chemically implausible
    # DeltaG magnitudes -- e.g. rxn05017 (Group Contribution ~15,900 kcal/mol) -- that
    # would otherwise crush the entire rest of the distribution onto a few pixels on a
    # linear axis. The bulk distribution is heavy-tailed (lots of near-zero reactions
    # plus a long tail of biologically ordinary large reactions out to a few hundred
    # kcal/mol), so a statistical rule (percentile trim, IQR fence) flags hundreds of
    # perfectly ordinary points as "outliers." Instead use a fixed, chemically-motivated
    # cutoff: single/few-step reaction DeltaG'deg rarely exceeds ~1,000 kcal/mol in
    # magnitude in this database; a handful of aggregate/polymer reactions with large
    # stoichiometric coefficients do. Flag only those (documented, not silently dropped).
    EXTREME_CUTOFF = 1500.0
    in_range = (np.abs(xs) <= EXTREME_CUTOFF) & (np.abs(ys) <= EXTREME_CUTOFF)
    off_ids = [rid for rid, keep in zip(common, in_range) if not keep]

    xs_r, ys_r = xs[in_range], ys[in_range]
    r_robust = np.corrcoef(xs_r, ys_r)[0, 1] if in_range.sum() > 1 else float("nan")
    lo = min(xs_r.min(), ys_r.min())
    hi = max(xs_r.max(), ys_r.max())
    pad = 0.06 * (hi - lo)
    lo, hi = lo - pad, hi + pad

    fig, ax = plt.subplots(figsize=(7, 6.2), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    ax.plot([lo, hi], [lo, hi], linestyle="--",
             linewidth=1.2, color=BASELINE, zorder=1, label="_nolegend_")

    for cat in CATEGORY_ORDER:
        idx = [i for i, (c, keep) in enumerate(zip(cats, in_range)) if c == cat and keep]
        label = f"{CATEGORY_LEGEND[cat]} ({cat_counts.get(cat, 0):,})"
        bulk = cat == "No change"   # now the majority; keep it recessive
        ax.scatter(xs[idx], ys[idx], s=9 if bulk else 16, linewidths=0,
                    color=CATEGORY_COLOR[cat], alpha=0.35 if bulk else 0.7,
                    label=label, zorder=2 if bulk else 3)

    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel(AXIS_TITLE[src_a], color=INK_SECONDARY, fontsize=11)
    ax.set_ylabel(AXIS_TITLE[src_b], color=INK_SECONDARY, fontsize=11)
    ax.set_title(
        f"{SOURCE_LABEL[src_a]} vs {SOURCE_LABEL[src_b]}{subset_note}\n"
        f"n = {len(common):,} reactions with both sources' data · Pearson r = {r_all:.2f}",
        color=INK_PRIMARY, fontsize=12, pad=12,
    )
    if off_ids:
        import textwrap
        id_lines = textwrap.wrap(", ".join(off_ids), width=70)
        note = (f"axis zoomed to |ΔG'°| ≤ {EXTREME_CUTOFF:.0f} kcal/mol; "
                 f"{len(off_ids)} reaction(s) off-scale (shown-range r = {r_robust:.2f}):\n"
                 + "\n".join(id_lines))
        ax.text(0.5, -0.16, note, transform=ax.transAxes, ha="center", va="top",
                fontsize=7.5, color=INK_MUTED, linespacing=1.5)
    else:
        ax.text(0.5, -0.16, f"shown-range r = {r_robust:.2f} (no off-scale points)",
                transform=ax.transAxes, ha="center", va="top",
                fontsize=7.5, color=INK_MUTED, linespacing=1.5)
    ax.grid(True, color=GRIDLINE, linewidth=0.8, zorder=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(BASELINE)
    ax.tick_params(colors=INK_MUTED, labelsize=9)

    legend = ax.legend(loc="upper left", frameon=False, fontsize=9,
                        labelcolor=INK_PRIMARY, markerscale=1.8,
                        title="Reversibility transition", title_fontsize=9)
    legend.get_title().set_color(INK_SECONDARY)

    fig.tight_layout()
    if off_ids:
        n_lines = note.count("\n") + 1
        fig.subplots_adjust(bottom=0.14 + 0.028 * n_lines)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"dg_scatter_{src_a}_vs_{src_b}.png"
    fig.savefig(out_path, facecolor=SURFACE)
    plt.close(fig)

    # Machine-readable counts alongside the PNG: the legend is not a data source,
    # and prose elsewhere quotes these numbers.
    stats_path = out_dir / "category_counts.tsv"
    header = ("pair\tn\tpearson_r\tshown_range_r\t"
              + "\t".join(CATEGORY_ORDER) + "\n")
    line = (f"{src_a}_vs_{src_b}\t{len(common)}\t{r_all:.4f}\t{r_robust:.4f}\t"
            + "\t".join(str(cat_counts.get(c, 0)) for c in CATEGORY_ORDER) + "\n")
    if stats_path.exists() and stats_path.read_text().startswith(header):
        with open(stats_path, "a") as fh:
            fh.write(line)
    else:
        stats_path.write_text(header + line)

    print(f"wrote {out_path} ({len(common)} points)")
    return out_path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--subset", type=Path, default=None,
                    help="JSON list of reaction ids to restrict the plots to")
    ap.add_argument("--out-subdir", default=None,
                    help="write PNGs to reports/figures/<NAME>/ instead of the default")
    ap.add_argument("--no-dgp-mask", action="store_true",
                    help="do NOT drop dGPredictor values whose staged KEGG reaction id "
                         "is unbacked by a ModelSEED alias; reproduces the pre-mask "
                         "figures")
    ap.add_argument("--subset-label", default="",
                    help="short label appended to each plot title, e.g. ' — core-model reactions'")
    args = ap.parse_args()

    out_dir = OUT_DIR if args.out_subdir is None else (
        ANALYSIS_DIR / "reports" / "thermoComparison" / "figures" / args.out_subdir)

    print("loading reactions from live ModelSEEDDatabase checkout...")
    reactions = Reactions().loadReactions()
    print(f"  {len(reactions)} reactions loaded")

    dgp_mask = set() if args.no_dgp_mask else load_mask()
    data = load_source_data(reactions, dgp_mask=dgp_mask)
    for src, label in SOURCE_LABEL.items():
        print(f"  {label}: {len(data[src])} reactions with a usable DeltaG")

    if args.subset is not None:
        import json
        keep = set(json.loads(args.subset.read_text()))
        print(f"\nrestricting to subset of {len(keep)} reaction ids from {args.subset}")
        missing = keep - set(reactions)
        if missing:
            print(f"  warning: {len(missing)} subset id(s) absent from MSDB: "
                  f"{sorted(missing)[:8]}")
        data = {src: {r: v for r, v in d.items() if r in keep} for src, d in data.items()}
        for src, label in SOURCE_LABEL.items():
            print(f"  {label}: {len(data[src])} of {len(keep)} subset reactions "
                  f"have a usable DeltaG")

    (out_dir / "category_counts.tsv").unlink(missing_ok=True)
    for src_a, src_b in PAIRS:
        plot_pair(src_a, src_b, data, out_dir=out_dir, subset_note=args.subset_label)
    print(f"wrote {out_dir / 'category_counts.tsv'}")


if __name__ == "__main__":
    main()
