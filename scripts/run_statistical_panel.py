#!/usr/bin/env python3
"""Monte-Carlo + analytic-P propagation of CC posterior through the panel.

For each ``ReversibilityConfig`` variant in ``--variants`` (default:
``baseline``, ``3.5``, ``H4``, plus the new posterior-rule variants
``pforward_50`` / ``pforward_95``), we do two things:

  1. **Analytic per-reaction P(direction)** -- one row per panel reaction
     with ``P(forward)``, ``P(reverse)``, ``P(reversible)`` from the
     marginal CC normal on ΔG′° at the cascade's concentration term.

  2. **Monte-Carlo cascade replay** -- N samples; each sample resamples
     each reaction's ΔG′° from ``N(deltag, deltagerr)``, replays the full
     cascade, then runs the 100-model panel FBA.  We aggregate
     5/50/95% growth-flux quantiles + ``P(grows)`` per panel model.

Outputs (under ``results/statistical_panel/``):

  - ``panel_distribution__{variant}__N{n}.csv``  per-model flux quantiles
  - ``p_direction__{variant}.csv``               per-panel-rxn P(direction)
  - ``summary.csv``                              one row per variant
  - ``pivot.md``                                 human-readable interpretation

The MC sampler assumes per-reaction independence under the marginal CC
normal -- §3.4 of the review wants the full covariance matrix
(``standard_dg_prime_multi``), but that covariance is not on disk.
Documented in pivot.md.

All heavy state is cached in the kbcache shared with notebook 06.
Cold full run (3 configs × 50 samples) takes ~20-30 min on a workstation;
re-run is a no-op.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import os
import statistics
import sys
import time
from pathlib import Path
from typing import Optional

SCRIPTS = Path(__file__).resolve().parent
ANALYSIS_ROOT = SCRIPTS.parent
MSDB_ROOT = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase"))

sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(MSDB_ROOT / "Libs" / "Python"))

import reversibility_lib as lib       # noqa: E402
import reversibility_stats as rstats  # noqa: E402
import growth_heuristics as gh        # noqa: E402
import variant_catalog                # noqa: E402

PANEL_FILE = ANALYSIS_ROOT / "results" / "selected_ids.txt"
OUT_DIR = ANALYSIS_ROOT / "results" / "statistical_panel"


# Variants exercised by default.  Extend with the two new posterior-rule
# configs that exist only here (not in variant_catalog -- they are
# experiments of the statistical implementation, not heuristic-review items
# the catalog promises to ship).
def _cfg_pforward(threshold: float):
    def _build():
        return lib.ReversibilityConfig(p_forward_threshold=threshold)
    return _build


EXTRA_VARIANTS = [
    {"tag": "pforward_50", "title": "P(direction) >= 0.50 rule (§3.6)",
     "section": "§ 3.6 (new)", "cfg": _cfg_pforward(0.50)},
    {"tag": "pforward_95", "title": "P(direction) >= 0.95 rule (§3.6)",
     "section": "§ 3.6 (new)", "cfg": _cfg_pforward(0.95)},
]

DEFAULT_VARIANT_TAGS = ["baseline", "3.5", "H4", "pforward_50", "pforward_95"]


# ---------------------------------------------------------------------------
def load_session():
    from kbutillib.notebook import NotebookSession
    return NotebookSession.for_notebook(
        notebook_file=str(ANALYSIS_ROOT / "notebooks" /
                          "06_ReactionReversibilityHeuristics.ipynb"),
        project_name="core_models_analysis",
    )


def resolve_variants(tags):
    """Look up tags in variant_catalog + EXTRA_VARIANTS."""
    extras = {v["tag"]: v for v in EXTRA_VARIANTS}
    out = []
    for t in tags:
        if t in extras:
            out.append(extras[t])
        else:
            out.append(variant_catalog.variant_by_tag(t))
    return out


# ---------------------------------------------------------------------------
def mc_cascades_key(tag: str, n: int) -> str:
    return f"mc_cascades__{tag}__N{n}_v1"


def mc_fba_key(tag: str, n: int) -> str:
    return f"mc_fba__{tag}__N{n}_v1"


def ensure_mc(session, rxns, variant, n_samples: int,
              panel_ids, n_workers: int, seed: int = 0,
              verbose: bool = False):
    """Compute (and cache) the MC cascades + per-sample panel FBA.

    Returns ``(samples_revmaps, samples_fba)``.
    """
    tag = variant["tag"]
    cfg = variant["cfg"]()

    cascades = None
    fbas = None
    try:
        cascades = session.cache.load(mc_cascades_key(tag, n_samples))["samples"]
        print(f"  [mc] cascade cache hit ({len(cascades)} samples)")
    except KeyError:
        print(f"  [mc] resampling cascade, N={n_samples} ...")
        t = time.time()
        # sample_cascade mutates ``rxns`` in place and restores at exit,
        # so we pass the cached dict directly without an outer deepcopy --
        # otherwise we pay ~40 s per call for the deepcopy alone.
        cascades = rstats.sample_cascade(rxns, cfg=cfg, n_samples=n_samples,
                                         seed=seed, verbose=verbose)
        # kbcache only knows dict / json / dataframe -- wrap the list.
        session.cache.save(mc_cascades_key(tag, n_samples),
                           {"samples": cascades}, type_hint="dict",
                           metadata={"variant": tag, "n_samples": n_samples,
                                     "seed": seed})
        print(f"  [mc] cascade samples in {time.time()-t:.1f}s")

    try:
        fbas = session.cache.load(mc_fba_key(tag, n_samples))["samples"]
        print(f"  [mc] FBA cache hit ({len(fbas)} samples)")
    except KeyError:
        print(f"  [mc] panel FBA per sample (N={len(cascades)}, "
              f"{n_workers} workers per sample) ...")
        t = time.time()
        fbas = rstats.panel_fba_for_samples(panel_ids, cascades,
                                            n_workers=n_workers,
                                            verbose=verbose)
        session.cache.save(mc_fba_key(tag, n_samples),
                           {"samples": fbas}, type_hint="dict",
                           metadata={"variant": tag, "n_samples": n_samples,
                                     "seed": seed})
        print(f"  [mc] FBA samples in {time.time()-t:.1f}s")
    return cascades, fbas


# ---------------------------------------------------------------------------
def compute_point_estimate_fba(session, rxns, variant, panel_ids, n_workers):
    """Recompute (or load) the point-estimate panel FBA for a given variant."""
    tag = variant["tag"]
    cfg = variant["cfg"]()
    # Re-use notebook-06 cache where applicable.
    try:
        if tag == "baseline":
            return session.cache.load("fba_heuristic_baseline_v1")
        ts = tag.replace(".", "_")
        return session.cache.load(f"fba_variant_{ts}_v2")
    except KeyError:
        pass
    # Build the rev_map fresh.
    rev_map = lib.run_cascade(copy.deepcopy(rxns), db_level="EQ", cfg=cfg,
                              gc_first=True)
    rev_map = {r: rev for r, (_, rev) in rev_map.items()}
    return gh.run_panel(panel_ids, reversibility_map=rev_map,
                        baseline_map=None, n_workers=n_workers)


# ---------------------------------------------------------------------------
def main(argv: Optional[list] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--n-samples", type=int, default=50,
                    help="number of MC samples per variant (default 50)")
    ap.add_argument("--variants", action="append", default=None,
                    help="restrict to these variant tags (repeatable). "
                         "Default: " + ", ".join(DEFAULT_VARIANT_TAGS))
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args(argv)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    session = load_session()

    panel_ids = PANEL_FILE.read_text().split()
    print(f"panel: {len(panel_ids)} models, n_samples={args.n_samples}")

    msdb_rxns = session.cache.load("msdb_reactions_v1")
    # Load panel reaction union -- used to size the per-rxn P(direction) table.
    panel_rxnsets_path = ANALYSIS_ROOT / "site" / "data" / "panel_rxnsets.json"
    panel_rxn = set()
    if panel_rxnsets_path.exists():
        panel = json.loads(panel_rxnsets_path.read_text())
        for s in panel.values():
            panel_rxn.update(s)
    print(f"MSDB: {len(msdb_rxns)} reactions | panel union: {len(panel_rxn)}")

    tags = args.variants or DEFAULT_VARIANT_TAGS
    variants = resolve_variants(tags)
    print(f"variants: {[v['tag'] for v in variants]}")

    summary_rows = []
    for v in variants:
        tag = v["tag"]
        print(f"\n=== variant {tag!r}: {v['title']}")
        cfg = v["cfg"]()

        # --- 1. Analytic P(direction) on panel reactions --------------------
        agg = rstats.aggregate_p_direction(msdb_rxns, cfg=cfg,
                                            restrict_to=panel_rxn or None)
        c = agg["counts"]
        print(f"  P(direction) on {c['total_with_dg']} panel-rxn ΔG′° calls: "
              f"confident_fwd={c['confident_forward']}, "
              f"confident_rev={c['confident_reverse']}, "
              f"ambiguous={c['ambiguous']}")
        pdir_path = OUT_DIR / f"p_direction__{tag}.csv"
        with open(pdir_path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["rxn_id", "p_forward", "p_reverse", "p_reversible"])
            for rid, (pf, pr, prv) in sorted(agg["per_rxn"].items()):
                w.writerow([rid, f"{pf:.6f}", f"{pr:.6f}", f"{prv:.6f}"])

        # --- 2. MC cascade + panel FBA --------------------------------------
        cascades, fbas = ensure_mc(session, msdb_rxns, v, args.n_samples,
                                    panel_ids, args.workers, seed=args.seed,
                                    verbose=args.verbose)

        # --- 3. Per-model flux distribution ---------------------------------
        dist = rstats.summarize_panel_distribution(fbas,
                                                    quantiles=(0.05, 0.5, 0.95))
        # Also: point-estimate FBA for this variant (centered ΔG)
        pe = compute_point_estimate_fba(session, msdb_rxns, v, panel_ids,
                                         args.workers)
        pe_by_id = {r["model_id"]: r for r in pe}

        dist_path = OUT_DIR / f"panel_distribution__{tag}__N{args.n_samples}.csv"
        with open(dist_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=[
                "model_id", "n_samples", "mean_flux",
                "q05", "q50", "q95", "p_grows",
                "point_estimate_flux", "point_in_ci95",
            ])
            w.writeheader()
            for mid, st in sorted(dist.items()):
                pe_row = pe_by_id.get(mid, {})
                pe_flux = float(pe_row.get("growth_flux", 0.0))
                in_ci = st["q05"] <= pe_flux <= st["q95"]
                w.writerow({
                    "model_id":            mid,
                    "n_samples":           st["n_samples"],
                    "mean_flux":           f"{st['mean_flux']:.6f}",
                    "q05":                 f"{st['q05']:.6f}",
                    "q50":                 f"{st['q50']:.6f}",
                    "q95":                 f"{st['q95']:.6f}",
                    "p_grows":             f"{st['p_grows']:.4f}",
                    "point_estimate_flux": f"{pe_flux:.6f}",
                    "point_in_ci95":       int(in_ci),
                })

        # --- 4. Variant-level summary --------------------------------------
        flux_widths = [st["q95"] - st["q05"] for st in dist.values()]
        flux_widths_grow = [st["q95"] - st["q05"] for st in dist.values() if st["p_grows"] > 0]
        p_grows = [st["p_grows"] for st in dist.values()]
        n_models = len(dist)
        n_always_grow = sum(1 for p in p_grows if p == 1.0)
        n_never_grow = sum(1 for p in p_grows if p == 0.0)
        n_uncertain = n_models - n_always_grow - n_never_grow
        # point-estimate coverage
        n_in_ci = sum(1 for mid, st in dist.items()
                      if st["q05"] <= float(pe_by_id.get(mid, {}).get("growth_flux", 0.0)) <= st["q95"])
        # rev_map variance: how many panel reactions had > 1 unique rev across samples
        rxn_diversity = 0
        if panel_rxn and cascades:
            for rid in panel_rxn:
                vals = {s.get(rid) for s in cascades}
                vals.discard(None)
                if len(vals) > 1:
                    rxn_diversity += 1
        summary_rows.append({
            "variant":               tag,
            "title":                 v["title"],
            "n_samples":             args.n_samples,
            "n_panel_models":        n_models,
            "n_always_grow":         n_always_grow,
            "n_never_grow":          n_never_grow,
            "n_uncertain":           n_uncertain,
            "mean_p_grows":          (sum(p_grows) / n_models) if n_models else 0,
            "median_flux_width_ci90": statistics.median(flux_widths_grow) if flux_widths_grow else 0.0,
            "max_flux_width_ci90":   max(flux_widths_grow) if flux_widths_grow else 0.0,
            "n_models_point_in_ci95": n_in_ci,
            "n_panel_rxns_with_direction_variance": rxn_diversity,
        })
        print(f"  panel models: always-grow={n_always_grow}, never-grow={n_never_grow}, "
              f"uncertain (0<P(grows)<1)={n_uncertain}")
        print(f"  point estimate inside CI95 for {n_in_ci}/{n_models} models")
        print(f"  panel rxns with sample-direction variance: {rxn_diversity}")

    # ---- summary table + pivot.md ----
    sum_path = OUT_DIR / "summary.csv"
    with open(sum_path, "w", newline="") as fh:
        if summary_rows:
            w = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()))
            w.writeheader()
            w.writerows(summary_rows)

    md = []
    md.append("# Statistical propagation of CC ΔG′° uncertainty through the panel")
    md.append("")
    md.append(f"Panel size: **{len(panel_ids)} models** | "
              f"MC sample size: **N = {args.n_samples}** per variant | "
              f"Seed: {args.seed}")
    md.append("")
    md.append("Two methods from `Reaction_Reversibility_Heuristics_Review.md`:")
    md.append("")
    md.append("- **§ 3.6 — analytic P(direction).**  Per reaction, compute "
              "`P(mMdeltaG < 0)`, `P(mMdeltaG > 0)`, "
              "`P(|mMdeltaG| ≤ 0)` from the marginal CC normal "
              "`ΔG′° ~ N(deltag, deltagerr)` (treating the concentration term "
              "`RT·Σνᵢ ln cᵢ` as fixed).  Per-panel-rxn CSVs in "
              "`p_direction__{variant}.csv`.")
    md.append("- **§ 2.5 — Monte-Carlo cascade + FBA.**  Resample each "
              "reaction's ΔG′° from its marginal normal; replay the full "
              "cascade; run FBA on the panel; aggregate 5/50/95% growth-flux "
              "quantiles + `P(grows)` per model.  Per-model rows in "
              "`panel_distribution__{variant}__N{n}.csv`.")
    md.append("")
    md.append("**Independence caveat.** §3.4 of the review asks for the full "
              "CC covariance matrix (`standard_dg_prime_multi`).  That covariance "
              "is not on disk -- only the marginal σ per reaction is -- so the "
              "MC sampler resamples each ΔG′° independently.  Correlations "
              "between reactions that share component-contribution groups are "
              "therefore *not* propagated; the resulting flux distributions "
              "are wider than the true posterior would give.")
    md.append("")
    md.append("## Per-variant summary")
    md.append("")
    md.append("| variant | title | always-grow | never-grow | uncertain | mean P(grows) | median CI90 width | point ∈ CI95 | rxns w/ sample-direction variance |")
    md.append("|---------|-------|------------:|-----------:|----------:|--------------:|------------------:|-------------:|----------------------------------:|")
    for r in summary_rows:
        md.append(f"| `{r['variant']}` | {r['title']} | "
                  f"{r['n_always_grow']} | {r['n_never_grow']} | {r['n_uncertain']} | "
                  f"{r['mean_p_grows']:.3f} | {r['median_flux_width_ci90']:.3f} | "
                  f"{r['n_models_point_in_ci95']}/{r['n_panel_models']} | "
                  f"{r['n_panel_rxns_with_direction_variance']} |")
    md.append("")
    md.append("**Column glossary.**")
    md.append("- `always-grow` / `never-grow` / `uncertain` — partitions the "
              "100 panel models by `P(grows)` across the MC samples.")
    md.append("- `median CI90 width` — among grow-capable models, the median "
              "of `q95 − q05` of the growth-flux distribution. Bigger = more "
              "flux uncertainty propagated by the variant.")
    md.append("- `point ∈ CI95` — how often the point-estimate FBA (at "
              "centered ΔG′°) falls inside the MC's 5–95% interval. "
              "A high count means the point estimate is well-calibrated.")
    md.append("- `rxns w/ sample-direction variance` — panel reactions whose "
              "MC samples disagreed on direction.  These are the reactions "
              "actually driving flux variance.")
    md.append("")
    md.append("## How to interpret the new `pforward_*` variants")
    md.append("")
    md.append("The two `pforward_50` / `pforward_95` variants set the cascade's "
              "new `cfg.p_forward_threshold` knob (in `reversibility_lib`).  "
              "At `0.95`, the cascade drops both the ±2 kcal mMdeltaG band "
              "and the low-energy-points heuristic in favor of a single "
              "posterior-probability rule: a reaction is `>` if "
              "`P(mMdeltaG < 0) ≥ 0.95`, `<` if `P(mMdeltaG > 0) ≥ 0.95`, "
              "and reversible otherwise.  At `0.50` the rule reduces to "
              "the sign of the centered mMdeltaG.")

    pivot_path = OUT_DIR / "pivot.md"
    pivot_path.write_text("\n".join(md) + "\n")
    print(f"\nwrote {sum_path}")
    print(f"wrote {pivot_path}")
    print("\n" + "\n".join(md))


if __name__ == "__main__":
    main()
