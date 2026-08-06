#!/usr/bin/env python3
"""Why do Group Contribution and eQuilibrator track each other (r = 0.91) while
dGPredictor correlates with neither (r = 0.08 / 0.17)?

Runs the diagnostics behind the "Pairwise DeltaG agreement across sources" section
of reports/thermoComparison/THERMO_SOURCE_FBA_PIPELINE.md. Four tests:

  1. ADDITIVITY -- reconstruct each source's stored reaction DeltaG as
     sum(nu_i * dGf_i) from the *compound* table. GC and eQuilibrator both
     reconstruct exactly (r = 1.0000); dGPredictor has no compound-level
     formation energies in the database at all, so it structurally cannot.
     This is the root cause: GC and eQ are two linear maps of the SAME
     stoichiometry vector, so they are correlated by construction.

  2. SIZE SCALING -- because GC/eQ sum over stoichiometry, |DeltaG| must grow
     with reaction size. A direct reaction-level regressor has no such
     mechanism. Reported as r(|dG|, sum|coeff|) and the median-|dG| ratio
     between large- and small-stoichiometry reactions.

  3. OUTPUT RANGE -- dGPredictor's outputs are effectively bounded (~+/-400
     kcal/mol) while GC/eQ run to +/-16,000 on aggregate/polymer reactions.

  4. RANGE RESTRICTION -- Pearson r for each pair as a function of the |dG|
     window, plus the same statistics restricted to the 239 reactions that
     actually occur across the combined core models. Shows the headline r
     values are dominated by a handful of extreme-leverage reactions.

Usage:
    python3 scripts/analyze_thermo_source_dg_agreement.py [--subset PATH]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
MSDB_ROOT = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase"))
sys.path.insert(0, str(MSDB_ROOT / "Libs" / "Python"))
sys.path.insert(0, str(MSDB_ROOT / "Scripts" / "Thermodynamics"))

from BiochemPy import Compounds, Reactions  # noqa: E402
from reversibility_heuristics import (  # noqa: E402
    DEFAULT_HEURISTICS, run_reversibility, per_source_energy,
)

SOURCE_LABEL = {
    "group_contribution": "Group contribution",
    "equilibrator": "eQuilibrator",
    "dgpredictor": "dGPredictor",
}
SHORT = {"group_contribution": "GC", "equilibrator": "eQ", "dgpredictor": "dGP"}
PAIRS = [("group_contribution", "equilibrator"),
         ("group_contribution", "dgpredictor"),
         ("equilibrator", "dgpredictor")]
SENTINEL = 1e6


def load_source_data(reactions: dict) -> dict:
    """{source: {rxn_id: (dg, operator)}} -- same selection rule as the scatter script."""
    out = {src: {} for src in SOURCE_LABEL}
    for rxn_id, entry in reactions.items():
        if entry.get("status") == "EMPTY":
            continue
        for src, label in SOURCE_LABEL.items():
            _, operator, _ = run_reversibility(
                entry, per_source_energy(label), DEFAULT_HEURISTICS)
            if operator is None:
                continue
            pair = (entry.get("thermodynamics") or {}).get(label)
            stored_op = pair[2] if pair and len(pair) > 2 else None
            if stored_op in (None, "?"):
                continue
            out[src][rxn_id] = (float(pair[0]), operator)
    return out


def compound_dgf(cpds: dict, cid: str, label: str):
    entry = cpds.get(cid)
    if not entry:
        return None
    trip = (entry.get("thermodynamics") or {}).get(label)
    if not trip:
        return None
    val = float(trip[0])
    return None if abs(val) >= SENTINEL else val


def test_additivity(data, reactions, cpds) -> None:
    print("=" * 78)
    print("1. ADDITIVITY -- can the stored reaction dG be rebuilt as sum(nu_i * dGf_i)?")
    print("=" * 78)
    for src in ("group_contribution", "equilibrator"):
        label = SOURCE_LABEL[src]
        recon, stored = [], []
        for rid, (dg, _) in data[src].items():
            stoich = reactions[rid].get("stoichiometry")
            if not isinstance(stoich, list):
                continue
            total, ok = 0.0, True
            for part in stoich:
                val = compound_dgf(cpds, part.get("compound"), label)
                if val is None:
                    ok = False
                    break
                total += float(part.get("coefficient", 0) or 0) * val
            if ok:
                recon.append(total)
                stored.append(dg)
        recon, stored = np.array(recon), np.array(stored)
        resid = np.abs(recon - stored)
        print(f"  {label:<20} n={len(recon):>6}  r(recon, stored)={np.corrcoef(recon, stored)[0, 1]:.4f}"
              f"  median|resid|={np.median(resid):.3f}"
              f"  within 1 kcal={100 * (resid < 1).mean():.1f}%")
    print(f"  {'dGPredictor':<20} NO compound-level formation energies exist in the DB")
    print("  -> GC and eQ are two linear maps of the SAME stoichiometry vector;")
    print("     dGPredictor is a direct reaction-level regressor with no additive structure.\n")


def test_size_scaling(data, reactions) -> None:
    print("=" * 78)
    print("2. SIZE SCALING -- does |dG| grow with reaction size?")
    print("=" * 78)
    size = {}
    for rid, entry in reactions.items():
        stoich = entry.get("stoichiometry")
        if isinstance(stoich, list):
            size[rid] = sum(abs(float(p.get("coefficient", 0) or 0)) for p in stoich)
    print(f"  {'source':<20}{'r(|dG|, sum|coeff|)':>22}{'median|dG| big/small':>24}")
    for src, label in SOURCE_LABEL.items():
        ids = [r for r in data[src] if r in size]
        sz = np.array([size[r] for r in ids])
        val = np.abs(np.array([data[src][r][0] for r in ids]))
        big = val[sz >= 20]
        small = val[sz < 20]
        ratio = np.median(big) / max(np.median(small), 1e-9)
        print(f"  {label:<20}{np.corrcoef(sz, val)[0, 1]:>22.3f}"
              f"{f'{np.median(big):.1f} / {np.median(small):.1f}  = {ratio:.1f}x':>24}")
    print("  (big = sum|coeff| >= 20)\n")


def test_output_range(data) -> None:
    print("=" * 78)
    print("3. OUTPUT RANGE -- dG distribution per source (kcal/mol)")
    print("=" * 78)
    print(f"  {'source':<20}{'n':>8}{'std':>9}{'IQR':>8}{'min':>11}{'max':>11}")
    for src, label in SOURCE_LABEL.items():
        v = np.array([x[0] for x in data[src].values()])
        q1, q3 = np.percentile(v, [25, 75])
        print(f"  {label:<20}{len(v):>8}{v.std():>9.1f}{q3 - q1:>8.1f}{v.min():>11.1f}{v.max():>11.1f}")
    print("  Note the IQRs are comparable -- the compression is entirely in the TAILS.\n")


def _pair_stats(data, a, b, subset=None):
    common = sorted(set(data[a]) & set(data[b]))
    if subset is not None:
        common = [r for r in common if r in subset]
    x = np.array([data[a][r][0] for r in common])
    y = np.array([data[b][r][0] for r in common])
    rx, ry = np.argsort(np.argsort(x)), np.argsort(np.argsort(y))
    agree = np.mean([(data[a][r][1] == "=") == (data[b][r][1] == "=") for r in common])
    return dict(n=len(common), pearson=np.corrcoef(x, y)[0, 1],
                spearman=np.corrcoef(rx, ry)[0, 1], slope=np.polyfit(x, y, 1)[0],
                mad=np.abs(x - y).mean(), rev_agree=100 * agree, x=x, y=y)


def test_range_restriction(data, subset) -> None:
    print("=" * 78)
    print("4. RANGE RESTRICTION -- r as a function of the |dG| window")
    print("=" * 78)
    for a, b in PAIRS:
        common = sorted(set(data[a]) & set(data[b]))
        x = np.array([data[a][r][0] for r in common])
        y = np.array([data[b][r][0] for r in common])
        cells = []
        for cut in (50, 100, 200, 500, np.inf):
            m = (np.abs(x) <= cut) & (np.abs(y) <= cut)
            tag = "all" if np.isinf(cut) else f"<={cut:g}"
            cells.append(f"{tag}: {np.corrcoef(x[m], y[m])[0, 1]:.3f}")
        print(f"  {SHORT[a]}-{SHORT[b]:<4} " + "   ".join(cells))
    print()

    if subset is None:
        return
    print("=" * 78)
    print(f"   ALL ModelSEED reactions  vs  the {len(subset)} core-model reactions")
    print("=" * 78)
    print(f"  {'pair':<10}{'scope':<11}{'n':>6}{'Pearson':>9}{'Spearman':>10}"
          f"{'slope':>8}{'MAD':>8}{'rev-agree':>11}")
    for a, b in PAIRS:
        for scope, sub in (("all", None), ("core-239", subset)):
            s = _pair_stats(data, a, b, sub)
            print(f"  {SHORT[a]}-{SHORT[b]:<7}{scope:<11}{s['n']:>6}{s['pearson']:>9.3f}"
                  f"{s['spearman']:>10.3f}{s['slope']:>8.3f}{s['mad']:>8.1f}{s['rev_agree']:>10.1f}%")
    print("\n  Control -- is the subset result just range restriction?")
    print("  (all-reaction r recomputed inside the same |dG| window the subset occupies)")
    for a, b in PAIRS:
        s = _pair_stats(data, a, b, subset)
        lim = max(np.abs(s["x"]).max(), np.abs(s["y"]).max())
        f = _pair_stats(data, a, b, None)
        m = (np.abs(f["x"]) <= lim) & (np.abs(f["y"]) <= lim)
        print(f"    {SHORT[a]}-{SHORT[b]:<4} subset r={s['pearson']:.3f} (|dG|<={lim:.0f})"
              f"   | ALL rxns, same window r={np.corrcoef(f['x'][m], f['y'][m])[0, 1]:.3f}"
              f" (n={m.sum()})")
    print()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--subset", type=Path,
                    default=ANALYSIS_DIR / "results" / "core_models_unique_reactions.json",
                    help="JSON list of reaction ids for the subset comparison")
    args = ap.parse_args()

    print("loading live ModelSEEDDatabase checkout...")
    reactions = Reactions().loadReactions()
    cpds = Compounds().loadCompounds()
    print(f"  {len(reactions)} reactions, {len(cpds)} compounds\n")

    data = load_source_data(reactions)
    for src, label in SOURCE_LABEL.items():
        print(f"  {label}: {len(data[src])} reactions with a usable DeltaG")
    print()

    subset = None
    if args.subset and args.subset.exists():
        subset = set(json.loads(args.subset.read_text()))

    test_additivity(data, reactions, cpds)
    test_size_scaling(data, reactions)
    test_output_range(data)
    test_range_restriction(data, subset)


if __name__ == "__main__":
    main()
