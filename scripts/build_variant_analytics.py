#!/usr/bin/env python3
"""Precompute cross-variant analytics for the site's Analytics tab.

Reads the per-variant direction reports (thermo_variants/<tag>/...EQ.txt) and the
all-models FBA outputs and emits one compact JSON the browser renders as several
intense overview visualizations:

  - agreement[][]      pairwise directional agreement between every pair of
                       variants, over reactions BOTH call directional/reversible
                       (i.e. excluding reactions either leaves '?'), so the
                       heatmap reflects real disagreement rather than shared '?'.
  - agreement_n[][]    number of co-decided reactions behind each cell.
  - direction_dist     per-variant counts of '>','<','=','?'.
  - delta_hist         per-variant histogram of Δ growth flux over the models the
                       variant actually moves (|Δ|>1e-6), on a shared bin grid.
  - baseline_flux_hist all-models baseline growth-flux distribution (reference).

Output: site/data/variant_analytics.json. Run after the variant reports +
build_all_models_impact exist.
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
ANALYSIS_ROOT = SCRIPTS.parent
SITE_DATA = ANALYSIS_ROOT / "site" / "data"
TV = ANALYSIS_ROOT / "thermo_variants"
OUT = SITE_DATA / "variant_analytics.json"

DIRS = (">", "<", "=", "?")


def load_dir_map(tag: str) -> dict:
    out = {}
    rep = TV / tag / "Estimated_Reaction_Reversibility_Report_EQ.txt"
    with open(rep) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4:
                out[p[0]] = p[3]
    return out


def histogram(values, lo, hi, nbins):
    """Counts of values into nbins equal bins on [lo,hi]; outliers clamped to ends."""
    counts = [0] * nbins
    if hi <= lo:
        hi = lo + 1.0
    width = (hi - lo) / nbins
    for v in values:
        idx = int((v - lo) / width)
        if idx < 0:
            idx = 0
        elif idx >= nbins:
            idx = nbins - 1
        counts[idx] += 1
    edges = [round(lo + i * width, 3) for i in range(nbins + 1)]
    return {"bins": edges, "counts": counts}


def main() -> None:
    manifest = json.loads((SITE_DATA / "manifest.json").read_text())
    tags = [v["tag"] for v in manifest["variants"]]
    print(f"[analytics] {len(tags)} variants")

    maps = {t: load_dir_map(t) for t in tags}

    # --- direction distribution per variant ---
    direction_dist = {t: {d: 0 for d in DIRS} for t in tags}
    for t in tags:
        c = Counter(maps[t].values())
        for d in DIRS:
            direction_dist[t][d] = int(c.get(d, 0))

    # --- pairwise agreement over co-decided (both non-'?') reactions ---
    rxn_ids = set(maps[tags[0]])
    decided = {t: {r: v for r, v in maps[t].items() if v in (">", "<", "=")} for t in tags}
    n = len(tags)
    agreement = [[0.0] * n for _ in range(n)]
    agreement_n = [[0] * n for _ in range(n)]
    for i in range(n):
        di = decided[tags[i]]
        for j in range(i, n):
            dj = decided[tags[j]]
            common = di.keys() & dj.keys()
            if common:
                same = sum(1 for r in common if di[r] == dj[r])
                frac = same / len(common)
            else:
                frac = 1.0
            agreement[i][j] = agreement[j][i] = round(frac, 4)
            agreement_n[i][j] = agreement_n[j][i] = len(common)
    print(f"[analytics] agreement matrix {n}x{n} done")

    # --- per-variant Δ growth histograms (models the variant actually moves) ---
    # shared bin range from the global spread of non-trivial deltas
    all_deltas = {}
    gmax = 1.0
    for t in tags:
        if t == "baseline":
            continue
        f = SITE_DATA / f"all_models_variant_fba__{t}.json"
        if not f.exists():
            continue
        rows = json.loads(f.read_text())
        ds = [float(r.get("delta_flux", 0.0)) for r in rows
              if abs(float(r.get("delta_flux", 0.0))) > 1e-6]
        all_deltas[t] = ds
        for d in ds:
            gmax = max(gmax, abs(d))
    # cap the shared range so the bulk is visible (clamp extreme tails)
    rng = min(gmax, 80.0)
    delta_hist = {t: {**histogram(ds, -rng, rng, 40), "n_moved": len(ds)}
                  for t, ds in all_deltas.items()}

    # --- baseline overall growth-flux distribution ---
    base_path = SITE_DATA / "all_models_baseline_fba.json"
    baseline_flux_hist = None
    if base_path.exists():
        rows = json.loads(base_path.read_text())
        fluxes = [float(r.get("growth_flux", 0.0)) for r in rows]
        fmax = max(fluxes) if fluxes else 1.0
        baseline_flux_hist = histogram(fluxes, 0.0, fmax, 40)

    OUT.write_text(json.dumps({
        "tags": tags,
        "n_all_models": manifest.get("n_all_models") or len(json.loads(base_path.read_text())) if base_path.exists() else None,
        "agreement": agreement,
        "agreement_n": agreement_n,
        "direction_dist": direction_dist,
        "delta_hist": delta_hist,
        "baseline_flux_hist": baseline_flux_hist,
    }, separators=(",", ":")))
    print(f"[analytics] wrote {OUT.name} ({OUT.stat().st_size/1024:.0f} KB)")


if __name__ == "__main__":
    main()
