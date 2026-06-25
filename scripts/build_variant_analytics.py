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

    # --- load each variant's all-models FBA once: Δ vectors + non-trivial deltas ---
    delta_vec = {}     # tag -> {model_id: delta_flux}
    n_eval = {}        # tag -> models evaluated (rows in that variant's FBA file)
    for t in tags:
        if t == "baseline":
            continue
        f = SITE_DATA / f"all_models_variant_fba__{t}.json"
        if not f.exists():
            continue
        rows = json.loads(f.read_text())
        n_eval[t] = len(rows)
        delta_vec[t] = {r["model_id"]: float(r.get("delta_flux", 0.0)) for r in rows}

    # --- per-variant Δ growth histograms (own symmetric range; no silent clamping) ---
    delta_hist = {}
    for t, vec in delta_vec.items():
        ds = [d for d in vec.values() if abs(d) > 1e-6]
        if not ds:
            continue  # variant moves nothing (e.g. no panel/DB intersection) -> omit
        rng = max(abs(min(ds)), abs(max(ds)), 1e-6)
        delta_hist[t] = {**histogram(ds, -rng, rng, 41),
                         "n_moved": len(ds), "n_evaluated": n_eval.get(t, len(ds))}

    # --- effect-similarity: Pearson corr of per-model Δgrowth vectors (over the DB) ---
    def _pearson(a, b):
        n = len(a)
        ma = sum(a) / n
        mb = sum(b) / n
        cov = va = vb = 0.0
        for ai, bi in zip(a, b):
            da, db = ai - ma, bi - mb
            cov += da * db
            va += da * da
            vb += db * db
        return cov / ((va * vb) ** 0.5) if va > 0 and vb > 0 else 0.0

    corr_tags = [t for t in tags if t in delta_vec
                 and sum(1 for d in delta_vec[t].values() if abs(d) > 1e-6) >= 2]
    universe = sorted(set().union(*[set(delta_vec[t]) for t in corr_tags])) if corr_tags else []
    vecs = {t: [delta_vec[t].get(m, 0.0) for m in universe] for t in corr_tags}
    nc = len(corr_tags)
    effect_corr = [[0.0] * nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i, nc):
            c = 1.0 if i == j else round(_pearson(vecs[corr_tags[i]], vecs[corr_tags[j]]), 4)
            effect_corr[i][j] = effect_corr[j][i] = c
    print(f"[analytics] effect-correlation matrix {nc}x{nc} over {len(universe)} models")

    base_path = SITE_DATA / "all_models_baseline_fba.json"
    n_all = manifest.get("n_all_models") or (max(n_eval.values()) if n_eval else None)

    OUT.write_text(json.dumps({
        "tags": tags,
        "n_all_models": n_all,
        "agreement": agreement,
        "agreement_n": agreement_n,
        "direction_dist": direction_dist,
        "delta_hist": delta_hist,
        "corr_tags": corr_tags,
        "effect_corr": effect_corr,
    }, separators=(",", ":")))
    print(f"[analytics] wrote {OUT.name} ({OUT.stat().st_size/1024:.0f} KB)")


if __name__ == "__main__":
    main()
