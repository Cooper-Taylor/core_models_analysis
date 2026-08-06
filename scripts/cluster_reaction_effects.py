#!/usr/bin/env python3
"""Starter clustering / analysis over the all-models reaction-effects dataset.

Reads results/reaction_effects_all/ (built by build_reaction_effects_all.py),
builds a per-(model, reaction) feature vector from the 4-direction sweep, and
KMeans-clusters reactions into influence archetypes (essential / reversibility-
sensitive / EGC-inducing / inert / ...).

Features per (model, rxn), for growing models (base_flux > tol):
  off_frac   growth kept when knocked out (0=lethal, ~1=dispensable)
  fwd_frac   growth forced forward     / base
  rev_frac   growth forced reverse     / base
  all_frac   growth forced reversible  / base
  dir_spread (max-min growth over the 4 options) / base   (direction sensitivity)
  gain       best-direction growth - default-direction growth, / base (upside)
  egc_atp    max ATP-drain rate across the 4 directions (>0 => can seed an EGC)
  egc_flag   1 if any direction induces an ATP EGC

Writes results/reaction_effects_all/reaction_features.parquet (feature matrix +
cluster label) for downstream work. Run: python3 scripts/cluster_reaction_effects.py
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

# Cap BLAS/OpenMP threads BEFORE importing numpy/sklearn — on high-core machines
# OpenBLAS otherwise over-subscribes threads and segfaults during KMeans.
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "8")

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parent.parent / "results" / "reaction_effects_all"
TOL = 1e-6


def build_features() -> pd.DataFrame:
    eff = ds.dataset(ROOT / "effects").to_table(
        columns=["model_id", "rxn", "base", "default_dir", "dir",
                 "growth", "n_active", "egc_atp"]).to_pandas()
    base_flux = {json.loads(l)["model_id"]: json.loads(l)["base_flux"]
                 for l in open(ROOT / "model_flux_loops.jsonl")}
    eff["base_flux"] = eff["model_id"].map(base_flux)
    eff = eff[eff["base_flux"] > TOL].copy()          # growing models only

    g = eff.pivot_table(index=["model_id", "rxn", "base", "default_dir", "base_flux"],
                        columns="dir", values="growth").reset_index()
    egc = (eff.groupby(["model_id", "rxn"])["egc_atp"].max()
           .rename("egc_atp").reset_index())
    g = g.merge(egc, on=["model_id", "rxn"], how="left")
    g = g.rename(columns={"<": "g_lt", ">": "g_gt", "=": "g_eq", "?": "g_off"})
    for c in ("g_lt", "g_gt", "g_eq", "g_off"):
        g[c] = g[c].fillna(0.0)
    bf = g["base_flux"]
    opts = g[["g_lt", "g_gt", "g_eq", "g_off"]].values
    # default-direction growth for the "gain" feature (vectorized)
    dir_col = {"<": "g_lt", ">": "g_gt", "=": "g_eq", "?": "g_off"}
    g_default = g["g_gt"].values.copy()  # fallback for default_dir not in the map (e.g. "0")
    dd = g["default_dir"].values
    for d, col in dir_col.items():
        mask = dd == d
        g_default[mask] = g[col].values[mask]

    feat = pd.DataFrame({
        "model_id": g["model_id"], "rxn": g["rxn"], "base": g["base"],
        "off_frac": (g["g_off"] / bf).clip(0, 5),
        "fwd_frac": (g["g_gt"] / bf).clip(0, 5),
        "rev_frac": (g["g_lt"] / bf).clip(0, 5),
        "all_frac": (g["g_eq"] / bf).clip(0, 5),
        "dir_spread": ((opts.max(1) - opts.min(1)) / bf).clip(0, 5),
        "gain": ((opts.max(1) - g_default) / bf).clip(0, 5),
        "egc_atp": g["egc_atp"],
        "egc_flag": (g["egc_atp"] > TOL).astype(int),
    })
    return feat


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("-k", type=int, default=6, help="number of clusters")
    ap.add_argument("--sample", type=int, default=0, help="subsample rows for KMeans fit")
    args = ap.parse_args(argv)

    feat = build_features()
    print(f"feature matrix: {len(feat):,} (model,reaction) rows over "
          f"{feat['model_id'].nunique()} growing models, {feat['base'].nunique()} unique base reactions")

    cols = ["off_frac", "fwd_frac", "rev_frac", "all_frac", "dir_spread", "gain", "egc_atp", "egc_flag"]
    X = StandardScaler().fit_transform(feat[cols].values)
    fit_X = X
    if args.sample and len(X) > args.sample:
        idx = np.random.default_rng(0).choice(len(X), args.sample, replace=False)
        fit_X = X[idx]
    km = KMeans(n_clusters=args.k, random_state=0, n_init=4).fit(fit_X)
    feat["cluster"] = km.predict(X)

    print(f"\nKMeans k={args.k} cluster profiles (mean feature values):")
    prof = feat.groupby("cluster")[cols].mean()
    prof["n"] = feat.groupby("cluster").size()
    prof["pct"] = (100 * prof["n"] / len(feat)).round(1)
    pd.set_option("display.width", 200, "display.max_columns", 20)
    print(prof.round(3).to_string())

    out = ROOT / "reaction_features.parquet"
    feat.to_parquet(out, compression="zstd", index=False)
    print(f"\nwrote {out} ({out.stat().st_size/1e6:.1f} MB) — feature matrix + cluster labels")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
