#!/usr/bin/env python3
"""HDBSCAN clustering of the 14,578 ModelSEED reactions that carry a DeltaG'deg
from all three thermodynamic sources, in the 3D space

    (Group Contribution, eQuilibrator, dGPredictor)   [kcal/mol]

and a chemical characterisation of whatever clusters fall out.

Two representations are clustered, because the answer depends on the metric and
saying so is part of the result:

  raw      -- Euclidean distance in kcal/mol. This is the space the scatter
              plots are drawn in, so clusters here are directly readable, but
              the distribution is heavy-tailed and a handful of reactions run to
              10^3-10^4 kcal/mol.
  symlog   -- sgn(x) * log10(1 + |x|) applied per axis. Preserves sign and
              order, compresses the tail, and gives near-zero reactions (the
              bulk of metabolism) room to separate.

For every cluster the script reports size, the centroid and spread on each
axis, the cross-source agreement inside it, and the chemical features that are
enriched relative to the whole three-source set (cofactors, EC class, KEGG
mapping provenance, element content, functional-group changes).

Writes results/thermo_agreement/clusters_3d_{raw,symlog}.tsv (per-reaction
labels) and cluster_profiles_{raw,symlog}.tsv (per-cluster profile).
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import HDBSCAN

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
DATA_DIR = ANALYSIS_DIR / "results" / "thermo_agreement"

AXES = ["dg_gc", "dg_eq", "dg_dgp"]
MIN_CLUSTER_SIZE = 30
MIN_SAMPLES = 10

# Binary chemistry features tested for enrichment inside each cluster.
FLAG_FEATURES = [
    "kegg_vouched", "is_transport", "is_obsolete", "all_have_smiles",
    "has_S", "has_P", "has_N", "has_halogen", "has_metal",
    "cof_atp", "cof_adp", "cof_amp", "cof_nad", "cof_nadh", "cof_nadp",
    "cof_nadph", "cof_fad", "cof_coa", "cof_o2", "cof_h2o2", "cof_co2",
    "cof_nh3", "cof_pi", "cof_ppi", "cof_sam", "cof_sah", "cof_glu",
    "cof_gln", "cof_thf", "cof_h2o", "cof_h",
    "d_phosphoanhydride_nz", "d_thioester_nz", "d_aldehyde_nz",
    "d_ketone_nz", "d_amide_nz", "d_carboxylate_nz", "d_hydroxyl_nz",
    "d_aromatic_c_nz", "d_alkene_nz",
    "status_ok", "multi_kegg", "kegg_heavily_reused", "generic_formula",
    "ec1", "ec2", "ec3", "ec4", "ec5", "ec6",
]


def symlog(a: np.ndarray) -> np.ndarray:
    return np.sign(a) * np.log10(1.0 + np.abs(a))


def add_flags(t: pd.DataFrame) -> pd.DataFrame:
    t = t.copy()
    for col in ("phosphoanhydride", "thioester", "aldehyde", "ketone", "amide",
                "carboxylate", "hydroxyl", "aromatic_c", "alkene"):
        t[f"d_{col}_nz"] = (t[f"d_{col}"] != 0).astype(int)
    t["status_ok"] = (t["status"] == "OK").astype(int)
    t["multi_kegg"] = (t["n_kegg_with_dg"] > 1).astype(int)
    t["kegg_heavily_reused"] = (t["kegg_max_reuse"] > 20).astype(int)
    t["generic_formula"] = (t["n_generic_formula"] > 0).astype(int)
    ecs = t["ec_class"].fillna("none").astype(str)
    for d in "123456":
        t[f"ec{d}"] = ecs.str.split(";").apply(lambda p, dd=d: int(dd in p))
    return t


def profile(t: pd.DataFrame, labels: np.ndarray, space: str) -> pd.DataFrame:
    base = {f: t[f].mean() for f in FLAG_FEATURES}
    rows = []
    for lab in sorted(set(labels)):
        sub = t[labels == lab]
        gc, eq, dgp = (sub[a].to_numpy(float) for a in AXES)
        rec = {
            "space": space,
            "cluster": int(lab),
            "label": "noise" if lab == -1 else f"C{lab}",
            "n": len(sub),
            "gc_median": float(np.median(gc)),
            "eq_median": float(np.median(eq)),
            "dgp_median": float(np.median(dgp)),
            "gc_iqr": float(np.subtract(*np.percentile(gc, [75, 25]))),
            "dgp_iqr": float(np.subtract(*np.percentile(dgp, [75, 25]))),
            "mad_gc_eq": float(np.median(np.abs(gc - eq))),
            "mad_gc_dgp": float(np.median(np.abs(gc - dgp))),
            "r_gc_dgp": (float(np.corrcoef(gc, dgp)[0, 1])
                         if len(sub) > 8 and gc.std() > 0 and dgp.std() > 0 else np.nan),
            "frac_kegg_vouched": float(sub["kegg_vouched"].mean()),
        }
        # top enriched binary features (ratio vs the whole three-source set)
        enr = []
        for f in FLAG_FEATURES:
            rate, b = sub[f].mean(), base[f]
            if rate >= 0.35 and rate >= 2.0 * max(b, 0.01) and sub[f].sum() >= 8:
                enr.append((rate / max(b, 0.01), f, rate))
        enr.sort(reverse=True)
        rec["enriched"] = "; ".join(f"{f} {rate:.0%} ({ratio:.1f}x)"
                                    for ratio, f, rate in enr[:6]) or "-"
        names = sub["name"].dropna().astype(str)
        rec["example_rxns"] = ";".join(sub["rxn"].head(4))
        rec["example_names"] = " | ".join(n[:38] for n in names.head(3))
        rec["top_ec_class"] = (sub["ec_class"].mode().iat[0]
                               if not sub["ec_class"].dropna().empty else "")
        rows.append(rec)
    return pd.DataFrame(rows).sort_values("n", ascending=False)


def run_space(t: pd.DataFrame, space: str) -> tuple[np.ndarray, pd.DataFrame]:
    X = t[AXES].to_numpy(float)
    if space == "symlog":
        X = symlog(X)
    # Standardise so no single axis dominates the metric purely through spread.
    X = (X - X.mean(axis=0)) / X.std(axis=0)
    model = HDBSCAN(min_cluster_size=MIN_CLUSTER_SIZE, min_samples=MIN_SAMPLES,
                    cluster_selection_method="eom")
    labels = model.fit_predict(X)
    n_clust = len(set(labels)) - (1 if -1 in labels else 0)
    noise = float((labels == -1).mean())
    print(f"\n=== {space}: {n_clust} clusters, noise {noise:.1%} "
          f"({int((labels == -1).sum())} reactions) ===")
    return labels, profile(t, labels, space)


def main() -> None:
    df = pd.read_csv(DATA_DIR / "reaction_features.tsv", sep="\t", low_memory=False)
    t = add_flags(df[df["n_sources"] == 3].reset_index(drop=True))
    print(f"{len(t)} reactions with all three sources")

    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 90)
    for space in ("raw", "symlog"):
        labels, prof = run_space(t, space)
        out = t[["rxn", "name", "ec_class"] + AXES].copy()
        out["cluster"] = labels
        out.to_csv(DATA_DIR / f"clusters_3d_{space}.tsv", sep="\t",
                   index=False, float_format="%.3f")
        prof.to_csv(DATA_DIR / f"cluster_profiles_{space}.tsv", sep="\t",
                    index=False, float_format="%.4f")
        show = ["label", "n", "gc_median", "eq_median", "dgp_median",
                "mad_gc_dgp", "r_gc_dgp", "frac_kegg_vouched", "enriched"]
        print(prof[show].head(22).to_string(index=False))


if __name__ == "__main__":
    main()
