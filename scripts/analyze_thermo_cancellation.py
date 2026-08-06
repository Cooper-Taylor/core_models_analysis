#!/usr/bin/env python3
"""Test the proposed mechanism behind the concordant / discordant reaction sets.

Hypothesis: the families where all three sources agree are those where the
reaction DeltaG is a *large fraction* of the formation-energy turnover, and the
families where they disagree are group-transfer reactions where DeltaG is a
small difference between two large numbers, so every source's per-compound
error is amplified.

Cancellation ratio for a reaction:

    turnover  = sum_i |nu_i * dGf_i|        (Group Contribution compound energies)
    kappa     = |sum_i nu_i * dGf_i| / turnover = |dG_rxn| / turnover

kappa near 1 => no cancellation (the reaction really does change the total
formation energy by that much). kappa near 0 => near-total cancellation: the
two sides have almost the same formation energy and DeltaG is the residue.

CAVEAT: ``turnover`` depends on where the compound-energy scale puts its zero,
so kappa is a descriptive heuristic, not an invariant. The convention-free
companion measure is the propagated signal-to-noise ratio, built from
ModelSEED's own per-compound uncertainty ``deltagerr``:

    sigma_rxn = sqrt( sum_i (nu_i * err_i)^2 )
    snr       = |dG_rxn| / sigma_rxn

snr is unchanged by shifting every compound energy by a constant, and answers
the same question directly: is this reaction's DeltaG large or small compared
with the error already accumulated in the compound energies it is built from?

The prediction is that cross-source disagreement (median |dG_a - dG_b|,
expressed relative to the compound-energy scale of the reaction) rises as kappa
falls, and that the discordant families identified by the family scan sit at
low kappa.

Restricted throughout to reactions whose staged KEGG id is vouched for by a
ModelSEED alias, so the reaction-identity mis-mapping artifact does not
contaminate the chemistry signal.

Writes results/thermo_agreement/cancellation_by_family.tsv.
"""
from __future__ import annotations

import glob
import json
import os
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

MSDB_ROOT = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase"))
ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
OUT_DIR = ANALYSIS_DIR / "results" / "thermo_agreement"
BIOCHEM = MSDB_ROOT / "Biochemistry"


def main() -> None:
    cpd_dg: dict[str, float] = {}
    cpd_err: dict[str, float] = {}
    for path in sorted(glob.glob(str(BIOCHEM / "compound_*.json"))):
        for entry in json.load(open(path)):
            val = entry.get("deltag")
            # 10000000 is the ModelSEED sentinel for "no estimate"
            if isinstance(val, (int, float)) and abs(val) < 1e6:
                cpd_dg[entry["id"]] = float(val)
                err = entry.get("deltagerr")
                cpd_err[entry["id"]] = (float(err)
                                        if isinstance(err, (int, float)) and abs(err) < 1e6
                                        else 0.0)

    stoich: dict[str, dict[str, float]] = {}
    for path in sorted(glob.glob(str(BIOCHEM / "reaction_*.json"))):
        for entry in json.load(open(path)):
            vec: dict[str, float] = defaultdict(float)
            for item in entry.get("stoichiometry") or []:
                coeff = float(item.get("coefficient", 0) or 0)
                if coeff:
                    vec[item["compound"]] += coeff
            stoich[entry["id"]] = {k: v for k, v in vec.items() if v != 0}

    df = pd.read_csv(OUT_DIR / "reaction_features.tsv", sep="\t", low_memory=False)
    t = df[(df["n_sources"] == 3) & (df["kegg_vouched"] == 1)].copy()
    print(f"{len(t)} three-source reactions with a ModelSEED-vouched KEGG mapping")

    turnover, kappa, sigma, snr = [], [], [], []
    for rid in t["rxn"]:
        vec = stoich.get(rid, {})
        if not vec or any(c not in cpd_dg for c in vec):
            turnover.append(np.nan); kappa.append(np.nan)
            sigma.append(np.nan); snr.append(np.nan)
            continue
        terms = np.array([v * cpd_dg[c] for c, v in vec.items()], float)
        tot = float(np.abs(terms).sum())
        net = float(terms.sum())
        turnover.append(tot)
        kappa.append(abs(net) / tot if tot > 0 else np.nan)
        sig = float(np.sqrt(sum((v * cpd_err[c]) ** 2 for c, v in vec.items())))
        sigma.append(sig)
        snr.append(abs(net) / sig if sig > 0 else np.nan)
    t["turnover"] = turnover
    t["kappa"] = kappa
    t["sigma_rxn"] = sigma
    t["snr"] = snr
    t = t[np.isfinite(t["kappa"])].copy()
    print(f"{len(t)} with complete Group-Contribution compound energies")

    for a, b in (("gc", "eq"), ("gc", "dgp"), ("eq", "dgp")):
        t[f"absdiff_{a}_{b}"] = (t[f"dg_{a}"] - t[f"dg_{b}"]).abs()

    # ---- overall: does disagreement track cancellation?
    print("\nSpearman rho( kappa , cross-source |difference| ) -- negative means "
          "more cancellation goes with more disagreement:")
    for a, b in (("gc", "eq"), ("gc", "dgp"), ("eq", "dgp")):
        rho = stats.spearmanr(t["kappa"], t[f"absdiff_{a}_{b}"]).statistic
        print(f"  {a}-{b}: rho = {rho:+.3f}")

    print("\nSpearman rho( snr , cross-source |difference| ) -- the "
          "convention-free version of the same test:")
    ts = t[np.isfinite(t["snr"])]
    for a, b in (("gc", "eq"), ("gc", "dgp"), ("eq", "dgp")):
        rho = stats.spearmanr(ts["snr"], ts[f"absdiff_{a}_{b}"]).statistic
        print(f"  {a}-{b}: rho = {rho:+.3f}   (n = {len(ts)})")

    print("\nby propagated signal-to-noise decile (|dG| / sigma_rxn):")
    ts = ts.copy()
    ts["sbin"] = pd.qcut(ts["snr"], 10, labels=False, duplicates="drop")
    rows = []
    for k, grp in ts.groupby("sbin"):
        rows.append({
            "snr_decile": int(k) + 1,
            "snr_range": f"{grp['snr'].min():.2f}-{grp['snr'].max():.2f}",
            "n": len(grp),
            "median_abs_dg": float(grp["dg_gc"].abs().median()),
            "median_sigma": float(grp["sigma_rxn"].median()),
            "mad_gc_eq": float(grp["absdiff_gc_eq"].median()),
            "mad_gc_dgp": float(grp["absdiff_gc_dgp"].median()),
            "r_gc_dgp": float(np.corrcoef(grp["dg_gc"], grp["dg_dgp"])[0, 1]),
            "r_gc_eq": float(np.corrcoef(grp["dg_gc"], grp["dg_eq"])[0, 1]),
        })
    pd.set_option("display.width", 200)
    print(pd.DataFrame(rows).to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    print("\nby kappa decile:")
    t["kbin"] = pd.qcut(t["kappa"], 10, labels=False, duplicates="drop")
    rows = []
    for k, grp in t.groupby("kbin"):
        rows.append({
            "kappa_decile": int(k) + 1,
            "kappa_range": f"{grp['kappa'].min():.3f}-{grp['kappa'].max():.3f}",
            "n": len(grp),
            "median_abs_dg": float(grp["dg_gc"].abs().median()),
            "median_turnover": float(grp["turnover"].median()),
            "mad_gc_eq": float(grp["absdiff_gc_eq"].median()),
            "mad_gc_dgp": float(grp["absdiff_gc_dgp"].median()),
            "r_gc_dgp": float(np.corrcoef(grp["dg_gc"], grp["dg_dgp"])[0, 1]),
            "r_gc_eq": float(np.corrcoef(grp["dg_gc"], grp["dg_eq"])[0, 1]),
        })
    dec = pd.DataFrame(rows)
    pd.set_option("display.width", 200)
    print(dec.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    # ---- per family
    fams = {
        "redox: NAD(P)(H)": t[["cof_nad", "cof_nadh", "cof_nadp", "cof_nadph"]].max(axis=1) == 1,
        "EC 1 oxidoreductase": t["ec_class"].fillna("none").str.contains("1"),
        "aldehyde group change": t["d_aldehyde"] != 0,
        "ketone group change": t["d_ketone"] != 0,
        "small (<=6 C)": t["max_carbon"] <= 6,
        "1 group type changes": t["n_groups_changed"] <= 1,
        "phosphoryl transfer (ATP/ADP/AMP)": t[["cof_atp", "cof_adp", "cof_amp"]].max(axis=1) == 1,
        "phosphoanhydride change": t["d_phosphoanhydride"] != 0,
        "methyl transfer (SAM/SAH)": t[["cof_sam", "cof_sah"]].max(axis=1) == 1,
        "acyl transfer (thioester/CoA)": t["d_thioester"] != 0,
        "amide-N transfer (Gln/Glu)": t[["cof_gln", "cof_glu"]].max(axis=1) == 1,
        "amide group change": t["d_amide"] != 0,
        "EC 2 transferase": t["ec_class"].fillna("none").str.contains("2"),
        "EC 3 hydrolase": t["ec_class"].fillna("none").str.contains("3"),
        "EC 6 ligase": t["ec_class"].fillna("none").str.contains("6"),
    }
    rows = []
    for name, mask in fams.items():
        g = t[mask.to_numpy()]
        if len(g) < 30:
            continue
        rows.append({
            "family": name, "n": len(g),
            "median_kappa": float(g["kappa"].median()),
            "median_snr": float(g["snr"].median()),
            "median_sigma": float(g["sigma_rxn"].median()),
            "median_turnover": float(g["turnover"].median()),
            "median_abs_dg": float(g["dg_gc"].abs().median()),
            "mad_gc_eq": float(g["absdiff_gc_eq"].median()),
            "mad_gc_dgp": float(g["absdiff_gc_dgp"].median()),
            "mad_eq_dgp": float(g["absdiff_eq_dgp"].median()),
            "r_gc_eq": float(np.corrcoef(g["dg_gc"], g["dg_eq"])[0, 1]),
            "r_gc_dgp": float(np.corrcoef(g["dg_gc"], g["dg_dgp"])[0, 1]),
            "r_eq_dgp": float(np.corrcoef(g["dg_eq"], g["dg_dgp"])[0, 1]),
        })
    fam_df = pd.DataFrame(rows).sort_values("median_snr", ascending=False)
    fam_df.to_csv(OUT_DIR / "cancellation_by_family.tsv", sep="\t", index=False,
                  float_format="%.4f")
    print("\nper family, sorted by propagated signal-to-noise (high = dG is large "
          "relative to accumulated compound-energy error):")
    print(fam_df.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    for key in ("median_kappa", "median_snr"):
        r1 = stats.spearmanr(fam_df[key], fam_df["r_gc_dgp"]).statistic
        r2 = stats.spearmanr(fam_df[key], fam_df["r_gc_eq"]).statistic
        print(f"\nacross families: rho({key}, r_gc_dgp) = {r1:+.3f}"
              f"   rho({key}, r_gc_eq) = {r2:+.3f}")


if __name__ == "__main__":
    main()
