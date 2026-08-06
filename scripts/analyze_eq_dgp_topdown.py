#!/usr/bin/env python3
"""Top-down eQuilibrator vs dGPredictor-ModelSEED: reaction chemistry first,
then enzymes, then metabolites.

Layer 1  ORGANIC TRANSFORMATION CLASS. Every key-subset reaction is assigned one
         chemistry class (organic_reaction_types.py) and we ask which classes
         carry the disagreement, relative to how common they are.

Layer 2  ENZYME. Inside the worst classes, which EC subfamilies and enzyme names
         concentrate?

Layer 3  METABOLITE. Which compounds do the two methods value differently --
         reported ONLY with the model-free validation, because raw fitted
         offsets are not identifiable (see gauge note below).

On the two quantities that need defining
----------------------------------------
sigma  = dGPredictor-ModelSEED's OWN reported uncertainty for that reaction. It
         is the BayesianRidge posterior standard deviation, staged in kJ/mol as
         ``dG_uncer`` and stored as element [1] of the reaction's
         ``thermodynamics["dGPredictor-ModelSEED"]`` triple in kcal/mol. Nothing
         here computes it; we only test whether it is informative.

fitted offset
       = NOT a regression of one method's dG on the other. For every reaction r
         take the disagreement d_r = dG_eq(r) - dG_dgp(r) and solve, by least
         squares over the stoichiometric matrix,
             d_r  ~  sum_i  nu_ir * x_i
         where nu_ir is compound i's coefficient in reaction r. x_i is the
         per-compound offset. This is legitimate because both methods are
         additive over compound formation energies, so their difference is too.

THE GAUGE PROBLEM. The stoichiometric matrix S has a non-trivial null space: any
z with S z = 0 can be added to a solution x without changing a single predicted
reaction value. Element conservation supplies such z for free -- adding a fixed
amount per carbon atom to every compound leaves every balanced reaction exactly
unchanged (verified numerically in the output below). So an individual x_i is
NOT identifiable, and a large fitted offset can be pure bookkeeping. Every
compound claim here is therefore validated against a model-free quantity: the
observed median |dG_eq - dG_dgp| over the reactions that actually contain that
compound, compared with the subset-wide baseline.

Outputs (results/eq_vs_dgpms/):
  reaction_class_breakdown.tsv   layer 1
  enzyme_breakdown.tsv           layer 2
  metabolite_validated.tsv       layer 3, with the model-free check
  gauge_demo.tsv                 the null-space demonstration
"""
from __future__ import annotations

import glob
import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.sparse.linalg import lsqr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from organic_reaction_types import classify, CLASS_ORDER  # noqa: E402

MSDB_ROOT = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/tmp/devsnap"))
ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
OUT_DIR = Path(os.environ.get("EQDGP_OUT",
                              str(ANALYSIS_DIR / "results" / "eq_vs_dgpms")))
BIOCHEM = MSDB_ROOT / "Biochemistry"
RNG = np.random.default_rng(20260806)
# Threshold for calling a reaction "discordant". This is a CHOSEN round number
# (~4x the subset baseline of 3.44), not a derived one -- it isolates a clear
# tail without being sensitive to where exactly it sits: >4 -> 45.9% of
# reactions, >10 -> 26.4%, >15 -> 17.6%, >20 -> 12.3%, >30 -> 8.4%, and the
# class ranking is stable across that range. For reference the cascade's own
# reversible band is +/-2.0 kcal/mol on mMdeltaG (reversibility_heuristics.py
# line 327), so a disagreement well below 15 can already flip a direction call;
# 15 is deliberately conservative.
DISCORDANT = 15.0   # kcal/mol


def load(pattern):
    out = {}
    for p in sorted(glob.glob(str(BIOCHEM / pattern))):
        for e in json.load(open(p)):
            out[e["id"]] = e
    return out


def main() -> None:
    key = pd.read_csv(OUT_DIR / "key_subset.tsv", sep="\t", low_memory=False)
    key["absdiff"] = (key["dg_eq"] - key["dg_dgp"]).abs()
    key["diff"] = key["dg_eq"] - key["dg_dgp"]
    rxns, cpds = load("reaction_*.json"), load("compound_*.json")
    stoich = {r: {i["compound"]: float(i.get("coefficient", 0) or 0)
                  for i in (rxns.get(r, {}).get("stoichiometry") or [])}
              for r in key["rxn"]}
    baseline = key["absdiff"].median()
    print(f"key subset n = {len(key)},  baseline median |eQ - dGP| = {baseline:.2f} kcal/mol")
    print(f"discordant threshold = {DISCORDANT:g} kcal/mol "
          f"({(key['absdiff'] > DISCORDANT).sum()} reactions, "
          f"{(key['absdiff'] > DISCORDANT).mean():.1%})\n")

    # ------------------------------------------------ layer 1: chemistry class
    labels = [classify(r, stoich[r["rxn"]], cpds) for _, r in key.iterrows()]
    key["chem_class"] = [a for a, _ in labels]
    key["chem_subclass"] = [b for _, b in labels]

    rows = []
    for cls in CLASS_ORDER:
        sub = key[key["chem_class"] == cls]
        if len(sub) < 15:
            continue
        rows.append({
            "class": cls, "n": len(sub),
            "share_of_subset": len(sub) / len(key),
            "median_absdiff": sub["absdiff"].median(),
            "vs_baseline": sub["absdiff"].median() / baseline,
            "frac_discordant": (sub["absdiff"] > DISCORDANT).mean(),
            "enrichment_in_discordant": ((sub["absdiff"] > DISCORDANT).mean()
                                         / (key["absdiff"] > DISCORDANT).mean()),
            "median_abs_dg_eq": sub["dg_eq"].abs().median(),
            "rel_error": sub["absdiff"].median() / max(sub["dg_eq"].abs().median(), 1e-9),
            "median_sigma": sub["dgp_uncertainty"].median(),
            "r": np.corrcoef(sub["dg_eq"], sub["dg_dgp"])[0, 1] if len(sub) > 8 else np.nan,
        })
    cls_df = pd.DataFrame(rows).sort_values("median_absdiff", ascending=False)
    cls_df.to_csv(OUT_DIR / "reaction_class_breakdown.tsv", sep="\t", index=False,
                  float_format="%.4f")
    pd.set_option("display.width", 230)
    print("=== LAYER 1: organic transformation class, worst first ===")
    print(cls_df[["class", "n", "median_absdiff", "vs_baseline", "rel_error",
                  "frac_discordant", "enrichment_in_discordant", "median_sigma", "r"]]
          .to_string(index=False, float_format=lambda v: f"{v:.2f}"))

    print("\n  subclass detail for the three worst classes:")
    for cls in cls_df["class"].head(3):
        for sc, g in key[key["chem_class"] == cls].groupby("chem_subclass"):
            if len(g) >= 10:
                print(f"    {cls[:34]:36s} | {sc[:32]:34s} n={len(g):5d} "
                      f"median|d|={g['absdiff'].median():7.2f}")

    # ------------------------------------------------------- layer 2: enzymes
    disc = key[key["absdiff"] > DISCORDANT]
    rows = []
    for ec3, g in key.assign(
            ec3=key["ec"].fillna("").astype(str).str.split(";").str[0]
            .str.split(".").str[:3].str.join(".")).groupby("ec3"):
        if len(g) < 20 or ec3 in ("", "nan"):
            continue
        rows.append({"ec3": ec3, "n": len(g),
                     "median_absdiff": g["absdiff"].median(),
                     "frac_discordant": (g["absdiff"] > DISCORDANT).mean(),
                     "example_enzyme": g["name"].iloc[0][:52],
                     "dominant_class": g["chem_class"].mode().iat[0]})
    enz = pd.DataFrame(rows).sort_values("median_absdiff", ascending=False)
    enz.to_csv(OUT_DIR / "enzyme_breakdown.tsv", sep="\t", index=False,
               float_format="%.4f")
    print("\n=== LAYER 2: EC subfamilies (>=20 reactions), worst first ===")
    print(enz.head(14).to_string(index=False, float_format=lambda v: f"{v:.2f}"))

    # ------------------------------------- gauge demonstration (for the report)
    ids = list(key["rxn"])
    all_c = sorted({c for r in ids for c in stoich[r]})
    idx = {c: i for i, c in enumerate(all_c)}
    rr, cc, vv = [], [], []
    for i, r in enumerate(ids):
        for c, v in stoich[r].items():
            rr.append(i); cc.append(idx[c]); vv.append(v)
    S = sparse.csr_matrix((vv, (rr, cc)), shape=(len(ids), len(all_c)))
    import re as _re
    gauge_rows = []
    for el in ("C", "N", "O", "P", "S"):
        z = np.zeros(len(all_c))
        for c, i in idx.items():
            m = _re.search(rf"{el}(\d*)(?![a-z])", str(cpds.get(c, {}).get("formula") or ""))
            if m:
                z[i] = int(m.group(1) or 1)
        if z.any():
            sz = np.abs(S @ z)
            gauge_rows.append({"element": el, "max_abs_change_kcal": sz.max(),
                               "frac_reactions_unchanged": float((sz < 1e-9).mean())})
    gd = pd.DataFrame(gauge_rows)
    gd.to_csv(OUT_DIR / "gauge_demo.tsv", sep="\t", index=False, float_format="%.6g")
    print("\n=== GAUGE: adding a fixed amount per atom of an element changes nothing ===")
    print(gd.to_string(index=False))

    # ---------------------------------------------------- layer 3: metabolites
    counts = Counter(c for r in ids for c in stoich[r])
    keep = sorted(c for c, n in counts.items() if n >= 5)
    kidx = {c: i for i, c in enumerate(keep)}
    usable = [r for r in ids if stoich[r] and all(c in kidx for c in stoich[r])]
    rr, cc, vv = [], [], []
    for i, r in enumerate(usable):
        for c, v in stoich[r].items():
            rr.append(i); cc.append(kidx[c]); vv.append(v)
    M = sparse.csr_matrix((vv, (rr, cc)), shape=(len(usable), len(keep)))
    tgt = key.set_index("rxn").loc[usable, "diff"].to_numpy(float)
    sol = lsqr(M, tgt, damp=1e-2, iter_lim=8000)[0]

    by_rxn = key.set_index("rxn")
    rows = []
    for c, i in kidx.items():
        if counts[c] < 15:
            continue
        members = [r for r in ids if c in stoich[r]]
        obs = by_rxn.loc[members, "absdiff"].median()
        rows.append({
            "compound": c, "name": cpds.get(c, {}).get("name", ""),
            "formula": cpds.get(c, {}).get("formula", ""),
            "n_reactions": len(members),
            "fitted_offset": sol[i],
            "observed_median_absdiff": obs,
            "ratio_vs_baseline": obs / baseline,
            "verdict": ("REAL" if obs / baseline >= 2 else
                        "gauge artifact" if abs(sol[i]) > 15 and obs / baseline < 1.5
                        else "unremarkable"),
            "median_sigma": by_rxn.loc[members, "dgp_uncertainty"].median(),
        })
    met = pd.DataFrame(rows).sort_values("observed_median_absdiff", ascending=False)
    met.to_csv(OUT_DIR / "metabolite_validated.tsv", sep="\t", index=False,
               float_format="%.4f")
    print("\n=== LAYER 3: metabolites, ranked by OBSERVED disagreement (not fitted) ===")
    print(met.head(16)[["name", "formula", "n_reactions", "fitted_offset",
                        "observed_median_absdiff", "ratio_vs_baseline",
                        "median_sigma", "verdict"]]
          .to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    print("\n  largest FITTED offsets that are NOT real (the gauge trap):")
    trap = met[met["verdict"] == "gauge artifact"].reindex(
        met["fitted_offset"].abs().sort_values(ascending=False).index).dropna(how="all")
    print(trap.head(8)[["name", "n_reactions", "fitted_offset",
                        "observed_median_absdiff", "ratio_vs_baseline"]]
          .to_string(index=False, float_format=lambda v: f"{v:.2f}"))

    key.to_csv(OUT_DIR / "key_subset_classified.tsv", sep="\t", index=False,
               float_format="%.4f")
    print(f"\nwrote outputs to {OUT_DIR}")


if __name__ == "__main__":
    main()
