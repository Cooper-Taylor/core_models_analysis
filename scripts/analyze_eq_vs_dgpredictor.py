#!/usr/bin/env python3
"""eQuilibrator vs dGPredictor, reconciled on ModelSEED reaction identity.

The two sources reach a ModelSEED reaction by completely different routes, and
the comparison is only meaningful once both are pinned to the *same* reaction:

  eQuilibrator  ModelSEED cpd --InChIKey--> MetaNetX id, then the ModelSEED
                stoichiometry is rebuilt in MetaNetX ids and handed to
                ComponentContribution.standard_dg_prime() at pH 7.0, I = 0.25 M,
                298.15 K (Retrieve_eQuilibrator_Reactions_Energies.py). The
                *reaction* is ModelSEED's own; only the compound identities are
                translated. Reactions are computed only when EVERY reagent maps
                ("EQC"); partial coverage ("EQP") yields nothing.

  dGPredictor   ModelSEED rxn --> a KEGG reaction id --> dGPredictor's own
                prediction for the KEGG reaction. The reaction itself is KEGG's,
                so a wrong id silently substitutes different chemistry. That is
                what build_dgpredictor_kegg_mask.py withholds.

So the two need different audits. dGPredictor's is the KEGG mask (already
applied upstream of this script, via reaction_features.tsv). eQuilibrator's is
here: its InChIKey match has three fallback tiers, and the two loose ones lose
information that matters thermodynamically.

  tier 1  full InChIKey                    exact
  tier 2  first two blocks                 protonation-state-blind
  tier 3  first block (connectivity) only  STEREO-BLIND -- conflates anomers,
                                           D/L pairs, cis/trans

Tier 3 makes stereoisomers indistinguishable, so a ModelSEED reaction that
interconverts them is handed to eQuilibrator as A = A. Worse, the retrieval
script accumulates each side into a dict keyed by MetaNetX id
(``lhs[mnx_id] = |coeff|``), so when two ModelSEED compounds in one reaction
collapse onto the same MetaNetX id the second silently OVERWRITES the first
rather than summing -- the reaction eQuilibrator scored is then not the
reaction ModelSEED wrote.

This script reconstructs that map, flags both failure modes, defines the
high-confidence subset where each source is on firm ground, and compares them.

Outputs (results/eq_vs_dgp/):
  reconciliation.tsv        per-reaction join + provenance flags for both sources
  key_subset.tsv            the high-confidence comparison set
  concordant.tsv / discordant.tsv
  compound_offsets.tsv      per-metabolite eQ-minus-dGP offsets (least squares)
  class_breakdown.tsv       agreement by reaction class
  mechanism_tests.tsv       the four method-level hypotheses, tested

Requires the ``eq3`` env (rdkit) only for the optional structure checks; runs
without it otherwise.
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

MSDB_ROOT = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase"))
ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
BIOCHEM = MSDB_ROOT / "Biochemistry"
FEATURES = Path(os.environ.get(
    "EQDGP_FEATURES",
    str(ANALYSIS_DIR / "results" / "thermo_agreement" / "reaction_features.tsv")))
OUT_DIR = Path(os.environ.get("EQDGP_OUT", str(ANALYSIS_DIR / "results" / "eq_vs_dgp")))
# Which stored dGPredictor record fills dg_dgp in FEATURES. Selects the KEGG-id
# filters (only meaningful for the original) and the uncertainty tiering (only
# meaningful for the retrain, whose reported error is calibrated).
DGP_LABEL = os.environ.get("EQDGP_DGP_LABEL", "dGPredictor")
IS_RETRAIN = DGP_LABEL == "dGPredictor-ModelSEED"

EQ_TBL = BIOCHEM / "Thermodynamics" / "eQuilibrator" / "MetaNetX_Reaction_Energies.tbl"
EQ_CPD_TBL = BIOCHEM / "Thermodynamics" / "eQuilibrator" / "MetaNetX_Compound_Energies.tbl"
MNX_STRUCT = BIOCHEM / "Structures" / "MetaNetX" / "Structures_in_ModelSEED_and_eQuilibrator.txt"

RNG = np.random.default_rng(20260806)

# Both sources are reported in kcal/mol in ModelSEED. Method-level constants
# that differ between them and that the mechanism tests below probe:
EQ_IONIC_STRENGTH = 0.25   # M, set in Retrieve_eQuilibrator_Reactions_Energies.py
DGP_IONIC_STRENGTH = 0.10  # M, hardcoded in dGPredictor's get_dG0_prime()
PROTON = "cpd00067"

# eQuilibrator marks compounds it cannot estimate by inflating their variance
# (the ``1e6 * sigmas_inf @ sigmas_inf.T`` term in standard_dg_formation), which
# propagates to a reaction uncertainty three orders of magnitude above any real
# one. Across the 19,510 stored eQuilibrator reaction records the distribution
# is strictly bimodal: genuine uncertainties top out at 64.3 kcal/mol, sentinels
# start at 1,000, and NOTHING falls between 100 and 1,000. So this cut sits in an
# empty gap rather than being a tuned threshold. 2,748 records (14.1%) are
# sentinels -- eQuilibrator declaring it has no estimate. ModelSEED stores the
# uncertainty faithfully but nothing downstream reads it, so those energies are
# currently used as though they were measurements.
EQ_SENTINEL_UNCERTAINTY = 100.0  # kcal/mol


# --------------------------------------------------------------- loading MSDB
def load_json_dir(pattern: str) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for path in sorted(glob.glob(str(BIOCHEM / pattern))):
        for entry in json.load(open(path)):
            out[entry["id"]] = entry
    return out


def build_seed_to_mnx() -> tuple[dict[str, str], dict[str, int]]:
    """Reproduce Retrieve_eQuilibrator_Reactions_Energies.py's compound map.

    Returns (cpd -> mnx, cpd -> tier) where tier is 1/2/3 as documented above.
    """
    problem: set[str] = set()
    for line in open(EQ_CPD_TBL):
        parts = line.rstrip("\n").split("\t")
        if len(parts) > 1 and ("energy" in parts[1] or parts[1] == "nan"):
            problem.add(parts[0])

    full: dict[str, str] = {}
    two: dict[str, str] = {}
    one: dict[str, str] = {}
    for line in open(MNX_STRUCT):
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 2:
            continue
        mnx, key = parts
        if mnx in problem:
            continue
        full.setdefault(key, mnx)
        two.setdefault("-".join(key.split("-")[:2]), mnx)
        one.setdefault(key.split("-")[0], mnx)

    cpds = load_json_dir("compound_*.json")
    seed_mnx: dict[str, str] = {}
    tier: dict[str, int] = {}
    for cid, entry in cpds.items():
        key = entry.get("inchikey")
        if not key:
            continue
        if key in full:
            seed_mnx[cid], tier[cid] = full[key], 1
        elif "-".join(key.split("-")[:2]) in two:
            seed_mnx[cid], tier[cid] = two["-".join(key.split("-")[:2])], 2
        elif key.split("-")[0] in one:
            seed_mnx[cid], tier[cid] = one[key.split("-")[0]], 3
    return seed_mnx, tier


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("loading ModelSEED + feature table ...")
    feat = pd.read_csv(FEATURES, sep="\t", low_memory=False)
    rxns = load_json_dir("reaction_*.json")
    cpds = load_json_dir("compound_*.json")

    # ---------------------------------------------------------- (2) reconcile
    both = feat[feat["dg_eq"].notna() & feat["dg_dgp"].notna()].copy()
    print(f"\nreactions with BOTH eQuilibrator and dGPredictor (post KEGG mask): {len(both)}")

    seed_mnx, tier = build_seed_to_mnx()
    print(f"  ModelSEED compounds mapped to MetaNetX: {len(seed_mnx)}"
          f"  (tier1 {sum(1 for v in tier.values() if v == 1)},"
          f" tier2 {sum(1 for v in tier.values() if v == 2)},"
          f" tier3 {sum(1 for v in tier.values() if v == 3)})")

    eq_err: dict[str, float] = {}
    for line in open(EQ_TBL):
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 3:
            try:
                eq_err[parts[0]] = float(parts[2])
            except ValueError:
                continue

    # dGPredictor's own reported uncertainty. For the retrain this is calibrated
    # (rho = +0.67 against |eQ - dGP|), so it is a usable quality signal in its
    # own right; for the original it is uniformly ~0.35 and uninformative.
    dgp_err = {}
    for rid, entry in rxns.items():
        t = (entry.get("thermodynamics") or {}).get(DGP_LABEL)
        if t and len(t) > 1 and isinstance(t[1], (int, float)):
            dgp_err[rid] = abs(float(t[1]))

    rows = []
    for rxn_id in both["rxn"]:
        entry = rxns.get(rxn_id, {})
        stoich = entry.get("stoichiometry") or []
        tiers, mnx_seen, collide = [], defaultdict(list), False
        n_proton = 0.0
        charge_sub = charge_prod = 0.0
        for item in stoich:
            cid = item["compound"]
            coeff = float(item.get("coefficient", 0) or 0)
            if cid == PROTON:
                n_proton += coeff
            ch = cpds.get(cid, {}).get("charge")
            if isinstance(ch, (int, float)):
                if coeff < 0:
                    charge_sub += abs(coeff) * ch * ch
                else:
                    charge_prod += coeff * ch * ch
            if cid in seed_mnx:
                tiers.append(tier[cid])
                mnx_seen[seed_mnx[cid]].append((cid, coeff))
        # collision = two DIFFERENT ModelSEED compounds sharing one MetaNetX id
        for mnx, members in mnx_seen.items():
            if len({m[0] for m in members}) > 1:
                collide = True
        rows.append({
            "rxn": rxn_id,
            "eq_worst_tier": max(tiers) if tiers else 0,
            "eq_all_tier1": int(bool(tiers) and max(tiers) == 1),
            "eq_stereo_blind": int(bool(tiers) and max(tiers) == 3),
            "eq_mnx_collision": int(collide),
            "net_proton": n_proton,
            "d_charge_sq": charge_prod - charge_sub,
            "eq_uncertainty": eq_err.get(rxn_id, np.nan),
            "eq_sentinel": int(eq_err.get(rxn_id, 0.0) > EQ_SENTINEL_UNCERTAINTY),
            "dgp_uncertainty": dgp_err.get(rxn_id, np.nan),
        })
    prov = pd.DataFrame(rows)
    both = both.merge(prov, on="rxn", how="left")

    print("\neQuilibrator compound-mapping provenance across those reactions:")
    print(f"  all reagents matched on the full InChIKey (tier 1) : "
          f"{int(both['eq_all_tier1'].sum()):6d}")
    print(f"  worst reagent matched protonation-blind (tier 2)   : "
          f"{int((both['eq_worst_tier'] == 2).sum()):6d}")
    print(f"  worst reagent matched STEREO-BLIND (tier 3)        : "
          f"{int((both['eq_worst_tier'] == 3).sum()):6d}")
    print(f"  two ModelSEED compounds collapse to one MetaNetX id: "
          f"{int(both['eq_mnx_collision'].sum()):6d}  <-- overwritten in the retrieval dict")
    print(f"  eQuilibrator declares NO estimate (sentinel uncertainty > "
          f"{EQ_SENTINEL_UNCERTAINTY:g} kcal/mol): {int(both['eq_sentinel'].sum()):6d}")

    both.to_csv(OUT_DIR / "reconciliation.tsv", sep="\t", index=False,
                float_format="%.4f")

    # -------------------------------------------------------- (3) key subset
    key = both[
        (both["eq_all_tier1"] == 1)
        & (both["eq_mnx_collision"] == 0)
        & (both["status"] == "OK")
        & (both["is_transport"] == 0)
        & (both["all_have_smiles"] == 1)
        & (both["n_generic_formula"] == 0)
        & (both["eq_sentinel"] == 0)
    ].copy()
    if not IS_RETRAIN:
        # Only the KEGG-keyed original can be corrupted by multi-id averaging.
        key = key[key["n_kegg_with_dg"] == 1].copy()
    print(f"\nKEY SUBSET (both sources on firm ground): {len(key)} reactions")
    for label, mask in [
        ("eQ tier-1 mapping only", both["eq_all_tier1"] == 1),
        ("no MetaNetX collision", both["eq_mnx_collision"] == 0),
        ("status OK", both["status"] == "OK"),
        ("non-transport", both["is_transport"] == 0),
        ("all reagents have a structure", both["all_have_smiles"] == 1),
        ("no generic R/polymer formula", both["n_generic_formula"] == 0),
        ("eQuilibrator has a real estimate", both["eq_sentinel"] == 0),
    ] + ([] if IS_RETRAIN else [("exactly one KEGG id", both["n_kegg_with_dg"] == 1)]):
        print(f"    {label:34s} keeps {int(mask.sum()):6d}")

    key["diff"] = key["dg_eq"] - key["dg_dgp"]
    key["absdiff"] = key["diff"].abs()
    key.to_csv(OUT_DIR / "key_subset.tsv", sep="\t", index=False, float_format="%.4f")

    # ------------------------------------------------------------- (4) compare
    x, y = key["dg_eq"].to_numpy(float), key["dg_dgp"].to_numpy(float)
    print(f"\n--- eQuilibrator vs dGPredictor on the key subset (n = {len(key)}) ---")
    print(f"  Pearson r        {np.corrcoef(x, y)[0, 1]:+.3f}")
    print(f"  Spearman rho     {stats.spearmanr(x, y).statistic:+.3f}")
    print(f"  median |diff|    {np.median(np.abs(x - y)):.2f} kcal/mol")
    print(f"  median  diff     {np.median(x - y):+.2f} kcal/mol   (eQ minus dGP)")
    print(f"  IQR of diff      {np.percentile(x - y, 75) - np.percentile(x - y, 25):.2f}")
    for thr in (1, 2, 5, 10, 30):
        print(f"  within {thr:2d} kcal/mol: {(key['absdiff'] <= thr).mean():6.1%}")

    if IS_RETRAIN:
        print("\n--- agreement by the retrain's OWN reported uncertainty ---")
        tiers = [("high conf   (u <= 3)", key["dgp_uncertainty"] <= 3),
                 ("medium     (3 < u <= 30)", key["dgp_uncertainty"].between(3, 30, "right")),
                 ("low        (u > 30)", key["dgp_uncertainty"] > 30)]
        for label, m in tiers:
            sub = key[m.to_numpy()]
            if len(sub) < 10:
                continue
            xa, ya = sub["dg_eq"].to_numpy(float), sub["dg_dgp"].to_numpy(float)
            print(f"  {label:26s} n={len(sub):6d}  r={np.corrcoef(xa, ya)[0, 1]:6.3f}  "
                  f"median|d|={np.median(np.abs(xa - ya)):6.2f}  "
                  f"within2={np.mean(np.abs(xa - ya) <= 2):6.1%}")

    conc = key[key["absdiff"] <= 2].sort_values("absdiff")
    disc = key[key["absdiff"] > 15].sort_values("absdiff", ascending=False)
    print(f"\n  concordant (<= 2 kcal/mol): {len(conc)}")
    print(f"  discordant (> 15 kcal/mol): {len(disc)}")
    cols = ["rxn", "name", "ec", "dg_eq", "dg_dgp", "diff", "absdiff", "definition",
            "cofactors", "net_proton", "d_charge_sq", "n_participants", "max_carbon",
            "dgp_uncertainty", "eq_uncertainty"]
    conc[cols].to_csv(OUT_DIR / "concordant.tsv", sep="\t", index=False, float_format="%.3f")
    disc[cols].to_csv(OUT_DIR / "discordant.tsv", sep="\t", index=False, float_format="%.3f")

    # ------------------------------------- (5a) class breakdown
    classes = {
        "EC 1 oxidoreductase": key["ec_class"].fillna("none").str.contains("1"),
        "EC 2 transferase": key["ec_class"].fillna("none").str.contains("2"),
        "EC 3 hydrolase": key["ec_class"].fillna("none").str.contains("3"),
        "EC 4 lyase": key["ec_class"].fillna("none").str.contains("4"),
        "EC 5 isomerase": key["ec_class"].fillna("none").str.contains("5"),
        "EC 6 ligase": key["ec_class"].fillna("none").str.contains("6"),
        "NAD(P)(H) redox": key[["cof_nad", "cof_nadh", "cof_nadp", "cof_nadph"]].max(axis=1) == 1,
        "O2 involved": key["cof_o2"] == 1,
        "ATP/ADP/AMP": key[["cof_atp", "cof_adp", "cof_amp"]].max(axis=1) == 1,
        "phosphoanhydride change": key["d_phosphoanhydride"] != 0,
        "SAM/SAH methyl transfer": key[["cof_sam", "cof_sah"]].max(axis=1) == 1,
        "CoA thioester change": key["d_thioester"] != 0,
        "CO2": key["cof_co2"] == 1,
        "net proton != 0": key["net_proton"] != 0,
        "net proton == 0": key["net_proton"] == 0,
        "aromatic ring present": key["total_arom_rings"] > 0,
        "contains S": key["has_S"] == 1,
        "contains P": key["has_P"] == 1,
    }
    crows = []
    for name, mask in classes.items():
        sub = key[mask.to_numpy()]
        if len(sub) < 20:
            continue
        xa, ya = sub["dg_eq"].to_numpy(float), sub["dg_dgp"].to_numpy(float)
        crows.append({
            "class": name, "n": len(sub),
            "r": float(np.corrcoef(xa, ya)[0, 1]),
            "rho": float(stats.spearmanr(xa, ya).statistic),
            "median_absdiff": float(np.median(np.abs(xa - ya))),
            "median_diff": float(np.median(xa - ya)),
            "frac_within_2": float((np.abs(xa - ya) <= 2).mean()),
            "frac_over_15": float((np.abs(xa - ya) > 15).mean()),
            "median_abs_dg_eq": float(np.median(np.abs(xa))),
        })
    cls = pd.DataFrame(crows).sort_values("median_absdiff")
    cls.to_csv(OUT_DIR / "class_breakdown.tsv", sep="\t", index=False, float_format="%.4f")
    pd.set_option("display.width", 200)
    print("\n--- agreement by reaction class (sorted best to worst) ---")
    print(cls.to_string(index=False))

    # ------------------------------------- (5b) per-metabolite attribution
    #
    # eQuilibrator's reaction value is exactly additive over its compound
    # formation energies. If dGPredictor were too, (eQ - dGP) would itself be a
    # sum of per-compound offsets. Fitting that by least squares over the
    # stoichiometry matrix therefore asks: which metabolites carry the
    # disagreement? A compound with a large fitted offset and many supporting
    # reactions is one the two methods value differently.
    print("\n--- per-metabolite attribution of (eQ - dGP) ---")
    ids = list(key["rxn"])
    stoich = {}
    for rxn_id in ids:
        vec = defaultdict(float)
        for item in (rxns.get(rxn_id, {}).get("stoichiometry") or []):
            c = float(item.get("coefficient", 0) or 0)
            if c:
                vec[item["compound"]] += c
        stoich[rxn_id] = {k: v for k, v in vec.items() if v != 0}
    counts = Counter(c for r in ids for c in stoich[r])
    keep_cpd = sorted(c for c, n in counts.items() if n >= 5)
    idx = {c: i for i, c in enumerate(keep_cpd)}
    usable = [r for r in ids if stoich[r] and all(c in idx for c in stoich[r])]
    print(f"  fitting {len(usable)} reactions over {len(keep_cpd)} compounds (>=5 reactions each)")

    rr, cc, vv = [], [], []
    for i, r in enumerate(usable):
        for c, v in stoich[r].items():
            rr.append(i); cc.append(idx[c]); vv.append(v)
    S = sparse.csr_matrix((vv, (rr, cc)), shape=(len(usable), len(keep_cpd)))
    target = key.set_index("rxn").loc[usable, "diff"].to_numpy(float)
    perm = RNG.permutation(len(usable))
    n_test = max(1, len(usable) // 5)
    test_i, train_i = perm[:n_test], perm[n_test:]
    sol = lsqr(S[train_i], target[train_i], damp=1e-2, iter_lim=8000)[0]
    pred = S[test_i] @ sol
    resid = target[test_i] - pred
    ss_tot = float(((target[test_i] - target[test_i].mean()) ** 2).sum())
    print(f"  held-out R^2 of an additive per-compound model of the disagreement: "
          f"{1 - float(resid @ resid) / ss_tot:.3f}")
    print("  (high => the disagreement is per-metabolite and systematic;"
          " low => it is reaction-specific)")

    off = pd.DataFrame({
        "compound": keep_cpd,
        "name": [cpds.get(c, {}).get("name", "") for c in keep_cpd],
        "formula": [cpds.get(c, {}).get("formula", "") for c in keep_cpd],
        "charge": [cpds.get(c, {}).get("charge") for c in keep_cpd],
        "n_reactions": [counts[c] for c in keep_cpd],
        "offset_eq_minus_dgp": sol,
    }).sort_values("offset_eq_minus_dgp", key=np.abs, ascending=False)
    off.to_csv(OUT_DIR / "compound_offsets.tsv", sep="\t", index=False, float_format="%.3f")
    print("\n  metabolites the two methods value most differently "
          "(fitted offset, kcal/mol, >=15 reactions):")
    print(off[off["n_reactions"] >= 15].head(18).to_string(index=False))

    # ------------------------------------- (5c) method-level mechanism tests
    print("\n--- method-level hypotheses ---")
    mrows = []

    # H1 ionic strength: eQ at 0.25 M, dGPredictor at 0.10 M. The Debye-Huckel
    # term scales with the change in sum(charge^2) across the reaction, so a
    # real ionic-strength mismatch shows up as (eQ - dGP) correlating with
    # d_charge_sq.
    m = key["d_charge_sq"].notna()
    rho = stats.spearmanr(key.loc[m, "d_charge_sq"], key.loc[m, "diff"])
    mrows.append({"hypothesis": "H1 ionic-strength mismatch (0.25 vs 0.10 M)",
                  "test": "rho(diff, delta sum charge^2)", "n": int(m.sum()),
                  "statistic": rho.statistic, "p": rho.pvalue})

    # H2 proton handling: dGPredictor drops H+ (KEGG C00080) from the group
    # change vector entirely, then applies a per-compound Legendre transform.
    # If that leaves a residue, reactions with a non-zero net proton count
    # should be offset relative to those without.
    a = key.loc[key["net_proton"] != 0, "diff"]
    b = key.loc[key["net_proton"] == 0, "diff"]
    u = stats.mannwhitneyu(a, b, alternative="two-sided")
    mrows.append({"hypothesis": "H2 net-proton handling",
                  "test": f"median diff {a.median():+.2f} (H+ != 0, n={len(a)}) vs "
                          f"{b.median():+.2f} (H+ == 0, n={len(b)})",
                  "n": len(key), "statistic": float(u.statistic), "p": float(u.pvalue)})
    rho2 = stats.spearmanr(key["net_proton"], key["diff"])
    mrows.append({"hypothesis": "H2b net-proton dose response",
                  "test": "rho(diff, net H+ coefficient)", "n": len(key),
                  "statistic": rho2.statistic, "p": rho2.pvalue})

    # H3 stereo blindness: dGPredictor's active decomposition path is the
    # no-stereo variant, so a pure stereochemical interconversion has an
    # all-zero group-change vector and must predict dG = 0.
    iso = key[key["ec_class"].fillna("none").str.contains("5")]
    if len(iso) >= 10:
        near_zero = (iso["dg_dgp"].abs() < 0.5).mean()
        base = (key["dg_dgp"].abs() < 0.5).mean()
        mrows.append({"hypothesis": "H3 stereo/isomer blindness",
                      "test": f"EC 5 isomerases with |dGP| < 0.5: {near_zero:.1%} "
                              f"vs {base:.1%} baseline",
                      "n": len(iso), "statistic": float(near_zero - base), "p": np.nan})

    # H4 reactant-contribution anchoring: eQuilibrator uses measured reactants
    # directly where it has them, dGPredictor interpolates everything from
    # radius-1 fragments. Reactions built only from very common metabolites
    # should therefore agree better.
    common = {c for c, n in counts.items() if n >= 50}
    key["_all_common"] = [all(c in common for c in stoich[r]) for r in key["rxn"]]
    a2 = key.loc[key["_all_common"], "absdiff"]
    b2 = key.loc[~key["_all_common"], "absdiff"]
    if len(a2) >= 10 and len(b2) >= 10:
        u2 = stats.mannwhitneyu(a2, b2, alternative="two-sided")
        mrows.append({"hypothesis": "H4 common-metabolite anchoring",
                      "test": f"median |diff| {a2.median():.2f} (all reagents common, "
                              f"n={len(a2)}) vs {b2.median():.2f} (n={len(b2)})",
                      "n": len(key), "statistic": float(u2.statistic), "p": float(u2.pvalue)})

    mech = pd.DataFrame(mrows)
    mech.to_csv(OUT_DIR / "mechanism_tests.tsv", sep="\t", index=False, float_format="%.4g")
    print(mech.to_string(index=False))

    print(f"\nwrote outputs to {OUT_DIR}")


if __name__ == "__main__":
    main()
