#!/usr/bin/env python3
"""Find *sets* of ModelSEED reactions where all three thermodynamic sources
agree, and sets where they do not, then test which chemical / structural
features separate them.

Reads results/thermo_agreement/reaction_features.tsv (built by
build_thermo_agreement_features.py) and writes:

  results/thermo_agreement/family_stats.tsv         -- one row per candidate
      reaction family: n, the three pairwise Pearson/Spearman r, median |diff|,
      and a DeltaG-matched null control.
  results/thermo_agreement/concordance_enrichment.tsv -- per-feature
      enrichment in the concordant vs. discordant tails.
  results/thermo_agreement/concordant_set.tsv / discordant_set.tsv

Why the matched null matters
----------------------------
Pearson r between two sources over a set of reactions is driven as much by the
*spread* of DeltaG within that set as by how well the sources track each other:
a family whose reactions all sit near 0 kcal/mol will show a low r even if the
sources agree to within 1 kcal/mol, and a family containing a few +/-500
kcal/mol reactions will show a high r even if everything near zero is noise.
So every family r is compared against the r of random reaction sets drawn to
match that family's own Group-Contribution DeltaG distribution (decile-
stratified resampling). ``excess_*`` = family r minus the matched-null mean;
that is the part of the correlation the family's chemistry actually explains.

Median |pairwise difference| (``mad_*``) is reported alongside because it is
immune to both effects and answers "do these sources agree in kcal/mol".
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
OUT_DIR = ANALYSIS_DIR / "results" / "thermo_agreement"

PAIRS = [("gc", "eq"), ("gc", "dgp"), ("eq", "dgp")]
MIN_FAMILY = 40
N_NULL = 400
RNG = np.random.default_rng(20260804)


# ---------------------------------------------------------------- statistics
def pair_stats(sub: pd.DataFrame, a: str, b: str) -> dict:
    x = sub[f"dg_{a}"].to_numpy(float)
    y = sub[f"dg_{b}"].to_numpy(float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 8 or np.std(x) == 0 or np.std(y) == 0:
        return {f"r_{a}_{b}": np.nan, f"rho_{a}_{b}": np.nan,
                f"mad_{a}_{b}": np.nan, f"n_{a}_{b}": len(x)}
    return {
        f"r_{a}_{b}": float(np.corrcoef(x, y)[0, 1]),
        f"rho_{a}_{b}": float(stats.spearmanr(x, y).statistic),
        f"mad_{a}_{b}": float(np.median(np.abs(x - y))),
        f"n_{a}_{b}": int(len(x)),
    }


def matched_null(pool: pd.DataFrame, family: pd.DataFrame, a: str, b: str,
                 n_null: int = N_NULL) -> tuple[float, float, float]:
    """Mean/sd of r over random pools matched to the family's dg_gc deciles.

    Returns (null_mean_r, null_sd_r, empirical_p) where p is the two-sided
    fraction of null draws at least as extreme as the observed r.
    """
    fam = family.dropna(subset=[f"dg_{a}", f"dg_{b}", "dg_gc"])
    if len(fam) < 12:
        return (np.nan, np.nan, np.nan)
    obs = np.corrcoef(fam[f"dg_{a}"], fam[f"dg_{b}"])[0, 1]

    # Decile-stratify the pool on dg_gc using the family's own quantiles, then
    # draw the same per-stratum counts the family has.
    edges = np.unique(np.quantile(fam["dg_gc"], np.linspace(0, 1, 11)))
    if len(edges) < 3:
        return (np.nan, np.nan, np.nan)
    edges[0], edges[-1] = -np.inf, np.inf
    fam_bins = np.digitize(fam["dg_gc"], edges[1:-1])
    pool_ok = pool.dropna(subset=[f"dg_{a}", f"dg_{b}", "dg_gc"])
    pool_bins = np.digitize(pool_ok["dg_gc"], edges[1:-1])
    idx_by_bin = {b_: np.flatnonzero(pool_bins == b_) for b_ in np.unique(fam_bins)}
    counts = pd.Series(fam_bins).value_counts().to_dict()
    if any(len(idx_by_bin.get(b_, [])) < 2 for b_ in counts):
        return (np.nan, np.nan, np.nan)

    xs = pool_ok[f"dg_{a}"].to_numpy(float)
    ys = pool_ok[f"dg_{b}"].to_numpy(float)
    draws = []
    for _ in range(n_null):
        pick = np.concatenate([RNG.choice(idx_by_bin[b_], size=c, replace=True)
                               for b_, c in counts.items()])
        xv, yv = xs[pick], ys[pick]
        if np.std(xv) == 0 or np.std(yv) == 0:
            continue
        draws.append(np.corrcoef(xv, yv)[0, 1])
    if len(draws) < 50:
        return (np.nan, np.nan, np.nan)
    draws = np.asarray(draws)
    p = float((np.abs(draws - draws.mean()) >= abs(obs - draws.mean())).mean())
    return (float(draws.mean()), float(draws.std()), p)


# ------------------------------------------------------------------ families
def build_families(df: pd.DataFrame) -> dict[str, pd.Series]:
    """name -> boolean mask over df. Families are deliberately overlapping."""
    fam: dict[str, pd.Series] = {}

    fam["ALL (3-source)"] = pd.Series(True, index=df.index)

    # --- provenance / bookkeeping
    fam["KEGG map vouched by MSDB alias"] = df["kegg_vouched"] == 1
    fam["KEGG map inferred (no MSDB alias)"] = df["kegg_vouched"] == 0
    fam["KEGG id unique to this reaction"] = df["kegg_max_reuse"] <= 1
    fam["KEGG id shared by >50 reactions"] = df["kegg_max_reuse"] > 50
    fam["KEGG: single id"] = df["n_kegg_with_dg"] == 1
    fam["KEGG: multiple ids (averaged)"] = df["n_kegg_with_dg"] > 1
    fam["structures: all participants have SMILES"] = df["all_have_smiles"] == 1
    fam["structures: >=1 participant lacks SMILES"] = df["all_have_smiles"] == 0
    fam["formula: generic R/X/polymer participant"] = df["n_generic_formula"] > 0
    fam["formula: fully specified"] = df["n_generic_formula"] == 0
    fam["transport reaction"] = df["is_transport"] == 1
    fam["non-transport"] = df["is_transport"] == 0
    fam["status OK"] = df["status"] == "OK"
    fam["status not OK (unbalanced)"] = df["status"] != "OK"

    # --- EC class
    ec_names = {"1": "1 oxidoreductase", "2": "2 transferase", "3": "3 hydrolase",
                "4": "4 lyase", "5": "5 isomerase", "6": "6 ligase",
                "7": "7 translocase"}
    for digit, label in ec_names.items():
        fam[f"EC {label}"] = df["ec_class"].fillna("none").str.split(";").apply(
            lambda parts, d=digit: d in parts)
    fam["EC none"] = df["ec_class"].fillna("none") == "none"

    # --- size / stoichiometry
    fam["participants <= 4"] = df["n_participants"] <= 4
    fam["participants 5-6"] = df["n_participants"].between(5, 6)
    fam["participants >= 7"] = df["n_participants"] >= 7
    fam["sum|coeff| <= 6"] = df["sum_abs_coeff"] <= 6
    fam["sum|coeff| > 12"] = df["sum_abs_coeff"] > 12
    fam["max|coeff| > 2"] = df["max_abs_coeff"] > 2

    # --- molecular size
    fam["max mass < 200 Da"] = df["max_mass"] < 200
    fam["max mass 200-500 Da"] = df["max_mass"].between(200, 500)
    fam["max mass > 500 Da"] = df["max_mass"] > 500
    fam["max carbon <= 6"] = df["max_carbon"] <= 6
    fam["max carbon 7-20"] = df["max_carbon"].between(7, 20)
    fam["max carbon > 20"] = df["max_carbon"] > 20

    # --- elements
    for el in ("S", "P", "N", "halogen", "metal"):
        fam[f"contains {el}"] = df[f"has_{el}"] == 1
        fam[f"no {el}"] = df[f"has_{el}"] == 0
    fam["C/H/O only"] = (df["has_N"] == 0) & (df["has_P"] == 0) & (df["has_S"] == 0) & \
                        (df["has_metal"] == 0) & (df["has_halogen"] == 0)

    # --- rings
    fam["no rings"] = df["total_rings"] == 0
    fam["has aromatic ring"] = df["total_arom_rings"] > 0
    fam["aliphatic ring only"] = (df["total_rings"] > 0) & (df["total_arom_rings"] == 0)

    # --- cofactors
    for col in [c for c in df.columns if c.startswith("cof_")]:
        name = col[4:]
        mask = df[col] == 1
        if mask.sum() >= MIN_FAMILY:
            fam[f"cofactor: {name}"] = mask
    fam["cofactor: none of tracked set"] = df["n_cofactors"] == 0
    fam["redox pair (NAD/NADP/FAD)"] = df[["cof_nad", "cof_nadh", "cof_nadp",
                                            "cof_nadph", "cof_fad", "cof_fadh2"]].max(axis=1) == 1
    fam["phosphoryl transfer (ATP/ADP/AMP)"] = df[["cof_atp", "cof_adp", "cof_amp"]].max(axis=1) == 1

    # --- functional group presence / change
    for col in [c for c in df.columns if c.startswith("g_")]:
        name = col[2:]
        mask = df[col] > 0
        if MIN_FAMILY <= mask.sum() <= len(df) - MIN_FAMILY:
            fam[f"group present: {name}"] = mask
    for col in [c for c in df.columns if c.startswith("d_") and c != "d_gc_eq"]:
        name = col[2:]
        mask = df[col] != 0
        if MIN_FAMILY <= mask.sum() <= len(df) - MIN_FAMILY:
            fam[f"group changes: {name}"] = mask

    fam["1 group type changes"] = df["n_groups_changed"] <= 1
    fam["2-3 group types change"] = df["n_groups_changed"].between(2, 3)
    fam[">=4 group types change"] = df["n_groups_changed"] >= 4

    return fam


# ----------------------------------------------------------------- main flow
def main() -> None:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--trusted-only", action="store_true",
                    help="restrict every analysis to reactions whose staged KEGG id is "
                         "one ModelSEED itself lists as an alias, so that chemistry "
                         "effects are not confounded with reaction-identity mis-mapping")
    ap.add_argument("--suffix", default="",
                    help="suffix appended to output filenames (e.g. _trusted)")
    args = ap.parse_args()
    sfx = args.suffix or ("_trusted" if args.trusted_only else "")

    df = pd.read_csv(OUT_DIR / "reaction_features.tsv", sep="\t", low_memory=False)
    three = df[df["n_sources"] == 3].reset_index(drop=True).copy()
    if args.trusted_only:
        three = three[three["kegg_vouched"] == 1].reset_index(drop=True).copy()
        print("restricted to reactions with a ModelSEED-vouched KEGG mapping")
    print(f"{len(three)} reactions with all three sources")

    # ---- 1. family scan
    families = build_families(three)
    print(f"scanning {len(families)} candidate families ...")
    rows = []
    for name, mask in families.items():
        sub = three[mask.to_numpy()]
        if len(sub) < MIN_FAMILY:
            continue
        row = {"family": name, "n": len(sub)}
        for a, b in PAIRS:
            row.update(pair_stats(sub, a, b))
            nm, nsd, p = matched_null(three, sub, a, b)
            row[f"null_r_{a}_{b}"] = nm
            row[f"excess_{a}_{b}"] = (row[f"r_{a}_{b}"] - nm) if np.isfinite(nm) else np.nan
            row[f"p_{a}_{b}"] = p
        row["median_abs_dg_gc"] = float(np.median(np.abs(sub["dg_gc"])))
        row["iqr_dg_gc"] = float(np.subtract(*np.percentile(sub["dg_gc"], [75, 25])))
        rows.append(row)
    fam_df = pd.DataFrame(rows).sort_values("r_gc_dgp", ascending=False)
    fam_df.to_csv(OUT_DIR / f"family_stats{sfx}.tsv", sep="\t", index=False, float_format="%.4f")
    print(f"wrote family_stats{sfx}.tsv ({len(fam_df)} families)")

    # ---- 2. per-reaction concordance classes
    #
    # "Concordant" = the three sources agree to within 5 kcal/mol pairwise. 5 is
    # about the quoted accuracy floor of group-contribution style estimators, so
    # anything inside it is agreement-within-method-error rather than a real
    # difference of opinion. "Discordant" = >30 kcal/mol pairwise spread, i.e.
    # far enough apart to flip a reversibility call under any threshold in use.
    three["concordant"] = three["max_abs_pairdiff"] <= 5
    three["discordant"] = three["max_abs_pairdiff"] > 30
    print(f"  concordant (<=5 kcal/mol): {three['concordant'].sum()}")
    print(f"  discordant (>30 kcal/mol): {three['discordant'].sum()}")
    print(f"  middle:                    {(~three['concordant'] & ~three['discordant']).sum()}")

    keep = ["rxn", "name", "ec", "status", "equation", "definition",
            "dg_gc", "dg_eq", "dg_dgp", "max_abs_pairdiff", "dg_range3",
            "n_participants", "sum_abs_coeff", "max_mass", "max_carbon",
            "elements", "cofactors", "n_kegg_with_dg", "kegg_spread_kcal",
            "all_have_smiles", "n_generic_formula", "total_rings",
            "total_arom_rings", "n_groups_changed", "is_transport"]
    three[three["concordant"]].sort_values("max_abs_pairdiff")[keep].to_csv(
        OUT_DIR / f"concordant_set{sfx}.tsv", sep="\t", index=False, float_format="%.3f")
    three[three["discordant"]].sort_values("max_abs_pairdiff", ascending=False)[keep].to_csv(
        OUT_DIR / f"discordant_set{sfx}.tsv", sep="\t", index=False, float_format="%.3f")

    # ---- 3. feature enrichment, concordant vs discordant
    #
    # Reported both raw and stratified by |dG| decile: a feature that merely
    # tracks reaction magnitude would otherwise look like a chemistry effect,
    # since large-|dG| reactions are mechanically more likely to disagree by
    # >30 kcal/mol.
    conc = three[three["concordant"]]
    disc = three[three["discordant"]]
    three["_dgbin"] = pd.qcut(np.abs(three["dg_mean3"]), 10, labels=False, duplicates="drop")

    bin_feats = [c for c in three.columns
                 if c.startswith(("cof_", "has_")) or c in
                 ("all_have_smiles", "is_transport", "is_obsolete")]
    bin_feats += ["_multi_kegg", "_generic", "_status_ok", "_arom", "_ring"]
    three["_multi_kegg"] = (three["n_kegg_with_dg"] > 1).astype(int)
    three["_generic"] = (three["n_generic_formula"] > 0).astype(int)
    three["_status_ok"] = (three["status"] == "OK").astype(int)
    three["_arom"] = (three["total_arom_rings"] > 0).astype(int)
    three["_ring"] = (three["total_rings"] > 0).astype(int)
    conc = three[three["concordant"]]
    disc = three[three["discordant"]]

    enr = []
    for feat in bin_feats:
        if feat not in three.columns:
            continue
        a = int(conc[feat].sum()); b = len(conc) - a
        c = int(disc[feat].sum()); d = len(disc) - c
        if min(a + c, b + d) < 20:
            continue
        odds, p = stats.fisher_exact([[a, b], [c, d]])
        # magnitude-stratified: mean within-decile difference in feature rate
        strat = []
        for _, grp in three.groupby("_dgbin"):
            gc_, gd_ = grp[grp["concordant"]], grp[grp["discordant"]]
            if len(gc_) >= 10 and len(gd_) >= 10:
                strat.append(gc_[feat].mean() - gd_[feat].mean())
        enr.append({
            "feature": feat,
            "rate_concordant": a / max(len(conc), 1),
            "rate_discordant": c / max(len(disc), 1),
            "odds_ratio": odds,
            "p_fisher": p,
            "strat_delta": float(np.mean(strat)) if strat else np.nan,
            "n_strata": len(strat),
        })

    num_feats = ["n_participants", "sum_abs_coeff", "max_abs_coeff", "max_mass",
                 "max_carbon", "max_heavy_atoms", "total_rings", "total_arom_rings",
                 "n_groups_changed", "n_elements", "n_cofactors", "n_kegg_with_dg",
                 "n_missing_smiles", "n_generic_formula"]
    for feat in num_feats:
        u = stats.mannwhitneyu(conc[feat].dropna(), disc[feat].dropna(),
                               alternative="two-sided")
        enr.append({
            "feature": feat,
            "rate_concordant": float(conc[feat].median()),
            "rate_discordant": float(disc[feat].median()),
            "odds_ratio": np.nan,
            "p_fisher": float(u.pvalue),
            "strat_delta": np.nan,
            "n_strata": np.nan,
        })
    enr_df = pd.DataFrame(enr).sort_values("p_fisher")
    enr_df.to_csv(OUT_DIR / f"concordance_enrichment{sfx}.tsv", sep="\t",
                  index=False, float_format="%.5g")
    print("wrote concordance_enrichment.tsv")

    # ---- 4. console summary
    pd.set_option("display.width", 200)
    show = ["family", "n", "r_gc_eq", "r_gc_dgp", "r_eq_dgp",
            "excess_gc_dgp", "mad_gc_dgp", "mad_gc_eq", "median_abs_dg_gc"]
    print("\n=== families where GC-dGPredictor correlates BEST ===")
    print(fam_df.nlargest(18, "r_gc_dgp")[show].to_string(index=False))
    print("\n=== families where GC-dGPredictor correlates WORST ===")
    print(fam_df.nsmallest(18, "r_gc_dgp")[show].to_string(index=False))
    print("\n=== biggest positive EXCESS over dG-matched null (gc vs dgp) ===")
    print(fam_df.nlargest(15, "excess_gc_dgp")[show].to_string(index=False))
    print("\n=== biggest negative excess (gc vs dgp) ===")
    print(fam_df.nsmallest(15, "excess_gc_dgp")[show].to_string(index=False))


if __name__ == "__main__":
    main()
