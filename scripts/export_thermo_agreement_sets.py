#!/usr/bin/env python3
"""Export the four named reaction sets that came out of the three-source
DeltaG agreement analysis, as TSVs under results/thermo_agreement/sets/.

  set_A_three_way_concordant.tsv
      All three sources within 5 kcal/mol of each other, KEGG mapping vouched.
      The reactions whose thermodynamics can be treated as settled.

  set_B_group_transfer_discordant.tsv
      KEGG mapping vouched (so this is not a mapping artifact) but the sources
      spread by more than 20 kcal/mol, restricted to group-transfer chemistry
      -- phosphoryl, methyl, acyl, amide-N. The genuinely hard chemistry.

  set_C_kegg_mismapped_suspects.tsv
      The actionable data defect: dGPredictor was run on a KEGG reaction id
      that ModelSEED does not list as an alias of this reaction, that same KEGG
      id is claimed by many other ModelSEED reactions, and dGPredictor lands far
      from Group Contribution and eQuilibrator *while those two agree with each
      other*. Ranked worst first.

  set_D_all_three_disagree.tsv
      Group Contribution and eQuilibrator themselves differ by more than 15
      kcal/mol. Not a dGPredictor problem -- reactions where no source agrees
      and the underlying compound energies should be distrusted.

Also writes sets/SUMMARY.md with counts and the top rows of each.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
DATA_DIR = ANALYSIS_DIR / "results" / "thermo_agreement"
OUT_DIR = DATA_DIR / "sets"

COLS = ["rxn", "name", "ec", "status", "definition",
        "dg_gc", "dg_eq", "dg_dgp", "max_abs_pairdiff",
        "cofactors", "n_participants", "max_carbon", "elements",
        "kegg_ids_staged", "kegg_ids_alias", "kegg_vouched", "kegg_max_reuse",
        "is_transport", "total_rings", "n_groups_changed"]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(DATA_DIR / "reaction_features.tsv", sep="\t", low_memory=False)
    t = df[df["n_sources"] == 3].copy()
    t["absdiff_gc_eq"] = (t["dg_gc"] - t["dg_eq"]).abs()
    t["absdiff_gc_dgp"] = (t["dg_gc"] - t["dg_dgp"]).abs()
    t["absdiff_eq_dgp"] = (t["dg_eq"] - t["dg_dgp"]).abs()
    vouched = t["kegg_vouched"] == 1

    redox = (t[["cof_nad", "cof_nadh", "cof_nadp", "cof_nadph"]].max(axis=1) == 1) | \
            t["ec_class"].fillna("none").str.contains("1")
    transfer = (t[["cof_atp", "cof_adp", "cof_amp", "cof_sam", "cof_sah",
                   "cof_glu", "cof_gln"]].max(axis=1) == 1) | \
               (t["d_phosphoanhydride"] != 0) | (t["d_thioester"] != 0)

    sets: dict[str, tuple[pd.DataFrame, str]] = {}

    # A large block of perfect agreement is trivial: isomerases and racemases
    # whose substrate and product share a group decomposition, so every method
    # returns exactly 0.00 by construction. Those are excluded -- agreement on a
    # reaction where nobody had to estimate anything is not evidence about the
    # methods. |dG| > 2 kcal/mol from at least one source is the filter.
    dg_scale = t[["dg_gc", "dg_eq", "dg_dgp"]].abs().max(axis=1)
    t["dg_scale"] = dg_scale
    nontrivial = dg_scale > 2

    a = t[vouched & nontrivial & (t["max_abs_pairdiff"] <= 5)].sort_values(
        "dg_scale", ascending=False)
    a = a.assign(is_redox=redox.reindex(a.index).astype(int))
    sets["set_A_three_way_concordant"] = (
        a, "all three sources within 5 kcal/mol and at least one |ΔG′°| > 2 "
           "kcal/mol (trivial isomerase zeros excluded), KEGG mapping vouched; "
           "sorted largest-ΔG first")

    b = t[vouched & transfer & (t["max_abs_pairdiff"] > 20)].sort_values(
        "max_abs_pairdiff", ascending=False)
    sets["set_B_group_transfer_discordant"] = (
        b, "group-transfer chemistry, vouched mapping, sources spread >20 kcal/mol")

    # Same chemistry, but restricted to ordinary single-step magnitudes so the
    # list is not swamped by prenyl/polymer reactions running to 10^3 kcal/mol.
    b2 = b[b["dg_scale"] <= 100]
    sets["set_B2_group_transfer_discordant_ordinary_scale"] = (
        b2, "as set B but restricted to reactions where every source is inside "
            "±100 kcal/mol -- ordinary metabolic chemistry, not aggregate/polymer "
            "reactions")

    # Set C describes what the mask REMOVED, so it is built from the mask file
    # rather than from the (already filtered) feature table -- in the filtered
    # table these reactions no longer have a dGPredictor value at all.
    mask_tsv = DATA_DIR / "dgpredictor_kegg_mask.tsv"
    if mask_tsv.exists():
        m = pd.read_csv(mask_tsv, sep="\t")
        c = m[m["keep"] == 0].sort_values(
            ["kegg_max_reuse", "rxn"], ascending=[False, True]).copy()
        c["dg_gc"] = c["rxn"].map(dict(zip(t["rxn"], t["dg_gc"])))
        c["dg_eq"] = c["rxn"].map(dict(zip(t["rxn"], t["dg_eq"])))
        c = c.rename(columns={"dg_dgpredictor": "dg_dgp"})
    else:
        c = pd.DataFrame(columns=["rxn", "name", "dg_gc", "dg_eq", "dg_dgp"])
    sets["set_C_kegg_mismapped_withheld"] = (
        c, "the reactions the mask withholds: their stored dGPredictor ΔG′° was "
           "predicted from a KEGG reaction ModelSEED does not list for them. "
           "Sorted by how many ModelSEED reactions share that KEGG id. "
           "`dg_dgp` is the withheld value; `dg_gc`/`dg_eq` are blank where the "
           "reaction has no Group-Contribution / eQuilibrator value either.")

    d = t[t["absdiff_gc_eq"] > 15].sort_values("absdiff_gc_eq", ascending=False)
    sets["set_D_all_three_disagree"] = (
        d, "Group Contribution and eQuilibrator themselves differ by >15 kcal/mol")

    lines = ["# Three-source ΔG agreement: named reaction sets", "",
             f"Built from `results/thermo_agreement/reaction_features.tsv` "
             f"({len(t):,} reactions carrying all three sources' ΔG′°).", ""]
    for name, (frame, desc) in sets.items():
        cols = [c_ for c_ in COLS if c_ in frame.columns]
        extra = [c_ for c_ in ("is_redox", "dg_scale", "absdiff_gc_eq",
                               "absdiff_gc_dgp", "absdiff_eq_dgp")
                 if c_ in frame.columns]
        frame[cols + extra].to_csv(OUT_DIR / f"{name}.tsv", sep="\t",
                                   index=False, float_format="%.3f")
        print(f"{name}.tsv  n={len(frame)}")
        lines += [f"## `{name}.tsv` — {len(frame):,} reactions", "", desc, ""]
        head = frame.head(8)[["rxn", "name", "dg_gc", "dg_eq", "dg_dgp"]]
        lines += ["| reaction | name | GC | eQ | dGP |",
                  "|---|---|---:|---:|---:|"]
        for _, row in head.iterrows():
            nm = str(row["name"])[:58]
            lines.append(f"| {row['rxn']} | {nm} | {row['dg_gc']:.2f} | "
                         f"{row['dg_eq']:.2f} | {row['dg_dgp']:.2f} |")
        lines.append("")

    (OUT_DIR / "SUMMARY.md").write_text("\n".join(lines))
    print(f"wrote {OUT_DIR / 'SUMMARY.md'}")


if __name__ == "__main__":
    main()
