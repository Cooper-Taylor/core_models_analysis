#!/usr/bin/env python3
"""Why do eQuilibrator and dGPredictor value specific metabolites differently?

Follow-on to analyze_eq_vs_dgpredictor.py, which showed that 80% of the
reaction-level disagreement is explained by a fixed per-metabolite offset. This
script asks what those metabolites have in common, and tests the explanation
against what each method actually does.

The two methods fail in different regimes, and the regimes are predictable from
their construction:

  eQuilibrator = component contribution. Two layers. The *reactant
  contribution* layer regresses directly on measured reactions, so any compound
  that appears in the TECRDB training set is anchored to data and carries a
  small uncertainty. Only the component orthogonal to that span falls back to
  group contribution, and those compounds carry a visibly larger uncertainty.
  eQuilibrator publishes that uncertainty per compound, so we can read off which
  layer a compound was served by.

  dGPredictor = one layer. Every compound is decomposed into radius-1
  atom-centred fragments (Chem.FindAtomEnvironmentOfRadiusN, radius=1) and the
  reaction is the net fragment-count change; an OLS fit with no intercept turns
  fragment counts into energies. There is no reactant-contribution shortcut, so
  a heavily-measured metabolite gets no special treatment. Radius 1 also means
  the descriptor sees each atom's immediate neighbours only -- it cannot
  represent conjugation, aromatic stabilisation, ring strain, or stereochemistry
  (the active decomposition path is explicitly the no-stereo variant).

Predictions tested here:
  P1  |offset| rises with eQuilibrator's own uncertainty for that compound
      (i.e. the disagreement concentrates where eQuilibrator itself fell back to
      group contribution).
  P2  |offset| rises with the extent of conjugation / aromatic + polyene
      systems, which radius-1 fragments cannot see.
  P3  |offset| rises with molecular size (more fragments, more accumulated
      error, and further outside the training distribution).
  P4  |offset| rises with ring count / fused-ring systems (strain and
      stabilisation are non-local).

Outputs (results/eq_vs_dgp/):
  metabolite_profile.tsv     per-compound offset + structural + uncertainty features
  metabolite_predictors.tsv  the four predictions, tested

Requires the ``eq3`` env (rdkit).
"""
from __future__ import annotations

import glob
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

MSDB_ROOT = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase"))
ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
BIOCHEM = MSDB_ROOT / "Biochemistry"
OUT_DIR = Path(os.environ.get("EQDGP_OUT", str(ANALYSIS_DIR / "results" / "eq_vs_dgp")))
EQ_CPD_TBL = BIOCHEM / "Thermodynamics" / "eQuilibrator" / "MetaNetX_Compound_Energies.tbl"
MNX_STRUCT = BIOCHEM / "Structures" / "MetaNetX" / "Structures_in_ModelSEED_and_eQuilibrator.txt"


def main() -> None:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")

    off = pd.read_csv(OUT_DIR / "compound_offsets.tsv", sep="\t")
    print(f"{len(off)} compounds with a fitted eQ-minus-dGP offset")

    cpds = {}
    for path in sorted(glob.glob(str(BIOCHEM / "compound_*.json"))):
        for entry in json.load(open(path)):
            cpds[entry["id"]] = entry

    # eQuilibrator's own per-compound uncertainty, via the MetaNetX id it used.
    # A small error means the compound sat in the reactant-contribution span
    # (anchored to measured reactions); a large one means it fell through to the
    # group-contribution layer.
    eq_err: dict[str, float] = {}
    for line in open(EQ_CPD_TBL):
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 3:
            try:
                eq_err[parts[0]] = float(parts[2])
            except ValueError:
                continue
    key2mnx: dict[str, str] = {}
    for line in open(MNX_STRUCT):
        parts = line.rstrip("\n").split("\t")
        if len(parts) == 2:
            key2mnx.setdefault(parts[1], parts[0])

    rows = []
    for _, r in off.iterrows():
        cid = r["compound"]
        entry = cpds.get(cid, {})
        smi = entry.get("smiles")
        rec = {
            "compound": cid, "name": r["name"], "formula": r["formula"],
            "charge": r["charge"], "n_reactions": r["n_reactions"],
            "offset": r["offset_eq_minus_dgp"],
            "abs_offset": abs(r["offset_eq_minus_dgp"]),
            "mass": entry.get("mass"),
        }
        ik = entry.get("inchikey")
        mnx = key2mnx.get(ik) if ik else None
        rec["eq_uncertainty"] = eq_err.get(mnx) if mnx else np.nan

        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is not None:
            rec["heavy_atoms"] = mol.GetNumHeavyAtoms()
            rec["n_rings"] = mol.GetRingInfo().NumRings()
            rec["n_arom_atoms"] = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
            # conjugation: bonds RDKit marks conjugated, i.e. delocalised
            # electron systems that a radius-1 fragment cannot represent
            rec["n_conj_bonds"] = sum(1 for b in mol.GetBonds() if b.GetIsConjugated())
            rec["frac_conj"] = (rec["n_conj_bonds"] / mol.GetNumBonds()
                                if mol.GetNumBonds() else 0.0)
            rec["n_stereo"] = len(Chem.FindMolChiralCenters(
                mol, includeUnassigned=True, useLegacyImplementation=False))
            rec["n_rotatable"] = rdMolDescriptors.CalcNumRotatableBonds(mol)
            rec["logp"] = Descriptors.MolLogP(mol)
            # distinct radius-1 environments = how many descriptor columns this
            # compound occupies in dGPredictor's design matrix
            envs = set()
            for atom in mol.GetAtoms():
                env = Chem.FindAtomEnvironmentOfRadiusN(mol, 1, atom.GetIdx())
                atoms = set()
                for bidx in env:
                    b = mol.GetBondWithIdx(bidx)
                    atoms.add(b.GetBeginAtomIdx()); atoms.add(b.GetEndAtomIdx())
                if not atoms:
                    atoms = {atom.GetIdx()}
                try:
                    envs.add(Chem.MolFragmentToSmiles(mol, atomsToUse=list(atoms),
                                                      bondsToUse=list(env),
                                                      canonical=True))
                except RuntimeError:
                    # RDKit refuses to canonicalise a fragment that cuts through
                    # a stereo double bond. dGPredictor hits the same wall and
                    # drops the compound as non-decomposable; count it so the
                    # frequency is visible rather than silently skipped.
                    envs.add(f"__uncanonicalisable_{atom.GetIdx()}")
                    rec["frag_failures"] = rec.get("frag_failures", 0) + 1
            rec["n_distinct_r1_frags"] = len(envs)
            rec.setdefault("frag_failures", 0)
        rows.append(rec)

    prof = pd.DataFrame(rows)
    prof.to_csv(OUT_DIR / "metabolite_profile.tsv", sep="\t", index=False,
                float_format="%.4f")

    # Weight by evidence: a compound seen in 5 reactions has a noisier fitted
    # offset than one seen in 400. Restrict the correlations to the better-
    # supported half so the tests are not dominated by fit noise.
    solid = prof[prof["n_reactions"] >= 10].dropna(subset=["heavy_atoms"])
    print(f"{len(solid)} compounds with >=10 supporting reactions and a structure")

    tests = []
    for label, col in [
        ("P1 eQuilibrator's own uncertainty", "eq_uncertainty"),
        ("P2 conjugated bonds", "n_conj_bonds"),
        ("P2b fraction of bonds conjugated", "frac_conj"),
        ("P2c aromatic atoms", "n_arom_atoms"),
        ("P3 heavy atoms", "heavy_atoms"),
        ("P3b mass", "mass"),
        ("P3c distinct radius-1 fragments", "n_distinct_r1_frags"),
        ("P4 rings", "n_rings"),
        ("-- stereocentres", "n_stereo"),
        ("-- rotatable bonds", "n_rotatable"),
        ("-- |charge|", "charge"),
    ]:
        sub = solid.dropna(subset=[col])
        if len(sub) < 20:
            continue
        v = sub[col].abs() if col == "charge" else sub[col]
        rho = stats.spearmanr(v, sub["abs_offset"])
        tests.append({"prediction": label, "n": len(sub),
                      "rho_vs_abs_offset": rho.statistic, "p": rho.pvalue})
    tdf = pd.DataFrame(tests).sort_values("rho_vs_abs_offset", ascending=False)
    tdf.to_csv(OUT_DIR / "metabolite_predictors.tsv", sep="\t", index=False,
               float_format="%.4g")
    pd.set_option("display.width", 200)
    print("\n--- what predicts a large eQ-vs-dGPredictor metabolite offset? ---")
    print(tdf.to_string(index=False))

    # eQuilibrator layer split, read off its published uncertainty
    u = solid.dropna(subset=["eq_uncertainty"])
    if len(u) > 40:
        lo = u[u["eq_uncertainty"] <= u["eq_uncertainty"].median()]
        hi = u[u["eq_uncertainty"] > u["eq_uncertainty"].median()]
        mw = stats.mannwhitneyu(lo["abs_offset"], hi["abs_offset"],
                                alternative="two-sided")
        print(f"\neQuilibrator-anchored vs group-estimated compounds "
              f"(split at its median published uncertainty):")
        print(f"  low  uncertainty (RC-anchored, n={len(lo)}): "
              f"median |offset| = {lo['abs_offset'].median():.2f} kcal/mol")
        print(f"  high uncertainty (GC fallback,  n={len(hi)}): "
              f"median |offset| = {hi['abs_offset'].median():.2f} kcal/mol")
        print(f"  Mann-Whitney p = {mw.pvalue:.3g}")

    print("\n--- most-implicated metabolites, with structural context ---")
    show = ["name", "formula", "charge", "n_reactions", "offset", "eq_uncertainty",
            "heavy_atoms", "n_rings", "n_arom_atoms", "n_conj_bonds", "n_distinct_r1_frags"]
    print(prof[prof["n_reactions"] >= 15].nlargest(15, "abs_offset")[show]
          .to_string(index=False))
    print("\n--- best-agreed metabolites (>=50 reactions) ---")
    print(prof[prof["n_reactions"] >= 50].nsmallest(12, "abs_offset")[show]
          .to_string(index=False))


if __name__ == "__main__":
    main()
