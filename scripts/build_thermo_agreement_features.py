#!/usr/bin/env python3
"""Build a per-reaction feature table for the three ModelSEED thermodynamic
sources (Group Contribution, eQuilibrator, dGPredictor) so that *sets* of
reactions with high vs. low cross-source agreement can be characterised
chemically and structurally.

Output: results/thermo_agreement/reaction_features.tsv  (one row per reaction
that has a stored, non-"?" DeltaG from at least two of the three sources).

Columns fall into four blocks:

  * ``dg_*`` / ``d_*``  -- the three source DeltaG'deg values (kcal/mol) and
    their pairwise signed differences.
  * stoichiometry / element / mass block -- computed from the ModelSEED
    reaction ``stoichiometry`` list joined to the compound records.
  * provenance block -- KEGG mapping multiplicity (dGPredictor is keyed by
    KEGG reaction id and MSDB stores the *mean* over all KEGG ids mapped to
    one ModelSEED reaction), EC number, status, transport flag.
  * RDKit structural block -- functional-group counts summed over the
    substrate and product sides, plus the net change across the reaction.
    Requires the ``eq3`` conda env (rdkit 2026.03.4).

Nothing in ModelSEEDDatabase is modified; this script only reads.
"""
from __future__ import annotations

import glob
import json
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

MSDB_ROOT = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase"))
ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
OUT_DIR = ANALYSIS_DIR / "results" / "thermo_agreement"
BIOCHEM = MSDB_ROOT / "Biochemistry"
DGP_JSON = BIOCHEM / "Thermodynamics" / "dGPredictor" / "json_files"

# The "dgp" slot is parametrised: ModelSEED carries two dGPredictor records.
#   "dGPredictor"            -- the original, predicted from a KEGG reaction and
#                               therefore exposed to the KEGG mis-mapping.
#   "dGPredictor-ModelSEED"  -- Freiburger's retrain on ModelSEED structures
#                               (radius 1 & 2, BayesianRidge, same 4,001
#                               measurements, pH 7 / I 0.25 M), staged keyed by
#                               ModelSEED reaction id, so no KEGG hop exists.
# --dgp-label selects which one fills dg_dgp; everything downstream is unchanged.
SOURCES = {
    "gc": "Group contribution",
    "eq": "eQuilibrator",
    "dgp": "dGPredictor",
}

# Cofactor / small-molecule participants worth flagging individually. These are
# the compounds whose own formation energies dominate a reaction DeltaG in the
# additive (GC / eQuilibrator) formulation, so they are the natural first place
# to look for systematic cross-source offsets.
COFACTORS = {
    "cpd00002": "atp", "cpd00008": "adp", "cpd00018": "amp",
    "cpd00003": "nad", "cpd00004": "nadh",
    "cpd00006": "nadp", "cpd00005": "nadph",
    "cpd00010": "coa", "cpd00015": "fad", "cpd00982": "fadh2",
    "cpd00050": "fmn", "cpd01270": "fmnh2",
    "cpd00016": "pydx5p", "cpd00087": "thf", "cpd00345": "mthf",
    "cpd00017": "sam", "cpd00019": "sah",
    "cpd00007": "o2", "cpd00025": "h2o2", "cpd00011": "co2",
    "cpd00013": "nh3", "cpd00067": "h", "cpd00001": "h2o",
    "cpd00009": "pi", "cpd00012": "ppi",
    "cpd00042": "gsh", "cpd00111": "gssg",
    "cpd11621": "fdx_ox", "cpd11620": "fdx_red",
    "cpd15500": "menaquinone8", "cpd15499": "menaquinol8",
    "cpd15560": "ubiquinone8", "cpd15561": "ubiquinol8",
    "cpd00109": "cytc_ox", "cpd00110": "cytc_red",
    "cpd11493": "acp", "cpd00023": "glu", "cpd00053": "gln",
}

METALS = {"Fe", "Cu", "Mn", "Zn", "Co", "Ni", "Mo", "Mg", "Ca", "Na", "K",
          "W", "V", "Cd", "Hg", "Ag", "Al"}
HALOGENS = {"F", "Cl", "Br", "I"}

# SMARTS for functional groups. Counts are summed over each side of the
# reaction weighted by stoichiometric coefficient, so ``d_<group>`` is the net
# number of that group created (positive) or destroyed (negative).
SMARTS = {
    "carboxylate": "[CX3](=O)[OX1H0-,OX2H1]",
    "phosphate": "[PX4](=O)([OX1,OX2])([OX1,OX2])[OX1,OX2]",
    "phosphoanhydride": "[PX4](=O)[OX2][PX4](=O)",
    "thioester": "[#6][CX3](=O)[SX2][#6]",
    "thiol": "[#6][SX2H1]",
    "disulfide": "[SX2][SX2]",
    "sulfate_ester": "[SX4](=O)(=O)([OX2])[OX1,OX2]",
    "sulfonate": "[SX4](=O)(=O)[OX1H0-,OX2H1]",
    "primary_amine": "[NX3;H2;!$(NC=O)]",
    "amide": "[NX3][CX3](=[OX1])",
    "aldehyde": "[CX3H1](=O)[#6]",
    "ketone": "[#6][CX3](=O)[#6]",
    "ester": "[#6][CX3](=O)[OX2H0][#6]",
    "hydroxyl": "[OX2H][#6]",
    "ether": "[OX2]([#6])[#6]",
    "alkene": "[CX3]=[CX3]",
    "alkyne": "[CX2]#[CX2]",
    "nitro": "[NX3](=O)=O",
    "nitrile": "[NX1]#[CX2]",
    "aromatic_c": "c",
    "aromatic_n": "n",
    "halide": "[F,Cl,Br,I]",
    "imine": "[CX3]=[NX2]",
    "guanidine": "[NX3][CX3](=[NX2])[NX3]",
    "quaternary_n": "[NX4+]",
    "acetal": "[CX4H1]([OX2])[OX2]",
    "epoxide": "[OX2r3]",
}

FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


def parse_formula(formula: str | None) -> dict[str, int] | None:
    """Element -> count. Returns None for missing / non-atomic formulas.

    ModelSEED uses ``R``/``X`` placeholders for generic side chains and
    ``*``/``n`` for polymers; those are reported so they can be flagged rather
    than silently parsed as real elements.
    """
    if not formula or formula in ("null", "noformula"):
        return None
    counts: dict[str, int] = {}
    pos = 0
    for m in FORMULA_TOKEN.finditer(formula):
        if m.start() != pos:
            return None  # unparseable chunk (e.g. '(' in a polymer formula)
        pos = m.end()
        el, num = m.group(1), m.group(2)
        counts[el] = counts.get(el, 0) + (int(num) if num else 1)
    if pos != len(formula):
        return None
    return counts


def load_compounds() -> dict[str, dict]:
    cpds: dict[str, dict] = {}
    for path in sorted(glob.glob(str(BIOCHEM / "compound_*.json"))):
        for entry in json.load(open(path)):
            cpds[entry["id"]] = entry
    return cpds


def load_reactions() -> dict[str, dict]:
    rxns: dict[str, dict] = {}
    for path in sorted(glob.glob(str(BIOCHEM / "reaction_*.json"))):
        for entry in json.load(open(path)):
            rxns[entry["id"]] = entry
    return rxns


def load_dgp_raw() -> dict[str, dict]:
    """rxn -> {n_kegg_with_dg, kegg_spread_kcal, n_kegg_no_dg, kegg_ids}.

    ``kegg_ids`` is the set of KEGG reaction ids that dGPredictor was actually
    run on for this ModelSEED reaction. This matters: dGPredictor predicts from
    the *KEGG* reaction, and Update_Reaction_dGPredictor_Energies.py stores the
    mean over all of them, so the stored value is only about the ModelSEED
    reaction to the extent that the KEGG mapping is right.
    """
    out: dict[str, dict] = {}
    for path in sorted(glob.glob(str(DGP_JSON / "reaction_*_dG.json"))):
        for rxn, payload in json.load(open(path)).items():
            if not isinstance(payload, dict):
                out[rxn] = {"n_kegg_with_dg": 0, "kegg_spread_kcal": None,
                            "n_kegg_no_dg": 0, "kegg_ids": ""}
                continue
            vals, nodg, used = [], 0, []
            for kegg, sub in payload.items():
                if isinstance(sub, dict) and "dG_mean" in sub:
                    vals.append(sub["dG_mean"])
                    used.append(kegg)
                else:
                    nodg += 1
            out[rxn] = {
                "n_kegg_with_dg": len(vals),
                "kegg_spread_kcal": round((max(vals) - min(vals)) / 4.184, 3) if len(vals) > 1 else 0.0,
                "n_kegg_no_dg": nodg,
                "kegg_ids": ";".join(used),
            }
    return out


def kegg_alias_ids(rxn: dict) -> set[str]:
    """KEGG reaction ids ModelSEED itself lists as aliases of this reaction."""
    ids: set[str] = set()
    for alias in rxn.get("aliases") or []:
        if alias.startswith("KEGG:"):
            ids |= {tok.strip() for tok in alias.split(":", 1)[1].split(";") if tok.strip()}
    return ids


def rdkit_mols(cpds: dict[str, dict]):
    """cpd_id -> (Mol|None, {group: count}). Silences RDKit parse chatter."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    patterns = {name: Chem.MolFromSmarts(sma) for name, sma in SMARTS.items()}
    bad = [n for n, p in patterns.items() if p is None]
    if bad:
        raise SystemExit(f"bad SMARTS: {bad}")

    info: dict[str, dict] = {}
    n_ok = n_fail = n_nosmiles = 0
    for cid, entry in cpds.items():
        smi = entry.get("smiles")
        if not smi:
            info[cid] = {"ok": False, "groups": {}, "rings": None,
                         "arom_rings": None, "heavy": None}
            n_nosmiles += 1
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            info[cid] = {"ok": False, "groups": {}, "rings": None,
                         "arom_rings": None, "heavy": None}
            n_fail += 1
            continue
        groups = {name: len(mol.GetSubstructMatches(pat))
                  for name, pat in patterns.items()}
        ri = mol.GetRingInfo()
        info[cid] = {
            "ok": True,
            "groups": groups,
            "rings": ri.NumRings(),
            "arom_rings": sum(1 for ring in ri.AtomRings()
                              if all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring)),
            "heavy": mol.GetNumHeavyAtoms(),
        }
        n_ok += 1
    print(f"  rdkit: {n_ok} parsed, {n_fail} failed to parse, {n_nosmiles} without SMILES")
    return info


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--no-dgp-mask", action="store_true",
                    help="keep dGPredictor values whose staged KEGG id is unbacked by a "
                         "ModelSEED alias, and write to reaction_features_unmasked.tsv. "
                         "Used to reproduce the evidence figure that motivates the mask.")
    ap.add_argument("--dgp-label", default="dGPredictor",
                    choices=["dGPredictor", "dGPredictor-ModelSEED"],
                    help="which stored dGPredictor record fills the dg_dgp column")
    ap.add_argument("--out", default=None, help="output filename under results/thermo_agreement/")
    cli = ap.parse_args()

    SOURCES["dgp"] = cli.dgp_label
    # The KEGG mask only describes the original KEGG-keyed record. The retrain is
    # staged per ModelSEED reaction id, so no mapping step exists to be wrong and
    # the mask must NOT be applied to it.
    if cli.dgp_label != "dGPredictor":
        cli.no_dgp_mask = True

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from build_dgpredictor_kegg_mask import load_mask
    dgp_mask = set() if cli.no_dgp_mask else load_mask()
    print(f"dgp slot = {cli.dgp_label!r}")
    print(f"dGPredictor KEGG mask: withholding {len(dgp_mask)} reaction(s)"
          + ("  [not applicable / disabled]" if cli.no_dgp_mask else ""))
    print("loading ModelSEED biochemistry ...")
    cpds = load_compounds()
    rxns = load_reactions()
    dgp_raw = load_dgp_raw()
    print(f"  {len(cpds)} compounds, {len(rxns)} reactions, {len(dgp_raw)} dGPredictor entries")
    print("decomposing compound structures ...")
    struct = rdkit_mols(cpds)

    # How many ModelSEED reactions share each staged KEGG id? A KEGG id claimed
    # by hundreds of ModelSEED reactions cannot be describing all of them.
    kegg_reuse: Counter[str] = Counter()
    for info in dgp_raw.values():
        for kid in filter(None, info.get("kegg_ids", "").split(";")):
            kegg_reuse[kid] += 1

    group_names = sorted(SMARTS)
    rows = []
    for rxn_id, rxn in sorted(rxns.items()):
        if rxn.get("status") == "EMPTY":
            continue
        thermo = rxn.get("thermodynamics") or {}
        dg = {}
        for key, label in SOURCES.items():
            # Reactions whose stored dGPredictor value was predicted from a KEGG
            # reaction that is not theirs are dropped from the dGPredictor
            # series only; their GC / eQuilibrator values are unaffected by the
            # KEGG mapping and are kept.
            if key == "dgp" and rxn_id in dgp_mask:
                continue
            triple = thermo.get(label)
            if not triple or len(triple) < 3:
                continue
            if triple[2] in (None, "?"):
                continue
            try:
                dg[key] = float(triple[0])
            except (TypeError, ValueError):
                continue
        if len(dg) < 2:
            continue

        stoich = rxn.get("stoichiometry") or []
        if not stoich:
            continue

        # --- stoichiometry / composition -----------------------------------
        elements: Counter[str] = Counter()
        n_sub = n_prod = 0
        sum_abs = 0.0
        max_abs = 0.0
        masses, carbons = [], []
        n_missing_formula = n_missing_smiles = n_generic = 0
        cof_present = set()
        heavy_atoms = []
        rings = arom_rings = 0
        side_groups = {"sub": Counter(), "prod": Counter()}

        for item in stoich:
            cid = item.get("compound")
            coeff = float(item.get("coefficient", 0) or 0)
            if coeff == 0:
                continue
            if coeff < 0:
                n_sub += 1
            else:
                n_prod += 1
            sum_abs += abs(coeff)
            max_abs = max(max_abs, abs(coeff))
            if cid in COFACTORS:
                cof_present.add(COFACTORS[cid])

            centry = cpds.get(cid, {})
            formula = item.get("formula") or centry.get("formula")
            if formula and re.search(r"[RX*]|\)n|\dn\b", str(formula)):
                n_generic += 1
            parsed = parse_formula(formula)
            if parsed is None:
                n_missing_formula += 1
            else:
                for el, cnt in parsed.items():
                    elements[el] += cnt * abs(coeff)
                carbons.append(parsed.get("C", 0))
            mass = centry.get("mass")
            if isinstance(mass, (int, float)):
                masses.append(mass)

            sinfo = struct.get(cid)
            if not sinfo or not sinfo["ok"]:
                n_missing_smiles += 1
            else:
                heavy_atoms.append(sinfo["heavy"])
                rings += sinfo["rings"]
                arom_rings += sinfo["arom_rings"]
                side = "sub" if coeff < 0 else "prod"
                for gname, gcount in sinfo["groups"].items():
                    side_groups[side][gname] += gcount * abs(coeff)

        el_set = set(elements)
        ec = rxn.get("ec_numbers") or []
        ec_classes = sorted({e.split(".")[0] for e in ec if e and e[0].isdigit()})
        aliases = rxn.get("aliases") or []
        n_kegg_alias = sum(len(a.split(":", 1)[1].split(";"))
                           for a in aliases if a.startswith("KEGG:"))
        raw = dgp_raw.get(rxn_id, {})
        staged_kegg = set(filter(None, raw.get("kegg_ids", "").split(";")))
        alias_kegg = kegg_alias_ids(rxn)
        # "vouched": every KEGG id dGPredictor was run on is one ModelSEED
        # itself lists for this reaction. Where ModelSEED lists no KEGG alias at
        # all but a KEGG id was staged anyway, the mapping came from somewhere
        # outside the curated aliases and cannot be checked against them.
        kegg_vouched = int(bool(alias_kegg) and staged_kegg <= alias_kegg)
        kegg_max_reuse = max((kegg_reuse[k] for k in staged_kegg), default=0)

        row = {
            "rxn": rxn_id,
            "name": rxn.get("name", ""),
            "status": rxn.get("status", ""),
            "is_transport": rxn.get("is_transport", 0),
            "is_obsolete": rxn.get("is_obsolete", 0),
            "stored_reversibility": rxn.get("reversibility", ""),
            "equation": rxn.get("equation", ""),
            "definition": rxn.get("definition", ""),
            "ec": ";".join(ec),
            "ec_class": ";".join(ec_classes) if ec_classes else "none",
            "n_ec": len(ec),
            "dg_gc": dg.get("gc"),
            "dg_eq": dg.get("eq"),
            "dg_dgp": dg.get("dgp"),
            "n_sources": len(dg),
            "n_kegg_alias": n_kegg_alias,
            "kegg_ids_staged": raw.get("kegg_ids", ""),
            "kegg_ids_alias": ";".join(sorted(alias_kegg)),
            "kegg_vouched": kegg_vouched,
            "kegg_max_reuse": kegg_max_reuse,
            "n_kegg_with_dg": raw.get("n_kegg_with_dg", 0),
            "n_kegg_no_dg": raw.get("n_kegg_no_dg", 0),
            "kegg_spread_kcal": raw.get("kegg_spread_kcal"),
            "n_sub": n_sub,
            "n_prod": n_prod,
            "n_participants": n_sub + n_prod,
            "sum_abs_coeff": sum_abs,
            "max_abs_coeff": max_abs,
            "n_elements": len(el_set),
            "elements": ";".join(sorted(el_set)),
            "has_S": int("S" in el_set),
            "has_P": int("P" in el_set),
            "has_N": int("N" in el_set),
            "has_halogen": int(bool(el_set & HALOGENS)),
            "has_metal": int(bool(el_set & METALS)),
            "has_Se": int("Se" in el_set),
            "n_generic_formula": n_generic,
            "n_missing_formula": n_missing_formula,
            "n_missing_smiles": n_missing_smiles,
            "all_have_smiles": int(n_missing_smiles == 0),
            "max_carbon": max(carbons) if carbons else 0,
            "max_mass": max(masses) if masses else 0.0,
            "max_heavy_atoms": max(heavy_atoms) if heavy_atoms else 0,
            "total_rings": rings,
            "total_arom_rings": arom_rings,
            "cofactors": ";".join(sorted(cof_present)) if cof_present else "",
            "n_cofactors": len(cof_present),
        }
        for cof in sorted(set(COFACTORS.values())):
            row[f"cof_{cof}"] = int(cof in cof_present)
        for gname in group_names:
            s, p = side_groups["sub"][gname], side_groups["prod"][gname]
            row[f"g_{gname}"] = s + p
            row[f"d_{gname}"] = p - s
        row["n_groups_changed"] = sum(
            1 for gname in group_names
            if side_groups["prod"][gname] != side_groups["sub"][gname])

        if "gc" in dg and "eq" in dg:
            row["diff_gc_eq"] = dg["gc"] - dg["eq"]
        if "gc" in dg and "dgp" in dg:
            row["diff_gc_dgp"] = dg["gc"] - dg["dgp"]
        if "eq" in dg and "dgp" in dg:
            row["diff_eq_dgp"] = dg["eq"] - dg["dgp"]
        if len(dg) == 3:
            vals = [dg["gc"], dg["eq"], dg["dgp"]]
            row["dg_range3"] = max(vals) - min(vals)
            row["dg_mean3"] = sum(vals) / 3
            row["max_abs_pairdiff"] = max(abs(row["diff_gc_eq"]),
                                          abs(row["diff_gc_dgp"]),
                                          abs(row["diff_eq_dgp"]))
        rows.append(row)

    import pandas as pd
    df = pd.DataFrame(rows)
    out = OUT_DIR / (cli.out or ("reaction_features_unmasked.tsv" if cli.no_dgp_mask
                                 else "reaction_features.tsv"))
    df.to_csv(out, sep="\t", index=False)
    print(f"wrote {out}  ({len(df)} reactions x {len(df.columns)} columns)")
    print(f"  with all 3 sources: {(df['n_sources'] == 3).sum()}")
    print(f"  gc+eq only:         {((df['n_sources'] == 2) & df['dg_dgp'].isna()).sum()}")


if __name__ == "__main__":
    main()
