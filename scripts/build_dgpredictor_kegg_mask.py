#!/usr/bin/env python3
"""Build the canonical mask of ModelSEED reactions whose stored dGPredictor
DeltaG'deg should NOT be used, because the KEGG reaction it was predicted from
is not that reaction.

Background. dGPredictor predicts from a KEGG reaction.
``Biochemistry/Thermodynamics/dGPredictor/json_files/`` stages predictions as
``ModelSEED-rxn -> KEGG-R-id -> {dG_mean, dG_uncer}`` and
``Update_Reaction_dGPredictor_Energies.py`` stores the mean over those ids. The
staged files assign a KEGG id to ~50% of the database, but ModelSEED itself
lists a KEGG alias for only ~25%. Where ModelSEED does list one the staged id
always matches; the surplus is unbacked, heavily reused (one KEGG id can be
claimed by 861 ModelSEED reactions) and, checked against KEGG's own reaction
definitions, usually chemically unrelated.

Criterion (``--strict``, the default): drop the dGPredictor value for any
reaction where the staged KEGG id is not one ModelSEED lists as an alias of
that reaction. This is a single crisp rule reproducible from ModelSEEDDatabase
alone.

``--lenient`` additionally *rescues* unvouched reactions that share at least one
non-ubiquitous participant name with their staged KEGG reaction's definition --
these are mostly legacy 1:1 correspondences from the original KEGG-seeded
ModelSEED build (rxn00019 -> R00024 and similar) that were later dropped from
the alias list. Rescuing on a name-string match is weaker evidence than the
alias, so it is not the default; the flag exists so the sensitivity of any
downstream result to that choice can be measured.

This script only READS ModelSEEDDatabase. The mask is written into
core_models_analysis; nothing in the database is modified, and the stored
dGPredictor records stay where they are. Consumers apply the mask at read time.

Output: results/thermo_agreement/dgpredictor_kegg_mask.tsv
        results/thermo_agreement/dgpredictor_kegg_mask.json  (list of dropped ids)
"""
from __future__ import annotations

import argparse
import glob
import json
import os
from collections import Counter, defaultdict
from pathlib import Path

MSDB_ROOT = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase"))
ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
BIOCHEM = MSDB_ROOT / "Biochemistry"
DGP_JSON = BIOCHEM / "Thermodynamics" / "dGPredictor" / "json_files"
KEGG_DEFS = MSDB_ROOT / "Scripts" / "Release" / "archived" / "kegg_reactions.txt"
OUT_DIR = ANALYSIS_DIR / "results" / "thermo_agreement"

MASK_TSV = OUT_DIR / "dgpredictor_kegg_mask.tsv"
MASK_JSON = OUT_DIR / "dgpredictor_kegg_mask.json"

# Participants too common to carry any evidence that two reactions are the same.
UBIQUITOUS = {"h+", "h2o", "co2", "phosphate", "ppi", "o2", "nad", "nadh",
              "nadp", "nadph", "atp", "adp", "amp", "coa", "nh3"}


def load_mask(path: Path | None = None) -> set[str]:
    """Reaction ids whose dGPredictor value must be ignored. Consumer entry point.

    Returns an empty set (with a warning) if the mask has not been built, so a
    consumer script still runs -- but unfiltered, which is the state the mask
    exists to fix, so the warning is loud.
    """
    path = path or MASK_JSON
    if not path.exists():
        print(f"  WARNING: dGPredictor KEGG mask not found at {path}; "
              f"running UNFILTERED. Run build_dgpredictor_kegg_mask.py first.")
        return set()
    obj = json.loads(path.read_text())
    return set(obj["dropped_reactions"])


def kegg_alias_ids(rxn: dict) -> set[str]:
    ids: set[str] = set()
    for alias in rxn.get("aliases") or []:
        if alias.startswith("KEGG:"):
            ids |= {t.strip() for t in alias.split(":", 1)[1].split(";") if t.strip()}
    return ids


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--strict", action="store_true", default=True,
                   help="drop every unvouched reaction (default)")
    g.add_argument("--lenient", action="store_true",
                   help="rescue unvouched reactions whose participants appear in "
                        "the staged KEGG reaction's definition")
    args = ap.parse_args()
    lenient = args.lenient

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("loading ModelSEED reactions ...")
    rxns: dict[str, dict] = {}
    for path in sorted(glob.glob(str(BIOCHEM / "reaction_*.json"))):
        for entry in json.load(open(path)):
            rxns[entry["id"]] = entry

    print("loading staged dGPredictor KEGG assignments ...")
    staged: dict[str, list[str]] = {}
    reuse: Counter[str] = Counter()
    for path in sorted(glob.glob(str(DGP_JSON / "reaction_*_dG.json"))):
        for rxn, payload in json.load(open(path)).items():
            if not isinstance(payload, dict):
                continue
            used = sorted(k for k, v in payload.items()
                          if isinstance(v, dict) and "dG_mean" in v)
            if used:
                staged[rxn] = used
                for k in used:
                    reuse[k] += 1

    kegg_def: dict[str, str] = {}
    if KEGG_DEFS.exists():
        for line in open(KEGG_DEFS):
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 1 and ";" in parts[1]:
                kegg_def[parts[0].replace("rn:", "")] = parts[1].split(";", 1)[1].strip().lower()
    print(f"  {len(staged)} staged reactions, {len(kegg_def)} KEGG definitions available")

    rows = []
    dropped: list[str] = []
    for rxn_id, kegg_ids in sorted(staged.items()):
        entry = rxns.get(rxn_id)
        if entry is None or entry.get("status") == "EMPTY":
            continue
        # Only reactions that actually carry a stored dGPredictor record matter.
        thermo = (entry.get("thermodynamics") or {}).get("dGPredictor")
        if not thermo:
            continue

        alias = kegg_alias_ids(entry)
        vouched = bool(alias) and set(kegg_ids) <= alias

        # participant-name overlap with the staged KEGG reaction's own definition
        names = [s["name"].lower() for s in (entry.get("stoichiometry") or [])]
        informative = [n for n in names if n not in UBIQUITOUS]
        overlap = None
        for kid in kegg_ids:
            kd = kegg_def.get(kid)
            if kd is None or not informative:
                continue
            hits = sum(1 for n in informative if n in kd)
            frac = hits / len(informative)
            overlap = frac if overlap is None else max(overlap, frac)

        if vouched:
            keep, reason = True, "kegg_id_is_modelseed_alias"
        elif lenient and overlap is not None and overlap > 0:
            keep, reason = True, "unvouched_but_participants_match_kegg_definition"
        elif overlap is not None and overlap == 0:
            keep, reason = False, "unvouched_and_no_shared_participant"
        elif overlap is None:
            keep, reason = False, "unvouched_and_unverifiable"
        else:
            keep, reason = False, "unvouched"

        if not keep:
            dropped.append(rxn_id)
        rows.append({
            "rxn": rxn_id,
            "name": entry.get("name", ""),
            "kegg_ids_staged": ";".join(kegg_ids),
            "kegg_ids_alias": ";".join(sorted(alias)),
            "kegg_max_reuse": max((reuse[k] for k in kegg_ids), default=0),
            "vouched": int(vouched),
            "participant_overlap": "" if overlap is None else f"{overlap:.3f}",
            "dg_dgpredictor": thermo[0],
            "keep": int(keep),
            "reason": reason,
        })

    import pandas as pd
    df = pd.DataFrame(rows)
    df.to_csv(MASK_TSV, sep="\t", index=False)
    MASK_JSON.write_text(json.dumps({
        "mode": "lenient" if lenient else "strict",
        "n_with_stored_dgpredictor": len(df),
        "n_dropped": len(dropped),
        "dropped_reactions": sorted(dropped),
    }, indent=1))

    print(f"\nmode: {'lenient' if lenient else 'strict'}")
    print(f"reactions with a stored dGPredictor value : {len(df):6d}")
    print(f"  kept                                    : {int(df['keep'].sum()):6d}")
    print(f"  dropped                                 : {len(dropped):6d}")
    print("\nbreakdown by reason:")
    for reason, grp in df.groupby("reason"):
        print(f"  {'KEEP' if grp['keep'].iat[0] else 'DROP'}  {len(grp):6d}  {reason}")
    print(f"\nwrote {MASK_TSV}")
    print(f"wrote {MASK_JSON}")


if __name__ == "__main__":
    main()
