#!/usr/bin/env python3
"""Where the original dGPredictor's ModelSEED->KEGG reaction mapping came from,
and why 17,784 of them are wrong.

The staged predictions ModelSEED dev consumes live in
``Biochemistry/Thermodynamics/dGPredictor/json_files/reaction_*_dG.json`` and are
keyed ``rxn##### -> KEGG R-id -> {dG_mean, dG_uncer}``.
``Update_Reaction_dGPredictor_Energies.py`` just averages over the KEGG ids it
finds there, so the mapping was fixed upstream, in the freiburgermsu repo.

This script establishes, from the data alone:
  1. how many staged KEGG ids are backed by a KEGG alias on the ModelSEED
     reaction record (and how many conflict -- none do);
  2. that the unbacked ones are NOT explained by any alias file on disk;
  3. that they ARE explained, at 100%, by a carry-forward: each alias-less
     reaction was handed the KEGG id of the nearest PRECEDING reaction in file
     order that had one.

(3) is the signature of the loop in ``dG_prediction_modelseed_dev_branch.ipynb``,
which sets ``kegg_id_str`` inside an inner ``for`` without ever initialising or
resetting it -- so a reaction whose alias list contains no KEGG entry silently
keeps the previous reaction's id, raises nothing, and never reaches the
``except: 'No KEGG id'`` branch. The ``.py`` port in the repo has since added
``kegg_id_str = None`` per iteration, but the staged data predates that.

Finally it reports what the retrain does with the affected reactions.
"""
from __future__ import annotations

import collections
import glob
import json
from pathlib import Path

DEV = Path("/scratch/ctaylor/tmp/devsnap2")
REPO = Path("/scratch/ctaylor/dgpredictor_repo")
THERMO = DEV / "Biochemistry" / "Thermodynamics" / "dGPredictor"
OUT = Path(__file__).resolve().parents[1] / "results"


def kegg_aliases(entry: dict) -> set[str]:
    ids: set[str] = set()
    for a in entry.get("aliases") or []:
        if a.startswith("KEGG:"):
            ids |= {t.strip() for t in a.split(":", 1)[1].split(";") if t.strip()}
    return ids


def load_alias_tsv(path: Path, kegg_only: bool = True) -> dict[str, set[str]]:
    m: dict[str, set[str]] = collections.defaultdict(set)
    with open(path) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            if kegg_only and p[2] != "KEGG":
                continue
            m[p[0]].add(p[1])
    return dict(m)


def main() -> None:
    staged: dict[str, set[str]] = {}
    for p in sorted(glob.glob(str(THERMO / "json_files" / "reaction_*_dG.json"))):
        for rxn, kegg_map in json.load(open(p)).items():
            if isinstance(kegg_map, dict) and kegg_map:
                staged[rxn] = set(kegg_map)

    order: list[str] = []
    alias: dict[str, set[str]] = {}
    thermo: dict[str, dict] = {}
    for p in sorted(glob.glob(str(DEV / "Biochemistry" / "reaction_*.json"))):
        for e in json.load(open(p)):
            order.append(e["id"])
            alias[e["id"]] = kegg_aliases(e)
            thermo[e["id"]] = e.get("thermodynamics") or {}

    vouched = sum(1 for r, s in staged.items() if alias.get(r) and s & alias[r])
    conflict = sum(1 for r, s in staged.items() if alias.get(r) and not s & alias[r])
    unbacked = [r for r in staged if not alias.get(r)]

    # (2) no alias file on disk accounts for the unbacked ids
    files = {
        "dev Aliases/Unique_ModelSEED_Reaction_Aliases.txt":
            load_alias_tsv(DEV / "Biochemistry" / "Aliases"
                           / "Unique_ModelSEED_Reaction_Aliases.txt"),
        "repo data/Unique_ModelSEED_Reaction_Aliases.txt":
            load_alias_tsv(REPO / "data" / "Unique_ModelSEED_Reaction_Aliases.txt"),
    }
    from_files = {name: sum(1 for r in unbacked if staged[r] & m.get(r, set()))
                  for name, m in files.items()}

    # (3) carry-forward reconstruction
    carry: set[str] | None = None
    predicted: dict[str, set[str] | None] = {}
    for r in order:
        if alias[r]:
            carry = alias[r]
        predicted[r] = carry
    exact = sum(1 for r in unbacked if predicted.get(r) == staged[r])

    reuse = collections.Counter(s for r in unbacked for s in staged[r])

    retrained = set(json.load(open(THERMO / "modelseed_retrained_dG.json")))
    covered = [r for r in unbacked if r in retrained]

    res = {
        "staged_reactions": len(staged),
        "staged_id_is_a_modelseed_alias": vouched,
        "staged_id_conflicts_with_alias": conflict,
        "reaction_has_no_kegg_alias_at_all": len(unbacked),
        "explained_by_an_alias_file_on_disk": from_files,
        "explained_by_carry_forward_from_previous_reaction": {
            "n": exact, "fraction": round(exact / len(unbacked), 6)},
        "top_reused_ids": reuse.most_common(8),
        "affected_reactions_carrying_a_stored_dGPredictor_record":
            sum(1 for r in unbacked if thermo.get(r, {}).get("dGPredictor")),
        "retrain_behaviour_on_affected_reactions": {
            "gives_its_own_value": len(covered),
            "declines_no_complete_structure": len(unbacked) - len(covered),
        },
    }
    (OUT / "kegg_mapping_trace.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
