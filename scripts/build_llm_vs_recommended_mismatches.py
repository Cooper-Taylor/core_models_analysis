#!/usr/bin/env python3
"""Cross-tabulate the recommended direction against the LLM ensemble's call.

ModelSEED reactions now carry up to four thermodynamic sources in
`thermodynamics`: "Group contribution", "dGPredictor", "eQuilibrator", and
"LLMs" (a council of LLMs voting a direction with no energy value -- see
ModelSEEDDatabase/Scripts/Thermodynamics/Add_LLM_Direction_Calls.py). The
top-level `reversibility` field is a SEPARATE recommendation computed by
Apply_Evidence_Grades_And_Recommendation.py with fixed precedence
eQuilibrator > dGPredictor > Group contribution -- "LLMs" never enters that
precedence, so it is possible, and apparently common, for the recommendation
and the LLM ensemble to disagree.

This script finds every reaction where both calls exist and are directional
(">" / "<" / "=", excluding "?" abstentions and reactions with no LLM entry at
all), and splits the disagreements into the six possible mismatch types:

    reversible_vs_forward    recommended "="  LLM ">"
    reversible_vs_reverse    recommended "="  LLM "<"
    forward_vs_reversible    recommended ">"  LLM "="
    forward_vs_reverse       recommended ">"  LLM "<"
    reverse_vs_reversible    recommended "<"  LLM "="
    reverse_vs_forward       recommended "<"  LLM ">"

Obsolete reactions (`is_obsolete`) are excluded throughout -- they move the
mismatch total by 1 reaction out of 8,467 and are not part of the live
database.

Data: ModelSEED dev @ 078a395f (/scratch/ctaylor/tmp/devsnap_078a395f),
fetched 2026-10-07.

Outputs:
    results/llm_vs_recommended/mismatches_by_group.json
        {group_key: [reaction record, ...]}  -- the 8,466 disagreements, full
        detail, keyed exactly by the six names above for later lookup.
    results/llm_vs_recommended/summary_counts.json
        the 3x3 recommended x LLM contingency table plus coverage counts
        (how many reactions were excluded for "?" / missing, and why), read
        by scripts/plot_llm_vs_recommended_confusion.py.
"""
from __future__ import annotations

import glob
import json
import os
from collections import Counter
from pathlib import Path

MSDB_ROOT = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/tmp/devsnap_078a395f"))
MSDB_SHA = "078a395f"
ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
OUT_DIR = ANALYSIS_DIR / "results" / "llm_vs_recommended"
BIOCHEM = MSDB_ROOT / "Biochemistry"

DIRECTIONS = (">", "=", "<")
LABEL = {">": "forward", "=": "reversible", "<": "reverse"}

# recommended_op, llm_op -> group key, in the order the user asked for them
GROUPS = {
    ("=", ">"): "reversible_vs_forward",
    ("=", "<"): "reversible_vs_reverse",
    (">", "="): "forward_vs_reversible",
    (">", "<"): "forward_vs_reverse",
    ("<", "="): "reverse_vs_reversible",
    ("<", ">"): "reverse_vs_forward",
}


def reaction_record(rxn: dict, llm_dir: str) -> dict:
    th = rxn.get("thermodynamics") or {}
    return {
        "id": rxn["id"],
        "name": rxn.get("name"),
        "abbreviation": rxn.get("abbreviation"),
        "equation": rxn.get("equation"),
        "definition": rxn.get("definition"),
        "ec_numbers": rxn.get("ec_numbers"),
        "is_transport": bool(rxn.get("is_transport")),
        "recommended_reversibility": rxn.get("reversibility"),
        "llm_direction": llm_dir,
        "thermo_evidence": rxn.get("thermo-evidence"),
        "thermodynamics": {
            src: vals for src, vals in th.items()
        },
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    groups: dict[str, list[dict]] = {name: [] for name in GROUPS.values()}
    matrix: Counter[tuple[str, str]] = Counter()
    n_total = n_obsolete = n_rec_unknown = n_llm_missing = n_llm_abstain = 0
    n_comparable = 0

    for path in sorted(glob.glob(str(BIOCHEM / "reaction_*.json"))):
        for rxn in json.load(open(path)):
            n_total += 1
            if rxn.get("is_obsolete"):
                n_obsolete += 1
                continue

            rec = rxn.get("reversibility")
            th = rxn.get("thermodynamics") or {}
            llm = th.get("LLMs")

            if llm is None:
                n_llm_missing += 1
                continue
            llm_dir = llm[2]
            if llm_dir == "?":
                n_llm_abstain += 1
                continue
            if rec not in DIRECTIONS:
                n_rec_unknown += 1
                continue

            n_comparable += 1
            matrix[(rec, llm_dir)] += 1
            if rec != llm_dir:
                key = GROUPS[(rec, llm_dir)]
                groups[key].append(reaction_record(rxn, llm_dir))

    with open(OUT_DIR / "mismatches_by_group.json", "w") as fh:
        json.dump(groups, fh, indent=2)

    summary = {
        "provenance": {"msdb_dev_sha": MSDB_SHA, "snapshot": str(MSDB_ROOT)},
        "n_reactions_total": n_total,
        "n_obsolete_excluded": n_obsolete,
        "n_recommended_unknown_excluded": n_rec_unknown,
        "n_llm_missing_excluded": n_llm_missing,
        "n_llm_abstained_excluded": n_llm_abstain,
        "n_comparable": n_comparable,
        "matrix": {f"{rec}|{llm}": count for (rec, llm), count in matrix.items()},
        "group_counts": {name: len(recs) for name, recs in groups.items()},
        "n_agree": sum(count for (rec, llm), count in matrix.items() if rec == llm),
        "n_mismatch": sum(len(recs) for recs in groups.values()),
    }
    with open(OUT_DIR / "summary_counts.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"reactions scanned: {n_total}  (obsolete excluded: {n_obsolete})")
    print(f"excluded -- recommended '?': {n_rec_unknown}, no LLM entry: {n_llm_missing}, "
          f"LLM abstained '?': {n_llm_abstain}")
    print(f"comparable (both directional): {n_comparable}")
    print()
    header = "recommended \\ LLM".ljust(18) + "".join(f"{LABEL[d] + ' (' + d + ')':>16}" for d in DIRECTIONS)
    print(header)
    for rec in DIRECTIONS:
        row = f"{LABEL[rec] + ' (' + rec + ')':<18}"
        for llm in DIRECTIONS:
            row += f"{matrix[(rec, llm)]:>16,}"
        print(row)
    print()
    for (rec, llm), key in GROUPS.items():
        print(f"  {key:<24} recommended {rec!r} vs LLM {llm!r}: {len(groups[key]):>6,}")
    print(f"\n  total mismatches: {summary['n_mismatch']:,}   total agreements: {summary['n_agree']:,}")
    print(f"\nwrote {OUT_DIR / 'mismatches_by_group.json'}")
    print(f"wrote {OUT_DIR / 'summary_counts.json'}")


if __name__ == "__main__":
    main()
