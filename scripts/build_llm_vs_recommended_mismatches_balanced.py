#!/usr/bin/env python3
"""Cross-tabulate the recommended direction against the LLM ensemble's call,
restricted to reactions whose mass/charge balance `status` tokenises to
exactly {"OK"} or exactly {"CI"} (charge-imbalanced, mass already balanced).

This is a filtered variant of build_llm_vs_recommended_mismatches.py -- see
that script for the full background on the `reversibility` (recommendation)
vs `thermodynamics.LLMs` (council vote) distinction. The only change here is
an extra population filter applied before the recommendation/LLM comparison:
a reaction's `status` field (computed by BiochemPy.Reactions.balanceReaction(),
see Libs/Python/BiochemPy/Reactions.py, and tokenised the same way
Scripts/Validation/Validate_Reactions.py does -- confirmed in
reports/reactionBalanceFlags/REACTION_BALANCE_FLAGS.md) must be exactly "OK"
or exactly "CI". Reactions carrying a mass imbalance (`MI`, alone or combined
with `CI`), a formula error (`CPDFORMERROR`), a combinatorial-expansion flag
(`CK`), or an `EMPTY` equation are excluded, so that equation-quality noise
cannot leak into a thermodynamics-direction comparison.

Same frozen snapshot as the unfiltered matrix, for direct comparability:
ModelSEED dev @ 078a395f (/scratch/ctaylor/tmp/devsnap_078a395f), fetched
2026-10-07.

Outputs (results/llm_vs_recommended_balanced/):
    mismatches_by_group.json   same shape as the unfiltered version
    summary_counts.json        same shape, plus n_status_excluded_by_group
                                broken out by the status group that excluded
                                the reaction (read by
                                plot_llm_vs_recommended_confusion_balanced.py)
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
OUT_DIR = ANALYSIS_DIR / "results" / "llm_vs_recommended_balanced"
BIOCHEM = MSDB_ROOT / "Biochemistry"

DIRECTIONS = (">", "=", "<")
LABEL = {">": "forward", "=": "reversible", "<": "reverse"}

# recommended_op, llm_op -> group key, same six as the unfiltered script
GROUPS = {
    ("=", ">"): "reversible_vs_forward",
    ("=", "<"): "reversible_vs_reverse",
    (">", "="): "forward_vs_reversible",
    (">", "<"): "forward_vs_reverse",
    ("<", "="): "reverse_vs_reversible",
    ("<", ">"): "reverse_vs_forward",
}

KEPT_STATUS = ({"OK"}, {"CI"})


def status_tokens(status: str | None) -> set[str]:
    """Exact parsing convention from Scripts/Validation/Validate_Reactions.py."""
    status = status or ""
    return {f[:f.find(':')] if ':' in f else f for f in status.split('|')}


def status_group(toks: set[str]) -> str:
    """Exclusion-funnel label only -- matches analyze_reaction_balance_flags.py's
    grouping of the 8 combinations actually present in the database."""
    if toks == {"CPDFORMERROR"}:
        return "CPDFORMERROR"
    if toks == {"MI"}:
        return "MI_only"
    if toks == {"CI", "MI"}:
        return "CI+MI"
    if "CK" in toks:
        return "CK_combo"
    if toks == {"EMPTY"}:
        return "EMPTY"
    return "other"


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
        "status": rxn.get("status"),
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
    status_excluded: Counter[str] = Counter()

    for path in sorted(glob.glob(str(BIOCHEM / "reaction_*.json"))):
        for rxn in json.load(open(path)):
            n_total += 1
            if rxn.get("is_obsolete"):
                n_obsolete += 1
                continue

            toks = status_tokens(rxn.get("status"))
            if toks not in KEPT_STATUS:
                status_excluded[status_group(toks)] += 1
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
        "provenance": {
            "msdb_dev_sha": MSDB_SHA,
            "snapshot": str(MSDB_ROOT),
            "status_filter": "status tokens == {'OK'} or == {'CI'} only",
        },
        "n_reactions_total": n_total,
        "n_obsolete_excluded": n_obsolete,
        "n_status_excluded": sum(status_excluded.values()),
        "n_status_excluded_by_group": dict(status_excluded),
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
    print(f"excluded by status filter (kept OK + CI-only): {sum(status_excluded.values()):,}")
    for grp, n in status_excluded.most_common():
        print(f"    {grp:<14} {n:>6,}")
    print(f"excluded -- recommended '?': {n_rec_unknown}, no LLM entry: {n_llm_missing}, "
          f"LLM abstained '?': {n_llm_abstain}")
    print(f"comparable (both directional, OK/CI-only): {n_comparable}")
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
