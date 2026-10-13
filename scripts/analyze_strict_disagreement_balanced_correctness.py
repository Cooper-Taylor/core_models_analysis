#!/usr/bin/env python3
"""Re-derive the strict-disagreement "which side is correct" percentages
(reports/strictDisagreementEnzymes/STRICT_DISAGREEMENT_ENZYMES.md, section 1)
restricted to the same OK/CI-only population used by
build_llm_vs_recommended_mismatches_balanced.py -- i.e. excluding reactions
whose ModelSEED `status` carries a mass imbalance (`MI` alone or combined with
`CI`). No new research: every verdict (`direction_conflict.enzymology_supports`)
was already collected by the literature subagents behind
results/strict_disagreement_enzymes/strict_disagreement_enzymes.json; this
script only re-tabulates it over a smaller population.

Of the 1,370 strict disagreements, the per-reaction `modelseed_status` field
only ever takes four combinations here (checked against the full 8-combination
database-wide set in reports/reactionBalanceFlags/ -- CPDFORMERROR/CK/EMPTY
simply don't occur in this subset):

    OK        660   kept
    CI        219   kept   (charge-imbalanced, mass already balanced)
    MI        150   excluded (mass-imbalanced)
    CI+MI     341   excluded (mass-and-charge-imbalanced)

Outputs (results/strict_disagreement_enzymes/):
    correctness_balanced.json   verdict_by_group (+ pooled) for the kept 879,
                                 alongside the original 1,370 figures for
                                 side-by-side comparison
    correctness_balanced.tsv    the same numbers, flat, for the figure script
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

ANALYSIS_DIR = Path("/scratch/ctaylor/core_models_analysis")
DELIVERABLE = ANALYSIS_DIR / "results" / "strict_disagreement_enzymes" / "strict_disagreement_enzymes.json"
SUMMARY_IN = ANALYSIS_DIR / "results" / "strict_disagreement_enzymes" / "summary.json"
OUT_DIR = ANALYSIS_DIR / "results" / "strict_disagreement_enzymes"

GROUPS = ("forward_vs_reverse", "reverse_vs_forward")
VERDICTS = ("llm", "modelseed", "neither", "unclear")
KEPT_STATUS = ({"OK"}, {"CI"})


def status_tokens(status: str | None) -> set[str]:
    status = status or ""
    return {f[:f.find(':')] if ':' in f else f for f in status.split('|')}


def tabulate(reactions: list[dict]) -> dict[str, Counter]:
    by_group: dict[str, Counter] = {g: Counter() for g in GROUPS}
    for r in reactions:
        dc = r["direction_conflict"]
        v = dc.get("enzymology_supports")
        if not v:
            continue
        by_group[dc["mismatch_group"]][v] += 1
    return by_group


def main() -> None:
    deliverable = json.load(open(DELIVERABLE))
    reactions = deliverable["reactions"]
    original_summary = json.load(open(SUMMARY_IN))

    status_combo_counts: Counter = Counter()
    kept, excluded = [], []
    for r in reactions:
        toks = status_tokens(r.get("modelseed_status"))
        status_combo_counts['+'.join(sorted(toks))] += 1
        (kept if toks in KEPT_STATUS else excluded).append(r)

    by_group_kept = tabulate(kept)
    by_group_all = {g: Counter(c) for g, c in original_summary["verdict_by_group"].items()}

    def group_block(counter: Counter) -> dict:
        n = sum(counter.values())
        return {
            "n": n,
            "counts": {v: counter.get(v, 0) for v in VERDICTS},
            "pct": {v: round(100.0 * counter.get(v, 0) / n, 1) if n else 0.0 for v in VERDICTS},
        }

    pooled_kept = Counter()
    for g in GROUPS:
        pooled_kept.update(by_group_kept[g])
    pooled_all = Counter()
    for g in GROUPS:
        pooled_all.update(by_group_all[g])

    result = {
        "description": (
            "Strict-disagreement enzymology verdicts (which side the literature backs), "
            "re-tabulated after excluding mass-imbalanced reactions (status MI or CI+MI), "
            "matching the OK/CI-only filter used for the confusion-matrix rebuild."
        ),
        "status_filter": "modelseed_status tokens == {'OK'} or == {'CI'} only",
        "status_combinations_present": dict(status_combo_counts),
        "n_total": len(reactions),
        "n_kept": len(kept),
        "n_excluded": len(excluded),
        "by_group": {
            "OK_CI_only": {g: group_block(by_group_kept[g]) for g in GROUPS},
            "all_1370": {g: group_block(by_group_all[g]) for g in GROUPS},
        },
        "pooled": {
            "OK_CI_only": group_block(pooled_kept),
            "all_1370": group_block(pooled_all),
        },
    }

    with open(OUT_DIR / "correctness_balanced.json", "w") as fh:
        json.dump(result, fh, indent=2)

    with open(OUT_DIR / "correctness_balanced.tsv", "w") as fh:
        fh.write("population\tgroup\tn\tllm\tmodelseed\tneither\tunclear\t"
                  "llm_pct\tmodelseed_pct\tneither_pct\tunclear_pct\n")
        for pop_key, pop_label in (("OK_CI_only", "ok_ci_only"), ("all_1370", "all_1370")):
            for g in GROUPS:
                b = result["by_group"][pop_key][g]
                fh.write(f"{pop_label}\t{g}\t{b['n']}\t" +
                         "\t".join(str(b['counts'][v]) for v in VERDICTS) + "\t" +
                         "\t".join(str(b['pct'][v]) for v in VERDICTS) + "\n")
            b = result["pooled"][pop_key]
            fh.write(f"{pop_label}\tpooled\t{b['n']}\t" +
                     "\t".join(str(b['counts'][v]) for v in VERDICTS) + "\t" +
                     "\t".join(str(b['pct'][v]) for v in VERDICTS) + "\n")

    print(f"status combinations in the 1,370: {dict(status_combo_counts)}")
    print(f"kept (OK/CI-only): {len(kept)}   excluded (MI / CI+MI): {len(excluded)}")
    print()
    for pop_key, pop_label in (("all_1370", "ALL 1,370 (original)"), ("OK_CI_only", "OK/CI-ONLY (879)")):
        print(f"=== {pop_label} ===")
        for g in GROUPS:
            b = result["by_group"][pop_key][g]
            print(f"  {g:<20} n={b['n']:>5}  " +
                  "  ".join(f"{v}={b['counts'][v]} ({b['pct'][v]}%)" for v in VERDICTS))
        b = result["pooled"][pop_key]
        print(f"  {'pooled':<20} n={b['n']:>5}  " +
              "  ".join(f"{v}={b['counts'][v]} ({b['pct'][v]}%)" for v in VERDICTS))
        print()

    print(f"wrote {OUT_DIR / 'correctness_balanced.json'}")
    print(f"wrote {OUT_DIR / 'correctness_balanced.tsv'}")


if __name__ == "__main__":
    main()
