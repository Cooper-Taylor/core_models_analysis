#!/usr/bin/env python3
"""Derive the user-facing deliverables from classifications.json +
reaction_ec_map.json. Re-run after every `enzyme_functionality_batch.py ingest`.

Outputs (results/enzyme_functionality/):
    reaction_classification.json
        {"monofunctional": [rxn ids], "bifunctional": [...],
         "polyfunctional": [...], "incomplete_ec": [...]}
        Reactions whose multiple ECs disagree in category take the
        max-multiplicity label (polyfunctional > bifunctional >
        monofunctional); these are listed separately in mixed_ec_reactions.
    category_counts.tsv        category, n_reactions, n_ecs  (figure source data)
        Includes a "pending" row (reactions with >=1 complete EC where NONE
        of them are classified yet) so the true 20,490-reaction denominator
        is always recoverable from this file alone, even mid-pipeline.
    ec_frequency_coverage.tsv  ec_rank, cumulative_reactions_covered, is_classified
        The theoretical coverage curve over all 3,883 ECs in frequency order
        (independent of how many are actually classified yet) -- shows why
        classifying the most-cited ECs first front-loads reaction coverage.
        is_classified marks how far actual progress has reached.
"""
from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from ec_functionality_lib import CATEGORIES, CATEGORY_RANK, ensure_out_dir


def main() -> None:
    out_dir = ensure_out_dir()
    reaction_ec_map = json.load(open(out_dir / "reaction_ec_map.json"))
    ec_freq = json.load(open(out_dir / "ec_frequency.json"))
    classifications = json.load(open(out_dir / "classifications.json"))

    buckets: dict[str, list[str]] = {c: [] for c in CATEGORIES}
    mixed_ec_reactions = []
    ec_n_classified_ecs = Counter()

    for rxn in reaction_ec_map["reactions"]:
        rid = rxn["id"]
        if rxn["status"] == "incomplete":
            buckets["incomplete_ec"].append(rid)
            continue

        cats = set()
        unclassified = []
        for ec in rxn["complete_ecs"]:
            rec = classifications.get(ec)
            if rec is None:
                unclassified.append(ec)
            else:
                cats.add(rec["category"])

        if unclassified and not cats:
            continue  # not yet classified at all -- leave out of every bucket for now
        if unclassified and cats:
            # partially classified reaction: still bucket by what's known so far,
            # but flag it since the label may change once the rest are classified.
            mixed_ec_reactions.append({"id": rid, "reason": "partially_classified",
                                        "classified_cats": sorted(cats),
                                        "unclassified_ecs": unclassified})

        if len(cats) > 1:
            mixed_ec_reactions.append({"id": rid, "reason": "disagreeing_ecs",
                                        "classified_cats": sorted(cats)})

        label = max(cats, key=lambda c: CATEGORY_RANK[c])
        buckets[label].append(rid)

    for cat in CATEGORIES:
        buckets[cat].sort()

    with open(out_dir / "reaction_classification.json", "w") as fh:
        json.dump({
            **{c: buckets[c] for c in CATEGORIES},
            "mixed_ec_reactions": mixed_ec_reactions,
        }, fh, indent=2)

    for ec, rec in classifications.items():
        ec_n_classified_ecs[rec["category"]] += 1

    n_total_reactions = reaction_ec_map["n_reactions"]
    n_pending = n_total_reactions - sum(len(buckets[c]) for c in CATEGORIES)
    with open(out_dir / "category_counts.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["category", "n_reactions", "n_ecs"])
        for cat in CATEGORIES:
            n_ecs = ec_n_classified_ecs.get(cat, 0) if cat != "incomplete_ec" else ""
            w.writerow([cat, len(buckets[cat]), n_ecs])
        w.writerow(["pending", n_pending, ""])

    ecs_by_freq = ec_freq["ecs"]  # already sorted descending by n_reactions
    with open(out_dir / "ec_frequency_coverage.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["ec_rank", "cumulative_unique_reactions_covered", "is_classified"])
        covered: set[str] = set()
        for rank, r in enumerate(ecs_by_freq, start=1):
            covered.update(r["reaction_ids"])
            w.writerow([rank, len(covered), int(r["ec"] in classifications)])

    n_classified_reactions = sum(len(buckets[c]) for c in CATEGORIES)
    print(f"ECs classified: {len(classifications)}/{ec_freq['n_unique_ecs']}")
    print(f"reactions assigned a final category: {n_classified_reactions}/{n_total_reactions}")
    for cat in CATEGORIES:
        print(f"  {cat:<16} {len(buckets[cat]):>6} reactions")
    print(f"  {'pending':<16} {n_pending:>6} reactions (>=1 complete EC, none classified yet)")
    if mixed_ec_reactions:
        print(f"flagged (disagreeing or partially-classified EC sets): {len(mixed_ec_reactions)}")
    print()
    print(f"wrote {out_dir / 'reaction_classification.json'}")
    print(f"wrote {out_dir / 'category_counts.tsv'}")
    print(f"wrote {out_dir / 'ec_frequency_coverage.tsv'}")


if __name__ == "__main__":
    main()
