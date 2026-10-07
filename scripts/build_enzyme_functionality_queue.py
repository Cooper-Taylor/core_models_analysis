#!/usr/bin/env python3
"""One-time, idempotent extraction for the enzyme-functionality classification
pipeline (see reports/enzymeFunctionality/ENZYME_FUNCTIONALITY_CLASSIFICATION.md).

Walks the same 20,490-reaction "comparable" set as
build_llm_vs_recommended_mismatches.py and splits it by EC completeness:

  - reactions with >=1 complete EC (`\\d+\\.\\d+\\.\\d+\\.\\d+`) need a literature
    search per unique EC to determine monofunctional/bifunctional/polyfunctional
  - reactions with no ec_numbers at all, or where every ec_numbers entry is
    incomplete (`x.x.x.-`), go straight into the 4th category -- no search
    needed.

Outputs (results/enzyme_functionality/):
    reaction_ec_map.json   per-reaction {id, name, ec_numbers, complete_ecs, status}
    ec_frequency.json      unique complete ECs, sorted descending by how many
                            reactions in the comparable set cite them -- this
                            IS the work queue for enzyme_functionality_batch.py
    progress.json          batch checkpoint (created only if absent)
    classifications.json   EC -> classification (created empty only if absent)

Safe to re-run: reaction_ec_map.json and ec_frequency.json are always
regenerated from the snapshot, but progress.json / classifications.json are
only *initialized* if missing -- re-running this script never erases batch
progress already made.
"""
from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from ec_functionality_lib import (
    MSDB_ROOT, MSDB_SHA, ensure_out_dir, is_complete_ec, iter_comparable_reactions,
)


def main() -> None:
    out_dir = ensure_out_dir()

    reaction_records = []
    ec_reactions: dict[str, list[str]] = defaultdict(list)
    ec_example_name: dict[str, str] = {}
    n_total = 0
    n_no_ec = 0
    n_incomplete_only = 0
    n_has_complete = 0

    for rxn in iter_comparable_reactions():
        n_total += 1
        rid = rxn["id"]
        name = rxn.get("name")
        ecs = rxn.get("ec_numbers") or []
        complete = [e for e in ecs if is_complete_ec(e)]
        if not ecs:
            status = "incomplete"
            n_no_ec += 1
        elif not complete:
            status = "incomplete"
            n_incomplete_only += 1
        else:
            status = "complete"
            n_has_complete += 1
            for e in complete:
                ec_reactions[e].append(rid)
                ec_example_name.setdefault(e, name)

        reaction_records.append({
            "id": rid,
            "name": name,
            "ec_numbers": ecs,
            "complete_ecs": complete,
            "status": status,
        })

    reaction_records.sort(key=lambda r: r["id"])
    with open(out_dir / "reaction_ec_map.json", "w") as fh:
        json.dump({
            "provenance": {"msdb_dev_sha": MSDB_SHA, "snapshot": str(MSDB_ROOT)},
            "n_reactions": n_total,
            "reactions": reaction_records,
        }, fh, indent=2)

    ec_freq = sorted(
        (
            {
                "ec": ec,
                "n_reactions": len(rids),
                "reaction_ids": sorted(rids),
                "example_reaction_name": ec_example_name[ec],
            }
            for ec, rids in ec_reactions.items()
        ),
        key=lambda r: (-r["n_reactions"], r["ec"]),
    )
    with open(out_dir / "ec_frequency.json", "w") as fh:
        json.dump({
            "provenance": {"msdb_dev_sha": MSDB_SHA, "snapshot": str(MSDB_ROOT)},
            "n_unique_ecs": len(ec_freq),
            "ecs": ec_freq,
        }, fh, indent=2)

    progress_path = out_dir / "progress.json"
    if not progress_path.exists():
        with open(progress_path, "w") as fh:
            json.dump({
                "batch_size": 500,
                "n_unique_ecs": len(ec_freq),
                "ec_order": [r["ec"] for r in ec_freq],
                "batches": [],
            }, fh, indent=2)
        print(f"initialized {progress_path}")
    else:
        print(f"{progress_path} already exists -- left untouched")

    classifications_path = out_dir / "classifications.json"
    if not classifications_path.exists():
        with open(classifications_path, "w") as fh:
            json.dump({}, fh, indent=2)
        print(f"initialized {classifications_path}")
    else:
        print(f"{classifications_path} already exists -- left untouched")

    print()
    print(f"comparable reactions: {n_total:,}")
    print(f"  no ec_numbers field at all:            {n_no_ec:,}")
    print(f"  ec_numbers present but all incomplete: {n_incomplete_only:,}")
    print(f"  -> incomplete_ec (4th category, no search needed): {n_no_ec + n_incomplete_only:,}")
    print(f"  >=1 complete EC (needs literature search): {n_has_complete:,}")
    print(f"unique complete EC numbers to classify: {len(ec_freq):,}")
    print()
    print("top 10 most-cited ECs:")
    for r in ec_freq[:10]:
        print(f"  {r['ec']:<12} {r['n_reactions']:>4} reactions   e.g. {r['example_reaction_name']}")
    print()
    print(f"wrote {out_dir / 'reaction_ec_map.json'}")
    print(f"wrote {out_dir / 'ec_frequency.json'}")


if __name__ == "__main__":
    main()
