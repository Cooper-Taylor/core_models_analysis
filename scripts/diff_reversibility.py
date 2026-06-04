#!/usr/bin/env python3
"""Diff reversibility (column 9, header 'reversibility') between dev and claude-changes
branches of /scratch/ctaylor/ModelSEEDDatabase, across the 61 Biochemistry/reaction_NN.tsv
shards. Reads via `git show <branch>:<path>` so no checkout is needed.

Writes outputs to /scratch/ctaylor/core_models_analysis/results/:
  - rev_map_dev.json
  - rev_map_claude.json
  - rev_diff_dev_vs_claude.csv  (rxn_id, dev_rev, claude_rev, dev_direction,
                                 claude_direction, dev_deltag, claude_deltag)

The schema for reaction TSVs has no separate 'direction' column; reversibility is the
single direction-indicator column. We still emit dev_direction / claude_direction columns
(equal to the reversibility) for downstream compatibility with the task spec.
"""

import csv
import io
import json
import subprocess
import sys
from collections import Counter
from pathlib import Path
import os

REPO = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase"))
OUT = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR", "/scratch/ctaylor/core_models_analysis") + "/results")
SHARD_FMT = "Biochemistry/reaction_{:02d}.tsv"
N_SHARDS = 61
BRANCHES = ("dev", "claude-changes")


def read_shard(branch: str, shard_path: str) -> str:
    """Return the raw contents of <branch>:<shard_path> via git show."""
    result = subprocess.run(
        ["git", "show", f"{branch}:{shard_path}"],
        cwd=REPO,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout


def parse_branch(branch: str):
    """Return dict rxn_id -> {'reversibility', 'direction', 'deltag'}."""
    out = {}
    for i in range(N_SHARDS):
        path = SHARD_FMT.format(i)
        text = read_shard(branch, path)
        reader = csv.DictReader(io.StringIO(text), delimiter="\t")
        for row in reader:
            rxn_id = row["id"]
            out[rxn_id] = {
                "reversibility": row.get("reversibility", ""),
                # No separate direction column exists in the schema; mirror reversibility
                # to satisfy the requested CSV columns.
                "direction": row.get("reversibility", ""),
                "deltag": row.get("deltag", ""),
            }
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    print("Parsing dev ...", file=sys.stderr)
    dev = parse_branch("dev")
    print(f"  {len(dev)} reactions", file=sys.stderr)

    print("Parsing claude-changes ...", file=sys.stderr)
    claude = parse_branch("claude-changes")
    print(f"  {len(claude)} reactions", file=sys.stderr)

    dev_rev = {k: v["reversibility"] for k, v in dev.items()}
    claude_rev = {k: v["reversibility"] for k, v in claude.items()}

    (OUT / "rev_map_dev.json").write_text(json.dumps(dev_rev, indent=0, sort_keys=True))
    (OUT / "rev_map_claude.json").write_text(json.dumps(claude_rev, indent=0, sort_keys=True))

    dev_ids = set(dev.keys())
    claude_ids = set(claude.keys())
    only_dev = sorted(dev_ids - claude_ids)
    only_claude = sorted(claude_ids - dev_ids)
    common = sorted(dev_ids & claude_ids)

    changed = []
    transitions = Counter()
    for rid in common:
        d = dev[rid]["reversibility"]
        c = claude[rid]["reversibility"]
        if d != c:
            changed.append(rid)
            transitions[f"{d}->{c}"] += 1

    diff_csv_path = OUT / "rev_diff_dev_vs_claude.csv"
    with diff_csv_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "rxn_id", "dev_rev", "claude_rev",
            "dev_direction", "claude_direction",
            "dev_deltag", "claude_deltag",
        ])
        for rid in changed:
            writer.writerow([
                rid,
                dev[rid]["reversibility"],
                claude[rid]["reversibility"],
                dev[rid]["direction"],
                claude[rid]["direction"],
                dev[rid]["deltag"],
                claude[rid]["deltag"],
            ])

    summary = {
        "n_total_reactions_dev": len(dev),
        "n_total_reactions_claude": len(claude),
        "n_common": len(common),
        "n_changed": len(changed),
        "by_transition": dict(transitions),
        "only_in_dev_count": len(only_dev),
        "only_in_claude_count": len(only_claude),
        "only_in_dev_sample": only_dev[:20],
        "only_in_claude_sample": only_claude[:20],
        "sample_changed": [
            {
                "rxn_id": rid,
                "dev_rev": dev[rid]["reversibility"],
                "claude_rev": claude[rid]["reversibility"],
                "dev_deltag": dev[rid]["deltag"],
                "claude_deltag": claude[rid]["deltag"],
            }
            for rid in changed[:10]
        ],
    }
    (OUT / "rev_diff_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
