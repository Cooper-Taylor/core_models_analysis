#!/usr/bin/env python3
"""Resumable batch driver for the enzyme-functionality classification pipeline.

Two subcommands:

  plan [--batch-size N] [--chunk-size N]
      Looks at progress.json + ec_frequency.json + classifications.json,
      takes the next up-to-batch-size *unclassified* ECs in frequency order,
      splits them into chunks, writes
      results/enzyme_functionality/batch_<NN>_manifest.json, and prints a
      cost estimate for that EXACT manifest. Does not call any model or MCP
      tool -- this is the checkpoint the user wants to see before any spend.

  ingest <chunk_result.json> [<chunk_result.json> ...] [--batch NN] [--force]
      Merges worker-subagent output files into classifications.json.
      Each input file must be a JSON object {ec: {category, confidence,
      rationale, evidence}, ...}. Validates category in {monofunctional,
      bifunctional, polyfunctional} and EC format before writing. Refuses to
      overwrite an existing classification unless --force is passed. Marks
      the named batch as done in progress.json (if --batch given).

Cost model (Sonnet 5 pricing, $2/MTok in $10/MTok out -- see
reports/enzymeFunctionality/ENZYME_FUNCTIONALITY_CLASSIFICATION.md for the
assumptions). RECALIBRATED after batch 0 (500 ECs, 20 general-purpose
subagents): actual usage averaged ~5,481 total tokens/EC and ~1.38 tool
calls/EC -- the tool-call count matched the original ~1-2-searches-per-EC
assumption almost exactly, but total tokens ran ~2.6x the pre-batch estimate,
most likely because a general-purpose subagent's full tool catalog (not just
the 2-3 pubmed tools actually used) is loaded every turn, and a ~25-enzyme,
30+-turn conversation resends its growing transcript each turn. The constants
below are scaled by that observed ~2.6x factor; re-check them against actual
usage after future batches and adjust again if the gap reopens.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from ec_functionality_lib import CATEGORIES, OUT_DIR, ensure_out_dir, is_complete_ec

PER_ENZYME_INPUT_TOKENS = 3900
PER_ENZYME_OUTPUT_TOKENS = 1235
PER_CHUNK_OVERHEAD_INPUT_TOKENS = 9100
SONNET5_INPUT_PER_MTOK = 2.00
SONNET5_OUTPUT_PER_MTOK = 10.00
VALID_CATEGORIES = {"monofunctional", "bifunctional", "polyfunctional"}


def load(path: Path, default):
    if not path.exists():
        return default
    with open(path) as fh:
        return json.load(fh)


def estimate_cost(n_ecs: int, n_chunks: int) -> dict:
    input_tokens = n_ecs * PER_ENZYME_INPUT_TOKENS + n_chunks * PER_CHUNK_OVERHEAD_INPUT_TOKENS
    output_tokens = n_ecs * PER_ENZYME_OUTPUT_TOKENS
    low = (input_tokens / 1e6 * SONNET5_INPUT_PER_MTOK
           + output_tokens / 1e6 * SONNET5_OUTPUT_PER_MTOK)
    # ambiguous enzymes needing a second search/fetch can roughly double
    # per-enzyme token spend -- quote a realistic range, not a false-precision point.
    high = low * 1.8
    return {
        "n_ecs": n_ecs, "n_chunks": n_chunks,
        "est_input_tokens": input_tokens, "est_output_tokens": output_tokens,
        "est_cost_low_usd": round(low, 2), "est_cost_high_usd": round(high, 2),
    }


def cmd_plan(args: argparse.Namespace) -> None:
    out_dir = ensure_out_dir()
    progress = load(out_dir / "progress.json", None)
    freq = load(out_dir / "ec_frequency.json", None)
    classifications = load(out_dir / "classifications.json", {})
    if progress is None or freq is None:
        sys.exit("progress.json / ec_frequency.json missing -- run "
                 "build_enzyme_functionality_queue.py first")

    ec_order = progress["ec_order"]
    remaining = [ec for ec in ec_order if ec not in classifications]
    if not remaining:
        print("all ECs already classified -- nothing to plan")
        return

    batch_size = args.batch_size
    chunk_size = args.chunk_size
    batch_ecs = remaining[:batch_size]

    freq_by_ec = {r["ec"]: r for r in freq["ecs"]}
    chunks = [batch_ecs[i:i + chunk_size] for i in range(0, len(batch_ecs), chunk_size)]
    manifest_chunks = []
    for chunk in chunks:
        manifest_chunks.append([
            {
                "ec": ec,
                "n_reactions": freq_by_ec[ec]["n_reactions"],
                "example_reaction_name": freq_by_ec[ec]["example_reaction_name"],
            }
            for ec in chunk
        ])

    batch_index = len(progress["batches"])
    manifest_path = out_dir / f"batch_{batch_index:02d}_manifest.json"
    cost = estimate_cost(len(batch_ecs), len(chunks))
    manifest = {
        "batch_index": batch_index,
        "n_ecs": len(batch_ecs),
        "n_chunks": len(chunks),
        "chunks": manifest_chunks,
        "cost_estimate": cost,
    }
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)

    n_total = progress["n_unique_ecs"]
    n_done = n_total - len(remaining)
    print(f"batch {batch_index}: {len(batch_ecs)} ECs across {len(chunks)} chunks "
          f"(chunk size {chunk_size})")
    print(f"progress before this batch: {n_done}/{n_total} ECs classified")
    print()
    print(f"COST ESTIMATE for this batch (Sonnet 5, ${SONNET5_INPUT_PER_MTOK}/${SONNET5_OUTPUT_PER_MTOK} per MTok):")
    print(f"  ~{cost['est_input_tokens']:,} input tokens, ~{cost['est_output_tokens']:,} output tokens")
    print(f"  estimated cost: ${cost['est_cost_low_usd']:.2f} - ${cost['est_cost_high_usd']:.2f}")
    print()
    print(f"wrote {manifest_path}")


def cmd_ingest(args: argparse.Namespace) -> None:
    out_dir = ensure_out_dir()
    classifications = load(out_dir / "classifications.json", {})
    progress = load(out_dir / "progress.json", None)
    if progress is None:
        sys.exit("progress.json missing -- run build_enzyme_functionality_queue.py first")

    n_added = 0
    for path_str in args.files:
        path = Path(path_str)
        data = json.loads(path.read_text())
        for ec, rec in data.items():
            if not is_complete_ec(ec):
                sys.exit(f"{path}: {ec!r} is not a complete EC number")
            category = rec.get("category")
            if category not in VALID_CATEGORIES:
                sys.exit(f"{path}: {ec} has invalid category {category!r} "
                         f"(must be one of {sorted(VALID_CATEGORIES)})")
            if ec in classifications and not args.force:
                sys.exit(f"{ec} already classified -- pass --force to overwrite")
            classifications[ec] = rec
            n_added += 1

    with open(out_dir / "classifications.json", "w") as fh:
        json.dump(classifications, fh, indent=2, sort_keys=True)

    if args.batch is not None:
        progress["batches"].append({"batch_index": args.batch, "n_ecs_ingested": n_added})
        with open(out_dir / "progress.json", "w") as fh:
            json.dump(progress, fh, indent=2)

    print(f"ingested {n_added} classifications from {len(args.files)} file(s)")
    print(f"total classified: {len(classifications)}/{progress['n_unique_ecs']}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p_plan = sub.add_parser("plan", help="plan the next batch and print its cost estimate")
    p_plan.add_argument("--batch-size", type=int, default=500)
    p_plan.add_argument("--chunk-size", type=int, default=25)
    p_plan.set_defaults(func=cmd_plan)

    p_ingest = sub.add_parser("ingest", help="merge worker results into classifications.json")
    p_ingest.add_argument("files", nargs="+")
    p_ingest.add_argument("--batch", type=int, default=None,
                           help="mark this batch index done in progress.json")
    p_ingest.add_argument("--force", action="store_true")
    p_ingest.set_defaults(func=cmd_ingest)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
