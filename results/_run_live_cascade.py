"""Run the live reversibility cascade and diff vs stored TSV reversibility."""
import csv
import json
import sys

sys.path.insert(0, '/scratch/ctaylor/core_models_analysis/scripts')
sys.path.insert(0, '/scratch/ctaylor/ModelSEEDDatabase/Libs/Python')

import reversibility_lib as lib
from BiochemPy import Reactions

RESULTS = "/scratch/ctaylor/core_models_analysis/results"
LIVE_CSV = f"{RESULTS}/rxn_directions_cascade_live.csv"
LIVE_JSON = f"{RESULTS}/rxn_directions_cascade_live.json"
DIFF_CSV = f"{RESULTS}/rev_diff_stored_vs_cascade_live.csv"
STORED_MAP = f"{RESULTS}/rev_map_dev.json"
SUMMARY_JSON = f"{RESULTS}/_live_cascade_summary.json"


def main() -> None:
    print("Loading reactions ...", flush=True)
    rxns = Reactions().loadReactions()
    print(f"Loaded {len(rxns)} reactions", flush=True)

    print("Running cascade (EQ, gc_first=True, default ReversibilityConfig) ...", flush=True)
    cascade = lib.run_cascade(
        rxns, db_level="EQ", cfg=lib.ReversibilityConfig(), gc_first=True
    )
    print(f"Cascade produced {len(cascade)} entries", flush=True)

    # Capture deltag / deltagerr for the diff sampling (post-cascade; GC pass
    # mutates only 'reversibility' on the rxn dicts, deltag/deltagerr are not
    # touched).
    dg = {rid: rxns[rid].get("deltag") for rid in cascade}
    dge = {rid: rxns[rid].get("deltagerr") for rid in cascade}

    # --- Live outputs ------------------------------------------------------
    live_flag = {rid: rev for rid, (_status, rev) in cascade.items()}
    live_status = {rid: status for rid, (status, _rev) in cascade.items()}

    with open(LIVE_CSV, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["rxn_id", "reversibility", "status"])
        for rid in sorted(cascade):
            w.writerow([rid, live_flag[rid], live_status[rid]])

    with open(LIVE_JSON, "w") as fh:
        json.dump(live_flag, fh, sort_keys=True)
    print(f"Wrote {LIVE_CSV} and {LIVE_JSON}", flush=True)

    # --- Stored map --------------------------------------------------------
    with open(STORED_MAP) as fh:
        stored = json.load(fh)
    print(f"Loaded {len(stored)} stored entries from {STORED_MAP}", flush=True)

    # --- Diff --------------------------------------------------------------
    all_ids = set(live_flag) | set(stored)
    n_match = n_diff = n_only_live = n_only_stored = 0
    transitions: dict[str, int] = {}
    diff_rows: list[dict] = []
    sample_diffs: list[dict] = []

    for rid in sorted(all_ids):
        in_live = rid in live_flag
        in_stored = rid in stored
        lv = live_flag.get(rid)
        sv = stored.get(rid)
        if in_live and in_stored:
            if lv == sv:
                n_match += 1
                continue
            n_diff += 1
            key = f"{sv}->{lv}"
            transitions[key] = transitions.get(key, 0) + 1
            diff_rows.append({
                "rxn_id": rid,
                "stored": sv,
                "live": lv,
                "deltag": dg.get(rid, ""),
                "deltagerr": dge.get(rid, ""),
                "status_live": live_status.get(rid, ""),
            })
        elif in_live:
            n_only_live += 1
            diff_rows.append({
                "rxn_id": rid,
                "stored": "",
                "live": lv,
                "deltag": dg.get(rid, ""),
                "deltagerr": dge.get(rid, ""),
                "status_live": live_status.get(rid, ""),
            })
        else:
            n_only_stored += 1
            diff_rows.append({
                "rxn_id": rid,
                "stored": sv,
                "live": "",
                "deltag": "",
                "deltagerr": "",
                "status_live": "",
            })

    # Pick 5 sample diffs: prefer mismatched-but-present-in-both rows
    paired_diffs = [r for r in diff_rows if r["stored"] and r["live"]]
    sample_pool = paired_diffs if paired_diffs else diff_rows
    sample_diffs = sample_pool[:5]

    with open(DIFF_CSV, "w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["rxn_id", "stored", "live", "deltag", "deltagerr", "status_live"],
        )
        w.writeheader()
        for row in diff_rows:
            w.writerow(row)
    print(f"Wrote diff CSV {DIFF_CSV} ({len(diff_rows)} rows)", flush=True)

    summary = {
        "n_total": len(all_ids),
        "n_live": len(live_flag),
        "n_stored": len(stored),
        "n_match_stored": n_match,
        "n_diff_stored": n_diff,
        "n_only_in_live": n_only_live,
        "n_only_in_stored": n_only_stored,
        "transitions_vs_stored": transitions,
        "sample_diffs": sample_diffs,
    }
    with open(SUMMARY_JSON, "w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True, default=str)
    print(f"Wrote summary {SUMMARY_JSON}", flush=True)

    print("\n=== SUMMARY ===")
    for k in ("n_total", "n_live", "n_stored", "n_match_stored", "n_diff_stored",
              "n_only_in_live", "n_only_in_stored"):
        print(f"  {k}: {summary[k]}")
    print("  top transitions:")
    for k, v in sorted(transitions.items(), key=lambda kv: -kv[1])[:15]:
        print(f"    {k}: {v}")
    print("  sample diffs:")
    for s in sample_diffs:
        print(f"    {s}")


if __name__ == "__main__":
    main()
