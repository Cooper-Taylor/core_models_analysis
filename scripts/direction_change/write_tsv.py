#!/usr/bin/env python3
"""Write the direction-change target reactions out as TSV.

Two files, because "which reaction ids were in each model" has two readings:

``direction_change_targets_by_reaction.tsv``
    one row per target reaction (76), with per-panel present/active counts and
    flags. This is the summary that answers "which of the 76".

``direction_change_targets_by_model.tsv``
    one row per (panel, model, reaction) where the reaction is present, with a
    flag for whether it carried flux. This is the ground truth the summary is
    aggregated from.
"""
from __future__ import annotations

import csv
import glob
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from cma import paths  # noqa: E402

RESULTS = paths.results("direction_change_targets")
MSDB = paths.msdb("live")

PANELS = {
    "core_kegg2": "core",
    "ms2_gsm": "gsm",
    "ms2_gsm_auxo": "gsm_auxo",
}


def msdb_names(ids: set[str]) -> dict[str, str]:
    out = {}
    for f in sorted(glob.glob(str(MSDB / "Biochemistry" / "reaction_*.json"))):
        for r in json.load(open(f)):
            if r["id"] in ids:
                out[r["id"]] = (r.get("name") or "").strip()
    return out


def main() -> None:
    d = json.load(open(RESULTS / "presence_activity.json"))
    targets: dict[str, str] = d["targets"]
    panels: dict[str, dict] = d["panels"]
    names = msdb_names(set(targets))

    # ---- per-reaction summary ---------------------------------------------
    path = RESULTS / "direction_change_targets_by_reaction.tsv"
    cols = ["reaction_id", "change", "name"]
    for short in PANELS.values():
        cols += [f"present_{short}", f"growing_{short}", f"active_{short}"]
    cols += ["present_total", "growing_total", "active_total",
             "present_any", "active_any", "n_panels_present", "n_panels_active"]

    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for rid in sorted(targets):
            row = [rid, targets[rid], names.get(rid, "")]
            np_, na = 0, 0
            pa = aa = False
            tot_p = tot_g = tot_a = 0
            for key, short in PANELS.items():
                v = panels[key]["per_reaction"][rid]
                row += [v["present"], v.get("present_growing", ""), v["active_pfba"]]
                np_ += 1 if v["present"] else 0
                na += 1 if v["active_pfba"] else 0
                pa = pa or bool(v["present"])
                aa = aa or bool(v["active_pfba"])
                tot_p += v["present"]
                tot_g += v.get("present_growing", 0)
                tot_a += v["active_pfba"]
            row += [tot_p, tot_g, tot_a, int(pa), int(aa), np_, na]
            w.writerow(row)
    print(f"wrote {path}  ({len(targets)} rows)")

    # ---- per-model detail --------------------------------------------------
    detail = RESULTS / "per_model_records.json"
    if not detail.exists():
        print(f"NOTE: {detail.name} not found; per-model TSV skipped "
              f"(re-run presence_and_activity.py with --dump-records)")
        return
    recs = json.load(open(detail))
    path2 = RESULTS / "direction_change_targets_by_model.tsv"
    n = 0
    with open(path2, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["panel", "model_id", "grows", "reaction_id", "change",
                    "present", "active_pfba", "active_fba"])
        for panel_key, rows in recs.items():
            short = PANELS[panel_key]
            for rec in rows:
                act_p = set(rec["active_pfba"])
                act_f = set(rec["active_fba"])
                for rid in rec["present"]:
                    w.writerow([short, rec["model_id"], int(rec["grows"]), rid,
                                targets[rid], 1, int(rid in act_p), int(rid in act_f)])
                    n += 1
    print(f"wrote {path2}  ({n} rows)")


if __name__ == "__main__":
    main()
