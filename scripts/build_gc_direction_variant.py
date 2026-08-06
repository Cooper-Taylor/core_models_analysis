#!/usr/bin/env python3
"""Build a Group Contribution reaction-direction *variant* for the site.

Group Contribution (Jankowski 2008) is the original MSDB thermodynamic
default -- unlike dGPredictor and eQuilibrator-3.0, it was never promoted
past the 100-model descriptive panel (see results/thermo_sources/) to a
full 5,683-model site variant. This closes that gap using the exact same
overlay pattern as build_eq3_direction_variants.py / build_ai_direction_variant.py:
a reaction keeps its baseline cascade direction unless the Group Contribution
operator table has an opinion on it, so the diff / FBA effect isolates GC's
calls specifically.

Source: results/rxn_directions_group-contribution.csv (cols: rxn_id, operator,
dg, dge), produced by direction_pipeline.snapshot_msdb_per_source(source=
"Group contribution") reading Biochemistry/reaction_NN.json shards directly
from ModelSEEDDatabase (git show origin/dev:...), 25,812 reactions with a
non-sentinel dG value out of ~56,000 GC-bearing entries.

Writes thermo_variants/group_contribution/ reports + cfg.json and idempotently
appends "group_contribution" to thermo_variants/manifest.json; everything
downstream (build_site_data.py / build_all_models_impact.py) then picks it up
like any other variant.

IMPORTANT: like the AI and eq3 variants, this is NOT added to
variant_catalog.VARIANTS (it has no ReversibilityConfig cascade).
"""

from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path
from typing import Optional

SCRIPTS = Path(__file__).resolve().parent
ANALYSIS_ROOT = SCRIPTS.parent
sys.path.insert(0, str(SCRIPTS))

from build_ai_direction_variant import parse_report, write_report, count_new_rev

OUT_ROOT = ANALYSIS_ROOT / "thermo_variants"
DEFAULT_BASELINE_REPORT = (
    OUT_ROOT / "baseline" / "Estimated_Reaction_Reversibility_Report_EQ.txt"
)
DEFAULT_GC_CSV = ANALYSIS_ROOT / "results" / "rxn_directions_group-contribution.csv"


def load_gc(csv_path: Path) -> dict:
    """Return {rxn: operator} from the Group Contribution per-reaction table."""
    directions = {}
    with open(csv_path) as fh:
        for row in csv.DictReader(fh):
            rid = row["rxn_id"]
            d = (row.get("operator") or "").strip()
            if d in (">", "<", "=", "?"):
                directions[rid] = d
    return directions


def build_one(tag, dir_map, base_rows, summary_extra):
    dest = OUT_ROOT / tag
    dest.mkdir(parents=True, exist_ok=True)
    rows, in_base, changed = {}, 0, 0
    for rxn, (status, _old, base_nr) in base_rows.items():
        new_rev, st = base_nr, status
        if rxn in dir_map:
            in_base += 1
            new_rev = dir_map[rxn]
            st = f"gc: {new_rev}"
            if new_rev != base_nr:
                changed += 1
        rows[rxn] = (st, base_nr, new_rev)
    outside = sum(1 for r in dir_map if r not in base_rows)
    write_report(rows, dest / "Estimated_Reaction_Reversibility_Report_EQ.txt")
    write_report(rows, dest / "Estimated_Reaction_Reversibility_Report.txt")
    write_report(rows, dest / "Estimated_Reaction_Reversibility_Report_GC.txt", drop_old_rev=True)
    counts_eq = {"total": len(rows), **{f"new_rev={k}": v for k, v in count_new_rev(rows).items()}}
    summary = {
        **summary_extra,
        "tag": tag,
        "cfg": {
            "_summary": "external Group Contribution direction map (no local cascade)",
            "source_csv": str(DEFAULT_GC_CSV),
            "overlay_on": "baseline",
            "n_gc_total": len(dir_map),
            "n_gc_in_baseline": in_base,
            "n_gc_outside_baseline": outside,
            "n_changed_vs_baseline": changed,
        },
        "files": {
            "EQ": "Estimated_Reaction_Reversibility_Report_EQ.txt",
            "GC": "Estimated_Reaction_Reversibility_Report_GC.txt",
            "unfiltered": "Estimated_Reaction_Reversibility_Report.txt",
        },
        "counts": {"EQ": counts_eq, "GC": dict(counts_eq), "unfiltered": dict(counts_eq)},
    }
    with open(dest / "cfg.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    print(f"[{tag}] in_baseline={in_base} outside={outside} changed_vs_baseline={changed}")
    return summary


def append_manifest(summaries: list, t0: float) -> None:
    manifest_path = OUT_ROOT / "manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {"variants": []}
    manifest.setdefault("variants", [])
    tags = {s["tag"] for s in summaries}
    manifest["variants"] = [v for v in manifest["variants"] if v.get("tag") not in tags]
    for s in summaries:
        entry = dict(s)
        entry["elapsed_s"] = round(time.time() - t0, 2)
        manifest["variants"].append(entry)
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2, default=str)
    print(f"updated {manifest_path} (+{[s['tag'] for s in summaries]}; "
          f"{len(manifest['variants'])} variants total)")


def main(argv: Optional[list] = None) -> None:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--csv", default=str(DEFAULT_GC_CSV))
    ap.add_argument("--baseline-report", default=str(DEFAULT_BASELINE_REPORT))
    args = ap.parse_args(argv)

    t0 = time.time()
    base_rows = parse_report(Path(args.baseline_report))
    print(f"baseline reactions: {len(base_rows)}")
    directions = load_gc(Path(args.csv))
    print(f"gc rows: directions={len(directions)}")

    summaries = [
        build_one(
            "group_contribution", directions, base_rows,
            {
                "title": "Group Contribution directions (Jankowski 2008)",
                "apt_title": "Group-contribution method reaction directionality "
                             "(Jankowski 2008), MSDB default thermodynamic source",
                "description": (
                    "Reaction directions from the Group Contribution method "
                    "(Jankowski 2008), the thermodynamic source MSDB shipped "
                    "with prior to eQuilibrator/dGPredictor integration. Read "
                    "directly from each reaction's bundled 'Group contribution' "
                    "dG/operator triple (results/rxn_directions_group-contribution.csv, "
                    "25,812 of ~56,000 reactions with a non-sentinel estimate). "
                    "Overlaid on baseline."
                ),
                "citations": ["Jankowski 2008"],
                "section": "§ Thermodynamic source comparison",
            },
        ),
    ]
    append_manifest(summaries, t0)


if __name__ == "__main__":
    main()
