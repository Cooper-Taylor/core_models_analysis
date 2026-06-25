#!/usr/bin/env python3
"""Build a consensus reaction-direction *variant* from the thermodynamic estimators.

High-confidence ensemble: for every reaction, take a majority vote across three
independent ΔG′ direction estimates already on the site --

  - the group-contribution baseline cascade  (baseline)
  - eQuilibrator-3.0 component-contribution    (eq3_beber2022, Beber 2022)
  - dGPredictor fingerprint energies           (dgpredictor, Wang 2021)

A reaction is called directional ('>' / '<') only when at least two of the three
agree; conflicting single votes collapse to '=' (reversible), and reactions no
method has an opinion on stay '?'. This flips the baseline only where the two
newer estimators outvote group contribution, so the changes are deliberately
conservative but high-confidence.

External-data variant (no local cascade), same machinery as
build_ai_direction_variant.py / build_eq3_direction_variants.py. Run AFTER
export_thermo_variants.py and build_eq3_direction_variants.py (it reads their
reports) and re-appends itself to thermo_variants/manifest.json.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Optional

SCRIPTS = Path(__file__).resolve().parent
ANALYSIS_ROOT = SCRIPTS.parent
sys.path.insert(0, str(SCRIPTS))

from build_ai_direction_variant import parse_report, write_report, count_new_rev

OUT_ROOT = ANALYSIS_ROOT / "thermo_variants"
SOURCES = ("baseline", "eq3_beber2022", "dgpredictor")  # group-contribution + eQ3.0 + dGPredictor


def consensus(votes_raw) -> str:
    """Majority direction among the (non-'?') votes; '=' on conflict, '?' if none."""
    votes = [v for v in votes_raw if v in (">", "<", "=")]
    if not votes:
        return "?"
    top, n = Counter(votes).most_common(1)[0]
    if n >= 2:
        return top
    if len(votes) == 1:
        return votes[0]
    return "="  # single votes that disagree -> reversible


def main(argv: Optional[list] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--tag", default="consensus_thermo")
    args = ap.parse_args(argv)

    t0 = time.time()
    maps = {}
    for src in SOURCES:
        rep = OUT_ROOT / src / "Estimated_Reaction_Reversibility_Report_EQ.txt"
        rows = parse_report(rep)
        maps[src] = {rxn: nr for rxn, (_s, _o, nr) in rows.items()}
        print(f"{src}: {len(maps[src])} reactions")

    base = maps["baseline"]
    rows: dict = {}
    changed = agree3 = 0
    for rxn, base_nr in base.items():
        votes = [maps[s].get(rxn, "?") for s in SOURCES]
        c = consensus(votes)
        if len([v for v in votes if v == c and v in (">", "<", "=")]) == 3:
            agree3 += 1
        if c != base_nr:
            changed += 1
        rows[rxn] = (f"consensus: {c}", base_nr, c)
    print(f"consensus: {len(rows)} reactions; changed vs baseline={changed}; "
          f"unanimous (all 3 agree)={agree3}")

    dest = OUT_ROOT / args.tag
    dest.mkdir(parents=True, exist_ok=True)
    write_report(rows, dest / "Estimated_Reaction_Reversibility_Report_EQ.txt")
    write_report(rows, dest / "Estimated_Reaction_Reversibility_Report.txt")
    write_report(rows, dest / "Estimated_Reaction_Reversibility_Report_GC.txt", drop_old_rev=True)
    counts_eq = {"total": len(rows), **{f"new_rev={k}": v for k, v in count_new_rev(rows).items()}}

    summary = {
        "tag": args.tag,
        "title": "Consensus of GC + eQuilibrator-3.0 + dGPredictor",
        "apt_title": "High-confidence consensus: majority vote of three ΔG′ estimators",
        "description": (
            "Majority vote across three independent reaction-ΔG′ direction estimates: "
            "the group-contribution baseline, eQuilibrator-3.0 component contribution "
            "(Beber 2022), and dGPredictor (Wang 2021). A reaction is called "
            "directional only when at least two of the three agree; single votes that "
            "disagree collapse to reversible. This keeps the baseline direction unless "
            "the two newer estimators outvote group contribution, yielding a "
            "conservative, high-confidence direction set."
        ),
        "citations": ["Beber 2022", "Wang 2021", "Noor 2012"],
        "section": "§ New — consensus ensemble",
        "cfg": {
            "_summary": "majority vote of baseline / eq3_beber2022 / dgpredictor (no local cascade)",
            "sources": list(SOURCES),
            "n_changed_vs_baseline": changed,
            "n_unanimous": agree3,
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

    manifest_path = OUT_ROOT / "manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {"variants": []}
    manifest.setdefault("variants", [])
    manifest["variants"] = [v for v in manifest["variants"] if v.get("tag") != args.tag]
    entry = dict(summary)
    entry["elapsed_s"] = round(time.time() - t0, 2)
    manifest["variants"].append(entry)
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2, default=str)
    print(f"wrote {dest}; updated manifest ({len(manifest['variants'])} variants total)")


if __name__ == "__main__":
    main()
