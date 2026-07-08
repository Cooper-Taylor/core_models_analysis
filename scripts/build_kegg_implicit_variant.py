#!/usr/bin/env python3
"""Build the "implicit KEGG model directions" reaction-direction *variant*.

The 5,683 KEGG core models already encode a direction for every reaction in their
on-disk COBRA flux bounds (lower_bound/upper_bound). This turns that as-built
direction into a reversibility "variant" the existing pipeline can display --
WITHOUT running the thermodynamic cascade -- exactly like
build_ai_direction_variant.py / build_consensus_variant.py.

Direction source: ``build_method_comparison.load_kegg_default()`` -- the
per-reaction majority on-disk bound direction across every model that contains
the reaction ((0,1000)=>, (-1000,0)=<, (-1000,1000)==; blocked ignored). Models
write reactions in MSDB-canonical orientation (verified: 0 flips), so the bound
direction maps straight to >/</= and is directly comparable to the cascade.

Like the other external overlays it is **baseline-overlaid**: a reaction keeps its
baseline cascade direction unless the models express one (237 in-model reactions).
This isolates the models' own directionality -- only those reactions diff vs
baseline / perturb FBA -- and keeps it apples-to-apples with the cascade baseline
the Variant Browser diffs against.

It writes ``thermo_variants/kegg_implicit/`` report files + ``cfg.json`` and
idempotently appends a ``kegg_implicit`` entry to ``thermo_variants/manifest.json``.
Everything downstream is then unchanged:

    build_site_data.py          -> site/data/manifest.json + variants/kegg_implicit.json
    build_all_models_impact.py  -> site/data/all_models_variant_fba__kegg_implicit.json
    site/serve.py + app.js      -> renders the variant like any other

IMPORTANT -- do NOT add this variant to ``variant_catalog.VARIANTS``: several
scripts call ``variant["cfg"]()`` unconditionally, and this variant has no
``ReversibilityConfig`` cascade. The standalone-script + manifest-append design
keeps the catalog pure.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Optional

SCRIPTS = Path(__file__).resolve().parent
ANALYSIS_ROOT = SCRIPTS.parent
sys.path.insert(0, str(SCRIPTS))

from build_ai_direction_variant import parse_report, write_report, count_new_rev
from build_method_comparison import load_kegg_default  # {seed_id: '>'/'<'/'='} from model bounds

OUT_ROOT = ANALYSIS_ROOT / "thermo_variants"
DEFAULT_BASELINE_REPORT = (
    OUT_ROOT / "baseline" / "Estimated_Reaction_Reversibility_Report_EQ.txt"
)


def main(argv: Optional[list] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--tag", default="kegg_implicit")
    ap.add_argument("--baseline-report", default=str(DEFAULT_BASELINE_REPORT),
                    help="baseline EQ report to overlay onto (default: %(default)s)")
    ap.add_argument("--out", default=str(OUT_ROOT),
                    help="thermo_variants output root (default: %(default)s)")
    args = ap.parse_args(argv)

    t0 = time.time()
    out_root = Path(args.out)
    dest = out_root / args.tag
    dest.mkdir(parents=True, exist_ok=True)

    base_rows = parse_report(Path(args.baseline_report))
    baseline_map = {rxn: nr for rxn, (_s, _o, nr) in base_rows.items()}
    print(f"baseline reactions: {len(baseline_map)}")

    kegg_map = load_kegg_default()  # majority on-disk bound direction per seed id
    print(f"KEGG implicit reactions: {len(kegg_map)}")

    # Overlay the model-built directions onto baseline, but only for reactions the
    # baseline (and thus the models/site) knows about.
    rows: dict = {}
    in_base = changed = 0
    for rxn, (status, _old, base_nr) in base_rows.items():
        new_rev, st = base_nr, status
        if rxn in kegg_map:
            in_base += 1
            new_rev = kegg_map[rxn]
            st = f"kegg: {new_rev}"
            if new_rev != base_nr:
                changed += 1
        # old_rev column carries the baseline direction (decorative for the site)
        rows[rxn] = (st, base_nr, new_rev)
    outside = sum(1 for r in kegg_map if r not in baseline_map)
    print(f"KEGG within baseline set: {in_base} | outside (ignored): {outside} "
          f"| direction changed vs baseline: {changed}")

    # Three reports (build_site_data reads EQ; GC/unfiltered written for parity).
    write_report(rows, dest / "Estimated_Reaction_Reversibility_Report_EQ.txt")
    write_report(rows, dest / "Estimated_Reaction_Reversibility_Report.txt")
    write_report(rows, dest / "Estimated_Reaction_Reversibility_Report_GC.txt",
                 drop_old_rev=True)

    counts_eq = {"total": len(rows),
                 **{f"new_rev={k}": v for k, v in count_new_rev(rows).items()}}

    summary = {
        "tag": args.tag,
        "title": "KEGG implicit directions (as-built model bounds)",
        "apt_title": ("The reaction directions the KEGG core models ship with "
                      "(on-disk flux bounds), overlaid on the baseline cascade"),
        "description": (
            "The direction each reaction was actually built with in the 5,683 KEGG "
            "core models, read straight from their on-disk COBRA flux bounds "
            "(lower_bound/upper_bound): (0,1000)->'>', (-1000,0)->'<', "
            "(-1000,1000)->'='. Taken as the per-reaction majority across every model "
            "containing the reaction (blocked bounds ignored; models are written in "
            "MSDB-canonical orientation, verified 0 flips). Applied as an overlay on "
            "the baseline cascade: a reaction keeps its baseline direction unless the "
            "models express one, so the diff and FBA effect isolate the models' own "
            "directionality."
        ),
        "citations": ["ModelSEED / KBase core models — on-disk flux bounds"],
        "section": "Model as-built",
        # Descriptive only -- NOT a ReversibilityConfig. build_site_data copies
        # this into the payload but never calls it.
        "cfg": {
            "_summary": "external KEGG on-disk direction map (no thermodynamic cascade)",
            "source": "core_models_kegg2 on-disk lower_bound/upper_bound",
            "overlay_on": "baseline",
            "n_kegg_total": len(kegg_map),
            "n_kegg_in_baseline": in_base,
            "n_kegg_outside_baseline": outside,
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

    # Idempotently append/replace the entry in manifest.json (create if absent).
    manifest_path = out_root / "manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
    else:
        manifest = {"variants": []}
    manifest.setdefault("variants", [])
    manifest["variants"] = [v for v in manifest["variants"] if v.get("tag") != args.tag]
    entry = dict(summary)
    entry["elapsed_s"] = round(time.time() - t0, 2)
    manifest["variants"].append(entry)
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2, default=str)

    print(f"wrote {dest}")
    print(f"updated {manifest_path} (tag={args.tag!r}; "
          f"{len(manifest['variants'])} variants total)")
    print(f"counts EQ: {counts_eq}")


if __name__ == "__main__":
    main()
