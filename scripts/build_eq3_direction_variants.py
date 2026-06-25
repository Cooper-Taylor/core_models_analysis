#!/usr/bin/env python3
"""Build eQuilibrator-3.0 (Beber 2022) reaction-direction *variants* for the site.

Two external-data variants, overlaid on the baseline cascade (same pattern as
``build_ai_direction_variant.py`` -- a reaction keeps its baseline direction
unless eQ-3.0 has an opinion, so the diff / FBA effect isolates eQ-3.0's calls):

  - ``eq3_beber2022`` -- eQuilibrator-3.0 component-contribution directions
    (Beber 2022) at the standard reversibility-index cutoff |log10 gamma| = 3,
    read straight from the ``Beber_2022`` column.
  - ``eq3_gamma1``    -- the same eQ-3.0 energies but with an aggressive cutoff
    |log10 gamma| = 1, recomputed from the per-reaction ``ln_RI`` column
    (Noor 2012/2013): '>' if ln_RI < -ln(10), '<' if > +ln(10), else '='.

Source: ``results/rxn_directions_eq3_2022.tsv`` (cols: rxn_id, Beber_2022, ln_RI)
produced earlier by ``estimate_directions_eq3.py`` using equilibrator-api 0.7.0
(= eQuilibrator 3.0 / component-contribution, Beber 2022).

Writes ``thermo_variants/<tag>/`` reports + cfg.json and idempotently appends each
``<tag>`` to ``thermo_variants/manifest.json``; everything downstream
(build_site_data / build_all_models_impact / build_panel_rxn_pipeline) then
picks them up like any other variant.

IMPORTANT: like the AI variant, these are NOT added to ``variant_catalog.VARIANTS``
(they have no ReversibilityConfig cascade). Run this AFTER
``export_thermo_variants.py`` (which rewrites the manifest from the catalog only).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from pathlib import Path
from typing import Optional

SCRIPTS = Path(__file__).resolve().parent
ANALYSIS_ROOT = SCRIPTS.parent
sys.path.insert(0, str(SCRIPTS))

# Reuse the AI-variant overlay machinery (identical report/manifest format).
from build_ai_direction_variant import parse_report, write_report, count_new_rev

OUT_ROOT = ANALYSIS_ROOT / "thermo_variants"
DEFAULT_BASELINE_REPORT = (
    OUT_ROOT / "baseline" / "Estimated_Reaction_Reversibility_Report_EQ.txt"
)
EQ3_TSV = ANALYSIS_ROOT / "results" / "rxn_directions_eq3_2022.tsv"

LN10 = math.log(10.0)  # |log10 gamma| = 1 threshold, in nat units


def load_eq3(tsv: Path) -> tuple[dict, dict]:
    """Return ({rxn: Beber_2022 direction}, {rxn: ln_RI})."""
    directions, lnri = {}, {}
    with open(tsv) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            rid = row["rxn_id"]
            d = (row.get("Beber_2022") or "").strip()
            if d in (">", "<", "=", "?"):
                directions[rid] = d
            raw = (row.get("ln_RI") or "").strip()
            try:
                lnri[rid] = float(raw)
            except (TypeError, ValueError):
                pass
    return directions, lnri


def dir_from_lnri(ln_ri: float, thresh: float) -> str:
    """eQuilibrator reversibility-index direction call (Noor 2012/2013).

    ln gamma carries the SAME sign as ΔG′ₘ, so ln_RI < -thresh -> forward,
    > +thresh -> reverse, else reversible."""
    if ln_ri < -thresh:
        return ">"
    if ln_ri > thresh:
        return "<"
    return "="


def build_one(tag, dir_map, base_rows, summary_extra):
    """Overlay {rxn: dir} onto baseline; write reports + cfg.json; return summary."""
    dest = OUT_ROOT / tag
    dest.mkdir(parents=True, exist_ok=True)
    rows, in_base, changed = {}, 0, 0
    for rxn, (status, _old, base_nr) in base_rows.items():
        new_rev, st = base_nr, status
        if rxn in dir_map:
            in_base += 1
            new_rev = dir_map[rxn]
            st = f"eq3: {new_rev}"
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
            "_summary": "external eQuilibrator-3.0 direction map (no local cascade)",
            "source_tsv": str(EQ3_TSV),
            "overlay_on": "baseline",
            "n_eq3_total": len(dir_map),
            "n_eq3_in_baseline": in_base,
            "n_eq3_outside_baseline": outside,
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
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--tsv", default=str(EQ3_TSV))
    ap.add_argument("--baseline-report", default=str(DEFAULT_BASELINE_REPORT))
    args = ap.parse_args(argv)

    t0 = time.time()
    base_rows = parse_report(Path(args.baseline_report))
    print(f"baseline reactions: {len(base_rows)}")
    directions, lnri = load_eq3(Path(args.tsv))
    print(f"eq3 rows: directions={len(directions)} ln_RI={len(lnri)}")

    gamma1 = {rid: dir_from_lnri(v, LN10) for rid, v in lnri.items()}

    summaries = [
        build_one(
            "eq3_beber2022", directions, base_rows,
            {
                "title": "eQuilibrator 3.0 directions (Beber 2022)",
                "apt_title": "eQuilibrator-3.0 component-contribution directionality "
                             "(Beber 2022), |log10 gamma| = 3",
                "description": (
                    "Reaction directions from eQuilibrator 3.0's component-contribution "
                    "energies (Beber 2022) using its reversibility index at the standard "
                    "|log10 gamma| = 3 cutoff (Noor 2012/2013). eQuilibrator 3.0 is the "
                    "current state-of-the-art ΔG′ estimator (full covariance, Mg2+/pH "
                    "corrections), so this swaps the baseline's bundled group-contribution "
                    "energies for the newest data. Overlaid on baseline."
                ),
                "citations": ["Beber 2022", "Noor 2012", "Noor 2013"],
                "section": "§ New — eQuilibrator 3.0",
            },
        ),
        build_one(
            "eq3_gamma1", gamma1, base_rows,
            {
                "title": "eQuilibrator 3.0 directions, |log10 gamma| = 1 (aggressive)",
                "apt_title": "eQuilibrator-3.0 energies at an aggressive 10-fold cutoff "
                             "(Beber 2022 + Noor 2013)",
                "description": (
                    "The same eQuilibrator-3.0 component-contribution energies as "
                    "eq3_beber2022, but with an aggressive reversibility-index cutoff "
                    "|log10 gamma| = 1 (a reaction is directional once a ~10-fold "
                    "concentration shift would reverse it). Combines the newest energy "
                    "data with the most aggressive threshold, so it makes the largest "
                    "direction change of the new set. Overlaid on baseline."
                ),
                "citations": ["Beber 2022", "Noor 2012", "Noor 2013"],
                "section": "§ New — eQuilibrator 3.0",
            },
        ),
    ]
    append_manifest(summaries, t0)


if __name__ == "__main__":
    main()
