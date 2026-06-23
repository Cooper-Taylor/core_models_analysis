#!/usr/bin/env python3
"""Build an AI-predicted reaction-directionality *variant* for the site.

Standalone producer that turns an ``AICurationCacheReactionDirectionality.json``
(reaction directions predicted by Claude Opus 4.8 via Argo) into a reversibility
"variant" the existing pipeline can display -- WITHOUT running the thermodynamic
cascade.

It writes ``thermo_variants/<tag>/`` report files in the exact format produced by
``export_thermo_variants.write_report`` (sourcing ``new_rev`` from the AI map),
plus ``cfg.json``, and idempotently appends a ``<tag>`` entry to
``thermo_variants/manifest.json``. Everything downstream is then unchanged:

    build_site_data.py          -> site/data/manifest.json + variants/<tag>.json
    build_all_models_impact.py  -> site/data/all_models_variant_fba__<tag>.json
    site/serve.py + app.js      -> renders the variant like any other

The AI variant is **baseline-overlaid**: every reaction keeps its baseline
direction unless the AI expressed an opinion on it. This isolates the AI's
effect (only AI-driven changes diff vs baseline / perturb FBA), matching the
``override_bounds(..., baseline_map=...)`` semantics used elsewhere.

IMPORTANT -- do NOT add this variant to ``variant_catalog.VARIANTS``: six scripts
(export_thermo_variants, eval_defaults_panel, run_statistical_panel [x3],
run_variant_source_panel) call ``variant["cfg"]()`` unconditionally, and an AI
variant has no ``ReversibilityConfig`` cascade. The standalone-script +
manifest-append design keeps the catalog pure and those scripts working.

Direction mapping (via direction_pipeline._normalize_direction):
    forward -> '>'   reverse -> '<'   reversible -> '='   uncertain -> '?'
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

from direction_pipeline import _normalize_direction  # forward/reverse/... -> >/</=/?

OUT_ROOT = ANALYSIS_ROOT / "thermo_variants"
DEFAULT_BASELINE_REPORT = (
    OUT_ROOT / "baseline" / "Estimated_Reaction_Reversibility_Report_EQ.txt"
)
DEFAULT_AI_JSON = (
    ANALYSIS_ROOT / "data" / "ai_curation" / "all_modelseed"
    / "AICurationCacheReactionDirectionality.json"
)


# ---------------------------------------------------------------------------
def parse_report(path: Path) -> dict:
    """Parse a 4-col MSDB EQ report -> {rxn: (status, old_rev, new_rev)}."""
    out = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4:
                out[parts[0]] = (parts[1], parts[2], parts[3])
    return out


def load_ai_map(ai_json_path: Path, drop_errors: bool = False) -> tuple[dict, int]:
    """Load {rxn: direction-char} from an AICurationCache directionality JSON."""
    data = json.loads(Path(ai_json_path).read_text())
    ai_map: dict = {}
    skipped_err = 0
    for rxn, entry in data.items():
        if not isinstance(entry, dict):
            continue
        if drop_errors and entry.get("errors"):
            skipped_err += 1
            continue
        ai_map[rxn] = _normalize_direction(entry.get("directionality"))
    return ai_map, skipped_err


def write_report(rows: dict, path: Path, drop_old_rev: bool = False) -> None:
    """Write rows {rxn: (status, old_rev, new_rev)} in MSDB report format.

    Mirrors export_thermo_variants.write_report: EQ + unfiltered keep 3 fields
    (status, old_rev, new_rev); GC drops old_rev. Sorted, tab-separated.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for rxn in sorted(rows):
            status, old_rev, new_rev = rows[rxn]
            row = [status, old_rev, new_rev]
            if drop_old_rev:
                del row[1]
            fh.write(rxn + "\t" + "\t".join(row) + "\n")


def count_new_rev(rows: dict) -> dict:
    c: dict = {}
    for _status, _old, new_rev in rows.values():
        c[new_rev] = c.get(new_rev, 0) + 1
    return c


# ---------------------------------------------------------------------------
def main(argv: Optional[list] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ai-json", default=str(DEFAULT_AI_JSON),
                    help="AICurationCacheReactionDirectionality.json (default: %(default)s)")
    ap.add_argument("--tag", default="ai_opus48")
    ap.add_argument("--baseline-report", default=str(DEFAULT_BASELINE_REPORT),
                    help="baseline EQ report to overlay onto (default: %(default)s)")
    ap.add_argument("--out", default=str(OUT_ROOT),
                    help="thermo_variants output root (default: %(default)s)")
    ap.add_argument("--title", default="AI-predicted directions (Claude Opus 4.8)")
    ap.add_argument("--apt-title",
                    default="Reaction directions predicted by Claude Opus 4.8 (via Argo), "
                            "overlaid on the MSDB baseline")
    ap.add_argument("--description",
                    default="Reaction directions predicted by the Claude Opus 4.8 model "
                            "(through the Argo gateway) from each reaction's stoichiometry. "
                            "Applied as an overlay on the baseline cascade: a reaction keeps "
                            "its baseline direction unless the AI expressed an opinion, so the "
                            "diff and FBA effect isolate the AI's calls. "
                            "forward->'>', reverse->'<', reversible->'=', uncertain->'?'.")
    ap.add_argument("--section", default="AI / LLM")
    ap.add_argument("--citation", action="append", default=None,
                    help="repeatable; defaults to a single Claude Opus 4.8 citation")
    ap.add_argument("--drop-errors", action="store_true",
                    help="skip AI entries whose 'errors' list is non-empty")
    args = ap.parse_args(argv)

    t0 = time.time()
    out_root = Path(args.out)
    dest = out_root / args.tag
    dest.mkdir(parents=True, exist_ok=True)

    base_rows = parse_report(Path(args.baseline_report))
    baseline_map = {rxn: nr for rxn, (_s, _o, nr) in base_rows.items()}
    print(f"baseline reactions: {len(baseline_map)}")

    ai_map, skipped_err = load_ai_map(Path(args.ai_json), drop_errors=args.drop_errors)
    print(f"AI reactions: {len(ai_map)} (errors skipped: {skipped_err})")

    # Overlay AI directions onto baseline, but only for reactions the baseline
    # (and thus the models/site) knows about.
    rows: dict = {}
    in_base = changed = 0
    for rxn, (status, _old, base_nr) in base_rows.items():
        new_rev, st = base_nr, status
        if rxn in ai_map:
            in_base += 1
            new_rev = ai_map[rxn]
            st = f"AI: {new_rev}"
            if new_rev != base_nr:
                changed += 1
        # old_rev column carries the baseline direction (decorative for the site)
        rows[rxn] = (st, base_nr, new_rev)
    outside = sum(1 for r in ai_map if r not in baseline_map)
    print(f"AI within baseline set: {in_base} | outside (ignored): {outside} "
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
        "title": args.title,
        "apt_title": args.apt_title,
        "description": args.description,
        "citations": args.citation or ["Claude Opus 4.8 (Anthropic), via the Argo gateway"],
        "section": args.section,
        # Descriptive only -- NOT a ReversibilityConfig. build_site_data copies
        # this into the payload but never calls it.
        "cfg": {
            "_summary": "external AI direction map (no thermodynamic cascade)",
            "source_json": str(args.ai_json),
            "model": "claudeopus48 (Claude Opus 4.8 via Argo)",
            "overlay_on": "baseline",
            "n_ai_total": len(ai_map),
            "n_ai_in_baseline": in_base,
            "n_ai_outside_baseline": outside,
            "n_changed_vs_baseline": changed,
            "drop_errors": bool(args.drop_errors),
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
