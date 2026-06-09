#!/usr/bin/env python3
"""Export every ReversibilityConfig variant in MSDB report format.

For each variant in ``variant_catalog.VARIANTS`` we write three files
that mirror the layout produced by
``ModelSEEDDatabase/Scripts/Thermodynamics/Estimate_Reaction_Reversibility.py``:

    Estimated_Reaction_Reversibility_Report_EQ.txt   (4 cols, EQ-level)
    Estimated_Reaction_Reversibility_Report_GC.txt   (3 cols, GC-level)
    Estimated_Reaction_Reversibility_Report.txt      (4 cols, unfiltered)

Plus a small ``cfg.json`` capturing the knob settings, and a top-level
``manifest.json`` indexing all variants.

Outputs live under
``core_models_analysis/thermo_variants/{tag}/`` -- *never* under
``ModelSEEDDatabase/`` (the source of truth is left untouched per the
project policy).

The script reads the cached MSDB reactions blob produced by notebook 06
(``msdb_reactions_v1``) so a cold start without any cache still works
(it loads via BiochemPy in that case, ~45s).
"""

from __future__ import annotations

import argparse
import copy
import dataclasses
import json
import os
import sys
import time
from pathlib import Path
from typing import Optional

SCRIPTS = Path(__file__).resolve().parent
ANALYSIS_ROOT = SCRIPTS.parent
MSDB_ROOT = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase"))

sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(MSDB_ROOT / "Libs" / "Python"))

import reversibility_lib as lib
import variant_catalog


OUT_ROOT = ANALYSIS_ROOT / "thermo_variants"


# ---------------------------------------------------------------------------
# Reaction database load (cached if notebook 06 was ever run)
# ---------------------------------------------------------------------------
def load_msdb_reactions(use_cache: bool = True) -> dict:
    """Load MSDB reactions, preferring the notebook 06 cache."""
    if use_cache:
        try:
            from kbutillib.notebook import NotebookSession
            session = NotebookSession.for_notebook(
                notebook_file=str(ANALYSIS_ROOT / "notebooks" /
                                  "06_ReactionReversibilityHeuristics.ipynb"),
                project_name="core_models_analysis",
            )
            return session.cache.load("msdb_reactions_v1")
        except Exception as exc:
            print(f"[load_msdb_reactions] cache miss ({type(exc).__name__}): "
                  f"falling back to BiochemPy load")

    from BiochemPy import Reactions
    return Reactions().loadReactions()


# ---------------------------------------------------------------------------
# Cascade runner that captures both the GC-pass and the requested-level
# pass, so we can write the ``old_rev`` column of the EQ / unfiltered
# reports exactly the way ``Estimate_Reaction_Reversibility.main`` does.
# ---------------------------------------------------------------------------
def run_one_level(rxns: dict, cfg: lib.ReversibilityConfig, db_level: str,
                  gc_first: bool = True) -> dict:
    """Return ``{rxn_id: (status, old_rev, new_rev, source_label)}``.

    Faithfully mirrors ``Estimate_Reaction_Reversibility.main``:
      * for db_level == 'EQ' with gc_first, the GC pass runs first and
        mutates ``rxn_entry['reversibility']`` so the EQ pass's
        ``Incomplete (GCC)`` fallback returns the right value, AND so the
        report's third column carries the *prior* GC-derived direction;
      * for db_level == 'GC', old_rev comes from the original JSON's
        ``reversibility`` field (which is also what the MSDB GC report
        encodes, modulo the writer dropping that column);
      * for db_level == '' (unfiltered), no gc_first preroll: old_rev is
        the original JSON's ``reversibility``.
    """
    work = copy.deepcopy(rxns)

    if db_level == "EQ" and gc_first:
        # GC pre-roll: estimate each reaction under GC and overwrite
        # work[rxn]['reversibility'] so the EQ pass's Incomplete-GCC
        # fallback can read it.  This is exactly what the cached
        # baseline_cascade did -- we replay it here so the report rows
        # come out the same.
        for rxn_id in sorted(work):
            entry = work[rxn_id]
            _, rev, _ = lib.estimate_one(entry, "GC", cfg)
            entry["reversibility"] = rev

    out = {}
    for rxn_id in sorted(work):
        entry = work[rxn_id]
        old_rev = entry["reversibility"]
        status, new_rev, src = lib.estimate_one(entry, db_level, cfg)
        out[rxn_id] = (status, old_rev, new_rev, src)
    return out


def write_report(rows: dict, path: Path, drop_old_rev: bool = False) -> None:
    """Write ``rows`` in the MSDB report format (sorted, tab-separated).

    Matches ``Estimate_Reaction_Reversibility._write_report``:
      * GC report has 3 columns (status, new_rev) -- old_rev dropped;
      * EQ + unfiltered have 4 columns (status, old_rev, new_rev).
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for rxn in sorted(rows):
            status, old_rev, new_rev, _src = rows[rxn]
            row = [status, old_rev, new_rev]
            if drop_old_rev:
                del row[1]
            fh.write(rxn + "\t" + "\t".join(row) + "\n")


def cfg_to_dict(cfg: lib.ReversibilityConfig) -> dict:
    """Lossy but human-readable dump of a ReversibilityConfig.

    Collections (per-met concentration tables, ln_RI map) are summarized
    as their length to keep the JSON small.
    """
    d = {}
    for f in dataclasses.fields(cfg):
        v = getattr(cfg, f.name)
        if isinstance(v, dict):
            d[f.name] = {"_summary": "dict", "n_entries": len(v)}
        elif isinstance(v, (list, tuple)) and len(v) > 16:
            d[f.name] = {"_summary": "sequence", "n_entries": len(v)}
        else:
            d[f.name] = v
    return d


# ---------------------------------------------------------------------------
# Per-variant export driver
# ---------------------------------------------------------------------------
def export_variant(rxns: dict, variant: dict, out_root: Path,
                   levels: tuple = ("EQ", "GC", "")) -> dict:
    """Write the three MSDB-format reports + cfg.json for one variant.

    Returns a summary dict suitable for ``manifest.json``.
    """
    tag = variant["tag"]
    cfg = variant["cfg"]()
    dest = out_root / tag
    dest.mkdir(parents=True, exist_ok=True)

    summary = {
        "tag": tag,
        "title": variant["title"],
        "apt_title": variant.get("apt_title", variant["title"]),
        "description": variant.get("description", ""),
        "citations": variant.get("citations", []),
        "section": variant["section"],
        "cfg": cfg_to_dict(cfg),
        "files": {},
        "counts": {},
    }
    t0 = time.time()
    for level in levels:
        rows = run_one_level(rxns, cfg, level, gc_first=(level == "EQ"))
        if level == "":
            fname = "Estimated_Reaction_Reversibility_Report.txt"
            drop = False
            label = "unfiltered"
        elif level == "GC":
            fname = "Estimated_Reaction_Reversibility_Report_GC.txt"
            drop = True
            label = "GC"
        else:
            fname = f"Estimated_Reaction_Reversibility_Report_{level}.txt"
            drop = False
            label = level
        write_report(rows, dest / fname, drop_old_rev=drop)
        summary["files"][label] = fname
        summary["counts"][label] = {
            "total": len(rows),
            **{f"new_rev={k}": v for k, v in _count_field(rows, 2).items()},
        }

    with open(dest / "cfg.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=str)

    summary["elapsed_s"] = round(time.time() - t0, 2)
    return summary


def _count_field(rows: dict, field_idx: int) -> dict:
    counts = {}
    for r in rows.values():
        v = r[field_idx]
        counts[v] = counts.get(v, 0) + 1
    return counts


# ---------------------------------------------------------------------------
def main(argv: Optional[list] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(OUT_ROOT),
                    help="output root (default: %(default)s)")
    ap.add_argument("--only", action="append", default=None,
                    help="only export variant with this tag (repeatable)")
    ap.add_argument("--no-cache", action="store_true",
                    help="bypass the notebook 06 kbcache; reload via BiochemPy")
    args = ap.parse_args(argv)

    out_root = Path(args.out)
    out_root.mkdir(parents=True, exist_ok=True)

    print(f"loading MSDB reactions (cache: {not args.no_cache})...")
    t = time.time()
    rxns = load_msdb_reactions(use_cache=not args.no_cache)
    print(f"  {len(rxns)} reactions in {time.time()-t:.1f}s")

    selected = (variant_catalog.VARIANTS if args.only is None
                else [variant_catalog.variant_by_tag(t) for t in args.only])
    manifest = {"variants": []}
    for v in selected:
        print(f"--- variant {v['tag']:>12} : {v['title']}")
        s = export_variant(rxns, v, out_root)
        for label, c in s["counts"].items():
            extras = {k: v for k, v in c.items() if k != "total"}
            print(f"      {label:>10}  total={c['total']}  {extras}")
        print(f"      elapsed: {s['elapsed_s']}s")
        manifest["variants"].append(s)

    with open(out_root / "manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2, default=str)
    print(f"wrote {out_root / 'manifest.json'}  "
          f"({len(manifest['variants'])} variants)")


if __name__ == "__main__":
    main()
