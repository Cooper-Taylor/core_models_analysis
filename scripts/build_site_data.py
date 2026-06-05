#!/usr/bin/env python3
"""Build static JSON snapshots for the heuristics-explorer website.

Reads:
  - ``thermo_variants/manifest.json`` + per-variant report files
    (produced by ``export_thermo_variants.py``)
  - ``notebooks/.kbcache/`` blobs (msdb_reactions_v1,
    fba_*, reversibility_variant_*, rxnsets_by_model)
  - ``data/core_models_kegg2/{model_id}.json`` (cobra panel models)

Writes ``site/data/``:
  - ``manifest.json``                    -- variants + panel summary
  - ``baseline.json``                    -- baseline cascade (status + rev)
  - ``variants/{tag}.json``              -- diff vs baseline + per-panel FBA
  - ``reactions.json``                   -- compact rxn index (name, eq, ec,
                                            which variants change it, which
                                            panel models contain it)
  - ``panel.json``                       -- model ids and per-model summary

After running, ``site/serve.py`` exposes the data through the
visualization frontend plus a live FBA endpoint for per-reaction
"what if I flip this off / reverse it" exploration.
"""

from __future__ import annotations

import argparse
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

THERMO_VARIANTS_DIR = ANALYSIS_ROOT / "thermo_variants"
SITE_DATA = ANALYSIS_ROOT / "site" / "data"
MODELS_DIR = ANALYSIS_ROOT / "data" / "core_models_kegg2"
PANEL_FILE = ANALYSIS_ROOT / "results" / "selected_ids.txt"
STATS_DIR = ANALYSIS_ROOT / "results" / "statistical_panel"


# ---------------------------------------------------------------------------
def load_session():
    from kbutillib.notebook import NotebookSession
    return NotebookSession.for_notebook(
        notebook_file=str(ANALYSIS_ROOT / "notebooks" /
                          "06_ReactionReversibilityHeuristics.ipynb"),
        project_name="core_models_analysis",
    )


def parse_report(path: Path, drop_old_rev: bool) -> dict:
    """Parse a MSDB-format report.

    Returns ``{rxn_id: {"status": str, "old_rev": str | None, "new_rev": str}}``.
    GC reports lack the old_rev column so it's None there.
    """
    out = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if drop_old_rev:
                rxn, status, new = parts[0], parts[1], parts[2]
                out[rxn] = {"status": status, "old_rev": None, "new_rev": new}
            else:
                rxn, status, old, new = parts[0], parts[1], parts[2], parts[3]
                out[rxn] = {"status": status, "old_rev": old, "new_rev": new}
    return out


def variant_diff(baseline_map: dict, variant_map: dict) -> dict:
    """``baseline_map`` and ``variant_map`` are ``{rxn: new_rev}``."""
    diffs = []
    transitions = {}
    for rxn in baseline_map:
        b = baseline_map[rxn]
        v = variant_map.get(rxn)
        if v is not None and b != v:
            diffs.append({"rxn": rxn, "base": b, "new": v})
            key = f"{b}->{v}"
            transitions[key] = transitions.get(key, 0) + 1
    return {"diffs": diffs, "transitions": transitions, "n_changed": len(diffs)}


# ---------------------------------------------------------------------------
def build_reactions_index(msdb_rxns: dict, panel_rxn_sets: dict,
                          changed_by_variant: dict) -> dict:
    """Build the per-reaction lookup used by the site reaction explorer.

    Keys: msdb reaction ids that appear *either* in any panel model *or*
    in at least one variant's diff vs baseline.  Reactions absent from
    both lists are excluded -- they wouldn't surface in the UI anyway and
    they would bloat the payload by an order of magnitude.
    """
    panel_rxn_union = set()
    for s in panel_rxn_sets.values():
        panel_rxn_union.update(s)

    surfacing = set(panel_rxn_union)
    for diffs in changed_by_variant.values():
        for d in diffs:
            surfacing.add(d["rxn"])

    out = {}
    for rxn_id in sorted(surfacing):
        r = msdb_rxns.get(rxn_id)
        if r is None:
            out[rxn_id] = {
                "id": rxn_id,
                "name": "(unknown)",
                "equation": "",
                "definition": "",
                "in_panel": rxn_id in panel_rxn_union,
                "panel_freq": 0,
                "is_transport": None,
                "deltag": None,
                "deltagerr": None,
                "ec_numbers": [],
                "pathways": [],
                "stoichiometry": [],
                "changed_by": [],
            }
            continue
        out[rxn_id] = {
            "id": rxn_id,
            "name": r.get("name"),
            "equation": r.get("equation"),
            "definition": r.get("definition"),
            "in_panel": rxn_id in panel_rxn_union,
            "panel_freq": sum(1 for s in panel_rxn_sets.values() if rxn_id in s),
            "is_transport": int(r.get("is_transport") or 0),
            "deltag": r.get("deltag"),
            "deltagerr": r.get("deltagerr"),
            "ec_numbers": r.get("ec_numbers") or [],
            "pathways": r.get("pathways") or [],
            "stoichiometry": [
                {"cpd": s["compound"], "name": s.get("name"),
                 "coef": s["coefficient"], "cpt": s.get("compartment", 0),
                 "formula": s.get("formula")}
                for s in (r.get("stoichiometry") or [])
            ],
            "changed_by": [],  # filled in below
        }

    for tag, diffs in changed_by_variant.items():
        for d in diffs:
            if d["rxn"] in out:
                out[d["rxn"]]["changed_by"].append(
                    {"variant": tag, "base": d["base"], "new": d["new"]})
    return out


# ---------------------------------------------------------------------------
def build_panel_info(panel_ids: list, ondisk_fba: list,
                     heur_baseline_fba: list, panel_rxn_sets: dict) -> dict:
    """Per-model summary table for the panel browser."""
    odi = {r["model_id"]: r for r in ondisk_fba}
    hbi = {r["model_id"]: r for r in heur_baseline_fba}
    out = []
    for mid in panel_ids:
        n_rxns = len(panel_rxn_sets.get(mid, []))
        o = odi.get(mid, {})
        h = hbi.get(mid, {})
        out.append({
            "model_id": mid,
            "n_rxns": n_rxns,
            "ondisk_grows": bool(o.get("grows")),
            "ondisk_flux": float(o.get("growth_flux", 0.0)),
            "heur_baseline_grows": bool(h.get("grows")),
            "heur_baseline_flux": float(h.get("growth_flux", 0.0)),
        })
    return out


# ---------------------------------------------------------------------------
def main(argv: Optional[list] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(SITE_DATA))
    args = ap.parse_args(argv)

    out_root = Path(args.out)
    (out_root / "variants").mkdir(parents=True, exist_ok=True)

    print("loading kbcache...")
    s = load_session()
    msdb_rxns = s.cache.load("msdb_reactions_v1")
    print(f"  {len(msdb_rxns)} reactions")
    panel_rxn_sets = s.cache.load("rxnsets_by_model")
    print(f"  {len(panel_rxn_sets)} model rxn-sets")

    panel_ids = PANEL_FILE.read_text().split()
    print(f"  panel: {len(panel_ids)} models")
    ondisk_fba = s.cache.load("fba_ondisk_v1")
    heur_baseline_fba = s.cache.load("fba_heuristic_baseline_v1")

    print("loading thermo_variants reports...")
    manifest_in = json.loads((THERMO_VARIANTS_DIR / "manifest.json").read_text())
    baseline_eq = parse_report(
        THERMO_VARIANTS_DIR / "baseline" / "Estimated_Reaction_Reversibility_Report_EQ.txt",
        drop_old_rev=False)
    baseline_map = {k: v["new_rev"] for k, v in baseline_eq.items()}

    # Per-variant: load report, diff vs baseline, attach FBA panel results
    print("computing per-variant diffs + assembling FBA payloads...")
    variant_payloads = {}
    changed_by_variant = {}
    summary_rows = []
    for v in manifest_in["variants"]:
        tag = v["tag"]
        if tag == "baseline":
            summary_rows.append({
                "tag": "baseline", "title": v["title"], "section": v["section"],
                "n_changed_vs_baseline": 0,
                "n_models_flip": 0, "n_models_flux_change": 0,
                "rev_counts": v["counts"]["EQ"],
            })
            continue
        eq = parse_report(THERMO_VARIANTS_DIR / tag /
                          "Estimated_Reaction_Reversibility_Report_EQ.txt",
                          drop_old_rev=False)
        vmap = {k: v_["new_rev"] for k, v_ in eq.items()}
        diff = variant_diff(baseline_map, vmap)
        changed_by_variant[tag] = diff["diffs"]

        # Variant FBA from cache (notebook 06 ran them)
        try:
            vfba = s.cache.load(f"fba_variant_{tag.replace('.', '_')}_v2")
        except KeyError:
            vfba = []

        # Panel grow-status / flux diff vs heuristic baseline
        odi = {r["model_id"]: r for r in heur_baseline_fba}
        rows = []
        n_grow = 0
        n_flux = 0
        for r in vfba:
            b = odi.get(r["model_id"], {})
            d = {
                "model_id": r["model_id"],
                "baseline_grows": bool(b.get("grows", False)),
                "variant_grows": bool(r.get("grows", False)),
                "baseline_flux": float(b.get("growth_flux", 0.0)),
                "variant_flux": float(r.get("growth_flux", 0.0)),
                "delta_flux": float(r.get("growth_flux", 0.0)) -
                              float(b.get("growth_flux", 0.0)),
                "n_overrides": int(r.get("n_overrides", 0)),
            }
            if d["baseline_grows"] != d["variant_grows"]:
                n_grow += 1
            if abs(d["delta_flux"]) > 1e-6:
                n_flux += 1
            rows.append(d)

        payload = {
            "tag": tag,
            "title": v["title"],
            "section": v["section"],
            "cfg": v["cfg"],
            "rev_counts": v["counts"]["EQ"],
            "n_changed": diff["n_changed"],
            "transitions": diff["transitions"],
            "diffs": diff["diffs"],
            "panel_fba": rows,
            "n_models_flip": n_grow,
            "n_models_flux_change": n_flux,
        }
        variant_payloads[tag] = payload
        summary_rows.append({
            "tag": tag, "title": v["title"], "section": v["section"],
            "n_changed_vs_baseline": diff["n_changed"],
            "n_models_flip": n_grow,
            "n_models_flux_change": n_flux,
            "rev_counts": v["counts"]["EQ"],
        })
        print(f"  {tag:>12}  changed: {diff['n_changed']:>5}  "
              f"models flip: {n_grow:>3}  flux change: {n_flux:>3}")

    # Reactions index, split into panel-full and others-compact for snappy load.
    print("building reactions index...")
    reactions_index = build_reactions_index(msdb_rxns, panel_rxn_sets, changed_by_variant)
    n_panel = sum(1 for r in reactions_index.values() if r["in_panel"])
    n_other = len(reactions_index) - n_panel
    print(f"  {len(reactions_index)} reactions in index "
          f"(panel union: {n_panel}, changed-only: {n_other})")

    reactions_panel = {rid: r for rid, r in reactions_index.items() if r["in_panel"]}

    # ---- attach analytic P(direction) per variant where available --------
    # Pulled from results/statistical_panel/p_direction__{tag}.csv (written by
    # scripts/run_statistical_panel.py).  We only embed in reactions_panel
    # because the website's reaction detail loads that file first.
    if STATS_DIR.exists():
        import csv
        for csv_path in sorted(STATS_DIR.glob("p_direction__*.csv")):
            tag = csv_path.stem.replace("p_direction__", "")
            with csv_path.open() as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    rid = row["rxn_id"]
                    if rid not in reactions_panel:
                        continue
                    reactions_panel[rid].setdefault("p_direction", {})[tag] = {
                        "p_forward":    float(row["p_forward"]),
                        "p_reverse":    float(row["p_reverse"]),
                        "p_reversible": float(row["p_reversible"]),
                    }
            print(f"  attached P(direction) for variant {tag!r} "
                  f"to {sum(1 for r in reactions_panel.values() if tag in r.get('p_direction', {}))} rxns")

    # ---- attach per-model flux distributions where available -------------
    flux_dist = {}
    if STATS_DIR.exists():
        import csv
        for csv_path in sorted(STATS_DIR.glob("panel_distribution__*.csv")):
            stem = csv_path.stem.replace("panel_distribution__", "")
            # stem looks like 'baseline__N50' -> split out tag and n_samples
            if "__N" in stem:
                tag, n_part = stem.split("__N", 1)
            else:
                tag, n_part = stem, "?"
            with csv_path.open() as fh:
                for row in csv.DictReader(fh):
                    mid = row["model_id"]
                    flux_dist.setdefault(mid, {})[tag] = {
                        "n_samples":           int(row["n_samples"]),
                        "mean_flux":           float(row["mean_flux"]),
                        "q05":                 float(row["q05"]),
                        "q50":                 float(row["q50"]),
                        "q95":                 float(row["q95"]),
                        "p_grows":             float(row["p_grows"]),
                        "point_estimate_flux": float(row["point_estimate_flux"]),
                        "point_in_ci95":       bool(int(row["point_in_ci95"])),
                    }
            print(f"  attached panel-flux distribution for variant {tag!r} "
                  f"(N={n_part}) for {len(flux_dist)} models")
    # Compact: keep just the searchable fields and the changed_by list.
    reactions_other = {
        rid: {
            "id": r["id"],
            "name": r["name"],
            "definition": r["definition"],
            "is_transport": r["is_transport"],
            "deltag": r["deltag"],
            "deltagerr": r["deltagerr"],
            "ec_numbers": r["ec_numbers"],
            "changed_by": r["changed_by"],
        }
        for rid, r in reactions_index.items() if not r["in_panel"]
    }

    # Panel info
    panel_info = build_panel_info(panel_ids, ondisk_fba, heur_baseline_fba, panel_rxn_sets)

    # ---- emit files
    print(f"writing JSON to {out_root}/ ...")
    with open(out_root / "manifest.json", "w") as fh:
        json.dump({
            "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "panel_size": len(panel_ids),
            "n_msdb_rxns": len(msdb_rxns),
            "n_panel_rxns": len(reactions_index),
            "variants": summary_rows,
        }, fh, indent=2, default=str)

    with open(out_root / "baseline.json", "w") as fh:
        json.dump({
            "tag": "baseline",
            "n_rxns": len(baseline_map),
            "by_rev": _count_field(baseline_map),
            "map": baseline_map,
        }, fh, separators=(",", ":"), default=str)

    for tag, payload in variant_payloads.items():
        with open(out_root / "variants" / f"{tag}.json", "w") as fh:
            json.dump(payload, fh, separators=(",", ":"), default=str)

    with open(out_root / "reactions_panel.json", "w") as fh:
        json.dump(reactions_panel, fh, separators=(",", ":"), default=str)
    with open(out_root / "reactions_other.json", "w") as fh:
        json.dump(reactions_other, fh, separators=(",", ":"), default=str)

    # Merge flux distributions into panel model rows.
    if flux_dist:
        for m in panel_info:
            d = flux_dist.get(m["model_id"])
            if d:
                m["flux_distribution"] = d
    with open(out_root / "panel.json", "w") as fh:
        json.dump({
            "models": panel_info,
        }, fh, indent=2, default=str)

    # Compact per-panel "rxnset" payload for the network/flux explorer
    with open(out_root / "panel_rxnsets.json", "w") as fh:
        json.dump({mid: sorted(panel_rxn_sets.get(mid, []))
                   for mid in panel_ids}, fh, separators=(",", ":"))

    print("done.")


def _count_field(m: dict) -> dict:
    out = {}
    for v in m.values():
        out[v] = out.get(v, 0) + 1
    return out


if __name__ == "__main__":
    main()
