#!/usr/bin/env python3
"""Build per-variant impact data across the FULL 5,683-model database.

The site's "Variant Browser" currently summarizes impact only against the
100-model descriptive growth panel.  This script extends that with the
analogous numbers across every core-models JSON under
``data/core_models_kegg2/``.

Outputs (under ``site/data/``):

  - ``all_models_rxnsets_summary.json``
      {model_id: n_seed_rxns}  -- one int per model (used to derive
      "models that contain at least one changed rxn" without shipping
      the full 5,683-model rxn lists to the browser).

  - ``all_models_baseline_fba.json``
      [{model_id, status, grows, growth_flux, n_overrides}, ...]
      heuristic-baseline (cascade default) rebound FBA for every model.

  - ``all_models_variants.json``
      Per-variant summary statistics:
        {tag: {
            n_models_containing_changed_rxn,
            n_rxn_instances_touched,
            max_increase: {model_id, delta_flux, baseline_flux, variant_flux},
            max_decrease: {model_id, ...},
            mean_flux, median_flux, std_flux,         -- across ALL 5,683 models
            n_models_grow,                            -- under variant
            n_grew_default_now_not, n_not_default_now_grew,
            n_models_flip,                            -- either direction
            n_flux_change,
        }}

  - ``all_models_variant_fba__{tag}.json``
      [{model_id, baseline_flux, variant_flux, delta_flux,
        baseline_grows, variant_grows, n_overrides}, ...]
      Only includes models whose reaction set intersected the variant's
      changed-reaction list.  For models with no intersection, the
      summary assumes variant_flux == baseline_flux (cascading
      reactions outside the model are inert).

Run order: this script does (rxnsets index) + (baseline FBA) once,
then iterates the variants from ``site/data/manifest.json``.  Use
``--skip-baseline`` to reuse a prior baseline pass.

Performance:
  * rxnset index: ~5 minutes (raw JSON parse, no cobra)
  * baseline FBA: ~20-40 minutes (depends on n_workers and CPU count)
  * each variant : ~5-30 minutes (depends on intersection size)

The script is incremental: existing per-variant FBA files are reused
unless ``--force-variant TAG`` or ``--force-all-variants`` is passed.
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import statistics
import sys
import time
from pathlib import Path
from typing import Optional

SCRIPTS = Path(__file__).resolve().parent
ANALYSIS_ROOT = SCRIPTS.parent
MODELS_DIR = ANALYSIS_ROOT / "data" / "core_models_kegg2"
SITE_DATA = ANALYSIS_ROOT / "site" / "data"
MSDB_ROOT = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase"))

sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(MSDB_ROOT / "Libs" / "Python"))


# ---------------------------------------------------------------------------
def list_model_ids() -> list:
    ids = []
    for p in sorted(MODELS_DIR.glob("*.json")):
        ids.append(p.stem)
    return ids


def _extract_rxnset(model_id: str) -> "tuple[str, list]":
    """Parse a model JSON and return its SEED reaction ID set.

    Raw JSON only — cobra-free for speed (the heuristic baseline FBA
    below loads via cobra anyway, but the rxnset lookup is the
    bottleneck if we use cobra for it too).
    """
    p = MODELS_DIR / f"{model_id}.json"
    d = json.loads(p.read_text())
    seeds = set()
    for r in d.get("reactions", []):
        anno = (r.get("annotation") or {}).get("seed.reaction")
        if anno:
            seeds.add(anno)
    return model_id, sorted(seeds)


def build_rxnsets(model_ids: list, n_workers: int) -> dict:
    """{model_id: [seed_rxn_id, ...]} for every model."""
    print(f"[rxnsets] indexing {len(model_ids)} models with {n_workers} workers...")
    t = time.time()
    if n_workers == 1:
        out = dict(_extract_rxnset(mid) for mid in model_ids)
    else:
        with mp.Pool(n_workers) as pool:
            out = dict(pool.imap_unordered(_extract_rxnset, model_ids, chunksize=64))
    print(f"[rxnsets] done in {time.time()-t:.1f}s "
          f"({sum(len(v) for v in out.values()):,} (model,rxn) pairs)")
    return out


# ---------------------------------------------------------------------------
def load_variant_maps(baseline_map: dict) -> dict:
    """Read each variant payload and apply its diff to the baseline map.

    Returns {tag: {rxn: rev}} including the synthetic 'baseline' entry.
    """
    out = {"baseline": dict(baseline_map)}
    for vfile in sorted((SITE_DATA / "variants").glob("*.json")):
        tag = vfile.stem
        payload = json.loads(vfile.read_text())
        vmap = dict(baseline_map)
        for d in payload["diffs"]:
            vmap[d["rxn"]] = d["new"]
        out[tag] = vmap
    return out


def variant_changed_rxns(baseline_map: dict, variant_map: dict) -> set:
    out = set()
    for r, b in baseline_map.items():
        v = variant_map.get(r)
        if v is not None and v != b:
            out.add(r)
    return out


# ---------------------------------------------------------------------------
def run_panel_fba_full(model_ids: list, reversibility_map: dict,
                       n_workers: int) -> list:
    """Wrapper around growth_heuristics.run_panel with a progress timer."""
    import growth_heuristics as gh
    t = time.time()
    print(f"[fba] running {len(model_ids)} models with {n_workers} workers...",
          flush=True)
    out = gh.run_panel(model_ids, reversibility_map=reversibility_map,
                       baseline_map=None, n_workers=n_workers)
    dt = time.time() - t
    n_ok = sum(1 for r in out if not r.get("error"))
    n_grow = sum(1 for r in out if r.get("grows"))
    print(f"[fba] done in {dt:.1f}s ({n_ok}/{len(out)} ok, {n_grow} growers)",
          flush=True)
    return out


# ---------------------------------------------------------------------------
def summarize_variant(tag: str, baseline_rows: list, variant_rows: list,
                      models_with_change: set, n_all_models: int,
                      n_rxn_instances_touched: int) -> dict:
    """Build the per-variant summary block from the FBA rows."""
    bbi = {r["model_id"]: r for r in baseline_rows}
    # For models not in variant_rows, variant_flux == baseline_flux.
    vbi = {r["model_id"]: r for r in variant_rows}

    # Build the comprehensive per-model variant flux series:
    deltas = []           # (model_id, baseline, variant, delta)
    flux_now = []
    n_grow = 0
    n_grew_then_not = 0   # was grower in baseline, now not under variant
    n_not_then_grew = 0   # vice versa
    n_flip = 0
    n_flux_change = 0
    for mid, b in bbi.items():
        b_flux = float(b.get("growth_flux", 0.0))
        b_grew = bool(b.get("grows"))
        v = vbi.get(mid)
        if v is not None:
            v_flux = float(v.get("growth_flux", 0.0))
            v_grew = bool(v.get("grows"))
        else:
            v_flux = b_flux
            v_grew = b_grew
        delta = v_flux - b_flux
        deltas.append({"model_id": mid, "baseline_flux": b_flux,
                       "variant_flux": v_flux, "delta_flux": delta,
                       "baseline_grows": b_grew, "variant_grows": v_grew})
        flux_now.append(v_flux)
        if v_grew:
            n_grow += 1
        if b_grew and not v_grew:
            n_grew_then_not += 1
        if (not b_grew) and v_grew:
            n_not_then_grew += 1
        if b_grew != v_grew:
            n_flip += 1
        if abs(delta) > 1e-6:
            n_flux_change += 1

    # Greatest increase / decrease
    deltas.sort(key=lambda d: d["delta_flux"])
    max_decrease = deltas[0] if deltas else None
    max_increase = deltas[-1] if deltas else None

    mean_flux = (sum(flux_now) / len(flux_now)) if flux_now else 0.0
    med_flux = statistics.median(flux_now) if flux_now else 0.0
    std_flux = statistics.pstdev(flux_now) if len(flux_now) > 1 else 0.0

    return {
        "tag": tag,
        "n_all_models": n_all_models,
        "n_models_containing_changed_rxn": len(models_with_change),
        "n_rxn_instances_touched": n_rxn_instances_touched,
        "max_increase": max_increase,
        "max_decrease": max_decrease,
        "mean_flux": mean_flux,
        "median_flux": med_flux,
        "std_flux": std_flux,
        "n_models_grow": n_grow,
        "n_grew_default_now_not": n_grew_then_not,
        "n_not_default_now_grew": n_not_then_grew,
        "n_models_flip": n_flip,
        "n_models_flux_change": n_flux_change,
    }


# ---------------------------------------------------------------------------
def main(argv: Optional[list] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--workers", type=int,
                    default=max(1, min(mp.cpu_count() - 1, 16)))
    ap.add_argument("--skip-rxnsets", action="store_true",
                    help="Reuse existing all_models_rxnsets_summary.json + the per-model rxnset cache")
    ap.add_argument("--skip-baseline", action="store_true",
                    help="Reuse existing all_models_baseline_fba.json")
    ap.add_argument("--variants", nargs="*", default=None,
                    help="Only process these variant tags (default: all)")
    ap.add_argument("--force-all-variants", action="store_true",
                    help="Recompute per-variant FBA even if a cached file exists")
    ap.add_argument("--smoke", type=int, default=0,
                    help="Smoke run: use only the first N models (default: all)")
    args = ap.parse_args(argv)

    SITE_DATA.mkdir(parents=True, exist_ok=True)
    rxnsets_path = SITE_DATA / "all_models_rxnsets.json"
    rxnsets_summary_path = SITE_DATA / "all_models_rxnsets_summary.json"
    baseline_fba_path = SITE_DATA / "all_models_baseline_fba.json"

    # ----- 1. rxnsets index ------------------------------------------------
    model_ids = list_model_ids()
    if args.smoke:
        model_ids = model_ids[:args.smoke]
        print(f"[smoke] truncated to {len(model_ids)} models", flush=True)

    if args.skip_rxnsets and rxnsets_path.exists():
        rxnsets = json.loads(rxnsets_path.read_text())
        # Filter to the model set we're actually running this pass.
        rxnsets = {mid: rxnsets[mid] for mid in model_ids if mid in rxnsets}
        print(f"[rxnsets] reused cache ({len(rxnsets)} models)")
    else:
        rxnsets = build_rxnsets(model_ids, args.workers)
        # Persist the full lists (one-time JSON, ~5,000 keys * ~100 strs)
        rxnsets_path.write_text(json.dumps(
            {mid: sorted(s) for mid, s in rxnsets.items()},
            separators=(",", ":"),
        ))
        # Persist a compact summary the browser actually loads.
        rxnsets_summary_path.write_text(json.dumps(
            {mid: len(s) for mid, s in rxnsets.items()},
            separators=(",", ":"),
        ))
        print(f"[rxnsets] wrote {rxnsets_path.name} + {rxnsets_summary_path.name}")

    # ----- 2. load baseline map + variant maps ----------------------------
    print("[maps] loading baseline + variant maps...", flush=True)
    baseline_payload = json.loads((SITE_DATA / "baseline.json").read_text())
    baseline_map = baseline_payload["map"]
    variant_maps = load_variant_maps(baseline_map)
    print(f"[maps] {len(variant_maps)} variants loaded")

    # ----- 3. baseline FBA across all models ------------------------------
    if args.skip_baseline and baseline_fba_path.exists():
        baseline_rows = json.loads(baseline_fba_path.read_text())
        # Restrict to our model_ids.
        baseline_rows = [r for r in baseline_rows if r["model_id"] in set(model_ids)]
        print(f"[baseline] reused {baseline_fba_path.name} ({len(baseline_rows)} rows)")
    else:
        baseline_rows = run_panel_fba_full(model_ids, baseline_map, args.workers)
        baseline_fba_path.write_text(json.dumps(baseline_rows, separators=(",", ":"),
                                                default=str))
        print(f"[baseline] wrote {baseline_fba_path.name}")

    # ----- 4. per-variant FBA ---------------------------------------------
    summary = {}
    tags = args.variants or [t for t in variant_maps if t != "baseline"]
    for tag in tags:
        if tag == "baseline":
            continue
        vmap = variant_maps.get(tag)
        if vmap is None:
            print(f"[skip] unknown variant: {tag}")
            continue

        changed_rxns = variant_changed_rxns(baseline_map, vmap)
        # Models whose reaction set intersects the variant's changed rxns
        models_with_change = []
        n_rxn_instances_touched = 0
        ch_set = set(changed_rxns)
        for mid, rxs in rxnsets.items():
            inter = ch_set & set(rxs)
            if inter:
                models_with_change.append(mid)
                n_rxn_instances_touched += len(inter)

        per_variant_path = SITE_DATA / f"all_models_variant_fba__{tag}.json"
        if per_variant_path.exists() and not args.force_all_variants:
            variant_rows = json.loads(per_variant_path.read_text())
            variant_rows = [r for r in variant_rows if r["model_id"] in set(model_ids)]
            print(f"[variant {tag}] reused {per_variant_path.name} ({len(variant_rows)} rows)")
        else:
            if not models_with_change:
                variant_rows = []
                print(f"[variant {tag}] no models contain any changed rxn; skipping FBA")
            else:
                print(f"[variant {tag}] {len(changed_rxns)} changed rxns; "
                      f"{len(models_with_change)} models affected; running FBA...", flush=True)
                variant_rows = run_panel_fba_full(
                    models_with_change, vmap, args.workers)
            # Enrich each row with baseline / variant / delta fields so the
            # UI does not have to join against all_models_baseline_fba.json
            # client-side. Keeps the cobra-output fields too.
            bbi = {r["model_id"]: r for r in baseline_rows}
            enriched = []
            for r in variant_rows:
                b = bbi.get(r["model_id"], {})
                b_flux = float(b.get("growth_flux", 0.0))
                v_flux = float(r.get("growth_flux", 0.0))
                enriched.append({
                    **r,
                    "baseline_flux": b_flux,
                    "variant_flux": v_flux,
                    "delta_flux": v_flux - b_flux,
                    "baseline_grows": bool(b.get("grows", False)),
                    "variant_grows": bool(r.get("grows", False)),
                })
            per_variant_path.write_text(json.dumps(enriched, separators=(",", ":"),
                                                   default=str))
            variant_rows = enriched
            print(f"[variant {tag}] wrote {per_variant_path.name} (enriched)")

        block = summarize_variant(
            tag, baseline_rows, variant_rows, set(models_with_change),
            len(model_ids), n_rxn_instances_touched,
        )
        summary[tag] = block
        print(f"[variant {tag}] mean_flux={block['mean_flux']:.4f}  "
              f"flip={block['n_models_flip']}  "
              f"flux_change={block['n_models_flux_change']}  "
              f"max_inc={block['max_increase']['delta_flux']:+.4f} "
              f"({block['max_increase']['model_id']})  "
              f"max_dec={block['max_decrease']['delta_flux']:+.4f} "
              f"({block['max_decrease']['model_id']})", flush=True)

    # ----- 5. write summary -----------------------------------------------
    out_path = SITE_DATA / "all_models_variants.json"
    # Preserve existing variants that weren't recomputed this pass.
    existing = {}
    if out_path.exists():
        try:
            existing_doc = json.loads(out_path.read_text())
            # File format is {generated_at, n_all_models, variants: {...}};
            # we only want to carry forward the inner per-variant blocks.
            existing = existing_doc.get("variants", {}) or {}
        except json.JSONDecodeError:
            pass
    # Add the synthetic baseline summary so the UI can show "all_models" for it too.
    b_flux = [float(r.get("growth_flux", 0.0)) for r in baseline_rows]
    baseline_block = {
        "tag": "baseline",
        "n_all_models": len(model_ids),
        "n_models_containing_changed_rxn": 0,
        "n_rxn_instances_touched": 0,
        "max_increase": None,
        "max_decrease": None,
        "mean_flux": (sum(b_flux) / len(b_flux)) if b_flux else 0.0,
        "median_flux": statistics.median(b_flux) if b_flux else 0.0,
        "std_flux": statistics.pstdev(b_flux) if len(b_flux) > 1 else 0.0,
        "n_models_grow": sum(1 for r in baseline_rows if r.get("grows")),
        "n_grew_default_now_not": 0,
        "n_not_default_now_grew": 0,
        "n_models_flip": 0,
        "n_models_flux_change": 0,
    }
    merged = dict(existing)
    merged["baseline"] = baseline_block
    merged.update(summary)
    out_path.write_text(json.dumps({
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "n_all_models": len(model_ids),
        "variants": merged,
    }, indent=2, default=str))
    print(f"[summary] wrote {out_path}")


if __name__ == "__main__":
    main()
