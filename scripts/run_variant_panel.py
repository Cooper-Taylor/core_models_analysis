"""Execute panel FBA for SEVERAL reaction-direction variants.

Writes:
  - results/variant_panel_fba.csv      (long form: variant, model_id, status, growth_flux, grows, n_overrides)
  - results/variant_panel_summary.csv  (per-variant: n_panel, n_growers, frac_growers, mean_flux, median_flux, n_status_optimal)
  - results/variant_panel_summary.json (compact JSON for the harness return)
"""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, "/scratch/ctaylor/core_models_analysis/scripts")

import pandas as pd

import direction_pipeline as dp
import growth_heuristics as gh

ANALYSIS_DIR = Path("/scratch/ctaylor/core_models_analysis")
RESULTS_DIR = ANALYSIS_DIR / "results"
N_WORKERS = 4

VARIANT_ORDER = [
    "on_disk",
    "msdb_dev",
    "msdb_claude",
    "all_reversible",
    "all_forward",
    "flip_eq_to_gt",
    "branch_diff_only",
]


def _ensure_msdb_direction_csv(branch: str) -> Path:
    """Make sure rxn_directions_msdb_<branch>.csv exists; snapshot via git if not."""
    p = RESULTS_DIR / f"rxn_directions_msdb_{branch}.csv"
    if p.exists():
        return p
    print(f"  snapshotting MSDB branch '{branch}' to {p.name} ...", flush=True)
    dp.snapshot_msdb(branch=branch, out_path=p)
    return p


def _load_or_snapshot(branch: str) -> dict:
    """Prefer existing rev_map_<branch>.json (fast); fall back to MSDB snapshot."""
    json_p = RESULTS_DIR / f"rev_map_{branch}.json"
    if json_p.exists():
        return dp.load_directions_from_path(json_p)
    csv_p = _ensure_msdb_direction_csv(branch)
    return dp.load_directions_from_path(csv_p)


def _ensure_csv_from_map(branch: str, dmap: dict) -> Path:
    """Persist a JSON-backed map to rxn_directions_msdb_<branch>.csv if missing."""
    csv_p = RESULTS_DIR / f"rxn_directions_msdb_{branch}.csv"
    if not csv_p.exists():
        dp.save_directions_to_path(dmap, csv_p)
    return csv_p


def run_all(ids: list) -> dict:
    """Run all variants. Returns {variant_name: results_list}."""
    failed: list = []
    results: dict = {}

    # Prepare direction maps once.
    print("[load] msdb_dev direction map ...", flush=True)
    dev_map = _load_or_snapshot("dev")
    _ensure_csv_from_map("dev", dev_map)
    print(f"  dev map size: {len(dev_map)}", flush=True)

    print("[load] msdb_claude direction map ...", flush=True)
    claude_map = _load_or_snapshot("claude")
    _ensure_csv_from_map("claude", claude_map)
    print(f"  claude map size: {len(claude_map)}", flush=True)

    # Derived maps.
    all_forward_map = {rid: ">" for rid in dev_map}
    flip_eq_to_gt_map = {
        rid: (">" if rev == "=" else rev) for rid, rev in dev_map.items()
    }
    # branch_diff_only: start from dev; for rxns that DIFFER in claude, use claude's value.
    branch_diff_only_map = dict(dev_map)
    for rid, claude_rev in claude_map.items():
        if dev_map.get(rid) != claude_rev:
            branch_diff_only_map[rid] = claude_rev
    n_branch_diff = sum(
        1 for rid in branch_diff_only_map if branch_diff_only_map[rid] != dev_map.get(rid)
    )
    print(f"  branch_diff_only differences applied: {n_branch_diff}", flush=True)

    plans = [
        ("on_disk", dict(reversibility_map=None)),
        # baseline_map=None to force full rebind from the supplied map
        ("msdb_dev", dict(reversibility_map=dev_map, baseline_map=None)),
        ("msdb_claude", dict(reversibility_map=claude_map, baseline_map=None)),
        ("all_reversible", dict(reversibility_map=None, ignore_bounds=True)),
        ("all_forward", dict(reversibility_map=all_forward_map, baseline_map=None)),
        ("flip_eq_to_gt", dict(reversibility_map=flip_eq_to_gt_map, baseline_map=None)),
        ("branch_diff_only", dict(reversibility_map=branch_diff_only_map, baseline_map=None)),
    ]

    for name, kwargs in plans:
        t0 = time.time()
        print(f"[run] variant={name} ...", flush=True)
        try:
            res = gh.run_panel(ids, n_workers=N_WORKERS, **kwargs)
            results[name] = res
            n_grow = sum(1 for r in res if r.get("grows"))
            mean_flux = (
                sum(float(r.get("growth_flux", 0.0)) for r in res) / max(1, len(res))
            )
            dt = time.time() - t0
            print(
                f"  done variant={name} models={len(res)} growers={n_grow} "
                f"mean_flux={mean_flux:.4g} elapsed={dt:.1f}s",
                flush=True,
            )
        except Exception as e:
            print(f"  FAILED variant={name}: {type(e).__name__}: {e}", flush=True)
            failed.append(name)
    return results, failed


def write_long_csv(results: dict, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["variant", "model_id", "status", "growth_flux", "grows", "n_overrides"])
        for variant in VARIANT_ORDER:
            res = results.get(variant)
            if not res:
                continue
            for r in res:
                writer.writerow([
                    variant,
                    r.get("model_id", ""),
                    r.get("status", ""),
                    float(r.get("growth_flux", 0.0) or 0.0),
                    bool(r.get("grows", False)),
                    int(r.get("n_overrides", 0) or 0),
                ])


def write_summary_csv(results: dict, out_path: Path) -> dict:
    rows = []
    for variant in VARIANT_ORDER:
        res = results.get(variant)
        if not res:
            continue
        n = len(res)
        n_grow = sum(1 for r in res if r.get("grows"))
        fluxes = [float(r.get("growth_flux", 0.0) or 0.0) for r in res]
        n_optimal = sum(1 for r in res if r.get("status") == "optimal")
        mean_flux = sum(fluxes) / n if n else 0.0
        median_flux = float(pd.Series(fluxes).median()) if fluxes else 0.0
        rows.append({
            "variant": variant,
            "n_panel": n,
            "n_growers": n_grow,
            "frac_growers": (n_grow / n) if n else 0.0,
            "mean_flux": mean_flux,
            "median_flux": median_flux,
            "n_status_optimal": n_optimal,
        })
    pd.DataFrame(rows).to_csv(out_path, index=False)
    return {
        "grow_counts": {row["variant"]: row["n_growers"] for row in rows},
        "mean_flux": {row["variant"]: row["mean_flux"] for row in rows},
        "variants_run": [row["variant"] for row in rows],
    }


def main() -> None:
    ids = (RESULTS_DIR / "selected_ids.txt").read_text().split()
    print(f"[init] panel size: {len(ids)}", flush=True)
    results, failed = run_all(ids)

    long_path = RESULTS_DIR / "variant_panel_fba.csv"
    summary_path = RESULTS_DIR / "variant_panel_summary.csv"
    summary_json = RESULTS_DIR / "variant_panel_summary.json"

    write_long_csv(results, long_path)
    aux = write_summary_csv(results, summary_path)

    payload = {
        "variants": aux["variants_run"],
        "failed_variants": failed,
        "panel_csv_path": str(long_path),
        "summary_csv_path": str(summary_path),
        "grow_counts": aux["grow_counts"],
        "mean_flux": aux["mean_flux"],
    }
    summary_json.write_text(json.dumps(payload, indent=2))
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
