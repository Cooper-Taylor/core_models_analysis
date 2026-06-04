"""Append the LIVE cascade reaction-direction variant to the panel FBA tables.

Idempotent: re-runs drop any prior 'cascade_live' rows and re-append.

Writes (in place):
  - results/variant_panel_fba.csv      (long form: variant, model_id, status, growth_flux, grows, n_overrides)
  - results/variant_panel_summary.csv  (per-variant: n_panel, n_growers, frac_growers, mean_flux, median_flux, n_status_optimal)

Mirrors run_variant_panel.py style. No CLI args.
"""
from __future__ import annotations
import os

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR", "/scratch/ctaylor/core_models_analysis")) / "scripts"))

import pandas as pd

import direction_pipeline as dp  # noqa: F401  (parity with sibling driver)
import growth_heuristics as gh

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR", "/scratch/ctaylor/core_models_analysis"))
RESULTS_DIR = ANALYSIS_DIR / "results"
N_WORKERS = 4

# Sort order for summary: keep cascade_live next to msdb_dev for readability.
SUMMARY_ORDER = [
    "on_disk",
    "msdb_dev",
    "cascade_live",
    "msdb_claude",
    "all_reversible",
    "all_forward",
    "flip_eq_to_gt",
    "branch_diff_only",
]


def _order_key(variant: str) -> tuple:
    try:
        return (0, SUMMARY_ORDER.index(variant))
    except ValueError:
        return (1, variant)


def run_cascade_live(ids: list) -> list:
    cascade_map_path = RESULTS_DIR / "rxn_directions_cascade_live.json"
    cascade_map = json.loads(cascade_map_path.read_text())
    print(f"[load] cascade_live map size: {len(cascade_map)}", flush=True)

    t0 = time.time()
    print(f"[run] variant=cascade_live n_models={len(ids)} workers={N_WORKERS} ...", flush=True)
    # baseline_map=None forces full rebind so cascade_map is authoritative.
    res = gh.run_panel(
        ids,
        reversibility_map=cascade_map,
        baseline_map=None,
        n_workers=N_WORKERS,
    )
    n_grow = sum(1 for r in res if r.get("grows"))
    mean_flux = sum(float(r.get("growth_flux", 0.0) or 0.0) for r in res) / max(1, len(res))
    print(
        f"  done variant=cascade_live models={len(res)} growers={n_grow} "
        f"mean_flux={mean_flux:.4g} elapsed={time.time() - t0:.1f}s",
        flush=True,
    )
    return res


def append_long_csv(results: list, out_path: Path) -> None:
    df = pd.read_csv(out_path)
    df = df[df["variant"] != "cascade_live"].copy()
    new_rows = pd.DataFrame(
        [
            {
                "variant": "cascade_live",
                "model_id": r.get("model_id", ""),
                "status": r.get("status", ""),
                "growth_flux": float(r.get("growth_flux", 0.0) or 0.0),
                "grows": bool(r.get("grows", False)),
                "n_overrides": int(r.get("n_overrides", 0) or 0),
            }
            for r in results
        ]
    )
    out = pd.concat([df, new_rows], ignore_index=True)
    out.to_csv(out_path, index=False)


def append_summary_csv(results: list, out_path: Path) -> None:
    df = pd.read_csv(out_path)
    df = df[df["variant"] != "cascade_live"].copy()
    n = len(results)
    n_grow = sum(1 for r in results if r.get("grows"))
    fluxes = [float(r.get("growth_flux", 0.0) or 0.0) for r in results]
    n_optimal = sum(1 for r in results if r.get("status") == "optimal")
    new_row = pd.DataFrame(
        [
            {
                "variant": "cascade_live",
                "n_panel": n,
                "n_growers": n_grow,
                "frac_growers": (n_grow / n) if n else 0.0,
                "mean_flux": (sum(fluxes) / n) if n else 0.0,
                "median_flux": float(pd.Series(fluxes).median()) if fluxes else 0.0,
                "n_status_optimal": n_optimal,
            }
        ]
    )
    out = pd.concat([df, new_row], ignore_index=True)
    out["__order__"] = out["variant"].map(_order_key)
    out = out.sort_values("__order__", kind="stable").drop(columns="__order__").reset_index(drop=True)
    out.to_csv(out_path, index=False)


def compute_delta_vs_msdb_dev(long_path: Path) -> dict:
    df = pd.read_csv(long_path)
    cas = df[df["variant"] == "cascade_live"].set_index("model_id")
    dev = df[df["variant"] == "msdb_dev"].set_index("model_id")
    common = cas.index.intersection(dev.index)
    cas = cas.loc[common]
    dev = dev.loc[common]

    # Normalize grows -> bool (CSV gives strings)
    def _to_bool(x):
        if isinstance(x, bool):
            return x
        return str(x).strip().lower() == "true"

    cas_grow = cas["grows"].map(_to_bool)
    dev_grow = dev["grows"].map(_to_bool)
    models_changed_growth = int((cas_grow != dev_grow).sum())

    delta = (cas["growth_flux"].astype(float) - dev["growth_flux"].astype(float)).abs()
    mean_abs_delta = float(delta.mean()) if len(delta) else 0.0
    max_abs_delta = float(delta.max()) if len(delta) else 0.0
    count_diff = int((delta > 1e-6).sum())

    # Direction of growth changes (cascade vs dev)
    cascade_gained = int(((~dev_grow) & cas_grow).sum())
    cascade_lost = int((dev_grow & (~cas_grow)).sum())

    return {
        "models_changed_growth": models_changed_growth,
        "cascade_gained_growth_vs_dev": cascade_gained,
        "cascade_lost_growth_vs_dev": cascade_lost,
        "mean_abs_delta": mean_abs_delta,
        "max_abs_delta": max_abs_delta,
        "count_diff": count_diff,
        "n_models_compared": int(len(common)),
    }


def main() -> dict:
    ids = (RESULTS_DIR / "selected_ids.txt").read_text().split()
    print(f"[init] panel size: {len(ids)}", flush=True)

    long_path = RESULTS_DIR / "variant_panel_fba.csv"
    summary_path = RESULTS_DIR / "variant_panel_summary.csv"

    results = run_cascade_live(ids)
    append_long_csv(results, long_path)
    append_summary_csv(results, summary_path)

    delta = compute_delta_vs_msdb_dev(long_path)

    summary_df = pd.read_csv(summary_path)
    grow_counts = {row["variant"]: int(row["n_growers"]) for _, row in summary_df.iterrows()}
    mean_flux = {row["variant"]: float(row["mean_flux"]) for _, row in summary_df.iterrows()}

    payload = {
        "variant_csv_path": str(long_path),
        "summary_csv_path": str(summary_path),
        "grow_counts": grow_counts,
        "mean_flux": mean_flux,
        "delta_vs_msdb_dev": delta,
    }
    print(json.dumps(payload, indent=2))
    return payload


if __name__ == "__main__":
    main()
