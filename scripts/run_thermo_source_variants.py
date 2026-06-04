"""Run the four panel-FBA variants and emit long-form + summary CSVs.

Variants:
  - kbase_baseline : pure on-disk JSON bounds (reversibility_map=None)
  - gc             : per-source override using Group contribution operators
  - eq             : per-source override using eQuilibrator operators
  - dgp            : per-source override using dGPredictor operators

Outputs (always overwritten):
  - results/thermo_sources/panel_fba_long.csv
  - results/thermo_sources/panel_fba_summary.csv

Idempotent.  No CLI arguments -- importing or executing this file as
``python scripts/run_thermo_source_variants.py`` rebuilds both artifacts.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

# Make sibling scripts/ importable when run directly.
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import direction_pipeline as dp  # noqa: E402
import growth_heuristics as gh  # noqa: E402

ANALYSIS_DIR = Path("/scratch/ctaylor/core_models_analysis")
RESULTS_DIR = ANALYSIS_DIR / "results"
THERMO_DIR = RESULTS_DIR / "thermo_sources"
SELECTED_IDS_FILE = RESULTS_DIR / "selected_ids.txt"

LONG_CSV = THERMO_DIR / "panel_fba_long.csv"
SUMMARY_CSV = THERMO_DIR / "panel_fba_summary.csv"

SOURCE_SPECS = [
    ("gc", "Group contribution"),
    ("eq", "eQuilibrator"),
    ("dgp", "dGPredictor"),
]

# Canonical ordering for the variant axis -- shared between the CLI and the
# notebook so panel_fba_long.csv / panel_fba_summary.csv come out byte-identical
# regardless of producer.
VARIANT_ORDER = ["kbase_baseline", "gc", "eq", "dgp"]

# Human-readable label for each variant (mirrors what the notebook stashes in
# ``runs[variant]['label']``).
VARIANT_LABELS = {
    "kbase_baseline": "KBase baseline (on-disk bounds)",
    "gc":             "Group contribution",
    "eq":             "eQuilibrator",
    "dgp":            "dGPredictor",
}

# Canonical column orders -- imported by build_thermo_source_comparison_notebook
# so both producers emit the same header in the same order.
LONG_COLS = [
    "variant",
    "label",
    "model_id",
    "status",
    "growth_flux",
    "grows",
    "n_overrides",
]

SUMMARY_COLS = [
    "variant",
    "label",
    "n_panel",
    "n_growers",
    "frac_growers",
    "mean_flux",
    "median_flux",
    "n_status_optimal",
    "mean_overrides_per_model",
]

GROWTH_THRESHOLD = gh.GROWTH_THRESHOLD


def _ensure_source_snapshots() -> dict:
    """Make sure the per-source operator CSVs exist; snapshot if missing."""
    THERMO_DIR.mkdir(parents=True, exist_ok=True)
    paths: dict = {}
    for slug, source in SOURCE_SPECS:
        csv_path = THERMO_DIR / f"rxn_directions_{slug}.csv"
        json_path = THERMO_DIR / f"rxn_directions_{slug}.json"
        if not (csv_path.exists() and json_path.exists()):
            dp.snapshot_msdb_per_source(
                branch="origin/dev",
                source=source,
                out_csv=csv_path,
                out_json=json_path,
            )
        paths[slug] = csv_path
    return paths


def _model_overrides_expected(model_id: str, source_map: dict) -> int:
    """Count reactions in this model whose seed.reaction id is in source_map."""
    model_path = gh.MODELS_DIR / f"{model_id}.json"
    if not model_path.exists():
        return 0
    with model_path.open() as fh:
        model = json.load(fh)
    skip_prefixes = ("EX_", "SK_", "DM_", "bio")
    n = 0
    for rxn in model.get("reactions", []):
        rid = rxn.get("id") or ""
        if rid.startswith(skip_prefixes):
            continue
        anno = rxn.get("annotation") or {}
        seed = anno.get("seed.reaction")
        if isinstance(seed, list):
            seed = seed[0] if seed else None
        if seed and seed in source_map:
            n += 1
    return n


def _run_variants(ids: list) -> dict:
    """Run all four variants; return {variant: results_list}."""
    paths = _ensure_source_snapshots()

    results: dict = {}

    # 1. Baseline: on-disk bounds, no rewrite.
    print(f"[run] kbase_baseline ({len(ids)} models, 4 workers)")
    results["kbase_baseline"] = gh.run_panel(ids, reversibility_map=None, n_workers=4)

    # 2-4. Per-source variants: full rebind (baseline_map=None).
    for slug, source in SOURCE_SPECS:
        print(f"[run] {slug} from {paths[slug].name}")
        results[slug] = dp.panel_growth_with_source(
            ids,
            paths[slug],
            baseline_map=None,
            n_workers=4,
        )

    return results, paths


def _assert_overrides_match_coverage(
    variant_results: list,
    source_csv: Path,
    variant_label: str,
) -> list:
    """For each model, check n_overrides ~= count of model rxns whose seed is in source.

    Returns a list of (model_id, n_overrides, n_expected) for models where
    n_overrides == 0 (informative -- source has no opinion on any of the
    model's reactions).
    """
    source_map = dp.load_per_source_operators(source_csv)
    zero_override_models = []
    for r in variant_results:
        mid = r["model_id"]
        n_over = int(r.get("n_overrides", 0) or 0)
        n_exp = _model_overrides_expected(mid, source_map)
        # n_overrides should equal n_expected exactly (every covered rxn gets
        # rewritten when baseline_map=None and the source operator is non-empty).
        if n_over != n_exp:
            # Tolerate small mismatches but flag any nontrivial deviation.
            if abs(n_over - n_exp) > 1:
                print(
                    f"[warn] {variant_label} {mid}: n_overrides={n_over} "
                    f"!= expected={n_exp}"
                )
        if n_over == 0:
            zero_override_models.append((mid, n_over, n_exp))
            print(
                f"[info] {variant_label} {mid}: n_overrides=0 "
                f"(source covers none of this model's seed.reaction ids)"
            )
    return zero_override_models


def _build_long_df(results_by_variant: dict) -> pd.DataFrame:
    rows = []
    for variant, results in results_by_variant.items():
        label = VARIANT_LABELS.get(variant, variant)
        for r in results:
            status = r.get("status", "") or ""
            flux = float(r.get("growth_flux", 0.0) or 0.0)
            grows = bool(r.get("grows", False))
            rows.append({
                "variant": variant,
                "label": label,
                "model_id": r.get("model_id"),
                "status": status,
                "growth_flux": flux,
                "grows": grows,
                "n_overrides": int(r.get("n_overrides", 0) or 0),
            })
    df = pd.DataFrame(rows, columns=LONG_COLS)
    # Sort deterministically: variant in VARIANT_ORDER, then model_id ascending.
    # Anything outside VARIANT_ORDER falls back to alphabetical, after the known
    # variants.
    df = sort_long_df(df)
    return df


def sort_long_df(df: pd.DataFrame) -> pd.DataFrame:
    """Sort panel_fba_long rows by (variant in VARIANT_ORDER, model_id ascending).

    Exposed so the notebook (and any other consumer) can apply the same canonical
    ordering before writing the CSV.
    """
    order_map = {v: i for i, v in enumerate(VARIANT_ORDER)}
    fallback = len(VARIANT_ORDER)
    df = df.copy()
    df["__variant_order"] = df["variant"].map(lambda v: order_map.get(v, fallback))
    df = df.sort_values(
        ["__variant_order", "variant", "model_id"],
        kind="mergesort",
    ).drop(columns="__variant_order").reset_index(drop=True)
    return df


def _build_summary_df(long_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    # Preserve VARIANT_ORDER for the summary as well so the rows always line up.
    order_map = {v: i for i, v in enumerate(VARIANT_ORDER)}
    variants_present = sorted(
        long_df["variant"].unique(),
        key=lambda v: (order_map.get(v, len(VARIANT_ORDER)), v),
    )
    for variant in variants_present:
        sub = long_df[long_df["variant"] == variant]
        n_panel = int(len(sub))
        n_growers = int(sub["grows"].sum())
        rows.append({
            "variant": variant,
            "label": VARIANT_LABELS.get(variant, variant),
            "n_panel": n_panel,
            "n_growers": n_growers,
            "frac_growers": (n_growers / n_panel) if n_panel else 0.0,
            "mean_flux": float(sub["growth_flux"].mean()) if n_panel else 0.0,
            "median_flux": float(sub["growth_flux"].median()) if n_panel else 0.0,
            "n_status_optimal": int((sub["status"] == "optimal").sum()),
            "mean_overrides_per_model": float(sub["n_overrides"].mean()) if n_panel else 0.0,
        })
    return pd.DataFrame(rows, columns=SUMMARY_COLS)


def _compute_deltas(long_df: pd.DataFrame) -> dict:
    """For each non-baseline variant, compute deltas vs baseline."""
    base = long_df[long_df["variant"] == "kbase_baseline"].set_index("model_id")
    out: dict = {}
    for variant in long_df["variant"].unique():
        if variant == "kbase_baseline":
            continue
        var_sub = long_df[long_df["variant"] == variant].set_index("model_id")
        joined = base.join(var_sub, lsuffix="_base", rsuffix="_var", how="inner")
        delta = (joined["growth_flux_var"] - joined["growth_flux_base"]).abs()
        changed = int((joined["grows_base"] != joined["grows_var"]).sum())
        count_diff = int((delta > 1e-6).sum())
        out[variant] = {
            "models_changed_growth": changed,
            "mean_abs_delta_flux": float(delta.mean()),
            "max_abs_delta_flux": float(delta.max()),
            "count_flux_diff": count_diff,
        }
    return out


def main() -> dict:
    THERMO_DIR.mkdir(parents=True, exist_ok=True)
    ids = SELECTED_IDS_FILE.read_text().split()
    ids = [i for i in ids if i]
    print(f"[main] panel size = {len(ids)}")

    results_by_variant, paths = _run_variants(ids)

    # Coverage assertions per non-baseline variant.
    zero_override_summary: dict = {}
    for slug, _ in SOURCE_SPECS:
        zeros = _assert_overrides_match_coverage(
            results_by_variant[slug], paths[slug], slug
        )
        zero_override_summary[slug] = [m for (m, _, _) in zeros]
        # Also assert that at least one model in the panel has overrides > 0.
        nonzero = [
            r for r in results_by_variant[slug]
            if int(r.get("n_overrides", 0) or 0) > 0
        ]
        assert len(nonzero) > 0, (
            f"variant {slug}: no model had n_overrides > 0 -- "
            "source map appears empty or all panel models lack seed annos"
        )

    long_df = _build_long_df(results_by_variant)
    summary_df = _build_summary_df(long_df)
    deltas = _compute_deltas(long_df)

    # lineterminator='\n' is explicit so the file is identical regardless of
    # OS / pandas default (pandas writes LF by default but pinning it makes the
    # contract obvious to readers).
    long_df.to_csv(LONG_CSV, index=False, lineterminator="\n")
    summary_df.to_csv(SUMMARY_CSV, index=False, lineterminator="\n")

    print(f"[write] {LONG_CSV} ({len(long_df)} rows)")
    print(f"[write] {SUMMARY_CSV} ({len(summary_df)} rows)")

    return {
        "long_csv": str(LONG_CSV),
        "summary_csv": str(SUMMARY_CSV),
        "summary": summary_df.to_dict(orient="records"),
        "deltas_vs_baseline": deltas,
        "zero_override_models": zero_override_summary,
    }


if __name__ == "__main__":
    result = main()
    print(json.dumps(result, indent=2, default=str))
