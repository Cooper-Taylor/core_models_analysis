"""CLI wrapper -- per-source network-level analysis CSVs.

For each thermodynamics source (Group contribution / eQuilibrator / dGPredictor)
walks the 100 panel models and produces two artifacts under
``results/thermo_sources/``:

  - ``coverage_<slug>.csv``  -- per-model coverage stats (wide schema)
  - ``overrides_<slug>.csv`` -- per-model bound-class transitions (wide schema)

All heavy lifting lives in :mod:`direction_pipeline`: this file just
loads the panel ids, calls the helpers, writes the CSVs, and prints a
small summary dict.

ModelSEEDDatabase is *never* touched -- the only MSDB-derived input we
read is the pre-computed snapshot at
``results/rxn_directions_msdb_dev.csv`` (produced read-only via
``git show`` by ``direction_pipeline.snapshot_msdb``).
"""

from __future__ import annotations
import os

import json
import sys
from pathlib import Path
from typing import Iterable, Mapping

import pandas as pd

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR", "/scratch/ctaylor/core_models_analysis"))
SCRIPTS_DIR = ANALYSIS_DIR / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import direction_pipeline as dp  # noqa: E402

RESULTS_DIR = ANALYSIS_DIR / "results"
THERMO_DIR = RESULTS_DIR / "thermo_sources"
PANEL_IDS_PATH = RESULTS_DIR / "selected_ids.txt"
MSDB_SNAPSHOT_PATH = RESULTS_DIR / "rxn_directions_msdb_dev.csv"

# slug -> per-source CSV under thermo_sources/
SOURCES: dict = {
    "gc":  THERMO_DIR / "rxn_directions_gc.csv",
    "eq":  THERMO_DIR / "rxn_directions_eq.csv",
    "dgp": THERMO_DIR / "rxn_directions_dgp.csv",
}

COVERAGE_COLS = [
    "model_id",
    "n_reactions",
    "n_non_exchange",
    "n_with_seed_anno",
    "n_in_msdb",
    "n_covered_by_source",
    "frac_covered",
    "n_uncovered_by_source",
]


def _load_panel_ids(path: Path = PANEL_IDS_PATH) -> list:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def _write_csv(path: Path, cols: list, rows: Iterable[dict]) -> None:
    """Write rows to ``path`` as a CSV with LF line endings.

    Uses pandas DataFrame.to_csv (LF by default, pinned with ``lineterminator``
    for clarity) so output is byte-identical to what the notebook produces with
    ``df.to_csv(...)``.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(list(rows), columns=cols)
    df.to_csv(path, index=False, lineterminator="\n")


def _coverage_rows(panel_ids: list, source_csv: Path, msdb_map: dict) -> list:
    """Build wide-schema coverage rows by augmenting dp.source_coverage.

    dp.source_coverage already gives the canonical
    (n_reactions, n_with_seed_anno, n_covered_by_source, frac_covered) tuple
    used by the figures.  This helper adds the legacy
    (n_non_exchange, n_in_msdb, n_uncovered_by_source) columns by replaying
    the per-model walk for those counts only.
    """
    import growth_heuristics as gh  # local import to keep the module light

    cov = dp.source_coverage(panel_ids, source_csv)
    skip = ("EX_", "SK_", "DM_", "bio")
    rows = []
    for mid in panel_ids:
        base = cov[mid]
        # Recompute the two legacy columns (n_non_exchange, n_in_msdb).
        n_non_exchange = 0
        n_in_msdb = 0
        model_path = gh.MODELS_DIR / f"{mid}.json"
        if model_path.exists():
            with model_path.open() as fh:
                model = json.load(fh)
            for rxn in model.get("reactions", []):
                rid = rxn.get("id") or ""
                if rid.startswith(skip):
                    continue
                n_non_exchange += 1
                anno = rxn.get("annotation") or {}
                seed = anno.get("seed.reaction")
                if isinstance(seed, list):
                    seed = seed[0] if seed else None
                if seed and seed in msdb_map:
                    n_in_msdb += 1
        rows.append({
            "model_id":              mid,
            "n_reactions":           base["n_reactions"],
            "n_non_exchange":        n_non_exchange,
            "n_with_seed_anno":      base["n_with_seed_anno"],
            "n_in_msdb":             n_in_msdb,
            "n_covered_by_source":   base["n_covered_by_source"],
            "frac_covered":          base["frac_covered"],
            "n_uncovered_by_source": base["n_with_seed_anno"] - base["n_covered_by_source"],
        })
    return rows


def build_tables(
    panel_ids_path: Path = PANEL_IDS_PATH,
    sources: Mapping[str, Path] = SOURCES,
    msdb_snapshot_path: Path = MSDB_SNAPSHOT_PATH,
    out_dir: Path = THERMO_DIR,
) -> dict:
    """Build per-source coverage + overrides CSVs and return path / stats dicts."""
    panel_ids = _load_panel_ids(panel_ids_path)
    msdb_map = dp.load_msdb_reversibility_map(msdb_snapshot_path)

    coverage_paths: dict = {}
    overrides_paths: dict = {}
    stats: dict = {}

    for slug, source_csv in sources.items():
        cov_rows = _coverage_rows(panel_ids, source_csv, msdb_map)
        ovr_rows = dp.override_transitions(
            panel_ids, source_csv, msdb_snapshot_path, slug
        )

        cov_path = out_dir / f"coverage_{slug}.csv"
        ovr_path = out_dir / f"overrides_{slug}.csv"
        _write_csv(cov_path, COVERAGE_COLS, cov_rows)
        _write_csv(ovr_path, list(dp.OVERRIDE_TRANSITION_COLS), ovr_rows)
        coverage_paths[slug] = str(cov_path)
        overrides_paths[slug] = str(ovr_path)

        total_seen = sum(r["n_with_seed_anno"] for r in cov_rows)
        total_covered = sum(r["n_covered_by_source"] for r in cov_rows)
        frac_acc = sum(r["frac_covered"] for r in cov_rows)
        zero_cov = sum(1 for r in cov_rows if r["n_covered_by_source"] == 0)
        stats[slug] = {
            "total_panel_reactions_seen": int(total_seen),
            "total_covered":              int(total_covered),
            "mean_frac_covered":          (frac_acc / len(panel_ids)) if panel_ids else 0.0,
            "models_with_zero_coverage":  int(zero_cov),
        }

    return {
        "coverage_paths":            coverage_paths,
        "overrides_paths":           overrides_paths,
        "per_source_coverage_stats": stats,
    }


if __name__ == "__main__":
    out = build_tables()
    print(json.dumps(out, indent=2, sort_keys=True))
