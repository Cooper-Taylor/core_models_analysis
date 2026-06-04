"""
Helpers for the reaction-direction pipeline notebook (09).

This module sits *next to* ``growth_heuristics`` and lets the notebook:

  - load a ``{rxn_id: reversibility}`` map from a TSV / CSV / JSON file
  - persist such a map back to CSV (the canonical on-disk form for this
    pipeline)
  - snapshot a chosen MSDB branch's reaction shards into a CSV under
    ``results/`` **without ever touching the MSDB working tree**
  - run the panel FBA with the supplied direction map and return both
    the raw results and a small DataFrame-shaped summary
  - line several named runs up for plotting (a long-form DataFrame keyed
    by variant name)

Everything reads MSDB via ``git show <branch>:<shard>`` so the on-disk
working tree of ``/scratch/ctaylor/ModelSEEDDatabase`` is never mutated.
All output paths live under
``/scratch/ctaylor/core_models_analysis/results/``.
"""

from __future__ import annotations
import os

import csv
import io
import json
import subprocess
from pathlib import Path
from typing import Iterable, Mapping, Optional

import pandas as pd

import growth_heuristics as gh

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR", "/scratch/ctaylor/core_models_analysis"))
RESULTS_DIR = ANALYSIS_DIR / "results"
PROJECT_ROOT = ANALYSIS_DIR  # alias used by per-source helpers' path guards
MSDB_DEFAULT_ROOT = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase"))
MSDB_SHARD_FMT = "Biochemistry/reaction_{:02d}.tsv"
MSDB_JSON_SHARD_FMT = "Biochemistry/reaction_{:02d}.json"
MSDB_N_SHARDS = 61

VALID_DIRECTIONS = {">", "<", "=", "?"}

# Sentinel energy ModelSEED writes when no estimate is available.
_THERMO_SENTINEL = 10000000

# Canonical per-source labels (the keys used inside the JSON
# ``thermodynamics`` dict in PR #263's new format).
PER_SOURCE_LABELS = ("Group contribution", "eQuilibrator", "dGPredictor")

# Lower-cased, dash-joined slugs used for output filenames.
_SOURCE_SLUGS = {
    "Group contribution": "group-contribution",
    "eQuilibrator": "equilibrator",
    "dGPredictor": "dgpredictor",
}


def _source_slug(source: str) -> str:
    """Return the canonical filename slug for a per-source label."""
    if source in _SOURCE_SLUGS:
        return _SOURCE_SLUGS[source]
    return source.strip().lower().replace(" ", "-")


# ---------------------------------------------------------------------------
# Direction map I/O
# ---------------------------------------------------------------------------
def _normalize_direction(value: str) -> str:
    """Coerce common spellings of a direction marker to ``>``/``<``/``=``/``?``."""
    if value is None:
        return "?"
    v = str(value).strip()
    if v in VALID_DIRECTIONS:
        return v
    low = v.lower()
    mapping = {
        "forward": ">",
        "fwd": ">",
        "f": ">",
        "reverse": "<",
        "rev": "<",
        "r": "<",
        "reversible": "=",
        "both": "=",
        "bidirectional": "=",
        "unknown": "?",
        "": "?",
        "none": "?",
        "null": "?",
    }
    if low in mapping:
        return mapping[low]
    # Last-chance: take the first character if it's already one of ours.
    if v and v[0] in VALID_DIRECTIONS:
        return v[0]
    return "?"


def _row_to_pair(row: Mapping[str, str]) -> Optional[tuple]:
    """Pull (rxn_id, direction) out of a CSV/TSV row using flexible keys."""
    rid = (
        row.get("rxn_id")
        or row.get("id")
        or row.get("reaction_id")
        or row.get("reaction")
    )
    if not rid:
        return None
    raw = (
        row.get("reversibility")
        if row.get("reversibility") not in (None, "")
        else row.get("direction")
    )
    if raw is None:
        raw = row.get("rev")
    return (str(rid).strip(), _normalize_direction(raw))


def load_directions_from_path(path: Path) -> dict:
    """Load a ``{rxn_id: reversibility}`` map from CSV / TSV / JSON.

    JSON shape can be either:
      - ``{rxn_id: direction, ...}``  (direct map)
      - ``[{"rxn_id": "...", "reversibility": ">"}, ...]``  (records)

    CSV/TSV files must have a header row.  Recognised id columns:
    ``rxn_id``, ``id``, ``reaction_id``, ``reaction``.  Recognised
    direction columns: ``reversibility``, ``direction``, ``rev``.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Direction source not found: {p}")
    suffix = p.suffix.lower()

    if suffix == ".json":
        data = json.loads(p.read_text())
        if isinstance(data, dict):
            return {str(k): _normalize_direction(v) for k, v in data.items()}
        if isinstance(data, list):
            out = {}
            for row in data:
                pair = _row_to_pair(row) if isinstance(row, dict) else None
                if pair is not None:
                    out[pair[0]] = pair[1]
            return out
        raise ValueError(f"Unexpected JSON layout in {p}")

    delim = "\t" if suffix in (".tsv", ".tab") else ","
    out = {}
    with p.open() as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        for row in reader:
            pair = _row_to_pair(row)
            if pair is not None:
                out[pair[0]] = pair[1]
    return out


def save_directions_to_path(
    direction_map: Mapping[str, str], path: Path
) -> Path:
    """Persist ``direction_map`` as a two-column CSV (``rxn_id,reversibility``)."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["rxn_id", "reversibility"])
        for rid in sorted(direction_map):
            writer.writerow([rid, direction_map[rid]])
    return p


# ---------------------------------------------------------------------------
# MSDB snapshot (read-only via git)
# ---------------------------------------------------------------------------
def _git_show(branch: str, path: str, repo: Path) -> str:
    """Return the contents of ``<branch>:<path>`` from ``repo`` via git show."""
    result = subprocess.run(
        ["git", "show", f"{branch}:{path}"],
        cwd=str(repo),
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout


def snapshot_msdb(
    branch: str = "dev",
    out_path: Optional[Path] = None,
    repo: Path = MSDB_DEFAULT_ROOT,
    n_shards: int = MSDB_N_SHARDS,
    extra_columns: Iterable[str] = ("deltag", "deltagerr", "status", "is_transport"),
) -> Path:
    """Snapshot the reaction-direction column of ``branch`` to CSV.

    Reads each of the 61 ``Biochemistry/reaction_NN.tsv`` shards via
    ``git show`` -- the working tree is never touched -- and writes a CSV
    of ``rxn_id, reversibility[, deltag, deltagerr, status, is_transport]``
    under ``results/``.

    Pass ``out_path=None`` to default to
    ``results/rxn_directions_msdb_<branch>.csv``.
    """
    if out_path is None:
        safe_branch = branch.replace("/", "_")
        out_path = RESULTS_DIR / f"rxn_directions_msdb_{safe_branch}.csv"
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    cols = ["rxn_id", "reversibility", *extra_columns]
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols)
        writer.writeheader()
        for i in range(n_shards):
            shard = MSDB_SHARD_FMT.format(i)
            text = _git_show(branch, shard, repo)
            reader = csv.DictReader(io.StringIO(text), delimiter="\t")
            for row in reader:
                rid = row.get("id")
                if not rid:
                    continue
                out_row = {
                    "rxn_id": rid,
                    "reversibility": _normalize_direction(row.get("reversibility")),
                }
                for col in extra_columns:
                    out_row[col] = row.get(col, "")
                writer.writerow(out_row)
    return out_path


# ---------------------------------------------------------------------------
# Panel growth wrapper
# ---------------------------------------------------------------------------
def panel_growth(
    model_ids: Iterable[str],
    direction_map: Mapping[str, str],
    baseline_map: Optional[Mapping[str, str]] = None,
    n_workers: Optional[int] = None,
) -> dict:
    """Run the FBA panel using ``direction_map`` and return results + summary.

    Returns ``{'results': [...], 'summary': DataFrame, 'totals': dict}``
    where ``totals`` collapses the panel to ``{n_models, n_grow,
    mean_flux, median_flux, max_flux}``.
    """
    ids = list(model_ids)
    results = gh.run_panel(
        ids,
        reversibility_map=dict(direction_map),
        baseline_map=dict(baseline_map) if baseline_map is not None else None,
        n_workers=n_workers,
    )
    summary_df = pd.DataFrame(results)
    if not summary_df.empty:
        totals = {
            "n_models": int(len(summary_df)),
            "n_grow": int(summary_df["grows"].sum()),
            "mean_flux": float(summary_df["growth_flux"].mean()),
            "median_flux": float(summary_df["growth_flux"].median()),
            "max_flux": float(summary_df["growth_flux"].max()),
        }
    else:
        totals = {"n_models": 0, "n_grow": 0, "mean_flux": 0.0,
                  "median_flux": 0.0, "max_flux": 0.0}
    return {"results": results, "summary": summary_df, "totals": totals}


# ---------------------------------------------------------------------------
# Cross-run helpers
# ---------------------------------------------------------------------------
def compare_runs(name_to_results: Mapping[str, Iterable[dict]]) -> pd.DataFrame:
    """Pivot several ``panel_growth['results']`` lists into long-form.

    Output columns: ``variant``, ``model_id``, ``grows``, ``growth_flux``,
    ``status``, ``n_overrides``.  Suitable for groupby/plotting.
    """
    rows = []
    for name, results in name_to_results.items():
        for r in results:
            rows.append({
                "variant": name,
                "model_id": r.get("model_id"),
                "grows": bool(r.get("grows")),
                "growth_flux": float(r.get("growth_flux", 0.0) or 0.0),
                "status": r.get("status", ""),
                "n_overrides": int(r.get("n_overrides", 0) or 0),
            })
    return pd.DataFrame(rows)


def variant_totals(name_to_results: Mapping[str, Iterable[dict]]) -> pd.DataFrame:
    """Collapse :func:`compare_runs` output to one row per variant."""
    df = compare_runs(name_to_results)
    if df.empty:
        return df
    agg = (
        df.groupby("variant", sort=False)
        .agg(
            n_models=("model_id", "count"),
            n_grow=("grows", "sum"),
            mean_flux=("growth_flux", "mean"),
            median_flux=("growth_flux", "median"),
            max_flux=("growth_flux", "max"),
        )
        .reset_index()
    )
    agg["n_grow"] = agg["n_grow"].astype(int)
    return agg


# ---------------------------------------------------------------------------
# Live cascade re-run (reversibility_lib + BiochemPy)
# ---------------------------------------------------------------------------
# Where BiochemPy lives in the upstream MSDB checkout.  This is on the import
# path only when ``run_cascade_live`` is invoked, so the rest of the module
# (and its existing call sites) never pay the load cost.
_MSDB_LIBS_PATH = MSDB_DEFAULT_ROOT / "Libs" / "Python"


def _assert_under_results(path: Path, label: str) -> None:
    """Refuse paths that would escape ``results/`` -- MSDB stays read-only."""
    try:
        Path(path).resolve().relative_to(RESULTS_DIR.resolve())
    except ValueError as exc:
        raise ValueError(
            f"{label} must live under {RESULTS_DIR}; got {path}"
        ) from exc


def run_cascade_live(
    out_csv: Optional[Path] = None,
    out_json: Optional[Path] = None,
    cfg=None,
) -> dict:
    """Re-run the MSDB cascade against the current ModelSEEDDatabase JSON data.

    Returns a ``{rxn_id: reversibility}`` map produced by ``reversibility_lib``'s
    port of ``Estimate_Reaction_Reversibility.py`` (EQ pass, ``gc_first=True``).
    Optionally persists the map to CSV (two-column ``rxn_id,reversibility`` plus
    a ``status`` column) and JSON under ``results/``.

    Parameters
    ----------
    out_csv : Path, optional
        Where to write the CSV.  Defaults to
        ``results/rxn_directions_cascade_live.csv``.  Pass an explicit
        ``None`` only via the keyword if you want to skip writing CSV (the
        default kicks in when the argument is not supplied at all).  Must
        live under ``results/`` when supplied.
    out_json : Path, optional
        Where to write the JSON map.  Defaults to
        ``results/rxn_directions_cascade_live.json``.  Same path constraint.
    cfg : ``reversibility_lib.ReversibilityConfig``, optional
        Knob-set for the cascade.  Default = byte-for-byte upstream baseline
        (``ReversibilityConfig()``).

    Notes
    -----
    Imports of ``reversibility_lib`` and ``BiochemPy`` are deferred to this
    call so existing notebook cells that only need the lighter helpers in
    this module don't pay the load cost.  ModelSEEDDatabase is treated as
    read-only -- ``BiochemPy.Reactions().loadReactions()`` only reads JSON.
    """
    import sys

    # Make the MSDB BiochemPy importable on first call only.
    msdb_libs = str(_MSDB_LIBS_PATH)
    if msdb_libs not in sys.path:
        sys.path.insert(0, msdb_libs)

    import reversibility_lib as lib  # lazy: heavy import path warm-up
    from BiochemPy import Reactions  # lazy: pulls the MSDB JSON parser

    if cfg is None:
        cfg = lib.ReversibilityConfig()

    # Resolve defaults -- if the caller passed ``None`` (or omitted), use the
    # canonical paths under ``results/``.  When a path is supplied, enforce
    # it lives under ``results/`` so this helper cannot scribble elsewhere.
    if out_csv is None:
        out_csv = RESULTS_DIR / "rxn_directions_cascade_live.csv"
    if out_json is None:
        out_json = RESULTS_DIR / "rxn_directions_cascade_live.json"
    out_csv = Path(out_csv)
    out_json = Path(out_json)
    assert out_csv is not None, "out_csv must be a path"
    assert out_json is not None, "out_json must be a path"
    _assert_under_results(out_csv, "out_csv")
    _assert_under_results(out_json, "out_json")

    # Load reactions (read-only) and run the cascade.
    rxns = Reactions().loadReactions()
    cascade = lib.run_cascade(rxns, db_level="EQ", cfg=cfg, gc_first=True)

    direction_map = {rid: rev for rid, (_status, rev) in cascade.items()}
    status_map = {rid: status for rid, (status, _rev) in cascade.items()}

    # Persist.  Parent dir is always ``results/`` thanks to the assertion.
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["rxn_id", "reversibility", "status"])
        for rid in sorted(direction_map):
            writer.writerow([rid, direction_map[rid], status_map[rid]])

    out_json.parent.mkdir(parents=True, exist_ok=True)
    with out_json.open("w") as fh:
        json.dump(direction_map, fh, sort_keys=True)

    return direction_map


# ---------------------------------------------------------------------------
# Per-source direction snapshot (PR #263's new ``thermodynamics`` layout)
# ---------------------------------------------------------------------------
def snapshot_msdb_per_source(
    branch: str = "origin/dev",
    source: str = "Group contribution",
    out_csv: Optional[Path] = None,
    out_json: Optional[Path] = None,
    repo: Path = MSDB_DEFAULT_ROOT,
    n_shards: int = MSDB_N_SHARDS,
) -> dict:
    """Extract a per-source ``{rxn_id: operator}`` map from MSDB JSON shards.

    Reads ``Biochemistry/reaction_NN.json`` for ``NN in 0..n_shards-1`` from
    ``branch`` via ``git show`` (the MSDB working tree is never touched) and
    looks at ``rxn['thermodynamics'][source]`` for each reaction.  The new
    PR #263 format puts ``[dg, dge, operator]`` triples under each source
    label inside ``thermodynamics``.

    Skips reactions whose entry is missing, has a sentinel energy
    (``sublist[0] == 10000000`` or ``None``), or is shorter than length 3.
    The captured operator (``sublist[2]``) is normally one of
    ``'>' / '<' / '=' / '?'``.

    Optionally persists the result as CSV with columns
    ``rxn_id, operator, dg, dge`` (path defaults to
    ``results/rxn_directions_{slug}.csv``) and as JSON
    (``{rxn_id: operator}``, default ``results/rxn_directions_{slug}.json``).
    Any caller-supplied paths are asserted to live under
    ``PROJECT_ROOT / 'results'`` so MSDB stays read-only.

    Returns the ``{rxn_id: operator}`` dict.
    """
    slug = _source_slug(source)
    if out_csv is None:
        out_csv = RESULTS_DIR / f"rxn_directions_{slug}.csv"
    if out_json is None:
        out_json = RESULTS_DIR / f"rxn_directions_{slug}.json"
    out_csv = Path(out_csv)
    out_json = Path(out_json)
    _assert_under_results(out_csv, "out_csv")
    _assert_under_results(out_json, "out_json")

    operator_map: dict = {}
    rows: list = []  # (rxn_id, operator, dg, dge), held for CSV write
    for i in range(n_shards):
        shard = MSDB_JSON_SHARD_FMT.format(i)
        text = _git_show(branch, shard, repo)
        data = json.loads(text)
        for rxn in data:
            rid = rxn.get("id")
            if not rid:
                continue
            thermo = rxn.get("thermodynamics") or {}
            sub = thermo.get(source)
            if sub is None:
                continue
            if len(sub) < 3:
                continue
            dg = sub[0]
            if dg == _THERMO_SENTINEL or dg is None:
                continue
            dge = sub[1]
            operator = sub[2]
            operator_map[rid] = operator
            rows.append((rid, operator, dg, dge))

    # Sort once so on-disk artifacts are stable across runs.
    rows.sort(key=lambda r: r[0])

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["rxn_id", "operator", "dg", "dge"])
        for r in rows:
            writer.writerow(r)

    out_json.parent.mkdir(parents=True, exist_ok=True)
    with out_json.open("w") as fh:
        json.dump(operator_map, fh, sort_keys=True)

    return operator_map


def load_per_source_operators(csv_path: Path) -> dict:
    """Load a ``{rxn_id: operator}`` map from a per-source CSV.

    Reads CSVs produced by :func:`snapshot_msdb_per_source` (columns
    ``rxn_id, operator, dg, dge``).  Rows whose ``operator`` is empty
    are skipped so callers don't accidentally rewrite a reaction's bounds
    to ``(-1000, 1000)`` based on a blank.
    """
    p = Path(csv_path)
    if not p.exists():
        raise FileNotFoundError(f"Per-source CSV not found: {p}")
    out: dict = {}
    with p.open() as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rid = (row.get("rxn_id") or "").strip()
            op = (row.get("operator") or "").strip()
            if not rid or not op:
                continue
            out[rid] = op
    return out


def panel_growth_with_source(
    panel_ids: Iterable[str],
    source_csv_path: Path,
    baseline_map: Optional[Mapping[str, str]] = None,
    n_workers: int = 4,
) -> list:
    """Run the FBA panel using a per-source direction map.

    Loads the operator map via :func:`load_per_source_operators`, then
    delegates to :func:`growth_heuristics.run_panel`.  The contract is
    intentional: reactions *not* in the source map are left at their
    on-disk KBase bounds because ``growth_heuristics.override_bounds``
    only touches reactions whose ``seed.reaction`` annotation is in the
    supplied map (see ``override_bounds`` -- ``reversibility_map.get(seed)``
    returns ``None`` and the loop ``continue``s).

    Returns the raw ``run_panel`` results list (one dict per model).
    """
    source_map = load_per_source_operators(source_csv_path)
    return gh.run_panel(
        list(panel_ids),
        reversibility_map=source_map,
        baseline_map=dict(baseline_map) if baseline_map is not None else None,
        n_workers=n_workers,
    )


# ---------------------------------------------------------------------------
# Override-transition bookkeeping (wide schema)
# ---------------------------------------------------------------------------
# Skip the same exchange/sink/demand/biomass prefixes that growth_heuristics
# uses everywhere else; this is the canonical "real" reaction filter.
_OVERRIDE_SKIP_PREFIXES = ("EX_", "SK_", "DM_", "bio")

OVERRIDE_TRANSITION_COLS = (
    "model_id",
    "source",
    "n_overrides_applied",
    "n_fwd_to_rev",
    "n_rev_to_fwd",
    "n_fwd_to_reversible",
    "n_reversible_to_fwd",
    "n_rev_to_reversible",
    "n_reversible_to_rev",
    "n_unchanged",
    "n_other",
)

_TRANSITION_BUCKETS = {
    ("fwd", "rev"):         "n_fwd_to_rev",
    ("rev", "fwd"):         "n_rev_to_fwd",
    ("fwd", "reversible"):  "n_fwd_to_reversible",
    ("reversible", "fwd"):  "n_reversible_to_fwd",
    ("rev", "reversible"):  "n_rev_to_reversible",
    ("reversible", "rev"):  "n_reversible_to_rev",
}


def _bound_class(rev: Optional[str]) -> str:
    """Reduce a reversibility flag to {'fwd', 'rev', 'reversible'}.

    Mirrors ``gh._bounds_for_rev`` so the labels we emit match what the
    override step would actually do at the bounds level:
        '>' -> (0, +B)         -> 'fwd'
        '<' -> (-B, 0)         -> 'rev'
        anything else ('?', '=', '', None) -> (-B, +B) -> 'reversible'
    """
    lb, ub = gh._bounds_for_rev(rev or "")
    if ub > 0 and lb >= 0:
        return "fwd"
    if lb < 0 and ub <= 0:
        return "rev"
    return "reversible"


def load_msdb_reversibility_map(snapshot_csv_path: Path) -> dict:
    """Load ``{rxn_id: reversibility}`` from a CSV produced by snapshot_msdb.

    Expects a header row with columns ``rxn_id`` and ``reversibility``.
    Empty reversibilities are kept as-is so ``_bound_class`` can decide
    how to coerce them (it treats unknowns as ``reversible``).
    """
    p = Path(snapshot_csv_path)
    if not p.exists():
        raise FileNotFoundError(f"MSDB snapshot CSV not found: {p}")
    out: dict = {}
    with p.open() as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rid = (row.get("rxn_id") or "").strip()
            rev = (row.get("reversibility") or "").strip()
            if rid:
                out[rid] = rev
    return out


def override_transitions(
    panel_ids: Iterable[str],
    source_csv_path: Path,
    msdb_snapshot_csv_path: Path,
    source_slug: str,
) -> list:
    """Per-model bound-class transition tally for one per-source variant.

    For each model, walks reactions whose ``annotation['seed.reaction']`` id
    is present in the per-source operator map.  Compares the **baseline**
    bound class (derived from the MSDB-snapshot reversibility -- the same
    classification the on-disk KBase bounds encode) against the
    **variant** bound class (derived from the source operator) and tallies
    transitions into the six wide ``n_<before>_to_<after>`` columns.

    Returns
    -------
    list[dict]
        One row per model with the wide schema declared in
        :data:`OVERRIDE_TRANSITION_COLS`.  Same shape the standalone
        ``build_thermo_source_network_tables.py`` CLI used to emit.
    """
    source_map = load_per_source_operators(source_csv_path)
    msdb_map = load_msdb_reversibility_map(msdb_snapshot_csv_path)
    rows = []
    for mid in panel_ids:
        row = {col: 0 for col in OVERRIDE_TRANSITION_COLS}
        row["model_id"] = mid
        row["source"] = source_slug
        model_path = gh.MODELS_DIR / f"{mid}.json"
        if not model_path.exists():
            rows.append(row)
            continue
        with model_path.open() as fh:
            model = json.load(fh)
        n_overrides_applied = 0
        for rxn in model.get("reactions", []):
            rid = rxn.get("id") or ""
            if rid.startswith(_OVERRIDE_SKIP_PREFIXES):
                continue
            anno = rxn.get("annotation") or {}
            seed = anno.get("seed.reaction")
            if isinstance(seed, list):
                seed = seed[0] if seed else None
            if not seed:
                continue
            new_op = source_map.get(seed)
            if new_op is None:
                continue
            n_overrides_applied += 1
            base_rev = msdb_map.get(seed)
            before = _bound_class(base_rev) if base_rev else "reversible"
            after = _bound_class(new_op)
            if before == after:
                row["n_unchanged"] += 1
                continue
            bucket = _TRANSITION_BUCKETS.get((before, after))
            if bucket is None:
                row["n_other"] += 1
            else:
                row[bucket] += 1
        row["n_overrides_applied"] = n_overrides_applied
        rows.append(row)
    return rows


def source_coverage(
    panel_ids: Iterable[str],
    source_csv_path: Path,
) -> dict:
    """Per-model coverage stats for a per-source direction map.

    For each model in ``panel_ids`` loads its JSON, walks every reaction
    whose id does **not** start with ``EX_ / SK_ / DM_ / bio``, and counts:

      - ``n_reactions``:        total non-(exchange/sink/demand/biomass)
        reactions
      - ``n_with_seed_anno``:   how many of those carry an
        ``annotation['seed.reaction']`` id
      - ``n_covered_by_source``: how many of *those* ids are in the source's
        operator map
      - ``frac_covered``:        ``n_covered_by_source / max(1, n_with_seed_anno)``

    The biomass/sink/demand/exchange skip list mirrors the same prefix set
    used by ``growth_heuristics.fba_one`` when stripping bounds.
    """
    source_map = load_per_source_operators(source_csv_path)
    out: dict = {}
    skip_prefixes = ("EX_", "SK_", "DM_", "bio")
    for mid in panel_ids:
        model_path = gh.MODELS_DIR / f"{mid}.json"
        if not model_path.exists():
            out[mid] = {
                "n_reactions": 0,
                "n_with_seed_anno": 0,
                "n_covered_by_source": 0,
                "frac_covered": 0.0,
            }
            continue
        with model_path.open() as fh:
            model = json.load(fh)
        n_rxn = 0
        n_anno = 0
        n_cov = 0
        for rxn in model.get("reactions", []):
            rid = rxn.get("id") or ""
            if rid.startswith(skip_prefixes):
                continue
            n_rxn += 1
            anno = rxn.get("annotation") or {}
            seed = anno.get("seed.reaction")
            # KBase models sometimes store a list under the annotation key.
            if isinstance(seed, list):
                seed = seed[0] if seed else None
            if not seed:
                continue
            n_anno += 1
            if seed in source_map:
                n_cov += 1
        out[mid] = {
            "n_reactions": n_rxn,
            "n_with_seed_anno": n_anno,
            "n_covered_by_source": n_cov,
            "frac_covered": (n_cov / n_anno) if n_anno else 0.0,
        }
    return out
