#!/usr/bin/env python3
"""
Build notebooks/09_ReactionDirectionPipeline.ipynb.

The notebook re-runs panel growth (the 100 representative core models) for
several sources of reaction direction:

  (a) directions loaded from a user-specified path (RXN_DIRECTION_SOURCE)
      -- defaults to a local CSV under ``results/`` seeded on first run
      from the dev-branch MSDB snapshot;
  (b) an editable overlay CSV that the user mutates in-notebook
      (``rxn_directions_overlay.csv``);
  (c) an in-notebook patch function (``patch_directions``) that rewrites
      the local CSV and re-runs FBA.

MSDB is read-only -- snapshots happen via ``git show`` and every output
lands under ``core_models_analysis/results/`` or ``.kbcache/``.

Generated, not hand-edited.  Run this script to rebuild.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import nbformat as nbf

ROOT = Path("/scratch/ctaylor/core_models_analysis")
NOTEBOOK_PATH = ROOT / "notebooks" / "09_ReactionDirectionPipeline.ipynb"


def md(text: str) -> nbf.NotebookNode:
    return nbf.v4.new_markdown_cell(textwrap.dedent(text).strip("\n"))


def code(src: str) -> nbf.NotebookNode:
    return nbf.v4.new_code_cell(textwrap.dedent(src).strip("\n"))


# ---------------------------------------------------------------------------
# Cell sources
# ---------------------------------------------------------------------------
INTRO_MD = """
# 09 -- Reaction Direction Pipeline

This notebook re-runs **panel FBA growth** for the 100 representative core
models across several sources of reaction direction (reversibility):

1. **(a) From a user-specified source path** -- ``RXN_DIRECTION_SOURCE``
   points at a CSV / TSV / JSON file keyed by ``rxn_id``.  On first run
   it is seeded from the dev-branch MSDB snapshot.
2. **(b) From an editable overlay** -- the notebook copies the current
   map into ``rxn_directions_overlay.csv``, applies a worked example
   mutation (every reaction with ``deltag > +5`` forced to ``<``), and
   re-runs FBA on the mutated map.
3. **(c) From an in-notebook patch function** -- a user-editable
   ``patch_directions(direction_map) -> dict`` rewrites the local CSV
   in place and re-runs FBA.
4. **(d) From a live re-run of the MSDB cascade** -- the
   ``Snapshot cascade-live`` cell re-runs ``reversibility_lib``'s port of
   ``Estimate_Reaction_Reversibility.py`` against the current MSDB JSON
   data and persists the resulting map to
   ``results/rxn_directions_cascade_live.csv`` / ``.json``.  The
   cross-variant summary picks it up automatically as
   ``RUNS['cascade_live']``.

### Constraint -- MSDB is read-only

This pipeline **never modifies anything under**
``/scratch/ctaylor/ModelSEEDDatabase``.  Direction snapshots are taken
with ``git show <branch>:<shard>`` so the MSDB working tree is never
touched.  All working data lives under
``/scratch/ctaylor/core_models_analysis/results/`` (durable outputs) or
the per-notebook ``.kbcache/`` (memoised intermediates).
"""

SETUP_CELL = """
from pathlib import Path
import sys, json, copy, shutil, time, csv

PROJECT_ROOT = Path('/scratch/ctaylor/core_models_analysis')
RESULTS      = PROJECT_ROOT / 'results'
SCRIPTS      = PROJECT_ROOT / 'scripts'

# Make project helpers importable.
sys.path.insert(0, str(SCRIPTS))

# KBUtils NotebookSession -- caches heavy intermediates next to the .ipynb.
from kbutillib.notebook import NotebookSession
session = NotebookSession.for_notebook(
    notebook_file=__file__ if '__file__' in dir() else None,
    project_name='core_models_analysis',
)
print('Cache directory:', session.kbcache_dir)
print('Notebook name :', session.notebook_name)

import direction_pipeline as dp
import growth_heuristics as gh
"""

PARAMS_CELL = """
# === Parameters =====================================================
# Edit these to change which source feeds the pipeline.  The defaults
# point at locally cached CSVs -- nothing here writes back to MSDB.

# (a) Source-of-truth for the 'from_source' run.  Override to point at
# any TSV / CSV / JSON keyed by rxn_id with a `reversibility` (or
# `direction`) column / key.  On first run, this file is seeded from the
# MSDB dev-branch snapshot below.
RXN_DIRECTION_SOURCE = Path('/scratch/ctaylor/core_models_analysis/results/rxn_directions_local.csv')

# Working overlay file -- (b) reads from / writes to this path.
OVERLAY_PATH         = Path('/scratch/ctaylor/core_models_analysis/results/rxn_directions_overlay.csv')

# (b)'s sample mutation threshold: any reaction with deltag > this value
# is forced to '<' (reverse-only).
OVERLAY_DELTAG_THRESHOLD = 5.0

# Panel of model ids to grow.
PANEL_IDS_PATH       = Path('/scratch/ctaylor/core_models_analysis/results/selected_ids.txt')

# Read-only MSDB checkout -- we only ever invoke `git show <branch>:...`
# against this path.
MSDB_ROOT            = Path('/scratch/ctaylor/ModelSEEDDatabase')

# Snapshot CSVs written by the next cell.
MSDB_DEV_SNAPSHOT    = RESULTS / 'rxn_directions_msdb_dev.csv'
MSDB_CLAUDE_SNAPSHOT = RESULTS / 'rxn_directions_msdb_claude.csv'

# FBA worker count -- keep modest; the pool spawns one process per worker.
N_WORKERS = 4

# Load the panel ids.
PANEL_IDS = PANEL_IDS_PATH.read_text().split()
print(f'Panel: {len(PANEL_IDS)} models from {PANEL_IDS_PATH}')

# Registry used by the cross-variant summary at the bottom.
RUNS = {}
"""

SNAPSHOT_CELL = """
# === Snapshot MSDB direction columns (read-only) =====================
# Uses `git show <branch>:Biochemistry/reaction_NN.tsv` -- the MSDB
# working tree is *never* modified.  Idempotent: skips work if the
# snapshot CSV already exists.

def _ensure_snapshot(branch: str, out_path: Path):
    if out_path.exists():
        n = sum(1 for _ in open(out_path)) - 1
        print(f'  exists  -> {out_path}  ({n} rows)')
        return out_path
    t = time.time()
    p = dp.snapshot_msdb(branch=branch, out_path=out_path, repo=MSDB_ROOT)
    n = sum(1 for _ in open(p)) - 1
    print(f'  wrote   -> {p}  ({n} rows in {time.time()-t:.1f}s)')
    return p

print('Snapshotting MSDB branches (read-only via git show):')
_ensure_snapshot('dev',            MSDB_DEV_SNAPSHOT)
_ensure_snapshot('claude-changes', MSDB_CLAUDE_SNAPSHOT)

# Convenience maps for downstream cells.
DEV_MAP    = dp.load_directions_from_path(MSDB_DEV_SNAPSHOT)
CLAUDE_MAP = dp.load_directions_from_path(MSDB_CLAUDE_SNAPSHOT)
print(f'dev   : {len(DEV_MAP)} reactions')
print(f'claude: {len(CLAUDE_MAP)} reactions')
"""

CASCADE_LIVE_CELL = """
# === Snapshot cascade-live (re-run the MSDB cascade on JSON data) =====
# Calls dp.run_cascade_live(...) -- a lazy port of the upstream
# Estimate_Reaction_Reversibility.py cascade -- against the live MSDB
# JSON.  Writes:
#
#   results/rxn_directions_cascade_live.csv
#   results/rxn_directions_cascade_live.json
#
# MSDB is read-only: the helper only loads JSON via BiochemPy and the
# output paths are asserted to live under results/.

CASCADE_LIVE_CSV  = RESULTS / 'rxn_directions_cascade_live.csv'
CASCADE_LIVE_JSON = RESULTS / 'rxn_directions_cascade_live.json'
STORED_REV_JSON   = RESULTS / 'rev_map_dev.json'

t = time.time()
cascade_live_map = dp.run_cascade_live(
    out_csv=CASCADE_LIVE_CSV,
    out_json=CASCADE_LIVE_JSON,
)
print(f'cascade_live: {len(cascade_live_map)} reactions  ({time.time()-t:.1f}s)')
print(f'  csv  -> {CASCADE_LIVE_CSV}')
print(f'  json -> {CASCADE_LIVE_JSON}')

# Compare vs the stored rev map snapshot for the dev branch.
if STORED_REV_JSON.exists():
    stored_rev = json.loads(STORED_REV_JSON.read_text())
    shared = set(cascade_live_map) & set(stored_rev)
    n_match = sum(1 for rid in shared if cascade_live_map[rid] == stored_rev[rid])
    print(
        f'  vs {STORED_REV_JSON.name}: '
        f'shared={len(shared)}  match={n_match}  diff={len(shared) - n_match}  '
        f'only_live={len(set(cascade_live_map) - set(stored_rev))}  '
        f'only_stored={len(set(stored_rev) - set(cascade_live_map))}'
    )
else:
    print(f'  (stored rev map not found at {STORED_REV_JSON}; skipping compare)')
"""

REQ_A_CELL = """
# === (a) Re-run growth from the SPECIFIED source path ================
# RXN_DIRECTION_SOURCE may point anywhere -- a CSV, TSV, or JSON keyed
# by rxn_id.  On first run we seed the local default by copying the
# dev-branch MSDB snapshot, so the pipeline is runnable end-to-end out
# of the box.

if not RXN_DIRECTION_SOURCE.exists():
    print(f'Seeding {RXN_DIRECTION_SOURCE} from dev snapshot ...')
    RXN_DIRECTION_SOURCE.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(MSDB_DEV_SNAPSHOT, RXN_DIRECTION_SOURCE)

source_map = dp.load_directions_from_path(RXN_DIRECTION_SOURCE)
print(f'Loaded {len(source_map)} directions from {RXN_DIRECTION_SOURCE}')

t = time.time()
source_run = dp.panel_growth(PANEL_IDS, source_map, n_workers=N_WORKERS)
print(f'panel FBA done in {time.time()-t:.1f}s -- totals:', source_run['totals'])

RUNS['from_source'] = source_run
"""

REQ_B_CELL = """
# === (b) Edit a TEMPORARY overlay, then re-run =======================
# Loads the source map, mutates it (worked example: force every
# reaction with deltag > OVERLAY_DELTAG_THRESHOLD to '<'), saves the
# result to OVERLAY_PATH, and re-runs FBA against the mutated map.
#
# Edit the body of `_mutate(...)` below to try other temporary changes.

# Pull deltag from the dev snapshot so the mutation has thermodynamic
# context regardless of which source the user pointed at.
def _load_deltag_table(snapshot_path: Path) -> dict:
    out = {}
    with open(snapshot_path) as fh:
        for row in csv.DictReader(fh):
            try:
                out[row['rxn_id']] = float(row.get('deltag') or 'nan')
            except ValueError:
                out[row['rxn_id']] = float('nan')
    return out

deltag_by_rxn = _load_deltag_table(MSDB_DEV_SNAPSHOT)

def _mutate(direction_map: dict) -> dict:
    out = dict(direction_map)
    n_flipped = 0
    for rid, rev in list(out.items()):
        dg = deltag_by_rxn.get(rid)
        if dg is not None and dg == dg and dg > OVERLAY_DELTAG_THRESHOLD:
            if out[rid] != '<':
                out[rid] = '<'
                n_flipped += 1
    print(f'  mutation: forced {n_flipped} rxns (deltag > {OVERLAY_DELTAG_THRESHOLD}) to "<"')
    return out

overlay_map = _mutate(source_map)
dp.save_directions_to_path(overlay_map, OVERLAY_PATH)
print(f'overlay written to {OVERLAY_PATH}')

t = time.time()
overlay_run = dp.panel_growth(PANEL_IDS, overlay_map, n_workers=N_WORKERS)
print(f'panel FBA done in {time.time()-t:.1f}s -- totals:', overlay_run['totals'])

RUNS['from_overlay'] = overlay_run
"""

REQ_C_CELL = """
# === (c) IN-NOTEBOOK PATCH -- modify the local CSV from script =======
# Edit the body of `patch_directions` below to alter the directions
# however you like.  Running this cell:
#
#   1. loads the current local map (RXN_DIRECTION_SOURCE),
#   2. passes it through `patch_directions`,
#   3. rewrites RXN_DIRECTION_SOURCE in place with the patched map,
#   4. re-runs panel FBA against the patched map.
#
# The worked example below flips every '=' direction to '>'.

def patch_directions(direction_map: dict) -> dict:
    '''User-editable: take the current direction map and return a new one.

    Default behaviour: rewrite every fully-reversible ('=') reaction
    as forward-only ('>').  Replace the body with whatever experiment
    you want to try.
    '''
    out = dict(direction_map)
    flipped = 0
    for rid, rev in list(out.items()):
        if rev == '=':
            out[rid] = '>'
            flipped += 1
    print(f'  patch_directions: flipped {flipped} "=" rxns to ">"')
    return out

current_map = dp.load_directions_from_path(RXN_DIRECTION_SOURCE)
patched_map = patch_directions(current_map)

# Step 3: write the patched map back to the local CSV.
dp.save_directions_to_path(patched_map, RXN_DIRECTION_SOURCE)
print(f'patched local CSV: {RXN_DIRECTION_SOURCE} ({len(patched_map)} rxns)')

# Step 4: re-run panel FBA against the patched map.
t = time.time()
patched_run = dp.panel_growth(PANEL_IDS, patched_map, n_workers=N_WORKERS)
print(f'panel FBA done in {time.time()-t:.1f}s -- totals:', patched_run['totals'])

RUNS['from_patch'] = patched_run
"""

SUMMARY_CELL = """
# === Cross-variant summary ===========================================
# One row per variant with grow-count and mean / median flux.  Also
# emits a long-form DataFrame (variant x model_id) for plotting.

# Add the cascade_live variant to the registry by loading the JSON map
# written by the snapshot cell and running the FBA panel against it.
# baseline_map=None forces a full rebind so the cascade map is
# authoritative.
if 'cascade_live' not in RUNS:
    if not CASCADE_LIVE_JSON.exists():
        raise FileNotFoundError(
            f'cascade-live JSON missing: {CASCADE_LIVE_JSON}.  '
            'Run the "Snapshot cascade-live" cell above first.'
        )
    cascade_live_loaded = json.loads(CASCADE_LIVE_JSON.read_text())
    print(f'cascade_live: panel FBA on {len(cascade_live_loaded)} directions ...')
    t = time.time()
    cascade_live_results = gh.run_panel(
        PANEL_IDS,
        reversibility_map=cascade_live_loaded,
        baseline_map=None,
        n_workers=N_WORKERS,
    )
    cascade_live_totals = {
        'n_models': len(cascade_live_results),
        'n_grow': sum(1 for r in cascade_live_results if r.get('grows')),
    }
    print(f'  done in {time.time()-t:.1f}s -- totals:', cascade_live_totals)
    RUNS['cascade_live'] = {
        'results': cascade_live_results,
        'totals': cascade_live_totals,
    }

name_to_results = {name: run['results'] for name, run in RUNS.items()}
totals_df = dp.variant_totals(name_to_results)
print('Per-variant totals:')
print(totals_df.to_string(index=False))

long_df = dp.compare_runs(name_to_results)
print(f'\\nLong-form table: {len(long_df)} rows  ({long_df["variant"].nunique()} variants)')

# Persist to results/ + cache so downstream visualization phase can pick
# them up without rerunning FBA.
SUMMARY_CSV = RESULTS / 'direction_pipeline_summary.csv'
LONG_CSV    = RESULTS / 'direction_pipeline_long.csv'
totals_df.to_csv(SUMMARY_CSV, index=False)
long_df.to_csv(LONG_CSV, index=False)
print(f'wrote {SUMMARY_CSV}')
print(f'wrote {LONG_CSV}')

try:
    session.cache.save('direction_pipeline_runs_v1', RUNS, type_hint='dict',
                       metadata={'panel': str(PANEL_IDS_PATH)})
except Exception as exc:
    print(f'(cache save skipped: {exc})')
"""

FIGURES_MD = """
## Figures

The visualization phase (Phase 4) writes its figures under
`reports/figures/direction_pipeline/`.  Once that phase has run, the
expected outputs are:

- `reports/figures/direction_pipeline/grow_counts_by_variant.png`
- `reports/figures/direction_pipeline/flux_distribution_by_variant.png`
- `reports/figures/direction_pipeline/per_model_heatmap.png`

The summary tables this notebook writes
(`results/direction_pipeline_summary.csv` and
`results/direction_pipeline_long.csv`) are the canonical input for those
figures.
"""


# ---------------------------------------------------------------------------
# Notebook assembly
# ---------------------------------------------------------------------------
def build_notebook() -> nbf.NotebookNode:
    cells = [
        md(INTRO_MD),                 # cell 0 -- markdown / intro + constraint
        code(SETUP_CELL),             # cell 1 -- imports + KBUtils session
        code(PARAMS_CELL),            # cell 2 -- parameter cell
        code(SNAPSHOT_CELL),          # cell 3 -- MSDB snapshots (idempotent)
        code(CASCADE_LIVE_CELL),      # cell 4 -- (d) live cascade snapshot
        code(REQ_A_CELL),             # cell 5 -- (a) source-driven run
        code(REQ_B_CELL),             # cell 6 -- (b) overlay-driven run
        code(REQ_C_CELL),             # cell 7 -- (c) in-notebook patch
        code(SUMMARY_CELL),           # cell 8 -- cross-variant table (adds cascade_live)
        md(FIGURES_MD),               # cell 9 -- figure links
    ]
    nb = nbf.v4.new_notebook()
    nb.cells = cells
    nb.metadata = {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python"},
    }
    return nb


def main() -> None:
    nb = build_notebook()
    NOTEBOOK_PATH.parent.mkdir(parents=True, exist_ok=True)
    nbf.write(nb, NOTEBOOK_PATH)
    print(f"wrote {NOTEBOOK_PATH}  ({len(nb.cells)} cells)")


if __name__ == "__main__":
    main()
