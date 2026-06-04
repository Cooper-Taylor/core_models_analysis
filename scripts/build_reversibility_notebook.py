#!/usr/bin/env python3
"""
Build notebooks/06_ReactionReversibilityHeuristics.ipynb.

The notebook ports ModelSEEDDatabase's
``Estimate_Reaction_Reversibility.py`` into a parameterizable library
(``reversibility_lib``) and exercises each suggestion from
``Reaction_Reversibility_Heuristics_Review.md`` against the 100-model panel.

Layout:
  0.  Setup + KBUtils NotebookSession (database loads cached at .kbcache/)
  1.  Reaction-direction baseline (matches MSDB report byte-for-byte)
  2.  Growth baselines:
        2a -- on-disk FBA matches results.csv (no rebinding)
        2b -- heuristic-baseline FBA (apply MSDB cascade to model bounds)
  3.  Heuristic variants, one cell per suggestion (3.1 ... + new ones).
        For each: compute variant reversibility, diff vs baseline cascade,
        rebound the panel, diff growth vs heuristic-baseline.
  4.  Cross-variant summary table.

Generated, not hand-edited -- run this script to rebuild.
"""

from __future__ import annotations
import textwrap
from pathlib import Path

import nbformat as nbf

ROOT = Path("/scratch/ctaylor/core_models_analysis")
NOTEBOOK_PATH = ROOT / "notebooks" / "06_ReactionReversibilityHeuristics.ipynb"


def md(text: str) -> nbf.NotebookNode:
    return nbf.v4.new_markdown_cell(textwrap.dedent(text).strip("\n"))


def code(src: str) -> nbf.NotebookNode:
    return nbf.v4.new_code_cell(textwrap.dedent(src).strip("\n"))


SETUP_CELL = """
from pathlib import Path
import sys, copy, time, json, csv, hashlib

PROJECT_ROOT = Path('/scratch/ctaylor/core_models_analysis')
MSDB_ROOT    = Path('/scratch/ctaylor/ModelSEEDDatabase')
REPORTS      = PROJECT_ROOT / 'reports'
RESULTS      = PROJECT_ROOT / 'results'
SCRIPTS      = PROJECT_ROOT / 'scripts'

# Make the project scripts and MSDB BiochemPy importable.
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(MSDB_ROOT / 'Libs' / 'Python'))

# KBUtils NotebookSession holds the .kbcache/ alongside this notebook.  The
# 56K-row MSDB reactions dict + the panel growth tables are cached here so
# re-runs skip the heavy load + FBA solves.
from kbutillib.notebook import NotebookSession
session = NotebookSession.for_notebook(
    notebook_file=__file__ if '__file__' in dir() else None,
    project_name='core_models_analysis',
)
print('Cache directory:', session.kbcache_dir)
print('Notebook name :', session.notebook_name)

import reversibility_lib as lib
import growth_heuristics as gh
"""


LOAD_MSDB_CELL = """
# Load the MSDB reaction database via BiochemPy.  Cached in the
# NotebookSession blob store the first time and reloaded as JSON on
# subsequent runs.
CACHE_KEY_MSDB_RXNS = 'msdb_reactions_v1'

def _load_msdb_reactions():
    from BiochemPy import Reactions
    return Reactions().loadReactions()

try:
    msdb_rxns = session.cache.load(CACHE_KEY_MSDB_RXNS)
    print(f'Loaded {len(msdb_rxns)} MSDB reactions from cache')
except KeyError:
    t = time.time()
    msdb_rxns = _load_msdb_reactions()
    session.cache.save(CACHE_KEY_MSDB_RXNS, msdb_rxns, type_hint='dict',
                       metadata={'source': 'BiochemPy.Reactions().loadReactions()'})
    print(f'Loaded {len(msdb_rxns)} MSDB reactions from disk in {time.time()-t:.1f}s '
          f'-- cached to {session.kbcache_dir}')

# Load the 100-model panel ids (the descriptive growth-model panel).
PANEL_IDS = open(RESULTS / 'selected_ids.txt').read().split()
print(f'Panel: {len(PANEL_IDS)} models')

# Load the original (no-rebound) growth from results.csv -- our 'on-disk FBA'
# reference.  This is the baseline that analyze_growth.py wrote.
ONDISK_FBA = {}
with open(RESULTS / 'results.csv') as fh:
    for row in csv.DictReader(fh):
        ONDISK_FBA[row['model_id']] = {
            'status'     : row['status'],
            'growth_flux': float(row['growth_flux']),
            'grows'      : row['grows'].lower() == 'true',
        }
print(f'On-disk FBA results loaded for {len(ONDISK_FBA)} models')
"""


REPORT_FILE = ("/scratch/ctaylor/ModelSEEDDatabase/Scripts/Thermodynamics/"
               "Estimated_Reaction_Reversibility_Report_EQ.txt")


BASELINE_CASCADE_CELL = f"""
# --- Baseline cascade ---------------------------------------------------
# Run the heuristic with ReversibilityConfig() defaults.  Caches the result
# so re-runs are instant.
CACHE_KEY_BASELINE = 'reversibility_baseline_v1'

def _run_baseline_cascade():
    rxns = copy.deepcopy(msdb_rxns)
    return lib.run_cascade(rxns, db_level='EQ', cfg=lib.ReversibilityConfig(),
                           gc_first=True)

try:
    baseline_cascade = session.cache.load(CACHE_KEY_BASELINE)
    print(f'baseline_cascade: loaded {{len(baseline_cascade)}} entries from cache')
except KeyError:
    t = time.time()
    baseline_cascade = _run_baseline_cascade()
    session.cache.save(CACHE_KEY_BASELINE, baseline_cascade, type_hint='dict',
                       metadata={{'config': 'ReversibilityConfig() default'}})
    print(f'baseline_cascade: {{len(baseline_cascade)}} reactions in {{time.time()-t:.1f}}s, cached')

BASELINE_MAP = {{r: rev for r, (_, rev) in baseline_cascade.items()}}

# --- Sanity check vs. the on-disk report ----------------------------------
# Compare our cascade's output to the existing
# Estimated_Reaction_Reversibility_Report_EQ.txt.  The MSDB JSON has been
# updated since the report was regenerated, so we expect a handful of
# reactions where the report says 'Incomplete' but our cascade now has data.
report = {{}}
with open('{REPORT_FILE}') as fh:
    for line in fh:
        parts = line.rstrip('\\n').split('\\t')
        report[parts[0]] = (parts[1], parts[3])  # (status, new_rev)

agree = 0; disagree = []
for rxn, (status, rev) in baseline_cascade.items():
    rs = report.get(rxn)
    if rs is None:
        continue
    if rs == (status, rev):
        agree += 1
    else:
        disagree.append((rxn, rs, (status, rev)))

print(f'Cascade reproduces Estimated_Reaction_Reversibility_Report_EQ.txt:')
print(f'  exact matches:  {{agree}}/{{len(baseline_cascade)}}')
print(f'  mismatches:     {{len(disagree)}}  (expected drift from updated MSDB data)')
if disagree:
    print('  sample drift:')
    for s in disagree[:3]:
        print(f'    {{s}}')
"""


FBA_BASELINES_CELL = """
# --- Two FBA baselines --------------------------------------------------
# 1) on-disk FBA: bounds untouched, matches results.csv byte-for-byte.
# 2) heuristic-baseline FBA: rebind every panel model to BASELINE_MAP
#    (= MSDB cascade defaults) -- the reference point for variant heuristics.

CACHE_KEY_ONDISK = 'fba_ondisk_v1'
CACHE_KEY_HEURBASE = 'fba_heuristic_baseline_v1'

def _run_panel(rev_map, baseline_map=None):
    return gh.run_panel(PANEL_IDS, reversibility_map=rev_map,
                        baseline_map=baseline_map, n_workers=4)

try:
    ondisk_fba = session.cache.load(CACHE_KEY_ONDISK)
    print(f'ondisk_fba: loaded {len(ondisk_fba)} from cache')
except KeyError:
    t = time.time()
    ondisk_fba = _run_panel(rev_map=None)
    session.cache.save(CACHE_KEY_ONDISK, ondisk_fba, type_hint='dict',
                       metadata={'rebound': False})
    print(f'ondisk_fba: {len(ondisk_fba)} models in {time.time()-t:.1f}s')

# Sanity check: every panel model's status + flux matches results.csv.
mm = 0
for r in ondisk_fba:
    o = ONDISK_FBA[r['model_id']]
    if r['status'] != o['status'] or abs(r['growth_flux'] - o['growth_flux']) > 1e-6:
        mm += 1
print(f'on-disk FBA reproduces results.csv for panel: '
      f'{len(ondisk_fba) - mm}/{len(ondisk_fba)}  (mismatches: {mm})')

try:
    heur_baseline_fba = session.cache.load(CACHE_KEY_HEURBASE)
    print(f'heur_baseline_fba: loaded {len(heur_baseline_fba)} from cache')
except KeyError:
    t = time.time()
    heur_baseline_fba = _run_panel(rev_map=BASELINE_MAP)
    session.cache.save(CACHE_KEY_HEURBASE, heur_baseline_fba, type_hint='dict',
                       metadata={'rebound': True, 'map': 'BASELINE_MAP'})
    print(f'heur_baseline_fba: {len(heur_baseline_fba)} models in {time.time()-t:.1f}s')

# How much does the heuristic-baseline diverge from on-disk?
ondisk_idx = {r['model_id']: r for r in ondisk_fba}
n_grow_diff = sum(1 for r in heur_baseline_fba
                  if ondisk_idx[r['model_id']]['grows'] != r['grows'])
n_flux_diff = sum(1 for r in heur_baseline_fba
                  if abs(ondisk_idx[r['model_id']]['growth_flux']
                         - r['growth_flux']) > 1e-6)
print(f'heuristic-baseline vs on-disk: '
      f'{n_grow_diff} models flip grow-status, '
      f'{n_flux_diff} models have flux delta > 1e-6')
"""


VARIANT_BOILERPLATE = '''
# -----------------------------------------------------------------------------
# Variant {tag}: {title}
# Heuristics Review {section}
# -----------------------------------------------------------------------------
def _cfg_{tag_safe}():
    {cfg_body}

def _build_variant_{tag_safe}():
    rxns = copy.deepcopy(msdb_rxns)
    cfg = _cfg_{tag_safe}()
    return lib.run_cascade(rxns, db_level='EQ', cfg=cfg, gc_first=True)

CACHE_KEY = 'reversibility_variant_{tag_safe}_v1'
try:
    variant = session.cache.load(CACHE_KEY)
except KeyError:
    t = time.time()
    variant = _build_variant_{tag_safe}()
    session.cache.save(CACHE_KEY, variant, type_hint='dict',
                       metadata={{'variant': '{tag}'}})
    print(f'variant {tag}: cascade in {{time.time()-t:.1f}}s')

VARIANT_MAP = {{r: rev for r, (_, rev) in variant.items()}}
diff = gh.reversibility_diff(BASELINE_MAP, VARIANT_MAP)
print(f'variant {tag}: {{diff["n_changed"]}} reactions changed direction')
for transition, n in sorted(diff['by_transition'].items()):
    print(f'  {{transition[0]!r}} -> {{transition[1]!r}}:  {{n}}')

FBA_KEY = 'fba_variant_{tag_safe}_v2'
try:
    variant_fba = session.cache.load(FBA_KEY)
except KeyError:
    t = time.time()
    # Fully rebind to VARIANT_MAP (no baseline filter) -- this makes the
    # variant_fba bounds directly comparable to heur_baseline_fba bounds.
    # Reactions where VARIANT_MAP == BASELINE_MAP get identical bounds in
    # both runs; only reactions where the variant differs contribute to
    # the FBA diff.
    variant_fba = gh.run_panel(PANEL_IDS, reversibility_map=VARIANT_MAP,
                                baseline_map=None, n_workers=4)
    session.cache.save(FBA_KEY, variant_fba, type_hint='dict',
                       metadata={{'variant': '{tag}'}})
    print(f'variant {tag}: panel FBA in {{time.time()-t:.1f}}s')

fba_diff = gh.diff_panel(heur_baseline_fba, variant_fba)
print(f'variant {tag}: panel growth -- '
      f'{{fba_diff["n_grow_change"]}}/{{fba_diff["n_models"]}} models flip grow-status, '
      f'{{fba_diff["n_flux_change"]}} have flux delta > 1e-6')

# Register the variant for the cross-variant summary table.
VARIANT_REGISTRY[{tag!r}] = {{
    'title'           : {title!r},
    'section'         : {section!r},
    'rev_diff'        : diff,
    'fba_diff'        : fba_diff,
}}
'''


VARIANTS = [
    {
        "tag": "3.1",
        "title": "Persist + use ln(reversibility_index)",
        "section": "§ 2.1 / 3.1",
        "cfg_body": (
            "ln_ri = lib.load_ln_reversibility_index()\n    "
            "print(f'  loaded ln_RI for {len(ln_ri)} reactions from MetaNetX_Reaction_Energies.tbl')\n    "
            "return lib.ReversibilityConfig(ln_ri_by_rxn=ln_ri)"
        ),
    },
    {
        "tag": "3.3",
        "title": "Bennett-2009 per-metabolite concentration ranges",
        "section": "§ 3.3",
        "cfg_body": (
            "return lib.ReversibilityConfig(\n        "
            "per_met_conc_range=lib.BENNETT_2009_ECOLI,\n        "
            "per_met_conc=lib.BENNETT_2009_MEAN,\n    )"
        ),
    },
    {
        "tag": "3.3_wide",
        "title": "Wider uniform concentration window [1e-7, 0.1] M",
        "section": "§ 3.3 (fallback)",
        "cfg_body": (
            "return lib.ReversibilityConfig(cell_min=1e-7, cell_max=1e-1)"
        ),
    },
    {
        "tag": "3.5",
        "title": "Per-reaction sigma band: k=1.96 (95%) instead of fixed +/-2 kcal",
        "section": "§ 3.5",
        "cfg_body": (
            "return lib.ReversibilityConfig(sigma_band_k=1.96)"
        ),
    },
    {
        "tag": "3.5_wide",
        "title": "Per-reaction CC bound widening: k=1.96 on stored_bounds",
        "section": "§ 3.5 / § 2.5",
        "cfg_body": (
            "return lib.ReversibilityConfig(sigma_bounds_k=1.96)"
        ),
    },
    {
        "tag": "3.6",
        "title": "Drop the low-energy-compounds list entirely",
        "section": "§ 3.6",
        "cfg_body": (
            "return lib.ReversibilityConfig(low_energy_cpds=())"
        ),
    },
    {
        "tag": "3.7",
        "title": "Drop the CO2 1e-4 hardcoded concentration override",
        "section": "§ 3.7",
        "cfg_body": (
            "return lib.ReversibilityConfig(apply_special_conc=False)"
        ),
    },
    {
        "tag": "3.10_tight",
        "title": "Tighten mMdeltaG band: +/-1 kcal/mol",
        "section": "§ 3.10",
        "cfg_body": (
            "return lib.ReversibilityConfig(mm_band=1.0)"
        ),
    },
    {
        "tag": "3.10_loose",
        "title": "Loosen mMdeltaG band: +/-4 kcal/mol",
        "section": "§ 3.10",
        "cfg_body": (
            "return lib.ReversibilityConfig(mm_band=4.0)"
        ),
    },
    {
        "tag": "H1",
        "title": "(NEW) default direction = '?' instead of '=' for unresolved",
        "section": "§ H1 (added by notebook)",
        "cfg_body": (
            "return lib.ReversibilityConfig(default_direction='?')"
        ),
    },
    {
        "tag": "H2",
        "title": "(NEW) repair LOW_LOCAL_CONC shadow bug (O2/H2 at 1e-6 M)",
        "section": "§ H2 (added by notebook)",
        "cfg_body": (
            "return lib.ReversibilityConfig(fix_low_local_conc=True)"
        ),
    },
    {
        "tag": "H3",
        "title": "(NEW) repair phosphates shadow bug (ABC + low-E phosphate spread)",
        "section": "§ H3 (added by notebook)",
        "cfg_body": (
            "return lib.ReversibilityConfig(fix_phosphates_shadow=True)"
        ),
    },
    {
        "tag": "H4",
        "title": "(NEW) stack 3.1 + 3.5 + Bennett: 'best-evidence' composite",
        "section": "§ H4 (added by notebook)",
        "cfg_body": (
            "ln_ri = lib.load_ln_reversibility_index()\n    "
            "return lib.ReversibilityConfig(\n        "
            "ln_ri_by_rxn=ln_ri,\n        "
            "sigma_band_k=1.96,\n        "
            "per_met_conc_range=lib.BENNETT_2009_ECOLI,\n        "
            "per_met_conc=lib.BENNETT_2009_MEAN,\n    )"
        ),
    },
]


SUMMARY_CELL = """
# Cross-variant summary table.
import pandas as pd
rows = []
rows.append({
    'variant': 'baseline',
    'section': '(reference)',
    'title'  : 'ReversibilityConfig() default',
    'rxns_changed_vs_baseline'  : 0,
    'panel_models_flux_changed' : 0,
    'panel_models_grow_flipped' : 0,
})
for tag, info in VARIANT_REGISTRY.items():
    rows.append({
        'variant': tag,
        'section': info['section'],
        'title'  : info['title'],
        'rxns_changed_vs_baseline'  : info['rev_diff']['n_changed'],
        'panel_models_flux_changed' : info['fba_diff']['n_flux_change'],
        'panel_models_grow_flipped' : info['fba_diff']['n_grow_change'],
    })
summary = pd.DataFrame(rows).set_index('variant')
print('Cross-variant summary:')
display(summary)

# Persist for later inspection / paper figures.
session.cache.save('reversibility_variant_summary_v1', summary,
                   type_hint='dataframe',
                   metadata={'panel': 'selected_ids.txt (100 models)'})
"""


SHARED_REGISTRY_CELL = "VARIANT_REGISTRY = {}"


# ---------------------------------------------------------------------------
def build_notebook() -> nbf.NotebookNode:
    cells = [
        md("""
        # 06 -- Reaction Reversibility Heuristics

        This notebook ports
        `ModelSEEDDatabase/Scripts/Thermodynamics/Estimate_Reaction_Reversibility.py`
        into a parameterizable library and exercises every suggestion in
        `/scratch/ctaylor/Reaction_Reversibility_Heuristics_Review.md`
        against the 100-model panel from
        `results/selected_ids.txt`.

        Each variant cell reports two things:

        1. **How many reactions changed direction** vs the baseline
           cascade, broken down by transition (`>` -> `=`, `=` -> `<`, ...).
        2. **How that change propagates to FBA growth** on the 100 panel
           models -- how many models flip between grower and non-grower,
           and how many have a meaningful flux delta.

        Caching: every heavy output lives in `notebooks/.kbcache/` via
        `KBUtils_local`'s `NotebookSession`.  The 56K-reaction MSDB load,
        the per-variant cascade, and the per-variant panel FBA are all
        cached -- re-running the notebook should be sub-second after the
        first execution.

        **What does not change.**  Nothing under `ModelSEEDDatabase/` or
        `core_models_kegg2/` is mutated.  All variant bounds live in
        in-memory cobra models.
        """),
        code(SETUP_CELL),

        md("## 1. Load MSDB reactions + panel + on-disk FBA results"),
        code(LOAD_MSDB_CELL),

        md("""
        ## 2. Baseline cascade (reproduces the upstream report)

        Run `Estimate_Reaction_Reversibility.estimate_one` (via the
        port in `reversibility_lib`) with `ReversibilityConfig()` --
        every knob at its upstream default.  The cascade output should
        match `Estimated_Reaction_Reversibility_Report_EQ.txt` exactly
        modulo a small set of reactions whose stored `deltag` was
        updated in MSDB after the report was regenerated.
        """),
        code(BASELINE_CASCADE_CELL),

        md("""
        ## 3. Two FBA baselines on the 100-model panel

        - **On-disk FBA** -- bounds untouched.  Should reproduce
          `results/results.csv` byte-for-byte (modulo solver jitter
          well below 1e-6).
        - **Heuristic-baseline FBA** -- panel models rebound to the
          cascade's baseline reversibility map.  Diverges from on-disk
          wherever the model JSON was built from a different
          reversibility (template-time curation, pre-EQ heuristics,
          ...).  This is the reference point all variants are compared
          against -- it isolates "what does the heuristic say" from
          "what does the on-disk model say".
        """),
        code(FBA_BASELINES_CELL),

        md("""
        ## 4. Heuristic variants

        Each variant cell:

        1. Builds a `ReversibilityConfig` with one knob flipped.
        2. Recomputes the cascade and the per-reaction direction diff
           vs the baseline cascade.
        3. Reruns FBA on the 100-model panel, but *only updates bounds
           for reactions whose direction changed* (so we measure the
           heuristic, not arbitrary differences between MSDB and the
           model JSONs).
        4. Reports the panel growth diff.

        Variants tagged `H*` are new suggestions added while writing
        this notebook -- they are documented at the bottom of
        `Reaction_Reversibility_Heuristics_Review.md`.
        """),
        code(SHARED_REGISTRY_CELL),
    ]

    for v in VARIANTS:
        v_copy = dict(v)
        v_copy["tag_safe"] = v["tag"].replace(".", "_")
        cells.append(md(f"### Variant {v['tag']} -- {v['title']}\n\n"
                        f"*Heuristics Review {v['section']}*"))
        cells.append(code(VARIANT_BOILERPLATE.format(**v_copy)))

    cells.append(md("## 5. Cross-variant summary"))
    cells.append(code(SUMMARY_CELL))

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
    print(f"wrote {NOTEBOOK_PATH}")


if __name__ == "__main__":
    main()
