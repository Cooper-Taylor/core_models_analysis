#!/usr/bin/env python3
"""
Build (or rebuild) the Jupyter notebooks under notebooks/.

Each notebook follows the same shape:
  1. Setup: import KBUtils, open a NotebookSession (cache lives next to the
     notebook in notebooks/.kbcache/).
  2. Render the matching reports/*.md verbatim as a markdown cell.
  3. Re-run the analysis interactively. Heavy outputs are cached via
     session.cache.save/.load so subsequent runs skip the recompute.
  4. Display the key tables / scalars produced by the analysis.

Run from anywhere:
    python3 scripts/build_notebooks.py
"""

from __future__ import annotations
import os

import json
import textwrap
from pathlib import Path

import nbformat as nbf

ROOT = Path(__file__).resolve().parent.parent
NOTEBOOK_DIR = ROOT / "notebooks"
REPORTS_DIR = ROOT / "reports"
NOTEBOOK_DIR.mkdir(exist_ok=True)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def md(text: str) -> nbf.NotebookNode:
    return nbf.v4.new_markdown_cell(textwrap.dedent(text).strip("\n"))


def code(src: str) -> nbf.NotebookNode:
    return nbf.v4.new_code_cell(textwrap.dedent(src).strip("\n"))


def render_report(path: Path) -> nbf.NotebookNode:
    """Embed an existing reports/*.md file verbatim as a markdown cell."""
    if not path.exists():
        return md(f"_Missing report file: `{path.name}` — run the script first._")
    body = path.read_text()
    return md(body)


SETUP_CELL = """
from pathlib import Path
import sys

# Find the project root (parent of notebooks/) so absolute paths work
# from anywhere the notebook is launched.
PROJECT_ROOT = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR", "/scratch/ctaylor/core_models_analysis"))
REPORTS = PROJECT_ROOT / 'reports'
RESULTS = PROJECT_ROOT / 'results'
LOGS    = PROJECT_ROOT / 'logs'
DATA    = PROJECT_ROOT / 'data'
SCRIPTS = PROJECT_ROOT / 'scripts'
sys.path.insert(0, str(SCRIPTS))

# KBUtils NotebookSession: opens .kbcache/ alongside this notebook so
# heavy outputs can be saved/loaded with provenance.
from kbutillib.notebook import NotebookSession
session = NotebookSession.for_notebook(notebook_file=__file__ if '__file__' in dir() else None,
                                       project_name='core_models_analysis')
print('Cache directory:', session.kbcache_dir)
print('Notebook name :', session.notebook_name)
"""


def new_notebook(title_md: str, cells: list) -> nbf.NotebookNode:
    nb = nbf.v4.new_notebook()
    nb.cells = [md(title_md), code(SETUP_CELL)] + cells
    nb.metadata = {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python"},
    }
    return nb


def write(nb: nbf.NotebookNode, name: str) -> None:
    path = NOTEBOOK_DIR / name
    nbf.write(nb, path)
    print(f"wrote {path}")


# ---------------------------------------------------------------------------
# 00 — Index
# ---------------------------------------------------------------------------
def build_index() -> nbf.NotebookNode:
    body_cells = [
        md("""
        # Core Models KEGG2 — Analysis Index

        This project asks: of the 5,683 ModelSEED metabolic models in
        `data/core_models_kegg2/`, which can produce biomass under complete
        media, why do the others fail, and which subset of growers best
        represents the diversity of the whole?

        The work is split into five notebooks. Each notebook re-runs its
        stage (cached via KBUtils `NotebookSession`) and renders the matching
        `reports/*.md` writeup verbatim, so this index is the only file you
        need to start from.
        """),
        md("""
        ## Notebook map

        | # | Notebook | Stage | Report |
        |---|---|---|---|
        | 01 | `01_GrowthFBA_Pipeline.ipynb` | FBA biomass solve over 5,683 models on complete media | `SUMMARY.md` |
        | 02 | `02_CharacteristicsAnalysis.ipynb` | Grower vs non-grower model-size + flux-distribution comparison | `CHARACTERISTICS.md` |
        | 03 | `03_GapAnalysis.ipynb` | Per-non-grower biomass-precursor reachability + the annotated interpretation | `GAP_ANALYSIS.md` + `INTERPRETATION.md` |
        | 04 | `04_ReactionPrevalence.ipynb` | Reactions enriched in growers vs non-growers | `REACTION_PREVALENCE.md` |
        | 05 | `05_DiversePanelSelection.ipynb` | 100-model panel that spans the 3,461 growers (descriptive test set) | `DIVERSE_SELECTION.md` |
        """),
        md("""
        ## Project layout

        ```
        core_models_analysis/
        ├── data/      core_models_kegg2/ → symlink to the model JSONs (5,683 files)
        ├── scripts/   analyze_growth.py · deeper_analysis.py · select_diverse.py · annotate.py · summarize.py · build_notebooks.py
        ├── notebooks/ this index + 5 stage notebooks · .kbcache/ (provenanced cache)
        ├── reports/   markdown writeups, rendered inside the notebooks
        ├── results/   results.csv · growers.csv · non_growers.csv · gap_per_model.json · selected_models.{csv,json} · selected_ids.txt
        └── logs/      run.log · deeper.log · failures.log
        ```

        **Quick access to the diverse-panel IDs (descriptive test set):**
        ```python
        from pathlib import Path
        ids = Path(os.environ.get('CORE_MODELS_ANALYSIS_DIR', '/scratch/ctaylor/core_models_analysis'), 'results', 'selected_ids.txt').read_text().split()
        # → 100 model IDs spanning the 3,461 growers
        ```
        """),
        code("""
        # Confirm everything resolves
        from pathlib import Path
        ROOT = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR", "/scratch/ctaylor/core_models_analysis"))
        for sub in ('data', 'scripts', 'notebooks', 'reports', 'results', 'logs'):
            d = ROOT / sub
            n = sum(1 for _ in d.iterdir()) if d.exists() else 0
            print(f'{sub:9s} {"✓" if d.exists() else "✗"}  {n:4d} entries')
        """),
        code("""
        # Show the 5 reports at a glance
        for r in sorted((ROOT / 'reports').glob('*.md')):
            print(f'• {r.name:25s} {r.stat().st_size:7d} bytes')
        """),
    ]
    nb = new_notebook("# 00 — Index", body_cells)
    return nb


# ---------------------------------------------------------------------------
# 01 — Growth FBA Pipeline
# ---------------------------------------------------------------------------
def build_01_growth_pipeline() -> nbf.NotebookNode:
    cells = [
        md("""
        ## What this notebook does

        Loads every model in `data/core_models_kegg2/`, restricts exchange
        reactions to the ModelSEED `KBaseMedia.cpd` complete media (347
        compounds), and runs FBA on `bio1`. A model "grows" if the optimum
        biomass flux exceeds 1e-6.

        The full multiprocess sweep is implemented in
        `scripts/analyze_growth.py`. Here we read the cached `results.csv`
        and surface the same numbers `summarize.py` writes into
        `SUMMARY.md`.
        """),
        code("""
        import csv, pandas as pd
        from collections import Counter
        from statistics import mean, median

        df = pd.read_csv(RESULTS / 'results.csv')
        df['grows'] = df['grows'].astype(str).str.lower() == 'true'
        print(f'{len(df)} models')
        df.head()
        """),
        code("""
        # Outcome breakdown — same shape as SUMMARY.md
        status_counts = Counter(df['status'])
        n_grow = int(df['grows'].sum())
        n_zero = int(((df['status'] == 'optimal') & (~df['grows'])).sum())
        n_other = len(df) - n_grow - n_zero
        print(f'  grows (flux > 1e-6):           {n_grow:5d}   ({100*n_grow/len(df):.1f}%)')
        print(f'  optimal but zero biomass:      {n_zero:5d}   ({100*n_zero/len(df):.1f}%)')
        print(f'  other (infeasible/error/none): {n_other:5d}')
        print('  solver-status histogram:', dict(status_counts))
        """),
        code("""
        # Cache the grower frame via KBUtils NotebookSession
        growers = df[df['grows']].reset_index(drop=True)
        session.cache.save('growers_frame', growers, type_hint='dataframe',
                           metadata={'n_growers': len(growers)})
        print(f'cached {len(growers)} grower rows under name "growers_frame"')
        print(session.cache.load('growers_frame').describe(include='all').iloc[:5])
        """),
        code("""
        # Distributions: biomass flux + reaction counts
        flux = growers['growth_flux']
        print(f'biomass flux  — min {flux.min():.2f}  median {flux.median():.2f}  '
              f'mean {flux.mean():.2f}  max {flux.max():.2f}')
        print(f'n_reactions   — min {df["n_reactions"].min()}  median {df["n_reactions"].median():.0f}  max {df["n_reactions"].max()}')
        """),
        md("---\n## Report: `reports/SUMMARY.md`"),
        render_report(REPORTS_DIR / "SUMMARY.md"),
    ]
    return new_notebook("# 01 — Growth FBA Pipeline", cells)


# ---------------------------------------------------------------------------
# 02 — Characteristics
# ---------------------------------------------------------------------------
def build_02_characteristics() -> nbf.NotebookNode:
    cells = [
        md("""
        ## What this notebook does

        Compares grower vs non-grower distributions across the size metrics
        (metabolites, reactions, genes, exchanges) and the biomass-flux
        bucket histogram, mirroring `scripts/deeper_analysis.py` part 1.
        """),
        code("""
        import pandas as pd
        from collections import Counter
        df = pd.read_csv(RESULTS / 'results.csv')
        df['grows'] = df['grows'].astype(str).str.lower() == 'true'
        keys = ['n_metabolites', 'n_reactions', 'n_genes',
                'n_exchanges_total', 'n_exchanges_open']
        tbl = df.groupby('grows')[keys].agg(['median', 'mean']).T
        tbl.columns = ['non-grower', 'grower']
        tbl
        """),
        code("""
        # Cache the characteristics table
        session.cache.save('characteristics_table', tbl.reset_index(),
                           type_hint='dataframe',
                           metadata={'description': 'grower vs non-grower size stats'})
        print('cached as characteristics_table')
        """),
        code("""
        # Flux-bucket histogram
        bins = [0, 1, 5, 10, 25, 50, 75, 100]
        growers = df[df['grows']]
        hist = []
        prev = 0
        for hi in bins[1:]:
            n = int(((growers['growth_flux'] > prev) & (growers['growth_flux'] <= hi)).sum())
            hist.append((f'({prev}, {hi}]', n))
            prev = hi
        pd.DataFrame(hist, columns=['flux range', 'n_growers'])
        """),
        md("---\n## Report: `reports/CHARACTERISTICS.md`"),
        render_report(REPORTS_DIR / "CHARACTERISTICS.md"),
    ]
    return new_notebook("# 02 — Characteristics: Grower vs Non-grower", cells)


# ---------------------------------------------------------------------------
# 03 — Gap analysis
# ---------------------------------------------------------------------------
def build_03_gap_analysis() -> nbf.NotebookNode:
    cells = [
        md("""
        ## What this notebook does

        Loads the precomputed per-non-grower gap analysis from
        `results/gap_per_model.json` (produced by
        `scripts/deeper_analysis.py` part 2 — for each non-grower it probes
        every biomass precursor with a demand reaction under complete
        media). Tallies the most-frequently-blocked precursors and shows
        the distribution of "how many precursors are blocked per model".
        """),
        code("""
        import json, pandas as pd
        from collections import Counter
        gap = json.loads((RESULTS / 'gap_per_model.json').read_text())
        print(f'{len(gap)} non-growers analysed')

        blocked = Counter()
        per_model_blocked = []
        errors = 0
        for mid, rec in gap.items():
            if 'error' in rec and rec.get('error'):
                errors += 1
                continue
            for met in rec.get('blocked', []):
                blocked[met.split('_')[0]] += 1
            per_model_blocked.append(len(rec.get('blocked', [])))
        print(f'  errors: {errors}')
        top = pd.DataFrame(blocked.most_common(15), columns=['cpd', 'n_blocked'])
        top['pct'] = (100 * top['n_blocked'] / len(gap)).round(1)
        top
        """),
        code("""
        # Cache the blocked-precursor counter
        session.cache.save('blocked_precursor_counts', dict(blocked),
                           type_hint='json',
                           metadata={'n_nongrowers': len(gap)})
        # Distribution: how many precursors are blocked per non-grower
        s = pd.Series(per_model_blocked)
        print(f'precursors blocked per non-grower — min {s.min()}  '
              f'median {int(s.median())}  mean {s.mean():.2f}  max {s.max()}')
        """),
        md("---\n## Report: `reports/GAP_ANALYSIS.md`"),
        render_report(REPORTS_DIR / "GAP_ANALYSIS.md"),
        md("---\n## Report: `reports/INTERPRETATION.md`  (annotated reading)"),
        render_report(REPORTS_DIR / "INTERPRETATION.md"),
    ]
    return new_notebook("# 03 — Gap Analysis: What Blocks Non-growers", cells)


# ---------------------------------------------------------------------------
# 04 — Reaction prevalence
# ---------------------------------------------------------------------------
def build_04_reaction_prevalence() -> nbf.NotebookNode:
    cells = [
        md("""
        ## What this notebook does

        Reads every model JSON, extracts the set of `seed.reaction`
        annotations, and computes the prevalence delta between growers and
        non-growers — mirroring `scripts/deeper_analysis.py` part 3. The
        per-model reaction set is cached so the diversity-selection notebook
        can reuse it without re-reading 5,683 files.
        """),
        code("""
        import json, pandas as pd, multiprocessing as mp
        from collections import Counter
        from pathlib import Path

        df = pd.read_csv(RESULTS / 'results.csv')
        df['grows'] = df['grows'].astype(str).str.lower() == 'true'
        MODELS = DATA / 'core_models_kegg2'

        def extract(mid):
            with open(MODELS / f'{mid}.json') as f:
                m = json.load(f)
            return mid, {r['annotation'].get('seed.reaction')
                         for r in m['reactions']
                         if r.get('annotation', {}).get('seed.reaction')}

        # Reuse cache if we ran this before
        try:
            cached = session.cache.load('rxnsets_by_model')
            print(f'loaded {len(cached)} reaction sets from cache')
            rxnsets = {mid: set(rxns) for mid, rxns in cached.items()}
        except KeyError:
            print('cache miss — extracting from JSON (may take ~30s)…')
            with mp.Pool(16) as pool:
                rxnsets = dict(pool.imap_unordered(extract, df['model_id'], chunksize=8))
            # save as plain dict of lists for JSON serializer
            session.cache.save('rxnsets_by_model',
                               {mid: sorted(s) for mid, s in rxnsets.items()},
                               type_hint='json',
                               metadata={'n_models': len(rxnsets)})
            print(f'cached {len(rxnsets)} reaction sets')
        """),
        code("""
        growers = set(df.loc[df['grows'], 'model_id'])
        nongrow = set(df.loc[~df['grows'], 'model_id'])
        gC, nC = Counter(), Counter()
        for mid, rs in rxnsets.items():
            tgt = gC if mid in growers else nC
            tgt.update(rs)

        all_r = set(gC) | set(nC)
        rows = []
        for r in all_r:
            gf, nf = gC[r] / len(growers), nC[r] / len(nongrow)
            rows.append((r, gC[r], 100*gf, nC[r], 100*nf, 100*(gf - nf)))
        prev_df = pd.DataFrame(rows, columns=['rxn', 'g_count', 'g_pct',
                                              'n_count', 'n_pct', 'delta_pct'])
        prev_df.sort_values('delta_pct', ascending=False).head(15)
        """),
        code("""
        prev_df.sort_values('delta_pct').head(15)   # most enriched in NON-growers
        """),
        md("---\n## Report: `reports/REACTION_PREVALENCE.md`"),
        render_report(REPORTS_DIR / "REACTION_PREVALENCE.md"),
    ]
    return new_notebook("# 04 — Reaction Prevalence: Growers vs Non-growers", cells)


# ---------------------------------------------------------------------------
# 05 — Diverse panel selection
# ---------------------------------------------------------------------------
def build_05_diverse_selection() -> nbf.NotebookNode:
    cells = [
        md("""
        ## What this notebook does

        Picks **100 models** from the 3,461 growers that together span the
        diversity of growth-capable networks. Selection uses four passes
        (greedy reaction coverage → metabolite coverage → forced extremes
        across multiple axes → farthest-point Jaccard sampling) so each
        pick carries an explicit reason. See `scripts/select_diverse.py`
        for the implementation.

        The notebook loads the cached results, validates panel
        representativeness, and exposes the 100 IDs for downstream
        descriptive testing of growing models.
        """),
        code("""
        import json, pandas as pd
        sel = json.loads((RESULTS / 'selected_models.json').read_text())
        print(f'panel size: {len(sel["ids"])}')
        print(f'reactions covered:    {sel["coverage"]["reactions_covered"]} / {sel["coverage"]["reactions_universe"]}')
        print(f'metabolites covered:  {sel["coverage"]["metabolites_covered"]} / {sel["coverage"]["metabolites_universe"]}')
        """),
        code("""
        # Selection-reason breakdown
        from collections import Counter
        kind = Counter()
        for r in sel['records']:
            tag = r['reason'].split(':')[0]
            kind[tag] += 1
        for k, v in kind.most_common():
            print(f'{v:3d}  {k}')
        """),
        code("""
        # Per-model summary (first 25)
        panel = pd.read_csv(RESULTS / 'selected_models.csv')
        panel.head(25)
        """),
        code("""
        # Cache the panel as a typed object for downstream notebooks
        session.cache.save('descriptive_test_panel',
                           {'ids': sel['ids'],
                            'coverage': sel['coverage']},
                           type_hint='json',
                           metadata={'n_models': len(sel['ids']),
                                     'role': 'descriptive test set for growing models'})
        ids = sel['ids']
        print(f'{len(ids)} IDs available — first 10:')
        for i in ids[:10]:
            print(' ', i)
        """),
        code("""
        # Programmatic access for downstream scripts
        from pathlib import Path
        print('From a fresh script:')
        print('  ids = Path(\"/scratch/ctaylor/core_models_analysis/results/selected_ids.txt\").read_text().split()')
        print()
        print('Or, in a notebook that opens its own NotebookSession:')
        print('  panel = session.cache.load(\"descriptive_test_panel\")')
        print('  ids = panel[\"ids\"]')
        """),
        code("""
        # Quick validation: distance from non-panel growers to panel
        import json, statistics
        df = pd.read_csv(RESULTS / 'results.csv')
        df['grows'] = df['grows'].astype(str).str.lower() == 'true'

        try:
            cached = session.cache.load('rxnsets_by_model')
            fps = {mid: set(rxns) for mid, rxns in cached.items()
                   if mid in set(df.loc[df['grows'], 'model_id'])}
        except KeyError:
            print('rxnsets_by_model cache not present — run notebook 04 first.')
            fps = None

        if fps:
            sel_set = [fps[i] for i in sel['ids']]
            def jacc(a, b):
                if not a and not b: return 0.0
                return 1.0 - len(a & b) / len(a | b)
            dists = [min(jacc(fps[m], s) for s in sel_set)
                     for m in fps if m not in set(sel['ids'])]
            print(f'non-panel growers: {len(dists)}')
            print(f'  median min-dist to panel: {statistics.median(dists):.3f}')
            print(f'  mean: {statistics.mean(dists):.3f}')
            print(f'  max (worst-represented):  {max(dists):.3f}')
        """),
        md("---\n## Report: `reports/DIVERSE_SELECTION.md`"),
        render_report(REPORTS_DIR / "DIVERSE_SELECTION.md"),
    ]
    return new_notebook("# 05 — Diverse Panel Selection (Descriptive Test Set)", cells)


# ---------------------------------------------------------------------------
# Drive
# ---------------------------------------------------------------------------
def main():
    write(build_index(), "00_Index.ipynb")
    write(build_01_growth_pipeline(), "01_GrowthFBA_Pipeline.ipynb")
    write(build_02_characteristics(), "02_CharacteristicsAnalysis.ipynb")
    write(build_03_gap_analysis(), "03_GapAnalysis.ipynb")
    write(build_04_reaction_prevalence(), "04_ReactionPrevalence.ipynb")
    write(build_05_diverse_selection(), "05_DiversePanelSelection.ipynb")


if __name__ == "__main__":
    main()
