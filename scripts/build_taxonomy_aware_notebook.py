#!/usr/bin/env python3
"""
Build notebook 08 — Taxonomy-Aware Diverse Panel Selection.

Mirrors the cell-construction idioms of ``scripts/build_notebooks.py``
(md / code / new_notebook helpers, KBUtils setup cell) but writes ONLY
``notebooks/08_TaxonomyAwareSelection.ipynb`` so re-running this script
does not touch notebooks 00–07.

The notebook walks the new 6-pass selection step-by-step so the
reviewer can see exactly what each pass contributes, then compares the
new panel to the original ``results/selected_*.json`` panel.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import nbformat as nbf

ROOT = Path(__file__).resolve().parent.parent
NOTEBOOK_DIR = ROOT / 'notebooks'
REPORTS_DIR = ROOT / 'reports'


def md(text: str) -> nbf.NotebookNode:
    return nbf.v4.new_markdown_cell(textwrap.dedent(text).strip('\n'))


def code(src: str) -> nbf.NotebookNode:
    return nbf.v4.new_code_cell(textwrap.dedent(src).strip('\n'))


def render_report(path: Path) -> nbf.NotebookNode:
    if not path.exists():
        return md(f'_Missing report file: `{path.name}` — '
                  f'run cell that writes it (Section 11) first._')
    return md(path.read_text())


SETUP_CELL = """
from pathlib import Path
import sys

# Find the project root (parent of notebooks/) so absolute paths work
# from anywhere the notebook is launched.
PROJECT_ROOT = Path('/scratch/ctaylor/core_models_analysis')
REPORTS = PROJECT_ROOT / 'reports'
RESULTS = PROJECT_ROOT / 'results'
LOGS    = PROJECT_ROOT / 'logs'
DATA    = PROJECT_ROOT / 'data'
SCRIPTS = PROJECT_ROOT / 'scripts'
sys.path.insert(0, str(SCRIPTS))

# KBUtils NotebookSession: opens .kbcache/ alongside this notebook so
# heavy outputs (the 3,461 reaction fingerprints) can be saved/loaded
# with provenance.
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
        'kernelspec': {'display_name': 'Python 3',
                       'language': 'python', 'name': 'python3'},
        'language_info': {'name': 'python'},
    }
    return nb


def build_08_taxonomy_aware_selection() -> nbf.NotebookNode:
    cells = [
        md("""
        ## What this notebook does

        Picks a **new 100-model panel** from the 3,461 growers that
        deliberately spans **NCBI taxonomy** in addition to network
        structure.  The original panel (`scripts/select_diverse.py`,
        notebook `05`, `results/selected_*`) covered `seed.reaction`
        and `cpd` sets very well, but only hit **19 of ~33
        grower-bearing phyla** — missing *Myxococcota* (21 growers),
        *Sphingobacteriia* (20), the entire *Micrococcales* order (122
        growers) and big industrial / clinical genera like
        *Pseudomonas* (84 growers) and *Burkholderia* (81).

        The taxonomy data was added in notebook `07` and lives at
        `results/ncbi_taxonomy.csv` (3,379 / 3,461 growers resolved;
        the rest are bucketed as `Unknown`).  This notebook layers
        that information into a 7-pass selection algorithm:

        1. **Phylum medoids** — one anchor per grower-bearing phylum.
        2. **Reaction-coverage core** — greedy max-coverage on
           `seed.reaction` IDs (12 picks).
        3. **Taxonomic-novelty fill** — greedy farthest-point on
           lineage rank-distance (12 picks, closes class/order gaps).
        4. **Metabolite-coverage layer** — greedy max-coverage on
           cpd IDs (8 picks).
        5. **Constrained extremes** — original axis extremes with
           a class-saturation guard (8 picks).
        6. **Hot-taxon medoids** — picks medoids for the top missing
           grower-heavy orders / genera / families (14 picks; closes
           broad orders like *Micrococcales* (122 growers),
           *Alteromonadales* (102), *Rhodobacterales* (95) that
           Pass 3 misses because their constituent genera are
           individually small).
        7. **Farthest-point Jaccard fill** — fills the rest by
           maximizing min Jaccard reaction distance.

        Each pick is tagged with the pass that selected it and a
        per-pass reason string.  All new artifacts use the `_tax`
        suffix so the original panel is left untouched for
        comparison:

        - `results/selected_ids_tax.txt`
        - `results/selected_models_tax.csv`
        - `results/selected_models_tax.json`
        - `reports/TAXONOMY_AWARE_SELECTION.md`

        Implementation lives in `scripts/select_diverse_tax.py`;
        this notebook calls into it and renders intermediate state.
        """),

        # ------------------------------------------------------------------
        md('## 1. Load inputs'),
        md("""
        We need three things:

        - **growers** — the 3,461 rows from `results.csv` where `grows=True`,
        - **reaction & metabolite fingerprints** — `seed.reaction` and
          `cpd` sets per model (cached as `rxnsets_by_model` /
          `cpdsets_by_model` to skip the JSON re-parse on subsequent
          runs),
        - **NCBI taxonomy** — `results/ncbi_taxonomy.csv` from notebook 07.
        """),
        code("""
        import json, pandas as pd, multiprocessing as mp
        from collections import Counter, defaultdict
        from pathlib import Path

        import select_diverse_tax as sdt

        growers = sdt.load_growers()
        growers_by_id = {g['model_id']: g for g in growers}
        grower_ids = list(growers_by_id)
        print(f'{len(growers)} growers loaded')
        """),
        code("""
        # Reaction + metabolite fingerprints.
        # Notebook 04 caches 'rxnsets_by_model' (all 5,683 models).  We
        # also cache the cpd sets here so this notebook is self-contained.
        MODELS = DATA / 'core_models_kegg2'

        def _extract(mid):
            with open(MODELS / f'{mid}.json') as f:
                m = json.load(f)
            rxns = {r['annotation'].get('seed.reaction')
                    for r in m['reactions']
                    if r.get('annotation', {}).get('seed.reaction')}
            cpds = {met['id'].split('_')[0] for met in m['metabolites']}
            return mid, sorted(rxns), sorted(cpds)

        try:
            cached_rxns = session.cache.load('rxnsets_by_model')
            rxns_by_id = {mid: set(rxns) for mid, rxns in cached_rxns.items()
                          if mid in growers_by_id}
            print(f'rxnsets_by_model cache hit: {len(rxns_by_id)} growers')
        except KeyError:
            print('rxnsets_by_model cache miss — extracting reactions from JSON ...')
            rxns_by_id = {}
            with mp.Pool(16) as pool:
                for mid, rxns, _ in pool.imap_unordered(_extract, grower_ids, chunksize=8):
                    rxns_by_id[mid] = set(rxns)
            session.cache.save('rxnsets_by_model',
                               {mid: sorted(s) for mid, s in rxns_by_id.items()},
                               type_hint='json',
                               metadata={'n_models': len(rxns_by_id),
                                         'source': 'notebook 08 extraction'})

        try:
            cached_cpds = session.cache.load('cpdsets_by_model')
            cpds_by_id = {mid: set(cpds) for mid, cpds in cached_cpds.items()
                          if mid in growers_by_id}
            print(f'cpdsets_by_model cache hit: {len(cpds_by_id)} growers')
        except KeyError:
            print('cpdsets_by_model cache miss — extracting metabolites from JSON ...')
            cpds_by_id = {}
            with mp.Pool(16) as pool:
                for mid, _, cpds in pool.imap_unordered(_extract, grower_ids, chunksize=8):
                    cpds_by_id[mid] = set(cpds)
            session.cache.save('cpdsets_by_model',
                               {mid: sorted(s) for mid, s in cpds_by_id.items()},
                               type_hint='json',
                               metadata={'n_models': len(cpds_by_id),
                                         'source': 'notebook 08 extraction'})

        print(f'reactions per grower: median '
              f'{int(pd.Series([len(s) for s in rxns_by_id.values()]).median())}; '
              f'metabolites per grower: median '
              f'{int(pd.Series([len(s) for s in cpds_by_id.values()]).median())}')
        """),
        code("""
        # NCBI taxonomy: ranks superkingdom -> species, plus organism_name.
        # Unresolved rows ('missing' status) get all ranks = 'Unknown' so they
        # cluster into a single bucket in Pass 1 rather than scattering.
        lineage = sdt.load_taxonomy()

        tax_df = pd.read_csv(RESULTS / 'ncbi_taxonomy.csv').fillna('Unknown')
        tax_grower = tax_df[tax_df['assembly_accession'].isin(growers_by_id)].copy()
        n_resolved = (tax_grower['status'] == 'resolved').sum()
        print(f'growers with taxonomy row: {len(tax_grower)} / {len(growers)}')
        print(f'  resolved: {n_resolved}; bucketed Unknown: {len(growers) - n_resolved}')
        tax_grower['phylum'] = tax_grower['phylum'].replace('Unknown', 'Unknown')
        # Quick phylum tally for context (top 10).
        print('\\nGrower phyla (top 10):')
        for phy, n in tax_grower['phylum'].value_counts().head(10).items():
            print(f'  {phy:30s} {n:>5}')
        """),

        # ------------------------------------------------------------------
        md('## 2. Audit the *original* panel against the grower taxonomy'),
        md("""
        Before we re-select, quantify the gaps the new algorithm needs to
        close. We load `results/selected_models.json` (the original 100)
        and tally how often each phylum / class / order / genus appears
        in the original panel vs. the full grower set.
        """),
        code("""
        original = json.loads((RESULTS / 'selected_models.json').read_text())
        original_ids = set(original['ids'])
        orig_lin = pd.DataFrame([
            {'model_id': mid, **{r: lineage[mid][r] for r in sdt.LINEAGE_RANKS}}
            for mid in original['ids']
        ])
        print(f'original panel size: {len(original_ids)}')
        print(f'original panel with resolved taxonomy: '
              f'{(orig_lin["phylum"] != "Unknown").sum()}')

        def rank_coverage(panel_lin, grower_lin, rank):
            panel_set = set(panel_lin[rank]) - {'Unknown'}
            grower_set = set(grower_lin[rank]) - {'Unknown'}
            return len(panel_set), len(grower_set), grower_set - panel_set

        rows = []
        for r in sdt.LINEAGE_RANKS:
            panel_n, grower_n, missing = rank_coverage(orig_lin, tax_grower, r)
            rows.append((r, panel_n, grower_n,
                         f'{100*panel_n/grower_n:.1f}%' if grower_n else '—',
                         len(missing)))
        pd.DataFrame(rows, columns=['rank', 'panel_distinct',
                                    'grower_distinct', 'coverage_pct',
                                    'taxa_missing_from_panel'])
        """),
        code("""
        # Which grower-bearing phyla / classes / orders does the original
        # panel skip entirely? Sorted by grower count desc -> the biggest
        # gaps the new algorithm should close.
        def missing_table(rank, top=10):
            panel_set = set(orig_lin[rank]) - {'Unknown'}
            counts = tax_grower[rank].value_counts()
            counts = counts[~counts.index.isin(panel_set | {'Unknown'})]
            return counts.head(top).rename_axis(rank).reset_index(name='grower_count')

        print('PHYLA the original panel misses (top 10 by grower count):')
        print(missing_table('phylum'))
        print('\\nCLASSES the original panel misses (top 10):')
        print(missing_table('class'))
        print('\\nORDERS the original panel misses (top 10):')
        print(missing_table('order'))
        print('\\nGENERA the original panel misses (top 10):')
        print(missing_table('genus'))
        """),

        # ------------------------------------------------------------------
        md('## 3. Pass 1 — phylum medoids'),
        md("""
        For every grower-bearing phylum, pick the model whose reaction
        fingerprint is the most central (minimum sum-Jaccard distance to
        its phylum-mates).  Phyla with only 1–2 growers get their lone
        member; phyla with > 50 growers are subsampled to 50
        deterministically (stride sampling on sorted IDs) to keep the
        all-pairs Jaccard tractable.

        This guarantees that **every phylum present in the growers** has
        at least one panel anchor — closing the largest original gap.
        """),
        code("""
        p1 = sdt.pass1_phylum_medoids(grower_ids, rxns_by_id, lineage)
        print(f'Pass 1 picked {len(p1)} medoids (one per grower-bearing phylum)')
        # Show first 10
        print('\\n# | model_id        | phylum                  | sum_jaccard')
        for i, (mid, phy, sj) in enumerate(p1[:15], 1):
            print(f'{i:2d} | {mid:15s} | {phy:24s} | {sj}')
        """),

        # ------------------------------------------------------------------
        md('## 4. Pass 2 — reaction-coverage core'),
        md("""
        Standard greedy max-coverage on `seed.reaction` IDs, **seeded
        with the Pass 1 picks**.  Each new pick adds the most
        previously-uncovered reactions.  Capped at 12 — fewer than the
        original's 20 because Pass 1 has already pulled in a lot of
        reaction diversity through the phylum anchors.
        """),
        code("""
        picked_after_p1 = {mid for mid, _, _ in p1}
        p2, covered_after_p2 = sdt.greedy_max_coverage(
            grower_ids, rxns_by_id, sdt.PASS2_RXN_CORE, picked_after_p1)
        print(f'Pass 2 picked {len(p2)} models')
        for i, (mid, gain) in enumerate(p2, 1):
            lin = lineage[mid]
            print(f'{i:2d} | {mid:15s} | +{gain:>3d} novel rxns | '
                  f'{lin["phylum"]:18s} {lin["class"]:18s} {lin["genus"]}')
        """),

        # ------------------------------------------------------------------
        md('## 5. Pass 3 — taxonomic-novelty fill'),
        md("""
        Greedy **farthest-point on lineage rank-distance**:
        7 = different superkingdom, 6 = different phylum, …, 1 = different
        species, 0 = identical lineage.  Tie-break with Jaccard reaction
        distance.  This is the pass that explicitly fills holes at the
        class / order / genus levels — for example, Pass 1 anchors
        Bacteroidota but doesn't necessarily hit *Sphingobacteriia*; Pass
        3 will pick the most-distant unrepresented class first.
        """),
        code("""
        picked_after_p2 = picked_after_p1 | {mid for mid, _ in p2}
        p3 = sdt.pass3_taxonomic_novelty(grower_ids, rxns_by_id, lineage,
                                         sdt.PASS3_TAX_NOVELTY,
                                         picked_after_p2)
        print(f'Pass 3 picked {len(p3)} models')
        print('# | model_id        | rank-d | tiebreak Jaccard | lineage (phylum/class/order/genus)')
        for i, (mid, rd, tj) in enumerate(p3, 1):
            lin = lineage[mid]
            print(f'{i:2d} | {mid:15s} | {rd:>5d}  | {tj:>14.3f}   | '
                  f'{lin["phylum"]}/{lin["class"]}/{lin["order"]}/{lin["genus"]}')
        """),

        # ------------------------------------------------------------------
        md('## 6. Pass 4 — metabolite-coverage layer'),
        md("""
        Same greedy idea as Pass 2, but on `cpd` IDs and seeded with
        every pick so far.  Catches transporter / cofactor diversity
        that the reaction-set passes miss.
        """),
        code("""
        picked_after_p3 = picked_after_p2 | {mid for mid, _, _ in p3}
        p4, _ = sdt.greedy_max_coverage(grower_ids, cpds_by_id,
                                        sdt.PASS4_CPD_LAYER, picked_after_p3)
        print(f'Pass 4 picked {len(p4)} models')
        for i, (mid, gain) in enumerate(p4, 1):
            lin = lineage[mid]
            print(f'{i:2d} | {mid:15s} | +{gain:>3d} novel cpds | '
                  f'{lin["phylum"]:18s} {lin["genus"]}')
        """),

        # ------------------------------------------------------------------
        md('## 7. Pass 5 — constrained extremes'),
        md("""
        Original 12 axis extremes (min/max reactions, metabolites, genes,
        flux, exchanges, fingerprint sizes) + 3 rare-reaction champions,
        but with a **class over-representation guard**: if an axis
        winner's class already has > 2× its expected panel share, walk
        down the sorted ranking up to 20 candidates to find one whose
        class is not yet saturated.  Cap at 10 new picks.
        """),
        code("""
        # Build the universe-prevalence / rare set before calling Pass 5.
        prevalence = Counter()
        for s in rxns_by_id.values():
            prevalence.update(s)
        rare_set = {r for r, c in prevalence.items()
                    if c / len(rxns_by_id) < 0.01}
        print(f'rare seed.reactions (<1% of growers): {len(rare_set)}')

        picked_after_p4 = picked_after_p3 | {mid for mid, _ in p4}
        p5 = sdt.pass5_constrained_extremes(growers_by_id, rxns_by_id,
                                            cpds_by_id, rare_set, lineage,
                                            picked_after_p4,
                                            cap=sdt.PASS5_EXTREMES)
        print(f'\\nPass 5 picked {sum(1 for mid, _ in p5 if mid not in picked_after_p4)} new models '
              f'(some axes may resolve to a model already picked)')
        for mid, reason in p5:
            lin = lineage[mid]
            tag = '(NEW)' if mid not in picked_after_p4 else '(dup)'
            print(f'  {tag} {mid:15s} | {lin["class"]:22s} | {reason}')
        """),

        # ------------------------------------------------------------------
        md('## 8. Pass 6 — hot-taxon medoids'),
        md("""
        Closes within-phylum gaps that Pass 3 misses because it
        prioritizes higher-rank distance.  For each rank in
        `['order', 'genus', 'family']` (broadest first), enumerate
        grower-bearing taxa with ≥20 growers; for each taxon NOT yet
        in the panel, add its medoid.  Cap 14 picks.

        This is what brings in *Micrococcales* (122 growers),
        *Alteromonadales* (102), *Rhodobacterales* (95),
        *Vibrionales*, *Sphingomonadales* — broad orders whose
        phylum is anchored but whose specific order had no
        representative.  Many also bring in big industrial / clinical
        genera (*Pseudomonas*, *Streptomyces*, *Vibrio*, *Salmonella*,
        *Acinetobacter*) as a side-effect of the order medoid.
        """),
        code("""
        picked_after_p5 = picked_after_p4 | {mid for mid, _ in p5}
        p6 = sdt.pass6_hot_taxon_medoids(grower_ids, rxns_by_id, lineage,
                                          sdt.PASS6_HOT_TAXA,
                                          picked_after_p5)
        print(f'Pass 6 picked {len(p6)} hot-taxon medoids')
        for i, (mid, rank, taxon, nm) in enumerate(p6, 1):
            lin = lineage[mid]
            print(f'{i:2d} | {mid:15s} | {rank:7s}={taxon:22s} | '
                  f'{nm:>4d} growers | {lin["phylum"]}/{lin["class"]}')
        """),

        # ------------------------------------------------------------------
        md('## 9. Pass 7 — farthest-point Jaccard fill'),
        md("""
        Fill the remaining slots up to 100 by repeatedly adding the
        grower whose reaction-set Jaccard distance to the panel is
        maximal.  Taxonomy-blind, on purpose — this pass is a pure
        diversity completion layer and lets the algorithm find odd
        outliers that didn't qualify under any other axis.
        """),
        code("""
        picked_after_p6 = picked_after_p5 | {mid for mid, _, _, _ in p6}
        remaining = 100 - len(picked_after_p6)
        print(f'panel size before Pass 7: {len(picked_after_p6)}; '
              f'fill {remaining} more')

        p7 = sdt.pass7_farthest_point(grower_ids, rxns_by_id,
                                      remaining, picked_after_p6)
        for i, (mid, dist) in enumerate(p7[:15], 1):
            lin = lineage[mid]
            print(f'{i:2d} | {mid:15s} | min-Jaccard {dist:.3f} | '
                  f'{lin["phylum"]}/{lin["class"]}/{lin["genus"]}')
        if len(p7) > 15:
            print(f'... {len(p7) - 15} more picks')
        """),

        # ------------------------------------------------------------------
        md('## 10. Run the full selection end-to-end & display the panel'),
        md("""
        Re-run the algorithm through `select_diverse_tax.run_selection`
        so the per-model records carry the merged `reason` strings,
        lineage columns, fingerprint sizes, and rare-reaction counts.
        """),
        code("""
        selection, coverage, tax_cov, rare_set, prevalence, total_cpds = \\
            sdt.run_selection(growers, rxns_by_id, cpds_by_id, lineage)

        pass_dist = Counter(s['pass_origin'] for s in selection)
        print('Final pass distribution:')
        for k in sorted(pass_dist):
            print(f'  pass {k}: {pass_dist[k]} picks')
        print(f'\\nReaction coverage:    {coverage["reactions_covered"]} / '
              f'{coverage["reactions_universe"]}')
        print(f'Metabolite coverage:  {coverage["metabolites_covered"]} / '
              f'{coverage["metabolites_universe"]}')
        """),
        code("""
        # First 25 rows for inspection.
        panel_df = pd.DataFrame(selection)
        panel_df[['model_id', 'pass_origin', 'organism_name', 'phylum',
                  'class', 'order', 'genus', 'n_reactions',
                  'n_rare_rxns', 'growth_flux']].head(25)
        """),

        # ------------------------------------------------------------------
        md('## 11. Write outputs & cache panel for downstream notebooks'),
        md("""
        Write the new panel to `results/selected_models_tax.{csv,json}`,
        `results/selected_ids_tax.txt`, and
        `reports/TAXONOMY_AWARE_SELECTION.md`.  The original
        `selected_models.*` artifacts are LEFT IN PLACE so both panels
        coexist.
        """),
        code("""
        sdt.write_outputs(selection, coverage, tax_cov, suffix='_tax')

        # Cache as a typed object for downstream notebooks
        session.cache.save('tax_aware_panel_v1',
                           {'ids': [s['model_id'] for s in selection],
                            'coverage': coverage,
                            'taxonomy_coverage': tax_cov,
                            'pass_distribution': dict(pass_dist)},
                           type_hint='json',
                           metadata={'n_models': len(selection),
                                     'algorithm': 'tax_aware_v1',
                                     'distinct_phyla': tax_cov['phylum_covered'],
                                     'distinct_classes': tax_cov['class_covered'],
                                     'role': 'taxonomy-aware descriptive test panel'})
        print('cached tax_aware_panel_v1')
        """),

        # ------------------------------------------------------------------
        md('## 12. Side-by-side comparison with the original panel'),
        md("""
        How does the new panel stack up against `results/selected_*.json`
        on the four axes that matter:

        1. **Taxonomic coverage** at each rank (phylum → genus)
        2. **Set coverage** (reactions, metabolites)
        3. **Representativeness** (Jaccard min-distance from every
           non-panel grower to its nearest panel member)
        4. **Per-phylum / per-class panel counts** (visual)
        """),
        code("""
        new = json.loads((RESULTS / 'selected_models_tax.json').read_text())
        # 1. Taxonomic coverage
        new_lin = pd.DataFrame([
            {'model_id': r['model_id'],
             **{rank: r[rank] for rank in sdt.LINEAGE_RANKS}}
            for r in new['records']
        ])

        cmp_rows = []
        for r in sdt.LINEAGE_RANKS:
            orig_set = set(orig_lin[r]) - {'Unknown'}
            new_set = set(new_lin[r]) - {'Unknown'}
            grower_set = set(tax_grower[r]) - {'Unknown'}
            cmp_rows.append({
                'rank': r,
                'grower_universe': len(grower_set),
                'original_panel': len(orig_set),
                'new_panel': len(new_set),
                'delta': len(new_set) - len(orig_set),
            })
        tax_cmp = pd.DataFrame(cmp_rows)
        print('TAXONOMIC COVERAGE:')
        print(tax_cmp.to_string(index=False))
        """),
        code("""
        # 2. Set coverage comparison
        set_rows = [
            ('reactions_covered',
             original['coverage']['reactions_covered'],
             new['coverage']['reactions_covered'],
             original['coverage']['reactions_universe']),
            ('metabolites_covered',
             original['coverage']['metabolites_covered'],
             new['coverage']['metabolites_covered'],
             original['coverage']['metabolites_universe']),
        ]
        set_cmp = pd.DataFrame(set_rows, columns=['metric', 'original', 'new',
                                                  'universe'])
        set_cmp['orig_pct'] = (100 * set_cmp['original'] / set_cmp['universe']).round(1)
        set_cmp['new_pct'] = (100 * set_cmp['new'] / set_cmp['universe']).round(1)
        set_cmp['delta_pct'] = (set_cmp['new_pct'] - set_cmp['orig_pct']).round(1)
        print('SET COVERAGE:')
        print(set_cmp.to_string(index=False))
        """),
        code("""
        # 3. Representativeness: min Jaccard distance from non-panel grower
        # to nearest panel member, for BOTH panels.
        import statistics

        def jacc(a, b):
            if not a and not b: return 0.0
            return 1.0 - len(a & b) / len(a | b)

        def panel_min_dists(panel_ids):
            sel_sets = [rxns_by_id[i] for i in panel_ids]
            non_panel = [m for m in rxns_by_id if m not in set(panel_ids)]
            dists = [min(jacc(rxns_by_id[m], s) for s in sel_sets)
                     for m in non_panel]
            return dists

        dists_orig = panel_min_dists(original['ids'])
        dists_new  = panel_min_dists(new['ids'])

        def stats(ds):
            ds = sorted(ds)
            n = len(ds)
            return {
                'n_non_panel': n,
                'median': statistics.median(ds),
                'mean': sum(ds) / n,
                'p90': ds[int(0.9 * (n - 1))],
                'max': ds[-1],
            }

        jstats = {'original': stats(dists_orig), 'new': stats(dists_new)}
        rep_cmp = pd.DataFrame(jstats).T
        print('REPRESENTATIVENESS (min Jaccard distance to nearest panel member):')
        print(rep_cmp.round(3).to_string())
        """),
        code("""
        # 4a. Per-phylum panel counts (top 20 phyla by grower count).
        top_phyla = (tax_grower['phylum'].value_counts()
                     .drop('Unknown', errors='ignore')
                     .head(20).index.tolist())

        def phylum_counts(panel_lin, phyla):
            c = panel_lin['phylum'].value_counts()
            return [int(c.get(p, 0)) for p in phyla]

        per_phylum = pd.DataFrame({
            'grower_count': [int(tax_grower['phylum'].value_counts().get(p, 0))
                             for p in top_phyla],
            'original_panel': phylum_counts(orig_lin, top_phyla),
            'new_panel': phylum_counts(new_lin, top_phyla),
        }, index=top_phyla)
        per_phylum['delta'] = per_phylum['new_panel'] - per_phylum['original_panel']
        print('Per-phylum panel composition (top 20 grower phyla):')
        print(per_phylum.to_string())
        """),
        code("""
        # 4b. Per-class panel counts (top 15 classes by grower count).
        top_classes = (tax_grower['class'].value_counts()
                       .drop('Unknown', errors='ignore')
                       .head(15).index.tolist())

        def class_counts(panel_lin, classes):
            c = panel_lin['class'].value_counts()
            return [int(c.get(cls, 0)) for cls in classes]

        per_class = pd.DataFrame({
            'grower_count': [int(tax_grower['class'].value_counts().get(c, 0))
                             for c in top_classes],
            'original_panel': class_counts(orig_lin, top_classes),
            'new_panel': class_counts(new_lin, top_classes),
        }, index=top_classes)
        per_class['delta'] = per_class['new_panel'] - per_class['original_panel']
        print('Per-class panel composition (top 15 grower classes):')
        print(per_class.to_string())
        """),
        code("""
        # 4c. Jaccard min-distance histogram (ASCII).
        def bucketize(dists, edges):
            buckets = [0] * (len(edges) - 1)
            for d in dists:
                for i in range(len(edges) - 1):
                    if edges[i] <= d < edges[i + 1]:
                        buckets[i] += 1
                        break
                else:
                    buckets[-1] += 1
            return buckets

        edges = [i / 20 for i in range(11)]  # 0.00, 0.05, 0.10, ..., 0.50
        b_orig = bucketize(dists_orig, edges)
        b_new  = bucketize(dists_new,  edges)

        hist = pd.DataFrame({
            'jaccard_bin': [f'[{edges[i]:.2f}, {edges[i+1]:.2f})'
                            for i in range(len(edges) - 1)],
            'original': b_orig,
            'new':      b_new,
            'delta':    [n - o for o, n in zip(b_orig, b_new)],
        })
        print('Min-Jaccard-distance distribution '
              '(every non-panel grower\\'s distance to its nearest panel member):')
        print(hist.to_string(index=False))
        print(f'\\nOriginal panel worst-case: {max(dists_orig):.3f}')
        print(f'New panel worst-case:      {max(dists_new):.3f}')
        """),

        # ------------------------------------------------------------------
        md('## 13. Taxonomic gap closure'),
        md("""
        Which phyla / classes / orders / genera that were **absent
        from the original panel** are now covered by the new panel?
        And which big-population gaps remain?
        """),
        code("""
        def gap_closure(rank, top=15):
            orig_set = set(orig_lin[rank]) - {'Unknown'}
            new_set  = set(new_lin[rank])  - {'Unknown'}
            grower_counts = tax_grower[rank].value_counts()
            closed = sorted(new_set - orig_set,
                            key=lambda t: -grower_counts.get(t, 0))
            still_missing = sorted(set(grower_counts.index) - new_set - {'Unknown'},
                                   key=lambda t: -grower_counts.get(t, 0))
            return [(t, int(grower_counts.get(t, 0))) for t in closed][:top], \\
                   [(t, int(grower_counts.get(t, 0))) for t in still_missing][:top]

        for r in ['phylum', 'class', 'order', 'genus']:
            closed, missing = gap_closure(r)
            print(f'\\n=== {r.upper()} ===')
            print('NEWLY COVERED by taxonomy-aware panel:')
            for t, n in closed:
                print(f'  + {t:30s} ({n} growers)')
            if not closed:
                print('  (none — same set of taxa as original at this rank)')
            print('STILL MISSING (top remaining gaps):')
            for t, n in missing[:5]:
                print(f'  - {t:30s} ({n} growers)')
        """),
        code("""
        # Persist Jaccard stats + gap closures into the report file.
        gap_for_report = {}
        for r in ['phylum', 'class', 'order', 'genus']:
            closed, _ = gap_closure(r, top=15)
            gap_for_report[r] = closed

        sdt.write_report(selection, coverage, tax_cov,
                         jaccard_stats={'original': stats(dists_orig),
                                        'new':      stats(dists_new)},
                         gap_closures=gap_for_report)
        """),

        # ------------------------------------------------------------------
        md('## 14. Programmatic access'),
        code("""
        from pathlib import Path
        print('From a fresh script:')
        print('  ids = Path(\"/scratch/ctaylor/core_models_analysis/results/selected_ids_tax.txt\").read_text().split()')
        print()
        print('Or, in a notebook that opens its own NotebookSession:')
        print('  panel = session.cache.load(\"tax_aware_panel_v1\")')
        print('  ids = panel[\"ids\"]')
        print('  print(panel[\"taxonomy_coverage\"])')
        """),

        md('---\n## Report: `reports/TAXONOMY_AWARE_SELECTION.md`'),
        render_report(REPORTS_DIR / 'TAXONOMY_AWARE_SELECTION.md'),
    ]
    return new_notebook(
        '# 08 — Taxonomy-Aware Diverse Panel Selection', cells)


def main():
    nb = build_08_taxonomy_aware_selection()
    out = NOTEBOOK_DIR / '08_TaxonomyAwareSelection.ipynb'
    nbf.write(nb, out)
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
