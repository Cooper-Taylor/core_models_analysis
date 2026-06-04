#!/usr/bin/env python3
"""
Select 100 diverse growth models that span the ~3,461 growers.

Approach
--------
A grower is a fingerprint: (set of seed.reaction IDs, set of cpd IDs,
n_reactions, n_metabolites, n_genes, growth_flux, exchanges).

We assemble the 100-model panel in three passes so different axes of
diversity are explicitly represented and easy to justify per pick:

  1. **Reaction-coverage core**  (~20 picks)
     Greedy max-coverage on seed.reaction IDs. Each pick is the model that
     adds the most previously-uncovered reactions. Captures rare enzymes.

  2. **Metabolite-coverage layer** (~15 picks)
     Same greedy idea on cpd IDs. Catches models that share enzymes with
     the core picks but carry distinct metabolite repertoires (different
     transporters / cofactors).

  3. **Farthest-point Jaccard sampling** (~50 picks)
     Repeatedly add the grower whose reaction-set Jaccard distance to the
     current panel is maximal. Spreads picks across the long tail of
     similar-but-not-identical networks.

  4. **Forced extremes** (~15 picks)
     Smallest / largest by reactions, metabolites, genes; lowest / highest
     growth flux; widest / narrowest exchange repertoire; rare-reaction
     champions (highest count of seed.reactions present in <1% of growers).

Each pick is tagged with the dimension(s) that made it stand out.
"""

import csv
import json
import multiprocessing as mp
from collections import Counter, defaultdict
from pathlib import Path

ANALYSIS_DIR = Path('/scratch/ctaylor/core_models_analysis')
MODELS_DIR = ANALYSIS_DIR / 'data' / 'core_models_kegg2'
RESULTS_CSV = ANALYSIS_DIR / 'results' / 'results.csv'
RESULTS_DIR = ANALYSIS_DIR / 'results'
REPORTS_DIR = ANALYSIS_DIR / 'reports'

TARGET_SIZE = 100
RXN_CORE_TARGET = 20
CPD_LAYER_TARGET = 15
EXTREMES_TARGET = 15
# remaining slots filled by farthest-point Jaccard sampling


def load_growers():
    growers = []
    with open(RESULTS_CSV) as f:
        for r in csv.DictReader(f):
            if r['grows'] != 'True':
                continue
            growers.append({
                'model_id': r['model_id'],
                'n_metabolites': int(r['n_metabolites']),
                'n_reactions': int(r['n_reactions']),
                'n_genes': int(r['n_genes']),
                'n_exchanges_total': int(r['n_exchanges_total']),
                'n_exchanges_open': int(r['n_exchanges_open']),
                'growth_flux': float(r['growth_flux']),
            })
    return growers


def extract_one(mid):
    with open(MODELS_DIR / f'{mid}.json') as f:
        m = json.load(f)
    rxns = set()
    for r in m['reactions']:
        sr = r.get('annotation', {}).get('seed.reaction')
        if sr:
            rxns.add(sr)
    cpds = {met['id'].split('_')[0] for met in m['metabolites']}
    return mid, rxns, cpds


def extract_features(growers):
    ids = [g['model_id'] for g in growers]
    rxns_by_id = {}
    cpds_by_id = {}
    with mp.Pool(16) as pool:
        for i, (mid, rxns, cpds) in enumerate(
                pool.imap_unordered(extract_one, ids, chunksize=8), 1):
            rxns_by_id[mid] = rxns
            cpds_by_id[mid] = cpds
            if i % 500 == 0:
                print(f'  extracted {i}/{len(ids)}')
    return rxns_by_id, cpds_by_id


def greedy_max_coverage(ids, sets_by_id, k, already_picked=None):
    """Standard greedy max-coverage. Returns [(id, novelty_count), ...]."""
    picked = list(already_picked or [])
    covered = set()
    for pid in picked:
        covered |= sets_by_id[pid]
    novelty = {}
    for _ in range(k):
        best_id, best_gain = None, -1
        for mid in ids:
            if mid in picked:
                continue
            gain = len(sets_by_id[mid] - covered)
            if gain > best_gain:
                best_id, best_gain = mid, gain
        if best_id is None or best_gain == 0:
            break
        picked.append(best_id)
        novelty[best_id] = best_gain
        covered |= sets_by_id[best_id]
    return picked, covered, novelty


def jaccard(a, b):
    if not a and not b:
        return 0.0
    return 1.0 - len(a & b) / len(a | b)


def farthest_point(ids, rxns_by_id, k, seed_ids):
    """
    Greedy farthest-point: each pick maximises the minimum Jaccard distance
    to any already-selected model. Returns the picks in order.
    """
    picked = list(seed_ids)
    # min_dist[i] = min Jaccard distance from ids[i] to any picked
    id_to_idx = {mid: i for i, mid in enumerate(ids)}
    min_dist = [float('inf')] * len(ids)
    for sid in seed_ids:
        s_set = rxns_by_id[sid]
        for i, mid in enumerate(ids):
            d = jaccard(s_set, rxns_by_id[mid])
            if d < min_dist[i]:
                min_dist[i] = d
    new_picks = []
    for _ in range(k):
        # pick the unpicked id with max min_dist
        best_idx, best_d = -1, -1.0
        for i, mid in enumerate(ids):
            if mid in picked:
                continue
            if min_dist[i] > best_d:
                best_idx, best_d = i, min_dist[i]
        if best_idx < 0:
            break
        new_id = ids[best_idx]
        picked.append(new_id)
        new_picks.append((new_id, best_d))
        # update min_dist with the new pick
        s_set = rxns_by_id[new_id]
        for i, mid in enumerate(ids):
            d = jaccard(s_set, rxns_by_id[mid])
            if d < min_dist[i]:
                min_dist[i] = d
    return new_picks


def pick_extremes(growers_by_id, rxns_by_id, cpds_by_id, rare_set,
                  already_picked):
    """
    Force-include extremes across several axes. Skips any already picked.
    Returns list of (model_id, reason).
    """
    picks = []
    g = list(growers_by_id.values())

    def add(model_id, reason):
        if model_id in already_picked or any(p[0] == model_id for p in picks):
            # update reason on an existing extreme pick
            for i, (pid, r) in enumerate(picks):
                if pid == model_id:
                    picks[i] = (pid, f'{r}; {reason}')
                    return
            return
        picks.append((model_id, reason))

    axes = [
        ('largest by reactions', lambda x: -x['n_reactions']),
        ('smallest by reactions', lambda x: x['n_reactions']),
        ('largest by metabolites', lambda x: -x['n_metabolites']),
        ('smallest by metabolites', lambda x: x['n_metabolites']),
        ('largest by genes', lambda x: -x['n_genes']),
        ('smallest by genes', lambda x: x['n_genes']),
        ('highest growth flux', lambda x: -x['growth_flux']),
        ('lowest growth flux (still growing)', lambda x: x['growth_flux']),
        ('widest open exchange repertoire', lambda x: -x['n_exchanges_open']),
        ('narrowest open exchange repertoire', lambda x: x['n_exchanges_open']),
        ('largest reaction fingerprint',
            lambda x: -len(rxns_by_id[x['model_id']])),
        ('largest metabolite fingerprint',
            lambda x: -len(cpds_by_id[x['model_id']])),
    ]
    for reason, key in axes:
        winner = min(g, key=key)
        add(winner['model_id'], reason)

    # Rare-reaction champions: most rare seed.reactions present
    rare_count = {gid: len(rxns_by_id[gid] & rare_set) for gid in growers_by_id}
    for mid, _ in sorted(rare_count.items(), key=lambda kv: -kv[1])[:3]:
        add(mid, 'carries many rare reactions (<1% of growers)')

    return picks


def main():
    print('Loading grower results...')
    growers = load_growers()
    print(f'  {len(growers)} growers')

    print('Extracting reaction & metabolite fingerprints...')
    rxns_by_id, cpds_by_id = extract_features(growers)
    growers_by_id = {g['model_id']: g for g in growers}
    ids = list(growers_by_id)

    # Global reaction prevalence (for "rare" definition)
    prevalence = Counter()
    for s in rxns_by_id.values():
        prevalence.update(s)
    total = len(ids)
    rare_set = {r for r, c in prevalence.items() if c / total < 0.01}
    print(f'  universe: {len(prevalence)} unique seed.reactions; '
          f'{len(rare_set)} appear in <1% of growers')

    # Per-pick metadata
    selection = []  # list of dicts: model_id, reason, novelty, ...
    picked_ids = set()

    def record(model_id, reason, **extra):
        if model_id in picked_ids:
            for s in selection:
                if s['model_id'] == model_id:
                    s['reason'] = f'{s["reason"]}; {reason}'
                    s.update(extra)
                    return
            return
        picked_ids.add(model_id)
        row = {'model_id': model_id, 'reason': reason, **extra}
        selection.append(row)

    # 1. Reaction-coverage core
    print(f'\nPass 1: greedy reaction max-coverage ({RXN_CORE_TARGET} picks)')
    picks1, covered_rxns, novelty1 = greedy_max_coverage(
        ids, rxns_by_id, RXN_CORE_TARGET)
    for i, pid in enumerate(picks1):
        record(pid, f'reaction-coverage rank {i+1}: '
                    f'adds {novelty1[pid]} previously-uncovered reactions',
               novel_rxns_added=novelty1[pid])
    print(f'  covered {len(covered_rxns)}/{len(prevalence)} reactions')

    # 2. Metabolite-coverage layer
    print(f'\nPass 2: greedy metabolite max-coverage '
          f'({CPD_LAYER_TARGET} picks, seeded with reaction core)')
    picks2, covered_cpds, novelty2 = greedy_max_coverage(
        ids, cpds_by_id, CPD_LAYER_TARGET, already_picked=picks1)
    new_in_2 = [p for p in picks2 if p not in picks1]
    for i, pid in enumerate(new_in_2):
        record(pid, f'metabolite-coverage rank {i+1}: '
                    f'adds {novelty2[pid]} previously-uncovered metabolites',
               novel_cpds_added=novelty2[pid])
    total_cpds = len({c for s in cpds_by_id.values() for c in s})
    print(f'  covered {len(covered_cpds)}/{total_cpds} metabolites')

    # 3. Forced extremes
    print(f'\nPass 3: extreme-axis picks')
    extremes = pick_extremes(growers_by_id, rxns_by_id, cpds_by_id,
                              rare_set, picked_ids)
    for pid, reason in extremes:
        record(pid, f'extreme: {reason}')
    print(f'  {len(extremes)} new extreme picks (selection now {len(selection)})')

    # 4. Farthest-point Jaccard sampling for the rest
    remaining = TARGET_SIZE - len(selection)
    print(f'\nPass 4: farthest-point Jaccard sampling for {remaining} more picks')
    new_picks = farthest_point(ids, rxns_by_id, remaining,
                                seed_ids=list(picked_ids))
    for pid, dist in new_picks:
        record(pid, f'farthest-point: min Jaccard distance to '
                    f'already-selected = {dist:.3f}',
               min_jaccard_to_selected=round(dist, 4))

    # Trim if we somehow overshot
    selection = selection[:TARGET_SIZE]
    print(f'\nFinal panel size: {len(selection)}')

    # Coverage check
    final_rxn_coverage = set()
    final_cpd_coverage = set()
    for s in selection:
        final_rxn_coverage |= rxns_by_id[s['model_id']]
        final_cpd_coverage |= cpds_by_id[s['model_id']]
    print(f'Panel covers {len(final_rxn_coverage)}/{len(prevalence)} '
          f'unique reactions and {len(final_cpd_coverage)}/{total_cpds} '
          f'unique metabolites across all growers.')

    # ------------------------------------------------------------------
    # Augment selection rows with the model stats and write outputs
    for s in selection:
        g = growers_by_id[s['model_id']]
        s['n_reactions'] = g['n_reactions']
        s['n_metabolites'] = g['n_metabolites']
        s['n_genes'] = g['n_genes']
        s['n_exchanges_open'] = g['n_exchanges_open']
        s['growth_flux'] = round(g['growth_flux'], 3)
        s['n_unique_rxns'] = len(rxns_by_id[s['model_id']])
        s['n_unique_cpds'] = len(cpds_by_id[s['model_id']])
        s['n_rare_rxns'] = len(rxns_by_id[s['model_id']] & rare_set)

    out_csv = RESULTS_DIR / 'selected_models.csv'
    fields = ['model_id', 'reason', 'n_reactions', 'n_metabolites', 'n_genes',
              'n_exchanges_open', 'growth_flux', 'n_unique_rxns',
              'n_unique_cpds', 'n_rare_rxns', 'novel_rxns_added',
              'novel_cpds_added', 'min_jaccard_to_selected']
    with open(out_csv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for s in selection:
            w.writerow({k: s.get(k, '') for k in fields})
    print(f'\nWrote {out_csv}')

    out_ids = RESULTS_DIR / 'selected_ids.txt'
    with open(out_ids, 'w') as f:
        for s in selection:
            f.write(s['model_id'] + '\n')
    print(f'Wrote {out_ids}')

    out_json = RESULTS_DIR / 'selected_models.json'
    with open(out_json, 'w') as f:
        json.dump({'ids': [s['model_id'] for s in selection],
                   'records': selection,
                   'coverage': {
                       'reactions_covered': len(final_rxn_coverage),
                       'reactions_universe': len(prevalence),
                       'metabolites_covered': len(final_cpd_coverage),
                       'metabolites_universe': total_cpds,
                   }}, f, indent=2)
    print(f'Wrote {out_json}')

    # ------------------------------------------------------------------
    # Markdown writeup
    md = REPORTS_DIR / 'DIVERSE_SELECTION.md'
    lines = []
    lines.append('# 100 Diverse Growth Models — Selection Report\n')
    lines.append(f'Selected from **{len(growers)} growing models** '
                 f'(of 5,683 total in `core_models_kegg2`).\n')
    lines.append('## Methodology\n')
    lines.append('Four passes, each tagged onto every pick so the role '
                 'of the model in the panel is explicit:\n')
    lines.append('1. **Reaction-coverage core** — greedy max-coverage on '
                 '`seed.reaction` IDs across growers. Each pick contributes '
                 'the most previously-uncovered reactions.')
    lines.append('2. **Metabolite-coverage layer** — same greedy idea on '
                 'compound IDs, seeded with the reaction core. Catches '
                 'transporter / cofactor diversity that the reaction-set '
                 'pass missed.')
    lines.append('3. **Forced extremes** — smallest/largest by reactions, '
                 'metabolites, genes; lowest/highest growth flux; '
                 'widest/narrowest open-exchange repertoire; '
                 'rare-reaction champions.')
    lines.append('4. **Farthest-point Jaccard sampling** — fills the rest '
                 'of the 100 slots by repeatedly adding the grower whose '
                 'reaction-set Jaccard distance to the already-picked '
                 'panel is maximal.\n')

    lines.append('## Coverage achieved by the 100-model panel\n')
    lines.append(f'- Reactions covered: **{len(final_rxn_coverage)} / '
                 f'{len(prevalence)}** unique `seed.reaction` IDs '
                 f'({100*len(final_rxn_coverage)/len(prevalence):.1f}%)')
    lines.append(f'- Metabolites covered: **{len(final_cpd_coverage)} / '
                 f'{total_cpds}** unique cpd IDs '
                 f'({100*len(final_cpd_coverage)/total_cpds:.1f}%)')
    lines.append(f'- Rare reactions (<1% prevalence) hit by panel: see '
                 f'`n_rare_rxns` column in CSV.\n')

    lines.append('## Files\n')
    lines.append('- `selected_ids.txt` — one model ID per line, ready for `for id in $(cat selected_ids.txt); do ...; done`')
    lines.append('- `selected_models.csv` — per-model metrics + selection reason')
    lines.append('- `selected_models.json` — same data, plus coverage stats, for programmatic use\n')

    lines.append('## Panel members (in the order they joined)\n')
    lines.append('| # | model_id | reaction count | growth flux | selection reason |')
    lines.append('|---|---|---|---|---|')
    for i, s in enumerate(selection, 1):
        lines.append(f'| {i} | `{s["model_id"]}` | {s["n_reactions"]} | '
                     f'{s["growth_flux"]} | {s["reason"]} |')
    md.write_text('\n'.join(lines) + '\n')
    print(f'Wrote {md}')


if __name__ == '__main__':
    main()
