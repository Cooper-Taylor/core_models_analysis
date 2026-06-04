#!/usr/bin/env python3
"""
Taxonomy-aware diverse panel selection: 100 representative growth models
drawn from the 3,461 growers, balancing network-set coverage with NCBI
phylogenetic spread.

Motivation
----------
The original ``select_diverse.py`` panel (see ``DIVERSE_SELECTION.md``)
is excellent at covering reactions and metabolites but only hits **19 of
~33 grower-bearing phyla**.  It misses Myxococcota (21 growers),
Sphingobacteriia (20), Micrococcales (122 growers in the panel-absent
order alone), and big industrial/clinical genera like *Pseudomonas* (84
growers) and *Burkholderia* (81).  This script reweights the four
original passes around an explicit taxonomic stratification so every
grower-bearing phylum gets at least one anchor and the panel spreads
deliberately across class / order / genus levels too.

Seven passes (target = 100)
---------------------------
1. **Phylum medoids** — one model per grower-bearing phylum
   (the model with min sum-Jaccard distance to its phylum-mates on
   reaction fingerprints).  Anchors the panel in every phylum.
2. **Reaction-coverage core** — greedy max-coverage on
   ``seed.reaction`` IDs, seeded with Pass 1.  Captures rare enzymes.
   12 picks.
3. **Taxonomic-novelty fill** — greedy farthest-point on
   *lineage rank-distance* (7=different superkingdom, 6=different
   phylum, …, 1=different species).  Tie-break with Jaccard reaction
   distance.  Closes class / order gaps. 12 picks.
4. **Metabolite-coverage layer** — greedy max-coverage on cpd IDs,
   seeded with all prior picks. 8 picks.
5. **Constrained extremes** — the original extreme axes
   (min/max reactions, metabolites, genes, flux, exchanges,
   fingerprints) PLUS top-3 rare-reaction carriers, with a class
   over-representation guard (skip a candidate whose class is already
   >2× its expected share; walk down the ranking until one passes).
   8 picks.
6. **Hot-taxon medoids** — for the top missing grower-heavy genera /
   orders / families (those with the most growers but no panel
   representative yet), add the medoid of each.  Closes within-phylum
   genus gaps that Pass 3 misses because it prioritizes higher-rank
   distance.  10 picks.
7. **Farthest-point Jaccard fill** — fills the remaining slots
   by maximizing min Jaccard reaction-set distance to the panel.
   Pure diversity completion.

Each pick is tagged with its pass and a reason string.  Records are
merged when a model qualifies under multiple passes; ``pass_origin``
records the FIRST pass that picked it.

Outputs
-------
- ``results/selected_ids_tax.txt``      — 100 model IDs, one per line
- ``results/selected_models_tax.csv``   — per-model metrics + reason + lineage
- ``results/selected_models_tax.json``  — same data + coverage + taxonomy stats
- ``reports/TAXONOMY_AWARE_SELECTION.md`` — narrative + panel table

The original ``selected_*`` artifacts are LEFT IN PLACE so both panels
can coexist for comparison.
"""

from __future__ import annotations
import os

import csv
import json
import multiprocessing as mp
from collections import Counter, defaultdict
from pathlib import Path

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR", "/scratch/ctaylor/core_models_analysis"))
MODELS_DIR = ANALYSIS_DIR / 'data' / 'core_models_kegg2'
RESULTS_CSV = ANALYSIS_DIR / 'results' / 'results.csv'
TAXONOMY_CSV = ANALYSIS_DIR / 'results' / 'ncbi_taxonomy.csv'
RESULTS_DIR = ANALYSIS_DIR / 'results'
REPORTS_DIR = ANALYSIS_DIR / 'reports'

TARGET_SIZE = 100
PASS2_RXN_CORE = 12
PASS3_TAX_NOVELTY = 12
PASS4_CPD_LAYER = 8
PASS5_EXTREMES = 8
PASS6_HOT_TAXA = 14
# Pass 1 takes whatever the phylum count allows, Pass 7 fills the rest.

# Ranks (broadest-first) that Pass 6 sweeps to find big missing taxa.
# Orders first so broad missing clades like Rhodobacterales (95 growers,
# no individual genus above threshold) get covered; then genera so big
# industrial / clinical bugs (Burkholderia, Escherichia, Streptomyces,
# Vibrio) that sit inside already-anchored orders are pulled in too.
HOT_TAXA_RANKS = ['order', 'genus', 'family']
# Skip taxa with fewer than this many growers — anything thinner gets
# picked up by Pass 7 (farthest-point) anyway.
HOT_TAXA_MIN_GROWERS = 20

LINEAGE_RANKS = ['superkingdom', 'phylum', 'class',
                 'order', 'family', 'genus', 'species']

# Pass 1 caps the medoid computation per phylum at this many candidates
# (sampled deterministically by sorted model_id).  Pseudomonadota has
# ~1,909 growers; full pairwise on that is 1.8M Jaccards.  50 keeps it
# fast and is plenty for picking a central representative.
PHYLUM_MEDOID_SAMPLE = 50

# Pass 5 walks down a sorted axis up to this many candidates before
# accepting the unconstrained extreme.
EXTREMES_WALK_LIMIT = 20

# Pass 5 over-representation threshold: skip a candidate whose class is
# already represented > THIS_MULT times its expected share.
EXTREMES_CLASS_MULT = 2.0


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------
def load_growers():
    """Return a list of grower dicts pulled from results.csv."""
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


def load_taxonomy():
    """Return dict[model_id] -> lineage dict.

    Missing rows are bucketed as 'Unknown' at every rank.  We do NOT
    propagate higher ranks for unresolved organisms; the goal is for
    Pass 1 to anchor the Unknown bucket as a single 'phylum'.
    """
    tax = {}
    with open(TAXONOMY_CSV) as f:
        for r in csv.DictReader(f):
            mid = r['assembly_accession']
            if r.get('status') == 'resolved':
                lin = {k: (r.get(k) or 'Unknown') for k in LINEAGE_RANKS}
            else:
                lin = {k: 'Unknown' for k in LINEAGE_RANKS}
            lin['organism_name'] = r.get('organism_name') or ''
            lin['tax_id'] = r.get('tax_id') or ''
            tax[mid] = lin
    return tax


def _extract_one(mid):
    with open(MODELS_DIR / f'{mid}.json') as f:
        m = json.load(f)
    rxns = set()
    for r in m['reactions']:
        sr = r.get('annotation', {}).get('seed.reaction')
        if sr:
            rxns.add(sr)
    cpds = {met['id'].split('_')[0] for met in m['metabolites']}
    return mid, rxns, cpds


def extract_features(growers, *, processes=16):
    """Extract reaction and metabolite fingerprints for every grower."""
    ids = [g['model_id'] for g in growers]
    rxns_by_id, cpds_by_id = {}, {}
    with mp.Pool(processes) as pool:
        for i, (mid, rxns, cpds) in enumerate(
                pool.imap_unordered(_extract_one, ids, chunksize=8), 1):
            rxns_by_id[mid] = rxns
            cpds_by_id[mid] = cpds
            if i % 500 == 0:
                print(f'  extracted {i}/{len(ids)}')
    return rxns_by_id, cpds_by_id


# ---------------------------------------------------------------------------
# Distance primitives
# ---------------------------------------------------------------------------
def jaccard(a, b):
    if not a and not b:
        return 0.0
    return 1.0 - len(a & b) / len(a | b)


def lineage_distance(lin_a, lin_b):
    """Return integer 0-7 capturing how distant two NCBI lineages are.

    7 = different superkingdom (or one is Unknown at superkingdom)
    6 = same superkingdom, different phylum
    5 = same phylum, different class
    4 = same class, different order
    3 = same order, different family
    2 = same family, different genus
    1 = same genus, different species
    0 = identical lineage down through species

    Both lineages treat 'Unknown' as a sentinel value:
    Unknown-vs-Unknown matches (so the Unknown bucket clusters); but
    Unknown-vs-resolved differs at whatever rank Unknown shows up.
    """
    for i, rank in enumerate(LINEAGE_RANKS):
        if lin_a.get(rank) != lin_b.get(rank):
            return 7 - i
    return 0


# ---------------------------------------------------------------------------
# Pass 1 — phylum medoids
# ---------------------------------------------------------------------------
def pass1_phylum_medoids(ids, rxns_by_id, lineage):
    """Pick one medoid per grower-bearing phylum.

    For each phylum P:
        candidates = sorted ids in P (cap sample at PHYLUM_MEDOID_SAMPLE)
        medoid = argmin over candidates of sum_{other in candidates}
                 jaccard(rxns[c], rxns[other])
        for singletons / pairs, just take the first.

    Returns list[(model_id, phylum, sum_jaccard)] in deterministic order
    (largest phylum first by grower count).
    """
    by_phylum = defaultdict(list)
    for mid in ids:
        by_phylum[lineage[mid]['phylum']].append(mid)
    # Largest phyla first so the panel order is meaningful.
    sorted_phyla = sorted(by_phylum.items(),
                          key=lambda kv: (-len(kv[1]), kv[0]))

    picks = []
    for phylum, members in sorted_phyla:
        members = sorted(members)
        if len(members) <= 2:
            mid = members[0]
            picks.append((mid, phylum, 0.0))
            continue
        sample = members
        if len(sample) > PHYLUM_MEDOID_SAMPLE:
            stride = len(sample) / PHYLUM_MEDOID_SAMPLE
            sample = [members[int(i * stride)]
                      for i in range(PHYLUM_MEDOID_SAMPLE)]
        # Sum Jaccard from each candidate to all others in the sample.
        sums = {}
        for c in sample:
            s = 0.0
            cs = rxns_by_id[c]
            for o in sample:
                if o == c:
                    continue
                s += jaccard(cs, rxns_by_id[o])
            sums[c] = s
        medoid = min(sums, key=sums.get)
        picks.append((medoid, phylum, round(sums[medoid], 3)))
    return picks


# ---------------------------------------------------------------------------
# Passes 2 & 4 — greedy set max-coverage
# ---------------------------------------------------------------------------
def greedy_max_coverage(ids, sets_by_id, k, already_picked):
    """Standard greedy. Returns [(id, novelty_count), ...]."""
    picked = list(already_picked)
    covered = set()
    for pid in picked:
        covered |= sets_by_id.get(pid, set())
    novelty = {}
    new_picks = []
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
        new_picks.append((best_id, best_gain))
    return new_picks, covered


# ---------------------------------------------------------------------------
# Pass 3 — taxonomic-novelty fill
# ---------------------------------------------------------------------------
def pass3_taxonomic_novelty(ids, rxns_by_id, lineage, k, already_picked):
    """Greedy farthest-point on lineage rank-distance.

    Per candidate, score = min lineage_distance to any picked.
    Tie-break: max Jaccard reaction distance to nearest picked.

    Returns list[(model_id, min_rank_distance, tiebreak_jaccard)].
    """
    picked = list(already_picked)
    # Precompute initial min_rank and min_jacc for each id.
    min_rank = {}
    min_jacc = {}
    for mid in ids:
        if mid in picked:
            continue
        mr, mj = 8, 1.0
        for sid in picked:
            d = lineage_distance(lineage[mid], lineage[sid])
            if d < mr:
                mr = d
                mj = jaccard(rxns_by_id[mid], rxns_by_id[sid])
            elif d == mr:
                jd = jaccard(rxns_by_id[mid], rxns_by_id[sid])
                if jd < mj:
                    mj = jd
        min_rank[mid] = mr
        min_jacc[mid] = mj
    new_picks = []
    for _ in range(k):
        best = None
        best_key = (-1, -1.0)
        for mid, mr in min_rank.items():
            if mid in picked:
                continue
            mj = min_jacc[mid]
            key = (mr, mj)
            if key > best_key:
                best, best_key = mid, key
        if best is None or best_key[0] <= 0:
            break
        picked.append(best)
        new_picks.append((best, best_key[0], round(best_key[1], 3)))
        # Update min_rank/min_jacc with the new pick.
        new_lin = lineage[best]
        new_rxn = rxns_by_id[best]
        for mid in list(min_rank):
            if mid in picked:
                continue
            d = lineage_distance(lineage[mid], new_lin)
            if d < min_rank[mid]:
                min_rank[mid] = d
                min_jacc[mid] = jaccard(rxns_by_id[mid], new_rxn)
            elif d == min_rank[mid]:
                jd = jaccard(rxns_by_id[mid], new_rxn)
                if jd < min_jacc[mid]:
                    min_jacc[mid] = jd
    return new_picks


# ---------------------------------------------------------------------------
# Pass 5 — constrained extremes
# ---------------------------------------------------------------------------
def pass5_constrained_extremes(growers_by_id, rxns_by_id, cpds_by_id,
                                rare_set, lineage, already_picked,
                                cap=PASS5_EXTREMES):
    """Force-include axis-extreme models with a class-saturation guard.

    Each of the 12 axes is sorted; we walk the top EXTREMES_WALK_LIMIT
    candidates and pick the first whose class is at most
    EXTREMES_CLASS_MULT × its expected share among the panel so far.
    """
    g = list(growers_by_id.values())
    n_growers = len(g)
    grower_class_n = Counter(lineage[gg['model_id']]['class'] for gg in g)

    def panel_class_share(picked):
        cnt = Counter(lineage[p]['class'] for p in picked)
        total = max(len(picked), 1)
        return {c: n / total for c, n in cnt.items()}, total

    picks = []
    picked_ids = set(already_picked)

    def add(model_id, reason):
        if model_id in picked_ids:
            for i, (pid, r) in enumerate(picks):
                if pid == model_id:
                    picks[i] = (pid, f'{r}; {reason}')
                    return
            return
        picks.append((model_id, reason))
        picked_ids.add(model_id)

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
        if len(picks) >= cap:
            break
        ranked = sorted(g, key=key)
        share, total = panel_class_share(picked_ids)
        chosen = None
        for cand in ranked[:EXTREMES_WALK_LIMIT]:
            cls = lineage[cand['model_id']]['class']
            expected = grower_class_n.get(cls, 1) / n_growers
            current = share.get(cls, 0.0)
            # No saturation issue if the class has no panel members yet.
            if total == 0 or current <= EXTREMES_CLASS_MULT * expected:
                chosen = cand
                break
        if chosen is None:
            chosen = ranked[0]
        add(chosen['model_id'], f'extreme: {reason}')

    # Rare-reaction champions (top 3)
    if len(picks) < cap:
        rare_count = {gid: len(rxns_by_id[gid] & rare_set) for gid in growers_by_id}
        for mid, _ in sorted(rare_count.items(), key=lambda kv: -kv[1]):
            if len(picks) >= cap:
                break
            add(mid, 'extreme: carries many rare reactions (<1% of growers)')
    return picks


# ---------------------------------------------------------------------------
# Pass 6 — hot-taxon medoids
# ---------------------------------------------------------------------------
def pass6_hot_taxon_medoids(ids, rxns_by_id, lineage, k, already_picked,
                             ranks=HOT_TAXA_RANKS,
                             min_growers=HOT_TAXA_MIN_GROWERS):
    """Add medoids for the top missing high-prevalence taxa.

    For each rank in ``ranks`` (deepest first), enumerate taxa sorted
    by grower count descending.  Skip taxa already represented in the
    panel.  For each missing taxon, pick its medoid (min sum-Jaccard
    on reaction fingerprints, sample-capped at PHYLUM_MEDOID_SAMPLE).
    Stops at ``k`` total picks or when every rank is exhausted.
    """
    by_rank_taxon = {r: defaultdict(list) for r in ranks}
    for mid in ids:
        for r in ranks:
            t = lineage[mid].get(r, 'Unknown')
            if t == 'Unknown':
                continue
            by_rank_taxon[r][t].append(mid)

    panel_present = {r: {lineage[p].get(r) for p in already_picked}
                     for r in ranks}

    new_picks = []
    picked_now = set(already_picked)
    for r in ranks:
        # Sorted by count desc, then alphabetical (deterministic).
        candidates = sorted(by_rank_taxon[r].items(),
                            key=lambda kv: (-len(kv[1]), kv[0]))
        for taxon, members in candidates:
            if len(new_picks) >= k:
                break
            if taxon in panel_present[r]:
                continue
            if len(members) < min_growers:
                # Remaining taxa at this rank are all too small — move on.
                break
            members = sorted(members)
            if len(members) <= 2:
                medoid = members[0]
            else:
                sample = members
                if len(sample) > PHYLUM_MEDOID_SAMPLE:
                    stride = len(sample) / PHYLUM_MEDOID_SAMPLE
                    sample = [members[int(i * stride)]
                              for i in range(PHYLUM_MEDOID_SAMPLE)]
                sums = {}
                for c in sample:
                    cs = rxns_by_id[c]
                    s = 0.0
                    for o in sample:
                        if o == c:
                            continue
                        s += jaccard(cs, rxns_by_id[o])
                    sums[c] = s
                medoid = min(sums, key=sums.get)
            if medoid in picked_now:
                # Medoid happens to already be in panel under another tag.
                panel_present[r].add(taxon)
                continue
            new_picks.append((medoid, r, taxon, len(members)))
            picked_now.add(medoid)
            # Update panel_present across all ranks since this pick lands
            # in concrete lineages we don't want to double-count.
            for rr in ranks:
                panel_present[rr].add(lineage[medoid].get(rr))
        if len(new_picks) >= k:
            break
    return new_picks


# ---------------------------------------------------------------------------
# Pass 7 — farthest-point Jaccard fill
# ---------------------------------------------------------------------------
def pass7_farthest_point(ids, rxns_by_id, k, already_picked):
    picked = list(already_picked)
    min_dist = {mid: float('inf') for mid in ids if mid not in picked}
    for sid in picked:
        s_set = rxns_by_id[sid]
        for mid in min_dist:
            d = jaccard(s_set, rxns_by_id[mid])
            if d < min_dist[mid]:
                min_dist[mid] = d
    new_picks = []
    for _ in range(k):
        if not min_dist:
            break
        best_id = max(min_dist, key=min_dist.get)
        d = min_dist.pop(best_id)
        new_picks.append((best_id, round(d, 4)))
        s_set = rxns_by_id[best_id]
        for mid in min_dist:
            nd = jaccard(s_set, rxns_by_id[mid])
            if nd < min_dist[mid]:
                min_dist[mid] = nd
    return new_picks


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def run_selection(growers, rxns_by_id, cpds_by_id, lineage,
                  target_size=TARGET_SIZE):
    """Run the full 6-pass algorithm and return a list of selection
    dicts plus aggregate stats.
    """
    growers_by_id = {g['model_id']: g for g in growers}
    ids = list(growers_by_id)

    # Global reaction prevalence (for rare-set definition).
    prevalence = Counter()
    for s in rxns_by_id.values():
        prevalence.update(s)
    total_models = len(ids)
    rare_set = {r for r, c in prevalence.items() if c / total_models < 0.01}
    total_cpds = len({c for s in cpds_by_id.values() for c in s})

    selection = []
    picked_ids = set()
    pass_origin = {}

    def record(model_id, reason, pass_num, **extra):
        if model_id in picked_ids:
            for s in selection:
                if s['model_id'] == model_id:
                    s['reason'] = f'{s["reason"]}; {reason}'
                    s.update({k: v for k, v in extra.items()
                              if s.get(k) in (None, '')})
                    return
            return
        picked_ids.add(model_id)
        pass_origin[model_id] = pass_num
        selection.append({'model_id': model_id, 'reason': reason,
                          'pass_origin': pass_num, **extra})

    # Pass 1 — phylum medoids
    print('Pass 1: phylum medoids')
    p1 = pass1_phylum_medoids(ids, rxns_by_id, lineage)
    for mid, phylum, sj in p1:
        record(mid, f'phylum-medoid: {phylum} (sum_jaccard={sj})', 1,
               medoid_phylum=phylum, medoid_sum_jaccard=sj)
    print(f'  picked {len(p1)} medoids (one per grower-bearing phylum); '
          f'panel size now {len(selection)}')

    # Pass 2 — reaction coverage core
    print(f'Pass 2: greedy reaction max-coverage ({PASS2_RXN_CORE} picks)')
    p2, _ = greedy_max_coverage(ids, rxns_by_id, PASS2_RXN_CORE, picked_ids)
    for i, (mid, gain) in enumerate(p2, 1):
        record(mid,
               f'reaction-coverage rank {i}: adds {gain} previously-uncovered reactions',
               2, novel_rxns_added=gain)
    print(f'  added {len(p2)} picks (panel size {len(selection)})')

    # Pass 3 — taxonomic novelty fill
    print(f'Pass 3: taxonomic-novelty greedy ({PASS3_TAX_NOVELTY} picks)')
    p3 = pass3_taxonomic_novelty(ids, rxns_by_id, lineage,
                                 PASS3_TAX_NOVELTY, picked_ids)
    for mid, rank_d, tjac in p3:
        lin = lineage[mid]
        record(mid,
               f'taxonomic-novelty: rank-distance {rank_d} '
               f'({lin["phylum"]}/{lin["class"]}/{lin["order"]}/{lin["genus"]})',
               3, tax_rank_distance=rank_d, tax_tiebreak_jaccard=tjac)
    print(f'  added {len(p3)} picks (panel size {len(selection)})')

    # Pass 4 — metabolite coverage layer
    print(f'Pass 4: greedy metabolite max-coverage ({PASS4_CPD_LAYER} picks)')
    p4, _ = greedy_max_coverage(ids, cpds_by_id, PASS4_CPD_LAYER, picked_ids)
    for i, (mid, gain) in enumerate(p4, 1):
        record(mid,
               f'metabolite-coverage rank {i}: adds {gain} previously-uncovered metabolites',
               4, novel_cpds_added=gain)
    print(f'  added {len(p4)} picks (panel size {len(selection)})')

    # Pass 5 — constrained extremes
    print(f'Pass 5: constrained extremes (cap {PASS5_EXTREMES})')
    p5 = pass5_constrained_extremes(growers_by_id, rxns_by_id, cpds_by_id,
                                     rare_set, lineage, picked_ids,
                                     cap=PASS5_EXTREMES)
    for mid, reason in p5:
        record(mid, reason, 5)
    new_p5 = [(mid, r) for mid, r in p5 if pass_origin.get(mid) == 5]
    print(f'  added {len(new_p5)} new extreme picks (panel size {len(selection)})')

    # Pass 6 — hot-taxon medoids
    print(f'Pass 6: hot-taxon medoids (cap {PASS6_HOT_TAXA})')
    p6 = pass6_hot_taxon_medoids(ids, rxns_by_id, lineage,
                                  PASS6_HOT_TAXA, picked_ids)
    for mid, rank, taxon, n_members in p6:
        record(mid,
               f'hot-taxon medoid: {rank}={taxon} ({n_members} growers)',
               6, hot_taxon_rank=rank, hot_taxon=taxon,
               hot_taxon_grower_count=n_members)
    print(f'  added {len(p6)} picks (panel size {len(selection)})')

    # Pass 7 — farthest-point Jaccard fill
    remaining = target_size - len(selection)
    print(f'Pass 7: farthest-point Jaccard fill ({remaining} more picks)')
    p7 = pass7_farthest_point(ids, rxns_by_id, remaining, picked_ids)
    for mid, dist in p7:
        record(mid,
               f'farthest-point: min Jaccard distance to selected = {dist:.3f}',
               7, min_jaccard_to_selected=dist)
    print(f'  added {len(p7)} picks (final panel size {len(selection)})')

    selection = selection[:target_size]

    # Augment with model stats + fingerprint sizes + lineage.
    for s in selection:
        g = growers_by_id[s['model_id']]
        lin = lineage[s['model_id']]
        s['n_reactions'] = g['n_reactions']
        s['n_metabolites'] = g['n_metabolites']
        s['n_genes'] = g['n_genes']
        s['n_exchanges_open'] = g['n_exchanges_open']
        s['growth_flux'] = round(g['growth_flux'], 3)
        s['n_unique_rxns'] = len(rxns_by_id[s['model_id']])
        s['n_unique_cpds'] = len(cpds_by_id[s['model_id']])
        s['n_rare_rxns'] = len(rxns_by_id[s['model_id']] & rare_set)
        s['organism_name'] = lin.get('organism_name', '')
        s['tax_id'] = lin.get('tax_id', '')
        for r in LINEAGE_RANKS:
            s[r] = lin.get(r, 'Unknown')

    # Coverage summary
    final_rxn = set()
    final_cpd = set()
    for s in selection:
        final_rxn |= rxns_by_id[s['model_id']]
        final_cpd |= cpds_by_id[s['model_id']]
    coverage = {
        'reactions_covered': len(final_rxn),
        'reactions_universe': len(prevalence),
        'metabolites_covered': len(final_cpd),
        'metabolites_universe': total_cpds,
    }

    # Taxonomy coverage summary
    grower_distinct = {r: len({lineage[m][r] for m in ids})
                       for r in LINEAGE_RANKS}
    panel_distinct = {r: len({s[r] for s in selection})
                      for r in LINEAGE_RANKS}
    tax_coverage = {f'{r}_covered': panel_distinct[r] for r in LINEAGE_RANKS}
    tax_coverage.update({f'{r}_universe': grower_distinct[r]
                         for r in LINEAGE_RANKS})

    return selection, coverage, tax_coverage, rare_set, prevalence, total_cpds


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------
PANEL_CSV_FIELDS = [
    'model_id', 'pass_origin', 'reason',
    'organism_name', 'tax_id',
    'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species',
    'n_reactions', 'n_metabolites', 'n_genes', 'n_exchanges_open',
    'growth_flux', 'n_unique_rxns', 'n_unique_cpds', 'n_rare_rxns',
    'medoid_phylum', 'medoid_sum_jaccard',
    'novel_rxns_added', 'novel_cpds_added',
    'tax_rank_distance', 'tax_tiebreak_jaccard',
    'hot_taxon_rank', 'hot_taxon', 'hot_taxon_grower_count',
    'min_jaccard_to_selected',
]


def write_outputs(selection, coverage, tax_coverage, *, suffix='_tax'):
    out_csv = RESULTS_DIR / f'selected_models{suffix}.csv'
    with open(out_csv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=PANEL_CSV_FIELDS)
        w.writeheader()
        for s in selection:
            w.writerow({k: s.get(k, '') for k in PANEL_CSV_FIELDS})
    print(f'wrote {out_csv}')

    out_ids = RESULTS_DIR / f'selected_ids{suffix}.txt'
    with open(out_ids, 'w') as f:
        for s in selection:
            f.write(s['model_id'] + '\n')
    print(f'wrote {out_ids}')

    out_json = RESULTS_DIR / f'selected_models{suffix}.json'
    pass_dist = Counter(s['pass_origin'] for s in selection)
    with open(out_json, 'w') as f:
        json.dump({
            'ids': [s['model_id'] for s in selection],
            'records': selection,
            'coverage': coverage,
            'taxonomy_coverage': tax_coverage,
            'pass_distribution': dict(pass_dist),
            'algorithm': 'tax_aware_v1',
        }, f, indent=2)
    print(f'wrote {out_json}')

    return out_csv, out_ids, out_json


def write_report(selection, coverage, tax_coverage, suffix='_tax',
                 jaccard_stats=None, gap_closures=None):
    md = REPORTS_DIR / 'TAXONOMY_AWARE_SELECTION.md'
    lines = []
    lines.append('# Taxonomy-Aware Diverse Panel (100 models)\n')
    lines.append('A taxonomy-aware re-selection of 100 representative growers, '
                 'derived from the same 3,461 growers used by '
                 '`scripts/select_diverse.py`. The original panel '
                 'covered network sets well but missed entire phyla and big '
                 'industrial/clinical genera; this panel adds an explicit '
                 'taxonomic stratification pass and a phylogenetic-novelty '
                 'pass that close those gaps.\n')

    lines.append('## Methodology\n')
    lines.append('Seven passes; each pick is tagged with the pass that selected it '
                 '(`pass_origin` column) plus a per-pass reason.\n')
    lines.append('1. **Phylum medoids** — one medoid per grower-bearing phylum '
                 '(model with min sum-Jaccard reaction distance to its phylum-mates; '
                 'sample capped at 50 per phylum). Guarantees every phylum is anchored.')
    lines.append(f'2. **Reaction-coverage core** — greedy max-coverage on '
                 f'`seed.reaction` IDs, seeded with Pass 1. {PASS2_RXN_CORE} picks.')
    lines.append(f'3. **Taxonomic-novelty fill** — greedy farthest-point on '
                 f'lineage rank-distance (7=different superkingdom, 1=different '
                 f'species). Tie-break with Jaccard reaction distance. '
                 f'{PASS3_TAX_NOVELTY} picks; closes class/order gaps.')
    lines.append(f'4. **Metabolite-coverage layer** — greedy max-coverage on '
                 f'cpd IDs, seeded with all prior picks. {PASS4_CPD_LAYER} picks.')
    lines.append(f'5. **Constrained extremes** — 12 axis extremes + top-3 '
                 f'rare-reaction carriers, but skip any candidate whose class '
                 f'is already represented > {EXTREMES_CLASS_MULT}× its expected '
                 f'share (walk down the ranking up to {EXTREMES_WALK_LIMIT} '
                 f'before accepting). Cap {PASS5_EXTREMES}.')
    lines.append(f'6. **Hot-taxon medoids** — fills genus/family/order gaps that '
                 f'Pass 3 misses by sweeping the top grower-heavy missing taxa '
                 f'(threshold {HOT_TAXA_MIN_GROWERS} growers) and adding their '
                 f'medoid. Cap {PASS6_HOT_TAXA}; closes Pseudomonas, '
                 f'Burkholderia, etc. that sit inside already-anchored phyla.')
    lines.append(f'7. **Farthest-point Jaccard fill** — fills the remaining '
                 f'slots by maximizing min reaction-set Jaccard distance to '
                 f'the panel. Pure diversity completion.\n')

    pass_dist = Counter(s['pass_origin'] for s in selection)
    lines.append('## Pass distribution\n')
    lines.append('| Pass | Description | Picks |')
    lines.append('|---|---|---|')
    descs = {
        1: 'phylum medoids', 2: 'reaction coverage',
        3: 'taxonomic novelty', 4: 'metabolite coverage',
        5: 'constrained extremes', 6: 'hot-taxon medoids',
        7: 'farthest-point Jaccard',
    }
    for k in sorted(pass_dist):
        lines.append(f'| {k} | {descs[k]} | {pass_dist[k]} |')
    lines.append('')

    lines.append('## Coverage achieved by the 100-model panel\n')
    lines.append(f'- Reactions covered: **{coverage["reactions_covered"]} / '
                 f'{coverage["reactions_universe"]}** unique `seed.reaction` IDs '
                 f'({100*coverage["reactions_covered"]/coverage["reactions_universe"]:.1f}%)')
    lines.append(f'- Metabolites covered: **{coverage["metabolites_covered"]} / '
                 f'{coverage["metabolites_universe"]}** unique cpd IDs '
                 f'({100*coverage["metabolites_covered"]/coverage["metabolites_universe"]:.1f}%)')
    lines.append('')

    lines.append('## Taxonomic coverage\n')
    lines.append('| Rank | distinct in panel | distinct among growers | coverage |')
    lines.append('|---|---|---|---|')
    for r in LINEAGE_RANKS:
        cov = tax_coverage[f'{r}_covered']
        uni = tax_coverage[f'{r}_universe']
        pct = 100 * cov / uni if uni else 0.0
        lines.append(f'| {r} | {cov} | {uni} | {pct:.1f}% |')
    lines.append('')

    if jaccard_stats:
        lines.append('## Representativeness (Jaccard min-distance to nearest panel member)\n')
        lines.append('| panel | median | mean | p90 | max |')
        lines.append('|---|---|---|---|---|')
        for label, st in jaccard_stats.items():
            lines.append(f'| {label} | {st["median"]:.3f} | {st["mean"]:.3f} | '
                         f'{st["p90"]:.3f} | {st["max"]:.3f} |')
        lines.append('')

    if gap_closures:
        lines.append('## Gap closure vs original panel\n')
        for rank, closed in gap_closures.items():
            if not closed:
                continue
            lines.append(f'### {rank} (newly covered)\n')
            lines.append(f'| {rank} | grower count |')
            lines.append('|---|---|')
            for taxon, n in closed:
                lines.append(f'| {taxon} | {n} |')
            lines.append('')

    lines.append('## Files\n')
    lines.append('- `results/selected_ids_tax.txt` — 100 model IDs')
    lines.append('- `results/selected_models_tax.csv` — per-model metrics, '
                 'lineage, and selection reason')
    lines.append('- `results/selected_models_tax.json` — same data + coverage stats')
    lines.append('- `scripts/select_diverse_tax.py` — reproducible selection')
    lines.append('- `notebooks/08_TaxonomyAwareSelection.ipynb` — interactive walkthrough\n')

    lines.append('## Panel members (in selection order)\n')
    lines.append('| # | model_id | pass | phylum | class | order | genus | reason |')
    lines.append('|---|---|---|---|---|---|---|---|')
    for i, s in enumerate(selection, 1):
        lines.append(f'| {i} | `{s["model_id"]}` | {s["pass_origin"]} | '
                     f'{s.get("phylum","")} | {s.get("class","")} | '
                     f'{s.get("order","")} | {s.get("genus","")} | {s["reason"]} |')
    md.write_text('\n'.join(lines) + '\n')
    print(f'wrote {md}')


# ---------------------------------------------------------------------------
def main():
    print('Loading growers...')
    growers = load_growers()
    print(f'  {len(growers)} growers')

    print('Loading NCBI taxonomy...')
    lineage = load_taxonomy()
    n_resolved = sum(1 for g in growers
                     if lineage.get(g['model_id'], {}).get('phylum') != 'Unknown')
    print(f'  {len(lineage)} accessions in taxonomy table; '
          f'{n_resolved}/{len(growers)} growers have resolved taxonomy')

    print('Extracting reaction & metabolite fingerprints from JSON...')
    rxns_by_id, cpds_by_id = extract_features(growers)

    selection, coverage, tax_coverage, rare_set, prevalence, total_cpds = \
        run_selection(growers, rxns_by_id, cpds_by_id, lineage)

    write_outputs(selection, coverage, tax_coverage)
    write_report(selection, coverage, tax_coverage)


if __name__ == '__main__':
    main()
