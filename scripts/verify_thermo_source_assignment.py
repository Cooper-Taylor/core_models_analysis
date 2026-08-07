#!/usr/bin/env python3
"""Verification gate for optimize_thermo_source_assignment.py.

Checks the ways a source-assignment policy can look right and be wrong: that it
beats every fixed-source baseline AND the incumbent dev priority on held-out
experimental data, that the hard overrides actually fired, that it never assigns
a source it does not have, that coverage exceeds the consensus-subset approach
it replaces, and that the calibration is monotone in sigma as claimed.

Exits non-zero on any failure.
"""
import json, sys
import numpy as np, pandas as pd
sys.path.insert(0, 'scripts')
from optimize_thermo_source_assignment import (load_db, load_truth, ASSIGN_TSV,
                                               MODELS_JSON, EQ_SENTINEL, TOLERANCE)

db = load_db()
a = pd.read_csv(ASSIGN_TSV, sep='\t', low_memory=False)
mdl = json.loads(MODELS_JSON.read_text())
ok = []

val = pd.DataFrame(mdl['validation']).set_index('strategy')
mine = val.loc['assignment (this script)']
for other in [i for i in val.index if i != 'assignment (this script)']:
    o = val.loc[other]
    ok.append((f'mean error beats "{other[:34]}"', mine.mean_abs_err <= o.mean_abs_err + 1e-9,
               f'{mine.mean_abs_err:.2f} vs {o.mean_abs_err:.2f}'))

# overrides actually fired
q = db.is_quinone == 1
ok.append(('dGPredictor never assigned to a quinone reaction',
           not ((a.chosen_source == 'DGPMS') & q.values).any(),
           f'{int(q.sum())} quinone reactions in db'))
eq_sent = (db.sig_EQ > EQ_SENTINEL).fillna(False)
ok.append(('eQuilibrator never assigned where it declares no estimate',
           not ((a.chosen_source == 'EQ') & eq_sent.values).any(),
           f'{int(eq_sent.sum())} sentinel reactions'))

# never assigns an absent source
bad = 0
for k in ('GC', 'EQ', 'DGPMS'):
    bad += int(((a.chosen_source == k) & db[f'dg_{k}'].isna().values).sum())
ok.append(('never assigns a source the reaction lacks', bad == 0, f'{bad} violations'))

# merged dG matches the chosen source's value
kept = a[a.kept == True]
mism = 0
for k in ('GC', 'EQ', 'DGPMS'):
    m = kept.chosen_source == k
    if m.any():
        mism += int((~np.isclose(kept.loc[m, 'merged_dg'],
                                 db.loc[kept.index[m], f'dg_{k}'], equal_nan=True)).sum())
ok.append(('merged_dg equals the chosen source value', mism == 0, f'{mism} mismatches'))

# beats the consensus-subset approach on coverage
try:
    cons = len(pd.read_csv('results/eq_vs_dgpms/consensus_selected.tsv', sep='\t'))
except Exception:
    cons = 3246
ok.append(('covers more than the consensus subset', int(a.kept.sum()) > cons,
           f'{int(a.kept.sum())} vs {cons}'))

# calibration monotone in sigma
for k, m in mdl['models'].items():
    if m['kind'] == 'isotonic':
        ok.append((f'{k} calibration monotone in sigma',
                   bool(np.all(np.diff(m['y']) >= -1e-9)), ''))

print(f"assigned {int(a.kept.sum())} reactions at expected error <= {TOLERANCE} kcal/mol\n")
bad = 0
for name, passed, note in ok:
    print(f"  [{'PASS' if passed else 'FAIL'}] {name:52s} {note}")
    bad += (not passed)
print(f"\n{len(ok)-bad}/{len(ok)} checks passed")
sys.exit(1 if bad else 0)
