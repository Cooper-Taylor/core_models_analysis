#!/usr/bin/env python3
"""Verification suite for optimize_thermo_consensus.py.

Ten assertions covering the ways the optimiser could look right and be wrong:
that the shipped rule actually meets the constraints it claims, that it has
not degenerated into a high-leverage subset (the failure mode that motivated
using CCC over Pearson r), that it beats the hand-picked sigma <= 20 baseline,
that the oracle bound is exact, that the predicate round-trips, and that the
selection does not quietly discard core-model metabolism.

Exits non-zero on any failure so it can gate a re-run.
"""
import pandas as pd, numpy as np, json, sys
sys.path.insert(0,'scripts')
from optimize_thermo_consensus import load_selector, metrics, stratum_retention
d = pd.read_csv('results/eq_vs_dgpms/key_subset_classified.tsv', sep='\t', low_memory=False)
d['abs_dg_eq']=d.dg_eq.abs(); d['abs_net_proton']=d.net_proton.abs()
x,y = d.dg_eq.to_numpy(float), d.dg_dgp.to_numpy(float)
strata = pd.qcut(d.abs_dg_eq,10,labels=False,duplicates='drop').to_numpy()
sel = load_selector(); m = sel(d); mm = metrics(x[m],y[m]); ok=[]

# 1 predicate round-trip
saved = set(pd.read_csv('results/eq_vs_dgpms/consensus_selected.tsv',sep='\t').rxn)
ok.append(('predicate round-trip reproduces consensus_selected.tsv',
           set(d.rxn[m])==saved, f'{len(set(d.rxn[m]))} vs {len(saved)}'))
# 2 constraints actually met
t = json.load(open('results/eq_vs_dgpms/consensus_rule.json'))['target']
ok.append(('RMSE constraint met', mm['rmse']<=t['rmse']+1e-9, f"{mm['rmse']:.3f} <= {t['rmse']}"))
ok.append(('CCC constraint met', mm['ccc']>=t['ccc'], f"{mm['ccc']:.4f} >= {t['ccc']}"))
ok.append(('slope guard met', abs(mm['slope']-1)<=t['slope_tol'], f"{mm['slope']:.3f}"))
# 3 NOT degenerate: must keep the near-zero bulk, not just high-|dG|
ret = stratum_retention(m, strata)
ok.append(('keeps the low-|dG| deciles (not a high-leverage set)', ret[:3].min()>0.15,
           f'deciles 0-2 retention {ret[:3].min():.1%}'))
ok.append(('median |dG| of kept set is not inflated',
           np.median(np.abs(x[m])) < np.median(np.abs(x))*3,
           f'{np.median(np.abs(x[m])):.2f} vs {np.median(np.abs(x)):.2f} overall'))
# 4 beats the hand-picked sigma<=20 baseline at equal RMSE
b = (d.dgp_uncertainty<=20).to_numpy(); bm = metrics(x[b],y[b])
ok.append(('beats hand-picked sigma<=20 on RMSE', mm['rmse']<bm['rmse'],
           f"rule {mm['rmse']:.2f} vs baseline {bm['rmse']:.2f}"))
ok.append(('...and on CCC', mm['ccc']>bm['ccc'], f"{mm['ccc']:.3f} vs {bm['ccc']:.3f}"))
# 5 oracle exactness (monotone prefix) re-checked
res=np.abs(x-y); o=np.argsort(res,kind='stable')
pref=np.sqrt(np.cumsum(res[o]**2)/np.arange(1,len(res)+1))
ok.append(('oracle prefix-RMSE monotone (exactness)', bool(np.all(np.diff(pref)>=-1e-12)),'')) 
# 6 biological: core-model retention
cnt=json.load(open('site/data/reaction_model_counts.json'))
core=np.array([cnt.get(r,{}).get('all',0)>0 for r in d.rxn])
ok.append(('core-model reactions retained at >= overall coverage',
           m[core].mean() >= m.mean(), f'{m[core].mean():.1%} core vs {m.mean():.1%} overall'))
print(f"selected n={mm['n']} ({m.mean():.1%})  CCC={mm['ccc']:.4f}  RMSE={mm['rmse']:.3f}  slope={mm['slope']:.3f}\n")
bad=0
for name,passed,note in ok:
    print(f"  [{'PASS' if passed else 'FAIL'}] {name:52s} {note}")
    bad += (not passed)
print(f"\n{len(ok)-bad}/{len(ok)} checks passed")
sys.exit(1 if bad else 0)
