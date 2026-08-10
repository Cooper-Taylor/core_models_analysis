#!/usr/bin/env python3
"""How close is the PREDICTED error (ehat) to the TRUE error?

optimize_thermo_source_assignment.py validates that the ASSIGNMENT beats the
baselines. This asks the narrower, more sceptical question: taken on its own, is
ehat a good predictor of a source's actual error?

Repeated 70/30 splits of the TECRDB gold set, refitting per split, with held-out
reactions excluded from the silver tier too so nothing leaks. Reports:
  * conservatism    -- median predicted vs median actual
  * P(actual <= ehat)  -- how often it is a genuine upper bound
  * rho(ehat, actual)  -- whether it RANKS reactions by error
  * promise-keeping -- P(actual <= 2 | ehat <= 2), the number to quote

Headline: ehat is a conservative THRESHOLD (promise kept 87-91% of the time) but
a weak RANKING (rho ~0.11-0.24) inside the TECRDB range -- because TECRDB holds
only well-measured central metabolism and none of the extreme reactions the
filter exists to exclude. See section 5b of THERMO_SOURCE_ASSIGNMENT.md.
"""
import sys, numpy as np, pandas as pd
sys.path.insert(0,'scripts')
from optimize_thermo_source_assignment import (load_db, load_truth, fit_error_models,
                                               predict_error, EQ_SENTINEL)
db = load_db(); truth = load_truth(db)
rng = np.random.default_rng(11)
rows=[]
for rep in range(20):                      # repeated splits, so this isn't one lucky draw
    idx = rng.permutation(len(truth)); cut = int(len(truth)*0.7)
    tr, te = truth.iloc[idx[:cut]], truth.iloc[idx[cut:]]
    mdl = fit_error_models(tr, db[~db.rxn.isin(te.rxn)])
    eh = predict_error(te, mdl)
    for k in ('GC','EQ','DGPMS'):
        m = te[f'dg_{k}'].notna() & eh[f'ehat_{k}'].notna()
        if m.sum() < 20: continue
        pred = eh.loc[m, f'ehat_{k}'].to_numpy(float)
        act  = (te.loc[m, f'dg_{k}'] - te.loc[m,'tecrdb']).abs().to_numpy(float)
        rows.append({'rep':rep,'source':k,'n':int(m.sum()),
                     'spearman':pd.Series(pred).corr(pd.Series(act),method='spearman'),
                     'median_pred':np.median(pred),'median_act':np.median(act),
                     'cover_pred_ge_act':float((act<=pred).mean()),
                     'p_act_le_2_given_pred_le_2': float((act[pred<=2]<=2).mean()) if (pred<=2).any() else np.nan,
                     'n_pred_le_2': int((pred<=2).sum())})
r = pd.DataFrame(rows)
print('Predicted vs TRUE error on held-out TECRDB, mean over 20 random splits:\n')
print(f'{"source":8s} {"n/split":>8} {"rho(pred,act)":>14} {"med pred":>9} {"med actual":>11} {"P(act<=pred)":>13}')
for k,g in r.groupby('source'):
    print(f'{k:8s} {g.n.mean():8.0f} {g.spearman.mean():14.3f} {g.median_pred.mean():9.2f} '
          f'{g.median_act.mean():11.2f} {g.cover_pred_ge_act.mean():13.1%}')
print('\nTHE PROMISE-KEEPING RATE -- when ehat says "<= 2 kcal/mol", how often is it true?')
for k,g in r.groupby('source'):
    print(f'  {k:8s} n~{g.n_pred_le_2.mean():5.0f}/split   actual error <= 2 in '
          f'{g.p_act_le_2_given_pred_le_2.mean():.1%} of cases')
