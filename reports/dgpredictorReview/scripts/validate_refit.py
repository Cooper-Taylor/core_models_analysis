#!/usr/bin/env python3
"""Check that the locally refit BayesianRidge reproduces the fine-tuned model's
shipped predictions, then measure how much of each predicted reaction's group
change falls on groups the model never saw (coefficient exactly zero).

Reference: data/modelseed_all_reaction_dG_predictions.json in the repo, which
carries dG_model_only (the raw BayesianRidge output, before the Legendre pH
correction) for every reaction the fine-tuned model could predict.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from scipy.sparse import csr_matrix
from sklearn.linear_model import BayesianRidge

REPO = Path("/scratch/ctaylor/dgpredictor_repo")
DATA = REPO / "data"
OUT = Path(__file__).resolve().parents[1] / "results"

SKIP_CPDS = {"cpd00067", "cpd11640"}
PAD = 44


def build_vec(stoich, dec_r1, dec_r2, idx_r1, idx_r2, n_r1):
    acc: dict[int, float] = {}
    for cpd, s in stoich.items():
        if cpd in SKIP_CPDS:
            continue
        for grp, cnt in dec_r1.get(cpd, {}).items():
            k = idx_r1.get(grp)
            if k is not None:
                acc[k] = acc.get(k, 0.0) + cnt * s
        for grp, cnt in dec_r2.get(cpd, {}).items():
            k = idx_r2.get(grp)
            if k is not None:
                acc[n_r1 + PAD + k] = acc.get(n_r1 + PAD + k, 0.0) + cnt * s
    return acc


def main() -> None:
    dec_r1 = json.loads((DATA / "modelseed_decompose_r1.json").read_text())
    dec_r2 = json.loads((DATA / "modelseed_decompose_r2.json").read_text())
    names_r1 = (DATA / "modelseed_group_names_r1.txt").read_text().split("\n")
    names_r2 = (DATA / "modelseed_group_names_r2.txt").read_text().split("\n")
    idx_r1 = {n: i for i, n in enumerate(names_r1)}
    idx_r2 = {n: i for i, n in enumerate(names_r2)}
    n_r1, n_r2 = len(names_r1), len(names_r2)
    n_feat = n_r1 + PAD + n_r2 + PAD

    kegg2ms = json.loads((DATA / "kegg_to_modelseed_compound_map.json").read_text())
    cc = loadmat(DATA / "component_contribution_python.mat")
    train_S, b = cc["train_S"], cc["b"].flatten()
    cids = [str(c).strip() for c in cc["train_cids"].flatten()]

    reactions, y_all = [], []
    for j in range(train_S.shape[1]):
        col = train_S[:, j]
        nz = np.nonzero(col)[0]
        if nz.size == 0:
            continue
        st, ok = {}, True
        for i in nz:
            ms = kegg2ms.get(cids[i])
            if ms is None:
                ok = False
                break
            st[ms] = st.get(ms, 0.0) + float(col[i])
        if ok and st:
            reactions.append(st)
            y_all.append(b[j])

    decomposable = set(dec_r1) | SKIP_CPDS
    rows, cols, vals, y = [], [], [], []
    for st, yy in zip(reactions, y_all):
        if any(c not in decomposable for c in st if c not in SKIP_CPDS):
            continue
        acc = build_vec(st, dec_r1, dec_r2, idx_r1, idx_r2, n_r1)
        ri = len(y)
        for k, v in acc.items():
            if v != 0.0:
                rows.append(ri); cols.append(k); vals.append(v)
        y.append(yy)
    X = csr_matrix((vals, (rows, cols)), shape=(len(y), n_feat))
    used_idx = np.flatnonzero(np.asarray((X != 0).sum(axis=0)).ravel() > 0)
    model = BayesianRidge(tol=1e-6, fit_intercept=False, compute_score=True)
    model.fit(np.asarray(X[:, used_idx].todense()), np.asarray(y))
    coef = np.zeros(n_feat)
    coef[used_idx] = model.coef_
    zero_mask = coef == 0.0

    # ---- compare against shipped predictions --------------------------------
    preds = json.loads((DATA / "modelseed_all_reaction_dG_predictions.json").read_text())
    stoich_all = json.loads((DATA / "modelseed_reaction_stoich.json").read_text())
    if isinstance(preds, dict) and "reactions" in preds:
        preds = preds["reactions"]

    diffs, mine, theirs = [], [], []
    zero_share, n_unseen_groups = [], []
    for rxn, rec in preds.items():
        ref = rec.get("dG_model_only")
        st = stoich_all.get(rxn)
        if ref is None or st is None:
            continue
        acc = build_vec(st, dec_r1, dec_r2, idx_r1, idx_r2, n_r1)
        if not acc:
            continue
        k = np.fromiter(acc.keys(), dtype=int)
        v = np.fromiter(acc.values(), dtype=float)
        pred = float(np.dot(coef[k], v))
        diffs.append(pred - float(ref)); mine.append(pred); theirs.append(float(ref))
        tot = np.abs(v).sum()
        if tot > 0:
            zero_share.append(float(np.abs(v[zero_mask[k]]).sum() / tot))
            n_unseen_groups.append(int(zero_mask[k].sum()))

    diffs = np.asarray(diffs); mine = np.asarray(mine); theirs = np.asarray(theirs)
    zero_share = np.asarray(zero_share); n_unseen_groups = np.asarray(n_unseen_groups)

    res = {
        "n_reactions_compared": int(diffs.size),
        "refit_vs_shipped_dG_model_only": {
            "pearson_r": float(np.corrcoef(mine, theirs)[0, 1]),
            "median_abs_diff_kJ": float(np.median(np.abs(diffs))),
            "p95_abs_diff_kJ": float(np.percentile(np.abs(diffs), 95)),
            "max_abs_diff_kJ": float(np.abs(diffs).max()),
            "frac_within_1_kJ": float((np.abs(diffs) < 1).mean()),
        },
        "unseen_group_exposure": {
            "frac_reactions_touching_a_zero_coef_group":
                float((n_unseen_groups > 0).mean()),
            "median_frac_of_|group change| on zero-coef groups":
                float(np.median(zero_share)),
            "mean_frac_of_|group change| on zero-coef groups":
                float(zero_share.mean()),
            "frac_reactions_with_>50pct_on_zero_coef":
                float((zero_share > 0.5).mean()),
            "frac_reactions_with_100pct_on_zero_coef":
                float((zero_share > 0.999).mean()),
        },
    }
    (OUT / "refit_validation.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
