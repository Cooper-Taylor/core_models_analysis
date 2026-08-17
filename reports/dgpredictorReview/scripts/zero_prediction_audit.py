#!/usr/bin/env python3
"""Audit the reactions the fine-tuned dGPredictor returns dG_model_only = 0 for.

There are two very different reasons a group-contribution model can return zero,
and the shipped prediction file does not distinguish them:

  BLIND   the reaction's group-change vector is identically zero -- every
          fragment created is also destroyed. Isomerisations, racemisations,
          some intramolecular rearrangements. The decomposition literally cannot
          see the reaction. Radius-2 fragments help but do not cure this.

  UNSEEN  the group-change vector is non-zero, but every fragment it touches has
          coefficient exactly zero because no training reaction ever contained
          it. The model has no information and says so by predicting zero.

Both come back as "0.00 +/- 3.31 kJ/mol", i.e. tagged with sqrt(1/alpha_), the
smallest uncertainty the model can emit. That is the failure mode worth naming:
maximum stated confidence on the reactions the model understands least.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from scipy.sparse import csr_matrix
from sklearn.linear_model import BayesianRidge

DATA = Path("/scratch/ctaylor/dgpredictor_repo/data")
OUT = Path(__file__).resolve().parents[1] / "results"
SKIP_CPDS = {"cpd00067", "cpd11640"}
PAD = 44


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
    S, b = cc["train_S"], cc["b"].flatten()
    cids = [str(c).strip() for c in cc["train_cids"].flatten()]

    def vec(st):
        acc = {}
        for cpd, s in st.items():
            if cpd in SKIP_CPDS:
                continue
            for g, c in dec_r1.get(cpd, {}).items():
                k = idx_r1.get(g)
                if k is not None:
                    acc[k] = acc.get(k, 0.0) + c * s
            for g, c in dec_r2.get(cpd, {}).items():
                k = idx_r2.get(g)
                if k is not None:
                    acc[n_r1 + PAD + k] = acc.get(n_r1 + PAD + k, 0.0) + c * s
        return {k: v for k, v in acc.items() if abs(v) > 1e-12}

    rxns, ys = [], []
    for j in range(S.shape[1]):
        nz = np.nonzero(S[:, j])[0]
        if nz.size == 0:
            continue
        st, ok = {}, True
        for i in nz:
            ms = kegg2ms.get(cids[i])
            if ms is None:
                ok = False
                break
            st[ms] = st.get(ms, 0.0) + float(S[i, j])
        if ok and st:
            rxns.append(st); ys.append(b[j])

    decomposable = set(dec_r1) | SKIP_CPDS
    rows, cols, vals, y = [], [], [], []
    for st, yy in zip(rxns, ys):
        if any(c not in decomposable for c in st if c not in SKIP_CPDS):
            continue
        acc = vec(st)
        ri = len(y)
        for k, v in acc.items():
            rows.append(ri); cols.append(k); vals.append(v)
        y.append(yy)
    X = csr_matrix((vals, (rows, cols)), shape=(len(y), n_feat))
    used = np.flatnonzero(np.asarray((X != 0).sum(axis=0)).ravel() > 0)
    m = BayesianRidge(tol=1e-6, fit_intercept=False, compute_score=True)
    m.fit(np.asarray(X[:, used].todense()), np.asarray(y))
    coef = np.zeros(n_feat)
    coef[used] = m.coef_
    zero = coef == 0.0
    floor = float(np.sqrt(1.0 / m.alpha_))

    preds = json.loads((DATA / "modelseed_all_reaction_dG_predictions.json").read_text())
    stoich = json.loads((DATA / "modelseed_reaction_stoich.json").read_text())

    blind, unseen, mixed, other = [], [], [], []
    for rxn, rec in preds.items():
        if abs(rec["dG_model_only"]) > 1e-9:
            continue
        st = stoich.get(rxn)
        if st is None:
            other.append(rxn); continue
        acc = vec(st)
        if not acc:
            blind.append((rxn, rec["dG_std"]))
        else:
            k = np.fromiter(acc.keys(), int)
            if zero[k].all():
                unseen.append((rxn, rec["dG_std"]))
            else:
                mixed.append((rxn, rec["dG_std"]))

    def stat(rows_):
        if not rows_:
            return {"n": 0}
        s = np.array([r[1] for r in rows_])
        return {"n": len(rows_),
                "n_at_sigma_floor": int(np.isclose(s, floor, atol=1e-3).sum()),
                "median_sigma_kJ": round(float(np.median(s)), 2),
                "examples": [r[0] for r in rows_[:8]]}

    res = {
        "sigma_floor_kJ": round(floor, 3),
        "n_predictions_total": len(preds),
        "n_with_dG_model_only_exactly_zero":
            int(sum(1 for r in preds.values() if abs(r["dG_model_only"]) < 1e-9)),
        "BLIND_group_change_vector_is_identically_zero": stat(blind),
        "UNSEEN_all_touched_groups_have_zero_coefficient": stat(unseen),
        "MIXED_nonzero_vector_that_happens_to_cancel_in_energy": stat(mixed),
        "no_stoichiometry_record": len(other),
    }
    (OUT / "zero_prediction_audit.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
