#!/usr/bin/env python3
"""How much of each predicted reaction lies OUTSIDE the span of the training set.

Component-contribution answers this question explicitly: it decomposes every
query into a reactant-contribution component, a group-contribution component,
and an orthogonal remainder, and assigns the remainder RMSE_inf -- an
intentionally enormous uncertainty that says "this part of your reaction is
outside anything I was trained on". The projection matrices for that split
(P_R_rc, P_N_rc, P_R_gc, P_N_gc) are literally stored in the parameter file.

dGPredictor's BayesianRidge has no equivalent. This script computes the same
quantity for it: the fraction of each reaction's group-change vector that lies
in the null space of the training feature matrix, i.e. the part of the reaction
for which the reported coefficients are the L2 prior talking, not the data.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from scipy.sparse import csr_matrix

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

    kegg2ms = json.loads((DATA / "kegg_to_modelseed_compound_map.json").read_text())
    cc = loadmat(DATA / "component_contribution_python.mat")
    S = cc["train_S"]
    cids = [str(c).strip() for c in cc["train_cids"].flatten()]
    decomposable = set(dec_r1) | SKIP_CPDS

    rows, cols, vals, nrow = [], [], [], 0
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
        if not ok or not st:
            continue
        if any(c not in decomposable for c in st if c not in SKIP_CPDS):
            continue
        for k, v in vec(st).items():
            rows.append(nrow); cols.append(k); vals.append(v)
        nrow += 1
    X = csr_matrix((vals, (rows, cols)), shape=(nrow, n_feat))
    used = np.flatnonzero(np.asarray((X != 0).sum(axis=0)).ravel() > 0)
    Xd = np.asarray(X[:, used].todense())

    # orthonormal basis of the training ROW space
    _, sv, Vt = np.linalg.svd(Xd, full_matrices=False)
    tol = max(Xd.shape) * np.finfo(float).eps * sv[0]
    rank = int((sv > tol).sum())
    V = Vt[:rank]                       # (rank, n_used)
    pos = {int(f): i for i, f in enumerate(used)}
    print(f"training feature matrix {Xd.shape}, rank {rank}")

    preds = json.loads((DATA / "modelseed_all_reaction_dG_predictions.json").read_text())
    stoich = json.loads((DATA / "modelseed_reaction_stoich.json").read_text())

    shares, ids = [], []
    for rxn in preds:
        st = stoich.get(rxn)
        if st is None:
            continue
        acc = vec(st)
        if not acc:
            continue
        x = np.zeros(used.size)
        outside_vocab = 0.0
        for k, v in acc.items():
            i = pos.get(k)
            if i is None:
                outside_vocab += v * v      # group with no training support at all
            else:
                x[i] += v
        n2 = float(x @ x) + outside_vocab
        if n2 <= 0:
            continue
        proj = V @ x
        in2 = float(proj @ proj)
        shares.append(float(np.sqrt(max(n2 - in2, 0.0) / n2)))
        ids.append(rxn)
    shares = np.asarray(shares)

    res = {
        "training_feature_matrix": {"shape": list(Xd.shape), "rank": rank,
                                    "n_nonzero_coefficients": int(used.size),
                                    "n_prior_determined_directions": int(used.size - rank)},
        "n_reactions_scored": int(shares.size),
        "out_of_span_norm_fraction": {
            "median": round(float(np.median(shares)), 4),
            "mean": round(float(shares.mean()), 4),
            "frac_below_0.1": round(float((shares < 0.1).mean()), 4),
            "frac_above_0.5": round(float((shares > 0.5).mean()), 4),
            "frac_above_0.9": round(float((shares > 0.9).mean()), 4),
        },
    }
    (OUT / "rowspace_coverage.json").write_text(json.dumps(res, indent=2))
    np.save(OUT / "out_of_span_share.npy", shares)
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
