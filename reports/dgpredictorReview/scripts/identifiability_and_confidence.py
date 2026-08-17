#!/usr/bin/env python3
"""Two structural questions about the fine-tuned dGPredictor's fingerprint basis,
answered against the three group vocabularies.

1. IDENTIFIABILITY -- how many independent degrees of freedom does each
   decomposition actually give the regression over its own training set?
   Reported as rank(G) / #groups / #observations.

2. WHAT HAPPENS OUTSIDE THE TRAINED SPAN -- component-contribution splits the
   training span into a reactant-contribution part, a group-contribution part,
   and an orthogonal remainder it assigns RMSE_inf (a deliberately enormous
   uncertainty). BayesianRidge has no such construct: a feature direction never
   seen in training gets coefficient exactly zero AND contributes nothing to the
   posterior variance, so the model reports its *smallest* uncertainty,
   sqrt(1/alpha_), precisely where it knows least. This script measures how
   often that happens on real ModelSEED reactions.
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
        nz = np.nonzero(train_S[:, j])[0]
        if nz.size == 0:
            continue
        st, ok = {}, True
        for i in nz:
            ms = kegg2ms.get(cids[i])
            if ms is None:
                ok = False
                break
            st[ms] = st.get(ms, 0.0) + float(train_S[i, j])
        if ok and st:
            reactions.append(st)
            y_all.append(b[j])

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
                    kk = n_r1 + PAD + k
                    acc[kk] = acc.get(kk, 0.0) + c * s
        return acc

    decomposable = set(dec_r1) | SKIP_CPDS
    rows, cols, vals, y = [], [], [], []
    train_cpds = set()
    for st, yy in zip(reactions, y_all):
        if any(c not in decomposable for c in st if c not in SKIP_CPDS):
            continue
        acc = vec(st)
        ri = len(y)
        for k, v in acc.items():
            if v != 0.0:
                rows.append(ri); cols.append(k); vals.append(v)
        y.append(yy)
        train_cpds |= set(st)
    X = csr_matrix((vals, (rows, cols)), shape=(len(y), n_feat))
    used_idx = np.flatnonzero(np.asarray((X != 0).sum(axis=0)).ravel() > 0)
    Xd = np.asarray(X[:, used_idx].todense())

    model = BayesianRidge(tol=1e-6, fit_intercept=False, compute_score=True)
    model.fit(Xd, np.asarray(y))
    coef = np.zeros(n_feat)
    coef[used_idx] = model.coef_
    zero_mask = coef == 0.0
    sigma_floor = float(np.sqrt(1.0 / model.alpha_))

    # ---- 1. identifiability -------------------------------------------------
    # compound x learned-group incidence, over the compounds in the training set
    tc = sorted(c for c in train_cpds if c not in SKIP_CPDS)
    Gd = np.zeros((len(tc), used_idx.size))
    pos = {f: i for i, f in enumerate(used_idx)}
    for r, cpd in enumerate(tc):
        for g, c in dec_r1.get(cpd, {}).items():
            k = idx_r1.get(g)
            if k is not None and k in pos:
                Gd[r, pos[k]] += c
        for g, c in dec_r2.get(cpd, {}).items():
            k = idx_r2.get(g)
            if k is not None and (n_r1 + PAD + k) in pos:
                Gd[r, pos[n_r1 + PAD + k]] += c

    identifiability = {
        "finetuned_dgpredictor": {
            "n_training_observations": int(len(y)),
            "n_training_compounds": len(tc),
            "n_groups_declared": n_r1 + n_r2,
            "n_groups_with_any_training_support": int(used_idx.size),
            "rank_of_compound_x_group_matrix": int(np.linalg.matrix_rank(Gd)),
            "rank_of_reaction_feature_matrix": int(np.linalg.matrix_rank(Xd)),
            "observations_per_free_parameter":
                round(len(y) / int(np.linalg.matrix_rank(Xd)), 2),
        },
    }

    # ---- 2. behaviour outside the trained span ------------------------------
    preds = json.loads((DATA / "modelseed_all_reaction_dG_predictions.json").read_text())
    stoich_all = json.loads((DATA / "modelseed_reaction_stoich.json").read_text())

    share, sig, dgm = [], [], []
    for rxn, rec in preds.items():
        st = stoich_all.get(rxn)
        if st is None or rec.get("dG_std") is None:
            continue
        acc = vec(st)
        if not acc:
            continue
        k = np.fromiter(acc.keys(), int)
        v = np.fromiter(acc.values(), float)
        tot = np.abs(v).sum()
        if tot <= 0:
            continue
        share.append(float(np.abs(v[zero_mask[k]]).sum() / tot))
        sig.append(float(rec["dG_std"]))
        dgm.append(float(rec["dG_model_only"]))
    share, sig, dgm = np.asarray(share), np.asarray(sig), np.asarray(dgm)

    bins = [(0.0, 1e-9), (1e-9, 0.1), (0.1, 0.25), (0.25, 0.5), (0.5, 0.999), (0.999, 1.01)]
    tiers = []
    for lo, hi in bins:
        m = (share >= lo) & (share < hi)
        if not m.any():
            continue
        tiers.append({
            "unseen_group_share": f"[{lo:g}, {hi:g})",
            "n_reactions": int(m.sum()),
            "median_reported_sigma_kJ": round(float(np.median(sig[m])), 2),
            "median_abs_dG_model_kJ": round(float(np.median(np.abs(dgm[m]))), 2),
        })

    outside = {
        "bayesianridge_sigma_floor_kJ": round(sigma_floor, 3),
        "note": ("sqrt(1/alpha_): the uncertainty BayesianRidge reports when a "
                 "reaction's feature vector lies entirely outside the trained "
                 "span, i.e. its most confident possible answer"),
        "spearman_rho_unseen_share_vs_reported_sigma":
            round(float(_spearman(share, sig)), 4),
        "tiers": tiers,
        "n_reactions_predicted_exactly_zero_because_all_groups_unseen":
            int(((share > 0.999) & (np.abs(dgm) < 1e-9)).sum()),
    }

    res = {"identifiability": identifiability, "outside_trained_span": outside}
    (OUT / "identifiability_and_confidence.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    return float(np.corrcoef(ra, rb)[0, 1])


if __name__ == "__main__":
    main()
