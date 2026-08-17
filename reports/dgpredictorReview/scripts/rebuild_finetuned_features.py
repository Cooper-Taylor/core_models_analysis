#!/usr/bin/env python3
"""Rebuild the fine-tuned (ModelSEED) dGPredictor training feature matrix and
refit BayesianRidge, so the per-fingerprint coefficients can be read out.

Why rebuild instead of loading the shipped model?
    freiburgermsu/dGPredictor gitignores ``*/*.pkl`` and ``*/*.mat``, so
    ``model/modelseed_M12_model_BR.pkl`` is NOT in the repository. Only
    ``data/modelseed_training.mat.zip`` (1.2 GB unzipped) ships, which does not
    fit in this machine's RAM.

    Everything needed to reproduce the fit exactly IS in the repo, though:
      * data/component_contribution_python.mat  -> train_S, train_cids, b
      * data/kegg_to_modelseed_compound_map.json
      * data/modelseed_decompose_r{1,2}.json
      * data/modelseed_group_names_r{1,2}.txt
    This script replays retrain_modelseed.py steps 4-5 verbatim, but keeps the
    feature matrix sparse.

The zero-column trick
    BayesianRidge(fit_intercept=False) gives coef_ = V^T D V X^T y. A feature
    column that is identically zero across the training set contributes a zero
    row/column to V and a zero entry to X^T y, so its coefficient is *exactly*
    zero and it has no influence on the coefficients of any other feature.
    Dropping those columns before fitting therefore yields coefficients
    numerically identical to the full-width fit, at a fraction of the memory.

Outputs (results/):
    finetuned_group_coefficients.tsv   one row per group in the full vocabulary
    finetuned_fit_summary.json         fit metrics + vocabulary accounting
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from scipy.sparse import csr_matrix
from sklearn.linear_model import BayesianRidge
from sklearn.metrics import mean_squared_error, r2_score

REPO = Path("/scratch/ctaylor/dgpredictor_repo")
DATA = REPO / "data"
OUT = Path(__file__).resolve().parents[1] / "results"
OUT.mkdir(parents=True, exist_ok=True)

# retrain_modelseed.py: SKIP_IN_FINGERPRINT = {'cpd00067', 'cpd11640'}  (H+, H2)
SKIP_CPDS = {"cpd00067", "cpd11640"}
PAD = 44  # retrain_modelseed.py pads [r1 | 44 zeros | r2 | 44 zeros]


def load_group_names(path: Path) -> list[str]:
    return path.read_text().split("\n")


def main() -> None:
    print("loading decompositions ...", flush=True)
    dec_r1 = json.loads((DATA / "modelseed_decompose_r1.json").read_text())
    dec_r2 = json.loads((DATA / "modelseed_decompose_r2.json").read_text())
    names_r1 = load_group_names(DATA / "modelseed_group_names_r1.txt")
    names_r2 = load_group_names(DATA / "modelseed_group_names_r2.txt")
    idx_r1 = {n: i for i, n in enumerate(names_r1)}
    idx_r2 = {n: i for i, n in enumerate(names_r2)}
    n_r1, n_r2 = len(names_r1), len(names_r2)
    n_feat = n_r1 + PAD + n_r2 + PAD
    print(f"  vocabulary: r1={n_r1}  r2={n_r2}  total feature width={n_feat}")

    kegg2ms = json.loads((DATA / "kegg_to_modelseed_compound_map.json").read_text())

    print("loading component_contribution_python.mat ...", flush=True)
    cc = loadmat(DATA / "component_contribution_python.mat")
    train_S = cc["train_S"]
    cids = [str(c).strip() for c in cc["train_cids"].flatten()]
    b = cc["b"].flatten()

    # ---- retrain_modelseed.load_training_reactions_from_mat -----------------
    reactions, y_all, n_unmapped = [], [], 0
    for j in range(train_S.shape[1]):
        col = train_S[:, j]
        nz = np.nonzero(col)[0]
        if nz.size == 0:
            continue
        stoich, ok = {}, True
        for i in nz:
            ms = kegg2ms.get(cids[i])
            if ms is None:
                ok = False
                break
            stoich[ms] = stoich.get(ms, 0.0) + float(col[i])
        if ok and stoich:
            reactions.append(stoich)
            y_all.append(b[j])
        else:
            n_unmapped += 1
    y_all = np.asarray(y_all)
    print(f"  training reactions mapped to ModelSEED: {len(reactions)} "
          f"({n_unmapped} unmapped)")

    # ---- retrain_modelseed.build_feature_matrix (sparse) --------------------
    decomposable = set(dec_r1) | SKIP_CPDS
    rows, cols, vals, y = [], [], [], []
    n_skipped = 0
    r = 0
    for stoich in reactions:
        if any(c not in decomposable for c in stoich if c not in SKIP_CPDS):
            n_skipped += 1
            r += 1
            continue
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
                    kk = n_r1 + PAD + k
                    acc[kk] = acc.get(kk, 0.0) + cnt * s
        ri = len(y)
        for k, v in acc.items():
            if v != 0.0:
                rows.append(ri)
                cols.append(k)
                vals.append(v)
        y.append(y_all[r])
        r += 1
    y = np.asarray(y)
    X = csr_matrix((vals, (rows, cols)), shape=(len(y), n_feat))
    print(f"  feature matrix: {X.shape}  nnz={X.nnz}  ({n_skipped} reactions skipped)")

    # ---- drop identically-zero columns --------------------------------------
    used = np.asarray((X != 0).sum(axis=0)).ravel() > 0
    used_idx = np.flatnonzero(used)
    Xd = np.asarray(X[:, used_idx].todense())
    print(f"  columns ever non-zero in training: {used_idx.size} / {n_feat} "
          f"({100.0 * used_idx.size / n_feat:.2f}%)")
    print(f"  dense reduced matrix: {Xd.shape} = {Xd.nbytes / 1e6:.0f} MB", flush=True)

    # ---- refit ---------------------------------------------------------------
    print("fitting BayesianRidge(tol=1e-6, fit_intercept=False) ...", flush=True)
    model = BayesianRidge(tol=1e-6, fit_intercept=False, compute_score=True)
    model.fit(Xd, y)
    pred = model.predict(Xd)
    mse = float(mean_squared_error(y, pred))
    r2 = float(r2_score(y, pred))
    print(f"  MSE={mse:.2f} kJ^2/mol^2   R2={r2:.4f}")

    coef_full = np.zeros(n_feat)
    coef_full[used_idx] = model.coef_

    # ---- write out -----------------------------------------------------------
    # how many distinct training compounds carry each group
    train_cpds = sorted({c for rx in reactions for c in rx})
    supp_r1: dict[str, int] = {}
    supp_r2: dict[str, int] = {}
    for c in train_cpds:
        for g in dec_r1.get(c, {}):
            supp_r1[g] = supp_r1.get(g, 0) + 1
        for g in dec_r2.get(c, {}):
            supp_r2[g] = supp_r2.get(g, 0) + 1

    lines = ["radius\tgroup\tfeature_index\tcoefficient_kJ_per_mol\t"
             "used_in_training\tn_training_compounds\tn_all_modelseed_compounds"]
    all_r1: dict[str, int] = {}
    all_r2: dict[str, int] = {}
    for d, acc in ((dec_r1, all_r1), (dec_r2, all_r2)):
        for gs in d.values():
            for g in gs:
                acc[g] = acc.get(g, 0) + 1
    for rad, names, off, supp, alln in ((1, names_r1, 0, supp_r1, all_r1),
                                        (2, names_r2, n_r1 + PAD, supp_r2, all_r2)):
        for i, name in enumerate(names):
            fi = off + i
            lines.append(f"{rad}\t{name}\t{fi}\t{float(coef_full[fi]):.10g}\t"
                         f"{int(used[fi])}\t{supp.get(name, 0)}\t{alln.get(name, 0)}")
    (OUT / "finetuned_group_coefficients.tsv").write_text("\n".join(lines) + "\n")

    summary = {
        "source_repo": str(REPO),
        "n_training_reactions_in_cc_mat": int(train_S.shape[1]),
        "n_training_reactions_mapped": len(reactions),
        "n_training_reactions_featurized": int(len(y)),
        "n_training_compounds": len(train_cpds),
        "vocabulary": {
            "r1": n_r1, "r2": n_r2, "feature_width": n_feat,
            "r1_used_in_training": int(used[:n_r1].sum()),
            "r2_used_in_training": int(used[n_r1 + PAD:n_r1 + PAD + n_r2].sum()),
            "total_used_in_training": int(used_idx.size),
        },
        "fit": {"mse_kJ2": mse, "r2": r2,
                "alpha_": float(model.alpha_), "lambda_": float(model.lambda_)},
        "coefficients": {
            "n_nonzero": int((coef_full != 0).sum()),
            "max_abs_kJ": float(np.abs(coef_full).max()),
        },
    }
    (OUT / "finetuned_fit_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    sys.exit(main())
