#!/usr/bin/env python3
"""Test *why* dGPredictor decorrelates from Group Contribution / eQuilibrator on
some reaction families and tracks them on others.

Three experiments, all read-only over ModelSEEDDatabase:

E1. ADDITIVITY. Group Contribution and eQuilibrator both store per-compound
    formation energies, so their reaction DeltaG is exactly ``sum(nu_i * dGf_i)``.
    dGPredictor has no compound-level energies -- it regresses the reaction
    directly. Fit implied per-compound energies to each source by least squares
    over the stoichiometry matrix, with a held-out split, and compare
    out-of-sample R^2. High R^2 => the source behaves like a sum over
    compounds; low R^2 => it does not.

E2. SIGNATURE NULL SPACE. dGPredictor represents a reaction by the *change* in
    atom-centred molecular signatures (Morgan/ECFP-style environments),
    ``dSig = sum(nu_i * sig(compound_i))``. Any model of that form must return
    the same DeltaG for two reactions with the same dSig. Group reactions into
    exact-dSig equivalence classes and decompose the Group-Contribution DeltaG
    variance into between-class (learnable) and within-class (structurally
    invisible). The within-class share is an upper bound on how well *any*
    signature-difference model can do -- independent of training data.

E3. FAMILY BREAKDOWN. Repeat E2 per reaction family so the redox-vs-group-
    transfer split seen in the correlation scan can be attributed.

Writes results/thermo_agreement/{additivity_fit.tsv, signature_classes.tsv,
nullspace_by_family.tsv} and prints a summary.

Requires the ``eq3`` conda env (rdkit).
"""
from __future__ import annotations

import glob
import json
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import lsqr

MSDB_ROOT = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase"))
ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
OUT_DIR = ANALYSIS_DIR / "results" / "thermo_agreement"
BIOCHEM = MSDB_ROOT / "Biochemistry"
RNG = np.random.default_rng(20260804)

SOURCES = {"gc": "Group contribution", "eq": "eQuilibrator", "dgp": "dGPredictor"}


def load_json_dir(pattern: str) -> dict[str, dict]:
    out = {}
    for path in sorted(glob.glob(str(BIOCHEM / pattern))):
        for entry in json.load(open(path)):
            out[entry["id"]] = entry
    return out


def compound_signatures(cpds: dict, radius: int = 1) -> dict[str, Counter]:
    """cpd -> Counter of atom-centred signature ids (identity of the local
    environment), the descriptor family dGPredictor is built on.

    Compounds without a parseable SMILES get no entry, and any reaction that
    touches one is dropped from the signature experiments (it is exactly the
    case dGPredictor cannot decompose either).
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=2 ** 20)
    sigs: dict[str, Counter] = {}
    for cid, entry in cpds.items():
        smi = entry.get("smiles")
        if not smi:
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        fp = gen.GetSparseCountFingerprint(mol)
        sigs[cid] = Counter(dict(fp.GetNonzeroElements()))
    return sigs


def main() -> None:
    feat = pd.read_csv(OUT_DIR / "reaction_features.tsv", sep="\t", low_memory=False)
    three = feat[feat["n_sources"] == 3].copy()
    keep_ids = set(three["rxn"])
    print(f"{len(three)} reactions with all three sources")

    print("loading biochemistry ...")
    cpds = load_json_dir("compound_*.json")
    rxns = load_json_dir("reaction_*.json")

    # stoichiometry per reaction, restricted to the 3-source set
    stoich: dict[str, dict[str, float]] = {}
    for rid in keep_ids:
        entry = rxns.get(rid)
        if not entry:
            continue
        vec: dict[str, float] = defaultdict(float)
        for item in entry.get("stoichiometry") or []:
            coeff = float(item.get("coefficient", 0) or 0)
            if coeff:
                vec[item["compound"]] += coeff
        vec = {k: v for k, v in vec.items() if v != 0}
        if vec:
            stoich[rid] = vec
    print(f"  {len(stoich)} reactions with usable stoichiometry")

    dg = {src: dict(zip(three["rxn"], three[f"dg_{src}"])) for src in SOURCES}

    # ------------------------------------------------------------------ E1
    print("\nE1. additivity: fitting implied per-compound energies ...")
    rid_list = sorted(stoich)
    cpd_counts = Counter(c for rid in rid_list for c in stoich[rid])
    # keep compounds seen in >=3 reactions so the fit is not memorising singletons
    cpd_list = sorted(c for c, n in cpd_counts.items() if n >= 3)
    cpd_idx = {c: i for i, c in enumerate(cpd_list)}
    usable = [rid for rid in rid_list if all(c in cpd_idx for c in stoich[rid])]
    print(f"  {len(usable)} reactions over {len(cpd_list)} compounds "
          f"(compounds appearing in >=3 reactions)")

    rows, cols, vals = [], [], []
    for i, rid in enumerate(usable):
        for c, v in stoich[rid].items():
            rows.append(i); cols.append(cpd_idx[c]); vals.append(v)
    S = sparse.csr_matrix((vals, (rows, cols)), shape=(len(usable), len(cpd_list)))

    perm = RNG.permutation(len(usable))
    n_test = len(usable) // 5
    test_i, train_i = perm[:n_test], perm[n_test:]

    fit_rows = []
    for src in SOURCES:
        y = np.array([dg[src][rid] for rid in usable], float)
        sol = lsqr(S[train_i], y[train_i], damp=1e-3, atol=1e-10, btol=1e-10,
                   iter_lim=6000)[0]
        for split, idx in (("train", train_i), ("test", test_i)):
            pred = S[idx] @ sol
            resid = y[idx] - pred
            ss_res = float(resid @ resid)
            ss_tot = float(((y[idx] - y[idx].mean()) ** 2).sum())
            fit_rows.append({
                "source": SOURCES[src], "split": split, "n": len(idx),
                "r2": 1 - ss_res / ss_tot,
                "median_abs_resid": float(np.median(np.abs(resid))),
                "frac_within_1kcal": float((np.abs(resid) < 1).mean()),
                "frac_within_5kcal": float((np.abs(resid) < 5).mean()),
            })
    fit_df = pd.DataFrame(fit_rows)
    fit_df.to_csv(OUT_DIR / "additivity_fit.tsv", sep="\t", index=False,
                  float_format="%.4f")
    print(fit_df.to_string(index=False))

    # ------------------------------------------------------------------ E2
    print("\nE2. signature null space ...")
    sigs = compound_signatures(cpds, radius=1)
    print(f"  signatures for {len(sigs)} compounds")

    dsig_key: dict[str, tuple] = {}
    n_drop = 0
    for rid in rid_list:
        vec = stoich[rid]
        if any(c not in sigs for c in vec):
            n_drop += 1
            continue
        acc: Counter = Counter()
        for c, coeff in vec.items():
            for sid, cnt in sigs[c].items():
                acc[sid] += coeff * cnt
        acc = {k: v for k, v in acc.items() if v != 0}
        dsig_key[rid] = tuple(sorted(acc.items()))
    print(f"  {len(dsig_key)} reactions with a computable signature difference "
          f"({n_drop} dropped for missing structures)")

    classes: dict[tuple, list[str]] = defaultdict(list)
    for rid, key in dsig_key.items():
        classes[key].append(rid)
    multi = {k: v for k, v in classes.items() if len(v) > 1}
    n_in_multi = sum(len(v) for v in multi.values())
    print(f"  {len(classes)} distinct signature-difference classes; "
          f"{len(multi)} classes contain >1 reaction, covering {n_in_multi} reactions")

    # Empty dSig = reaction invisible to a signature-difference model entirely
    empty = [rid for rid, k in dsig_key.items() if len(k) == 0]
    print(f"  {len(empty)} reactions have an EMPTY signature difference "
          f"(structurally indistinguishable from a no-op)")

    # variance decomposition on reactions that sit in a multi-member class
    cls_rows = []
    for key, members in multi.items():
        rec = {"n": len(members), "rxns": ";".join(sorted(members)[:8])}
        for src in SOURCES:
            v = np.array([dg[src][r] for r in members], float)
            rec[f"sd_{src}"] = float(v.std())
            rec[f"range_{src}"] = float(v.max() - v.min())
        rec["empty_dsig"] = int(len(key) == 0)
        cls_rows.append(rec)
    cls_df = pd.DataFrame(cls_rows).sort_values("range_gc", ascending=False)
    cls_df.to_csv(OUT_DIR / "signature_classes.tsv", sep="\t", index=False,
                  float_format="%.3f")

    def var_decomp(rids: list[str], src: str) -> tuple[float, float]:
        """(within-class share of variance, total variance) for one source."""
        by_cls: dict[tuple, list[float]] = defaultdict(list)
        for r in rids:
            by_cls[dsig_key[r]].append(dg[src][r])
        allv = np.array([v for vs in by_cls.values() for v in vs], float)
        if len(allv) < 8:
            return (np.nan, np.nan)
        tot = float(((allv - allv.mean()) ** 2).sum())
        within = float(sum(((np.array(vs) - np.mean(vs)) ** 2).sum()
                           for vs in by_cls.values()))
        return (within / tot if tot > 0 else np.nan, tot / len(allv))

    all_rids = sorted(dsig_key)
    for src in SOURCES:
        share, _ = var_decomp(all_rids, src)
        print(f"  within-signature-class share of {SOURCES[src]} variance: {share:.3%}")
    print("    (dGPredictor's own share is the empirical check: a pure "
          "signature-difference model would score 0%)")

    # ------------------------------------------------------------------ E3
    print("\nE3. per-family breakdown ...")
    three_idx = three.set_index("rxn")
    fam_defs = {
        "EC 1 oxidoreductase": lambda d: d["ec_class"].fillna("none").str.contains("1"),
        "redox pair (NAD/NADP/FAD)": lambda d: d[["cof_nad", "cof_nadh", "cof_nadp",
                                                  "cof_nadph", "cof_fad"]].max(axis=1) == 1,
        "aldehyde group changes": lambda d: d["d_aldehyde"] != 0,
        "phosphoryl transfer (ATP/ADP/AMP)": lambda d: d[["cof_atp", "cof_adp",
                                                          "cof_amp"]].max(axis=1) == 1,
        "phosphoanhydride changes": lambda d: d["d_phosphoanhydride"] != 0,
        "methyl transfer (SAM/SAH)": lambda d: d[["cof_sam", "cof_sah"]].max(axis=1) == 1,
        "acetal / glycosidic present": lambda d: d["g_acetal"] > 0,
        "thioester changes (CoA acyl)": lambda d: d["d_thioester"] != 0,
        "ester changes": lambda d: d["d_ester"] != 0,
        "transport": lambda d: d["is_transport"] == 1,
        "EC 2 transferase": lambda d: d["ec_class"].fillna("none").str.contains("2"),
        "EC 3 hydrolase": lambda d: d["ec_class"].fillna("none").str.contains("3"),
        "EC 6 ligase": lambda d: d["ec_class"].fillna("none").str.contains("6"),
    }
    rows = []
    for name, fn in fam_defs.items():
        mask = fn(three_idx)
        rids = [r for r in three_idx.index[mask.to_numpy()] if r in dsig_key]
        if len(rids) < 30:
            continue
        rec = {"family": name, "n": len(rids)}
        for src in SOURCES:
            share, var = var_decomp(rids, src)
            rec[f"within_class_share_{src}"] = share
        # degeneracy: how often reactions in this family collide in dSig
        keys = [dsig_key[r] for r in rids]
        rec["n_distinct_dsig"] = len(set(keys))
        rec["degeneracy"] = len(rids) / max(len(set(keys)), 1)
        rec["frac_empty_dsig"] = float(np.mean([len(k) == 0 for k in keys]))
        # how compressed is dGPredictor's output relative to GC on this family
        gcv = np.array([dg["gc"][r] for r in rids], float)
        dpv = np.array([dg["dgp"][r] for r in rids], float)
        rec["sd_gc"] = float(gcv.std())
        rec["sd_dgp"] = float(dpv.std())
        rec["sd_ratio_dgp_gc"] = float(dpv.std() / gcv.std()) if gcv.std() else np.nan
        rec["n_distinct_dgp_values"] = int(len(np.unique(np.round(dpv, 2))))
        rec["distinct_dgp_frac"] = rec["n_distinct_dgp_values"] / len(rids)
        rows.append(rec)
    fam_df = pd.DataFrame(rows).sort_values("within_class_share_gc")
    fam_df.to_csv(OUT_DIR / "nullspace_by_family.tsv", sep="\t", index=False,
                  float_format="%.4f")
    pd.set_option("display.width", 220)
    print(fam_df.to_string(index=False))


if __name__ == "__main__":
    main()
