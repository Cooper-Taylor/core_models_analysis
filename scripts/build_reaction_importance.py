#!/usr/bin/env python3
"""Reaction direction-importance across ALL core models + model clustering.

Part 1 (compute, parallel, resumable)
  For every core model (~5,683) and every reaction in it, start from the model's
  DEFAULT bounds and set that one reaction to each of the four direction options,
  re-solving biomass growth each time (all other reactions left at default):

      "<"(-1000,0)   ">"(0,1000)   "="(-1000,1000)   "?"(0,0)=off/knockout

  Per (model, reaction) we derive flux-based importance features:
      spread = max(growth over the 4 options) - min(...)      (direction sensitivity)
      boost  = max(growth) - baseline                          (best achievable gain)
      ess    = 1 if baseline>0 and growth(off)==0             (knockout-essential)
  Written one JSON line per model to results/reaction_importance_raw.jsonl.

Part 2 (analyze)
  - Aggregate per base reaction across all models -> importance table
    (n_present, n_sensitive, mean/max spread, n_essential, mean boost, ...).
  - Cluster the models by WHICH reactions most influence their growth: build a
    models x top-K-reaction matrix of relative influence (spread/baseline), scale,
    KMeans (k chosen by silhouette), PCA to 2-D for a scatter, and extract each
    cluster's signature reactions.
  - Emit the small site JSON site/data/reaction_importance.json for the website.

Pure cobra FBA (no modelseedpy needed for the solves). Clustering uses scikit-learn.

    PY=/mnt/homes/ctaylor/conda/miniforge3/envs/core_models_analysis/bin/python
    $PY scripts/build_reaction_importance.py --workers 96          # compute (resumable) + analyze
    $PY scripts/build_reaction_importance.py --analyze-only        # re-run analysis on existing raw
"""

from __future__ import annotations

import argparse
import datetime
import json
import math
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from direction_change_template_eval import build_base_to_model_index, _find_biomass

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR", "/scratch/ctaylor/core_models_analysis"))
RESULTS_DIR = ANALYSIS_DIR / "results"
SITE_DATA = ANALYSIS_DIR / "site" / "data"
MODELS_DIR = Path("/scratch/ctaylor/core_models_kegg2")
RAW_PATH = RESULTS_DIR / "reaction_importance_raw.jsonl"
SITE_PATH = SITE_DATA / "reaction_importance.json"

OPTION_BOUNDS = {"<": (-1000.0, 0.0), ">": (0.0, 1000.0), "=": (-1000.0, 1000.0), "?": (0.0, 0.0)}
TOL = 1e-6


# ---------------------------------------------------------------------------
# Part 1 — compute
# ---------------------------------------------------------------------------
def _init(_):
    import cobra
    cobra.Configuration().solver = "glpk"


def eval_model(model_id):
    try:
        import cobra
        model = cobra.io.load_json_model(str(MODELS_DIR / f"{model_id}.json"))
        bio = _find_biomass(model)
        if bio is None:
            return {"model_id": model_id, "error": "no biomass"}
        model.objective = bio
        model.objective_direction = "max"
        v = model.slim_optimize()
        baseline = round(float(v), 6) if (v is not None and v == v) else 0.0

        base2mdl = build_base_to_model_index(model)
        rows = []
        for base, mids in base2mdl.items():
            for rid in mids:
                r = model.reactions.get_by_id(rid)
                gv = []
                goff = None
                for opt, (lb, ub) in OPTION_BOUNDS.items():
                    try:
                        with model:
                            r.bounds = (lb, ub)
                            x = model.slim_optimize()
                        x = float(x) if (x is not None and x == x) else None
                    except Exception:
                        x = None
                    if opt == "?":
                        goff = x
                    if x is not None:
                        gv.append(x)
                if not gv:
                    continue
                spread = max(gv) - min(gv)
                boost = max(gv) - baseline
                ess = 1 if (baseline > TOL and goff is not None and goff < TOL) else 0
                rows.append([base, round(spread, 6), round(boost, 6), ess])
        return {"model_id": model_id, "baseline": baseline, "rxn": rows}
    except Exception as exc:
        return {"model_id": model_id, "error": f"{type(exc).__name__}: {exc}"}


def run_compute(workers, limit):
    ids = sorted(p.stem for p in MODELS_DIR.glob("*.json"))
    if limit:
        ids = ids[:limit]
    done = set()
    if RAW_PATH.exists():
        for line in RAW_PATH.read_text().splitlines():
            try:
                done.add(json.loads(line)["model_id"])
            except Exception:
                pass
    todo = [m for m in ids if m not in done]
    print(f"compute: total={len(ids)} done={len(done)} todo={len(todo)} workers={workers}", flush=True)
    if not todo:
        return
    import multiprocessing as mp
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    n = 0
    with RAW_PATH.open("a") as out, \
            mp.Pool(workers, initializer=_init, initargs=(None,), maxtasksperchild=40) as pool:
        for rec in pool.imap_unordered(eval_model, todo, chunksize=1):
            out.write(json.dumps(rec) + "\n")
            out.flush()
            n += 1
            if n % 200 == 0 or n == len(todo):
                rate = n / max(time.time() - t0, 1e-9)
                print(f"  {n}/{len(todo)}  {rate:.1f} models/s  ETA {(len(todo)-n)/max(rate,1e-9)/60:.1f} min",
                      flush=True)


# ---------------------------------------------------------------------------
# Part 2 — analyze + cluster
# ---------------------------------------------------------------------------
def _load_names():
    names = {}
    for f in ("reactions_panel.json", "reactions_other.json"):
        p = SITE_DATA / f
        if not p.exists():
            continue
        try:
            d = json.loads(p.read_text())
            for k, v in d.items():
                if isinstance(v, dict) and v.get("name"):
                    names.setdefault(k, v["name"])
        except Exception:
            pass
    return names


def run_analyze(topk, kmin, kmax):
    import numpy as np
    from sklearn.preprocessing import StandardScaler
    from sklearn.cluster import KMeans
    from sklearn.decomposition import PCA
    from sklearn.metrics import silhouette_score

    if not RAW_PATH.exists():
        raise SystemExit(f"no raw data at {RAW_PATH}; run compute first")

    # aggregate
    agg = {}  # base -> dict of counters
    model_feat = {}   # mid -> {base: relspread}
    model_meta = {}   # mid -> {baseline, top_base, top_spread}
    n_err = 0
    for line in RAW_PATH.read_text().splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        if "error" in rec:
            n_err += 1
            continue
        mid = rec["model_id"]
        base_flux = rec.get("baseline", 0.0) or 0.0
        rxns = rec.get("rxn", [])
        # per-model normalization: influence = spread / this model's max spread, so the
        # clustering captures the PATTERN of influential reactions, scale-free and robust
        # to non-growers (baseline ~ 0).
        model_max = max((sp for _, sp, _, _ in rxns), default=0.0)
        feat = {}
        top_base, top_spread = None, 0.0
        for base, spread, boost, ess in rxns:
            a = agg.setdefault(base, {"n_present": 0, "n_sensitive": 0, "sum_spread": 0.0,
                                      "max_spread": 0.0, "n_ess": 0, "sum_boost": 0.0, "n_boost": 0})
            a["n_present"] += 1
            if spread > TOL:
                a["n_sensitive"] += 1
                a["sum_spread"] += spread
                a["max_spread"] = max(a["max_spread"], spread)
                feat[base] = spread / model_max if model_max > TOL else 0.0   # 0..1 within-model
                if spread > top_spread:
                    top_spread, top_base = spread, base
            if ess:
                a["n_ess"] += 1
            if boost > TOL:
                a["n_boost"] += 1
                a["sum_boost"] += boost
        model_feat[mid] = feat
        model_meta[mid] = {"baseline": round(base_flux, 4), "top_base": top_base,
                           "top_spread": round(top_spread, 4)}

    names = _load_names()
    n_models = len(model_feat)

    # importance table (rank by n_sensitive, tiebreak mean spread)
    imp = []
    for base, a in agg.items():
        ns = a["n_sensitive"]
        mean_spread = a["sum_spread"] / ns if ns else 0.0
        imp.append({
            "base": base, "name": names.get(base, ""),
            "n_present": a["n_present"], "n_sensitive": ns,
            "frac_sensitive": round(ns / a["n_present"], 4) if a["n_present"] else 0.0,
            "mean_spread": round(mean_spread, 4), "max_spread": round(a["max_spread"], 4),
            "n_essential": a["n_ess"],
            "n_boost": a["n_boost"],
            "mean_boost": round(a["sum_boost"] / a["n_boost"], 4) if a["n_boost"] else 0.0,
        })
    imp.sort(key=lambda x: (x["n_sensitive"], x["mean_spread"]), reverse=True)

    # ---- clustering: models x top-K influential reactions ----
    feat_bases = [r["base"] for r in imp[:topk]]
    bidx = {b: j for j, b in enumerate(feat_bases)}
    mids = sorted(model_feat.keys())
    X = np.zeros((len(mids), len(feat_bases)), dtype=np.float32)
    for i, mid in enumerate(mids):
        for b, val in model_feat[mid].items():
            j = bidx.get(b)
            if j is not None:
                X[i, j] = val
    Xs = StandardScaler().fit_transform(X)

    # choose k by silhouette (subsample for speed)
    rng_idx = np.arange(len(mids))
    sample = rng_idx if len(mids) <= 3000 else rng_idx[:: max(1, len(mids) // 3000)]
    best_k, best_s, best_labels, best_km = None, -1, None, None
    for k in range(kmin, kmax + 1):
        km = KMeans(n_clusters=k, n_init=5, random_state=0)
        labels = km.fit_predict(Xs)
        try:
            s = silhouette_score(Xs[sample], labels[sample])
        except Exception:
            s = -1
        print(f"  k={k} silhouette={s:.3f}", flush=True)
        if s > best_s:
            best_k, best_s, best_labels, best_km = k, s, labels, km
    labels = best_labels
    print(f"chosen k={best_k} (silhouette {best_s:.3f})", flush=True)

    # 2-D PCA for the scatter
    pca = PCA(n_components=2, random_state=0)
    XY = pca.fit_transform(Xs)

    # per-cluster signatures: reactions with the highest mean relative influence in-cluster
    clusters_info = []
    for c in range(best_k):
        members = [i for i in range(len(mids)) if labels[i] == c]
        if not members:
            clusters_info.append({"id": c, "size": 0, "signature": [], "mean_baseline": 0})
            continue
        colmean = X[members].mean(axis=0)
        order = np.argsort(colmean)[::-1][:8]
        sig = [{"base": feat_bases[j], "name": names.get(feat_bases[j], ""),
                "mean_influence": round(float(colmean[j]), 4)} for j in order if colmean[j] > 0]
        mb = float(np.mean([model_meta[mids[i]]["baseline"] for i in members]))
        clusters_info.append({"id": c, "size": len(members), "mean_baseline": round(mb, 3),
                              "signature": sig})

    models_out = []
    for i, mid in enumerate(mids):
        models_out.append({
            "id": mid, "x": round(float(XY[i, 0]), 3), "y": round(float(XY[i, 1]), 3),
            "c": int(labels[i]),
            "top": model_meta[mid]["top_base"], "baseline": model_meta[mid]["baseline"],
        })

    out = {
        "generated": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "n_models": n_models, "n_errors": n_err, "tol": TOL,
        "option_note": "growth tested per reaction under <(-1000,0) >(0,1000) =(-1000,1000) ?=off(0,0), one change at a time vs default",
        "importance_top": imp[:400],
        "clustering": {
            "k": best_k, "silhouette": round(best_s, 3),
            "feature_reactions": feat_bases,
            "pca_explained": [round(float(x), 3) for x in pca.explained_variance_ratio_],
            "clusters": clusters_info,
            "models": models_out,
        },
    }
    SITE_DATA.mkdir(parents=True, exist_ok=True)
    SITE_PATH.write_text(json.dumps(out))
    print(f"DONE analyze: {n_models} models, {len(imp)} reactions, k={best_k} -> {SITE_PATH} "
          f"({SITE_PATH.stat().st_size/1e6:.1f} MB)", flush=True)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--workers", type=int, default=96)
    p.add_argument("--limit", type=int, default=0)
    p.add_argument("--compute-only", action="store_true")
    p.add_argument("--analyze-only", action="store_true")
    p.add_argument("--topk", type=int, default=150, help="# top reactions used as clustering features")
    p.add_argument("--kmin", type=int, default=4)
    p.add_argument("--kmax", type=int, default=10)
    args = p.parse_args(argv)

    if not args.analyze_only:
        run_compute(args.workers, args.limit)
    if not args.compute_only:
        run_analyze(args.topk, args.kmin, args.kmax)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
