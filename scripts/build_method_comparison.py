#!/usr/bin/env python3
"""Compare reaction-direction methods against the KEGG core-model default.

Reference (the direction the 5,683 KEGG core models were actually built with):
  - KEGG_default   : per-reaction majority of the on-disk COBRA flux bounds
                     (lower_bound/upper_bound) across every core model that
                     contains the reaction -- (0,1000)=>, (-1000,0)=<,
                     (-1000,1000)==, (0,0)=blocked (ignored). Models write each
                     reaction in MSDB's canonical orientation (verified: 0 flips),
                     so the bound direction is directly comparable to the heuristics.

Heuristics / predictions (all computed from ModelSEEDDatabase-stored values):
  - Jankowski_2007 : group-contribution ΔG' with the Henry-2007 feasibility
                     window  (MSDB thermodynamics["Group contribution"]).
  - Flamholz_2012  : eQuilibrator reversibility index, |ln γ| > ln(1000) rule,
                     from MSDB's stored eQ energies
                     (Biochemistry/Thermodynamics/eQuilibrator/
                      MetaNetX_Reaction_Energies.tbl  ln_RI column).
  - Opus_4.8       : Claude Opus 4.8 predictions
                     (data/ai_curation/all_modelseed/AICurationCacheReactionDirectionality.json).

Four reaction scopes are emitted as `modes`:
  - models / models_no_transport : the ~237 reactions present in the core models,
        with KEGG_default as the first (reference) method, then the heuristics.
  - all / no_transport           : the original all-MSDB symmetric view over the
        three heuristics (no KEGG_default -- the default only exists in-model).

For each scope it emits, over the union of reactions each pair co-covers:
  - per-method direction distribution (>, <, =, ?),
  - an N×N pairwise agreement matrix (over reactions both methods call directionally),
  - the pairwise confusion matrices (direction x direction counts),
ordered so the reference sits first (model scopes) or by rank similarity (wide scopes).

Output: site/data/method_comparison.json.
"""

from __future__ import annotations

import glob
import itertools
import json
import os
import sys
from collections import Counter
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
ROOT = SCRIPTS.parent
MSDB = os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase")
# Core-model JSONs (read-only). The repo symlinks data/core_models_kegg2 at the
# unpacked 5,683-model directory; KEGG_MODELS_DIR overrides for other layouts.
KEGG_MODELS_DIR = os.environ.get(
    "KEGG_MODELS_DIR",
    str(ROOT / "data" / "core_models_kegg2") if (ROOT / "data" / "core_models_kegg2").exists()
    else "/scratch/ctaylor/core_models_kegg2")
OUT = ROOT / "site" / "data" / "method_comparison.json"
sys.path.insert(0, str(SCRIPTS))

from reversibility_lib import ReversibilityConfig, _walk_stoichiometry, load_ln_reversibility_index
from estimate_directions_literature import source_pair, feasibility_dir, LNRI_THRESH, GC_LABEL, LLM_JSON
from direction_pipeline import _normalize_direction
from seed_annotation import normalize_seed_id

REFERENCE = "KEGG_default"                                  # core-model-built direction
HEURISTICS = ["Jankowski_2007", "Flamholz_2012", "Opus_4.8"]
METHODS_MODEL = [REFERENCE] + HEURISTICS                    # model scopes (reference first)
METHODS = HEURISTICS                                        # wide all-MSDB scopes
DIRS = [">", "<", "=", "?"]
DIRECTIONAL = (">", "<", "=")


def _bound_dir(lb, ub):
    """COBRA flux bounds -> direction call ('>', '<', '=', or 'blocked')."""
    if lb < 0 and ub > 0:
        return "="
    if lb >= 0 and ub > 0:
        return ">"
    if lb < 0 and ub <= 0:
        return "<"
    return "blocked"


def load_kegg_default():
    """Per-reaction direction the KEGG core models were built with.

    Tally the on-disk bound direction over every core model containing each
    reaction (keyed by its normalized ``seed.reaction`` annotation) and take the
    majority directional call, ignoring 'blocked' (0,0) instances. Reactions that
    are blocked in every model contribute no default. Models write reactions in
    MSDB's canonical orientation (0 flips verified), so no re-orientation is needed.
    """
    votes = {}
    files = sorted(glob.glob(os.path.join(KEGG_MODELS_DIR, "GCF_*.json")))
    for f in files:
        for rxn in json.load(open(f)).get("reactions", []):
            anno = rxn.get("annotation") or {}
            seed = anno.get("seed.reaction")
            if isinstance(seed, list):
                seed = seed[0] if seed else None
            if not seed:
                continue
            seed = normalize_seed_id(seed)
            d = _bound_dir(rxn.get("lower_bound", 0), rxn.get("upper_bound", 0))
            votes.setdefault(seed, Counter())[d] += 1
    kegg = {}
    for sid, c in votes.items():
        directional = Counter({k: v for k, v in c.items() if k in DIRECTIONAL})
        if directional:
            kegg[sid] = directional.most_common(1)[0][0]
    print(f"[method-cmp] KEGG_default: {len(files)} models -> {len(votes)} reactions "
          f"({len(kegg)} with a directional default)")
    return kegg


def load_methods():
    """Return ({method: {rxn: dir}}, {rxn: is_transport}).

    KEGG_default comes from the core-model bounds; the heuristics from MSDB.
    """
    cfg = ReversibilityConfig()
    jank, is_transport = {}, {}
    for f in sorted(glob.glob(os.path.join(MSDB, "Biochemistry", "reaction_*.json"))):
        for rxn in json.load(open(f)):
            rid = rxn["id"]
            is_transport[rid] = 1 if rxn.get("is_transport") in (1, "1", True) else 0
            gcp = source_pair(rxn.get("thermodynamics"), GC_LABEL)
            stoich = rxn.get("stoichiometry")
            if gcp is not None and stoich:
                terms = _walk_stoichiometry(stoich, cfg)
                jank[rid] = feasibility_dir(gcp[0], gcp[1], terms, cfg)
    flam = {}
    for rid, lnri in load_ln_reversibility_index().items():  # MSDB MetaNetX_Reaction_Energies.tbl
        flam[rid] = ">" if lnri < -LNRI_THRESH else "<" if lnri > LNRI_THRESH else "="
    opus = {}
    for rid, e in json.load(open(LLM_JSON)).items():
        if isinstance(e, dict) and e.get("directionality"):
            opus[rid] = _normalize_direction(e["directionality"])
    return {
        REFERENCE: load_kegg_default(),
        "Jankowski_2007": jank,
        "Flamholz_2012": flam,
        "Opus_4.8": opus,
    }, is_transport


def seriate(sim):
    """Order indices 0..n-1 to maximize summed adjacent similarity (brute force; n<=4)."""
    n = len(sim)
    best, best_score = list(range(n)), float("-inf")
    for perm in itertools.permutations(range(n)):
        s = sum(sim[perm[k]][perm[k + 1]] for k in range(n - 1))
        if s > best_score:
            best_score, best = s, list(perm)
    return best


def analyze(maps, keep, methods, reference=None):
    """`keep` = set of rxn ids in scope. Returns the comparison payload.

    `methods` is the ordered method list for this scope. When `reference` is set
    (the KEGG_default-anchored model scopes) it is pinned first and the remaining
    methods are ordered by descending agreement with it, so the matrices read as
    default-vs-heuristics; otherwise the method axis is rank-similarity seriated.
    Confusion matrices are emitted in `methods` order, so reference-anchored pairs
    come first with the reference on the rows.
    """
    n = len(methods)
    sub = {m: {r: d for r, d in maps[m].items() if r in keep} for m in methods}
    dist = {m: {d: 0 for d in DIRS} for m in methods}
    for m in methods:
        for d in sub[m].values():
            dist[m][d] = dist[m].get(d, 0) + 1

    # pairwise agreement over co-decided (both directional) reactions
    agr_matrix = [[0.0] * n for _ in range(n)]
    agree_n = [[0] * n for _ in range(n)]
    dec = {m: {r: d for r, d in sub[m].items() if d in DIRECTIONAL} for m in methods}
    for i in range(n):
        for j in range(n):
            ci, cj = dec[methods[i]], dec[methods[j]]
            common = ci.keys() & cj.keys()
            agree_n[i][j] = len(common)
            agr_matrix[i][j] = (sum(1 for r in common if ci[r] == cj[r]) / len(common)) if common else 1.0
    if reference is not None and reference in methods:
        r0 = methods.index(reference)
        rest = sorted((k for k in range(n) if k != r0), key=lambda k: -agr_matrix[r0][k])
        morder = [r0] + rest
    else:
        morder = seriate(agr_matrix)

    # pairwise confusion matrices (direction x direction, includes '?')
    confusions = []
    for a, b in itertools.combinations(range(n), 2):
        ma, mb = methods[a], methods[b]
        common = sub[ma].keys() & sub[mb].keys()
        idx = {d: k for k, d in enumerate(DIRS)}
        cm = [[0] * 4 for _ in range(4)]
        for r in common:
            cm[idx[sub[ma][r]]][idx[sub[mb][r]]] += 1
        # rank-based similarity seriation of the 4 categories (cluster confused ones)
        sim = [[0.0] * 4 for _ in range(4)]
        for x in range(4):
            for y in range(4):
                if x != y:
                    sim[x][y] = cm[x][y] + cm[y][x]
        order = seriate(sim)
        cats = [DIRS[k] for k in order]
        cm_s = [[cm[order[x]][order[y]] for y in range(4)] for x in range(4)]
        diag = sum(cm[k][k] for k in range(4))
        confusions.append({
            "a": ma, "b": mb, "n": len(common),
            "agree": round(diag / len(common), 4) if common else 0.0,
            "cats": cats, "matrix": cm_s,
        })

    return {
        "methods": [methods[k] for k in morder],
        "agreement": [[round(agr_matrix[morder[i]][morder[j]], 4) for j in range(n)] for i in range(n)],
        "agreement_n": [[agree_n[morder[i]][morder[j]] for j in range(n)] for i in range(n)],
        "dist": dist,
        "confusion": confusions,
        "n_reactions": len(keep),
    }


def main():
    maps, is_transport = load_methods()
    all_rids = set(is_transport) | set().union(*[set(maps[m]) for m in HEURISTICS])
    n_tx = sum(is_transport.get(r, 0) for r in all_rids)
    # Core-model scope: the reactions for which the models encode a default.
    model_rids = set(maps[REFERENCE])
    print(f"[method-cmp] wide: {len(all_rids)} reactions ({n_tx} transport); "
          f"Jankowski={len(maps['Jankowski_2007'])} Flamholz={len(maps['Flamholz_2012'])} Opus={len(maps['Opus_4.8'])}")
    print(f"[method-cmp] core-model scope: {len(model_rids)} reactions with a KEGG_default")

    keep_all = all_rids
    keep_notx = {r for r in all_rids if not is_transport.get(r)}
    keep_models = model_rids
    keep_models_notx = {r for r in model_rids if not is_transport.get(r)}
    out = {
        "methods_all": METHODS_MODEL,
        "reference": REFERENCE,
        "counts": {
            "models": len(keep_models), "models_no_transport": len(keep_models_notx),
            "all": len(keep_all), "no_transport": len(keep_notx), "transport": n_tx,
        },
        "modes": {
            "models": analyze(maps, keep_models, METHODS_MODEL, reference=REFERENCE),
            "models_no_transport": analyze(maps, keep_models_notx, METHODS_MODEL, reference=REFERENCE),
            "all": analyze(maps, keep_all, METHODS),
            "no_transport": analyze(maps, keep_notx, METHODS),
        },
    }
    OUT.write_text(json.dumps(out, separators=(",", ":")))
    print(f"[method-cmp] wrote {OUT.name} ({OUT.stat().st_size/1024:.0f} KB)")
    for mode in ("models", "models_no_transport", "all", "no_transport"):
        m = out["modes"][mode]
        print(f"  [{mode}] method order={m['methods']}")
        for c in m["confusion"]:
            print(f"     {c['a']} vs {c['b']}: agree={c['agree']} over n={c['n']} (cats {c['cats']})")


if __name__ == "__main__":
    main()
