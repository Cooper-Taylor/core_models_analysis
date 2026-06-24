#!/usr/bin/env python3
"""Estimate MSDB reaction directions from literature thermodynamic methods and
compare them to the Claude-Opus-4.8 LLM predictions.

Columns produced (one TSV row per MSDB reaction with >=1 prediction):
    rxn_id  Jankowski_2008  Flamholz_2012  LLM_Opus_4.8

Directionality rule -- EACH method uses its OWN approach (do not mix them):

  Jankowski_2008 (group contribution) -- Henry-2007 delta-G' feasibility range
  (Jankowski 2008 relies on Henry 2007 for directionality):
      terms = _walk_stoichiometry(stoich, cfg)          # conc. window 1e-5..0.02 M
      stored_max, stored_min = _stored_bounds(dG, dGerr, terms, cfg)
      ">" if stored_max < 0 ; "<" if stored_min > 0 ; else "="
  (feasibility branch only -- not the full MSDB cascade.)

  Flamholz_2012 (eQuilibrator) -- eQuilibrator's OWN reversibility index
  ln(gamma) = (2/N)*(dG'm/RT) (Noor 2012, reported by eQuilibrator), taken from
  MSDB's precomputed ln_RI column:
      ">" if ln_RI < -ln(1000) ; "<" if ln_RI > ln(1000) ; else "="
  (Deliberately NOT the Henry-2007 rule -- that would apply 2007/2008 heuristics
  to the 2012 method.)

Energy / index sources, read from MSDB:
    Jankowski_2008  <- rxn["thermodynamics"]["Group contribution"]  [dG, dGerr] kcal/mol
    Flamholz_2012   <- Biochemistry/Thermodynamics/eQuilibrator/MetaNetX_Reaction_Energies.tbl
                       (ln_RI column; MSDB's bundled eQuilibrator data)
    LLM_Opus_4.8    <- data/ai_curation/all_modelseed/AICurationCacheReactionDirectionality.json
        directionality forward/reverse/reversible/uncertain -> >/</=/?

The eQ-3.0 (Beber_2022) column is produced separately by estimate_directions_eq3.py
(eQuilibrator 3.0's own ln_reversibility_index) and merged in.

Missing energy / missing reaction -> "NA".
"""
from __future__ import annotations

import csv
import glob
import json
import os
import sys
from collections import Counter

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(SCRIPTS)
MSDB = os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase")
sys.path.insert(0, SCRIPTS)

from reversibility_lib import (
    ReversibilityConfig, _walk_stoichiometry, _stored_bounds, SENTINEL_DG,
    load_ln_reversibility_index,
)
from direction_pipeline import _normalize_direction

LNRI_THRESH = 6.9077553  # ln(1000) == |log10 gamma| of 3 (Noor 2012 reversibility cutoff)

LLM_JSON = os.path.join(
    ROOT, "data", "ai_curation", "all_modelseed",
    "AICurationCacheReactionDirectionality.json")
OUT_TSV = os.path.join(ROOT, "results", "reaction_directions_literature_vs_llm.tsv")
OUT_SUMMARY = os.path.join(ROOT, "results", "reaction_directions_literature_summary.txt")

GC_LABEL = "Group contribution"
EQ_LABEL = "eQuilibrator"


def source_pair(thermo, label):
    """Return (dG, dGerr) for a source, or None if missing/sentinel."""
    if not isinstance(thermo, dict):
        return None
    p = thermo.get(label)
    if not p or p[0] is None:
        return None
    try:
        dg = float(p[0])
    except (TypeError, ValueError):
        return None
    if dg == SENTINEL_DG:
        return None
    try:
        dge = float(p[1])
    except (TypeError, ValueError):
        dge = 0.0
    return dg, dge


def feasibility_dir(dg, dge, terms, cfg):
    """delta-G' feasibility-range direction (Henry 2007)."""
    smax, smin = _stored_bounds(dg, dge, terms, cfg)
    if smax < 0:
        return ">"
    if smin > 0:
        return "<"
    return "="


def main():
    cfg = ReversibilityConfig()  # MSDB baseline: 1e-5..0.02 M, RT=0.5921 kcal/mol

    # --- LLM predictions ---
    llm_raw = json.load(open(LLM_JSON))
    llm = {}
    for rid, e in llm_raw.items():
        if isinstance(e, dict) and e.get("directionality"):
            llm[rid] = _normalize_direction(e["directionality"])
    print(f"[lit] LLM reactions: {len(llm)}")

    # --- Jankowski_2008: group-contribution energies + Henry-2007 feasibility range
    gc = {}
    n_rxn = 0
    for f in sorted(glob.glob(os.path.join(MSDB, "Biochemistry", "reaction_*.json"))):
        for rxn in json.load(open(f)):
            rid = rxn["id"]
            n_rxn += 1
            gcp = source_pair(rxn.get("thermodynamics"), GC_LABEL)
            stoich = rxn.get("stoichiometry")
            if gcp is None or not stoich:
                continue
            terms = _walk_stoichiometry(stoich, cfg)
            gc[rid] = feasibility_dir(gcp[0], gcp[1], terms, cfg)

    # --- Flamholz_2012: eQuilibrator's OWN directionality (reversibility index)
    #     on MSDB's bundled eQ energies (MetaNetX_Reaction_Energies.tbl ln_RI col).
    #     |ln gamma| < ln(1000) -> reversible; ln gamma < 0 -> forward; > 0 -> reverse.
    #     (NOT the Henry-2007 feasibility rule used for the GC column.)
    eq = {}
    for rid, lnri in load_ln_reversibility_index().items():
        eq[rid] = ">" if lnri < -LNRI_THRESH else "<" if lnri > LNRI_THRESH else "="

    print(f"[lit] MSDB reactions scanned: {n_rxn}")
    print(f"[lit] Jankowski_2008 (GC feasibility) directions: {len(gc)}")
    print(f"[lit] Flamholz_2012 (eQ reversibility index) directions: {len(eq)}")

    # --- union of reactions with >=1 prediction ---
    ids = sorted(set(gc) | set(eq) | set(llm))
    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["rxn_id", "Jankowski_2008", "Flamholz_2012", "LLM_Opus_4.8"])
        for rid in ids:
            w.writerow([rid, gc.get(rid, "NA"), eq.get(rid, "NA"), llm.get(rid, "NA")])
    print(f"[lit] wrote {OUT_TSV} ({len(ids)} rows)")

    # --- summary: distributions + pairwise agreement on co-present reactions ---
    def dist(d):
        return dict(Counter(d.values()))

    def agree(a, b):
        common = set(a) & set(b)
        if not common:
            return 0, 0, 0.0
        same = sum(1 for r in common if a[r] == b[r])
        return same, len(common), 100.0 * same / len(common)

    lines = []
    lines.append(f"reactions in file (union): {len(ids)}")
    lines.append(f"Jankowski_2008 (GC) non-NA: {len(gc)}   dist: {dist(gc)}")
    lines.append(f"Flamholz_2012  (eQ) non-NA: {len(eq)}   dist: {dist(eq)}")
    lines.append(f"LLM_Opus_4.8        non-NA: {len(llm)}   dist: {dist(llm)}")
    lines.append("")
    for name, a, b in [("GC vs eQ", gc, eq), ("GC vs LLM", gc, llm),
                       ("eQ vs LLM", eq, llm)]:
        s, n, pct = agree(a, b)
        lines.append(f"{name:10s}: {s}/{n} agree ({pct:.1f}%) on co-present reactions")
    summary = "\n".join(lines)
    open(OUT_SUMMARY, "w").write(summary + "\n")
    print("\n" + summary)
    print(f"\n[lit] wrote {OUT_SUMMARY}")


if __name__ == "__main__":
    main()
