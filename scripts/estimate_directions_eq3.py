#!/usr/bin/env python3
"""Estimate MSDB reaction directions from eQuilibrator 3.0 (Beber et al. 2022) energies.

Run with the `eq3` env (equilibrator-api 0.7.0 == component-contribution / eQ 3.0):
    XDG_CACHE_HOME=/scratch/ctaylor/eq_cache \
      /mnt/homes/ctaylor/conda/miniforge3/envs/eq3/bin/python scripts/estimate_directions_eq3.py

For each MSDB reaction we map every metabolite to an eQuilibrator compound
(KEGG -> ChEBI -> BiGG -> MetaNetX -> InChIKey, first hit wins), build the
reaction, and compute the standard transformed Gibbs energy with eQ 3.0
(component contribution; eQ defaults pH 7, I=0.25 M, no Mg, 298.15 K).

Direction uses eQuilibrator's OWN method -- the reversibility index
ln(gamma) = (2/N)*(dG'm/RT) (Noor 2012, reported by eQuilibrator), via
cc.ln_reversibility_index(): reversible if |log10 gamma| < 3
(|ln gamma| < ln(1000) = 6.908), else forward (ln gamma < 0) / reverse
(ln gamma > 0).  This is deliberately NOT the Henry-2007 feasibility rule used
for the group-contribution (Jankowski) column.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import re
import sys
import time

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(SCRIPTS)
MSDB = os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase")
sys.path.insert(0, SCRIPTS)

from equilibrator_api import ComponentContribution, Reaction

KCAL = 4.184  # kJ per kcal
LNRI_THRESH = 6.9077553  # ln(1000) == |log10 gamma| of 3 (Noor 2012 reversibility cutoff)
_KEGG = re.compile(r"KEGG:\s*(C\d{5})")
_CHEBI = re.compile(r"ChEBI:\s*(\d+)")
_BIGG = re.compile(r"BiGG:\s*([^\s;]+)")
_MNX = re.compile(r"MetaNetX:\s*(MNXM\d+)")


def compound_accessions(c):
    """Ordered eQuilibrator accessions for a ModelSEED compound (best first)."""
    al = c.get("aliases") or []
    text = al if isinstance(al, str) else "\n".join(al)
    accs = []
    m = _KEGG.search(text)
    if m:
        accs.append("kegg:" + m.group(1))
    m = _CHEBI.search(text)
    if m:
        accs.append("chebi:CHEBI:" + m.group(1))
    m = _BIGG.search(text)
    if m:
        accs.append("bigg.metabolite:" + m.group(1))
    m = _MNX.search(text)
    if m:
        accs.append("metanetx.chemical:" + m.group(1))
    ik = c.get("inchikey")
    if ik:
        accs.append(ik)  # bare InChIKey resolves via cc.get_compound
    return accs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="0 = all reactions")
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "rxn_directions_eq3_2022.tsv"))
    args = ap.parse_args()

    print("[eq3] loading ComponentContribution (eQ 3.0 cache) ...", flush=True)
    cc = ComponentContribution()

    cpd_accs = {}
    for f in glob.glob(os.path.join(MSDB, "Biochemistry", "compound_*.json")):
        for c in json.load(open(f)):
            cpd_accs[c["id"]] = compound_accessions(c)

    rxns = []
    for f in sorted(glob.glob(os.path.join(MSDB, "Biochemistry", "reaction_*.json"))):
        rxns.extend(json.load(open(f)))
    if args.limit:
        rxns = rxns[: args.limit]
    print(f"[eq3] compounds: {len(cpd_accs)} | reactions to attempt: {len(rxns)}", flush=True)

    resolved = {}  # cpd -> Compound or None (cached)

    def resolve(cpd):
        if cpd in resolved:
            return resolved[cpd]
        obj = None
        for acc in cpd_accs.get(cpd, []):
            try:
                obj = cc.get_compound(acc)
            except Exception:
                obj = None
            if obj is not None:
                break
        resolved[cpd] = obj
        return obj

    out = {}
    n_ok = n_unmapped = n_unbal = n_err = 0
    t0 = time.time()
    for i, rxn in enumerate(rxns):
        rid = rxn["id"]
        stoich = rxn.get("stoichiometry")
        if not stoich:
            continue
        sd, ok = {}, True
        for s in stoich:
            obj = resolve(s["compound"])
            if obj is None:
                ok = False
                break
            sd[obj] = sd.get(obj, 0.0) + float(s["coefficient"])
        if not ok:
            n_unmapped += 1
            continue
        sd = {k: v for k, v in sd.items() if v != 0}
        if not sd:
            continue
        try:
            r = Reaction(sd)
            if not r.is_balanced():
                n_unbal += 1
                continue
            # eQuilibrator's OWN directionality: the reversibility index
            # ln(gamma) = (2/N)*(dG'm/RT) (Noor 2012), as eQ reports it.
            lnri = float(cc.ln_reversibility_index(r).value.m_as(""))
        except Exception:
            n_err += 1
            continue
        d = ">" if lnri < -LNRI_THRESH else "<" if lnri > LNRI_THRESH else "="
        out[rid] = (d, lnri)
        n_ok += 1
        if (i + 1) % 2000 == 0:
            print(f"  {i+1}/{len(rxns)} | ok {n_ok} unmapped {n_unmapped} "
                  f"unbal {n_unbal} err {n_err} | {time.time()-t0:.0f}s", flush=True)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["rxn_id", "Beber_2022", "ln_RI"])
        for rid in sorted(out):
            d, lnri = out[rid]
            w.writerow([rid, d, f"{lnri:.3f}"])
    print(f"[eq3] ok {n_ok} | unmapped {n_unmapped} | unbalanced {n_unbal} "
          f"| err {n_err} | {time.time()-t0:.0f}s")
    print(f"[eq3] wrote {args.out} ({len(out)} reactions)")


if __name__ == "__main__":
    main()
