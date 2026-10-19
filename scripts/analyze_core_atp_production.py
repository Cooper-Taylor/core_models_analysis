#!/usr/bin/env python3
"""Core-model ATP production (``bio2``) beside core-model biomass (``bio1``).

Every core model in ``data/core_models_kegg2/`` carries TWO reactions whose
name is "biomass", and they are not the same kind of thing:

  bio1   the actual core biomass: 12 central-carbon precursors drawn against
         41.26 ATP and the NAD(P)H pools, producing cpd11416_c0 (Biomass).
  bio2   ATP + H2O -> ADP + Pi + H+. Pure ATP hydrolysis. Despite the name it
         is a *demand* reaction, not a biomass -- it is the test reaction for
         the ModelSEED v2 criterion "Production of ATP from ADP and phosphate
         in at least one media is verified for the reconstructed core model"
         (Faria et al., bioRxiv 2023.10.04.556561). Its flux is bounded only
         by the 1000 upper bound, so it saturates whenever it is feasible at
         all; only the binary feasible/infeasible is meaningful.

This matters because ``results/results.csv`` -- the baseline every other
analysis reads -- optimises **bio1**, and on a single complete medium that
puts 0 of the 299 archaeal core models at non-zero biomass. Reporting that
number without ``bio2`` beside it invites the reading that the archaeal core
models are broken.

They are not, because **ModelSEED never asks a core model to grow.** The
core-stage criterion is ATP production -- "Production of ATP from ADP and
phosphate in at least one media is verified for the reconstructed core model"
-- tested across 54 minimal media formulations (Supplementary Table 1), and
core-level gap-filling is performed only for ATP: "If gap-filling is necessary
for ATP production, we gap-fill the core model in all ATP biosynthesis media
formulations." Biomass gap-filling happens only after the model is expanded to
the genome scale, against auxotrophy or glucose minimal media. So a core model
with zero biomass flux has not failed anything ModelSEED asked of it, and the
grower fraction derived from bio1 is a convention of this repo rather than a
pipeline pass rate.

This script therefore re-solves every core model twice on the same complete
medium the baseline uses, once per objective, and reports the pass rates split
by domain. It writes ``results/core_atp_production.tsv`` (per model) and
prints the contingency table.

It is deliberately NOT part of build_baseline_biomass_taxonomy.py: that script
only reads cached results and solves nothing, and this one solves ~11k LPs.

Usage::

    python3 scripts/analyze_core_atp_production.py --jobs 32
"""
from __future__ import annotations

import argparse
import csv
import logging
import multiprocessing as mp
import os
from collections import Counter
from pathlib import Path

logging.getLogger("cobra").setLevel(logging.ERROR)

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
MODELS_DIR = Path(os.environ.get("CORE_MODELS_DIR",
                                 "/scratch/ctaylor/core_models_kegg2"))
MEDIA_CPD = Path(os.environ.get(
    "COMPLETE_MEDIA", "/scratch/ctaylor/ModelSEEDDatabase/Media/KBaseMedia.cpd"))
OUT_TSV = ANALYSIS_DIR / "results" / "core_atp_production.tsv"

GROWTH_THRESHOLD = 1e-6
_MEDIA: set[str] = set()


def _init() -> None:
    global _MEDIA
    _MEDIA = {line.strip() for line in open(MEDIA_CPD) if line.strip()}


def _solve(model_id: str) -> dict:
    from cobra.io import load_json_model
    row = {"model_id": model_id, "bio1": 0.0, "bio2": 0.0,
           "status1": "", "status2": "", "error": ""}
    try:
        m = load_json_model(str(MODELS_DIR / f"{model_id}.json"))
        # Same media gate as analyze_growth.py: an exchange is opened only if
        # its compound is in the complete-media list; everything else is shut.
        for ex in m.exchanges:
            cpd = list(ex.metabolites)[0].id.rsplit("_", 1)[0]
            ex.lower_bound = -1000.0 if cpd in _MEDIA else 0.0
        ids = {r.id for r in m.reactions}
        for key, obj in (("1", "bio1"), ("2", "bio2")):
            if obj not in ids:
                row[f"status{key}"] = "absent"
                continue
            m.objective = m.reactions.get_by_id(obj)
            sol = m.optimize()
            row[f"status{key}"] = sol.status
            row[obj] = float(sol.objective_value or 0.0)
    except Exception as exc:  # noqa: BLE001
        row["error"] = f"{type(exc).__name__}: {exc}"
    return row


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=32)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    ids = sorted(p.stem for p in MODELS_DIR.glob("GCF_*.json"))
    if args.limit:
        ids = ids[: args.limit]
    print(f"solving {len(ids)} core models x 2 objectives on "
          f"{MEDIA_CPD.name} ...")

    with mp.get_context("fork").Pool(args.jobs, initializer=_init) as pool:
        rows = list(pool.imap_unordered(_solve, ids, chunksize=16))
    rows.sort(key=lambda r: r["model_id"])

    # join taxonomy for the domain split
    tax: dict[str, dict] = {}
    for name in ("ncbi_taxonomy.csv", "ncbi_taxonomy_fill.csv"):
        path = ANALYSIS_DIR / "results" / name
        if not path.exists():
            continue
        for r in csv.DictReader(open(path)):
            acc = r["assembly_accession"]
            if acc not in tax or tax[acc].get("status") != "resolved":
                tax[acc] = r

    for r in rows:
        lin = tax.get(r["model_id"], {})
        r["superkingdom"] = lin.get("superkingdom", "") or "Unassigned"
        r["phylum"] = lin.get("phylum", "") or "Unassigned"
        r["organism_name"] = lin.get("organism_name", "")
        r["biomass_ok"] = r["bio1"] > GROWTH_THRESHOLD
        r["atp_ok"] = r["bio2"] > GROWTH_THRESHOLD

    cols = ["model_id", "superkingdom", "phylum", "organism_name", "bio1",
            "bio2", "biomass_ok", "atp_ok", "status1", "status2", "error"]
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT_TSV}")

    # ---- report ----------------------------------------------------------
    print()
    print(f"{'domain':14s} {'n':>6s} {'biomass (bio1)':>16s} {'ATP (bio2)':>14s}")
    by_dom: dict[str, list] = {}
    for r in rows:
        by_dom.setdefault(r["superkingdom"], []).append(r)
    for dom, grp in sorted(by_dom.items(), key=lambda kv: -len(kv[1])):
        nb = sum(1 for r in grp if r["biomass_ok"])
        na = sum(1 for r in grp if r["atp_ok"])
        print(f"{dom:14s} {len(grp):6d} {nb:7d} ({100*nb/len(grp):5.1f}%) "
              f"{na:6d} ({100*na/len(grp):5.1f}%)")
    print()
    print("contingency (all models):")
    cont = Counter((r["biomass_ok"], r["atp_ok"]) for r in rows)
    for (b, a), n in sorted(cont.items(), key=lambda kv: -kv[1]):
        print(f"  biomass={'yes' if b else ' no'}  ATP={'yes' if a else ' no'}  {n:6d}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
