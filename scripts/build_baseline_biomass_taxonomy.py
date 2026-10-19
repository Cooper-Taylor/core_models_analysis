#!/usr/bin/env python3
"""Baseline biomass distribution and genome taxonomy for the ModelSEED v2 model sets.

Two model sets from the ModelSEED v2 manuscript (Faria et al., bioRxiv
2023.10.04.556561) are summarised side by side:

  core      5,683 core models      /scratch/ctaylor/core_models_kegg2/
  full      5,420 genome-scale     /scratch/ctaylor/modelseed2_gs_models/

"Baseline" here means strictly the models AS THEY SIT ON DISK: each model's own
stored lower/upper bounds, its own biomass reaction, the media it was gap-filled
against, and nothing else. No direction map, no thermodynamic source, no
heuristic cascade and no gap-filling is applied. Nothing in this script solves
an LP -- it reads FBA results that were already computed and cached:

  core      results/results.csv
            ModelSEED ``KBaseMedia.cpd`` complete media (347 compounds),
            written by scripts/analyze_growth.py.
  full      results/influence_v201/{gsm_gmm,gsm_auxo}/growth.json, the
            ``on_disk`` entry of each per-model ``runs`` block -- the pipeline's
            reserved label for "no direction overrides at all" (an empty map,
            so no reaction bound is touched before the media is applied).

THE THREE PANELS ARE NOT ON A COMMON MEDIA, so their absolute biomass fluxes
are not comparable and this script never pools them. The core panel runs on a
347-compound complete media with the models' own (largely unconstrained)
exchange bounds, which puts its fluxes two orders of magnitude above the
genome-scale panels; those run on the defined media the manuscript gap-filled
them against, at a -10 mmol/gDW/h default uptake. Compare the SHAPE of each
distribution and the grower fraction across panels, not the flux values.

Taxonomy is joined per genome from ``results/ncbi_taxonomy.csv``, with the 122
suppressed RefSeq accessions and the one genome-scale-only genome recovered
into ``results/ncbi_taxonomy_fill.csv`` by scripts/resolve_missing_taxonomy.py.

Outputs land in ``results/baseline_biomass_taxonomy/``:

  per_model.tsv          one row per (panel, model): growth, grows, lineage
  biomass_summary.tsv    per panel: grower counts and the flux quantiles
  biomass_deciles.tsv    per panel: decile edges of the grower flux
  taxonomy_<rank>.tsv    per rank: genome counts and grower rates per taxon
  biomass_by_phylum.tsv  per panel x phylum: n, growers, median/mean flux
  summary.json           everything above, machine-readable

Usage::

    python3 scripts/build_baseline_biomass_taxonomy.py
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import statistics as stats
import sys
from collections import Counter, defaultdict
from pathlib import Path

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
MS2_DIR = Path(os.environ.get("MS2_GSM_DIR", "/scratch/ctaylor/modelseed2_gs_models"))
RESULTS = ANALYSIS_DIR / "results"
OUT_DIR = RESULTS / "baseline_biomass_taxonomy"

#: FBA objective above which a model counts as growing. Matches the threshold
#: used by analyze_growth.py and by cma.stages (GROWTH_THRESHOLD).
GROWTH_THRESHOLD = 1e-6

RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]

#: panel key -> (model set label, media label, provenance note)
PANELS = {
    "core": ("Core models (5,683)",
             "ModelSEED complete media (KBaseMedia.cpd, 347 compounds)",
             "results/results.csv"),
    "full_gmm": ("Genome-scale models (5,420)",
                 "Glucose minimal media (Carbon-D-Glucose, -10 mmol/gDW/h)",
                 "results/influence_v201/gsm_gmm/growth.json [on_disk]"),
    "full_auxo": ("Genome-scale models (5,420)",
                  "Auxotrophy media (-10 mmol/gDW/h)",
                  "results/influence_v201/gsm_auxo/growth.json [on_disk]"),
}

# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------


def load_taxonomy() -> dict[str, dict]:
    """assembly accession -> lineage dict, base table overlaid with the fill."""
    tax: dict[str, dict] = {}
    with open(RESULTS / "ncbi_taxonomy.csv") as fh:
        for row in csv.DictReader(fh):
            tax[row["assembly_accession"]] = row
    fill = RESULTS / "ncbi_taxonomy_fill.csv"
    if fill.exists():
        with open(fill) as fh:
            for row in csv.DictReader(fh):
                acc = row["assembly_accession"]
                # The fill only ever adds information: it either introduces a
                # genome the base table never had, or replaces a row the base
                # table left blank. A resolved base row is never overwritten.
                if acc not in tax or tax[acc].get("status") != "resolved":
                    tax[acc] = row
    else:
        print(f"WARNING: {fill} missing -- run scripts/resolve_missing_taxonomy.py "
              "first; 122 genomes will be unclassified", file=sys.stderr)
    return tax


def genome_id(model_id: str) -> str:
    """Strip the KBase model suffix to recover the assembly accession.

    The two media carry different suffixes -- ``.RAST.GMM`` for glucose minimal
    and ``.RAST.auxo`` for auxotrophy -- so everything from ``.RAST.`` onward
    goes, rather than matching the media tags by name.

    >>> genome_id("GCF_000005825.2.RAST.GMM")
    'GCF_000005825.2'
    >>> genome_id("GCF_000005825.2.RAST.auxo")
    'GCF_000005825.2'
    """
    return re.sub(r"\.RAST\..*$", "", model_id)


def load_core() -> list[dict]:
    rows = []
    with open(RESULTS / "results.csv") as fh:
        for r in csv.DictReader(fh):
            flux = float(r["growth_flux"] or 0.0)
            rows.append({
                "model_id": r["model_id"],
                "genome_id": r["model_id"],
                "status": r["status"],
                "growth": flux,
                "grows": flux > GROWTH_THRESHOLD and r["status"] == "optimal",
                "n_reactions": int(r["n_reactions"]),
                "n_metabolites": int(r["n_metabolites"]),
                "n_genes": int(r["n_genes"]),
            })
    return rows


def load_full(tag: str) -> list[dict]:
    """Per-model ``on_disk`` growth out of an influence-pipeline growth.json."""
    path = RESULTS / "influence_v201" / tag / "growth.json"
    with open(path) as fh:
        payload = json.load(fh)
    sizes = load_full_sizes()
    rows = []
    for entry in payload["per_model"]:
        run = entry.get("runs", {}).get("on_disk")
        if run is None:
            continue
        gid = genome_id(entry["model_id"])
        size = sizes.get(entry["model_id"], {})
        rows.append({
            "model_id": entry["model_id"],
            "genome_id": gid,
            "status": run["status"],
            "growth": float(run["growth"]),
            "grows": bool(run["grows"]),
            "n_reactions": size.get("n_reactions", ""),
            "n_metabolites": size.get("n_metabolites", ""),
            "n_genes": size.get("n_genes", ""),
        })
    return rows


def load_full_sizes() -> dict[str, dict]:
    """Model sizes from the two cobra conversion manifests, keyed by model id."""
    out: dict[str, dict] = {}
    for sub in ("gmm", "auxotrophy"):
        man = MS2_DIR / sub / "cobra_manifest.tsv"
        if not man.exists():
            continue
        with open(man) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("ok") != "True":
                    continue
                out[r["model_id"]] = {
                    "n_reactions": int(r["n_reactions"]),
                    "n_metabolites": int(r["n_metabolites"]),
                    "n_genes": int(r["n_genes"]),
                }
    return out


# ---------------------------------------------------------------------------
# summarising
# ---------------------------------------------------------------------------


def quantile(sorted_vals: list[float], q: float) -> float:
    """Linear-interpolated quantile; ``sorted_vals`` must be sorted ascending."""
    if not sorted_vals:
        return float("nan")
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = q * (len(sorted_vals) - 1)
    lo = int(pos)
    hi = min(lo + 1, len(sorted_vals) - 1)
    frac = pos - lo
    return sorted_vals[lo] * (1 - frac) + sorted_vals[hi] * frac


def describe(values: list[float]) -> dict:
    vals = sorted(values)
    if not vals:
        return {k: None for k in ("n", "min", "q1", "median", "mean", "q3",
                                  "max", "sd", "cv", "iqr")}
    mean = stats.fmean(vals)
    q1, med, q3 = quantile(vals, .25), quantile(vals, .5), quantile(vals, .75)
    sd = stats.pstdev(vals) if len(vals) > 1 else 0.0
    return {"n": len(vals), "min": vals[0], "q1": q1, "median": med,
            "mean": mean, "q3": q3, "max": vals[-1], "sd": sd,
            "cv": (sd / mean) if mean else None, "iqr": q3 - q1}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=OUT_DIR,
                    help="output directory (default results/baseline_biomass_taxonomy)")
    args = ap.parse_args()
    out_dir: Path = args.out
    out_dir.mkdir(parents=True, exist_ok=True)

    tax = load_taxonomy()
    panels = {"core": load_core(),
              "full_gmm": load_full("gsm_gmm"),
              "full_auxo": load_full("gsm_auxo")}

    # ---- per-model table -------------------------------------------------
    per_model_path = out_dir / "per_model.tsv"
    cols = (["panel", "model_id", "genome_id", "media", "status", "growth",
             "grows", "n_reactions", "n_metabolites", "n_genes", "organism_name",
             "tax_id"] + RANKS)
    unmatched: Counter = Counter()
    with open(per_model_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for panel, rows in panels.items():
            media = PANELS[panel][1]
            for r in rows:
                lineage = tax.get(r["genome_id"])
                if lineage is None:
                    unmatched[panel] += 1
                    lineage = {}
                rec = {"panel": panel, "media": media,
                       "organism_name": lineage.get("organism_name", ""),
                       "tax_id": lineage.get("tax_id", "")}
                rec.update({k: r[k] for k in ("model_id", "genome_id", "status",
                                              "growth", "grows", "n_reactions",
                                              "n_metabolites", "n_genes")})
                for rank in RANKS:
                    rec[rank] = lineage.get(rank, "") or "Unassigned"
                w.writerow(rec)
    if unmatched:
        print(f"WARNING: genomes with no taxonomy row: {dict(unmatched)}",
              file=sys.stderr)

    # ---- biomass distribution per panel ----------------------------------
    summary: dict = {"growth_threshold": GROWTH_THRESHOLD, "panels": {}}
    sum_rows, dec_rows = [], []
    for panel, rows in panels.items():
        label, media, prov = PANELS[panel]
        growers = [r["growth"] for r in rows if r["grows"]]
        desc = describe(growers)
        nonopt = sum(1 for r in rows if r["status"] != "optimal")
        rec = {"panel": panel, "model_set": label, "media": media,
               "provenance": prov, "n_models": len(rows),
               "n_growers": len(growers),
               "frac_growers": round(len(growers) / len(rows), 4) if rows else None,
               "n_zero_growth": len(rows) - len(growers) - nonopt,
               "n_non_optimal": nonopt, **desc}
        summary["panels"][panel] = rec
        sum_rows.append(rec)
        dec = {"panel": panel, "media": media}
        for i in range(11):
            dec[f"p{i * 10}"] = quantile(sorted(growers), i / 10)
        dec_rows.append(dec)

    with open(out_dir / "biomass_summary.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(sum_rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(sum_rows)
    with open(out_dir / "biomass_deciles.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(dec_rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(dec_rows)

    # ---- taxonomy distribution per rank ----------------------------------
    # Counted over GENOMES, not models: the core and genome-scale sets are two
    # model families over (almost) the same genome list, so a per-model count
    # would double-count every genome that appears in both.
    genomes_by_panel = {p: {r["genome_id"] for r in rows}
                        for p, rows in panels.items()}
    genomes_by_panel["full"] = genomes_by_panel["full_gmm"] | genomes_by_panel["full_auxo"]
    grows_by = {p: {r["genome_id"] for r in rows if r["grows"]}
                for p, rows in panels.items()}

    summary["genome_sets"] = {
        "core": len(genomes_by_panel["core"]),
        "full": len(genomes_by_panel["full"]),
        "shared": len(genomes_by_panel["core"] & genomes_by_panel["full"]),
        "core_only": len(genomes_by_panel["core"] - genomes_by_panel["full"]),
        "full_only": len(genomes_by_panel["full"] - genomes_by_panel["core"]),
        "union": len(genomes_by_panel["core"] | genomes_by_panel["full"]),
    }

    summary["taxonomy"] = {}
    for rank in RANKS:
        counts: dict[str, dict] = defaultdict(lambda: {
            "core": 0, "full": 0, "union": 0,
            "core_growers": 0, "full_gmm_growers": 0, "full_auxo_growers": 0})
        for gid in genomes_by_panel["core"] | genomes_by_panel["full"]:
            taxon = (tax.get(gid, {}).get(rank) or "").strip() or "Unassigned"
            c = counts[taxon]
            c["union"] += 1
            if gid in genomes_by_panel["core"]:
                c["core"] += 1
                if gid in grows_by["core"]:
                    c["core_growers"] += 1
            if gid in genomes_by_panel["full"]:
                c["full"] += 1
                if gid in grows_by["full_gmm"]:
                    c["full_gmm_growers"] += 1
                if gid in grows_by["full_auxo"]:
                    c["full_auxo_growers"] += 1
        ordered = sorted(counts.items(), key=lambda kv: (-kv[1]["union"], kv[0]))
        path = out_dir / f"taxonomy_{rank}.tsv"
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow([rank, "n_genomes_union", "n_core", "n_full",
                        "core_growers", "core_grower_frac",
                        "full_gmm_growers", "full_gmm_grower_frac",
                        "full_auxo_growers", "full_auxo_grower_frac"])
            for taxon, c in ordered:
                w.writerow([
                    taxon, c["union"], c["core"], c["full"],
                    c["core_growers"],
                    round(c["core_growers"] / c["core"], 4) if c["core"] else "",
                    c["full_gmm_growers"],
                    round(c["full_gmm_growers"] / c["full"], 4) if c["full"] else "",
                    c["full_auxo_growers"],
                    round(c["full_auxo_growers"] / c["full"], 4) if c["full"] else "",
                ])
        summary["taxonomy"][rank] = {
            "n_distinct": len([t for t, _ in ordered if t != "Unassigned"]),
            "n_unassigned_genomes": counts.get("Unassigned", {}).get("union", 0),
            "top": [{"taxon": t, **c} for t, c in ordered[:25]],
        }

    # ---- biomass by phylum ------------------------------------------------
    by_phylum: dict[tuple, list] = defaultdict(list)
    tot: dict[tuple, int] = Counter()
    for panel, rows in panels.items():
        for r in rows:
            ph = (tax.get(r["genome_id"], {}).get("phylum") or "").strip() or "Unassigned"
            tot[(panel, ph)] += 1
            if r["grows"]:
                by_phylum[(panel, ph)].append(r["growth"])
    with open(out_dir / "biomass_by_phylum.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["panel", "phylum", "n_models", "n_growers", "grower_frac",
                    "min", "q1", "median", "mean", "q3", "max", "sd"])
        for (panel, ph), n in sorted(tot.items(), key=lambda kv: (kv[0][0], -kv[1])):
            vals = by_phylum.get((panel, ph), [])
            d = describe(vals)
            w.writerow([panel, ph, n, len(vals),
                        round(len(vals) / n, 4) if n else "",
                        *[("" if d[k] is None else round(d[k], 6))
                          for k in ("min", "q1", "median", "mean", "q3", "max", "sd")]])

    with open(out_dir / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=1)

    # ---- console report ---------------------------------------------------
    print(f"wrote {out_dir}/")
    print()
    print(f"{'panel':11s} {'n':>6s} {'growers':>8s} {'%':>6s} "
          f"{'min':>9s} {'q1':>9s} {'median':>9s} {'mean':>9s} {'q3':>9s} {'max':>9s}")
    for panel, rec in summary["panels"].items():
        print(f"{panel:11s} {rec['n_models']:6d} {rec['n_growers']:8d} "
              f"{100 * rec['frac_growers']:5.1f}% "
              f"{rec['min']:9.3f} {rec['q1']:9.3f} {rec['median']:9.3f} "
              f"{rec['mean']:9.3f} {rec['q3']:9.3f} {rec['max']:9.3f}")
    print()
    g = summary["genome_sets"]
    print(f"genomes: core {g['core']}, full {g['full']}, shared {g['shared']}, "
          f"core-only {g['core_only']}, full-only {g['full_only']}, union {g['union']}")
    for rank in ("superkingdom", "phylum", "genus", "species"):
        t = summary["taxonomy"][rank]
        print(f"  {rank:13s} {t['n_distinct']:5d} distinct"
              f"   (unassigned genomes: {t['n_unassigned_genomes']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
