#!/usr/bin/env python3
"""Fill the taxonomy gaps left by ``results/ncbi_taxonomy.csv``.

That table was built against the NCBI *datasets* API, which serves only
assemblies that are still live in RefSeq. 122 of the 5,683 core-model
accessions have since been suppressed or replaced (E. coli O157:H7 Sakai,
Dechloromonas aromatica RCB, ...), so they came back blank. One further
accession -- GCF_000009085.1, Campylobacter jejuni NCTC 11168 -- is present in
the ModelSEED v2 genome-scale set but absent from the core set, so it has no
row in that table at all.

Entrez (eutils) keeps suppressed assemblies, so this script goes back through
``esearch``/``esummary`` on ``db=assembly`` to recover the organism name and
tax id, then ``efetch`` on ``db=taxonomy`` for the ranked lineage. Output is a
SEPARATE file, ``results/ncbi_taxonomy_fill.csv``, with the same columns as the
original -- the original is left untouched because several other scripts read
it and expect its exact contents.

The run is resumable: accessions already present in the fill file are skipped,
so an interrupted run can simply be repeated.

Usage::

    python3 scripts/resolve_missing_taxonomy.py            # resolve what is missing
    python3 scripts/resolve_missing_taxonomy.py --force    # ignore existing fill file
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
TAXONOMY_CSV = ANALYSIS_DIR / "results" / "ncbi_taxonomy.csv"
FILL_CSV = ANALYSIS_DIR / "results" / "ncbi_taxonomy_fill.csv"
GSM_MANIFEST = Path(os.environ.get(
    "MS2_GSM_DIR", "/scratch/ctaylor/modelseed2_gs_models")) / "gmm" / "manifest.tsv"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
# Entrez allows 3 requests/sec without an API key, 10/sec with one.
API_KEY = os.environ.get("NCBI_API_KEY", "")
DELAY = 0.12 if API_KEY else 0.36

FIELDS = ["assembly_accession", "current_accession", "organism_name", "tax_id",
          "assembly_name", "status", "superkingdom", "kingdom", "phylum",
          "class", "order", "family", "genus", "species"]
RANKS = ["superkingdom", "kingdom", "phylum", "class", "order", "family",
         "genus", "species"]


def _get(endpoint: str, params: dict) -> bytes:
    if API_KEY:
        params = {**params, "api_key": API_KEY}
    url = f"{EUTILS}/{endpoint}?{urllib.parse.urlencode(params)}"
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=60) as fh:
                data = fh.read()
            time.sleep(DELAY)
            return data
        except Exception as exc:  # noqa: BLE001 -- transient Entrez failures
            if attempt == 4:
                raise
            print(f"  retry {attempt + 1}/4 after {type(exc).__name__}: {exc}",
                  file=sys.stderr)
            time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")


def assembly_uids(accessions: list[str]) -> dict[str, str]:
    """accession -> assembly UID, via one esearch per accession.

    Batching with OR loses the accession->uid correspondence, and 120-odd
    single queries cost under a minute at the unauthenticated rate limit.
    """
    out: dict[str, str] = {}
    for i, acc in enumerate(accessions, 1):
        raw = _get("esearch.fcgi", {"db": "assembly", "term": acc,
                                    "retmode": "json"})
        import json as _json
        ids = _json.loads(raw)["esearchresult"].get("idlist", [])
        if ids:
            out[acc] = ids[0]
        print(f"  [{i}/{len(accessions)}] {acc} -> "
              f"{ids[0] if ids else 'NOT FOUND'}", file=sys.stderr)
    return out


def assembly_summaries(uid_by_acc: dict[str, str]) -> dict[str, dict]:
    """assembly UID -> {organism_name, tax_id, assembly_name, current_accession}."""
    uids = list(uid_by_acc.values())
    out: dict[str, dict] = {}
    for start in range(0, len(uids), 100):
        chunk = uids[start:start + 100]
        raw = _get("esummary.fcgi", {"db": "assembly", "id": ",".join(chunk)})
        root = ET.fromstring(raw)
        for doc in root.iter("DocumentSummary"):
            uid = doc.get("uid") or ""
            out[uid] = {
                "organism_name": (doc.findtext("Organism") or "").strip(),
                "tax_id": (doc.findtext("Taxid") or "").strip(),
                "assembly_name": (doc.findtext("AssemblyName") or "").strip(),
                "current_accession": (doc.findtext("AssemblyAccession") or "").strip(),
            }
    return out


def taxonomy_lineages(tax_ids: list[str]) -> dict[str, dict]:
    """tax id -> {rank: name} for the eight ranks we keep."""
    out: dict[str, dict] = {}
    uniq = sorted({t for t in tax_ids if t})
    for start in range(0, len(uniq), 100):
        chunk = uniq[start:start + 100]
        raw = _get("efetch.fcgi", {"db": "taxonomy", "id": ",".join(chunk),
                                   "retmode": "xml"})
        root = ET.fromstring(raw)
        for taxon in root.findall("Taxon"):
            tid = taxon.findtext("TaxId") or ""
            lineage = {}
            for node in taxon.findall("./LineageEx/Taxon"):
                rank = (node.findtext("Rank") or "").strip()
                if rank in RANKS:
                    lineage[rank] = (node.findtext("ScientificName") or "").strip()
            own_rank = (taxon.findtext("Rank") or "").strip()
            if own_rank in RANKS:
                lineage[own_rank] = (taxon.findtext("ScientificName") or "").strip()
            # NCBI renamed superkingdom -> domain in 2024; keep both populated.
            if "superkingdom" not in lineage:
                for node in taxon.findall("./LineageEx/Taxon"):
                    if (node.findtext("Rank") or "") == "domain":
                        lineage["superkingdom"] = (
                            node.findtext("ScientificName") or "").strip()
            out[tid] = lineage
    return out


def targets() -> list[str]:
    """Accessions needing a lineage: unresolved core rows + GSM-only genomes."""
    need: list[str] = []
    seen: set[str] = set()
    with open(TAXONOMY_CSV) as fh:
        for row in csv.DictReader(fh):
            seen.add(row["assembly_accession"])
            if row["status"] != "resolved" or not row["phylum"].strip():
                need.append(row["assembly_accession"])
    if GSM_MANIFEST.exists():
        with open(GSM_MANIFEST) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                gid = row.get("genome_id", "").strip()
                if gid and gid not in seen:
                    need.append(gid)
                    seen.add(gid)
    return sorted(set(need))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true",
                    help="re-resolve everything, ignoring the existing fill file")
    args = ap.parse_args()

    want = targets()
    done: dict[str, dict] = {}
    if FILL_CSV.exists() and not args.force:
        with open(FILL_CSV) as fh:
            done = {r["assembly_accession"]: r for r in csv.DictReader(fh)}
    todo = [a for a in want if a not in done]
    print(f"{len(want)} accessions need a lineage; {len(done)} already in "
          f"{FILL_CSV.name}; resolving {len(todo)}", file=sys.stderr)
    if not todo:
        print("nothing to do")
        return 0

    uid_by_acc = assembly_uids(todo)
    summaries = assembly_summaries(uid_by_acc)
    lineages = taxonomy_lineages([s["tax_id"] for s in summaries.values()])

    for acc in todo:
        uid = uid_by_acc.get(acc)
        summ = summaries.get(uid, {}) if uid else {}
        lin = lineages.get(summ.get("tax_id", ""), {})
        row = {f: "" for f in FIELDS}
        row["assembly_accession"] = acc
        row["current_accession"] = summ.get("current_accession", "")
        row["organism_name"] = summ.get("organism_name", "")
        row["tax_id"] = summ.get("tax_id", "")
        row["assembly_name"] = summ.get("assembly_name", "")
        row["status"] = "resolved" if lin else "missing"
        for rank in RANKS:
            row[rank] = lin.get(rank, "")
        done[acc] = row

    FILL_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(FILL_CSV, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS)
        writer.writeheader()
        for acc in sorted(done):
            writer.writerow({f: done[acc].get(f, "") for f in FIELDS})

    n_ok = sum(1 for r in done.values() if r["status"] == "resolved")
    print(f"wrote {FILL_CSV}  ({n_ok}/{len(done)} resolved)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
