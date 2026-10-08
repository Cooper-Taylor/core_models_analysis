#!/usr/bin/env python3
"""Build the per-reaction literature-research queue for the strict-disagreement set.

`reports/llmVsRecommended/LLM_VS_RECOMMENDED_DIRECTION.md` splits the 8,466
recommended-vs-LLM direction mismatches into six groups. Two of them are
*strict* disagreements -- the recommendation and the LLM council pick opposite
sides rather than merely disagreeing about whether the reaction is
bidirectional:

    forward_vs_reverse    recommended ">"  LLM "<"      62
    reverse_vs_forward    recommended "<"  LLM ">"   1,308
                                                    ------
                                                     1,370

This script assembles those 1,370 reactions into a research queue. The unit of
work is deliberately the REACTION, not the EC number: a separate EC-keyed
pipeline already exists (`build_enzyme_functionality_queue.py` and the
`results/enzyme_functionality/` tree) and is NOT reused here. Each reaction is
researched on its own -- which enzyme actually catalyses it, and is that enzyme
bifunctional or polyfunctional -- so that a reaction can be matched to the
specific enzyme that runs it rather than inheriting a label from an EC class
shared with unrelated chemistry.

Every queue entry carries the identifiers a literature search can work from:
the ModelSEED name and any `Name:` synonyms, the human-readable `definition`,
the raw `ec_numbers` (complete or partial -- context only, never the key), and
the KEGG / MetaCyc / BiGG / Rhea cross-references, which are the main handle
for the many reactions whose ModelSEED `name` is just a bare KEGG R-number.

Reactions are tiered only by how much identifying information exists, so the
research prompt can be adapted and so coverage is honestly reportable:

    named          ModelSEED name or synonym names a real activity
    xref_only      no usable name, but a KEGG / MetaCyc / BiGG / Rhea id to chase
    opaque         neither -- nothing to search on

Data: ModelSEED dev @ 078a395f (/scratch/ctaylor/tmp/devsnap_078a395f).

Outputs:
    results/strict_disagreement_enzymes/research_queue.json
        {"reactions": [entry, ...], "counts": {...}} -- all 1,370, in group
        then id order, each with every identifier found.
    results/strict_disagreement_enzymes/chunks/chunk_NN.json
        the same entries split into fixed-size chunks, one per research
        subagent, so a re-run hands out identical work.
"""
from __future__ import annotations

import glob
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SNAPSHOT = "/scratch/ctaylor/tmp/devsnap_078a395f/Biochemistry"
MISMATCHES = os.path.join(ROOT, "results", "llm_vs_recommended", "mismatches_by_group.json")
OUTDIR = os.path.join(ROOT, "results", "strict_disagreement_enzymes")

STRICT_GROUPS = ("forward_vs_reverse", "reverse_vs_forward")
CHUNK_SIZE = 25

# A ModelSEED `name` that is really just an identifier echoed back: a bare KEGG
# R-number, the reaction's own rxn id, a MetaCyc RXN- handle (optionally with a
# ".c" compartment or ".metaexp.*" suffix), or empty.
PLACEHOLDER_NAME = re.compile(
    r"^(?:\s*|R\d{5,}|rxn\d{5,}|RXN[-0-9A-Za-z]*(?:\.[0-9A-Za-z._]+)?)$"
)
XREF_PREFIXES = ("KEGG:", "MetaCyc:", "BiGG:", "Rhea:", "AraCyc:", "MetaNetX:", "BrachyCyc:")


def load_snapshot():
    """id -> full dev-snapshot reaction record."""
    snap = {}
    for path in sorted(glob.glob(os.path.join(SNAPSHOT, "reaction_*.json"))):
        with open(path) as fh:
            for rec in json.load(fh):
                snap[rec["id"]] = rec
    return snap


def synonyms(record):
    """The `Name:` alias line, split into individual synonym strings."""
    for alias in record.get("aliases") or []:
        if alias.startswith("Name:"):
            return [s.strip() for s in alias[5:].split(";") if s.strip()]
    return []


def xrefs(record):
    """{'KEGG': ['R00004'], 'MetaCyc': [...], ...} from the alias block."""
    out = {}
    for alias in record.get("aliases") or []:
        for prefix in XREF_PREFIXES:
            if alias.startswith(prefix):
                db = prefix.rstrip(":")
                ids = [s.strip() for s in alias[len(prefix):].split(";") if s.strip()]
                if ids:
                    out.setdefault(db, []).extend(ids)
    return out


def best_name(mismatch, record):
    """The most informative human-readable name, or '' if only placeholders."""
    name = (mismatch.get("name") or "").strip()
    if name and not PLACEHOLDER_NAME.match(name):
        return name
    for syn in synonyms(record):
        if not PLACEHOLDER_NAME.match(syn):
            return syn
    return ""


def tier_of(name, xref):
    if name:
        return "named"
    return "xref_only" if xref else "opaque"


def main():
    snap = load_snapshot()
    with open(MISMATCHES) as fh:
        groups = json.load(fh)

    entries = []
    for group in STRICT_GROUPS:
        for mismatch in sorted(groups[group], key=lambda r: r["id"]):
            record = snap.get(mismatch["id"], {})
            name = best_name(mismatch, record)
            xref = xrefs(record)
            syns = synonyms(record)
            entries.append({
                "id": mismatch["id"],
                "mismatch_group": group,
                "research_tier": tier_of(name, xref),
                "modelseed_name": (mismatch.get("name") or "").strip(),
                "search_name": name,
                "synonyms": [s for s in syns if s != name][:8],
                "abbreviation": mismatch.get("abbreviation") or "",
                "equation": mismatch.get("equation") or "",
                "definition": mismatch.get("definition") or "",
                "ec_numbers": mismatch.get("ec_numbers") or [],
                "xrefs": xref,
                "pathways": (record.get("pathways") or [])[:4],
                "is_transport": bool(mismatch.get("is_transport")),
                "recommended_reversibility": mismatch.get("recommended_reversibility"),
                "llm_direction": mismatch.get("llm_direction"),
                "thermo_evidence": mismatch.get("thermo_evidence") or {},
                "thermodynamics": mismatch.get("thermodynamics") or {},
            })

    counts = {
        "total": len(entries),
        "by_group": {g: sum(1 for e in entries if e["mismatch_group"] == g) for g in STRICT_GROUPS},
        "by_tier": {
            t: sum(1 for e in entries if e["research_tier"] == t)
            for t in ("named", "xref_only", "opaque")
        },
        "with_any_ec": sum(1 for e in entries if e["ec_numbers"]),
        "transport": sum(1 for e in entries if e["is_transport"]),
    }

    os.makedirs(os.path.join(OUTDIR, "chunks"), exist_ok=True)
    with open(os.path.join(OUTDIR, "research_queue.json"), "w") as fh:
        json.dump({"counts": counts, "reactions": entries}, fh, indent=1)

    # Clear stale chunks so a re-run after a queue change can't leave orphans
    # that a later ingest would silently pick up.
    for stale in glob.glob(os.path.join(OUTDIR, "chunks", "chunk_*.json")):
        os.remove(stale)
    n_chunks = 0
    for i in range(0, len(entries), CHUNK_SIZE):
        chunk = entries[i:i + CHUNK_SIZE]
        path = os.path.join(OUTDIR, "chunks", f"chunk_{i // CHUNK_SIZE:02d}.json")
        with open(path, "w") as fh:
            json.dump(chunk, fh, indent=1)
        n_chunks += 1

    print(json.dumps(counts, indent=1))
    print(f"\nwrote {n_chunks} chunks of <= {CHUNK_SIZE} to {OUTDIR}/chunks/")


if __name__ == "__main__":
    main()
