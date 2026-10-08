#!/usr/bin/env python3
"""Merge the per-chunk literature findings into the strict-disagreement deliverable.

Reads the research queue written by `build_strict_disagreement_queue.py` and
every `findings/chunk_NN.json` returned by the research subagents, joins them
on reaction id, and emits the deliverable plus summary tables.

The join is strict on purpose: a finding whose id is not in the queue is an
error (a subagent invented or mistyped an id), and a queue reaction with no
finding is reported as `missing` rather than silently dropped, so partial runs
are always visible as partial.

Outputs:
    results/strict_disagreement_enzymes/strict_disagreement_enzymes.json
        the deliverable -- every one of the 1,370 reactions with its
        ModelSEED identity, the reaction itself, the enzyme behind it, the
        mono/bi/polyfunctional call with evidence, and the direction conflict.
    results/strict_disagreement_enzymes/summary.json
        counts by functionality, confidence, tier, and mismatch group.
    results/strict_disagreement_enzymes/functionality_counts.tsv
        figure-ready counts.
"""
from __future__ import annotations

import collections
import csv
import glob
import re
import json
import os
import sys

SNAPSHOT = "/scratch/ctaylor/tmp/devsnap_078a395f/Biochemistry"

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUTDIR = os.path.join(ROOT, "results", "strict_disagreement_enzymes")
QUEUE = os.path.join(OUTDIR, "research_queue.json")
FINDINGS = os.path.join(OUTDIR, "findings")

CATEGORIES = ("monofunctional", "bifunctional", "polyfunctional", "unknown")
VERDICTS = ("llm", "modelseed", "neither", "unclear")
DEFECTS = ("none", "cofactors_stripped", "mass_unbalanced",
           "lumped_multienzyme", "source_disowned")
DIRECTION_WORD = {">": "forward", "<": "reverse", "=": "reversible", "?": "unspecified"}


def _snapshot_records():
    """id -> the full snapshot reaction record."""
    out = {}
    for path in sorted(glob.glob(os.path.join(SNAPSHOT, "reaction_*.json"))):
        with open(path) as fh:
            for rec in json.load(fh):
                out[rec["id"]] = rec
    return out


def _obsolete_ids():
    """Reaction ids the snapshot marks obsolete (excluded from baselines)."""
    out = set()
    for path in sorted(glob.glob(os.path.join(SNAPSHOT, "reaction_*.json"))):
        with open(path) as fh:
            for rec in json.load(fh):
                if rec.get("is_obsolete"):
                    out.add(rec["id"])
    return out


def load_status():
    """id -> ModelSEED's own `status` string for the reaction.

    ModelSEED already audits its equations and records the verdict here:
    "OK", or codes such as "MI:..." (mass imbalanced) and "CI:..." (charge
    imbalanced). Carrying it into the deliverable gives the literature-derived
    `equation_defect` an independent, deterministic check from the database
    itself rather than resting on the research agents' reading alone.
    """
    status = {}
    for path in sorted(glob.glob(os.path.join(SNAPSHOT, "reaction_*.json"))):
        with open(path) as fh:
            for rec in json.load(fh):
                status[rec["id"]] = rec.get("status") or ""
    return status


def load_backfill():
    """id -> {direction_verdict, equation_defect} for chunks researched before
    those two fields were added to the protocol.

    Chunks 00-10 were dispatched under the original protocol, which asked only
    for a free-text `direction_note`. Rather than re-run ~$19 of literature
    search, a later pass re-read those findings and classified the text they
    already contained into the two structured fields; the result lands in
    `findings/backfill_*.json` and is applied here only where the research
    chunk itself left the field empty.
    """
    back = {}
    for path in sorted(glob.glob(os.path.join(FINDINGS, "backfill_*.json"))):
        with open(path) as fh:
            for rec in json.load(fh):
                if rec.get("id"):
                    back[rec["id"]] = rec
    return back


def load_findings():
    """id -> finding dict, plus the list of problems found while loading."""
    found, problems = {}, []
    for path in sorted(glob.glob(os.path.join(FINDINGS, "chunk_*.json"))):
        try:
            with open(path) as fh:
                records = json.load(fh)
        except json.JSONDecodeError as exc:
            problems.append(f"{os.path.basename(path)}: unparseable JSON ({exc})")
            continue
        if not isinstance(records, list):
            problems.append(f"{os.path.basename(path)}: top level is not a list")
            continue
        for rec in records:
            rid = rec.get("id")
            if not rid:
                problems.append(f"{os.path.basename(path)}: record with no id")
                continue
            if rid in found:
                problems.append(f"{rid}: duplicate finding (later copy ignored)")
                continue
            cat = rec.get("functionality")
            if cat not in CATEGORIES:
                problems.append(f"{rid}: bad functionality {cat!r} -> coerced to 'unknown'")
                rec["functionality"] = "unknown"
            v = rec.get("direction_verdict")
            if v is not None and v not in VERDICTS:
                problems.append(f"{rid}: bad direction_verdict {v!r} -> cleared")
                rec["direction_verdict"] = ""
            d = rec.get("equation_defect")
            if d is not None and d not in DEFECTS:
                problems.append(f"{rid}: bad equation_defect {d!r} -> cleared")
                rec["equation_defect"] = ""
            found[rid] = rec
    return found, problems


def parse_mi(status):
    """The MI: atom deltas of a ModelSEED status string, as {element: delta}."""
    for part in (status or "").split("|"):
        if part.startswith("MI:"):
            out = {}
            for tok in part[3:].split("/"):
                if ":" in tok:
                    el, val = tok.rsplit(":", 1)
                    try:
                        out[el] = int(val)
                    except ValueError:
                        return {}
            return out
    return {}


def orphan_adp_multiple(status):
    """How many whole ADP molecules are missing from one side, or 0.

    A ModelSEED equation that lost an ADP leaves the exact signature
    C:-10k / N:-5k / O:-10k / P:-2k (H varies because the charge balancer
    adds free protons). This is a mechanically detectable curation defect
    rather than an enzymological disagreement, and it inflates the reaction's
    estimated dG by hundreds of kJ/mol -- which is what flips the direction.
    """
    d = parse_mi(status)
    if not d:
        return 0
    for k in range(1, 9):
        if (d.get("C"), d.get("N"), d.get("O"), d.get("P")) == (-10 * k, -5 * k, -10 * k, -2 * k):
            return k
    return 0


EXPANSION_TAG = re.compile(r"\.([a-z]+exp)\.")


def expansion_tag(record):
    """The combinatorial-expansion namespace of a reaction, e.g. 'metaexp'.

    ModelSEED instantiates class-level source reactions ("an aldehyde + NAD+
    = a carboxylate + NADH") by crossing the member lists, and several of the
    resulting entries pair a substrate with a product that is not its own
    reaction partner. The namespace is recorded in the abbreviation as
    ".<source>exp.".
    """
    text = f"{record.get('abbreviation') or ''} {record.get('name') or ''}"
    m = EXPANSION_TAG.search(text)
    return m.group(1) if m else None


def defect_signatures(status_by_id, strict_ids, obsolete, records_by_id=None):
    """Mechanical, database-native defect classes in the strict set."""
    nonobs = [r for r in status_by_id if r not in obsolete]
    strict_adp = [r for r in strict_ids if orphan_adp_multiple(status_by_id.get(r, ""))]
    base_adp = [r for r in nonobs if orphan_adp_multiple(status_by_id.get(r, ""))]
    sig = collections.Counter()
    for rid in strict_ids:
        d = parse_mi(status_by_id.get(rid, ""))
        if d:
            sig[", ".join(f"{k}:{v}" for k, v in sorted(d.items()))] += 1
    exp_strict = collections.Counter()
    exp_base = collections.Counter()
    if records_by_id:
        for rid in strict_ids:
            t = expansion_tag(records_by_id.get(rid, {}))
            if t:
                exp_strict[t] += 1
        for rid, rec in records_by_id.items():
            if rid in obsolete:
                continue
            t = expansion_tag(rec)
            if t:
                exp_base[t] += 1
    n_exp, n_exp_base = sum(exp_strict.values()), sum(exp_base.values())
    return {
        "combinatorial_expansions": {
            "in_strict_set": n_exp,
            "by_namespace": dict(exp_strict.most_common()),
            "in_whole_database": n_exp_base,
            "pct_of_strict_set": round(100.0 * n_exp / len(strict_ids), 1) if strict_ids else 0.0,
            "pct_of_database": round(100.0 * n_exp_base / max(1, len(status_by_id) - len(obsolete)), 2),
            # Only mildly enriched: expansions are not themselves a cause of
            # strict disagreement, but the ones that land here are reliably
            # defective -- see the per-chunk findings.
            "enrichment_vs_baseline": round(
                (n_exp / len(strict_ids)) /
                (n_exp_base / max(1, len(status_by_id) - len(obsolete))), 1)
            if n_exp_base and strict_ids else 0.0,
        },
        "orphan_adp": {
            "in_strict_set": len(strict_adp),
            "in_whole_database": len(base_adp),
            "share_of_all_database_cases_that_are_in_this_set": (
                round(100.0 * len(strict_adp) / len(base_adp), 1) if base_adp else 0.0),
            "enrichment_vs_baseline": (
                round((len(strict_adp) / len(strict_ids)) / (len(base_adp) / len(nonobs)), 1)
                if base_adp and nonobs else 0.0),
            "reaction_ids": sorted(strict_adp),
        },
        "most_common_mass_imbalance_signatures": dict(sig.most_common(8)),
    }


def identify_dropped_compounds(strict_ids, status_by_id):
    """For each mass-imbalanced reaction, name the single compound that is missing.

    A ModelSEED `status` of e.g. MI:C:-20/H:-26/N:-6/O:-12/S:-2 says the
    equation is short exactly C20H26N6O12S2 -- which is GSSG, a compound that
    already exists in the database as cpd00111. These are not estimation
    failures or hard curation problems: a product was dropped on import, and
    the counterpart is sitting in the compound table. Matching the residual
    formula against every compound turns "this equation is broken" into a
    specific, mechanically checkable repair.

    Only exact whole-compound matches are reported, and H is ignored in the
    match because ModelSEED's balancer adds free protons to patch charge.
    """
    formulas = collections.defaultdict(list)
    for path in sorted(glob.glob(os.path.join(SNAPSHOT, "compound_*.json"))):
        with open(path) as fh:
            for rec in json.load(fh):
                f = rec.get("formula")
                if f:
                    formulas[f].append((rec["id"], rec.get("name", "")))

    def canon(formula):
        """{element: count} from a formula string like 'C20H26N6O12S2'."""
        out = {}
        for el, num in re.findall(r"([A-Z][a-z]?)(\d*)", formula or ""):
            if el:
                out[el] = out.get(el, 0) + (int(num) if num else 1)
        out.pop("H", None)
        return tuple(sorted(out.items()))

    by_atoms = collections.defaultdict(list)
    for f, items in formulas.items():
        by_atoms[canon(f)].extend(items)

    found = collections.Counter()
    detail = {}
    for rid in strict_ids:
        d = parse_mi(status_by_id.get(rid, ""))
        if not d:
            continue
        neg = {el: -v for el, v in d.items() if el != "H" and v < 0}
        if not neg or any(v < 0 for v in d.values() if v > 0):
            pass
        if not neg:
            continue
        key = tuple(sorted(neg.items()))
        hits = by_atoms.get(key)
        if hits:
            cid, name = sorted(hits)[0]
            found[f"{cid} ({name})"] += 1
            detail[rid] = {"missing_compound": cid, "name": name}
    return {
        "reactions_with_an_identifiable_missing_compound": len(detail),
        "most_commonly_dropped": dict(found.most_common(10)),
        "per_reaction": detail,
    }


SAM_CPD = "cpd00017"  # S-Adenosyl-L-methionine


def sam_methyltransferase_sign(entries):
    """Is dGPredictor systematically wrong-signed on SAM methyl transfer?

    SAM -> SAH methyl transfer is strongly exergonic (order -60 kJ/mol): it is
    one of the most reliably one-way steps in metabolism, which is why nearly
    every methyltransferase in this set was judged "llm" on enzymology. If the
    estimator returns a positive dG for these, the resulting "<" recommendation
    is an artifact of the estimate, not a claim about the chemistry.

    The repeated-value tally matters as much as the sign: group-contribution
    style estimators reuse one class-level number across every reaction
    sharing a fragment pattern, so a single mis-signed group value flips a
    whole family of reactions at once rather than one.
    """
    hits = [e for e in entries if SAM_CPD in (e.get("equation") or "")]
    vals = []
    for e in hits:
        t = (e.get("thermodynamics") or {}).get("dGPredictor")
        if t and isinstance(t[0], (int, float)):
            vals.append(round(float(t[0]), 2))
    if not vals:
        return {"sam_reactions": len(hits), "with_estimate": 0}
    positive = [v for v in vals if v > 0]
    repeats = {str(v): n for v, n in collections.Counter(vals).most_common() if n > 1}
    return {
        "sam_reactions": len(hits),
        "with_estimate": len(vals),
        "positive_uphill": len(positive),
        "pct_positive": round(100.0 * len(positive) / len(vals), 1),
        "median_dg_kj_mol": round(sorted(vals)[len(vals) // 2], 2),
        "expected_sign": "negative (SAM->SAH methyl transfer is strongly exergonic)",
        "identical_values_reused": repeats,
    }


def metacyc_source_direction(entries):
    """Does ModelSEED's recommended direction contradict its own source?

    The pinned snapshot ships a MetaCyc dump under Provenance/MetaCyc/, whose
    `equation` column carries the arrow MetaCyc itself asserts: "=>" for an
    irreversible forward reaction, "<=>" for a reversible one. For every
    strict-disagreement reaction with a MetaCyc (or AraCyc) cross-reference we
    can therefore ask a question that needs no literature and no model
    judgment: when ModelSEED recommends "<", what did the source say?

    This is the independent check on the literature finding. The research
    agents concluded from enzymology that the LLM council is usually right;
    this asks the primary database the same question mechanically.
    """
    path = os.path.join(SNAPSHOT, "Provenance", "MetaCyc", "MetaCyc_reactions.tsv")
    if not os.path.exists(path):
        return {"available": False}
    eq = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rid = (row.get("id") or "").strip()
            if rid:
                eq[rid] = row.get("equation") or ""

    def arrow(text):
        if "<=>" in text:
            return "="
        if "=>" in text:
            return ">"
        if "<=" in text:
            return "<"
        return None

    table = collections.Counter()
    matched = 0
    for e in entries:
        xref = e.get("xrefs") or {}
        src = None
        for key in ("MetaCyc", "AraCyc"):
            for ident in xref.get(key) or []:
                if ident in eq:
                    src = arrow(eq[ident])
                    break
            if src:
                break
        if src is None:
            continue
        matched += 1
        table[(src, e["recommended_reversibility"], e["llm_direction"])] += 1

    if not matched:
        return {"available": True, "matched": 0}
    forward_overridden = sum(v for (src, rec, _), v in table.items()
                             if src == ">" and rec == "<")
    source_agrees_llm = sum(v for (src, _, llm), v in table.items() if src == llm)
    source_agrees_rec = sum(v for (src, rec, _), v in table.items() if src == rec)
    return {
        "available": True,
        "matched": matched,
        "of_total": len(entries),
        "source_writes_forward_modelseed_recommends_reverse": forward_overridden,
        "pct_forward_overridden": round(100.0 * forward_overridden / matched, 1),
        "source_direction_agrees_with_llm": source_agrees_llm,
        "pct_source_agrees_with_llm": round(100.0 * source_agrees_llm / matched, 1),
        "source_direction_agrees_with_recommendation": source_agrees_rec,
        "pct_source_agrees_with_recommendation": round(100.0 * source_agrees_rec / matched, 1),
        "breakdown_source_recommended_llm": {
            f"{src}|{rec}|{llm}": v for (src, rec, llm), v in table.most_common()
        },
    }


def imbalance_enrichment(status_by_id, strict_ids):
    """Is the strict-disagreement set enriched for equations ModelSEED itself
    flags as mass/charge imbalanced?

    This is the one measure in the report that needs no literature and no model
    judgment: ModelSEED audits its own equations and stores the verdict in
    `status`. Comparing the strict set against the whole non-obsolete database
    -- and against the wider recommended-vs-LLM mismatch set -- says whether
    broken stoichiometry travels with strict direction conflicts specifically,
    or is just background noise.
    """
    snap_obsolete = set()
    for path in sorted(glob.glob(os.path.join(SNAPSHOT, "reaction_*.json"))):
        with open(path) as fh:
            for rec in json.load(fh):
                if rec.get("is_obsolete"):
                    snap_obsolete.add(rec["id"])

    def flagged(rid):
        return bool(status_flags(status_by_id.get(rid, "")))

    nonobs = [r for r in status_by_id if r not in snap_obsolete]
    mismatch_path = os.path.join(ROOT, "results", "llm_vs_recommended",
                                 "mismatches_by_group.json")
    with open(mismatch_path) as fh:
        groups = json.load(fh)
    all_mismatch = {x["id"] for g in groups.values() for x in g}

    def rate(ids):
        ids = [i for i in ids if i in status_by_id]
        hit = sum(1 for i in ids if flagged(i))
        return {"n": len(ids), "flagged": hit,
                "pct": round(100.0 * hit / len(ids), 1) if ids else 0.0}

    baseline = rate(nonobs)
    strict = rate(strict_ids)
    out = {
        "all_non_obsolete": baseline,
        "all_recommended_vs_llm_mismatches": rate(all_mismatch),
        "strict_disagreements": strict,
        "enrichment_vs_baseline": round(
            (strict["pct"] / baseline["pct"]) if baseline["pct"] else 0.0, 2),
        "flag_codes_in_strict_set": dict(collections.Counter(
            code for i in strict_ids for code in status_flags(status_by_id.get(i, ""))
        )),
    }
    return out


def status_flags(status):
    """The imbalance codes in a ModelSEED status string, e.g. ['MI', 'CI']."""
    if not status or status == "OK":
        return []
    return sorted({part.split(":", 1)[0].strip()
                   for part in status.split("|") if ":" in part})


def merge(entry, finding, status=""):
    """One deliverable record: ModelSEED identity + reaction + enzyme + conflict."""
    f = finding or {}
    rec = finding.get("functionality") if finding else None
    return {
        # --- identity -------------------------------------------------
        "modelseed_id": entry["id"],
        "name": entry["modelseed_name"] or entry["search_name"] or "(unnamed in ModelSEED)",
        "synonyms": entry["synonyms"],
        "abbreviation": entry["abbreviation"],
        # --- the reaction ---------------------------------------------
        "reaction": entry["definition"],          # human-readable compound names
        "reaction_equation": entry["equation"],   # cpd ids
        "ec_numbers": entry["ec_numbers"],
        "xrefs": entry["xrefs"],
        "pathways": entry["pathways"],
        "is_transport": entry["is_transport"],
        "modelseed_status": status,
        "modelseed_imbalance_flags": status_flags(status),
        # --- the enzyme -----------------------------------------------
        "enzyme_name": f.get("enzyme_name", ""),
        "gene_symbols": f.get("gene_symbols", []),
        "functionality": rec or "not_researched",
        "distinct_activities": f.get("distinct_activities", []),
        "organism_context": f.get("organism_context", ""),
        "confidence": f.get("confidence", ""),
        "rationale": f.get("rationale", ""),
        "evidence_pmids": f.get("evidence_pmids", []),
        # --- why this reaction is in the set --------------------------
        "direction_conflict": {
            "mismatch_group": entry["mismatch_group"],
            "modelseed_recommended": DIRECTION_WORD.get(
                entry["recommended_reversibility"], entry["recommended_reversibility"]
            ),
            "llm_ensemble": DIRECTION_WORD.get(entry["llm_direction"], entry["llm_direction"]),
            "thermo_evidence": entry["thermo_evidence"],
            "thermodynamics": entry["thermodynamics"],
            "enzymology_supports": f.get("direction_verdict", ""),
            "equation_defect": f.get("equation_defect", ""),
            "note": f.get("direction_note", ""),
        },
        "research_tier": entry["research_tier"],
    }


def main():
    with open(QUEUE) as fh:
        queue = json.load(fh)
    entries = queue["reactions"]
    by_id = {e["id"]: e for e in entries}

    findings, problems = load_findings()
    backfill = load_backfill()
    n_filled = 0
    for rid, rec in findings.items():
        extra = backfill.get(rid)
        if not extra:
            continue
        for field in ("direction_verdict", "equation_defect"):
            if not rec.get(field) and extra.get(field):
                rec[field] = extra[field]
                n_filled += 1
    if backfill:
        problems.append(
            f"note: backfilled {n_filled} field(s) across {len(backfill)} "
            f"pre-protocol reactions from findings/backfill_*.json"
        )
    orphans = sorted(set(findings) - set(by_id))
    for rid in orphans:
        problems.append(f"{rid}: finding has no matching queue reaction (ignored)")

    status_by_id = load_status()
    records = [merge(e, findings.get(e["id"]), status_by_id.get(e["id"], "")) for e in entries]
    missing = [r["modelseed_id"] for r in records if r["functionality"] == "not_researched"]

    researched = [r for r in records if r["functionality"] != "not_researched"]
    summary = {
        "total_reactions": len(records),
        "researched": len(researched),
        "missing": len(missing),
        "by_functionality": dict(collections.Counter(r["functionality"] for r in records)),
        "by_confidence": dict(collections.Counter(r["confidence"] for r in researched if r["confidence"])),
        "by_tier": dict(collections.Counter(r["research_tier"] for r in records)),
        "by_mismatch_group": dict(
            collections.Counter(r["direction_conflict"]["mismatch_group"] for r in records)
        ),
        "with_direction_note": sum(1 for r in records if r["direction_conflict"]["note"]),
        "by_enzymology_supports": dict(collections.Counter(
            r["direction_conflict"]["enzymology_supports"] for r in researched
            if r["direction_conflict"]["enzymology_supports"]
        )),
        "by_equation_defect": dict(collections.Counter(
            r["direction_conflict"]["equation_defect"] for r in researched
            if r["direction_conflict"]["equation_defect"]
        )),
        "modelseed_flags_imbalance": sum(1 for r in records if r["modelseed_imbalance_flags"]),
        # Does the literature-derived defect call line up with ModelSEED's own audit?
        "defect_vs_modelseed_status": {
            d: {
                "modelseed_flags_imbalance": sum(
                    1 for r in researched
                    if r["direction_conflict"]["equation_defect"] == d and r["modelseed_imbalance_flags"]
                ),
                "modelseed_says_ok": sum(
                    1 for r in researched
                    if r["direction_conflict"]["equation_defect"] == d and not r["modelseed_imbalance_flags"]
                ),
            }
            for d in DEFECTS
        },
        # How often does ModelSEED hand a CONFIDENT thermodynamic grade to a
        # reaction whose own equation it has already flagged as imbalanced?
        # The grade describes the number, not the chemistry, so a gold grade on
        # a mass/charge-broken equation is confidence in an answer to the wrong
        # question. Deterministic -- both fields are ModelSEED's own.
        "grade_vs_imbalance": {
            grade: {
                "total": sum(1 for r in records
                             if (r["direction_conflict"]["thermo_evidence"] or {}).get("grade") == grade),
                "modelseed_flags_imbalance": sum(
                    1 for r in records
                    if (r["direction_conflict"]["thermo_evidence"] or {}).get("grade") == grade
                    and r["modelseed_imbalance_flags"]),
            }
            for grade in ("gold", "silver", "bronze")
        },
        # The two strict groups are NOT interchangeable. `reverse_vs_forward`
        # (1,308) and `forward_vs_reverse` (62) have different causes, so a
        # pooled "the LLM is right N% of the time" would be misleading --
        # always read the verdict split per group.
        "verdict_by_group": {
            g: dict(collections.Counter(
                r["direction_conflict"]["enzymology_supports"] for r in researched
                if r["direction_conflict"]["mismatch_group"] == g
                and r["direction_conflict"]["enzymology_supports"]
            ))
            for g in ("forward_vs_reverse", "reverse_vs_forward")
        },
        "defect_by_group": {
            g: dict(collections.Counter(
                r["direction_conflict"]["equation_defect"] for r in researched
                if r["direction_conflict"]["mismatch_group"] == g
                and r["direction_conflict"]["equation_defect"]
            ))
            for g in ("forward_vs_reverse", "reverse_vs_forward")
        },
        "modelseed_imbalance_enrichment": imbalance_enrichment(
            status_by_id, [e["id"] for e in entries]),
        "metacyc_source_direction": metacyc_source_direction(entries),
        "sam_methyltransferase_sign": sam_methyltransferase_sign(entries),
        "dropped_compounds": identify_dropped_compounds(
            [e["id"] for e in entries], status_by_id),
        "mechanical_defect_signatures": defect_signatures(
            status_by_id, [e["id"] for e in entries], _obsolete_ids(),
            _snapshot_records()),
        "verdict_backfill_needed": sorted(
            r["modelseed_id"] for r in researched
            if not r["direction_conflict"]["enzymology_supports"]
        ),
        "multifunctional_total": sum(
            1 for r in records if r["functionality"] in ("bifunctional", "polyfunctional")
        ),
        "problems": problems,
    }
    # Cross-tab: does being multifunctional track with which way the conflict runs?
    summary["functionality_by_group"] = {
        g: dict(collections.Counter(
            r["functionality"] for r in records if r["direction_conflict"]["mismatch_group"] == g
        ))
        for g in ("forward_vs_reverse", "reverse_vs_forward")
    }

    deliverable = {
        "description": (
            "Enzyme functionality (mono/bi/polyfunctional) behind the 1,370 ModelSEED "
            "reactions where the recommended reversibility and the LLM ensemble's "
            "direction call point in strictly opposite directions."
        ),
        "source_report": "reports/llmVsRecommended/LLM_VS_RECOMMENDED_DIRECTION.md",
        "modelseed_snapshot": "dev @ 078a395f",
        "mismatch_groups_included": ["forward_vs_reverse", "reverse_vs_forward"],
        "summary": summary,
        "reactions": records,
    }

    with open(os.path.join(OUTDIR, "strict_disagreement_enzymes.json"), "w") as fh:
        json.dump(deliverable, fh, indent=1)
    with open(os.path.join(OUTDIR, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    with open(os.path.join(OUTDIR, "functionality_counts.tsv"), "w") as fh:
        fh.write("category\treactions\n")
        for cat in CATEGORIES + ("not_researched",):
            fh.write(f"{cat}\t{summary['by_functionality'].get(cat, 0)}\n")

    print(json.dumps({k: v for k, v in summary.items() if k != "problems"}, indent=1))
    if problems:
        print(f"\n{len(problems)} problem(s):", file=sys.stderr)
        for p in problems[:40]:
            print("  " + p, file=sys.stderr)
    if missing:
        print(f"\nstill unresearched: {len(missing)} reactions "
              f"(first few: {', '.join(missing[:8])})", file=sys.stderr)


if __name__ == "__main__":
    main()
