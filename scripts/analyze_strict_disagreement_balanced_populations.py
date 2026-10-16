#!/usr/bin/env python3
"""Classify WHY the strict disagreements happen, restricted to two
status-defined populations that cannot blame a ModelSEED-flagged mass
imbalance for the conflict:

    OK_CI_only (879)  modelseed_status tokens == {'OK'} or == {'CI'}
                       (charge-imbalanced but mass-balanced is kept; mass-
                       imbalanced statuses MI / CI+MI are excluded)
    OK_only    (660)  modelseed_status tokens == {'OK'} exactly (no flag
                       at all, not even charge imbalance)

This reuses every piece of the 504-equation-sound analysis
(analyze_strict_disagreement_sound_subset.py / _residual_pathways.py) --
mechanism-tagging on direction_conflict.note, the 5 systematic
dGPredictor/eQuilibrator clusters, and the KEGG/MetaCyc pathway bucketing of
whatever's left -- but reruns it on these two status-filtered populations
instead of the literature-derived equation_defect=="none" filter alone, so
the same reaction can land in a different place here than in
sound_subset_analysis.json (that 504 includes 18 reactions ModelSEED's own
status flags MI/CI+MI, and conversely 120 of these 879 carry a
literature-found "mass_unbalanced" equation_defect despite ModelSEED's
automated checker saying OK/CI -- both are reported as their own finding,
not reconciled away).

No new research: every field read here was already collected by the
literature subagents behind strict_disagreement_enzymes.json.

Outputs (results/strict_disagreement_enzymes/):
    balanced_population_analysis.json   full nested breakdown, both populations
    balanced_population_clusters.tsv    cluster sizes for all 3 populations
                                         (all_1370, ok_ci_only, ok_only), figure data
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.dirname(HERE)
OUTDIR = os.path.join(ROOT, "results", "strict_disagreement_enzymes")
DELIVERABLE = os.path.join(OUTDIR, "strict_disagreement_enzymes.json")

from analyze_strict_disagreement_sound_subset import (
    tags_of, thermo, claim_systematic_clusters,
)
from analyze_strict_disagreement_residual_pathways import (
    pathway_text, fallback_text, categorize,
)

GROUPS = ("forward_vs_reverse", "reverse_vs_forward")
VERDICTS = ("llm", "modelseed", "neither", "unclear")
DEFECTS = ("none", "cofactors_stripped", "mass_unbalanced", "source_disowned", "lumped_multienzyme")


def status_tokens(status: str | None) -> frozenset:
    status = status or ""
    return frozenset(f[:f.find(':')] if ':' in f else f for f in status.split('|'))


def counter_dict(c: Counter) -> dict:
    return dict(c.most_common())


def verdict_block(reactions: list[dict]) -> dict:
    by_group = {g: Counter() for g in GROUPS}
    for r in reactions:
        dc = r["direction_conflict"]
        v = dc.get("enzymology_supports")
        if v:
            by_group[dc["mismatch_group"]][v] += 1
    pooled = Counter()
    for g in GROUPS:
        pooled.update(by_group[g])

    def block(c: Counter) -> dict:
        n = sum(c.values())
        return {"n": n, "counts": {v: c.get(v, 0) for v in VERDICTS},
                "pct": {v: round(100.0 * c.get(v, 0) / n, 1) if n else 0.0 for v in VERDICTS}}

    return {"by_group": {g: block(by_group[g]) for g in GROUPS}, "pooled": block(pooled)}


def mechanism_table(reactions: list[dict]) -> dict:
    tag_counts = Counter()
    tag_examples: dict[str, list[str]] = {}
    n_tagged = 0
    for r in reactions:
        tags = tags_of(r["direction_conflict"].get("note"))
        if tags:
            n_tagged += 1
        for tag in tags:
            tag_counts[tag] += 1
            tag_examples.setdefault(tag, [])
            if len(tag_examples[tag]) < 5:
                tag_examples[tag].append(r["modelseed_id"])
    n = len(reactions) or 1
    return {
        "n": len(reactions), "tagged_at_least_one": n_tagged,
        "tags": [{"mechanism": t, "n": c, "share": round(c / n, 4), "examples": tag_examples[t]}
                 for t, c in tag_counts.most_common()],
    }


def residual_pathway_table(residual: list[dict]) -> dict:
    if not residual:
        return {"n_residual": 0, "categories": []}
    assignments = []
    for r in residual:
        ptext = pathway_text(r)
        cat = categorize(ptext) if ptext else None
        prov = "pathway-annotated" if cat else None
        if cat is None:
            cat = categorize(fallback_text(r))
            prov = "text-inferred" if cat else None
        if cat is None:
            cat, prov = "No pathway identifiable", "none"
        assignments.append({"modelseed_id": r["modelseed_id"], "category": cat, "provenance": prov})
    by_id = {a["modelseed_id"]: a for a in assignments}
    cat_counts = Counter(a["category"] for a in assignments)
    cat_source: dict[str, Counter] = {}
    for r in residual:
        cat = by_id[r["modelseed_id"]]["category"]
        cat_source.setdefault(cat, Counter())[thermo(r).get("source")] += 1
    return {
        "n_residual": len(residual),
        "n_pathway_annotated": sum(1 for a in assignments if a["provenance"] == "pathway-annotated"),
        "n_text_inferred": sum(1 for a in assignments if a["provenance"] == "text-inferred"),
        "n_no_pathway_identifiable": sum(1 for a in assignments if a["provenance"] == "none"),
        "categories": [{"category": c, "n": n, "share_of_residual": round(n / len(residual), 4),
                         "thermo_source": counter_dict(cat_source[c])}
                        for c, n in cat_counts.most_common()],
    }


def analyze_population(reactions: list[dict], label: str) -> dict:
    n = len(reactions)
    defect_counts = Counter(r["direction_conflict"].get("equation_defect") for r in reactions)

    sound = [r for r in reactions if r["direction_conflict"].get("equation_defect") == "none"]
    llm_sound = [r for r in sound if r["direction_conflict"].get("enzymology_supports") == "llm"]
    ms_sound = [r for r in sound if r["direction_conflict"].get("enzymology_supports") == "modelseed"]

    clusters, residual = claim_systematic_clusters(llm_sound)
    claimed = len(llm_sound) - len(residual)

    return {
        "label": label,
        "n": n,
        "equation_defect": {
            "counts": {d: defect_counts.get(d, 0) for d in DEFECTS},
            "share": {d: round(100.0 * defect_counts.get(d, 0) / n, 1) if n else 0.0 for d in DEFECTS},
        },
        "correctness_verdict": verdict_block(reactions),
        "n_equation_sound": len(sound),
        "share_equation_sound": round(len(sound) / n, 4) if n else 0.0,
        "sound_correctness_verdict": verdict_block(sound),
        "n_llm_sound": len(llm_sound),
        "n_modelseed_sound": len(ms_sound),
        "sound_thermo_source": counter_dict(Counter(thermo(r).get("source") for r in sound)),
        "sound_thermo_grade": counter_dict(Counter(thermo(r).get("grade") for r in sound)),
        "sound_thermo_assessment": counter_dict(Counter(thermo(r).get("assessment") for r in sound)),
        "llm_sound_systematic_clusters": clusters,
        "llm_sound_clustered_total": claimed,
        "llm_sound_clustered_share": round(claimed / len(llm_sound), 4) if llm_sound else 0.0,
        "llm_sound_residual": {
            "n": len(residual),
            "thermo_source": counter_dict(Counter(thermo(r).get("source") for r in residual)),
            "thermo_grade": counter_dict(Counter(thermo(r).get("grade") for r in residual)),
            "confidence": counter_dict(Counter(r.get("confidence") for r in residual)),
        },
        "direction_mechanism_tags_llm_sound": mechanism_table(llm_sound),
        "residual_pathways": residual_pathway_table(residual),
    }


def main() -> None:
    deliverable = json.load(open(DELIVERABLE))
    reactions = deliverable["reactions"]

    ok_ci_only = [r for r in reactions if status_tokens(r.get("modelseed_status")) in ({"OK"}, {"CI"})]
    ok_only = [r for r in reactions if status_tokens(r.get("modelseed_status")) == frozenset({"OK"})]

    out = {
        "description": (
            "Classification of WHY the strict direction disagreements happen, restricted "
            "to reactions ModelSEED's own status field does not flag as mass-imbalanced."
        ),
        "n_total_1370": len(reactions),
        "populations": {
            "ok_ci_only": analyze_population(ok_ci_only, "OK or CI-only (879)"),
            "ok_only": analyze_population(ok_only, "OK only, no flags (660)"),
        },
    }

    with open(os.path.join(OUTDIR, "balanced_population_analysis.json"), "w") as fh:
        json.dump(out, fh, indent=2)

    # also re-derive the all-1370 clusters for a 3-population comparison TSV
    all_sound = [r for r in reactions if r["direction_conflict"].get("equation_defect") == "none"]
    all_llm_sound = [r for r in all_sound if r["direction_conflict"].get("enzymology_supports") == "llm"]
    all_clusters, all_residual = claim_systematic_clusters(all_llm_sound)

    with open(os.path.join(OUTDIR, "balanced_population_clusters.tsv"), "w") as fh:
        fh.write("population\tcluster\tn\tn_llm_sound_total\n")
        for pop_label, clusters, llm_sound_n, residual_n in (
            ("all_1370", all_clusters, len(all_llm_sound), len(all_residual)),
            ("ok_ci_only", out["populations"]["ok_ci_only"]["llm_sound_systematic_clusters"],
             out["populations"]["ok_ci_only"]["n_llm_sound"], out["populations"]["ok_ci_only"]["llm_sound_residual"]["n"]),
            ("ok_only", out["populations"]["ok_only"]["llm_sound_systematic_clusters"],
             out["populations"]["ok_only"]["n_llm_sound"], out["populations"]["ok_only"]["llm_sound_residual"]["n"]),
        ):
            for c in clusters:
                fh.write(f"{pop_label}\t{c['cluster']}\t{c['n']}\t{llm_sound_n}\n")
            fh.write(f"{pop_label}\tResidual (case-specific enzymology)\t{residual_n}\t{llm_sound_n}\n")

    for pop_key in ("ok_ci_only", "ok_only"):
        p = out["populations"][pop_key]
        print(f"=== {p['label']} ===")
        print(f"  equation_defect: {p['equation_defect']['counts']}  (sound={p['n_equation_sound']}, {p['share_equation_sound']:.1%})")
        cv = p["correctness_verdict"]["pooled"]
        print(f"  correctness (all {p['n']}): " + "  ".join(f"{v}={cv['counts'][v]} ({cv['pct'][v]}%)" for v in VERDICTS))
        scv = p["sound_correctness_verdict"]["pooled"]
        print(f"  correctness (sound {p['n_equation_sound']}): " + "  ".join(f"{v}={scv['counts'][v]} ({scv['pct'][v]}%)" for v in VERDICTS))
        print(f"  llm_sound={p['n_llm_sound']}  clustered={p['llm_sound_clustered_total']} ({p['llm_sound_clustered_share']:.1%})  residual={p['llm_sound_residual']['n']}")
        for c in p["llm_sound_systematic_clusters"]:
            print(f"    {c['cluster']:<42} n={c['n']:>3}  dominant={c['dominant_source']} (purity {c['dominant_source_purity']:.0%})")
        print()

    print(f"wrote {os.path.join(OUTDIR, 'balanced_population_analysis.json')}")
    print(f"wrote {os.path.join(OUTDIR, 'balanced_population_clusters.tsv')}")


if __name__ == "__main__":
    main()
