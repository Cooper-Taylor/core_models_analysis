#!/usr/bin/env python3
"""Sub-cluster the "Residual (case-specific enzymology)" bucket of the
strict-disagreement analysis, and test whether that bucket is the SAME set of
reactions across the three populations in
strict_disagreement_balanced_clusters.png.

Two questions, both left open by the sections this builds on
(analyze_strict_disagreement_sound_subset.py -> the 5 systematic clusters +
a 248-reaction residual; analyze_strict_disagreement_balanced_populations.py
-> the same thing re-run on two modelseed_status-defined populations):

  1. SET IDENTITY. The balanced-populations figure shows a residual of
     248 / 242 / 233 across all_1370 / ok_ci_only / ok_only and reads the
     stability of the numbers as independent corroboration. But each
     population is a ROW-SUBSET of the one before it (660 OK-only rows are a
     subset of the 879 OK-or-CI rows, which are a subset of all 1,370), and
     every filter applied downstream -- equation_defect == "none",
     enzymology_supports == "llm", the 5 cluster predicates -- is a pure
     per-reaction predicate. So the residuals are nested by construction, and
     this script proves it rather than asserting it: it emits the actual ID
     sets, the subset relations, and the identity of every dropped reaction.
     (The report's "neither is a subset of the 504" is true of the 879/660
     POPULATIONS, but not of the sound/LLM-correct/residual sets actually
     compared in that figure -- 476 <= 486 <= 504, gap 10 and 18.)

  2. SUB-STRUCTURE. The residual was characterised as having "no further
     single-mechanism pattern" because no further THERMO-SOURCE-pure cluster
     was found at the scale of the big five. Searched at a finer grain --
     enzyme-name tokens appearing >=4 times, ranked by thermo-source purity
     rather than guessed at -- ten sub-clusters do emerge, almost all 100%
     single-source, covering ~31% of the residual. They are smaller than the
     big five (4-19 reactions, not 35-81), which is why a coarse pass missed
     them; they are not weaker (purity 83-100%).

     The headline chemistry: the residual's sub-clusters are HYDROLYTIC and
     RING-REARRANGEMENT chemistry on xenobiotic substrates, where the big
     five were activated-cofactor group transfer (SAM/PAPS/ATP) and two
     structural motifs (glycoside, quinone). Two consecutive steps of the
     SAME pathway -- chlorocatechol modified-ortho degradation -- appear as
     two separate sub-clusters under two DIFFERENT thermo sources
     (chloromuconate cycloisomerase, 10/10 dGPredictor; dienelactone
     hydrolase, 5/5 eQuilibrator), which is hard to explain as anything but
     a per-source coverage gap on that pathway's compounds.

Two hypotheses tested here that did NOT hold, reported because they were
checked (the residual is NOT near-equilibrium chemistry, and it is NOT
mostly duplicate entries of a few enzymes) -- see `rejected_hypotheses`.

No new research and no ModelSEED biochemistry read: every field used was
already collected by the literature subagents behind
strict_disagreement_enzymes.json.

Outputs (results/strict_disagreement_enzymes/):
    residual_subclusters.json   nesting proof, sub-cluster membership per
                                 population, rejected hypotheses, dG stats
    residual_subclusters.tsv    sub-cluster sizes for all 3 populations
                                 (figure source data)
"""
from __future__ import annotations

import json
import os
import re
import statistics as st
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.dirname(HERE)
OUTDIR = os.path.join(ROOT, "results", "strict_disagreement_enzymes")
DELIVERABLE = os.path.join(OUTDIR, "strict_disagreement_enzymes.json")

from analyze_strict_disagreement_sound_subset import (
    SYSTEMATIC_CLUSTERS, claim_systematic_clusters,
)
from analyze_strict_disagreement_balanced_populations import status_tokens
from analyze_strict_disagreement_residual_pathways import (
    pathway_text, fallback_text, categorize,
)

# Ten sub-clusters, claimed in this priority order so a reaction counts once.
# Matched on enzyme_name + name (the residual's EC coverage is too sparse to
# cluster on -- only 84/248 carry a complete 4-field EC). Each was found by
# ranking every enzyme-name token with n>=4 by thermo-source purity, NOT by
# guessing chemistry first; the labels below name the chemistry those tokens
# turned out to share. `defect_class` groups them the way the big five are
# grouped, by which source produced the wrong confident call.
SUBCLUSTERS = [
    ("Chloromuconate cycloisomerase", "dGPredictor",
     r"cycloisomerase|muconate lactonis|muconate lactoniz"),
    ("Dienelactone hydrolase", "eQuilibrator",
     r"dienelactone|carboxymethylenebutenolidase"),
    ("Ring aminohydrolase / hydrolytic deaminase", "dGPredictor",
     r"aminohydrolase|deaminase|chlorohydrolase"),
    ("Triple-bond hydratase (alkyne / nitrile)", "dGPredictor",
     r"acetylen|alkynoate|nitrile hydratase|isonitrile hydratase|nitrilase|hydrase"),
    ("cis-Dihydrodiol dehydrogenase", "dGPredictor",
     r"dihydrodiol dehydrogenase|dihydro-[\d,]*-diol dehydrogenase|diol dehydrogenase"),
    ("Epoxide hydrolase", "dGPredictor",
     r"epoxide hydrolase"),
    ("Flavonol glycoside acyl / glucosyl transfer", "eQuilibrator",
     r"flavonol|kaempferol|o-glucoside"),
    ("Hydrolytic / reductive dehalogenase", "dGPredictor",
     r"dehalogenase|alkylhalidase"),
    ("HCCA isomerase (naphthalene, NahD)", "eQuilibrator",
     r"chromene-2-carboxylate|chromene"),
    ("Acyl-CoA dehydrogenase (beta-oxidation)", "eQuilibrator",
     r"acyl-coa dehydrogenase"),
]

SRCMAP = {"dGP": "dGPredictor", "eQ": "eQuilibrator",
          "GC": "Group contribution", "opentecr": "openTECR"}


def src(r: dict) -> str | None:
    return (r["direction_conflict"].get("thermo_evidence") or {}).get("source")


def name_text(r: dict) -> str:
    return ((r.get("enzyme_name") or "") + " " + (r.get("name") or "")).lower()


def driving_dg(r: dict) -> tuple[float | None, float | None]:
    """ΔG'° and uncertainty from the source that actually drove the call."""
    s = SRCMAP.get(src(r))
    v = (r["direction_conflict"].get("thermodynamics") or {}).get(s)
    if not v:
        return None, None
    try:
        return float(v[0]), float(v[1])
    except (TypeError, ValueError):
        return None, None


def topology(r: dict) -> str:
    """Water/cofactor shape of the equation -- the chemistry axis the five
    big clusters are blind to (they key on cofactor identity, not on whether
    the reaction is a hydrolysis at all)."""
    parts = re.split(r"\s*(?:=>|<=>|<=)\s*", r.get("reaction") or "")
    if len(parts) != 2:
        return "unparsed"
    lhs, rhs = parts[0].lower(), parts[1].lower()
    if re.search(r"\bnad\b|nadh|nadp|\bfad\b|\bcoa\b|\batp\b|\badp\b|\bo2\b", lhs + " " + rhs):
        return "cofactor-coupled"
    water_l, water_r = bool(re.search(r"\bh2o\b", lhs)), bool(re.search(r"\bh2o\b", rhs))
    if water_l and water_r:
        return "H2O on both sides"
    if water_l:
        return "hydrolysis (H2O substrate, no cofactor)"
    if water_r:
        return "condensation / dehydration (H2O product)"
    return "intramolecular (no H2O, no cofactor)"


def pathway_of(r: dict) -> str:
    ptext = pathway_text(r)
    cat = categorize(ptext) if ptext else None
    if cat:
        return cat
    return categorize(fallback_text(r)) or "No pathway identifiable"


def counter_dict(c: Counter) -> dict:
    return dict(c.most_common())


def residual_of(reactions: list[dict]) -> tuple[list[dict], list[dict], list[dict]]:
    sound = [r for r in reactions if r["direction_conflict"].get("equation_defect") == "none"]
    llm = [r for r in sound if r["direction_conflict"].get("enzymology_supports") == "llm"]
    _, residual = claim_systematic_clusters(llm)
    return sound, llm, residual


def subcluster(residual: list[dict]) -> tuple[list[dict], list[dict]]:
    """Priority-ordered, non-overlapping. Returns (clusters, leftover)."""
    seen: set[str] = set()
    out = []
    for label, defect, pat in SUBCLUSTERS:
        rx = re.compile(pat, re.I)
        grp = [r for r in residual
               if r["modelseed_id"] not in seen and rx.search(name_text(r))]
        seen.update(r["modelseed_id"] for r in grp)
        sc = Counter(src(r) for r in grp)
        top, topn = sc.most_common(1)[0] if grp else (None, 0)
        out.append({
            "cluster": label,
            "defect_class": defect,
            "n": len(grp),
            "dominant_source": top,
            "dominant_source_purity": round(topn / len(grp), 4) if grp else 0.0,
            "thermo_source": counter_dict(sc),
            "pathways": counter_dict(Counter(pathway_of(r) for r in grp)),
            "modelseed_ids": sorted(r["modelseed_id"] for r in grp),
        })
    leftover = [r for r in residual if r["modelseed_id"] not in seen]
    return out, leftover


def dg_block(group: list[dict]) -> dict:
    vals = [(abs(d), u) for d, u in (driving_dg(r) for r in group) if d is not None]
    if not vals:
        return {"n": len(group), "n_with_dg": 0}
    dgs = sorted(v[0] for v in vals)
    sigmas = [v[1] for v in vals if v[1] is not None]
    return {
        "n": len(group), "n_with_dg": len(dgs),
        "median_abs_dg": round(st.median(dgs), 2),
        "q1_abs_dg": round(dgs[len(dgs) // 4], 2),
        "q3_abs_dg": round(dgs[3 * len(dgs) // 4], 2),
        "n_near_equilibrium_abs_dg_lt_10": sum(1 for x in dgs if x < 10),
        "share_near_equilibrium": round(sum(1 for x in dgs if x < 10) / len(dgs), 4),
        "median_reported_sigma": round(st.median(sigmas), 2) if sigmas else None,
    }


def main() -> None:
    reactions = json.load(open(DELIVERABLE))["reactions"]

    pops = {
        "all_1370": reactions,
        "ok_ci_only": [r for r in reactions
                       if status_tokens(r.get("modelseed_status")) in ({"OK"}, {"CI"})],
        "ok_only": [r for r in reactions
                    if status_tokens(r.get("modelseed_status")) == frozenset({"OK"})],
    }

    per_pop, sets = {}, {}
    for key, rows in pops.items():
        sound, llm, residual = residual_of(rows)
        clusters, leftover = subcluster(residual)
        sets[key] = {
            "rows": {r["modelseed_id"] for r in rows},
            "sound": {r["modelseed_id"] for r in sound},
            "llm_sound": {r["modelseed_id"] for r in llm},
            "residual": {r["modelseed_id"] for r in residual},
        }
        claimed = sum(c["n"] for c in clusters)
        per_pop[key] = {
            "n_rows": len(rows), "n_sound": len(sound), "n_llm_sound": len(llm),
            "n_residual": len(residual),
            "subclusters": clusters,
            "n_subclustered": claimed,
            "share_of_residual_subclustered": round(claimed / len(residual), 4) if residual else 0.0,
            "n_leftover": len(leftover),
            "leftover_thermo_source": counter_dict(Counter(src(r) for r in leftover)),
            "residual_topology": counter_dict(Counter(topology(r) for r in residual)),
            "residual_pathways": counter_dict(Counter(pathway_of(r) for r in residual)),
        }

    # ---- 1. nesting proof -------------------------------------------------
    order = ["all_1370", "ok_ci_only", "ok_only"]
    nesting = {}
    for level in ("rows", "sound", "llm_sound", "residual"):
        chain = []
        for a, b in zip(order, order[1:]):
            sa, sb = sets[a][level], sets[b][level]
            chain.append({
                "outer": a, "inner": b,
                "inner_is_subset_of_outer": sb <= sa,
                "n_outer": len(sa), "n_inner": len(sb),
                "n_dropped": len(sa - sb),
                "n_gained_by_inner": len(sb - sa),
                "dropped_ids": sorted(sa - sb) if level == "residual" else None,
            })
        nesting[level] = chain

    by_id = {r["modelseed_id"]: r for r in reactions}
    dropped_detail = []
    for step in nesting["residual"]:
        for rid in step["dropped_ids"] or []:
            r = by_id[rid]
            dropped_detail.append({
                "modelseed_id": rid, "dropped_between": f"{step['outer']} -> {step['inner']}",
                "modelseed_status": r.get("modelseed_status"),
                "imbalance_flags": r.get("modelseed_imbalance_flags"),
                "thermo_source": src(r), "enzyme_name": r.get("enzyme_name"),
            })

    # ---- 2. rejected hypotheses ------------------------------------------
    _, llm_all, residual_all = residual_of(reactions)
    resid_ids = {r["modelseed_id"] for r in residual_all}
    clustered_all = [r for r in llm_all if r["modelseed_id"] not in resid_ids]

    def norm_name(r: dict) -> str:
        s = re.sub(r"\(.*?\)", "", (r.get("enzyme_name") or "").lower())
        s = re.sub(r"[^a-z\s-]", " ", s)
        s = re.sub(r"\b(the|a|an|of|type|family|like|activity|unidentified|unknown|"
                   r"putative|and|or|in|on|acting|its|variant|domain|protein|enzyme)\b", " ", s)
        return re.sub(r"\s+", " ", s).strip()

    name_counts = Counter(norm_name(r) for r in residual_all if norm_name(r))
    n_shared = sum(v for v in name_counts.values() if v >= 2)

    rejected = {
        "near_equilibrium_chemistry": {
            "hypothesis": ("The residual is reactions whose true |dG'o| is near zero, so "
                           "any estimator noise flips the sign -- which would make the "
                           "residual a precision problem, not a coverage problem."),
            "verdict": "REJECTED",
            "residual_dg": dg_block(residual_all),
            "clustered_dg": dg_block(clustered_all),
            "note": ("Median |dG'o| is essentially the same in both (residual vs the five "
                     "clusters), and only ~10% of the residual is under 10 kJ/mol. The "
                     "residual is NOT near-equilibrium chemistry. The one real difference "
                     "runs the other way: the residual's driving estimates carry a much "
                     "SMALLER reported uncertainty than the clustered ones, i.e. they are "
                     "wrong while claiming to be precise."),
        },
        "duplicate_enzyme_entries": {
            "hypothesis": ("The residual is inflated by many near-duplicate ModelSEED "
                           "entries of a few enzymes, so '248 individually documented "
                           "enzymes' overstates its diversity."),
            "verdict": "MOSTLY REJECTED",
            "n_residual": len(residual_all),
            "n_distinct_normalised_enzyme_names": len(name_counts),
            "n_sharing_a_name_with_another": n_shared,
            "share_sharing_a_name": round(n_shared / len(residual_all), 4),
            "largest_repeats": [{"enzyme": k, "n": v}
                                for k, v in name_counts.most_common(5) if v >= 3],
            "note": ("216 distinct normalised enzyme names across 248 reactions -- the "
                     "'individually documented' framing is ~85% accurate. The one real "
                     "exception is chloromuconate cycloisomerase (8 entries by name, 10 by "
                     "the cluster predicate), which is why it is the single largest "
                     "sub-cluster below."),
        },
    }

    out = {
        "description": ("Sub-structure of the 'Residual (case-specific enzymology)' bucket, "
                        "and a proof that the residual is the SAME nested set of reactions "
                        "across all three populations of "
                        "strict_disagreement_balanced_clusters.png."),
        "nesting": nesting,
        "residual_dropped_reactions": dropped_detail,
        "populations": per_pop,
        "rejected_hypotheses": rejected,
    }

    with open(os.path.join(OUTDIR, "residual_subclusters.json"), "w") as fh:
        json.dump(out, fh, indent=2)

    with open(os.path.join(OUTDIR, "residual_subclusters.tsv"), "w") as fh:
        fh.write("population\tcluster\tdefect_class\tn\tn_residual_total\t"
                 "dominant_source\tdominant_source_purity\n")
        for key in order:
            p = per_pop[key]
            for c in p["subclusters"]:
                fh.write(f"{key}\t{c['cluster']}\t{c['defect_class']}\t{c['n']}\t"
                         f"{p['n_residual']}\t{c['dominant_source']}\t"
                         f"{c['dominant_source_purity']:.4f}\n")
            fh.write(f"{key}\tStill unexplained\tnone\t{p['n_leftover']}\t"
                     f"{p['n_residual']}\t\t0\n")

    # ---- report ----------------------------------------------------------
    print("=== 1. is the residual the same set across populations? ===")
    for level in ("rows", "sound", "llm_sound", "residual"):
        for step in nesting[level]:
            print(f"  {level:<10} {step['inner']:<11} subset of {step['outer']:<11} "
                  f"{str(step['inner_is_subset_of_outer']):<5} "
                  f"({step['n_outer']} -> {step['n_inner']}, dropped {step['n_dropped']}, "
                  f"gained {step['n_gained_by_inner']})")
    print()
    print("=== 2. residual sub-clusters ===")
    hdr = f"{'sub-cluster':<44}{'all':>6}{'ok_ci':>7}{'ok':>6}   source (all_1370)"
    print(hdr)
    print("-" * len(hdr))
    for i, (label, _, _) in enumerate(SUBCLUSTERS):
        cells = [per_pop[k]["subclusters"][i]["n"] for k in order]
        c = per_pop["all_1370"]["subclusters"][i]
        print(f"  {label:<42}{cells[0]:>6}{cells[1]:>7}{cells[2]:>6}   "
              f"{c['dominant_source']} {c['dominant_source_purity']:.0%}")
    print("-" * len(hdr))
    print(f"  {'SUB-CLUSTERED':<42}" + "".join(
        f"{per_pop[k]['n_subclustered']:>{w}}" for k, w in zip(order, (6, 7, 6))))
    print(f"  {'still unexplained':<42}" + "".join(
        f"{per_pop[k]['n_leftover']:>{w}}" for k, w in zip(order, (6, 7, 6))))
    print(f"  {'RESIDUAL TOTAL':<42}" + "".join(
        f"{per_pop[k]['n_residual']:>{w}}" for k, w in zip(order, (6, 7, 6))))
    print()
    print("=== 3. rejected hypotheses ===")
    for k, v in rejected.items():
        print(f"  {k}: {v['verdict']}")
    rd, cd = rejected["near_equilibrium_chemistry"]["residual_dg"], \
        rejected["near_equilibrium_chemistry"]["clustered_dg"]
    print(f"    median |dG| residual={rd['median_abs_dg']} vs clustered={cd['median_abs_dg']}; "
          f"near-equilibrium {rd['share_near_equilibrium']:.1%} vs {cd['share_near_equilibrium']:.1%}; "
          f"median sigma {rd['median_reported_sigma']} vs {cd['median_reported_sigma']}")
    print()
    print("=== 4. residual chemistry topology (all_1370) ===")
    for k, v in per_pop["all_1370"]["residual_topology"].items():
        print(f"  {k:<46} {v:>3}")
    print()
    print(f"wrote {os.path.join(OUTDIR, 'residual_subclusters.json')}")
    print(f"wrote {os.path.join(OUTDIR, 'residual_subclusters.tsv')}")


if __name__ == "__main__":
    main()
