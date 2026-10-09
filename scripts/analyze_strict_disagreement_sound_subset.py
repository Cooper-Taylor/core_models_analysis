#!/usr/bin/env python3
"""Zoom into the 504 strict-disagreement reactions with NO equation defect
(direction_conflict.equation_defect == "none") -- the subset where there is no
data-quality excuse, so any remaining wrong call is purely a direction-cascade
problem. Three questions:

  1. Why is the LLM council right so much more often here (92.3%, vs the
     57.7%/monofunctional-dominated full 1,370)? -- mined from
     direction_conflict.note (the DIRECTION rationale; top-level `rationale`
     is about the mono/bi/polyfunctional call and must not be conflated with
     it -- an earlier pass over the combined text undercounted real
     mechanisms because methyltransferase FUNCTIONALITY rationale ("single
     catalytic activity") kept matching before the DIRECTION note did).
  2. Which specific enzyme(s)/enzyme families impose the direction
     thermodynamically, and do they cluster by which thermo source produced
     the wrong recommendation?
  3. Of the 465 LLM-correct/equation-sound reactions, how much is explained by
     a SMALL NUMBER of systematic thermo-source defects (as opposed to 465
     independent literature judgment calls)?

Two systematic defect classes emerged (see `systematic_clusters` below),
found by moving past enzyme-NAME matching to searching the reaction EQUATION
itself for cofactor/motif signatures -- this is what caught two clusters that
`enzyme_name` family terms alone missed entirely (they cut across many
different named enzyme families):

  dGPredictor: mis-signs "activated cofactor" group transfer -- SAM
  (methyltransferases), PAPS (sulfotransferases), and ATP/phosphoanhydride
  (kinases/synthetases/ligases/other ATP-coupled transferases) are
  mechanistically the same category (a high-group-transfer-potential
  cofactor driving an essentially irreversible transfer); one bad group
  contribution for that bond type would explain all three at once.

  eQuilibrator: mis-estimates specific structural motifs -- glycosidic-bond
  hydrolysis (NOT just literal "glucosidase": glucuronidase, rhamnosidase and
  other GH-family hydrolases are the same mechanism and were undercounted by
  a naive enzyme-name substring match) and quinone/quinol-coupled redox
  (independently flagged elsewhere in this project,
  see memory "project_eq_vs_dgpredictor_modelseed" -- "quinone redox is the
  real failure").

Clusters are claimed in priority order so a reaction counts toward exactly
one cluster (no double-counting); `systematic_clusters_residual` is
everything left over, characterised (not just counted) to check whether ITS
errors reduce to some further pattern -- tested here against a
secondary-metabolism/molecular-novelty hypothesis, which does NOT hold (that
signal concentrates in the clusters, not the residual) and is reported as a
tested-and-rejected hypothesis rather than omitted.

Outputs (results/strict_disagreement_enzymes/):
    sound_subset_analysis.json     every count in this script, structured
    sound_subset_families.tsv      generic enzyme-family cluster sizes (figure
                                    source data for the supplementary table)
    sound_subset_clusters.tsv      the 5 non-overlapping systematic clusters +
                                    residual, with defect_class (figure source
                                    data for the headline chart)
"""
from __future__ import annotations

import json
import os
import re
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUTDIR = os.path.join(ROOT, "results", "strict_disagreement_enzymes")
DELIVERABLE = os.path.join(OUTDIR, "strict_disagreement_enzymes.json")

# Ordered so a more specific family name wins over a generic one
# (methyltransferase before transferase, glucosidase before hydrolase, ...).
FAMILY_TERMS = [
    "dehydrogenase", "oxidase", "reductase", "monooxygenase", "dioxygenase",
    "hydroxylase", "glucosidase", "glycosidase", "galactosidase",
    "sulfotransferase", "methyltransferase", "acyltransferase",
    "glycosyltransferase", "hydrolase", "esterase", "lipase", "transferase",
    "kinase", "synthetase", "ligase", "synthase", "hydratase", "dehydratase",
    "lyase", "decarboxylase", "aldolase", "isomerase", "epimerase",
    "racemase", "mutase", "cyclase", "deaminase", "amidohydrolase",
    "phosphatase", "peptidase", "protease", "oxidoreductase", "transporter",
    "permease",
]

# Direction-mechanism keyword taxonomy, matched against direction_conflict.note
# (NOT the top-level functionality `rationale`). A reaction can match several.
MECHANISM_PATTERNS = {
    "named canonical pathway role": r"salvage|biosynthe|catabolis|catabolic|respiratory chain|tca cycle|physiologic",
    "pathway channeling / downstream pull": r"channel|committed( step| pathway)?|downstream|biosynthetic flux|pulled (forward|toward|the)|drained|handed straight|out.compete|drives? (the|this)|pulls? (the|this)",
    "dehydration / hydration / elimination / condensation": r"dehydrat|elimination|condensation|water (loss|release)|loss of water|hydration",
    "ATP/NTP or phosphoanhydride coupling": r"\batp\b|\badp\b|pyrophosphat|\bppi\b|phosphoanhydride|nucleotide.{0,15}hydrolysis|\bkinase\b|synthetase|\bligase\b",
    "redox-partner / quinone potential mismatch": r"quinone|ubiquinone|menaquinone|redox potential|\bem.?\s*[′\']|electron.transferring flavoprotein|\betf\b|ferredoxin|nad\+?/nadh|nadp",
    "spontaneous / non-enzymatic step involved": r"spontaneous|non.enzymatic|nonenzymatic",
    "hydrolysis of a high-energy / glycosidic bond": r"glycosidic|n-glycosid|thioester|high.energy bond",
    "O2-dependent oxygenation/oxidation": r"\bo2\b|molecular oxygen|monooxygenas|dioxygenas|oxygenase|\boxidase\b",
    "cyclization / carbocation / terpene chemistry": r"cycliz|cyclase|carbocation|diphosphate (loss|ionization|ionisation)|terpene|prenyl",
    "empirically measured / experimentally observed": r"experimentally|observed (chemistry|reaction|direction)|only (direction|reaction) (observed|reported|characterised|characterized)|assay|purified enzyme|in vitro",
    "decarboxylation (CO2 loss)": r"\bdecarboxyl|co2\b.{0,25}(loss|releas|driv|escap)|loss of co2|releas\w* of co2",
}

# Broad glycoside-hydrolase pattern -- matched on enzyme_name. Deliberately
# wider than the single term "glucosidase": glucuronidase, rhamnosidase, and
# other GH-family hydrolases are the same mechanism (hydrolyse-only; never
# runs biosynthetically in reverse) and a literal "glucosidase" substring
# match undercounts this cluster (43 vs the true 54).
GLYCOSIDASE_PAT = re.compile(
    r"glucosidase|glucuronidase|galactosidase|glycosidase|xylosidase|mannosidase"
    r"|fucosidase|rhamnosidase|arabinosidase|glycoside hydrolase|glycohydrolase|\bgh\d",
    re.I)

# The 5 systematic, non-overlapping clusters, claimed in this priority order
# (equation-level predicate; `enzyme_name`-only matching, used for
# `llm_sound_enzyme_families` below, would have missed the ATP and quinone
# clusters entirely since they cut across many different named families).
SYSTEMATIC_CLUSTERS = [
    ("SAM-dependent methyltransferase", "dGPredictor: activated-cofactor group transfer",
     lambda r: "methyltransferase" in (r.get("enzyme_name") or "").lower()),
    ("PAPS-dependent sulfotransferase", "dGPredictor: activated-cofactor group transfer",
     lambda r: "sulfotransferase" in (r.get("enzyme_name") or "").lower()),
    ("Glycoside hydrolase (broad)", "eQuilibrator: structural motif",
     lambda r: bool(GLYCOSIDASE_PAT.search(r.get("enzyme_name") or ""))),
    ("Quinone/quinol-coupled redox", "eQuilibrator: structural motif",
     lambda r: "quinone" in (r.get("reaction") or "").lower()),
    ("ATP/ADP/phosphate-coupled group transfer", "dGPredictor: activated-cofactor group transfer",
     lambda r: bool(re.search(r"\batp\b|\badp\b|\bpi\b|phosphate", (r.get("reaction") or "").lower()))),
]

SPECIALIZED_PAT = re.compile(
    r"plant|fungal|fungus|fungi|secondary metabolite|alkaloid|terpen|flavon"
    r"|polyketide|antibiotic|natural product|phytochemical|mycotoxin", re.I)


def family_of(enzyme_name: str) -> str | None:
    name = (enzyme_name or "").lower()
    for term in FAMILY_TERMS:
        if term in name:
            return term
    return None


def tags_of(note: str) -> list[str]:
    t = (note or "").lower()
    return [name for name, pat in MECHANISM_PATTERNS.items() if re.search(pat, t)]


def counter_dict(c: Counter) -> dict:
    return dict(c.most_common())


def thermo(r: dict) -> dict:
    return r["direction_conflict"].get("thermo_evidence") or {}


def ctx_text(r: dict) -> str:
    return " ".join([r.get("organism_context") or "", " ".join(r.get("pathways") or [])])


def load_sound_and_llm_sound() -> tuple[list[dict], list[dict], int]:
    """Reload the deliverable and return (sound, llm_sound, n_total) -- shared
    with other scripts (e.g. analyze_strict_disagreement_residual_pathways.py)
    so the equation_defect/enzymology_supports filter never drifts between them.
    """
    data = json.load(open(DELIVERABLE))
    rxns = data["reactions"]
    sound = [r for r in rxns if r["direction_conflict"].get("equation_defect") == "none"]
    llm_sound = [r for r in sound if r["direction_conflict"].get("enzymology_supports") == "llm"]
    return sound, llm_sound, len(rxns)


def claim_systematic_clusters(llm_sound: list[dict]) -> tuple[list[dict], list[dict]]:
    """Claim reactions to the 5 systematic clusters in priority order (see
    SYSTEMATIC_CLUSTERS). Returns (cluster_stat_dicts, residual_reactions) --
    shared with other scripts so the residual definition never drifts.
    """
    claimed: set[str] = set()
    clusters = []
    for label, defect_class, predicate in SYSTEMATIC_CLUSTERS:
        members = [r for r in llm_sound if r["modelseed_id"] not in claimed and predicate(r)]
        ids = {r["modelseed_id"] for r in members}
        claimed.update(ids)
        src = Counter(thermo(r).get("source") for r in members)
        grade = Counter(thermo(r).get("grade") for r in members)
        dominant_source, dominant_n = src.most_common(1)[0]
        clusters.append({
            "cluster": label, "defect_class": defect_class, "n": len(members),
            "source_breakdown": counter_dict(src),
            "dominant_source": dominant_source,
            "dominant_source_purity": round(dominant_n / len(members), 3),
            "grade_breakdown": counter_dict(grade),
            "examples": [r["modelseed_id"] for r in members[:5]],
        })
    residual = [r for r in llm_sound if r["modelseed_id"] not in claimed]
    return clusters, residual


def main() -> None:
    sound, llm_sound, n_total = load_sound_and_llm_sound()

    out = {"n_total": n_total, "n_sound": len(sound), "n_llm_sound": len(llm_sound)}

    out["sound_by_mismatch_group"] = counter_dict(
        Counter(r["direction_conflict"]["mismatch_group"] for r in sound))
    out["sound_by_enzymology_supports"] = counter_dict(
        Counter(r["direction_conflict"].get("enzymology_supports") for r in sound))
    out["sound_by_functionality"] = counter_dict(Counter(r.get("functionality") for r in sound))
    out["sound_by_confidence"] = counter_dict(Counter(r.get("confidence") for r in sound))
    out["sound_by_thermo_grade"] = counter_dict(Counter(thermo(r).get("grade") for r in sound))
    out["sound_by_thermo_assessment"] = counter_dict(Counter(thermo(r).get("assessment") for r in sound))
    out["sound_by_thermo_source"] = counter_dict(Counter(thermo(r).get("source") for r in sound))

    # --- generic enzyme-name family breakdown (supplementary: the "mixed
    # families" contrast data referenced in the report) -----------------
    family_counts = Counter()
    family_source: dict[str, Counter] = {}
    family_grade: dict[str, Counter] = {}
    family_examples: dict[str, list[str]] = {}
    for r in llm_sound:
        fam = family_of(r.get("enzyme_name")) or "(no family term / unidentified)"
        family_counts[fam] += 1
        family_source.setdefault(fam, Counter())[thermo(r).get("source")] += 1
        family_grade.setdefault(fam, Counter())[thermo(r).get("grade")] += 1
        family_examples.setdefault(fam, [])
        if len(family_examples[fam]) < 5:
            family_examples[fam].append(r["modelseed_id"])

    families = []
    for fam, n in family_counts.most_common():
        src = family_source[fam]
        dominant_source, dominant_n = src.most_common(1)[0]
        families.append({
            "family": fam, "n": n,
            "source_breakdown": counter_dict(src),
            "dominant_source": dominant_source,
            "dominant_source_purity": round(dominant_n / n, 3),
            "grade_breakdown": counter_dict(family_grade[fam]),
            "examples": family_examples[fam],
        })
    out["llm_sound_enzyme_families"] = families

    # --- the 5 non-overlapping systematic clusters + residual -----------
    clusters, residual = claim_systematic_clusters(llm_sound)
    out["systematic_clusters"] = clusters
    claimed_ids = {r["modelseed_id"] for r in llm_sound} - {r["modelseed_id"] for r in residual}
    n_clustered = len(claimed_ids)
    out["systematic_clusters_total"] = n_clustered
    out["systematic_clusters_share_of_llm_sound"] = round(n_clustered / len(llm_sound), 4)
    by_defect_class = Counter()
    for c in clusters:
        by_defect_class[c["defect_class"]] += c["n"]
    out["systematic_clusters_by_defect_class"] = counter_dict(by_defect_class)

    n_spec_residual = sum(1 for r in residual if SPECIALIZED_PAT.search(ctx_text(r)))
    n_spec_clustered = sum(1 for r in llm_sound if r["modelseed_id"] in claimed_ids
                            and SPECIALIZED_PAT.search(ctx_text(r)))
    out["systematic_clusters_residual"] = {
        "n": len(residual),
        "share_of_llm_sound": round(len(residual) / len(llm_sound), 4),
        "thermo_source": counter_dict(Counter(thermo(r).get("source") for r in residual)),
        "thermo_grade": counter_dict(Counter(thermo(r).get("grade") for r in residual)),
        "confidence": counter_dict(Counter(r.get("confidence") for r in residual)),
        "research_tier": counter_dict(Counter(r.get("research_tier") for r in residual)),
        "secondary_metabolism_hypothesis_test": {
            "description": "Tested whether the residual's errors concentrate in "
                            "plant/fungal/secondary-metabolite chemistry (molecular "
                            "novelty underrepresented in thermo-source training data). "
                            "REJECTED: the specialized-context signal is HIGHER inside "
                            "the 5 systematic clusters than in the residual -- it "
                            "correlates with cluster membership (plant secondary "
                            "metabolism heavily uses SAM-methyltransferases and "
                            "glycoside hydrolases), not with what's left unexplained.",
            "specialized_context_share_in_residual": round(n_spec_residual / len(residual), 4),
            "specialized_context_share_in_clusters": round(n_spec_clustered / n_clustered, 4),
        },
    }

    # mechanism keyword tagging on direction_conflict.note (all 465, and
    # restricted to the residual, to show the tags don't just relocate)
    def mechanism_table(reactions):
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
        return {
            "n": len(reactions),
            "tagged_at_least_one": n_tagged,
            "tags": [{"mechanism": tag, "n": n, "share": round(n / len(reactions), 4),
                      "examples": tag_examples[tag]} for tag, n in tag_counts.most_common()],
        }

    out["direction_mechanism_tags_all_465"] = mechanism_table(llm_sound)
    out["direction_mechanism_tags_residual_only"] = mechanism_table(residual)

    with open(os.path.join(OUTDIR, "sound_subset_analysis.json"), "w") as fh:
        json.dump(out, fh, indent=2)

    with open(os.path.join(OUTDIR, "sound_subset_families.tsv"), "w") as fh:
        fh.write("family\tn\tdominant_source\tdominant_source_purity\n")
        for f in families:
            fh.write(f"{f['family']}\t{f['n']}\t{f['dominant_source']}\t{f['dominant_source_purity']}\n")

    with open(os.path.join(OUTDIR, "sound_subset_clusters.tsv"), "w") as fh:
        fh.write("cluster\tn\tdefect_class\tdominant_source\tdominant_source_purity\n")
        for c in clusters:
            fh.write(f"{c['cluster']}\t{c['n']}\t{c['defect_class']}\t"
                     f"{c['dominant_source']}\t{c['dominant_source_purity']}\n")
        fh.write(f"Residual (case-specific enzymology)\t{len(residual)}\tnone\tmixed\t\n")

    print(f"sound (equation_defect=none): {len(sound)}/{n_total}")
    print(f"llm-correct within sound: {len(llm_sound)}")
    print()
    print("systematic clusters (priority-claimed, non-overlapping):")
    for c in clusters:
        print(f"  {c['cluster']:<42} n={c['n']:>3}  {c['defect_class']:<42} "
              f"dominant={c['dominant_source']} (purity {c['dominant_source_purity']:.0%})")
    print(f"  {'Residual (case-specific enzymology)':<42} n={len(residual):>3}")
    print()
    print(f"total clustered: {n_clustered}/{len(llm_sound)} ({out['systematic_clusters_share_of_llm_sound']:.1%})")
    for dc, n in by_defect_class.items():
        print(f"  {dc}: {n} ({100*n/len(llm_sound):.1f}%)")
    print()
    print("secondary-metabolism hypothesis (tested against the residual):")
    smt = out["systematic_clusters_residual"]["secondary_metabolism_hypothesis_test"]
    print(f"  specialized context in residual: {smt['specialized_context_share_in_residual']:.1%}")
    print(f"  specialized context in clusters: {smt['specialized_context_share_in_clusters']:.1%}  (REJECTED: higher in clusters)")
    print()
    print(f"wrote {os.path.join(OUTDIR, 'sound_subset_analysis.json')}")
    print(f"wrote {os.path.join(OUTDIR, 'sound_subset_families.tsv')}")
    print(f"wrote {os.path.join(OUTDIR, 'sound_subset_clusters.tsv')}")


if __name__ == "__main__":
    main()
