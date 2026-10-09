#!/usr/bin/env python3
"""Which biological pathways/networks are the 248 "residual" strict
disagreements from? (the reactions NOT explained by either of the two
systematic thermo-source defect classes in
analyze_strict_disagreement_sound_subset.py -- see that script for
`load_sound_and_llm_sound` / `claim_systematic_clusters`, reused here so the
residual definition never drifts between the two analyses.)

Each reaction carries a `pathways` list (KEGG `rn#####` + MetaCyc PGDB
pathway/class strings, e.g. "KEGG: rn00625 (Chloroalkane and chloroalkene
degradation)"); 155/248 (62.5%) have at least one. Parenthetical names are
pulled out and bucketed into broad categories by keyword match (priority
order, first match wins, so a reaction counts once). For the 93 reactions
with no `pathways` entry, the same keyword set is applied to `enzyme_name` +
`rationale` + `organism_context` + `direction_conflict.note` as a fallback --
flagged separately (`source: "text-inferred"` vs `"pathway-annotated"`) since
it is a weaker signal than a curated KEGG/MetaCyc assignment.

Outputs (results/strict_disagreement_enzymes/):
    residual_pathways.json   category counts, per-category source/family
                              breakdown, examples; full per-reaction category
                              assignment
    residual_pathways.tsv    category, n, n_pathway_annotated, n_text_inferred
                              (figure source data)
"""
from __future__ import annotations

import json
import os
import re
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.dirname(HERE)
OUTDIR = os.path.join(ROOT, "results", "strict_disagreement_enzymes")

from analyze_strict_disagreement_sound_subset import load_sound_and_llm_sound, claim_systematic_clusters

# Priority-ordered: first category whose pattern matches wins, so each
# reaction is assigned to exactly one. Xenobiotic/pollutant degradation is
# checked first because its vocabulary ("degradation", "aromatic compound")
# would otherwise be swallowed by the generic MetaCyc superclass
# "Degradation/Utilization/Assimilation" without naming what's degraded.
CATEGORIES = [
    ("Xenobiotic / pollutant biodegradation",
     r"xenobiotic|biodegradation|detoxification|drug metabolism|dehalogenat|mercury"
     r"|chloro|dichlor|trichlor|fluoro\w*benzoate|atrazine|triazine|thiocyanate"
     r"|carbon disulfide|nicotine degradation"
     r"|\bddt\b|naphthalene|toluene|nitroaromatic|polycyclic aromatic|pesticide"
     r"|aromatic compound degradation|degradation of aromatic|chloroaromatic"
     r"|caprolactam|benzoate degradation"),
    ("Secondary metabolite biosynthesis (plant/microbial)",
     r"secondary metabolite|phenylpropanoid|flavonoid|flavonol|isoflavonoid|anthocyanin"
     r"|terpenoid|carotenoid|limonene|pinene|alkaloid|phytoalexin|antibiotic biosynthesis"
     r"|nitrogen-containing secondary|sulfur-containing secondary|polyketide|stilbenoid"
     r"|lignan|glucosinolate"),
    ("Steroid / sterol / bile-acid metabolism",
     r"steroid|cholesterol|bile acid|sterol"),
    ("Fatty acid / lipid metabolism",
     r"fatty acid|lipid degradation|lipid metabolism|glycerolipid|glycerophospholipid"
     r"|sphingolipid|beta-oxidation"),
    ("Amino acid metabolism",
     r"amino acid (metabolism|biosynthesis|degradation)|cysteine and methionine"
     r"|taurine|glycine, serine|lysine (biosynthesis|degradation)|arginine and proline"
     r"|tyrosine metabolism|tryptophan metabolism|phenylalanine metabolism"
     r"|histidine metabolism|valine, leucine|amine and polyamine|choline degradation"),
    ("Nucleotide metabolism",
     r"purine|pyrimidine|nucleotide|nucleoside|guanosine|adenosine|cytidine"),
    ("Cofactor / vitamin / coenzyme biosynthesis",
     r"cofactor|vitamin|coenzyme|nicotinate and nicotinamide|folate|biotin metabolism"
     r"|riboflavin|thiamine metabolism|pantothenate|porphyrin|heme biosynthesis"
     r"|molybdopterin|ubiquinone biosynthesis"),
    ("Carbohydrate / central carbon metabolism",
     r"carbohydrate metabolism|glycolysis|pentose phosphate|pyruvate metabolism"
     r"|butanoate metabolism|propanoate metabolism|glyoxylate|citrate cycle"),
    ("Energy / precursor-metabolite generation",
     r"generation of precursor metabolite|energy metabolism|methane metabolism"
     r"|oxidative phosphorylation|nitrogen metabolism|sulfur metabolism"),
    ("Glycan / cell-envelope metabolism",
     r"glycan|peptidoglycan|murein|lipopolysaccharide|cell wall"),
]


def pathway_text(r: dict) -> str:
    return " ".join(r.get("pathways") or [])


def fallback_text(r: dict) -> str:
    return " ".join([
        r.get("enzyme_name") or "", r.get("rationale") or "",
        r.get("organism_context") or "", r["direction_conflict"].get("note") or "",
    ])


def categorize(text: str) -> str | None:
    for label, pat in CATEGORIES:
        if re.search(pat, text, re.I):
            return label
    return None


def thermo(r: dict) -> dict:
    return r["direction_conflict"].get("thermo_evidence") or {}


def counter_dict(c: Counter) -> dict:
    return dict(c.most_common())


def main() -> None:
    _, llm_sound, _ = load_sound_and_llm_sound()
    _, residual = claim_systematic_clusters(llm_sound)

    assignments = []
    for r in residual:
        ptext = pathway_text(r)
        cat = categorize(ptext) if ptext else None
        provenance = "pathway-annotated" if cat else None
        if cat is None:
            cat = categorize(fallback_text(r))
            provenance = "text-inferred" if cat else None
        if cat is None:
            cat = "No pathway identifiable"
            provenance = "none"
        assignments.append({"modelseed_id": r["modelseed_id"], "category": cat, "provenance": provenance})

    by_id = {a["modelseed_id"]: a for a in assignments}
    cat_counts = Counter(a["category"] for a in assignments)
    cat_provenance = {}
    cat_examples: dict[str, list[str]] = {}
    cat_source: dict[str, Counter] = {}
    for r in residual:
        a = by_id[r["modelseed_id"]]
        cat = a["category"]
        cat_provenance.setdefault(cat, Counter())[a["provenance"]] += 1
        cat_source.setdefault(cat, Counter())[thermo(r).get("source")] += 1
        cat_examples.setdefault(cat, [])
        if len(cat_examples[cat]) < 6:
            cat_examples[cat].append(r["modelseed_id"])

    categories_out = []
    for cat, n in cat_counts.most_common():
        categories_out.append({
            "category": cat, "n": n,
            "share_of_residual": round(n / len(residual), 4),
            "provenance": counter_dict(cat_provenance[cat]),
            "thermo_source": counter_dict(cat_source[cat]),
            "examples": cat_examples[cat],
        })

    out = {
        "n_residual": len(residual),
        "n_pathway_annotated": sum(1 for a in assignments if a["provenance"] == "pathway-annotated"),
        "n_text_inferred": sum(1 for a in assignments if a["provenance"] == "text-inferred"),
        "n_no_pathway_identifiable": sum(1 for a in assignments if a["provenance"] == "none"),
        "categories": categories_out,
        "per_reaction": assignments,
    }

    with open(os.path.join(OUTDIR, "residual_pathways.json"), "w") as fh:
        json.dump(out, fh, indent=2)

    with open(os.path.join(OUTDIR, "residual_pathways.tsv"), "w") as fh:
        fh.write("category\tn\tn_pathway_annotated\tn_text_inferred\n")
        for c in categories_out:
            fh.write(f"{c['category']}\t{c['n']}\t{c['provenance'].get('pathway-annotated', 0)}\t"
                     f"{c['provenance'].get('text-inferred', 0)}\n")

    print(f"residual: {len(residual)}")
    print(f"pathway-annotated: {out['n_pathway_annotated']}, text-inferred: {out['n_text_inferred']}, "
          f"no pathway identifiable: {out['n_no_pathway_identifiable']}")
    print()
    for c in categories_out:
        print(f"  {c['category']:<50} n={c['n']:>3} ({c['share_of_residual']:.1%})  "
              f"source={c['thermo_source']}")
    print()
    print(f"wrote {os.path.join(OUTDIR, 'residual_pathways.json')}")
    print(f"wrote {os.path.join(OUTDIR, 'residual_pathways.tsv')}")


if __name__ == "__main__":
    main()
