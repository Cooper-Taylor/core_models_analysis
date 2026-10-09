#!/usr/bin/env python3
"""Database-wide analysis of ModelSEED's reaction `status` field (the mass/
charge balance check every reaction carries) -- scope is ALL 48,384
non-obsolete reactions, not any disagreement subset. Answers three questions:

  1. What flags exist, and how many reactions carry each (+ combinations)?
  2. What features/patterns characterise the reactions in each flag group?
  3. How tractable is each group to actually fix?

`status` is computed by BiochemPy.Reactions.balanceReaction() (confirmed by
reading Libs/Python/BiochemPy/Reactions.py in the full ModelSEEDDatabase
checkout, lines 316-442 -- not present in the data-only snapshot). It is a
`|`-joined string of codes, each optionally carrying a `:`-delimited payload;
the exact parsing convention (`fields = status.split('|')`; a field's type is
whatever precedes its first `:`) is taken from
Scripts/Validation/Validate_Reactions.py so this script tokenises status
strings exactly the way ModelSEED's own validator does.

Outputs (results/reaction_balance_flags/):
    flag_counts.tsv             (A) every flag token + combination, with counts
    cpdformerror_compounds.tsv  which missing-formula compounds drive CPDFORMERROR,
                                 ranked by how many reactions they appear in
    ci_only_analysis.json       charge-imbalance magnitude distribution,
                                 combinatorial-expansion-template overlap
    mi_only_analysis.json       mass-residual-vs-known-compound-formula matches
    both_analysis.json          EC-annotation and imbalance-multiplicity profile
                                 of the combined CI+MI group
    group_tractability.tsv      figure source data: identified-pattern vs
                                 residual share per group (the "ease of fixing" chart)

Data: ModelSEED dev @ 078a395f (/scratch/ctaylor/tmp/devsnap_078a395f),
confirmed still == origin/dev HEAD as of 2026-10-09 (unchanged since
2026-10-06). Scope: all 48,384 non-obsolete reactions, not a disagreement
subset -- a different population from reports/strictDisagreementEnzymes/.
"""
from __future__ import annotations

import glob
import json
import os
import re
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MSDB_SHA = "078a395f"
MSDB_ROOT = os.environ.get("MSDB_ROOT", "/scratch/ctaylor/tmp/devsnap_078a395f")
BIOCHEM = os.path.join(MSDB_ROOT, "Biochemistry")
OUTDIR = os.path.join(ROOT, "results", "reaction_balance_flags")

EXP_PAT = re.compile(r'\.(meta|chlamy|maize|plant|fungal|gram|eco)?exp\.', re.I)


def status_tokens(status: str | None) -> list[str]:
    """Exact parsing convention from Scripts/Validation/Validate_Reactions.py."""
    status = status or ""
    return [f[:f.find(':')] if ':' in f else f for f in status.split('|')]


def parse_formula(f: str | None) -> dict[str, float] | None:
    if not f:
        return None
    d: dict[str, float] = {}
    for elem, cnt in re.findall(r'([A-Z][a-z]?)(-?\d*\.?\d*)', f):
        if not elem:
            continue
        cnt = cnt if cnt not in ("", "-") else "1"
        try:
            d[elem] = d.get(elem, 0) + float(cnt)
        except ValueError:
            return None
    return d


def load_reactions() -> list[dict]:
    rxns = []
    for path in sorted(glob.glob(os.path.join(BIOCHEM, "reaction_*.json"))):
        for rxn in json.load(open(path)):
            if not rxn.get("is_obsolete"):
                rxns.append(rxn)
    return rxns


def load_compounds() -> tuple[dict[str, str], dict[tuple, list[str]]]:
    cpd_name: dict[str, str] = {}
    formula_index: dict[tuple, list[str]] = {}
    for path in sorted(glob.glob(os.path.join(BIOCHEM, "compound_*.json"))):
        for cpd in json.load(open(path)):
            cpd_name[cpd["id"]] = cpd.get("name")
            f = cpd.get("formula")
            if f:
                parsed = parse_formula(f)
                if parsed:
                    sig = tuple(sorted((e, round(c, 2)) for e, c in parsed.items() if e != "R"))
                    formula_index.setdefault(sig, []).append(cpd["id"])
    return cpd_name, formula_index


def group_by_tokens(rxns: list[dict]) -> dict[str, list[dict]]:
    groups: dict[str, list[dict]] = {"OK": [], "CPDFORMERROR": [], "CI_only": [],
                                      "MI_only": [], "both": [], "CK_combo": [], "EMPTY": []}
    for r in rxns:
        toks = set(status_tokens(r.get("status")))
        if toks == {"OK"}:
            groups["OK"].append(r)
        elif toks == {"CPDFORMERROR"}:
            groups["CPDFORMERROR"].append(r)
        elif toks == {"CI"}:
            groups["CI_only"].append(r)
        elif toks == {"MI"}:
            groups["MI_only"].append(r)
        elif toks == {"CI", "MI"}:
            groups["both"].append(r)
        elif "CK" in toks:
            groups["CK_combo"].append(r)
        elif toks == {"EMPTY"}:
            groups["EMPTY"].append(r)
    return groups


def analyze_cpdformerror(rxns: list[dict], cpd_name: dict[str, str]) -> dict:
    bad_cpd_hits: Counter = Counter()
    # a compound is "bad" if it has no formula; cross-check via name lookup absence of formula
    # is already implicit: reactions in this group only arise when >=1 reagent lacks a formula,
    # so every cpd id referenced that has no name->formula is a candidate; re-derive directly:
    cpd_has_formula = set()
    for path in sorted(glob.glob(os.path.join(BIOCHEM, "compound_*.json"))):
        for cpd in json.load(open(path)):
            if cpd.get("formula"):
                cpd_has_formula.add(cpd["id"])
    for r in rxns:
        cids = set(re.findall(r'(cpd\d{5})', r.get("equation", "") or ""))
        for cid in cids:
            if cid not in cpd_has_formula:
                bad_cpd_hits[cid] += 1
    top10_share = sum(n for _, n in bad_cpd_hits.most_common(10))
    return {
        "n": len(rxns),
        "n_distinct_bad_compounds": len(bad_cpd_hits),
        "top_compounds": [{"id": cid, "name": cpd_name.get(cid), "n_reactions": n}
                           for cid, n in bad_cpd_hits.most_common(25)],
        "top10_reaction_hits": top10_share,
    }, bad_cpd_hits


def analyze_ci_only(rxns: list[dict]) -> dict:
    mag_counts: Counter = Counter()
    for r in rxns:
        m = re.search(r'CI:(-?[\d.]+)', r["status"])
        if m:
            mag_counts[abs(float(m.group(1)))] += 1

    exp_rxns = [r for r in rxns if EXP_PAT.search(r.get("abbreviation", "") or "")]
    n_h_in_exp = sum(1 for r in exp_rxns if "cpd00067" in re.findall(r'(cpd\d{5})', r.get("equation", "") or ""))
    n_transport = sum(1 for r in rxns if r.get("is_transport"))

    return {
        "n": len(rxns),
        "magnitude_distribution": {str(k): v for k, v in sorted(mag_counts.items())},
        "n_transport": n_transport,
        "n_combinatorial_expansion": len(exp_rxns),
        "n_expansion_involving_H+": n_h_in_exp,
        "expansion_examples": [r["id"] for r in exp_rxns[:10]],
    }


def analyze_mi_only(rxns: list[dict], formula_index: dict, cpd_name: dict[str, str]) -> dict:
    match_counter: Counter = Counter()
    n_matched = 0
    for r in rxns:
        m = re.search(r'MI:([^|]+)', r["status"])
        if not m:
            continue
        resid: dict[str, float] = {}
        ok = True
        for p in m.group(1).split('/'):
            if ':' not in p:
                ok = False
                break
            elem, val = p.split(':')
            if elem == 'R':
                ok = False
                break
            try:
                resid[elem] = round(float(val), 2)
            except ValueError:
                ok = False
                break
        if not ok or not resid:
            continue
        for sign in (1, -1):
            sig = tuple(sorted((e, sign * v) for e, v in resid.items()))
            if sig in formula_index:
                match_counter[formula_index[sig][0]] += 1
                n_matched += 1
                break

    return {
        "n": len(rxns),
        "n_residual_matches_known_compound": n_matched,
        "top_matches": [{"id": cid, "name": cpd_name.get(cid), "n_reactions": n}
                         for cid, n in match_counter.most_common(20)],
    }


def analyze_both(rxns: list[dict]) -> dict:
    n_fractional = sum(1 for r in rxns
                        if re.search(r'\((\d*\.\d+)\)', r.get("equation", "") or ""))
    ec_counter: Counter = Counter()
    for r in rxns:
        ecs = sorted(set(e.split('.')[0] for e in (r.get('ec_numbers') or []) if e and e[0].isdigit()))
        if not ecs:
            ec_counter['no EC'] += 1
        else:
            for c in ecs:
                ec_counter[c] += 1
    n_no_ec = ec_counter.get('no EC', 0)

    n_elems_dist: Counter = Counter()
    n_multi = 0
    for r in rxns:
        m = re.search(r'MI:([^|]+)', r["status"])
        if m:
            n = len(m.group(1).split('/'))
            n_elems_dist[n] += 1
            if n >= 3:
                n_multi += 1

    n_transport = sum(1 for r in rxns if r.get("is_transport"))

    return {
        "n": len(rxns),
        "n_no_ec": n_no_ec,
        "n_fractional_coefficient": n_fractional,
        "n_transport": n_transport,
        "n_3plus_elements_imbalanced": n_multi,
        "elements_imbalanced_distribution": {str(k): v for k, v in sorted(n_elems_dist.items())},
    }


def main() -> None:
    os.makedirs(OUTDIR, exist_ok=True)
    rxns = load_reactions()
    cpd_name, formula_index = load_compounds()
    groups = group_by_tokens(rxns)

    n_total = len(rxns)
    print(f"non-obsolete reactions: {n_total}")

    # (A) flag counts, including sub-combinations
    combo_counter: Counter = Counter()
    for r in rxns:
        combo_counter[tuple(sorted(set(status_tokens(r.get("status")))))] += 1
    with open(os.path.join(OUTDIR, "flag_counts.tsv"), "w") as fh:
        fh.write("combination\tn\tshare\n")
        for combo, n in combo_counter.most_common():
            fh.write(f"{'+'.join(combo)}\t{n}\t{n/n_total:.4f}\n")

    cpdform_result, bad_cpd_hits = analyze_cpdformerror(groups["CPDFORMERROR"], cpd_name)
    with open(os.path.join(OUTDIR, "cpdformerror_compounds.tsv"), "w") as fh:
        fh.write("compound_id\tname\tn_reactions\n")
        for cid, n in bad_cpd_hits.most_common(50):
            fh.write(f"{cid}\t{cpd_name.get(cid)}\t{n}\n")
    with open(os.path.join(OUTDIR, "cpdformerror_analysis.json"), "w") as fh:
        json.dump(cpdform_result, fh, indent=2)

    ci_result = analyze_ci_only(groups["CI_only"])
    with open(os.path.join(OUTDIR, "ci_only_analysis.json"), "w") as fh:
        json.dump(ci_result, fh, indent=2)

    mi_result = analyze_mi_only(groups["MI_only"], formula_index, cpd_name)
    with open(os.path.join(OUTDIR, "mi_only_analysis.json"), "w") as fh:
        json.dump(mi_result, fh, indent=2)

    both_result = analyze_both(groups["both"])
    with open(os.path.join(OUTDIR, "both_analysis.json"), "w") as fh:
        json.dump(both_result, fh, indent=2)

    # group_tractability.tsv: identified-pattern vs residual, for the summary chart
    rows = [
        ("CPDFORMERROR", len(groups["CPDFORMERROR"]), cpdform_result["top10_reaction_hits"],
         "top-10 generic cofactor placeholders"),
        ("CI only", len(groups["CI_only"]), ci_result["n_combinatorial_expansion"],
         "combinatorial-expansion template bug"),
        ("MI only", len(groups["MI_only"]), mi_result["n_residual_matches_known_compound"],
         "residual matches a known compound (dropped product)"),
        ("CI+MI both", len(groups["both"]), len(groups["both"]) - both_result["n_3plus_elements_imbalanced"],
         "single/double-element imbalance (simpler cases)"),
    ]
    with open(os.path.join(OUTDIR, "group_tractability.tsv"), "w") as fh:
        fh.write("group\tn_total\tn_identified_pattern\tn_residual\tpattern_label\n")
        for label, n, n_identified, pattern_label in rows:
            fh.write(f"{label}\t{n}\t{n_identified}\t{n - n_identified}\t{pattern_label}\n")

    print()
    print("=== (A) flag combinations ===")
    for combo, n in combo_counter.most_common():
        print(f"  {'+'.join(combo):<20} {n:>6}  {100*n/n_total:.2f}%")
    print()
    print(f"CPDFORMERROR: {cpdform_result['n']} reactions, {cpdform_result['n_distinct_bad_compounds']} "
          f"distinct bad-formula compounds, top 10 cover {cpdform_result['top10_reaction_hits']} reaction-hits")
    print(f"CI-only: {ci_result['n']}, combinatorial-expansion-driven: {ci_result['n_combinatorial_expansion']} "
          f"({100*ci_result['n_combinatorial_expansion']/ci_result['n']:.1f}%)")
    print(f"MI-only: {mi_result['n']}, residual matches known compound: "
          f"{mi_result['n_residual_matches_known_compound']} "
          f"({100*mi_result['n_residual_matches_known_compound']/mi_result['n']:.1f}%)")
    print(f"Both: {both_result['n']}, no EC: {both_result['n_no_ec']} "
          f"({100*both_result['n_no_ec']/both_result['n']:.1f}%), "
          f"3+ elements imbalanced: {both_result['n_3plus_elements_imbalanced']} "
          f"({100*both_result['n_3plus_elements_imbalanced']/both_result['n']:.1f}%)")
    print()
    print(f"wrote outputs to {OUTDIR}/")


if __name__ == "__main__":
    main()
