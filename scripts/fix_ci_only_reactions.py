#!/usr/bin/env python3
"""Of the 2,824 charge-imbalance-only (`CI`, no `MI`) reactions characterised
in reports/reactionBalanceFlags/REACTION_BALANCE_FLAGS.md, find which ones
have a mechanically verifiable, zero-blast-radius fix, and compute it.

Background: that report guessed (but explicitly did NOT verify) that ~778
"combinatorial expansion" reactions shared one templating bug ("a proton
coefficient that wasn't recalculated correctly"). Digging into actual
equations shows that guess doesn't hold up -- e.g. rxn12950 is just
`NADH => NAD` with its real substrate pair missing entirely, not a proton
slip. What *does* hold up, verified independently below, is a different and
more precise mechanism:

For every CI-only reaction, ModelSEED's own balanceReaction() has already
confirmed mass is fully balanced (that's what "CI-only, no MI" means) -- so
whatever's wrong is a pure charge-accounting error, not a missing atom. Many
of these reactions carry a redox-cofactor couple (NAD/NADH, NADP/NADPH,
oxidized/reduced ferredoxin, cytochrome c2+/c3+, ...): two compounds with an
*identical* chemical formula but different formal charge. Rescaling such a
pair's shared coefficient by a single integer multiplier changes the
reaction's net charge but can never change its net mass (same formula on
both sides, any multiple of it still cancels). So for each CI-only reaction,
we test every reactant/product pair sharing an exact formula match: does
some integer rescaling of that pair's coefficients exist that *exactly*
zeroes the charge? If there's exactly one such pair-and-multiplier, that's a
verified, local, minimal fix -- it touches only that one reaction's existing
terms, never a compound definition, so it cannot regress any other reaction.

Three outcomes per reaction:
  - solved (n_candidates == 1): unique fix, verified end-to-end by rebuilding
    the full term list and recomputing mass+charge from scratch.
  - ambiguous (n_candidates > 1): more than one redox pair in the equation
    could independently explain the imbalance -- which one is actually wrong
    needs the reaction's real mechanism, not arithmetic. Not auto-fixed.
  - unresolved (n_candidates == 0): no same-formula pair lever exists at all,
    or none gives a chemically sane (positive, <= MAXCOEF) integer solution.
    These need the same per-reaction literature review as every other
    imbalance cluster -- nothing here claims otherwise.

Data: ModelSEED dev @ 078a395f (/scratch/ctaylor/tmp/devsnap_078a395f),
reconfirmed == origin/dev HEAD on 2026-10-12 (unchanged since 2026-10-06).

Outputs (results/ci_only_fixes/):
    solved.tsv            one row per uniquely-fixed reaction: id, old/new
                           equation, the pair and coefficient change
    ambiguous.tsv          the handful with >1 viable pair
    pair_family_counts.tsv how the 389 solved reactions break down by
                            cofactor-pair family (figure source)
    outcome_counts.tsv     solved / ambiguous / unresolved counts (figure
                            source)
    summary.json            headline numbers for the report
"""
from __future__ import annotations

import csv
import glob
import json
import os
import re
from collections import Counter
from pathlib import Path

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
MSDB_SHA = "078a395f"
MSDB_ROOT = os.environ.get("MSDB_ROOT", "/scratch/ctaylor/tmp/devsnap_078a395f")
BIOCHEM = os.path.join(MSDB_ROOT, "Biochemistry")
OUTDIR = ANALYSIS_DIR / "results" / "ci_only_fixes"

TERM_RE = re.compile(r'\(([\-\d.]+)\)\s*([A-Za-z0-9_]+)\[(\d+)\]')
OP_RE = re.compile(r'(<=>|<=|=>|<=?)')
MAXCOEF = 12  # reject "solutions" needing absurdly large stoichiometry
TOL = 1e-6


def status_tokens(status: str | None) -> list[str]:
    status = status or ""
    return [f[:f.find(':')] if ':' in f else f for f in status.split('|')]


def parse_formula(f: str | None) -> dict[str, int]:
    d: dict[str, int] = {}
    if not f:
        return d
    for elem, cnt in re.findall(r'([A-Z][a-z]?)(\d*)', f):
        if not elem:
            continue
        d[elem] = d.get(elem, 0) + (int(cnt) if cnt else 1)
    return d


def load_compounds() -> dict[str, dict]:
    cpds: dict[str, dict] = {}
    for path in sorted(glob.glob(os.path.join(BIOCHEM, "compound_*.json"))):
        for c in json.load(open(path)):
            cpds[c["id"]] = c
    return cpds


def load_reactions() -> list[dict]:
    rxns = []
    for path in sorted(glob.glob(os.path.join(BIOCHEM, "reaction_*.json"))):
        for r in json.load(open(path)):
            if not r.get("is_obsolete"):
                rxns.append(r)
    return rxns


def parse_equation(eq: str):
    parts = OP_RE.split(eq.replace(' ', ''))
    if len(parts) != 3:
        return None
    left, _op, right = parts

    def parse_side(s):
        terms: Counter = Counter()
        for coef, cid, _comp in TERM_RE.findall(s):
            terms[cid] += float(coef)
        return terms

    return parse_side(left), parse_side(right)


def net_balance(reactants: Counter, products: Counter, cpds: dict):
    net_elem: Counter = Counter()
    net_charge = 0.0
    bad = False
    for cid, coef in reactants.items():
        c = cpds.get(cid)
        if c is None or c.get("charge") is None or not c.get("formula"):
            bad = True
            continue
        net_charge -= coef * c["charge"]
        for e, n in parse_formula(c["formula"]).items():
            net_elem[e] -= coef * n
    for cid, coef in products.items():
        c = cpds.get(cid)
        if c is None or c.get("charge") is None or not c.get("formula"):
            bad = True
            continue
        net_charge += coef * c["charge"]
        for e, n in parse_formula(c["formula"]).items():
            net_elem[e] += coef * n
    return net_elem, net_charge, bad


def solve_pairs(reactants: Counter, products: Counter, cpds: dict, net_charge_total: float):
    """Same-formula reactant/product pairs, rescaled by one integer
    multiplier, that exactly zero net charge while leaving mass untouched
    (b' - a' is held equal to the pair's current b - a, which is exactly the
    condition under which the pair's own net-mass contribution is unchanged).
    """
    candidates = []
    for rcid, a in reactants.items():
        rc = cpds.get(rcid)
        if not rc or not rc.get("formula") or rc.get("charge") is None:
            continue
        rform = parse_formula(rc["formula"])
        for pcid, b in products.items():
            pc = cpds.get(pcid)
            if not pc or not pc.get("formula") or pc.get("charge") is None:
                continue
            if parse_formula(pc["formula"]) != rform:
                continue
            dq = pc["charge"] - rc["charge"]
            if abs(dq) < TOL:
                continue
            k = b - a
            old_contrib = -a * rc["charge"] + b * pc["charge"]
            req = old_contrib - net_charge_total
            a_prime = (req - k * pc["charge"]) / dq
            b_prime = a_prime + k
            a_r, b_r = round(a_prime), round(b_prime)
            if abs(a_prime - a_r) > TOL or abs(b_prime - b_r) > TOL:
                continue
            if a_r < 1 or b_r < 1 or a_r > MAXCOEF or b_r > MAXCOEF:
                continue  # reject zero/negative (= deleting the reagent) and absurd scale
            if a_r == a and b_r == b:
                continue
            candidates.append({
                "reactant": rcid, "reactant_name": rc.get("name"),
                "product": pcid, "product_name": pc.get("name"),
                "formula": rc["formula"],
                "charge_react": rc["charge"], "charge_prod": pc["charge"],
                "old_coef_react": a, "old_coef_prod": b,
                "new_coef_react": a_r, "new_coef_prod": b_r,
            })
    return candidates


def new_equation_str(eq: str, cand: dict) -> str:
    """Render the corrected equation by substituting the two coefficients."""
    def sub(match):
        coef, cid, comp = match.groups()
        if cid == cand["reactant"] and float(coef) == cand["old_coef_react"]:
            return f"({cand['new_coef_react']}) {cid}[{comp}]"
        if cid == cand["product"] and float(coef) == cand["old_coef_prod"]:
            return f"({cand['new_coef_prod']}) {cid}[{comp}]"
        return match.group(0)
    return TERM_RE.sub(sub, eq)


PAIR_FAMILY = [
    ("NAD / NADH", {"NAD", "NADH"}),
    ("NADP / NADPH", {"NADP", "NADPH"}),
    ("Ferredoxin", None),  # matched by substring below
    ("Cytochrome", None),
    ("Other electron carrier", None),
]


def pair_family(reactant_name: str, product_name: str) -> str:
    names = {reactant_name, product_name}
    if names == {"NAD", "NADH"}:
        return "NAD / NADH"
    if names == {"NADP", "NADPH"}:
        return "NADP / NADPH"
    joined = f"{reactant_name} {product_name}".lower()
    if "ferredoxin" in joined:
        return "Ferredoxin"
    if "cytochrome" in joined:
        return "Cytochrome"
    return "Other electron carrier"


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    cpds = load_compounds()
    rxns = load_reactions()
    ci_only = [r for r in rxns if set(status_tokens(r.get("status"))) == {"CI"}]
    print(f"CI-only reactions: {len(ci_only)}")

    solved_rows = []
    ambiguous_rows = []
    n_unresolved = 0

    for r in ci_only:
        parsed = parse_equation(r["equation"])
        if not parsed:
            n_unresolved += 1
            continue
        reactants, products = parsed
        net_elem, net_charge, bad = net_balance(reactants, products, cpds)
        mass_ok = (not bad) and all(abs(v) < TOL for v in net_elem.values())
        if not mass_ok:
            # shouldn't happen -- CI-only is defined as mass-balanced -- but
            # don't silently assume it; count as unresolved if it ever does.
            n_unresolved += 1
            continue
        if abs(net_charge) < TOL:
            n_unresolved += 1
            continue
        cands = solve_pairs(reactants, products, cpds, net_charge)
        if len(cands) == 1:
            c = cands[0]
            new_eq = new_equation_str(r["equation"], c)
            # end-to-end verification: rebuild the full term list and
            # recompute mass+charge from scratch (not just the pair algebra)
            new_reactants = Counter(reactants)
            new_products = Counter(products)
            new_reactants[c["reactant"]] = c["new_coef_react"]
            new_products[c["product"]] = c["new_coef_prod"]
            v_elem, v_charge, v_bad = net_balance(new_reactants, new_products, cpds)
            verified = (not v_bad) and all(abs(v) < TOL for v in v_elem.values()) and abs(v_charge) < TOL
            if not verified:
                n_unresolved += 1
                continue
            solved_rows.append({
                "id": r["id"], "old_status": r["status"],
                "old_equation": r["equation"], "new_equation": new_eq,
                "old_definition": r.get("definition"),
                "pair_family": pair_family(c["reactant_name"], c["product_name"]),
                "reactant": c["reactant"], "reactant_name": c["reactant_name"],
                "product": c["product"], "product_name": c["product_name"],
                "old_coef_react": c["old_coef_react"], "new_coef_react": c["new_coef_react"],
                "old_coef_prod": c["old_coef_prod"], "new_coef_prod": c["new_coef_prod"],
            })
        elif len(cands) > 1:
            for c in cands:
                ambiguous_rows.append({
                    "id": r["id"], "old_status": r["status"], "old_equation": r["equation"],
                    "reactant": c["reactant"], "reactant_name": c["reactant_name"],
                    "product": c["product"], "product_name": c["product_name"],
                    "old_coef_react": c["old_coef_react"], "new_coef_react": c["new_coef_react"],
                    "old_coef_prod": c["old_coef_prod"], "new_coef_prod": c["new_coef_prod"],
                })
        else:
            n_unresolved += 1

    n_ambiguous_rxns = len({row["id"] for row in ambiguous_rows})
    print(f"uniquely solved: {len(solved_rows)}")
    print(f"ambiguous (>1 candidate pair): {n_ambiguous_rxns}")
    print(f"unresolved: {n_unresolved}")
    assert len(solved_rows) + n_ambiguous_rxns + n_unresolved == len(ci_only)

    with open(OUTDIR / "solved.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(solved_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(solved_rows)

    if ambiguous_rows:
        with open(OUTDIR / "ambiguous.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(ambiguous_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(ambiguous_rows)

    family_counts = Counter(row["pair_family"] for row in solved_rows)
    with open(OUTDIR / "pair_family_counts.tsv", "w") as fh:
        fh.write("family\tn\n")
        for fam, n in family_counts.most_common():
            fh.write(f"{fam}\t{n}\n")

    with open(OUTDIR / "outcome_counts.tsv", "w") as fh:
        fh.write("outcome\tn\n")
        fh.write(f"solved\t{len(solved_rows)}\n")
        fh.write(f"ambiguous\t{n_ambiguous_rxns}\n")
        fh.write(f"unresolved\t{n_unresolved}\n")

    summary = {
        "msdb_sha": MSDB_SHA,
        "n_ci_only": len(ci_only),
        "n_solved": len(solved_rows),
        "n_ambiguous": n_ambiguous_rxns,
        "n_unresolved": n_unresolved,
        "pair_family_counts": dict(family_counts),
        "coefficient_change_examples": Counter(
            (row["old_coef_react"], row["new_coef_react"]) for row in solved_rows
        ).most_common(10),
    }
    with open(OUTDIR / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=str)

    print(f"\nwrote outputs to {OUTDIR}/")


if __name__ == "__main__":
    main()
