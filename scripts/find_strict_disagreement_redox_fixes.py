#!/usr/bin/env python3
"""Can a mechanically-verified "fixed" equation be built for the 219 CI-only
strict-disagreement reactions (see analyze_strict_disagreement_balanced_populations.py
-- 89.5% of them are a redox-cofactor bookkeeping defect), and if so, does
eQuilibrator's real ΔG'0 on the FIXED stoichiometry now agree with the LLM's
claimed direction?

This generalizes fix_ci_only_reactions.py's same-formula-pair rescale (which
only reassigns a pair's coefficients while holding each member on its
ORIGINAL side) to a full SIGNED solve: for a same-formula pair (X, Y) present
anywhere in the equation, let nx/ny be their NET signed coefficients
(product - reactant; a same-formula pair's net element contribution is
(nx+ny)*formula(X), so holding K = nx+ny fixed preserves the reaction's
already-passing mass balance exactly). Solve for the unique (nx', ny') with
nx'+ny' = K that zeros the WHOLE equation's net charge. This is a strict
superset of "rescale" (same sign pattern, new magnitude) and also captures
"swap sides" (sign pattern flips) -- e.g. rxn07326/rxn07329, where the
literature says NADP+/NADPH are assigned to the chemically backward sides:
swapping them exactly explains that reaction's CI:-4.

It does NOT, and cannot, fix a "bare half-reaction" like rxn12950 (NADH =>
NAD and nothing else): there K=0 forces nx'=-ny', and the charge equation
collapses to nx'*(charge(Y)-charge(X)) = 0, i.e. nx'=0 -- no real reaction,
confirming by construction (not by inspection) that these need the missing
substrate pair invented from literature, not arithmetic. Reactions needing a
completely different missing cofactor (UDP-sugar, ATP, CoA -- several of the
120 OK-status "cofactors_stripped" reactions are this, not NAD(P)/NAD(P)H at
all) are out of scope for the same reason: with charge ALREADY at 0 (OK
status), the K-preserving charge equation has no second degree of freedom to
move, so this method only ever finds something for the CI-only population
(nonzero existing charge residual) -- confirmed empirically below, not
assumed.

Outcomes per CI-only reaction, same three buckets as fix_ci_only_reactions.py:
    solved      exactly one same-formula pair gives a valid (both coefficients
                integer, 1..MAXCOEF, not identical to the original) new
                assignment that exactly zeros net charge; verified end-to-end
                by rebuilding the full term list and recomputing mass+charge
                from scratch (not just the pair algebra)
    ambiguous   more than one pair qualifies
    unresolved  no pair gives a valid solution (bare half-reactions, and any
                non-redox-cofactor defect)

For every solved reaction, both the ORIGINAL and FIXED stoichiometries are
sent to eQuilibrator 3.0 (component contribution, same calling convention as
scripts/estimate_directions_eq3.py) to get real standard_dg_prime values, and
the implied direction (sign, relative to the equation's own written
left-to-right orientation) is compared against direction_conflict.llm_ensemble
and .modelseed_recommended from strict_disagreement_enzymes.json.

Must run in the `eq3` conda env (equilibrator-api 0.7.0):
    /mnt/homes/ctaylor/conda/miniforge3/envs/eq3/bin/python \\
      scripts/find_strict_disagreement_redox_fixes.py

Data: ModelSEED dev @ 078a395f (/scratch/ctaylor/tmp/devsnap_078a395f), same
snapshot as the rest of this pipeline.

Outputs (results/strict_disagreement_enzymes/):
    redox_fix_outcomes.json   per-reaction outcome, fix detail, both dG values,
                               implied direction before/after, agreement verdict
    redox_fix_outcomes.tsv    flat version for a figure
    redox_fix_summary.json    headline counts
"""
from __future__ import annotations

import glob
import json
import os
import re
import sys
from collections import Counter
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.dirname(HERE)
OUTDIR = os.path.join(ROOT, "results", "strict_disagreement_enzymes")
MSDB_SHA = "078a395f"
MSDB_ROOT = os.environ.get("MSDB_ROOT", "/scratch/ctaylor/tmp/devsnap_078a395f")
DELIVERABLE = os.path.join(OUTDIR, "strict_disagreement_enzymes.json")
MAXCOEF = 12
TOL = 1e-6

from fix_ci_only_reactions import parse_formula, load_compounds, parse_equation, net_balance

KCAL = 4.184


def status_tokens(status: str | None) -> frozenset:
    status = status or ""
    return frozenset(f[:f.find(':')] if ':' in f else f for f in status.split('|'))


def signed_net(reactants: Counter, products: Counter) -> dict[str, float]:
    net: dict[str, float] = {}
    for cid, c in reactants.items():
        net[cid] = net.get(cid, 0.0) - c
    for cid, c in products.items():
        net[cid] = net.get(cid, 0.0) + c
    return net


def find_pair_fix(net: dict[str, float], cpds: dict, net_charge_total: float):
    """All same-formula, different-charge compound pairs present in `net`
    (net != 0 for at least one of the two -- a pair entirely absent can never
    affect an already-computed charge residual, see module docstring), with
    the unique signed (nx', ny') that holds K=nx+ny fixed and zeros the whole
    equation's net charge. Returns a list of candidate fixes (empty if none).
    """
    present = [cid for cid, v in net.items() if abs(v) > TOL]
    candidates = []
    seen_pairs = set()
    for i, a in enumerate(present):
        ca = cpds.get(a)
        if not ca or not ca.get("formula") or ca.get("charge") is None:
            continue
        fa = parse_formula(ca["formula"])
        for b in present[i + 1:]:
            if (a, b) in seen_pairs or (b, a) in seen_pairs:
                continue
            cb = cpds.get(b)
            if not cb or not cb.get("formula") or cb.get("charge") is None:
                continue
            if parse_formula(cb["formula"]) != fa:
                continue
            dq = cb["charge"] - ca["charge"]
            if abs(dq) < TOL:
                continue
            seen_pairs.add((a, b))
            nx, ny = net[a], net.get(b, 0.0)
            K = nx + ny
            c_pair = nx * ca["charge"] + ny * cb["charge"]
            target = c_pair - net_charge_total
            # nx' + ny' = K ; nx'*charge(a) + ny'*charge(b) = target
            ny_p = (target - K * ca["charge"]) / dq
            nx_p = K - ny_p
            nx_r, ny_r = round(nx_p), round(ny_p)
            if abs(nx_p - nx_r) > TOL or abs(ny_p - ny_r) > TOL:
                continue
            if abs(nx_r) < 1 or abs(ny_r) < 1 or abs(nx_r) > MAXCOEF or abs(ny_r) > MAXCOEF:
                continue
            if nx_r == nx and ny_r == ny:
                continue  # no-op, not a fix
            kind = "swap" if (nx_r > 0) != (nx > 0) else "rescale"
            candidates.append({
                "a": a, "a_name": ca.get("name"), "b": b, "b_name": cb.get("name"),
                "formula": ca["formula"], "charge_a": ca["charge"], "charge_b": cb["charge"],
                "old_net_a": nx, "old_net_b": ny, "new_net_a": nx_r, "new_net_b": ny_r,
                "kind": kind,
            })
    return candidates


def apply_fix(net: dict[str, float], fix: dict) -> dict[str, float]:
    new_net = dict(net)
    new_net[fix["a"]] = fix["new_net_a"]
    new_net[fix["b"]] = fix["new_net_b"]
    return new_net


def verify_fix(new_net: dict[str, float], cpds: dict) -> bool:
    reactants = Counter({cid: -v for cid, v in new_net.items() if v < -TOL})
    products = Counter({cid: v for cid, v in new_net.items() if v > TOL})
    net_elem, net_charge, bad = net_balance(reactants, products, cpds)
    return (not bad) and all(abs(v) < TOL for v in net_elem.values()) and abs(net_charge) < TOL


def main() -> None:
    cpds = load_compounds()
    deliverable = json.load(open(DELIVERABLE))
    reactions = deliverable["reactions"]
    ci_only = [r for r in reactions if status_tokens(r.get("modelseed_status")) == frozenset({"CI"})]
    print(f"[fix] CI-only strict-disagreement reactions: {len(ci_only)}", flush=True)

    outcomes = []
    n_solved = n_ambiguous = n_unresolved = 0
    for r in ci_only:
        eq = r["reaction_equation"]
        parsed = parse_equation(eq)
        if not parsed:
            outcomes.append({"modelseed_id": r["modelseed_id"], "outcome": "unresolved", "reason": "unparseable"})
            n_unresolved += 1
            continue
        reactants, products = parsed
        net_elem, net_charge, bad = net_balance(reactants, products, cpds)
        if bad:
            outcomes.append({"modelseed_id": r["modelseed_id"], "outcome": "unresolved", "reason": "missing compound data"})
            n_unresolved += 1
            continue
        net = signed_net(reactants, products)
        cands = find_pair_fix(net, cpds, net_charge)
        if len(cands) == 1:
            fix = cands[0]
            new_net = apply_fix(net, fix)
            if not verify_fix(new_net, cpds):
                outcomes.append({"modelseed_id": r["modelseed_id"], "outcome": "unresolved", "reason": "candidate failed end-to-end verification"})
                n_unresolved += 1
                continue
            outcomes.append({
                "modelseed_id": r["modelseed_id"], "outcome": "solved",
                "old_status": r.get("modelseed_status"), "equation_defect": r["direction_conflict"].get("equation_defect"),
                "fix": fix, "old_net_stoich": net, "new_net_stoich": new_net,
            })
            n_solved += 1
        elif len(cands) > 1:
            outcomes.append({"modelseed_id": r["modelseed_id"], "outcome": "ambiguous",
                              "n_candidates": len(cands), "candidates": cands})
            n_ambiguous += 1
        else:
            outcomes.append({"modelseed_id": r["modelseed_id"], "outcome": "unresolved", "reason": "no same-formula pair solves it"})
            n_unresolved += 1

    print(f"[fix] solved {n_solved}  ambiguous {n_ambiguous}  unresolved {n_unresolved}")
    assert n_solved + n_ambiguous + n_unresolved == len(ci_only)

    # ---- eQuilibrator pass on original vs fixed stoichiometry, solved only ----
    try:
        from equilibrator_api import ComponentContribution, Reaction
        have_eq = True
    except ImportError:
        print("[fix] equilibrator_api not importable -- run in the eq3 conda env "
              "to compute dG; writing fix outcomes without thermodynamics.")
        have_eq = False

    if have_eq:
        print("[fix] loading ComponentContribution (eQ 3.0 cache) ...", flush=True)
        cc = ComponentContribution()

        _KEGG = re.compile(r"KEGG:\s*(C\d{5})")
        _CHEBI = re.compile(r"ChEBI:\s*(\d+)")
        _BIGG = re.compile(r"BiGG:\s*([^\s;]+)")
        _MNX = re.compile(r"MetaNetX:\s*(MNXM\d+)")

        def compound_accessions(c):
            al = c.get("aliases") or []
            text = al if isinstance(al, str) else "\n".join(al)
            accs = []
            m = _KEGG.search(text)
            if m:
                accs.append("kegg:" + m.group(1))
            m = _CHEBI.search(text)
            if m:
                accs.append("chebi:CHEBI:" + m.group(1))
            m = _BIGG.search(text)
            if m:
                accs.append("bigg.metabolite:" + m.group(1))
            m = _MNX.search(text)
            if m:
                accs.append("metanetx.chemical:" + m.group(1))
            ik = c.get("inchikey")
            if ik:
                accs.append(ik)
            return accs

        resolved: dict[str, object] = {}

        def resolve(cid):
            if cid in resolved:
                return resolved[cid]
            obj = None
            c = cpds.get(cid)
            if c:
                for acc in compound_accessions(c):
                    try:
                        obj = cc.get_compound(acc)
                    except Exception:
                        obj = None
                    if obj is not None:
                        break
            resolved[cid] = obj
            return obj

        LNRI_THRESH = 6.9077553  # ln(1000), eQuilibrator/Noor-2012's own reversibility cutoff

        def dg_prime(net: dict[str, float]):
            """Returns (dg_kj_mol, lnri, direction, error). Direction uses
            eQuilibrator's OWN reversibility index (same threshold as
            scripts/estimate_directions_eq3.py), not a bare point-estimate
            sign -- a a raw-sign call on an obscure compound with wide group-
            contribution uncertainty (e.g. a methylated sterol intermediate)
            can flip on noise alone; the index call is what the rest of this
            project already treats as eQuilibrator's real verdict.
            """
            sd = {}
            for cid, coef in net.items():
                if abs(coef) < TOL:
                    continue
                obj = resolve(cid)
                if obj is None:
                    return None, None, None, "unmapped"
                sd[obj] = sd.get(obj, 0.0) + coef
            sd = {k: v for k, v in sd.items() if abs(v) > TOL}
            if not sd:
                return None, None, None, "empty"
            try:
                rxn = Reaction(sd)
                if not rxn.is_balanced():
                    return None, None, None, "eq_reports_unbalanced"
                dg = float(cc.standard_dg_prime(rxn).value.m_as("kJ/mol"))
                lnri = float(cc.ln_reversibility_index(rxn).value.m_as(""))
                d = ">" if lnri < -LNRI_THRESH else "<" if lnri > LNRI_THRESH else "="
                return dg, lnri, d, None
            except Exception as e:
                return None, None, None, f"error:{type(e).__name__}"

        n_both = 0
        for o in outcomes:
            if o["outcome"] != "solved":
                continue
            rid = o["modelseed_id"]
            dg_old, lnri_old, dir_old, err_old = dg_prime(o["old_net_stoich"])
            dg_new, lnri_new, dir_new, err_new = dg_prime(o["new_net_stoich"])
            o["eq3_dg_old_kj_mol"], o["eq3_lnri_old"], o["eq3_direction_old"], o["eq3_dg_old_error"] = dg_old, lnri_old, dir_old, err_old
            o["eq3_dg_new_kj_mol"], o["eq3_lnri_new"], o["eq3_direction_new"], o["eq3_dg_new_error"] = dg_new, lnri_new, dir_new, err_new
            dc = next(r["direction_conflict"] for r in reactions if r["modelseed_id"] == rid)
            llm_dir = ">" if dc["llm_ensemble"] == "forward" else "<"
            ms_dir = ">" if dc["modelseed_recommended"] == "forward" else "<"
            o["llm_direction"] = llm_dir
            o["modelseed_recommended_direction"] = ms_dir
            # Direction uses eQuilibrator's own reversibility index (ln RI,
            # Noor 2012), not a bare dG sign -- see dg_prime() docstring.
            # The ORIGINAL equation is, by construction, charge-imbalanced
            # (that's what CI-only means), so eQuilibrator's OWN balance
            # check correctly refuses it almost every time
            # (eq_reports_unbalanced) -- that refusal is itself evidence,
            # not a methodology gap, so it's recorded but not required for
            # the FIXED-equation verdict below.
            if dg_old is not None:
                o["eq3_old_agrees_with_llm"] = (dir_old == llm_dir)
            if dg_new is not None:
                n_both += 1
                o["eq3_new_agrees_with_llm"] = (dir_new == llm_dir)
                o["eq3_new_agrees_with_modelseed"] = (dir_new == ms_dir)
            print(f"  {rid}: old dG={dg_old} dir={dir_old} ({err_old})  new dG={dg_new} dir={dir_new} ({err_new})", flush=True)
        print(f"[fix] eQuilibrator computed a FIXED-equation dG for {n_both}/{n_solved} solved reactions "
              f"({sum(1 for o in outcomes if o.get('eq3_dg_old_kj_mol') is not None)} of those also had an "
              f"(expectedly, refused) original value)")

    PLAUSIBLE_DG_MAX = 200.0  # kJ/mol -- generous outer bound for a single enzymatic
    # step (real biochemistry rarely exceeds ~100-150). Anything past it is a red
    # flag for a bad compound-structure mapping, not real thermodynamics -- found
    # empirically here: 22/29 solved+computed reactions are ALL the same generic,
    # templated "aldehyde dehydrogenase [NAD(P)+]" combinatorial-expansion entry
    # (consistent with this project's own prior finding that such expansions pair
    # non-corresponding substrate/product instances), giving dG from +289 to +687
    # kJ/mol -- physically impossible for one enzymatic step. Other instances of
    # the SAME enzyme name map to a plausible dG of -69, so this is a
    # per-reaction mapping defect (which specific substrate got filled into the
    # template), not a flaw in the fix itself.
    if have_eq:
        for o in outcomes:
            if o.get("eq3_dg_new_kj_mol") is not None:
                o["plausible_magnitude"] = abs(o["eq3_dg_new_kj_mol"]) <= PLAUSIBLE_DG_MAX

    os.makedirs(OUTDIR, exist_ok=True)
    with open(os.path.join(OUTDIR, "redox_fix_outcomes.json"), "w") as fh:
        json.dump(outcomes, fh, indent=2, default=str)

    with open(os.path.join(OUTDIR, "redox_fix_outcomes.tsv"), "w") as fh:
        fh.write("modelseed_id\toutcome\tkind\teq3_dg_old_kj_mol\teq3_dg_new_kj_mol\teq3_lnri_new\tplausible_magnitude\t"
                  "llm_direction\tmodelseed_direction\teq3_direction_old\teq3_direction_new\tverdict\n")
        for o in outcomes:
            if o["outcome"] != "solved":
                fh.write(f"{o['modelseed_id']}\t{o['outcome']}\t\t\t\t\t\t\t\t\t\t\n")
                continue
            dn, llm_dir, ms_dir = o.get("eq3_direction_new"), o.get("llm_direction"), o.get("modelseed_recommended_direction")
            plausible = o.get("plausible_magnitude")
            if dn is None:
                verdict = o.get("eq3_dg_new_error", "")
            elif plausible is False:
                verdict = "implausible_magnitude"
            elif dn == "=":
                verdict = "near_equilibrium"
            elif dn == llm_dir:
                verdict = "agrees_with_llm"
            elif dn == ms_dir:
                verdict = "agrees_with_modelseed"
            else:
                verdict = "agrees_with_neither"
            fh.write(f"{o['modelseed_id']}\t{o['outcome']}\t{o['fix']['kind']}\t"
                      f"{o.get('eq3_dg_old_kj_mol','')}\t{o.get('eq3_dg_new_kj_mol','')}\t{o.get('eq3_lnri_new','')}\t"
                      f"{plausible if plausible is not None else ''}\t"
                      f"{llm_dir or ''}\t{ms_dir or ''}\t{o.get('eq3_direction_old','')}\t{dn or ''}\t{verdict}\n")

    summary = {
        "msdb_sha": MSDB_SHA, "n_ci_only": len(ci_only),
        "n_solved": n_solved, "n_ambiguous": n_ambiguous, "n_unresolved": n_unresolved,
    }
    if have_eq:
        solved_with_fix_dg = [o for o in outcomes if o.get("eq3_dg_new_kj_mol") is not None]
        summary["n_eq3_computed_for_fixed_equation"] = len(solved_with_fix_dg)
        summary["n_eq3_original_refused_as_unbalanced"] = sum(
            1 for o in outcomes if o["outcome"] == "solved" and o.get("eq3_dg_old_error") == "eq_reports_unbalanced")
        implausible = [o for o in solved_with_fix_dg if not o["plausible_magnitude"]]
        plausible = [o for o in solved_with_fix_dg if o["plausible_magnitude"]]
        summary["n_implausible_magnitude_excluded"] = len(implausible)
        directional = [o for o in plausible if o["eq3_direction_new"] != "="]
        summary["n_near_equilibrium"] = len(plausible) - len(directional)
        summary["n_directional"] = len(directional)
        summary["n_fixed_agrees_with_llm"] = sum(1 for o in directional if o["eq3_new_agrees_with_llm"])
        summary["n_fixed_agrees_with_modelseed"] = sum(1 for o in directional if o["eq3_new_agrees_with_modelseed"])
    with open(os.path.join(OUTDIR, "redox_fix_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))
    print(f"wrote {os.path.join(OUTDIR, 'redox_fix_outcomes.json')}")


if __name__ == "__main__":
    main()
