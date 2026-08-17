#!/usr/bin/env python3
"""Can the fine-tuned dGPredictor produce a per-MOLECULE energy?

Mechanically, obviously yes: its reaction prediction is exactly
    dG_model(reaction) = sum_i  nu_i * f(compound_i),      f(c) = sum over c's
                                                            fragments of coef
so f(c) is already sitting there as a well-defined number for every compound.
The question is whether f(c) means anything on its own.

The test: compute f(c) for every ModelSEED compound and compare it against the
formation energies ModelSEED already stores from Group Contribution and
eQuilibrator (compound_*.json -> thermodynamics, kcal/mol). If the model carries
a usable absolute scale, f should track them. Then repeat the comparison on
DIFFERENCES across balanced reactions, where any constant-per-atom offset
cancels, to show that the same coefficients are fine there.

That contrast is the answer: reaction energies are differences, and the training
data was made entirely of differences, so the absolute level was never pinned.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

REPO = Path("/scratch/ctaylor/dgpredictor_repo")
DATA = REPO / "data"
DEV = Path("/scratch/ctaylor/tmp/devsnap2")
RES = Path(__file__).resolve().parents[1] / "results"
KJ_PER_KCAL = 4.184


def main() -> None:
    rows = [ln.split("\t") for ln in
            (RES / "finetuned_group_coefficients.tsv").read_text().strip().split("\n")]
    hdr, rows = rows[0], rows[1:]
    c = {n: i for i, n in enumerate(hdr)}
    coef = {(r[c["radius"]], r[c["group"]]): float(r[c["coefficient_kJ_per_mol"]])
            for r in rows}

    dec_r1 = json.loads((DATA / "modelseed_decompose_r1.json").read_text())
    dec_r2 = json.loads((DATA / "modelseed_decompose_r2.json").read_text())

    def f(cpd: str) -> float | None:
        """The per-molecule number implied by the model, kJ/mol."""
        if cpd not in dec_r1:
            return None
        t = 0.0
        for g, n in dec_r1[cpd].items():
            t += n * coef.get(("1", g), 0.0)
        for g, n in dec_r2.get(cpd, {}).items():
            t += n * coef.get(("2", g), 0.0)
        return t

    # reference formation energies already in ModelSEED (kcal/mol -> kJ/mol)
    ref: dict[str, dict[str, float]] = {}
    import glob
    for p in sorted(glob.glob(str(DEV / "Biochemistry" / "compound_*.json"))):
        for e in json.load(open(p)):
            th = e.get("thermodynamics") or {}
            d = {}
            for src in ("Group contribution", "eQuilibrator"):
                v = th.get(src)
                if v and abs(float(v[0])) < 1e6:
                    d[src] = float(v[0]) * KJ_PER_KCAL
            if d:
                ref[e["id"]] = d

    out: dict = {}
    for src in ("Group contribution", "eQuilibrator"):
        pairs = [(f(k), v[src]) for k, v in ref.items()
                 if src in v and f(k) is not None]
        x = np.array([p[0] for p in pairs])
        y = np.array([p[1] for p in pairs])
        out[f"per_molecule_vs_{src}"] = {
            "n_compounds": int(len(x)),
            "pearson_r": round(float(np.corrcoef(x, y)[0, 1]), 3),
            "median_abs_diff_kJ": round(float(np.median(np.abs(x - y))), 1),
            "median_model_value_kJ": round(float(np.median(x)), 1),
            "median_reference_kJ": round(float(np.median(y)), 1),
        }

    # the same coefficients, used on DIFFERENCES across balanced reactions
    stoich = json.loads((DATA / "modelseed_reaction_stoich.json").read_text())
    preds = json.loads((DATA / "modelseed_all_reaction_dG_predictions.json").read_text())
    mine, theirs = [], []
    for rxn, rec in preds.items():
        st = stoich.get(rxn)
        if st is None:
            continue
        vals = [(f(cp), nu) for cp, nu in st.items()]
        if any(v is None for v, _ in vals if _ != 0):
            continue
        tot = sum((v or 0.0) * nu for v, nu in vals)
        mine.append(tot)
        theirs.append(rec["dG_model_only"])
    mine, theirs = np.asarray(mine), np.asarray(theirs)
    out["same_f_summed_over_a_balanced_reaction"] = {
        "n_reactions": int(mine.size),
        "max_abs_diff_vs_shipped_prediction_kJ": float(np.abs(mine - theirs).max()),
        "note": "f() is exactly what the model uses; summing it over a reaction "
                "reproduces the shipped prediction, so f is real -- it is only "
                "its ABSOLUTE level that is undetermined",
    }

    # a few familiar compounds, for the write-up
    demo = {}
    for cpd, name in [("cpd00001", "H2O"), ("cpd00002", "ATP"), ("cpd00008", "ADP"),
                      ("cpd00009", "Phosphate"), ("cpd00003", "NAD"),
                      ("cpd00020", "Pyruvate"), ("cpd00027", "D-Glucose")]:
        v = f(cpd)
        demo[name] = {
            "dGPredictor_per_molecule_kJ": None if v is None else round(v, 1),
            "Group_contribution_kJ": round(ref.get(cpd, {}).get("Group contribution", float("nan")), 1),
            "eQuilibrator_kJ": round(ref.get(cpd, {}).get("eQuilibrator", float("nan")), 1),
        }
    out["examples"] = demo

    (RES / "per_molecule_test.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
