#!/usr/bin/env python3
"""Three-way comparison of the "group" alphabets used by the fine-tuned
dGPredictor, the ModelSEED Group Contribution method, and eQuilibrator.

Fine-tuned dGPredictor
    Atom-centred RDKit environments (Morgan-style, radius 1 and radius 2)
    written out as canonical fragment SMILES. Vocabulary from
    ``data/modelseed_group_names_r{1,2}.txt`` in the freiburgermsu repo; the
    LEARNED subset and its per-group energies come from
    ``rebuild_finetuned_features.py`` (verified to reproduce the shipped
    predictions to 6e-7 kJ/mol).

ModelSEED Group Contribution  (source A)
    The MFAToolkit / Mavrovouniotis-Jankowski named groups, read out of
    ``Biochemistry/Thermodynamics/ModelSEED/*_MolAnalysis.tbl`` column 3, which
    is the literal decomposition ModelSEED stores per structure.

eQuilibrator  (source B)
    Two vintages, both present locally:
      * component-contribution 0.7 / eQuilibrator 3.0 -- ``train_G`` columns of
        the cached ``cc_params.npz``: 163 real groups + 50 one-hot placeholder
        columns for non-decomposable compounds.
      * the 2.x-vintage ``G`` matrix (673 x 241) shipped inside
        ``data/component_contribution_python.mat`` in the dGPredictor repo --
        i.e. the group matrix dGPredictor itself was benchmarked against.

Writes a TSV of the learned fingerprints, per-vocabulary summary JSON, and the
character breakdown used in the write-up.
"""
from __future__ import annotations

import collections
import glob
import json
import re
from pathlib import Path

import numpy as np

REPO = Path("/scratch/ctaylor/dgpredictor_repo")
DATA = REPO / "data"
MSDB = Path("/scratch/ctaylor/ModelSEEDDatabase")
CC_PARAMS = Path("/scratch/ctaylor/eq_cache/equilibrator/cc_params.npz")
EQ3_PYTHON = "/mnt/homes/ctaylor/conda/miniforge3/envs/eq3/bin/python"
OUT = Path(__file__).resolve().parents[1] / "results"
OUT.mkdir(parents=True, exist_ok=True)

HEAVY = re.compile(r"Cl|Br|[BCNOPSFIbcnops]|\[[^\]]+\]")


def n_heavy_atoms(frag: str) -> int:
    """Rough heavy-atom count of a fragment SMILES (bracket atoms count once)."""
    s = re.sub(r"\[[^\]]+\]", "X", frag)
    s = s.replace("Cl", "X").replace("Br", "X")
    return sum(1 for ch in s if ch in "BCNOPSFIbcnopsX")


def main() -> None:
    # ---------------- fine-tuned dGPredictor --------------------------------
    coef_rows = [ln.split("\t") for ln in
                 (OUT / "finetuned_group_coefficients.tsv").read_text().strip().split("\n")]
    header, coef_rows = coef_rows[0], coef_rows[1:]
    col = {n: i for i, n in enumerate(header)}
    learned = [r for r in coef_rows if r[col["used_in_training"]] == "1"]
    declared_r1 = [r for r in coef_rows if r[col["radius"]] == "1"]
    declared_r2 = [r for r in coef_rows if r[col["radius"]] == "2"]
    learned_r1 = [r for r in learned if r[col["radius"]] == "1"]
    learned_r2 = [r for r in learned if r[col["radius"]] == "2"]

    coefs = np.array([float(r[col["coefficient_kJ_per_mol"]]) for r in learned])
    frags = [r[col["group"]] for r in learned]

    charged = [f for f in frags if "+]" in f or "-]" in f]
    aromatic = [f for f in frags if any(c in f for c in "cnops")]
    sizes = [n_heavy_atoms(f) for f in frags]

    # top energies, for the write-up
    order = np.argsort(-np.abs(coefs))
    top = [{"group": frags[i], "radius": int(learned[i][col["radius"]]),
            "coef_kJ_per_mol": round(float(coefs[i]), 2),
            "n_training_compounds": int(learned[i][col["n_training_compounds"]]),
            "n_modelseed_compounds": int(learned[i][col["n_all_modelseed_compounds"]])}
           for i in order[:25]]

    # ---------------- ModelSEED Group Contribution --------------------------
    gc_groups: collections.Counter[str] = collections.Counter()
    gc_per_cpd: list[int] = []
    for f in sorted(glob.glob(str(MSDB / "Biochemistry" / "Thermodynamics"
                                  / "ModelSEED" / "*_MolAnalysis.tbl"))):
        for line in open(f, errors="replace"):
            p = line.rstrip("\n").split("\t")
            if len(p) < 3 or not p[2] or p[2] == "null":
                continue
            names = [e.rsplit(":", 1)[0] for e in p[2].split("|") if ":" in e]
            gc_per_cpd.append(len(names))
            gc_groups.update(set(names))

    # ---------------- eQuilibrator ------------------------------------------
    # cc_params.npz round-trips pandas DataFrames, so read it through the
    # package (in the one env that has it) rather than re-deriving its layout.
    eq = {}
    try:
        import subprocess
        out = subprocess.run(
            [EQ3_PYTHON, "-c",
             "import json, numpy as np;"
             "from component_contribution.parameters import CCModelParameters as P;"
             f"p=P.from_npz('{CC_PARAMS}');"
             "G=p.train_G.values; S=p.train_S.values; M=S.T@G;"
             "print(json.dumps({'cols':[str(c) for c in p.train_G.columns],"
             "'Nc':int(p.dimensions.at['Nc','number']),"
             "'Ng':int(p.dimensions.at['Ng','number']),"
             "'Ng_full':int(p.dimensions.at['Ng_full','number']),"
             "'median_groups_per_compound':float(np.median((G[:,:p.dimensions.at['Ng','number']]!=0).sum(axis=1))),"
             "'rank_of_group_matrix':int(np.linalg.matrix_rank(G[:,:p.dimensions.at['Ng','number']])),"
             "'reaction_feature_matrix_shape':list(M.shape),"
             "'reaction_feature_matrix_rank':int(np.linalg.matrix_rank(M)),"
             "'n_columns_used_in_training':int(((M!=0).sum(axis=0)>0).sum())}))"],
            capture_output=True, text=True, check=True)
        eqd = json.loads(out.stdout)
        eq = {
            "vintage": "component-contribution 0.7 / eQuilibrator 3.0",
            "n_real_groups": eqd["Ng"],
            "n_placeholder_columns": eqd["Ng_full"] - eqd["Ng"],
            "n_training_compounds": eqd["Nc"],
            "group_names": eqd["cols"][:eqd["Ng"]],
            "median_groups_per_compound": eqd["median_groups_per_compound"],
            "rank_of_group_matrix": eqd["rank_of_group_matrix"],
            "reaction_feature_matrix_shape": eqd["reaction_feature_matrix_shape"],
            "reaction_feature_matrix_rank": eqd["reaction_feature_matrix_rank"],
            "n_columns_used_in_training": eqd["n_columns_used_in_training"],
        }
    except Exception as exc:  # pragma: no cover - environment dependent
        eq = {"error": f"could not read {CC_PARAMS} via {EQ3_PYTHON}: {exc}"}

    # 2.x-vintage G from the dGPredictor repo itself
    from scipy.io import loadmat
    cc = loadmat(DATA / "component_contribution_python.mat")
    G2, S2 = cc["G"], cc["train_S"]
    M2 = S2.T @ G2
    eq2 = {
        "vintage": "component-contribution 2.x (shipped in dGPredictor repo)",
        "shape": list(G2.shape),
        "n_columns": int(G2.shape[1]),
        "n_training_compounds": int(G2.shape[0]),
        "median_groups_per_compound": float(np.median((G2 != 0).sum(axis=1))),
        "rank_of_group_matrix": int(np.linalg.matrix_rank(G2)),
        "reaction_feature_matrix_shape": list(M2.shape),
        "reaction_feature_matrix_rank": int(np.linalg.matrix_rank(M2)),
        "n_columns_used_in_training": int(((M2 != 0).sum(axis=0) > 0).sum()),
    }

    summary = {
        "finetuned_dgpredictor": {
            "construction": "RDKit atom-environment fragments (radius 1 and 2), "
                            "canonical fragment SMILES, one per heavy atom per radius",
            "declared_vocabulary": {"r1": len(declared_r1), "r2": len(declared_r2),
                                    "total": len(declared_r1) + len(declared_r2)},
            "learned_vocabulary": {"r1": len(learned_r1), "r2": len(learned_r2),
                                   "total": len(learned)},
            "frac_declared_that_is_learned":
                round(len(learned) / (len(declared_r1) + len(declared_r2)), 4),
            "character": {
                "carry_explicit_charge": len(charged),
                "carry_aromatic_atoms": len(aromatic),
                "heavy_atoms_per_fragment": {
                    "min": int(min(sizes)), "median": float(np.median(sizes)),
                    "max": int(max(sizes))},
            },
            "coefficients_kJ_per_mol": {
                "median_abs": round(float(np.median(np.abs(coefs))), 3),
                "max_abs": round(float(np.abs(coefs).max()), 2),
                "n_above_50": int((np.abs(coefs) > 50).sum()),
            },
            "top25_by_abs_energy": top,
        },
        "modelseed_group_contribution": {
            "construction": "hand-curated named chemical groups "
                            "(Mavrovouniotis / Jankowski), assigned by MFAToolkit",
            "vocabulary_size": len(gc_groups),
            "median_groups_per_compound": float(np.median(gc_per_cpd)),
            "structures_decomposed": len(gc_per_cpd),
            "has_origin_term": "Origin" in gc_groups,
            "has_undecomposable_marker": "NoGroup" in gc_groups,
            "n_structures_marked_NoGroup": gc_groups.get("NoGroup", 0),
            "group_names": sorted(gc_groups),
        },
        "equilibrator_current": eq,
        "equilibrator_in_dgpredictor_repo": eq2,
    }
    (OUT / "vocabulary_comparison.json").write_text(json.dumps(summary, indent=2))

    # learned fingerprints, sorted by |energy|
    lines = ["radius\tgroup\tcoef_kJ_per_mol\tn_heavy_atoms\tn_training_compounds\t"
             "n_modelseed_compounds"]
    for i in order:
        r = learned[i]
        lines.append(f"{r[col['radius']]}\t{frags[i]}\t{coefs[i]:.4f}\t"
                     f"{n_heavy_atoms(frags[i])}\t{r[col['n_training_compounds']]}\t"
                     f"{r[col['n_all_modelseed_compounds']]}")
    (OUT / "learned_fingerprints.tsv").write_text("\n".join(lines) + "\n")

    print(json.dumps({k: (v if k != "finetuned_dgpredictor" else
                          {kk: vv for kk, vv in v.items() if kk != "top25_by_abs_energy"})
                      for k, v in summary.items()
                      if k != "modelseed_group_contribution"}, indent=2)[:4000])
    g = summary["modelseed_group_contribution"]
    print(json.dumps({k: v for k, v in g.items() if k != "group_names"}, indent=2))
    print(f"\nwrote {OUT/'vocabulary_comparison.json'} and {OUT/'learned_fingerprints.tsv'}")


if __name__ == "__main__":
    main()
