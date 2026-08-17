#!/usr/bin/env python3
"""Side-by-side of how the three vocabularies encode the same chemistry.

Not a computed alignment -- the three alphabets are not in the same space, and
pretending otherwise would be the whole error. This just pulls, for a handful of
functional groups that matter thermodynamically, every entry each vocabulary
carries for them, so the structural differences (protonation states enumerated
vs not, Mg vs not, ring context vs not, whole-molecule entries vs not) are
visible as evidence rather than assertion.
"""
from __future__ import annotations

import json
import re
from pathlib import Path



RES = Path(__file__).resolve().parents[1] / "results"

# (theme, dGPredictor fragment-SMILES regex, GC name regex, eQ group-name regex)
THEMES = [
    ("phosphate / phosphoanhydride", r"P", r"(?i)pho|PO3|PO4|Itriphos", r"PO3|PO2|PO4"),
    ("carboxylate", r"^O=C\(\[?O-?\]?\)|C\(=O\)\[O-\]|^O=C\(O\)", r"(?i)COO|carbamate|acetate|formate|oxalate", r"COO|-C=O|>C=O"),
    ("thiol / thioester / disulfide", r"S", r"(?i)^W?S|thio|smide|triS|FeS", r"-S|C\(=O\)S|SO3|SOO"),
    ("amine / amide", r"N", r"(?i)NH|amide|urea|amine", r"-N|=N|N-CO|NC\(=N\)N"),
    ("aromatic / heteroaromatic", r"[cnops]", r"(?i)aromatic|Nap|trp", r"ring [a-z=<>-]|two fused rings"),
    ("magnesium", r"Mg", r"(?i)^Mg", r"Mg1"),
]


def main() -> None:
    voc = json.loads((RES / "vocabulary_comparison.json").read_text())
    gc_names = voc["modelseed_group_contribution"]["group_names"]
    eq_names = voc["equilibrator_current"]["group_names"]

    learned = [ln.split("\t") for ln in
               (RES / "learned_fingerprints.tsv").read_text().strip().split("\n")[1:]]
    dgp = [(r[1], float(r[2]), int(r[4])) for r in learned]

    out = ["# How each vocabulary encodes the same chemistry",
           "",
           "Counts are entries in that vocabulary matching the theme; the examples "
           "are verbatim group names.", ""]
    rows = []
    for theme, re_dgp, re_gc, re_eq in THEMES:
        d = [g for g, c, n in dgp if re.search(re_dgp, g)]
        g_ = [g for g in gc_names if re.search(re_gc, g)]
        e = [g for g in eq_names if re.search(re_eq, g)]
        rows.append((theme, len(d), len(g_), len(e)))
        out += [f"## {theme}", "",
                f"**fine-tuned dGPredictor** ({len(d)} learned fingerprints): "
                + ", ".join(f"`{x}`" for x in sorted(d, key=len)[:14])
                + (" …" if len(d) > 14 else ""),
                "",
                f"**ModelSEED Group Contribution** ({len(g_)} groups): "
                + ", ".join(f"`{x}`" for x in g_[:14]) + (" …" if len(g_) > 14 else ""),
                "",
                f"**eQuilibrator / component-contribution** ({len(e)} groups): "
                + ", ".join(f"`{x}`" for x in e[:14]) + (" …" if len(e) > 14 else ""),
                ""]

    out += ["## Summary counts", "",
            "| theme | fine-tuned dGPredictor | ModelSEED GC | eQuilibrator |",
            "|---|---:|---:|---:|"]
    for theme, a, b, c in rows:
        out += [f"| {theme} | {a} | {b} | {c} |"]

    # structural properties
    n_charged_eq = sum(1 for g in eq_names if re.search(r"Z-?[1-9]", g))
    n_mg_eq = sum(1 for g in eq_names if "Mg1" in g)
    n_ring_eq = sum(1 for g in eq_names if g.startswith("ring ") or "fused rings" in g)
    n_charged_dgp = sum(1 for g, _, _ in dgp if "+]" in g or "-]" in g)
    n_arom_dgp = sum(1 for g, _, _ in dgp if any(ch in g for ch in "cnops"))
    out += ["", "## Structural properties of the alphabets", "",
            "| property | fine-tuned dGPredictor | ModelSEED GC | eQuilibrator |",
            "|---|---|---|---|",
            "| unit | RDKit atom environment, canonical fragment SMILES | "
            "named chemical group | named chemical group |",
            "| built by | RDKit, automatically, per compound | MFAToolkit rules | "
            "curated SMARTS-style rules |",
            f"| protonation state in the label | no — implicit in the pH-7 SMILES "
            f"({n_charged_dgp} learned fragments carry an explicit charge) | no | "
            f"yes — every group is `[Hn Zq Mgm]` ({n_charged_eq} charged variants) |",
            f"| Mg binding in the label | no | one free-ion group (`Mg`) | "
            f"yes ({n_mg_eq} Mg-bound variants) |",
            f"| ring context in the label | implicit (aromatic lowercase atoms in "
            f"{n_arom_dgp} fragments) | yes (`RW…`, `T…`, `HeteroAromatic`) | "
            f"yes ({n_ring_eq} explicit `ring`/`fused rings` groups) |",
            "| whole-molecule entries | no | yes (`H2O`, `CO2`, `urea`, `acetate`, …) | "
            "yes (50 one-hot placeholder columns for non-decomposable compounds) |",
            "| origin / per-molecule constant | none | `Origin` | "
            "`Origin [H0 Z0 Mg0]` |",
            "| undecomposable marker | reaction is simply not predicted | `NoGroup` | "
            "one-hot placeholder column + RMSE_inf |",
            ]
    (RES / "vocabulary_crosswalk.md").write_text("\n".join(out) + "\n")
    print("\n".join(out[-40:]))
    print(f"\nwrote {RES / 'vocabulary_crosswalk.md'}")


if __name__ == "__main__":
    main()
