#!/usr/bin/env python3
"""Compare the thermo-source variants of the core-model sweep.

Consumes ``build_graded_direction_maps.py`` + ``run_graded_fba_all_models.py``
and answers three questions:

  1. GROWTH -- how many of the 5,683 core models grow under each variant, and
     which models change verdict. Growth counts are reported alongside a
     permissiveness measure, because a map that calls everything reversible
     will "win" on growth without being right about anything.

  2. DIRECTION ACCURACY -- for reactions with a TECRDB measurement, run the
     cascade on the EXPERIMENTAL dG'^0 to get a reference direction, then score
     each variant against it. ``graded``/``graded_trusted`` are excluded from
     this scoring as circular (they use TECRDB); ``graded_heldout`` is the
     variant that can be scored.

  3. THE CORE SET -- grade distribution and source mix over the 239 reactions
     that appear in at least one core model, which is where all of this
     actually bites.

Outputs (results/thermo_grades_fba/):
    variant_growth.tsv        per-variant growth totals + deltas vs implicit
    variant_agreement.tsv     pairwise growth-verdict disagreement
    direction_accuracy.tsv    vs the TECRDB reference, all + core + hard subset
    core_reaction_grades.tsv  one row per core reaction: grades, ops, pick
"""
from __future__ import annotations

import collections
import itertools
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
MSDB_CODE = Path(os.environ.get("MSDB_CODE", "/scratch/ctaylor/ModelSEEDDatabase"))
DATA = ANALYSIS_DIR / "results" / "thermo_grades_fba"
GRADES = ANALYSIS_DIR / "results" / "thermo_grades"
MODELS_DIR = ANALYSIS_DIR / "data" / "core_models_kegg2"
TECRDB_CSV = Path(os.environ.get(
    "TECRDB_COMPARISON",
    "/scratch/ctaylor/dgpredictor_tecrdb/results/tecrdb_vs_dgpredictor_modelseed.csv"))

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(MSDB_CODE / "Scripts" / "Thermodynamics"))

VARIANTS = ["implicit", "gc", "eq", "dgpms", "graded", "graded_trusted", "graded_heldout"]
MAPPED = VARIANTS[1:]
CIRCULAR = {"graded", "graded_trusted"}   # contain TECRDB, cannot be scored on it
PRETTY = {"implicit": "model's own bounds", "gc": "Group Contribution only",
          "eq": "eQuilibrator only", "dgpms": "dGPredictor-ModelSEED only",
          "graded": "graded (recommended)", "graded_trusted": "graded, SILVER floor",
          "graded_heldout": "graded, TECRDB held out"}


def core_reaction_set() -> collections.Counter:
    """{seed_rxn_id: number of core models containing it}, read straight from
    the model JSONs (cobra is not needed to count annotations)."""
    from seed_annotation import seed_id
    counts = collections.Counter()
    for path in sorted(MODELS_DIR.glob("*.json")):
        d = json.load(open(path))
        counts.update({s for s in (seed_id(r) for r in d.get("reactions", [])) if s})
    return counts


def tecrdb_reference() -> dict:
    """{rxn: (operator, match_tier)} -- the direction the cascade returns when
    fed the EXPERIMENTAL dG'^0 instead of a prediction."""
    from reversibility_heuristics import DEFAULT_HEURISTICS, explicit_energy, run_reversibility
    from build_graded_direction_maps import load_reactions
    rx = load_reactions()
    t = pd.read_csv(TECRDB_CSV)
    t["dg"] = t.tecrdb_dG_kJ / 4.184
    t["sd"] = t.tecrdb_dG_sd_kJ / 4.184
    t = t.sort_values("match_tier").groupby("modelseed_rxn", as_index=False).first()
    ref = {}
    for r in t.itertuples():
        e = rx.get(r.modelseed_rxn)
        if e is None or e.get("status") == "EMPTY":
            continue
        _, op, _ = run_reversibility(e, explicit_energy(float(r.dg), float(r.sd)),
                                     DEFAULT_HEURISTICS)
        if op:
            ref[r.modelseed_rxn] = (op, r.match_tier)
    return ref


def growth_tables(res: pd.DataFrame, cov: pd.DataFrame, core: set) -> tuple:
    grows = {v: res[f"fba_grows_{v}"].astype(bool) for v in VARIANTS}
    rows = []
    for v in VARIANTS:
        has = cov[f"has_{v}"] if v != "implicit" else None
        rev = np.nan if has is None else float((cov.loc[has, f"op_{v}"] == "=").mean())
        core_has = None if has is None else (has & cov.index.isin(core))
        rows.append({
            "variant": v, "label": PRETTY[v],
            "n_directions": 0 if has is None else int(has.sum()),
            "n_core_directions": 0 if has is None else int(core_has.sum()),
            "frac_reversible": rev,
            "frac_reversible_core": np.nan if has is None else
                float((cov.loc[core_has, f"op_{v}"] == "=").mean()),
            "n_grows": int(grows[v].sum()),
            "pct_grows": float(grows[v].mean()),
            "gained_vs_implicit": int((grows[v] & ~grows["implicit"]).sum()),
            "lost_vs_implicit": int((~grows[v] & grows["implicit"]).sum()),
            "median_bounds_changed": float(res[f"fba_n_bounds_changed_{v}"].median()),
        })
    pair = [{"a": a, "b": b, "n_differ": int((grows[a] != grows[b]).sum()),
             "a_grows_only": int((grows[a] & ~grows[b]).sum()),
             "b_grows_only": int((grows[b] & ~grows[a]).sum())}
            for a, b in itertools.combinations(VARIANTS, 2)]
    return pd.DataFrame(rows), pd.DataFrame(pair)


def direction_accuracy(cov: pd.DataFrame, ref: dict, core: set) -> pd.DataFrame:
    se = [r for r in ref if ref[r][1] == "stereo_exact" and r in cov.index]
    subsets = {
        "all_stereo_exact": se,
        "core_models_only": [r for r in se if r in core],
        "reference_directional": [r for r in se if ref[r][0] != "="],
        "reference_reversible": [r for r in se if ref[r][0] == "="],
    }
    rows = []
    for name, ids in subsets.items():
        for v in MAPPED:
            n = ok = 0
            for r in ids:
                if not cov.at[r, f"has_{v}"]:
                    continue
                n += 1
                ok += (cov.at[r, f"op_{v}"] == ref[r][0])
            rows.append({"subset": name, "variant": v, "label": PRETTY[v],
                         "n_scored": n, "n_correct": ok,
                         "accuracy": ok / n if n else np.nan,
                         "circular": v in CIRCULAR})
    return pd.DataFrame(rows)


def core_table(cov: pd.DataFrame, core_counts: collections.Counter, ref: dict) -> pd.DataFrame:
    wide = pd.read_csv(GRADES / "source_grades_wide.tsv", sep="\t", low_memory=False)
    wide = wide.set_index("rxn")
    rows = []
    for rxn, n_models in sorted(core_counts.items(), key=lambda x: -x[1]):
        row = {"rxn": rxn, "n_core_models": n_models}
        if rxn in wide.index:
            w = wide.loc[rxn]
            row.update({"name": w["name"], "ec": w.get("ec"),
                        "grade_GC": w.get("grade_GC"), "grade_EQ": w.get("grade_EQ"),
                        "grade_DGPMS": w.get("grade_DGPMS"),
                        "grade_TECRDB": w.get("grade_TECRDB"),
                        "best_grade": w.get("best_grade"),
                        "best_source": w.get("best_source"), "birge": w.get("birge")})
        if rxn in cov.index:
            for v in MAPPED:
                row[f"op_{v}"] = cov.at[rxn, f"op_{v}"] if cov.at[rxn, f"has_{v}"] else ""
            row["graded_pick"] = cov.at[rxn, "src_graded"]
        if rxn in ref:
            row["tecrdb_direction"] = ref[rxn][0]
            row["tecrdb_match_tier"] = ref[rxn][1]
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    res = pd.read_csv(DATA / "model_results.csv")
    cov = pd.read_csv(DATA / "rxn_source_coverage.csv", low_memory=False).set_index("rxn_id")
    core_counts = core_reaction_set()
    core = set(core_counts)
    print(f"{len(res)} models, {len(core)} distinct core reactions "
          f"({len(core & set(cov.index))} of them non-EMPTY in the snapshot)")

    growth, pair = growth_tables(res, cov, core)
    print("\n=== growth ===")
    print(growth[["label", "n_directions", "n_core_directions", "frac_reversible_core",
                  "n_grows", "pct_grows", "gained_vs_implicit",
                  "lost_vs_implicit"]].to_string(index=False))

    ref = tecrdb_reference()
    acc = direction_accuracy(cov, ref, core)
    print("\n=== direction accuracy vs the experimental reference ===")
    for subset, g in acc.groupby("subset", sort=False):
        print(f"  {subset}:")
        for r in g.itertuples():
            flag = "   [circular -- uses TECRDB]" if r.circular else ""
            print(f"     {r.label:32s} {r.n_correct:4d}/{r.n_scored:4d}  "
                  f"{r.accuracy:6.1%}{flag}")

    ct = core_table(cov, core_counts, ref)
    print("\n=== the 239 core reactions ===")
    print("  best grade:", ct.best_grade.value_counts().to_dict(),
          " no source:", int(ct.best_grade.isna().sum()))
    print("  graded map picked:", ct.loc[ct.graded_pick.notna()
                                         & (ct.graded_pick != ""), "graded_pick"]
          .value_counts().to_dict())
    print("  with a TECRDB measurement:", int(ct.tecrdb_direction.notna().sum()))

    growth.to_csv(DATA / "variant_growth.tsv", sep="\t", index=False, float_format="%.4f")
    pair.to_csv(DATA / "variant_agreement.tsv", sep="\t", index=False)
    acc.to_csv(DATA / "direction_accuracy.tsv", sep="\t", index=False, float_format="%.4f")
    ct.to_csv(DATA / "core_reaction_grades.tsv", sep="\t", index=False, float_format="%.4f")
    print(f"\nwrote variant_growth.tsv, variant_agreement.tsv, direction_accuracy.tsv, "
          f"core_reaction_grades.tsv under {DATA}")


if __name__ == "__main__":
    main()
