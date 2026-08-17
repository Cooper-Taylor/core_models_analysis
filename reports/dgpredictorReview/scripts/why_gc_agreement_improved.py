#!/usr/bin/env python3
"""Why does the fine-tuned dGPredictor correlate so much better with Group
Contribution than the original does?

Pooled on ModelSEED dev, r(GC, dGPredictor) = 0.22 and r(GC, dGPredictor-ModelSEED)
= 0.80. This script decomposes that gap into candidate causes and rules three of
them out.

Candidates tested
  coverage      the two variants cover different reaction sets, so the pooled r's
                are computed on different populations
  mis-mapping   the ORIGINAL predictor's stored value is often another reaction's
                dG (a KEGG id ModelSEED does not list as an alias)
  zero outputs  the retrain returns exactly 0.00 for 3,673 reactions, which might
                pile onto the origin and manufacture agreement
  leverage      Pearson r on these sources is tail-dominated; both GC and the
                retrain are strictly additive over the SAME ModelSEED structures,
                so both scale with molecule size and agree trivially on big
                reactions

Everything is reported as r AND Spearman rho AND median |delta|, because the
three disagree here and only quoting r would give the wrong answer.

Reads results/_devsnap2_thermo.pkl (built by the snippet in REVIEW.md, or
regenerate: one row per non-EMPTY dev reaction with dg_/sig_ per source, eQ
sentinels sigma > 100 dropped).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
sys.path.insert(0, "/scratch/ctaylor/core_models_analysis/scripts")
from build_dgpredictor_kegg_mask import load_mask  # noqa: E402

CUT = 1500.0        # same implausible-magnitude cutoff as the scatter scripts
PRED = Path("/scratch/ctaylor/dgpredictor_repo/data"
            "/modelseed_all_reaction_dG_predictions.json")


def stats(sub: pd.DataFrame, a: str, b: str) -> dict | None:
    x = sub[f"dg_{a}"].to_numpy(float)
    y = sub[f"dg_{b}"].to_numpy(float)
    if len(x) < 3:
        return None
    return {
        "n": int(len(x)),
        "r": round(float(np.corrcoef(x, y)[0, 1]), 3),
        "rho": round(float(pd.Series(x).corr(pd.Series(y), method="spearman")), 3),
        "median_abs_delta": round(float(np.median(np.abs(y - x))), 2),
        "frac_sign_flip": round(float(np.mean(np.sign(x) * np.sign(y) < 0)), 3),
    }


def main() -> None:
    df = pd.read_pickle(RES / "_devsnap2_thermo.pkl")
    mask = load_mask()
    zero = {k for k, v in json.loads(PRED.read_text()).items()
            if abs(v["dG_model_only"]) < 1e-9}

    def pair(a: str, b: str, idx=None) -> pd.DataFrame:
        m = df[f"dg_{a}"].notna() & df[f"dg_{b}"].notna()
        s = df[m]
        s = s[(s[f"dg_{a}"].abs() <= CUT) & (s[f"dg_{b}"].abs() <= CUT)]
        return s if idx is None else s[s.index.isin(idx)]

    out: dict = {}

    # -- the claim, on each variant's own coverage --------------------------
    out["own_coverage"] = {f"{ref}_vs_{v}": stats(pair(ref, v), ref, v)
                           for ref in ("GC", "EQ") for v in ("BASE", "FT")}

    # -- rule out coverage: identical reaction sets --------------------------
    common = {ref: set(pair(ref, "BASE").index) & set(pair(ref, "FT").index)
              for ref in ("GC", "EQ")}
    out["like_for_like"] = {
        f"{ref}_vs_{v}": stats(pair(ref, v, common[ref]), ref, v)
        for ref in ("GC", "EQ") for v in ("BASE", "FT")}

    # -- the mis-mapping split ----------------------------------------------
    cm = common["GC"]
    split = {}
    for tag, keep in (("vouched", lambda s: ~s.index.isin(mask)),
                      ("mismapped", lambda s: s.index.isin(mask))):
        for v in ("BASE", "FT"):
            s = pair("GC", v, cm)
            split[f"{tag}_GC_vs_{v}"] = stats(s[keep(s)], "GC", v)
    # the mechanism: a near-constant series cannot correlate with anything
    for tag, keep in (("vouched", ~pair("GC", "BASE", cm).index.isin(mask)),
                      ("mismapped", pair("GC", "BASE", cm).index.isin(mask))):
        s = pair("GC", "BASE", cm)[keep]
        split[f"{tag}_value_reuse"] = {
            "n": int(len(s)),
            "distinct_BASE_values": int(s["dg_BASE"].nunique()),
            "distinct_BASE_frac": round(float(s["dg_BASE"].nunique() / len(s)), 4),
            "distinct_FT_values": int(df.loc[s.index, "dg_FT"].nunique()),
            "median_abs_GC": round(float(s["dg_GC"].abs().median()), 2),
        }
    out["kegg_mismapping_split"] = split

    # -- rule out the retrain's zero outputs ---------------------------------
    z = {}
    for v in ("BASE", "FT"):
        s = pair("GC", v, cm)
        z[f"all_GC_vs_{v}"] = stats(s, "GC", v)
        z[f"drop_FT_zero_GC_vs_{v}"] = stats(s[~s.index.isin(zero)], "GC", v)
    out["zero_output_control"] = z

    # -- leverage ladder ------------------------------------------------------
    lev, lev_v = {}, {}
    for lim in (None, 300, 100, 50, 20):
        for v in ("BASE", "FT"):
            s = pair("GC", v, cm)
            sv = s[~s.index.isin(mask)]
            if lim:
                s = s[s["dg_GC"].abs() < lim]
                sv = sv[sv["dg_GC"].abs() < lim]
            lev[f"absGC_lt_{lim}_GC_vs_{v}"] = stats(s, "GC", v)
            lev_v[f"absGC_lt_{lim}_GC_vs_{v}"] = stats(sv, "GC", v)
    out["leverage_ladder_all"] = lev
    out["leverage_ladder_vouched_only"] = lev_v

    (RES / "why_gc_agreement_improved.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
