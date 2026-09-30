#!/usr/bin/env python3
"""Which reactions a direction tier would modify, split by gapfill provenance.

    python3 scripts/analyze_gapfill_direction_impact.py
    python3 scripts/analyze_gapfill_direction_impact.py --limit 200   # smoke test

The ModelSEED v2 genome-scale models exist twice: the same 5,419 genomes
reconstructed once and gap-filled against glucose minimal media and again
against auxotrophy media. Comparing the pair separates what the genome
supports from what gap-filling added.

    original   present in BOTH media models and gap-filled in neither.
               The genome-derived core: what annotation alone supports.
    gapfilled  carries gapfill_data in either model.
    media-only present in one media model and absent from the other. Not
               flagged as gap-filled, but conditional on the medium, so it is
               reported separately rather than folded into either class.

For each class the script reports how many reactions a direction tier
(gold / silver / bronze, exclusive and cumulative) would call differently from
the bounds the model file already carries, at both reaction and occurrence
level, and classifies the modified reactions by the cofactors they use.

Outputs results/gapfill_direction_impact/{summary,by_cofactor,reactions}.{json,tsv}.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
from cma import paths  # noqa: E402
from cma.directions import DirectionMap  # noqa: E402

GS = Path(os.environ.get("MS2_GS_MODELS_DIR", "/scratch/ctaylor/modelseed2_gs_models"))

TIERS = ["v201_gold_only", "v201_silver_only", "v201_bronze_only",
         "v201_gold", "v201_gold_silver", "v201_gold_silver_bronze"]
TIER_LABEL = {"v201_gold_only": "gold", "v201_silver_only": "silver",
              "v201_bronze_only": "bronze", "v201_gold": "gold (cumulative)",
              "v201_gold_silver": "gold + silver",
              "v201_gold_silver_bronze": "gold + silver + bronze"}

#: Cofactor families, by ModelSEED compound id. A reaction is tagged with every
#: family it touches, so the tags overlap by design: an ATP-dependent
#: dehydrogenase is both "ATP/ADP" and "NAD(H)".
COFACTORS = {
    "ATP/ADP/AMP": {"cpd00002", "cpd00008", "cpd00018"},
    "NAD(H)": {"cpd00003", "cpd00004"},
    "NADP(H)": {"cpd00006", "cpd00005"},
    "CoA": {"cpd00010"},
    "FAD(H2)": {"cpd00015", "cpd00982"},
    "phosphate": {"cpd00009"},
    "pyrophosphate": {"cpd00012"},
    "quinone": {"cpd15560", "cpd15561", "cpd15499", "cpd15500", "cpd29674", "cpd15352"},
    "THF / folate": {"cpd00087", "cpd00016", "cpd00345", "cpd00141"},
    "SAM": {"cpd00017", "cpd00019"},
    "glutathione": {"cpd00042", "cpd00111"},
    "ferredoxin": {"cpd11621", "cpd11620"},
    "CO2": {"cpd00011"},
    "ammonia": {"cpd00013"},
    "O2": {"cpd00007"},
}
#: Present in almost every reaction; tagging on them would say nothing.
UBIQUITOUS = {"cpd00001", "cpd00067"}


def dir_from_bounds(lb, ub) -> str:
    if lb < 0 < ub:
        return "="
    if ub > 0:
        return ">"
    if lb < 0:
        return "<"
    return "0"


# ---------------------------------------------------------------------------
def _scan_pair(genome: str) -> tuple | None:
    """Classify every reaction of one genome across its two media models."""
    out = {}
    for slug, tag in (("gmm", "GMM"), ("auxotrophy", "auxo")):
        p = GS / slug / "models" / f"{genome}.RAST.{tag}.mdl.json.gz"
        if not p.exists():
            return None
        try:
            with gzip.open(p, "rt", encoding="utf-8") as fh:
                d = json.load(fh)
        except Exception:
            return None
        m = {}
        for r in d.get("modelreactions", []):
            base = r["id"].rsplit("_", 1)[0]
            if not base.startswith("rxn"):
                continue
            direction = r.get("direction", "=")
            maxf = float(r.get("maxforflux", 1000) or 0)
            maxr = float(r.get("maxrevflux", 1000) or 0)
            lb, ub = ((0.0, maxf) if direction == ">" else
                      (-maxr, 0.0) if direction == "<" else (-maxr, maxf))
            m[base] = (bool(r.get("gapfill_data")), dir_from_bounds(lb, ub))
        out[slug] = m

    g, a = out["gmm"], out["auxotrophy"]
    shared = set(g) & set(a)
    rows = []
    for rxn in shared:
        gf = g[rxn][0] or a[rxn][0]
        rows.append((rxn, "gapfilled" if gf else "original", g[rxn][1]))
    for rxn in set(g) - shared:
        rows.append((rxn, "media-only", g[rxn][1]))
    for rxn in set(a) - shared:
        rows.append((rxn, "media-only", a[rxn][1]))
    return genome, rows


def load_stoichiometry(ref: str) -> dict:
    """{rxn_id: set(compound ids)} from a ModelSEED release."""
    repo = paths.msdb("live")
    out = {}
    for i in range(61):
        blob = subprocess.run(
            ["git", "-C", str(repo), "show", f"{ref}:Biochemistry/reaction_{i:02d}.json"],
            capture_output=True).stdout
        if not blob:
            continue
        for r in json.loads(blob):
            # v2.0.1 stores stoichiometry as a list of dicts; older releases
            # used a ';'-delimited string, and the equation is the fallback.
            cpds = set()
            st = r.get("stoichiometry")
            if isinstance(st, list):
                cpds = {t.get("compound") for t in st if isinstance(t, dict)
                        and t.get("compound")}
            elif isinstance(st, str):
                cpds = set(re.findall(r"cpd\d{5}", st))
            if not cpds:
                cpds = set(re.findall(r"cpd\d{5}", r.get("equation") or ""))
            out[r["id"]] = cpds
    return out


def cofactors_of(cpds: set) -> list:
    real = cpds - UBIQUITOUS
    tags = [name for name, ids in COFACTORS.items() if real & ids]
    return tags or ["(none of the tracked families)"]


# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--limit", type=int, default=None, help="cap paired genomes")
    ap.add_argument("--jobs", type=int, default=48)
    ap.add_argument("--ref", default="v2.0.1")
    ap.add_argument("--out", type=Path,
                    default=paths.results("gapfill_direction_impact"))
    args = ap.parse_args()

    gmm = {p.name.replace(".RAST.GMM.mdl.json.gz", "")
           for p in (GS / "gmm" / "models").glob("*.json.gz")}
    aux = {p.name.replace(".RAST.auxo.mdl.json.gz", "")
           for p in (GS / "auxotrophy" / "models").glob("*.json.gz")}
    genomes = sorted(gmm & aux)
    if args.limit:
        genomes = genomes[: args.limit]
    print(f"paired genomes: {len(genomes):,}"
          f"  (unpaired: {len(gmm ^ aux)})", flush=True)

    # class -> rxn -> how many genomes carry it that way
    freq = {c: collections.Counter() for c in ("original", "gapfilled", "media-only")}
    # rxn -> on-disk direction -> occurrences (the direction is near-constant per rxn)
    ondisk = collections.defaultdict(collections.Counter)
    n_done = 0
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        for res in ex.map(_scan_pair, genomes, chunksize=16):
            if res is None:
                continue
            n_done += 1
            for rxn, cls, direction in res[1]:
                freq[cls][rxn] += 1
                ondisk[rxn][direction] += 1
    print(f"scanned {n_done:,} genome pairs", flush=True)

    # A reaction is "original" if it is ever genome-derived; gapfilled only if
    # it is never genome-derived in any model. Otherwise one gap-filled instance
    # anywhere would reclassify a reaction the genome plainly encodes.
    klass = {}
    for rxn in set(freq["original"]) | set(freq["gapfilled"]) | set(freq["media-only"]):
        if freq["original"][rxn]:
            klass[rxn] = "original"
        elif freq["gapfilled"][rxn]:
            klass[rxn] = "gapfilled"
        else:
            klass[rxn] = "media-only"

    maps = {}
    for t in TIERS:
        try:
            maps[t] = DirectionMap.load(t).ops
        except Exception as exc:  # noqa: BLE001
            print(f"  skipping {t}: {exc}")
    stoich = load_stoichiometry(args.ref)

    def dominant(rxn):
        c = ondisk[rxn]
        return c.most_common(1)[0][0] if c else "?"

    summary, per_rxn = [], []
    cof = collections.defaultdict(collections.Counter)   # modified
    bg = collections.defaultdict(collections.Counter)    # all called
    n_called = {}
    for t, ops in maps.items():
        for cls in ("original", "gapfilled", "media-only"):
            rxns = [r for r in klass if klass[r] == cls]
            called = [r for r in rxns if r in ops]
            modified = [r for r in called if ops[r] != dominant(r)]
            inst_all = sum(freq[cls][r] for r in rxns)
            inst_mod = sum(freq[cls][r] for r in modified)
            summary.append({
                "tier": TIER_LABEL[t], "tier_key": t, "class": cls,
                "reactions_in_class": len(rxns),
                "reactions_called_by_tier": len(called),
                "reactions_modified": len(modified),
                "pct_of_class_modified": round(100 * len(modified) / len(rxns), 2) if rxns else 0,
                "pct_of_called_modified": round(100 * len(modified) / len(called), 2) if called else 0,
                "occurrences_in_class": inst_all,
                "occurrences_modified": inst_mod,
                "pct_occurrences_modified": round(100 * inst_mod / inst_all, 2) if inst_all else 0,
            })
            if t == "v201_gold_silver_bronze":
                n_called[cls] = (len(called), len(modified))
                # Enrichment, not raw share. "29% of modified reactions use ATP"
                # means nothing on its own if 29% of all reactions use ATP; what
                # matters is the rate among reactions the tier CALLED, modified
                # against not-modified.
                modset = set(modified)
                for r in called:
                    for fam in cofactors_of(stoich.get(r, set())):
                        bg[cls][fam] += 1
                        if r in modset:
                            cof[cls][fam] += 1
                    per_rxn.append({
                        "rxn": r, "class": cls, "on_disk": dominant(r),
                        "tier_call": ops[r], "n_genomes": freq[cls][r],
                        "cofactors": ";".join(cofactors_of(stoich.get(r, set()))),
                    })

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2))
    with (args.out / "summary.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0]), delimiter="\t")
        w.writeheader(); w.writerows(summary)
    with (args.out / "reactions.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["rxn", "class", "on_disk", "tier_call",
                                           "n_genomes", "cofactors"], delimiter="\t")
        w.writeheader(); w.writerows(sorted(per_rxn, key=lambda x: -x["n_genomes"]))
    enrich = {}
    for cls in cof:
        ncall, nmod = n_called.get(cls, (0, 0))
        base = nmod / ncall if ncall else 0
        rows = []
        for fam, n_all in bg[cls].most_common():
            n_mod = cof[cls][fam]
            rate = n_mod / n_all if n_all else 0
            rows.append({"cofactor": fam, "called": n_all, "modified": n_mod,
                         "modified_rate": round(100 * rate, 1),
                         "enrichment_vs_class": round(rate / base, 2) if base else None})
        enrich[cls] = {"class_modified_rate": round(100 * base, 1),
                       "n_called": ncall, "n_modified": nmod, "families": rows}
    (args.out / "by_cofactor.json").write_text(json.dumps(enrich, indent=2))

    # ---- console report ----
    print(f"\n=== reaction provenance across {n_done:,} paired genomes ===")
    for cls in ("original", "gapfilled", "media-only"):
        n = sum(1 for r in klass.values() if r == cls)
        inst = sum(freq[cls][r] for r in klass if klass[r] == cls)
        print(f"  {cls:<11} {n:>6,} distinct reactions   {inst:>10,} occurrences")

    print(f"\n=== how many would a tier call differently from the model file? ===")
    print(f"{'tier':<24} {'class':<11} {'in class':>9} {'called':>8} "
          f"{'modified':>9} {'% class':>8} {'% called':>9}")
    for s in summary:
        print(f"{s['tier']:<24} {s['class']:<11} {s['reactions_in_class']:>9,} "
              f"{s['reactions_called_by_tier']:>8,} {s['reactions_modified']:>9,} "
              f"{s['pct_of_class_modified']:>7.1f}% {s['pct_of_called_modified']:>8.1f}%")

    print(f"\n=== cofactor families of modified reactions (gold+silver+bronze) ===")
    print("    enrichment > 1 means the tier changes that family more often than")
    print("    it changes the class as a whole\n")
    for cls in ("original", "gapfilled"):
        e = enrich.get(cls)
        if not e:
            continue
        print(f"  {cls}: {e['n_modified']:,} of {e['n_called']:,} called reactions "
              f"modified ({e['class_modified_rate']}%)")
        print(f"     {'cofactor':<32} {'called':>7} {'modified':>9} {'rate':>7} {'enrich':>7}")
        for r in sorted(e["families"], key=lambda x: -(x["enrichment_vs_class"] or 0)):
            if r["called"] < 15:
                continue
            print(f"     {r['cofactor']:<32} {r['called']:>7,} {r['modified']:>9,} "
                  f"{r['modified_rate']:>6.1f}% {r['enrichment_vs_class']:>7.2f}")
        print()
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
