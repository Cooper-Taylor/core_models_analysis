#!/usr/bin/env python3
"""How much of each model set can each direction source actually reach?

Growth effects are bounded by a prior question: of the reactions a model set
contains, how many does a direction source have an opinion about? A source that
covers 40,000 MSDB reactions is still inert on a model whose reactions it never
mentions.

Writes results/influence_v201/coverage.json:

  * the reaction universe of each model set (distinct MSDB ids, and how often
    each appears across models)
  * per source: how many of that universe it calls, weighted by how many models
    each reaction appears in -- the quantity that predicts FBA impact
  * how many of each model's reactions a source would rewrite, and how many of
    those differ from the bounds already on disk

    python3 scripts/analyze_direction_coverage.py
"""

from __future__ import annotations

import argparse
import collections
import json
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))

from cma import directions as dirs_reg, models as models_reg, paths  # noqa: E402
from cma.directions import DirectionMap  # noqa: E402

SETS = ["v201_group_contribution", "v201_dgpredictor", "v201_equilibrator",
        "v201_llm_council", "claude_opus48", "v201_gold_only", "v201_silver_only",
        "v201_bronze_only", "v201_gold", "v201_gold_silver", "v201_gold_silver_bronze"]


def dir_from_bounds(lb, ub) -> str:
    if lb < 0 and ub > 0:
        return "="
    if ub > 0:
        return ">"
    if lb < 0:
        return "<"
    return "0"


def _scan_one(args):
    """Read one model's reactions straight from JSON -- no cobra, no solver.

    Coverage needs only reaction ids and their stored bounds, so building a
    cobra model per file costs about a second each and buys nothing. Reading the
    JSON directly and fanning out across processes turns a ~70 minute scan into
    under a minute.
    """
    import gzip as _gz
    path, loader, anno_key, strip_sfx, fallback = args
    try:
        # Detect by extension, not by loader name: a model set may declare a
        # `module:function` loader (the genome-scale sets do), in which case a
        # name test silently falls through to the text branch and every gzipped
        # model fails to parse.
        if str(path).endswith(".gz"):
            with _gz.open(path, "rt", encoding="utf-8") as fh:
                d = json.load(fh)
        else:
            d = json.loads(Path(path).read_text())
    except Exception as exc:  # noqa: BLE001
        # Return the error rather than None: a silent skip turned a total
        # parse failure into "0 models scanned", which reads like a result.
        return {"error": f"{Path(path).name}: {type(exc).__name__}: {exc}"}
    out = []
    for r in d.get("reactions", []):
        rid = r.get("id", "")
        if rid.startswith(("EX_", "SK_", "DM_", "bio")):
            continue
        raw = (r.get("annotation") or {}).get(anno_key)
        if isinstance(raw, (list, tuple)):
            raw = raw[0] if raw else None
        if not raw and fallback == "id_prefix":
            cand = re.sub(r"_[a-z]\d*$", "", rid)
            raw = cand if cand.startswith("rxn") else None
        if not raw:
            continue
        if strip_sfx:
            raw = re.sub(r"_[a-z]$", "", raw)
        out.append((raw, dir_from_bounds(r.get("lower_bound", 0), r.get("upper_bound", 0))))
    return out


def scan(model_set_name: str, limit, sample_step: int, jobs: int = 48) -> dict:
    """Reaction universe of a model set, plus the on-disk direction of each."""
    ms = models_reg.get(model_set_name)
    ids = ms.model_ids()
    if sample_step > 1:
        ids = ids[::sample_step]
    if limit:
        ids = ids[:limit]
    sk = ms.seed_key
    tasks = [(str(ms.path_for(i)), ms.loader, sk.annotation,
              sk.strip_compartment_suffix, sk.fallback) for i in ids]

    freq = collections.Counter()
    ondisk = collections.defaultdict(collections.Counter)
    per_model_rxn = []
    with ProcessPoolExecutor(max_workers=jobs) as ex:
        errors = []
        for res in ex.map(_scan_one, tasks, chunksize=16):
            if isinstance(res, dict):
                errors.append(res["error"])
                continue
            if res is None:
                continue
            seen = set()
            for seed, op in res:
                ondisk[seed][op] += 1
                seen.add(seed)
            freq.update(seen)
            per_model_rxn.append(len(seen))
    if errors:
        print(f"  WARNING: {len(errors)} model(s) unreadable, e.g. {errors[0]}")
        if len(errors) == len(tasks):
            raise SystemExit(f"every model in {model_set_name!r} failed to parse -- "
                             "refusing to report a coverage table of zeros")
    return {"model_set": ms.name, "n_models": len(per_model_rxn), "freq": freq,
            "ondisk": ondisk, "per_model_rxn": per_model_rxn}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sets", nargs="*", default=["core_kegg2", "ms2_gsm"])
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--jobs", type=int, default=48)
    ap.add_argument("--sample-step", type=int, default=1,
                    help="scan every Nth model (1 = all)")
    ap.add_argument("--out", type=Path, default=paths.results("influence_v201"))
    args = ap.parse_args()

    maps = {}
    for key in SETS:
        try:
            maps[key] = DirectionMap.load(key).ops
        except Exception as exc:  # noqa: BLE001
            print(f"  skipping {key}: {exc}")

    out = {}
    for name in args.sets:
        print(f"scanning {name} ...", flush=True)
        sc = scan(name, args.limit, args.sample_step, args.jobs)
        freq, ondisk = sc["freq"], sc["ondisk"]
        universe = set(freq)
        total_instances = sum(freq.values())
        print(f"  {sc['n_models']} models, {len(universe)} distinct MSDB reactions, "
              f"{total_instances} reaction instances")

        per_source = {}
        for key, ops in maps.items():
            hit = universe & set(ops)
            # weight by how many models each reaction appears in: this is what
            # decides whether a source moves the panel or just one model
            inst = sum(freq[r] for r in hit)
            # of those, how many differ from the direction already on disk
            changed_inst = 0
            changed_rxn = 0
            for r in hit:
                dominant = ondisk[r].most_common(1)[0][0] if ondisk[r] else "0"
                if ops[r] != dominant:
                    changed_rxn += 1
                    changed_inst += freq[r]
            per_source[key] = {
                "reactions_called": len(hit),
                "pct_of_universe": round(100.0 * len(hit) / len(universe), 2) if universe else 0,
                "reaction_instances_called": inst,
                "pct_of_instances": round(100.0 * inst / total_instances, 2) if total_instances else 0,
                "reactions_differing_from_on_disk": changed_rxn,
                "instances_differing_from_on_disk": changed_inst,
                "pct_instances_differing": round(100.0 * changed_inst / total_instances, 2)
                if total_instances else 0,
            }
        prm = sorted(sc["per_model_rxn"])
        out[name] = {
            "n_models": sc["n_models"],
            "distinct_msdb_reactions": len(universe),
            "reaction_instances": total_instances,
            "reactions_per_model": {"min": prm[0], "median": prm[len(prm) // 2],
                                    "max": prm[-1]} if prm else {},
            "on_disk_direction_mix": dict(collections.Counter(
                ondisk[r].most_common(1)[0][0] for r in universe if ondisk[r])),
            "per_source": per_source,
        }

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "coverage.json").write_text(json.dumps(out, indent=2))
    print(f"\nwrote {args.out / 'coverage.json'}")

    for name, d in out.items():
        print(f"\n=== {name}: {d['n_models']} models, "
              f"{d['distinct_msdb_reactions']} distinct reactions ===")
        print(f"{'source':<28} {'called':>7} {'%univ':>7} {'%inst':>7} {'%inst≠disk':>11}")
        for k, v in d["per_source"].items():
            print(f"{k:<28} {v['reactions_called']:>7} {v['pct_of_universe']:>7} "
                  f"{v['pct_of_instances']:>7} {v['pct_instances_differing']:>11}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
