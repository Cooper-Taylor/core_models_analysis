#!/usr/bin/env python3
"""Extract every reaction-direction set from ModelSEED Biochemistry v2.0.1.

    python3 scripts/build_direction_sets_v201.py
    python3 scripts/build_direction_sets_v201.py --ref v2.0.1 --out results/direction_sets_v201

Produces eleven `{reaction_id: > < =}` tables under ``--out``, in three families.

PER-SOURCE -- each source's own direction call
    v201_group_contribution   thermodynamics["Group contribution"][2]
    v201_dgpredictor          thermodynamics["dGPredictor"][2]
    v201_equilibrator         thermodynamics["eQuilibrator"][2]
    v201_llm_council          thermodynamics["LLMs"][2]
    claude_opus48             a standalone single-model run, NOT from the release

GRADE-EXCLUSIVE -- the recommended direction, one evidence tier only
    v201_gold_only            reactions graded gold
    v201_silver_only          reactions graded silver
    v201_bronze_only          reactions graded bronze

GRADE-CUMULATIVE -- the same, accumulating tiers
    v201_gold                 gold
    v201_gold_silver          gold + silver
    v201_gold_silver_bronze   gold + silver + bronze

Three things about this data decide the shape of the output.

**A grade belongs to a reaction, not to a source.** ``SourceGrading/.gitignore``
states it outright: "a reaction is graded on its collective evidence, so the
release must not also carry a per-source verdict the scheme does not make." The
per-source table is internal. So a grade-filtered set is *the recommended
direction restricted to reactions at that tier*, which is the only reading the
scheme supports.

**The recommended direction is ``reversibility``**, assigned by
``Apply_Evidence_Grades_And_Recommendation.py`` with a fixed precedence of
eQuilibrator > dGPredictor > Group contribution. Every ungraded reaction carries
``?``, so ``reversibility`` is populated exactly on the graded set.

**``?`` is dropped, never emitted.** A direction map is applied by opening a
reaction to whatever it says, and ``?`` is applied as fully reversible
(-1000, 1000). Emitting ``?`` for a reaction that was graded but could not be
called would therefore *remove* a constraint the model already had, in the name
of evidence. Absent means "no opinion", and the model keeps its own bounds. The
counts of dropped ``?`` are reported per set.

The LLM council is in the release but is deliberately outside the grading: it
reasons from the reaction rather than from an energy, carries no dG, and cannot
become the recommended direction.
"""

from __future__ import annotations

import argparse
import collections
import json
import subprocess
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))

from cma import paths  # noqa: E402

CALLABLE = ("<", ">", "=")

#: label in the release -> output slug
SOURCES = {
    "Group contribution": "group_contribution",
    "dGPredictor": "dgpredictor",
    "eQuilibrator": "equilibrator",
    "LLMs": "llm_council",
}

GRADES = ("gold", "silver", "bronze")

#: A standalone single-model Claude run. Not part of the release: the council
#: entry in v2.0.1 is an ensemble whose members include this model, so the two
#: are related but not the same thing.
CLAUDE_JSON = (paths.analysis_dir() / "data" / "ai_curation" / "all_modelseed"
               / "AICurationCacheReactionDirectionality.json")

CLAUDE_WORDS = {"forward": ">", "reverse": "<", "reversible": "=",
                "uncertain": "?", "unknown": "?"}


def load_release(ref: str) -> list:
    """Every reaction record of ``ref``, read with ``git show`` (tree untouched)."""
    repo = paths.msdb("live")
    out = []
    for i in range(61):
        shard = f"Biochemistry/reaction_{i:02d}.json"
        try:
            blob = subprocess.run(["git", "-C", str(repo), "show", f"{ref}:{shard}"],
                                  capture_output=True, check=True).stdout
        except subprocess.CalledProcessError:
            continue
        out.extend(json.loads(blob))
    if not out:
        raise SystemExit(f"no reaction shards found at {ref} in {repo}")
    return out


def load_claude() -> dict:
    if not CLAUDE_JSON.exists():
        print(f"  claude: {CLAUDE_JSON} absent -- skipping")
        return {}
    raw = json.loads(CLAUDE_JSON.read_text())
    out = {}
    for rxn, val in raw.items():
        if isinstance(val, dict):
            val = (val.get("direction") or val.get("directionality")
                   or val.get("value") or "")
        v = str(val).strip().lower()
        op = CLAUDE_WORDS.get(v, val if val in CALLABLE else "?")
        if op in CALLABLE:
            out[rxn.split("_")[0]] = op
    return out


def build(records: list) -> dict:
    """Every direction set, as ``{set_name: {rxn: op}}``, plus a stats table."""
    sets: dict = {}
    stats: dict = {}

    # --- per source ------------------------------------------------------
    for label, slug in SOURCES.items():
        ops, dropped = {}, 0
        for r in records:
            t = r.get("thermodynamics")
            if not isinstance(t, dict) or label not in t:
                continue
            v = t[label]
            op = str(v[2]).strip() if len(v) > 2 else ""
            if op in CALLABLE:
                ops[r["id"]] = op
            else:
                dropped += 1
        sets[f"v201_{slug}"] = ops
        stats[f"v201_{slug}"] = {"family": "source", "source": label,
                                 "n_callable": len(ops), "n_dropped_uncallable": dropped}

    # --- grade families --------------------------------------------------
    by_grade: dict = {g: {} for g in GRADES}
    dropped_by_grade = collections.Counter()
    for r in records:
        te = r.get("thermo-evidence")
        if not isinstance(te, dict):
            continue
        g = (te.get("grade") or "").lower()
        if g not in by_grade:
            continue
        op = str(r.get("reversibility") or "").strip()
        if op in CALLABLE:
            by_grade[g][r["id"]] = op
        else:
            dropped_by_grade[g] += 1

    for g in GRADES:
        sets[f"v201_{g}_only"] = dict(by_grade[g])
        stats[f"v201_{g}_only"] = {
            "family": "grade-exclusive", "grades": [g],
            "n_callable": len(by_grade[g]),
            "n_dropped_uncallable": dropped_by_grade[g],
        }

    cumulative: dict = {}
    for i, g in enumerate(GRADES):
        cumulative.update(by_grade[g])
        name = "v201_" + "_".join(GRADES[: i + 1])
        sets[name] = dict(cumulative)
        stats[name] = {
            "family": "grade-cumulative", "grades": list(GRADES[: i + 1]),
            "n_callable": len(cumulative),
            "n_dropped_uncallable": sum(dropped_by_grade[x] for x in GRADES[: i + 1]),
        }

    return sets, stats


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", default="v2.0.1", help="git ref of the release")
    ap.add_argument("--out", type=Path,
                    default=paths.results("direction_sets_v201"))
    args = ap.parse_args()

    print(f"reading {args.ref} from {paths.msdb('live')} ...", flush=True)
    records = load_release(args.ref)
    print(f"  {len(records)} reactions")

    sets, stats = build(records)

    claude = load_claude()
    if claude:
        sets["claude_opus48"] = claude
        stats["claude_opus48"] = {"family": "source", "source": "Claude Opus 4.8",
                                  "n_callable": len(claude), "n_dropped_uncallable": 0,
                                  "note": "standalone single-model run, not part of the release"}

    args.out.mkdir(parents=True, exist_ok=True)
    order = [f"v201_{s}" for s in SOURCES.values()] + ["claude_opus48"] \
        + [f"v201_{g}_only" for g in GRADES] \
        + ["v201_gold", "v201_gold_silver", "v201_gold_silver_bronze"]

    print(f"\n{'set':<28} {'family':<18} {'calls':>7} {'>':>6} {'<':>6} {'=':>6} {'dropped ?':>10}")
    rows = []
    for name in order:
        ops = sets.get(name)
        if ops is None:
            continue
        (args.out / f"{name}.json").write_text(
            json.dumps(dict(sorted(ops.items())), separators=(",", ":")))
        c = collections.Counter(ops.values())
        st = stats[name]
        print(f"{name:<28} {st['family']:<18} {len(ops):>7} {c['>']:>6} {c['<']:>6} "
              f"{c['=']:>6} {st['n_dropped_uncallable']:>10}")
        rows.append({"name": name, **st, "operators": dict(c),
                     "path": str((args.out / f"{name}.json").relative_to(paths.analysis_dir()))})

    manifest = {"release": args.ref, "n_reactions_in_release": len(records), "sets": rows}
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"\nwrote {len(rows)} sets + manifest.json to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
