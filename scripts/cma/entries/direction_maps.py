"""Direction maps: the ``{rxn_id: '>'|'<'|'='|'?'}`` tables.

Loaders are lazy -- a spec registers a *callable*, never a loaded map, so
importing the registry stays cheap and a missing artifact fails where it is
used rather than at import.

``coverage`` is declared. ``complete`` means the map has an opinion about every
MSDB reaction (possibly ``?``); ``partial`` means an absent key is *no opinion*.
The distinction changes published growth numbers substantially, so it is
carried rather than normalized away.
"""

from __future__ import annotations

from pathlib import Path

from .. import paths
from ..directions import register
from ..kinds import DirectionMapSpec


def _json_map(path: Path):
    def loader():
        import json

        return json.loads(Path(path).read_text())

    return loader


def _csv_map(path: Path):
    def loader():
        import csv

        out = {}
        with Path(path).open() as fh:
            for row in csv.DictReader(fh):
                rid = row.get("rxn_id") or row.get("id") or row.get("reaction_id")
                op = row.get("reversibility") or row.get("operator") or row.get("direction")
                if rid and op:
                    out[rid] = op.strip()
        return out

    return loader


_R = paths.results

# --- per-thermo-source maps, produced by build_thermo_source_direction_maps --
for _key, _slug, _src in (
    ("gc", "group-contribution", "gc"),
    ("eq", "equilibrator", "eq"),
    ("dgp", "dgpredictor", "dgp"),
):
    register(DirectionMapSpec(
        key=f"thermo_{_key}",
        loader=_json_map(_R(f"rxn_directions_{_slug}.json")),
        coverage="partial",
        source=_src,
        description=f"Per-reaction operators from the {_src} ΔG′° source.",
        artifacts={"json": _R(f"rxn_directions_{_slug}.json"),
                   "csv": _R(f"rxn_directions_{_slug}.csv")},
    ))

# --- MSDB branch snapshots --------------------------------------------------
register(DirectionMapSpec(
    key="msdb_dev",
    loader=_json_map(_R("rev_map_dev.json")),
    coverage="complete",
    description="MSDB dev branch stored reversibility.",
    artifacts={"json": _R("rev_map_dev.json")},
))
register(DirectionMapSpec(
    key="msdb_claude",
    loader=_json_map(_R("rev_map_claude.json")),
    coverage="complete",
    description="MSDB claude-changes branch stored reversibility.",
    artifacts={"json": _R("rev_map_claude.json")},
))

# --- the live cascade -------------------------------------------------------
register(DirectionMapSpec(
    key="cascade_live",
    loader=_csv_map(_R("rxn_directions_cascade_live.csv")),
    coverage="complete",
    description="Directions from a live run of the ported cascade.",
    artifacts={"csv": _R("rxn_directions_cascade_live.csv"),
               "json": _R("rxn_directions_cascade_live.json")},
))

register(DirectionMapSpec(
    key="dgpredictor_energies",
    loader=_json_map(_R("rxn_directions_dgpredictor.json")),
    coverage="partial",
    source="dgp",
    description="Cascade run with dGPredictor energies substituted.",
    artifacts={"json": _R("rxn_directions_dgpredictor.json")},
))


# --- overlay variants: directions parsed from their MSDB-format reports -----
# Each variant under thermo_variants/<tag>/ writes the upstream report shape
#   rxn_id \t status \t old_rev \t new_rev
# and the *new_rev* column is the map. These six were produced by standalone
# builder scripts and had no registry entry of any kind before.
def _report_map(tag: str, level: str = "EQ"):
    suffix = f"_{level}" if level else ""

    def loader():
        p = paths.thermo_variants(tag, f"Estimated_Reaction_Reversibility_Report{suffix}.txt")
        out = {}
        for line in p.read_text().splitlines():
            if not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) >= 4:
                out[cols[0]] = cols[3].strip()
            elif len(cols) == 3:          # the GC report drops the old_rev column
                out[cols[0]] = cols[2].strip()
        return out

    return loader


for _tag, _src, _desc in (
    ("ai_opus48", None, "Per-reaction directions predicted by Claude Opus 4.8."),
    ("eq3_beber2022", "eq", "eQuilibrator 3.0 (Beber 2022) directions."),
    ("eq3_gamma1", "eq", "eQuilibrator 3.0 at a 10-fold reversibility-index cutoff."),
    ("kegg_implicit", None, "Directions implied by the on-disk KEGG model bounds."),
    ("group_contribution", "gc", "Group-contribution source directions, as exported."),
    ("consensus_thermo", None, "Majority vote across the thermodynamic sources."),
):
    register(DirectionMapSpec(
        key=_tag,
        loader=_report_map(_tag),
        coverage="complete",
        source=_src,
        description=_desc,
        artifacts={"report_eq": paths.thermo_variants(_tag, "Estimated_Reaction_Reversibility_Report_EQ.txt")},
    ))


# ---------------------------------------------------------------------------
# ModelSEED Biochemistry v2.0.1 (the 2026 update, tagged 2026-09-15)
# ---------------------------------------------------------------------------
# Built by scripts/build_direction_sets_v201.py, which reads the release with
# `git show v2.0.1:...` and never touches the MSDB working tree.
#
# Three families:
#   source            each estimator's own direction call
#   grade-exclusive   the recommended direction, one evidence tier only
#   grade-cumulative  the same, accumulating tiers
#
# A grade belongs to a REACTION, not to a source -- the release deliberately
# ships one grade per reaction and keeps the per-source table internal. So a
# grade set is the recommended direction (precedence eQuilibrator > dGPredictor
# > Group contribution) restricted to reactions at that tier.
#
# None of these maps contains '?'. A direction map is applied by rewriting
# bounds, and '?' is applied as fully reversible, so emitting it for a reaction
# that was graded but could not be called would REMOVE a constraint in the name
# of evidence. Absent means no opinion, and the model keeps its own bounds.
_V201 = paths.results("direction_sets_v201")

for _key, _coverage, _src, _desc in (
    # --- per source ---
    ("v201_group_contribution", "partial", "gc",
     "v2.0.1 Group contribution direction calls (19,246 callable)."),
    ("v201_dgpredictor", "partial", "dgp",
     "v2.0.1 dGPredictor direction calls (12,177 callable)."),
    ("v201_equilibrator", "partial", "eq",
     "v2.0.1 eQuilibrator direction calls (21,218 callable)."),
    ("v201_llm_council", "partial", None,
     "v2.0.1 LLM council: three models predict, a fourth audits, a fifth "
     "adjudicates. Carries no dG and is deliberately excluded from the "
     "evidence grading (44,850 callable, 87% forward)."),
    ("claude_opus48", "partial", None,
     "Claude Opus 4.8 single-model direction calls (51,515). A standalone run, "
     "not part of the release; the council's members include this model, so "
     "the two are related but not the same."),
    # --- grade-exclusive ---
    ("v201_gold_only", "partial", None,
     "Recommended direction for gold-graded reactions only (3,089 callable "
     "of 3,434 graded)."),
    ("v201_silver_only", "partial", None,
     "Recommended direction for silver-graded reactions only (16,150 of 18,388)."),
    ("v201_bronze_only", "partial", None,
     "Recommended direction for bronze-graded reactions only (8,875 of 11,277)."),
    # --- grade-cumulative ---
    ("v201_gold", "partial", None,
     "Cumulative tier 1: gold (3,089 callable)."),
    ("v201_gold_silver", "partial", None,
     "Cumulative tier 2: gold + silver (19,239 callable)."),
    ("v201_gold_silver_bronze", "partial", None,
     "Cumulative tier 3: gold + silver + bronze (28,114 callable). Equals the "
     "release's own `reversibility` field wherever it is callable."),
):
    register(DirectionMapSpec(
        key=_key,
        loader=_json_map(_V201 / f"{_key}.json"),
        coverage=_coverage,
        source=_src,
        data_ref="msdb:v2.0.1",
        description=_desc,
        artifacts={"json": _V201 / f"{_key}.json"},
    ))
