"""The thermodynamic sources.

Transcribes today's reality exactly: the four ΔG′° labels that live in MSDB's
per-reaction ``thermodynamics`` dict, plus their spellings as actually used.
Counts measured over all 56,012 reactions of the live clone:

    Group contribution      56,002
    dGPredictor-ModelSEED   31,924
    dGPredictor             27,715
    eQuilibrator            25,028

``dGPredictor-ModelSEED`` is present in the data but was never a first-class
source in any script -- it is registered here so promoting it is a one-line
`layers` change rather than 20 literal edits.
"""

from __future__ import annotations

from ..kinds import ThermoSource
from ..sources import register

register(ThermoSource(
    key="gc",
    display="Group contribution",
    label="Group contribution",
    slugs={"short": "gc", "upper": "GC", "long": "group_contribution",
           "dash": "group-contribution"},
    kind="msdb_sublist",
    db_level="GC",
    color="#2a78d6",
    color_families=("bars", "scatter"),
    axis_title="Group contribution ΔG′° (kcal/mol)",
    section="§ 2.2",
    citations=("Jankowski2008",),
))

register(ThermoSource(
    key="eq",
    display="eQuilibrator",
    label="eQuilibrator",
    slugs={"short": "eq", "upper": "EQ", "long": "equilibrator", "dash": "equilibrator"},
    kind="msdb_sublist",
    db_level="EQ",
    # The EQ level falls back to group contribution when a reaction is
    # incomplete -- upstream behaviour, recorded as data rather than as an
    # `if` buried in the cascade.
    incomplete_fallback="GC",
    color="#e8833a",
    color_families=("bars", "scatter"),
    axis_title="eQuilibrator ΔG′° (kcal/mol)",
    section="§ 2.3",
    citations=("Beber2022", "Noor2013"),
))

register(ThermoSource(
    key="dgp",
    display="dGPredictor",
    label="dGPredictor",
    slugs={"short": "dgp", "upper": "DGP", "long": "dgpredictor", "dash": "dgpredictor"},
    kind="msdb_sublist",
    db_level="DGP",
    color="#3aa76d",
    color_families=("bars", "scatter"),
    axis_title="dGPredictor ΔG′° (kcal/mol)",
    section="§ New — dGPredictor energies",
    citations=("Wang2021",),
))

register(ThermoSource(
    key="dgpms",
    display="dGPredictor-ModelSEED",
    label="dGPredictor-ModelSEED",
    slugs={"short": "dgpms", "upper": "DGPMS", "long": "dgpredictor_modelseed",
           "dash": "dgpredictor-modelseed"},
    kind="msdb_sublist",
    db_level=None,
    # Present in MSDB (31,924 reactions) but not yet promoted into the panel /
    # all-models sweeps. Add "panel" here to promote it.
    layers=(),
    color="#9a6fb0",
    color_families=("bars",),
    axis_title="dGPredictor-ModelSEED ΔG′° (kcal/mol)",
    section="§ New — retrained dGPredictor",
    citations=("Freiburger2025",),
))
