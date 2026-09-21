"""Frozen pre-refactor literals, and the assertion that the registry matches them.

This module is scaffolding with a declared end of life. It holds a snapshot of
the literals that were scattered across the pipeline before the registry
existed, so ``cma_check.py --parity`` can prove that every projection the
registry offers is byte-identical to what the old code hardcoded. Once every
consumer reads the registry, the frozen copies here are deleted.
"""

from __future__ import annotations

FROZEN = {
    # scripts/direction_pipeline.py:58
    "PER_SOURCE_LABELS": ("Group contribution", "eQuilibrator", "dGPredictor"),
    # scripts/direction_pipeline.py:61
    "SOURCE_SLUGS": {
        "Group contribution": "group-contribution",
        "eQuilibrator": "equilibrator",
        "dGPredictor": "dgpredictor",
    },
    # scripts/run_thermo_source_variants.py:43
    "SOURCE_SPECS": [("gc", "Group contribution"), ("eq", "eQuilibrator"), ("dgp", "dGPredictor")],
    # scripts/run_thermo_source_variants.py:52
    "VARIANT_ORDER": ["kbase_baseline", "gc", "eq", "dgp"],
    # scripts/run_variant_source_panel.py:53
    "SOURCES_UPPER": ("GC", "EQ", "DGP"),
    # scripts/analyze_thermo_signature_nullspace.py:53
    "SOURCES_SHORT": {"gc": "Group contribution", "eq": "eQuilibrator", "dgp": "dGPredictor"},
    # the 14 cascade tags, in catalog order
    "CASCADE_TAGS": [
        "baseline", "3.1", "3.3", "3.3_wide", "3.5", "3.5_wide", "3.6", "3.7",
        "3.10_tight", "3.10_loose", "H4", "ri_gamma1", "ri_gamma2", "dgpredictor",
    ],
    # the six that existed only in thermo_variants/manifest.json
    "OVERLAY_TAGS": [
        "ai_opus48", "eq3_beber2022", "eq3_gamma1", "consensus_thermo",
        "kegg_implicit", "group_contribution",
    ],
    # scripts/growth_heuristics.py:29-30
    "MODELS_DIRNAME": "core_models_kegg2",
    "MEDIA_FILE_REL": "Media/KBaseMedia.cpd",
}


def per_source_labels() -> tuple:
    from . import sources

    return tuple(s.display for s in sources.all("panel") if s.label)


def source_slugs() -> dict:
    from . import sources

    return {s.display: s.slug("dash") for s in sources.all("panel") if s.label}


def source_specs() -> list:
    from . import sources

    return [(s.slug("short"), s.display) for s in sources.all("panel")]


def sources_upper() -> tuple:
    from . import sources

    return tuple(s.slug("upper") for s in sources.all("panel"))


def sources_short() -> dict:
    from . import sources

    return {s.slug("short"): s.display for s in sources.all("panel")}


def cascade_tags() -> list:
    from . import variants

    return variants.tags("cascade")


def overlay_tags() -> list:
    from . import variants

    return [v.tag for v in variants.all() if v.kind in ("overlay", "derived")]


PROJECTIONS = {
    "PER_SOURCE_LABELS": per_source_labels,
    "SOURCE_SLUGS": source_slugs,
    "SOURCE_SPECS": source_specs,
    "SOURCES_UPPER": sources_upper,
    "SOURCES_SHORT": sources_short,
    "CASCADE_TAGS": cascade_tags,
}


def assert_all() -> list[str]:
    """Return a list of mismatches; empty means the registry reproduces the literals."""
    problems = []
    for name, fn in PROJECTIONS.items():
        want, got = FROZEN[name], fn()
        if isinstance(want, list) and isinstance(got, tuple):
            got = list(got)
        if isinstance(want, tuple) and isinstance(got, list):
            got = tuple(got)
        if want != got:
            problems.append(f"{name}:\n    frozen   = {want!r}\n    registry = {got!r}")
    # the overlay tags are order-insensitive (they come from a dict-backed manifest)
    if sorted(overlay_tags()) != sorted(FROZEN["OVERLAY_TAGS"]):
        problems.append(
            f"OVERLAY_TAGS:\n    frozen   = {sorted(FROZEN['OVERLAY_TAGS'])!r}"
            f"\n    registry = {sorted(overlay_tags())!r}"
        )
    return problems
