"""Shared helpers for reading ``annotation['seed.reaction']`` from cobra model JSONs.

Background
----------
Most cobra reactions in ``core_models_kegg2/*.json`` carry a clean SEED id
in their annotation::

    "annotation": {"sbo": "SBO:0000176", "seed.reaction": "rxn00549"}

A small set of transport reactions (17 distinct SEED ids — empirically all
``is_transport=1`` in MSDB) instead carry the SEED id with a stray ``_c``
suffix::

    "annotation": {"seed.reaction": "rxn11322_c", ...}

That suffix is not part of the ModelSEED identifier — ``rxn11322`` is the
real MSDB record. Reading the annotation verbatim then comparing against
the cascade's ``{rxn_id: reversibility}`` map silently fails for these
reactions, splitting prevalence counts into a ``rxn11322`` bucket and a
``rxn11322_c`` bucket and skipping the bound override during rebound FBA.

This helper normalizes the annotation at read time so downstream code
operates on the canonical MSDB id.

See ``reports/DUPLICATE_REACTIONS_INVESTIGATION.md`` for the full
investigation, decision, and impact.
"""

from __future__ import annotations

import re
from typing import Any, Optional

# Compartment-letter suffix on a SEED id, e.g. ``_c`` in ``rxn11322_c``.
# We strip a single trailing ``_<letter>`` (with no digits — the bug
# always uses the letter-only form). We do NOT strip ``_c0`` / ``_e0``
# from cobra reaction IDs — those are the legitimate compartment marker
# on the cobra-side id and are not present in seed.reaction annotations.
_SEED_COMPARTMENT_SUFFIX = re.compile(r"_[a-z]$")


def normalize_seed_id(seed: Optional[str]) -> Optional[str]:
    """Strip a stray compartment-letter suffix from a SEED reaction id.

    >>> normalize_seed_id("rxn00549")
    'rxn00549'
    >>> normalize_seed_id("rxn11322_c")
    'rxn11322'
    >>> normalize_seed_id(None)  # returns None unchanged
    """
    if not seed:
        return seed
    return _SEED_COMPARTMENT_SUFFIX.sub("", seed)


def seed_id(reaction: Any) -> Optional[str]:
    """Return the normalized SEED id for ``reaction`` (cobra rxn or raw JSON dict).

    Accepts either a cobra ``Reaction`` object (uses ``.annotation``) or a
    raw JSON-loaded reaction dict (uses ``["annotation"]``). Returns None
    if the reaction has no SEED annotation.
    """
    if reaction is None:
        return None
    anno = getattr(reaction, "annotation", None)
    if anno is None and isinstance(reaction, dict):
        anno = reaction.get("annotation")
    if not anno:
        return None
    raw = anno.get("seed.reaction") if hasattr(anno, "get") else None
    return normalize_seed_id(raw)
