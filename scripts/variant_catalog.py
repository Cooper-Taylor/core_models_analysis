"""Canonical list of ReversibilityConfig variants -- now a view over ``cma``.

This module used to hold the variant definitions. They live in
``cma/entries/variants.py``; this file is the backwards-compatible projection,
so every existing importer keeps working unchanged::

    import variant_catalog as vc
    for v in vc.VARIANTS:        # 14 cascade variants, same order, same 7 keys
        cfg = v["cfg"]()         # still a zero-arg ReversibilityConfig factory

Why the move
------------
``VARIANTS`` could only ever describe a *cascade* variant, because ``cfg`` was a
mandatory ``ReversibilityConfig`` factory and six call sites invoke it
unconditionally. Six other variants -- ``ai_opus48``, ``consensus_thermo``,
``eq3_beber2022``, ``eq3_gamma1``, ``group_contribution``, ``kegg_implicit`` --
carry a precomputed direction map instead, so they could not be registered at
all. They lived only in ``thermo_variants/manifest.json``, which
``export_thermo_variants.py`` rebuilt from this list on every run and therefore
erased.

The registry records both, discriminated by ``Variant.kind``. ``VARIANTS`` below
deliberately projects **only** the cascade ones, so the unconditional
``v["cfg"]()`` callers cannot be handed an entry they would choke on.

New code should prefer the registry, which can see all twenty::

    from cma import variants
    variants.all()               # every variant
    variants.all("overlay")      # just the direction-map ones
    variants.get("ai_opus48")
"""

from __future__ import annotations

from cma import variants as _registry

#: The 14 cascade variants, in their historical order, as 7-key dicts.
VARIANTS: list[dict] = _registry.legacy("cascade")


def variant_by_tag(tag: str) -> dict:
    """Look up one cascade variant by tag (unchanged signature and return)."""
    return _registry.by_tag(tag)


__all__ = ["VARIANTS", "variant_by_tag"]
