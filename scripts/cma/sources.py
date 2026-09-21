"""Thermodynamic source registry.

Replaces ~30 hand-maintained literal lists. Before this module the same source
appeared under six mutually incompatible spellings with no mapping table, so a
partial edit silently dropped a source from one artifact and nothing failed.
"""

from __future__ import annotations

from .kinds import ThermoSource
from .registry import Registry, ensure_loaded

_reg: Registry[ThermoSource] = Registry("thermo source")

register = _reg.register


def get(key: str) -> ThermoSource:
    ensure_loaded()
    return _reg.get(key)


def all(layer: str | None = None) -> list[ThermoSource]:
    """Every registered source, optionally restricted to one layer."""
    ensure_loaded()
    if layer is None:
        return _reg.all()
    return _reg.where(lambda s: layer in s.layers)


def keys() -> list[str]:
    ensure_loaded()
    return _reg.keys()


def resolve(name: str) -> ThermoSource:
    """Look up a source by its key or by ANY of its declared spellings.

    This is the function that lets old code keep passing ``'gc'``,
    ``'Group contribution'`` or ``'group-contribution'`` interchangeably.
    """
    ensure_loaded()
    hit = _reg.maybe(name)
    if hit is not None:
        return hit
    for s in _reg.all():
        if name == s.display or name == s.label or name in s.slugs.values():
            return s
    raise KeyError(
        f"no thermo source spelled {name!r}. Known keys: {_reg.keys()}; "
        "add a spelling to that source's `slugs` rather than inventing one."
    )


def slug(name: str, style: str = "long") -> str:
    return resolve(name).slug(style)


def labels() -> dict[str, str]:
    """``{key: exact MSDB thermodynamics label}`` for MSDB-resident sources."""
    return {s.key: s.label for s in all() if s.label}


def display_map(style: str = "long") -> dict[str, str]:
    return {s.slug(style): s.display for s in all()}


def palette(family: str = "bars") -> dict[str, str]:
    """``{key: colour}`` restricted to sources that opt into this family.

    Slot caps are data: a scatter palette that only has three distinguishable
    inks says so, instead of silently reusing one when a fourth source appears.
    """
    return {s.key: s.color for s in all() if s.color and family in s.color_families}


def pairs(style: str = "long", layer: str | None = None) -> list[tuple[str, str]]:
    """Every unordered pair of sources, in registration order.

    Three pairs today; six the moment a fourth source registers, with no edit
    to any plotting script.
    """
    ss = [s.slug(style) for s in all(layer)]
    return [(a, b) for i, a in enumerate(ss) for b in ss[i + 1:]]


def operators(key: str) -> dict:
    """Resolve a source's ``{rxn_id: operator}`` map, applying its mask."""
    s = resolve(key)
    if s.operators is None:
        raise ValueError(f"thermo source {s.key!r} declares no operators loader")
    ops = s.operators()
    if s.mask is not None:
        masked = s.mask()
        ops = {k: v for k, v in ops.items() if k not in masked}
    return ops


def energy_fn(key: str):
    s = resolve(key)
    if s.energy is None:
        raise ValueError(f"thermo source {s.key!r} declares no energy loader")
    return s.energy()
