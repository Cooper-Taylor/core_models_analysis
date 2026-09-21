"""Variant registry: cascade configs, direction overlays, and composites.

The pre-registry catalog typed ``cfg`` as a *mandatory* zero-arg
``ReversibilityConfig`` factory, which meant a variant whose directions came
from anywhere else could not be registered. Six such variants existed and lived
only in ``thermo_variants/manifest.json`` -- which ``export_thermo_variants.py``
rebuilt from the catalog on every run, erasing them.

``kind`` is the discriminator that fixes it. ``legacy()`` projects only the
cascade variants, in their original order, as the original seven-key dicts, so
the six call sites that do ``variant['cfg']()`` unconditionally keep working and
cannot see an entry they would choke on.
"""

from __future__ import annotations

from .kinds import Variant
from .registry import Registry, RegistryError, ensure_loaded

RESERVED_TAGS = frozenset({"baseline"})
KINDS = ("cascade", "overlay", "derived")

_reg: Registry[Variant] = Registry("variant", key_attr="tag")


def register(v: Variant) -> Variant:
    if v.kind not in KINDS:
        raise RegistryError(f"variant {v.tag!r}: kind must be one of {KINDS}, got {v.kind!r}")
    if v.kind == "cascade" and v.cfg is None:
        raise RegistryError(f"variant {v.tag!r}: kind='cascade' requires a cfg callable")
    if v.kind == "overlay" and not v.direction_map:
        raise RegistryError(f"variant {v.tag!r}: kind='overlay' requires direction_map")
    if v.kind == "derived" and v.combine is None:
        raise RegistryError(f"variant {v.tag!r}: kind='derived' requires a combine callable")
    if not tag_safe(v.tag):
        raise RegistryError(
            f"variant tag {v.tag!r} is not filesystem/URL safe -- it becomes a "
            "directory under thermo_variants/ and a key in site JSON"
        )
    return _reg.register(v)


_TAG_OK = frozenset("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-")


def tag_safe(tag: str) -> bool:
    # NB: this module defines a function named `all`, which shadows the builtin
    # at module scope -- so use a set test rather than `all(...)` here.
    return bool(tag) and not (set(tag) - _TAG_OK)


def get(tag: str) -> Variant:
    ensure_loaded()
    return _reg.get(tag)


def all(kind: str | None = None) -> list[Variant]:
    ensure_loaded()
    items = _reg.all() if kind is None else _reg.where(lambda v: v.kind == kind)
    return sorted(items, key=lambda v: (v.order, v.tag))


def tags(kind: str | None = None) -> list[str]:
    return [v.tag for v in all(kind)]


def legacy(kind: str = "cascade") -> list[dict]:
    """The pre-registry ``VARIANTS`` list: same entries, order and seven keys."""
    return [v.legacy_dict() for v in all(kind)]


def by_tag(tag: str) -> dict:
    """Drop-in for the old ``variant_catalog.variant_by_tag``."""
    return get(tag).legacy_dict()


def in_notebook() -> list[Variant]:
    return [v for v in all() if v.in_notebook]


def in_stats() -> list[Variant]:
    return [v for v in all() if v.in_stats]


def in_presentation() -> list[Variant]:
    return [v for v in all() if v.in_presentation]
