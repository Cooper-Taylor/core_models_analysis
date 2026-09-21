"""Figure-set registry. The TSV at ``scripts/figures.tsv`` remains the source."""

from __future__ import annotations

from .kinds import FigureSet
from .registry import Registry, ensure_loaded

_reg: Registry[FigureSet] = Registry("figure set", key_attr="name")
register = _reg.register


def get(name: str) -> FigureSet:
    ensure_loaded()
    return _reg.get(name)


def all(tag: str | None = None) -> list[FigureSet]:
    ensure_loaded()
    return _reg.all() if tag is None else _reg.where(lambda f: tag in f.tags)


def keys() -> list[str]:
    ensure_loaded()
    return _reg.keys()
