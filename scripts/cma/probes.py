"""Probe registry. See :class:`cma.kinds.Probe`."""

from __future__ import annotations

from .kinds import Probe
from .registry import Registry, ensure_loaded

_reg: Registry[Probe] = Registry("probe", key_attr="name")
register = _reg.register


def get(name: str) -> Probe:
    ensure_loaded()
    return _reg.get(name)


def all() -> list[Probe]:
    ensure_loaded()
    return _reg.all()


def keys() -> list[str]:
    ensure_loaded()
    return _reg.keys()
