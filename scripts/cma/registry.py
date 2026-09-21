"""A tiny ordered, duplicate-rejecting registry.

Every plugin kind in :mod:`cma` is an instance of :class:`Registry`. The
contract is deliberately small:

* insertion order is preserved (it is what drives figure axes and table rows);
* registering the same key twice raises, rather than silently winning;
* ``load()`` imports a fixed, named list of entry modules -- never a directory
  scan. Ten modules under ``scripts/`` call ``sys.exit`` at import time, so a
  glob-and-import would take the process down.
"""

from __future__ import annotations

import threading
from collections.abc import Callable, Iterator
from typing import Generic, TypeVar

T = TypeVar("T")


class RegistryError(RuntimeError):
    """Raised for duplicate keys, unknown keys, or post-freeze mutation."""


class Registry(Generic[T]):
    def __init__(self, kind: str, key_attr: str = "key") -> None:
        self.kind = kind
        self._key_attr = key_attr
        self._items: dict[str, T] = {}
        self._frozen = False

    # -- registration -------------------------------------------------------
    def register(self, item: T) -> T:
        key = getattr(item, self._key_attr)
        if self._frozen:
            raise RegistryError(
                f"{self.kind} registry is frozen; cannot register {key!r}. "
                "Register in scripts/cma/entries/, which loads before freeze."
            )
        if key in self._items:
            raise RegistryError(
                f"duplicate {self.kind} key {key!r} -- already registered by "
                f"{type(self._items[key]).__name__}"
            )
        self._items[key] = item
        return item

    def freeze(self) -> None:
        self._frozen = True

    # -- lookup -------------------------------------------------------------
    def get(self, key: str) -> T:
        try:
            return self._items[key]
        except KeyError:
            raise RegistryError(
                f"unknown {self.kind} {key!r}. Registered: {', '.join(self._items) or '(none)'}"
            ) from None

    def maybe(self, key: str) -> T | None:
        return self._items.get(key)

    def all(self) -> list[T]:
        return list(self._items.values())

    def keys(self) -> list[str]:
        return list(self._items)

    def where(self, pred: Callable[[T], bool]) -> list[T]:
        return [v for v in self._items.values() if pred(v)]

    def __contains__(self, key: object) -> bool:
        return key in self._items

    def __iter__(self) -> Iterator[T]:
        return iter(self._items.values())

    def __len__(self) -> int:
        return len(self._items)

    def __repr__(self) -> str:
        return f"<Registry {self.kind}: {len(self._items)} entries>"


# ---------------------------------------------------------------------------
# One-shot loading of the entry modules.
# ---------------------------------------------------------------------------
# Fixed order. thermo_sources before variants because an overlay variant may
# name a source; model_sets before probes because a probe may name a media.
ENTRY_MODULES = (
    "cma.entries.model_sets",
    "cma.entries.thermo_sources",
    "cma.entries.direction_maps",
    "cma.entries.variants",
    "cma.entries.probes",
    "cma.entries.figure_sets",
)

_loaded = False
_lock = threading.Lock()


def load(force: bool = False) -> None:
    """Import every entry module exactly once, then freeze the registries."""
    global _loaded
    with _lock:
        if _loaded and not force:
            return
        import importlib

        for name in ENTRY_MODULES:
            importlib.import_module(name)
        _loaded = True


def ensure_loaded() -> None:
    if not _loaded:
        load()
