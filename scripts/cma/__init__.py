"""``cma`` -- the registry core for core_models_analysis.

Adding a thermo source, a heuristic variant, a direction map, a model set, a
figure or a probe is one entry in ``cma/entries/``. Everything downstream reads
the registry instead of its own copy of the literal.

    from cma import paths, sources, variants, models, directions

Import is deliberately cheap: no ``os.environ`` is read at import time and
neither cobra, pandas, numpy nor matplotlib is imported at module scope
(``site/serve.py`` is stdlib-only by design). Submodules are resolved lazily.

Note: this package shadows the PyPI distribution ``cma`` (CMA-ES), which this
project does not use. Do not add it to ``requirements.txt``.
"""

from __future__ import annotations

import importlib
from typing import TYPE_CHECKING

__all__ = [
    "paths", "kinds", "registry", "sources", "variants", "directions",
    "models", "probes", "figures", "manifest",
]

if TYPE_CHECKING:  # pragma: no cover
    from . import (  # noqa: F401
        directions,
        figures,
        kinds,
        manifest,
        models,
        paths,
        probes,
        registry,
        sources,
        variants,
    )


def __getattr__(name: str):
    if name in __all__:
        mod = importlib.import_module(f".{name}", __name__)
        globals()[name] = mod
        return mod
    raise AttributeError(f"module 'cma' has no attribute {name!r}")


def __dir__():
    return sorted(__all__)
