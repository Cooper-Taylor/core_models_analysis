"""Figure sets.

Mirrors ``scripts/figures.tsv``, which stays the hand-edited source of truth --
the ``make-plot`` and ``change-plot`` skills depend on that, and
``regen_figures.py`` keeps working unchanged. Registering the rows here adds
what the TSV cannot express: the ``needs`` edges into the stage DAG and the
palette family a figure draws from, so a source added to the registry appears
in the right plots without a literal edit.
"""

from __future__ import annotations

import csv

from .. import paths
from ..figures import register
from ..kinds import FigureSet


def _rows():
    tsv = paths.scripts("figures.tsv")
    if not tsv.exists():
        return []
    lines = [ln for ln in tsv.read_text().splitlines() if not ln.startswith(">")]
    return [r for r in csv.DictReader(lines, delimiter="\t") if r.get("name")]


for _r in _rows():
    register(FigureSet(
        name=_r["name"],
        script=_r.get("script", ""),
        args=_r.get("args") or "",
        out_dir=_r.get("out_dir") or "",
        tags=tuple(t.strip() for t in (_r.get("tags") or "").split(",") if t.strip()),
        values_sidecar=_r.get("values") or None,
    ))
