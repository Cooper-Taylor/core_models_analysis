"""``thermo_variants/manifest.json``: merge, never rewrite.

The bug this fixes
------------------
``export_thermo_variants.py`` opened ``manifest = {"variants": []}``, filled it
from ``variant_catalog.VARIANTS`` (or from ``--only``), and dumped it wholesale.
The manifest holds **20** variants; the catalog holds **14**. Any full export
therefore erased the six variants produced by the standalone overlay builders
(``ai_opus48``, ``consensus_thermo``, ``eq3_beber2022``, ``eq3_gamma1``,
``group_contribution``, ``kegg_implicit``), and an ``--only`` run reduced the
manifest to just the tags named on the command line. The directories survived
on disk; the index the site reads did not.

Five of the overlay builders then each carried their own hand-copied
``append_manifest`` body. This module is the single implementation.
"""

from __future__ import annotations

import json
import os
import tempfile
from collections.abc import Iterable, Mapping
from pathlib import Path

from . import paths


def manifest_path(root: Path | None = None) -> Path:
    return (Path(root) if root else paths.thermo_variants()) / "manifest.json"


def read(root: Path | None = None) -> dict:
    p = manifest_path(root)
    if not p.exists():
        return {"variants": []}
    try:
        data = json.loads(p.read_text())
    except json.JSONDecodeError as exc:
        raise ValueError(f"{p} is not valid JSON: {exc}") from None
    data.setdefault("variants", [])
    return data


def entry_for(tag: str, root: Path | None = None) -> dict | None:
    for v in read(root).get("variants", []):
        if v.get("tag") == tag:
            return v
    return None


def _atomic_write(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=str(path.parent), suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as fh:
            json.dump(payload, fh, indent=2, default=str)
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def merge(summaries: Iterable[Mapping], root: Path | None = None) -> dict:
    """Upsert ``summaries`` into the manifest, preserving every other entry.

    Entries are matched on ``tag``. An existing tag is replaced in place, so
    ordering is stable across runs; a new tag is appended. Nothing is ever
    dropped -- removing an entry is :func:`prune`, which requires confirmation.
    """
    data = read(root)
    variants = list(data.get("variants", []))
    index = {v.get("tag"): i for i, v in enumerate(variants) if v.get("tag")}
    for s in summaries:
        tag = s.get("tag")
        if not tag:
            raise ValueError(f"manifest entry has no 'tag': {sorted(s)}")
        if tag in index:
            variants[index[tag]] = dict(s)
        else:
            index[tag] = len(variants)
            variants.append(dict(s))
    data["variants"] = variants
    _atomic_write(manifest_path(root), data)
    return data


def prune(tags: Iterable[str], root: Path | None = None, confirm: bool = False) -> dict:
    """Remove entries. Refuses unless ``confirm=True`` -- removal is the rare case."""
    tags = set(tags)
    if not confirm:
        raise ValueError(
            f"prune() would remove {sorted(tags)} from the manifest; "
            "pass confirm=True if that is intended"
        )
    data = read(root)
    data["variants"] = [v for v in data.get("variants", []) if v.get("tag") not in tags]
    _atomic_write(manifest_path(root), data)
    return data


def validate(root: Path | None = None) -> list[str]:
    """Report inconsistencies between the manifest, the registry and the disk."""
    from . import variants as variant_reg

    problems: list[str] = []
    data = read(root)
    listed = [v.get("tag") for v in data.get("variants", [])]
    if len(listed) != len(set(listed)):
        dupes = {t for t in listed if listed.count(t) > 1}
        problems.append(f"duplicate tags in manifest: {sorted(dupes)}")

    tv_root = Path(root) if root else paths.thermo_variants()
    on_disk = {p.name for p in tv_root.iterdir() if p.is_dir()} if tv_root.exists() else set()
    for tag in on_disk - set(listed):
        problems.append(f"{tag}: directory on disk but absent from manifest")
    for tag in set(listed) - on_disk:
        problems.append(f"{tag}: in manifest but no directory on disk")

    try:
        registered = set(variant_reg.tags())
    except Exception:  # registry not loadable in this context
        registered = set()
    if registered:
        for tag in set(listed) - registered:
            problems.append(f"{tag}: in manifest but not registered in cma/entries/variants.py")
    return problems
