"""Model sets, media and panels.

``growth_heuristics.MODELS_DIR`` was a module constant, and 19 other modules
declared their own copy (5 of them as absolute paths with no env-var escape).
The literal ``core_models_kegg2`` appears in 43 lines. A second collection of
models therefore could not coexist with the first -- you overwrote
``results/`` or you forked the script.

A model set declares what used to be assumed: where the files are, how to load
them, which medium, how to find biomass, and how a model reaction maps to an
MSDB reaction id. ``results_namespace=None`` means "write exactly where the
pipeline already writes", which is what makes adoption a no-op for the core set.
"""

from __future__ import annotations

import os

from .kinds import Media, ModelSet, Panel
from .registry import Registry, ensure_loaded

_sets: Registry[ModelSet] = Registry("model set", key_attr="name")
_media: Registry[Media] = Registry("media", key_attr="name")
_panels: Registry[Panel] = Registry("panel", key_attr="name")

register = _sets.register
register_media = _media.register
register_panel = _panels.register

DEFAULT_SET = "core_kegg2"


# --- lookup ----------------------------------------------------------------
def get(name: str | None = None) -> ModelSet:
    ensure_loaded()
    return _sets.get(name or active_name())


def all() -> list[ModelSet]:
    ensure_loaded()
    return _sets.all()


def keys() -> list[str]:
    ensure_loaded()
    return _sets.keys()


def media(name: str) -> Media:
    ensure_loaded()
    return _media.get(name)


def panel(name: str) -> Panel:
    ensure_loaded()
    return _panels.get(name)


def panels(model_set: str | None = None) -> list[Panel]:
    ensure_loaded()
    if model_set is None:
        return _panels.all()
    return _panels.where(lambda p: p.model_set == model_set)


# --- the active set --------------------------------------------------------
def active_name() -> str:
    """Which model set is in force.

    Read from the environment on *every* call rather than cached at import, so
    a spawn child re-reads it and a parent that switches sets after import is
    not left on a stale path. A mismatch there surfaces downstream as
    ``n_overrides=0``, which reads like a successful run.
    """
    return os.environ.get("CMA_MODEL_SET", DEFAULT_SET)


def activate(name: str) -> ModelSet:
    """Make ``name`` the active set for this process and its children."""
    ms = get(name)
    os.environ["CMA_MODEL_SET"] = ms.name
    return ms


def active() -> ModelSet:
    return get(active_name())


# --- helpers used by the runner and by growth_heuristics -------------------
def load_model(model_set: ModelSet, model_id: str):
    """Load one model, honouring the set's declared loader."""
    path = model_set.path_for(model_id)
    kind = model_set.loader
    if kind == "cobra_json":
        from cobra.io import load_json_model

        return load_json_model(str(path))
    if kind == "cobra_json_gz":
        import gzip
        import json

        from cobra.io import model_from_dict

        with gzip.open(path, "rt", encoding="utf-8") as fh:
            return model_from_dict(json.load(fh))
    if ":" in kind:
        import importlib

        mod_name, fn_name = kind.split(":", 1)
        return getattr(importlib.import_module(mod_name), fn_name)(path)
    raise ValueError(f"model set {model_set.name!r}: unknown loader {kind!r}")


def seed_id(reaction, policy) -> str | None:
    """Map a cobra reaction to its MSDB reaction id under a SeedKeyPolicy."""
    import re

    anno = getattr(reaction, "annotation", None)
    if anno is None and isinstance(reaction, dict):
        anno = reaction.get("annotation")
    raw = anno.get(policy.annotation) if anno else None
    if isinstance(raw, (list, tuple)):
        if not policy.allow_list_values:
            raise TypeError(
                f"{policy.annotation} is a list on {getattr(reaction, 'id', '?')}; "
                "set allow_list_values=True on this model set's SeedKeyPolicy"
            )
        raw = raw[0] if raw else None
    if not raw and policy.fallback == "id_prefix":
        raw = getattr(reaction, "id", None) or (reaction.get("id") if isinstance(reaction, dict) else None)
        if raw:
            raw = re.sub(r"_[a-z]\d*$", "", raw)
        if raw and not raw.startswith("rxn"):
            raw = None
    if not raw:
        return None
    if policy.strip_compartment_suffix:
        raw = re.sub(r"_[a-z]$", "", raw)
    return raw


def find_biomass(model, policy):
    """Biomass reaction under a BiomassPolicy (mirrors the legacy precedence)."""
    for rid in (*policy.prefer, *policy.then):
        if rid in model.reactions:
            return model.reactions.get_by_id(rid)
    for r in model.reactions:
        if r.id.lower().startswith(policy.fallback_prefix) and not r.id.startswith(policy.exclude_prefix):
            return r
    return None


def apply_media(model, media_spec: Media) -> int:
    """Restrict uptake to the medium. Returns the number of exchanges opened."""
    table = media_spec.loader()
    opened = 0
    for rxn in model.reactions:
        if not rxn.id.startswith("EX_"):
            continue
        mets = list(rxn.metabolites)
        if len(mets) != 1:
            continue
        cpd = mets[0].id.split("_")[0]
        if cpd in table:
            rxn.lower_bound = float(table[cpd])
            opened += 1
        else:
            rxn.lower_bound = 0.0
        if rxn.upper_bound < 1000.0:
            rxn.upper_bound = 1000.0
    return opened
