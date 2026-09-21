"""Turn a command-line string into a live object.

Every ``beginPipeline`` argument accepts either a **registry key** or a
**filesystem path**, and this module is the only place that decides which. That
is what lets you point the pipeline at a brand-new directory of models or a
direction table you produced five minutes ago, without editing any Python.

    resolve_model_set("ms2_gsm")                 # a registered set
    resolve_model_set("/data/my_models")         # a directory, inspected
    resolve_model_set("/data/models/*.json.gz")  # a glob
    resolve_model_set("ids.txt")                 # a file listing ids or paths

    resolve_directions("thermo_gc")              # registered
    resolve_directions("out/my_dirs.json")       # {rxn: op} JSON
    resolve_directions("out/my_dirs.csv")        # rxn_id,reversibility
    resolve_directions("thermo_variants/H4/Estimated_..._EQ.txt")   # MSDB report

    resolve_heuristic("H4")                      # registered variant
    resolve_heuristic("my_knobs.yaml")           # ReversibilityConfig knobs
    resolve_heuristic("sigma_band_k=1.96,mm_band=3")   # inline

    resolve_thermo("gc")                         # registered source
    resolve_thermo("my_energies.tsv")            # rxn_id, dg, dge

Ambiguity rule: a string that names a registry entry wins, unless it also
exists on disk, in which case the path wins and a note is emitted. Names and
paths rarely collide, and silently preferring the wrong one is the failure mode
worth making loud.
"""

from __future__ import annotations

import csv
import glob as _glob
import json
import os
from dataclasses import replace
from pathlib import Path

from . import directions as directions_reg
from . import models as models_reg
from . import sources as sources_reg
from . import variants as variants_reg
from .kinds import BiomassPolicy, Media, ModelSet, SeedKeyPolicy


class ResolveError(ValueError):
    """A command-line argument could not be turned into an object."""


def _is_pathish(s: str) -> bool:
    return ("/" in s or "\\" in s or s.startswith(".")
            or Path(s).exists() or bool(_glob.glob(s)))


#: Directory names that describe a *format*, not a dataset, so they make a
#: useless run label. ``.../modelseed2_gs_models/gmm/cobra`` should be called
#: ``gmm``, not ``cobra``.
_GENERIC_DIRNAMES = {"cobra", "models", "model", "json", "data", "out", "output", "sbml"}


def _label_for(d: Path) -> str:
    parts = [p for p in d.resolve().parts if p not in ("/", "")]
    for name in reversed(parts[-3:]):
        if name.lower() not in _GENERIC_DIRNAMES:
            return name
    return d.name


def _load_structured(path: Path) -> dict:
    """Read a .json / .yaml / .yml / .toml mapping."""
    suffix = path.suffix.lower()
    text = path.read_text()
    if suffix == ".json":
        return json.loads(text)
    if suffix in (".yaml", ".yml"):
        import yaml

        return yaml.safe_load(text) or {}
    if suffix == ".toml":
        import tomllib

        return tomllib.loads(text)
    raise ResolveError(
        f"{path}: unsupported config format {suffix!r}; use .json, .yaml or .toml"
    )


#: Column names the pipeline reserves for itself. A user file whose stem
#: collides with one would silently replace that column -- ``on_disk`` is the
#: unmodified-bounds baseline every diff is measured against.
RESERVED_LABELS = frozenset({"on_disk"})


def disambiguate(labels: list, sources: list) -> list:
    """Make labels unique, keeping them readable.

    Two ``--directions`` files called ``d.json`` in different directories would
    otherwise collapse into one column, silently discarding a whole input.
    Collisions gain their parent directory; if that is still not enough they
    gain an index.
    """
    out, seen = [], {}
    for lab, src in zip(labels, sources, strict=True):
        if lab in RESERVED_LABELS:
            lab = f"{lab}_input"
        base = lab
        if labels.count(base) > 1:
            parent = Path(src).parent.name
            lab = f"{parent}/{base}" if parent else base
        n = seen.get(lab, 0)
        seen[lab] = n + 1
        out.append(lab if n == 0 else f"{lab}#{n + 1}")
    return out


# ---------------------------------------------------------------------------
# Model sets
# ---------------------------------------------------------------------------
_LOADER_FOR_SUFFIX = {".json": "cobra_json", ".json.gz": "cobra_json_gz",
                      ".xml": "sbml", ".sbml": "sbml"}


def _detect_models(d: Path) -> tuple[str, str]:
    """(id_glob, loader) for a directory, by majority file type."""
    counts: dict[str, int] = {}
    for p in list(d.iterdir())[:400]:
        if not p.is_file():
            continue
        name = p.name
        for sfx in (".json.gz", ".json", ".xml", ".sbml", ".xml.gz"):
            if name.endswith(sfx):
                counts[sfx] = counts.get(sfx, 0) + 1
                break
    if not counts:
        raise ResolveError(
            f"{d}: no model files found (looked for *.json, *.json.gz, *.xml, *.sbml)"
        )
    if len(counts) > 1:
        shown = ", ".join(f"{n} x {s}" for s, n in sorted(counts.items(), key=lambda kv: -kv[1]))
        raise ResolveError(
            f"{d}: mixed model formats ({shown}). Picking the majority would "
            f"silently ignore the rest. Pass a glob that selects one format, "
            f"e.g. --models '{d}/*{max(counts, key=counts.get)}'."
        )
    sfx = next(iter(counts))
    loader = _LOADER_FOR_SUFFIX.get(sfx)
    if loader is None or loader == "sbml":
        raise ResolveError(
            f"{d}: {sfx} models are recognised but there is no loader for them yet. "
            "Convert to cobra JSON, or register a model set with "
            "loader='module:function'."
        )
    return f"*{sfx}", loader


def resolve_model_set(spec: str, media: str | None = None,
                      name: str | None = None) -> ModelSet:
    """A registry key, a directory, a glob, or a file listing ids/paths."""
    if not _is_pathish(spec) and spec in models_reg.keys():
        ms = models_reg.get(spec)
        # A registry entry only says where the models *should* be. Verify they
        # are actually there: otherwise a run resolves cleanly, reports "0 of 0
        # models", and fails later somewhere that does not mention the data.
        d = ms.models_dir()
        if not d.exists():
            raise ResolveError(
                f"model set {spec!r} is registered but its directory does not exist:\n"
                f"  {d}\n"
                "  Fetch the data with `python3 scripts/fetch_data.py --all`, or point at\n"
                "  an existing copy (see `beginPipeline --doctor`)."
            )
        if not any(d.glob(ms.id_glob)):
            raise ResolveError(
                f"model set {spec!r} points at {d}, which contains no files matching "
                f"{ms.id_glob!r}."
            )
        return replace(ms, media=media) if media else ms

    p = Path(spec).expanduser()

    # a glob
    if any(ch in spec for ch in "*?[") and not p.exists():
        matches = sorted(_glob.glob(spec))
        if not matches:
            raise ResolveError(f"{spec!r} matched no files")
        root = Path(os.path.commonpath([str(Path(m).parent) for m in matches]))
        pattern = "*" + "".join(Path(matches[0]).suffixes)
        _, loader = _detect_models(root)
        return _adhoc_set(name or _label_for(root), root, pattern, loader, media)

    if not p.exists():
        known = models_reg.keys()
        raise ResolveError(
            f"{spec!r} is neither a registered model set {known} nor an existing path"
        )

    if p.is_dir():
        pattern, loader = _detect_models(p)
        return _adhoc_set(name or _label_for(p), p, pattern, loader, media)

    # a single model file
    if "".join(p.suffixes) in _LOADER_FOR_SUFFIX or p.suffix in _LOADER_FOR_SUFFIX:
        sfx = ".json.gz" if p.name.endswith(".json.gz") else p.suffix
        return _adhoc_set(name or p.stem, p.parent, p.name, _LOADER_FOR_SUFFIX[sfx], media)

    # a text file listing ids or paths
    entries = [ln.strip() for ln in p.read_text().splitlines()
               if ln.strip() and not ln.startswith("#")]
    if not entries:
        raise ResolveError(f"{p}: empty id list")

    if "/" in entries[0] or Path(entries[0]).is_absolute():
        root = Path(entries[0]).parent
        pattern, loader = _detect_models(root)
    else:
        # Bare ids. They carry no directory, so find the registered model set
        # that actually contains them rather than guessing the list file's own
        # directory -- which is almost never where the models live.
        root = pattern = loader = None
        for cand in models_reg.all():
            d = cand.models_dir()
            if not d.exists():
                continue
            sfx = "".join(Path(cand.id_glob).suffixes)
            if (d / f"{entries[0]}{sfx}").exists():
                root, pattern, loader = d, cand.id_glob, cand.loader
                name = name or f"{cand.name}:{p.stem}"
                break
        if root is None:
            searched = [str(c.models_dir()) for c in models_reg.all()]
            raise ResolveError(
                f"{p}: lists bare model ids (first: {entries[0]!r}) but no registered "
                f"model set contains them.\nSearched: {', '.join(searched)}\n"
                "Either write full paths in the list, or pass the directory with "
                "--models <dir> and the id list with --panel <path>."
            )
    ms = _adhoc_set(name or p.stem, root, pattern, loader, media)
    # Strip ONLY a real model-file extension. Path.suffixes treats an assembly
    # version as one, so "GCF_000005825.2" used to become "GCF_000005825" and
    # every resolved path pointed at a file that does not exist.
    model_sfx = "".join(Path(pattern).suffixes)
    ids = []
    for e in entries:
        nm = Path(e).name
        for sfx in (model_sfx, ".json.gz", ".json", ".xml", ".sbml"):
            if sfx and nm.endswith(sfx):
                nm = nm[: -len(sfx)]
                break
        ids.append(nm)
    object.__setattr__(ms, "_explicit_ids", ids)
    return ms


def _adhoc_set(name: str, root: Path, pattern: str, loader: str,
               media: str | None) -> ModelSet:
    """A ModelSet for a path the registry has never heard of.

    Conservative defaults: the widest biomass precedence, and a SEED-key policy
    that falls back to the reaction-id prefix so a model whose reactions carry
    no ``seed.reaction`` annotation still matches a direction map.
    """
    if not media:
        print(f"note: model set {name!r} has no declared medium; using "
              "'kbase_complete' (347 compounds, uptake -1000). If these models "
              "were gap-filled on a defined medium, pass --media <file> -- the "
              "difference is roughly 200x in growth flux.")
    return ModelSet(
        name=name,
        models_dir=lambda root=root: root,
        id_glob=pattern,
        loader=loader,
        media=media or "kbase_complete",
        biomass=BiomassPolicy(prefer=("bio1",), then=("bio2", "biomass", "Biomass")),
        seed_key=SeedKeyPolicy(fallback="id_prefix", require_nonzero_overrides=False),
        results_namespace=name,
        description=f"ad-hoc model set from {root}",
    )


def model_ids(ms: ModelSet) -> list[str]:
    explicit = getattr(ms, "_explicit_ids", None)
    return list(explicit) if explicit else ms.model_ids()


# ---------------------------------------------------------------------------
# Media
# ---------------------------------------------------------------------------
def resolve_media(spec: str) -> Media:
    """A registered media name, or a ``.cpd`` / ``.json`` / ``.tsv`` file."""
    if not _is_pathish(spec):
        return models_reg.media(spec)
    p = Path(spec).expanduser()
    if not p.exists():
        raise ResolveError(f"media {spec!r} is neither registered nor a path")

    def loader(p=p):
        if p.suffix == ".json":
            d = json.loads(p.read_text())
            if isinstance(d, dict) and "mediacompounds" in d:
                out = {}
                for mc in d["mediacompounds"]:
                    cid = mc.get("id") or mc.get("compound_ref", "").rsplit("/", 1)[-1]
                    if cid:
                        out[cid] = -float(abs(mc.get("maxFlux", 100)))
                return out
            return {k: float(v) for k, v in d.items()}
        out = {}
        for line in p.read_text().splitlines():
            parts = line.replace(",", "\t").split("\t")
            cid = parts[0].strip()
            if cid and not cid.startswith("#"):
                out[cid] = float(parts[1]) if len(parts) > 1 else default_uptake
        return out

    # A bare .cpd list carries no flux column. Defaulting it to -1000 makes a
    # defined minimal medium behave like a rich one: the same 20-compound file
    # gave ~200x the growth of the registered version, which sets -10.
    # A small list is a defined medium; a large one is a complete medium.
    n_lines = sum(1 for ln in p.read_text().splitlines()
                  if ln.strip() and not ln.startswith("#"))
    default_uptake = -1000.0 if n_lines > 100 else -10.0
    return Media(name=p.stem, loader=loader, default_uptake=default_uptake,
                 description=(f"media from {p}; {n_lines} compounds, uptake "
                              f"{default_uptake} where the file gives no flux column"))


# ---------------------------------------------------------------------------
# Direction maps
# ---------------------------------------------------------------------------
def resolve_directions(spec: str) -> tuple[str, dict]:
    """``(label, {rxn_id: operator})`` from a registry key or any table on disk."""
    from .directions import DirectionMap, normalize_operator

    if not _is_pathish(spec) and spec in directions_reg.keys():
        return spec, DirectionMap.load(spec).ops

    p = Path(spec).expanduser()
    if not p.exists():
        raise ResolveError(
            f"{spec!r} is neither a registered direction map {directions_reg.keys()} "
            "nor an existing path"
        )

    raw: dict = {}
    if p.suffix == ".json":
        d = json.loads(p.read_text())
        raw = {k: v for k, v in d.items()} if isinstance(d, dict) else {}
    elif p.suffix in (".csv", ".tsv"):
        delim = "," if p.suffix == ".csv" else "\t"
        with p.open() as fh:
            for row in csv.DictReader(fh, delimiter=delim):
                rid = (row.get("rxn_id") or row.get("id") or row.get("reaction_id")
                       or row.get("reaction"))
                op = (row.get("reversibility") or row.get("operator")
                      or row.get("direction") or row.get("new_rev"))
                if rid and op:
                    raw[rid] = op
    else:
        # MSDB-format report: rxn \t status \t [old_rev \t] new_rev
        for line in p.read_text().splitlines():
            if not line.strip():
                continue
            c = line.split("\t")
            if len(c) >= 3:
                raw[c[0]] = c[-1]
    if not raw:
        raise ResolveError(f"{p}: parsed zero direction entries")
    ops = {k: normalize_operator(v) for k, v in raw.items()}
    label = p.stem
    if label in RESERVED_LABELS:
        label = f"{label}_input"
    DirectionMap(key=label, ops=ops).validate()
    return label, ops


# ---------------------------------------------------------------------------
# Heuristics (ReversibilityConfig knob-sets)
# ---------------------------------------------------------------------------
def heuristic_knobs() -> dict:
    """``{knob: type}`` for every field of ReversibilityConfig -- the documentation."""
    import dataclasses

    import reversibility_lib as lib

    return {f.name: f.type for f in dataclasses.fields(lib.ReversibilityConfig)}


def resolve_heuristic(spec: str):
    """``(label, ReversibilityConfig)`` from a variant tag, a config file, or inline knobs."""
    import reversibility_lib as lib

    if not _is_pathish(spec) and "=" not in spec:
        if spec in variants_reg.tags("cascade"):
            return spec, variants_reg.get(spec).cfg()
        raise ResolveError(
            f"{spec!r} is not a registered cascade variant "
            f"{variants_reg.tags('cascade')}, a config file, or inline knobs"
        )

    if "=" in spec and not Path(spec).exists():
        knobs = {}
        for part in spec.split(","):
            if not part.strip():
                continue
            k, _, v = part.partition("=")
            knobs[k.strip()] = _coerce(v.strip())
        return f"inline({spec})", _build_cfg(lib, knobs, spec)

    p = Path(spec).expanduser()
    if not p.exists():
        raise ResolveError(f"heuristic {spec!r} is neither a variant tag nor a path")
    data = _load_structured(p)
    knobs = data.get("knobs", data)
    label = data.get("tag") or data.get("name") or p.stem
    return label, _build_cfg(lib, knobs, str(p))


def _coerce(v: str):
    low = v.lower()
    if low in ("none", "null", ""):
        return None
    if low in ("true", "false"):
        return low == "true"
    try:
        return int(v)
    except ValueError:
        pass
    try:
        return float(v)
    except ValueError:
        return v


# Knobs whose value is a data table rather than a scalar; a config file names a
# loader or a path and this expands it.
_TABLE_KNOBS = {
    "ln_ri_by_rxn": "load_ln_reversibility_index",
    "energy_override_by_rxn": "load_dgpredictor_energies",
    "per_met_conc_range": None,
    "per_met_conc": None,
}


def _build_cfg(lib, knobs: dict, origin: str):
    import dataclasses

    valid = {f.name for f in dataclasses.fields(lib.ReversibilityConfig)}
    clean = {}
    for k, v in knobs.items():
        if k not in valid:
            raise ResolveError(
                f"{origin}: unknown heuristic knob {k!r}.\n"
                f"Valid knobs: {', '.join(sorted(valid))}"
            )
        if k in _TABLE_KNOBS and isinstance(v, str):
            # Three ways to name a data table, tried in order:
            #   'auto'                      -> the loader reversibility_lib provides
            #   a name in reversibility_lib -> that module-level table
            #   anything else               -> a path on disk
            if v == "auto":
                loader = _TABLE_KNOBS[k]
                if not loader:
                    raise ResolveError(
                        f"{origin}: {k}='auto' has no default loader; "
                        "name a table in reversibility_lib or give a path"
                    )
                v = getattr(lib, loader)()
            elif getattr(lib, v, None) is not None and not Path(v).expanduser().exists():
                v = getattr(lib, v)
            else:
                pth = Path(v).expanduser()
                if not pth.exists():
                    named = [n for n in dir(lib) if n.isupper() and isinstance(getattr(lib, n), dict)]
                    raise ResolveError(
                        f"{origin}: {k}={v!r} is neither 'auto', a table in "
                        f"reversibility_lib ({', '.join(named) or 'none'}), nor a readable path"
                    )
                if k == "ln_ri_by_rxn":
                    v = lib.load_ln_reversibility_index(str(pth))
                elif pth.suffix in (".json", ".yaml", ".yml", ".toml"):
                    v = _load_structured(pth)
                else:
                    v = _two_col_table(pth)
        clean[k] = v
    return lib.ReversibilityConfig(**clean)


def _two_col_table(p: Path) -> dict:
    """``{id: (value, uncertainty)}`` from a 3-column CSV/TSV, sentinels dropped."""
    out = {}
    for line in p.read_text().splitlines():
        c = line.replace(",", "\t").split("\t")
        if len(c) >= 3 and not c[0].startswith("#"):
            try:
                a, b = float(c[1]), float(c[2])
            except ValueError:
                continue
            if abs(a) >= THERMO_SENTINEL or abs(b) >= THERMO_SENTINEL:
                continue
            out[c[0].strip()] = (a, b)
    return out


# ---------------------------------------------------------------------------
# Thermo sources
# ---------------------------------------------------------------------------
def resolve_thermo(spec: str) -> tuple[str, dict]:
    """``(label, {rxn_id: (dg, dge)})`` from a registered source or a table."""
    if not _is_pathish(spec):
        src = sources_reg.resolve(spec)
        return src.key, _energies_from_msdb(src)
    p = Path(spec).expanduser()
    if not p.exists():
        raise ResolveError(f"thermo source {spec!r} is neither registered nor a path")
    if p.suffix in (".json", ".yaml", ".yml", ".toml"):
        d = _load_structured(p)
        table = {k: (float(v[0]), float(v[1])) for k, v in d.items()}
    else:
        table = _two_col_table(p)
    if not table:
        raise ResolveError(
            f"{p}: parsed zero usable rows. Expected three columns "
            "(reaction id, dG, uncertainty), comma- or tab-separated. "
            "Rows carrying ModelSEED's 10000000 sentinel are dropped as "
            "'no estimate'."
        )
    return p.stem, table


#: ModelSEED writes this in place of a ΔG′° it could not estimate. Treating it
#: as a number makes every downstream statistic meaningless -- the median of the
#: group-contribution table comes out as 10000000 rather than about -2 kcal/mol.
THERMO_SENTINEL = 10000000


def _energies_from_msdb(src) -> dict:
    """``{rxn: (dg, dge)}`` for one source, sentinel rows dropped.

    A reaction whose ΔG′° or uncertainty is the sentinel has no estimate from
    this source, so it is absent from the table rather than present with a
    placeholder value.
    """
    from . import paths

    out = {}
    if not src.label:
        return out
    for shard in sorted((paths.msdb("live") / "Biochemistry").glob("reaction_[0-9][0-9].json")):
        for r in json.loads(shard.read_text()):
            t = r.get("thermodynamics")
            if isinstance(t, dict) and src.label in t:
                pair = t[src.label]
                try:
                    dg, dge = float(pair[0]), float(pair[1])
                except (TypeError, ValueError, IndexError):
                    continue
                if abs(dg) >= THERMO_SENTINEL or abs(dge) >= THERMO_SENTINEL:
                    continue
                out[r["id"]] = (dg, dge)
    return out
