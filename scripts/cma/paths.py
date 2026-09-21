"""Every path the analysis reads or writes, resolved in one place.

Resolution order, for each root, first hit wins:

1. an environment variable (``MSDB_ROOT``, ``CORE_MODELS_ANALYSIS_DIR``, ...)
2. a key in a config file -- ``cma.toml`` beside the repo, or
   ``$XDG_CONFIG_HOME/cma/config.toml`` (default ``~/.config/cma/config.toml``)
3. a search of conventional locations relative to the repo
4. for optional roots, ``None``; for required ones, an error naming what to set

Nothing is hardcoded to one machine. ``analysis_dir()`` is derived from where
this file lives, so a clone works wherever it is put, and the external data
roots are searched for next to the repo before the user is asked.

Two roots that used to share a name are separate here. ``MSDB_ROOT`` is the
live ModelSEEDDatabase clone; ``MSDB_SNAPSHOT_ROOT`` is a pinned snapshot that
part of the grading pipeline reads. They held different biochemistry, so
exporting one variable used to silently repoint half the scripts at the other's
data.

    from cma import paths
    paths.analysis_dir()          # this checkout
    paths.msdb()                  # the live MSDB clone, or a clear error
    paths.msdb(required=False)    # None if it is not installed
    paths.doctor()                # what is present, what is missing, what to do
"""

from __future__ import annotations

import os
import shutil
import sys
from functools import lru_cache
from pathlib import Path

# scripts/cma/paths.py -> repo root is two levels up
_PKG = Path(__file__).resolve().parent
_REPO = _PKG.parent.parent


class MissingRoot(RuntimeError):
    """A required external data root could not be located."""


# ---------------------------------------------------------------------------
# config file
# ---------------------------------------------------------------------------
def config_paths() -> list[Path]:
    xdg = os.environ.get("XDG_CONFIG_HOME")
    home_cfg = Path(xdg) / "cma" / "config.toml" if xdg else Path.home() / ".config" / "cma" / "config.toml"
    return [_REPO / "cma.toml", _REPO.parent / "cma.toml", home_cfg]


@lru_cache(maxsize=1)
def _config() -> dict:
    for p in config_paths():
        if p.is_file():
            try:
                import tomllib

                data = tomllib.loads(p.read_text())
            except Exception:  # noqa: BLE001 -- a broken config must not be fatal
                continue
            out = dict(data.get("paths", {}))
            out["_config_file"] = str(p)
            return out
    return {}


def config_file() -> str | None:
    return _config().get("_config_file")


# ---------------------------------------------------------------------------
# root resolution
# ---------------------------------------------------------------------------
#: name -> (env var, config key,candidate relative locations, what it is)
_ROOTS = {
    "msdb": ("MSDB_ROOT", "msdb",
             ("../ModelSEEDDatabase", "ModelSEEDDatabase", "data/ModelSEEDDatabase",
              "../../ModelSEEDDatabase"),
             "the ModelSEEDDatabase clone (biochemistry, media, thermodynamics)"),
    # Deliberately does NOT fall back to the live clone. The two held different
    # biochemistry, and quietly substituting one for the other is the bug this
    # separation exists to prevent: callers pass required=False and handle None.
    "msdb_snapshot": ("MSDB_SNAPSHOT_ROOT", "msdb_snapshot",
                      ("../tmp/devsnap2", "../devsnap2", "data/devsnap2"),
                      "a pinned MSDB snapshot, read by the grading/scatter scripts",
                      "Biochemistry/reaction_00.json"),
    "core_models": ("CORE_MODELS_DIR", "core_models",
                    ("data/core_models_kegg2", "../core_models_kegg2"),
                    "the 5,683 core model JSONs", "*.json"),
    "gs_models": ("MS2_GS_MODELS_DIR", "gs_models",
                  ("data/modelseed2_gs_models", "../modelseed2_gs_models"),
                  "the ModelSEED v2 genome-scale models", "gmm"),
}


def _looks_right(p: Path, marker: str) -> bool:
    """Is this directory actually the thing we are looking for?"""
    if not p.is_dir():
        return False
    if "*" in marker:
        return any(p.glob(marker))
    return (p / marker).exists()


def _resolve(name: str, required: bool) -> Path | None:
    env, key, candidates, what, marker = _ROOTS[name]
    raw = os.environ.get(env)
    if raw:
        p = Path(raw).expanduser()
        if not p.exists():
            raise MissingRoot(f"{env}={raw!r} does not exist (expected {what})")
        if not _looks_right(p, marker):
            raise MissingRoot(
                f"{env}={raw!r} exists but does not look like {what}: "
                f"expected to find {marker!r} inside it"
            )
        return p

    cfg = _config().get(key)
    if cfg:
        p = Path(cfg).expanduser()
        if not p.exists():
            raise MissingRoot(f"{config_file()} sets paths.{key}={cfg!r}, which does not exist")
        if not _looks_right(p, marker):
            raise MissingRoot(
                f"{config_file()} sets paths.{key}={cfg!r}, which does not look like "
                f"{what}: expected {marker!r} inside it"
            )
        return p

    for rel in candidates:
        p = (_REPO / rel).resolve()
        if _looks_right(p, marker):
            return p

    if not required:
        return None
    tried = ", ".join(str((_REPO / c).resolve()) for c in candidates)
    raise MissingRoot(
        f"cannot find {what}.\n"
        f"  Set {env}=/path/to/it, or add it to cma.toml as paths.{key}.\n"
        f"  Looked in: {tried}\n"
        f"  Run `beginPipeline --doctor` for the full picture."
    )


def analysis_dir() -> Path:
    """This checkout. Derived from the package location, never hardcoded."""
    raw = os.environ.get("CORE_MODELS_ANALYSIS_DIR")
    return Path(raw).expanduser() if raw else _REPO


def msdb(which: str = "live", required: bool = True) -> Path | None:
    if which == "live":
        return _resolve("msdb", required)
    if which == "snapshot":
        return _resolve("msdb_snapshot", required)
    raise ValueError(f"msdb(which=) must be 'live' or 'snapshot', got {which!r}")


def msdb_snapshot(required: bool = True) -> Path | None:
    return _resolve("msdb_snapshot", required)


def core_models_dir(required: bool = True) -> Path | None:
    return _resolve("core_models", required)


def gs_models_dir(required: bool = True) -> Path | None:
    return _resolve("gs_models", required)


def msdb_code() -> Path:
    return msdb("live") / "Libs" / "Python"


def add_msdb_to_syspath() -> None:
    p = str(msdb_code())
    if p not in sys.path:
        sys.path.insert(0, p)


def python_exe() -> str:
    """Interpreter used when a driver forks a child job.

    Defaults to the interpreter running right now, which is correct in a venv,
    a conda env or a bare system Python alike. ``CMA_PYTHON`` overrides it.
    """
    return os.environ.get("CMA_PYTHON") or sys.executable


# --- output trees ----------------------------------------------------------
def results(*parts: str, model_set: str | None = None) -> Path:
    base = analysis_dir() / "results"
    if model_set:
        base = base / model_set
    return base.joinpath(*parts) if parts else base


def reports(*parts: str) -> Path:
    return (analysis_dir() / "reports").joinpath(*parts)


def figures(*parts: str) -> Path:
    return (analysis_dir() / "reports" / "figures").joinpath(*parts)


def thermo_variants(*parts: str) -> Path:
    return (analysis_dir() / "thermo_variants").joinpath(*parts)


def site_data(*parts: str) -> Path:
    return (analysis_dir() / "site" / "data").joinpath(*parts)


def scripts(*parts: str) -> Path:
    return (analysis_dir() / "scripts").joinpath(*parts)


def kbcache_root(stage: str = "notebooks") -> Path:
    return analysis_dir() / ("notebooks/.kbcache" if stage == "notebooks" else ".kbcache")


def assert_under_results(path: Path, label: str = "output") -> Path:
    rp, base = Path(path).resolve(), results().resolve()
    if not str(rp).startswith(str(base)):
        raise ValueError(f"{label} must live under {base}, got {rp}")
    return rp


# ---------------------------------------------------------------------------
# diagnostics
# ---------------------------------------------------------------------------
def describe() -> dict:
    out = {"analysis_dir": str(analysis_dir()), "python": python_exe(),
           "config_file": config_file() or "(none)"}
    for name in _ROOTS:
        try:
            p = _resolve(name, required=False)
            out[name] = str(p) if p else "(not found)"
        except MissingRoot as exc:
            out[name] = f"(error: {exc})"
    return out


def doctor() -> tuple[int, str]:
    """Human-readable environment report. Returns (exit_code, text)."""
    L = ["cma environment", ""]
    L.append(f"  repo             {analysis_dir()}")
    L.append(f"  python           {python_exe()}")
    L.append(f"  python version   {sys.version.split()[0]}")
    L.append(f"  config file      {config_file() or '(none; optional)'}")
    L.append("")

    rc = 0
    missing_data = []
    L.append("data roots")
    for name, (env, key, _cands, what, _marker) in _ROOTS.items():
        try:
            p = _resolve(name, required=False)
        except MissingRoot as exc:
            p = None
            L.append(f"  {name:<16} ERROR  {exc}")
            rc = 1
            continue
        if p:
            L.append(f"  {name:<16} ok     {p}")
        else:
            L.append(f"  {name:<16} absent -- {what}")
            L.append(f"  {'':<16}        set {env}, or paths.{key} in cma.toml")
            if name != "msdb_snapshot":     # only the grading scripts need that one
                missing_data.append(name)
    L.append("")

    L.append("python packages")
    for mod, why in (("cobra", "required: FBA"), ("pandas", "required: tables"),
                     ("numpy", "required"), ("matplotlib", "figures"),
                     ("plotly", "interactive figures"), ("pytest", "tests"),
                     ("yaml", "heuristic .yaml configs"),
                     ("kbutillib", "optional: notebook caching")):
        try:
            __import__(mod)
            L.append(f"  {mod:<16} ok     ({why})")
        except ImportError:
            hard = why.startswith("required")
            L.append(f"  {mod:<16} {'MISSING' if hard else 'absent '} ({why})")
            if hard:
                rc = 1
    L.append("")

    L.append("external tools")
    for tool, why in (("git", "required: reads MSDB releases with `git show`"),):
        found = shutil.which(tool)
        L.append(f"  {tool:<16} {'ok     ' + found if found else 'MISSING'} ({why})")
        if not found:
            rc = 1
    L.append("")

    if rc:
        L.append("Python requirements are missing. Run:")
        L.append("  python3 -m pip install -e '.[all]'")
    if missing_data:
        rc = rc or 1
        L.append(f"Data is missing ({', '.join(missing_data)}), so most of the pipeline")
        L.append("cannot run yet. Fetch it with:")
        L.append("  python3 scripts/fetch_data.py --all")
        L.append("or point at existing copies with the environment variables above.")
    if not rc:
        L.append("Everything required is present.")
    return rc, "\n".join(L)
