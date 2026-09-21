"""The pipeline: a registry of output stages, and the runner that walks them.

Each thing the pipeline can produce is a :class:`Stage` -- a name, a one-line
description, its dependencies, and a function. ``beginPipeline`` selects stages
by name (``--only growth,loops``), and the default is every stage. Adding a new
output is one ``@stage`` decorator, not a new driver script.

A stage receives a :class:`Ctx` carrying the resolved inputs (model set,
direction maps, heuristics, thermo sources, media) and returns a JSON-able
dict. The runner handles dependency order, timing, parallelism, writing
``<out>/<stage>.json``, and the run manifest.
"""

from __future__ import annotations

import json
import time
import traceback
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

from .registry import Registry


# ---------------------------------------------------------------------------
@dataclass
class Ctx:
    """Everything a stage may read. Populated by ``beginPipeline``."""

    model_set: object = None
    model_ids: list = field(default_factory=list)
    media: object = None
    #: ``{label: {rxn_id: operator}}`` -- one entry per --directions argument
    direction_maps: dict = field(default_factory=dict)
    #: ``{label: ReversibilityConfig}`` -- one per --heuristics argument
    heuristics: dict = field(default_factory=dict)
    #: ``{label: {rxn_id: (dg, dge)}}`` -- one per --thermo argument
    thermo: dict = field(default_factory=dict)
    out_dir: Path = Path(".")
    workers: int = 8
    limit: int | None = None
    baseline: str | None = None
    #: results of already-completed stages, keyed by stage name
    results: dict = field(default_factory=dict)
    options: dict = field(default_factory=dict)

    def log(self, msg: str) -> None:
        print(f"    {msg}", flush=True)


@dataclass(frozen=True)
class Stage:
    name: str
    fn: Callable
    description: str = ""
    needs: tuple[str, ...] = ()
    #: inputs without which the stage is skipped rather than failed (ALL needed)
    requires: tuple[str, ...] = ()
    #: skip unless at least ONE of these inputs is present
    requires_any: tuple[str, ...] = ()
    heavy: bool = False
    tags: tuple[str, ...] = ()

    @property
    def key(self) -> str:
        return self.name


_reg: Registry[Stage] = Registry("stage", key_attr="name")


def stage(name: str, description: str = "", needs: tuple = (),
          requires: tuple = (), requires_any: tuple = (),
          heavy: bool = False, tags: tuple = ()):
    """Register a pipeline stage."""

    def deco(fn):
        _reg.register(Stage(name=name, fn=fn, description=description,
                            needs=tuple(needs), requires=tuple(requires),
                            requires_any=tuple(requires_any),
                            heavy=heavy, tags=tuple(tags)))
        return fn

    return deco


def get(name: str) -> Stage:
    _load()
    return _reg.get(name)


def all_stages() -> list[Stage]:
    _load()
    return _reg.all()


def names() -> list[str]:
    _load()
    return _reg.keys()


_loaded = False


def _load() -> None:
    global _loaded
    if not _loaded:
        _loaded = True
        from . import stages  # noqa: F401  (registers via @stage)


# ---------------------------------------------------------------------------
def select(requested: list[str] | None, exclude: list[str] | None = None,
           include_heavy: bool = True) -> list[Stage]:
    """Resolve ``--only`` / ``--exclude`` / ``--no-heavy`` into an ordered list.

    Exclusion is applied **after** dependency expansion, not before. The first
    version expanded dependencies last, so ``--no-heavy`` and
    ``--exclude growth`` were silently undone: ``growth_diff`` pulled ``growth``
    back in and the heavy FBA sweep ran anyway, which is the opposite of what
    the flag says.

    A stage whose dependency has been excluded is dropped too, with a note,
    rather than running against missing inputs.
    """
    _load()
    banned = set(exclude or ())
    if not include_heavy:
        banned |= {s.name for s in all_stages() if s.heavy}

    if not requested or "all" in requested:
        chosen = [s.name for s in all_stages()]
    else:
        chosen = []
        for r in requested:
            if r.startswith("tag:"):
                tag = r[4:]
                hit = [s.name for s in all_stages() if tag in s.tags]
                if not hit:
                    known = sorted({t for s in all_stages() for t in s.tags})
                    raise ValueError(
                        f"no stage carries the tag {tag!r}. Known tags: {', '.join(known)}"
                    )
                chosen += hit
            else:
                get(r)          # raises, listing the valid names
                chosen.append(r)

    # expand dependencies first ...
    wanted, queue = set(), list(chosen)
    while queue:
        n = queue.pop()
        if n in wanted:
            continue
        wanted.add(n)
        queue.extend(get(n).needs)

    # ... then apply exclusion, and drop anything orphaned by it
    wanted -= banned
    changed = True
    while changed:
        changed = False
        for n in sorted(wanted):
            missing_dep = [d for d in get(n).needs if d not in wanted]
            if missing_dep:
                print(f"note: skipping {n!r} -- it depends on "
                      f"{', '.join(missing_dep)}, which is excluded")
                wanted.discard(n)
                changed = True

    ordered, seen = [], set()

    def visit(n: str, trail=()):
        if n in seen:
            return
        if n in trail:
            raise ValueError(f"circular stage dependency: {' -> '.join((*trail, n))}")
        for d in get(n).needs:
            if d in wanted:
                visit(d, (*trail, n))
        seen.add(n)
        ordered.append(get(n))

    for s in all_stages():
        if s.name in wanted:
            visit(s.name)
    return ordered


def _missing(stage_obj: Stage, ctx: Ctx) -> list[str]:
    """Inputs that are absent, phrased the way the user would supply them."""
    out = [req for req in stage_obj.requires if not getattr(ctx, req, None)]
    if stage_obj.requires_any:
        if not any(getattr(ctx, r, None) for r in stage_obj.requires_any):
            out.append(" or ".join(stage_obj.requires_any))
    return out


def run(stages_to_run: list[Stage], ctx: Ctx, keep_going: bool = True) -> dict:
    """Execute stages in order, writing each result under ``ctx.out_dir``."""
    ctx.out_dir.mkdir(parents=True, exist_ok=True)
    summary = {"stages": {}, "order": [s.name for s in stages_to_run]}

    for s in stages_to_run:
        missing = _missing(s, ctx)
        if missing:
            print(f"[skip] {s.name}: needs {', '.join(missing)}", flush=True)
            summary["stages"][s.name] = {"status": "skipped", "missing": missing}
            continue
        print(f"[run ] {s.name}  -- {s.description}", flush=True)
        t0 = time.time()
        try:
            result = s.fn(ctx)
            elapsed = round(time.time() - t0, 2)
            ctx.results[s.name] = result
            path = ctx.out_dir / f"{s.name}.json"
            path.write_text(json.dumps(result, indent=2, default=str))
            n = result.get("_n") if isinstance(result, dict) else None
            print(f"[ok  ] {s.name}  {elapsed}s -> {path.name}"
                  + (f"  ({n} rows)" if n is not None else ""), flush=True)
            summary["stages"][s.name] = {"status": "ok", "elapsed_s": elapsed,
                                         "output": path.name}
        except Exception as exc:  # noqa: BLE001
            elapsed = round(time.time() - t0, 2)
            tb = traceback.format_exc(limit=6)
            print(f"[FAIL] {s.name}  {elapsed}s  {type(exc).__name__}: {exc}", flush=True)
            summary["stages"][s.name] = {"status": "failed", "elapsed_s": elapsed,
                                         "error": f"{type(exc).__name__}: {exc}",
                                         "traceback": tb}
            if not keep_going:
                break
    return summary
