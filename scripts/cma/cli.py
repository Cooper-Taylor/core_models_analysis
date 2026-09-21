"""Command-line entry points.

``beginPipeline``  the pipeline itself (:func:`main`)
``cma-check``      validate the registries (:func:`main_check`)
``cma-new``        scaffold a new registry entry (:func:`main_new`)

All three are installed as console scripts by ``pip install -e .``; the
``beginPipeline`` file at the repo root is a shim so the command also works
from a plain checkout with nothing installed.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

from . import directions as directions_reg
from . import models as models_reg
from . import paths, pipeline, resolve
from . import sources as sources_reg
from . import variants as variants_reg


# ---------------------------------------------------------------------------
def cmd_list() -> int:
    print("STAGES  (--only NAME[,NAME...]   default: all)\n")
    width = max(len(s.name) for s in pipeline.all_stages())
    for s in pipeline.all_stages():
        flags = []
        if s.heavy:
            flags.append("heavy")
        if s.needs:
            flags.append("needs " + ",".join(s.needs))
        if s.requires:
            flags.append("requires --" + "/--".join(s.requires))
        if s.requires_any:
            flags.append("requires one of --" + "/--".join(s.requires_any))
        tail = f"   [{'; '.join(flags)}]" if flags else ""
        print(f"  {s.name:<{width}}  {s.description}{tail}")
    tags = sorted({t for s in pipeline.all_stages() for t in s.tags})
    print(f"\n  select a group with --only tag:NAME   tags: {', '.join(tags)}")

    print("\nREGISTERED MODEL SETS  (--models)\n")
    for ms in models_reg.all():
        print(f"  {ms.name:<14} {ms.description}")
    print("\nREGISTERED DIRECTION MAPS  (--directions)\n")
    print("  " + ", ".join(directions_reg.keys()))
    print("\nREGISTERED HEURISTIC VARIANTS  (--heuristics)\n")
    print("  " + ", ".join(variants_reg.tags("cascade")))
    print("\nREGISTERED THERMO SOURCES  (--thermo)\n")
    print("  " + ", ".join(sources_reg.keys()))
    print("\nREGISTERED MEDIA  (--media)\n")
    print("  " + ", ".join(m.name for m in models_reg._media.all()))
    print("\nHEURISTIC KNOBS  (for --heuristics 'k=v,...' or a .yaml/.json/.toml file)\n")
    for k, t in resolve.heuristic_knobs().items():
        print(f"  {k:<26} {t}")
    return 0


def cmd_describe(name: str) -> int:
    s = pipeline.get(name)
    print(f"{s.name}\n{'-' * len(s.name)}\n{s.description}\n")
    print(f"  depends on : {', '.join(s.needs) or '(nothing)'}")
    print(f"  requires   : {', '.join('--' + r for r in s.requires) or '(no extra inputs)'}")
    if s.requires_any:
        print(f"  needs one of: {', '.join('--' + r for r in s.requires_any)}")
    print(f"  heavy      : {s.heavy}")
    print(f"  tags       : {', '.join(s.tags) or '(none)'}")
    print(f"  writes     : <out>/{s.name}.json")
    if s.fn.__doc__:
        print("\n" + s.fn.__doc__.strip())
    return 0


# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="beginPipeline", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    g = p.add_argument_group("inputs (each accepts a PATH or a registered KEY; repeatable)")
    g.add_argument("--models", "-m", action="append", metavar="PATH|KEY",
                   help="directory of models, a glob, a file listing ids, or a registered "
                        "model set. Repeat to pool several.")
    g.add_argument("--directions", "-d", action="append", default=[], metavar="PATH|KEY",
                   help="a {rxn_id: > < = ?} table: .json, .csv/.tsv, an MSDB report .txt, "
                        "or a registered direction map.")
    g.add_argument("--heuristics", "-H", action="append", default=[], metavar="PATH|TAG|KNOBS",
                   help="a registered cascade variant, a .yaml/.json/.toml of "
                        "ReversibilityConfig knobs, or inline 'k=v,k=v'.")
    g.add_argument("--thermo", "-t", action="append", default=[], metavar="PATH|KEY",
                   help="a ΔG′° table (rxn_id, dg, dge) or a registered thermo source.")
    g.add_argument("--panel", metavar="PATH|NAME",
                   help="restrict the model set to the ids in this file (one per line), "
                        "or a registered panel name. Combines with any --models.")
    g.add_argument("--media", metavar="PATH|NAME",
                   help="a .cpd/.json/.tsv medium or a registered media name. "
                        "Default: the model set's own.")

    g = p.add_argument_group("what to produce")
    g.add_argument("--only", "-o", metavar="A,B",
                   help="comma-separated stage names, or tag:NAME. Default: all.")
    g.add_argument("--exclude", "-x", metavar="A,B", help="stages to skip.")
    g.add_argument("--no-heavy", action="store_true",
                   help="skip stages marked heavy (the FBA sweeps).")

    g = p.add_argument_group("run control")
    g.add_argument("--out", type=Path, metavar="DIR",
                   help="output directory. Default: results/pipeline/<model set>.")
    g.add_argument("--limit", type=int, metavar="N", help="use only the first N models.")
    g.add_argument("--workers", "-j", type=int, default=max(1, min(16, (os.cpu_count() or 4) - 2)))
    g.add_argument("--baseline", metavar="LABEL",
                   help="direction set that growth_diff compares against. Default: on_disk.")
    g.add_argument("--db-level", default="EQ", choices=("EQ", "GC", ""),
                   help="MSDB energy level for the heuristic cascade. Default: EQ.")
    g.add_argument("--option", action="append", default=[], metavar="K=V",
                   help="extra stage option, e.g. --option top_n=50")
    g.add_argument("--stop-on-error", action="store_true",
                   help="halt at the first failing stage instead of continuing.")

    g = p.add_argument_group("discovery")
    g.add_argument("--list", "-l", action="store_true",
                   help="list stages and everything registered, then exit.")
    g.add_argument("--describe", metavar="STAGE", help="explain one stage, then exit.")
    g.add_argument("--dry-run", action="store_true",
                   help="resolve inputs and print the plan without running anything.")
    g.add_argument("--doctor", action="store_true",
                   help="report what is installed, what data roots were found, "
                        "and what is missing; then exit.")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    if args.doctor:
        rc, text = paths.doctor()
        print(text)
        return rc
    if args.list:
        return cmd_list()
    if args.describe:
        return cmd_describe(args.describe)
    if not args.models:
        build_parser().print_usage()
        print("\nbeginPipeline: --models is required. "
              "Try `python3 beginPipeline --list` to see what is available.")
        return 2

    # ---- resolve every input -------------------------------------------
    try:
        model_sets = [resolve.resolve_model_set(s, media=args.media) for s in args.models]
        ms = model_sets[0]
        ids = []
        for s in model_sets:
            ids.extend(resolve.model_ids(s))

        if args.panel:
            pth = Path(args.panel).expanduser()
            if pth.exists():
                wanted = [ln.strip() for ln in pth.read_text().splitlines()
                          if ln.strip() and not ln.startswith("#")]
            else:
                wanted = models_reg.panel(args.panel).ids()
            available = set(ids)
            keep = [w for w in wanted if w in available]
            missing = [w for w in wanted if w not in available]
            if not keep:
                print(f"beginPipeline: --panel {args.panel} selected none of the "
                      f"{len(ids)} models in {ms.name!r}. First panel id: "
                      f"{wanted[0] if wanted else '(empty)'}")
                return 2
            if missing:
                print(f"note: --panel listed {len(missing)} id(s) absent from "
                      f"{ms.name!r}, e.g. {missing[:3]}")
            ids = keep
        if len(model_sets) > 1:
            print(f"note: pooling {len(model_sets)} model sets under the loader of "
                  f"{ms.name!r}; pass one --models if they differ in format.")

        media = (resolve.resolve_media(args.media) if args.media
                 else models_reg.media(ms.media))

        raw_dirs = [resolve.resolve_directions(spec) for spec in args.directions]
        dir_labels = resolve.disambiguate([lab for lab, _ in raw_dirs], args.directions)
        dmaps = {lab: ops for lab, (_, ops) in zip(dir_labels, raw_dirs, strict=True)}

        raw_heur = [resolve.resolve_heuristic(spec) for spec in args.heuristics]
        heur_labels = resolve.disambiguate([lab for lab, _ in raw_heur], args.heuristics)
        heur = {lab: cfg for lab, (_, cfg) in zip(heur_labels, raw_heur, strict=True)}

        raw_thermo = [resolve.resolve_thermo(spec) for spec in args.thermo]
        thermo_labels = resolve.disambiguate([lab for lab, _ in raw_thermo], args.thermo)
        thermo = {lab: t for lab, (_, t) in zip(thermo_labels, raw_thermo, strict=True)}
    except resolve.ResolveError as exc:
        print(f"beginPipeline: {exc}")
        return 2

    options = {}
    for kv in args.option:
        k, _, v = kv.partition("=")
        try:
            options[k.strip()] = int(v)
        except ValueError:
            options[k.strip()] = v
    options["db_level"] = args.db_level

    out_dir = args.out or paths.results("pipeline", ms.name)

    ctx = pipeline.Ctx(
        model_set=ms, model_ids=ids, media=media, direction_maps=dmaps,
        heuristics=heur, thermo=thermo, out_dir=Path(out_dir), workers=args.workers,
        limit=args.limit, baseline=args.baseline, options=options,
    )

    try:
        chosen = pipeline.select(
            [s.strip() for s in args.only.split(",")] if args.only else None,
            [s.strip() for s in args.exclude.split(",")] if args.exclude else None,
            include_heavy=not args.no_heavy,
        )
    except Exception as exc:
        print(f"beginPipeline: {exc}")
        print("\nRun `python3 beginPipeline --list` to see the available stages.")
        return 2

    n_sel = len(ids[: args.limit] if args.limit else ids)
    print(f"model set   : {ms.name}  ({n_sel} of {len(ids)} models, {ms.loader})")
    print(f"media       : {media.name}")
    print(f"directions  : {', '.join(dmaps) or '(none — on-disk bounds only)'}")
    print(f"heuristics  : {', '.join(heur) or '(none)'}")
    print(f"thermo      : {', '.join(thermo) or '(none)'}")
    print(f"stages      : {', '.join(s.name for s in chosen)}")
    print(f"output      : {out_dir}")

    if args.dry_run:
        print("\n--dry-run: inputs resolved, nothing executed.")
        return 0

    print()
    t0 = time.time()
    summary = pipeline.run(chosen, ctx, keep_going=not args.stop_on_error)
    summary["elapsed_s"] = round(time.time() - t0, 2)
    # Everything needed to reproduce the run. The first version omitted
    # --panel/--only/--exclude/--baseline/--no-heavy/--out, so run.json
    # described a different run from the one that happened.
    summary["inputs"] = {
        "models": args.models, "directions": args.directions,
        "heuristics": args.heuristics, "thermo": args.thermo,
        "media": args.media, "panel": args.panel, "limit": args.limit,
        "only": args.only, "exclude": args.exclude, "no_heavy": args.no_heavy,
        "baseline": args.baseline, "db_level": args.db_level,
        "workers": args.workers, "option": args.option, "out": str(out_dir),
        "stop_on_error": args.stop_on_error,
    }
    summary["resolved"] = {
        "model_set": ms.name, "models_dir": str(ms.models_dir()),
        "loader": ms.loader, "id_glob": ms.id_glob,
        "n_models": n_sel, "n_models_available": len(ids),
        "media": media.name, "media_default_uptake": media.default_uptake,
        "direction_labels": list(dmaps), "heuristic_labels": list(heur),
        "thermo_labels": list(thermo),
        "stages": [st.name for st in chosen],
    }
    summary["command"] = "python3 beginPipeline " + " ".join(
        __import__("shlex").quote(a) for a in (sys.argv[1:] or []))
    summary["paths"] = paths.describe()
    (Path(out_dir) / "run.json").write_text(json.dumps(summary, indent=2, default=str))

    ok = sum(1 for v in summary["stages"].values() if v["status"] == "ok")
    failed = [k for k, v in summary["stages"].items() if v["status"] == "failed"]
    skipped = [k for k, v in summary["stages"].items() if v["status"] == "skipped"]
    print(f"\n{ok} ok, {len(failed)} failed, {len(skipped)} skipped "
          f"in {summary['elapsed_s']}s -> {out_dir}")
    if failed:
        print(f"failed: {', '.join(failed)}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())


def main_check(argv=None) -> int:
    """``cma-check`` -- validate registries, parity and per-kind conformance."""
    import runpy

    target = paths.scripts("cma_check.py")
    if argv is not None:
        sys.argv = ["cma_check.py", *argv]
    try:
        runpy.run_path(str(target), run_name="__main__")
    except SystemExit as exc:
        return int(exc.code or 0)
    return 0


def main_new(argv=None) -> int:
    """``cma-new`` -- scaffold a registry entry."""
    import runpy

    target = paths.scripts("cma_new.py")
    if argv is not None:
        sys.argv = ["cma_new.py", *argv]
    try:
        runpy.run_path(str(target), run_name="__main__")
    except SystemExit as exc:
        return int(exc.code or 0)
    return 0
