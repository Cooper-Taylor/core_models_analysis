#!/usr/bin/env python3
"""Validate the cma registry: parity, conformance, and per-kind pre-flight.

    python3 scripts/cma_check.py                      # everything
    python3 scripts/cma_check.py --parity             # registry == frozen literals
    python3 scripts/cma_check.py --manifest           # thermo_variants/manifest.json
    python3 scripts/cma_check.py --kind modelset --key ms2_gsm   # load 3 models, report

Exit codes: 0 pass, 1 a real failure, 2 a data/reference mismatch that means
"refuse to judge" rather than "this is broken".
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from cma import compat, directions, figures, manifest, models, paths, probes, sources, variants  # noqa: E402


def _head(title: str) -> None:
    print(f"\n=== {title} ===")


def check_parity() -> int:
    _head("parity: registry projections vs frozen pre-refactor literals")
    problems = compat.assert_all()
    if not problems:
        print(f"  OK  {len(compat.PROJECTIONS) + 1} projections reproduce their literals")
        return 0
    for p in problems:
        print(f"  FAIL {p}")
    return 1


def check_registries() -> int:
    _head("registries")
    rc = 0
    print(f"  thermo sources : {len(sources.all())} "
          f"({len(sources.all('panel'))} in the panel layer) {sources.keys()}")
    v_all = variants.all()
    print(f"  variants       : {len(v_all)} = "
          f"{len(variants.all('cascade'))} cascade, {len(variants.all('overlay'))} overlay, "
          f"{len(variants.all('derived'))} derived")
    print(f"  direction maps : {len(directions.all())} {directions.keys()}")
    print(f"  model sets     : {len(models.all())} {models.keys()}")
    print(f"  panels         : {[p.name for p in models.panels()]}")
    print(f"  probes         : {len(probes.all())} {probes.keys()}")
    print(f"  figure sets    : {len(figures.all())}")

    # cross-references must resolve
    for v in v_all:
        if v.kind == "overlay" and v.direction_map not in directions.keys():
            print(f"  FAIL variant {v.tag!r} names direction map {v.direction_map!r}, "
                  f"which is not registered")
            rc = 1
    for ms in models.all():
        try:
            models.media(ms.media)
        except Exception as exc:
            print(f"  FAIL model set {ms.name!r} names media {ms.media!r}: {exc}")
            rc = 1
    for p in models.panels():
        if p.model_set not in models.keys():
            print(f"  FAIL panel {p.name!r} names unknown model set {p.model_set!r}")
            rc = 1
    if rc == 0:
        print("  OK  every cross-reference resolves")
    return rc


def check_manifest() -> int:
    _head("thermo_variants/manifest.json")
    problems = manifest.validate()
    data = manifest.read()
    print(f"  {len(data.get('variants', []))} entries on disk")
    if not problems:
        print("  OK  manifest, registry and directories agree")
        return 0
    for p in problems:
        print(f"  WARN {p}")
    return 0  # advisory: the manifest can legitimately lead the registry


def check_modelset(key: str, n: int = 3) -> int:
    _head(f"model set {key!r}")
    try:
        ms = models.get(key)
    except Exception as exc:
        print(f"  FAIL {exc}")
        return 1
    d = ms.models_dir()
    print(f"  dir        : {d}  (exists: {d.exists()})")
    if not d.exists():
        print("  FAIL models_dir does not exist")
        return 1
    ids = ms.model_ids()
    print(f"  models     : {len(ids)}  glob={ms.id_glob!r} loader={ms.loader!r}")
    print(f"  media      : {ms.media}")
    print(f"  namespace  : {ms.results_namespace or '(bare results/)'}")
    if not ids:
        print("  FAIL no models matched the glob")
        return 1

    med = models.media(ms.media)
    try:
        table = med.loader()
        print(f"  media cpds : {len(table)} (uptake {med.default_uptake})")
    except Exception as exc:
        print(f"  FAIL media {ms.media!r} will not load: {exc}")
        return 1

    rc = 0
    for mid in ids[:n]:
        try:
            m = models.load_model(ms, mid)
        except Exception as exc:
            print(f"  FAIL {mid}: load failed: {type(exc).__name__}: {exc}")
            rc = 1
            continue
        seeds = [models.seed_id(r, ms.seed_key) for r in m.reactions]
        n_seed = sum(1 for s in seeds if s)
        bio = models.find_biomass(m, ms.biomass)
        n_big_bio = sum(1 for r in m.reactions
                        if r.id.lower().startswith(ms.biomass.fallback_prefix)
                        and len(r.metabolites) >= ms.biomass.real_biomass_min_mets)
        n_small_bio = sum(1 for r in m.reactions
                          if r.id.lower().startswith(ms.biomass.fallback_prefix)
                          and len(r.metabolites) < ms.biomass.real_biomass_min_mets)
        opened = models.apply_media(m, med)
        cpts = sorted({mm.compartment for mm in m.metabolites})
        print(f"  {mid[:38]:<38} rxn={len(m.reactions):>5} seed={n_seed:>5} "
              f"bio={bio.id if bio else 'NONE':<6} EX_open={opened:>4} cpts={cpts}")
        if bio is None:
            print(f"    FAIL no biomass reaction found under {ms.biomass}")
            rc = 1
        if n_seed == 0 and ms.seed_key.require_nonzero_overrides:
            print("    FAIL zero reactions resolved to a SEED id -- a sweep would "
                  "silently override nothing")
            rc = 1
        if set(cpts) - set(ms.compartments):
            print(f"    WARN compartments {sorted(set(cpts) - set(ms.compartments))} "
                  f"not declared in the model set")
        if ms.biomass.atp_probe is None and n_small_bio == 0 and n_big_bio > 0:
            print(f"    WARN flux_loops would close all {n_big_bio} bio* reactions and "
                  "leave no ATP probe open -> zero energy-generating cycles, silently. "
                  "Declare BiomassPolicy.atp_probe or inject one.")
    if rc == 0:
        print("  OK")
    return rc


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--parity", action="store_true")
    ap.add_argument("--registries", action="store_true")
    ap.add_argument("--manifest", action="store_true")
    ap.add_argument("--kind", choices=("modelset",))
    ap.add_argument("--key")
    ap.add_argument("-n", type=int, default=3, help="models to probe for --kind modelset")
    args = ap.parse_args()

    if args.kind == "modelset":
        return check_modelset(args.key or models.active_name(), args.n)

    selected = args.parity or args.registries or args.manifest
    rc = 0
    if not selected or args.registries:
        rc |= check_registries()
    if not selected or args.parity:
        rc |= check_parity()
    if not selected or args.manifest:
        rc |= check_manifest()
    print()
    print("RESULT:", "pass" if rc == 0 else "FAIL")
    return rc


if __name__ == "__main__":
    sys.exit(main())
