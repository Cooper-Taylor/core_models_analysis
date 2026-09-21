#!/usr/bin/env python3
"""Scaffold a new registry entry: the "just tell me what to write" command.

    python3 scripts/cma_new.py thermo-source --key eq3 --display "eQuilibrator 3.0"
    python3 scripts/cma_new.py variant  --key my_variant --kind cascade
    python3 scripts/cma_new.py variant  --key llm_gpt5  --kind overlay
    python3 scripts/cma_new.py model-set --key ecoli_gsm
    python3 scripts/cma_new.py direction-map --key tecrdb
    python3 scripts/cma_new.py probe    --key atp_yield

Prints a valid, commented entry to paste into the matching file under
``scripts/cma/entries/``, plus the command that validates it. Pass ``--append``
to write it there directly.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from cma import paths  # noqa: E402

TARGET = {
    "thermo-source": "thermo_sources.py",
    "variant": "variants.py",
    "direction-map": "direction_maps.py",
    "model-set": "model_sets.py",
    "probe": "probes.py",
    "figure-set": "figure_sets.py",
}


def thermo_source(key: str, display: str) -> str:
    return f'''
register(ThermoSource(
    key={key!r},
    display={display!r},
    # The EXACT key inside MSDB's per-reaction `thermodynamics` dict, or None
    # for a source that does not live in MSDB.
    label={display!r},
    # Every spelling this source may be referred to by. Add one here rather
    # than inventing it at a call site -- `sources.resolve()` accepts them all.
    slugs={{"short": {key!r}, "upper": {key.upper()!r},
           "long": {key!r}, "dash": {key!r}}},
    # msdb_sublist | msdb_toplevel_gated | external_energy
    #   | external_operators | derived
    kind="msdb_sublist",
    db_level=None,
    color="#888888",
    color_families=("bars",),
    axis_title="{display} ΔG′° (kcal/mol)",
    # ("panel", "all_models") puts it in the sweeps; () registers it without
    # promoting it, which is how dgpms is parked today.
    layers=("panel", "all_models"),
    section="§ New",
    citations=(),
))
'''


def variant(key: str, kind: str) -> str:
    if kind == "cascade":
        return f'''
def _{key}_cfg():
    # Any knob on ReversibilityConfig. Defaults reproduce MSDB exactly, so
    # change only what the variant is about.
    return lib.ReversibilityConfig(sigma_band_k=1.96)


register(Variant(
    tag={key!r},
    kind="cascade",
    cfg=_{key}_cfg,
    title="One-line legacy label",
    apt_title="A descriptive title shown to website users",
    description=("Two to four sentences: what the baseline does, what this "
                 "variant changes, and why. Self-contained for a reader who "
                 "has not read the heuristics review."),
    citations=(),
    section="§ New",
    order=900,
    in_notebook=True,
))
'''
    return f'''
register(Variant(
    tag={key!r},
    kind="overlay",
    # A key registered in cma/entries/direction_maps.py.
    direction_map={key!r},
    overlay_on="baseline",
    status_prefix={key!r},
    title="One-line legacy label",
    apt_title="A descriptive title shown to website users",
    description="What produced these directions and why they differ from baseline.",
    citations=(),
    section="§ New — external direction map",
    order=900,
))
'''


def direction_map(key: str) -> str:
    return f'''
register(DirectionMapSpec(
    key={key!r},
    # LAZY: register a loader, never a loaded map.
    loader=_json_map(_R("rxn_directions_{key}.json")),
    # 'complete' = an opinion (possibly '?') for every MSDB reaction.
    # 'partial'  = an absent key means NO opinion; the consumer falls back.
    coverage="partial",
    source=None,
    description="Where these directions came from.",
    artifacts={{"json": _R("rxn_directions_{key}.json")}},
))
'''


def model_set(key: str) -> str:
    return f'''
register(ModelSet(
    name={key!r},
    models_dir=lambda: Path(os.environ.get(
        "{key.upper()}_DIR", "/scratch/ctaylor/{key}")),
    id_glob="*.json",                 # "*.json.gz" if gzipped
    loader="cobra_json",              # cobra_json | cobra_json_gz | module:fn
    media="kbase_complete",
    biomass=BiomassPolicy(
        prefer=("bio1",), then=("bio2",),
        real_biomass_min_mets=10,
        # None means there is NO maintenance reaction to use as an ATP probe.
        # flux_loops then closes every biomass and silently finds zero cycles.
        atp_probe=None,
    ),
    seed_key=SeedKeyPolicy(
        annotation="seed.reaction",
        strip_compartment_suffix=True,
        fallback="id_prefix",
        require_nonzero_overrides=True,
    ),
    compartments=("c0", "e0"),
    cytosol="c0",
    id_to_accession=lambda mid: mid,
    # None writes to the bare results/ tree (only the core set should).
    results_namespace={key!r},
    default_panel=None,
    description="What this collection is.",
))
'''


def probe(key: str) -> str:
    return f'''
def _{key}_one(model_id, ctx):
    # ctx carries: direction_map, baseline_map, media, biomass, flux()
    gh = _gh()
    return gh.fba_one(model_id, reversibility_map=ctx.get("direction_map"))


register(Probe(
    name={key!r},
    fn=_{key}_one,
    needs_media=True,
    needs_biomass=True,
    # True if the probe parallelizes internally (cobra's FVA does): the runner
    # then refuses to nest it inside its own pool.
    internal_pool=False,
    out_template="site/data/{{ns}}panel_{key}.json",
    description="What this measures.",
))
'''


BUILDERS = {
    "thermo-source": lambda a: thermo_source(a.key, a.display or a.key),
    "variant": lambda a: variant(a.key, a.kind),
    "direction-map": lambda a: direction_map(a.key),
    "model-set": lambda a: model_set(a.key),
    "probe": lambda a: probe(a.key),
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("entry_kind", metavar="kind", choices=sorted(BUILDERS))
    ap.add_argument("--key", required=True)
    ap.add_argument("--display")
    ap.add_argument("--kind", dest="variant_kind", default="cascade",
                    choices=("cascade", "overlay", "derived"),
                    help="for `variant` only: which flavour of variant to scaffold")
    ap.add_argument("--append", action="store_true",
                    help="write the stub into the entries file instead of printing it")
    args = ap.parse_args()

    entry_kind = args.entry_kind
    stub = BUILDERS[entry_kind](
        argparse.Namespace(key=args.key, display=args.display, kind=args.variant_kind))

    target = paths.scripts("cma", "entries", TARGET[entry_kind])
    if args.append:
        with target.open("a") as fh:
            fh.write("\n" + stub.strip() + "\n")
        print(f"appended to {target}")
    else:
        print(f"# --- paste into {target} ---")
        print(stub.strip())
    print()
    print("# then validate:")
    if entry_kind == "model-set":
        print(f"#   python3 scripts/cma_check.py --kind modelset --key {args.key}")
    else:
        print("#   python3 scripts/cma_check.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
