# `cma` — the registry core

Adding a **thermo source**, a **heuristic variant**, a **direction map**, a
**model set**, a **figure** or a **probe** is one entry in
[`entries/`](entries/). Everything downstream reads the registry instead of
keeping its own copy of the literal.

```python
from cma import paths, sources, variants, directions, models, probes
```

## Why

| concept | before | now |
|---|---|---|
| a thermo source | re-declared in ~30 scripts under 6 incompatible spellings | one `ThermoSource` |
| a variant | `cfg` had to be a `ReversibilityConfig` factory, so 6 of 20 variants could not be registered at all | `Variant.kind` ∈ `cascade` / `overlay` / `derived` |
| a direction map | 5 on-disk shapes, no schema check anywhere | one `DirectionMap` with `validate()` and `content_hash` |
| a model set | `MODELS_DIR` constant in 19 files, 5 of them env-proof absolutes | one `ModelSet` |
| paths | 11 env vars, 6 idioms, `MSDB_ROOT` with two different defaults | `paths.py`, all functions |

## Add something

```bash
python3 scripts/cma_new.py thermo-source --key eq3 --display "eQuilibrator 3.0"
python3 scripts/cma_new.py variant   --key llm_gpt5 --kind overlay
python3 scripts/cma_new.py model-set --key ecoli_gsm
python3 scripts/cma_check.py                              # validate
python3 scripts/cma_check.py --kind modelset --key ms2_gsm
```

`cma_new.py` prints a commented, valid stub for the matching file under
`entries/`; `--append` writes it there.

## Guarantees

`cma_check.py --parity` asserts that every projection the registry offers is
identical to the literal the old code hardcoded, using the frozen copies in
[`compat.py`](compat.py). `variant_catalog.VARIANTS` is now a view over the
registry and was verified to produce the same 14 entries, in the same order,
with `ReversibilityConfig` objects equal field-for-field, and the same cascade
output across all 56,012 reactions.

`compat.py` is scaffolding with a declared end of life: once every consumer
reads the registry, the frozen copies are deleted.

## Layout

```
registry.py   ordered, duplicate-rejecting Registry[T]; loads a FIXED list of
              entry modules (never a directory scan — ten scripts/ modules
              call sys.exit at import)
kinds.py      the six frozen dataclasses. Read this first.
paths.py      every root, as functions. Splits MSDB_ROOT (live clone) from
              MSDB_SNAPSHOT_ROOT (pinned snapshot) — they were conflated, so
              exporting MSDB_ROOT silently repointed 9 scripts at other data.
sources.py    ThermoSource registry: resolve(), slug(), pairs(), palette()
variants.py   Variant registry + the legacy() projection
directions.py DirectionMap value object + registry
models.py     ModelSet / Media / Panel, load_model(), seed_id(), apply_media()
probes.py     Probe registry
figures.py    FigureSet registry (scripts/figures.tsv stays the source of truth)
manifest.py   thermo_variants/manifest.json: merge, never rewrite
compat.py     frozen pre-refactor literals + the parity assertion
entries/      the six files you actually edit
```

## Import weight

Importing `cma` reads no environment variable and imports neither cobra,
pandas, numpy nor matplotlib — `site/serve.py` is deliberately stdlib-only.
Submodules resolve lazily through `__getattr__`.

This package shadows the PyPI distribution `cma` (CMA-ES), which this project
does not use. Do not add it to `requirements.txt`.
