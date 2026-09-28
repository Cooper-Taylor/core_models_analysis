# `beginPipeline`

One command runs the whole reaction-direction analysis.

```bash
python3 beginPipeline --models <path|key> [--directions ...] [--heuristics ...] [--thermo ...]
```

Every input accepts **either a filesystem path or the key of something already
registered**, so you can point it at models, direction tables, heuristic
configs or ΔG tables that nothing knows about yet, without editing any Python.
All input flags repeat, and the pipeline runs the cross-product: every model is
solved under every direction set you supply.

By default it produces **everything it knows how to produce**. Narrow that with
`--only`.

Run it with an interpreter that has `cobra` and `modelseedpy`. On this machine
that is `/mnt/homes/ctaylor/conda/miniforge3/envs/core_models_analysis/bin/python`;
`python3` works if that environment is active.

```bash
python3 beginPipeline --list              # every stage + everything registered
python3 beginPipeline --describe growth   # explain one stage
python3 beginPipeline --dry-run --models /data/my_models    # resolve, run nothing
```

---

## The five inputs

### `--models` / `-m` — required

| you pass | what happens |
|---|---|
| `core_kegg2` | a registered model set |
| `/data/my_models` | a directory; the file type is detected and the loader chosen |
| `'/data/models/*.json.gz'` | a glob (quote it so the shell does not expand it) |
| `/data/one_model.json` | a single model |
| `ids.txt` | a file of model ids or paths, one per line |

Repeatable. `.json`, `.json.gz` and SBML are recognised. An ad-hoc directory
gets conservative defaults: the widest biomass precedence, and a SEED-key
policy that falls back to the reaction-id prefix, so models with no
`seed.reaction` annotation still match a direction map.

### `--directions` / `-d`

A `{reaction: > < = ?}` table. Accepts a registered direction-map key, a
`.json` object, a `.csv`/`.tsv` with a reaction column and a direction column,
or an MSDB-format report (`rxn \t status \t old \t new`). Operators are
normalised, so `forward` / `fwd` / `>` all mean the same thing, and the result
is schema-checked before anything runs.

### `--heuristics` / `-H`

Three forms:

```bash
--heuristics H4                              # a registered cascade variant
--heuristics my_knobs.yaml                   # a config file
--heuristics 'sigma_band_k=1.96,mm_band=3'   # inline
```

A config file is `.yaml`, `.json` or `.toml`:

```yaml
tag: sigma95
knobs:
  sigma_band_k: 1.96
  sigma_bounds_k: 1.96
```

Every key under `knobs` is a field of `ReversibilityConfig`; `--list` prints
them all with their types. Anything you leave out keeps the ModelSEED default,
so a spec states only what it changes.

Four knobs take a data table rather than a number:

| knob | accepts |
|---|---|
| `ln_ri_by_rxn` | `auto`, or a path to a MetaNetX-style energies table |
| `energy_override_by_rxn` | `auto`, or a path to `id, dG, uncertainty` |
| `per_met_conc_range` | the name of a table in `reversibility_lib` (e.g. `BENNETT_2009_ECOLI`), or a path |
| `per_met_conc` | the name of a table (e.g. `BENNETT_2009_MEAN`), or a path |

Worked examples are in [`examples/`](examples/).

Each heuristic is run through the cascade over all 56,012 MSDB reactions and
becomes another direction set, so heuristics and `--directions` are compared
side by side.

### `--thermo` / `-t`

A ΔG′° table: a registered source key, or a file of `rxn_id, dg, dge`.

### `--media`

A registered media name, or a `.cpd` / `.json` / `.tsv` file. Defaults to the
model set's own medium. This matters: the ModelSEED v2 genome-scale models were
gap-filled on a 20-compound defined medium, not the 347-compound complete one.

### `--panel`

Restrict the models to the ids in a file, or to a registered panel name. Works
with any `--models`.

---

## Choosing what comes out

| flag | effect |
|---|---|
| *(nothing)* | every stage — the default |
| `--only growth,growth_diff` | just those, plus whatever they depend on |
| `--only tag:inventory` | every stage carrying that tag |
| `--exclude growth` | everything but that, and anything depending on it |
| `--no-heavy` | skip the FBA sweeps |

Dependencies are pulled in automatically and ordered. A stage whose inputs are
absent is **skipped and reported**, not failed: asking for `thermo_summary`
without `--thermo` prints `[skip] thermo_summary: needs thermo`.

### Stages

| stage | produces |
|---|---|
| `inputs` | every resolved input, echoed back |
| `model_summary` | size distribution, compartments, SEED-id coverage |
| `direction_summary` | operator composition and pairwise agreement of the maps |
| `heuristics` | runs each heuristic's cascade over MSDB into a direction map |
| `thermo_summary` | coverage and ΔG′° distribution per source |
| `growth` | FBA biomass flux per model per direction set *(heavy)* |
| `growth_diff` | which models gain or lose growth, against a baseline |
| `direction_changes` | reactions whose direction differs between the supplied sets |
| `report` | `REPORT.md`, a human-readable summary of everything that ran |

---

### Other run control

| flag | effect |
|---|---|
| `--limit N` | use only the first N models. Applied after `--panel`. |
| `--workers N` / `-j` | parallel FBA workers. Default: CPUs minus two, capped at 16. |
| `--baseline LABEL` | which direction set `growth_diff` compares against. Must name a set in this run, or the stage fails rather than diffing against something else. Default `on_disk`. |
| `--db-level EQ\|GC` | MSDB energy level for the heuristic cascade. Default `EQ`. |
| `--option K=V` | extra stage option, e.g. `--option summary_sample=1000`, `--option top_n=50`. |
| `--stop-on-error` | halt at the first failing stage instead of continuing. |

### How a direction becomes a bound

| operator | bounds applied |
|---|---|
| `>` | `(0, 1000)` |
| `<` | `(-1000, 0)` |
| `=` | `(-1000, 1000)` |
| `?` | `(-1000, 1000)` — same as reversible |

`?` means *no opinion*, and it is applied permissively: a map full of `?`
**removes** constraints rather than adding them. Word and arrow spellings are
accepted (`forward`, `<=>`, `-->`); anything unrecognised becomes `?`, so check
`direction_summary`'s operator counts if a map behaves unexpectedly.

A direction set that matches **zero** reactions in your models is reported as a
warning during `growth`, because its numbers would otherwise be identical to
`on_disk` and look like a real result.

## Output

```
<out>/<stage>.json     one file per stage
<out>/REPORT.md        human-readable summary
<out>/run.json         what ran, how long, and every resolved input
```

Default `--out` is `results/pipeline/<model set>`, so two runs of the same
model set overwrite each other unless you pass `--out`. `run.json` records the
exact command, every flag, and everything they resolved to.

Exit codes: **0** all stages ok, **1** a stage failed, **2** an input could not
be resolved. A failing stage does not stop the others unless you pass
`--stop-on-error`.

---

### Media, and a trap worth knowing

If you point `--models` at a directory, the pipeline has no way to know what
medium those models were built for, so it falls back to the 347-compound
complete medium and says so. For models gap-filled on a defined medium that is
roughly a **200-fold** difference in growth flux. Pass `--media` explicitly:

```bash
--media /scratch/ctaylor/modelseed2_gs_models/gmm/media/Carbon-D-Glucose.cpd
```

A `.cpd` file has no flux column, so uptake is inferred: a list of more than
100 compounds is treated as a complete medium (`-1000`), a shorter one as a
defined medium (`-10`).

## Worked examples

```bash
# everything, on a registered set
python3 beginPipeline --models core_kegg2

# a new directory of models and two direction tables of your own
python3 beginPipeline --models /data/my_models \
                      --directions /data/gc.json --directions /data/eq.csv

# compare a registered heuristic against one of your own, on genome-scale
# models, using the medium they were gap-filled on
python3 beginPipeline --models ms2_gsm \
                      --media /scratch/ctaylor/modelseed2_gs_models/gmm/media/Carbon-D-Glucose.cpd \
                      --heuristics H4 --heuristics examples/heuristic_sigma95.yaml \
                      --limit 100

# just the cheap inventory, no FBA
python3 beginPipeline --models ms2_gsm --only tag:inventory

# the 100-model panel, against the group-contribution directions
python3 beginPipeline --models core_kegg2 --panel core_100 --directions thermo_gc
```

---

## The written analysis

The direction-influence study is typeset as a PDF from the same result files
the pipeline writes, so the document cannot drift from the data:

```bash
python3 scripts/regen_figures.py direction_influence   # the four figures
python3 scripts/build_influence_report.py              # PDF + Markdown
python3 scripts/build_influence_report.py --format md  # just one
```

Outputs `reports/DIRECTION_INFLUENCE_V201.pdf` and `.md`. Both are renderings
of a single content definition in `document()`, so they cannot disagree with
each other, and every number is read at build time from
`results/influence_v201/` and `results/direction_sets_v201/`, so neither can
disagree with the data. The earlier hand-maintained Markdown drifted from the
results on sixteen separate claims before a recount caught it.

Figures are registered as `direction_influence` in `scripts/figures.tsv`.

## Adding a new stage

One decorated function in [`scripts/cma/stages.py`](scripts/cma/stages.py). No
driver script, no argparse, no multiprocessing scaffold:

```python
@stage("atp_yield", "ATP produced per unit glucose", needs=("growth",), heavy=True)
def _atp_yield(ctx):
    ...
    return {"_n": len(rows), "per_model": rows}
```

It appears in `--list`, is selectable with `--only atp_yield`, and writes
`<out>/atp_yield.json` automatically.

To make a model set, direction map, heuristic, thermo source or figure
permanently available by name instead of by path, register it once — see
[`scripts/cma/README.md`](scripts/cma/README.md) and `scripts/cma_new.py`.
