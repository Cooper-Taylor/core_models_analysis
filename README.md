# core_models_analysis

End-to-end analysis of the 5,683 [ModelSEED](https://modelseed.org/)
core metabolic models in `core_models_kegg2`. The bulk of the work
**compares biological models** — first by who-grows-and-who-doesn't,
then by how reaction directionality (forward / reverse / reversible)
propagates from the ModelSEED database (MSDB) into FBA growth outcomes
across multiple thermodynamic sources.

The repository is organized as a runnable, notebook-driven pipeline:
each notebook reads the artifacts produced by an earlier stage, runs
its own analysis, caches heavy intermediates, and embeds the matching
markdown report.

---

## Quickstart

```bash
git clone https://github.com/Cooper-Taylor/core_models_analysis
cd core_models_analysis
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[all]"

beginPipeline --doctor                    # what is installed, what data is missing
python3 scripts/fetch_data.py --all       # clone MSDB, fetch the models
beginPipeline --list                      # every stage and everything registered
beginPipeline --models core_kegg2 --limit 20 --directions v201_gold
```

`--doctor` is the one to run first: it reports every data root it found or
could not find, every Python package, and exactly which variable to set for
anything missing.

Requires Python 3.11+ and `git`. The heavy inputs (ModelSEEDDatabase, the model
JSONs) are not vendored; `fetch_data.py` gets them and puts them where the
package already looks, so there is usually nothing to configure. If your copies
live elsewhere, set `MSDB_ROOT` and friends, or copy `cma.toml.example` to
`cma.toml`.

### The one command

```bash
beginPipeline --models <path|name> [--directions ...] [--heuristics ...] [--only ...]
```

Every input takes a filesystem path or the name of something registered, so you
can point it at models or direction tables nothing knows about yet without
editing any Python. By default it produces everything it can; `--only` narrows.
Full guide in [PIPELINE.md](PIPELINE.md).

---

## What's in here

| Stage | Notebook | What it does | Embedded report |
|---|---|---|---|
| 00 | `00_Index.ipynb` | Project index + quick access to the descriptive test set | — |
| 01 | `01_GrowthFBA_Pipeline.ipynb` | FBA biomass solve over all 5,683 models on the ModelSEED complete media | `SUMMARY.md` |
| 02 | `02_CharacteristicsAnalysis.ipynb` | Grower vs non-grower model-size + flux distributions | `CHARACTERISTICS.md` |
| 03 | `03_GapAnalysis.ipynb` | Per-non-grower precursor reachability + annotated reading | `GAP_ANALYSIS.md`, `INTERPRETATION.md` |
| 04 | `04_ReactionPrevalence.ipynb` | Reactions enriched in growers vs non-growers | `REACTION_PREVALENCE.md` |
| 05 | `05_DiversePanelSelection.ipynb` | 100-model panel that spans the 3,461 growers | `DIVERSE_SELECTION.md` |
| 06 | `06_ReactionReversibilityHeuristics.ipynb` | Parameterizable port of MSDB's `Estimate_Reaction_Reversibility.py` exercised against every `Reaction_Reversibility_Heuristics_Review.md` suggestion on the 100-model panel | — |
| 07 | `07_NCBITaxonomy.ipynb` | NCBI taxonomy lookup for the 3,461 growers | — |
| 08 | `08_TaxonomyAwareSelection.ipynb` | Taxonomy-aware diverse panel (alternative to 05) | `TAXONOMY_AWARE_SELECTION.md` |
| 09 | `09_ReactionDirectionPipeline.ipynb` | Reaction-direction-driven growth pipeline: audits whether on-disk bounds track MSDB; diffs MSDB branches; reruns FBA under arbitrary direction sources | `REACTION_DIRECTION_PIPELINE.md` |
| 10 | `10_ThermoSourceComparison.ipynb` | Cross-source comparison on the 100-model panel: KBase baseline vs MSDB group-contribution / eQuilibrator / dGPredictor | — |

Notebooks 06–10 are the comparison-of-biological-models core of the
project; 01–05 produce the upstream artifacts they consume.

## Layout

```
core_models_analysis/
├── README.md                          this file
├── requirements.txt                   Python deps
├── notebooks/                         interactive walkthroughs (KBUtils-backed)
│   ├── 00_Index.ipynb … 10_ThermoSourceComparison.ipynb
│   └── README.md                      per-notebook detail
├── scripts/                           regenerable pipeline + notebook builders
│   ├── analyze_growth.py              FBA over all 5,683 models
│   ├── summarize.py                   SUMMARY.md
│   ├── deeper_analysis.py             CHARACTERISTICS, GAP_ANALYSIS, REACTION_PREVALENCE
│   ├── annotate.py                    INTERPRETATION
│   ├── select_diverse.py              100-model panel
│   ├── select_diverse_tax.py          taxonomy-aware panel
│   ├── reversibility_lib.py           parameterizable port of MSDB Estimate_Reaction_Reversibility.py
│   ├── growth_heuristics.py           panel-rebound + FBA driver
│   ├── direction_pipeline.py          reaction-direction pipeline helpers
│   ├── run_thermo_source_variants.py  per-source variant runner
│   ├── thermo_source_figures.py       per-source comparison figures
│   ├── build_*.py                     notebook builders (source of truth — edit these, not the .ipynb)
│   └── …
├── reports/                           markdown writeups, rendered inside each notebook
│   ├── SUMMARY.md, CHARACTERISTICS.md, GAP_ANALYSIS.md, REACTION_PREVALENCE.md,
│   ├── INTERPRETATION.md, DIVERSE_SELECTION.md, TAXONOMY_AWARE_SELECTION.md,
│   ├── REACTION_DIRECTION_PIPELINE.md
│   └── figures/                       PNGs embedded by the notebooks
└── results/                           data artifacts (CSV / JSON)
    ├── results.csv, growers.csv, non_growers.csv, gap_per_model.json
    ├── selected_ids*.txt, selected_models*.{csv,json}
    ├── ncbi_taxonomy.{csv,json}
    ├── rxn_directions_*.{csv,json}    per-source direction maps
    ├── rev_map_{dev,claude}.json      MSDB branch snapshots
    └── thermo_sources/                per-source coverage + FBA tables
```

`data/`, `logs/`, and `notebooks/.kbcache/` are excluded from the repo
(see [`.gitignore`](.gitignore)) — see the next section for what you
need to supply locally.

---

## External dependencies

Four things this repo deliberately does **not** vendor:

### 1. ModelSEEDDatabase (required)

Many scripts and notebooks read MSDB reaction shards, the `KBaseMedia.cpd`
complete media, and per-source thermodynamics dicts. The notebooks treat
the MSDB working tree as **read-only** and snapshot specific branches via
`git show <branch>:Biochemistry/reaction_NN.json`.

```bash
git clone https://github.com/ModelSEED/ModelSEEDDatabase /scratch/ctaylor/ModelSEEDDatabase
# Fetch the dev branch used by notebooks 09 / 10 (per-source thermo dicts).
git -C /scratch/ctaylor/ModelSEEDDatabase fetch origin dev:dev
```

The default path baked into the scripts is
`/scratch/ctaylor/ModelSEEDDatabase`. If you put it elsewhere, see
[Adapting paths](#adapting-paths-to-your-environment) below.

### 2. KBUtils_Local (`kbutillib`) (required for notebooks 06–10)

The `kbutillib.notebook.NotebookSession` helper provides the
content-addressed cache under `notebooks/.kbcache/` that lets the heavy
intermediates (the 56K-reaction MSDB load, per-variant FBA results,
panel descriptors) survive across kernel restarts.

```bash
git clone <your-KBUtils_Local-url> /scratch/ctaylor/KBUtils_Local
pip install -e /scratch/ctaylor/KBUtils_Local
```

Notebooks 00–05 will run without it (they only persist data via CSV /
JSON under `results/`), but the comparison-heavy notebooks (06+) assume
it is importable.

### 3. core_models_kegg2 (5,683 model JSONs)

The input dataset itself. The repo expects the symlink
`data/core_models_kegg2` to point at the unpacked
`core_models_kegg2/` directory of 5,683 KBase JSON models (and
optionally `data/core_models.tar` for the source tarball). These are
not redistributed here.

```bash
mkdir -p data
ln -s /path/to/core_models_kegg2 data/core_models_kegg2
ln -s /path/to/core_models.tar   data/core_models.tar   # optional
```

### 4. ModelSEED v2 genome-scale models (optional)

The 5,420 gap-filled genome-scale models from the ModelSEED v2 manuscript
(Faria et al. 2023), in two media conditions, downloaded from the public KBase
workspaces 155807 and 155808. They need no KBase auth token -- both workspaces
are world-readable.

```bash
python3 /scratch/ctaylor/modelseed2_gs_models/scripts/download_ms2_models.py
python3 /scratch/ctaylor/modelseed2_gs_models/scripts/convert_to_cobra.py --jobs 32
```

They are registered as the model sets `ms2_gsm` (glucose minimal media) and
`ms2_gsm_auxo` (auxotrophy media) in
[`scripts/cma/entries/model_sets.py`](scripts/cma/entries/model_sets.py); point
`MS2_GS_MODELS_DIR` elsewhere if you unpack them somewhere other than
`/scratch/ctaylor/modelseed2_gs_models`. Validate with:

```bash
python3 scripts/cma_check.py --kind modelset --key ms2_gsm
```

They matter because they touch **3,769** distinct ModelSEED reactions against
the core panel's **239** -- a fifteenfold expansion of the reaction surface a
direction map or thermodynamic source can move. See
`/scratch/ctaylor/modelseed2_gs_models/README.md`.

---

## Adding a source, a variant, a model set: the `cma` registry

Thermo sources, heuristic variants, direction maps, model sets, figures and
probes are declared once in [`scripts/cma/entries/`](scripts/cma/entries/) and
read from there by everything downstream, instead of being re-declared as
literals across the pipeline.

```bash
python3 scripts/cma_new.py thermo-source --key eq3 --display "eQuilibrator 3.0"
python3 scripts/cma_new.py variant --key llm_gpt5 --kind overlay
python3 scripts/cma_new.py model-set --key ecoli_gsm
python3 scripts/cma_check.py          # validate registries + parity with the old literals
python3 scripts/check_goldens.py --all   # prove the numbers did not move
```

See [`scripts/cma/README.md`](scripts/cma/README.md).
---

## Setup

```bash
pip install -e ".[all]"        # or ".[dev]" for just the tests
```

Extras: `figures` (matplotlib, plotly, kaleido), `notebooks` (jupyter),
`dev` (pytest, ruff), `all`.

Installing puts three commands on PATH -- `beginPipeline`, `cma-check` and
`cma-new`. They also work straight from a checkout with nothing installed, via
the `beginPipeline` shim at the repo root.

### Where the data goes

`cma.paths` resolves each external root in this order, first hit wins:

1. an environment variable
2. a key in `cma.toml` (beside the repo, or `~/.config/cma/config.toml`)
3. conventional locations next to the repo
4. otherwise an error naming the variable to set

| root | env var | what it is |
|---|---|---|
| MSDB | `MSDB_ROOT` | the ModelSEEDDatabase clone; required by almost everything |
| MSDB snapshot | `MSDB_SNAPSHOT_ROOT` | a pinned snapshot the grading scripts read |
| core models | `CORE_MODELS_DIR` | the 5,683 core model JSONs |
| genome-scale models | `MS2_GS_MODELS_DIR` | the ModelSEED v2 models |

The two MSDB roots are deliberately separate variables. They have held
different biochemistry, so a single shared name silently repointed part of the
pipeline at the other one's data.

`beginPipeline --doctor` prints what each resolved to.

## Running the pipeline

### From scratch (regenerate everything)

```bash
# Stage 1 — growth + characteristics + gaps + prevalence + diverse panel
python3 scripts/analyze_growth.py
python3 scripts/summarize.py
python3 scripts/deeper_analysis.py
python3 scripts/annotate.py
python3 scripts/select_diverse.py

# Stage 2 — taxonomy-aware panel
python3 scripts/select_diverse_tax.py

# Stage 3 — per-source thermodynamics variants (notebook 10 backing data)
python3 scripts/run_thermo_source_variants.py
python3 scripts/build_thermo_source_figures.py

# Rebuild every notebook from its builder
python3 scripts/build_notebooks.py
python3 scripts/build_reversibility_notebook.py
python3 scripts/build_taxonomy_aware_notebook.py
python3 scripts/build_direction_pipeline_notebook.py
python3 scripts/build_thermo_source_comparison_notebook.py

# Re-execute the notebooks (populates .kbcache/)
cd notebooks && jupyter execute --inplace *.ipynb
```

### Presentation figures (thermodynamics → ATP/growth impact)

Six interactive figures that summarize how the reversibility heuristics reshape
core-model growth/ATP flux. Reads existing artifacts only (variant diff JSONs, the
N=50 statistical panel, the live cascade) — no FBA is recomputed.

```bash
python3 scripts/build_presentation_figures.py            # all 6, HTML + PNG + site export
python3 scripts/build_presentation_figures.py --no-png   # HTML only (skip kaleido/Chrome)
python3 scripts/build_presentation_figures.py --figures 1,3,6
```

Outputs: `reports/presentation/index.html` (scrollable dashboard), `figN_*.html`
(interactive) and `png/figN.png` (static, for slides); plus `site/data/figures/*.html`
+ `manifest.json` and a `site/figures.html` viewer for embedding in the static site.
PNG export needs Chrome — install once with `plotly_get_chrome` (or pass `--no-png`).

### Re-run a single notebook interactively

```bash
cd notebooks
jupyter lab 10_ThermoSourceComparison.ipynb
```

First execution of notebook 06 / 10 takes ~45 s (loads the 56K-reaction
MSDB, runs every variant); subsequent runs hit `.kbcache/` and finish
under 5 s.

### Inspect cached intermediates

```bash
sqlite3 notebooks/.kbcache/catalog.sqlite 'select id, type, n_bytes from cache_objects;'
```

```python
from kbutillib.notebook import NotebookSession
session = NotebookSession.for_notebook(project_name='core_models_analysis')
session.cache.load('descriptive_test_panel')        # → {'ids': [...], 'coverage': {...}}
session.cache.load('msdb_reactions_v1')             # → full MSDB reactions dict (~56K)
```

---

## The 100-model descriptive test set

`results/selected_ids.txt` is the 100-model panel used by every
comparison notebook from 06 onward:

```python
ids = open('results/selected_ids.txt').read().split()
# 100 model IDs spanning the 3,461 growers
```

Methodology, per-model selection reason, and Jaccard-coverage validation
live in `reports/DIVERSE_SELECTION.md` (genome-similarity-based) and
`reports/TAXONOMY_AWARE_SELECTION.md` (taxonomy-aware alternative).

---

## Where to start reading

- **Just the headlines** → `reports/SUMMARY.md` then
  `reports/REACTION_DIRECTION_PIPELINE.md`.
- **The pipeline end-to-end** → notebooks 01 → 05 → 06 → 09 → 10.
- **The comparison work specifically** → notebooks 06, 09, 10 and the
  corresponding scripts (`reversibility_lib.py`, `direction_pipeline.py`,
  `run_thermo_source_variants.py`, `thermo_source_figures.py`).
