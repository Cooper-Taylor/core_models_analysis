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

Three things this repo deliberately does **not** vendor:

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

---

## Setup

```bash
# 1. Clone
git clone <this-repo-url> core_models_analysis
cd core_models_analysis

# 2. Create a conda / venv environment with Python 3.12
conda create -n core_models_analysis python=3.12 -y
conda activate core_models_analysis

# 3. Install dependencies
pip install -r requirements.txt

# 4. Install KBUtils_Local in editable mode (see above)
pip install -e /path/to/KBUtils_Local

# 5. Point data/ at your local copies (see external dependencies)
mkdir -p data
ln -s /path/to/core_models_kegg2     data/core_models_kegg2
ln -s /path/to/ModelSEEDDatabase     /scratch/ctaylor/ModelSEEDDatabase  # or adapt paths
```

### Adapting paths to your environment

Every script reads its two root paths from environment variables and
falls back to the original defaults if they are unset:

| Env var | Default | Meaning |
|---|---|---|
| `CORE_MODELS_ANALYSIS_DIR` | `/scratch/ctaylor/core_models_analysis` | this repo on disk |
| `MSDB_ROOT` | `/scratch/ctaylor/ModelSEEDDatabase` | ModelSEEDDatabase clone |

If your clones live at the defaults, do nothing. Otherwise export the
overrides before running anything:

```bash
export CORE_MODELS_ANALYSIS_DIR=/path/to/core_models_analysis
export MSDB_ROOT=/path/to/ModelSEEDDatabase
```

Add these to your shell profile (or a project-local `.envrc` if you use
[direnv](https://direnv.net/)) so every shell, notebook kernel, and
`jupyter execute` invocation inherits them.

---

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
