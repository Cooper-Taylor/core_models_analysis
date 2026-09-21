# Contributing

## Setup

```bash
git clone https://github.com/Cooper-Taylor/core_models_analysis
cd core_models_analysis
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[all]"
beginPipeline --doctor          # says what external data is still missing
python3 scripts/fetch_data.py --all
```

## Adding things

Thermo sources, heuristic variants, direction maps, model sets, figures and
probes are **registry entries**, not new scripts. One entry in
`scripts/cma/entries/`, and everything downstream picks it up.

```bash
cma-new thermo-source --key eq3 --display "eQuilibrator 3.0"
cma-new model-set --key my_models
cma-check                       # registries, cross-references, parity
```

A new pipeline **output** is one decorated function in `scripts/cma/stages.py`.
It appears in `beginPipeline --list` and is selectable with `--only` with no
other change. See `PIPELINE.md`.

## Before you push

```bash
pytest -q                       # no external data needed
ruff check scripts/cma tests
cma-check                       # registry parity against the frozen literals
python3 scripts/check_goldens.py --all   # needs ModelSEEDDatabase; ~2 min
```

`check_goldens.py` is the one that matters for anything touching the cascade.
It regenerates each variant and compares the **direction** column separately
from the status text, because status strings embed formatted floats and would
otherwise produce diffs that are not behaviour changes. Exit 2 means the pinned
MSDB reference does not match the one on disk, so it refuses to judge rather
than failing for data reasons.

## Conventions

- **Never modify** `ModelSEEDDatabase/` or `core_models_kegg2/`. Releases are read
  with `git show <tag>:<path>`; the working tree is left alone.
- **Every figure gets a row** in `scripts/figures.tsv` so it can be regenerated
  with `python3 scripts/regen_figures.py <name>`.
- **Paths come from `cma.paths`**, which are functions, not module constants.
  A constant baked into a default argument freezes at import and makes the
  environment variable inert.
- `cma` must stay cheap to import: no `os.environ` reads at import time, and no
  cobra, pandas, numpy or matplotlib at module scope. CI checks this.
