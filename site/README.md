# Reaction-Reversibility Heuristics Explorer

Web UI for browsing every `ReversibilityConfig` variant from notebook
06 against the 100-model descriptive growth panel, drilling into
individual reactions, and (optionally) sandboxing per-reaction FBA
effects.

## Quickstart (one command, no install)

```sh
cd /path/to/core_models_analysis
python3 site/serve.py
# open http://127.0.0.1:8765/ in a browser
```

That is the whole quickstart. No `pip install`, no virtualenv, no env
vars. Requires Python 3.9 or newer (tested on 3.10 / 3.11 / 3.12 /
3.13). The `site/` directory ships with every JSON snapshot the static
UI needs.

If port 8765 is already in use, the server tries 8766..8775
automatically and prints the chosen port. To pin a specific port, pass
`--port N` (an explicit port will not auto-fall-back; it fails fast
with a clear error).

## Static vs live mode

`serve.py` runs in one of two modes. The mode is announced on stdout
at startup (`[mode] static (FBA disabled): ...` or `[mode] live (FBA
enabled): ...`) and is also exposed via `/api/health` so the UI can
gate the FBA controls.

| UI surface                                  | Static (default) | Live FBA (`--live`) |
|---------------------------------------------|------------------|---------------------|
| Variant Browser tab                         | works            | works               |
| Reaction Explorer tab (search, detail)      | works            | works               |
| Reaction Explorer: per-mode panel sweep     | hidden + notice  | works               |
| Flux Sandbox tab                            | hidden           | works               |
| `/api/panel_fba`, `/api/reaction_impact`    | 503 with banner  | works               |

### How the mode is chosen

Default is **static** — even on a workstation that has every live-mode
prerequisite installed. To opt into live mode you must pass `--live`
explicitly. This means a fresh clone always works with the literal
quickstart command above.

Auto-detection checks all six prerequisites and reports them; `--live`
errors out if any are missing:

1. `importlib.util.find_spec('cobra')` returns a real spec (i.e. cobra
   is importable);
2. `<repo>/../ModelSEEDDatabase/Libs/Python/` exists (sibling of
   `core_models_analysis/`);
3. `notebooks/.kbcache/` exists (built by notebook 06);
4. `results/selected_ids.txt` exists (the 100-model panel);
5. `site/data/baseline.json` exists (built by `scripts/build_site_data.py`);
6. `site/data/panel_rxnsets.json` exists (same script).

| Flag               | Effect                                                                                      |
|--------------------|---------------------------------------------------------------------------------------------|
| (none)             | Static mode. FBA is disabled regardless of prerequisites.                                   |
| `--static`         | Force static mode. Useful in tests/CI.                                                      |
| `--live`           | Force live mode. Errors out at startup if any of (1)–(6) is missing, listing the offender.  |

### What ships in the repo (and what doesn't)

`site/` is self-contained for static use — you can zip just that
directory, drop it elsewhere, run `python3 site/serve.py`, and the
static UI works. The data files shipped in git:

- `site/data/manifest.json`              (~6 KB)
- `site/data/reactions_panel.json`       (~380 KB)
- `site/data/reactions_other.json`       (~4.2 MB — non-panel reactions, lazy-loaded)
- `site/data/variants/*.json`            (~1.2 MB total — per-variant cascade diffs)

NOT shipped (regenerable by `scripts/build_site_data.py`):

- `site/data/baseline.json`              (FBA baseline cascade map)
- `site/data/panel.json`                 (FBA-only)
- `site/data/panel_rxnsets.json`         (per-model reaction sets)

The `.gitignore` uses an explicit-deny list, so any NEW file written
under `site/data/` by a future build will be committed by default. If
you add an FBA-only output, also add it to the .gitignore block.

## Enabling live FBA

Live mode requires the full upstream pipeline:

1. Python deps: `pip install -r requirements.txt` (cobra, etc.). Or
   activate an existing env with cobra installed.
2. ModelSEEDDatabase cloned as a sibling of `core_models_analysis/`:
   ```sh
   cd <parent-of-core_models_analysis>
   git clone https://github.com/ModelSEED/ModelSEEDDatabase.git
   ```
3. The 100-model panel id list at `results/selected_ids.txt` (already
   committed in this repo).
4. The notebook-06 kbcache at `notebooks/.kbcache/` (run notebook 06
   once to populate, or copy from a teammate).
5. The FBA snapshots at `site/data/baseline.json` and
   `site/data/panel_rxnsets.json` — regenerate with
   ```sh
   cd /path/to/core_models_analysis
   python3 scripts/export_thermo_variants.py   # ~2-3 min
   python3 scripts/build_site_data.py          # ~20s
   ```
6. Launch with `--live`:
   ```sh
   python3 site/serve.py --live
   ```
   The first `/api/panel_fba` call spawns multiprocessing workers and
   takes 5–30 s; subsequent calls within the same run are faster.

## Pieces

```
site/
├── README.md            this file
├── serve.py             stdlib HTTP server + auto static/live mode
├── static/
│   ├── index.html       single-page UI
│   ├── app.js           UI logic + fetch wrappers
│   └── style.css        styling
└── data/                JSON snapshots (see "What ships" above)
    ├── manifest.json
    ├── baseline.json            (FBA-only, gitignored)
    ├── panel.json               (FBA-only, gitignored)
    ├── panel_rxnsets.json       (FBA-only, gitignored)
    ├── reactions_panel.json
    ├── reactions_other.json
    └── variants/{tag}.json
```

## Rebuilding the data

The two-script pipeline (live-mode users only — static-mode users
never need to run these because the static-mode files ship in git):

1. `scripts/export_thermo_variants.py` — writes
   `core_models_analysis/thermo_variants/{tag}/Estimated_Reaction_Reversibility_Report{,_EQ,_GC}.txt`
   for every variant, mirroring the layout MSDB writes under
   `ModelSEEDDatabase/Scripts/Thermodynamics/`.  Re-run when
   `variant_catalog.py` or `reversibility_lib.py` changes.

2. `scripts/build_site_data.py` — reads the per-variant reports + the
   notebook 06 kbcache (variant cascade + FBA results) and writes the
   JSON snapshots above.  Re-run after `export_thermo_variants.py` or
   when notebook 06 caches change.

```sh
cd /path/to/core_models_analysis
python3 scripts/export_thermo_variants.py        # ~2-3 min
python3 scripts/build_site_data.py               # ~20s
```

## What the three tabs do

1. **Variant Browser** — one row per `ReversibilityConfig` variant
   with: reactions changed vs MSDB baseline, panel models that flip
   grower/non-grower, panel models whose biomass flux moved.  Clicking
   a row shows the transition matrix (`> -> =`, `= -> <`, ...) and the
   top changed reactions; clicking a rxn ID jumps to the Reaction
   Explorer.

2. **Reaction Explorer** — searchable list of every reaction that
   either appears in any panel model or has its direction changed by
   some variant.  Click any reaction to see its full stoichiometry,
   ΔG′°, which variants change its direction. In **live mode** there
   is also a **per-mode panel sweep** widget that does live FBA under
   `off / forward / reverse / free / as_is` for that one reaction
   across panel models that contain it.

3. **Flux Sandbox** (live mode only) — pick any variant, add an
   arbitrary set of per-reaction overrides (`forward`, `reverse`,
   `free`, `off`), and run FBA across the full panel or a chosen
   subset.  Backed by `/api/panel_fba`, which uses
   `growth_heuristics.run_panel` with a small monkey-patch that adds
   the `off` mode (`lb=ub=0`).

## API

| route                       | method | body / params                                        | returns                                          |
|----------------------------|--------|-----------------------------------------------------|--------------------------------------------------|
| `/api/health`              | GET    | --                                                  | `{ok, static_mode, n_variants}` (never raises)   |
| `/api/rxn/<rxn_id>`        | GET    | --                                                  | per-rxn metadata + `changed_by` list             |
| `/api/panel_fba`           | POST   | `{variant, overrides, models?, n_workers?}`         | `{n_models, elapsed_s, results: [...]}` (503 in static mode) |
| `/api/reaction_impact`     | POST   | `{rxn_id, variant?, modes?, models?}`               | `{baseline, by_mode: {mode: {model: {...}}}}` (503 in static mode) |
| `/static/*`, `/data/*`     | GET    | --                                                  | files                                            |

In static mode the two POST endpoints return HTTP 503 with body
`{error: 'FBA disabled in static mode — ...', static_mode: true}`.
The frontend uses `/api/health.static_mode` to hide the UI that
depends on them.

## Advanced launch options

```sh
# Pin the port:
python3 site/serve.py --port 8769

# Preload everything on startup (slightly faster first click; same total work):
python3 site/serve.py --preload

# Force live mode (will exit with an error if prerequisites are missing):
python3 site/serve.py --live

# Allow LAN access. NOTE: on macOS this triggers a firewall prompt for
# the `python3` binary on first launch.
python3 site/serve.py --host 0.0.0.0
```

The server is stdlib-only (`http.server`, `socketserver`,
`threading`).  No Flask / FastAPI / Node required.  `multiprocessing`
is used only by live-mode FBA via `growth_heuristics.run_panel`.

## Constraints honored

- Nothing under `ModelSEEDDatabase/` or `core_models_kegg2/` is
  modified.  All variant reports live under
  `core_models_analysis/thermo_variants/{tag}/`; the live FBA endpoint
  rebinds cobra models *in memory* only. Static mode honors this
  trivially — it never touches either tree.
- The live FBA path shares notebook 06's kbcache at
  `core_models_analysis/notebooks/.kbcache/` — a cold first build
  (before notebook 06 was ever run) will load the 56K MSDB reactions
  through BiochemPy, which adds ~45 s.
