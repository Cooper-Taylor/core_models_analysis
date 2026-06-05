# Reaction-Reversibility Heuristics Explorer

Web UI for browsing every `ReversibilityConfig` variant from notebook
06 against the 100-model descriptive growth panel, drilling into
individual reactions, and sandboxing per-reaction FBA effects.

## Pieces

```
site/
├── README.md            this file
├── serve.py             stdlib HTTP server + live FBA endpoint
├── static/
│   ├── index.html       single-page UI (Variant Browser, Reaction Explorer, Flux Sandbox)
│   ├── app.js           UI logic + fetch wrappers
│   └── style.css        styling
└── data/                JSON snapshots (rebuilt by scripts/build_site_data.py)
    ├── manifest.json
    ├── baseline.json
    ├── panel.json
    ├── panel_rxnsets.json
    ├── reactions_panel.json
    ├── reactions_other.json
    └── variants/{tag}.json
```

## Rebuilding the data

Two scripts:

1. `scripts/export_thermo_variants.py` --- writes
   `core_models_analysis/thermo_variants/{tag}/Estimated_Reaction_Reversibility_Report{,_EQ,_GC}.txt`
   for every variant, mirroring the layout MSDB writes under
   `ModelSEEDDatabase/Scripts/Thermodynamics/`.  Re-run this when
   `variant_catalog.py` or `reversibility_lib.py` changes.

2. `scripts/build_site_data.py` --- reads the per-variant reports + the
   notebook 06 kbcache (variant cascade + FBA results) and writes the
   JSON snapshots above.  Re-run this after `export_thermo_variants.py`
   or when notebook 06 caches change.

```sh
cd /scratch/ctaylor/core_models_analysis
python3 scripts/export_thermo_variants.py        # ~2-3 min
python3 scripts/build_site_data.py               # ~20s
```

## Launching the server

```sh
cd /scratch/ctaylor/core_models_analysis
python3 site/serve.py --port 8769 --preload
# -> open http://127.0.0.1:8769/ in a browser
```

Flags:

- `--host` (default `127.0.0.1`) --- pass `0.0.0.0` for LAN access.
- `--port` (default `8765`) --- the example above uses `8769` because
  some other process was holding `8765` on the workstation; pick any
  free port.
- `--preload` --- load all per-variant maps + reactions index up front
  rather than on the first request (1s warm-up vs slightly faster first
  click).

The server is stdlib-only (`http.server`, `socketserver`,
`multiprocessing`).  No Flask / FastAPI / Node required.

## What the three tabs do

1. **Variant Browser** --- one row per `ReversibilityConfig` variant
   with: reactions changed vs MSDB baseline, panel models that flip
   grower/non-grower, panel models whose biomass flux moved.  Clicking
   a row shows the transition matrix (`> -> =`, `= -> <`, ...) and the
   top changed reactions; clicking a rxn ID jumps to the Reaction
   Explorer.

2. **Reaction Explorer** --- searchable list of every reaction that
   either appears in any panel model or has its direction changed by
   some variant.  Click any reaction to see its full stoichiometry,
   ΔG′°, which variants change its direction, and run a
   **per-mode panel sweep** that does live FBA under
   `off / forward / reverse / free / as_is` for that one reaction across
   panel models that contain it.  Useful for "would this one heuristic
   flip actually matter for biology?"

3. **Flux Sandbox** --- pick any variant, add an arbitrary set of
   per-reaction overrides (`forward`, `reverse`, `free`, `off`), and
   run FBA across the full panel or a chosen subset.  Backed by the
   `/api/panel_fba` endpoint, which uses
   `growth_heuristics.run_panel` with a small monkey-patch that adds
   the `off` mode (`lb=ub=0`).

## API

| route                       | method | body / params                                              | returns                                       |
|----------------------------|--------|------------------------------------------------------------|-----------------------------------------------|
| `/api/health`              | GET    | --                                                         | `{ok, n_variants}`                            |
| `/api/rxn/<rxn_id>`        | GET    | --                                                         | per-rxn metadata + `changed_by` list          |
| `/api/panel_fba`           | POST   | `{variant, overrides, models?, n_workers?}`                | `{n_models, elapsed_s, results: [...]}`       |
| `/api/reaction_impact`     | POST   | `{rxn_id, variant?, modes?, models?}`                      | `{baseline, by_mode: {mode: {model: {...}}}}` |
| `/static/*`, `/data/*`     | GET    | --                                                         | files                                         |

## Constraints honored

- Nothing under `ModelSEEDDatabase/` or `core_models_kegg2/` is
  modified.  All variant reports live under
  `core_models_analysis/thermo_variants/{tag}/`; the live FBA endpoint
  rebinds cobra models *in memory* only.
- The site shares notebook 06's kbcache at
  `core_models_analysis/notebooks/.kbcache/` --- a cold first build
  (before notebook 06 was ever run) will load the 56K MSDB reactions
  through BiochemPy, which adds ~45 s.
