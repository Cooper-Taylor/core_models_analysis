# Per-reaction direction sensitivity + heuristic calls (Panel Models chart)

**Date:** 2026-07-06
**Compute:** [`scripts/build_reaction_direction_effects.py`](../scripts/build_reaction_direction_effects.py)
**Data:** `site/data/reaction_direction_effects_panel.json` (100 panel models, 2.1 MB)
**UI:** new **“Reaction-direction heuristics — growth under <, >, =, ?”** section inside the
site’s **Panel Models** tab (select a model → the chart renders in its detail pane).

---

## What it answers

Starting from each panel model’s **default** bounds, every reaction is changed **one
at a time** to each of the four ModelSEED direction options, and the model’s biomass
growth (max flux) is re-solved — so each number is the *isolated* growth effect of that
single reaction’s direction. For the same reaction we also show **where the four
heuristics send it**.

Four options (per the user’s choice, `?` = off/knockout):

| option | bounds | meaning |
|---|---|---|
| `<` | (−1000, 0) | reverse-only |
| `>` | (0, 1000) | forward-only |
| `=` | (−1000, 1000) | reversible |
| `?` | (0, 0) | **unknown → blocked / knocked out** |

Four heuristics (“where the reaction is sent”):

| letter | heuristic | source |
|---|---|---|
| **D** | Default | the reaction’s own in-place bounds |
| **J** | Jankowski (group contribution, Henry-2007 rule) | `Jankowski_2008` |
| **F** | Flamholz 2012 (eQuilibrator) | `Flamholz_2012` |
| **O** | Claude Opus 4.8 | `LLM_Opus_4.8` |

(J/F/O from `results/reaction_directions_literature_vs_llm.tsv`.)

## The chart

For the selected model, a scrollable grid — one row per reaction — showing:

- **heuristic calls**: four badges `D:· J:· F:· O:·` giving the direction each heuristic
  assigns (`NA` = no call, `?` = uncertain);
- **four growth cells** (`<`, `>`, `=`, `?`), heat-coloured by growth relative to the
  model’s baseline (≈baseline green · boosted blue · reduced yellow · strongly-down
  orange · ~dead red · × infeasible grey);
- inside each cell, the **letters of the heuristics that send the reaction to that
  option** — so you can read straight across: “Opus sends this reaction to `>`, which
  gives growth X; Jankowski/Flamholz send it to `<`, growth Y.”

Controls: filter box, **only growth-sensitive** (default on — reactions whose growth
actually varies across options), **only heuristic disagreements**. Rows are ranked by
direction sensitivity (spread of growth across the four options). A collapsible
**Cross-panel patterns** block lists the reactions where heuristics disagree in the
most models and the reactions most often essential-when-off.

### Example (E. coli-like panel model)
- `rxn08173`: heuristics **disagree** — J/F → `<`, O → `>`, D → `=`; growth `<`=66,
  `>`=87.2, `=`=87.2, off=66 → Opus’s call gives higher growth than the thermo methods.
- `rxn00777`: `>`=0 and off=0 but `<`/`=`=87.2 → forcing it forward (or off) kills
  growth; all heuristics agree `=`. (Essential-when-off in 66/100 panel models.)

## How to (re)generate

```bash
PY=/mnt/homes/ctaylor/conda/miniforge3/envs/core_models_analysis/bin/python
$PY scripts/build_reaction_direction_effects.py --workers 48   # 100 models in ~3 s
```
Pure cobra FBA (no modelseedpy needed) — for each (model, reaction) it toggles bounds
inside a `with model:` context and re-solves biomass; ~72k solves total, 0 errors.

## How to view

```bash
cd site && $PY serve.py --static     # http://localhost:8080  → Panel Models → pick a model
```
Works in static mode (precomputed JSON; no `--live`/KBase/FBA needed).

## Files changed

- **new** `scripts/build_reaction_direction_effects.py`, `site/data/reaction_direction_effects_panel.json`.
- `site/static/app.js` — additive: a `STATE.panelRxnDirEffects` loader, a section in the
  Panel Models detail template, and `renderPmRxnDirEffects` / `renderRdeTable` /
  `renderRdeGlobal` (mirrors the existing `renderPm*` precomputed-chart pattern; no
  existing behaviour changed; `node --check` clean).
- `site/static/style.css` — appended `.rde-*` styles.

## Caveats

- Growth = uncapped max biomass flux under the model’s exchange bounds (values can
  exceed nominal rates for very open models); use it **comparatively** across options,
  not as an absolute rate.
- Only reactions with a base id present in the model are probed; a heuristic with `NA`
  for a reaction simply shows no letter in any option cell.
- `=` and `?` differ here only because `?` is defined as off (0,0); if a heuristic emits
  `?` (Opus uncertainty), that call maps to the off column.
