# Cumulative direction-change growth trajectory (Panel Models chart)

**Date:** 2026-07-06
**Compute:** [`scripts/build_cumulative_direction_growth.py`](../scripts/build_cumulative_direction_growth.py)
**Data:** `site/data/cumulative_direction_growth_panel.json` (100 panel models, ~1.0 MB)
**UI:** **“Cumulative direction changes — growth trajectory”** section in the site’s
**Panel Models** tab (select a model → the line chart renders in its detail pane).

---

## What it answers

For each panel model and each heuristic (**Jankowski** group contribution,
**Flamholz 2012** eQuilibrator, **Claude Opus 4.8**):

1. **Rank.** Starting from the model’s **default** bounds, take every reaction whose
   heuristic direction differs from the model’s own direction, and measure each one’s
   *individual* effect on biomass growth (that single change vs the default baseline).
   Sort the reactions by that individual effect, largest first.
2. **Accumulate.** Apply the changes **one on another** in that sorted order, recording
   growth after each prefix: `[baseline, r1, r1+r2, r1+r2+r3, …]`.
3. **Graph.** Plot growth vs the number of reactions applied — one line per heuristic,
   with the model’s default growth as a dashed reference.

Both steps are driven by the KBUtilLib feature we’ve been using,
**`MSTemplateUtils.diff_template_evaluation`**:
- the ranking uses **`mode="independent"`** (each change vs the shared default baseline);
- the accumulation uses **`mode="cumulative"`** (each change vs the previous state).
The offline runner is subclassed so `_evaluate_model_quality` reports biomass flux (a
single `slim_optimize`), and the growth of the baseline + every step is read out of the
sequence of evaluations `diff_template_evaluation` performs. (No FVA battery here, so it
is fast and GLPK-stable; the `_apply_perturbation` bounds edits come from the library.)

## The chart

For the selected model: an SVG line chart, x = # reactions applied, y = biomass growth,
one colored line per heuristic (toggle via the legend), plus a dashed **default**
baseline. Hovering a point shows the reaction just added and its individual Δ. Below the
chart, a collapsible per-heuristic table lists the reactions **in applied order** with
their individual growth Δ.

Reading it: a line that climbs then dips means the highest-individual-effect changes come
first and later (small- or negative-effect) changes pull growth back down — and because
of network interactions the cumulative curve is **not** the sum of the individual effects.

### Example (`GCF_000164985.3`, baseline growth 67.75)
- **Jankowski**: 57 changes; climbs 67.8 → 85.9 → 90.4 → … to a **peak 147.7 at step 51**,
  settling at **135.3** — later changes reduce growth from the peak.
- **Flamholz 2012**: 42 changes; peak 129.1 at step 27.
- **Claude Opus 4.8**: 45 changes but nearly flat (peak 76.0) — Opus’s calls are close to
  the model’s own directions, so they barely move growth.

### Panel-wide (mean peak lift over baseline)
Jankowski **+54.0**, Flamholz **+33.9**, Opus **+11.1** — the thermodynamic methods
relax many reactions and lift growth substantially; Opus is conservative.

## How to (re)generate

```bash
PY=/mnt/homes/ctaylor/conda/miniforge3/envs/core_models_analysis/bin/python
$PY scripts/build_cumulative_direction_growth.py --workers 48   # 100 models in ~22 s
```
Parallel across models; 0 errors on the full panel. J/F/O direction calls come from
`results/reaction_directions_literature_vs_llm.tsv`; only concrete calls (`<`/`>`/`=`)
that differ from the model’s own direction are applied (`NA`/`?` skipped).

## How to view

```bash
cd site && $PY serve.py --static     # http://localhost:8080 → Panel Models → pick a model
```
Static (precomputed JSON); no `--live`/KBase/FBA needed.

## Files

- **new** `scripts/build_cumulative_direction_growth.py`, `site/data/cumulative_direction_growth_panel.json`.
- `site/static/app.js` — additive: `STATE.panelCumDirGrowth` loader, a section in the
  Panel Models detail template, and `renderPmCumDirGrowth` / `drawCdg` (SVG line chart).
  `node --check` clean; headless-render tested (3 polylines + legend + ordered lists).
- `site/static/style.css` — appended `.cdg-*` styles.

## Caveats

- Growth = uncapped max biomass flux; compare across steps/heuristics, not as an
  absolute rate. Ranking is by **signed** individual Δ (largest first), so growth-reducing
  changes are applied last.
- Interactions: the cumulative curve reflects real re-solved growth at each prefix, which
  differs from summing individual effects (hence peaks/dips).
