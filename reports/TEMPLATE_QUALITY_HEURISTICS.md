# Template Quality under reaction-direction heuristics (website feature)

**Date:** 2026-07-06
**Script:** [`scripts/template_quality_heuristics.py`](../scripts/template_quality_heuristics.py)
**Data:** `site/data/template_quality_all.json` (5,683 models) + `…_all.jsonl` (raw, resumable)
**UI:** `site/static/template_quality.html` → served at **`/static/template_quality.html`**
(linked from the main explorer nav as **“Template Quality ↗”**)
**Depends on:** [`scripts/direction_change_template_eval.py`](../scripts/direction_change_template_eval.py)
(the `OfflineTemplateEval` runner + shims) and
[`reports/TEMPLATE_DIRECTION_EVAL.md`](TEMPLATE_DIRECTION_EVAL.md);
function-level defects in the KBUtilLib suite are in
`KBUtils_Local/agent-io/audits/2026-07-06-ms-template-utils-review.md`.

---

## What it does

For **every core model (all 5,683)** it runs the KBUtilLib `MSTemplateUtils`
evaluation battery — FVA reaction classification (dead / forward-only /
reverse-only / reversible / essential), closed-mode loop detection, and
producible/consumable-metabolite sweeps — under **four reaction-direction
schemes**, then measures how those qualities change and surfaces the comparison
in the website.

| scheme | direction source | provenance |
|---|---|---|
| **Default** | the model's own in-place bounds | `core_models_kegg2/*.json` |
| **Jankowski** | group contribution (Henry-2007 feasibility rule) | `Jankowski_2008` column |
| **Flamholz 2012** | eQuilibrator reversibility index | `Flamholz_2012` column |
| **Claude Opus 4.8** | LLM-predicted directionality | `LLM_Opus_4.8` column |

Sources (2)–(4) are the columns of
`results/reaction_directions_literature_vs_llm.tsv` (built by
`estimate_directions_literature.py`) — each method uses **its own** directionality
rule, not the MSDB cascade. A direction is applied by rewriting the matching
in-model reaction's bounds (`>`→(0,1000), `<`→(-1000,0), `=`→(-1000,1000));
`NA`/`?` (no call / uncertain) leave the model's bound unchanged. Growth is the
model's actual max biomass flux (an honest metric, not the suite's essential-count
proxy — see audit #4/#10).

> Note on labels: the user's "2007 Jankowski" maps to the group-contribution
> method (column `Jankowski_2008`; group contribution relies on Henry 2007 for
> directionality). "2012 Flamholtz" maps to eQuilibrator (`Flamholz_2012`).

## Headline result (mean per model, across all 5,683)

| scheme | biomass flux | closed-mode loops | dead | reversible | rxns overridden | models growth ↑ / ↓ |
|---|--:|--:|--:|--:|--:|--:|
| Default | 34.7 | 9.3 | 63.0 | 24.7 | — | — |
| **Jankowski (group contribution)** | **77.4** | **49.7** | 56.9 | 60.2 | 42.9 | **4023 / 84** |
| **Flamholz 2012 (eQuilibrator)** | 47.5 | 36.1 | 60.7 | 46.2 | 29.4 | 3278 / 607 |
| **Claude Opus 4.8** | 29.6 | 24.0 | 68.1 | 26.4 | 35.1 | 1295 / 2428 |

**Interpretation.** The two thermodynamic methods are **permissive** — they
relax many reactions to reversible, which raises growth (group contribution lifts
4,023 models) but **inflates thermodynamically-infeasible loops ~4–5×** (9 → 50 / 36
closed-mode reactions per model). **Claude Opus 4.8 is conservative** — its
reversibility mix stays close to the model's own bounds, it introduces the fewest
loops of the three override schemes, and it tends to *reduce* growth (2,428 models
lose growth vs 1,295 that gain). In short: the LLM directions are the most
constraint-preserving; the thermo methods trade loop-inflation for growth. This is
exactly the quality/feasibility trade-off the template suite is meant to expose,
now visible per-model and in aggregate for the whole core-model collection.

## How to view

```bash
cd site
/mnt/homes/ctaylor/conda/miniforge3/envs/core_models_analysis/bin/python serve.py --static
# open http://localhost:8080/  ->  click "Template Quality ↗"  (or go straight to
#      http://localhost:8080/static/template_quality.html )
```
(The page is pure static + JSON — no `--live`/KBase/FBA needed to browse.) It shows
(a) the cross-model averages table with a **panel-only (100)** toggle, and (b) a
per-model detail table (search any model id; Δ vs Default is colour-coded by whether
the change is favourable). Seeded with E. coli `GCF_000005845.2`.

## How to (re)generate

```bash
PY=/mnt/homes/ctaylor/conda/miniforge3/envs/core_models_analysis/bin/python
$PY scripts/template_quality_heuristics.py --panel                # 100-model panel
$PY scripts/template_quality_heuristics.py --all --workers 96     # all 5,683 (~7 min on this box)
$PY scripts/template_quality_heuristics.py --all --resume         # continue an interrupted run
```
Parallel (multiprocessing), resumable (writes `…_all.jsonl` incrementally then
aggregates), and robust (per-model try/except; the full run completed **0 errors**).
~13 models/s across 96 workers.

## Requirements & caveats

- Needs the conda env with cobra + **modelseedpy 0.4.2** + `kbutillib` (the pulled
  `ms_template_utils.py`). The `MSTemplateUtils` suite is driven through
  `direction_change_template_eval.OfflineTemplateEval`, which supplies the three
  documented shims (offline construction, cobra objective in place of the private
  `ObjectivePkg`, per-stage `model.copy()` isolation to avoid the GLPK
  process-abort). See `TEMPLATE_DIRECTION_EVAL.md` §5.
- Media is `None` in the battery, so the `rich` and `minimal` reaction-class passes
  are identical (audit #9/#19); the UI reports the `rich` set. `functional_biolog_media`
  is empty here (`MSGrowthPhenotypes.from_dict` is absent in modelseedpy 0.4.2), so
  it is not shown.
- Growth = max biomass flux under the model's exchange bounds (uncapped, so values
  reach the 1000 flux ceiling for very open models); use it comparatively across
  schemes, not as an absolute rate.
- `NA`/`?` directions are skipped (bound left unchanged), so a scheme only perturbs
  reactions it has a concrete call for (mean 29–43 reactions/model).
```
