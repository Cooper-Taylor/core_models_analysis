# Template-direction evaluation: wiring reaction-direction changes into KBUtilLib's `diff_template_evaluation`

**Date:** 2026-07-06
**Script:** [`scripts/direction_change_template_eval.py`](../scripts/direction_change_template_eval.py)
**Outputs:** `results/template_direction_eval/`
**Companion audit** (defects + suggestions in the KBUtilLib functions):
`KBUtils_Local/agent-io/audits/2026-07-06-ms-template-utils-review.md`

---

## 1. What this connects

KBUtilLib recently gained `MSTemplateUtils` (`kbutillib/ms_template_utils.py`, upstream
commit `d1fbe88`), a "Template Evaluation Suite" with four public functions:

| function | role |
|---|---|
| `build_full_template_model(template)` | build a full cobra model from a ModelSEED template (needs modelseedpy + KBase) |
| `evaluate_template_quality(template, …)` | run the full FBA/FVA battery → one structured report (reaction classes, closed-mode loops, Biolog growth, producible/consumable metabolites) |
| `render_template_report(report)` | report dict → markdown |
| **`diff_template_evaluation(model, perturbations, mode=…)`** | apply model-level edits (add/remove/**modify bounds**) and report, per edit, what changed across every report category |

The **reaction-direction work** in this repo produces per-source direction maps in
`results/rxn_directions_*.json` (`{rxn_id: ">"|"<"|"="|"?"}`) from the ModelSEED
thermodynamics cascade and from alternate sources (group-contribution, eQuilibrator,
dGPredictor). **A reaction-direction change is exactly a bounds change**, so it maps
directly onto a `diff_template_evaluation` `modify` perturbation. That is what the
script does: pick a *baseline* and a *new* direction source, turn every reaction whose
direction differs (and that exists in the model) into a `modify` perturbation, and let
`diff_template_evaluation` attribute the functional consequences.

```
>  forward-only  ->  (lower=0,     upper=+1000)
<  reverse-only  ->  (lower=-1000, upper=0)
=  reversible    ->  (lower=-1000, upper=+1000)
?  unknown       ->  SKIPPED   (no thermodynamic call is not a bound)
```

This reuses the repo's existing conventions: `_bounds_for_rev` magnitude (1000) from
`growth_heuristics.py`, and `seed_annotation.normalize_seed_id` to map a map key
`rxn00549` to the model reaction id `rxn00549_c0` (handling the stray-`_c` annotation bug).

## 2. How to run

Use the conda env that has cobra + modelseedpy + kbutillib:
`/mnt/homes/ctaylor/conda/miniforge3/envs/core_models_analysis/bin/python`

```bash
PY=/mnt/homes/ctaylor/conda/miniforge3/envs/core_models_analysis/bin/python

# 1. Offline logic check (no modelseedpy needed):
$PY scripts/direction_change_template_eval.py --self-test

# 2. Build the perturbation set only, no FBA (no modelseedpy needed):
$PY scripts/direction_change_template_eval.py --dry-run \
      --model-id GCF_000005845.2 --baseline cascade_live --new group-contribution

# 3. Full live diff (default = batch: one combined edit, ~2 evals):
$PY scripts/direction_change_template_eval.py \
      --model-id GCF_000005845.2 --baseline cascade_live --new group-contribution --mode batch

# 4. Per-edit attribution (independent mode, ~5 s/edit — cap with --limit):
$PY scripts/direction_change_template_eval.py --mode independent --limit 8
```

**Key flags:** `--baseline {cascade_live|group-contribution|equilibrator|dgpredictor|msdb_dev|msdb_claude|model|<path>}`
(`model` = derive the baseline direction from the model's *own* bounds, which avoids
emitting no-op modifies), `--new <source>`, `--mode {independent|cumulative|batch}`,
`--limit N`, `--solver`, `--verbose`.

**Outputs** (in `results/template_direction_eval/`): `<tag>.diff.json` (full diff),
`<tag>.summary.csv` (one row per perturbation: provenance + changed-category counts +
growth), `<tag>.summary.md`, and `<tag>.growth.csv` (per-reaction honest biomass-flux
growth sensitivity). Dry-run writes `<tag>.perturbations.json`.

## 3. Demonstrated result (E. coli core model `GCF_000005845.2`)

**Pairing:** baseline = `cascade_live`, new = `group-contribution`. Of the 4,720
map-wide direction changes, **31** touch reactions present in this 213-reaction model
(16 `=`→`>`, 13 `>`→`=`, 2 `<`→`=`).

**`batch` mode — applying all 31 flips at once** (`…__batch.diff.json`) reports, vs the
baseline:

- **`closed_mode_reactions`: +54** — switching to group-contribution directions opens
  54 new reactions that can carry flux with all exchanges closed, i.e. new
  thermodynamically-infeasible loops. This is the most striking signal.
- **reaction classes:** +13 reversible, −11 forward_only, −2 dead; essential set shifts
  (bio1 +7/−4, bio2 +6/−5).
- **metabolites:** `consumable.complete` −1; `producible` +1/−1.
- **shipped `growth_change` proxy:** reports `interpretation: "gained"` (bio1 essential
  14→17). **This is misleading** (see §5): the *actual* max biomass flux is unchanged
  (87.21 → 87.21) for every single flip — see `…__batch.growth.csv`.

**`independent` mode** attributes these to individual flips: e.g. `rxn00122` (`>`→`=`)
alone changes 6 categories (+30/−14 members); most single flips change nothing, and
none change max biomass flux.

Interpretation: for a well-connected core model, re-sourcing directions from
group-contribution mainly **relaxes constraints and creates loops** (closed-mode +54)
rather than changing what the model can grow on — a concrete, quantitative version of
"what did this direction change actually do?"

## 4. Environment changes made to get a live run

The user opted to install dependencies and run live. Changes made:

**`KBUtils_Local`** (git working tree; uncommitted `ai_curation_utils.py` WIP left untouched):
- `git checkout origin/main --` pulled **`src/kbutillib/ms_fba_utils.py`** (Phase‑1
  battery), **`src/kbutillib/ms_template_utils.py`** (new), and
  **`src/kbutillib/data/biolog_phenotypes.json`** (Biolog stash). These are needed for
  `from kbutillib.ms_template_utils import MSTemplateUtils` to resolve (the package is
  an editable install pointing at `KBUtils_Local/src`).

**conda env `core_models_analysis`:**
- `pip install modelseedpy` → **modelseedpy 0.4.2** (+ deps: scipy 1.18, scikit-learn,
  networkx, chemw, chemicals, fluids, pubchempy, sigfig, joblib, narwhals — **no
  downgrade** of cobra 0.31.1 / numpy 2.4 / pandas 2.3).
- `pip install highspy` (1.15.1) — installed while diagnosing solver crashes; **not
  ultimately required** (GLPK works once evaluation stages are isolated), but harmless.

**`core_models_analysis`:** added `scripts/direction_change_template_eval.py`, this
report, and `results/template_direction_eval/` outputs.

## 5. Compatibility shims (why the script has an `OfflineTemplateEval` subclass)

`MSTemplateUtils` was written against an internal ModelSEEDpy build; three
narrowly-scoped, documented shims (all in the script) let it run on the public stack
without changing what `diff_template_evaluation` computes:

- **S1 — offline construction.** `KBModelUtils.__init__` demands a KBase token and
  builds a KBaseAPI client; the diff path needs only `self.MSModelUtil`. Built via
  `__new__` + wiring the modelseedpy classes.
- **S2 — objective shim.** `set_objective_from_string` calls
  `pkgmgr.getpkg("ObjectivePkg")`, which exists only in the private modelseedpy fork
  (absent from PyPI 0.4.2 *and* ModelSEEDpy GitHub `main`). Overridden to set the cobra
  objective directly. (`ObjConstPkg`, used for the fraction-of-optimum pass, is public
  and left as-is.)
- **S3 — stage isolation.** `run_fva`'s growth-forced pass adds an `ObjConstPkg`
  "biomass ≥ fraction·optimum" constraint that is **never removed**. The shipped
  `_evaluate_model_quality` runs all stages on one model, so
  `find_closed_mode_reactions` (which zeroes every exchange) inherits that constraint,
  becomes infeasible, and **GLPK aborts the process (SIGABRT)**. The override runs each
  stage on a fresh `model.copy()`. This is a genuine upstream bug (audit #5/#5b).

Also: `simulate_biolog` no-ops here (`MSGrowthPhenotypes.from_dict` is absent in
modelseedpy 0.4.2), so `functional_biolog_media` comes back empty — handled, not fatal.

Because the shipped `growth_change` is an essential-reaction *count*, not biomass flux
(audit #4/#10), the script **also** computes an honest per-reaction biomass-flux growth
delta (`…growth.csv`).

## 6. Findings on the KBUtilLib functions

The four functions are correct in intent and shape, but the review found **22 verified
defects/gaps** (3 candidate findings refuted). The execution blockers (S2/S3 above),
the fake `growth_change` metric, and the media-`None` rich≡minimal duplication are the
most consequential. Full list, severities, and suggested fixes:
`KBUtils_Local/agent-io/audits/2026-07-06-ms-template-utils-review.md`.

## 7. Caveats / next steps

- Direction maps key on **base** ids; the script intersects with the loaded model, so
  only in-model reactions are perturbed (31 for E. coli; the map-wide count is far
  larger). Choose `--new dgpredictor` for a larger per-model edit set (dGPredictor
  diverges most from the cascade).
- `--baseline cascade_live` (map-vs-map) can emit rows where the model's *on-disk*
  bounds already equal the target (visible as `#cats=0`, `Δgrowth=0`); this reflects a
  real divergence between on-disk core-model bounds and the cascade map (cf. notebook
  09). Use `--baseline model` to attribute strictly against the model's actual bounds
  (57 edits for E. coli).
- To evaluate a full ModelSEED **template** (thousands of reactions, both GP/GN
  biomass, real Biolog panels) rather than a core model, `build_full_template_model` /
  `evaluate_template_quality` need modelseedpy **with** `ObjectivePkg` + a KBase token —
  out of scope here; the core model is the runnable, demo-sized target.
