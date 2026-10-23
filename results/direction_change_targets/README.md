# Direction-change target reactions: presence and FBA activity

76 ModelSEED reactions whose direction changes — 56 with an explicit change
(`= -> >` or `= -> <`) and 20 labelled "no change" in the source list but which
change because of thermodynamic uncertainty. All 76 are currently reversible
(`=`), so every change is reversible → irreversible.

Measured across all three model panels with the repo's own `cma` model-set API
(model sets, media, biomass policy, seed-key mapping). **No ModelSEEDpy fbapkg
code is involved.**

## Headline

| panel | models | growers | targets present | targets active |
|---|---:|---:|---:|---:|
| core (`core_kegg2`, kbase_complete) | 5,683 | 3,461 | **5 / 76** | **3 / 76** |
| full / GSM (`ms2_gsm`, glucose minimal) | 5,420 | 4,347 | **36 / 76** | **27 / 76** |
| full / GSM (`ms2_gsm_auxo`, auxotrophy) | 5,420 | 5,397 | **35 / 76** | **27 / 76** |
| **union** | 16,523 | | **36 / 76** | **30 / 76** |

40 of the 76 appear in no model at all. 6 more are present but never carry
flux. Split by category: explicit change 27/56 present and 22/56 active;
"no change" 9/20 present and 8/20 active.

## Activity is measured on pFBA, not plain FBA

Plain FBA picks an arbitrary vertex among alternate optima, so whether a given
reaction carries flux in it is not reproducible. pFBA is used throughout, the
same convention as `build_flux_distributions.py`, with `ACTIVE_TOL = 1e-9`.
The plain-FBA answer is recorded alongside in `active_fba`: it gives 31 rather
than 30, the extra one (rxn15230) being an alternate-optimum artefact.

## Files

| file | rows | what |
|---|---:|---|
| `direction_change_targets_by_reaction.tsv` | 76 | one row per target reaction, per-panel present/growing/active counts and flags |
| `figures/direction_change_funnel.png` | — | the three-stage drop, split by category |
| `figures/direction_change_by_reaction.png` | — | per-reaction model counts, three facets |
| `direction_change_targets_by_model.tsv` | 154,936 | one row per (panel, model, reaction) where the reaction is present, flagged for flux |
| `presence_activity.json` | — | the raw aggregate the TSVs are built from |

`per_model_records.json` is the intermediate the by-model TSV is built from.
It is gitignored; regenerate with `--dump-records`.

### `direction_change_targets_by_reaction.tsv`

`reaction_id`, `change`, `name` (from MSDB), then a
`present_<panel>` / `growing_<panel>` / `active_<panel>` triple for each of
`core`, `gsm`, `gsm_auxo`, then `present_total`, `growing_total`,
`active_total` summed over panels, then `present_any`, `active_any` (0/1),
`n_panels_present`, `n_panels_active`.

The three counts nest: **active ≤ growing ≤ present**, checked on every row.

* `present_<panel>` — models that contain the reaction at all.
* `growing_<panel>` — of those, how many have non-zero biomass. A non-grower
  has no flux solution, so this is the denominator `active` should be read
  against, not `present`. The gap matters: rxn00011 is in 5,352 core models
  but only 3,394 of them grow.
* `active_<panel>` — of those, how many carry flux through it under pFBA.

### `direction_change_targets_by_model.tsv`

`panel`, `model_id`, `grows`, `reaction_id`, `change`, `present`,
`active_pfba`, `active_fba`. Only rows where the reaction is present in that
model are emitted.

## Figures

Registered in `scripts/figures.tsv` under tag `direction_change`:

```bash
python3 scripts/regen_figures.py --tag direction_change
```

* `direction_change_funnel` — all 76 → present in ≥1 model (36) → active in
  ≥1 pFBA solution (30), each stage split into the 56 explicit changes and the
  20 labelled "no change".
* `direction_change_by_reaction` — one row per reaction, three facets; pale bar
  is models containing it, solid bar the subset running it.

Each writes a `.tsv` of the numbers it asserts beside the image.

## Reproduce

```bash
python3 scripts/direction_change/presence_and_activity.py -j 32 --dump-records
python3 scripts/direction_change/write_tsv.py
python3 scripts/direction_change/summarise.py
```

Roughly 5 minutes for all 16,523 models at `-j 32`. Deterministic: a second
run reproduced every count.

## Validation

- The `cma` seed-key mapping (model reaction id → MSDB `rxn` id) was checked
  against a raw-JSON scan of the same 40 core models and agreed exactly.
- Counts agree with the precomputed reaction universes in
  `modelseed2_gs_models/`: 5/76 against `reaction_universe_core.txt` (233
  reactions), 36/76 against `reaction_universe.txt` (3,769).
- Reaction orientation was checked against MSDB's canonical stoichiometry in
  150 core and 60 GSM models: **0 flipped**. This matters because the `= -> >`
  labels are written in MSDB's frame.
- 0 solver errors across all 16,523 models.

## Caveat: "active" is an upper bound on impact

Because all 76 go reversible → irreversible, a reaction is only genuinely
constrained if it carries flux **in the direction being forbidden**. A reaction
already running the permitted way is unaffected. The orientation check above is
the prerequisite for that calculation, but it has not been run — so treat
30/76 as the set of candidates, not the set of reactions that will change flux.
