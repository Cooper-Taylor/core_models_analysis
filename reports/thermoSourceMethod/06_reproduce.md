# 6. Reproducing, provenance, and file inventory

## 6.1 Pipeline

Run in order from `core_models_analysis/scripts/`. Each stage writes the inputs
the next one reads.

```bash
cd /scratch/ctaylor/core_models_analysis/scripts

python3 grade_thermo_sources.py            # ~4 min  → results/thermo_grades/
python3 build_graded_direction_maps.py     # ~2 min  → results/thermo_grades_fba/
python3 run_graded_fba_all_models.py --workers 32   # ~4 min, 39,781 LP solves
python3 analyze_graded_fba.py              # ~3 min
python3 analyze_implicit_directions.py     # ~2 min  (the implicit baseline, §5.3)
python3 plot_graded_fba.py                 # ~20 s

python3 recommend_thermo_source.py         # ~8 min  → results/thermo_recommendation/
```

Useful flags:

```bash
python3 grade_thermo_sources.py --tecrdb-skeleton-gold      # skeleton matches GOLD, not SILVER
python3 run_graded_fba_all_models.py --limit 40 --workers 8 # smoke test
python3 recommend_thermo_source.py --target direction --tolerance 0.2
```

## 6.2 Environment

| variable | default | what it selects |
|---|---|---|
| `MSDB_ROOT` | `/scratch/ctaylor/tmp/devsnap2` | the **data** snapshot |
| `MSDB_CODE` | `/scratch/ctaylor/ModelSEEDDatabase` | the **cascade code** checkout |
| `CORE_MODELS_ANALYSIS_DIR` | `/scratch/ctaylor/core_models_analysis` | output root |
| `TECRDB_COMPARISON` | `/scratch/ctaylor/dgpredictor_tecrdb/results/tecrdb_vs_dgpredictor_modelseed.csv` | the experimental reference |
| `GRADES_OUT`, `RECOMMEND_OUT` | under `results/` | where each stage writes |

Python: `numpy`, `pandas`, `scikit-learn` (isotonic regression), `cobra`
(GLPK), `matplotlib`. The default `core_models_analysis` env has all of these.
RDKit is **not** needed — the structure matching was done upstream, in the
`eq3` env, and enters here as a finished CSV.

`build_graded_direction_maps.load_reactions` re-implements the only transform
`BiochemPy.Reactions.loadReactions` applies (`None` → the string `"null"`),
because the dev archive ships no `Libs/Python`. If that ever diverges, the
cascade's `notes` handling is what breaks first.

## 6.3 Provenance

| item | value |
|---|---|
| reaction/compound data | ModelSEED `dev` @ **49563c6f**, archived to `/scratch/ctaylor/tmp/devsnap2` |
| Group Contribution vintage | Convention A rebuild, dev `ad34d6ab` (2026-08-07) — 27,313 non-sentinel, median σ 10.28 |
| dGPredictor-ModelSEED vintage | the ModelSEED-retrained model, 31,924 reactions, keyed by `rxnNNNNN` |
| cascade code | local checkout, `DEFAULT_HEURISTICS` verified identical in order to the snapshot's |
| KEGG mask | **not applied** — the legacy KEGG-keyed `dGPredictor` label is never read |
| TECRDB | eQuilibrator Zenodo doi:10.5281/zenodo.3978440, 4,544 rows → 802 stereo-exact + 748 skeleton matches |
| core models | 5,683 Kegg2 models, `data/core_models_kegg2/` |
| media | `ModelSEEDDatabase/Media/KBaseMedia.cpd`, KBase complete |
| run date | 2026-08-12 |

Each output stamps its own provenance:
`direction_maps_summary.json` → `msdb_data` / `msdb_code`;
`grade_calibration.json` → `msdb_root`;
`manifest.json` → model count, worker count, elapsed time.

**Nothing in `ModelSEEDDatabase` or `core_models_kegg2` is modified by any of
these scripts.** Every filter is applied at read time.

## 6.4 Outputs

### `results/thermo_grades/` — the grades

| file | rows | contents |
|---|---:|---|
| `source_grades.tsv` | 80,335 | long form, one row per (reaction × source): `rxn, name, ec, source, dg, sigma, operator, ehat, p_ok, z, birge, n_src, struct_zero, grade, reason, n_anchor, n_proxy` |
| `source_grades_heldout.tsv` | 78,785 | same, TECRDB removed and Rule 1 disabled |
| `source_grades_wide.tsv` | 56,002 | one row per reaction, one grade column per source, plus `best_grade` / `best_source` |
| `grade_calibration.json` | — | fitted ĥ curves, thresholds, veto counts, the §3.5 validation table |
| `grade_frontier.tsv` | 4 | grade and reason counts per source |

### `results/thermo_grades_fba/` — the simulations

| file | rows | contents |
|---|---:|---|
| `rxn_directions_<variant>.json` | — | `{rxn_id: operator}` for each of the six mapped variants |
| `rxn_source_coverage.csv` | 56,002 | per reaction: `has_/op_/status_` for every variant, plus which source each graded variant picked |
| `direction_maps_summary.json` | — | coverage, operator mix, graded source mix, provenance |
| `model_results.csv` | 5,683 | **one row per model**: inventory, per-variant direction coverage, overrides, bounds changed, FBA status, growth flux, grow/no-grow |
| `summary_stats.json` | — | combined unions, growth totals, median bounds changed |
| `manifest.json` | — | run metadata |
| `variant_growth.tsv` | 7 | §5.4, one row per variant |
| `variant_agreement.tsv` | 21 | §5.5, all pairs |
| `direction_accuracy.tsv` | 24 | §5.6, 6 variants × 4 subsets |
| `core_reaction_grades.tsv` | 239 | §5.7, the per-core-reaction table |
| `implicit_directions.tsv` | 239 | §5.3, the models' native direction per core reaction vs every variant |
| `implicit_summary.json` | — | §5.3, the aggregate implicit tables |

### `results/thermo_recommendation/` — the recommendations

| file | rows | contents |
|---|---:|---|
| `recommendation_direction.tsv` | 56,002 | per-source ΔG and risk, the choice, `recommended_dg`, `recommended_sigma`, `risk`, `n_feasible`, `kept` |
| `recommendation_magnitude.tsv` | 56,002 | same, under the argmin-ê rule |
| `recommendation_models.json` | — | the τ calibration (*k*ₛ, coverage before/after) and the full §4.5 ablation |

### This folder

`tables/` holds copies of every summary table above (including
`implicit_directions.tsv` and `implicit_summary.json`), so it can be read without
the results tree. `figures/` holds the three figures. The per-model
`model_results.csv` (5,683 rows) and the long `source_grades.tsv` (80,335 rows)
are left in `results/` rather than copied, being too large to be useful here.

## 6.5 Consumer entry points

```python
from grade_thermo_sources import load_grades, recommended_energy_map
from recommend_thermo_source import load_recommendation

load_grades()                                  # long form
load_grades(wide=True)                         # one row per reaction
recommended_energy_map(min_grade="SILVER")     # {rxn: (dg, sigma, source)} for a direction map
recommended_energy_map(heldout=True)           # no TECRDB, for unbiased scoring
load_recommendation("direction")               # the §4.6 rule's output
```

`recommended_energy_map` returns exactly what `build_graded_direction_maps.py`
feeds the cascade, so a new variant is a two-line change.

## 6.6 Related reports

| report | what it adds |
|---|---|
| `reports/thermoComparison/THERMO_SOURCE_GRADING_PROPOSAL.md` | the grading design as originally proposed, with the rejected alternatives |
| `reports/thermoComparison/THERMO_SOURCE_RECOMMENDER.md` | the recommender in narrative form |
| `reports/thermoComparison/GRADED_SOURCE_CORE_MODEL_ANALYSIS.md` | the simulations in narrative form |
| `reports/thermoComparison/THERMO_SOURCE_ASSIGNMENT.md` | the earlier single-source assignment (ê), which the magnitude target reuses |
| `reports/thermoComparison/EQUILIBRATOR_VS_DGPREDICTOR_MODELSEED.md` | where the σ-calibration and quinone findings come from |
| `reports/thermoComparison/THERMO_SOURCE_AGREEMENT_STRUCTURE.md` | the KEGG mis-mapping defect and the chemistry of source disagreement |
| `reports/thermoComparison/THERMO_SOURCE_FBA_PIPELINE.md` | the 2026-08-03 sweep this one supersedes |
