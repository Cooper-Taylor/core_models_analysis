# 5. The core-model simulations

Flux-balance analysis over all 5,683 Kegg2 core models under each thermodynamic
source in isolation and under the graded source. Complete results.

Scripts: `build_graded_direction_maps.py`, `run_graded_fba_all_models.py`,
`analyze_graded_fba.py`, `analyze_implicit_directions.py`, `plot_graded_fba.py`.
Data: `results/thermo_grades_fba/`; copies in [`tables/`](tables/).

---

## 5.1 Design

Seven variants over the same 5,683 models, same media, same objective, same
overlay policy. **Only the ΔG′° feeding the cascade varies.**

| variant | direction from |
|---|---|
| `implicit` | no override at all — the bounds baked into the model file |
| `gc` | Group Contribution's own ΔG′°, nothing else |
| `eq` | eQuilibrator's own ΔG′°, nothing else |
| `dgpms` | dGPredictor-ModelSEED's own ΔG′°, nothing else |
| `graded` | per reaction, the ΔG′° of the best-graded source (any grade) |
| `graded_trusted` | same, but reactions whose best grade is BRONZE get no call |
| `graded_heldout` | same as `graded` with TECRDB removed and Rule 1 disabled — the only graded variant scoreable against TECRDB without circularity |

**Held constant.**

- *The cascade.* `DEFAULT_HEURISTICS`, unmodified, for all seven. No LLM
  heuristics and no eQuilibrator reversibility index.
- *The energy accessor.* `per_source_energy(label)` — that source's own ΔG′° —
  not `top_level_energy`, which evaluates the shared canonical ΔG′° gated only
  by source eligibility and, in its eQuilibrator mode, silently falls back to
  Group Contribution's reversibility. Either would blend sources and break
  "direction from source X alone".
- *The overlay policy.* `growth_heuristics.override_bounds` rewrites bounds only
  for reactions the variant has an opinion about; every other reaction keeps its
  native on-disk bound. This is what makes single-source variants meaningful —
  no other source is substituted into the gaps.
- *The FBA.* `apply_media` restricts exchange uptake to KBase complete media;
  `find_biomass_reaction` + `model.optimize()`, cobra's default GLPK solver.
  Each model is loaded once and copied per variant, so the seven runs never
  interfere.

5,683 models × 7 variants = **39,781 LP solves**, 32 workers, **0 errors**, all
39,781 returning solver status `optimal`.

## 5.2 Coverage

### Database-wide direction maps

| variant | reactions with a direction | `>` | `<` | `=` |
|---|---:|---:|---:|---:|
| Group Contribution | 27,313 | 9,579 | 1,413 | 16,321 |
| eQuilibrator | 25,028 | 10,583 | 1,590 | 12,855 |
| dGPredictor-ModelSEED | 31,924 | 9,922 | 1,231 | 20,771 |
| graded | **33,289** | 12,671 | 1,862 | 18,756 |
| graded, SILVER floor | 23,103 | 8,405 | 948 | 13,750 |
| graded, TECRDB held out | 33,289 | 12,705 | 1,894 | 18,690 |

### Which source the graded map used

| | database-wide | on the 239 core reactions |
|---|---:|---:|
| Group Contribution | 12,031 | 46 |
| eQuilibrator | 10,128 | 27 |
| dGPredictor-ModelSEED | 9,842 | 69 |
| **TECRDB** | 1,288 | **67** |

Roughly a third of the graded map's core-reaction calls rest on experimental
data, which is not true anywhere else in the database (§1.5).

### Per model

| variant | min | median | mean | max |
|---|---:|---:|---:|---:|
| unique reactions in the model | 20 | 128 | 123.1 | 187 |
| unique compounds in the model | 41 | 124 | 119.2 | 163 |
| — with a direction under `gc` | 18 | 114 | 109.9 | 165 |
| — under `eq` | 2 | 92 | 88.2 | 141 |
| — under `dgpms` | 19 | 117 | 112.8 | 172 |
| — under `graded` | 19 | 118 | 113.8 | 173 |
| — under `graded_trusted` | 19 | 107 | 102.5 | 151 |

### Combined across all models, of the 239 core reactions

| variant | core reactions with a direction |
|---|---:|
| Group Contribution | 200 |
| eQuilibrator | 173 |
| dGPredictor-ModelSEED | 208 |
| graded | **209** |
| graded, SILVER floor | 179 |
| graded, TECRDB held out | 209 |

### Bounds actually changed, per model

Distinct from "touched": the count of reactions whose (lower, upper) bound pair
differs from the model's native one after the overlay.

| variant | min | median | mean | max |
|---|---:|---:|---:|---:|
| `implicit` | 0 | 0 | 0.0 | 0 |
| `gc` | 5 | 39 | 37.8 | 62 |
| `eq` | 0 | 21 | 20.0 | 39 |
| `dgpms` | 5 | 40 | 38.8 | 66 |
| `graded` | 5 | 37 | 36.0 | 61 |
| `graded_trusted` | 5 | 28 | 27.2 | 44 |
| `graded_heldout` | 5 | 38 | 36.3 | 61 |

## 5.3 The implicit baseline — what the models ship with

`implicit` applies no override at all: each model runs on the reaction bounds
baked in when it was built and gap-filled. It is the reference every other
variant is measured against, so it is worth knowing what those bounds actually
say rather than treating them as a black box.

Script: `analyze_implicit_directions.py`. Output:
`results/thermo_grades_fba/implicit_directions.tsv`, `implicit_summary.json`.
A reaction's native operator is read from its bounds — `lb < 0 < ub` → `=`,
`lb ≥ 0 < ub` → `>`, `lb < 0, ub ≤ 0` → `<`, otherwise blocked.

### The native bounds are a single global map

Across all 699,422 (model, reaction) pairs:

| native operator | pairs | share |
|---|---:|---:|
| `=` reversible | 358,244 | 51.2% |
| `>` forward | 272,917 | 39.0% |
| `<` reverse | 62,515 | 8.9% |
| blocked | 5,746 | 0.8% |

**All 239 core reactions carry a unanimous native direction across all 5,683
models** — no reaction is forward in one model and reversible in another. The
implicit bounds are therefore not 5,683 independent decisions but one map
applied everywhere, presumably inherited from the template the models were built
from. That is what makes it scoreable like any other direction source.

Per reaction, that map is:

| modal native operator | core reactions | share |
|---|---:|---:|
| `=` reversible | 108 | 45.2% |
| `>` forward | 104 | 43.5% |
| `<` reverse | 25 | 10.5% |
| blocked | 2 | 0.8% |

So the shipped bounds are *less* permissive than every thermodynamic variant
except eQuilibrator: 45.2% reversible against Group Contribution's 82.5% and the
graded map's 76.6%.

### The native bounds are the least accurate direction source tested

Scored against the experimental reference Λ\* exactly as in §5.6, on the 65 core
reactions with a `stereo_exact` TECRDB match:

| direction source | correct | accuracy |
|---|---:|---:|
| **implicit (native bounds)** | **44 / 65** | **67.7%** |
| Group Contribution | 59 / 65 | 90.8% |
| dGPredictor-ModelSEED | 63 / 65 | 96.9% |
| graded, TECRDB held out | 63 / 65 | 96.9% |
| eQuilibrator | 64 / 65 | 98.5% |

Every thermodynamic source beats the bounds the models were shipped with, by 23
to 31 points.

**And the errors are one-sided.** Of the 21 mismatches:

| native says | experiment says | n |
|---|---|---:|
| `>` forward | `=` reversible | **14** |
| `<` reverse | `=` reversible | **5** |
| `=` reversible | `<` reverse | 2 |

19 of 21 are the model **over-constraining** — forcing a direction on a reaction
the thermodynamics call reversible — against 2 in the opposite direction. That
asymmetry is the mechanism behind §5.4: every thermodynamic variant grows *more*
models than `implicit` (+133 to +279) because it is mostly relaxing constraints
that were not thermodynamically justified in the first place.

### Where the variants move it

Agreement with the native bound, over core reactions each variant covers:

| variant | agrees | of | share |
|---|---:|---:|---:|
| eQuilibrator | 127 | 173 | 73.4% |
| graded, SILVER floor | 126 | 177 | 71.2% |
| graded | 136 | 207 | 65.7% |
| Group Contribution | 125 | 198 | 63.1% |
| dGPredictor-ModelSEED | 130 | 206 | 63.1% |

Native → graded transitions on the core set:

| transition | n | |
|---|---:|---|
| `=` → `=` | 93 | agreed reversible |
| `>` → `>` | 32 | agreed forward |
| `<` → `<` | 11 | agreed reverse |
| **`>` → `=`** | **54** | **relaxed** — the dominant change |
| **`<` → `=`** | **11** | **relaxed** |
| `=` → `<` | 3 | tightened |
| `=` → `>` | 1 | tightened |
| `<` → `>` | 2 | reversed |

65 relaxations against 4 tightenings and 2 reversals. The graded map's effect on
a core model is overwhelmingly to *remove* a direction constraint, which is both
why growth goes up and why growth is a poor way to judge it (§5.4, below).

## 5.4 Growth

![growth](figures/fig1_growth.png)

| variant | models growing | % of 5,683 | gained vs `implicit` | lost vs `implicit` |
|---|---:|---:|---:|---:|
| model's own bounds | 3,461 | 60.9% | — | — |
| Group Contribution | 3,656 | 64.3% | +206 | −11 |
| eQuilibrator | 3,570 | 62.8% | +133 | −24 |
| dGPredictor-ModelSEED | 3,717 | 65.4% | +279 | −23 |
| **graded** | **3,715** | **65.4%** | +277 | −23 |
| graded, SILVER floor | 3,689 | 64.9% | +251 | −23 |
| graded, TECRDB held out | 3,715 | 65.4% | +277 | −23 |

Every thermodynamic variant grows more models than the models' own shipped
bounds. The total spread across all seven is 256 models, 4.5 percentage points.

### Growth flux

Non-growers contribute a flux of 0 to the "all models" columns.

| variant | median (all) | mean (all) | median (growers) | mean (growers) |
|---|---:|---:|---:|---:|
| implicit | 32.18 | 31.91 | 52.30 | 52.39 |
| Group Contribution | 55.80 | 51.20 | 85.22 | 79.59 |
| eQuilibrator | 31.39 | 25.86 | 38.17 | 41.17 |
| dGPredictor-ModelSEED | 35.98 | 33.31 | 49.45 | 50.93 |
| graded | 34.56 | 30.67 | 40.63 | 46.91 |
| graded, SILVER floor | 34.29 | 30.45 | 40.63 | 46.90 |
| graded, TECRDB held out | 34.56 | 30.66 | 40.63 | 46.90 |

Group Contribution's growers carry roughly twice eQuilibrator's flux, which is
the same phenomenon as the growth count — see below.

### Growth counts do not rank the sources

| variant | share of core direction calls that are reversible `=` |
|---|---:|
| Group Contribution | 82.5% |
| dGPredictor-ModelSEED | 78.4% |
| graded | 76.6% |
| graded, SILVER floor | 76.0% |
| graded, TECRDB held out | 75.6% |
| eQuilibrator | 63.6% |

A map that calls more reactions reversible removes more constraints and grows
more models whether or not it is right. **eQuilibrator grows the fewest models
and is the most accurate against experiment (§5.6), so on this data growth count
and correctness point in opposite directions.** Any claim of the form "source X
is better because more models grow under it" is unsupported. Report the
permissiveness column alongside, always.

## 5.5 Pairwise agreement

Number of models whose grow/no-grow verdict differs, and in whose favour. All 21
pairs; full table in [`tables/variant_agreement.tsv`](tables/variant_agreement.tsv).

| a | b | differ | a grows only | b grows only |
|---|---|---:|---:|---:|
| implicit | gc | 217 | 11 | 206 |
| implicit | eq | 157 | 24 | 133 |
| implicit | dgpms | 302 | 23 | 279 |
| implicit | graded | 300 | 23 | 277 |
| implicit | graded_trusted | 274 | 23 | 251 |
| implicit | graded_heldout | 300 | 23 | 277 |
| gc | eq | 86 | **86** | **0** |
| gc | dgpms | 95 | 17 | 78 |
| gc | graded | 93 | 17 | 76 |
| gc | graded_trusted | 115 | 41 | 74 |
| gc | graded_heldout | 93 | 17 | 76 |
| eq | dgpms | 147 | **0** | **147** |
| eq | graded | 145 | **0** | **145** |
| eq | graded_trusted | 119 | 0 | 119 |
| eq | graded_heldout | 145 | 0 | 145 |
| dgpms | graded | **2** | 2 | 0 |
| dgpms | graded_trusted | 28 | 28 | 0 |
| dgpms | graded_heldout | 2 | 2 | 0 |
| graded | graded_trusted | 26 | 26 | 0 |
| graded | graded_heldout | **0** | 0 | 0 |
| graded_trusted | graded_heldout | 26 | 0 | 26 |

Three things stand out.

**`graded` and `dgpms` differ on 2 of 5,683 models.** dGPredictor-ModelSEED is
GOLD on 78 core reactions — more than any other predictor — so it wins the
graded pick most often, and where it loses the winner usually agrees with it
anyway (10 direction disagreements over 208 co-covered core reactions).

**Growth differences are nested, not scattered.** `eq` vs `graded` differ on 145
models and `graded` grows all 145; `gc` vs `eq` differ on 86 and `gc` grows all
86. One variant's growers are essentially a superset of the other's, which is the
signature of a permissiveness difference rather than a substantive disagreement.

**`graded` and `graded_heldout` differ on 0 models**, despite differing on 1,288
reactions' worth of source assignment. Removing the experimental data changes
which source is used but almost never changes the resulting direction, because
where TECRDB exists the predictors usually agree with it.

## 5.6 Direction accuracy against experiment

![direction accuracy](figures/fig2_direction_accuracy.png)

**Method.** For each reaction with a `stereo_exact` TECRDB match, run the *same
cascade* on the *experimental* ΔG′° to obtain the reference direction Λ\*
([§2.4](02_notation.md#24-the-cascade)) — same heuristics, same concentration
assumptions, differing only in that the energy is measured. Then score each
variant's operator against it.

`graded` and `graded_trusted` are excluded as **circular**: they use TECRDB, so
they reproduce Λ\* perfectly by construction (802/802). `graded_heldout` exists
precisely so this comparison can be made.

| variant | all 802 | reference DIRECTIONAL (155) | reference `=` (647) | core reactions (65) |
|---|---:|---:|---:|---:|
| Group Contribution | 85.5% | **51.0%** | 93.8% | 90.8% |
| **eQuilibrator** | **95.5%** | **93.5%** | **95.9%** | **98.5%** |
| dGPredictor-ModelSEED | 91.8% | 85.8% | 93.2% | 96.9% |
| graded, TECRDB held out | 93.4% | 89.0% | 94.4% | 96.9% |
| graded *(circular)* | 100% | 100% | 100% | 100% |
| graded, SILVER floor *(circular)* | 100% | 100% | 100% | 100% |

Full table: [`tables/direction_accuracy.tsv`](tables/direction_accuracy.tsv).

**Group Contribution's 51.0% on directional reactions is the finding to carry
forward.** On reactions the experiment says are one-way, it is at chance. Its
93.8% on reversible reactions is largely the mMΔG band absorbing a wrong number —
returning `=` for the wrong reason. This is consistent with its grade
distribution (23 GOLD of 239 core reactions) and with σ_GC correlating with real
error at only ρ = +0.176.

The reference's own operator mix is 647 `=`, 124 `>`, 31 `<`, which is why the
directional subset is the discriminating one: a variant can score in the low 90s
overall while being at chance on the reactions that actually constrain a model.

## 5.7 The core reaction set

![core grades](figures/fig3_core_grades.png)

239 distinct reactions appear in at least one core model; 233 of them are
non-EMPTY in the snapshot. Grades:

| source | GOLD | SILVER | BRONZE | no data |
|---|---:|---:|---:|---:|
| TECRDB | 65 | 4 | 0 | 170 |
| eQuilibrator | 61 | 74 | 34 | 70 |
| dGPredictor-ModelSEED | 78 | 84 | 20 | 57 |
| Group Contribution | 23 | 88 | 89 | 39 |

Per reaction, the best grade available: **103 GOLD, 76 SILVER, 30 BRONZE**, and
30 core reactions with no thermodynamic source at all.

Direction disagreements between the graded map and each single source, over
reactions both cover:

| pair | disagreements | co-covered |
|---|---:|---:|
| graded vs Group Contribution | 19 | 200 |
| graded vs eQuilibrator | 19 | 173 |
| graded vs dGPredictor-ModelSEED | 10 | 208 |

The complete per-reaction table —
[`tables/core_reaction_grades.tsv`](tables/core_reaction_grades.tsv), 239 rows —
carries for each core reaction: how many models contain it, name and EC, the
grade of all four sources, the best grade and best source, the Birge ratio, the
operator under every variant, which source the graded map picked, and the TECRDB
direction and match tier where one exists. That is the table to use for
case-by-case curation.

## 5.8 What to conclude

1. **Coverage is where the graded map wins.** 33,289 reactions and 209 of 239
   core reactions, against eQuilibrator's 25,028 and 173 — and every call
   arrives with a grade attached, so a consumer can decline the bad ones.
2. **It does not win on per-reaction direction accuracy.** Held out, 93.4%
   against eQuilibrator's 95.5%. §4.2 explains why, and §4.6 is the rule that
   fixes it.
3. **Do not use Group Contribution alone to set direction.** 51% on directional
   reactions is chance.
4. **Do not read growth counts as a quality ranking** without §5.4's
   permissiveness column beside them.
5. **Use `graded_trusted` when a wrong direction costs more than a missing one.**
   It gives up 26 growers and 10,186 reactions to exclude the BRONZE tier.
