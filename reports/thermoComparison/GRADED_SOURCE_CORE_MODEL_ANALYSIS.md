# Core-model analysis, redone per thermodynamic source and under the graded source

*Reaction direction determined from Group Contribution alone, eQuilibrator alone,
dGPredictor-ModelSEED alone, and from the gold/silver/bronze graded pick — then
FBA over all 5,683 Kegg2 core models for each. 2026-08-12.*

Scripts: [`grade_thermo_sources.py`](/scratch/ctaylor/core_models_analysis/scripts/grade_thermo_sources.py),
[`build_graded_direction_maps.py`](/scratch/ctaylor/core_models_analysis/scripts/build_graded_direction_maps.py),
[`run_graded_fba_all_models.py`](/scratch/ctaylor/core_models_analysis/scripts/run_graded_fba_all_models.py),
[`analyze_graded_fba.py`](/scratch/ctaylor/core_models_analysis/scripts/analyze_graded_fba.py),
[`plot_graded_fba.py`](/scratch/ctaylor/core_models_analysis/scripts/plot_graded_fba.py).
Data: `results/thermo_grades/`, `results/thermo_grades_fba/`.
Figures: `reports/thermoComparison/figures/graded_fba/`.
Grading method: [`THERMO_SOURCE_GRADING_PROPOSAL.md`](THERMO_SOURCE_GRADING_PROPOSAL.md).

---

## Headline

**On the reactions where an experiment can check us, no arbitration scheme beats
"use eQuilibrator when it has an answer."** eQuilibrator alone reproduces the
experimental direction on 95.5% of 802 matched reactions; the graded pick with
the experiment held out gets 93.4%, dGPredictor-ModelSEED 91.8%, Group
Contribution 85.5%. On the hard subset where the experiment says the reaction is
*directional*, the spread widens: 93.5% / 89.0% / 85.8% / **51.0%**.

What the graded source buys is not per-reaction accuracy but **coverage carried
with a quality label**: 33,289 reactions get a direction versus eQuilibrator's
25,028, 209 of the 239 core reactions versus 173, and every one of them arrives
tagged GOLD, SILVER or BRONZE so a downstream consumer can decline the bad ones.
It also uses the measured value outright on the 1,288 reactions that have one.

Growth counts separate the variants by only 4.5 percentage points and, as §4
shows, track permissiveness more than correctness. They are reported because they
were asked for, not because they rank the sources.

---

## 1. What was run, and what was held constant

Six FBA variants over the same 5,683 models, same media, same biomass objective,
same overlay policy:

| variant | direction from |
|---|---|
| `implicit` | no override — the bounds baked into the model file |
| `gc` | **Group Contribution only** |
| `eq` | **eQuilibrator only** |
| `dgpms` | **dGPredictor-ModelSEED only** |
| `graded` | **the best-graded source per reaction** (any grade) — the recommended map |
| `graded_trusted` | same, but reactions whose best grade is BRONZE get no call |
| `graded_heldout` | same as `graded` with TECRDB removed and the measurement-derived grade disabled — the only graded variant that can be scored against TECRDB without circularity |

**Held constant: the cascade.** All six use `DEFAULT_HEURISTICS` from
`reversibility_heuristics.py` unmodified — ATP synthase → ABC transporter →
stored ΔG bounds → mMΔG band → low-energy points → default. The only thing that
varies is which ΔG′° feeds in. No LLM heuristics, no eQuilibrator reversibility
index.

**Held constant: the overlay policy.** `growth_heuristics.override_bounds`
rewrites bounds only for reactions the variant has an opinion about; everything
else keeps its native on-disk bound. That is what makes "direction from source X
alone" meaningful — no other source is silently substituted into the gaps.

### Snapshot change from the 2026-08-03 sweep — read this before comparing numbers

The earlier sweep read the live `/scratch/ctaylor/ModelSEEDDatabase` working
tree. This one reads **dev @ 49563c6f** (`/scratch/ctaylor/tmp/devsnap2`), because
the live tree has neither of the two things this analysis needs:

- **`dGPredictor-ModelSEED`**, the retrain keyed directly by ModelSEED reaction
  id. The live tree only has the legacy KEGG-keyed `dGPredictor`, 62% of whose
  records were predicted from a KEGG reaction ModelSEED does not list for them.
  The retrain is structurally immune to that, so **no KEGG mask is applied here** —
  it is not needed. Coverage goes from 141 core reactions (masked legacy source)
  to 208.
- **Group Contribution rebuilt under Convention A** (dev `ad34d6ab`): 53% of
  values changed, coverage 25,812 → 27,313, σ roughly doubled.

Cascade code still comes from the local checkout; its `DEFAULT_HEURISTICS` order
is byte-identical to the snapshot's, verified before the run. Old-run growth
totals were GC 3,642 / eQ 3,570 / legacy-dGP 3,546 / implicit 3,461, so
eQuilibrator and implicit are unchanged and the two moved sources moved for the
documented reasons.

---

## 2. Coverage

| variant | reactions with a direction | of the 239 core reactions | median bounds changed per model |
|---|---:|---:|---:|
| Group Contribution | 27,313 | 200 | 39 |
| eQuilibrator | 25,028 | 173 | 21 |
| dGPredictor-ModelSEED | 31,924 | 208 | 40 |
| **graded (recommended)** | **33,289** | **209** | 37 |
| graded, SILVER floor | 23,103 | 179 | 28 |

The graded map is the union, so it covers more than any single source. The SILVER
floor deliberately gives back 10,186 reactions — those whose best available
grade is BRONZE — on the grounds that a bad direction is worse than none.

**Which source the graded map actually used**, database-wide: Group Contribution
12,031, eQuilibrator 10,128, dGPredictor-ModelSEED 9,842, TECRDB 1,288. On the
core set: dGPredictor-ModelSEED 69, **TECRDB 67**, Group Contribution 46,
eQuilibrator 27.

**The core set is far better measured than the database as a whole.** 69 of the
239 core reactions (29%) carry a TECRDB measurement, against 1,550 of 56,002
(2.8%) database-wide — a 10× enrichment, because central metabolism is exactly
what NIST measured. So on core models the graded map is standing on experimental
data for roughly a third of its calls, which is not true anywhere else.

Grades on the core set (`fig3_core_grades.png`):

| source | GOLD | SILVER | BRONZE | no data |
|---|---:|---:|---:|---:|
| TECRDB | 65 | 4 | 0 | 170 |
| eQuilibrator | 61 | 74 | 34 | 70 |
| dGPredictor-ModelSEED | 78 | 84 | 20 | 57 |
| Group Contribution | 23 | 88 | 89 | 39 |

Per reaction, the best grade available is GOLD for 103, SILVER for 76, BRONZE for
30, and 30 core reactions have no thermodynamic source at all.

---

## 2b. The implicit baseline — what the models ship with

`implicit` runs each model on its own on-disk bounds. Worth knowing what those
say, since everything else is measured against them
(`analyze_implicit_directions.py` → `results/thermo_grades_fba/implicit_*`).

**All 239 core reactions carry a unanimous native direction across all 5,683
models**, so the shipped bounds are one global map, not 5,683 decisions. That
map is 45.2% reversible / 43.5% forward / 10.5% reverse / 0.8% blocked — less
permissive than every variant except eQuilibrator.

Scored against the experimental reference on the 65 scoreable core reactions,
**it is the least accurate direction source tested**:

| direction source | accuracy |
|---|---:|
| **implicit (native bounds)** | **67.7%** (44/65) |
| Group Contribution | 90.8% |
| dGPredictor-ModelSEED / graded held out | 96.9% |
| eQuilibrator | 98.5% |

**The errors are one-sided:** of 21 mismatches, 14 are native `>` where the
experiment says `=` and 5 are native `<` → `=`, against 2 the other way. 19 of
21 are the model over-constraining a reaction the thermodynamics call
reversible. That is the mechanism behind §3 — every thermodynamic variant grows
more models than `implicit` because it is mostly relaxing constraints that were
never thermodynamically justified. On the core set the graded map makes 65
relaxations against 4 tightenings and 2 reversals.

## 3. Growth

![growth](figures/graded_fba/fig1_growth.png)

| variant | models growing | % of 5,683 | gained vs implicit | lost vs implicit |
|---|---:|---:|---:|---:|
| model's own bounds | 3,461 | 60.9% | — | — |
| Group Contribution | 3,656 | 64.3% | +206 | −11 |
| eQuilibrator | 3,570 | 62.8% | +133 | −24 |
| dGPredictor-ModelSEED | 3,717 | 65.4% | +279 | −23 |
| **graded (recommended)** | **3,715** | **65.4%** | +277 | −23 |
| graded, SILVER floor | 3,689 | 64.9% | +251 | −23 |

Every thermodynamic variant grows more models than the models' own shipped
bounds, and the total spread across all six is 256 models (4.5 points).

Where they disagree:

| pair | models differing | |
|---|---:|---|
| graded vs eQuilibrator | 145 | graded grows all 145, eQuilibrator none |
| graded vs Group Contribution | 93 | graded grows 76, GC grows 17 |
| graded vs dGPredictor-ModelSEED | 2 | dGPredictor grows 2, graded none |
| graded vs graded (SILVER floor) | 26 | dropping BRONZE calls loses 26 growers |

**graded and dgpms are nearly the same map on core models** — they differ on 2 of
5,683. That is not a coincidence: dGPredictor-ModelSEED is graded GOLD on 78 core
reactions, more than any predictor, so it wins the pick most often, and where it
loses the winner usually agrees with it anyway (10 direction disagreements over
208 co-covered core reactions).

### Growth counts do not rank the sources

The right panel of the figure is the reason. The share of core direction calls
that are reversible `=` runs Group Contribution 82.5%, dGPredictor-ModelSEED
78.4%, graded 76.6%, eQuilibrator 63.6%. A map that calls more reactions
reversible removes more constraints and therefore grows more models regardless of
whether it is right. eQuilibrator grows the fewest models *and* is the most
accurate against experiment — so on this data growth count and correctness point
in opposite directions. Any claim of the form "source X is better because more
models grow" is unsupported.

---

## 4. Direction accuracy against experiment — the metric that does rank them

![direction accuracy](figures/graded_fba/fig2_direction_accuracy.png)

Method: for each reaction with a TECRDB match, run the **same cascade** on the
**experimental** ΔG′° instead of a predicted one. That yields a reference
direction on the same footing as every variant — same heuristics, same
concentration assumptions — differing only in that the energy is measured. Then
score each variant's operator against it. 802 stereo-exact matches.

`graded` and `graded_trusted` are excluded: they *use* TECRDB, so they reproduce
it perfectly by construction. `graded_heldout` is the scoreable one.

| variant | all 802 | reference DIRECTIONAL (155) | reference `=` (647) |
|---|---:|---:|---:|
| Group Contribution | 85.5% | **51.0%** | 93.8% |
| **eQuilibrator** | **95.5%** | **93.5%** | **95.9%** |
| dGPredictor-ModelSEED | 91.8% | 85.8% | 93.2% |
| graded, TECRDB held out | 93.4% | 89.0% | 94.4% |

Restricted to core-model reactions (65 with a stereo-exact match): eQuilibrator
98.5%, dGPredictor-ModelSEED 96.9%, graded-heldout 96.9%, Group Contribution
90.8%.

**Group Contribution's 51% on directional reactions is the finding to carry
forward.** On reactions the experiment says are one-way, it is at chance. Its
93.8% on reversible reactions is mostly the mMΔG band absorbing a wrong number —
getting `=` right for the wrong reason. This is consistent with its grade
distribution (23 GOLD of 239 core reactions) and with the earlier measurement
that its σ correlates with its real error at only ρ = +0.176.

---

## 5. Why the graded pick loses to eQuilibrator, and what that means

Held-out graded picks dGPredictor-ModelSEED on 521 of the 802 matched reactions
and eQuilibrator on only 274 — yet eQuilibrator is right 98.9% of the times it
*is* picked, against dGPredictor-ModelSEED's 90.4%. The ranking is choosing the
weaker source on this set.

The cause is a real dissociation, not a bug in the implementation:

> **The grade is calibrated on ΔG′° magnitude; direction is a different target.**
> `p_ok = P(|ΔG_s − ΔG*| ≤ 2 kcal/mol | σ)`. On the anchor set,
> dGPredictor-ModelSEED is genuinely the *better magnitude* predictor in its GOLD
> tier (98% within 2 kcal/mol vs eQuilibrator's 94%) — and the *worse direction*
> predictor. Direction errors concentrate near ΔG′° ≈ 0, where magnitude error is
> small, so optimising one does not optimise the other.

Testing the obvious fix — refit the calibration on
`1[direction(source) = direction(experiment)]` instead of on magnitude, 20 random
70/30 splits of the 802:

| ranking strategy | held-out direction accuracy |
|---|---|
| direction-calibrated grade | 95.2% ± 1.4% |
| magnitude-calibrated grade (shipped) | 93.3% ± 1.4% |
| eQuilibrator alone | 95.4% ± 1.0% |
| eQuilibrator first, graded as fallback | **95.5% ± 1.0%** |
| dGPredictor-ModelSEED alone | 91.7% ± 1.7% |

Direction-calibrating recovers most of the gap, and then the three leading
strategies are within each other's noise. So the honest reading is: on the
measured regime, arbitration does not beat the best single source, but it does
not have to lose to it either — and it covers 8,261 more reactions.

**What this does not show.** TECRDB is central metabolism, which is where all
three sources are at their best and where eQuilibrator's component-contribution
layer is anchored to measured data by construction. The 30,000-odd reactions
where the sources disagree wildly are not in it. The grading's BRONZE tier is aimed at
exactly those — its BRONZE tier is validated at 3.3 (eQuilibrator), 8.7 (Group
Contribution) and 20.8 (dGPredictor-ModelSEED) kcal/mol median error against 0.3
for GOLD — and that benefit is invisible in this table.

---

## 6. Recommendation

1. **For reaction direction specifically, prefer eQuilibrator where it has a
   non-sentinel value**, and fall back to the graded pick elsewhere. That is the
   best-scoring strategy tested (95.5%) and it retains full graded coverage.
   Concretely: ship `graded` with eQuilibrator promoted ahead of the p_ok
   tie-break within a grade tier.
2. **Add a direction-calibrated grade alongside the magnitude one.** Same
   machinery, different target; it closes most of the gap on its own and it is
   the correct calibration for a reversibility pipeline. Keep the magnitude grade
   for ΔG-consumers.
3. **Use `graded_trusted` when a wrong direction is costlier than a missing one.**
   It gives up 26 growers and 10,186 reactions to exclude the BRONZE tier.
4. **Do not use Group Contribution alone to set direction.** 51% on directional
   reactions is chance.
5. **Do not report growth counts as a quality ranking** without the permissiveness
   figure beside them.

---

## 7. Reproducing

```bash
cd /scratch/ctaylor/core_models_analysis/scripts
python3 grade_thermo_sources.py          # ~4 min  -> results/thermo_grades/
python3 build_graded_direction_maps.py   # ~2 min  -> results/thermo_grades_fba/
python3 run_graded_fba_all_models.py --workers 32   # ~4 min, 39,781 LP solves
python3 analyze_graded_fba.py            # ~3 min
python3 plot_graded_fba.py
```

`MSDB_ROOT` selects the data snapshot (default `devsnap2`), `MSDB_CODE` the
cascade checkout. 5,683 models × 7 variants, 0 errors.

---

## 8. Limitations

1. **The reference is the cascade on measured energies, not a measured
   direction.** If a heuristic is wrong, the reference inherits it. What the
   comparison isolates is the contribution of the *energy*, which is the question
   asked — but it cannot vindicate the cascade itself.
2. **802 reactions, all central metabolism.** Every accuracy number above is
   conditioned on the easy, well-measured part of the database. Group
   Contribution's 51% is a lower bound on how bad it gets, not an average.
3. **69 core reactions have measurements; the other 170 do not.** Core-set
   conclusions rest on those 65–69 scoreable reactions.
4. **Growth is a coarse readout.** A direction error only shows up if it happens
   to gate biomass. Reaction-level metrics (§4) are the sensitive ones; the
   per-reaction table `core_reaction_grades.tsv` is there for case-by-case work.
5. **Snapshot-specific.** Convention A changed 53% of Group Contribution values;
   these numbers move if the snapshot moves. `direction_maps_summary.json` stamps
   the commit.
