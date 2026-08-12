# 4. The recommendation algorithm — choosing which source to use

Notation: [`02_notation.md`](02_notation.md). Script:
`scripts/recommend_thermo_source.py`. Output: `results/thermo_recommendation/`.

---

## 4.1 The decision problem

For each reaction *i*, choose one source *s*\*(*i*) ∈ F(*i*), or abstain. The
general form is a minimum-risk decision:

> *s*\*(*i*) = argmin over *s* ∈ F(*i*) of **𝓛_T( *s*, *i* )**
>
> use *s*\*(*i*) if 𝓛_T( *s*\*, *i* ) ≤ *E*\*, otherwise return nothing

where 𝓛_T is the expected loss under the **target** *T* — what the caller is
going to do with the number. Two targets are supported because they turn out to
disagree:

| target *T* | loss | 𝓛_T( *s*, *i* ) |
|---|---|---|
| **magnitude** | \|ΔG*ₛ* − ΔG\*\| | ê*ₛ*(*i*) — the calibrated expected absolute error |
| **direction** | 𝟙[ Λ*ₛ* ≠ Λ\* ] | ρ*ₛ*(*i*) = 1 − P( Λ\*(*i*) = Λ*ₛ*(*i*) ) |

For the magnitude target the argmin rule is correct and is what ships. For the
direction target it is not, and §4.5–4.6 replace it.

## 4.2 Why the grade cannot serve as the selector

The grade of §3 is calibrated on **magnitude** — *p*ₛ = P(\|ΔG*ₛ* − ΔG\*\| ≤ 2).
On that target dGPredictor-ModelSEED is genuinely the better source: 98% of its
GOLD tier lands within 2 kcal/mol against eQuilibrator's 94%.

So a grade-ranked pick chooses dGPredictor-ModelSEED on **521 of the 802** anchor
reactions and eQuilibrator on 274. But eQuilibrator's *direction* call is right
**98.9%** of the times it is picked, and dGPredictor-ModelSEED's only **90.4%**.
The ranking systematically prefers the weaker source for this job.

The cause is structural. **Direction errors concentrate where magnitude error is
smallest**: near ΔG′° ≈ 0, which is exactly where the cascade's ±2 kcal/mol band
edge sits. Optimising E\|error\| does not optimise P(right side of the band).
A single statistic cannot serve both targets, which is why *T* is an input.

## 4.3 Calibrating σ to a usable standard deviation

The direction risk needs a probability distribution over the truth, not a point
estimate of error. Starting from ê and the Gaussian identity E\|X\| = τ√(2/π):

> **τ*ₛ*(*i*) = *k*ₛ · ê*ₛ*(*i*) / √(2/π)**,   floored at 0.05 kcal/mol

*k*ₛ is one scalar per source, fitted so that the interval ΔG*ₛ* ± τ*ₛ* actually
covers the nominal 68.3% of measured errors on the anchor. Without it the
intervals are too wide, because ê is fitted partly on the proxy tier, which is an
upper bound rather than a measurement.

| source | *k*ₛ | anchor n | coverage before | coverage after | target |
|---|---:|---:|---:|---:|---:|
| Group Contribution | 0.539 | 802 | 0.762 | 0.683 | 0.683 |
| eQuilibrator | 0.598 | 794 | 0.866 | 0.694 | 0.683 |
| dGPredictor-ModelSEED | 0.648 | 802 | 0.813 | 0.687 | 0.683 |

All three needed shrinking by roughly a third — i.e. ê is conservative, as §5b of
the assignment report found independently. The before/after coverage is printed
on every run and stored in `tables/recommendation_models.json`, so the
calibration is checkable rather than asserted.

## 4.4 The direction risk, computed exactly

**The key observation:** holding stoichiometry fixed, the cascade's operator is
a **piecewise-constant function of ΔG′°**. Every heuristic that reads the energy
compares it against a threshold, so 𝒞(·, e; i) can only change value at finitely
many points, all in closed form:

| breakpoint in ΔG′° | which comparison it comes from |
|---|---|
| − *e* − *a*(*i*) | stored-bounds maximum, ΔG + *e* + *a* = 0 |
| + *e* − *b*(*i*) | stored-bounds minimum, ΔG − *e* + *b* = 0 |
| − 2 − *c*(*i*) | lower edge of the mMΔG band, mMΔG = −2 |
| + 2 − *c*(*i*) | upper edge of the mMΔG band, mMΔG = +2 |
| − *c*(*i*) | mMΔG sign flip, which selects the low-energy branch |
| 2 / *P*(*i*) − *c*(*i*) | the low-energy threshold, *P*·mMΔG = 2 (omitted when *P* = 0) |

Sort these into an ordered partition of the real line,

> −∞ = *t*₀ < *t*₁ < … < *t*_m < *t*_{m+1} = +∞,   *I*_j = (*t*_j, *t*_{j+1})

and evaluate the **real cascade** once at an interior point of each interval,
ω_j = 𝒞( midpoint(*I*_j), σ\* ; *i* ). Because the function is constant on each
interval, ω_j is the operator on the whole of *I*_j — this is exact, not a
discretisation. Evaluating the actual cascade rather than a re-derivation also
means the calculation cannot drift from the cascade's behaviour if a heuristic
changes.

Then, modelling the truth as ΔG\* ~ N( ΔG*ₛ*(*i*), τ*ₛ*(*i*)² ):

> **P( Λ\*(*i*) = Λ*ₛ*(*i*) ) = Σ_j : ω_j = Λ*ₛ*(*i*) [ Φ( (t_{j+1} − ΔG*ₛ*) / τ*ₛ* ) − Φ( (t_j − ΔG*ₛ*) / τ*ₛ* ) ]**
>
> **ρ*ₛ*(*i*) = 1 − P( Λ\*(*i*) = Λ*ₛ*(*i*) )**

with Φ the standard normal CDF. A closed-form sum of normal CDFs — no Monte
Carlo, no sampling error, and exact with respect to the full cascade rather than
an approximation of one heuristic.

The inner *e* used when evaluating the reference operator is σ\* = 0.15 kcal/mol,
the median TECRDB experimental standard deviation, since Λ\* is defined as the
cascade fed a *measurement*.

**Reactions the cascade decides without reading ΔG′° at all** — ATP synthase and
ABC transporters, which match at heuristics 1 and 2 — yield a single interval,
and ρ*ₛ* = 0 exactly for every source.

Resulting distribution:

| source | n | median ρ | share with ρ ≤ 0.05 |
|---|---:|---:|---:|
| eQuilibrator | 20,059 | 0.000 | 76.6% |
| Group Contribution | 27,313 | 0.048 | 50.3% |
| dGPredictor-ModelSEED | 31,413 | 0.081 | 47.4% |

## 4.5 Validation, and the negative result

Held-out direction accuracy on the anchor set, 20 random 70/30 splits. Coverage
is the share of held-out reactions the strategy will answer for. The only fitted
parameters the split protects are the three *k*ₛ; the risk model has no
per-reaction free parameters.

| strategy | accuracy | coverage |
|---|---:|---:|
| **priority EQ > DG > GC** | **95.9% ± 1.1** | **100%** |
| eQuilibrator only | 95.9% ± 1.2 | 98.3% |
| priority + risk veto at ρ > 0.20 | 94.2% ± 1.1 | 100% |
| priority + risk veto at ρ > 0.05 | 93.8% ± 1.4 | 100% |
| argmin τ*ₛ* | 93.7% ± 1.1 | 100% |
| argmin ê*ₛ* (the magnitude rule) | 93.7% ± 1.1 | 100% |
| priority + risk veto at ρ > 0.02 | 93.0% ± 1.3 | 100% |
| dGPredictor-ModelSEED only | 91.6% ± 1.1 | 100% |
| **argmin ρ*ₛ* (the risk rule)** | **90.9% ± 1.2** | 100% |
| Group Contribution only | 85.1% ± 2.3 | 100% |

Two things to read off this table.

**Every uncertainty-based arbitration lost to a fixed priority order.** argmin ρ
came second-to-last, below simply always using dGPredictor-ModelSEED.

**Layering a risk veto on top of priority made it monotonically worse** —
94.2 → 93.8 → 93.0 as the veto tightens. A veto can only ever move a reaction
*off* the better source and *onto* a worse one, so tightening it trades accuracy
away.

### Why argmin ρ fails

ρ*ₛ* is P(this source's own call is overturned by this source's own
uncertainty). That is **precision, not accuracy**. The integral in §4.4 is
centred on ΔG*ₛ*, so it measures how far that point estimate sits from a
breakpoint relative to its own noise — and it is structurally blind to the point
estimate being displaced from the truth in the first place.

Concretely: a source that is confidently wrong — small τ*ₛ*, far from any
breakpoint, wrong region — scores ρ ≈ 0 and wins the argmin against a source
that is right but happens to sit near a band edge. The quantity is well-defined
and correctly computed. It is simply not the quantity that selects a source.

## 4.6 The algorithm as shipped

```
RECOMMEND( reaction i, target T ):

  Step 0  MEASUREMENT
          if i has a TECRDB match, return (ΔG_TEC, σ_TEC, "TECRDB").

  Step 1  FEASIBILITY                                          [uses uncertainty]
          F(i) = A(i) \ V(i)
          V removes: eQuilibrator σ > 100 (sentinel, 4,934)
                     eQuilibrator MetaNetX collision (35)
                     dGPredictor-ModelSEED on a quinone (511)
          if F(i) = ∅ : return nothing.

  Step 2  CALIBRATE                                            [uses uncertainty]
          τ_s(i) = k_s · ê_s(i) / √(2/π)   for each s ∈ F(i)

  Step 3  SELECT
          T = magnitude :  s* = argmin_{s ∈ F(i)} ê_s(i)       [uses uncertainty]
          T = direction :  s* = first source in (EQ, DG, GC) that lies in F(i)
                                                               [does NOT]

  Step 4  RISK
          compute ρ_{s*}(i) by the interval integral of §4.4.

  Step 5  ABSTAIN                                              [uses uncertainty]
          return nothing if the risk exceeds the tolerance:
              direction : ρ_{s*}(i) > 0.35
              magnitude : ê_{s*}(i) > 2.0 kcal/mol
          otherwise return (ΔG_{s*}, σ_{s*}, s*, risk).
```

**The priority order is empirical, not a preference.** It is the measured
accuracy ordering on the anchor: eQuilibrator 95.5% > dGPredictor-ModelSEED
91.8% > Group Contribution 85.5% for direction. §4.7 is the caveat on it.

## 4.7 What the uncertainty is and is not for

The clean summary of this whole line of work:

> **The reported uncertainties are usable *within* a source and not *between*
> sources.**

Three jobs they do, all validated:

1. **Feasibility.** The eQuilibrator sentinel is the single most valuable
   uncertainty signal in the database — 4,934 reactions where the source states
   outright that it has no estimate. Reading it costs nothing and prevents a
   meaningless number entering a model.
2. **Abstention.** Within a source, σ is informative. eQuilibrator's direction
   accuracy by its own σ quartile:

   | quartile | σ range | n | eQuilibrator | dGPredictor-MS | Group Contribution |
   |---|---|---:|---:|---:|---:|
   | Q1 | 0.00–0.17 | 205 | **100.0%** | 98.5% | 83.9% |
   | Q2 | 0.17–0.36 | 208 | 98.6% | 91.8% | 88.5% |
   | Q3 | 0.36–0.57 | 209 | 91.9% | 92.3% | 87.1% |
   | Q4 | 0.57–1.12 | 200 | 92.0% | 84.5% | 83.5% |

   So "should I answer at all" is a question σ can answer even though "which
   source" is not.
3. **Magnitude arbitration.** For ΔG′° itself, argmin ê is right and beats
   always-eQuilibrator on mean error, 1.03 vs 1.75 kcal/mol.

### The caveat that limits the priority order

**eQuilibrator is fitted on TECRDB.** Its reactant-contribution layer is anchored
to exactly these measurements, and dGPredictor was trained on 4,001 of them. Only
Group Contribution is genuinely out-of-sample. The benchmark that produced the
priority order is therefore partly in-sample for two of the three contestants.

The quartile table above is consistent with that reading: eQuilibrator's
advantage is concentrated where its own σ is smallest — 100.0% in Q1 — and by Q3
it is level with dGPredictor-ModelSEED (91.9% vs 92.3%, n = 209, a
one-reaction difference). That is what partial memorisation would look like.

Treat "prefer eQuilibrator" as the best rule available on the evidence in hand,
not as a settled fact about the methods. The clean test is a held-out measurement
set neither model was fitted on, and it has not been run.

## 4.8 Results

### Which source gets recommended, across all covered reactions

The universe is the **33,289 reactions with at least one feasible source** —
i.e. everything the recommender could in principle answer for. The remaining
22,713 of the 56,002 non-EMPTY reactions have no thermodynamic source at all and
are outside the question. Shares are of 33,289 and sum to 100%.

**Target = direction** (priority rule, abstain at ρ > 0.35):

| source recommended | reactions | share of 33,289 | share of the 27,240 answered |
|---|---:|---:|---:|
| eQuilibrator | 17,574 | **52.79%** | 64.52% |
| dGPredictor-ModelSEED | 7,259 | **21.81%** | 26.65% |
| TECRDB | 1,550 | **4.66%** | 5.69% |
| Group contribution | 857 | **2.57%** | 3.15% |
| *(abstained — risk above tolerance)* | 6,049 | **18.17%** | — |
| **total** | **33,289** | **100.00%** | **100.00%** |

**Target = magnitude** (argmin ê, abstain at ê > 2.0 kcal/mol):

| source recommended | reactions | share of 33,289 | share of the 11,841 answered |
|---|---:|---:|---:|
| dGPredictor-ModelSEED | 5,499 | **16.52%** | 46.44% |
| eQuilibrator | 4,792 | **14.40%** | 40.47% |
| TECRDB | 1,550 | **4.66%** | 13.09% |
| Group contribution | 0 | **0.00%** | 0.00% |
| *(abstained — ê above tolerance)* | 21,448 | **64.43%** | — |
| **total** | **33,289** | **100.00%** | **100.00%** |

Four things this table shows.

**The two targets disagree about the winner.** eQuilibrator takes 64.5% of the
direction answers; dGPredictor-ModelSEED takes 46.4% of the magnitude answers
and edges eQuilibrator there. That inversion is §4.2 made operational — it is
the same dissociation, counted.

**Group Contribution is recommended for almost nothing**: 2.57% under direction,
and *never* under magnitude, because its calibrated ê never falls below 2
kcal/mol anywhere in the database (its floor is 3.04). Under direction it wins
only the 857 reactions where it is the sole feasible source. Given its 51%
accuracy on directional reactions (§5.6 of the simulations write-up), that is
the right outcome — but note it is reached by the priority order, not by any
uncertainty-based judgement.

**The magnitude target abstains on 64% of what it could answer.** That is not a
failure; ê ≤ 2 kcal/mol is a demanding bar and the tolerance is a knob. Loosening
it trades coverage for stated accuracy along the frontier in
`THERMO_SOURCE_ASSIGNMENT.md` §7.

**The direction target abstains on 18%.** These are the 6,049 reactions where
even the best feasible source has ρ > 0.35 — the source's own calibrated
uncertainty is wide enough, relative to the nearest cascade breakpoint, that its
call is close to a coin flip. This is the abstention job of §4.7 doing work.

For context, of the 33,289 covered reactions, 17,389 have all three predictors
feasible, 10,718 have two, and 5,182 have exactly one.

### Files

| target | tolerance | answered | file |
|---|---|---:|---|
| direction | ρ ≤ 0.35 | 27,240 | `recommendation_direction.tsv` |
| magnitude | ê ≤ 2.0 | 11,841 | `recommendation_magnitude.tsv` |

Columns in `recommendation_<target>.tsv`: `rxn, name, ec, status`, per-source
`dg_*` and `risk_*`, then `chosen_source, chosen_label, recommended_dg,
recommended_sigma, risk, n_feasible, kept`.

The two targets disagree about which source to use on a large fraction of the
reactions they both cover — that disagreement *is* the §4.2 finding, made
operational.
