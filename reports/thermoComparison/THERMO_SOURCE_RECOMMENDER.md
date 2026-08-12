# Recommending a thermodynamic source per reaction

*What the gold/silver/bronze grade is for, what it is not for, and the decision
rule that replaces it for source selection. 2026-08-12.*

Script: [`recommend_thermo_source.py`](/scratch/ctaylor/core_models_analysis/scripts/recommend_thermo_source.py).
Data: `results/thermo_recommendation/`.
Related: [`THERMO_SOURCE_GRADING_PROPOSAL.md`](THERMO_SOURCE_GRADING_PROPOSAL.md) (the grade),
[`THERMO_SOURCE_ASSIGNMENT.md`](THERMO_SOURCE_ASSIGNMENT.md) (the magnitude rule),
[`GRADED_SOURCE_CORE_MODEL_ANALYSIS.md`](GRADED_SOURCE_CORE_MODEL_ANALYSIS.md) (where the problem surfaced).

---

## 1. What the grade is, in one page

For each reaction, each of the four sources gets an independent label.

**TECRDB** is always GOLD — it is a measurement. (Skeleton-tier matches, which
are blind to stereochemistry, are capped at SILVER: the measurement is gold, the
match to a ModelSEED reaction is not.)

The three predictors are labelled by a cascade:

1. **Measured** — where TECRDB exists, compare directly. \|error\| ≤ 1 → GOLD,
   ≤ 3 → SILVER, above → BRONZE. Per source, so on one reaction eQuilibrator can
   be GOLD and Group Contribution BRONZE.
2. **Self-confidence** — otherwise, `p_ok = P(|ΔG_s − ΔG*| ≤ 2 kcal/mol | σ_s)`,
   fitted per source by isotonic regression against measured error. ≥ 0.90 GOLD,
   ≥ 0.70 SILVER, else BRONZE. This is what converts three incomparable σ scales
   into one comparable number.
3. **Corroboration, asymmetrically** — agreement with other sources can lift a
   BRONZE to SILVER but never make a GOLD; being the outlier in a discrepant set
   costs one tier.
4. **Vetoes** → UNGRADED: eQuilibrator sentinels, MetaNetX collisions,
   dGPredictor-ModelSEED on quinones.

**What it is good for, demonstrated:** the tiers separate real error cleanly.
Held out from the measurement, median \|error\| runs GOLD 0.32 → SILVER 0.46 →
BRONZE 3.33 kcal/mol for eQuilibrator, and 0.32 → 0.55 → **20.78** for
dGPredictor-ModelSEED. As a *trust label* on a number you already have, it works.

**What it is not good for:** choosing between sources. That was never validated,
and when tested it failed — §2.

---

## 2. Why the grade is the wrong tool for selection

The grade's target is ΔG′° **magnitude**. On that target dGPredictor-ModelSEED
is genuinely the better source (98% of its GOLD tier within 2 kcal/mol against
eQuilibrator's 94%). So a grade-ranked pick chooses dGPredictor-ModelSEED on 521
of the 802 measured reactions — and gets the **direction** wrong more often than
just always using eQuilibrator, because eQuilibrator is right 98.9% of the times
it is picked and dGPredictor-ModelSEED only 90.4%.

The reason is structural, not a bug: **direction errors concentrate where
magnitude error is smallest**, at ΔG′° ≈ 0, right where the cascade's ±2 kcal/mol
band decides. Optimising magnitude does not optimise direction.

So the recommender takes the target as an input.

---

## 3. The algorithm

```
RECOMMEND(reaction i, target T ∈ {direction, magnitude}):

  0.  MEASUREMENT      if TECRDB has i, return the experimental ΔG′°.

  1.  FEASIBILITY      drop any source vetoed by its own uncertainty:
                         eQuilibrator σ > 100          (4,934 — the source
                                                        disclaiming the reaction)
                         eQuilibrator MetaNetX collision (35)
                         dGPredictor-ModelSEED on a quinone (511)

  2.  CALIBRATE        τ_s(i) = k_s · ê_s(i) / √(2/π)
                       ê from the two-tier isotonic fit; k_s one scalar per
                       source, fitted so ±τ covers 68.3% of measured error.
                       Fitted values: GC 0.539, eQ 0.598, dGPMS 0.648;
                       coverage 0.762 → 0.683, 0.866 → 0.694, 0.813 → 0.687.

  3.  SELECT           T = magnitude:  s* = argmin_s ê_s(i)
                       T = direction:  s* = first feasible in EQ > DGPMS > GC

  4.  RISK             compute the chosen source's direction risk (§4).

  5.  ABSTAIN          return nothing if risk > tolerance
                       (default 0.35 for direction, ê ≤ 2 kcal/mol for magnitude).
```

Step 3 is the part that surprised me, and §5 is the evidence for it.

---

## 4. The direction risk — exact, not sampled

The uncertainty still needs to be turned into something meaningful about
*direction*. The trick is that the cascade's operator is a **piecewise-constant
function of ΔG′°**. Holding stoichiometry fixed, it can only change at:

| breakpoint | which heuristic |
|---|---|
| −σ − RT(pdt_max + rct_min) | stored-bounds max crosses 0 |
| +σ − RT(pdt_min + rct_max) | stored-bounds min crosses 0 |
| −2 − RT·rgt_sum, +2 − RT·rgt_sum | the mMΔG ±2 band edges |
| −RT·rgt_sum | mMΔG sign flip |
| 2/points − RT·rgt_sum | low-energy-points threshold |

Evaluate the **real cascade** once inside each interval, then integrate
N(ΔG_s, τ_s²) over the intervals whose operator matches the source's own call.
That is a closed-form sum of normal CDFs — no Monte Carlo, and exact with
respect to the cascade rather than an approximation of one heuristic. ATP
synthase and ABC transporters short-circuit before any ΔG is read, so their risk
is exactly 0.

Resulting risk: eQuilibrator median 0.000 (76.6% at ≤ 0.05), Group Contribution
0.048 (50.3%), dGPredictor-ModelSEED 0.081 (47.4%).

---

## 5. Validation, including the negative result

Held-out direction accuracy on the 802 TECRDB stereo-exact anchors, 20 × 70/30
splits. Coverage is the share of held-out reactions the strategy will answer for.

| strategy | accuracy | coverage |
|---|---:|---:|
| **priority EQ > DGPMS > GC** | **95.9% ± 1.1** | **100%** |
| eQuilibrator only | 95.9% ± 1.2 | 98.3% |
| priority + risk veto at 0.20 | 94.2% ± 1.1 | 100% |
| priority + risk veto at 0.05 | 93.8% ± 1.4 | 100% |
| argmin calibrated τ | 93.7% ± 1.1 | 100% |
| argmin ê (the magnitude rule) | 93.7% ± 1.1 | 100% |
| priority + risk veto at 0.02 | 93.0% ± 1.3 | 100% |
| dGPredictor-ModelSEED only | 91.6% ± 1.1 | 100% |
| **argmin direction risk** | **90.9% ± 1.2** | 100% |
| Group Contribution only | 85.1% ± 2.3 | 100% |

**Every uncertainty-based arbitration lost to a fixed priority order**, and
layering a risk veto on top of priority made it monotonically worse — the veto
only ever pushes a reaction off the better source onto a worse one.

**Why argmin-risk fails, and it is worth being precise about this.** The risk is
P(this source's own call is overturned by this source's own uncertainty). That is
**precision, not accuracy**. A source that is confidently wrong — small τ, far
from any breakpoint, wrong region — scores risk ≈ 0 and beats a source that is
right but sits near a band edge. An integral centred on a biased point estimate
cannot see the bias. The quantity is well-defined and correctly computed; it is
just not the quantity that selects a source.

---

## 6. What the uncertainty is still doing

Three jobs, all validated, none of them arbitration between sources:

1. **Feasibility.** The sentinel is the single most valuable uncertainty signal
   in the database — 4,934 reactions where eQuilibrator explicitly says it has no
   estimate. Reading it costs nothing and prevents a nonsense value.
2. **Abstention.** Within a source, uncertainty is informative. eQuilibrator's
   direction accuracy by its own σ quartile: **100.0% / 98.6% / 91.9% / 92.0%**.
   So "should I answer at all" is a question σ can answer even though "which
   source" is not.
3. **Magnitude arbitration.** For ΔG′° itself, argmin ê is the right rule and
   beats always-eQuilibrator on mean error (1.03 vs 1.75 kcal/mol).

That is the clean summary of this whole line of work: **the reported
uncertainties are usable within a source and not between sources.**

---

## 7. Output

| target | reactions recommended | mix |
|---|---:|---|
| direction (risk ≤ 0.35) | 27,240 | eQuilibrator 17,574 · dGPredictor-MS 7,259 · TECRDB 1,550 · Group Contribution 857 |
| magnitude (ê ≤ 2) | 11,841 | dGPredictor-MS 5,499 · eQuilibrator 4,792 · TECRDB 1,550 |

6,049 reactions have a feasible source but are abstained on (risk > 0.35);
22,713 have no feasible source at all. Columns: per-source ΔG and risk, the
choice, the recommended ΔG′° and σ, and `kept`.

```bash
python3 scripts/recommend_thermo_source.py                    # both targets
python3 scripts/recommend_thermo_source.py --target direction --tolerance 0.2
```

---

## 8. The caveat that limits all of this

**eQuilibrator is fitted on TECRDB.** Its reactant-contribution layer is anchored
to exactly these measurements, and dGPredictor was trained on 4,001 of them; only
Group Contribution is genuinely out-of-sample. So the priority order is derived
from a benchmark that two of the three contestants have partly seen.

The stratification is consistent with that: eQuilibrator's advantage is
concentrated where its own σ is smallest — 100% in its lowest σ quartile, and by
the third quartile it is level with dGPredictor-ModelSEED (91.9% vs 92.3%, n=209,
a one-reaction difference). That is what partial memorisation looks like.

So treat "prefer eQuilibrator" as the best rule available on the evidence we
have, not as a settled fact about the methods. The clean test would be a held-out
measurement set neither model was fitted on. Two other standing limits also
apply: 802 reactions, all central metabolism; and the reference direction is the
cascade run on measured energies, so it inherits any error in the cascade itself.
