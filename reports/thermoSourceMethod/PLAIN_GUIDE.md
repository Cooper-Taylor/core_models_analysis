# A plain-language guide: the steps, and how we check the grades are meaningful

Companion to [`THERMO_SOURCE_METHOD.md`](THERMO_SOURCE_METHOD.md), which has the
mathematics, the full tables and the provenance. This document has neither. It
answers two questions:

- **[Part A](#part-a--the-steps)** — what does each step actually do?
- **[Part B](#part-b--how-we-evaluate-the-grading-criteria)** — how do we know
  the gold/silver/bronze labels mean anything?

---

# Part A — the steps

## The situation

A metabolic model needs to know which way each reaction can run. That comes from
ΔG′°, the reaction's free energy change. ModelSEED has up to four values for it
per reaction, from sources that disagree, and each source also states how
confident it is.

The naive plan — believe whichever source claims to be most confident — fails,
because the sources are wrong about themselves in different directions. Group
Contribution claims **4.4× more** error than it actually has; eQuilibrator
claims **1.6× less**. Comparing their confidence numbers directly just rewards
whichever source is most modest.

Everything below exists to fix that and then use the result.

---

## Step 1 — Collect the four sources

| | what it is | reactions covered |
|---|---|---:|
| **Group Contribution** | adds up energy contributions of chemical groups | 27,313 |
| **eQuilibrator** | component contribution, anchored on measured reactions | 25,028 |
| **dGPredictor-ModelSEED** | machine learning over molecular fragments | 31,924 |
| **TECRDB** | actual laboratory measurements | 1,550 |

Of 56,002 reactions, **33,337 have at least one source and 22,665 have none**.
That last number is a hard ceiling — nothing downstream improves it.

TECRDB is different in kind from the other three: it is data, not a prediction.
It gets used twice — once as the yardstick everything is measured against, and
once as a source in its own right.

---

## Step 2 — Throw out the values that are not really values

Some stored numbers are not noisy estimates; they are non-estimates wearing the
costume of one. These are removed outright rather than being modelled as
"uncertain":

| removed | count | why |
|---|---:|---|
| eQuilibrator sentinels | 4,934 | the source inflates its uncertainty by 10⁶ to say "I cannot estimate this". Real values stop at σ 65.35; sentinels start at 7,504.61. Nothing lives in between, so the cutoff is not a judgement call. |
| eQuilibrator ID collisions | 35 | a retrieval bug lets two compounds sharing an external ID overwrite each other |
| dGPredictor on quinones | 511 | the model regressed badly on this one chemistry — it disagrees with eQuilibrator on the *sign* 52.8% of the time |

**This step is the single highest-value use of the reported uncertainties in the
whole pipeline**, and it is almost free: a source explicitly declaring "no
estimate" is worth more than any inference we can make about it.

---

## Step 3 — Learn what each source's confidence actually means

**The problem:** "σ = 1" means something different for each source.

**What we do:** for each source separately, take the 802 reactions where we have
a lab measurement, plot its claimed σ against its actual error, and fit a curve.
The curve is forced to be *monotone* — a source claiming more uncertainty can
never be predicted to be more accurate — but otherwise no shape is imposed,
because there is no theory saying what shape it should be.

**Two outputs** from the same fit, differing only in the question asked:

- **ê** — "how far off is this number likely to be?" (in kcal/mol)
- **p** — "what is the chance this number is within 2 kcal/mol of the truth?"

2 kcal/mol is not an arbitrary threshold: it is the width of the band the
direction cascade itself treats as "reversible", so **p** reads directly as
*"is this number good enough to call the direction?"*

**One wrinkle worth knowing.** The 802 measured reactions are all well-studied
central metabolism, which means they are all *low-σ*. But three quarters of the
database sits at σ values far beyond anything the measurements cover. Fitting
only on measurements and extrapolating flat would mean being most optimistic
exactly where a source is least reliable — backwards for a safety check. So the
fit is padded with a second, weaker tier: reactions where a *trusted* source
(one inside its own well-measured σ range) can stand in as a stand-in reference.
Those points carry ⅓ the weight, because a stand-in gives an upper bound on the
error, not a measurement of it.

**What comes out:** eQuilibrator and dGPredictor get informative curves. Group
Contribution's is nearly flat — 0.845 down to 0.511 across the whole database.
That is not a broken fit; it is the finding that **GC's confidence barely
predicts its accuracy at all**. The consequence is structural: GC can never earn
a top grade from its own confidence, and it never does.

---

## Step 4 — Let the sources check each other

Now that the three confidences are on one scale, they can be combined.

Each source is weighted by how precise we now *know* it to be, and we ask: **is
the spread between these sources what their stated uncertainties predict?**

- Spread as expected → they corroborate each other.
- Spread much bigger → at least one is wrong, and the maths also says *which*
  one is the outlier.

This is a standard technique from particle physics (the "Birge ratio", used when
independent labs report conflicting measurements of the same constant).

**Two traps this step has to avoid, both real:**

1. **Free agreement.** 28.5% of the reactions where all sources agree are cases
   where every source says "about zero" — transport reactions and reactions
   whose chemistry cancels out. The agreement is imposed by the stoichiometry,
   not earned. These are flagged and excluded from counting as corroboration.
2. **Agreement bought with vagueness.** A source with a huge uncertainty is
   trivially "consistent" with everything. In one real example the two sources
   differ by 7 kcal/mol and still look consistent, because the third source's
   uncertainty is 9 kcal/mol wide.

---

## Step 5 — Label every number: gold, silver or bronze

Each source on each reaction gets its own label, independently. On the same
reaction, eQuilibrator can be gold and Group Contribution bronze — and that
happens: on `rxn00001`, Group Contribution has the wrong *sign*.

In order:

1. **If we have a lab measurement, just compare to it.** Within 1 kcal/mol →
   gold. Within 3 → silver. Worse → bronze. This overrides everything else.
2. **Otherwise, use the source's own calibrated confidence.** p ≥ 0.90 → gold.
   p ≥ 0.70 → silver. Else bronze.
3. **A bronze can be rescued to silver if other sources corroborate it** — but
   corroboration can *never* create a gold.
4. **A source outvoted by the others drops one level.**

**Why rule 3 is capped, in one sentence:** two fallible predictors agreeing is
weak evidence — they can share a blind spot — whereas a lab measurement is
strong evidence, so agreement should be able to raise a floor but not confer top
marks. This is not a preference; it was tested, and letting corroboration create
golds diluted what gold means (Part B, method 5).

---

## Step 6 — Pick one source to actually use

A label is not a choice. Asked "which of these three should I use?", the
answer depends on **what you are going to do with the number**:

- **If you need the energy itself** → use the source with the smallest expected
  error. The uncertainties work well for this.
- **If you need the direction** → use a fixed priority order: eQuilibrator,
  then dGPredictor, then Group Contribution.

The second one is deliberately *not* uncertainty-based, and that is the most
surprising result in this work. Every scheme we tried that used the
uncertainties to choose between sources performed **worse** than the flat
priority list. Details and the reason in Part B, method 8.

The uncertainties are still used here for two jobs: removing broken sources
(Step 2) and deciding when to give **no answer at all** rather than a coin flip.

---

## Step 7 — Turn an energy into a direction

The chosen ΔG′° is fed through ModelSEED's existing reversibility cascade,
completely unmodified. Six different direction maps were built this way, varying
*only* which energy goes in, so any difference downstream is attributable to the
thermodynamic source and nothing else.

---

## Step 8 — Run the metabolic models

All 5,683 core models, seven times each — once per variant plus the models'
own original settings. 39,781 optimisations, no failures.

**The headline is a caution, not a ranking.** Growth counts differ by only 4.5
percentage points across all seven variants, and they track how *permissive* a
variant is rather than how *correct*: eQuilibrator grows the fewest models and
is the most accurate. Anyone reporting "source X is better because more models
grow under it" has measured permissiveness.

The genuinely useful finding from this step is about the models themselves: the
reaction directions they ship with are the **least accurate** of everything
tested (67.7% correct against experiment, versus 90.8–98.5% for the
thermodynamic sources), and 19 of their 21 errors are the same mistake —
forcing a direction on a reaction that is actually reversible.

---

# Part B — how we evaluate the grading criteria

A grade is a claim: *"gold means you can trust this number."* That claim has to
be tested, and it can be tested in more than one way. Nine methods are used
below. Each says what it can show and — as importantly — what it cannot.

## The core difficulty: we mostly cannot see the truth

We have lab measurements for **802** of 56,002 reactions. Every direct
evaluation happens on those 802, and they are not a random sample — they are
well-studied central metabolism, i.e. the *easy* part of the database. So every
number in this section is conditioned on the easy half, and the tiers designed
to catch exotic chemistry are the least well tested.

That constraint shapes all nine methods.

---

## Method 1 — Hold out the answer, then grade, then check

**The primary test.** The grading cascade's first rule is "if we have a
measurement, use it" — so on the 802 measured reactions, the grade is trivially
correct and proves nothing.

So we **switch that rule off**, regrade those 802 reactions using only the
source's calibrated confidence and cross-source agreement, and *then* compare
against the measurement the grader was not allowed to see.

**What it should show if the grades are real:** error rising monotonically from
gold to silver to bronze, for every source.

| source | grade | n | median error | within 2 kcal/mol |
|---|---|---:|---:|---:|
| eQuilibrator | GOLD | 246 | **0.32** | 94% |
| | SILVER | 529 | 0.46 | 85% |
| | BRONZE | 14 | **3.33** | **0%** |
| dGPredictor-ModelSEED | GOLD | 184 | **0.32** | 98% |
| | SILVER | 608 | 0.55 | 82% |
| | BRONZE | 10 | **20.78** | **0%** |
| Group Contribution | GOLD | 0 | — | — |
| | SILVER | 517 | 1.28 | 69% |
| | BRONZE | 285 | **8.68** | 34% |

Monotone in every column for every source. Gold-to-bronze separates by 10× for
eQuilibrator and 65× for dGPredictor.

**What it cannot show:** whether bronze works on the exotic chemistry it exists
for. Those bronze rows are n = 10 to 285, all central metabolism.

---

## Method 2 — Require monotonicity, not just separation

A grade that put gold and bronze far apart but shuffled silver would be a bad
grade. So the test above is read column by column — median, mean, within-1,
within-2, p90 — and each must be ordered.

It passes, with **one honest exception** worth recording rather than hiding:
Group Contribution's bronze tier has a *higher* within-1-kcal/mol rate (33%)
than its median suggests. That tier is bimodal — a third of it is nearly exact
and the rest is badly wrong. The lesson is to read the median and the p90 rather
than the mean for that source.

---

## Method 3 — Guard against circularity explicitly

Any evaluation that uses the same measurements the grader used is worthless. Two
guards:

1. **Rule 1 is disabled** during evaluation (Method 1).
2. For the downstream direction test, a **separate variant is built that has
   never seen TECRDB at all** — no measurement override, and TECRDB removed as
   a source. Without it the graded map scores a meaningless 100%, because it
   contains the answer.

That second guard matters: the plain graded map scores **802/802** against
experiment. The held-out version scores **93.4%**. Reporting the first number
would be nonsense, and the trap is easy to fall into.

---

## Method 4 — Turn each rule off and measure what changes

A cascade with four rules invites the assumption that all four matter. Turning
them on one at a time says otherwise. Bronze-tier median error:

| source | confidence only | + corroboration lift | + demotion (shipped) |
|---|---|---|---|
| Group Contribution | 1.66 (n=779) | **8.68** (n=285) | 8.68 (n=285) |
| eQuilibrator | 2.26 (n=21) | 3.47 (n=11) | 3.33 (n=14) |
| dGPredictor-ModelSEED | 2.57 (n=5) | 16.23 (n=2) | **20.78** (n=10) |

**The lift does nearly all the work.** By moving corroborated reactions *out* of
bronze it concentrates the genuinely bad ones there.

**The demotion is mostly decorative.** It fires on 4,137 rows but changes a
grade on only 298 of them — **7%** — because the rest were already bronze. Its
one real contribution is to dGPredictor, adding 8 bad reactions to bronze.

*(An earlier version of this analysis credited the demotion with the
1.66 → 8.68 improvement. That was wrong; it is the lift. The ablation is what
caught it.)*

---

## Method 5 — Test the alternative you rejected

The design says corroboration may raise a bronze to silver but never create a
gold. That is a choice, so the other choice was implemented and measured.

Letting corroboration promote all the way to gold:

| | golds | measured guarantee (within 2 kcal/mol) |
|---|---:|---:|
| eQuilibrator, capped (shipped) | 2,443 | **94%** |
| eQuilibrator, uncapped | 9,157 | 90% |
| dGPredictor, capped (shipped) | 5,808 | **98%** |
| dGPredictor, uncapped | — | 91% |

Nearly 4× more golds, and gold means measurably less. The cap stays.

---

## Method 6 — Audit which rule actually fired

Every graded row records *why* it got its grade. Counting those reasons tests
whether the design behaves as described:

| reason | GOLD | SILVER | BRONZE |
|---|---:|---:|---:|
| own confidence ≥ 0.90 | **7,119** | — | — |
| measured against lab data | **1,421** | 574 | 398 |
| own confidence ≥ 0.70 | — | 11,301 | — |
| corroborated by others | — | **26,984** | — |
| outvoted by others | — | — | 4,086 |
| nothing supported it | — | — | **26,902** |

Three checks this passes or reveals:

- **Corroboration produces zero golds** — the asymmetry is real in the output,
  not just claimed in the design.
- **Group Contribution has no "own confidence ≥ 0.90" rows anywhere.** All 309
  of its golds come from lab measurements. This is Step 3's flat curve showing
  up as a hard structural fact.
- **Corroboration is the most-used rule in the system** (26,984 silvers, more
  than double the next). Worth flagging as a caution rather than a success:
  the middle tier rests mostly on the *weaker* kind of evidence.

---

## Method 7 — Check the calibration is honest, not just fitted

The confidence curves make a probabilistic claim, so that claim is checked
directly: if we say a value is within ±τ, is it, at the stated rate?

| source | coverage before correction | after | target |
|---|---:|---:|---:|
| Group Contribution | 0.762 | 0.683 | 0.683 |
| eQuilibrator | 0.866 | 0.694 | 0.683 |
| dGPredictor-ModelSEED | 0.813 | 0.687 | 0.683 |

All three were over-cautious and needed shrinking by about a third. This is
printed on every run, so the calibration is checkable rather than asserted.

---

## Method 8 — Test the grade at the job it was *not* designed for

The grade sorts numbers by trustworthiness. Does that make it a good way to
*choose* a source? Tested directly, on held-out data, against every alternative
including trivial ones:

| strategy | direction accuracy |
|---|---:|
| fixed priority: eQuilibrator → dGPredictor → Group Contribution | **95.9%** |
| always eQuilibrator | 95.9% |
| pick by smallest calibrated uncertainty | 93.7% |
| pick by smallest expected error | 93.7% |
| always dGPredictor | 91.6% |
| **pick by lowest direction risk** | **90.9%** |
| always Group Contribution | 85.1% |

**Every uncertainty-based choice lost to the flat priority list**, and the most
sophisticated one came second-to-last.

**Why**, and it is worth understanding rather than memorising: the risk measures
*"how likely is this source's answer to survive its own uncertainty"* — which is
**precision, not accuracy**. A source that is confidently wrong scores
beautifully. The calculation is centred on the source's own estimate, so it
cannot see that estimate being displaced from the truth.

This is the strongest evidence for the document's central claim: **the reported
uncertainties are usable within a source and not between sources.**

---

## Method 9 — Check it downstream, where it is actually used

Grades that look good in a table but change nothing in practice would not be
worth shipping. Two end-use checks:

- **Direction accuracy.** The held-out graded map gets 93.4% of reaction
  directions right, against eQuilibrator alone at 95.5% — so on the measured
  set the grading does not beat the best single source. Its advantage is
  **coverage**: 33,289 reactions versus 25,028, each carrying a label.
- **Applying a quality floor.** Dropping every bronze-graded reaction removes
  10,186 reactions and changes the growth verdict on 26 of 5,683 models. Small,
  but it is the intended knob: use it when a wrong direction costs more than a
  missing one.

---

## What would change our minds

Stated in advance, so the grading is falsifiable rather than merely defended:

| observation | what it would mean |
|---|---|
| Error not monotone across tiers on a fresh measurement set | the grade does not measure trustworthiness |
| Gold's within-2 rate falling well below ~90% out of sample | the gold threshold is too loose |
| Bronze performing as well as silver on exotic chemistry | the tiers only separate on easy reactions |
| A source's confidence curve inverting on new data | the calibration does not transfer |

## The three limitations that bound all nine methods

1. **802 reactions, all central metabolism.** Everything here is conditioned on
   the easy part of the database.
2. **The benchmark is partly in-sample.** eQuilibrator is *fitted* on these
   measurements and dGPredictor was trained on 4,001 of them; only Group
   Contribution is genuinely out-of-sample. So the accuracy ordering between
   sources is not fully independent evidence — the fair reading is that
   "prefer eQuilibrator" is the best rule available, not a settled fact.
3. **The direction reference is the cascade run on measured energies**, not a
   curated direction. It isolates the contribution of the energy, which is the
   question being asked, but it cannot validate the cascade itself.
