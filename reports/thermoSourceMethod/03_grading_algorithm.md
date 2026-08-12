# 3. The grading algorithm — assigning gold / silver / bronze per source

Notation: [`02_notation.md`](02_notation.md). Script:
`scripts/grade_thermo_sources.py`. Output: `results/thermo_grades/`.

**What this produces.** For each reaction *i* and each source *s* ∈ 𝓢⁺, a label
G*ₛ*(*i*) ∈ {GOLD, SILVER, BRONZE, UNGRADED} answering *how much should I trust
this particular number*. It does **not** answer *which number should I use* —
that is [`04`](04_recommendation_algorithm.md), and §04.2 shows why the two must
be separated.

---

## 3.1 Why σ cannot be used raw

If σ*ₛ* were an honest Gaussian standard deviation, no calibration would be
needed: for X ~ N(0, σ²), E\|X\| = σ√(2/π) ≈ 0.798 σ, so the expected error
would follow directly from the reported number. Measured against the anchor set,
none of the three obeys this:

| source | median σ | median ε | ratio ε/σ | vs the Gaussian 0.798 |
|---|---:|---:|---:|---|
| Group Contribution | 4.35 | 1.57 | 0.368 | σ **overstates** error 2.2× |
| dGPredictor-ModelSEED | 0.91 | 0.47 | 0.522 | σ **overstates** error 1.5× |
| eQuilibrator | 0.36 | 0.45 | 1.260 | σ **understates** error 1.6× |

Two sources are pessimistic and one is optimistic, by different factors. A rule
that compares raw σ across sources therefore rewards whichever source
self-reports most pessimistically, independent of accuracy. *(This is not
hypothetical: it is the mechanism by which dev's
`Promote_Reaction_Thermodynamics_to_Canonical.py`, which breaks ties on
min \|reported error\|, prefers the KEGG-mis-mapped legacy dGPredictor over its
correctly-keyed replacement on 95.3% of the reactions carrying both.)*

Calibrating **per source** puts all three on the single scale of *expected error
against truth*.

## 3.2 Calibration

### The estimator

Both calibrations are **isotonic regressions** — monotone and non-parametric.
Monotone because a source reporting more uncertainty must not be predicted to be
*more* accurate; non-parametric because there is no theory fixing the functional
form.

> ĝ*ₛ* = argmin over non-**decreasing** *f* of Σ_j w_j ( ε_j − *f*(σ_j) )²
>
> ĥ*ₛ* = argmin over non-**increasing** *f* of Σ_j w_j ( 𝟙[ε_j ≤ τ] − *f*(σ_j) )²

ĝ*ₛ* predicts the *magnitude* of the error; ĥ*ₛ* predicts the *probability* of
being within τ. Same data, same fitting machinery, different response variable.

### Two tiers of fitting data, and why one is not enough

TECRDB covers well-measured central metabolism, which is exactly the **low-σ**
regime. It cannot constrain the range the model must actually work over:

| source | anchor σ p50 | anchor σ p90 | database σ p50 | database σ p90 |
|---|---:|---:|---:|---:|
| dGPredictor-ModelSEED | 0.91 | 1.22 | 21.17 | 52.89 |
| eQuilibrator | 0.36 | 0.70 | 0.78 | — |
| Group Contribution | — | 13.06 | 10.28 | — |

75.6% of database reactions for dGPredictor-ModelSEED (43.4% eQuilibrator, 27.8%
Group Contribution) lie beyond the anchor's σ p90. Fitting on the anchor alone
and clipping would assign those the error learned at σ ≈ 1.2 — underestimating
error precisely where a source is least reliable, which is backwards for a
safety filter.

So each fit uses two tiers:

| tier | response ε_j | weight | rationale |
|---|---|---:|---|
| **anchor** | \|ΔG*ₛ* − ΔG\*\| on 𝓐 | 3 | a measurement |
| **proxy** | \|ΔG*ₛ* − ΔG_ref\| where the reference source is inside its trusted-σ band | 1 | an *upper bound* on the error, not a measurement — hence the lower weight |

The proxy reference is eQuilibrator for GC and DG, and dGPredictor-ModelSEED for
EQ. "Trusted σ band" means σ_EQ ≤ 0.70 or σ_DG ≤ 1.22 — each source's anchor p90,
i.e. the range where the measurements actually constrain it. TECRDB establishes
that eQuilibrator below σ 0.70 is accurate to a median 0.45 kcal/mol; that is
what earns it the right to stand in as a reference where gold data runs out.

### Fitted ĥ

| source | anchor n | proxy n | *p* at min σ | *p* at max σ | knots | anchor fraction within τ |
|---|---:|---:|---:|---:|---:|---:|
| eQuilibrator | 794 | 4,011 | 0.999 | 0.000 | 24 | 0.856 |
| dGPredictor-ModelSEED | 802 | 11,183 | 0.988 | 0.000 | 44 | 0.848 |
| Group Contribution | 802 | 10,025 | 0.845 | 0.511 | 18 | 0.566 |

Group Contribution's curve is nearly flat, spanning only 0.845 → 0.511. **That is
the honest output, not a failure of the fit**: GC's σ correlates with its
measured error at only ρ = +0.176. The consequence is structural — GC can never
reach GOLD from its own σ, and in the results it never does.

*(Using ê instead of p would be worse still. On the Convention A rebuild GC's ê
spans only 3.04 → 5.70 kcal/mol across the entire database — a curve that cannot
rank anything, and one that puts GC below any ê threshold under 3 for the wrong
reason.)*

## 3.3 Cross-source consistency

Sources are treated as independent measurements of the same quantity and
combined by precision weighting, with **ê rather than σ** as the scale — that is
what makes them commensurable. Definitions of *w*ₛ, ΔḠ, χ², R and *z*ₛ are in
[§2.6](02_notation.md#26-fusion-and-consistency-statistics).

R is the **PDG scale factor** for discrepant measurements: R ≈ 1 means the spread
among sources is what their stated uncertainties predict; R ≫ 1 means at least
one is wrong, and *z*ₛ identifies which.

Validated against the anchor — R ranks true error where ê does not:

| R | n | median \|ΔḠ − ΔG\*\| | mean | fraction ≤ 2 |
|---|---:|---:|---:|---:|
| R ≤ 1 | 568 | **0.36** | 0.66 | 93% |
| 1 < R ≤ 2 | 157 | 0.64 | 1.07 | 84% |
| 2 < R ≤ 5 | 59 | **3.18** | 3.04 | 42% |
| R > 5 | 11 | **5.69** | 17.73 | **0%** |

A 16× monotone spread in median error. For contrast, the ê accept/reject split
separates dGPredictor's accepted from rejected reactions by 0.46 vs 0.60 —
essentially not at all.

Database-wide, of the 28,107 reactions with *n* ≥ 2: 13,071 at R ≤ 1, 6,729 at
1 < R ≤ 2, 4,286 at 2 < R ≤ 5, 890 discrepant at R > 5.

### Corroboration is used asymmetrically — this is the key design decision

> **Agreement lifts BRONZE to SILVER and never creates GOLD.
> Being outvoted costs one tier.**

Agreement between two fallible predictors is weak evidence: eQuilibrator and
Group Contribution share group-contribution lineage so they can be wrong the same
way, and 11% of the R ≤ 1 set are structural zeros (Z = 1) where agreement is
imposed by the stoichiometry. Disagreement is strong evidence: someone is
definitely wrong, and *z*ₛ names them.

This was tested, not assumed. Allowing corroboration to promote all the way to
GOLD grew eQuilibrator's GOLD column from 2,443 to 9,157 reactions but diluted
its measured guarantee from 94% to 90% within 2 kcal/mol (dGPredictor-ModelSEED
98% → 91%). Meanwhile the demotion half carries most of the discriminating
power: without it, Group Contribution's BRONZE tier has a median error of 1.66
kcal/mol — indistinguishable from its SILVER — and with it, 8.68.

## 3.4 The decision cascade

Applied independently to each *s* ∈ 𝓢. Most-specific-first; each rule writes a
`reason` slug so every label is auditable back to the rule that produced it.

```
Rule 0  UNGRADED   s ∉ F(i)          — no value, or vetoed (§1.4)

Rule 1  MEASURED   i ∈ 𝓐 :   ε_s(i) ≤ 1  → GOLD
                              ε_s(i) ≤ 3  → SILVER
                              otherwise   → BRONZE
                   terminal — a measurement outranks every inference below

Rule 2  BASE       p_s(i) ≥ 0.90 → GOLD    "self-certain"
                   p_s(i) ≥ 0.70 → SILVER  "self-confident"
                   otherwise     → BRONZE  "uncorroborated"

Rule 3  FLOOR      if currently BRONZE and n(i) ≥ 2 and Z(i) = 0
                     and R(i) ≤ 1.5 and z_s(i) ≤ 1
                   → SILVER  "corroborated"        (never above SILVER)

Rule 4  DEMOTE     if n(i) ≥ 2 and R(i) > 2 and z_s(i) > 3
                   → one tier down  "outvoted"
```

TECRDB is graded separately and trivially: GOLD wherever it exists, except that
`skeleton`-tier matches are capped at SILVER. The measurement is gold in both
tiers; the *match* to a ModelSEED reaction is not, and a skeleton match is blind
to stereochemistry — the mechanism by which hexokinase/aldose cases were
conflated to glucose data. `--tecrdb-skeleton-gold` disables the cap.

Thresholds as shipped: `p_gold 0.90, p_silver 0.70, r_corrob 1.5, z_corrob 1.0,
r_outvote 2.0, z_outvote 3.0, meas_gold 1.0, meas_silver 3.0`
(`tables/grade_calibration.json`).

## 3.5 Validation

Grades recomputed on the anchor with **Rule 1 disabled**, so the label is
inferred from *p*ₛ and the consistency statistics only, then scored against the
measurement it was not allowed to see.

| source | grade | n | median \|ε\| | mean | within 1 | within 2 | p90 |
|---|---|---:|---:|---:|---:|---:|---:|
| eQuilibrator | GOLD | 246 | **0.32** | 0.56 | 87% | 94% | 1.45 |
| | SILVER | 529 | 0.46 | 0.95 | 68% | 85% | 2.54 |
| | BRONZE | 14 | **3.33** | 3.54 | 0% | **0%** | 4.83 |
| dGPredictor-ModelSEED | GOLD | 184 | **0.32** | 0.56 | 78% | 98% | 1.36 |
| | SILVER | 608 | 0.55 | 1.35 | 65% | 82% | 3.09 |
| | BRONZE | 10 | **20.78** | 18.80 | 0% | **0%** | 21.45 |
| Group Contribution | GOLD | 0 | — | — | — | — | — |
| | SILVER | 517 | 1.28 | 1.62 | 42% | 69% | 3.59 |
| | BRONZE | 285 | **8.68** | 7.73 | 33% | 34% | 15.23 |

*(Verbatim from `tables/grade_calibration.json` → `validation`.)*

Monotone in every column, for every source, on data withheld from the label.
GOLD → BRONZE separates by 10× for eQuilibrator and 65× for
dGPredictor-ModelSEED. Group Contribution never reaches GOLD by inference, as
§3.2 predicted.

The one non-monotonicity worth naming: Group Contribution's BRONZE tier has a
*higher* within-1 rate (33%) than its own median suggests, because that tier is
bimodal — a third of it is nearly exact and the rest is badly wrong. The median
and the p90 tell the story the mean hides.

**As a trust label, the grade works.** That claim is what §3.5 supports, and it
is the only claim made for it.

## 3.6 Results

Database-wide, over 56,002 non-EMPTY reactions:

| source | GOLD | SILVER | BRONZE | UNGRADED |
|---|---:|---:|---:|---:|
| TECRDB | 802 | 748 *(skeleton match)* | 0 | 54,452 |
| eQuilibrator | 2,423 | 13,514 | 4,122 | 35,943 |
| dGPredictor-ModelSEED | 5,808 | 12,404 | 13,201 | 24,589 |
| Group Contribution | 309 | 12,941 | 14,063 | 28,689 |

Reason breakdown for dGPredictor-ModelSEED, as an example: 5,271 self-certain,
2,094 self-confident, 10,123 corroborated, 11,175 uncorroborated, 1,948 outvoted,
802 measured.

Per reaction, taking the best grade available: **6,771 GOLD, 16,332 SILVER,
10,186 BRONZE**, and 22,713 with no source at all.

Full counts and reason breakdowns: [`tables/grade_frontier.tsv`](tables/grade_frontier.tsv).
Fitted curves, thresholds and the validation table above:
[`tables/grade_calibration.json`](tables/grade_calibration.json).

## 3.7 What the grade is used for downstream

1. **A quality floor on a direction map.** `graded_trusted` (§5) drops every
   reaction whose best grade is BRONZE — 10,186 reactions — on the principle
   that a wrong direction is worse than a missing one.
2. **Per-reaction triage.** `core_reaction_grades.tsv` carries the grade of every
   source on every core reaction, for case-by-case curation.
3. **Not for source selection.** See [`04`](04_recommendation_algorithm.md).
