# Choosing a thermodynamic source per reaction, across all of ModelSEED

**What this is.** An algorithm that decides, for every ModelSEED reaction, which
of the three thermodynamic sources to use for its ΔG′° — or that none of them is
good enough. It is calibrated against experimental measurements and validated on
data held back from that calibration.

**Result in one line.** At an expected-error tolerance of 2 kcal/mol it assigns a
source to **11,576 reactions** with a held-out mean error of **1.03 kcal/mol**,
against **1.75** for the priority rule ModelSEED `dev` currently uses.

Companion to `EQUILIBRATOR_VS_DGPREDICTOR_MODELSEED.md`, which established the
inputs this depends on. Code: `scripts/optimize_thermo_source_assignment.py`.
Figure: `figures/eq_vs_dgpms/fig8_source_assignment.png`.

---

## 1. The problem

Three methods estimate reaction free energies. Each covers a different,
overlapping slice of the database, and where they overlap they often disagree —
sometimes enough to reverse the reaction's predicted direction.

| source | reactions covered | share of database |
|---|---:|---:|
| dGPredictor-ModelSEED | 31,924 | 57.0% |
| Group Contribution | 25,812 | 46.1% |
| eQuilibrator | 25,028 | 44.7% |
| **at least one source** | **32,466** | **58.0%** |
| all three | 20,672 | 36.9% |

Crucially, **no source is uniformly best.** Against experimental values the
median absolute errors are eQuilibrator 0.45, dGPredictor-ModelSEED 0.47, Group
Contribution 1.60 kcal/mol. The first two are near-tied overall, but not on the
same reactions — so picking per reaction is worth doing.

## 2. The algorithm in plain terms

Think of three property appraisers. Each has valued a different overlapping set
of houses, and each states how confident they are in each valuation. You want
the best number for every house.

1. **Find houses that actually sold.** For 802 reactions there are real
   laboratory measurements. That is ground truth: you can see exactly how wrong
   each method was.
2. **Learn what each appraiser's confidence is worth.** Not all "I'm confident"
   claims are equal. For each method, learn the curve *when this method reports
   this much uncertainty, how far off does it typically turn out to be?*
3. **Apply that curve everywhere.** Every reaction now has a predicted error for
   each method that covers it — including reactions with no experimental data,
   because the prediction comes from the method's own stated confidence.
4. **Pick, then filter.** Use whichever method is predicted to be least wrong.
   If even that one is above tolerance, assign nothing.

Step 4 is the part that is easy to drop and shouldn't be — see §6.

## 3. The algorithm precisely

### 3.1 Notation

| symbol | meaning | units | observable? |
|---|---|---|---|
| *i* | a ModelSEED reaction | — | — |
| *s* | a source ∈ {EQ, DGPMS, GC} | — | — |
| ΔG*s*(*i*) | source *s*'s estimate of the reaction's ΔG′° | kcal/mol | yes, stored |
| **σ*s*(*i*)** | source *s*'s **own reported standard deviation** for that estimate | kcal/mol | yes, stored |
| ΔG\*(*i*) | the **true** ΔG′° | kcal/mol | only for 802 reactions |
| ε*s*(*i*) = \|ΔG*s*(*i*) − ΔG\*(*i*)\| | the **true absolute error** | kcal/mol | only for those 802 |
| **ê(*i*,*s*)** | **estimate of E[ε*s*(*i*) \| σ*s*(*i*)]** — the predicted error | kcal/mol | computed for all |
| ĝ*s* | the fitted calibration curve, so ê(*i*,*s*) = ĝ*s*(σ*s*(*i*)) | — | fitted |
| *E*\* | the error tolerance the user chooses | kcal/mol | chosen |

The whole method is the gap between rows 2 and 3 of the bolded entries: **σ is a
claim the source makes about itself; ε is the truth; ê is a calibrated bridge
from the claim to the truth.**

### 3.2 What σ actually is, per source

σ is *not* a measured error. It is each method's internal estimate of the
dispersion of its own answer, computed before anyone compares it to anything.
ModelSEED stores it as element [1] of the `thermodynamics` triple
`[ΔG, σ, operator]`, in kcal/mol.

| source | how σ is produced |
|---|---|
| **eQuilibrator** | Propagated from the component-contribution covariance matrix. For a reaction with stoichiometry vector **ν**, σ² = **ν**ᵀ **Σ** **ν**, where **Σ** is the covariance of the compound formation energies. Compounds it cannot estimate enter through a separate term inflated by 10⁶, which is the sentinel. |
| **dGPredictor-ModelSEED** | The BayesianRidge **posterior predictive standard deviation**, σ = √Var[ΔG \| **x**, 𝒟], for that reaction's group-change feature vector **x**. |
| **Group Contribution** | Propagated uncertainty of the fitted group energies. |

These are three different quantities that happen to share a unit. There is no
reason a priori that "σ = 1" means the same thing to any two of them — and it
does not, which §3.3 measures.

### 3.3 Why ê is needed, and why it is not just σ

If σ were an honest Gaussian standard deviation, the expected absolute error
would follow directly:

> E\|X − μ\| = σ·√(2/π) ≈ **0.798 σ**

So a perfectly calibrated source would need no calibration at all — you could
set ê = 0.798σ and stop. Measured against TECRDB, none of the three obeys this:

| source | median σ | median true error ε | ratio ε/σ | vs 0.798 |
|---|---:|---:|---:|---|
| Group Contribution | 4.35 | 1.60 | **0.368** | σ **overstates** error 2.2× |
| dGPredictor-MS | 0.91 | 0.47 | **0.522** | σ **overstates** error 1.5× |
| eQuilibrator | 0.36 | 0.45 | **1.260** | σ **understates** error 1.6× |

**This is the justification for the whole calibration step.** The three sources'
confidence scales are not comparable — one is optimistic, two are pessimistic, by
different factors. Comparing raw σ across sources to decide which to trust would
systematically favour whichever source happens to be most pessimistic in its
self-reporting, regardless of accuracy. ĝ*s* is fitted **per source**, which puts
all three on the common scale of *expected error against truth*.

(Practical consequence, quantified in §6: on the 7,720 reactions where both
eQuilibrator and dGPredictor are usable, choosing by ê disagrees with choosing by
raw σ on **28.7%** of them.)

### 3.4 The decision rule, term by term

For each reaction *i*, define three nested sets:

```
  A(i) = { s : ΔG_s(i) is defined }                    availability
  V(i) = { s : s is vetoed on i }                      hard overrides, §4
  F(i) = A(i) \ V(i)                                   feasible sources
```

then

```
  s*(i)  = argmin              ê(i, s)                 (1) pick
             s ∈ F(i)

  keep(i) = [ F(i) ≠ ∅ ]  ∧  [ ê(i, s*(i)) ≤ E* ]      (2) filter

  ΔG_merged(i) = ΔG_{s*(i)}(i)          if keep(i)     (3) emit
                 unassigned             otherwise
```

and the objective being maximised is simply

```
  coverage(E*) = | { i : keep(i) } |
```

Reading each line:

- **(1) pick** — among sources that exist and are not vetoed, take the one whose
  *predicted error* is smallest. Not the one with the smallest σ (§3.3), and not
  a fixed priority order.
- **(2) filter** — a reaction is only assigned if even the best available source
  clears the tolerance. This is the step that makes single-source reactions
  useful rather than dangerous (§6), and the one most easily dropped by mistake.
- **(3) emit** — the merged value is *that source's* number, unmodified. Nothing
  is averaged, reconciled or adjusted; the algorithm selects, it does not
  correct.

Because `keep` and `s*` depend only on reaction *i*, **reactions do not
interact**. There is no combinatorial search: given ê the optimum is obtained by
inspection in one pass, and it is exactly optimal, not approximate.

**Why ê cannot be derived from cross-source disagreement alone.** For two
sources *a*, *b*, the observable |ΔG*a* − ΔG*b*| ≤ ε*a* + ε*b* by the triangle
inequality. It bounds the *sum* of their errors and cannot attribute the blame:
a 20 kcal/mol disagreement is consistent with (20, 0), (0, 20) or (10, 10).
Splitting it requires an external ΔG\*, which is what TECRDB supplies.

## 4. Calibration

### Ground truth

TECRDB (NIST, *Thermodynamics of Enzyme-Catalyzed Reactions*), matched to
ModelSEED reactions by the SMILES→InChIKey multiset pipeline in
`/scratch/ctaylor/dgpredictor_tecrdb`. Only the **802 `stereo_exact`** matches
are used — the tier that distinguishes anomers and D/L pairs.

### Why one calibration tier was not enough

The obvious approach — regress each source's observed error on its own reported
σ, using TECRDB — fails, and the failure is instructive.

TECRDB covers well-studied central metabolism, which is exactly the **low-σ**
regime:

| | TECRDB p50 / p90 / max | whole database p50 / p90 / max |
|---|---|---|
| dGPredictor-MS σ | 0.91 / 1.22 / 21.6 | **21.17 / 52.89 / 2039** |
| eQuilibrator σ | 0.36 / 0.70 / 1.1 | 0.59 / 1.58 / 65.3 |
| Group Contribution σ | 4.35 / 6.53 / 11.5 | 5.06 / 9.71 / 387 |

**75.6% of dGPredictor's database reactions lie beyond the TECRDB p90** (43.4%
eQuilibrator, 27.8% Group Contribution). A curve fitted on gold data alone, then
clipped at its edge, assigns all of them the error learned at σ ≈ 1.2 — it is
*most optimistic exactly where the source is least reliable*, which is the
opposite of what a safety filter must do. Measured, that fit gives a Spearman
correlation between σ and true error of **−0.066** for dGPredictor: worse than
useless.

### The two-tier fit

For each source *s* we build a set of calibration points (*x*ₖ, *y*ₖ, *w*ₖ) — an
input σ, a target error, and a weight — and fit ĝ*s* through them.

**The two kinds of target.**

```
  gold    x = σ_s(i),  y = |ΔG_s(i) − ΔG_TECRDB(i)|,       w = 3
          i ranges over the 802 stereo-exact TECRDB matches.
          y is the TRUE error ε_s(i). This is a measurement.

  silver  x = σ_s(i),  y = |ΔG_s(i) − ΔG_ref(i)|,          w = 1
          i ranges over every reaction where source s and the reference ref
          are both present and  σ_ref(i) ≤ τ_ref.
          y BOUNDS ε_s(i); it is not a measurement of it.

          ref(EQ) = DGPMS,  τ_DGPMS = 1.22
          ref(GC) = ref(DGPMS) = EQ,  τ_EQ = 0.70
```

Each τ is that reference's TECRDB σ-p90 — the edge of the range where gold data
actually constrains it. Within it, eQuilibrator is accurate to a median 0.45
kcal/mol against experiment, which is what licenses it to stand in as a
reference at all. Silver only bounds the error (by the triangle-inequality
argument in §3.4), so it carries ⅓ the weight of gold and is accounted
separately in the table below.

**The fit.** ĝ*s* is the weighted **isotonic regression** of *y* on *x*:

```
  ĝ_s  =  argmin        Σ_k  w_k · ( y_k − f(x_k) )²
           f nondecreasing
```

solved exactly by pool-adjacent-violators (PAVA), then evaluated by linear
interpolation and clamped at the fitted endpoints:

```
  ê(i, s) = ĝ_s( σ_s(i) )
```

Two properties motivate this choice over, say, a linear or log-linear fit:

- **Monotone** — the constraint *f* nondecreasing encodes the one thing we are
  certain of: a source reporting *more* uncertainty must never be predicted to be
  *more* accurate. A least-squares line has no such guarantee and can slope the
  wrong way on noisy data, which is exactly what the gold-only fit did
  (ρ = −0.066).
- **Non-parametric** — no functional form is assumed for σ → ε, because there is
  no theory saying it should be linear, quadratic or anything else, and §3.3
  shows the constant of proportionality is not even the same across sources.

| source | gold n | silver n | gold median err | ρ(σ, **fitting target**) |
|---|---:|---:|---:|---:|
| dGPredictor-MS | 802 | 11,183 | 0.475 | **+0.612** |
| eQuilibrator | 794 | 4,011 | 0.454 | **+0.354** |
| Group Contribution | 802 | 9,808 | 1.600 | **−0.082** |

**Read that last column carefully.** It is the correlation against the
*combined fitting target*, which silver dominates (11,183 vs 802 points for
dGPredictor). Silver is `|source − eQuilibrator|`, so +0.612 mostly says *σ
predicts cross-source disagreement well* — which is what makes the extrapolation
behave, but it is **not** the same as "σ predicts true error". §5b measures that
separately, and the number is much smaller.

Group Contribution stays negative: **its self-reported uncertainty carries no
information**, so its `ê` is effectively a constant. That is a finding in its
own right.

### Two hard overrides

Failures a σ-only model cannot see, both established in the companion report:

- **eQuilibrator sentinels** (4,934 reactions) — the source explicitly declaring
  it has no estimate, by inflating its variance. Never assigned.
- **dGPredictor on the quinone/quinol couple** (1,028 reactions) — 52.8%
  sign-wrong. Never assigned.

## 5. Validation

Fitted on a 70% split of TECRDB, scored on the held-out 30% (n = 241), against
every fixed-source policy and the incumbent — `dev`'s
`Promote_Reaction_Thermodynamics_to_Canonical.py`, which takes eQuilibrator then
Group Contribution, then the ML tier, lowest reported error within a tier.

| strategy | median | **mean** | p90 |
|---|---:|---:|---:|
| **this algorithm** | **0.45** | **1.03** | 2.23 |
| always eQuilibrator | 0.47 | 1.74 | 2.23 |
| always dGPredictor-MS | 0.52 | 1.27 | 2.96 |
| dev priority (EQ > GC, then ML) | 0.48 | 1.75 | 2.23 |
| always Group Contribution | 1.42 | 3.59 | 10.08 |

Medians are near-tied. **The gain is in the mean — 41% below the incumbent —
which means the algorithm is avoiding catastrophic cases, not improving typical
ones.** For direction calls that is the failure mode that matters.

## 5b. How close is ê to the true error?

The validation above shows the *assignment* beats the baselines. This asks a
narrower and more sceptical question: as a predictor of a source's actual error,
how good is ê? Measured on held-out TECRDB over 20 random 70/30 splits.

| source | median ê predicted | median actual error | P(actual ≤ ê) | ρ(ê, actual) |
|---|---:|---:|---:|---:|
| dGPredictor-MS | 1.99 | **0.50** | 77.7% | **0.111** |
| eQuilibrator | 1.83 | **0.45** | 76.9% | **0.124** |
| Group Contribution | 5.79 | 1.65 | 75.6% | 0.236 |

Three things follow, and they are not all favourable.

**ê is conservative, by roughly 4×.** It predicts ~2 kcal/mol where the truth is
~0.5, and it is an upper bound on the actual error about 77% of the time. For a
safety filter that is the right direction to err.

**The promise it makes is kept.** When ê says "≤ 2 kcal/mol", the actual error is
≤ 2 in **86.8%** of cases for dGPredictor and **91.1%** for eQuilibrator. This is
the number to quote for "how much can I trust the tolerance".

**But ê does not finely RANK reactions by error inside the well-measured
regime.** ρ(ê, actual) is only 0.11–0.24, and at the decision boundary the
separation is weak:

| source | group | median actual | mean actual | % over 5 kcal/mol |
|---|---|---:|---:|---:|
| dGPredictor-MS | accepted (ê ≤ 2) | 0.46 | 1.25 | 5.4% |
| dGPredictor-MS | rejected (ê > 2) | 0.60 | 1.44 | 6.1% |
| eQuilibrator | accepted (ê ≤ 2) | 0.43 | **0.80** | 0.9% |
| eQuilibrator | rejected (ê > 2) | 0.56 | **2.03** | 2.8% |

For eQuilibrator the filter separates on the mean (0.80 vs 2.03, a 2.5× gap) and
on the tail. For dGPredictor, within TECRDB, it barely separates at all.

**Why, and why this is expected rather than a defect.** TECRDB contains only
well-studied central metabolism, so even its *rejected* reactions are accurate —
median 0.60 kcal/mol. The reactions the filter exists to catch, high-σ exotic
chemistry with errors of tens of kcal/mol, are simply **not in TECRDB**. The
filter cannot be shown to work on a test set that contains none of the cases it
targets.

So the evidence that the filter does its job is indirect, and comes from two
other places: σ tracks cross-source disagreement at ρ = 0.61 across the *full*
σ range (companion report §4), and the assignment beats every baseline on
held-out **mean** error, 1.03 vs 1.75 (§5) — a mean-not-median gain, i.e. it is
removing catastrophic cases.

**The honest one-line summary:** ê is a well-behaved conservative *threshold* —
trust "≤ 2 kcal/mol" at about 87–91% — but it is a weak *ranking*, and its
performance on the extreme reactions it is designed to exclude is inferred, not
measured.


## 6. What the filter does, and why "only one source" is not a free pass

A natural simplification is *if only one source covers the reaction, just use
it; only arbitrate when there are two or more.* The data says otherwise.

**4,644 reactions have exactly one usable source. The algorithm rejects 87.7% of
them:**

| only source available | reactions | kept | rejected | median predicted error |
|---|---:|---:|---:|---:|
| dGPredictor only | 4,232 | 530 | **3,702** | **8.99** kcal/mol |
| eQuilibrator only | 270 | 42 | 228 | 2.57 |
| Group Contribution only | 142 | **0** | 142 | 5.70 |

Auto-accepting these would import ~4,000 reactions at a median expected error of
9 kcal/mol — about 4× the tolerance, and enough to flip a direction call on its
own. Group Contribution alone never clears the bar: 0 of 142.

The reason is intuitive once stated: **a reaction only one method can reach is
usually unusual chemistry, and unusual chemistry is where every method is least
reliable.** "It is the only number I have" is closer to a warning than a
reassurance.

So the algorithm runs uniformly — pick, then filter. With one source the pick is
trivial, but the filter still applies, and for single-source reactions the filter
*is* the entire value.

**Where there is a real choice, the calibration earns its keep.** Among the
7,720 kept reactions where both eQuilibrator and dGPredictor are usable, `ê`
disagrees with the naive "take whichever reports the smaller σ" on **2,215
(28.7%)** — because the two methods' confidence scales mean different things,
which is precisely what TECRDB calibration corrects.

## 7. Coverage: the frontier

| tolerance ê ≤ | reactions | of database | eQuilibrator | dGPredictor-MS | Group Contribution |
|---:|---:|---:|---:|---:|---:|
| 0.5 | 1,537 | 2.7% | 1,537 | 0 | 0 |
| 1 | 6,072 | 10.8% | 1,904 | 4,168 | 0 |
| 1.5 | 8,045 | 14.4% | 3,877 | 4,168 | 0 |
| **2** | **11,576** | **20.7%** | 5,370 | 6,206 | 0 |
| 3 | 19,310 | 34.5% | 12,377 | 6,898 | 35 |
| 5 | 21,638 | 38.6% | 14,004 | 7,447 | 187 |
| 7.5 | 29,136 | 52.0% | 15,322 | 7,989 | 5,825 |
| 10 | 30,685 | 54.8% | 15,497 | 9,363 | 5,825 |
| no limit | 32,343 | 57.8% | 15,497 | 11,021 | 5,825 |

Two things to read off it:

- **Group Contribution is never chosen below ê ≤ 3.** It only picks up
  assignments at loose tolerances, where nothing better is available.
- **The unassigned reactions are not unreachable**, just not reachable at that
  tolerance. This is why the frontier is published instead of a single cut — the
  operating point is a decision about how much error your application tolerates,
  not a fact about the data.

The shipped default is ê ≤ 2 kcal/mol, chosen as roughly the cascade's own ±2.0
kcal/mol reversible band, so a selected value will rarely flip a direction call
through estimation error alone.

## 8. Using it

```python
from optimize_thermo_source_assignment import load_assignment

a = load_assignment()
# columns: rxn, chosen_source, merged_dg, merged_operator, ehat, kept
#          plus ehat_GC / ehat_EQ / ehat_DGPMS per-source predictions
usable = a[a.kept]
```

| file | contents |
|---|---|
| `results/eq_vs_dgpms/source_assignment.tsv` | per reaction: chosen source, merged ΔG′° and operator, predicted error |
| `results/eq_vs_dgpms/source_assignment_frontier.tsv` | the coverage table above |
| `results/eq_vs_dgpms/source_assignment_models.json` | fitted calibration curves + validation |

Reproduce:

```bash
python scripts/optimize_thermo_source_assignment.py   # fit + assign
python scripts/verify_thermo_source_assignment.py     # 12 gating assertions
python scripts/plot_thermo_source_assignment.py       # fig8
python scripts/verify_ehat_calibration.py             # section 5b: is ehat a good error predictor?
```

`verify_thermo_source_assignment.py` exits non-zero if the assignment fails to
beat any of the four baselines on mean error, if either override did not fire,
if a source was assigned to a reaction that lacks it, if the merged ΔG does not
match the chosen source, or if any calibration curve is non-monotone. Currently
12/12.

## 9. Limitations

- **`ê` is only as good as TECRDB's reach.** 802 gold reactions, all low-σ
  central metabolism. Everything beyond is calibrated against a proxy reference
  that *bounds* a source's error rather than measuring it. More experimental
  coverage of unusual chemistry would improve this more than any algorithmic
  change.
- **`ê` is a threshold, not a ranking** (§5b). It keeps its "≤ 2 kcal/mol"
  promise 87–91% of the time, but correlates only ρ ≈ 0.11–0.24 with true error
  within TECRDB, and for dGPredictor barely separates accepted from rejected
  there. That is because TECRDB holds none of the extreme reactions the filter
  targets — so its benefit on those is inferred from cross-source behaviour, not
  directly measured.
- **It selects, it does not correct.** Nothing is reconciled or repaired; the
  algorithm declines to use a source where it expects it to be wrong.
- **Agreement is not correctness.** Where sources agree they can still be jointly
  wrong; the quinone override is the one place fault is positively attributed.
- **Group Contribution's `ê` is near-constant**, since its σ carries no signal.
  It is effectively ranked by its global median error rather than per reaction.
- **The overrides are hand-specified**, derived from the companion report rather
  than learned. A different snapshot or a retrained model would need them
  revisited.
- Fitted against ModelSEED `dev` @ 34992d39. Refit if either the database or the
  dGPredictor build changes.
