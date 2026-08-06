# eQuilibrator vs dGPredictor-ModelSEED

Where do the two thermodynamic estimators disagree, and why? This reads
top-down: **what kind of chemistry** disagrees, then **which enzymes**, then
**which metabolites**. The metabolite layer comes last deliberately — it is the
layer most easily misread, and two of its quantities need defining before any
claim from it is believable.

**Sources.** `dGPredictor-ModelSEED` (the ModelSEED-retrained variant) and
eQuilibrator, both from upstream `origin/dev` @ **34992d39** (2026-08-04),
extracted read-only to `/scratch/ctaylor/tmp/devsnap`. Nothing modified.
Key subset **n = 11,097**; baseline median |eQ − dGP| = **3.44 kcal/mol**;
17.6% of reactions (1,948) disagree by more than 15 kcal/mol.

---

## 0. Definitions

Everything below is built from five quantities. They are defined here so no
section depends on outside context.

### The disagreement, Δ

Per reaction, `Δ = ΔG′°(eQuilibrator) − ΔG′°(dGPredictor-ModelSEED)`, both in
kcal/mol, both read from the same reaction's `thermodynamics` dict in ModelSEED.
Class and enzyme tables report **median |Δ|** over the reactions in the group —
median rather than mean because the distribution is heavy-tailed (subset mean
7.4 vs median 3.4 kcal/mol).

### The baseline, 3.44 kcal/mol

The same statistic over all 11,097 key-subset reactions at once. It is a
reference point, not a threshold. **"× baseline" = class median ÷ 3.44**, so
1.0× means the class disagrees exactly as much as the database average does.

### Relative error

`class median |Δ| ÷ class median |ΔG′°(eQuilibrator)|` — the disagreement
expressed against the size of the energies in that class.

Two properties to keep in mind. It is a **ratio of two medians**, not the median
of per-reaction ratios; that avoids blow-up from reactions with ΔG ≈ 0 in the
denominator, but it describes a typical class member rather than a typical
reaction. And it takes eQuilibrator's ΔG as the reference scale, which presumes
eQuilibrator is the more reliable of the two — defensible given §2, but assumed.

It exists because absolute error alone is misleading for large reactions: a
class running at |ΔG| ≈ 100 kcal/mol will show a big absolute Δ while being
proportionally accurate.

### "% discordant" and the 15 kcal/mol cut

Fraction of reactions in a group with **|Δ| > 15 kcal/mol**. This is a *chosen*
round number — roughly 4× baseline — not one derived from the direction
heuristics. For reference the cascade's own reversible band is **±2.0 kcal/mol
on mMdeltaG** (`reversibility_heuristics.py:327`), so a disagreement well under
15 can already flip a direction call; 15 is deliberately conservative and
isolates a clear tail.

Nothing hinges on where exactly it sits — the class ranking is stable across the
whole range:

| cut | reactions | share of subset |
|---|---:|---:|
| > 4 kcal/mol | 5,093 | 45.9% |
| > 10 | 2,925 | 26.4% |
| **> 15** | **1,948** | **17.6%** |
| > 20 | 1,370 | 12.3% |
| > 30 | 930 | 8.4% |

### σ — the model's own error bar, not something computed here

`dGPredictor-ModelSEED` reports a **BayesianRidge posterior standard deviation**
for every prediction. It is staged in kJ/mol as `dG_uncer` and stored as element
[1] of the reaction's triple:

```
rxn00001   staged:  {"dG_mean": -15.7618, "dG_uncer": 3.6286}   kJ/mol
           stored:  [-3.77, 0.87, ">"]                          kcal/mol
                      ΔG     σ   operator      3.6286 / 4.184 = 0.87 ✓
```

So "σ ≤ 3 kcal/mol" means *the model itself says it is confident here*. The only
question this report asks about σ is whether that self-report is trustworthy
(§4: it is).

### "Fitted offset" — an inference, because the values do not exist

This one needs the most care, because the obvious reading is wrong twice over.

**First: dGPredictor has no per-compound energies to read.** It emits a reaction
ΔG and nothing else. Checked directly against the compound records:

| source | compounds with a stored formation energy |
|---|---:|
| eQuilibrator | 30,607 |
| Group contribution | 7,021 |
| **dGPredictor (either variant)** | **0** |

Its regression coefficients are per-*fragment* formation energies internally, so
in principle a compound value could be assembled as
`Σ(fragment counts × fragment energies)` — but ModelSEED stores only the
reaction-level predictions, the trained coefficients are not in the repository,
and those fragment energies carry the same identifiability problem one level
down. So per-compound numbers must be **inferred**, not extracted.

**Second: the inference is not a regression of one method's ΔG on the other.**
That would be a single line through the scatter. Instead, for every reaction *r*
take the disagreement `d_r` and solve, by least squares over the
**stoichiometric matrix** S:

```
for each reaction r:   d_r  =  ΔG_eq(r) − ΔG_dgp(r)      (known; 11,097 of them)
solve:                 d_r  ≈  Σ_i  ν_ir · x_i           (6,767 unknowns x_i)
```

`ν_ir` is compound *i*'s stoichiometric coefficient in reaction *r*; `x_i` is
that compound's offset. Solved with `scipy.sparse.linalg.lsqr` on the sparse
11,097 × 6,767 matrix with light ridge damping.

The justification: eQuilibrator's reaction ΔG is *exactly* `Σν·ΔGf`, and
dGPredictor behaves additively enough that their difference should decompose the
same way — which is why the fit reaches held-out **R² = 0.977** at the reaction
level.

`x_i` therefore estimates **how differently the two methods implicitly value
metabolite *i***. It is not dGPredictor's formation energy for that compound.

### …and why a fitted offset can still be meaningless

**The solution is not unique.** S has a null space: any `z` with `S·z = 0` can be
added to `x` and *every predicted reaction value is unchanged*. Element
conservation supplies such `z` for free. Verified directly on this matrix:

| element | max change over all 11,097 reactions | reactions unchanged |
|---|---:|---:|
| C | 0.0 kcal/mol | 100% |
| N | 0.0 | 100% |
| O | 0.0 | 100% |
| P | 0.0 | 100% |
| S | 0.0 | 100% |

Add one kcal/mol per carbon atom to every compound in the database and nothing
observable moves. An individual offset can therefore drift by tens of kcal/mol
with no consequence — that freedom is the **gauge**.

**Consequence:** every metabolite claim in §3 is validated against a model-free
quantity — the *observed* median |Δ| over the reactions actually containing that
compound, versus the 3.44 baseline. The fitted offset alone is never sufficient.
`metabolite_validated.tsv` is the file to use; `compound_offsets.tsv` holds the
raw fit and should not be quoted on its own. Figure: `fig4_gauge.png`.

---

## 1. Which chemistry disagrees

Each reaction is assigned exactly one **organic transformation class** by a
priority cascade over what bonds change (`organic_reaction_types.py`) — what
carrier takes the electrons, whether O₂ is incorporated or reduced, which bond
is made or broken. This is deliberately not EC class: EC 1 lumps hydride
transfer, O₂ insertion, disulfide exchange and quinone reduction together, and
those are four different problems for a fragment-based estimator.

Columns are as defined in §0: median |Δ| is the median absolute difference
between the two ΔG′° values over the class; × baseline compares that with the
subset-wide 3.44 kcal/mol; % discordant is the share above 15 kcal/mol; relative
error is the class's median |Δ| divided by its own median |ΔG′°|.
Figure: `fig1_reaction_class.png`.

| transformation class | n | median \|Δ\| | × baseline | % discordant | rel. error |
|---|---:|---:|---:|---:|---:|
| **Redox: quinone / quinol** | 562 | **50.79** | **14.8×** | **88%** | 181% |
| Group transfer: methyl / C1 | 42 | 13.09 | 3.8× | 43% | 178% |
| Oxygenation: O₂ incorporated | 1,159 | 7.72 | 2.2× | 37% | **8%** |
| Group transfer: glycosyl | 95 | 4.44 | 1.3× | 22% | 158% |
| Addition / elimination | 449 | 3.25 | 0.9× | 12% | 140% |
| Redox: carbon, hydride transfer | 2,719 | 3.18 | 0.9× | 12% | 65% |
| Hydrolysis | 1,515 | 3.05 | 0.9× | 12% | 76% |
| Redox: O₂ as terminal acceptor | 396 | 2.79 | 0.8× | 6% | **11%** |
| Redox: heteroatom (S) | 1,220 | 2.54 | 0.7× | 5% | 54% |
| C–C: carboxylation / decarboxylation | 417 | 2.43 | 0.7× | 7% | 56% |
| Group transfer: amino / amido | 195 | 2.22 | 0.6× | 6% | 109% |
| C–C: aldol / Claisen | 149 | 1.93 | 0.6× | 9% | 58% |
| Isomerisation / rearrangement | 430 | 1.86 | 0.5× | 7% | (373%)* |
| Group transfer: phosphoryl | 73 | **1.74** | 0.5× | **0%** | 35% |

\* Isomerisations have a median |ΔG′°| of 0.5 kcal/mol, so relative error is
meaningless for them — the absolute agreement (1.86) is what counts.

**Read the last two columns together.** Absolute and relative error give
different rankings, and the difference is the whole story for O₂:

- **Quinone/quinol is catastrophic on both** — 50.79 kcal/mol absolute, and
  181% of the reaction's own ΔG. It is not a scaling issue; the answer is wrong.
- **Oxygenation looks bad but isn't.** 7.72 kcal/mol absolute is 2.2× baseline,
  but those reactions have a median |ΔG′°| near 101 kcal/mol, so the relative
  error is **8%** — the best in the table alongside O₂-as-acceptor at 11%. O₂
  chemistry ranks high on absolute error only because O₂ chemistry is large.
- **Phosphoryl transfer is the cleanest**: 1.74 kcal/mol and *zero* reactions
  discordant by more than 15.

Everything below glycosyl transfer sits at or under the baseline.

### What this table does *not* say

It says no class stands out from the others. That is not the same as saying the
classes agree. The test compares each class against the rest of the subset, so
if disagreement were spread evenly, every class would score "no difference" —
not because they agree, but because they are uniformly mediocre. That is close
to what happens. Removing the quinone class entirely barely moves the totals:

| | whole subset | quinones removed |
|---|---:|---:|
| n | 11,097 | 10,535 |
| median \|Δ\| | 3.44 | **3.13** |
| within 2 kcal/mol | 43.2% | **37.8%** |
| disagree > 15 kcal/mol | 17.6% | **13.8%** |
| disagree on **sign** | — | **22.6%** |

Quinones are 5.1% of reactions and 29.3% of the total disagreement — a real,
concentrated lump, but **71% of the disagreement is everywhere else**, diffuse
and not attributable to any transformation class. More than one reaction in five
still disagrees on *direction*, which is what matters for reversibility work.

The structure the chemistry classification fails to find, σ finds (§4), and it
survives removing quinones: with the class excluded, σ ≤ 3 gives median |Δ| 0.78
and **zero** reactions above 15 kcal/mol, against 12.06 and 40.6% for σ > 30.

---

## 2. The dominant failure: the quinone / hydroquinone couple

562 reactions, 5.0× enriched in the discordant tail, and 88% of them disagree by
more than 15 kcal/mol. Figure: `fig2_quinone.png`.

Two subclasses, both the same underlying chemistry — a two-electron, two-proton
aromatic redox couple:

| subclass | n | median \|Δ\| |
|---|---:|---:|
| prenyl-quinone carriers (ubiquinone, menaquinone, phylloquinone) | 291 | **86.54** |
| catecholic dihydroxyarenes (catechol, protocatechuate, homogentisate) | 271 | 37.04 |

### It is a direction error, not a magnitude error

For phenylacetyl-CoA dehydrogenase (`rxn35639`):

```
(1) H2O + (1) Phenylacetyl-CoA + (2) Ubiquinone-8  =>  (1) Phenylglyoxylyl-CoA + (2) Ubiquinol-8

eQuilibrator            −45.71 kcal/mol
dGPredictor-ModelSEED  +141.13 kcal/mol   (σ = 163)
```

A dehydrogenation coupled to quinone reduction is strongly exergonic, so a large
positive value has the **sign** wrong. Across the class, **52.8% of quinone
reactions disagree on the sign of ΔG′°**, against 22.6% for everything else. **dGPredictor-ModelSEED is the source at
fault here**, not eQuilibrator. This is the same species of localised redox
regression as its known disulfide/glutathione-reductase failure against TECRDB.

### The enzyme layer confirms it independently

Same Δ as §1 — median `|ΔG_eq − ΔG_dgp|` — grouped by EC subfamily rather than
by chemistry class, with % discordant again the share above 15 kcal/mol. The
grouping uses the **first** EC number listed for a reaction, truncated to three
fields; 1,353 of 11,097 reactions list more than one EC and are assigned by list
order rather than adjudicated.

Ranking EC subfamilies (≥20 reactions) by disagreement:

| EC | n | median \|Δ\| | % discordant | example |
|---|---:|---:|---:|---|
| **1.3.5** | 24 | **88.04** | **100%** | fumarate reductase |
| **1.1.5** | 81 | 79.34 | 96% | glycerol-3-phosphate dehydrogenase |
| **1.6.5** | 139 | 50.47 | 99% | NADPH:*p*-benzoquinone oxidoreductase |
| 1.13.11 | 119 | 41.50 | 72% | β-carotene 15,15′-dioxygenase |
| 5.5.1 | 36 | 17.75 | 64% | flavanone lyase (decyclizing) |
| 1.14.12 | 47 | 16.41 | 66% | anthranilate 1,2-dioxygenase |

The top three are **EC 1.x.5.x** — and in the EC nomenclature the third digit
`5` *means* "with a quinone or similar compound as acceptor". The classification
arrived at from bond changes and the classification arrived at from enzyme
naming agree without being told about each other.

Below the quinone families, the next tier is ring-cleaving dioxygenases
(1.13.11, 41.50) and aromatic hydroxylases — still aromatic redox chemistry.

### Why this class specifically

Both estimators decompose molecules into local fragments; the retrain uses
radius-1 **and** radius-2 atom environments. Quinone ⇌ hydroquinone is a change
in **delocalisation across the whole ring**, not a change in any local
environment. The aromatic stabilisation energy that dominates the couple is
invisible to a descriptor with a two-bond horizon. eQuilibrator sidesteps this
because ubiquinone and its quinol appear in measured reactions, so its reactant-
contribution layer anchors them to data rather than reconstructing them.

**But the model flags it itself.** These reactions carry a median σ of **80.3
kcal/mol**, and **99.6% land in the low-confidence tier (σ > 30); 0% in the
high-confidence tier.** The failure is real and self-quarantining.

---

## 3. Metabolites — validated, not fitted

Figures: `fig4_gauge.png`, `fig5_metabolites.png`. The fitted offset, why it has
to be inferred rather than read off, and why it can be meaningless are all in
§0 — the short version is that dGPredictor stores **no** compound-level
energies, so per-compound values come from a least-squares solve over the
stoichiometric matrix whose solution is only unique up to element conservation.

Ranked by **observed** disagreement over the reactions containing each compound
(baseline 3.44 kcal/mol), with the fitted offset shown for comparison:

| metabolite | n | fitted offset | **observed** \|Δ\| | × baseline | σ |
|---|---:|---:|---:|---:|---:|
| ubiquinone-8/9/10 + quinols (6 species) | 44–55 | ±45.9 | **86.91** | 25.3× | 80.3 |
| 3-dehydroshikimate | 22 | −10.82 | 63.58 | 18.5× | 64.8 |
| catechol | 21 | −48.94 | 59.22 | 17.2× | 36.8 |
| 2,3-dihydroxybenzoate | 17 | −42.12 | 59.01 | 17.2× | 38.8 |
| 5-dehydroquinate | 23 | −13.70 | 54.98 | 16.0× | 61.9 |
| demethylphylloquinone / phytyl-naphthoquinol | 16 | ±27.53 | 50.47 | 14.7× | 61.9 |
| 2-demethylmenaquinone-8 / -quinol | 29 | ±27.53 | 50.26 | 14.6× | 61.9 |
| phytonadiol (vitamin K₁) | 26 | −26.41 | 47.82 | 13.9× | 68.5 |

Every entry is a quinone, a quinol, or a dihydroxyarene — the metabolite layer
recovers §1's class with no additional input. The shikimate-pathway pair
(3-dehydroshikimate, 5-dehydroquinate) enters because those reactions feed
protocatechuate/catechol chemistry.

### The gauge trap, concretely

These have the **largest fitted offsets in the whole table** and represent **no
real disagreement at all**:

| metabolite | n | fitted offset | observed \|Δ\| | × baseline |
|---|---:|---:|---:|---:|
| arachidonyl-CoA | 20 | **+56.58** | **0.80** | 0.23× |
| arachidonate | 34 | +52.51 | 2.41 | 0.70× |
| **NAD** | 1,799 | **+47.12** | **3.03** | 0.88× |
| NADP | 2,016 | +44.63 | 3.92 | 1.14× |
| NADH | 1,789 | +34.11 | 3.04 | 0.88× |
| dopamine | 17 | −39.92 | 4.50 | 1.31× |

Arachidonyl-CoA has the largest fitted offset of any metabolite and the
*smallest* observed disagreement in the table — the reactions containing it
agree four times better than average. Its C41 skeleton and its appearance in
only 20 reactions, mostly alongside arachidonate and CoA, leave it weakly
constrained, so it absorbs gauge freedom.

NAD is the case worth internalising: fitted +47.12 across 1,799 reactions looks
authoritative, and reporting it would have produced the headline *"the two
methods disagree by 47 kcal/mol on NAD."* They do not. NAD and NADH almost
always appear together, so only their difference is pinned, and reactions
containing both disagree by a median **3.03 kcal/mol** — better than baseline.

---

## 4. σ is calibrated, and that is the practical result

Figure: `fig3_sigma.png`.

The original KEGG-based dGPredictor reported a median ±0.35 on everything,
including reactions where it had scored entirely different chemistry — its error
bar carried no information. The retrain's does: **ρ = +0.672** between reported
σ and observed disagreement (n = 24,542), monotone across the whole range.

| reported σ | observed median \|Δ\| |
|---:|---:|
| 0.88 | 0.48 |
| 10.93 | 2.37 |
| 18.91 | 4.50 |
| 26.31 | 10.03 |
| 37.02 | 13.83 |
| 56.93 | 37.42 |

σ over-states the error consistently (it is conservative), but it orders
correctly, which is all a filter needs. So the model can be tiered on its own
output with no external evidence:

| tier | n | r | median \|Δ\| | within 2 kcal/mol |
|---|---:|---:|---:|---:|
| **high (σ ≤ 3)** | 2,272 | **0.993** | **0.79** | **77.4%** |
| medium (3 < σ ≤ 30) | 6,338 | 0.955 | 3.59 | 31.5% |
| low (σ > 30) | 2,487 | 0.769 | 14.98 | 9.6% |

**On its high-confidence fifth, the two methods are effectively
interchangeable** (r = 0.993) — and not one of those 2,272 reactions disagrees
by more than 15 kcal/mol (§5). Pooled over all tiers: r = 0.857, median 3.44,
median signed −0.42 (no bias).

---

## 5. Does the disagreement land on reactions that matter?

§1 and §4 describe the database. The operational question is narrower: of the
reactions actually present in the 5,683 Kegg2 core models, are they in the
confident half of this comparison or the doubtful half?

"Biologically significant" is operationalised two ways, both from this repo
rather than asserted:

- **prevalence** — how many of the 5,683 core models contain the reaction
  (`site/data/reaction_model_counts.json`, field `all`);
- **direction-sensitivity of growth** — in what fraction of models containing it
  does swapping its bound direction move the FBA objective, aggregated across
  all 5,683 per-model sweeps in `results/reaction_effects_all/`
  (cached to `results/eq_vs_dgpms/rxn_growth_sensitivity.tsv`).

Figure: `fig6_biological_significance.png`.

### By prevalence

| core models containing it | n | median \|Δ\| | in the σ ≤ 3 tier |
|---|---:|---:|---:|
| none | 10,979 | 3.46 | 19.9% |
| 1–499 | 15 | 1.26 | ~53% |
| 500–1,999 | 29 | 0.85 | **75.9%** |
| **≥ 2,000** | 74 | **0.63** | **73.0%** |

### By growth direction-sensitivity

Among the 119 core-model reactions that survive into the key subset:

| direction moves growth in… | n | median \|Δ\| | in the σ ≤ 3 tier |
|---|---:|---:|---:|
| < 25% of models | 86 | 0.86 | 65.1% |
| 25–60% | 20 | 0.59 | 85.0% |
| **> 60% of models** | 13 | **0.15** | **92.3%** |

Core-model reactions overall are **71.4%** high-confidence against **20.5%** for
the subset at large, with a median |Δ| of **0.70 vs 3.47 kcal/mol**.

### Reading

The two measures agree and point the same way: the more a reaction actually
matters — more models carry it, and its direction more often changes predicted
growth — the more confident dGPredictor is and the closer the two methods sit.
Reactions where they differ by tens of kcal/mol are overwhelmingly ones that
appear in **no** core model.

This reframes §1 without contradicting it. The disagreement is real and
widespread across the *database*, but it is concentrated in metabolism that is
not being modelled. For the reactions the reversibility work actually touches,
the two sources are close to interchangeable — and the residual risk there is
direction near ΔG ≈ 0 (§4), not magnitude.

Two limits. The rank correlations are weak (ρ = −0.080 for prevalence vs σ,
ρ = −0.063 vs |Δ|) because 10,979 of 11,097 reactions sit in zero core models,
so the distribution is too skewed for a correlation to mean much — the binned
comparison is the meaningful view. And only 119 of the 239 core reactions reach
the key subset, so this describes those, not all 239.

### Using σ as a filter: how well does it actually work?

Treating "σ above a cut" as a detector for "|Δ| > 15 kcal/mol":

| cut | reactions kept | large differences removed | still > 15 among kept | median \|Δ\| kept |
|---|---:|---:|---:|---:|
| σ ≤ 3 | 20.5% | **100.0%** | **0** | 0.79 |
| σ ≤ 10 | 25.7% | 99.9% | 2 | 1.06 |
| **σ ≤ 20** | **52.6%** | **91.2%** | 171 (2.9%) | 1.76 |
| σ ≤ 30 | 77.6% | 63.8% | 706 (8.2%) | 2.50 |
| σ ≤ 50 | 93.3% | 32.6% | 1,313 (12.7%) | 3.04 |

**Not one of the 2,272 reactions with σ ≤ 3 disagrees by more than 15 kcal/mol.**
Zero. On this dataset σ is not merely correlated with disagreement at the low
end — it is an exclusion guarantee.

Three things follow.

**Discarding only the low tier is not enough.** σ > 30 catches 63.8% of the
large differences; the medium tier still holds 706 of them. A cut near σ ≤ 20
is the better operating point: half the data retained, 91% of large differences
gone.

**The filter is imprecise in the discard direction.** Of the 8,825 reactions
σ ≤ 3 throws away, **50.2% actually agreed within 5 kcal/mol**. You buy the
guarantee by discarding a great deal of perfectly good agreement. It also
*declines to use* dGPredictor rather than reconciling it with eQuilibrator —
nothing is repaired, only withheld.

**On core-model reactions the cost is small and the benefit is total.**

| cut | core reactions kept | all reactions kept | quinone class kept |
|---|---:|---:|---:|
| σ ≤ 3 | **71.2%** (84/118) | 20.5% | 0.5% |
| σ ≤ 20 | **78.8%** (93/118) | 52.6% | 4.3% |
| σ ≤ 30 | 82.2% (97/118) | 77.6% | 22.8% |
| unfiltered | 100% | 100% | 100% |

The filter is far gentler on core-model reactions than on the database at large,
and it removes exactly the right ones. Unfiltered, those 118 reactions have a
maximum |Δ| of **97.58 kcal/mol with 21 above 15**; at any cut from σ ≤ 3 to
σ ≤ 30 the maximum falls to **11.37 with none above 15**.

## 6. Head-to-head against the original, and coverage

On the 7,871 reactions covered by eQuilibrator **and both** dGPredictor
variants, with KEGG-mismapped and eQ-sentinel rows removed from both:

| | r | median \|Δ\| | within 2 |
|---|---:|---:|---:|
| original (KEGG-based) | 0.766 | **2.88** | **42.4%** |
| dGPredictor-ModelSEED | **0.922** | 3.47 | 37.8% |

**The retrain is not uniformly more accurate.** It wins decisively on
correlation — it suppresses the original's catastrophic outliers — and its
typical error is slightly worse.

The case for switching is elsewhere:

| | original | retrain |
|---|---|---|
| what it scores | a **KEGG** reaction | the **ModelSEED** reaction |
| mis-mapping failure mode | 17,271 reactions affected | **impossible by construction** |
| coverage | 27,715 | **31,924** (+~11,400 with no KEGG counterpart) |
| co-covered with eQuilibrator | 15,300 | **24,542** |
| key subset | 5,292 | **11,097** |
| ionic strength | I = 0.10 M | I = 0.25 M (**matches eQuilibrator**) |
| usable confidence signal | none | **yes** |

---

## 7. Method-level hypotheses

| | result |
|---|---|
| **H1** ionic-strength mismatch | **moot** — both now at I = 0.25 M. A residual ρ = 0.120 vs Δ∑charge² persists, so something charge-related remains, but it is not the setting. |
| **H2** proton handling | small: median Δ −0.30 (H⁺ ≠ 0) vs −0.54 (H⁺ = 0); dose-response ρ = 0.143 |
| **H3** stereo blindness | **halved but alive**: 17.2% of EC 5 isomerases return \|ΔG\| < 0.5 vs 4.7% baseline (3.7×), down from 29.3% vs 6.3% with radius-1 only. Radius-2 fragments help; they do not cure it. |
| **H4** measured-chemistry anchoring | **confirmed, stronger**: median \|Δ\| 0.65 where every reagent is a common metabolite (n = 260) vs 3.52 otherwise (n = 10,837), p = 6e-31 |

---

## Recommendations

1. **Adopt `dGPredictor-ModelSEED`; retire the KEGG-based `dGPredictor`.**
   The mis-mapping failure mode disappears, comparable coverage more than
   doubles, and you gain a confidence signal.
2. **Use σ as the quality gate, cutting near σ ≤ 20.** That keeps 52.6% of
   reactions and removes 91.2% of the large differences; σ ≤ 3 removes 100% but
   retains only a fifth, and dropping merely σ > 30 catches just 63.8%. On
   core-model reactions a σ ≤ 20 cut keeps 78.8% while taking the worst |Δ| from
   97.58 down to 11.37 kcal/mol. It also auto-quarantines the quinone failure
   (95.7% of the class removed at σ ≤ 20, 99.5% at σ ≤ 3). Note the filter
   *withholds* rather than reconciles, and about half of what a σ ≤ 3 cut
   discards would have been fine — so prefer the loosest cut that meets your
   tolerance rather than the tightest.
3. **Treat quinone/quinol reactions from dGPredictor as unusable** regardless of
   tier — 562 reactions, sign frequently wrong. Prefer eQuilibrator there.
4. **Still filter eQuilibrator's sentinels** — 4,491 in the co-covered set,
   unchanged on dev — and fix the MetaNetX collision in
   `Retrieve_eQuilibrator_Reactions_Energies.py` (35 reactions; `lhs[mnx_id] =`
   should accumulate, not assign).
5. **Fix the upstream promotion order.** dev's
   `Promote_Reaction_Thermodynamics_to_Canonical.py` picks within its ML tier by
   smallest reported error, so the old overconfident model beats the retrain
   95.3% of the time. Combined with (1), the highest-value upstream fix.

## Reproducing

```bash
MSDB_ROOT=/scratch/ctaylor/tmp/devsnap \
  python scripts/build_thermo_agreement_features.py \
    --dgp-label dGPredictor-ModelSEED --out reaction_features_dgpms_dev.tsv   # eq3 env

export EQDGP_FEATURES=results/thermo_agreement/reaction_features_dgpms_dev.tsv
export EQDGP_OUT=results/eq_vs_dgpms EQDGP_DGP_LABEL=dGPredictor-ModelSEED
export EQDGP_FIGS=reports/thermoComparison/figures/eq_vs_dgpms MSDB_ROOT=/scratch/ctaylor/tmp/devsnap
python scripts/analyze_eq_vs_dgpredictor.py     # reconciliation, tiers, mechanisms
python scripts/analyze_eq_dgp_topdown.py        # layers 1-3 + the gauge demo
python scripts/plot_eq_dgp_topdown.py           # fig1-fig5
python scripts/plot_eq_dgp_biological_scatter.py  # fig6 (reads reaction_effects_all/)
```

## Caveats

- **Per-compound offsets are not identifiable, and are inferred rather than
  read.** dGPredictor stores no compound-level formation energies (0, against
  eQuilibrator's 30,607), so §3's numbers come from a least-squares attribution
  whose solution is unique only up to the stoichiometric null space. Every
  metabolite claim above is backed by the observed-disagreement column. Do not
  quote the fitted offsets alone — `compound_offsets.tsv` exists but
  `metabolite_validated.tsv` is the one to use.
- **§5 covers 119 of the 239 core reactions** — the rest lack one of the two
  sources or fall outside the key subset. Growth direction-sensitivity is a
  single-reaction FBA perturbation, so it measures marginal effect in
  isolation, not importance under simultaneous changes.
- **The 15 kcal/mol discordance cut is chosen, not derived** (~4× baseline). The
  cascade's own reversible band is ±2.0 kcal/mol on mMdeltaG, so real direction
  flips start well below it. Class rankings are stable from >4 to >30 (§0).
- **Layer 2 groups on the first-listed EC number only** (1,353 reactions list
  several).
- The transformation classifier is a rule cascade, so borderline reactions land
  by priority order rather than by chemical adjudication. 1,657 reactions
  (15%) fall to "Other / unclassified" and sit at 0.9× baseline, so nothing
  large is hiding there.
- "dGPredictor is the one at fault on quinones" is argued from the sign of a
  quinone-coupled dehydrogenation, not validated against a measurement.
- The σ filter is validated against eQuilibrator, so "removes the large
  differences" means "removes the reactions where the two sources disagree" —
  not, on its own, "removes the reactions where dGPredictor is wrong". §2 argues
  dGPredictor is the source at fault on quinones; elsewhere the attribution is
  not established.
- The zero-large-differences result at σ ≤ 3 is an empirical property of these
  11,097 reactions, not a guarantee the model carries to new ones.
- eQuilibrator here is dev's re-run, so the original-dGPredictor numbers quoted
  from the earlier report used an older eQuilibrator. §6 holds both constant on
  dev and is the fair comparison.
