# eQuilibrator vs dGPredictor, reconciled on ModelSEED reaction identity

Companion to `THERMO_SOURCE_AGREEMENT_STRUCTURE.md`, which handled all three
sources and found that most of dGPredictor's apparent disagreement was a KEGG
reaction-id mis-mapping. That mask is applied throughout here. This document
asks the narrower question: **on reactions where both eQuilibrator and
dGPredictor are on firm ground, do they agree, and where they don't, why?**

Nothing in `ModelSEEDDatabase/` is modified. All filtering is read-time.

---

## 1. Methodology

Both numbers are ΔG′° in kcal/mol sitting in the same `thermodynamics` dict, but
they are produced by methods that differ in almost every respect that matters.

### eQuilibrator — component contribution

Noor et al.'s component contribution, run through `equilibrator-api`. Two
estimators are combined:

- **Reactant contribution (RC).** Regresses directly on measured reactions
  (TECRDB equilibrium constants). Any compound inside the span of the training
  reactions is anchored to measurement.
- **Group contribution (GC).** Covers only the component orthogonal to that
  span — compounds the measurements cannot reach.

The split is visible in the output: eQuilibrator publishes a per-compound
uncertainty, small for RC-anchored compounds and large for GC fallbacks, with a
**non-diagonal covariance** so reactions sharing reactants have correlated
errors. Reaction ΔG′° is *exactly* additive over Legendre-transformed compound
formation energies (verified in the equilibrator docs; my own held-out additive
refit reproduces eQuilibrator to a median 0.010 kcal/mol).

How ModelSEED invokes it (`Retrieve_eQuilibrator_Reactions_Energies.py`):

1. Each ModelSEED compound → MetaNetX id **by InChIKey**, with three fallback
   tiers: full key → first two blocks (protonation-blind) → **first block only
   (stereo-blind)**.
2. Reactions are computed **only when every reagent maps** (noted `EQC`);
   partial coverage (`EQP`) yields nothing.
3. The **ModelSEED stoichiometry** is rebuilt in MetaNetX ids and handed to
   `ComponentContribution.standard_dg_prime()` at **pH 7.0, I = 0.25 M,
   298.15 K**.

The critical point: **the reaction eQuilibrator scores is ModelSEED's own.**
Only the compound identities are translated.

### dGPredictor — radius-1 molecular signatures

From the source (`decompose_groups.py`, `predict.py`):

- Every compound is decomposed by `Chem.FindAtomEnvironmentOfRadiusN(m, radius, i)`
  with **radius = 1** — each atom plus its immediate neighbours, serialised to a
  canonical fragment SMILES. A compound is a bag of fragment counts.
- A reaction is the **net fragment-count change**, `rule_df[rid] += molsigna_df[met] * stoic`.
- **Protons and electrons are dropped outright** (`if met == "C00080" or met == "C00282"`).
- The active decomposition path under `__main__` is the **no-stereo** variant
  (`get_rxn_rule_no_stero_remove_TECRDB_mets`).
- The executed fit is `LinearRegression(fit_intercept=False)` on
  `X = S.T @ G` against measured ΔG° — **ordinary least squares, no
  regularisation, no intercept**, with in-sample metrics only. Bayesian Ridge
  exists and yields `dG_std`, and that is the path that produced ModelSEED's
  stored `dG_uncer`.
- ΔG′° is obtained afterwards by adding a per-compound Legendre correction at
  **pH 7, I = 0.1 M, 298.15 K**.

The consequences that matter downstream:

| | eQuilibrator | dGPredictor |
|---|---|---|
| what is scored | the ModelSEED reaction | a **KEGG** reaction |
| well-measured compounds | anchored via RC | no special treatment |
| descriptor locality | reactant-level | **radius-1**, cannot see conjugation or ring strain |
| stereochemistry | InChIKey (tier 1) | **blind** (no-stereo path) |
| protons | pseudoisomer transform | **dropped**, then transformed |
| ionic strength | 0.25 M | 0.10 M |

---

## 2. Reconciliation: joining only on ModelSEED identity

8,768 reactions carry both an eQuilibrator and a (KEGG-mask-surviving)
dGPredictor value. Each source then needs its own provenance audit, because they
fail differently.

**dGPredictor** — the KEGG mask, already applied
(`dgpredictor_kegg_mask.json`, 17,271 reactions withheld).

**eQuilibrator** — two failure modes, found here:

| check | reactions |
|---|---:|
| all reagents matched on the full InChIKey (tier 1) | 6,501 |
| worst reagent matched protonation-blind (tier 2) | 1,961 |
| worst reagent matched **stereo-blind** (tier 3) | 303 |
| two ModelSEED compounds collapse onto one MetaNetX id | **15** |
| eQuilibrator declares **no estimate** (sentinel uncertainty) | **1,148** |

Two of these are defects rather than looseness:

- **MetaNetX collisions (15 reactions).** The retrieval script accumulates each
  side into a dict keyed by MetaNetX id (`lhs[mnx_id] = |coeff|`). When two
  distinct ModelSEED compounds collapse onto one MetaNetX id, the second
  **overwrites** the first instead of summing. The reaction eQuilibrator scored
  is then not the reaction ModelSEED wrote.
- **Sentinel uncertainties (1,148 reactions).** eQuilibrator flags compounds it
  cannot estimate by inflating their variance (the `1e6 * sigmas_inf` term).
  Across all 19,510 stored eQuilibrator records the distribution is **strictly
  bimodal**: real uncertainties top out at 64.3 kcal/mol, sentinels start at
  1,000, and **nothing falls between 100 and 1,000**. 2,748 records (14.1%) are
  sentinels. ModelSEED stores the uncertainty faithfully — but nothing
  downstream reads it, so these are currently used as if they were measurements.
  The largest is ±121,605 kcal/mol.

  **S-adenosyl-L-methionine is the most consequential case**: eQuilibrator
  publishes ±23,900 for it, and it appears in 239 reactions. This fully explains
  why SAM/SAH methyl transfer was the worst-agreeing family in the three-source
  analysis (median |Δ| 15.0 kcal/mol, ρ = 0.13) — eQuilibrator never had an
  estimate for it.

---

## 3. The key subset

Filters, all conjunctive: eQuilibrator tier-1 mapping · no MetaNetX collision ·
eQuilibrator has a real estimate · `status == OK` · non-transport · every
reagent has a structure · no generic R/polymer formula · exactly one KEGG id.

**5,292 reactions.** (`results/eq_vs_dgp/key_subset.tsv`)

---

## 4. How well they agree

| | value |
|---|---:|
| Pearson r | **0.925** |
| Spearman ρ | 0.735 |
| median \|eQ − dGP\| | **2.71 kcal/mol** |
| median (eQ − dGP) | −0.18 kcal/mol (no systematic bias) |
| within 1 / 2 / 5 / 10 kcal/mol | 30.7% / 43.2% / 66.5% / 83.1% |
| concordant (≤ 2 kcal/mol) | 2,287 |
| discordant (> 15 kcal/mol) | 548 |

For comparison, the same subset before the eQuilibrator sentinel filter gave
r = 0.810 and median |Δ| 3.28 — removing 779 reactions eQuilibrator had already
disclaimed accounts for most of the improvement.

### Two distinct failure modes

They should not be conflated, and they matter for different purposes:

- **Absolute error grows with |ΔG|** — from 1.7 kcal/mol at |ΔG| < 1 up to
  12.2 kcal/mol above 100. Matters for energetics.
- **Direction error collapses onto ΔG ≈ 0** — 48% of reactions with
  |ΔG| < 1 kcal/mol disagree on sign, falling to 3% above 100. Matters for
  reversibility calls, which is what this project uses these numbers for.

A reaction near equilibrium is where the two methods most often disagree about
which way it runs, even though they agree on the number to within ~1.7 kcal/mol.

---

## 5. Which reactions, which metabolites, and why

### By reaction class

| best agreement | n | median \|Δ\| | r |
|---|---:|---:|---:|
| phosphoanhydride change | 51 | **0.75** | 0.82 |
| ATP/ADP/AMP | 41 | 0.87 | 0.81 |
| EC 5 isomerase | 335 | 1.39 | 0.28 |
| net proton = 0 | 2,618 | 1.89 | 0.92 |
| EC 3 hydrolase | 956 | 2.02 | 0.43 |

| worst agreement | n | median \|Δ\| | r |
|---|---:|---:|---:|
| **O₂-involving** | 864 | **5.88** | 0.86 |
| net proton ≠ 0 | 2,674 | 3.36 | 0.93 |
| EC 1 oxidoreductase | 2,196 | 3.35 | 0.93 |
| NAD(P)(H) redox | 1,859 | 3.25 | 0.93 |

The **discordant tail is overwhelmingly EC 1** (59.9% of reactions with
|Δ| > 15, vs 41.5% of the subset) and specifically **O₂-dependent oxidases and
monooxygenases acting on plant secondary metabolites** — benzylisoquinoline
alkaloids (N-methylcoclaurine, scoulerine, cheilanthifoline), flavonoids,
phenolics. In the worst cases the two disagree by ~200 kcal/mol *and on sign*:
`rxn03585` vanillate oxidase, eQ −118.0 vs dGP +121.9.

### The O₂ effect is a fixed offset, not a breakdown

Tempting to read the examples above as "dGPredictor fails on O₂". The dose
response says otherwise:

| O₂ consumed | n | median (eQ − dGP) | median ΔG′° |
|---:|---:|---:|---:|
| 0 | 4,428 | −0.11 | −1.0 |
| 1 | 813 | −2.02 | −95.1 |
| 2 | 23 | −4.10 | −194.4 |

A clean **~2 kcal/mol per O₂**, consistent with a single per-compound formation-
energy difference for O₂, and matching the fitted O₂ offset of +1.39 from the
independent least-squares attribution. On a −95 kcal/mol reaction that is 2%
relative error, and **O₂ reactions agree on sign 98.2% of the time** — versus
24.4% sign disagreement for non-O₂ reactions. O₂ chemistry dominates the
*absolute*-error ranking only because those reactions are large. The ~200
kcal/mol alkaloid cases are a genuine but rare tail (1.8%).

### Per-metabolite attribution

eQuilibrator is exactly additive over compounds, so if dGPredictor were too, the
disagreement would decompose into fixed per-metabolite offsets. Fitting that by
least squares (held out a fifth):

**Held-out R² = 0.215.** So only about a fifth of the residual disagreement is
systematic per-metabolite; the rest is reaction-specific. (Before the
eQuilibrator sentinel filter this figure was 0.805 — but that was driven almost
entirely by squalene-cyclase reactions eQuilibrator had disclaimed, and it is
not a real result.)

Largest offsets (kcal/mol, eQ − dGP, ≥15 supporting reactions): palmitoyl-CoA
−27.0, homogentisate −24.9, L-lysine +22.7, prephenate +21.4, caffeoyl-CoA
−16.9, maltose +12.9, L-glutamine −12.6, urea −12.4, retinal −12.2, citrate
−11.7, L-arginine −11.4, L-aspartate −10.2, L-glutamate −10.0 (199 reactions).

Best-agreed, and this is the informative half: **NAD −0.20 (909 reactions),
FADH₂ −0.18, NADH −1.22 (894), NADP −2.08 (959)**, plus H₂O, O₂, pyruvate,
acetate, glucose.

### What explains it — and what doesn't

Testing the structural predictions against the per-metabolite offsets
(Spearman ρ vs |offset|, compounds with ≥10 supporting reactions):

| prediction | ρ | p |
|---|---:|---:|
| fraction of bonds conjugated | **+0.241** | 1.0e-4 |
| conjugated bond count | +0.149 | 0.017 |
| eQuilibrator's own uncertainty | +0.090 | 0.22 (ns) |
| heavy atoms | −0.058 | 0.35 (ns) |
| molecular mass | −0.094 | 0.13 (ns) |
| distinct radius-1 fragments | −0.140 | 0.025 |
| rings | −0.093 | 0.14 (ns) |

**Only conjugation *density* survives**, and modestly. This is the one result
that points straight at construction: radius-1 atom-centred fragments see each
atom's immediate neighbours and nothing more, so delocalisation — which is
non-local by definition — is invisible to the descriptor while eQuilibrator's
measured reactants carry it implicitly.

Three things I expected to matter and which **do not**:

- **Molecular size.** NAD, NADP, NADH and FADH₂ are 44–53 heavy atoms with
  20–26 conjugated bonds and agree to 0.2–2.1 kcal/mol across ~900 reactions
  each. Size is not the problem.
- **eQuilibrator's own RC-vs-GC layer.** Predictive before the sentinel filter
  (4.11 vs 7.35 kcal/mol, p = 6e-5), null after it (2.90 vs 2.93, p = 0.20).
  The earlier signal was the sentinels.
- **Metabolite centrality.** ρ(|offset|, n_reactions) = −0.101, p = 0.11, and
  the binned medians are non-monotonic (3.44 → 2.33 → 3.70 → 2.53). A
  "central metabolism agrees, periphery doesn't" story is *not* supported at the
  metabolite level.

At the **reaction** level, though, the analogous test is strongly positive:
reactions in which *every* reagent is a common metabolite (≥50 reactions) agree
to a median **0.53 kcal/mol vs 2.79** for the rest (p = 5e-17, n = 81 vs 5,211).
That is a conjunction over all reagents, not a property of any one — consistent
with errors accumulating across a reaction rather than concentrating in single
metabolites.

### Method-level tests

| hypothesis | result |
|---|---|
| **H1** ionic strength 0.25 vs 0.10 M | detectable, immaterial: ρ(Δ, Δ∑charge²) = 0.070, p = 3e-7 |
| **H2** dGPredictor drops H⁺ | detectable, immaterial: median Δ −0.47 (H⁺ ≠ 0) vs −0.04 (H⁺ = 0), p = 2e-10; dose response ρ = 0.052 |
| **H3** stereo blindness | **confirmed**: 29.3% of EC 5 isomerases have \|dGP\| < 0.5 kcal/mol vs 6.3% baseline (4.6×) |
| **H4** common-metabolite anchoring | **confirmed**: 0.53 vs 2.79 kcal/mol median \|Δ\|, p = 5e-17 |

H1 and H2 are real but tiny — statistically visible only because n = 5,292. The
ionic-strength mismatch is a genuine unfixed inconsistency in the database (the
two sources are not reported at the same conditions) but it is worth under
0.1 kcal/mol in practice.

**H3 is the sharpest methodological finding.** dGPredictor's active
decomposition path is stereo-free, so a pure stereochemical interconversion has
an all-zero fragment-change vector and *must* return ΔG = 0. Nearly a third of
EC 5 isomerases in the subset show exactly that. eQuilibrator's tier-3 InChIKey
fallback has the same blind spot for 303 reactions — so on epimerases and
racemases the two agree on 0 for the same wrong reason, and their apparent
agreement there should not be read as mutual confirmation.

---

## Recommendations

1. **Filter eQuilibrator's sentinel-uncertainty reactions** the way dGPredictor's
   mis-mapped ones are now filtered. 2,748 records database-wide; they are
   currently indistinguishable from measurements to every downstream consumer,
   and they carry 8× the disagreement.
2. **Fix the MetaNetX collision** in `Retrieve_eQuilibrator_Reactions_Energies.py`
   (`lhs[mnx_id] = |coeff|` should accumulate, not assign). 15 reactions here,
   likely more across the full eQuilibrator coverage.
3. **Treat EC 5 isomerase ΔG from dGPredictor as absent, not zero.**
4. **For reversibility work specifically**, the risk is not the large-|ΔG|
   reactions — it is the near-equilibrium band, where the two sources agree on
   magnitude to ~1.7 kcal/mol but disagree on direction 48% of the time.

## Reproducing

| script | env | output |
|---|---|---|
| `analyze_eq_vs_dgpredictor.py` | `core_models_analysis` | `results/eq_vs_dgp/{reconciliation,key_subset,concordant,discordant,class_breakdown,compound_offsets,mechanism_tests}.tsv` |
| `analyze_eq_dgp_metabolites.py` | `eq3` (rdkit) | `results/eq_vs_dgp/{metabolite_profile,metabolite_predictors}.tsv` |
| `plot_eq_vs_dgpredictor.py` | `core_models_analysis` | `reports/thermoComparison/figures/eq_vs_dgp/*.png` |

## Caveats

- The key subset is 5,292 of 8,768 co-covered reactions. Conclusions are about
  reactions where both methods are on firm ground, by construction — the
  excluded 3,476 are not necessarily wrong, just unattributable.
- Per-metabolite offsets are a least-squares attribution, not measurements. With
  held-out R² = 0.215 the individual values are indicative; the named worst
  offenders are worth checking one at a time before acting on any of them.
- The `--no-stereo` and OLS choices are read from dGPredictor's `__main__`
  block on the public repository. ModelSEED's staged predictions came from an
  opaque upload (commit 3870679) with no generator script in the repository, so
  I cannot verify that the staged values used exactly that path — the stored
  `dG_uncer` implies the Bayesian Ridge branch was used, not plain OLS.
- H3's 29.3% figure counts reactions where dGPredictor returns ≈0. Some
  isomerisations genuinely have ΔG ≈ 0, so that number is an upper bound on the
  blind-spot rate; the 4.6× enrichment over baseline is the defensible part.
