# 2. Notation

Every symbol used in `03` and `04`, defined once. Units are kcal/mol throughout
unless stated. "Observable" means it can be read or computed for every reaction;
quantities marked otherwise exist only on the 802-reaction anchor set.

## 2.1 Index sets

| symbol | meaning |
|---|---|
| 𝓡 | the set of non-EMPTY ModelSEED reactions, \|𝓡\| = 56,002 |
| *i* | a reaction, *i* ∈ 𝓡 |
| 𝓢 | the three **predictors**, 𝓢 = {GC, EQ, DG} |
| 𝓢⁺ | predictors plus the measurement, 𝓢⁺ = 𝓢 ∪ {TEC} |
| *s* | a source, *s* ∈ 𝓢 (or 𝓢⁺ where stated) |
| 𝓐 | the **anchor set** — reactions with a `stereo_exact` TECRDB match, \|𝓐\| = 802 |
| 𝓒 | the **core set** — reactions appearing in ≥ 1 Kegg2 core model, \|𝓒\| = 239 |
| 𝓜 | the core models, \|𝓜\| = 5,683 |

## 2.2 Observables

| symbol | meaning | observable? |
|---|---|---|
| ΔG*ₛ*(*i*) | source *s*'s estimate of the standard transformed Gibbs energy of reaction | yes, stored |
| σ*ₛ*(*i*) | source *s*'s **own reported** standard deviation for that estimate | yes, stored |
| **ν**(*i*) | the stoichiometry vector | yes |
| A(*i*) | availability, A(*i*) = { *s* ∈ 𝓢 : ΔG*ₛ*(*i*) is defined and non-sentinel } | yes |
| V(*i*) | the vetoed sources on *i* (§1.4) | yes |
| F(*i*) | the **feasible** sources, F(*i*) = A(*i*) \ V(*i*) | yes |

## 2.3 Truth and error

| symbol | meaning | available on |
|---|---|---|
| ΔG\*(*i*) | the true (measured) ΔG′° | 𝓐 only |
| σ\*(*i*) | the experimental standard deviation of that measurement | 𝓐 only |
| ε*ₛ*(*i*) = \|ΔG*ₛ*(*i*) − ΔG\*(*i*)\| | the **true absolute error** of source *s* | 𝓐 only |

The distinction between σ (a claim a source makes about itself, before anyone
compares it to anything) and ε (what is actually wrong with it) is the axis the
whole method turns on.

## 2.4 The cascade

The reversibility cascade is a deterministic function

> **𝒞( δ, e ; i ) ∈ { '>', '<', '=' }**

returning the direction operator for reaction *i* when supplied an energy δ and
a reported error *e*. Internally it is the ordered list
`atp_synthase → abc_transporter → stored_bounds → mmdeltag_band → low_energy →
default`, first match wins. `'>'` is forward-only, `'<'` reverse-only, `'='`
reversible.

| symbol | meaning |
|---|---|
| Λ*ₛ*(*i*) = 𝒞(ΔG*ₛ*(*i*), σ*ₛ*(*i*) ; *i*) | the direction source *s* implies |
| Λ\*(*i*) = 𝒞(ΔG\*(*i*), σ\*(*i*) ; *i*) | the **reference direction** — the same cascade fed the experiment. Defined on 𝓐 |

Λ\* is the yardstick for every direction-accuracy number in this folder. It is
deliberately *not* a literature-curated direction: holding the cascade fixed and
varying only the energy is what isolates the contribution of the thermodynamic
source, which is the question being asked.

Three quantities the cascade computes from stoichiometry alone, reused in §4:

| symbol | definition |
|---|---|
| *a*(*i*) | RT·(pdt_max + rct_min) — the concentration term in the stored-bounds **maximum** |
| *b*(*i*) | RT·(pdt_min + rct_max) — the same for the **minimum** |
| *c*(*i*) | RT·rgt_sum — the concentration term in mMΔG, so mMΔG = ΔG′° + *c* |
| *P*(*i*) | the low-energy **points** score (phosphate spread + low-energy compound coefficients); may be negative |

with RT = 0.0019858 · 298.15 kcal/mol and the concentration range 10⁻⁵–2×10⁻² M.

## 2.5 Calibrated quantities

Fitted maps from a source's own σ to something on a common, meaningful scale.

| symbol | meaning |
|---|---|
| ĝ*ₛ* | the **magnitude calibration**, a non-decreasing map σ ↦ expected \|error\| |
| ê*ₛ*(*i*) = ĝ*ₛ*(σ*ₛ*(*i*)) | predicted expected absolute error of source *s* on reaction *i* |
| ĥ*ₛ* | the **probability calibration**, a non-increasing map σ ↦ P(\|error\| ≤ τ) |
| *p*ₛ(*i*) = ĥ*ₛ*(σ*ₛ*(*i*)) | probability source *s* is within τ of the truth on reaction *i* |
| τ | the tolerance defining "close enough", **τ = 2.0 kcal/mol** |
| *k*ₛ | a single per-source scalar correcting the Gaussian scale (§4.3) |
| τ*ₛ*(*i*) = *k*ₛ · ê*ₛ*(*i*) / √(2/π) | the **calibrated uncertainty** — a standard deviation on a scale comparable across sources |

τ = 2.0 is not free: it is the half-width of the cascade's own reversible band
(`mmdeltag_band_heuristic`, `−2.0 ≤ mMΔG ≤ 2.0`), so *p*ₛ reads as "the
probability this number is good enough to call the direction".

The conversion ê → τ uses the Gaussian identity **E\|X\| = τ√(2/π) ≈ 0.798 τ**
for X ~ N(0, τ²); *k*ₛ then corrects the fact that the sources are not Gaussian
and ê is fitted partly on an upper bound.

## 2.6 Fusion and consistency statistics

Computed per reaction over the feasible sources, with *n*(*i*) = \|F(*i*)\|.

| symbol | definition |
|---|---|
| *w*ₛ(*i*) = 1 / ê*ₛ*(*i*)² | precision weight — the **calibrated** error, never the raw σ |
| ΔḠ(*i*) = Σₛ *w*ₛΔG*ₛ* / Σₛ *w*ₛ | the precision-weighted combination |
| χ²(*i*) = Σₛ *w*ₛ (ΔG*ₛ* − ΔḠ)² | the weighted dispersion, df = *n* − 1 |
| **R(*i*) = √(χ² / (n − 1))** | the **Birge ratio** (PDG scale factor). R ≈ 1 ⇒ the spread among sources is exactly what their stated uncertainties predict; R ≫ 1 ⇒ at least one is wrong |
| *z*ₛ(*i*) = \|ΔG*ₛ*(*i*) − ΔḠ(*i*)\| / ê*ₛ*(*i*) | the per-source standardised residual — which source is the outlier |
| Z(*i*) | the **structural-zero** flag: every feasible source reports \|ΔG′°\| < 0.5, or the reaction is a transport reaction |

ΔḠ is an internal construct — a reference point for computing *z*ₛ — not a
recommended value. Nothing in this folder ships a fused energy.

Z matters because agreement can be free rather than earned: 1,445 of the 13,071
reactions at R ≤ 1 are transport or net-cancelling stoichiometries where every
source reports ≈ 0 by construction.

## 2.7 Outputs

| symbol | meaning |
|---|---|
| G*ₛ*(*i*) ∈ {GOLD, SILVER, BRONZE, UNGRADED} | the **grade**: a trust label on ΔG*ₛ*(*i*) (§3) |
| ρ*ₛ*(*i*) = 1 − P( Λ\*(*i*) = Λ*ₛ*(*i*) ) | the **direction risk** of source *s* (§4.4) |
| *T* | the recommendation **target**, *T* ∈ {magnitude, direction} |
| *s*\*(*i*) | the **recommended** source under target *T* (§4) |
| *E*\* | the abstention tolerance; ê ≤ 2.0 for magnitude, ρ ≤ 0.35 for direction |

## 2.8 Naming collision, resolved

`optimize_thermo_source_assignment.py` already uses "gold" and "silver" for its
two **calibration-data tiers**. Those are fitting *inputs*. Throughout this
folder they are called **anchor** (TECRDB measurements, weight 3) and **proxy**
(a trusted-σ source standing in as a reference, weight 1); gold/silver/bronze
refer only to the emitted grade G*ₛ*.
