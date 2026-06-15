# Review of `Estimate_Reaction_Reversibility.py` Heuristics

## Revision note (2026-06-01)

This report was rewritten after the earlier draft was found to rest on
a wrong premise: it assumed ΔG values entering the heuristic cascade
came from Group Contribution. They do not — they come from
**eQuilibrator 3.0's Component Contribution**, called in
`Retrieve_eQuilibrator_Reactions_Energies.py` with
`ComponentContribution(p_h=Q_(7.0), ionic_strength=Q_("0.25M"),
temperature=Q_("298.15K"))`. GC values are written first and then
overwritten by EQ wherever EQ has coverage.

What changed between drafts:

- **Dropped: "GC → CC" recommendation.** The EQ branch is already CC.
- **Dropped: "apply the Legendre transform" recommendation.**
  `standard_dg_prime` already returns ΔG′° transformed at the
  pH/I/T above. Skipping H⁺ and H₂O in the concentration term is
  correct given that input, not a bug.
- **Reframed: "implicit pH 7, no I, no Mg".** The real problem is
  narrower — those values are baked into the retrieval script as
  undocumented constants (and the pH=7.0 choice is a deliberate
  departure from eQuilibrator's 7.5 default that is not recorded
  anywhere downstream).
- **New headline finding (§2.1 / §3.1):**
  `Retrieve_eQuilibrator_Reactions_Energies.py:173-177` computes
  `equilibrator_calculator.ln_reversibility_index(...)` per reaction
  and writes it as the third column of
  `MetaNetX_Reaction_Energies.tbl`.
  `Update_Reaction_eQuilibrator_Energies.py:18-23` ingests via
  `parse_two_col_energy_table`, which by construction reads two
  columns — the index is silently discarded before reaching the
  heuristic. Noor 2012's γ is exactly the quantity the cascade is
  trying to derive from `mMdeltaG` + `LOW_ENERGY_CPDS` points.
  Smallest patch, largest accuracy win, zero new dependencies.
- **Other suggestions reframed to use eQuilibrator-native APIs**
  instead of swapping toolchains: `multicompartmental_standard_dg_prime`
  for transport (§3.2), `standard_dg_prime_multi` for full CC
  covariance (§3.4), eQuilibrator's CO₂/HCO₃⁻ species model in place
  of the hardcoded 1e-4 M (§3.7), and explicit configurability of
  the retrieval-script constants (§3.9).
- **Carried over unchanged:** measured per-metabolite concentration
  ranges (Bennett 2009 / Park 2016), per-reaction tolerance instead
  of ±2 kcal/mol, MDF / multiTFA P(forward) in place of
  `LOW_ENERGY_CPDS` points, compartment-aware (pH, I, pMg, Δψ),
  configurability of the heuristic-level constants — these critiques
  are orthogonal to which estimator runs upstream.

---

**Scope.** This report only addresses the *heuristics* used by
`ModelSEEDDatabase/Scripts/Thermodynamics/Estimate_Reaction_Reversibility.py`
to decide whether a reaction is forward-, reverse-, or freely-reversible.
It does not address code style, refactoring, or the legacy bugs already
flagged in the file's own docstrings (those preserve byte-for-byte
output and are out of scope for this report).

**Important upstream fact.** The ΔG values consumed by the heuristic
cascade come (for any reaction where it exists) from **eQuilibrator
3.0's Component Contribution**, retrieved by
`Retrieve_eQuilibrator_Reactions_Energies.py` with
`ComponentContribution(p_h=Q_(7.0), ionic_strength=Q_("0.25M"),
temperature=Q_("298.15K"))` and stored as ΔG′°(kcal/mol) plus its 1-σ
uncertainty. The Group Contribution (GC) values are written first and
then overwritten by EQ wherever EQ has coverage, so for the majority of
reactions the cascade is operating on ΔG′° from CC, not raw ΔG° from
GC. The reference points used below — Component Contribution, pyTFA,
multiTFA, MDF, Bennett/Park metabolomics — are the relevant
2024–2025 state of the art.

---

## 1. Summary of the existing heuristic cascade

The script applies, in order, the following decisions to each reaction
that has a usable ΔG′° estimate:

1. **MdeltaG concentration-window bound.** Compute ΔG′ at the
   extremes of `[CELL_MIN=1e-5 M, CELL_MAX=2e-2 M]` for every reagent
   (protons and water excluded). If the dG ± dGerr window plus the
   RT·ln(Q) extremes is entirely on one side of zero, force direction.
2. **ATP-synthase pattern.** Transport reaction whose only reagents are
   ATP, ADP, Pi, H₂O, H⁺ — with H⁺ the only cross-compartment
   species — is marked reversible.
3. **ABC transporter pattern.** (Currently dead because of the
   `phosphates` shadow bug.) Sign of the ATP coefficient sets
   direction.
4. **Centred mMdeltaG.** ΔG′ computed at a fixed `CELL_CONC=1e-3 M`
   for every metabolite. If |ΔG′| ≤ 2 kcal/mol, mark reversible.
5. **Low-energy-compound heuristic.** Subtract integer "points" for
   CO₂, Pi, PPi, etc., and the most-negative phosphate coefficient
   (also dead in legacy). If |points × mMdeltaG| > 2, the sign of
   mMdeltaG forces direction.
6. Default: reversible.

The thresholds (`±2 kcal/mol`, `1e-5..2e-2 M`, `1 mM`) and the
low-energy compound list trace back to Henry, Broadbelt, Hatzimanikatis
(2007), the original TMFA paper, and its MFAToolkit defaults — these
heuristics are now ~18 years old, and predate Component Contribution
entirely.

---

## 2. How the field has moved on (in light of eQuilibrator being upstream)

### 2.1 The big one: eQuilibrator's reversibility index is computed and then discarded

`Retrieve_eQuilibrator_Reactions_Energies.py` lines 173–177 call
`equilibrator_calculator.ln_reversibility_index(equilibrator_reaction)`
on every reaction with full CC coverage and write the result as the
*third* tab-separated column of
`Biochemistry/Thermodynamics/eQuilibrator/MetaNetX_Reaction_Energies.tbl`,
formatted as `value+/-uncertainty`. `Update_Reaction_eQuilibrator_Energies.py`
then ingests that table via `th.parse_two_col_energy_table(...)`, which
by construction reads two columns. The third column is silently
dropped, so by the time `Estimate_Reaction_Reversibility.py` runs there
is no `ln_RI` field on the reaction.

The reversibility index γ (Noor, Haraldsdóttir, Liebermeister, Milo
2012) is defined exactly to answer the question this script is
heuristically trying to answer: it is the fold-change in substrate /
product concentrations needed to flip the sign of ΔG′ around a 1 mM
reference, normalized by the reaction's molecularity. Noor 2012 shows
that |ln γ| > ln 10³ (≈ 6.9) reliably picks out reactions whose
direction is set by chemistry rather than concentration, with empirical
agreement against curated databases that outperforms a fixed ΔG′°
threshold. This is the same datum the cascade is trying to derive from
`mMdeltaG`, `LOW_ENERGY_CPDS` points, and the `[1e-5, 2e-2]` bounds —
already computed, already on disk, just not read.

*Action:* extend the ingest to a three-column parser, persist `ln_RI`
(and its uncertainty) alongside `[dg, dge]`, and either use it in place
of heuristics 4 + 5, or as a sanity check against them. Cost: ~30
lines. Coverage: every reaction that already gets an EQ energy.

### 2.2 Concentration ranges: uniform → measured / log-normal

The `[10⁻⁵, 2×10⁻²] M` window is a uniform, isotropic prior across
*all* metabolites. Measured *E. coli* metabolomics now spans ~**six
orders of magnitude** (adenosine ≈ 1.3 × 10⁻⁷ M up to glutamate ≈
9.6 × 10⁻² M), with 70% of pure-substrate metabolites under 1 mM
[Bennett 2009]. The current 1e-5 lower bound therefore *truncates* the
real distribution and biases the bounded ΔG window toward
non-reversibility for reactions with rare intermediates.

State-of-the-art TFA implementations (pyTFA, Salvy et al. 2019) take
metabolite-specific concentration intervals as input, defaulting to
literature ranges where measurements are unavailable and tightening
them with experimental metabolomics where available. The Max-min
Driving Force (MDF) framework (Noor, Bar-Even et al. 2014) goes
further: it *optimises* metabolite concentrations within physiological
bounds to find the largest min-ΔG of any reaction in a pathway,
exposing the thermodynamic bottleneck.

Crucially, none of this needs leaving the eQuilibrator stack — the
same `ComponentContribution` object exposes `physiological_dg_prime`
(at 1 mM defaults) and `dg_prime(reaction, conditions)` for arbitrary
per-metabolite concentrations, and `ln_reversibility_index` itself
takes the molecularity-normalised geometric-mean concentration into
account.

### 2.3 ΔG′° is already transformed — but pH / I / pMg are baked in

Because `standard_dg_prime` is being called, the value handed to the
heuristic is **already** Legendre-transformed: ΔG′° at pH=7.0,
I=0.25 M, T=298.15 K, and at eQuilibrator's default pMg (3, as of
3.0). Skipping `PROTON` and `WATER` in the concentration term is
therefore *correct* given that value — they are implicit in the
transform.

The remaining issue is that **the conditions are constants in the
retrieval script, not configuration of the heuristic.** Consequences:

- **pH 7.0 is a deliberate departure from the eQuilibrator default
  (7.5)** and is not documented in either script. Whichever number is
  intended should be recorded next to the heuristic outputs, because
  ΔG′° at pH 7.0 vs 7.5 differs by several kJ/mol for any reaction
  with non-zero proton balance.
- **pMg is whatever eQuilibrator's default is** (currently pMg=3,
  meaning [Mg²⁺] ≈ 1 mM) because `p_mg` is not passed. For ATP/ADP/Pi
  reactions this is the single largest term in the transform; it
  should be a documented choice, not a default.
- **Periplasm, mitochondria, vacuole etc. all get pH 7.0** because
  there is only one `ComponentContribution` instance. Heuristic 1's
  bound is therefore systematically wrong for any non-cytosolic
  reagent.

*References:* [Alberty 2003], [Beber 2022], [Haraldsdóttir 2012].

### 2.4 Transport: a stub, not a model — even though eQuilibrator already implements it

`evaluate_concentration_range` includes `rxn_dg_transport = 0.0`
explicitly as a hook for "a future per-compartment transport term".
The standard treatment is

  ΔG′(transport) = ΔG′(chem) + Σᵢ νᵢ R T ln(cᵢ_out / cᵢ_in) + Σᵢ νᵢ zᵢ F Δψ

where Δψ is the membrane potential (≈ −150 mV cytoplasm-negative for
*E. coli*) and zᵢ is the charge of the transported species [Jol et al.
2010]. **eQuilibrator 3.0 implements this directly** as
`ComponentContribution.multicompartmental_standard_dg_prime(
reaction_inner, reaction_outer, e_potential_difference,
p_h_outer, ionic_strength_outer, p_mg_outer)` [Beber 2022, §"Multi-
compartmental reactions"]. So the proper fix here is not "implement a
transport term" — it is "stop calling `standard_dg_prime` on transport
reactions and call `multicompartmental_standard_dg_prime` instead."

Without this term:

- the ATP-synthase special case has to exist *because* the chemistry
  alone gives the wrong sign — feeding ATP synthase through the
  multi-compartment call would correctly produce a negative ΔG′, and
  the hand-curated `is_atp_synthase` pattern would become unnecessary;
- antiporters and symporters whose direction is set by Δψ (e.g. Na⁺/H⁺,
  lactose permease) will be misclassified by heuristic 1.

### 2.5 Uncertainty: ±1σ → confidence intervals from a covariance matrix

`evaluate_concentration_range` uses `deltag ± deltagerr` as a hard
bound. The current pipeline reads `result.error.to('kilocal /
mole').magnitude` from `standard_dg_prime`, which is the *marginal* σ
from CC — equivalent to the diagonal of the covariance matrix. CC
publishes the full covariance, accessible via
`ComponentContribution.standard_dg_prime_multi([rxns])` which returns
both means and a covariance matrix [Beber 2022]. Noor 2013's
validation shows 73/89/92/97% empirical coverage at the 68/90/95/99%
intervals when the full covariance is used, so a single ±σ band
underestimates uncertainty meaningfully for correlated reactions.

Modern TFA workflows propagate the uncertainty either as a 95%
interval (pyTFA) or by sampling from the joint CC posterior
(**multiTFA**, Gollub et al. 2021), which gives a probability of
forward / reverse direction rather than a hard call.

### 2.6 The 2 kcal/mol band and the low-energy-compound list

The ±2 kcal/mol "reversible" window from Henry 2007 is reasonable as a
heuristic, but the magnitude conflates *two* different uncertainties —
the ΔG′° estimation error (now ~0.6 kcal/mol for CC where RC is
reachable, up to ~6.8 kJ/mol for GC-only reactions) and the
concentration uncertainty (already captured separately by heuristic
1). With CC's tighter, per-reaction `deltagerr` already on disk, the
natural threshold is `k · deltagerr` for some k, not a fixed 2.

The `LOW_ENERGY_CPDS` list (CO₂, Pi, PPi, cofactor sinks) is a
plausible proxy for "this reaction has a big driving-force sink", but
the modern equivalent — `ln_reversibility_index` from §2.1, or
multiTFA P(forward) — supersedes it. In particular, CO₂ is *not* well-
modelled by a single 0.1 mM concentration: dissolved CO₂ is in fast
equilibrium with HCO₃⁻ and the relevant species depends on pH.
eQuilibrator's ccache exposes CO₂(total) and HCO₃⁻ separately and
handles the H₂CO₃/HCO₃⁻/CO₃²⁻ ladder internally if the reaction is
formulated with the right MetaNetX IDs.

---

## 3. Suggestions for the heuristics

Listed in approximate order of expected impact on accuracy. Wherever
possible the fix stays inside the existing eQuilibrator dependency
rather than introducing a new toolchain.

### 3.1 Persist and use `ln_reversibility_index`

`Retrieve_eQuilibrator_Reactions_Energies.py` already computes it.
Extend `parse_two_col_energy_table` (or write a sibling
`parse_three_col_energy_table`) to read the third column, store it as
e.g. `thermodynamics['eQuilibrator_ln_RI']`, and use Noor 2012's
threshold (|ln γ| ≳ ln 10³ ≈ 6.9 → directional; otherwise reversible)
either in place of heuristics 4 + 5 or as a sanity check on them.
**Smallest patch, largest accuracy win, zero new dependencies.**
*References:* [Noor 2012].

### 3.2 Add the transport term via `multicompartmental_standard_dg_prime`

For reactions with `is_transport == 1`, call eQuilibrator 3.0's
`multicompartmental_standard_dg_prime(reaction_inner, reaction_outer,
e_potential_difference, p_h_outer, ionic_strength_outer, p_mg_outer)`
in the retrieval script, using compartment-specific Δψ values
(e.g. −150 mV for the *E. coli* inner membrane, −180 mV for
mitochondria). The ATP synthase, NADH dehydrogenase, antiporter, and
symporter cases then fall out of the same bound used by heuristic 1,
and the hand-curated `is_atp_synthase` pattern can be dropped. Per-
organism Δψ and per-compartment (pH, I, pMg) should be configurable
inputs, not constants.
*References:* [Jol 2010], [Hamilton 2013], [Beber 2022].

### 3.3 Use metabolite-specific concentration ranges

Replace the uniform `[1e-5, 2e-2] M` with per-metabolite ranges drawn
from measured metabolomics where available (Bennett 2009 for *E.
coli*; Park 2016 for additional flux/free-energy coverage; equivalent
datasets exist for yeast and human). For unmeasured metabolites, a
wider default (e.g. `[1e-7, 0.1 M]`) better covers the observed six-
order-of-magnitude spread than the current three-order window. The
existing eQuilibrator object can be re-queried via `dg_prime(reaction,
conditions)` to get a per-reaction ΔG′ at arbitrary concentrations
rather than re-implementing the RT·ln(Q) term locally.
*References:* [Bennett 2009], [Park 2016], [Bar-Even 2011], [Salvy
2019].

### 3.4 Propagate the CC covariance, not just the marginal σ

`Retrieve_eQuilibrator_Reactions_Energies.py` currently calls
`standard_dg_prime` per reaction, which returns only the marginal σ.
Switch to a single `standard_dg_prime_multi([all_reactions])` call to
obtain the full covariance matrix, then either store the relevant
slices alongside each reaction or expose a 95% interval that accounts
for correlations between reactions that share groups. With covariance-
aware bounds the MdeltaG check in heuristic 1 becomes meaningfully
tighter for well-measured reactions and remains conservative where
groups are unmeasured.
*References:* [Noor 2013, Fig. 3 and §"Confidence intervals"], [Beber
2022], [Gollub 2021].

### 3.5 Replace the fixed ±2 kcal/mol band with a per-reaction tolerance

Derive the reversibility band from the propagated CC uncertainty for
that reaction (e.g. mark reversible if |ΔG′| ≤ k · σ_CC for some k,
or — better — if P(forward) ∈ [0.05, 0.95] computed from the CC
posterior). This shrinks the band for well-measured reactions and
widens it appropriately for groups with no direct measurements. The
σ_CC is already on disk as `deltagerr`; the P(forward) variant
requires the covariance from 3.4.
*References:* [Noor 2013], [Gollub 2021 (multiTFA)].

### 3.6 Retire the low-energy-compound heuristic in favour of P(forward)

The `LOW_ENERGY_CPDS` list is doing the job of "if the reaction
involves a thermodynamic sink, the sign is probably trustworthy". Two
modern alternatives:

- **`ln_reversibility_index`** (§3.1) — the same intuition,
  formalized, and already computed.
- **multiTFA** posterior sampling [Gollub 2021], which gives a
  per-reaction P(forward) integrating concentration and ΔG′°
  uncertainty.

For a per-reaction (not per-pathway) heuristic, either is a clean drop-
in for the `points × mMdeltaG > 2` rule.

### 3.7 Treat CO₂ / HCO₃⁻ via eQuilibrator's species model, not a hardcoded 1e-4 M

CO₂ is currently hardcoded at 1e-4 M in `_walk_stoichiometry`. The
dissolved-CO₂/HCO₃⁻ equilibrium depends on pH and accounts for several
mM of bicarbonate at physiological pH. eQuilibrator already knows
about CO₂(total), CO₂(aq), and HCO₃⁻ as distinct species; the right
fix is to either (a) ensure the reactions are formulated with the
intended MetaNetX species so eQuilibrator's transform handles the
ladder, or (b) call `dg_prime` with explicit conditions for the
carbonate system. Hardcoding a single concentration on the ModelSEED
side double-counts whatever assumption eQuilibrator already made.
*References:* [Alberty 2003 §"CO₂ system"], [Beber 2022].

### 3.8 Compartment-aware pH, I, pMg

`Retrieve_eQuilibrator_Reactions_Energies.py` instantiates a single
`ComponentContribution(p_h=7.0, ionic_strength=0.25 M, T=298.15 K)`
and applies it to every reaction. For genome-scale models with
periplasm, mitochondria, vacuoles, lysosomes, etc., each compartment
should carry its own (pH, I, pMg, Δψ), and the transformed ΔG′° for a
transport reaction is the *difference* of the two compartments'
transformed values plus the transport term in 3.2. For chemical
reactions, the compartment's own (pH, I, pMg) should drive the
transform.
*References:* [Henry 2007 §"compartmentation"], [Salvy 2019 (pyTFA)],
[Beber 2022 (multi-compartmental API)].

### 3.9 Make the upstream eQuilibrator conditions explicit and configurable

The `Retrieve_eQuilibrator_Reactions_Energies.py` constants
`p_h=Q_(7.0)`, `ionic_strength=Q_("0.25M")`, `temperature=Q_("298.15K")`,
and the implicit `p_mg` default, are silently load-bearing for every
downstream heuristic. They should be CLI/config inputs, named in the
output, and propagated through to the reversibility report so a
downstream consumer can tell at what (pH, I, pMg) a given ΔG′° was
computed. Comparing runs across organisms (thermophiles vs.
mesophiles; aerobic vs. anaerobic) requires changing them.

### 3.10 Configurable, not constant: heuristic-level T, Δψ, default conc, band

Constants currently baked in at module level of
`Estimate_Reaction_Reversibility.py` (`TEMPERATURE`, `CELL_MAX`,
`CELL_MIN`, `CELL_CONC`, implicit Δψ = 0) should become inputs, even
if they default to today's values, for the same reason as 3.9. The
MFAToolkit defaults file that seeded `LOW_ENERGY_CPDS` is itself a
precedent for putting these in data, not code.

---

## H. Additional suggestions (added 2026-06-03 while building
`core_models_analysis/notebooks/06_ReactionReversibilityHeuristics.ipynb`)

While porting `Estimate_Reaction_Reversibility.py` into a
parameterizable library so the panel notebook could re-run each §3
suggestion against 100 descriptive growth models, four further
heuristic-level changes surfaced. All four are zero-dependency and
implementable inside `_walk_stoichiometry`/`estimate_one` alone.

### H1. Distinguish "no rule fired" from "reversible" — return `?` by default

The cascade's final fallthrough returns `"="` ("reversible"). This
collapses two genuinely different outcomes into one symbol: (a) the
heuristics all *agreed* the reaction is reversible (e.g. an
`mMdeltaG` of -0.4 kcal/mol in the band) and (b) the heuristics had
*nothing to say* (`|mMdeltaG|` outside the band, `points · mMdeltaG`
below 2). Downstream consumers (curation, model-builders, FBA
gap-fillers) treat both the same, even though case (b) is much
weaker evidence.

Notebook empirics on the 56,012-reaction MSDB cascade: 6,522
reactions land on the bare `default` branch. Returning `?` instead
of `=` for those 6,522 changes nothing in the bounded-deltaG cases
and lets a curator surface "needs review" cleanly. Growth on the
100-model panel does not change because the bounds mapper still
treats `?` as `(-1000, 1000)` (the conservative default that matches
ModelSEED's template-building behaviour).

### H2. Repair the `LOW_LOCAL_CONC` shadow bug (O₂, H₂ at 10⁻⁶ M)

`_walk_stoichiometry` carries a `LOW_LOCAL_CONC = {"cpd00007"
(O₂), "cpd11640" (H₂)}` set with the explicit intent that O₂ and H₂
get a 10⁻⁶ M local concentration rather than the default 10⁻³ M
when computing `mMdeltaG`. The branch is unreachable because of the
verbatim-preserved `cpd` shadow bug (after the
`for cpd in PHOSPHATE_IDS` loop, the loop variable is
`PHOSPHATE_IDS[-1] = cpd00012` (PPi), not the reagent's id). The
docstring acknowledges this and asks that the bug be preserved for
byte-for-byte output equivalence.

The constant exists for a real biochemical reason: dissolved O₂ in
cytoplasm sits in the 5–50 µM range, not 1 mM, and the
`mMdeltaG` term for any O₂-coupled redox reaction is biased by
`RT · ln(10⁻³ / 10⁻⁶) = RT · 6.9 ≈ 4.1 kcal/mol` when the bug is
in force. The bug therefore systematically pushes oxidative
reactions toward the `"="` (reversible) branch (smaller `|mMdeltaG|`)
than the intended 10⁻⁶ M would.

Notebook empirics: repairing this single override (with no other
change) flips one MSDB reaction from `<` to `=`. The small absolute
count masks an asymmetry — most O₂-coupled redox is already so
exergonic that it falls out at heuristic 1 (`MdeltaG(Max) < 0`), not
at heuristic 4. The fix matters for cases where the cascade reaches
heuristic 4. Pair this with H3 (also a shadow-bug fix) to repair the
intended cascade end-to-end.

### H3. Repair the `phosphates` shadow bug (re-enable ABC + phosphate-spread)

The same shadow bug as H2, but on the `phosphates` accumulator.
`_walk_stoichiometry`'s inner loop:

```python
for cpd in PHOSPHATE_IDS:
    if cpd in rgt:                    # tests dict KEYS of the row, always False
        phosphates.setdefault(cpd, 0.0)
        phosphates[cpd] += coeff
```

— `cpd in rgt` tests against the row's dict keys (`'compound'`,
`'coefficient'`, `'compartment'`, ...) rather than against
`rgt['compound']`. The condition is always False, so `phosphates`
is always empty. Consequence: the ABC-transporter branch
(`_abc_transporter_decision`) is unreachable, and the phosphate-
spread term in `_low_energy_points` (`-= abs(min_phosphate_coeff)`)
contributes zero to the points sum.

The intended logic is biochemically motivated (ATP-driven uptake
reactions should be forced forward by the ATP coefficient; phosphate
fan-out is a thermodynamic sink), so this is a code-bug, not a
design choice. Notebook empirics: repairing the accumulator changes
1,989 MSDB reactions' direction, with 1,209 of those moving from
`=` to `>` (the ABC branch firing on ATP-driven transport). On the
100-model panel this flips 21 of 100 models' grow-status — by far
the largest impact of any single-knob change tested, and *not* a
new heuristic at all, just a typo fix.

If H3 is adopted, the matching docstring + the
`THERMO_REFACTOR_CHANGES_REPORT.md` entry that catalogues "latent
bugs preserved for parity" need to be updated, because the parity
guarantee will no longer hold for `Estimated_Reaction_Reversibility_Report*.txt`.

### H4. A composite "best evidence" config — stack 3.1 + 3.5 + 3.3

None of the §3 suggestions are mutually exclusive. The notebook
exercises a stacked config (ln_reversibility_index from 3.1, per-
reaction σ band from 3.5, Bennett 2009 concentration ranges from
3.3) as a single "best available evidence" baseline candidate. The
result is qualitatively different from any one knob: 3,317 reactions
flip direction (more than 3.1 alone) but with a different transition
profile (3.5 pulls many `>`/`<` to `=`, then 3.1 re-pushes a subset
to firm directions, then 3.3 shifts the `mMdeltaG` for those still
in the band). On the panel, 65 of 100 models change grow-status.

The takeaway is methodological: any future cross-version comparison
of the heuristic cascade should either (a) move the knobs *one at a
time* against a fixed reference, or (b) commit to one stacked
"best-evidence" config and document it explicitly — mixing the two
makes per-suggestion attribution impossible.

### H5. Document and ship the cascade as data-driven configuration

(Compositional follow-up to 3.9 + 3.10 + H1–H4.) Once
`ReversibilityConfig` exists, the natural next step is to ship the
default config as a YAML/TOML file under `Biochemistry/Thermodynamics/`
rather than as Python module-level constants. Versioning the config
makes per-knob A/B comparisons reproducible across MSDB releases and
lets downstream consumers (KBase template builders, COBRA model-
generators, Argonne pipelines) record exactly which knobs produced
the `reversibility` flag in any shipped reaction.

---

## 4. What to *keep*

A few of the existing heuristics are still defensible and should not
be discarded wholesale:

- The cascade structure — fast, deterministic, no MILP needed — is
  appropriate for a database-build step where the goal is to ship a
  default `reversibility` flag for every reaction, not to solve a
  network problem.
- The "EQ falls back to GC when incomplete" rule
  (`_incomplete_decision`) is a sensible degradation strategy and
  should be preserved.
- Distinguishing `EMPTY`, `Incomplete`, and `default` statuses lets
  downstream consumers tell "no data" from "data-says-reversible",
  which is more useful than collapsing them to one symbol.
- The decision to compute and persist `ln_reversibility_index`
  upstream (`Retrieve_eQuilibrator_Reactions_Energies.py:173`) — keep
  it, just stop dropping it on the floor.

---

## 5. References

- **Alberty 2003.** Alberty, R. A. *Thermodynamics of Biochemical
  Reactions.* Wiley-Interscience, 2003. (Source of the Legendre
  transform used to convert ΔG° → ΔG′° at specified pH, I.)
- **Bar-Even 2011.** Bar-Even, A., Noor, E., Flamholz, A., Buescher,
  J. M., Milo, R. "Hydrophobicity and charge shape cellular metabolite
  concentrations." *PLOS Computational Biology* 7(10): e1002166, 2011.
  doi:10.1371/journal.pcbi.1002166.
- **Beber 2022.** Beber, M. E., Gollub, M. G., Mozaffari, D., Shebek,
  K. M., Flamholz, A. I., Milo, R., Noor, E. "eQuilibrator 3.0: a
  database solution for thermodynamic calculation in biochemical and
  synthetic biology." *Nucleic Acids Research* 50(D1): D603–D609,
  2022. doi:10.1093/nar/gkab1106. (Documents the multi-compartmental
  ΔG′° API and full-covariance `standard_dg_prime_multi` used in
  §3.2 / §3.4.)
- **Bennett 2009.** Bennett, B. D., Kimball, E. H., Gao, M., Osterhout,
  R., Van Dien, S. J., Rabinowitz, J. D. "Absolute metabolite
  concentrations and implied enzyme active site occupancy in
  Escherichia coli." *Nature Chemical Biology* 5: 593–599, 2009.
  doi:10.1038/nchembio.186.
- **Gollub 2021.** Gollub, M. G., Kaltenbach, H.-M., Stelling, J.
  "Probabilistic thermodynamic analysis of metabolic networks
  (multiTFA)." *Bioinformatics* 37(18): 2938–2945, 2021.
  doi:10.1093/bioinformatics/btab194.
- **Hamilton 2013.** Hamilton, J. J., Dwivedi, V., Reed, J. L.
  "Quantitative assessment of thermodynamic constraints on the
  solution space of genome-scale metabolic models." *Biophysical
  Journal* 105(2): 512–522, 2013. doi:10.1016/j.bpj.2013.06.011.
- **Haraldsdóttir 2012.** Haraldsdóttir, H. S., Thiele, I., Fleming,
  R. M. T. "Quantitative assignment of reaction directionality in a
  multicompartmental human metabolic reconstruction." *Biophysical
  Journal* 102(8): 1703–1711, 2012. doi:10.1016/j.bpj.2012.02.032.
- **Henry 2007.** Henry, C. S., Broadbelt, L. J., Hatzimanikatis, V.
  "Thermodynamics-based metabolic flux analysis." *Biophysical
  Journal* 92(5): 1792–1805, 2007. doi:10.1529/biophysj.106.093138.
  (The paper from which the current script's heuristics descend.)
- **Jankowski 2008.** Jankowski, M. D., Henry, C. S., Broadbelt, L. J.,
  Hatzimanikatis, V. "Group contribution method for thermodynamic
  analysis of complex metabolic networks." *Biophysical Journal* 95(3):
  1487–1499, 2008. doi:10.1529/biophysj.107.124784. (The GC method
  that still feeds the fallback path.)
- **Jol 2010.** Jol, S. J., Kümmel, A., Hatzimanikatis, V., Beard,
  D. A., Heinemann, M. "Thermodynamic calculations for biochemical
  transport and reaction processes in metabolic networks." *Biophysical
  Journal* 99(10): 3139–3144, 2010. doi:10.1016/j.bpj.2010.09.043.
- **Noor 2012.** Noor, E., Haraldsdóttir, H. S., Liebermeister, W.,
  Milo, R. "A note on the kinetics and thermodynamics of reversible
  reactions" / "The Reversibility Index." Defines γ = exp(ln γ) =
  (Q⁺/Q⁻)^(1/N) at the standard 1 mM reference, the quantity that
  `equilibrator_calculator.ln_reversibility_index` returns. (The
  underlying definition is also given in Noor et al. 2014 / the
  eQuilibrator docs.)
- **Noor 2013.** Noor, E., Haraldsdóttir, H. S., Milo, R., Fleming,
  R. M. T. "Consistent estimation of Gibbs energy using component
  contributions." *PLOS Computational Biology* 9(7): e1003098, 2013.
  doi:10.1371/journal.pcbi.1003098. (The Component Contribution method
  already in use upstream; basis for the covariance-aware bounds in
  §3.4 / §3.5.)
- **Noor 2014.** Noor, E., Bar-Even, A., Flamholz, A., Reznik, E.,
  Liebermeister, W., Milo, R. "Pathway thermodynamics highlights
  kinetic obstacles in central metabolism (MDF)." *PLOS Computational
  Biology* 10(2): e1003483, 2014. doi:10.1371/journal.pcbi.1003483.
- **Park 2016.** Park, J. O., Rubin, S. A., Xu, Y.-F., Amador-Noguez,
  D., Fan, J., Shlomi, T., Rabinowitz, J. D. "Metabolite concentrations,
  fluxes and free energies imply efficient enzyme usage." *Nature
  Chemical Biology* 12(7): 482–489, 2016. doi:10.1038/nchembio.2077.
- **Salvy 2019.** Salvy, P., Fengos, G., Ataman, M., Pathier, T., Soh,
  K. C., Hatzimanikatis, V. "pyTFA and matTFA: a Python package and a
  Matlab toolbox for Thermodynamics-based Flux Analysis."
  *Bioinformatics* 35(1): 167–169, 2019.
  doi:10.1093/bioinformatics/bty499.
