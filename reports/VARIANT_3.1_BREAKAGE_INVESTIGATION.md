# Why variant 3.1 (eQuilibrator reversibility index) breaks so many models

**Variant:** 3.1 — replaces the points-and-bands heuristic with eQuilibrator's
`ln(reversibility_index)`, calling a reaction directional when |ln γ| > 6.9
(Noor 2012). See [REACTION_REVERSIBILITY_HEURISTICS_REVIEW.md](REACTION_REVERSIBILITY_HEURISTICS_REVIEW.md) §3.1.
**Question set:** (1) which reactions break the models, (2) does the panel's
~75% break-rate generalize to all ~3.5K growers, (3) does 3.1 ever rescue a
non-grower into growth.
**Compute:** `scripts/analyze_variant_31_breakage.py` →
`results/variant_31_breakage.json`. FBA on all 5,683 `core_models_kegg2`
models, heuristic-baseline rebinding (cascade defaults), single-reaction
knock-ins. No on-disk model JSON or MSDB file was modified.
**Generated:** 2026-06-16.

---

## ⚠ Correction (2026-06-16): the breakage was caused by a sign bug in our code, now fixed

After this investigation, a follow-up question — *is the rxn00199 flip even
biochemically viable?* — exposed the actual root cause: a **sign-inversion bug
in our cascade**, not a real property of the reversibility index.

`scripts/reversibility_lib.py` mapped the eQuilibrator `ln(reversibility_index)`
to a direction with the comparison **backwards**:

```python
direction = ">" if ln_ri > 0 else "<"      # WRONG
direction = ">" if ln_ri < 0 else "<"      # FIXED (2026-06-16)
```

eQuilibrator's `ln(reversibility_index)` carries the **same sign as ΔG′°**
(Noor 2012: ln γ̂ = (2/N)·ΔG′ₘ/RT), so a **negative** index means
forward-favorable — the same convention the cascade's own MdeltaG rule uses
(`stored_max < 0 → ">"`). Verified end-to-end:

- The upstream table `MetaNetX_Reaction_Energies.tbl` stores eQuilibrator's
  raw output faithfully (retrieval script applies no negation, same reaction
  orientation as the MSDB equation). The MSDB data is **correct**.
- Across all 16,742 reactions, `sign(ln_RI) == sign(ΔG′°)` for **96.8%**.
- Among the 9,076 strongly-directional reactions (|ln_RI| > 6.9), the buggy
  mapping agreed with the ΔG′° sign in **1** of them — i.e. it inverted
  **9,075 / 9,076**.
- rxn00199: ΔG′° = −3.33 kcal/mol, ln_RI = −8.40 (matches the hand-computed
  (2/N)·ΔG′ₘ/RT ≈ −8.3). Negative ⟹ forward. The buggy code called it reverse.

**What this means for the numbers below.** Everything in this report
(§Q1–Q3, the 2,377 broken growers, the rxn00199/rxn05145/rxn00256 culprits)
describes the **buggy** variant 3.1. It is a faithful diagnosis of *what the
bug did*, not a property of the reversibility-index method. With the sign
corrected, the strongly-directional reactions now match their ΔG′° (and the
MSDB baseline direction): rxn00199 → `>`, rxn00256 → `<`, rxn05145 → `>`, so
they no longer flip and the catastrophic breakage is expected to largely
disappear.

**Status: RE-RUN COMPLETE (2026-06-16).** The one-line fix is applied in
`scripts/reversibility_lib.py` and the whole pipeline was regenerated
(`export_thermo_variants.py` → `build_site_data.py` →
`build_all_models_impact.py` → `build_site_data.py`, then this analysis).
All shipped artifacts now reflect the corrected sign. The corrected results
are summarized in **§Post-fix results** immediately below; the original
Q1–Q3 analysis that follows is retained, clearly marked, as the record of
what the bug did.

Note: variant **H4** (the "best-evidence" stack) also sets `ln_ri_by_rxn`, so
it inherited the same inversion; it was regenerated too (see Post-fix
results). The default `ReversibilityConfig()` does **not** set `ln_ri_by_rxn`,
so the byte-for-byte MSDB baseline parity is unaffected by the fix.

---

## Post-fix results (corrected sign) — what variant 3.1 *actually* does

With the sign corrected, **3.1 is a mild, well-behaved thermodynamic
constraint, not a model-breaker.** Reference is the same heuristic baseline
(4,000 / 5,683 models grow; mean biomass flux 63.74).

| metric | buggy 3.1 | **corrected 3.1** |
|---|---:|---:|
| reactions changed vs baseline | 3,256 | **1,316** |
| transition types | `>`→`<` (1,874), `=`→dir (1,316), `<`→`>` (66) | **only `=`→`>` (1,072) and `=`→`<` (244)** |
| growers lost (grew→not), all 5,683 | 2,377 | **44** |
| non-growers rescued (→grew) | 59 | **0** |
| grow-status flips, panel (of 100) | 73 | **0** |
| grow-status flips, all (of 5,683) | 2,436 | **44** |
| growers under 3.1 (of baseline 4,000) | 1,682 | **3,956** |
| mean biomass flux | 40.57 | **48.15** |

**The corrected variant never flips an already-directional reaction.** Every
one of the 1,316 changes is a reaction MSDB left **reversible** (`=`) that the
reversibility index now resolves to a confident direction — exactly the
intended use of the index. The three reactions that drove the buggy breakage
(rxn00199 the IDH decarboxylation, rxn00256 citrate synthase, rxn05145 the
phosphate-ATPase) **no longer change at all**: their index agrees with their
existing `>`/`<`/`<` directions, so they stay put.

**Q1 (corrected) — which reactions still break anything?** Only two, and they
break very few. Single-reaction knock-ins over the corrected change set:

| reaction | base→new | breaks | name |
|---|---|---:|---|
| rxn01476 | `=`→`>` | 43 / 3,177 | 6-phospho-D-glucono-1,5-lactone lactonohydrolase |
| rxn00251 | `=`→`<` | 1 / 2,723 | phosphate:oxaloacetate carboxy-lyase (PEPCK) |

These account for **44/44 = 100%** of the (now tiny) breakage. They are
*legitimate* directional calls on reactions MSDB had left reversible — a real
modeling choice to evaluate, not an artifact. (rxn01476 forcing the
gluconolactonase forward removes a small reverse flux a few dozen models had
been relying on; worth a look, but it is 1.1% of growers, not 59%.)

**Q2 (corrected) — break-rate across the grower population.** 44 / 4,000 =
**1.1%** of all heuristic-baseline growers; **0 / 100** on the panel. The
"75%" was entirely the sign bug. (The earlier selection-bias analysis is now
moot — there is barely any breakage to be biased about.)

**Q3 (corrected) — non-grower rescues.** **0.** The 59 "rescues" in the buggy
run were the rxn05145 reversal artifact; with the correct sign rxn05145 stays
forward and nothing is spuriously rescued.

**Caveat — flux still drops without growth loss.** Corrected 3.1 changes
biomass *flux* in 3,590 of 5,683 models and lowers the mean from 63.74 to
48.15, with **no** model gaining flux (max Δ = +0.00). That is expected and
correct: adding directionality to 1,316 previously-reversible reactions can
only shrink the feasible space, trimming flux, but it rarely removes the last
path to biomass. So 3.1 is a genuine, defensible constraint with modest
biological cost — the opposite of the catastrophe the buggy run reported.

---

## TL;DR *(historical — describes the BUGGY run; see Post-fix results above for the corrected picture)*

- **One reaction explains 90% of the breakage.** Forcing **rxn00199**
  (oxalosuccinate → CO₂ + 2-oxoglutarate; the decarboxylation half of
  isocitrate dehydrogenase) from forward (`>`) to reverse (`<`) — which is
  what 3.1's reversibility index does — single-handedly kills growth in
  **2,138 of the 2,377** broken growers. Adding one more reaction
  (**rxn05145**, a phosphate-transporting ATPase) reaches **96%**; five
  reactions reach **98.6%**.
- **The "75%" does NOT generalize.** 3.1 breaks **59.4%** of all 4,000
  heuristic-baseline growers (**57.0%** of the 3,461 on-disk growers), not
  ~75%. The panel's 73% is a selection artifact: the diversity-maximizing
  panel is enriched for *fragile* models that lack a bypass around the
  broken TCA step (conditional break-rate 76.7% on the panel vs 59.5%
  DB-wide).
- **Yes, 3.1 rescues 59 non-growers into growth.** Almost all (54/59) are
  rescued by the *same* rxn05145 flip (`>`→`<`) that breaks others — in
  those models the forward phosphate-ATPase was a futile ATP drain;
  reversing it frees ATP for biomass.

---

> **The four sections below (Method, Q1, Q2, Q3) describe the BUGGY run**
> (sign-inverted variant 3.1). They are kept verbatim as the diagnostic
> record of what the bug did and how it was traced to rxn00199. For what
> corrected 3.1 actually does, see **§Post-fix results** above.

## Method

Each variant's reversibility map is applied in memory to every model by
rewriting cobra reaction bounds (`>`→(0,1000), `<`→(−1000,0), `=`/`?`→
(−1000,1000)), then FBA is run on biomass. The reference is the
**heuristic baseline** (every model rebound to the byte-for-byte MSDB
cascade defaults): 4,000 / 5,683 models grow. Under full 3.1, 1,682 grow —
a net loss of 2,318, decomposing into **2,377 growers lost** and **59
non-growers rescued**.

To attribute the 2,377 losses to specific reactions, only **22 of the
3,256** reactions 3.1 changes actually appear in any core model (the other
3,234 are in MSDB tails absent from these reconstructions). For each of
those 22, we build a map = baseline with **only that one reaction** set to
its 3.1 direction, run FBA across the growers that contain it, and record
which lose growth. A greedy union of these single-reaction break-sets shows
how much of the total breakage a handful of reactions explain.

---

## Q1 — Which reactions break the models?

### Single-reaction knock-in break counts

Each row: baseline cascade with **only this reaction** flipped to its 3.1
direction. "breaks" = baseline-growers containing it that lose growth.

| Reaction | base→new | growers w/ rxn | breaks | break-rate | name |
|---|---|---:|---:|---:|---|
| **rxn00199** | **>→<** | 3,905 | **2,139** | **54.8%** | oxalosuccinate carboxy-lyase (2-oxoglutarate-forming) |
| rxn00256 | <→> | 3,906 | 2,038 | 52.2% | acetyl-CoA:oxaloacetate C-acetyltransferase (citrate synthase, written reverse) |
| rxn05145 | >→< | 3,906 | 453 | 11.6% | phosphate-transporting ATPase |
| rxn00251 | =→> | 2,723 | 194 | 7.1% | phosphate:oxaloacetate carboxy-lyase (PEP carboxykinase) |
| rxn10122 | =→< | 2,916 | 62 | 2.1% | NADH dehydrogenase (ubiquinone-8, 3.5 H⁺) |
| rxn00001 | >→< | 3,929 | 45 | 1.1% | diphosphate phosphohydrolase (pyrophosphatase) |
| rxn01476 | =→< | 3,177 | 42 | 1.3% | 6-phospho-D-glucono-1,5-lactone lactonohydrolase |
| rxn08527 | =→< | 2,951 | 1 | 0.0% | fumarate reductase |

(The other 14 candidate reactions break 0 growers on their own.)

### Greedy cover — a couple of reactions explain almost everything

| step | reaction | base→new | marginal breaks | cumulative |
|---|---|---|---:|---:|
| 1 | **rxn00199** | >→< | **2,138** | 2,138 |
| 2 | rxn05145 | >→< | 154 | 2,292 |
| 3 | rxn00256 | <→> | 40 | 2,332 |
| 4 | rxn00251 | =→> | 8 | 2,340 |
| 5 | rxn01476 | =→< | 3 | 2,343 |

**2,343 / 2,377 = 98.6%** of the breakage is explained by single-reaction
knock-ins; only **34** models break solely from a *combination* of changes.

Note rxn00256 breaks 2,038 models *on its own* but contributes only **+40**
marginal after rxn00199 — the two hit the **same** models. Both collapse the
oxidative TCA branch:

- **rxn00199** `>→<`: the reaction makes 2-oxoglutarate (α-ketoglutarate)
  by decarboxylating oxalosuccinate (ΔG′° = −3.33 kcal/mol, firmly forward).
  Forcing it reverse removes the only route to α-ketoglutarate → no
  glutamate/glutamine/proline/arginine → no biomass.
- **rxn00256** `<→>`: this is citrate synthase written in the
  hydrolysis direction; MSDB runs it `<` (Acetyl-CoA + OAA → citrate). 3.1
  flips it to `>` (citrate → Acetyl-CoA + OAA), abolishing citrate synthesis
  and again starving the TCA cycle of the carbon that becomes α-ketoglutarate.

Either reversal independently severs the oxidative TCA cycle at the
α-ketoglutarate-supplying step, so they break the same ~2,000 growers;
rxn00199 is the dominant attributable cause because slightly more growers
depend on it without an alternative.

**Bottom line for Q1:** 3.1's breakage is not diffuse. It is essentially a
**single thermodynamic mis-call on rxn00199** (the IDH decarboxylation step),
forced reverse by the reversibility index despite a favorable forward ΔG′°,
plus a smaller contribution from a phosphate-ATPase (rxn05145).

---

## Q2 — Does the panel's ~75% break-rate generalize to all ~3.5K growers?

**No.** Measured on the same heuristic-baseline → 3.1 basis as the panel
flip metric:

| Grower population | n | break under 3.1 | break-rate |
|---|---:|---:|---:|
| **Panel (100 descriptive)** | 100 | 73 | **73.0%** |
| On-disk growers (results.csv) ∩ heuristic growers | 3,460 | 1,972 | **57.0%** |
| All heuristic-baseline growers | 4,000 | 2,377 | **59.4%** |

So DB-wide the rate is ~**59%**, roughly **14 points below** the panel's 73%.

### Why the panel over-states it — selection bias, not containment

The panel is built by greedy max-coverage on reactions/metabolites plus
rare-reaction champions ([DIVERSE_SELECTION.md](DIVERSE_SELECTION.md)), which
favors specialized / minimal models. Decomposing the gap on the dominant
reaction rxn00199:

| | contains rxn00199 | break-rate \| contains rxn00199 |
|---|---:|---:|
| All baseline growers | 3,905 / 4,000 = 97.6% | 2,323 / 3,905 = **59.5%** |
| Panel baseline growers | 86 / 100 = 86.0% | 66 / 86 = **76.7%** |

The panel actually contains rxn00199 **less** often (86% vs 97.6%), so the
gap is **not** higher exposure. It is **higher conditional fragility**: when
a panel model carries rxn00199 it breaks 76.7% of the time vs 59.5% across
the database. The diversity-maximizing selection picks models that more
often **lack a bypass** around the broken α-ketoglutarate step (e.g. a
glutamate/2-oxoglutarate transporter, an alternative aminotransferase, or a
reductive-TCA route). ~40% of the general grower population has such a
bypass and survives; the panel under-represents them.

**Bottom line for Q2:** 3.1 would change growth in a *majority* of growers
(~59%), but the panel's 73–75% over-states the population effect by ~14–16
points due to its enrichment for fragile/minimal models.

---

## Q3 — Does 3.1 flip any non-growers into growth?

**Yes — 59 models that do not grow under the baseline cascade grow under 3.1.**

| rescuer reaction | base→new | rescues | name / definition |
|---|---|---:|---|
| **rxn05145** | **>→<** | **54 / 59** | phosphate-transporting ATPase: `H₂O + ATP + Pi[periplasm] → ADP + 2 Pi[cyt] + H⁺` |
| rxn00199 | >→< | 1 | (incidental) |
| (combination only) | — | 4 | — |

The dominant rescuer is the **same rxn05145 flip (`>`→`<`)** that *breaks*
453 growers in Q1 — the effect's sign depends on network context:

- In the 54 rescued models, rxn05145 running **forward** (baseline) is a
  **futile ATP-hydrolysis / phosphate cycle**: it burns ATP (and shuttles
  phosphate across the membrane) in a loop the solver is forced to carry,
  leaving insufficient ATP for biomass, so baseline FBA returns no growth.
  3.1 forces it **reverse**, closing the drain and freeing ATP → growth
  appears.
- In the 453 *broken* models, the forward direction is on a productive
  pathway, so reversing it removes needed flux.

**Bottom line for Q3:** 3.1 rescues 59 non-growers, almost entirely by
reversing a phosphate-transporting ATPase that was acting as a futile ATP
sink. Net across the database 3.1 is strongly negative (−2,377 + 59), but it
is not purely destructive.

---

## Implications for adopting 3.1 *(updated post-fix)*

- **The "breaks 75% of models" headline was a sign bug, now fixed.** The
  rxn00199 "audit" this section originally called for resolved to the
  one-line sign inversion in §Correction; rxn00199 was never really called
  reverse by the index — our code read the index backwards. With the fix,
  rxn00199/rxn00256/rxn05145 don't change, and 3.1 breaks **1.1%** of
  growers (44/4,000), not ~59–75%.
- **Corrected 3.1 is adoptable as a mild constraint.** It only adds
  directionality to reactions MSDB left reversible (1,316 `=`→`>`/`<`),
  costs growth in 44 models, and rescues none. The two reactions behind the
  residual 44 (rxn01476, rxn00251) are legitimate `=`→directional calls
  worth a curation look, not artifacts.
- **It does trim flux broadly without killing growth** (mean 63.74→48.15;
  flux changes in 3,590 models; no model gains flux). That is the expected
  signature of a thermodynamic directionality constraint and is the real
  cost/benefit to weigh — far smaller than the buggy run implied.
- **Process lesson:** the panel's buggy "73/100 flip" looked like a strong
  biological signal; it was a code bug amplified by selection bias. Always
  sanity-check a directional call against ΔG′° sign (the cascade's own
  MdeltaG rule is the reference) before trusting a large impact number.

---

## Reproduction

```bash
cd core_models_analysis/
python3 scripts/analyze_variant_31_breakage.py   # ~2 min, 64 workers
# -> results/variant_31_breakage.json  (per-reaction break-sets, greedy cover)
```

Q2 population break-rates and Q3 rescue attribution are derived from
`site/data/all_models_baseline_fba.json`,
`site/data/all_models_variant_fba__3.1.json`,
`site/data/all_models_rxnsets.json`, and `results/results.csv` (on-disk
grower set). Reaction names/ΔG′° are from the notebook-06 kbcache
`msdb_reactions_v1`.
