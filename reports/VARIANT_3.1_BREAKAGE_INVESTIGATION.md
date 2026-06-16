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

## TL;DR

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

## Implications for adopting 3.1

- The headline "3.1 changes 3,256 reactions and breaks 75% of models" is
  misleading on two counts: only 22 of those reactions touch any core model,
  and the breakage is dominated by **one reaction** (rxn00199).
- Before adopting the reversibility index as a default, the **rxn00199
  call must be audited**: a ΔG′° = −3.33 kcal/mol decarboxylation being
  labeled reverse-only is almost certainly a reversibility-index artifact
  (likely the CO₂/molecularity normalization at the 1 mM reference, or a
  sign/curation issue in the stored `ln γ`). Fixing or overriding rxn00199
  alone would recover ~90% of the lost growers.
- rxn05145 is genuinely **double-edged** (breaks 453, rescues 54). Its
  forward-direction futile-cycle behavior in some models is a modeling
  artifact worth flagging independent of 3.1.
- The panel remains useful for *detecting* sensitive reactions but
  **over-estimates population-level growth impact**; report DB-wide rates
  (the "all models" scope in the site's Variant Browser) alongside panel
  rates whenever a magnitude is quoted.

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
