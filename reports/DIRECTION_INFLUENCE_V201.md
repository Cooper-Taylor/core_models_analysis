# How reaction-direction sources influence the core and genome-scale models

Eleven direction sources from ModelSEED Biochemistry v2.0.1, applied to 5,683 core models and 5,420 genome-scale models, measured against the bounds the model files already carry.

**Headline.** The same direction source has opposite effects on the two model sets. On the core models most sources *increase* growth; on the genome-scale models nearly all of them abolish it.

That is not a property of the sources. Two things produce it. The genome-scale networks are 15 times larger, so a source rewrites 3-31% of all reaction occurrences rather than touching a handful. And within that larger surface sit a few reactions that are thermodynamically uphill but biologically obligatory: one of them, the thiamine phosphomethylpyrimidine kinase, is present in **every** genome-scale model and is individually lethal in 72 of 73 tested. eQuilibrator's +15.59 kcal/mol for it is correct, and constraining the reaction to that direction still removes thiamine and with it biomass.

The practical consequence: **the core panel cannot validate a direction source**, because the reactions where thermodynamic and biological direction diverge are not in it.

---

## 1. Scope and method

**Baseline.** The bounds written in the model files themselves: for each reaction, forward if `lower_bound >= 0 < upper_bound`, reverse if `upper_bound <= 0 > lower_bound`, reversible if both. Nothing is recomputed. This is the `on_disk` row of every table below, and it is what the models were gap-filled and published with.

**Applying a source.** A reaction whose `annotation['seed.reaction']` appears in the source is rewritten: `>` to (0, 1000), `<` to (-1000, 0), `=` to (-1000, 1000). A reaction the source does not mention keeps its on-disk bounds. Media are applied afterwards, biomass is the objective, and a model counts as growing above 1e-6.

**Model sets.**

| | models | reactions/model (median) | distinct MSDB reactions | media |
|---|---:|---:|---:|---|
| core (`core_kegg2`) | 5,683 | 128 | 239 | 347-compound complete |
| genome-scale (`ms2_gsm`) | 5,420 | 1,085 | 3,634 | 20-compound glucose minimal |
| genome-scale (`ms2_gsm_auxo`) | 5,420 | 1,085 | 3,484 | 51-compound auxotrophy |

Both genome-scale sets come from the ModelSEED v2 manuscript's public KBase workspaces and were gap-filled on the media shown. The core models carry a **15-fold smaller reaction surface**, which is the single most important number in this document.

## 2. The eleven sources

Three thermodynamic estimators, two LLM routes, three exclusive evidence tiers and three cumulative ones. A grade belongs to a reaction, not to a source, so a tier is the release's recommended direction restricted to reactions at that tier.

| source | calls | directional | forward-biased? |
|---|---:|---:|---|
| Group contribution | 19,246 | 62% | 55% forward |
| dGPredictor | 12,177 | 73% | 56% forward |
| eQuilibrator | 21,218 | 48% | 43% forward |
| LLM council | 44,850 | 92% | 87% forward |
| Claude Opus 4.8 | 51,515 | 84% | 81% forward |
| gold only | 3,089 | 58% | 55% forward |
| silver only | 16,150 | 43% | 38% forward |
| bronze only | 8,875 | 62% | 47% forward |
| gold | 3,089 | 58% | 55% forward |
| gold + silver | 19,239 | 46% | 40% forward |
| gold + silver + bronze | 28,114 | 51% | 42% forward |

`directional` is the share of calls that are `>` or `<` rather than `=`. It is the constraint pressure a source applies, and it ranges three-fold: silver calls 43% of its reactions directional, the LLM council 92%.

## 3. Coverage: how much of each network a source can reach

Growth effects are downstream of a simpler question. Of the reactions a model set actually contains, how many does a source have an opinion about, and how many of those would it *change*? A source that covers 45,000 database reactions is still inert on a model whose reactions it never mentions.

The core models contain **239 distinct** ModelSEED reactions; the genome-scale models contain **3,634**. The last column is the one that predicts damage: the share of all reaction occurrences, across every model, whose direction the source would rewrite away from what the file already says.

| source | % of core reactions called | % of core occurrences changed | % of GSM reactions called | % of GSM occurrences changed |
|---|---:|---:|---:|---:|
| Group contribution | 44% | **11.5%** | 42% | **5.1%** |
| dGPredictor | 48% | **17.6%** | 33% | **12.5%** |
| eQuilibrator | 71% | **24.3%** | 60% | **18.5%** |
| LLM council | 85% | **18.3%** | 93% | **31.1%** |
| Claude Opus 4.8 | 97% | **28.5%** | 98% | **31.5%** |
| gold only | 23% | **7.2%** | 12% | **3.7%** |
| silver only | 27% | **12.2%** | 37% | **11.4%** |
| bronze only | 20% | **7.0%** | 16% | **3.3%** |
| gold + silver | 50% | **19.4%** | 49% | **15.0%** |
| gold + silver + bronze | 71% | **26.4%** | 65% | **18.3%** |

This column is the **best single predictor tested** on the genome-scale models, but it is a weak one and it does not generalise. Against grower count across the eleven sources: Spearman rho = -0.65 (p = 0.03) on glucose minimal, -0.45 (p = 0.17) on auxotrophy, and -0.03 on the core models, where it has essentially no ordering power. A different column, the share of the model set's distinct reactions a source calls at all, beats it on Pearson in both genome-scale sets. With eleven non-independent sets (gold is nested inside gold+silver, which is nested inside gold+silver+bronze) none of this supports a law.

The broad pattern is still visible: gold changes 3.7% of reaction occurrences and retains the most growth of any source, Claude and the council change 31% and retain almost none, and the cumulative tiers climb 3.7% -> 15.0% -> 18.3% while their auxotrophy growth falls 4,618 -> 1,866 -> 55. But there are clear counterexamples: bronze changes *less* than gold (3.3% against 3.7%) and retains far less growth, and dGPredictor changes only 12.5% yet leaves zero growers, the same as Claude at 31.5%.

Note also that the two LLM routes reach **84.5-97.6%** of the reactions in these networks, against 32.5-70.7% for any thermodynamic source. The separation is wide on the genome-scale models (93.0-97.6% against 32.5-60.0%) and much narrower on the core models, where the council reaches 84.5% and eQuilibrator 70.7%. The LLMs have an opinion about nearly everything, which is their value and, applied as hard constraints, their danger.

## 4. Core models: most sources help

All 5,683 core models, complete media. Baseline 3,461/5,683 growers.

| source | growers | vs baseline | mean flux of growers | reactions overridden/model |
|---|---:|---:|---:|---:|
| **on-disk baseline** | **3,461** | — | 52.3904 | 0 |
| Group contribution | 3,486 | +25 | 55.1213 | 59.2 |
| dGPredictor | 3,806 | +345 | 62.1426 | 61.8 |
| eQuilibrator | 3,948 | +487 | 70.1681 | 88.6 |
| LLM council | 3,300 | -161 | 49.9081 | 99.0 |
| Claude Opus 4.8 | 2,901 | -560 | 51.1276 | 122.0 |
| gold only | 3,339 | -122 | 41.8413 | 33.6 |
| silver only | 3,828 | +367 | 63.489 | 35.0 |
| bronze only | 3,483 | +22 | 51.1321 | 24.4 |
| gold + silver | 3,686 | +225 | 58.6591 | 68.6 |
| gold + silver + bronze | 3,556 | +95 | 48.3827 | 93.1 |

The ordering broadly follows constraint pressure, in the direction you would expect: eQuilibrator is the least directional source (48% of its calls are `>` or `<`) and gains the most models (+487); Claude is among the most directional (84%) and loses the most (-560); silver, the least directional tier (43%), gains 367 while gold at 58% loses 122.

Measured rather than eyeballed, the relationship is real but **negative and fragile**: more directional calls means worse growth, Spearman rho = -0.68 (p = 0.02, n = 11) using each set's directional fraction over the whole release. Computed over only the 239 reactions the core models actually contain -- the ones that can affect the FBA -- it falls to rho = -0.28 (p = 0.41) and is not significant. The eleven sets are also nested rather than independent. Treat it as a tendency, not a predictor.

**Gold alone costs growth while silver alone gains it.** That inverts the evidence hierarchy and is worth stating plainly: on this model set, better-graded directions are more often committal, and commitment removes degrees of freedom the models were using. The grade measures how well an energy is known, not how safe it is to impose.

## 5. Genome-scale models: nearly everything collapses

All 5,420 models, each on the medium it was gap-filled for.

| source | glucose minimal | auxotrophy | overrides/model |
|---|---:|---:|---:|
| **on-disk baseline** | **4,347** | **5,399** | 0 |
| Group contribution | 9 | 107 | 441 |
| dGPredictor | 0 | 0 | 368 |
| eQuilibrator | 0 | 515 | 662 |
| LLM council | 108 | 213 | 1030 |
| Claude Opus 4.8 | 0 | 0 | 1072 |
| gold only | 2,502 | 4,618 | 183 |
| silver only | 0 | 2,009 | 375 |
| bronze only | 153 | 178 | 129 |
| gold + silver | 0 | 1,866 | 559 |
| gold + silver + bronze | 0 | 55 | 688 |

On glucose minimal media six of the eleven sources leave **zero** growing models (dGPredictor, eQuilibrator, Claude, silver only, gold+silver, gold+silver+bronze); every one returns an optimal LP with biomass flux 0, so these are genuine no-growth solutions, not solver failures. The auxotrophy set, which starts from a near-complete baseline of 5,399, separates them better and shows a clean monotonic decay across the cumulative tiers: gold retains 4,618 models, gold+silver 1,866, gold+silver+bronze 55.

"Nearly everything collapses" is fair for glucose minimal and an overstatement for auxotrophy. There, only two of eleven sets reach zero; gold retains 85.5% of baseline and the silver-containing sets about a third. Six of eleven fall below 4%.

**The survivors grow faster.** Mean flux among growing models rises from 1.32 at baseline to 4.09 under gold and 6.52 under gold+silver. Directional constraints remove futile cycling and channel flux, so the models that tolerate the constraints benefit from them. The cost is concentrated entirely in the models that stop growing.

## 6. Why the genome-scale models collapse

This needed diagnosing rather than reporting, because a total collapse is as likely to be a bug as a result. It is not a bug. Three checks establish that.

**The conversion is faithful.** Direction round-trips exactly from the KBase originals: 327,872 reactions across a 300-model random sample, zero mismatches (`results/influence_v201/conversion_roundtrip.json`). The models also grow normally on their own bounds, 4,347 and 5,399 growers. Only the override breaks them.

**It is not the medium.** The collapse reproduces on the 347-compound complete medium, not just the 20-compound minimal one.

**It localises to individual reactions.** Taking one genome-scale model under eQuilibrator directions: 229 of its reactions differ from disk, 51 of those are constraining, and applying them one at a time shows that **3 reactions individually abolish growth**. The other 48 are harmless.

Repeating that scan over 73 models that all grow at baseline: every one of them has at least one individually lethal constraint, and 48 distinct reactions account for all of them. The distribution is extremely skewed.

| reaction | lethal in | grade | eQ | dGP | GC | council | Claude | name |
|---|---:|---|---|---|---|---|---|---|
| `rxn03108` | 72/73 | silver | < | – | – | > | > | ATP:4-amino-2-methyl-5-phosphomethylpyrimidine phosphotransferase |
| `rxn00359` | 18/73 | silver | > | > | > | > | > | 3'-phosphoadenylyl-sulfate sulfohydrolase |
| `rxn00763` | 12/73 | silver | < | – | – | = | = | Glycerol:NAD+ oxidoreductase |
| `rxn00379` | 9/73 | gold | < | < | – | = | < | ATP:sulfate adenylyltransferase |
| `rxn08801` | 9/73 | silver | > | – | – | > | > | Lysophospholipase L1 |
| `rxn05119` | 4/73 | silver | < | – | – | > | > | L-aspartate:NAD+ oxidoreductase (deaminating) |

**One reaction explains most of it.** `rxn03108`, the thiamine phosphomethylpyrimidine kinase, is present in 100% of the genome-scale models and lethal in 72 of 73. eQuilibrator computes +15.59 kcal/mol for it and therefore calls it reverse. That number is not wrong: the phosphoryl transfer is uphill in isolation. But the enzyme is an ATP-driven kinase in thiamine biosynthesis and must run forward, so constraining it to the thermodynamically favoured direction removes thiamine and with it biomass.

`rxn00379`, ATP sulfurylase, is the same failure with better evidence: graded **gold** from an openTECR *measurement*, +11 kcal/mol, called reverse by eQuilibrator, dGPredictor and Claude. In the cell it runs forward because pyrophosphatase removes the PPi product. The measurement is right and the direction call is right in isolation; the constraint is still wrong. It is a milder case than `rxn03108`: present in 66.6% of the models and individually lethal in roughly one in eight of those that carry it, against `rxn03108`'s 100% presence and near-universal lethality.

Removing just the top three from the eQuilibrator map restores 137 of 381 growers in a 480-model sample, from zero. The collapse is not diffuse over-constraint. It is a handful of obligatory, cofactor-driven reactions.

**Why the core models escape.** They contain 239 distinct reactions against 3,634. Thiamine and sulfate assimilation are simply absent from a core carbon-metabolism model, so the reactions that kill the genome-scale models are never touched. The core panel does not so much disagree with the genome-scale result as fail to test it.

## 7. What the sources agree about

- The graded tiers are largely **eQuilibrator**: 99.52% agreement over 20,868 shared reactions, 100.0% for bronze and 99.9% for silver. The recommendation precedence puts eQuilibrator first, so a tier behaves as a confidence filter on one source rather than as a consensus of several.

- The LLMs and the thermodynamics **disagree**: council against eQuilibrator 56.22% over 17,108 reactions.

- The LLMs **agree with each other**, 93.43% over 42,599. Claude is a member of the council, so that is shared lineage, not independent corroboration.

- dGPredictor agrees with the gold tier only 51.08% of the time over 2,596 reactions, barely above chance among three options.

On the lethal reactions the LLMs are often the ones that are biologically right: both call `rxn03108` forward where eQuilibrator calls it reverse, and the council calls `rxn00379` reversible where the thermodynamics call it reverse. They reason from the reaction rather than from an energy, which is exactly the case where an energy misleads. That advantage does not survive aggregation, because both LLM sets are 84-92% directional overall and lose growth everywhere else.

## 8. What this means

**A direction source cannot be validated on the core models.** They exercise 239 of the 56,012 reactions in the database and omit the pathways where thermodynamic direction and biological direction come apart. Every conclusion in this project's earlier work rests on that panel, and the genome-scale result does not contradict it so much as show that it was not a test.

**A grade is not a licence to constrain.** The grading scheme answers "how well is this energy known", and it answers it well: `rxn00379` is gold because openTECR measured it. It does not answer "is it safe to fix this reaction's direction in a model", and those two questions come apart precisely at the reactions that matter, the cofactor-driven ones. A tier used as an FBA constraint set needs a second filter for obligatory reactions.

**Reversibility is the safe default and the sources differ mainly in how often they use it.** Every set that is more than about 60% directional damages the genome-scale models badly. But no single summary statistic of a direction set reliably predicts its growth effect: of the candidates tested -- directional fraction, reactions called, occurrences changed, overrides per model -- none is significant on the core panel, and the best on the genome-scale panel reaches only rho = -0.65. What a source does depends on *which* reactions it constrains, not how many.

**Two candidate fixes**, neither tested here: exempt reactions whose stoichiometry contains ATP/PPi hydrolysis from directional constraint, which is the mechanism in every case examined; or intersect a tier with the LLM council and keep only reactions where both agree, trading coverage for safety.

## 9. Reproducing this

```bash
# the eleven direction sets, from the v2.0.1 tag
python3 scripts/build_direction_sets_v201.py

# core models, all 5,683, complete media          (~40 s at 64 workers)
python3 beginPipeline --models core_kegg2 --workers 64 \
    --directions v201_gold --directions v201_gold_silver \
    --directions v201_gold_silver_bronze --directions v201_equilibrator \
    --directions v201_dgpredictor --directions v201_group_contribution \
    --directions v201_llm_council --directions claude_opus48 \
    --only growth,growth_diff --out results/influence_v201/core_all

# genome-scale, all 5,420, on the medium they were gap-filled for  (~5 min)
python3 beginPipeline --models ms2_gsm --workers 64 \
    --media /scratch/ctaylor/modelseed2_gs_models/gmm/media/Carbon-D-Glucose.cpd \
    --directions ... --only growth,growth_diff --out results/influence_v201/gsm_gmm

# which reactions each source can reach in each model set
python3 scripts/analyze_direction_coverage.py
```

Outputs are under `results/influence_v201/` (`core_all/`, `gsm_gmm/`, `gsm_auxo/`), each with `growth.json`, `growth_diff.json` and a `run.json` recording the exact command. The direction sets and their agreement matrix are in `results/direction_sets_v201/`, described in [DIRECTION_SETS_V201.md](DIRECTION_SETS_V201.md).

## 10. Caveats

- **Growth is the only readout.** A direction set that preserves grower count may still distort flux distributions, and nothing here measures that. Loop and FVA stages exist in the pipeline but were not run at this scale.

- **Energy-generating-cycle detection does not work on these models.** The existing `flux_loops` probe closes every biomass reaction with 10 or more metabolites and relies on a smaller one staying open as an ATP probe. The genome-scale models have exactly one biomass and no maintenance reaction, so it would report zero cycles for all 5,420 with no error. The model-set checker warns about this; the probe needs an injected ATP sink before any loop claim can be made here.

- **The lethal-reaction scan covers 73 models**: every 60th of the set is 91 models, of which 73 grow at baseline and are therefore testable. The ranking of `rxn03108` is unambiguous and replicated on an independent sample; the long tail is not fully enumerated, and the per-reaction counts below `rxn03108` shift by a few percentage points between samples.

- **Claude Opus 4.8 is not part of the release.** It is a standalone single-model run held in this repo. The council includes this model among its five roles, so the two are not independent.

- **Grade tiers are largely one source.** Comparing a tier against eQuilibrator is close to comparing eQuilibrator against itself.

- **`?` is never applied.** Sets contain only `>`, `<` and `=`; a reaction a source declines to call keeps its on-disk bounds rather than being opened wide.

