# Reaction-direction sources and what they do to models

*Eleven direction sources from ModelSEED Biochemistry v2.0.1, applied to 5,683 core models and 5,420 genome-scale models*

## Summary

The same direction source has opposite effects on the two model sets. On the core models most sources **increase** the number that grow; on the genome-scale models most **abolish** growth entirely. That is not a property of the sources.

Two things produce it. The genome-scale networks are roughly fifteen times larger by distinct reaction count, so a source rewrites a far greater share of the network. And within that larger surface sit a few reactions that are thermodynamically uphill but biologically obligatory. One of them, the thiamine phosphomethylpyrimidine kinase, is present in every genome-scale model tested and on its own abolishes growth in 99% of them. The thermodynamic estimate for it is correct; the constraint derived from it is still wrong.

The practical consequence is that **the core panel cannot validate a direction source**. It exercises 239 of the 56,012 reactions in the database and omits the pathways where thermodynamic and biological direction diverge.

## 1. What was measured

**Baseline.** The bounds written in the model files themselves: forward, reverse or reversible according to each reaction's stored lower and upper bound. Nothing is recomputed. This is the *on-disk* row throughout, and it is what the models were gap-filled and published with.

**Applying a source.** A reaction whose SEED annotation appears in the source is rewritten to (0, 1000), (-1000, 0) or (-1000, 1000). A reaction the source does not mention keeps its on-disk bounds. Media are applied afterwards, biomass is the objective, and a model counts as growing above 1e-6.

| Model set | Models | Reactions / model | Distinct reactions | Media |
|---|---:|---:|---:|---|
| Core (core_kegg2) | 5,683 | 128 | 239 | 347-compound complete |
| Genome-scale (ms2_gsm) | 5,420 | 1,097 | 3,634 | 20-compound glucose minimal |
| Genome-scale (auxotrophy) | 5,420 | 1,085 | 3,484 | 51-compound auxotrophy |

Reactions per model is the median count of distinct ModelSEED reactions. The genome-scale models come from the public KBase workspaces of the ModelSEED v2 manuscript and were gap-filled on the media shown.

## 2. The eleven sources

Three thermodynamic estimators, two LLM routes, three evidence tiers and three cumulative tiers. A grade belongs to a *reaction*, not to a source: the release ships one grade per reaction and keeps the per-source table internal. A tier is therefore the release's recommended direction restricted to reactions at that tier, chosen by a fixed precedence of eQuilibrator over dGPredictor over group contribution.

| Source | Calls | Forward | Reverse | Reversible | Directional |
|---|---:|---:|---:|---:|---:|
| Group contribution | 19,246 | 10,514 | 1,423 | 7,309 | 62% |
| dGPredictor | 12,177 | 6,863 | 2,016 | 3,298 | 73% |
| eQuilibrator | 21,218 | 9,129 | 1,065 | 11,024 | 48% |
| LLM council | 44,850 | 38,945 | 2,350 | 3,555 | 92% |
| Claude Opus 4.8 | 51,515 | 41,909 | 1,428 | 8,178 | 84% |
| gold | 3,089 | 1,687 | 96 | 1,306 | 58% |
| silver | 16,150 | 6,065 | 922 | 9,163 | 43% |
| bronze | 8,875 | 4,179 | 1,300 | 3,396 | 62% |
| gold + silver | 19,239 | 7,752 | 1,018 | 10,469 | 46% |
| gold + silver + bronze | 28,114 | 11,931 | 2,318 | 13,865 | 51% |

No set contains an undecided call. A direction map is applied by rewriting bounds, and an undecided value would be applied as fully reversible, so emitting it for a reaction that was graded but could not be called would *remove* a constraint in the name of evidence. Absent means no opinion, and the model keeps its own bounds. **Directional** is the share of a source's calls that fix a direction rather than allow both: the constraint pressure it applies.

---

## 3. What each source does to growth

![Figure 1](figures/direction_influence/fig1_growth_by_panel.png)

**Figure 1.** Share of models that grow under each direction source, against the bounds already in the model files. The three panels share a source ordering; each has its own baseline.

| Source | Core | vs base | GSM glucose | GSM auxotrophy |
|---|---:|---:|---:|---:|
| **on-disk baseline** | **3,461** | -- | **4,347** | **5,399** |
| Group contribution | 3,486 | +25 | 9 | 107 |
| dGPredictor | 3,806 | +345 | 0 | 0 |
| eQuilibrator | 3,948 | +487 | 0 | 515 |
| LLM council | 3,300 | -161 | 108 | 213 |
| Claude Opus 4.8 | 2,901 | -560 | 0 | 0 |
| gold | 3,339 | -122 | 2,502 | 4,618 |
| silver | 3,828 | +367 | 0 | 2,009 |
| bronze | 3,483 | +22 | 153 | 178 |
| gold + silver | 3,686 | +225 | 0 | 1,866 |
| gold + silver + bronze | 3,556 | +95 | 0 | 55 |

Models that grow, of 5,683 core and 5,420 genome-scale. Zero models failed to load in any run.

### On the core models, most sources help

eQuilibrator is the least directional source and gains the most models; Claude is among the most directional and loses the most. **Gold alone costs growth while silver alone gains it**, which inverts the evidence hierarchy: 122 models stop growing under the best-graded directions and 367 start growing under the weaker ones. Gold calls 58% of its reactions directional against silver's 43%, and a directional call removes a degree of freedom the model was using. The grade measures how well an energy is known, not how safe it is to impose.

### On the genome-scale models, most collapse

On glucose minimal media **six of the eleven sources leave zero growing models**. Every one returns an optimal solution with biomass flux zero, so these are genuine no-growth results rather than solver failures. The auxotrophy set starts from a near-complete baseline and separates the sources better: gold retains 85% of baseline, the silver-containing sets about a third, and six of eleven fall below 4%. The cumulative tiers decay monotonically there, 4,618 to 1,866 to 55.

"Nearly everything collapses" is fair for glucose minimal and an overstatement for auxotrophy. One detail runs the other way: among the models that do survive, mean growth flux *rises*, from 1.32 at baseline to 4.09 under gold and 6.52 under gold plus silver. Directional constraints remove futile cycling and channel flux, so the cost is concentrated entirely in the models that stop growing.

## 4. The grading, against the original thermodynamic sources

![Figure 2](figures/direction_influence/fig2_grading_vs_thermo.png)

**Figure 2.** Left: how often an evidence tier and an estimator agree, over the reactions where both have an opinion. Right: the share of each source's calls that fix a direction rather than allow both.

**The tiers are largely eQuilibrator re-expressed.** The full graded set agrees with eQuilibrator on 99.52% of the 20,868 reactions they share, bronze at 100.0% and silver at 99.9%. That follows directly from the recommendation precedence putting eQuilibrator first. A tier behaves as a confidence filter on one source rather than as a consensus of several, so comparing a tier against eQuilibrator is close to comparing eQuilibrator against itself.

**The estimators disagree with each other more than the tiers suggest.** dGPredictor agrees with the gold tier on only 51.08% of the 2,596 reactions they share, barely above chance among three options. Gold is the tier where the evidence is strongest, so this is disagreement about well-measured reactions, not about hard ones.

**The LLM routes are a different kind of claim.** The council agrees with eQuilibrator on 56.22% of 17,108 shared reactions; it reasons from the reaction rather than from an energy, which is the release's own argument for keeping it out of the grading. The council and Claude agree with each other 93.43% of the time over 42,599 reactions, but Claude is one of the council's members, so that is shared lineage rather than independent corroboration.

The right panel of Figure 2 shows what the agreement numbers hide: the tiers are markedly *less committal* than the estimators they are drawn from. Silver allows both directions on 57% of its calls against dGPredictor's 27%. That is most of why the silver tier is gentler on the models than the source it mostly agrees with.

## 5. Why the genome-scale models collapse

A total collapse is as likely to be a defect as a result, so it was diagnosed rather than reported.

**The conversion is faithful.** Direction round-trips exactly from the KBase originals: 327,872 reactions across a 300-model random sample, 0 mismatches. The models also grow normally on their own bounds. Only the override breaks them.

**It is not the medium.** The collapse reproduces on the 347-compound complete medium, not only the 20-compound minimal one.

**It localises to individual reactions.** Taking one genome-scale model under eQuilibrator directions, 229 of its reactions differ from disk and 51 of those are constraining, but applying them one at a time shows that only three individually abolish growth. Repeating that scan across 73 models that all grow at baseline, every one has at least one individually lethal constraint, and the distribution is extremely skewed.

![Figure 3](figures/direction_influence/fig4_mechanism.png)

**Figure 3.** Share of tested genome-scale models in which constraining one reaction alone abolishes growth. Rates are from a 73-model scan and shift by a few points between samples; the ranking of the first does not.

**One reaction explains most of it.** The thiamine phosphomethylpyrimidine kinase (rxn03108) is present in every genome-scale model and lethal in 72 of 73. eQuilibrator computes +15.59 kcal/mol for it and therefore calls it reverse. That number is not wrong: the phosphoryl transfer is uphill in isolation. But the enzyme is an ATP-driven kinase in thiamine biosynthesis and must run forward, so constraining it to the thermodynamically favoured direction removes thiamine and with it biomass. Both LLM routes call this reaction forward.

ATP sulfurylase (rxn00379) is the same failure with better evidence: graded **gold** from an openTECR measurement, +11 kcal/mol, called reverse by eQuilibrator, dGPredictor and Claude. In the cell it runs forward because pyrophosphatase removes the pyrophosphate product. The measurement is right and the direction call is right in isolation; the constraint is still wrong. It is a milder case, present in 67% of models and lethal in about one in eight of those.

Removing the top three reactions from the eQuilibrator map restores 137 of 381 growers in a 480-model sample, from zero. The collapse is not diffuse over-constraint; it is a handful of obligatory, cofactor-driven reactions.

**Why the core models escape.** They contain 239 distinct reactions against 3,634. Thiamine biosynthesis and sulfate assimilation are simply absent from a core carbon-metabolism model, so the reactions that kill the genome-scale models are never touched.

## 6. Does anything predict what a source will do?

![Figure 4](figures/direction_influence/fig3_coverage_effect.png)

**Figure 4.** Share of reaction occurrences a source changes, against the share of models that grow. Note the independent vertical scales: the core panel spans 51-70%, the genome-scale panel 0-80%.

Not reliably. Of the candidate summary statistics -- directional fraction, reactions called, occurrences changed, overrides per model -- none is significant on the core panel, and the best on the genome-scale panel reaches only Spearman rho of -0.65. The share of occurrences a source changes gives rho = -0.65 on glucose minimal, -0.45 on auxotrophy and -0.03 on the core panel, where it has essentially no ordering power at all.

Directional fraction does better on the core panel, rho = -0.68, but only when computed over the whole release. Restricted to the 239 reactions the core models actually contain -- the ones that can affect the result -- it falls to -0.28 and is not significant. The eleven sets are also nested rather than independent: gold sits inside gold plus silver, which sits inside gold plus silver plus bronze.

Two counterexamples make the point concretely. Bronze changes *less* of the network than gold, 3.3% against 3.7%, and retains far less growth. dGPredictor changes only 12.5% and leaves zero growers, the same as Claude at 31.5%. **What a source does depends on which reactions it constrains, not how many.**

## 7. What follows

**A direction source cannot be validated on the core models.** They exercise 239 of 56,012 reactions and omit the pathways where thermodynamic and biological direction come apart. Earlier conclusions drawn from that panel are not contradicted so much as shown to be untested.

**A grade is not a licence to constrain.** The grading scheme answers how well an energy is known, and answers it well: ATP sulfurylase is gold because openTECR measured it. Whether it is safe to fix that reaction's direction in a model is a different question, and the two come apart precisely at the cofactor-driven reactions. A tier used as an FBA constraint set needs a second filter for obligatory reactions.

**Reversibility is the safe default.** Every set that is more than about 60% directional damages the genome-scale models badly. That is a tendency rather than a rule, but it is the most robust one available.

**Two candidate fixes, neither tested here.** Exempt reactions whose stoichiometry contains ATP or pyrophosphate hydrolysis from directional constraint, which is the mechanism in every case examined; or intersect a tier with the LLM council and keep only reactions where both agree, trading coverage for safety.

## 8. Caveats

**Growth is the only readout.** A direction set that preserves grower count may still distort flux distributions; nothing here measures that.

**Energy-generating-cycle detection does not work on these models.** The existing probe closes every biomass reaction above ten metabolites and relies on a smaller one staying open as an ATP probe. The genome-scale models have exactly one biomass and no maintenance reaction, so it would report zero cycles for all 5,420 with no error.

**The lethal-reaction scan covers 73 models**, the growers among every 60th of the set. The ranking of the first reaction is unambiguous and replicated on an independent sample; the long tail is not fully enumerated.

**Claude Opus 4.8 is not part of the release.** It is a standalone single-model run; the council includes this model among its five roles, so the two are not independent.

**Correlations rest on eleven non-independent sets.** No coefficient here exceeds 0.7 in magnitude and most are not significant.

## Reproducing this

```bash
python3 scripts/build_direction_sets_v201.py
python3 beginPipeline --models core_kegg2 --workers 64 \
    --directions v201_gold ... --only growth,growth_diff \
    --out results/influence_v201/core_all
python3 scripts/analyze_direction_coverage.py
python3 scripts/regen_figures.py direction_influence
python3 scripts/build_influence_report.py
```

Results are under results/influence_v201/ and results/direction_sets_v201/, each run with a run.json recording the exact command. Every number in this document is read from those files at build time, and the PDF and Markdown renderings come from one content definition in scripts/build_influence_report.py.

