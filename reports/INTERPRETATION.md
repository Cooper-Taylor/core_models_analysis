# Core Models KEGG2 — Biological Interpretation

5,683 ModelSEED metabolic models tested for biomass production under the ModelSEEDDatabase complete media (`KBaseMedia.cpd`, 347 compounds).

## Headline numbers
- **3,461 (60.9%)** produced biomass (FBA optimum > 1e-6 on `bio1`)
- **2,222 (39.1%)** solved to optimal but with zero biomass — networks with gaps
- 0 errors, 0 missing biomass reactions, 0 infeasible solves

## Growers are larger and gene-richer
| metric | grower median | non-grower median | gap |
|---|---|---|---|
| metabolites | 158 | 125 | +33 |
| reactions   | 175 | 118 | +57 |
| genes       | 189 | 101 | +88 |
| exchanges   | 27  | 23  | +4  |

The non-growers are nearly half-size by gene count — these are draft models with too little metabolism to close the biomass equation. The exchange-count gap is small, so it's not a transport problem; it's a cytoplasmic-pathway-completeness problem.

## Grower biomass flux is unimodal and high
- median 52.3, mean 52.4, max 87.2
- only 33 growers fall below flux 10; almost all of the 3,461 grow vigorously

## What's blocking non-growers — top biomass precursors that cannot be made
| compound id | name | non-growers blocked | % |
|---|---|---|---|
| cpd00002 | ATP | 2222 | 100.0% |
| cpd00003 | NAD | 2222 | 100.0% |
| cpd00005 | NADPH | 2222 | 100.0% |
| cpd00022 | Acetyl-CoA | 2222 | 100.0% |
| cpd00024 | 2-Oxoglutarate | 1850 | 83.3% |
| cpd00032 | Oxaloacetate | 1714 | 77.1% |
| cpd00101 | ribose-5-phosphate | 1339 | 60.3% |
| cpd00236 | D-Erythrose4-phosphate | 1257 | 56.6% |
| cpd00079 | D-glucose-6-phosphate | 1179 | 53.1% |
| cpd00072 | D-fructose-6-phosphate | 1162 | 52.3% |
| cpd00102 | Glyceraldehyde3-phosphate | 1068 | 48.1% |
| cpd00169 | 3-Phosphoglycerate | 1050 | 47.3% |
| cpd00061 | Phosphoenolpyruvate | 1029 | 46.3% |
| cpd00020 | Pyruvate | 440 | 19.8% |

**Reading the table:** four cofactors are universally blocked across all 2,222 non-growers — ATP, NAD+, NADPH, and Acetyl-CoA. This is the signature of a core-energy-metabolism gap: without a closed loop for regenerating these carriers, the biomass reaction (which consumes them stoichiometrically) cannot run at any rate.

Beyond cofactors, the next-most-blocked precursors are TCA-cycle and central-carbon intermediates (2-oxoglutarate, oxaloacetate, 3-phosphoglycerate, etc.), confirming that the missing capability is central carbon catabolism and oxidative phosphorylation, not biosynthesis of exotic biomass building blocks.

Per-model, non-growers are missing a median of 9 (of 14) biomass precursors — most of the biomass equation, not a single bottleneck.

## Which reactions correlate with growth?
Reaction prevalence compared between the two cohorts. Δ = grower% − non-grower%.

### Top reactions enriched in GROWERS (+Δ)
| seed.reaction | name | grower % | non-grower % | Δ |
|---|---|---|---|---|
| rxn05561 | Transport of dicarboxylates, extracellular | 80.2% | 17.5% | +62.7% |
| rxn00330 | acetyl-CoA:glyoxylate C-acetyltransferase (thioester-hydrolysing, carboxymethyl-forming) | 68.5% | 7.8% | +60.7% |
| rxn00006 | hydrogen-peroxide:hydrogen-peroxide oxidoreductase | 91.8% | 35.6% | +56.2% |
| rxn00336 | isocitrate glyoxylate-lyase (succinate-forming) | 63.2% | 7.9% | +55.3% |
| rxn02167 | (S)-3-Hydroxybutanoyl-CoA hydro-lyase | 82.6% | 27.5% | +55.1% |
| rxn03250 | (S)-Hydroxyhexanoyl-CoA hydro-lyase | 82.1% | 27.2% | +54.9% |
| rxn03240 | (S)-3-Hydroxyhexadecanoyl-CoA hydro-lyase | 82.1% | 27.2% | +54.9% |
| rxn02376 | R03316 | 87.6% | 34.0% | +53.6% |
| rxn00441 | R00621 | 87.6% | 34.0% | +53.6% |
| rxn01872 | succinyl-CoA:enzyme N6-(dihydrolipoyl)lysine S-succinyltransferase | 87.3% | 34.1% | +53.2% |
| rxn00379 | ATP:sulfate adenylyltransferase | 86.5% | 34.5% | +52.0% |
| rxn03249 | (S)-hydroxyhexanoyl-CoA:NAD+ oxidoreductase | 82.8% | 31.5% | +51.2% |
| rxn06777 | (S)-3-Hydroxytetradecanoyl-CoA:NAD+ oxidoreductase | 82.8% | 31.5% | +51.2% |
| rxn09003 | Nitrate reductase (Menaquinol-8) (periplasm) | 76.4% | 25.8% | +50.5% |
| rxn05937 | Ferredoxin:NADP+ oxidoreductase | 67.7% | 17.2% | +50.4% |
| rxn00175 | Acetate:CoA ligase (AMP-forming) | 89.8% | 39.7% | +50.1% |
| rxn01480 | 3-Hydroxy-2-methylpropanoate:NAD+ oxidoreductase | 76.0% | 26.2% | +49.7% |
| rxn09001 | Nitrate reductase (Ubiquinol-8) (periplasm) | 75.1% | 26.6% | +48.6% |
| rxn10126 | succinate dehyrdogenase | 76.4% | 28.2% | +48.2% |
| rxn14427 | Nitrate reductase cytochrome-c type (2 protons translocated) | 76.9% | 28.8% | +48.1% |

### Top reactions enriched in NON-GROWERS (−Δ)
| seed.reaction | name | grower % | non-grower % | Δ |
|---|---|---|---|---|
| rxn05466_c | TRANS-RXN-173.ce | 12.1% | 48.0% | -35.9% |
| rxn05319_c | TRANS-RXNBWI-115401.ce.maizeexp.OH_OH | 53.9% | 80.1% | -26.2% |
| rxn05759 | hydrogen:ferredoxin oxidoreductase | 2.9% | 17.8% | -14.9% |
| rxn05488_c | acetate reversible transport via proton symport | 81.8% | 94.2% | -12.4% |
| rxn13974 | pyruvate:ferredoxin 2-oxidoreductase (CoA-acetylating) | 27.9% | 40.1% | -12.2% |
| rxn00782 | D-glyceraldehyde-3-phosphate:NADP+ oxidoreductase (phosphorylating) | 1.6% | 13.7% | -12.1% |
| rxn02527 | R03544 | 7.5% | 19.5% | -12.0% |
| rxn05939 | 2-oxoglutarate synthase (rev) | 26.6% | 37.5% | -11.0% |
| rxn05559_c | formate transport in via proton symport | 70.9% | 81.2% | -10.3% |
| rxn40505 |  | 1.6% | 11.8% | -10.1% |
| rxn05469_c | pyruvate reversible transport via proton symport | 90.2% | 100.0% | -9.8% |
| rxn03644 | D-arabino-hex-3-ulose-6-phosphate isomerase | 9.2% | 16.8% | -7.6% |
| rxn03643 | D-arabino-hex-3-ulose-6-phosphate formaldehyde-lyase (D-ribulose-5-phosphate-forming) | 8.4% | 15.8% | -7.4% |
| rxn15962 | carbon monoxide dehydrogenase/acetyl-CoA synthase (CODH/ACS) | 0.3% | 7.4% | -7.1% |
| rxn00151 | ATP:pyruvate,phosphate phosphotransferase | 29.5% | 36.2% | -6.7% |

**Interpretation:** the grower-enriched reactions cluster on a few axes:
- **TCA cycle**: citrate synthase (rxn00256), isocitrate dehydrogenase (rxn01387), oxalosuccinate → 2-oxoglutarate (rxn00199), 2-oxoglutarate dehydrogenase E2 (rxn01872), 2-oxoglutarate decarboxylation via TPP (rxn00441), and the glyoxylate-shunt enzymes malate synthase (rxn00330) and isocitrate lyase (rxn00336).
- **Acetyl-CoA activation**: acetate:CoA ligase (rxn00175).
- **Fatty-acid β-oxidation**: hydroxyacyl-CoA hydro-lyases / dehydrogenases (rxn02167, rxn03249/03250, rxn03240, rxn06777).
- **Respiration / electron transport**: succinate dehydrogenase (rxn10126), ferredoxin-NADP+ oxidoreductase (rxn05937), and dissimilatory nitrate reductases (rxn09001/09003, rxn14412/14414/14427).
- **Sulfur assimilation**: ATP-sulfurylase (rxn00379), needed because the biomass equation pulls sulfur-containing building blocks.
- **Dicarboxylate uptake**: rxn05561 (transports succinate/malate-like substrates into the cell).

Non-growers are *enriched* in proton-symport transport reactions (rxn05319, rxn05466, rxn05469 pyruvate, rxn05488 acetate, rxn05559 formate) — they have the import machinery but not the downstream catabolism, so the imported carbon never reaches the biomass equation. They're also enriched in alternative/fermentative-style oxidoreductases (pyruvate:ferredoxin rxn13974, 2-oxoglutarate synthase rxn05939, CODH/ACS rxn15962), which are characteristic of obligate-anaerobe drafts where ModelSEED's default biomass — calibrated around aerobic central metabolism — doesn't close.

## Bottom line
Growth capability in this collection is essentially binary and tracks **model completeness**, not media composition. The 60.9% of models that grow have the central-carbon + TCA + energy backbone needed to feed the ModelSEED biomass equation; the 39.1% that don't are missing roughly two-thirds of that backbone (median 9 of 14 precursors unreachable, 100% missing the four energy/redox cofactors). On complete media these non-growers would need gap-filling — not nutrient addition — to grow.
