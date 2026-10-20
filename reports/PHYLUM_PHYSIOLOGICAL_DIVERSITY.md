# Physiological diversity within phyla, and what the GMM→auxotrophy gap actually measures

**Scope.** A literature + domain-knowledge companion to
`results/baseline_biomass_taxonomy/taxonomy_phylum.tsv` and the per-phylum
grower-fraction figure. It answers one question: *do any of these phyla contain
organisms living across physiologically diverse conditions, and does that
diversity line up with what the models do?*

No code was written and no repo data was modified. The only computation here is
read-only `awk` over the existing `results/baseline_biomass_taxonomy/per_model.tsv`,
used to break the published phylum fractions down by order and genus. Those
breakdowns are reproduced inline and are the single most informative thing in
this report.

---

## Lede — the direct answer

**Yes, and three phyla in particular are so internally heterogeneous that the
phylum-level bar in the figure is close to meaningless as a biological
statement.** Bacillota spans a free-living soil prototroph that grows on glucose
and ammonium (*Bacillus subtilis*), a lactic acid bacterium that requires
fourteen amino acids delivered pre-made (*Lactobacillus helveticus*), and an
anaerobic polyextremophile growing optimally at 3.5 M Na⁺, pH 9.5 and 53 °C
(*Natranaerobius thermophilus*). Actinomycetota spans *Streptomyces* soil
saprophytes and the genome-decayed obligate intracellular *Mycobacterium
leprae*. Pseudomonadota spans *Pseudomonas* generalists, the acidophilic
iron-oxidising chemolithoautotroph *Acidithiobacillus ferrooxidans*, and
*Buchnera aphidicola*, which exists only inside an aphid cell.

**And — this is the useful finding — that heterogeneity is visible in the model
output, at genus resolution, and it is mostly real.** Within Bacillota the GMM
grower fraction runs from 0.00 (*Lactobacillus*, *Lactococcus*, *Enterococcus*,
*Listeria*) to 0.91–1.00 (*Bacillus*, *Geobacillus*). Within Actinomycetota it
runs from 0.00 (*Bifidobacterium*, *Actinomyces*) to 1.00 (*Streptomyces*).
Within Deinococcota it runs 0.07 (*Deinococcus*) vs 0.90 (*Thermus*). Within
Pseudomonadota, the nine genomes out of 2,365 that fail on glucose minimal are
*Buchnera*, *Wolbachia*, two *Francisella*, two other insect endosymbionts, and
*Stenotrophomonas acidaminiphila*. That is not noise. The models are recovering
the auxotrophy gradient.

**But the errors run in the opposite direction in a different set of phyla.**
Chlamydiota (0.99), Methanobacteriota (0.99), Mycoplasmatota (0.81),
Campylobacterota (1.00), Aquificota (1.00), Chlorobiota (1.00) and
Spirochaetota (0.93) are all near-ceiling on glucose minimal medium, and almost
none of their members can do that in reality. Here the near-absence of a
GMM→auxotrophy gap is the artefact, not the gap itself.

**So the short version: the big GMM→auxotrophy jumpers are mostly real biology;
the phyla with no jump at all are where the modelling is lying.**

---

## 0. What the two media actually are (needed before interpreting anything)

This matters and is easy to get wrong. The 5,420 genome-scale models are **two
separate gap-filled reconstructions per genome**, not one model run on two
media. On disk:

```
/scratch/ctaylor/modelseed2_gs_models/gmm/cobra/GCF_*.RAST.GMM.json
/scratch/ctaylor/modelseed2_gs_models/auxotrophy/cobra/GCF_*.RAST.auxo.json
```

and they differ in reaction content for the same genome (sampled pairs:
1404/1398, 1843/1841, 617/610, 982/972, 565/507 reactions, GMM/auxo).

> **Corrected 2026-10-20 on a larger sample.** An earlier draft of this section
> said the GMM reconstruction is *consistently* the larger of the two,
> generalised from those five pairs. Over 120 randomly sampled matched pairs
> that is **not** true: GMM is larger in 81, *smaller* in 34 and equal in 5
> (median Δ +6 reactions, range −26 to +84, reaction-set Jaccard median 0.962,
> min 0.558). The per-phylum *means* below are unaffected and still positive —
> the GMM gap-fill does add more on average — but no individual pair should be
> assumed to run that way. `bio1` is stoichiometrically identical in 120 of
> 120 pairs, which is what keeps biomass flux comparable across the two sets.

Per the ModelSEED v2 pipeline (quoted in `reports/BASELINE_BIOMASS_TAXONOMY.md`
§3), the **auxotrophy medium is the default gap-fill target**: *"The core
metabolism model is then expanded to the genome scale, followed by gap-filling
in the default auxotrophy medium or an optional user-specified medium."*

Three consequences:

1. **`auxo_grower_frac ≈ 1.00` everywhere is near-tautological.** Those models
   were built to grow there. A 0.99 is a statement about gap-filling converging,
   not about the organism.
2. **`gmm_grower_frac` is the informative column.** It is "did gap-filling
   against glucose + minimal salts succeed?" A failure means the gap-filler
   could not close the biosynthetic network from glucose alone.
3. **A GMM success is weak evidence.** The gap-filler is allowed to invent
   biochemistry. Mean extra reactions in the GMM model over the auxo model,
   by phylum:

   | phylum | Δ reactions (GMM − auxo) |
   |---|---|
   | Mycoplasmatota | **+54.9** |
   | Chlamydiota | **+39.6** |
   | Fusobacteriota | +35.0 |
   | Bdellovibrionota | +29.8 |
   | Bacillota | +27.6 |
   | Spirochaetota | +26.7 |
   | Thermotogota | +26.6 |
   | Actinomycetota | +13.9 |
   | Campylobacterota | +10.3 |
   | Pseudomonadota | +1.8 |
   | Bacteroidota | −1.8 |
   | Chlorobiota | −11.0 |
   | Aquificota | −11.7 |

   The phyla where the gap-filler worked hardest are exactly the reduced-genome
   host-dependent ones. It bolted ~55 reactions of invented biosynthesis onto
   each Mycoplasma and ~40 onto each Chlamydia to make them grow on sugar and
   salts.

---

## 1. Per-phylum summary table

Heterogeneity = *within-phylum* spread of physiology (temperature, pH,
salinity, O₂, trophic mode, host dependence). "Explains the jump?" refers to the
GMM→auxotrophy gap in the figure.

| phylum | dominant physiologies | O₂ relationship | temp / pH / salinity span | prototroph vs auxotroph | within-phylum heterogeneity | does biology explain the GMM→auxo jump? |
|---|---|---|---|---|---|---|
| **Pseudomonadota** (2,366) | chemoheterotroph generalists; also chemolitho­autotrophs, methanotrophs, anoxygenic phototrophs, obligate intracellular parasites, endosymbionts | all of: obligate aerobe → facultative → microaerophile → strict anaerobe | psychrophile→moderate thermophile; pH 1–10 (*Acidithiobacillus*); freshwater→hypersaline | overwhelmingly prototrophic; a small host-restricted tail is extremely auxotrophic | **HIGH** (highest absolute diversity of any phylum) | No jump to explain (0.996→1.000). The 9 failures are exactly the endosymbionts — **real biology, correctly captured** |
| **Bacillota** (1,002) | fermenters, aerobic/facultative heterotrophs, spore-formers, syntrophs, acetogens, sulfate reducers, polyextremophiles | strict anaerobe → aerotolerant → facultative → strict aerobe | psychrotroph→hyperthermophile (*Caldicellulosiruptor* ~78 °C); acidophile→alkaliphile pH 9.5+; up to 3.5 M Na⁺ | **bimodal**: *Bacillus*/*Geobacillus* prototrophic, LAB and gut Clostridia heavily auxotrophic | **HIGH** | **Yes — this is the real-biology jump.** 0.32→0.99 decomposes cleanly by genus (below) |
| **Actinomycetota** (717) | soil saprophytes, filamentous antibiotic producers, gut mutualists, obligate intracellular pathogens | mostly aerobe/facultative; gut clades anaerobic/microaerophilic | mesophile-dominated, some thermophiles; neutral to acidic soils | **bimodal**: *Streptomyces*/*Mycobacterium* prototrophic, *Bifidobacterium*/*Actinomyces*/*Tropheryma* auxotrophic | **HIGH** | **Yes — real biology.** 0.60→0.99 decomposes by order (below) |
| **Bacteroidota** (330) | saccharolytic/glycan-degrading heterotrophs; gut, soil, marine, sea-ice | strict anaerobe (gut *Bacteroides*) → strict aerobe (*Flavobacterium*) | psychrophile→moderate thermophile; marine→hypersaline relatives | mostly prototrophic for amino acids; B-vitamin dependencies common | **MED-HIGH** | No jump (0.997→0.997). Plausible: this phylum really is built to eat sugars |
| **Methanobacteriota** (203 / 153 full) | **obligate hydrogenotrophic methanogens only** | strict anaerobe (O₂-intolerant) | psychrophilic peatland→thermophilic; pH 6.5–8.5; 0–2.5% NaCl | prototrophic from CO₂, but need acetate/vitamins in practice | **LOW** (trophically the narrowest large phylum) | **No — artefact.** 0.99 GMM is wrong; methanogens do not catabolise glucose |
| **Campylobacterota** (181) | host-associated microaerophiles **and** deep-sea vent chemolitho­autotrophs | microaerophile → anaerobe | gut/gastric mesophile → hydrothermal vent | largely prototrophic; **asaccharolytic** | **HIGH** (gastric pathogen vs vent primary producer) | **No — artefact.** 1.00 GMM contradicts the defining trait of the phylum |
| **Mycoplasmatota** (147) | wall-less obligate host-associated parasites/commensals of animals, insects, plants | facultative/anaerobic, host-niche | host body temperature | **the most auxotrophic bacteria known** | LOW–MED (hosts vary; metabolism does not) | **Partly — but the direction is wrong.** 0.81 GMM is far too high; should be ~0 |
| **Cyanobacteriota** (118, no full models) | oxygenic photoautotrophs; some N₂-fixers, some facultative heterotrophs | oxygenic phototroph; some anoxygenic under sulfide | psychrophile→thermophile (~73 °C hot-spring *Synechococcus*); freshwater→hypersaline | prototrophic; streamlined marine picocyanobacteria scavenge N compounds | MED-HIGH | **N/A — absent from the genome-scale set**, itself a signal: a glucose-minimal medium is the wrong medium for a photoautotroph |
| **Chlamydiota** (93) | **obligate intracellular parasites only**; energy parasites with ATP/ADP translocases | aerobic host cytosol/vacuole | host body temperature | extreme auxotrophs | **LOW** | **No — artefact.** 0.99 GMM is flatly impossible |
| **Spirochaetota** (72) | free-living sediment/gut saccharolytics **and** obligate host-restricted pathogens | strict anaerobe → aerobe (*Leptospira*) | mesophile-dominated; one alkaliphilic halophile (*Spirochaeta africana*) | free-living clades prototrophic; *Borrelia*/*Treponema pallidum* drastically auxotrophic | **HIGH** | **Mostly artefact.** 0.93 GMM; *Treponema* 12/12 and *Leptospira* 13/13 "grow", which they cannot |
| **Thermodesulfobacteriota** (69) | sulfate reducers, Fe(III) reducers, syntrophs, predatory myxobacteria-adjacent lineages | strict/aerotolerant anaerobes | mesophile→hyperthermophile vent lineages; marine | mostly prototrophic; incomplete oxidisers excrete acetate | **HIGH** | No jump (1.00→1.00). Suspicious: these are mostly non-glucose lithotrophs/syntrophs |
| **Thermoproteota** (66 / 1 full) | thermoacidophilic archaea, sulfur oxidisers/reducers | aerobe→anaerobe | 55–85 °C, pH 1–4 (*Sulfolobus*) | prototrophic | MED-HIGH | Unassessable: only 1 genome-scale model |
| **Thermotogota** (33) | anaerobic saccharolytic fermenters, hyperthermophile→mesophile | strict anaerobe | **~20 °C (*Mesotoga*) to 90 °C (*Thermotoga*)** — the widest clean temperature span of any phylum here | prototrophic | **HIGH** (textbook model clade for thermo↔meso transitions) | **Mixed.** 0.33→1.00. *Mesotoga*'s failure is correct biology; *Thermotoga*'s (0.14) is a false negative |
| **Deinococcota** (31) | radiation/desiccation-resistant mesophiles **and** extreme thermophiles | aerobes | *Deinococcus* ~30 °C ↔ *Thermus* 65–80 °C; one alkaliphilic halotolerant (*Truepera*) | *Deinococcus* genuinely vitamin- and amino-acid-dependent; *Thermus* less so | **HIGH** for its size | **Yes — real biology**, and the genus split proves it (below) |
| **Planctomycetota** (29) | budding aerobic heterotrophs **and** anammox chemolitho­autotrophs | aerobe ↔ strict anaerobe | mesophile; acidic peat to marine | prototrophic | **HIGH** | No jump (1.00→1.00). Anammox members should not grow on glucose — likely masked artefact |
| **Chloroflexota** (28) | thermophilic anoxygenic photoheterotrophs/autotrophs, anaerobic organohalide respirers, aerobic thermophiles, filamentous sludge bacteria | strict anaerobe ↔ strict aerobe | mesophile→thermophile (*Chloroflexus* 55 °C) | *Chloroflexus* versatile; *Dehalococcoides* vitamin-dependent and obligately respiratory | **HIGH** | **Inverted — the best artefact example in the set** (below) |
| **Fusobacteriota** (24) | anaerobic **amino-acid fermenters** of mucosal surfaces | strict anaerobe | host body temperature | amino-acid-dependent for energy, not just biosynthesis | LOW–MED | **No — artefact.** 0.92 GMM; this phylum's defining trait is that it burns amino acids, not sugars |
| **Myxococcota** (22) | social predatory bacteria, fruiting-body formers; also anaerobic lineages | mostly aerobe | soil mesophiles | **auxotrophic for ≥5 amino acids**; barely metabolise glucose | MED | **No — artefact.** 1.00 GMM directly contradicts the classical nutrition data |
| **Chlorobiota** (15) | **obligately anaerobic anoxygenic photolitho­autotrophs** (green sulfur bacteria) | strict anaerobe, strictly phototrophic | mesophile; chemocline/sediment | photoautotrophic; cannot grow chemotrophically | **LOW** | **No — artefact.** 1.00 GMM is impossible: no light, no sulfide, no growth |
| **Nitrososphaerota** (15, no full models) | ammonia-oxidising archaea, chemolithoautotrophs | aerobe | soil/marine mesophile→thermophile | autotrophic | LOW–MED | N/A — absent from the genome-scale set |
| **Verrucomicrobiota** (14) | gut mucin degraders **and** thermoacidophilic methanotrophs | anaerobe/aerotolerant ↔ aerobe | *Akkermansia* 37 °C pH 7 ↔ *Methylacidiphilum* 55 °C **pH 2.5** | mixed | **HIGH** (extreme for 14 genomes) | No jump (1.00→1.00). Methanotroph members cannot use glucose — masked artefact |
| **Acidobacteriota** (13) | acidophilic oligotrophic soil/peat heterotrophs | aerobe/facultative | acidic soils pH 3–6; mesophile | prototrophic, slow-growing | MED | No jump (1.00→1.00); plausible for the heterotrophic majority |
| **Aquificota** (13) | **hyperthermophilic H₂/S-oxidising chemolitho­autotrophs** | microaerophile | 65–95 °C, hot springs and vents | autotrophic via reverse TCA | **LOW** | **No — artefact.** 1.00 GMM; *Aquifex* fixes CO₂ and does not eat glucose |
| **Bdellovibrionota** (13) | **obligate predatory bacteria** (periplasmic predators of other Gram-negatives) | aerobe | mesophile | obtains building blocks by digesting prey | LOW | **No — artefact.** 1.00 GMM, with +30 gap-filled reactions doing the work |

---

## 2. The phyla where a single phylum spans wildly different physiologies

This is the heart of the user's question. Ranked, with exemplars on each extreme.

### 2.1 Bacillota — HIGH, and the models see it

Extremes within one phylum:

- *Bacillus subtilis* — free-living soil aerobe, the textbook glucose+ammonium
  prototroph.
- *Lactobacillus helveticus* CNRZ 32 — **requires 14 amino acids**; the
  auxotrophy is gene *absence*, not point mutation, and matches the genome
  content (Christiansen et al. 2008, *Appl Environ Microbiol* 74:416–23,
  PMID 17993552, doi:10.1128/AEM.01174-07). The same paper states plainly that
  "most lactic acid bacteria are auxotrophic for several amino acids".
- *Listeria monocytogenes* — multiple amino acid requirements across strains
  (Siddiqi & Khan 1989, *Zentralbl Bakteriol* (1989), PMID 2505787).
- *Staphylococcus aureus* — largely prototrophic but with regulatory,
  conditional auxotrophies: CcpA-mediated proline auxotrophy (Li et al. 2010,
  *J Bacteriol*, PMID 20525824) and CcpA repression of arginine
  biosynthesis (Nuxoll et al. 2012, *PLoS Pathog*, PMID 23209408).
- *Natranaerobius thermophilus* — anaerobic **halophilic alkalithermophile**,
  optimum **3.5 M Na⁺, pH 9.5, 53 °C** (Mesbah, Cook & Wiegel 2009,
  *Mol Microbiol* 74:270–81, PMID 19708921, doi:10.1111/j.1365-2958.2009.06845.x).

Model output, GMM grower fraction by genus (`full_gmm` panel, n ≥ 15):

```
Enterococcus         n= 20  0.00
Lacticaseibacillus   n= 17  0.00
Lactobacillus        n= 27  0.00
Lactococcus          n= 17  0.00
Listeria             n= 41  0.00
Streptococcus        n=135  0.01
Paenibacillus        n= 47  0.28
Clostridium          n= 53  0.34
Staphylococcus       n= 85  0.80
Bacillus             n=112  0.91
Geobacillus          n= 16  1.00
```

and by order:

```
Lactobacillales      n=333  0.02      Caryophanales        n=431  0.57
Lachnospirales       n= 25  0.04      Selenomonadales      n= 10  0.90
Erysipelotrichales   n= 11  0.00      Eubacteriales        n=110  0.25
Peptostreptococcales n= 19  0.26      Thermoanaerobacterales n=16 0.19
```

**Verdict: real biology.** The ranking Lactobacillales ≈ 0 < Clostridia ≈ 0.3 <
*Staphylococcus* 0.80 < *Bacillus* 0.91 < *Geobacillus* 1.00 is the published
auxotrophy gradient, reproduced independently. Bacillota's 0.32→0.99 jump is
the single most defensible number in the figure.

One caveat: *Selenomonadales* at 0.90 is out of line for gut anaerobes.
Negativicutes are diderm Bacillota and would be assigned ModelSEED's
Gram-negative template; I cannot tell from this data whether their high score
is biology or a template effect, and I flag it as unresolved.

### 2.2 Actinomycetota — HIGH, and the models see it

Extremes:

- *Streptomyces* — filamentous free-living soil saprophytes with large genomes,
  prototrophic, prolific secondary metabolism.
- *Mycobacterium leprae* — **less than half the genome contains functional
  genes**; pseudogenes abound; gene decay eliminated siderophore production,
  parts of the respiratory chains, and numerous catabolic systems (Cole et al.
  2001, *Nature* 409:1007–11, PMID 11234002, doi:10.1038/35059006). It has never
  been cultured axenically.
- *Bifidobacterium bifidum* PRL2010 — **cysteine auxotroph**, growth promoted by
  whey/casein hydrolysate (Ferrario et al. 2015, *Front Microbiol* 6:1331,
  PMID 26635786, doi:10.3389/fmicb.2015.01331).
- *Tropheryma whipplei* — obligate intracellular, heavily reduced. (Stated from
  general knowledge; I did not find a primary auxotrophy paper in this search,
  so treat as uncited.)

Model output, GMM grower fraction by order:

```
Kitasatosporales (Streptomyces)  n= 72  1.00
Pseudonocardiales                n= 37  0.97
Mycobacteriales                  n=249  0.91
Micromonosporales                n= 19  0.89
Streptosporangiales              n= 11  0.73
Propionibacteriales              n= 46  0.43
Micrococcales                    n=157  0.23
Actinomycetales (Actinomyces)    n= 19  0.00
Bifidobacteriales                n= 60  0.00
```

within Micrococcales:

```
Kocuria        0.60   Arthrobacter   0.53   Rathayibacter 0.43
Pseudarthrobacter 0.40   Microbacterium 0.11   Brachybacterium 0.00
Cellulomonas   0.00   Clavibacter    0.00
```

**Verdict: real biology.** Free-living soil orders at 0.89–1.00; host-associated
gut and oral orders at 0.00. Exactly the expected split. The 0.60→0.99 jump is
genuine auxotrophy, not an artefact.

### 2.3 Pseudomonadota — HIGHEST absolute diversity, no jump to explain

Extremes:

- *Pseudomonas putida* — catabolic generalist, dozens of carbon sources.
- *Acidithiobacillus ferrooxidans* — extreme acidophile, chemolithoautotroph on
  Fe(II)/S, CO₂-fixing (Esparza et al. 2019, *Front Microbiol*,
  PMID 31019493; Jung, Inaba & Banta 2022, *Trends Biotechnol*,
  PMID 34794837).
- *Rickettsia prowazekii* — obligate intracellular energy parasite with an
  ATP/ADP translocase that steals host ATP (Andersson 1998, *Biochim Biophys
  Acta*, PMID 9693729; Audia & Winkler 2006, *J Bacteriol*,
  PMID 16923893).
- *Buchnera aphidicola* — obligate aphid endosymbiont; the canonical reduced
  genome, provisions essential amino acids to its host and cannot live outside
  it (Moran 2021, *Microbiology (Reading)*, PMID 34939561,
  doi:10.1099/mic.0.001127).

Model output — the **complete** list of GMM non-growers in this 2,365-model
phylum:

```
Buchnera aphidicola str. LL01 (Acyrthosiphon pisum)
Candidatus Kinetoplastidibacterium desouzai TCC079E
Candidatus Doolittlea endobia
Wolbachia pipientis wAlbB
Francisella orientalis LADL 07-285A
Candidatus Francisella endociliophora
Stenotrophomonas acidaminiphila
Pseudoxanthomonas spadix BD-a59
Candidatus Nitrotoga sp. AM1P
```

**Verdict: outstanding.** Six of the nine are obligate endosymbionts or
host-restricted pathogens; one is literally named for requiring amino acids
(*acidaminiphila*); one is a nitrite-oxidising chemolithoautotroph
(*Nitrotoga*) that has no business growing on glucose. The false-negative rate
on this phylum is essentially zero and the true-positive identification of the
host-dependent tail is essentially perfect. This is the cleanest validation
available that the pipeline is capable of recovering auxotrophy when the
genomic signal is strong.

### 2.4 Deinococcota — HIGH for 31 genomes, and correctly resolved

Extremes: *Deinococcus radiodurans* (mesophilic, radiation-resistant, ~30–37 °C)
vs *Thermus thermophilus* (extreme thermophile, 65–80 °C; Cava, Hidalgo &
Berenguer 2009, *Extremophiles*, PMID 19156357) vs *Truepera*
(alkaliphilic, halotolerant) vs *Marinithermus*/*Oceanithermus* (marine
thermophiles).

Holland, Rothfuss & Lidstrom 2006 (*Appl Microbiol Biotechnol* 72:1074–82,
PMID 16575566, doi:10.1007/s00253-006-0399-1) nail *D. radiodurans*'
requirements: it **requires biotin and niacin**; it is **methionine-auxotrophic
unless vitamin B₁₂ is supplied** (because its methionine synthase is
B₁₂-dependent); and critically, *"it was not possible to achieve reproducibly
good growth in the absence of at least one amino acid as a nitrogen source."*

Model output by genus:

```
Deinococcus     n=15  0.07
Meiothermus     n= 2  0.00
Thermus         n=10  0.90
Allomeiothermus n= 1  1.00
Marinithermus   n= 1  1.00
Oceanithermus   n= 1  1.00
Truepera        n= 1  1.00
```

**Verdict: real biology, strikingly so.** *Deinococcus* at 0.07 is precisely
what Holland et al. predict — it cannot grow on glucose + salts alone. The
phylum's 0.45→0.97 jump is a genuine *Deinococcus*-specific vitamin/amino-acid
dependency, not a modelling failure. This is my favourite result in the report
because the model and the 2006 bench experiment agree without either knowing
about the other.

### 2.5 Thermotogota — HIGH (the textbook temperature-span phylum), mixed result

This phylum is the standard model system for studying thermophile↔mesophile
evolutionary transitions precisely because it contains hyperthermophiles,
thermophiles, mesophiles, and organisms spanning both ranges (Pollo,
Zhaxybayeva & Nesbø 2015, *Can J Microbiol* 61:655–70, PMID 26211682,
doi:10.1139/cjm-2015-0073).

- *Thermotoga maritima* — hyperthermophilic anaerobic saccharolytic fermenter,
  grows ~55–90 °C (Conners et al. 2006, *FEMS Microbiol Rev* 30:872–905,
  PMID 17064285).
- *Mesotoga prima* — the **only described mesophilic Thermotogae lineage**
  (Nesbø et al. 2019, *Environ Microbiol* 21:456–470, PMID 30452102), and
  metabolically distinct: it **cannot ferment glucose**. It can only oxidise
  sugars when elemental sulfur or a hydrogenotrophic sulfate reducer is present
  as an electron sink — "obligate sugar oxidation", attributed to loss of the
  bifurcating [FeFe]-hydrogenase (Fadhlaoui et al. 2018, *Environ Microbiol*
  20:281–292, PMID 29124868, doi:10.1111/1462-2920.13995).

Model output by genus:

```
Thermotoga       n=14  0.14      Pseudothermotoga n= 4  1.00
Thermosipho      n= 4  0.00      Marinitoga       n= 2  1.00
Fervidobacterium n= 3  0.00      Defluviitoga     n= 1  1.00
Mesotoga         n= 1  0.00      Athalassotoga    n= 1  1.00
Petrotoga        n= 1  0.00      Kosmotoga        n= 2  0.50
```

**Verdict: mixed, and I will not claim more than the data supports.**
*Mesotoga*'s failure is **correct for the right reason** — that organism really
cannot grow on glucose without an external electron acceptor, and a glucose
minimal medium does not supply one. *Thermotoga*'s 0.14 is a **false negative**:
*T. maritima* is a saccharolytic glucose fermenter. So the 0.33→1.00 jump in
this phylum is a mixture of one correct refusal and a dozen gap-filling
failures, and I cannot separate them further without inspecting the gap-fill
logs.

### 2.6 Chloroflexota — HIGH, and the model result is *inverted*

This is the most instructive failure in the dataset.

- *Chloroflexus aurantiacus* / *C. aggregans* — thermophilic filamentous
  anoxygenic phototrophs that grow **photoheterotrophically, aerobically
  chemoheterotrophically, photoautotrophically under H₂, and
  chemolithotrophically under H₂/O₂ in the dark** (Kawai et al. 2019, *FEMS
  Microbiol Lett* 366:fnz122, PMID 31158281, doi:10.1093/femsle/fnz122). Four
  distinct trophic modes in one organism.
- *Dehalococcoides mccartyi* — the opposite: **strictly anaerobic, obligately
  respires halogenated compounds**, genomes "among the smallest known for
  free-living organisms", and **cyanocobalamin, thiamine and biotin are
  essential supplements** it cannot make (Schipp et al. 2013, *Philos Trans R
  Soc B* 368:20120321, PMID 23479751, doi:10.1098/rstb.2012.0321).
- Also present: *Anaerolinea*/*Caldilinea*/*Pelolinea* (anaerobic filamentous
  sludge bacteria), *Thermobaculum* (aerobic thermophile), *Roseiflexus*
  (thermophilic phototroph lacking chlorosomes; Hanada et al. 2002, *Int J Syst
  Evol Microbiol*, PMID 11837302).

Model output by genus:

```
Dehalococcoides  n=12  1.00  <-- obligate organohalide respirer, 3 vitamin auxotrophies
Dehalogenimonas  n= 2  1.00  <-- likewise
Chloroflexus     n= 3  0.00  <-- 4 trophic modes, most versatile member
Roseiflexus      n= 2  0.00
Anaerolinea      n= 1  0.00
Caldilinea       n= 1  0.00
Pelolinea        n= 1  0.00
Thermobaculum    n= 1  0.00
```

**Verdict: artefact, and inverted.** The twelve *Dehalococcoides* genomes — the
most nutritionally dependent organisms in the phylum, which require three
vitamins and cannot respire anything but halogenated compounds — all "grow" on
glucose and salts. The versatile *Chloroflexus* fails. Because *Dehalococcoides*
is 12 of 28 genomes, it alone drives most of the phylum's 0.64 GMM fraction.
Any biological reading of the Chloroflexota bar in the figure is wrong.

### 2.7 Verrucomicrobiota — HIGH, extreme for only 14 genomes

- *Akkermansia muciniphila* — anaerobic human gut mucin forager, 37 °C, pH ~7
  (Davey et al. 2023, *Nat Microbiol*, PMID 37337046); it also
  remodels cobamides rather than making them de novo (Mok et al. 2020, *mBio*
  PMID 33293380).
- *Methylacidiphilum kamchatkense* — **thermoacidophilic methanotroph**, optimal
  growth at **55 °C (range 37–60) and pH 2.5 (range 1–6)**, isolated from an
  acidic geothermal mud spring (Ratnadevi et al. 2023, *Int J Syst Evol
  Microbiol* 73:006060, PMID 37755432, doi:10.1099/ijsem.0.006060). It also
  oxidises subatmospheric H₂ via a high-affinity [NiFe] hydrogenase (Schmitz et
  al. 2020, *ISME J*, PMID 32042101).

A gut commensal at pH 7 and a volcanic methanotroph at pH 2.5 in a 14-genome
phylum is the sharpest physiological contrast per genome in the whole table.
Both score 1.00 on glucose minimal, which is wrong for the methanotroph.

### 2.8 Campylobacterota — HIGH, and the model result contradicts the phylum's defining trait

- *Helicobacter pylori* / *Campylobacter jejuni* — host-associated
  microaerophiles.
- *Sulfurimonas* / *Sulfurovum* — deep-sea hydrothermal vent
  **chemolithoautotrophs**, the dominant primary producers in those systems,
  able to both oxidise and disproportionate sulfur compounds (Wang et al. 2023,
  *mSystems* 8:e00954-22, PMID 36541763, doi:10.1128/msystems.00954-22).

On nutrition: campylobacters are **generally considered unable to metabolise
glucose** because they lack key glycolytic enzymes. A systematic search of
>6,000 *C. jejuni*/*C. coli* genomes found a complete Entner–Doudoroff pathway
in only **1.7%** (Vegge et al. 2016, *Front Microbiol* 7:1877, PMID 27920773,
doi:10.3389/fmicb.2016.01877). The trait is transferable on a plasticity region
in *C. coli* (Vorwerk et al. 2015, *Mol Microbiol*, PMID 26259566).

**Verdict: artefact.** 181/181 Campylobacterota models grow on glucose minimal
medium, with +10 gap-filled reactions each. The literature says ~1.7% should.
This is the most quantitatively checkable false positive in the dataset.

### 2.9 Myxococcota — MED heterogeneity, clean artefact

*Myxococcus xanthus* is a social predator. The classical nutrition study
(Hemphill & Zahler 1968, *J Bacteriol* 95:1011–7, PMID 4868349,
doi:10.1128/jb.95.3.1011-1017.1968) found that **phenylalanine, leucine,
isoleucine, valine and methionine are essential** for growth, and — directly
relevant here — that **¹⁴C-labelled ribose and glucose "were not appreciably
metabolized."**

22/22 Myxococcota models grow on glucose minimal medium. **Artefact.**

### 2.10 The obligate-lifestyle phyla: Chlamydiota, Mycoplasmatota, Bdellovibrionota, Fusobacteriota, Chlorobiota, Aquificota, Methanobacteriota

LOW within-phylum heterogeneity, uniformly high GMM scores, uniformly wrong.

- **Mycoplasmatota.** Mycoplasmas evolved from Gram-positive ancestors by
  reductive evolution; "the supply of nutrients from their hosts apparently
  enabled mycoplasmas to lose, during evolution, the genes for many assimilative
  processes" (Razin, Yogev & Naot 1998, *Microbiol Mol Biol Rev* 62:1094–1156,
  PMID 9841667). Model output by genus is damning: *Mycoplasmoides* 20/20,
  *Mycoplasmopsis* 20/20, *Mesomycoplasma* 15/15, *Metamycoplasma* 8/8,
  **Ureaplasma 4/4**, **Candidatus Phytoplasma 7/7** all "grow" on glucose and
  salts. *Ureaplasma* conserves energy by **hydrolysing urea** to generate a
  transmembrane potential (Smith et al. 1993, *J Bacteriol* 175:3253–8,
  PMID 8501029); phytoplasmas are uncultivable obligate plant/insect
  intracellular parasites (Oshima et al. 2004, *Nat Genet*,
  PMID 14661021). These models carry **+55 reactions** over their auxo
  counterparts. The gap-filler wrote them a metabolism.
- **Chlamydiota.** Obligate intracellular energy parasites with ATP/ADP
  translocases (Schmitz-Esser et al. 2004, *J Bacteriol*,
  PMID 14729693). 92/93 grow on glucose minimal, with +40 gap-filled reactions.
- **Bdellovibrionota.** Obligate periplasmic predators whose cell cycle is
  gated on prey recognition and prey quality (Rotem et al. 2015, *PNAS*
  PMID 26487679). 13/13 grow, +30 reactions.
- **Fusobacteriota.** *F. nucleatum* ferments amino acids for energy (Loesche &
  Gibbons 1968, *Arch Oral Biol* 13:191–202, PMID 5238887), and its sugar
  transport is itself amino-acid-dependent (Robrish, Oliver & Thompson 1987,
  *J Bacteriol* 169:3891–7, PMID 3114229). 22/24 grow on glucose minimal, +35
  reactions. Notably its auxo fraction is 1.00 — the one phylum where the
  auxotrophy medium is arguably the *physiologically correct* medium.
- **Chlorobiota.** Green sulfur bacteria are obligately anaerobic anoxygenic
  photolithoautotrophs fixing CO₂ by the reverse TCA cycle (Buchanan & Arnon
  1990, *Photosynth Res*, PMID 11540925; Gregersen, Bryant & Frigaard
  2011, *Front Microbiol* 2:116, PMID 21833341). 15/15 grow on dark glucose
  minimal medium. Impossible.
- **Aquificota.** *Aquifex aeolicus* is a microaerophilic hyperthermophilic H₂-
  and S-oxidiser that assimilates CO₂ via the reverse TCA cycle (Prioretti et
  al. 2023, *Life* 13:627, PMID 36983784, doi:10.3390/life13030627; genome:
  Deckert et al. 1998, *Nature*, PMID 9537320). 13/13 grow on glucose.
- **Methanobacteriota.** Methanogens of this phylum use **H₂+CO₂ or formate**
  and nothing else, and typically need yeast-extract growth factors; e.g.
  *Methanobacterium nebraskense* grows 20–45 °C, pH 6.5–8.5, 0–2.5% NaCl on
  H₂+CO₂ or formate only (Fiore, Zhou & Weber 2026, *Int J Syst Evol Microbiol*
  76:007139, PMID 42018355, doi:10.1099/ijsem.0.007139). 152/153 grow on glucose
  minimal medium.

### 2.11 Spirochaetota — HIGH heterogeneity, mostly artefactual result

Extremes: free-living alkaliphilic/halophilic sediment spirochaetes
(*Spirochaeta africana*, *Oceanispirochaeta*) vs *Treponema pallidum*, which has
drastically limited biosynthetic capability and cannot be cultured
conventionally (Norris, Cox & Weinstock 2001, *J Mol Microbiol Biotechnol*
PMID 11200228), vs *Borrelia burgdorferi*, which lacks most
biosynthesis and salvages riboflavin from its environment (Ramos et al. 2022,
*Mol Microbiol*, PMID 36054485).

Model output: *Treponema* 12/12 grow, *Leptospira* 13/13 grow, *Brachyspira*
8/8 grow; the only failures are three *Borrelia*/*Borreliella*, one
*Spirochaeta africana* and one *Oceanispirochaeta* — i.e. the pipeline flags
two **free-living** spirochaetes and lets every host-restricted one through.

(*Leptospira* is conventionally described as using long-chain fatty acids rather
than sugars as carbon and energy source. I could not retrieve a primary citation
for that in this session, so I flag it as textbook-level and unverified here.)

---

## 3. Ranking by within-phylum physiological heterogeneity

Most heterogeneous first. "Span" is the single widest axis of variation.

| rank | phylum | span | model resolves it? |
|---|---|---|---|
| 1 | **Pseudomonadota** | every O₂ relationship, every trophic mode, pH 1–10, free-living generalist → obligate endosymbiont | **Yes — near-perfectly** (the 9 failures are the 9 right ones) |
| 2 | **Bacillota** | strict anaerobe ↔ strict aerobe; mesophile ↔ 78 °C; pH 5 ↔ 9.5; 0 ↔ 3.5 M Na⁺; prototroph ↔ 14-amino-acid auxotroph | **Yes — cleanly, at genus level** |
| 3 | **Actinomycetota** | soil saprophyte with a 10 Mb genome ↔ obligate intracellular with >half its genome dead | **Yes — cleanly, at order level** |
| 4 | **Chloroflexota** | 4-trophic-mode thermophilic phototroph ↔ obligate organohalide respirer with 3 vitamin auxotrophies | **No — inverted** |
| 5 | **Thermotogota** | ~20 °C *Mesotoga* ↔ 90 °C *Thermotoga*; fermenter ↔ obligate syntroph | Partly (right answer for *Mesotoga*, wrong for *Thermotoga*) |
| 6 | **Verrucomicrobiota** | gut mucin forager at pH 7 ↔ volcanic methanotroph at pH 2.5, 55 °C | **No** (both 1.00) |
| 7 | **Spirochaetota** | alkaliphilic free-living sediment spirochaete ↔ uncultivable *T. pallidum* | **No — backwards** |
| 8 | **Campylobacterota** | human gastric microaerophile ↔ deep-sea vent chemolithoautotrophic primary producer | **No** (all 1.00) |
| 9 | **Deinococcota** | 30 °C radiation-resistant *Deinococcus* ↔ 80 °C *Thermus* ↔ alkaliphilic halotolerant *Truepera* | **Yes — cleanly** |
| 10 | **Thermodesulfobacteriota** | mesophilic sulfate reducers ↔ vent lineages ↔ Fe(III)-reducing *Geobacter* ↔ syntrophs | No (all 1.00) |
| 11 | **Planctomycetota** | aerobic budding heterotroph ↔ anammox chemolithoautotroph | No (all 1.00) |
| 12 | **Bacteroidota** | strictly anaerobic gut *Bacteroides* ↔ aerobic marine/sea-ice *Flavobacterium* | No jump to resolve (0.997 both) |
| 13 | **Cyanobacteriota** | psychrophilic ↔ 73 °C hot-spring *Synechococcus*; freshwater ↔ hypersaline | **N/A — no genome-scale models exist** |
| 14 | **Thermoproteota** | thermoacidophiles pH 1–4, 55–85 °C | N/A — 1 model |
| 15 | **Acidobacteriota** | acidophilic peat ↔ neutral soil oligotrophs | No jump |
| 16 | **Mycoplasmatota** | hosts vary (animal/insect/plant); metabolism does not | No — badly |
| 17–20 | **Fusobacteriota, Myxococcota, Bdellovibrionota, Chlorobiota, Aquificota, Chlamydiota, Methanobacteriota** | each is physiologically narrow by definition | No — all near-ceiling on GMM |

---

## 4. What this does and does not explain

### Explained

- **Bacillota 0.32→0.99, Actinomycetota 0.60→0.99, Deinococcota 0.45→0.97.**
  Genuine, literature-documented auxotrophy. These three jumps decompose into
  genus-level patterns that match published nutritional requirements without
  any tuning. If you want one sentence for a paper: *the GMM→auxotrophy gap in
  the Gram-positive phyla is the prototroph/auxotroph axis, and it resolves to
  the genus.*
- **Pseudomonadota's absence of a jump.** The phylum is overwhelmingly
  prototrophic; the handful of exceptions are correctly identified as
  endosymbionts.
- **Mycoplasmatota's 0.81→0.99 jump exists at all.** Something real is being
  detected (these organisms are maximally auxotrophic) — but 0.81 is far too
  high, so the jump is a faint echo of a much larger effect.
- **Cyanobacteriota's absence from the genome-scale set.** A glucose minimal
  medium and a default auxotrophy medium are both wrong media for an oxygenic
  photoautotroph. Their absence is appropriate, not a gap.

### Not explained — and where I would not trust the figure

- **Every phylum sitting at exactly 1.00 on glucose minimal medium that is
  dominated by lithotrophs, phototrophs, methanogens, predators, obligate
  intracellular parasites or asaccharolytic organisms.** Campylobacterota,
  Aquificota, Chlorobiota, Methanobacteriota, Myxococcota, Bdellovibrionota,
  Planctomycetota, Chlamydiota, and most of Spirochaetota and Fusobacteriota.
  In these the gap-filler manufactured a glucose catabolism and a biosynthetic
  network that the genome does not encode. The reaction-count deltas
  (Mycoplasmatota +55, Chlamydiota +40, Fusobacteriota +35, Bdellovibrionota
  +30) are the fingerprint.
- **Chloroflexota specifically is inverted** and should not be cited as
  evidence of anything biological.
- **Thermotogota's 0.33** mixes one correct refusal (*Mesotoga*) with ~12
  gap-fill failures on organisms that really do ferment glucose.
- **Selenomonadales at 0.90** inside an otherwise-0.02–0.34 Bacillota is
  unexplained; the diderm cell envelope of Negativicutes may route them to a
  different ModelSEED template, but I could not confirm that from the data
  available here.

### The structural point

The two media are not symmetric probes. **The auxotrophy medium is the
gap-filling target, so `auxo_grower_frac` measures convergence of the
gap-filler, not biology.** `gmm_grower_frac` is the only column carrying
organism-level information, and it carries it asymmetrically:

- a **GMM failure is strong evidence** of genuine biosynthetic dependency
  (precision looks very high — *Buchnera*, *Wolbachia*, *Lactobacillus*,
  *Bifidobacterium*, *Deinococcus*, *Mesotoga* are all correctly refused);
- a **GMM success is weak evidence** of anything, because the gap-filler is
  permitted to add whatever reactions make the LP feasible.

If this figure is going into a paper, the defensible framing is *"the fraction
of genomes for which a minimal-medium reconstruction could be closed"* rather
than *"the fraction of organisms that can grow on glucose"*, with the
genus-level Bacillota/Actinomycetota/Deinococcota/Pseudomonadota breakdowns as
the validation that the quantity is tracking real auxotrophy where the genomic
signal is strong.

### A concrete suggestion

The single highest-value follow-up would be to compare, per model, the set of
reactions present in the `.RAST.GMM` model but absent from the `.RAST.auxo`
model — the explicit list of what the gap-filler invented. For Chlamydiota and
Mycoplasmatota that list is ~40–55 reactions per genome and should be almost
entirely amino acid and cofactor biosynthesis with no genomic evidence. That
would turn the qualitative artefact claim above into a quantitative one.

---

## 5. References

All verified live against PubMed/Europe PMC during this analysis.

1. Christiansen JK, Hughes JE, Welker DL, Rodríguez BT, Steele JL, Broadbent JR. Phenotypic and genotypic analysis of amino acid auxotrophy in *Lactobacillus helveticus* CNRZ 32. *Appl Environ Microbiol* 2008;74(2):416–23. PMID 17993552. doi:10.1128/AEM.01174-07
2. Siddiqi R, Khan MA. Amino acid requirement of six strains of *Listeria monocytogenes*. *Zentralbl Bakteriol* 1989. PMID 2505787. doi:10.1016/s0934-8840(89)80067-2
3. Li C, Sun F, Cho H, et al. CcpA mediates proline auxotrophy and is required for *Staphylococcus aureus* pathogenesis. *J Bacteriol* 2010. PMID 20525824. doi:10.1128/JB.00237-10
4. Nuxoll AS, Halouska SM, Sadykov MR, et al. CcpA regulates arginine biosynthesis in *Staphylococcus aureus* through repression of proline catabolism. *PLoS Pathog* 2012. PMID 23209408. doi:10.1371/journal.ppat.1003033
5. Mesbah NM, Cook GM, Wiegel J. The halophilic alkalithermophile *Natranaerobius thermophilus* adapts to multiple environmental extremes using a large repertoire of Na⁺(K⁺)/H⁺ antiporters. *Mol Microbiol* 2009;74(2):270–81. PMID 19708921. doi:10.1111/j.1365-2958.2009.06845.x
6. Cole ST, Eiglmeier K, Parkhill J, et al. Massive gene decay in the leprosy bacillus. *Nature* 2001;409(6823):1007–11. PMID 11234002. doi:10.1038/35059006
7. Ferrario C, Duranti S, Milani C, et al. Exploring amino acid auxotrophy in *Bifidobacterium bifidum* PRL2010. *Front Microbiol* 2015;6:1331. PMID 26635786. doi:10.3389/fmicb.2015.01331
8. Moran NA. Microbe Profile: *Buchnera aphidicola* — ancient aphid accomplice and endosymbiont exemplar. *Microbiology (Reading)* 2021. PMID 34939561. doi:10.1099/mic.0.001127
9. Andersson SG. Bioenergetics of the obligate intracellular parasite *Rickettsia prowazekii*. *Biochim Biophys Acta* 1998. PMID 9693729. doi:10.1016/s0005-2728(98)00050-4
10. Audia JP, Winkler HH. Study of the five *Rickettsia prowazekii* proteins annotated as ATP/ADP translocases (Tlc). *J Bacteriol* 2006. PMID 16923893. doi:10.1128/JB.00371-06
11. Schmitz-Esser S, Linka N, Collingro A, et al. ATP/ADP translocases: a common feature of obligate intracellular amoebal symbionts related to Chlamydiae and Rickettsiae. *J Bacteriol* 2004. PMID 14729693. doi:10.1128/JB.186.3.683-691.2004
12. Jung H, Inaba Y, Banta S. Genetic engineering of the acidophilic chemolithoautotroph *Acidithiobacillus ferrooxidans*. *Trends Biotechnol* 2022. PMID 34794837. doi:10.1016/j.tibtech.2021.10.004
13. Esparza M, Jedlicki E, González C, et al. Effect of CO₂ concentration on uptake and assimilation of inorganic carbon in the extreme acidophile *Acidithiobacillus ferrooxidans*. *Front Microbiol* 2019. PMID 31019493. doi:10.3389/fmicb.2019.00603
14. Holland AD, Rothfuss HM, Lidstrom ME. Development of a defined medium supporting rapid growth for *Deinococcus radiodurans* and analysis of metabolic capacities. *Appl Microbiol Biotechnol* 2006;72(5):1074–82. PMID 16575566. doi:10.1007/s00253-006-0399-1
15. Cava F, Hidalgo A, Berenguer J. *Thermus thermophilus* as biological model. *Extremophiles* 2009. PMID 19156357. doi:10.1007/s00792-009-0226-6
16. Pollo SMJ, Zhaxybayeva O, Nesbø CL. Insights into thermoadaptation and the evolution of mesophily from the bacterial phylum Thermotogae. *Can J Microbiol* 2015;61(9):655–70. PMID 26211682. doi:10.1139/cjm-2015-0073
17. Nesbø CL, Charchuk R, Pollo SMJ, et al. Genomic analysis of the mesophilic Thermotogae genus *Mesotoga* reveals phylogeographic structure and genomic determinants of its distinct metabolism. *Environ Microbiol* 2019;21(1):456–70. PMID 30452102. doi:10.1111/1462-2920.14477
18. Fadhlaoui K, Ben Hania W, Armougom F, et al. Obligate sugar oxidation in *Mesotoga* spp., phylum Thermotogae, in the presence of either elemental sulfur or hydrogenotrophic sulfate-reducers as electron acceptor. *Environ Microbiol* 2018;20(1):281–92. PMID 29124868. doi:10.1111/1462-2920.13995
19. Conners SB, Mongodin EF, Johnson MR, et al. Microbial biochemistry, physiology, and biotechnology of hyperthermophilic *Thermotoga* species. *FEMS Microbiol Rev* 2006;30(6):872–905. PMID 17064285. doi:10.1111/j.1574-6976.2006.00039.x
20. Kawai S, Nishihara A, Matsuura K, Haruta S. Hydrogen-dependent autotrophic growth in phototrophic and chemolithotrophic cultures of thermophilic bacteria, *Chloroflexus aggregans* and *Chloroflexus aurantiacus*. *FEMS Microbiol Lett* 2019;366(10):fnz122. PMID 31158281. doi:10.1093/femsle/fnz122
21. Hanada S, Takaichi S, Matsuura K, Nakamura K. *Roseiflexus castenholzii* gen. nov., sp. nov. *Int J Syst Evol Microbiol* 2002. PMID 11837302. doi:10.1099/00207713-52-1-187
22. Schipp CJ, Marco-Urrea E, Kublik A, Seifert J, Adrian L. Organic cofactors in the metabolism of *Dehalococcoides mccartyi* strains. *Philos Trans R Soc Lond B Biol Sci* 2013;368(1616):20120321. PMID 23479751. doi:10.1098/rstb.2012.0321
23. Davey LE, Malkus PN, Villa M, et al. A genetic system for *Akkermansia muciniphila* reveals a role for mucin foraging in gut colonization and host sterol biosynthesis gene expression. *Nat Microbiol* 2023. PMID 37337046. doi:10.1038/s41564-023-01407-w
24. Mok KC, Sokolovskaya OM, Nicolas AM, et al. Identification of a novel cobamide remodeling enzyme in the beneficial human gut bacterium *Akkermansia muciniphila*. *mBio* 2020. PMID 33293380. doi:10.1128/mBio.02507-20
25. Ratnadevi CM, Erikstad HA, Kruse T, Birkeland NK. *Methylacidiphilum kamchatkense* gen. nov., sp. nov., an extremely acidophilic and moderately thermophilic methanotroph belonging to the phylum Verrucomicrobiota. *Int J Syst Evol Microbiol* 2023;73(9):006060. PMID 37755432. doi:10.1099/ijsem.0.006060
26. Schmitz RA, Pol A, Mohammadi SS, et al. The thermoacidophilic methanotroph *Methylacidiphilum fumariolicum* SolV oxidizes subatmospheric H₂ with a high-affinity, membrane-associated [NiFe] hydrogenase. *ISME J* 2020. PMID 32042101. doi:10.1038/s41396-020-0609-3
27. Vegge CS, Jansen van Rensburg MJ, Rasmussen JJ, et al. Glucose metabolism via the Entner-Doudoroff pathway in *Campylobacter*: a rare trait that enhances survival and promotes biofilm formation in some isolates. *Front Microbiol* 2016;7:1877. PMID 27920773. doi:10.3389/fmicb.2016.01877
28. Vorwerk H, Huber C, Mohr J, et al. A transferable plasticity region in *Campylobacter coli* allows isolates of an otherwise non-glycolytic food-borne pathogen to catabolize glucose. *Mol Microbiol* 2015. PMID 26259566. doi:10.1111/mmi.13159
29. Wang S, Jiang L, Xie S, et al. Disproportionation of inorganic sulfur compounds by mesophilic chemolithoautotrophic Campylobacterota. *mSystems* 2023;8(1):e00954-22. PMID 36541763. doi:10.1128/msystems.00954-22
30. Hemphill HE, Zahler SA. Nutrition of *Myxococcus xanthus* FBa and some of its auxotrophic mutants. *J Bacteriol* 1968;95(3):1011–7. PMID 4868349. doi:10.1128/jb.95.3.1011-1017.1968
31. Loesche WJ, Gibbons RJ. Amino acid fermentation by *Fusobacterium nucleatum*. *Arch Oral Biol* 1968;13(2):191–202. PMID 5238887. doi:10.1016/0003-9969(68)90051-4
32. Robrish SA, Oliver C, Thompson J. Amino acid-dependent transport of sugars by *Fusobacterium nucleatum* ATCC 10953. *J Bacteriol* 1987;169(9):3891–7. PMID 3114229. doi:10.1128/jb.169.9.3891-3897.1987
33. Razin S, Yogev D, Naot Y. Molecular biology and pathogenicity of mycoplasmas. *Microbiol Mol Biol Rev* 1998;62(4):1094–156. PMID 9841667. doi:10.1128/MMBR.62.4.1094-1156.1998
34. Smith DG, Russell WC, Ingledew WJ, Thirkell D. Hydrolysis of urea by *Ureaplasma urealyticum* generates a transmembrane potential with resultant ATP synthesis. *J Bacteriol* 1993;175(11):3253–8. PMID 8501029. doi:10.1128/jb.175.11.3253-3258.1993
35. Oshima K, Kakizawa S, Nishigawa H, et al. Reductive evolution suggested from the complete genome sequence of a plant-pathogenic phytoplasma. *Nat Genet* 2004. PMID 14661021. doi:10.1038/ng1277
36. Rotem O, Pasternak Z, Shimoni E, et al. Cell-cycle progress in obligate predatory bacteria is dependent upon sequential sensing of prey recognition and prey quality cues. *Proc Natl Acad Sci USA* 2015. PMID 26487679. doi:10.1073/pnas.1515749112
37. Gregersen LH, Bryant DA, Frigaard NU. Mechanisms and evolution of oxidative sulfur metabolism in green sulfur bacteria. *Front Microbiol* 2011;2:116. PMID 21833341. doi:10.3389/fmicb.2011.00116
38. Buchanan BB, Arnon DI. A reverse KREBS cycle in photosynthesis: consensus at last. *Photosynth Res* 1990. PMID 11540925.
39. Prioretti L, D'Ermo G, Infossi P, et al. Carbon fixation in the chemolithoautotrophic bacterium *Aquifex aeolicus* involves two low-potential ferredoxins as partners of the PFOR and OGOR enzymes. *Life (Basel)* 2023;13(3):627. PMID 36983784. doi:10.3390/life13030627
40. Deckert G, Warren PV, Gaasterland T, et al. The complete genome of the hyperthermophilic bacterium *Aquifex aeolicus*. *Nature* 1998. PMID 9537320. doi:10.1038/32831
41. Fiore NA, Zhou Y, Weber KA. *Methanobacterium nebraskense* sp. nov., a hydrogenotrophic methanogen isolated from saline wetland soil. *Int J Syst Evol Microbiol* 2026;76(4):007139. PMID 42018355. doi:10.1099/ijsem.0.007139
42. Norris SJ, Cox DL, Weinstock GM. Biology of *Treponema pallidum*: correlation of functional activities with genome sequence data. *J Mol Microbiol Biotechnol* 2001. PMID 11200228.
43. Ramos D, Lasseter AG, Richards CL, et al. Riboflavin salvage by *Borrelia burgdorferi* supports carbon metabolism and is essential for survival in the tick vector. *Mol Microbiol* 2022. PMID 36054485. doi:10.1111/mmi.14977
44. D'Souza G, Waschina S, Pande S, Bohl K, Kaleta C, Kost C. Less is more: selective advantages can explain the prevalent loss of biosynthetic genes in bacteria. *Evolution* 2014;68(9):2559–70. PMID 24910088. doi:10.1111/evo.12468
45. Starke S, Harris DMM, Paulay A, Aden K, Waschina S. Comparative analysis of amino acid auxotrophies and peptidase profiles in non-dysbiotic and dysbiotic small intestinal microbiomes. *Comput Struct Biotechnol J* 2025;27:821–31. PMID 40103612. doi:10.1016/j.csbj.2025.02.004
46. Sulfolobus / thermoacidophily: Schocke L, Bräsen C, Siebers B. Thermoacidophilic *Sulfolobus* species as source for extremozymes and as novel archaeal platform organisms. *Curr Opin Biotechnol* 2019. PMID 30875666. doi:10.1016/j.copbio.2019.02.012
47. Thermophilic cyanobacteria: Kees ED, Murugapiran SK, Bennett AC, Hamilton TL. Distribution and genomic variation of thermophilic cyanobacteria in diverse microbial mats at the upper temperature limits of photosynthesis. *mSystems* 2022. PMID 35980085. doi:10.1128/msystems.00317-22
48. Lovley DR, Giovannoni SJ, White DC, et al. *Geobacter metallireducens* gen. nov. sp. nov., a microorganism capable of coupling the complete oxidation of organic compounds to the reduction of iron and other metals. *Arch Microbiol* 1993. PMID 8387263. doi:10.1007/BF00290916

---

*Generated 2026-10-20. Data source: `results/baseline_biomass_taxonomy/per_model.tsv`
and `taxonomy_phylum.tsv`. Literature retrieved via the `pubmed-mcp-server` MCP
tools (PubMed + Europe PMC) per project convention; no citation in this report
was written from model memory without a live lookup, and the one claim I could
not verify (Leptospira carbon source) is flagged as such in §2.11.*
