# Taxonomy-Aware Diverse Panel (100 models)

A taxonomy-aware re-selection of 100 representative growers, derived from the same 3,461 growers used by `scripts/select_diverse.py`. The original panel covered network sets well but missed entire phyla and big industrial/clinical genera; this panel adds an explicit taxonomic stratification pass and a phylogenetic-novelty pass that close those gaps.

## Methodology

Seven passes; each pick is tagged with the pass that selected it (`pass_origin` column) plus a per-pass reason.

1. **Phylum medoids** — one medoid per grower-bearing phylum (model with min sum-Jaccard reaction distance to its phylum-mates; sample capped at 50 per phylum). Guarantees every phylum is anchored.
2. **Reaction-coverage core** — greedy max-coverage on `seed.reaction` IDs, seeded with Pass 1. 12 picks.
3. **Taxonomic-novelty fill** — greedy farthest-point on lineage rank-distance (7=different superkingdom, 1=different species). Tie-break with Jaccard reaction distance. 12 picks; closes class/order gaps.
4. **Metabolite-coverage layer** — greedy max-coverage on cpd IDs, seeded with all prior picks. 8 picks.
5. **Constrained extremes** — 12 axis extremes + top-3 rare-reaction carriers, but skip any candidate whose class is already represented > 2.0× its expected share (walk down the ranking up to 20 before accepting). Cap 8.
6. **Hot-taxon medoids** — fills genus/family/order gaps that Pass 3 misses by sweeping the top grower-heavy missing taxa (threshold 20 growers) and adding their medoid. Cap 14; closes Pseudomonas, Burkholderia, etc. that sit inside already-anchored phyla.
7. **Farthest-point Jaccard fill** — fills the remaining slots by maximizing min reaction-set Jaccard distance to the panel. Pure diversity completion.

## Pass distribution

| Pass | Description | Picks |
|---|---|---|
| 1 | phylum medoids | 26 |
| 2 | reaction coverage | 8 |
| 3 | taxonomic novelty | 12 |
| 5 | constrained extremes | 8 |
| 6 | hot-taxon medoids | 14 |
| 7 | farthest-point Jaccard | 32 |

## Coverage achieved by the 100-model panel

- Reactions covered: **239 / 239** unique `seed.reaction` IDs (100.0%)
- Metabolites covered: **181 / 181** unique cpd IDs (100.0%)

## Taxonomic coverage

| Rank | distinct in panel | distinct among growers | coverage |
|---|---|---|---|
| superkingdom | 2 | 2 | 100.0% |
| phylum | 26 | 26 | 100.0% |
| class | 42 | 56 | 75.0% |
| order | 70 | 131 | 53.4% |
| family | 90 | 302 | 29.8% |
| genus | 97 | 870 | 11.1% |
| species | 100 | 2467 | 4.1% |

## Representativeness (Jaccard min-distance to nearest panel member)

| panel | median | mean | p90 | max |
|---|---|---|---|---|
| original | 0.192 | 0.181 | 0.255 | 0.304 |
| new | 0.149 | 0.147 | 0.237 | 0.309 |

## Gap closure vs original panel

### phylum (newly covered)

| phylum | grower count |
|---|---|
| Myxococcota | 21 |
| Rhodothermota | 3 |
| Gemmatimonadota | 3 |
| Calditrichota | 1 |
| Balneolota | 1 |
| Ignavibacteriota | 1 |

### class (newly covered)

| class | grower count |
|---|---|
| Myxococcia | 14 |
| Gemmatimonadia | 3 |
| Rhodothermia | 3 |
| Calditrichia | 1 |
| Ignavibacteria | 1 |
| Blastocatellia | 1 |
| Balneolia | 1 |
| Phycisphaerae | 1 |

### order (newly covered)

| order | grower count |
|---|---|
| Micrococcales | 122 |
| Rhodobacterales | 95 |
| Kitasatosporales | 70 |
| Vibrionales | 59 |
| Pseudonocardiales | 36 |
| Oceanospirillales | 35 |
| Synechococcales | 15 |
| Myxococcales | 14 |
| Deinococcales | 11 |
| Methylococcales | 5 |
| Gemmatimonadales | 3 |
| Planctomycetales | 3 |
| Rhodothermales | 3 |
| Geitlerinematales | 1 |
| Chloracidobacteriales | 1 |

### genus (newly covered)

| genus | grower count |
|---|---|
| Pseudomonas | 84 |
| Streptomyces | 63 |
| Vibrio | 49 |
| Salmonella | 47 |
| Acinetobacter | 43 |
| Xanthomonas | 37 |
| Francisella | 35 |
| Alteromonas | 20 |
| Microbacterium | 18 |
| Deinococcus | 11 |
| Gordonia | 11 |
| Flavobacterium | 11 |
| Amycolatopsis | 11 |
| Novosphingobium | 6 |
| Synechococcus | 5 |

## Files

- `results/selected_ids_tax.txt` — 100 model IDs
- `results/selected_models_tax.csv` — per-model metrics, lineage, and selection reason
- `results/selected_models_tax.json` — same data + coverage stats
- `scripts/select_diverse_tax.py` — reproducible selection
- `notebooks/08_TaxonomyAwareSelection.ipynb` — interactive walkthrough

## Panel members (in selection order)

| # | model_id | pass | phylum | class | order | genus | reason |
|---|---|---|---|---|---|---|---|
| 1 | `GCF_000746525.1` | 1 | Pseudomonadota | Gammaproteobacteria | Pseudomonadales | Pseudomonas | phylum-medoid: Pseudomonadota (sum_jaccard=11.409) |
| 2 | `GCF_000024785.1` | 1 | Actinomycetota | Actinomycetes | Mycobacteriales | Gordonia | phylum-medoid: Actinomycetota (sum_jaccard=10.912) |
| 3 | `GCF_004006435.1` | 1 | Bacillota | Bacilli | Caryophanales | Bacillus | phylum-medoid: Bacillota (sum_jaccard=12.441) |
| 4 | `GCF_003076455.1` | 1 | Bacteroidota | Flavobacteriia | Flavobacteriales | Flavobacterium | phylum-medoid: Bacteroidota (sum_jaccard=11.569) |
| 5 | `GCF_000304375.1` | 1 | Campylobacterota | Epsilonproteobacteria | Campylobacterales | Campylobacter | phylum-medoid: Campylobacterota (sum_jaccard=10.497) |
| 6 | `GCF_000016865.1` | 1 | Unknown | Unknown | Unknown | Unknown | phylum-medoid: Unknown (sum_jaccard=13.373) |
| 7 | `GCF_000317045.1` | 1 | Cyanobacteriota | Cyanophyceae | Geitlerinematales | Geitlerinema | phylum-medoid: Cyanobacteriota (sum_jaccard=8.379) |
| 8 | `GCF_001278055.1` | 1 | Thermodesulfobacteriota | Desulfuromonadia | Desulfuromonadales | Desulfuromonas | phylum-medoid: Thermodesulfobacteriota (sum_jaccard=6.84) |
| 9 | `GCF_009017495.1` | 1 | Deinococcota | Deinococci | Deinococcales | Deinococcus | phylum-medoid: Deinococcota (sum_jaccard=4.052) |
| 10 | `GCF_000280925.3` | 1 | Myxococcota | Myxococcia | Myxococcales | Pseudomyxococcus | phylum-medoid: Myxococcota (sum_jaccard=3.925) |
| 11 | `GCF_009720525.1` | 1 | Planctomycetota | Planctomycetia | Planctomycetales | Gimesia | phylum-medoid: Planctomycetota (sum_jaccard=3.317) |
| 12 | `GCF_000178955.2` | 1 | Acidobacteriota | Terriglobia | Terriglobales | Granulicella | phylum-medoid: Acidobacteriota (sum_jaccard=2.377) |
| 13 | `GCF_000317895.1` | 1 | Bdellovibrionota | Bdellovibrionia | Bdellovibrionales | Bdellovibrio | phylum-medoid: Bdellovibrionota (sum_jaccard=0.981) |
| 14 | `GCF_001747405.1` | 1 | Chlorobiota | Chlorobiia | Chlorobiales | Chlorobaculum | phylum-medoid: Chlorobiota (sum_jaccard=1.175) |
| 15 | `GCF_000018865.1` | 1 | Chloroflexota | Chloroflexia | Chloroflexales | Chloroflexus | phylum-medoid: Chloroflexota (sum_jaccard=0.803) |
| 16 | `GCF_000218625.1` | 1 | Deferribacterota | Deferribacteres | Deferribacterales | Flexistipes | phylum-medoid: Deferribacterota (sum_jaccard=0.885) |
| 17 | `GCF_000017605.1` | 1 | Spirochaetota | Leptospiria | Leptospirales | Leptospira | phylum-medoid: Spirochaetota (sum_jaccard=1.133) |
| 18 | `GCF_002310495.1` | 1 | Verrucomicrobiota | Opitutia | Opitutales | Nibricoccus | phylum-medoid: Verrucomicrobiota (sum_jaccard=0.916) |
| 19 | `GCF_000695095.2` | 1 | Gemmatimonadota | Gemmatimonadia | Gemmatimonadales | Gemmatimonas | phylum-medoid: Gemmatimonadota (sum_jaccard=0.314) |
| 20 | `GCF_001518995.2` | 1 | Rhodothermota | Rhodothermia | Rhodothermales | Unknown | phylum-medoid: Rhodothermota (sum_jaccard=0.446) |
| 21 | `GCF_003353065.1` | 1 | Balneolota | Balneolia | Balneolales | Cyclonatronum | phylum-medoid: Balneolota (sum_jaccard=0.0) |
| 22 | `GCF_001886815.1` | 1 | Calditrichota | Calditrichia | Calditrichales | Caldithrix | phylum-medoid: Calditrichota (sum_jaccard=0.0) |
| 23 | `GCF_000253035.1` | 1 | Chlamydiota | Chlamydiia | Parachlamydiales | Parachlamydia | phylum-medoid: Chlamydiota (sum_jaccard=0.0) |
| 24 | `GCF_000177635.2` | 1 | Chrysiogenota | Chrysiogenia | Chrysiogenales | Desulfurispirillum | phylum-medoid: Chrysiogenota (sum_jaccard=0.0) |
| 25 | `GCF_000279145.1` | 1 | Ignavibacteriota | Ignavibacteria | Ignavibacteriales | Melioribacter | phylum-medoid: Ignavibacteriota (sum_jaccard=0.0) |
| 26 | `GCF_000284315.1` | 1 | Nitrospirota | Nitrospiria | Nitrospirales | Leptospirillum | phylum-medoid: Nitrospirota (sum_jaccard=0.0) |
| 27 | `GCF_000018625.1` | 2 | Pseudomonadota | Gammaproteobacteria | Enterobacterales | Salmonella | reaction-coverage rank 1: adds 5 previously-uncovered reactions |
| 28 | `GCF_000025265.1` | 2 | Actinomycetota | Thermoleophilia | Solirubrobacterales | Conexibacter | reaction-coverage rank 2: adds 3 previously-uncovered reactions |
| 29 | `GCF_001266795.1` | 2 | Pseudomonadota | Gammaproteobacteria | Pseudomonadales | Marinobacter | reaction-coverage rank 3: adds 2 previously-uncovered reactions |
| 30 | `GCF_009688965.1` | 2 | Thermodesulfobacteriota | Desulfobacteria | Desulfobacterales | Desulfosarcina | reaction-coverage rank 4: adds 2 previously-uncovered reactions |
| 31 | `GCF_000008525.1` | 2 | Campylobacterota | Epsilonproteobacteria | Campylobacterales | Helicobacter | reaction-coverage rank 5: adds 1 previously-uncovered reactions |
| 32 | `GCF_000015745.1` | 2 | Bacillota | Bacilli | Caryophanales | Geobacillus | reaction-coverage rank 6: adds 1 previously-uncovered reactions |
| 33 | `GCF_000016745.1` | 2 | Thermodesulfobacteriota | Desulfuromonadia | Geobacterales | Geotalea | reaction-coverage rank 7: adds 1 previously-uncovered reactions |
| 34 | `GCF_000756615.1` | 2 | Bacillota | Bacilli | Caryophanales | Paenibacillus | reaction-coverage rank 8: adds 1 previously-uncovered reactions |
| 35 | `GCF_000015445.1` | 3 | Pseudomonadota | Alphaproteobacteria | Hyphomicrobiales | Bartonella | taxonomic-novelty: rank-distance 5 (Pseudomonadota/Alphaproteobacteria/Hyphomicrobiales/Bartonella) |
| 36 | `GCF_001189515.2` | 3 | Actinomycetota | Coriobacteriia | Coriobacteriales | Olsenella | taxonomic-novelty: rank-distance 5 (Actinomycetota/Coriobacteriia/Coriobacteriales/Olsenella) |
| 37 | `GCF_002082765.1` | 3 | Bacillota | Negativicutes | Veillonellales | Veillonella | taxonomic-novelty: rank-distance 5 (Bacillota/Negativicutes/Veillonellales/Veillonella) |
| 38 | `GCF_001688905.2` | 3 | Pseudomonadota | Betaproteobacteria | Burkholderiales | Unknown | taxonomic-novelty: rank-distance 5 (Pseudomonadota/Betaproteobacteria/Burkholderiales/Unknown) |
| 39 | `GCF_000194135.1` | 3 | Campylobacterota | Desulfurellia | Desulfurellales | Hippea | taxonomic-novelty: rank-distance 5 (Campylobacterota/Desulfurellia/Desulfurellales/Hippea) |
| 40 | `GCF_000507245.1` | 3 | Spirochaetota | Spirochaetia | Spirochaetales | Salinispira | taxonomic-novelty: rank-distance 5 (Spirochaetota/Spirochaetia/Spirochaetales/Salinispira) |
| 41 | `GCF_001688725.2` | 3 | Bacteroidota | Bacteroidia | Bacteroidales | Bacteroides | taxonomic-novelty: rank-distance 5 (Bacteroidota/Bacteroidia/Bacteroidales/Bacteroides) |
| 42 | `GCF_000021485.1` | 3 | Pseudomonadota | Acidithiobacillia | Acidithiobacillales | Acidithiobacillus | taxonomic-novelty: rank-distance 5 (Pseudomonadota/Acidithiobacillia/Acidithiobacillales/Acidithiobacillus) |
| 43 | `GCF_002005145.1` | 3 | Bacillota | Clostridia | Eubacteriales | Desulforamulus | taxonomic-novelty: rank-distance 5 (Bacillota/Clostridia/Eubacteriales/Desulforamulus) |
| 44 | `GCF_000226295.1` | 3 | Acidobacteriota | Blastocatellia | Chloracidobacteriales | Chloracidobacterium | taxonomic-novelty: rank-distance 5 (Acidobacteriota/Blastocatellia/Chloracidobacteriales/Chloracidobacterium) |
| 45 | `GCF_001659705.1` | 3 | Bacteroidota | Chitinophagia | Chitinophagales | Arachidicoccus | taxonomic-novelty: rank-distance 5 (Bacteroidota/Chitinophagia/Chitinophagales/Arachidicoccus) |
| 46 | `GCF_000025945.1` | 3 | Thermodesulfobacteriota | Desulfobulbia | Desulfobulbales | Desulfotalea | taxonomic-novelty: rank-distance 5 (Thermodesulfobacteriota/Desulfobulbia/Desulfobulbales/Desulfotalea) |
| 47 | `GCF_003261575.2` | 5 | Pseudomonadota | Gammaproteobacteria | Enterobacterales | Klebsiella | extreme: largest by reactions |
| 48 | `GCF_000283635.1` | 5 | Bacillota | Bacilli | Lactobacillales | Streptococcus | extreme: smallest by reactions |
| 49 | `GCF_000632475.1` | 5 | Pseudomonadota | Alphaproteobacteria | Rhodospirillales | Azospirillum | extreme: largest by metabolites |
| 50 | `GCF_000014405.1` | 5 | Bacillota | Bacilli | Lactobacillales | Lactobacillus | extreme: smallest by metabolites |
| 51 | `GCF_000599545.1` | 5 | Actinomycetota | Actinomycetes | Mycobacteriales | Rhodococcus | extreme: largest by genes |
| 52 | `GCF_000014425.1` | 5 | Bacillota | Bacilli | Lactobacillales | Lactobacillus | extreme: smallest by genes |
| 53 | `GCF_000021045.1` | 5 | Pseudomonadota | Gammaproteobacteria | Pseudomonadales | Azotobacter | extreme: highest growth flux |
| 54 | `GCF_000195855.1` | 5 | Actinomycetota | Actinomycetes | Mycobacteriales | Mycobacterium | extreme: lowest growth flux (still growing) |
| 55 | `GCF_003991855.1` | 6 | Actinomycetota | Actinomycetes | Micrococcales | Microbacterium | hot-taxon medoid: order=Micrococcales (122 growers) |
| 56 | `GCF_001562115.1` | 6 | Pseudomonadota | Gammaproteobacteria | Alteromonadales | Alteromonas | hot-taxon medoid: order=Alteromonadales (102 growers) |
| 57 | `GCF_006351965.1` | 6 | Pseudomonadota | Alphaproteobacteria | Rhodobacterales | Oceanicola | hot-taxon medoid: order=Rhodobacterales (95 growers) |
| 58 | `GCF_000816885.1` | 6 | Pseudomonadota | Gammaproteobacteria | Lysobacterales | Xanthomonas | hot-taxon medoid: order=Lysobacterales (86 growers) |
| 59 | `GCF_000013325.1` | 6 | Pseudomonadota | Alphaproteobacteria | Sphingomonadales | Novosphingobium | hot-taxon medoid: order=Sphingomonadales (81 growers) |
| 60 | `GCF_001278075.1` | 6 | Actinomycetota | Actinomycetes | Kitasatosporales | Streptomyces | hot-taxon medoid: order=Kitasatosporales (70 growers) |
| 61 | `GCF_000018445.1` | 6 | Pseudomonadota | Gammaproteobacteria | Moraxellales | Acinetobacter | hot-taxon medoid: order=Moraxellales (62 growers) |
| 62 | `GCF_000039765.1` | 6 | Pseudomonadota | Gammaproteobacteria | Vibrionales | Vibrio | hot-taxon medoid: order=Vibrionales (59 growers) |
| 63 | `GCF_000016105.1` | 6 | Pseudomonadota | Gammaproteobacteria | Thiotrichales | Francisella | hot-taxon medoid: order=Thiotrichales (44 growers) |
| 64 | `GCF_000010945.1` | 6 | Pseudomonadota | Alphaproteobacteria | Acetobacterales | Acetobacter | hot-taxon medoid: order=Acetobacterales (38 growers) |
| 65 | `GCF_000009105.1` | 6 | Pseudomonadota | Betaproteobacteria | Neisseriales | Neisseria | hot-taxon medoid: order=Neisseriales (36 growers) |
| 66 | `GCF_009429145.1` | 6 | Actinomycetota | Actinomycetes | Pseudonocardiales | Amycolatopsis | hot-taxon medoid: order=Pseudonocardiales (36 growers) |
| 67 | `GCF_009846525.1` | 6 | Pseudomonadota | Gammaproteobacteria | Oceanospirillales | Vreelandella | hot-taxon medoid: order=Oceanospirillales (35 growers) |
| 68 | `GCF_003260975.1` | 6 | Bacteroidota | Cytophagia | Cytophagales | Echinicola | hot-taxon medoid: order=Cytophagales (34 growers) |
| 69 | `GCF_000270285.1` | 7 | Actinomycetota | Coriobacteriia | Eggerthellales | Eggerthella | farthest-point: min Jaccard distance to selected = 0.413 |
| 70 | `GCF_002127965.1` | 7 | Pseudomonadota | Betaproteobacteria | Burkholderiales | Oxalobacter | farthest-point: min Jaccard distance to selected = 0.395 |
| 71 | `GCF_003073475.1` | 7 | Actinomycetota | Actinomycetes | Actinomycetales | Actinobaculum | farthest-point: min Jaccard distance to selected = 0.383 |
| 72 | `GCF_900186975.1` | 7 | Actinomycetota | Actinomycetes | Propionibacteriales | Cutibacterium | farthest-point: min Jaccard distance to selected = 0.382 |
| 73 | `GCF_008805035.1` | 7 | Pseudomonadota | Betaproteobacteria | Neisseriales | Eikenella | farthest-point: min Jaccard distance to selected = 0.379 |
| 74 | `GCF_000747315.1` | 7 | Actinomycetota | Actinomycetes | Mycobacteriales | Corynebacterium | farthest-point: min Jaccard distance to selected = 0.370 |
| 75 | `GCF_001888165.1` | 7 | Thermodesulfobacteriota | Desulfuromonadia | Desulfuromonadales | Syntrophotalea | farthest-point: min Jaccard distance to selected = 0.370 |
| 76 | `GCF_001042635.1` | 7 | Actinomycetota | Actinomycetes | Bifidobacteriales | Bifidobacterium | farthest-point: min Jaccard distance to selected = 0.369 |
| 77 | `GCF_000284095.1` | 7 | Bacillota | Negativicutes | Selenomonadales | Pseudoselenomonas | farthest-point: min Jaccard distance to selected = 0.364 |
| 78 | `GCF_009649955.1` | 7 | Bacillota | Clostridia | Eubacteriales | Heliorestis | farthest-point: min Jaccard distance to selected = 0.355 |
| 79 | `GCF_000828835.1` | 7 | Pseudomonadota | Gammaproteobacteria | Thiotrichales | Thioploca | farthest-point: min Jaccard distance to selected = 0.353 |
| 80 | `GCF_000298115.2` | 7 | Bacillota | Bacilli | Lactobacillales | Lentilactobacillus | farthest-point: min Jaccard distance to selected = 0.353 |
| 81 | `GCF_000260965.1` | 7 | Pseudomonadota | Gammaproteobacteria | Thiotrichales | Methylophaga | farthest-point: min Jaccard distance to selected = 0.352 |
| 82 | `GCF_000307165.1` | 7 | Bacillota | Bacilli | Caryophanales | Amphibacillus | farthest-point: min Jaccard distance to selected = 0.346 |
| 83 | `GCF_001262075.1` | 7 | Pseudomonadota | Betaproteobacteria | Burkholderiales | Ottowia | farthest-point: min Jaccard distance to selected = 0.344 |
| 84 | `GCF_000241025.1` | 7 | Pseudomonadota | Gammaproteobacteria | Pasteurellales | Aggregatibacter | farthest-point: min Jaccard distance to selected = 0.343 |
| 85 | `GCF_000599985.1` | 7 | Pseudomonadota | Gammaproteobacteria | Orbales | Gilliamella | farthest-point: min Jaccard distance to selected = 0.343 |
| 86 | `GCF_000014585.1` | 7 | Cyanobacteriota | Cyanophyceae | Synechococcales | Synechococcus | farthest-point: min Jaccard distance to selected = 0.341 |
| 87 | `GCF_000743945.1` | 7 | Pseudomonadota | Betaproteobacteria | Burkholderiales | Basilea | farthest-point: min Jaccard distance to selected = 0.333 |
| 88 | `GCF_000025545.1` | 7 | Pseudomonadota | Gammaproteobacteria | Chromatiales | Thioalkalivibrio | farthest-point: min Jaccard distance to selected = 0.331 |
| 89 | `GCF_000148405.1` | 7 | Pseudomonadota | Gammaproteobacteria | Lysobacterales | Xylella | farthest-point: min Jaccard distance to selected = 0.331 |
| 90 | `GCF_002162355.1` | 7 | Bacillota | Bacilli | Caryophanales | Tumebacillus | farthest-point: min Jaccard distance to selected = 0.329 |
| 91 | `GCF_000092365.1` | 7 | Actinomycetota | Actinomycetes | Actinomycetales | Arcanobacterium | farthest-point: min Jaccard distance to selected = 0.328 |
| 92 | `GCF_000012805.1` | 7 | Pseudomonadota | Gammaproteobacteria | Chromatiales | Nitrosococcus | farthest-point: min Jaccard distance to selected = 0.326 |
| 93 | `GCF_002302395.1` | 7 | Bacteroidota | Flavobacteriia | Flavobacteriales | Capnocytophaga | farthest-point: min Jaccard distance to selected = 0.325 |
| 94 | `GCF_000233715.2` | 7 | Bacillota | Clostridia | Eubacteriales | Desulfoscipio | farthest-point: min Jaccard distance to selected = 0.324 |
| 95 | `GCF_002443115.1` | 7 | Actinomycetota | Actinomycetes | Micrococcales | Dermabacter | farthest-point: min Jaccard distance to selected = 0.314 |
| 96 | `GCF_000284115.1` | 7 | Planctomycetota | Phycisphaerae | Phycisphaerales | Phycisphaera | farthest-point: min Jaccard distance to selected = 0.314 |
| 97 | `GCF_001644685.1` | 7 | Pseudomonadota | Gammaproteobacteria | Methylococcales | Methylomonas | farthest-point: min Jaccard distance to selected = 0.313 |
| 98 | `GCF_000020525.1` | 7 | Chlorobiota | Chlorobiia | Chlorobiales | Chloroherpeton | farthest-point: min Jaccard distance to selected = 0.312 |
| 99 | `GCF_003261295.1` | 7 | Pseudomonadota | Betaproteobacteria | Burkholderiales | Polynucleobacter | farthest-point: min Jaccard distance to selected = 0.310 |
| 100 | `GCF_000145255.1` | 7 | Pseudomonadota | Betaproteobacteria | Nitrosomonadales | Gallionella | farthest-point: min Jaccard distance to selected = 0.309 |
