# 100 Diverse Growth Models — Selection Report

Selected from **3393 growing models** (of 5,683 total in `core_models_kegg2`).

## Methodology

Four passes, each tagged onto every pick so the role of the model in the panel is explicit:

1. **Reaction-coverage core** — greedy max-coverage on `seed.reaction` IDs across growers. Each pick contributes the most previously-uncovered reactions.
2. **Metabolite-coverage layer** — same greedy idea on compound IDs, seeded with the reaction core. Catches transporter / cofactor diversity that the reaction-set pass missed.
3. **Forced extremes** — smallest/largest by reactions, metabolites, genes; lowest/highest growth flux; widest/narrowest open-exchange repertoire; rare-reaction champions.
4. **Farthest-point Jaccard sampling** — fills the rest of the 100 slots by repeatedly adding the grower whose reaction-set Jaccard distance to the already-picked panel is maximal.

## Coverage achieved by the 100-model panel

- Reactions covered: **234 / 234** unique `seed.reaction` IDs (100.0%)
- Metabolites covered: **181 / 181** unique cpd IDs (100.0%)
- Rare reactions (<1% prevalence) hit by panel: see `n_rare_rxns` column in CSV.

## Files

- `selected_ids.txt` — one model ID per line, ready for `for id in $(cat selected_ids.txt); do ...; done`
- `selected_models.csv` — per-model metrics + selection reason
- `selected_models.json` — same data, plus coverage stats, for programmatic use

## Panel members (in the order they joined)

| # | model_id | reaction count | growth flux | selection reason |
|---|---|---|---|---|
| 1 | `GCF_003261575.2` | 221 | 40.634 | reaction-coverage rank 1: adds 187 previously-uncovered reactions |
| 2 | `GCF_000632475.1` | 212 | 61.502 | reaction-coverage rank 2: adds 20 previously-uncovered reactions |
| 3 | `GCF_000164985.3` | 187 | 32.852 | reaction-coverage rank 3: adds 9 previously-uncovered reactions |
| 4 | `GCF_009688965.1` | 179 | 19.98 | reaction-coverage rank 4: adds 7 previously-uncovered reactions |
| 5 | `GCF_003441595.1` | 182 | 27.024 | reaction-coverage rank 5: adds 3 previously-uncovered reactions |
| 6 | `GCF_000010625.1` | 144 | 16.908 | reaction-coverage rank 6: adds 2 previously-uncovered reactions |
| 7 | `GCF_000011305.1` | 161 | 19.464 | reaction-coverage rank 7: adds 1 previously-uncovered reactions |
| 8 | `GCF_000015745.1` | 185 | 20.097 | reaction-coverage rank 8: adds 1 previously-uncovered reactions |
| 9 | `GCF_000016745.1` | 151 | 18.596 | reaction-coverage rank 9: adds 1 previously-uncovered reactions |
| 10 | `GCF_000025265.1` | 170 | 33.816 | reaction-coverage rank 10: adds 1 previously-uncovered reactions |
| 11 | `GCF_000200595.1` | 133 | 20.317 | reaction-coverage rank 11: adds 1 previously-uncovered reactions |
| 12 | `GCF_000756615.1` | 178 | 29.647 | reaction-coverage rank 12: adds 1 previously-uncovered reactions |
| 13 | `GCF_000283635.1` | 99 | 14.275 | extreme: smallest by reactions |
| 14 | `GCF_000014405.1` | 100 | 8.681 | extreme: smallest by metabolites |
| 15 | `GCF_000599545.1` | 197 | 57.355 | extreme: largest by genes |
| 16 | `GCF_000014425.1` | 106 | 8.657 | extreme: smallest by genes |
| 17 | `GCF_002591335.1` | 191 | 62.063 | extreme: highest growth flux |
| 18 | `GCF_000746585.1` | 112 | 1.674 | extreme: lowest growth flux (still growing) |
| 19 | `GCF_000008165.1` | 194 | 40.634 | extreme: widest open exchange repertoire |
| 20 | `GCF_000024425.1` | 137 | 2.902 | extreme: narrowest open exchange repertoire |
| 21 | `GCF_000018405.1` | 140 | 17.802 | extreme: carries many rare reactions (<1% of growers) |
| 22 | `GCF_000233715.2` | 156 | 26.203 | extreme: carries many rare reactions (<1% of growers) |
| 23 | `GCF_001688725.2` | 122 | 17.73 | farthest-point: min Jaccard distance to already-selected = 0.454 |
| 24 | `GCF_000270285.1` | 129 | 15.253 | farthest-point: min Jaccard distance to already-selected = 0.434 |
| 25 | `GCF_008805035.1` | 112 | 14.207 | farthest-point: min Jaccard distance to already-selected = 0.431 |
| 26 | `GCF_000237205.1` | 130 | 30.736 | farthest-point: min Jaccard distance to already-selected = 0.431 |
| 27 | `GCF_002082765.1` | 125 | 21.991 | farthest-point: min Jaccard distance to already-selected = 0.419 |
| 28 | `GCF_003352785.1` | 125 | 4.853 | farthest-point: min Jaccard distance to already-selected = 0.405 |
| 29 | `GCF_900186975.1` | 130 | 5.845 | farthest-point: min Jaccard distance to already-selected = 0.402 |
| 30 | `GCF_001688905.2` | 134 | 18.361 | farthest-point: min Jaccard distance to already-selected = 0.393 |
| 31 | `GCF_001262075.1` | 139 | 9.718 | farthest-point: min Jaccard distance to already-selected = 0.390 |
| 32 | `GCF_001020955.1` | 132 | 6.617 | farthest-point: min Jaccard distance to already-selected = 0.390 |
| 33 | `GCF_000284095.1` | 147 | 21.144 | farthest-point: min Jaccard distance to already-selected = 0.389 |
| 34 | `GCF_000020645.1` | 129 | 9.755 | farthest-point: min Jaccard distance to already-selected = 0.388 |
| 35 | `GCF_900232105.1` | 144 | 16.426 | farthest-point: min Jaccard distance to already-selected = 0.384 |
| 36 | `GCF_013377295.1` | 137 | 31.953 | farthest-point: min Jaccard distance to already-selected = 0.381 |
| 37 | `GCF_001746835.1` | 139 | 18.784 | farthest-point: min Jaccard distance to already-selected = 0.379 |
| 38 | `GCF_006874645.1` | 127 | 11.932 | farthest-point: min Jaccard distance to already-selected = 0.379 |
| 39 | `GCF_001659705.1` | 130 | 20.374 | farthest-point: min Jaccard distance to already-selected = 0.375 |
| 40 | `GCF_003355515.1` | 127 | 16.808 | farthest-point: min Jaccard distance to already-selected = 0.371 |
| 41 | `GCF_000179915.2` | 141 | 20.688 | farthest-point: min Jaccard distance to already-selected = 0.365 |
| 42 | `GCF_008118345.1` | 151 | 40.634 | farthest-point: min Jaccard distance to already-selected = 0.365 |
| 43 | `GCF_002005145.1` | 142 | 21.482 | farthest-point: min Jaccard distance to already-selected = 0.364 |
| 44 | `GCF_000507245.1` | 141 | 16.057 | farthest-point: min Jaccard distance to already-selected = 0.360 |
| 45 | `GCF_000284315.1` | 127 | 5.392 | farthest-point: min Jaccard distance to already-selected = 0.359 |
| 46 | `GCF_003143535.1` | 134 | 5.969 | farthest-point: min Jaccard distance to already-selected = 0.358 |
| 47 | `GCF_009649955.1` | 143 | 13.92 | farthest-point: min Jaccard distance to already-selected = 0.356 |
| 48 | `GCF_000014865.1` | 145 | 17.146 | farthest-point: min Jaccard distance to already-selected = 0.354 |
| 49 | `GCF_001888165.1` | 135 | 16.426 | farthest-point: min Jaccard distance to already-selected = 0.354 |
| 50 | `GCF_000008945.1` | 115 | 18.672 | farthest-point: min Jaccard distance to already-selected = 0.354 |
| 51 | `GCF_003966365.1` | 156 | 18.505 | farthest-point: min Jaccard distance to already-selected = 0.353 |
| 52 | `GCF_000025725.1` | 156 | 22.553 | farthest-point: min Jaccard distance to already-selected = 0.353 |
| 53 | `GCF_001636925.1` | 151 | 13.124 | farthest-point: min Jaccard distance to already-selected = 0.352 |
| 54 | `GCF_000196355.1` | 128 | 2.388 | farthest-point: min Jaccard distance to already-selected = 0.350 |
| 55 | `GCF_002302535.1` | 118 | 15.041 | farthest-point: min Jaccard distance to already-selected = 0.350 |
| 56 | `GCF_000014285.1` | 159 | 19.464 | farthest-point: min Jaccard distance to already-selected = 0.346 |
| 57 | `GCF_000260985.4` | 156 | 6.128 | farthest-point: min Jaccard distance to already-selected = 0.343 |
| 58 | `GCF_000298115.2` | 122 | 17.296 | farthest-point: min Jaccard distance to already-selected = 0.343 |
| 59 | `GCF_002127965.1` | 115 | 12.714 | farthest-point: min Jaccard distance to already-selected = 0.340 |
| 60 | `GCF_000599985.1` | 149 | 17.349 | farthest-point: min Jaccard distance to already-selected = 0.338 |
| 61 | `GCF_000143845.1` | 126 | 15.604 | farthest-point: min Jaccard distance to already-selected = 0.337 |
| 62 | `GCF_002951835.1` | 146 | 19.464 | farthest-point: min Jaccard distance to already-selected = 0.336 |
| 63 | `GCF_000600005.1` | 153 | 19.902 | farthest-point: min Jaccard distance to already-selected = 0.336 |
| 64 | `GCF_000266885.1` | 147 | 16.426 | farthest-point: min Jaccard distance to already-selected = 0.333 |
| 65 | `GCF_001610875.1` | 161 | 27.017 | farthest-point: min Jaccard distance to already-selected = 0.333 |
| 66 | `GCF_000195535.1` | 144 | 16.97 | farthest-point: min Jaccard distance to already-selected = 0.331 |
| 67 | `GCF_000224745.1` | 146 | 8.864 | farthest-point: min Jaccard distance to already-selected = 0.329 |
| 68 | `GCF_000025185.1` | 158 | 13.92 | farthest-point: min Jaccard distance to already-selected = 0.326 |
| 69 | `GCF_004571195.1` | 141 | 13.135 | farthest-point: min Jaccard distance to already-selected = 0.326 |
| 70 | `GCF_000007945.1` | 144 | 21.991 | farthest-point: min Jaccard distance to already-selected = 0.326 |
| 71 | `GCF_000009905.1` | 160 | 29.604 | farthest-point: min Jaccard distance to already-selected = 0.326 |
| 72 | `GCF_002288005.1` | 118 | 7.442 | farthest-point: min Jaccard distance to already-selected = 0.323 |
| 73 | `GCF_004114615.1` | 147 | 3.77 | farthest-point: min Jaccard distance to already-selected = 0.322 |
| 74 | `GCF_000973725.1` | 149 | 9.655 | farthest-point: min Jaccard distance to already-selected = 0.321 |
| 75 | `GCF_000307165.1` | 122 | 14.275 | farthest-point: min Jaccard distance to already-selected = 0.320 |
| 76 | `GCF_000724625.1` | 138 | 15.027 | farthest-point: min Jaccard distance to already-selected = 0.319 |
| 77 | `GCF_001632845.1` | 150 | 10.849 | farthest-point: min Jaccard distance to already-selected = 0.319 |
| 78 | `GCF_001443605.1` | 147 | 23.431 | farthest-point: min Jaccard distance to already-selected = 0.317 |
| 79 | `GCF_000092365.1` | 138 | 25.622 | farthest-point: min Jaccard distance to already-selected = 0.315 |
| 80 | `GCF_000046685.1` | 110 | 2.879 | farthest-point: min Jaccard distance to already-selected = 0.312 |
| 81 | `GCF_000968135.1` | 165 | 24.336 | farthest-point: min Jaccard distance to already-selected = 0.312 |
| 82 | `GCF_001021065.1` | 156 | 38.057 | farthest-point: min Jaccard distance to already-selected = 0.310 |
| 83 | `GCF_002005425.1` | 147 | 9.283 | farthest-point: min Jaccard distance to already-selected = 0.310 |
| 84 | `GCF_005280295.1` | 148 | 3.233 | farthest-point: min Jaccard distance to already-selected = 0.307 |
| 85 | `GCF_000016885.1` | 168 | 18.599 | farthest-point: min Jaccard distance to already-selected = 0.301 |
| 86 | `GCF_000521655.1` | 171 | 38.927 | farthest-point: min Jaccard distance to already-selected = 0.300 |
| 87 | `GCF_000023265.1` | 135 | 11.086 | farthest-point: min Jaccard distance to already-selected = 0.299 |
| 88 | `GCF_000828835.1` | 158 | 22.8 | farthest-point: min Jaccard distance to already-selected = 0.299 |
| 89 | `GCF_000010305.1` | 143 | 10.031 | farthest-point: min Jaccard distance to already-selected = 0.297 |
| 90 | `GCF_000253035.1` | 147 | 16.853 | farthest-point: min Jaccard distance to already-selected = 0.296 |
| 91 | `GCF_002302615.1` | 133 | 18.264 | farthest-point: min Jaccard distance to already-selected = 0.296 |
| 92 | `GCF_001544015.1` | 157 | 18.791 | farthest-point: min Jaccard distance to already-selected = 0.295 |
| 93 | `GCF_002007605.1` | 157 | 6.754 | farthest-point: min Jaccard distance to already-selected = 0.295 |
| 94 | `GCF_000284415.1` | 147 | 18.645 | farthest-point: min Jaccard distance to already-selected = 0.294 |
| 95 | `GCF_000725365.1` | 153 | 20.317 | farthest-point: min Jaccard distance to already-selected = 0.293 |
| 96 | `GCF_000016345.1` | 158 | 6.877 | farthest-point: min Jaccard distance to already-selected = 0.293 |
| 97 | `GCF_000968535.2` | 164 | 17.118 | farthest-point: min Jaccard distance to already-selected = 0.293 |
| 98 | `GCF_000164695.2` | 150 | 20.317 | farthest-point: min Jaccard distance to already-selected = 0.293 |
| 99 | `GCF_000219105.1` | 167 | 20.317 | farthest-point: min Jaccard distance to already-selected = 0.292 |
| 100 | `GCF_900066015.1` | 153 | 15.531 | farthest-point: min Jaccard distance to already-selected = 0.292 |
