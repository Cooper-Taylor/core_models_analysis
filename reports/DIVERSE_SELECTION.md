# 100 Diverse Growth Models — Selection Report

Selected from **3461 growing models** (of 5,683 total in `core_models_kegg2`).

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
| 1 | `GCF_003261575.2` | 221 | 81.095 | reaction-coverage rank 1: adds 187 previously-uncovered reactions |
| 2 | `GCF_000632475.1` | 212 | 82.456 | reaction-coverage rank 2: adds 20 previously-uncovered reactions |
| 3 | `GCF_000164985.3` | 187 | 67.751 | reaction-coverage rank 3: adds 9 previously-uncovered reactions |
| 4 | `GCF_009688965.1` | 179 | 32.565 | reaction-coverage rank 4: adds 7 previously-uncovered reactions |
| 5 | `GCF_003441595.1` | 182 | 77.258 | reaction-coverage rank 5: adds 3 previously-uncovered reactions |
| 6 | `GCF_000010625.1` | 144 | 30.826 | reaction-coverage rank 6: adds 2 previously-uncovered reactions |
| 7 | `GCF_000008525.1` | 117 | 13.7 | reaction-coverage rank 7: adds 1 previously-uncovered reactions |
| 8 | `GCF_000011305.1` | 161 | 49.465 | reaction-coverage rank 8: adds 1 previously-uncovered reactions |
| 9 | `GCF_000015745.1` | 185 | 60.29 | reaction-coverage rank 9: adds 1 previously-uncovered reactions |
| 10 | `GCF_000016745.1` | 151 | 35.323 | reaction-coverage rank 10: adds 1 previously-uncovered reactions |
| 11 | `GCF_000025265.1` | 170 | 55.173 | reaction-coverage rank 11: adds 1 previously-uncovered reactions |
| 12 | `GCF_000756615.1` | 178 | 62.063 | reaction-coverage rank 12: adds 1 previously-uncovered reactions |
| 13 | `GCF_000283635.1` | 99 | 29.597 | extreme: smallest by reactions |
| 14 | `GCF_000014405.1` | 100 | 8.422 | extreme: smallest by metabolites |
| 15 | `GCF_000599545.1` | 197 | 73.645 | extreme: largest by genes |
| 16 | `GCF_000014425.1` | 106 | 29.985 | extreme: smallest by genes |
| 17 | `GCF_000021045.1` | 199 | 87.211 | extreme: highest growth flux |
| 18 | `GCF_000195855.1` | 139 | 2.03 | extreme: lowest growth flux (still growing) |
| 19 | `GCF_000008165.1` | 194 | 81.095 | extreme: widest open exchange repertoire |
| 20 | `GCF_000007465.2` | 105 | 29.597 | extreme: narrowest open exchange repertoire |
| 21 | `GCF_000233715.2` | 156 | 36.095 | extreme: carries many rare reactions (<1% of growers) |
| 22 | `GCF_009688985.1` | 153 | 64.946 | extreme: carries many rare reactions (<1% of growers) |
| 23 | `GCF_001688725.2` | 122 | 36.465 | farthest-point: min Jaccard distance to already-selected = 0.454 |
| 24 | `GCF_002022605.1` | 111 | 23.104 | farthest-point: min Jaccard distance to already-selected = 0.442 |
| 25 | `GCF_000270285.1` | 129 | 27.496 | farthest-point: min Jaccard distance to already-selected = 0.437 |
| 26 | `GCF_008805035.1` | 112 | 13.597 | farthest-point: min Jaccard distance to already-selected = 0.431 |
| 27 | `GCF_000284095.1` | 147 | 36.947 | farthest-point: min Jaccard distance to already-selected = 0.426 |
| 28 | `GCF_000219045.1` | 132 | 34.545 | farthest-point: min Jaccard distance to already-selected = 0.415 |
| 29 | `GCF_900186975.1` | 130 | 41.125 | farthest-point: min Jaccard distance to already-selected = 0.412 |
| 30 | `GCF_000024945.1` | 132 | 56.251 | farthest-point: min Jaccard distance to already-selected = 0.411 |
| 31 | `GCF_001688905.2` | 134 | 36.723 | farthest-point: min Jaccard distance to already-selected = 0.397 |
| 32 | `GCF_002021985.1` | 132 | 16.931 | farthest-point: min Jaccard distance to already-selected = 0.397 |
| 33 | `GCF_003073475.1` | 135 | 33.713 | farthest-point: min Jaccard distance to already-selected = 0.391 |
| 34 | `GCF_000023745.1` | 138 | 31.184 | farthest-point: min Jaccard distance to already-selected = 0.390 |
| 35 | `GCF_001020955.1` | 132 | 23.205 | farthest-point: min Jaccard distance to already-selected = 0.390 |
| 36 | `GCF_000020645.1` | 129 | 38.401 | farthest-point: min Jaccard distance to already-selected = 0.388 |
| 37 | `GCF_001262075.1` | 139 | 46.889 | farthest-point: min Jaccard distance to already-selected = 0.381 |
| 38 | `GCF_900476035.1` | 148 | 61.856 | farthest-point: min Jaccard distance to already-selected = 0.377 |
| 39 | `GCF_000507245.1` | 141 | 33.964 | farthest-point: min Jaccard distance to already-selected = 0.376 |
| 40 | `GCF_001659705.1` | 130 | 50.161 | farthest-point: min Jaccard distance to already-selected = 0.375 |
| 41 | `GCF_001746835.1` | 139 | 33.052 | farthest-point: min Jaccard distance to already-selected = 0.375 |
| 42 | `GCF_000179915.2` | 141 | 35.018 | farthest-point: min Jaccard distance to already-selected = 0.368 |
| 43 | `GCF_000194135.1` | 120 | 16.772 | farthest-point: min Jaccard distance to already-selected = 0.368 |
| 44 | `GCF_000455605.1` | 128 | 19.715 | farthest-point: min Jaccard distance to already-selected = 0.364 |
| 45 | `GCF_002005145.1` | 142 | 58.204 | farthest-point: min Jaccard distance to already-selected = 0.364 |
| 46 | `GCF_004114615.1` | 147 | 28.044 | farthest-point: min Jaccard distance to already-selected = 0.359 |
| 47 | `GCF_009649955.1` | 143 | 15.052 | farthest-point: min Jaccard distance to already-selected = 0.356 |
| 48 | `GCF_001888165.1` | 135 | 24.196 | farthest-point: min Jaccard distance to already-selected = 0.356 |
| 49 | `GCF_001042635.1` | 120 | 15.53 | farthest-point: min Jaccard distance to already-selected = 0.353 |
| 50 | `GCF_000020525.1` | 128 | 34.347 | farthest-point: min Jaccard distance to already-selected = 0.352 |
| 51 | `GCF_000196355.1` | 128 | 26.033 | farthest-point: min Jaccard distance to already-selected = 0.352 |
| 52 | `GCF_008118345.1` | 151 | 64.906 | farthest-point: min Jaccard distance to already-selected = 0.349 |
| 53 | `GCF_000183405.1` | 133 | 32.453 | farthest-point: min Jaccard distance to already-selected = 0.347 |
| 54 | `GCF_000284315.1` | 127 | 39.868 | farthest-point: min Jaccard distance to already-selected = 0.346 |
| 55 | `GCF_000298115.2` | 122 | 33.904 | farthest-point: min Jaccard distance to already-selected = 0.343 |
| 56 | `GCF_002302395.1` | 131 | 30.415 | farthest-point: min Jaccard distance to already-selected = 0.341 |
| 57 | `GCF_002127965.1` | 115 | 25.954 | farthest-point: min Jaccard distance to already-selected = 0.340 |
| 58 | `GCF_000599985.1` | 149 | 58.281 | farthest-point: min Jaccard distance to already-selected = 0.338 |
| 59 | `GCF_000018865.1` | 148 | 6.828 | farthest-point: min Jaccard distance to already-selected = 0.337 |
| 60 | `GCF_000143845.1` | 126 | 23.453 | farthest-point: min Jaccard distance to already-selected = 0.337 |
| 61 | `GCF_000746585.1` | 112 | 21.756 | farthest-point: min Jaccard distance to already-selected = 0.336 |
| 62 | `GCF_002951835.1` | 146 | 47.811 | farthest-point: min Jaccard distance to already-selected = 0.336 |
| 63 | `GCF_001610875.1` | 161 | 56.556 | farthest-point: min Jaccard distance to already-selected = 0.333 |
| 64 | `GCF_000325705.1` | 163 | 36.843 | farthest-point: min Jaccard distance to already-selected = 0.331 |
| 65 | `GCF_000284415.1` | 147 | 52.298 | farthest-point: min Jaccard distance to already-selected = 0.331 |
| 66 | `GCF_004571195.1` | 141 | 31.682 | farthest-point: min Jaccard distance to already-selected = 0.331 |
| 67 | `GCF_013201825.1` | 148 | 33.809 | farthest-point: min Jaccard distance to already-selected = 0.328 |
| 68 | `GCF_000009905.1` | 160 | 40.251 | farthest-point: min Jaccard distance to already-selected = 0.327 |
| 69 | `GCF_001021065.1` | 156 | 76.113 | farthest-point: min Jaccard distance to already-selected = 0.327 |
| 70 | `GCF_000025185.1` | 158 | 24.609 | farthest-point: min Jaccard distance to already-selected = 0.326 |
| 71 | `GCF_000014285.2` | 162 | 45.35 | farthest-point: min Jaccard distance to already-selected = 0.324 |
| 72 | `GCF_002117445.1` | 139 | 26.061 | farthest-point: min Jaccard distance to already-selected = 0.323 |
| 73 | `GCF_000183745.1` | 153 | 16.111 | farthest-point: min Jaccard distance to already-selected = 0.321 |
| 74 | `GCF_000307165.1` | 122 | 29.597 | farthest-point: min Jaccard distance to already-selected = 0.320 |
| 75 | `GCF_001443605.1` | 147 | 25.956 | farthest-point: min Jaccard distance to already-selected = 0.319 |
| 76 | `GCF_000253035.1` | 147 | 41.181 | farthest-point: min Jaccard distance to already-selected = 0.318 |
| 77 | `GCF_000227745.2` | 167 | 46.894 | farthest-point: min Jaccard distance to already-selected = 0.317 |
| 78 | `GCF_000266885.1` | 147 | 39.299 | farthest-point: min Jaccard distance to already-selected = 0.316 |
| 79 | `GCF_000092365.1` | 138 | 55.135 | farthest-point: min Jaccard distance to already-selected = 0.315 |
| 80 | `GCF_009662475.1` | 134 | 17.849 | farthest-point: min Jaccard distance to already-selected = 0.312 |
| 81 | `GCF_000009985.1` | 170 | 55.581 | farthest-point: min Jaccard distance to already-selected = 0.312 |
| 82 | `GCF_003966365.1` | 156 | 51.116 | farthest-point: min Jaccard distance to already-selected = 0.312 |
| 83 | `GCF_002162355.1` | 147 | 17.387 | farthest-point: min Jaccard distance to already-selected = 0.310 |
| 84 | `GCF_010669225.1` | 137 | 42.963 | farthest-point: min Jaccard distance to already-selected = 0.308 |
| 85 | `GCF_003660165.1` | 153 | 27.636 | farthest-point: min Jaccard distance to already-selected = 0.306 |
| 86 | `GCF_000695095.2` | 151 | 41.448 | farthest-point: min Jaccard distance to already-selected = 0.306 |
| 87 | `GCF_000766665.1` | 144 | 14.938 | farthest-point: min Jaccard distance to already-selected = 0.304 |
| 88 | `GCF_000012725.1` | 159 | 23.177 | farthest-point: min Jaccard distance to already-selected = 0.303 |
| 89 | `GCF_009660225.1` | 151 | 48.366 | farthest-point: min Jaccard distance to already-selected = 0.302 |
| 90 | `GCF_000014005.1` | 162 | 49.382 | farthest-point: min Jaccard distance to already-selected = 0.301 |
| 91 | `GCF_000025725.1` | 156 | 59.863 | farthest-point: min Jaccard distance to already-selected = 0.299 |
| 92 | `GCF_000198775.1` | 151 | 24.5 | farthest-point: min Jaccard distance to already-selected = 0.298 |
| 93 | `GCF_000015565.1` | 164 | 35.323 | farthest-point: min Jaccard distance to already-selected = 0.298 |
| 94 | `GCF_000006725.1` | 128 | 34.246 | farthest-point: min Jaccard distance to already-selected = 0.297 |
| 95 | `GCF_900187045.1` | 153 | 75.289 | farthest-point: min Jaccard distance to already-selected = 0.295 |
| 96 | `GCF_002208135.1` | 161 | 41.603 | farthest-point: min Jaccard distance to already-selected = 0.295 |
| 97 | `GCF_000021905.1` | 161 | 60.377 | farthest-point: min Jaccard distance to already-selected = 0.294 |
| 98 | `GCF_000725365.1` | 153 | 60.951 | farthest-point: min Jaccard distance to already-selected = 0.293 |
| 99 | `GCF_000164695.2` | 150 | 56.142 | farthest-point: min Jaccard distance to already-selected = 0.293 |
| 100 | `GCF_000219105.1` | 167 | 48.68 | farthest-point: min Jaccard distance to already-selected = 0.292 |
