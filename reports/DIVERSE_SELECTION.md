# 100 Diverse Growth Models — Selection Report

Selected from **3461 growing models** (of 5,683 total in `core_models_kegg2`).

## Methodology

Four passes, each tagged onto every pick so the role of the model in the panel is explicit:

1. **Reaction-coverage core** — greedy max-coverage on `seed.reaction` IDs across growers. Each pick contributes the most previously-uncovered reactions.
2. **Metabolite-coverage layer** — same greedy idea on compound IDs, seeded with the reaction core. Catches transporter / cofactor diversity that the reaction-set pass missed.
3. **Forced extremes** — smallest/largest by reactions, metabolites, genes; lowest/highest growth flux; widest/narrowest open-exchange repertoire; rare-reaction champions.
4. **Farthest-point Jaccard sampling** — fills the rest of the 100 slots by repeatedly adding the grower whose reaction-set Jaccard distance to the already-picked panel is maximal.

## Coverage achieved by the 100-model panel

**Set coverage:**
- Reactions covered: **239 / 239** unique `seed.reaction` IDs (100.0%) — Pass 1 alone saturated after 12 picks
- Metabolites covered: **181 / 181** unique cpd IDs (100.0%) — Pass 2 added 0 extra picks (already covered)
- Rare reactions (<1% prevalence) hit by panel: see `n_rare_rxns` column in CSV

**Representativeness** (Jaccard distance from every non-panel grower to its nearest panel member):
- median: **0.192**, mean: 0.181, p90: 0.255, **max: 0.304**
- every one of the 3,361 non-panel growers is within Jaccard 0.30 of some panel member
- distribution:

  | Jaccard bin | growers in bin | % |
  |---|---|---|
  | [0.00, 0.05) |  116 |  3.5% |
  | [0.05, 0.10) |  302 |  9.0% |
  | [0.10, 0.15) |  545 | 16.2% |
  | [0.15, 0.20) |  922 | 27.4% |
  | [0.20, 0.25) | 1092 | 32.5% |
  | [0.25, 0.30) |  371 | 11.0% |
  | [0.30, 0.40) |   13 |  0.4% |

## Files

- `results/selected_ids.txt` — one model ID per line, ready for `for id in $(cat results/selected_ids.txt); do ...; done`
- `results/selected_models.csv` — per-model metrics + selection reason
- `results/selected_models.json` — same data, plus coverage stats, for programmatic use
- `scripts/select_diverse.py` — reproducible selection script
- `notebooks/05_DiversePanelSelection.ipynb` — interactive walkthrough

**Programmatic access:**
```python
ids = open('/scratch/ctaylor/core_models_analysis/results/selected_ids.txt').read().split()
# or:
import json
panel = json.load(open('/scratch/ctaylor/core_models_analysis/results/selected_models.json'))
panel['ids']                       # ordered list of 100 model IDs
panel['records']                   # per-model rows: id, reason, metrics
panel['coverage']                  # {reactions_covered, ..., metabolites_universe}
```

## Panel members (in the order they joined)

| # | model_id | reaction count | growth flux | selection reason |
|---|---|---|---|---|
| 1 | `GCF_003261575.2` | 221 | 81.095 | reaction-coverage rank 1: adds 187 previously-uncovered reactions |
| 2 | `GCF_009688965.1` | 179 | 32.565 | reaction-coverage rank 2: adds 23 previously-uncovered reactions |
| 3 | `GCF_000522545.2` | 210 | 83.513 | reaction-coverage rank 3: adds 12 previously-uncovered reactions |
| 4 | `GCF_000007845.1` | 192 | 75.345 | reaction-coverage rank 4: adds 4 previously-uncovered reactions |
| 5 | `GCF_003441595.1` | 182 | 77.258 | reaction-coverage rank 5: adds 4 previously-uncovered reactions |
| 6 | `GCF_000021805.1` | 143 | 32.215 | reaction-coverage rank 6: adds 3 previously-uncovered reactions |
| 7 | `GCF_000008525.1` | 117 | 13.7 | reaction-coverage rank 7: adds 1 previously-uncovered reactions |
| 8 | `GCF_000011305.1` | 161 | 49.465 | reaction-coverage rank 8: adds 1 previously-uncovered reactions |
| 9 | `GCF_000015745.1` | 185 | 60.29 | reaction-coverage rank 9: adds 1 previously-uncovered reactions |
| 10 | `GCF_000016745.1` | 151 | 35.323 | reaction-coverage rank 10: adds 1 previously-uncovered reactions |
| 11 | `GCF_000025265.1` | 170 | 55.173 | reaction-coverage rank 11: adds 1 previously-uncovered reactions |
| 12 | `GCF_000756615.1` | 178 | 62.063 | reaction-coverage rank 12: adds 1 previously-uncovered reactions |
| 13 | `GCF_000283635.1` | 99 | 29.597 | extreme: smallest by reactions |
| 14 | `GCF_000632475.1` | 212 | 82.456 | extreme: largest by metabolites; largest metabolite fingerprint |
| 15 | `GCF_000014405.1` | 100 | 8.422 | extreme: smallest by metabolites |
| 16 | `GCF_000599545.1` | 197 | 73.645 | extreme: largest by genes |
| 17 | `GCF_000014425.1` | 106 | 29.985 | extreme: smallest by genes |
| 18 | `GCF_000021045.1` | 199 | 87.211 | extreme: highest growth flux |
| 19 | `GCF_000195855.1` | 139 | 2.03 | extreme: lowest growth flux (still growing) |
| 20 | `GCF_000008165.1` | 194 | 81.095 | extreme: widest open exchange repertoire |
| 21 | `GCF_000007465.2` | 105 | 29.597 | extreme: narrowest open exchange repertoire |
| 22 | `GCF_000233715.2` | 156 | 36.095 | extreme: carries many rare reactions (<1% of growers) |
| 23 | `GCF_009688985.1` | 153 | 64.946 | extreme: carries many rare reactions (<1% of growers) |
| 24 | `GCF_000270285.1` | 129 | 27.496 | farthest-point: min Jaccard distance to already-selected = 0.458 |
| 25 | `GCF_008805035.1` | 112 | 13.597 | farthest-point: min Jaccard distance to already-selected = 0.447 |
| 26 | `GCF_002022605.1` | 111 | 23.104 | farthest-point: min Jaccard distance to already-selected = 0.442 |
| 27 | `GCF_012222825.1` | 135 | 42.464 | farthest-point: min Jaccard distance to already-selected = 0.434 |
| 28 | `GCF_000284095.1` | 147 | 36.947 | farthest-point: min Jaccard distance to already-selected = 0.427 |
| 29 | `GCF_000024945.1` | 132 | 56.251 | farthest-point: min Jaccard distance to already-selected = 0.423 |
| 30 | `GCF_002021985.1` | 132 | 16.931 | farthest-point: min Jaccard distance to already-selected = 0.419 |
| 31 | `GCF_001688905.2` | 134 | 36.723 | farthest-point: min Jaccard distance to already-selected = 0.413 |
| 32 | `GCF_900186975.1` | 130 | 41.125 | farthest-point: min Jaccard distance to already-selected = 0.412 |
| 33 | `GCF_000020525.1` | 128 | 34.347 | farthest-point: min Jaccard distance to already-selected = 0.411 |
| 34 | `GCF_001688725.2` | 122 | 36.465 | farthest-point: min Jaccard distance to already-selected = 0.406 |
| 35 | `GCF_000023745.1` | 138 | 31.184 | farthest-point: min Jaccard distance to already-selected = 0.390 |
| 36 | `GCF_001020955.1` | 132 | 23.205 | farthest-point: min Jaccard distance to already-selected = 0.390 |
| 37 | `GCF_002127965.1` | 115 | 25.954 | farthest-point: min Jaccard distance to already-selected = 0.385 |
| 38 | `GCF_001880285.1` | 146 | 37.102 | farthest-point: min Jaccard distance to already-selected = 0.385 |
| 39 | `GCF_000241025.1` | 157 | 53.458 | farthest-point: min Jaccard distance to already-selected = 0.384 |
| 40 | `GCF_000507245.1` | 141 | 33.964 | farthest-point: min Jaccard distance to already-selected = 0.382 |
| 41 | `GCF_001262075.1` | 139 | 46.889 | farthest-point: min Jaccard distance to already-selected = 0.381 |
| 42 | `GCF_001888165.1` | 135 | 24.196 | farthest-point: min Jaccard distance to already-selected = 0.380 |
| 43 | `GCF_003073475.1` | 135 | 33.713 | farthest-point: min Jaccard distance to already-selected = 0.380 |
| 44 | `GCF_000020545.1` | 140 | 27.519 | farthest-point: min Jaccard distance to already-selected = 0.378 |
| 45 | `GCF_000025945.1` | 133 | 33.543 | farthest-point: min Jaccard distance to already-selected = 0.371 |
| 46 | `GCF_005221305.1` | 143 | 44.538 | farthest-point: min Jaccard distance to already-selected = 0.371 |
| 47 | `GCF_000253035.1` | 147 | 41.181 | farthest-point: min Jaccard distance to already-selected = 0.371 |
| 48 | `GCF_001042635.1` | 120 | 15.53 | farthest-point: min Jaccard distance to already-selected = 0.369 |
| 49 | `GCF_000179915.2` | 141 | 35.018 | farthest-point: min Jaccard distance to already-selected = 0.368 |
| 50 | `GCF_000194135.1` | 120 | 16.772 | farthest-point: min Jaccard distance to already-selected = 0.368 |
| 51 | `GCF_001746835.1` | 139 | 33.052 | farthest-point: min Jaccard distance to already-selected = 0.368 |
| 52 | `GCF_008118345.1` | 151 | 64.906 | farthest-point: min Jaccard distance to already-selected = 0.367 |
| 53 | `GCF_000266885.1` | 147 | 39.299 | farthest-point: min Jaccard distance to already-selected = 0.366 |
| 54 | `GCF_000525675.1` | 120 | 13.201 | farthest-point: min Jaccard distance to already-selected = 0.362 |
| 55 | `GCF_000746585.1` | 112 | 21.756 | farthest-point: min Jaccard distance to already-selected = 0.361 |
| 56 | `GCF_000599985.1` | 149 | 58.281 | farthest-point: min Jaccard distance to already-selected = 0.361 |
| 57 | `GCF_004114615.1` | 147 | 28.044 | farthest-point: min Jaccard distance to already-selected = 0.359 |
| 58 | `GCF_009649955.1` | 143 | 15.052 | farthest-point: min Jaccard distance to already-selected = 0.356 |
| 59 | `GCF_000298115.2` | 122 | 33.904 | farthest-point: min Jaccard distance to already-selected = 0.353 |
| 60 | `GCF_000143845.1` | 126 | 23.453 | farthest-point: min Jaccard distance to already-selected = 0.352 |
| 61 | `GCF_000183405.1` | 133 | 32.453 | farthest-point: min Jaccard distance to already-selected = 0.351 |
| 62 | `GCF_001659705.1` | 130 | 50.161 | farthest-point: min Jaccard distance to already-selected = 0.351 |
| 63 | `GCF_001021065.1` | 156 | 76.113 | farthest-point: min Jaccard distance to already-selected = 0.348 |
| 64 | `GCF_000307165.1` | 122 | 29.597 | farthest-point: min Jaccard distance to already-selected = 0.346 |
| 65 | `GCF_000284315.1` | 127 | 39.868 | farthest-point: min Jaccard distance to already-selected = 0.346 |
| 66 | `GCF_003261295.1` | 133 | 21.118 | farthest-point: min Jaccard distance to already-selected = 0.346 |
| 67 | `GCF_001610875.1` | 161 | 56.556 | farthest-point: min Jaccard distance to already-selected = 0.345 |
| 68 | `GCF_002302395.1` | 131 | 30.415 | farthest-point: min Jaccard distance to already-selected = 0.341 |
| 69 | `GCF_003966365.1` | 156 | 51.116 | farthest-point: min Jaccard distance to already-selected = 0.341 |
| 70 | `GCF_000973725.1` | 149 | 9.583 | farthest-point: min Jaccard distance to already-selected = 0.340 |
| 71 | `GCF_000196135.1` | 126 | 18.422 | farthest-point: min Jaccard distance to already-selected = 0.339 |
| 72 | `GCF_000196355.1` | 128 | 26.033 | farthest-point: min Jaccard distance to already-selected = 0.339 |
| 73 | `GCF_000009905.1` | 160 | 40.251 | farthest-point: min Jaccard distance to already-selected = 0.337 |
| 74 | `GCF_000018865.1` | 148 | 6.828 | farthest-point: min Jaccard distance to already-selected = 0.336 |
| 75 | `GCF_000014005.1` | 162 | 49.382 | farthest-point: min Jaccard distance to already-selected = 0.333 |
| 76 | `GCF_002162355.1` | 147 | 17.387 | farthest-point: min Jaccard distance to already-selected = 0.333 |
| 77 | `GCF_001443605.1` | 147 | 25.956 | farthest-point: min Jaccard distance to already-selected = 0.331 |
| 78 | `GCF_000009985.1` | 170 | 55.581 | farthest-point: min Jaccard distance to already-selected = 0.330 |
| 79 | `GCF_000025185.1` | 158 | 24.609 | farthest-point: min Jaccard distance to already-selected = 0.329 |
| 80 | `GCF_000227745.2` | 167 | 46.894 | farthest-point: min Jaccard distance to already-selected = 0.329 |
| 81 | `GCF_000092365.1` | 138 | 55.135 | farthest-point: min Jaccard distance to already-selected = 0.328 |
| 82 | `GCF_003660165.1` | 153 | 27.636 | farthest-point: min Jaccard distance to already-selected = 0.327 |
| 83 | `GCF_000183745.1` | 153 | 16.111 | farthest-point: min Jaccard distance to already-selected = 0.321 |
| 84 | `GCF_004551665.1` | 149 | 33.599 | farthest-point: min Jaccard distance to already-selected = 0.321 |
| 85 | `GCF_003097575.1` | 141 | 20.711 | farthest-point: min Jaccard distance to already-selected = 0.320 |
| 86 | `GCF_003143555.1` | 154 | 35.762 | farthest-point: min Jaccard distance to already-selected = 0.318 |
| 87 | `GCF_000953635.1` | 151 | 33.949 | farthest-point: min Jaccard distance to already-selected = 0.315 |
| 88 | `GCF_000014285.1` | 159 | 45.35 | farthest-point: min Jaccard distance to already-selected = 0.315 |
| 89 | `GCF_000177635.2` | 144 | 33.574 | farthest-point: min Jaccard distance to already-selected = 0.314 |
| 90 | `GCF_001682195.2` | 165 | 29.07 | farthest-point: min Jaccard distance to already-selected = 0.314 |
| 91 | `GCF_000766665.1` | 144 | 14.938 | farthest-point: min Jaccard distance to already-selected = 0.313 |
| 92 | `GCF_000006725.1` | 128 | 34.246 | farthest-point: min Jaccard distance to already-selected = 0.312 |
| 93 | `GCF_009662475.1` | 134 | 17.849 | farthest-point: min Jaccard distance to already-selected = 0.312 |
| 94 | `GCF_002005425.1` | 147 | 33.255 | farthest-point: min Jaccard distance to already-selected = 0.312 |
| 95 | `GCF_010669225.1` | 137 | 42.963 | farthest-point: min Jaccard distance to already-selected = 0.308 |
| 96 | `GCF_000014965.1` | 161 | 52.588 | farthest-point: min Jaccard distance to already-selected = 0.307 |
| 97 | `GCF_000022565.1` | 147 | 57.085 | farthest-point: min Jaccard distance to already-selected = 0.305 |
| 98 | `GCF_000725345.1` | 180 | 56.251 | farthest-point: min Jaccard distance to already-selected = 0.305 |
| 99 | `GCF_000180175.2` | 161 | 53.963 | farthest-point: min Jaccard distance to already-selected = 0.305 |
| 100 | `GCF_004299785.2` | 121 | 33.809 | farthest-point: min Jaccard distance to already-selected = 0.304 |
