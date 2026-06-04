# Core Models KEGG2 — Growth Analysis Summary

- Models analyzed: **5683**
- Media: ModelSEEDDatabase `KBaseMedia.cpd` (complete media, 347 compounds)
- Growth criterion: FBA on biomass reaction (`bio1`) > 1e-6

## Outcomes
| Outcome | Count | % |
|---|---|---|
| **Grows (biomass flux > 1e-6)** | 3461 | 60.9% |
| Optimal solve but zero biomass flux | 2222 | 39.1% |
| Non-optimal solver status (infeasible/unbounded/etc.) | 0 | 0.0% |
| No biomass reaction found | 0 | 0.0% |
| Model load / processing error | 0 | 0.0% |

## Solver status distribution
- `optimal` : 5683

## Biomass flux among growers (n=3461)
- min: 2.03
- median: 52.3
- mean: 52.39
- max: 87.21

## Model sizes (n_reactions)
- min: 41  median: 156  max: 221  mean: 151.4

## Files
- `results/results.csv` — per-model row
- `results/growers.csv` — subset where biomass flux > 1e-6
- `results/non_growers.csv` — everything else (includes errors)
- `logs/failures.log` — solver-status / exception details
- `scripts/analyze_growth.py` — analysis script
- `scripts/summarize.py` — this summary script
- `notebooks/01_GrowthFBA_Pipeline.ipynb` — interactive walkthrough
