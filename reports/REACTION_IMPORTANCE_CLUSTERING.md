# Reaction direction-importance across all models + model clustering

**Date:** 2026-07-06
**Compute/analyze:** [`scripts/build_reaction_importance.py`](../scripts/build_reaction_importance.py)
**Raw (on disk, not shipped):** `results/reaction_importance_raw.jsonl` (~36 MB)
**Site data:** `site/data/reaction_importance.json` (~0.6 MB)
**UI:** standalone page `site/static/reaction_importance.html` → `/static/reaction_importance.html`
(nav link **“Reaction Importance ↗”**).

---

## What it does

For **all 5,683 core models**, every reaction is set — one at a time, from the model’s
default bounds — to each of the four options and biomass growth is re-solved:

    "<"(-1000,0)   ">"(0,1000)   "="(-1000,1000)   "?"(0,0)=off/knockout

That’s ~4.1M FBA solves (720 per model); it runs in ~30 s across 96 workers, **0 errors**.
Per (model, reaction) it derives:

- **direction-sensitivity** = spread of growth across the four options (how much the
  reaction’s direction can swing that model’s growth);
- **knockout-essential** = growth collapses to ~0 when the reaction is off (`?`);
- **boost** = best achievable growth gain over the default.

Then it (1) aggregates per reaction across all models into an importance ranking, and
(2) clusters the models by *which reactions* most influence their growth.

## Which reactions matter most (validation)

Ranked by number of models where the reaction is direction-sensitive, the top reactions
are the **central energy / glycolysis core** — a strong sanity check:

| rxn | name | models sensitive | knockout-essential | mean Δ growth |
|---|---|--:|--:|--:|
| rxn08173 | F(1)-ATPase | 3434 | 250 | 20.4 |
| rxn05145 | phosphate-transporting ATPase | 3423 | 1823 | 28.5 |
| rxn01100 | phosphoglycerate kinase | 3406 | 843 | 20.2 |
| rxn00459 | enolase | 3398 | 852 | 19.1 |
| rxn00781 | GAPDH | 3395 | 833 | 20.2 |
| rxn00558 | phosphoglucose isomerase | 3391 | 894 | 15.2 |
| rxn00777 | ribose-5-phosphate isomerase | 3375 | 2582 | 43.7 |

The page lets you re-rank the bar chart / table by any metric (models sensitive, mean Δ,
max Δ, knockout-essential, mean boost) and search/sort the full top-400 table.

## Model clusters

Each model is described by a vector of **within-model normalized influence**
(reaction’s growth-spread ÷ that model’s max spread) over the **top-150** globally
influential reactions — scale-free, so clustering captures the *pattern* of which
reactions drive a model, not its absolute growth. KMeans with **k chosen by silhouette
(k=6)**; a PCA projection gives the 2-D scatter.

Cluster sizes / mean baseline growth: the two dominant groups split roughly by growth
capability — `cluster 1` (2,841 models, mean growth 17.4) vs `cluster 3` (2,497, mean
growth 52.7) — plus smaller specialized clusters (14–230 models) with distinct signature
reactions. Click a cluster in the page to see its **signature reactions** (highest mean
within-model influence).

## How to view

```bash
cd site && python serve.py --static     # http://localhost:8080 → "Reaction Importance ↗"
```
Static (precomputed JSON); no `--live`/KBase/FBA needed. The scatter draws all 5,683
models; hovering a point shows its id, cluster, baseline growth, and top reaction.

## How to (re)generate

```bash
PY=/mnt/homes/ctaylor/conda/miniforge3/envs/core_models_analysis/bin/python
$PY scripts/build_reaction_importance.py --workers 96            # compute (resumable) + analyze
$PY scripts/build_reaction_importance.py --analyze-only          # re-cluster from existing raw
```
Compute is pure cobra FBA (resumable jsonl). Analysis uses scikit-learn (KMeans/PCA/
silhouette). `--topk` sets the # reactions used as clustering features (default 150),
`--kmin/--kmax` the silhouette scan range.

## Caveats

- Influence is **growth-flux-based** (single `slim_optimize` per test). The full FVA
  quality battery (closed-mode loops, producible/consumable) per reaction × direction ×
  5k models is computationally out of reach, so “other features” here means the
  flux-derived features (sensitivity, essentiality, boost), not the full battery.
- `?` = off/knockout (your earlier choice), so the four options are reverse / forward /
  reversible / off.
- Growth is uncapped biomass flux — comparative, not an absolute rate.
- PCA explains a modest fraction of variance and the silhouette is moderate (~0.24):
  the clusters are **soft** (the space is high-dimensional and continuous), so read the
  scatter as a projection and the clusters as broad groupings, not hard partitions.
- `results/reaction_importance_raw.jsonl` (~36 MB) is a regenerable intermediate — not
  committed / not shipped to the browser.
