# Per-source thermodynamic data: grading, recommendation, and the core-model simulations

A self-contained account of how each ModelSEED reaction's thermodynamic sources
are labelled and arbitrated, and of the flux-balance simulations run under each
one. Everything needed to read, check, or repeat the work is in this folder.

Built 2026-08-12 against ModelSEED `dev` @ **49563c6f**.

| file | contents |
|---|---|
| [`01_data.md`](01_data.md) | The inputs. Four sources, the snapshot, coverage, the experimental reference set, and the four known data defects that are vetoed. |
| [`02_notation.md`](02_notation.md) | Every symbol used anywhere in this folder, defined once. Read this before the two algorithm files. |
| [`03_grading_algorithm.md`](03_grading_algorithm.md) | The gold/silver/bronze **label**: calibration, the fusion statistics, the decision cascade, and its validation. |
| [`04_recommendation_algorithm.md`](04_recommendation_algorithm.md) | The **selection** rule: which source to actually use, per target. Includes the exact direction-risk integral, the full ablation, and the negative result that determined the design. |
| [`05_core_model_simulations.md`](05_core_model_simulations.md) | All core-model FBA data: inventory, coverage, growth, pairwise agreement, direction accuracy, and the per-reaction core table. |
| [`06_reproduce.md`](06_reproduce.md) | Exact commands, runtimes, environment, provenance, and file inventory. |
| [`tables/`](tables/) | The result tables themselves (TSV/JSON), copied here so the folder stands alone. |
| [`figures/`](figures/) | The three figures, referenced from `05`. |

---

## The one-paragraph summary

Three predictors (Group Contribution, eQuilibrator, dGPredictor-ModelSEED) and
one measurement set (TECRDB) supply ΔG′° for ModelSEED reactions. Each publishes
its own uncertainty, but on three scales that are not comparable, so the first
step is to calibrate each source's σ against measured error. From the calibrated
uncertainties we build two different things, for two different jobs. A **grade**
(gold/silver/bronze) labels how much a given number can be trusted — it works,
separating median error from 0.32 to 20.78 kcal/mol across tiers. A
**recommendation** picks which source to use — and here the same uncertainties
fail: every arbitration scheme built on them lost to a fixed source priority
when scored against experimental reaction directions. The shipped recommender
therefore uses uncertainty for feasibility vetoes, abstention, and magnitude
arbitration, but not for choosing between sources on direction. The core-model
simulations (5,683 models × 7 variants) quantify what each choice does to
metabolic models, and show that growth counts track how permissive a variant is
rather than how correct it is.

## Three results worth carrying out of this folder

1. **The reported uncertainties are usable *within* a source and not *between*
   sources.** Within a source they predict error well enough to grade and to
   abstain; between sources they pick the wrong one. §04 quantifies both halves.
2. **Growth counts do not rank thermodynamic sources.** The most permissive
   variant grows the most models; the most accurate variant grows the fewest.
   Always report permissiveness alongside. §05.3.
3. **The benchmark is partly in-sample.** eQuilibrator is fitted on TECRDB and
   dGPredictor trained on 4,001 of its measurements, so the accuracy ordering
   they win is not fully independent evidence. §04.7.
