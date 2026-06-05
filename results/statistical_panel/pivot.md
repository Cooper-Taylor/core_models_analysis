# Statistical propagation of CC ΔG′° uncertainty through the panel

Panel size: **100 models** | MC sample size: **N = 50** per variant | Seed: 0

Two methods from `Reaction_Reversibility_Heuristics_Review.md`:

- **§ 3.6 — analytic P(direction).**  Per reaction, compute `P(mMdeltaG < 0)`, `P(mMdeltaG > 0)`, `P(|mMdeltaG| ≤ 0)` from the marginal CC normal `ΔG′° ~ N(deltag, deltagerr)` (treating the concentration term `RT·Σνᵢ ln cᵢ` as fixed).  Per-panel-rxn CSVs in `p_direction__{variant}.csv`.
- **§ 2.5 — Monte-Carlo cascade + FBA.**  Resample each reaction's ΔG′° from its marginal normal; replay the full cascade; run FBA on the panel; aggregate 5/50/95% growth-flux quantiles + `P(grows)` per model.  Per-model rows in `panel_distribution__{variant}__N{n}.csv`.

**Independence caveat.** §3.4 of the review asks for the full CC covariance matrix (`standard_dg_prime_multi`).  That covariance is not on disk -- only the marginal σ per reaction is -- so the MC sampler resamples each ΔG′° independently.  Correlations between reactions that share component-contribution groups are therefore *not* propagated; the resulting flux distributions are wider than the true posterior would give.

## Per-variant summary

| variant | title | always-grow | never-grow | uncertain | mean P(grows) | median CI90 width | point ∈ CI95 | rxns w/ sample-direction variance |
|---------|-------|------------:|-----------:|----------:|--------------:|------------------:|-------------:|----------------------------------:|
| `baseline` | ReversibilityConfig() default (matches MSDB) | 100 | 0 | 0 | 1.000 | 18.802 | 99/100 | 24 |
| `3.5` | Per-reaction sigma band: k=1.96 (95%) replaces ±2 kcal | 100 | 0 | 0 | 1.000 | 18.802 | 99/100 | 24 |
| `H4` | (NEW) best-evidence composite: 3.1 + 3.5 + Bennett | 26 | 53 | 21 | 0.353 | 5.930 | 100/100 | 25 |
| `pforward_50` | P(direction) >= 0.50 rule (§3.6) | 8 | 71 | 21 | 0.158 | 35.636 | 100/100 | 58 |
| `pforward_95` | P(direction) >= 0.95 rule (§3.6) | 13 | 71 | 16 | 0.266 | 19.889 | 100/100 | 44 |

**Column glossary.**
- `always-grow` / `never-grow` / `uncertain` — partitions the 100 panel models by `P(grows)` across the MC samples.
- `median CI90 width` — among grow-capable models, the median of `q95 − q05` of the growth-flux distribution. Bigger = more flux uncertainty propagated by the variant.
- `point ∈ CI95` — how often the point-estimate FBA (at centered ΔG′°) falls inside the MC's 5–95% interval. A high count means the point estimate is well-calibrated.
- `rxns w/ sample-direction variance` — panel reactions whose MC samples disagreed on direction.  These are the reactions actually driving flux variance.

## How to interpret the new `pforward_*` variants

The two `pforward_50` / `pforward_95` variants set the cascade's new `cfg.p_forward_threshold` knob (in `reversibility_lib`).  At `0.95`, the cascade drops both the ±2 kcal mMdeltaG band and the low-energy-points heuristic in favor of a single posterior-probability rule: a reaction is `>` if `P(mMdeltaG < 0) ≥ 0.95`, `<` if `P(mMdeltaG > 0) ≥ 0.95`, and reversible otherwise.  At `0.50` the rule reduces to the sign of the centered mMdeltaG.
