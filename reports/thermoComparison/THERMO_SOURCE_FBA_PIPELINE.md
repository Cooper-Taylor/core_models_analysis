# Thermo-source FBA pipeline across all 5,683 Kegg2 core models

*How reaction direction was determined from four thermodynamic sources in
isolation, how that was applied to FBA across every core model, and the
compound/reaction coverage statistics that go with it. Date: 2026-08-03.*

Scripts: [`scripts/build_thermo_source_direction_maps.py`](/scratch/ctaylor/core_models_analysis/scripts/build_thermo_source_direction_maps.py),
[`scripts/run_thermo_source_fba_all_models.py`](/scratch/ctaylor/core_models_analysis/scripts/run_thermo_source_fba_all_models.py).
Outputs: [`results/thermo_source_fba_all_models/`](/scratch/ctaylor/core_models_analysis/results/thermo_source_fba_all_models/).

---

## 1. Goal

Run the full FBA pipeline over all 5,683 Kegg2 core models (`data/core_models_kegg2/`,
symlinked to `/scratch/ctaylor/core_models_kegg2`) four times, once per thermodynamic
source, where reaction direction in each run is determined **only** from that source's
own data:

1. thermodynamic data from the original ModelSEED database ("modelseed")
2. Group Contribution data only ("group_contribution")
3. eQuilibrator (2.0) data only ("equilibrator")
4. dGPredictor data only, fine-tuned on ModelSEED ("dgpredictor")

Apart from swapping which ΔG feeds in, no new heuristics were introduced — every
direction call uses the same, unmodified cascade already shipped in the local
ModelSEEDDatabase checkout.

## 2. The four thermodynamic sources

Every ModelSEED reaction that has thermodynamic data at all carries it as a
`{dg, dge}` pair, either at the canonical top level (`deltag`/`deltagerr`) or nested
under `thermodynamics[<source label>]`:

| Source | Field read | Coverage (of 56,002 non-EMPTY reactions) |
|---|---|---|
| Original ModelSEED | canonical top-level `deltag`/`deltagerr` | 24,910 |
| Group Contribution | `thermodynamics['Group contribution']` | 25,826 |
| eQuilibrator (2.0) | `thermodynamics['eQuilibrator']` | 19,498 |
| dGPredictor | `thermodynamics['dGPredictor']` | 27,715 |

All data was read directly from the **live, currently checked-out** ModelSEEDDatabase
working tree (`/scratch/ctaylor/ModelSEEDDatabase`, branch `claude-changes`) via
`BiochemPy.Reactions().loadReactions()` / `BiochemPy.Compounds().loadCompounds()` — the
same loader `Estimate_Reaction_Reversibility.py` itself uses — rather than any older
cached snapshot, so the numbers reflect what's actually on this machine right now.

### eQuilibrator version note

This repository has **two** eQuilibrator pathways, and they are not the same vintage:

- The static table wired into every reaction's `thermodynamics['eQuilibrator']` field
  (`Biochemistry/Thermodynamics/eQuilibrator/MetaNetX_Reaction_Energies.tbl`, generated
  once by `Retrieve_eQuilibrator_Reactions_Energies.py` some time ago). Elsewhere in this
  repo's reports this is labeled `Flamholz_2012` to distinguish it from the item below.
- A separate **live eQuilibrator 3.0** pathway (`equilibrator-api 0.7.0`, labeled
  `Beber_2022` elsewhere in this repo), used only for a different, optional
  reversibility-index heuristic (`thermo_variants/eq3_*`) — not for the raw per-reaction
  ΔG that feeds the cascade here.

Per instruction, this pipeline uses the **first** one (the static, already-wired-in table)
for "eQuilibrator (2.0)". If a comparison against the live eQuilibrator 3.0 pathway is
wanted later, that is a distinct, larger undertaking (per-reaction structure/identity
resolution against a live compound cache) and out of scope here.

### dGPredictor note

The `dGPredictor` thermodynamics entries in this checkout come from Andrew Freiburger's
ModelSEED-retrained dGPredictor branch (merged in via
`freiburger/dgpredictor-modelseed-retrained-energies`) — i.e. dGPredictor **fine-tuned on
ModelSEED**, as requested, not the original dGPredictor training set.

## 3. The heuristic cascade (unchanged)

Every direction call — for all four sources — runs through the exact same cascade,
`reversibility_heuristics.DEFAULT_HEURISTICS`
([`ModelSEEDDatabase/Scripts/Thermodynamics/reversibility_heuristics.py`](/scratch/ctaylor/ModelSEEDDatabase/Scripts/Thermodynamics/reversibility_heuristics.py)),
first-match-wins:

1. **ATP synthase heuristic** — transport reaction, >1 proton compartment, exactly the
   5 ATP-synthase reagents (ATP/ADP/Pi/H₂O/H⁺) → reversible (`=`).
2. **ABC transporter heuristic** — transport reaction with a nonzero ATP coefficient →
   direction follows the sign of that coefficient.
3. **Stored-ΔG-bounds heuristic** — ΔG′° bounded over the 1e-5–2e-2 M cell-concentration
   range stays entirely negative or entirely positive → forced direction.
4. **mMΔG band heuristic** — ΔG at fixed millimolar concentrations (with special low
   local concentrations for O₂/H₂/CO₂) falls in [-2, 2] kcal/mol → reversible (`=`).
5. **Low-energy-points heuristic** — phosphate-spread × low-energy-compound score →
   directional call.
6. **Default** — always fires if nothing else did → reversible (`=`).

No new or modified heuristics were added for this analysis, and no eQuilibrator
reversibility-index or LLM-based heuristics (used elsewhere in this repo) were invoked.

**Energy-source choice, and why it matters.** For Group Contribution / eQuilibrator /
dGPredictor, the cascade was fed `per_source_energy(label)` — that source's **own** ΔG —
not `top_level_energy(db_level)`. The latter is only *eligibility*-gated by source but
still evaluates the shared canonical top-level ΔG, and its `eQuilibrator` mode additionally
falls back to Group Contribution's own reversibility when eQuilibrator itself has no data.
Both behaviors would silently blend sources together, which conflicts with "determine
direction **only** from source X." `per_source_energy` avoids that: a reaction with no data
under a given source simply gets no direction call from it.

## 4. Applying a direction map to a model (FBA)

Reused unchanged from [`scripts/growth_heuristics.py`](/scratch/ctaylor/core_models_analysis/scripts/growth_heuristics.py):

- `override_bounds(model, direction_map)` rewrites bounds only for reactions present in
  the source's map (matched via `seed_annotation.seed_id`, which also normalizes 17
  known stray `_c`-suffixed SEED annotations). **Reactions the source has no data for
  keep the model's native, as-shipped bounds** — this "overlay, don't replace" policy
  matches every other per-source variant already in this repo, and is what makes
  "direction from source X alone" meaningful (no other source is silently substituted
  in for the gaps).
- `apply_media(model)` restricts exchange-reaction uptake to the standard KBase complete
  media (`ModelSEEDDatabase/Media/KBaseMedia.cpd`).
- `find_biomass_reaction(model)` + `model.optimize()` — standard single-objective FBA,
  cobra's default solver (GLPK).

Each model is loaded once (`cobra.io.load_json_model`) and copied per source before
overriding bounds, so the four FBA runs for a model never interfere with each other.
5,683 models × 4 sources = 22,732 LP solves, parallelized across 32 worker processes
(`scripts/growth_heuristics.py`'s multiprocessing pattern); the full sweep completed in
~29 seconds with 0 errors (`results/thermo_source_fba_all_models/manifest.json`).

## 5. Definitions used for the inventory/coverage statistics

- **Unique reaction** in a model = distinct `seed_annotation.seed_id(rxn)` (base
  `rxnNNNNN` ModelSEED ID, compartment-normalized). This deliberately excludes synthetic
  `EX_`/`SK_`/`DM_`/biomass pseudo-reactions, which carry no SEED annotation and have no
  thermodynamic data to begin with.
- **Unique compound** in a model = distinct base `cpdNNNNN` ID obtained by stripping the
  trailing `_<compartment>` off each metabolite id (e.g. `cpd00001_c0` → `cpd00001`), so a
  compound appearing in both cytosol and extracellular compartments counts once.
- **"Combined across all models"** = the union of these per-model ID sets over all 5,683
  models.
- **"Has a defined energy/direction under source X"** = the base ID appears in that
  source's coverage set (§2). This is a pure data-coverage statistic, independent of
  whether FBA actually changed that reaction's bound (many reactions the source *does*
  cover already had the same bound on disk).
- **dGPredictor compound-energy coverage is N/A (reported as 0).** dGPredictor is a
  direct reaction-ΔrG predictor; unlike Group Contribution (which estimates per-compound
  formation energies and sums them) or eQuilibrator (component contribution also yields
  compound energies), it has no compound-level formation-energy analog in this database.

## 6. Results

### (a) Unique compounds and reactions per model

Across all 5,683 models (`results/thermo_source_fba_all_models/model_results.csv`,
columns `n_unique_reactions` / `n_unique_compounds`):

| | min | median | mean | max |
|---|---|---|---|---|
| Unique reactions | 20 | 128 | 123.1 | 187 |
| Unique compounds | 41 | 124 | 119.2 | 163 |

Full per-model counts are in `model_results.csv` (5,683 rows — too many to reproduce
here). The core models are fairly homogeneous in composition (central-metabolism scope),
which is why the combined total in (b) is so much smaller than 5,683 × the per-model
average.

### (b) Total unique compounds/reactions across all combined core models

| | Combined unique count |
|---|---|
| Reactions | **239** |
| Compounds | **182** |

(`results/thermo_source_fba_all_models/summary_stats.json`, `combined_unique_reactions_all_models` /
`combined_unique_compounds_all_models`.)

### (c) Unique compounds with a defined energy under the thermodynamic source, per model

| Source | min | median | mean | max |
|---|---|---|---|---|
| Original ModelSEED | 38 | 110 | 105.9 | 139 |
| Group Contribution | 37 | 109 | 104.8 | 138 |
| eQuilibrator (2.0) | 35 | 100 | 96.9 | 127 |
| dGPredictor | — | — | — | N/A (0; no compound-level energies, §5) |

Per-model values: `model_results.csv` columns `n_compounds_with_energy_<source>`.

### (d) Unique reactions with a defined direction under the thermodynamic source, per model

| Source | min | median | mean | max |
|---|---|---|---|---|
| Original ModelSEED | 18 | 114 | 110.0 | 166 |
| Group Contribution | 18 | 114 | 109.9 | 165 |
| eQuilibrator (2.0) | 2 | 91 | 87.0 | 138 |
| dGPredictor | 11 | 95 | 91.7 | 141 |

Per-model values: `model_results.csv` columns `n_reactions_with_direction_<source>`.

### (e) Combined compound-energy / reaction-direction coverage across all models

| Source | Unique compounds w/ energy (of 182 combined) | Unique reactions w/ direction (of 239 combined) |
|---|---|---|
| Original ModelSEED | 156 | 202 |
| Group Contribution | 154 | 201 |
| eQuilibrator (2.0) | 143 | 170 |
| dGPredictor | 0 (N/A) | 176 |

(`results/thermo_source_fba_all_models/summary_stats.json`,
`combined_compounds_with_energy_by_source` / `combined_reactions_with_direction_by_source`.)

### Pairwise ΔG agreement across sources, all ModelSEED reactions

Script: [`scripts/plot_thermo_source_dg_scatter.py`](/scratch/ctaylor/core_models_analysis/scripts/plot_thermo_source_dg_scatter.py).
For every pair of {Group Contribution, eQuilibrator, dGPredictor}, plotting each
reaction's ΔG under both sources (intersection of coverage — not restricted to core-model
reactions) against each other:

> **Regenerated with the dGPredictor KEGG mask applied.** Every panel involving
> dGPredictor now excludes the 17,271 reactions whose stored dGPredictor ΔG′° was
> predicted from a KEGG reaction ModelSEED does not list for them — see
> `THERMO_SOURCE_AGREEMENT_STRUCTURE.md` §1 and
> `results/thermo_agreement/dgpredictor_kegg_mask.tsv`. That is why the
> dGPredictor n values are roughly half what they were and the correlations rose
> from 0.08/0.17 to 0.74/0.69. Group Contribution and eQuilibrator are untouched,
> so the GC-vs-eQuilibrator panel is unchanged. `--no-dgp-mask` reproduces the
> old figures.

- [`fig/thermo_source_dg_scatter/dg_scatter_group_contribution_vs_equilibrator.png`](/scratch/ctaylor/core_models_analysis/reports/thermoComparison/figures/thermo_source_dg_scatter/dg_scatter_group_contribution_vs_equilibrator.png) — n = 18,477, Pearson r = 0.91
- [`dg_scatter_group_contribution_vs_dgpredictor.png`](/scratch/ctaylor/core_models_analysis/reports/thermoComparison/figures/thermo_source_dg_scatter/dg_scatter_group_contribution_vs_dgpredictor.png) — n = 9,204, Pearson r = 0.74
- [`dg_scatter_equilibrator_vs_dgpredictor.png`](/scratch/ctaylor/core_models_analysis/reports/thermoComparison/figures/thermo_source_dg_scatter/dg_scatter_equilibrator_vs_dgpredictor.png) — n = 8,768, Pearson r = 0.69

Points are colored by the **reversibility transition** between the two sources (each
reaction's own operator from §3/§4, i.e. the unmodified cascade fed that source's own ΔG):
*No change* — the two sources make the identical call, either both reversible `=` or both
irreversible in the **same** direction (`>`→`>`, `<`→`<`); *Reversible → Irreversible*;
*Irreversible → Reversible*; or *Irreversible → Irreversible* — both irreversible in
**opposite** directions (`>`→`<` or `<`→`>`), which is the only genuine direction conflict
and the one that actually changes a flux model. Counts are written to
`category_counts.tsv` beside the PNGs; for Group Contribution vs eQuilibrator: No change
14,015; Reversible→Irreversible 2,146; Irreversible→Reversible 2,173;
Irreversible→Irreversible 143 (sums to n = 18,477).

> **Changed 2026-08-17.** The last category previously held *every* both-irreversible
> pair regardless of direction, which coloured perfect agreement identically to a reversal
> and made it the largest category in every panel (7,422 of 18,477 here, now 143). All
> transition-coloured figures were regenerated under the strict definition via
> `scripts/regen_figures.py --tag transition`. Group Contribution and eQuilibrator
track each other closely (r = 0.91). The compression described in the rest of this
section was measured *before* the KEGG mask; on the masked set dGPredictor's ΔG values
are far less compressed than stated here (out-of-sample additive residual 1.58 rather
than 9.21 kcal/mol). Historically, dGPredictor's ΔG values looked far more compressed
(mostly within ±200–300 kcal/mol regardless of what GC/eQuilibrator say) and barely
correlate with either — a real signal about the fine-tuned model's output distribution,
not a plotting artifact.

A handful of reactions (4–14 depending on the pair) carry chemically extreme
(non-sentinel, but implausibly large — into the thousands of kcal/mol, e.g.
`rxn05017` at ~15,900 kcal/mol under Group Contribution) ΔG values that would otherwise
crush the entire rest of the distribution onto a few pixels on a linear axis. Each plot's
axis is zoomed to |ΔG′°| ≤ 1,500 kcal/mol and explicitly lists which reactions were
excluded from view (not from the reported n or the headline Pearson r, which reflect all
data) plus the correlation recomputed on just the shown range, so nothing is silently
hidden.

#### Why GC and eQuilibrator track each other but dGPredictor does not

Script: [`scripts/analyze_thermo_source_dg_agreement.py`](/scratch/ctaylor/core_models_analysis/scripts/analyze_thermo_source_dg_agreement.py)
(four diagnostics; all numbers below are its output).

**1. GC and eQuilibrator are the same *kind* of estimator — additive over shared
compound formation energies.** Both store per-compound ΔG_f in the compound table, and
each source's stored *reaction* ΔG reconstructs **exactly** as Σ νᵢ·ΔG_f,ᵢ from its own
compound energies:

| Source | n | r(reconstructed, stored) | median abs. residual | within 1 kcal/mol |
|---|--:|--:|--:|--:|
| Group Contribution | 25,826 | **1.0000** | 0.000 | 100.0% |
| eQuilibrator | 16,227 | **1.0000** | 0.000 | 99.7% |
| dGPredictor | — | **N/A — no compound-level formation energies exist** | | |

So GC and eQuilibrator are two *linear maps of the same stoichiometry vector*, differing
only in the per-compound coefficients. They are correlated **by construction**.
dGPredictor is a direct reaction-level regressor (§5) with no additive structure — it
predicts ΔrG from a reaction fingerprint, so nothing forces it to agree.

**2. The additive structure makes GC/eQ scale with reaction size; dGPredictor does not.**

| Source | r(\|ΔG\|, Σ\|coeff\|) | median \|ΔG\|: large- vs small-stoichiometry |
|---|--:|--:|
| Group Contribution | 0.211 | 117.7 / 6.8 = **17.3×** |
| eQuilibrator | 0.315 | 120.7 / 6.5 = **18.7×** |
| dGPredictor | **0.006** | 12.4 / 7.3 = **1.7×** |

(large = Σ\|coeff\| ≥ 20.) Summing ~30 formation energies makes \|ΔG\| explode; a
fingerprint regressor has no such mechanism.

**3. dGPredictor's output is effectively bounded.** Across 27,715 reactions it never
leaves ≈ ±400 kcal/mol, while GC/eQ run to ±16,000 — but the **IQRs are comparable**, so
the compression is entirely in the *tails*, not the bulk:

| Source | n | std | IQR | min | max |
|---|--:|--:|--:|--:|--:|
| Group Contribution | 25,826 | 148.3 | 14.3 | −4,486.0 | 15,907.7 |
| eQuilibrator | 19,498 | 123.2 | 17.6 | −4,372.5 | 10,014.0 |
| dGPredictor | 27,715 | 48.7 | **18.9** | −392.7 | 324.7 |

**4. Consequently the headline r values are dominated by extreme-leverage reactions.**
Pearson r by \|ΔG\| window:

| Pair | ≤50 | ≤100 | ≤200 | ≤500 | all |
|---|--:|--:|--:|--:|--:|
| GC–eQ | 0.634 | 0.760 | 0.766 | 0.783 | **0.911** |
| GC–dGP | 0.269 | 0.405 | 0.415 | 0.263 | **0.079** |
| eQ–dGP | 0.279 | 0.421 | 0.422 | 0.322 | **0.167** |

A few dozen aggregate/polymer reactions (e.g. `rxn05017`: GC 15,908, eQ 10,014,
dGPredictor −1.7) simultaneously *inflate* GC–eQ toward 0.91 and *destroy* the
dGPredictor pairs. Rank correlation, which is immune to that leverage, tells a much
less extreme story: Spearman ρ = 0.744 (GC–eQ), 0.321 (GC–dGP), 0.443 (eQ–dGP).
**dGPredictor is not uncorrelated — it is correlated in rank but on a different,
bounded scale, and Pearson r on the raw values hides that.**

### Pairwise ΔG agreement restricted to the 239 core-model reactions

Same three plots, restricted to the **239 unique reactions that actually occur across the
combined 5,683 core models** (`results/core_models_unique_reactions.json`, regenerated as
the union of `annotation['seed.reaction']` over every model — matches the documented
`combined_unique_reactions_all_models` count; per-source coverage 201/170/176 reproduces
§(e) exactly). Command:

```
python3 scripts/plot_thermo_source_dg_scatter.py \
  --subset results/core_models_unique_reactions.json \
  --out-subdir thermo_source_dg_scatter_core239 \
  --subset-label "  ·  239 core-model reactions"
```

- [`fig/thermo_source_dg_scatter_core239/dg_scatter_group_contribution_vs_equilibrator.png`](/scratch/ctaylor/core_models_analysis/reports/thermoComparison/figures/thermo_source_dg_scatter_core239/dg_scatter_group_contribution_vs_equilibrator.png) — n = 169, r = 0.81
- [`dg_scatter_group_contribution_vs_dgpredictor.png`](/scratch/ctaylor/core_models_analysis/reports/thermoComparison/figures/thermo_source_dg_scatter_core239/dg_scatter_group_contribution_vs_dgpredictor.png) — n = 138, r = 0.70
- [`dg_scatter_equilibrator_vs_dgpredictor.png`](/scratch/ctaylor/core_models_analysis/reports/thermoComparison/figures/thermo_source_dg_scatter_core239/dg_scatter_equilibrator_vs_dgpredictor.png) — n = 136, r = 0.77

**The relationship does not hold.** The gap between the pairs largely closes:

| Pair | scope | n | Pearson | Spearman | slope | MAD (kcal/mol) | reversibility agreement |
|---|---|--:|--:|--:|--:|--:|--:|
| GC–eQ | all | 18,477 | 0.911 | 0.744 | 0.645 | 18.0 | 76.6% |
| GC–eQ | **core-239** | 169 | **0.806** | 0.639 | 0.869 | **6.6** | 73.4% |
| GC–dGP | all | 18,603 | 0.079 | 0.321 | 0.025 | 29.6 | 63.4% |
| GC–dGP | **core-239** | 169 | **0.453** | 0.422 | 0.654 | **10.7** | 69.8% |
| eQ–dGP | all | 15,300 | 0.167 | 0.443 | 0.067 | 28.5 | 72.6% |
| eQ–dGP | **core-239** | 157 | **0.507** | 0.679 | 0.761 | **7.4** | **91.7%** |

GC–eQ *falls* (0.91 → 0.81) while both dGPredictor pairs *rise sharply* (0.08 → 0.45,
0.17 → 0.51), and every pair's OLS slope moves toward 1 and its MAD roughly halves.

This is almost entirely a **range-restriction** effect, not a property of core metabolism
per se — the control recomputes the all-reaction r inside the same ΔG window the 239
occupy and lands within ~0.04 of the subset value every time:

| Pair | core-239 r (window) | all reactions, same window |
|---|--:|--:|
| GC–eQ | 0.806 (\|ΔG\| ≤ 101) | 0.763 (n = 16,426) |
| GC–dGP | 0.453 (\|ΔG\| ≤ 136) | 0.442 (n = 17,894) |
| eQ–dGP | 0.507 (\|ΔG\| ≤ 136) | 0.480 (n = 14,478) |

The core models contain only ordinary central-metabolism reactions (GC spread across the
239: std 15.3, range −45 to +101) — none of the large-stoichiometry aggregate reactions
that generate the extreme leverage. **For the reactions the core models actually use, all
three sources agree at a broadly comparable, moderate level**, and eQuilibrator/dGPredictor
in particular agree on reversibility for 91.7% of them — which matters more for FBA than
the ΔG values themselves.

### Two more sources: "current ModelSEED" and "implicit"

Two additional direction sources extend the same pipeline:

- **ModelSEED (current)** (`modelseed_current`) -- what `Estimate_Reaction_Reversibility.py`'s
  `EQ` mode actually computes against the live database today: eQuilibrator's
  own reversibility call (`top_level_energy('EQ')` through the unmodified
  cascade) when eQuilibrator has data for the reaction (19,083 reactions),
  falling back to the reaction's existing Group-Contribution-backed stored
  `reversibility` field when it doesn't (5,827 reactions) -- i.e. `estimate_one(rxn,
  'EQ')`, called directly rather than reimplemented
  (`scripts/build_modelseed_current_direction_map.py`). 24,910 reactions total,
  identical count to (and, on inspection, the same reactions as) the earlier
  canonical top-level `modelseed` source in this checkout -- expected, since
  the canonical `deltag` field was itself populated under the same
  eQuilibrator > Group Contribution > dGPredictor priority.
- **Implicit** (`implicit`) -- no direction map and no override at all: FBA is
  run on each model exactly as it ships, using whatever reaction bounds were
  implicitly baked in when that Kegg2 core model was built/gap-filled. Unlike
  every other source here this is inherently per-model, not a global
  `{rxn: operator}` map, so it's handled as an FBA-only special case in
  `run_thermo_source_fba_all_models.py` (`override_bounds` is simply skipped).

### FBA growth outcomes by source

Number of the 5,683 models that grow (optimal status, biomass flux > 1e-6) under each
source's direction map (all other reactions at native/on-disk bounds):

| Source | Models that grow |
|---|---|
| Original ModelSEED (canonical) | 3,610 |
| Group Contribution | 3,642 |
| eQuilibrator (2.0) | 3,570 |
| dGPredictor | 3,502 |
| ModelSEED (current) | 3,610 |
| Implicit (on-disk) | 3,461 |

Per-model growth flux and solver status for all six sources are in
`model_results.csv` (`fba_grows_<source>`, `fba_growth_flux_<source>`,
`fba_status_<source>`, `fba_n_overrides_<source>`).

### Growth-flux bar charts

Script: [`scripts/plot_thermo_source_growth_bars.py`](/scratch/ctaylor/core_models_analysis/scripts/plot_thermo_source_growth_bars.py).
Median ± standard deviation of biomass growth flux across the 5 sources used
for head-to-head comparison (Group Contribution, eQuilibrator (2.0),
dGPredictor, Implicit, ModelSEED (current)):

- [`fig/thermo_source_growth_bars/growth_flux_median_std_all_models.png`](/scratch/ctaylor/core_models_analysis/reports/thermoComparison/figures/thermo_source_growth_bars/growth_flux_median_std_all_models.png)
  -- over all 5,683 models (non-growing models contribute a flux of 0): medians
  30.7-47.8, largest spread under Group Contribution (std 39.9, driven by its
  higher reaction-direction coverage and larger share of growers).
- [`growth_flux_median_std_growing_only.png`](/scratch/ctaylor/core_models_analysis/reports/thermoComparison/figures/thermo_source_growth_bars/growth_flux_median_std_growing_only.png)
  -- over only the models that grow *under that source* (n differs per bar,
  3,502-3,642, annotated on each bar): Group Contribution grows noticeably
  faster on the models it does grow (median 77.8) than the other four sources
  (37.4-52.3), which cluster closer together.

Error bars are sample standard deviation (`ddof=1`); the lower whisker is
clipped at 0 (growth flux cannot be negative) rather than left to imply
negative growth.

**Outlier markers.** Each bar overlays open circles for individual models
whose growth flux falls beyond median ± 2 standard deviations (jittered
horizontally to reduce overlap). A literal ±1σ — the whisker actually drawn —
flags 30–60% of models per source on this skewed/bimodal data (many exact-zero
non-growers, a wide grower spread), so ±2σ was used instead to get a sparse,
plottable set: 0 outliers for most sources on the all-models chart (14 for
dGPredictor), 35–246 per source on the growing-only chart, with Group
Contribution's 246 the most (its long grower-side tail down toward ~3 flux
units, well below its much higher median). Markers use a single neutral dark
outline rather than each bar's own hue — most outliers land *below* the
bar's median, inside the solid fill, where a same-hue outline would be
invisible against the identically-colored bar behind it.

### Coverage bar charts (metabolites and reactions, per model)

Script: [`scripts/plot_thermo_source_coverage_bars.py`](/scratch/ctaylor/core_models_analysis/scripts/plot_thermo_source_coverage_bars.py).
Same median ± std dev / ±2σ-outlier convention as the growth-flux bar charts,
now applied to tasks 2(c)/(d)'s coverage counts, sampled per model (n = 5,683
per bar) in `reports/thermoComparison/figures/thermo_source_coverage_bars/`:

- `coverage_pct_reactions_per_model.png` / `coverage_abs_reactions_per_model.png`
  -- % / count of a model's unique reactions with a defined direction.
- `coverage_pct_compounds_per_model.png` / `coverage_abs_compounds_per_model.png`
  -- % / count of a model's unique compounds with a defined formation energy
  (dGPredictor is flat 0/0%, per §5).
- `coverage_pct_combined_all_models.png` -- a single-value (no per-model
  sampling, no error bars) grouped-bar chart: % of the 239 combined-unique
  reactions / 182 combined-unique compounds (§6b) covered by each source,
  reactions and metabolites side by side.

Here "ModelSEED" is the canonical top-level source (`modelseed` column) rather
than the eQuilibrator-with-GC-fallback `modelseed_current` used for the
growth-flux charts: the two are numerically identical for reaction coverage in
this checkout, but only `modelseed` has a compound-level coverage column (there
is no compound-level "eQuilibrator, fallback to Group Contribution" concept in
this database), so `modelseed` is the only source that covers both metabolites
and reactions consistently. It keeps the same magenta color used for
`modelseed_current` in the growth-flux charts, since both represent the same
conceptual "ModelSEED" entity.

Group Contribution and ModelSEED track each other closely on both metabolites
and reactions (84-89% per-model medians); eQuilibrator and dGPredictor sit
lower on reaction coverage (~71-74%), and eQuilibrator (~81%) covers more
compounds than dGPredictor (0%, no compound-level energies at all). The
per-model percentage distributions are tight (std ≈ 3-5 points) with a cluster
of higher-coverage models sitting just past the ±2σ line in most sources --
visible as the dense band of outlier circles just above each bar.

## 7. Reproduction

```bash
cd /scratch/ctaylor/core_models_analysis
python3 scripts/build_thermo_source_direction_maps.py   # ~10s, writes direction maps + coverage CSVs
python3 scripts/run_thermo_source_fba_all_models.py --workers 32   # ~30s, writes model_results.csv + summary_stats.json
```

## 8. Verification performed

- Reaction/compound coverage totals cross-checked against a hand-written independent
  count over the same MSDB checkout (agreed to within a handful of reactions attributable
  to the cascade's own edge-case handling, e.g. malformed `dge` values).
- 5 randomly sampled models' `n_unique_reactions`/`n_unique_compounds` and per-source
  coverage counts were independently recomputed from the raw cobra models and matched the
  CSV exactly.
- 3 randomly sampled (model, source) FBA runs were independently re-executed outside the
  pipeline (fresh `override_bounds` + `apply_media` + `optimize()`) and matched the
  recorded growth flux to full floating-point precision.
- `manifest.json` confirms `n_errors: 0` across all 5,683 models.

## 9. Caveats / scope notes

- This is single-objective (biomass) FBA per model per source — a growth/no-growth and
  flux-magnitude comparison across sources, **not** the exhaustive per-reaction ×
  per-direction sensitivity sweep that `results/reaction_effects_all/` already provides
  for the baseline cascade.
- "No data → keep native bounds" is a deliberate, documented policy (§4), not the only
  reasonable one — an alternative would be to leave uncovered reactions fully open
  (`?`/unconstrained). We matched this repo's existing convention for per-source variants.
- The eQuilibrator source used is the static, already-wired-in MSDB table, not a fresh
  live eQuilibrator 3.0 query (§2) — confirmed as the intended interpretation of
  "eQuilibrator (2.0)."
- Reaction/compound counts exclude boundary (`EX_`/`SK_`/`DM_`) and biomass
  pseudo-reactions, since those carry no ModelSEED thermodynamic data (§5).
