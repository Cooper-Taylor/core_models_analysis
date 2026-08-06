# ModelSEED Biochemistry Database — 2026 Update Paper Plan

**Draft plan • initial 2026-07-27 • updated 2026-07-29 (integrating notes from the paper working session led by Chris with Ray and Sam)**

> ---
> ### ⚑ Annotated working copy — 2026-08-04
>
> This is a **copy** of `ModelSEED/ModelSEEDDatabase@dev:Papers/NAR_Update_2026/PAPER_2026_PLAN.md`
> (133 lines, fetched 2026-08-04), annotated with the analysis completed in
> `/scratch/ctaylor/core_models_analysis`. The upstream file is unmodified.
>
> All additions are marked **[C.T. 2026-08]**. Original text is otherwise preserved verbatim
> so this copy can be diffed against upstream. Figures referenced here are assembled in
> [`figures/`](figures/) next to this file.
>
> **The three things a reader should take from this annotation:**
>
> 1. **The multi-source thermodynamics comparison is done** — at full corpus scale
>    (5,683 models, 6 direction sources, 22,732 LP solves, 0 errors) with figures ready.
> 2. **The headline ΔG-agreement statistic needs restating.** The naive "GC and eQuilibrator
>    agree (r = 0.91), dGPredictor does not (r = 0.08)" is *largely an artifact of a few dozen
>    extreme-leverage reactions*. On the reactions the models actually use, the three sources
>    agree at a broadly comparable, moderate level. Publishing r = 0.08 unqualified would be
>    a mischaracterisation of dGPredictor. See §Results ▸ Multi-source thermodynamics.
> 3. **⚠ The empirical hook's premise does not survive contact with the data, and the paper
>    should be re-framed around what does.** Swapping the *ΔG source* barely moves the
>    biology (3,461–3,642 growing models — a 3.2-point spread across all six sources).
>    Swapping the *direction-assignment rule* moves it enormously (100/100 → 26/100 models
>    reliably growing). **The sensitivity that matters is to the heuristic cascade, not to
>    the thermodynamics source.** Detailed in §Empirical study.
>
> Also flagged below: two of the four proposed per-model metrics **do not exist yet**, one is
> **blocked by an upstream dependency**, and two *unplanned* results (energy-generating cycles;
> Monte-Carlo direction uncertainty) are stronger than several planned ones.
>
> ---

## Snapshot

- **Prior paper:** Seaver et al. 2020, *The ModelSEED Biochemistry Database…*, Nucleic Acids Research 49(D1):D575-D588. DOI [10.1093/nar/gkaa746](https://doi.org/10.1093/nar/gkaa746). (Note: the DOI [10.1093/nar/gkaa1143](https://doi.org/10.1093/nar/gkaa1143) is the erratum, not the main paper. Cite gkaa746.)
- **Target venue:** NAR Database Issue (same as 2020 — keeps citation lineage; ~4000-word cap).
- **Draft window:** TBD by team.
- **Novel empirical hook:** systematic study of reaction-direction heuristic sensitivity across the ModelSEED v2 draft-model corpus ([10.1101/2023.10.04.556561](https://doi.org/10.1101/2023.10.04.556561)).
- **Working framing (from the 2026-07 working session):** the paper walks through three technical components in turn — compound **structures**, reaction **similarities**, and **dynamics** (the empirical study of how reaction-direction assignments from different thermodynamics sources impact the ModelSEED v2 core and full models — this *is* the empirical hook detailed below) — and then their integration into the repository, the website, and the manuscript.

## What the 2020 paper committed to (baseline for the update)

According to PubMed ([10.1093/nar/gkaa746](https://doi.org/10.1093/nar/gkaa746)), the 2020 paper's four positioning claims were:

1. **Rosetta-Stone integration** of KEGG + MetaCyc + BioCyc + 34 published models + BiGG/MetaNetX/Rhea aliases.
2. **Structure-first curation** — ChemAxon Marvin protonation @ pH 7; RDKit + OpenBabel for formula/charge derivation; mass-balance ledger via a `status` field.
3. **Thermodynamics** — eQuilibrator (primary) + Group Contribution (fallback), used to predict reaction reversibility.
4. **Community + Provenance** — GitHub PR workflow with Travis CI, quarterly releases.

**2020 baseline numbers to update:**

| | 2020 snapshot |
|---|---:|
| Compounds | 33,958 |
| Compounds with structures | 28,120 (83%) |
| Reactions | 36,193 |
| Mass+charge balanced | 25,457 (70%) |
| Functional in whole-DB FBA | 21,403 |
| Biolog conditions supported | 355 / 390 (91%) |
| Compounds with eQuilibrator ΔfG′ | 17,510 |
| Reactions with accepted ΔrG′ | 13,298 |

We'll refresh every row and add a `2026` column when the update numbers are in.

> **[C.T. 2026-08] Partial 2026 numbers now available (thermodynamics rows only).**
> Measured against the live `claude-changes` working tree of the MSDB checkout, which
> contains **56,002 non-EMPTY reactions**. These fill the bottom two rows; the structure /
> balance / Biolog rows still need the structure-curation stream.
>
> | | 2020 snapshot | 2026 (measured) | source |
> |---|---:|---:|---|
> | Reactions with ΔrG′, Group Contribution | — | **25,826** | `analyze_thermo_source_dg_agreement.py` |
> | Reactions with ΔrG′, eQuilibrator | 13,298 (accepted) | **19,498** | ″ |
> | Reactions with ΔrG′, dGPredictor | n/a (new) | **27,715** | ″ |
> | Reactions with ΔrG′, ModelSEED canonical | — | **24,910** | `THERMO_SOURCE_FBA_PIPELINE.md` |
> | Compounds with ΔfG′ | 17,510 (eQuilibrator) | *pending* | — |
>
> ⚠ **Snapshot sensitivity — reconcile before publishing.** The Group Contribution count is
> **25,812** when read from `origin/dev` but **25,826** on the local working tree. Whichever
> number goes in the paper must be quoted against the **`v2.0.0` tag** agreed in Open
> Decision 5, not against a working branch.

## Proposed section outline

### 1. Introduction

Frame: six years of community-driven biochemistry curation. Motivate: rise of atom mapping in metabolic-model applications, ML-based thermodynamics methods maturing, growing need for reaction-direction transparency and its downstream model impact.

### 2. Materials & Methods

- **Structure-curation pipeline** — PubChem validation stage (Ray's stereo-loss guard); RDKit canonicalization (one-time mass recanonicalization); curator override system (`Biochemistry/Curation/overrides/structure_picks/<curator>.tsv`); **provenance tracking** — every reconciled structure records *what* was chosen, *why*, and *who* decided (per-curator override files + rationale field), which distinguishes reconciled entries from silently overwritten ones; mass-balance exclusion mechanism (new file schema: `cpd_id, name, reason, date, curator`); per-source SMILES preserved; PubChem-alias mismap detection surfacing.
- **MetaCyc→KEGG structure alignment policy** — ModelSEED originated from KEGG conventions; when MetaCyc contains a structure absent from KEGG, the structure is adjusted so its formula matches what it would be under KEGG conventions. This prevents new MetaCyc-only structures from silently unbalancing reactions that already work under the KEGG-conforming set. Document this policy explicitly.
- **Formula-conflict resolver** — H-only vs skeleton conflicts; auto-match to `compound_*.json` for zero-disruption picks; explicit exclusion when no source matches (e.g., ascorbate radical, quercetin bissulfate).
- **ACP formula override framework** — pantetheine-inclusive formula overrides in `acps_formula_charge.tsv`.
- **AI-assisted database-wide ACP standardization** — acyl-**CoA** compounds are not affected (CoA is fully represented in the compound structure). Acyl-**ACP** compounds are different: the acyl group is bound to the **phosphopantetheine** cofactor of an acyl-carrier protein, and the phosphopantetheine moiety may or may not be included in the stored formula of a given acyl-ACP compound. That inconsistency generates false-positive mass-imbalance flags in reactions that produce or consume phosphopantetheine (biosynthetically) or that transfer acyl groups between the CoA-bound and ACP-bound pools. Methods paragraph should describe the standardization: for every acyl-ACP compound, include the phosphopantetheine formula in the compound record, so that (i) mass balance holds across the ACP-related reactions, and (ii) biosynthetic demand for phosphopantetheine becomes explicitly modelable in metabolic reconstructions. Applied AI-assisted across the full ACP family, extending the current `acps_formula_charge.tsv` framework.
- **Broader protein-carrier cofactor sweep** — the acyl-ACP / phosphopantetheine pattern (a protein-covalent cofactor carrying a reactive group between enzymes) recurs across several other classes. In each case the paper should scan the database for compounds that represent the loaded state, decide whether the carrier cofactor's atoms are included in the stored formula, and apply the same standardization to make the cofactor's biosynthetic demand modelable. Initial classes to survey:
    - **Biotinyl carriers** — biotin covalently attached (via amide bond to a lysine ε-amine) on biotin-carboxyl-carrier proteins (BCCP), used by pyruvate carboxylase, acetyl-CoA carboxylase, propionyl-CoA carboxylase, methylcrotonyl-CoA carboxylase and other carboxylases.
    - **Lipoyl carriers** — lipoic acid covalently attached (again via amide bond to a lysine ε-amine) on the H-protein of the glycine cleavage system and on the E2 subunits of the α-ketoacid dehydrogenase complexes (PDH, KGDH, BCKDH).
    - **Candidate extensions** to consider surveying but potentially defer: covalently bound **FMN/FAD** (e.g. the histidyl-FAD of succinate dehydrogenase), **molybdopterin** cofactors, **heme c** attached via thioether bonds in c-type cytochromes.

    The unifying claim is that biotinyl-carrier and lipoyl-carrier compounds should include the biotin and lipoyl moieties in their stored formulas for exactly the same reasons phosphopantetheine should be included in acyl-ACPs — the cofactor is biosynthesized, demand for it should be modelable, and reactions handling it should mass-balance. Report the classes surveyed and any deferred.
- **Multi-source thermodynamics** — rebuild of Group Contribution from MFAToolkit; dGpredictor retrained on ModelSEED; eQuilibrator refresh; OpenTECR integration. Per-reaction-class heuristics applied globally.
- **Reaction similarity** — describe pipeline for computing reaction similarity from updated compound structures; GPU regeneration is cheap (under ~10 minutes end-to-end per Ray), so the matrix is regenerated on demand as compound structures update. Include the selection procedure for the external reaction-embedding foundation model: candidates are evaluated for how discerning their embeddings are on ModelSEED reactions, and the paper defends the chosen model on scientific grounds rather than default adoption.
- **Atom mapping** — collaboration with the Nikoloski lab; describe methodology and coverage-to-date.
- **PR-time validation** — GitHub Actions replacing Travis (Travis is retired). New curated entries validated on submission via a workflow that runs the curation self-check.
- **Packaging** — PyPi distribution with a documented API.

### 3. Results

- **Growth statistics 2020→2026** (mirror Tables 2-3 of the 2020 paper).
- **Structure-curation improvements** — number of compounds with newly-picked structures, mass-balance exclusions applied, PubChem pipeline yield. Include the reduction in false-positive mass-imbalance flags attributable to protein-carrier-cofactor standardization (acyl-ACP phosphopantetheine + biotinyl-BCCP + lipoyl-carrier and whichever additional classes make the final cut), broken out per class: how many compounds required an override, how many reactions were rebalanced.
- **Reaction similarity** — the refreshed reaction-similarity matrix (regenerated against the updated compound structures) and the head-to-head comparison of external reaction-embedding foundation models, with the selection justified on discerning-power grounds.
- **Multi-source thermodynamics comparison** — heat map / correlation among the four ΔG sources; per-class heuristic direction assignments.

> **[C.T. 2026-08] ✅ DONE — figures and numbers ready.**
> Pipeline documented in [`reports/thermoComparison/THERMO_SOURCE_FBA_PIPELINE.md`](../THERMO_SOURCE_FBA_PIPELINE.md);
> agreement diagnostics in [`scripts/analyze_thermo_source_dg_agreement.py`](../../scripts/analyze_thermo_source_dg_agreement.py).
>
> **Figure 1 — pairwise ΔG′° agreement, all ModelSEED reactions.** Points coloured by
> reversibility transition between the two sources (No change / Rev→Irrev / Irrev→Rev /
> Irrev→Irrev), y = x reference, axis zoomed to |ΔG′°| ≤ 1,500 kcal/mol with off-scale
> reactions listed rather than silently dropped.
> · [F1a GC vs eQuilibrator](figures/F1a_dg_gc_vs_eq_all.png) — n = 18,477, r = 0.91
> · [F1b GC vs dGPredictor](figures/F1b_dg_gc_vs_dgp_all.png) — n = 18,603, r = 0.08
> · [F1c eQuilibrator vs dGPredictor](figures/F1c_dg_eq_vs_dgp_all.png) — n = 15,300, r = 0.17
>
> **Figure 2 — the same three panels restricted to the 239 reactions that actually occur
> across the 5,683 core models.** This is the comparison the paper should lead with, because
> it is the regime the models operate in.
> · [F2a GC vs eQ](figures/F2a_dg_gc_vs_eq_core239.png) — n = 169, r = **0.81**
> · [F2b GC vs dGP](figures/F2b_dg_gc_vs_dgp_core239.png) — n = 169, r = **0.45**
> · [F2c eQ vs dGP](figures/F2c_dg_eq_vs_dgp_core239.png) — n = 157, r = **0.51**
>
> **Figure 3 — coverage.** [F3 combined](figures/F3_coverage_combined.png) (% of the 239
> unique reactions / 182 unique compounds covered per source) and
> [F3b per-model](figures/F3b_coverage_rxn_per_model.png) (median ± σ over 5,683 models).
>
> **Figure 5 — ΔG distribution per source, stacked by assigned direction operator**
> ([F5](figures/F5_dg_distribution_per_source.png)).
> ⚠ *Axis label bug: this figure is labelled kJ/mol but ModelSEED stores kcal/mol
> (`thermo_source_figures.py:436`). Fix before submission.*
>
> **⚠ The headline correlation must be restated — r = 0.08 is misleading.** Two findings:
>
> *(a) GC and eQuilibrator are correlated by construction, not by independent agreement.*
> Both store per-compound ΔG_f, and each source's stored reaction ΔG reconstructs **exactly**
> as Σνᵢ·ΔG_f,ᵢ (GC r = 1.0000, n = 25,826, 100.0% within 1 kcal/mol; eQuilibrator r = 1.0000,
> n = 16,227, 99.7%). They are two linear maps of the *same* stoichiometry vector.
> dGPredictor has **no compound-level energies at all** — it is a direct reaction-level
> regressor. Presenting their agreement as mutual validation would be circular; the paper
> should say so explicitly.
>
> *(b) The r values are dominated by a few dozen extreme-leverage reactions.* Pearson r by
> |ΔG| window:
>
> | Pair | ≤50 | ≤100 | ≤200 | ≤500 | all |
> |---|--:|--:|--:|--:|--:|
> | GC–eQ | 0.634 | 0.760 | 0.766 | 0.783 | **0.911** |
> | GC–dGP | 0.269 | 0.405 | 0.415 | 0.263 | **0.079** |
> | eQ–dGP | 0.279 | 0.421 | 0.422 | 0.322 | **0.167** |
>
> Aggregate/polymer reactions (e.g. `rxn05017`: GC 15,908 · eQ 10,014 · dGPredictor −1.7)
> simultaneously inflate GC–eQ toward 0.91 and destroy the dGPredictor pairs. dGPredictor's
> output is effectively bounded to ±400 kcal/mol and does not scale with stoichiometry
> (r(|ΔG|, Σ|coeff|) = **0.006**, vs 0.211 GC / 0.315 eQ) — but its **IQR is the largest of
> the three** (18.9 vs 14.3 / 17.6), so the compression is purely in the tails.
> Rank correlation, immune to that leverage, is far less extreme: **Spearman ρ = 0.744 /
> 0.321 / 0.443**. Recommended phrasing: *dGPredictor is correlated in rank but on a
> different, bounded scale; Pearson r on raw ΔG hides this.*
>
> *Control — the core-239 result is range restriction, not a property of central metabolism.*
> Recomputing all-reaction r inside the same ΔG window the 239 occupy reproduces the subset
> value within 0.04 every time (0.763 / 0.442 / 0.480 vs 0.806 / 0.453 / 0.507).
>
> **For FBA, reversibility agreement matters more than ΔG agreement** — and it is high:
> eQuilibrator vs dGPredictor agree on direction for **91.7%** of the 239 core reactions.

- **NEW: reaction-direction heuristic sensitivity study** — systematic sweep over the ModelSEED v2 draft-model corpus ([10.1101/2023.10.04.556561](https://doi.org/10.1101/2023.10.04.556561)). Measures per-model: predicted growth rate, essential-gene set overlap, mass-balance survival rate, feasibility under Biolog conditions, as we swap between direction sources (eQuilibrator vs GC vs dGpredictor vs heuristic overlay). This is the empirical hook that lifts the paper above a numbers-refresh.

> **[C.T. 2026-08] ◐ PARTIALLY DONE — and the framing needs to change. See §Empirical study
> below for the full status matrix.** One of the four proposed metrics is complete at full
> corpus scale; two do not exist; one is blocked upstream. Two *unplanned* results are
> stronger than several planned ones.
>
> **Figure 4 — growth flux by direction source, 5,683 models** (median ± σ, ±2σ outliers):
> [F4a all models](figures/F4a_growth_all_models.png) ·
> [F4b growing models only](figures/F4b_growth_growing_only.png).
> **Figure 6 — bound-class transitions** induced by each source
> ([F6](figures/F6_override_transitions.png), 100-model panel).
- **Atom mapping coverage** — fraction of priority-scope reactions with mappings; example use cases in atom-tracking (C, P, S).

### 4. Discussion

- Where 6 years of community curation has taken the DB.
- Open problems: fragment-aware pKa/pKb (per Priorities list); EC number refresh beyond ExPASy name-matching (bound to have false / missing associations); obsolescence-leakage audit (obsolete compounds/reactions may be accidentally used); direction-field removal from the reaction schema now that thermoreversibility + template capture direction.

### 5. Data & software availability

- GitHub repo (dev/master branches, PR workflow, GitHub Actions CI).
- PyPi package + API.
- Updated website (with atom-mapping surfaced).
- Nikoloski-lab atom-mapping endpoint (as delivered).

## Empirical study — design sketch (reaction-direction sensitivity)

**Corpus:** ModelSEED v2 draft-model corpus published in [10.1101/2023.10.04.556561](https://doi.org/10.1101/2023.10.04.556561). Use the models as-published (no manual curation) to keep the study systematic.

**Direction sources to compare** (each is a full reversibility labeling over the ModelSEED reactions):

1. eQuilibrator-derived ΔrG′ + reversibility rules (this is the 2020 default).
2. Rebuilt Group Contribution ΔrG′ + reversibility rules.
3. dGpredictor (retrained on ModelSEED) ΔrG′ + reversibility rules.
4. Heuristic overlay on each of the above (per-reaction-class rules from the meeting notes — e.g., MetaCyc fatty-acyl ACP handling, sugar-isomer conventions).

**Per-model metrics:**

- Predicted growth rate under the model's original medium.
- Essential-gene set (single-gene knockout FBA) — Jaccard overlap across direction sources.
- Mass-balance survival rate (reactions removed as thermodynamically infeasible under each direction source).
- Feasibility of biomass production under the 390 Biolog conditions (mirroring the 2020 whole-DB FBA analysis).

**Reporting:** matrix of `(model × direction source × metric)` deltas; identify direction sources that systematically break vs preserve biology across the corpus.

> ## [C.T. 2026-08] Status of the empirical study — what exists, what doesn't
>
> Executed over the **Kegg2 core-model corpus: 5,683 models**, 6 direction sources,
> **22,732 LP solves, 32 workers, ~29 s wall, 0 errors**
> (`results/thermo_source_fba_all_models/`, driver
> [`run_thermo_source_fba_all_models.py`](../../scripts/run_thermo_source_fba_all_models.py)).
>
> ### Correction 1 — there are **six** direction sources implemented, not four
>
> The plan lists 4. The sweep runs 6, and the two extras are load-bearing controls:
>
> | # | Source | Note |
> |---|---|---|
> | 1 | eQuilibrator (2.0) | ⚠ the **static Flamholz-2012 MSDB table**, *not* live eQuilibrator 3.0. The paper must not imply otherwise. |
> | 2 | Group Contribution | rebuilt |
> | 3 | dGPredictor | ModelSEED-retrained (Freiburger branch) |
> | 4 | Original ModelSEED (canonical) | the 2020 default — the baseline the update is measured against |
> | 5 | ModelSEED (current) | eQuilibrator where available (19,083 rxns) + GC-backed fallback (5,827) = 24,910 |
> | 6 | **Implicit (on-disk)** | *no override at all* — models exactly as shipped. **This is the true null control** and the plan omits it. |
>
> *Note: sources 4 and 5 produce byte-identical direction maps and identical growth for this
> corpus — collapse them to one row in the paper.*
>
> ### Correction 2 — ⚠ the premise "direction **source** drives model behaviour" is not supported
>
> Growing models out of 5,683, by source:
>
> | Source | Growing | % |
> |---|--:|--:|
> | Group Contribution | 3,642 | 64.1% |
> | Original ModelSEED / ModelSEED (current) | 3,610 | 63.5% |
> | eQuilibrator (2.0) | 3,570 | 62.8% |
> | dGPredictor | 3,502 | 61.6% |
> | **Implicit (no override)** | 3,461 | 60.9% |
>
> **Total spread: 181 models — 3.2 percentage points.** Every thermodynamics source lands
> within ~3 points of *doing nothing at all*. A systematic sweep whose headline is "the
> sources broadly agree" is a much weaker hook than the plan assumes.
>
> Contrast the **direction-assignment rule** sweep (50-sample Monte-Carlo over ΔG uncertainty,
> 100-model panel, `results/statistical_panel/summary.csv`):
>
> | Variant | Always grow | Never grow | Uncertain | mean P(grows) |
> |---|--:|--:|--:|--:|
> | baseline (default cascade) | 100 | 0 | 0 | 1.000 |
> | 3.5 (σ band, k = 1.96) | 100 | 0 | 0 | 1.000 |
> | **H4 (best-evidence composite)** | **26** | 53 | 21 | 0.353 |
> | pforward ≥ 0.95 | 13 | 71 | 16 | 0.266 |
> | pforward ≥ 0.50 | 8 | 71 | 21 | 0.158 |
>
> 100/100 → 26/100 is a **74-point** swing, versus 3.2 points across all six ΔG sources.
>
> > **Recommended re-frame.** Keep the multi-source comparison as a *robustness* result —
> > "the choice of ΔG source is not what determines model behaviour" is a genuinely useful,
> > publishable negative finding, and it de-risks Open Decision 1. Then make the **heuristic
> > cascade** the empirical hook: it is where the sensitivity actually lives, it is
> > ModelSEED-specific (nobody else can run it), and 20 cascade variants are already built
> > (`thermo_variants/manifest.json`).
>
> ### Correction 3 — per-model metric status
>
> | Plan metric | Status | Evidence |
> |---|---|---|
> | **Growth rate** under original medium | ✅ **complete, full corpus** | 5,683 × 6, `model_results.csv` (`fba_growth_flux_*`, `fba_grows_*`, `fba_status_*`, `fba_n_overrides_*`). GC grows *faster where it grows*: median flux 77.8 vs ~37–52. |
> | **Essential-gene set overlap (Jaccard)** | ❌ **does not exist** | No gene-level analysis in the repo. Reaction-level single-KO essentiality exists (`build_growth_control.py`, 100-model panel, mean 29.0 essential of 117.8 tested) but runs **baseline map only** — no cross-source comparison, so no Jaccard is computable without a re-run. Either scope it in explicitly or drop it from the plan. |
> | **Mass-balance survival rate** | ❌ **does not exist** | Never computed or persisted. Nothing is currently removed as infeasible (0 infeasible solves, `n_errors: 0`). Needs defining before it can be measured. |
> | **Biolog feasibility (390 conditions)** | ⛔ **blocked upstream** | `functional_biolog_media` is `{}` in **all 5,683** records of `template_quality_all.jsonl`. Root cause: `MSGrowthPhenotypes.from_dict` absent in **modelseedpy 0.4.2**, so `simulate_biolog` no-ops (`reports/TEMPLATE_DIRECTION_EVAL.md:146`). Fails gracefully, produces nothing. **Blocking dependency — resolve or cut this row.** |
>
> ### Addition — two unplanned results that are stronger than several planned ones
>
> 1. **⭐ Energy-generating cycles.** **4,177 of 5,683 models (73.5%) carry ≥ 1 EGC**
>    (all ATP-group; `results/reaction_effects_all/model_flux_loops.jsonl`, 0 errors).
>    A thermodynamically-infeasible ATP-generating loop in ~3 of every 4 published draft
>    models is a **stronger and more actionable finding than any growth-rate delta**, and it
>    is exactly what reaction-direction assignment is supposed to prevent. Per-reaction ×
>    per-direction EGC probes exist for all 4 direction options across the whole corpus
>    (`effects/*.parquet`; 138,130,895 dual/shadow-price rows). Strong candidate for a
>    main-text figure — arguably *the* result of the paper.
> 2. **Direction uncertainty propagated by Monte-Carlo** (table above) — turns "which source
>    is right?" into a calibrated statement about how confident a direction call has to be
>    before the model's phenotype is stable.
>
> Also available if useful: **FVA** (100 models, 99.9% optimum: mean 42.3 blocked / 28.8
> flux-forced / 46.7 flexible) and **synthetic-lethal reaction pairs** (85/100 models have
> ≥1; 864 pairs). ⚠ Both are **baseline-map only** — no direction-sensitivity conclusion can
> be drawn from them as they stand.
>
> ### Scope caveat to state plainly in the paper
>
> The corpus swept here is the **Kegg2 *core* models** — central metabolism, median 128
> reactions and 124 compounds per model, **239 unique reactions and 182 unique compounds in
> total** across all 5,683. That is a much narrower slice than "the ModelSEED v2 draft-model
> corpus" implies, and it is *why* the ΔG sources agree so well: the core models contain none
> of the large-stoichiometry aggregate reactions where the sources diverge. Either extend the
> sweep to full genome-scale models or state the restriction explicitly — the current wording
> would over-claim.

**Explicitly out of scope:** ANIME-based heuristic-break detection. Even though ANIME appeared in the July 2024 meeting notes as a candidate, it is not part of this paper — save for a future methods paper if warranted.

## Open decisions to close with the team

Before the first full draft, the team should reach agreement on:

1. **Which direction sources make the final cut** for the empirical study — dropping any of the four saves complexity but weakens the systematic-sweep claim.

    > **[C.T. 2026-08] Evidence now available — recommend keeping all six, at no extra cost.**
    > All six are already computed; there is no complexity saving left to realise. Keeping
    > them is what *licenses* the robustness claim, and two are needed as controls:
    > **Implicit (no override)** is the null that makes "the sources barely differ"
    > interpretable, and **Original ModelSEED** is the 2020 baseline. Collapse *ModelSEED
    > (current)* into *Original ModelSEED* — byte-identical maps and identical growth on this
    > corpus. So: **5 rows, 6 sources.** Note also that the plan's item 4 ("heuristic overlay")
    > is not a ΔG source at all but a *rule* layer — it belongs on the cascade axis, where the
    > real sensitivity lives (20 variants built).
2. **Whether OpenTECR gets its own Methods paragraph or is folded into the eQuilibrator refresh.**
3. **How the Nikoloski atom-mapping work is credited** — co-authorship, collaboration acknowledgment, or joint methods citation depending on delivery timing.
4. **Timing of the direction-field removal** — do we describe it as done (schema-breaking release) or planned (Discussion)?
5. ~~**Cutoff date** for the compound / reaction growth statistics~~ — **decided 2026-07-29:** the paper snapshot is tagged **`v2.0.0`**, cut from `dev` after guide sections §5–§8 land. Until the tag is cut, the current `dev` HEAD is the provisional snapshot for tracking growth statistics. See `PAPER_2026_GUIDE.md` §1.
6. **Which external reaction-embedding foundation model** to use — requires a short evaluation (most-discerning embeddings on ModelSEED reactions) so the choice is defensible in the paper.
7. **Which protein-carrier-cofactor classes make the cut** — biotinyl and lipoyl are the strong candidates alongside acyl-ACP; whether covalent FMN/FAD, molybdopterin, and heme-c are included in this paper or deferred is a scope question that depends on how many compounds each class actually touches in the DB (a quick census would answer it).

## Explicitly out of scope for this paper

- ANIME-based heuristic-break detection.
- Kinetic-constant integration (mentioned in 2020 Introduction as a related need; still not the focus here).

## Infrastructure note (Travis → GitHub Actions)

The 2020 paper's Methods described "We utilize Travis CI along with scripts for testing data immediately, and reporting whether or not data in the pull request is valid." Travis is no longer in use. The replacement is a **GitHub Actions workflow that runs at PR submission time** to validate any newly curated entries against the existing schema and consistency checks. This should be a short Methods paragraph in the update paper, framed as an infrastructure upgrade rather than a novel contribution.

## Relationship to prior in-house work

- **Structure curation stream** (`MSD_Structures/`) — the PubChem pipeline, stereo-loss guard, curator override system, and formula-conflict resolver are all recent products of this working directory. The 2020 paper described mass-balance in one paragraph; the update should devote several paragraphs to the *pipeline* that produces mass-balanced compounds today.
- **Priority scope** — the ~9K reactions / ~6.5K compounds used in the v7.0 templates is the working slice for structure curation; the empirical study extends to the full ModelSEED v2 model corpus.

---

*Regenerate this file whenever paper scope shifts materially. Attribution for source paper information: PubMed ([10.1093/nar/gkaa746](https://doi.org/10.1093/nar/gkaa746)).*

---

## [C.T. 2026-08] Artifact index — what backs each annotation

Everything below is in `/scratch/ctaylor/core_models_analysis`. The ModelSEEDDatabase
checkout was **read only**; no MSDB or `core_models_kegg2` file was modified.

### Figures assembled for this paper ([`figures/`](figures/))

| Fig | File | Shows | Generator |
|---|---|---|---|
| 1a–c | `F1{a,b,c}_dg_*_all.png` | pairwise ΔG′° agreement, all MSDB reactions | `plot_thermo_source_dg_scatter.py` |
| 2a–c | `F2{a,b,c}_dg_*_core239.png` | same, restricted to the 239 core-model reactions | ″ `--subset` |
| 3 | `F3_coverage_combined.png` | % of 239 rxns / 182 cpds covered per source | `plot_thermo_source_coverage_bars.py` |
| 3b | `F3b_coverage_rxn_per_model.png` | per-model reaction coverage, median ± σ, n = 5,683 | ″ |
| 4a–b | `F4{a,b}_growth_*.png` | growth flux by source (all / growing-only) | `plot_thermo_source_growth_bars.py` |
| 5 | `F5_dg_distribution_per_source.png` | ΔG distribution stacked by direction operator ⚠ *axis mislabelled kJ/mol* | `thermo_source_figures.py` |
| 6 | `F6_override_transitions.png` | bound-class transitions per source (100-model panel) | ″ |

### Underlying reports and scripts

- [`reports/thermoComparison/THERMO_SOURCE_FBA_PIPELINE.md`](../THERMO_SOURCE_FBA_PIPELINE.md) — the 5,683-model × 6-source sweep, coverage, growth, and the ΔG-agreement analysis
- [`scripts/analyze_thermo_source_dg_agreement.py`](../../scripts/analyze_thermo_source_dg_agreement.py) — the four agreement diagnostics (additivity, size scaling, output range, range restriction)
- [`scripts/run_thermo_source_fba_all_models.py`](../../scripts/run_thermo_source_fba_all_models.py) — the FBA sweep
- [`scripts/build_flux_loops.py`](../../scripts/build_flux_loops.py) — energy-generating-cycle detection
- [`results/core_models_unique_reactions.json`](../../results/core_models_unique_reactions.json) — the 239-reaction union (regenerated from all 5,683 models; reproduces the documented count and per-source coverage exactly)
- `results/thermo_source_fba_all_models/`, `results/reaction_effects_all/`, `results/statistical_panel/`, `thermo_variants/manifest.json`

### Open items this annotation raises

1. **Essential-gene Jaccard** — does not exist; needs a gene-level KO sweep per direction source, or drop from the plan.
2. **Mass-balance survival rate** — undefined and uncomputed; needs a definition first.
3. **Biolog feasibility** — blocked on modelseedpy 0.4.2 (`MSGrowthPhenotypes.from_dict` missing).
4. **Corpus scope** — core models only (239 unique reactions); extend to genome-scale or state the restriction.
5. **Snapshot pinning** — regenerate all counts against tag `v2.0.0`, not a working branch (GC differs by 14 between `origin/dev` and the local tree).
6. **Figure 5 axis label** — kJ/mol → kcal/mol.
7. **Re-frame the empirical hook** from ΔG-source sensitivity (3.2 pt spread) to heuristic-cascade sensitivity (74 pt spread) + energy-generating cycles (73.5% of models).
