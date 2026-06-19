# Reaction-Reversibility Cascade — Defaults Decision (H1 rejected, H2+H3 adopted)

*Decision record for the change that adopts the H2 + H3 shadow-bug repairs as
the canonical reversibility cascade, rejects H1, re-derives the descriptive
100-model growth panel under the fixed cascade, and sets the new defaults.*

This supersedes the "evaluate as default" open questions in
[`REVERSIBILITY_HEURISTICS_IMPACT.md`](REVERSIBILITY_HEURISTICS_IMPACT.md) §6
for H1/H2/H3. The §3.x heuristics remain opt-in comparison variants.

---

## 1. Decisions

| Heuristic | Decision | Rationale |
|---|---|---|
| **H1** — return `?` for the bare-default branch | **Rejected; removed entirely** | The cascade's design intent is that the fall-through *is* "reversible" (`=`). At the FBA-bounds level `?` and `=` are identical (both → `(-1000, 1000)`), so H1 adds no thermodynamic information — it only relabels the "no rule fired" bucket. The `default_direction` knob, the `H1` variant, and all references were deleted from `reversibility_lib.py`, `variant_catalog.py`, `build_reversibility_notebook.py`, and the site data. |
| **H2** — O₂/H₂ `LOW_LOCAL_CONC` 1 µM override | **Adopted into baseline** | Shadow-bug repair, not a design choice. |
| **H3** — `phosphates` accumulator repair | **Adopted into baseline** | Shadow-bug repair; restores the ABC-transporter rule, the phosphate-spread term, the per-reagent proton/water skip, and the CO₂ override. |
| **§3.1, 3.3, 3.3_wide, 3.5, 3.5_wide, 3.6, 3.7, 3.10_*, H4** | **Remain opt-in** | Principled biology/statistics choices (not bug fixes). On the new panel they cost 0–60 growers; 3.1/3.3/H4 are flagged for experimental validation before adoption. |

**Why H2 and H3 are one change in the canonical script.** Both bugs live in the
same `cpd` variable shadow in `_walk_stoichiometry`: a loop reused the name
`cpd` (leaving it pinned to `PHOSPHATE_IDS[-1]` = PPi) *and* tested `cpd in rgt`
(the row's dict keys, never a compound id). The natural un-bugged loop —
`if cpd in PHOSPHATE_IDS:` on the real compound id — repairs both at once. The
port (`reversibility_lib.py`) previously exposed them as two independent knobs
(`fix_phosphates_shadow`, `fix_low_local_conc`); those knobs are now removed and
the fix is baked in so the port reproduces the fixed upstream byte-for-byte.

---

## 2. New default values

### Cascade configuration (`reversibility_lib.ReversibilityConfig()` default)
Reproduces the fixed `Estimate_Reaction_Reversibility.py` exactly:
`temperature=298.15`, `cell_min=1e-5`, `cell_max=2e-2`, `cell_conc=1e-3`,
`mm_band=2.0`, `apply_special_conc=True`, `co2_local_conc=1e-4`,
`low_energy_cpds=<MFAToolkit default 8>`, `ln_ri=None`, `per_met_conc*=None`,
`sigma_*=None`, `p_forward_threshold=None`. **No `default_direction` knob** (H1
removed); the fall-through is always `=`. H2 + H3 are baked into the
stoichiometry walk (not knobs).

### Canonical MSDB reaction directions (regenerated on `claude-changes`)
EQ-level `reversibility` after the fix (56,012 reactions):

| direction | count |
|---|---|
| `>` forward | 10,564 |
| `<` reverse | 2,078 |
| `=` reversible | 12,268 |
| `?` (no energy / incomplete) | 31,102 |

The fix changes **1,992 of 56,012** reaction directions vs the pre-fix cascade
(EQ): `=`→`>` 1,211, `>`→`=` 398, `=`→`<` 252, `<`→`=` 131 — the net forward
shift is ATP-driven uptake (ABC rule) and CO₂-sink reactions now resolved.
Parity check: the port baseline reproduces the fixed MSDB cascade with **0
mismatches across all 56,012 reactions at both GC and EQ**. Per-source
`thermodynamics` operators were refreshed for 5,899 entries.

### Descriptive 100-model growth panel (`results/selected_ids.txt`)
Re-derived from growers **under the adopted cascade** (see §3). 100 of 3,393
growers; full coverage (234/234 unique seed.reactions, 181/181 metabolites).
**59 of 100 models retained, 41 replaced** vs the previous (on-disk-bounds)
panel.

---

## 3. Methodology change — `results.csv` now reflects the adopted cascade

`select_diverse.py` selects the panel from the growers in `results/results.csv`.
Previously that file recorded FBA growth under each model's **on-disk
(template-time) bounds**, which are *not* synced with any current cascade. To
make the rerun meaningful, growth was recomputed with every reaction's bounds
**re-bound in memory from the H2+H3 cascade** (`scripts/regen_results_rebound.py`;
`core_models_kegg2/*.json` are never written). The grower universe moved
**3,461 → 3,393** (202 grew→not, 134 not→grew). An independent full-database FBA
in `build_all_models_impact.py` agrees exactly (3,393 growers), confirming the
rebind is consistent. The previous on-disk-bounds growth is preserved at
`results/results_ondisk_bounds.csv`; the previous panel at
`results/selected_ids_prev.txt`.

---

## 4. Evidence — variant impact vs the fixed baseline

All 100 panel models grow under the fixed baseline (mean biomass flux 18.66),
so any §3.x variant can only *remove* panel growth. `n_changed` is the all-DB
EQ direction-change count vs baseline; all-models grower counts start from
the 3,393 baseline.

| variant | Δrxns vs baseline | panel growers (Δ) | all-models growers | verdict |
|---|---:|---:|---:|---|
| **baseline (H2+H3)** | — | 100 | 3,393 | **default** |
| 3.1 reversibility index | 537 | 98 (−2) | 3,375 | opt-in |
| 3.3 Bennett conc. | 1,188 | 43 (−57) | 2,516 | opt-in |
| 3.3_wide conc. window | 1,767 | 100 (0) | 3,664 | opt-in |
| 3.5 1.96σ band | 246 | 100 (0) | 3,393 | opt-in (core-inert) |
| 3.5_wide 1.96σ bounds | 575 | 100 (0) | 3,393 | opt-in (core-inert) |
| 3.6 drop low-E list | 1,295 | 100 (0) | 3,887 | opt-in |
| 3.7 no special conc | 58 | 100 (0) | 3,746 | opt-in (now live post-fix) |
| 3.10_tight ±1 kcal | 15 | 100 (0) | 3,393 | opt-in |
| 3.10_loose ±4 kcal | 456 | 100 (0) | 3,777 | opt-in |
| H4 composite | 1,624 | 40 (−60) | 2,343 | opt-in |

Generated by `scripts/eval_defaults_panel.py`
(→ `results/defaults_panel_eval.json`) and `scripts/build_all_models_impact.py`
(→ `site/data/all_models_variants.json`). The aggressive directional variants
(3.3, H4) cost the most growth; the conservative bug-fixed baseline keeps the
panel intact, which is the desired default.

---

## 5. Reproduce

```bash
# MSDB (branch claude-changes): regenerate reversibility after the .py fix
cd ModelSEEDDatabase/Scripts/Thermodynamics
./Estimate_Reaction_Reversibility.py        # unfiltered report
./Estimate_Reaction_Reversibility.py GC
./Estimate_Reaction_Reversibility.py EQ     # canonical top-level reversibility
./Add_Reaction_Thermodynamics_Operators.py  # refresh per-source operators

# core_models_analysis (branch main)
python3 scripts/regen_results_rebound.py --workers 64   # results.csv under fix
python3 scripts/select_diverse.py                       # re-select 100 panel
python3 scripts/eval_defaults_panel.py                  # variant decision table
python3 scripts/export_thermo_variants.py --no-cache    # site cascade reports
python3 scripts/build_site_data.py                      # site JSON (bootstrap)
python3 scripts/build_all_models_impact.py --skip-rxnsets --force-all-variants --workers 64
python3 scripts/build_site_data.py                      # site JSON (final, fresh FBA)
```

## 6. Known follow-up

The statistical panel (`results/statistical_panel/`, P(direction) + Monte-Carlo
flux distributions surfaced as site badges) was computed under the **pre-fix**
baseline and the previous panel. Its Monte-Carlo cache (`mc_cascades__*` /
`mc_fba__*` in `notebooks/.kbcache/`) is keyed by variant/sample-count only, so
it does not invalidate on a baseline or panel change. Refreshing it correctly
requires clearing those cache entries and re-running
`scripts/run_statistical_panel.py`; left as a follow-up since it is independent
of the cascade-defaults decision recorded here.
