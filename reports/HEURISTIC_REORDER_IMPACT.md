# Heuristic-Reorder Impact: before vs after

**Change.** The reversibility cascade was reordered **most-specific → least-specific** so the structure-specific ATP-synthase and ABC-transporter heuristics run **before** the general MdeltaG stored-bounds rule:

```
before:  stored_bounds → atp_synthase → abc_transporter → (ln_RI) → mMdeltaG_band → low_energy → default
after:   atp_synthase → abc_transporter → stored_bounds → (ln_RI) → mMdeltaG_band → low_energy → default
```

Applied in MSDB (`Estimate_Reaction_Reversibility` `DEFAULT_HEURISTICS`, branch `claude-changes`, commit `873f5e4`) and in this repo's port (`scripts/reversibility_lib.py` `estimate_one`). The port (`ReversibilityConfig()`) reproduces MSDB's regenerated canonical reversibility **byte-for-byte (0/56,012 mismatches)**.

Reproduce: `python scripts/analyze_reorder_impact.py` → `results/reorder_impact.json`.

---

## 1. Direction-level change: 44 of 56,012 reactions

| transition | count | meaning |
|---|---:|---|
| `<` → `=` | 8 | ATP synthase: was forced reverse by ΔG bounds, now reversible |
| `>` → `=` | 6 | ATP synthase: was forced forward by ΔG bounds, now reversible |
| `<` → `>` | 30 | ATP-driven (ABC) transporter: now follows the ATP-coefficient sign |
| **total** | **44** | |

Deciding heuristic after the reorder: **14 ATP synthase** (→ `=`), **30 ABC transporter** (→ sign of the ATP coefficient). Every one of the 44 is a reaction the old order had handed to `stored_bounds` first, which forced a directional call before the structural heuristic could fire. No `deltag`/`deltagerr` energy values changed — this is purely a cascade-ordering effect.

The 44 are all ATP synthases (`rxn08173` F(1)-ATPase; eukaryotic/organellar ATP-synthase variants `rxn09528`–`rxn09530`, `rxn10042`, `rxn13766`, `rxn30646`, `Complex_V`, `ATPSYN-RXN.*`, …) and ATP-driven transporters (`rxn11221` choline ABC transporter; the Chlamydomonas `JM_Cre_*` series).

---

## 2. Model footprint: the impact reduces to ONE reaction

Of the 44 changed reactions, **only `rxn08173` (F(1)-ATPase) is present in any core model.** The other 43 (eukaryotic/organellar ATP-synthase isoforms and the `JM_Cre_*` ABC transporters) appear in **0** of the 5,683 `core_models_kegg2` models — they exist in MSDB but not in this bacterial KEGG2 panel.

| reaction | change | heuristic | # all models (of 5,683) | # panel (of 100) |
|---|---|---|---:|---:|
| `rxn08173` F(1)-ATPase | `<` → `=` | ATP synthase | **5,121** | **97** |
| all other 43 | — | ATPS / ABCT | 0 | 0 |

So **5,121 / 5,683 (90%)** of core models — and **97 / 100** panel models — contain the one reaction whose direction changed.

---

## 3. FBA impact on the 100-model panel

Panel growth FBA with the **old** baseline direction map vs the **new** one (full rebind, `baseline_map=None`, so only the 44 changed reactions can move flux):

| metric | value |
|---|---:|
| models with a growth-status flip (grow ↔ no-grow) | **0 / 100** |
| models with a biomass-flux change (>1e-6) | **91 / 100** |
| direction of flux change | **all 91 increases** |
| mean Δ flux (changed models) | **+18.8** |
| median Δ flux | +19.6 |
| range Δ flux | +2.7 … +28.8 |

Largest movers (old → new biomass flux): `GCF_000010305.1` 10.0 → 38.8 (+28.8), `GCF_001262075.1` 9.7 → 37.5 (+27.7), `GCF_002005425.1` 9.3 → 36.5 (+27.2), `GCF_013377295.1` 32.0 → 58.4 (+26.4).

**Interpretation.** Making F(1)-ATPase **reversible** (instead of forced reverse-only) removes a binding directional constraint: 91% of panel models can now route additional flux through ATP synthase, lifting biomass flux by ~50% on average. Crucially this is **growth-status-neutral** — no model crosses the grow/no-grow threshold — so the reorder improves the realism of the flux solution without changing which organisms grow. The 9 unchanged panel models are the 3 that lack `rxn08173` plus 6 where the F(1)-ATPase direction was not flux-limiting.

---

## 4. Site / data implications

The site **baseline is the port cascade output**, so the reorder propagates into every variant (each is a diff vs that baseline) and into the panel + all-models FBA. Because `rxn08173` is in 5,121 models, the all-models and panel FBA numbers move materially (flux up, growers unchanged). Rebuild path: `export_thermo_variants` → overlay variants (`build_ai/eq3/consensus/kegg_implicit`) → `build_all_models_impact` → `build_site_data` → downstream analytics. The 4-method direction comparison (`method_comparison.json`) is independent of the cascade order and does not change.
