# Variant × Source panel-FBA comparison

Reference: KBase baseline (on-disk core_models_kegg2/ bounds, no rebinding) -- 100/100 grow.

All cells: number of panel models out of 100 whose biomass flux changed vs KBase baseline, with grow-status flips shown in parentheses.

| variant | title | GC (flux Δ / grow-flip) | EQ (flux Δ / grow-flip) | DGP (flux Δ / grow-flip) |
|---------|-------|------------------------:|------------------------:|------------------------:|
| `baseline` | ReversibilityConfig() default (matches MSDB) | 98 / 0 | 98 / 0 | 99 / 0 |
| `3.1` | Persist + use ln(reversibility_index) (Noor 2012) | 100 / 75 | 100 / 75 | 99 / 0 |
| `3.3` | Bennett-2009 per-metabolite concentration ranges | 91 / 0 | 91 / 0 | 99 / 0 |
| `3.3_wide` | Wider uniform conc window [1e-7, 0.1] M | 98 / 0 | 98 / 0 | 99 / 0 |
| `3.5` | Per-reaction sigma band: k=1.96 (95%) replaces ±2 kcal | 98 / 0 | 98 / 0 | 99 / 0 |
| `3.5_wide` | Per-reaction CC bound widening: k=1.96 on stored_bounds | 98 / 0 | 98 / 0 | 99 / 0 |
| `3.6` | Drop the low-energy-compounds list entirely | 99 / 0 | 99 / 0 | 99 / 0 |
| `3.7` | Drop the CO2 1e-4 hardcoded concentration override | 98 / 0 | 98 / 0 | 99 / 0 |
| `3.10_tight` | Tighten mMdeltaG band: ±1 kcal/mol | 98 / 0 | 98 / 0 | 99 / 0 |
| `3.10_loose` | Loosen mMdeltaG band: ±4 kcal/mol | 98 / 0 | 98 / 0 | 99 / 0 |
| `H1` | (NEW) default direction = '?' for unresolved | 98 / 0 | 98 / 0 | 99 / 0 |
| `H2` | (NEW) repair LOW_LOCAL_CONC shadow bug (O2/H2 at 1e-6 M) | 98 / 0 | 98 / 0 | 99 / 0 |
| `H3` | (NEW) repair phosphates shadow bug (ABC + low-E phosphate spread) | 99 / 21 | 99 / 21 | 99 / 0 |
| `H4` | (NEW) best-evidence composite: 3.1 + 3.5 + Bennett | 100 / 65 | 100 / 65 | 99 / 0 |

## Mean |Δ flux| per (variant, source) -- among models with Δflux > 1e-6

| variant | GC | EQ | DGP |
|---------|---:|---:|----:|
| `baseline` | 30.184 | 30.184 | 66.682 |
| `3.1` | 41.907 | 41.907 | 66.682 |
| `3.3` | 14.316 | 14.316 | 66.682 |
| `3.3_wide` | 30.681 | 30.681 | 66.682 |
| `3.5` | 30.184 | 30.184 | 66.682 |
| `3.5_wide` | 30.184 | 30.184 | 66.682 |
| `3.6` | 62.125 | 62.125 | 66.682 |
| `3.7` | 30.184 | 30.184 | 66.682 |
| `3.10_tight` | 30.184 | 30.184 | 66.682 |
| `3.10_loose` | 31.531 | 31.531 | 66.682 |
| `H1` | 30.184 | 30.184 | 66.682 |
| `H2` | 30.184 | 30.184 | 66.682 |
| `H3` | 23.525 | 23.525 | 66.682 |
| `H4` | 34.178 | 34.178 | 66.682 |

## Panel-reaction coverage per source (baseline cascade)

`growth_heuristics._bounds_for_rev` maps `=` and `?` to the same `(-1000, 1000)` bounds.  Differences between sources are only visible to FBA when a reaction's bound class (`<` / `>` / free) changes.

| source | panel rxns with source data | baseline `>` in panel | baseline `<` in panel | baseline free in panel |
|--------|----------------------------:|----------------------:|----------------------:|------------------------:|
| GC | 185 / 239 | 34 | 9 | 196 |
| EQ | 186 / 239 | 35 | 9 | 195 |
| DGP | 0 / 239 | 0 | 0 | 239 |

## Interpretation

- **GC ≡ EQ on this panel.** GC and EQ cascades agree on every panel reaction up to the bounds-class collapse: at baseline they differ on a single panel reaction's bound class, and the FBA results are byte-identical across all 14 variants.  EQ's extra coverage (more `'>'` / `'<'` calls than GC) lands on reactions outside the panel.
- **DGP cascade is a no-op on this panel.** dGPredictor has energies for **0** of the 239 panel reactions, so every reaction comes back as Incomplete (`?`) and gets `(-1000, 1000)` bounds.  That's why the DGP column is identical across all 14 variants and why `mean|Δflux|` is much larger under DGP than GC/EQ — every model is fully rebound to reversibility, which differs from the KBase on-disk bounds in ~99/100 models.
- **The variant axis is real** but only shows in GC/EQ.  The biggest grow-flip variants (`3.1`: 75, `H4`: 65, `H3`: 21) all lose growers — never gain them — because forcing a reaction's direction can only ever remove flux capacity, never add it, on an FBA solve.
