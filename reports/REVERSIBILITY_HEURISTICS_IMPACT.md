# Reaction-Reversibility Heuristics: Impact Across the ModelSEED Database and Core-Model Panels

**Companion to:** [REACTION_REVERSIBILITY_HEURISTICS_REVIEW.md](REACTION_REVERSIBILITY_HEURISTICS_REVIEW.md)
**Pipeline:** notebook `06_ReactionReversibilityHeuristics.ipynb` (build script: `scripts/build_reversibility_notebook.py`) + `scripts/build_all_models_impact.py`
**Data:** `site/data/manifest.json`, `site/data/variants/*.json`, `site/data/all_models_variants.json`
**Generated:** 2026-06-15

---

## 1. Scope and methodology

The Heuristics Review proposes a set of one-knob configuration changes to the
ModelSEED cascade in
`ModelSEEDDatabase/Scripts/Thermodynamics/Estimate_Reaction_Reversibility.py`.
Each suggestion is encoded as a single entry in `scripts/variant_catalog.py`
and re-runs the cascade across all 56,012 MSDB reactions; the resulting
reversibility map is then applied (in memory) to two model populations:

| Population | Size | Source | What it tells us |
|---|---|---|---|
| **Descriptive panel** | 100 | `results/selected_ids.txt` — diverse growers from notebook 05 | Hand-picked reproducible coverage; FBA detail per model |
| **All core models** | 5,683 | `data/core_models_kegg2/*.json` | Database-wide impact; tail behavior; whole-pipeline robustness |

FBA in both populations is always *fully re-bound* against the chosen
reversibility map (no template-time bound leakage), so that any FBA delta
isolates the heuristic change itself.

**No on-disk artifacts under `ModelSEEDDatabase/` or `core_models_kegg2/` are
mutated.** Every variant runs the cascade in memory and overrides cobra
reaction bounds at solve time.

Each variant carries three independent impact axes:

1. **Database direction shift** — how many of the 56,012 MSDB reactions
   flip direction (`>` / `<` / `=` / `?`) versus the baseline cascade.
2. **Reaction footprint in the model population** — how many `(model,
   reaction)` instances those direction shifts touch when intersected with
   model reaction sets (244 distinct SEED reactions exist across all
   5,683 core models; many MSDB reactions never appear in a core model).
3. **Biological FBA impact** — biomass-flux change and grow-status flips
   when the variant's bounds are applied. Reported for both the 100-model
   descriptive panel and the full 5,683-model database.

The site at `core_models_analysis/site/` (run `python3 site/serve.py`)
exposes this same data interactively, including a scope toggle that
switches the Variant Browser between the two populations.

---

## 2. Database-wide direction shift

Counts are vs the byte-for-byte MSDB cascade baseline (`ReversibilityConfig()`
defaults), counted at the EQ level (i.e. reactions where eQuilibrator
Component Contribution provided ΔG′°). The `?` bucket (31,102 reactions in
every variant) reflects reactions with no usable ΔG′° from either CC or GC;
they never enter the heuristic cascade and so are inert under every knob
except H1.

| Variant | Section | MSDB rxns changed | New `>` | New `<` | New `=` | New `?` | Citation/source |
|---|---|---:|---:|---:|---:|---:|---|
| **baseline** | (reference) | 0 | 9,751 | 1,957 | 13,202 | 31,102 | Henry 2007 |
| **3.1** | § 3.1 | **3,256** | 8,187 | 4,837 | 11,886 | 31,102 | Noor 2012 — `ln_reversibility_index` |
| **3.3** | § 3.3 | 1,316 | 10,627 | 2,141 | 12,142 | 31,102 | Bennett 2009 per-metabolite ranges |
| **3.3_wide** | § 3.3 | 1,796 | 8,379 | 1,533 | 14,998 | 31,102 | Wider uniform `[1e-7, 0.1] M` |
| **3.5** | § 3.5 | 304 | 9,453 | 1,951 | 13,506 | 31,102 | 1.96·σ CC band on mMdeltaG |
| **3.5_wide** | § 3.5 | 546 | 9,311 | 1,851 | 13,748 | 31,102 | 1.96·σ CC band on stored bounds |
| **3.6** | § 3.6 | **3,070** | 6,869 | 1,769 | 16,272 | 31,102 | Drop low-energy-compounds list |
| **3.7** | § 3.7 | 0 | 9,751 | 1,957 | 13,202 | 31,102 | Disable CO₂/O₂ 1e-4 override (gated by H3) |
| **3.10_tight** | § 3.10 | 7 | 9,758 | 1,957 | 13,195 | 31,102 | mMdeltaG band ±1 kcal/mol |
| **3.10_loose** | § 3.10 | 696 | 9,137 | 1,875 | 13,898 | 31,102 | mMdeltaG band ±4 kcal/mol |
| **H1** | § H1 | **6,522** | 9,751 | 1,957 | 6,680 | **37,624** | Return `?` for bare default |
| **H2** | § H2 | 1 | 9,751 | 1,956 | 13,203 | 31,102 | Repair O₂/H₂ shadow bug |
| **H3** | § H3 | **1,989** | 10,562 | 2,079 | 12,269 | 31,102 | Repair `phosphates` shadow bug (ABC rule live) |
| **H4** | (composite) | **3,317** | 9,800 | 3,711 | 11,399 | 31,102 | Stack of 3.1 + 3.5 + Bennett |

**What the database-wide pattern says:**

- The two largest *direction-shift* footprints are **3.1** (the reversibility
  index from Noor 2012, 3,256 changes) and **3.6** (dropping the
  low-energy-compounds list, 3,070 changes) — independent of the model
  population. Both target the legacy `points × mMdeltaG > 2` heuristic
  by either replacing it with a principled quantity or removing it
  entirely; the residual changes are different sets of reactions.
- **3.7** alone is a no-op: the CO₂ override is unreachable without the
  H3 shadow-bug repair (and even with it, eQuilibrator's CO₂ species
  model would handle the carbonate ladder correctly without the override).
  This is a finding by itself — the heuristic was advertised as live, was
  dead, and disabling it has no effect until a *different* line is fixed.
- **H1** (return `?` instead of `=` for unresolved reactions) shifts 6,522
  reactions but **does not change any FBA bound** — `?` is mapped to
  `(-1000, 1000)` (the conservative ModelSEED default), same as `=`. It is
  a curation signal, not a biology signal.
- **H3** flips 1,989 reactions just by fixing a typo (`cpd in rgt` checks
  dict keys instead of the compound id, making the `phosphates` accumulator
  always empty); 1,209 of those move from `=` to `>` because ABC-driven
  uptake reactions become correctly directional.

---

## 3. Reaction footprint in model populations

A direction shift in MSDB only matters for biology if the affected reaction
actually appears in a model. The core-models population uses only **244
distinct SEED reactions** across all 5,683 models (typical model = 125–187
reactions), so many MSDB variants never touch the FBA bounds of any core
model even when they flip thousands of reactions in MSDB.

| Variant | DB rxns changed | Panel rxn-instances (100 mdl) | Panel models containing change | All-DB rxn-instances (5683 mdl) | All-DB models containing change |
|---|---:|---:|---:|---:|---:|
| 3.1 | 3,256 | 1,007 | 100 | 61,386 | 5,672 |
| 3.3 | 1,316 | 1,173 | 100 | 72,321 | 5,669 |
| 3.3_wide | 1,796 | 372 | 100 | 21,737 | 5,583 |
| 3.5 | 304 | **0** | 0 | **0** | 0 |
| 3.5_wide | 546 | **0** | 0 | 92 | 92 |
| 3.6 | 3,070 | 774 | 100 | 44,913 | 5,682 |
| 3.7 | 0 | 0 | 0 | 0 | 0 |
| 3.10_tight | 7 | **0** | 0 | **0** | 0 |
| 3.10_loose | 696 | 420 | 100 | 24,129 | 5,637 |
| H1 | 6,522 | 3,740 | 100 | 226,695 | 5,683 |
| H2 | 1 | 0 | 0 | 0 | 0 |
| H3 | 1,989 | 950 | 100 | 59,127 | 5,678 |
| H4 | 3,317 | 1,841 | 100 | 113,894 | 5,673 |

Take-aways:

- **3.5 and 3.5_wide are essentially core-model-inert.** The reactions
  whose CC σ pulls them across the reversibility band live in MSDB tails
  (cofactor-edge, secondary metabolism) absent from the curated core
  models. To see them you'd need to widen the population beyond core
  models (e.g. genome-scale reconstructions).
- **3.10_tight (±1 kcal/mol) changes 7 reactions in MSDB and zero in any
  core model.** Tightening the Henry-band has no biology consequence at
  this resolution. **3.10_loose (±4 kcal/mol) changes 696 MSDB reactions
  and 420 core-rxn-instances** — the asymmetry is real (more reactions
  sit just outside the existing band than just inside).
- **3.1 and 3.6 saturate the model population.** Almost every panel /
  all-models row contains at least one reaction the variant flipped,
  even though the DB-level change touches <6% of MSDB.
- **H1 has the largest reaction footprint of any variant** (226,695
  all-DB instances) and zero biology effect — `?` and `=` map to the
  same bounds. The number is a measure of how many core-model reactions
  are flagged "no rule fired" by the upstream cascade and would land on
  a curator's review queue if shipped.

---

## 4. Panel vs all-models FBA impact

The 100-model descriptive panel and the 5,683-model database give two
different lenses on the same heuristic change. Numbers below count models
whose biomass flux moves by more than `1e-6` (`flux-Δ`) or whose
grower/non-grower status flips (`flip`); both are evaluated against the
heuristic-baseline FBA (panel models rebound to the baseline cascade) so
that the comparison isolates the variant's effect.

| Variant | Panel flip (of 100) | Panel flux-Δ (of 100) | All-DB flip (of 5,683) | All-DB flux-Δ (of 5,683) | All-DB grew→not | All-DB →grew |
|---|---:|---:|---:|---:|---:|---:|
| 3.1 | **75** | 100 | **2,464** | 4,055 | 2,405 | 59 |
| 3.3 | 0 | 83 | 213 | 3,743 | 161 | 52 |
| 3.3_wide | 0 | 12 | 30 | 580 | 0 | 30 |
| 3.5 | 0 | 0 | 0 | 0 | 0 | 0 |
| 3.5_wide | 0 | 0 | 0 | 0 | 0 | 0 |
| 3.6 | 0 | 96 | 161 | 4,043 | 0 | 161 |
| 3.7 | 0 | 0 | 0 | 0 | 0 | 0 |
| 3.10_tight | 0 | 0 | 0 | 0 | 0 | 0 |
| 3.10_loose | 0 | 30 | 89 | 1,268 | 0 | 89 |
| H1 | 0 | 0 | 0 | 0 | 0 | 0 |
| H2 | 0 | 0 | 0 | 0 | 0 | 0 |
| H3 | **21** | 100 | **621** | 3,999 | 621 | 0 |
| H4 | **65** | 99 | **1,246** | 3,971 | 1,198 | 48 |

The all-DB column tracks per-variant mean biomass flux too (median, std,
and per-variant largest gainers / losers are in
`site/data/all_models_variants.json` and visible in the Variant Browser's
"all models" scope):

| Variant | All-DB mean flux | vs baseline (63.30) | All-DB median flux | All-DB std flux | All-DB n_grow (of 5683) | Largest gainer (Δ flux, model) | Largest loser (Δ flux, model) |
|---|---:|---|---:|---:|---:|---|---|
| baseline | 63.30 | — | 73.18 | 48.80 | 3,999 | — | — |
| 3.1 | 40.23 | **-23.07** | 0.00 | 66.81 | 1,653 | +114.50 (GCF_900116045.1) | -144.81 (GCF_002902865.1) |
| 3.3 | 49.49 | -13.81 | 52.32 | 40.58 | 3,890 | +94.48 (GCF_002849775.1) | -116.38 (GCF_002313045.1) |
| 3.3_wide | 64.08 | +0.78 | 74.18 | 48.75 | 4,029 | +120.16 (GCF_000767465.1) |  0.00 |
| 3.6 | **94.92** | **+31.62** | 111.99 | 68.72 | 4,160 | +177.28 (GCF_000022025.1) |  0.00 |
| 3.10_loose | 65.25 | +1.95 | 74.69 | 48.90 | 4,088 | +107.42 (GCF_001577265.1) |  0.00 |
| H3 | **15.18** | **-48.12** | 14.79 | 15.60 | 3,378 |  0.00 | -121.90 (GCF_003261575.2) |
| H4 | 52.70 | -10.60 | 22.12 | 56.56 | 2,849 | +122.49 (GCF_000264455.2) | -143.72 (GCF_002356295.1) |

(Variants with all-DB flip = 0 omitted — they match baseline exactly.)

---

## 5. Joint interpretation

**The descriptive panel under-represents directionalizing changes; the
all-models view captures their tail.** 3.6 (drop the low-energy-compounds
list) shows 0 panel grow-flips but 161 all-DB grow-flips — the LOW_ENERGY
rule was forcing direction on ~3,000 reactions, and removing it lets more
core models reach growth via newly-reversible cofactor flow. The panel
sees flux Δ on 96/100 models but no status changes because the panel was
chosen to be robust growers; the tail of the database is where you see
non-growers becoming growers (no `grew → not` for 3.6 anywhere). The site's
all-models scope surfaces these by listing largest gainers/losers per
variant.

**The two largest *directionalizing* changes (3.1 and H4) shrink the
grower set.** 3.1's reversibility index reclassifies ~3,300 reactions
including many that the heuristic had been treating as reversible; under
the tighter bounds, 2,405 / 5,683 models stop growing (and only 59 newly
do). H4 combines 3.1 + 3.5 + Bennett; the directional pressure adds up
the same way (1,198 lose growth, 48 gain). For users intending to use a
new heuristic as the *default*, these are the variants whose impact must
be understood before adoption — they will measurably reduce the count of
growing core models.

**H3 (the phosphates shadow-bug repair) is the single most disruptive
variant per character of code changed.** Three characters (the `cpd in
rgt` → `rgt['compound'] in PHOSPHATE_IDS` repair, or equivalent) flip
1,989 MSDB reactions and 621 core models from grower to non-grower.
This is a code-correctness fix, not a heuristic-design choice, and the
existing docstring in `Estimate_Reaction_Reversibility.py` describes
itself as preserving the bug for byte-for-byte parity. The all-models
impact column (the third-largest variant by all-DB grow-flip count)
quantifies what *should have* been happening all along under the
intended cascade.

**H1 is the inverse: the largest reaction footprint, zero biology
effect** — useful as a curation signal (which 6,522 reactions had no
rule fire?) but explicitly designed not to change any bound.

**3.7 and 3.10_tight are no-ops** at this resolution. The first is
unreachable due to a different shadow bug; the second tightens a band
no core-model reaction crosses.

**3.5 and 3.5_wide are MSDB-only changes** — measurable in the
MSDB direction-count table but invisible to any core-model FBA because
the affected reactions don't appear in core models. Their impact will
show up against a genome-scale model panel that includes the cofactor-
adjacent and secondary-metabolism reactions where CC σ is closest to
the band.

---

## 6. Recommendations by impact tier

### Adopt as code fix (no design choice)

- **H3** — `phosphates` shadow-bug repair. The cascade documentation
  itself describes this as latent-bug-preserved; the repair restores
  the intended ABC-transporter rule and is unambiguously a bug fix.
  Pair the merge with an update to the MSDB
  `THERMO_REFACTOR_CHANGES_REPORT.md` parity entry so consumers know
  the change occurred.
- **H2** — `LOW_LOCAL_CONC` shadow-bug repair. Identical bug pattern,
  near-zero impact (1 MSDB reaction flips, 0 core models), but pairs
  with H3 to restore the intended cascade end-to-end.

### Evaluate as default (large, principled biology change)

- **3.1** — Persist eQuilibrator's `ln_reversibility_index`. Largest
  reaction-direction footprint of any principled variant; replaces the
  ad-hoc points-and-bands heuristic with the formal quantity already
  computed upstream. ~30 lines to plumb through.
- **H4** — Composite of 3.1 + 3.5 + 3.3. The "best available evidence"
  stack. Its 1,198 grew→not all-DB count tells you the new default
  would shrink the growing-model set; that needs review against
  experimental panels before adoption.

### Adopt as curation surface

- **H1** — Distinguish "no rule fired" (`?`) from "agreed reversible"
  (`=`). Zero biology effect; large curation signal. Ship it as a
  metadata column rather than a bound override.

### Evaluate against a wider model population

- **3.5 / 3.5_wide** — CC σ band. Core-model-inert. Worth re-running
  the same pipeline against a genome-scale panel (e.g. iML1515,
  Recon3D-derived models) where the affected reactions live.

### Tertiary or contingent

- **3.3 / 3.3_wide** — concentration ranges. Real but smaller biology
  effect (3.3 has 161 grew→not, 3.3_wide has only 30 →grew). Useful as
  a Bennett 2009 sensitivity check rather than a default replacement.
- **3.6** — drop the low-energy-compounds list. Large positive biology
  shift (161 new growers, no losses) because relaxing forced direction
  opens cofactor flow. Pairs naturally with 3.1 (which restores principled
  directionality where it's warranted).
- **3.7** — disable CO₂ override. Currently a no-op (gated by H3).
  Re-evaluate after H3 lands.
- **3.10_tight / 3.10_loose** — band-width sensitivity. Educational
  knobs rather than recommendations.

---

## 7. Reproducing this report

```bash
cd core_models_analysis/
# 1. cascade outputs (depends on ModelSEEDDatabase + .kbcache built by notebook 06)
python3 scripts/export_thermo_variants.py
# 2. panel FBA aggregation (depends on cobra + .kbcache)
python3 scripts/build_site_data.py
# 3. all-models FBA aggregation (depends on cobra + core_models_kegg2/)
python3 scripts/build_all_models_impact.py --workers 64
# 4. browse interactively
python3 site/serve.py --port 8770 [--live]
```

Wall-clock on the workstation (336 cores, `--workers 64`):

- Step 1 (cascade): ~2 min for all 14 variants
- Step 2 (panel aggregation): ~30 s
- Step 3 (all-models FBA): ~3 min total (rxnsets index 0.5 s; baseline FBA 9 s; per-variant FBA 0–10 s each)

All-models variant FBA reuses cached files unless `--force-all-variants`
is passed.
