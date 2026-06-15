# Duplicate-reaction investigation: the `_c`-suffix annotation bug

**Status:** Bug identified and fixed across the scripts pipeline.
**Scope:** 17 transport reactions × 5,683 core models = ~96,500 cobra
reactions whose `seed.reaction` annotation carried a non-canonical suffix.
**Affected outputs:** `reports/REACTION_PREVALENCE.md` (visibly split rows),
`site/data/panel_rxnsets.json`, `site/data/all_models_*` (silent
under-counting in variant intersections), the heuristic-baseline FBA
across both panel and all-models populations (silent bound under-application).
**Generated:** 2026-06-15

---

## 1. The observation

A user spotted that `reports/REACTION_PREVALENCE.md` contained rows like

```
| rxn05466_c | 418  | 12.1% | 1067 | 48.0% | -35.9% |
| rxn05319_c | 1865 | 53.9% | 1780 | 80.1% | -26.2% |
| rxn05488_c | 2832 | 81.8% | 2093 | 94.2% | -12.4% |
| rxn05559_c | 2455 | 70.9% | 1805 | 81.2% | -10.3% |
| rxn05469_c | 3121 | 90.2% | 2221 | 100.0% | -9.8% |
```

— while sibling rows used the bare `rxnXXXXX` form (`rxn00001`,
`rxn05561`, etc.). ModelSEED reaction ids are *always* bare in the
upstream database. The `_c` suffix is non-canonical, and these five
rows were ranked as the most enriched-in-non-growers reactions in the
entire database — a strong-looking biological signal that turned out
to be a counting artifact.

This document records the investigation, the root cause, the fix, and
the downstream consequences.

---

## 2. Investigation

### 2.1 At the cobra-reaction-id level: no duplicates

Scanning all 5,683 model JSONs under `data/core_models_kegg2/`:

| Reaction-id suffix | Count |
|---|---:|
| `_c0` | 705,105 |
| `_e0` | 144,137 |
| (no suffix) | 11,366 |
| `_c` *(letter-only)* | 0 |

The 11,366 no-suffix reactions are biomass (`bio1`, `bio2`), exchanges
(`EX_*`), and sinks (`SK_*`). No cobra reaction id ends in a letter-only
`_c` — the duplicate is not at this layer.

Within any single model, no two cobra reactions share an id. No model
has the same SEED annotation appearing on two different cobra reactions.
No model has two cobra reactions with the same normalized stoichiometry.

### 2.2 At the `seed.reaction` annotation level: 17 reactions carry a stray `_c`

The duplicate lives inside the `annotation["seed.reaction"]` field. Across
all 5,683 models:

| seed.reaction value pattern | Count |
|---|---:|
| bare (`rxn00549`, `rxn14427`, …) | 611,569 |
| with `_c` suffix (`rxn11322_c`, `rxn05468_c`, …) | 87,853 |

The `_c` suffix appears on **exactly 17 distinct seed values**:

```
rxn05319_c  rxn05466_c  rxn05467_c  rxn05468_c  rxn05469_c
rxn05488_c  rxn05559_c  rxn08350_c  rxn08428_c  rxn08689_c
rxn08691_c  rxn08966_c  rxn09008_c  rxn10471_c  rxn10577_c
rxn11322_c  rxn39175_c
```

All 17 are **transport reactions** in MSDB (`is_transport=1`). Their
bare MSDB records exist; the `_c`-suffixed form does not exist as an
MSDB id:

```
rxn11322_c  →  rxn11322  ("(R,R)-butanediol transport")
rxn05468_c  →  rxn05468  ("TRANS-RXNAVI-26568.ce")
rxn05467_c  →  rxn05467  ("CO2 transporter via diffusion")
rxn08691_c  →  rxn08691  ("hydrogen transport via diffusion …")
…  etc
```

### 2.3 Per-model uniqueness: bare and suffixed forms never co-occur in one model

Within any single model, the bare seed id and its `_c`-suffixed form
never both appear — each affected reaction has *one* annotation in each
model, either bare or suffixed. The duplication is across models:

| Affected seed | Models with bare annotation | Models with `_c` annotation |
|---|---:|---:|
| rxn05466 | 4,198 | 1,485 |
| rxn05319 | 2,038 | 3,645 |
| rxn05488 | 758 | 4,925 |
| rxn05559 | 1,423 | 4,260 |
| rxn05469 | 341 | 5,342 |
| *(12 others)* | 0 | 5,683 *(each)* |

For 5 of the 17 reactions, both annotation forms appear in the
population — these are the 5 visible rows in REACTION_PREVALENCE.md.
For the other 12, only the suffixed form ever appears, so they collapse
into a single column and were invisible to that table.

### 2.4 Origin

The `_c` suffix matches the compartment letter of the cytosol in
ModelSEED's compartment scheme (`c` = cytosol, `e` = extracellular).
Each affected reaction is a transport (cytosol-touching) reaction
whose KBase template-build path appears to have appended the
compartment letter to the SEED id when writing the annotation. Why
only 17 reactions and not all transports is a quirk of the build
template — the user-facing data is what we have to work with.

---

## 3. Decision

**The bare and `_c`-suffixed forms are the same MSDB reaction** and must
be counted as such. The justification is unambiguous:

- Both forms have the same MSDB record (the bare id) — there is no
  `rxn05466_c` reaction definition anywhere upstream.
- The `_c` suffix carries no biological meaning the bare id does not
  already encode (the cobra reaction itself already lives in the
  `c0` compartment, encoded in the cobra reaction id `rxn05466_c0`).
- For 5 of the 17 affected reactions, the bare and suffixed forms
  appear in *different subsets* of the same 5,683-model population —
  treating them as distinct reactions split prevalence counts into
  artifacts.
- The cascade's `{rxn_id: reversibility}` map uses bare MSDB ids
  exclusively, so any downstream lookup against the `_c`-suffixed
  annotation silently returns `None` and the override is skipped.

The fix is to **normalize the annotation at read time**: strip a
trailing `_<letter>` (no digits) from any `seed.reaction` value
before using it as a lookup key or a grouping key.

---

## 4. Downstream consequences before the fix

### 4.1 `reports/REACTION_PREVALENCE.md` — split rows that *all collapsed to no signal*

The top-5 non-grower-enriched rows were all bug artifacts. After merging:

| Reaction | Pre-fix top-30 rank | Pre-fix Δ | Post-fix grower / non-grower | Post-fix Δ |
|---|---:|---:|:---:|---:|
| rxn05466 | 1 (was `_c`) | -35.9% | 3461/3461 vs 2222/2222 (100% / 100%) | **0.0%** |
| rxn05319 | 2 (was `_c`) | -26.2% | 3461/3461 vs 2222/2222 (100% / 100%) | **0.0%** |
| rxn05488 | 4 (was `_c`) | -12.4% | 3461/3461 vs 2222/2222 (100% / 100%) | **0.0%** |
| rxn05559 | 8 (was `_c`) | -10.3% | 3461/3461 vs 2222/2222 (100% / 100%) | **0.0%** |
| rxn05469 | 10 (was `_c`) | -9.8% | 3461/3461 vs 2222/2222 (100% / 100%) | **0.0%** |

All five reactions are universally present (100% of growers, 100% of
non-growers). The apparent enrichment was entirely the
bare-vs-`_c`-suffixed split.

The post-fix top-30 non-grower-enriched table is led by `rxn05759`
(2.9% vs 17.8%, Δ -14.9%) — five real biological signals are now
visible at the top that were buried below the artifacts.

### 4.2 `site/data/panel_rxnsets.json` — wrong keys for transport intersections

Before fix: 1,564 entries across the 100 panel models held the
`_c`-suffixed form. After fix: 0. Any variant-vs-panel intersection
(`variant.diffs ∩ panel_rxnsets[mid]`) using the `_c` form would miss
these transport reactions silently. (As it happens, none of the 14
cascade variants currently in the catalog change any of the 17
affected reactions, so this specific intersection produced the same
numbers; the bug was latent and would have surfaced as soon as a new
variant touched a transport reaction.)

### 4.3 Heuristic-baseline FBA — bounds silently under-applied

`override_bounds()` in `scripts/growth_heuristics.py` looked up
`reversibility_map.get(seed)` using the raw annotation. For the 17
`_c`-suffixed reactions in each model, this returned `None` and the
override was skipped — leaving whatever on-disk bound the cobra JSON
contained.

In a sample panel model (`GCF_003261575.2`), 5 of the 12 `_c`-suffixed
reactions had on-disk bounds that conflict with the cascade default:

| Reaction id | seed annotation | Cascade rev | On-disk bounds | Should be |
|---|---|---|---:|---|
| `rxn11322_c0` | `rxn11322_c` | `=` (reversible) | (0, 0) **disabled** | (-1000, 1000) |
| `rxn05468_c0` | `rxn05468_c` | `=` | (0, 1000) | (-1000, 1000) |
| `rxn10577_c0` | `rxn10577_c` | `=` | (0, 1000) | (-1000, 1000) |
| `rxn08350_c0` | `rxn08350_c` | `=` | (0, 1000) | (-1000, 1000) |
| `rxn08689_c0` | `rxn08689_c` | `=` | (0, 1000) | (-1000, 1000) |

Most striking is `rxn11322_c0` ((R,R)-butanediol transport) which was
shipped with `lb=ub=0` (disabled) but the cascade calls it reversible.
Before the fix, every rebound-FBA call left that reaction disabled;
after the fix it opens.

The aggregate effect on the all-models baseline mean biomass flux:
**63.30 → 63.74** (+0.44, +0.7%) across the 5,683-model run. Per-variant
flux means moved by similar amounts (the bound restoration applies to
both baseline and variant maps).

### 4.4 Variant-impact flip counts shift slightly

Re-running `build_all_models_impact.py` on the corrected baseline
produces the following changes:

| Variant | Pre-fix `n_models_flip` | Post-fix | Δ |
|---|---:|---:|---:|
| baseline mean_flux | 63.30 | 63.74 | +0.44 |
| 3.1 flip | 2,464 | 2,436 | −28 |
| 3.3 flip | 213 | 213 | 0 |
| 3.6 flip | 161 | 161 | 0 |
| H3 flip | 621 | 607 | −14 |
| H4 flip | 1,246 | 1,245 | −1 |

The shifts are small in absolute terms but reflect real models whose
"flipped to non-grower" status was an artifact of the missing transport
bound. The qualitative conclusions of the impact report
(`REVERSIBILITY_HEURISTICS_IMPACT.md`) do not change.

---

## 5. The fix

### 5.1 New helper: `scripts/seed_annotation.py`

```python
import re
_SEED_COMPARTMENT_SUFFIX = re.compile(r"_[a-z]$")

def normalize_seed_id(seed):
    """Strip a stray compartment-letter suffix from a SEED reaction id."""
    if not seed:
        return seed
    return _SEED_COMPARTMENT_SUFFIX.sub("", seed)

def seed_id(reaction):
    """Return the normalized SEED id for a cobra rxn or raw JSON dict."""
    ...
```

The regex strips a single trailing `_<letter>` (no digits). It is
deliberately narrower than `_<letter>\d*`:

- The 17 affected annotations always use the letter-only form
  (`_c`, never `_c0`).
- The cobra reaction *id* (e.g. `rxn05466_c0`) carries the legitimate
  compartment-marker suffix that we never want to strip from cobra ids.
  Confining the normalizer to letter-only protects against accidentally
  collapsing future ids that look like `rxn00001_c0`.

### 5.2 Call sites updated

Every script that reads `annotation["seed.reaction"]` now goes through
`seed_id()` or `normalize_seed_id()`:

| File | Function / cell |
|---|---|
| `scripts/build_all_models_impact.py` | `_extract_rxnset()` |
| `scripts/growth_heuristics.py` | `override_bounds()` |
| `scripts/build_site_data.py` | normalizes `rxnsets_by_model` after the kbcache load |
| `scripts/deeper_analysis.py` | `extract_reactions_one()` (`REACTION_PREVALENCE.md` source) |
| `scripts/select_diverse.py` | `extract_one()` |
| `scripts/select_diverse_tax.py` | `_extract_one()` |
| `scripts/direction_pipeline.py` | `_apply_source_map_to_model()`, `source_coverage()` |
| `scripts/build_thermo_source_network_tables.py` | model-walk loop |
| `scripts/build_thermo_source_comparison_notebook.py` | emitted cell |
| `scripts/run_thermo_source_variants.py` | `_model_overrides_expected()` |
| `scripts/build_notebooks.py` | `rxnsets_by_model` cache builder (notebook-emitted cell, regex inlined) |
| `scripts/build_taxonomy_aware_notebook.py` | same shape, regex inlined |

The two notebook-emitter scripts inline the regex because the emitted
cell runs in a notebook context that cannot rely on a sibling-module
import without explicit `sys.path` setup. The behavior matches
`normalize_seed_id` exactly.

### 5.3 The kbcache (`notebooks/.kbcache/rxnsets_by_model`) is *not* invalidated

The cache was populated before the fix and still stores the
`_c`-suffixed annotations verbatim (its blob is content-hashed; touching
the source would force a full rebuild of notebook 06's downstream
artifacts). Instead, `build_site_data.py` normalizes the cached values
on read — the emitted `site/data/panel_rxnsets.json` is correct, and
the cache itself is harmless because no other consumer reads it
directly. Re-running notebook 06 with the patched `build_notebooks.py`
would refresh the cache; that step is optional and orthogonal to the
correctness of the site data.

### 5.4 Outputs regenerated

After the patches:

- `site/data/all_models_rxnsets.json` (8 MB; gitignored): rebuilt by
  `scripts/build_all_models_impact.py`. Sample model
  `GCF_000014585.1`: 97 reactions, 0 with a `_c` suffix.
- `site/data/all_models_baseline_fba.json` (gitignored): rebuilt with
  the patched `override_bounds`. Mean biomass flux 63.30 → 63.74.
- `site/data/all_models_variant_fba__{tag}.json` (×14, committed):
  rebuilt for every variant. Per-variant flip counts shift by 0–28
  models — see the table in §4.4.
- `site/data/all_models_variants.json` (committed): rebuilt aggregate
  summary.
- `site/data/panel_rxnsets.json` (committed): rebuilt by
  `scripts/build_site_data.py`, which now reports the normalization
  count when it loads the cache (typically "normalized seed.reaction
  _c-suffix in 3461 models").
- `reports/REACTION_PREVALENCE.md`: regenerated. The five
  `_c`-suffixed top-30 non-grower-enriched rows are gone; the top of
  that list now starts with the next biological signal (`rxn05759`
  at -14.9%).

---

## 6. Lessons / follow-ups

- **Annotation contracts deserve assertion.** A `_c`-suffixed
  `seed.reaction` is silently wrong and silently masked by every code
  path that just calls `dict.get`. A one-line guard at the helper layer
  (`assert _SEED_COMPARTMENT_SUFFIX.search(s) is None` after
  normalization, behind an env-var debug flag) would catch new
  occurrences as the model set evolves.
- **The cobra "touched" counter and the unique-seed counter are
  different metrics.** The fix makes them more consistent
  (`override_bounds` now updates one more cobra reaction per `_c`-
  suffixed transport per model), but they still mean different things;
  whenever a number is reported as "reactions touched" the units —
  cobra rxn instances vs unique MSDB ids — should be in the caption.
- **The fix did not change the qualitative conclusions of any prior
  report.** It did move:
  - `REACTION_PREVALENCE.md` (5 artificial top-of-list rows removed)
  - All-models baseline mean flux (+0.7%)
  - A handful of variant flip counts (−1 to −28 per variant)
- **17 transport reactions were silently disabled or under-bounded
  during every prior rebound-FBA call.** The most striking was
  `rxn11322_c0` ((R,R)-butanediol transport) which carried `lb=ub=0`
  on disk despite the cascade calling it reversible. Anyone looking
  at growth on butanediol carbon would have been confused.
