# A gold / silver / bronze grade for each thermodynamic source, per reaction

**Status: shipped — this page is kept as the design record.** Numbers below are
from the working prototype run 2026-08-12 against `/scratch/ctaylor/tmp/devsnap2`
(ModelSEED `dev` @ 49563c6f; Group Contribution = the Convention A rebuild) and
differ from the shipped run by a few tens of reactions, mostly because the
MetaNetX-collision veto grew from 15 to 35 reactions once
`results/eq_vs_dgpms/reconciliation.tsv` was wired in. **For current numbers use
[`../thermoSourceMethod/THERMO_SOURCE_METHOD.md`](../thermoSourceMethod/THERMO_SOURCE_METHOD.md)**;
read this page for the *rationale* — where the problem came from (§0), why
corroboration is asymmetric (§2.3), and which simpler schemes were rejected (§4).

Companion to [`THERMO_SOURCE_ASSIGNMENT.md`](THERMO_SOURCE_ASSIGNMENT.md), which
picks *one winning source* per reaction. This does not pick a winner. It hands
back a **grade card**: for one ModelSEED reaction, every source that has a value
gets its own independent gold / silver / bronze label.

```
rxn00001   diphosphate phosphohydrolase                        EC 3.6.1.1

  source            ΔG′°     σ      p_ok    z     grade    reason
  TECRDB            —        —      —       —     —        not measured
  eQuilibrator      -4.07    0.28   0.96    0.31  GOLD     self-certain
  dGPredictor-MS    -3.12    1.90   0.88    0.62  SILVER   self-confident
  Group contrib.    -6.81    9.40   0.55    1.41  SILVER   corroborated
```

Four sources, four independent verdicts. TECRDB is graded GOLD wherever it
exists, because it is a measurement rather than a prediction.

---

## 0. Where this problem came from

The distinction the design turns on — a source's **reported** error (σ) versus
its **true** error (ε) — entered this project in four steps:

| when | where | what happened |
|---|---|---|
| **2026-07-21** | dgpredictor↔TECRDB work | First time the two were set against each other: *"exclude any reactions where the combined errors are larger than the difference between the energies predicted by DG and TECR."* |
| **2026-08-06** | eQ vs dGPredictor-ModelSEED | *"What is sigma here?"* → the retrained model's σ is **calibrated** (ρ = +0.672 vs \|eQ−dGP\|) where the original's was flat and useless; large disagreements can be filtered by self-reported confidence. |
| **2026-08-07 → 08-10** | `optimize_thermo_source_assignment.py` | The explicit ask: *"how exactly we are gathering the predictions between the error of the model, and how close it is to the true error?"* → ê, the isotonic bridge σ → E[ε], and §5b of the assignment report. |
| **2026-08-11** | σ/ê-filtered scatters | *"is sigma the model's error, or the energy?"* — plus the finding that Group Contribution's σ is only ρ = +0.176 against measured error. |

**§5b of the assignment report is why this proposal exists.** It found ê is a
good *threshold* (ê ≤ 2 ⇒ true error ≤ 2 about 87–91% of the time) but a poor
*ranking*: ρ(ê, actual) = 0.11–0.24, and dGPredictor's accepted vs rejected
groups differ by only 0.46 vs 0.60 kcal/mol. A three-bin grade cut straight from
ê would be a grade in name only. Worse, on the Convention A GC rebuild ê spans
just 3.04 → 5.70 kcal/mol across the entire database, so GC would be unrankable
and permanently non-gold for the wrong reason.

### Naming collision, resolved

`optimize_thermo_source_assignment.py` already uses "gold" and "silver" for its
two **calibration-data tiers** (TECRDB measurements vs the trusted-σ proxy
target). Those are fitting inputs. This document calls them **anchor** and
**proxy**, and reserves gold / silver / bronze for the emitted grade.

---

## 1. The four sources

| source | what it is | grading basis |
|---|---|---|
| **TECRDB** | NIST experimental ΔG′° = −RT ln K′, matched to ModelSEED by the SMILES→InChIKey multiset pipeline in `/scratch/ctaylor/dgpredictor_tecrdb` | **Always GOLD.** It is data, not a model. |
| **eQuilibrator** | component contribution, σ from **ν**ᵀ**Σν** | σ → calibrated, plus cross-source behaviour |
| **dGPredictor-ModelSEED** | BayesianRidge on radius-1/2 group changes, σ = posterior predictive sd | same |
| **Group contribution** | fitted group energies, propagated σ | same |

**One deliberate deviation on TECRDB.** The match to a ModelSEED reaction comes
in two tiers: `stereo_exact` (802 reactions, full InChIKey — distinguishes
anomers and D/L pairs) and `skeleton` (748 more, connectivity block only). The
*measurement* is gold in both, but a skeleton match may have attached the wrong
reaction's data — this is how hexokinase/aldose cases got conflated to glucose.
Proposal: `stereo_exact` → GOLD, `skeleton` → SILVER, with a
`--tecrdb-skeleton-gold` flag if you want the simpler rule. Nothing about
TECRDB is ever BRONZE.

TECRDB's own experimental scatter is small enough not to warrant further
grading: median sd 0.15 kcal/mol, p90 1.04, max 3.15. It is carried in the
output as `tecrdb_sd` and `n_measurements` (551 of the 1,550 rest on a single
measurement) so a consumer can apply its own floor.

---

## 2. What grades the three predictors

Each predictor is graded on its own merits. Three signals feed it, and they fail
in different ways — which is the point.

### 2.1 Its own calibrated confidence — `p_ok` (the primary signal)

Replace ê (an expected *magnitude*) with a calibrated *probability*:

> **p_ok(i, s) = P̂( ε_s(i) ≤ τ | σ_s(i) )**,  τ = 2.0 kcal/mol

τ = 2.0 is the reversible band the direction cascade actually uses
(`reversibility_heuristics.py:327`), so p_ok reads as *"the probability this
number is good enough to call the direction."*

Fitted exactly like ĝ_s in the assignment report — per-source isotonic
regression, decreasing, two-tier anchor (802 TECRDB stereo-exact, weight 3) +
proxy (|source − trusted-σ reference|, weight 1) — but on the indicator
`1[ε ≤ τ]` rather than on |ε|. Three reasons this is the better primitive for
grading:

1. **Comparable across sources.** ê inherits each σ's own scale, and §3.3 of the
   assignment report shows those scales differ by up to 3.4× in optimism (GC
   overstates error 2.2×, dGPredictor-MS 1.5×, eQuilibrator *understates* 1.6×).
   A probability does not have that problem, so "GOLD" means the same thing in
   the eQuilibrator column as in the GC column.
2. **Does not saturate.** GC's ê is effectively a constant; GC's p_ok spans
   0.85 → 0.51 — weak, but ordered and usable.
3. **Honest direction.** ê is conservative by ~4×; an isotonic-calibrated
   probability is unbiased by construction.

Fitted curves (prototype):

| source | anchor n | proxy n | p_ok at min σ | at median σ | at max σ | anchor frac(ε ≤ 2) |
|---|---:|---:|---:|---:|---:|---:|
| eQuilibrator | 794 | 4,011 | 0.999 | 0.811 | 0.000 | 0.856 |
| dGPredictor-MS | 802 | 11,183 | 0.988 | 0.383 | 0.000 | 0.848 |
| Group contribution | 802 | 10,025 | 0.845 | 0.522 | 0.511 | 0.566 |

GC's near-flatness is the honest output, not a bug: its σ carries real but weak
information (ρ = +0.176 against measured error). **GC therefore never reaches
GOLD on its own σ.** That is a finding, not a limitation of the scheme.

### 2.2 Direct measurement — TECRDB (overrides everything)

Where the reaction is in the stereo-exact tier, ε_s is *observed* for that
source, and nothing else needs to be inferred:

> ε_s ≤ 1 → **GOLD** ε_s ≤ 3 → **SILVER** ε_s > 3 → **BRONZE**

This is the most per-source signal there is — each predictor is scored against
the experiment individually, so on the same reaction eQuilibrator can be GOLD
while GC is BRONZE. 802 reactions × 3 sources.

### 2.3 Cross-source behaviour — used *asymmetrically*

Treat the sources as independent measurements of the same quantity, weighted by
their **calibrated** error scale ê_s (not raw σ — that is what makes them
commensurable):

```
w_s      = 1 / ê_s²
ΔG_fused = Σ w_s ΔG_s / Σ w_s
χ²       = Σ w_s (ΔG_s − ΔG_fused)²           df = n_src − 1
R        = sqrt( χ² / df )                     Birge ratio / PDG scale factor
z_s      = |ΔG_s − ΔG_fused| / ê_s             per-source residual  ← the per-source part
```

R ≈ 1 means the spread among sources is exactly what their stated uncertainties
predict; R ≫ 1 means at least one is wrong, and **z_s says which**. ΔG_fused is
an internal construct here — a reference point for computing z_s — not an output.

Measured against the 802 TECRDB reactions, R ranks true error where ê could not:

| Birge R | n | median \|fused − experiment\| | frac ≤ 2 |
|---|---:|---:|---:|
| R ≤ 1 | 568 | **0.36** | 93% |
| 1 < R ≤ 2 | 157 | 0.64 | 84% |
| 2 < R ≤ 5 | 59 | **3.18** | 42% |
| R > 5 | 11 | **5.69** | **0%** |

A 16× monotone spread, versus 0.46-vs-0.60 for the ê accept/reject split.

**But agreement and disagreement are not symmetric evidence, and the prototype
shows it.** Two fallible predictors agreeing is weak — they can be wrong the
same way (eQuilibrator and Group Contribution share group-contribution lineage),
and 11% of the R ≤ 1 set are structural zeros where agreement is imposed by the
stoichiometry rather than earned. Two predictors disagreeing is strong: someone
is definitely wrong, and z_s names them. So:

> **Corroboration can raise a BRONZE to SILVER. It can never create a GOLD.
> Being outvoted demotes one tier.**

Letting corroboration promote all the way to GOLD was tested and rejected: it
grew eQuilibrator's GOLD column from 2,443 to 9,157 but diluted its measured
guarantee from 94% to 90% within 2 kcal/mol, and dGPredictor's from 98% to 91%.
The asymmetric rule keeps GOLD pure and still gets the full benefit of the
demotion, which is where the discriminating power turned out to live.

---

## 3. The per-source cascade

Applied independently to each of GC, EQ, DGPMS. Most-specific-first, one `reason`
slug per row, so every label is auditable back to the rule that made it — the
same discipline as the reversibility cascade.

```
0. UNGRADED   no ΔG_s stored for this reaction
              OR eQuilibrator sentinel, σ > 100 (the source disclaims it)   4,934
              OR eQuilibrator MetaNetX collision (accumulate-vs-assign bug)    15
              OR dGPredictor-MS on a quinone/quinol couple, 52.8% sign-wrong  511
              OR legacy KEGG-keyed dGPredictor inside the 17,271-entry
                 mis-map mask, if that source is ever admitted

1. measured   reaction in TECRDB stereo_exact:
                 |ΔG_s − ΔG*| ≤ 1 → GOLD | ≤ 3 → SILVER | > 3 → BRONZE     (terminal)

2. base grade from the source's own calibrated confidence
                 p_ok ≥ 0.90 → GOLD    "self-certain"
                 p_ok ≥ 0.70 → SILVER  "self-confident"
                 else        → BRONZE  "uncorroborated"

3. floor      BRONZE and n_src ≥ 2 and R ≤ 1.5 and z_s ≤ 1
                 and not a structural zero (all-|ΔG| < 0.5, or transport)
                             → SILVER  "corroborated"      (never above SILVER)

4. demote     n_src ≥ 2 and R > 2 and z_s > 3
                             → one tier down   "outvoted"
```

### Validation — grades inferred, then scored against data they could not see

Rule 1 disabled, so the grade comes only from p_ok and cross-source behaviour;
then scored against the TECRDB measurement:

| source | grade | n | median \|ε\| | frac ≤ 2 | frac ≤ 1 |
|---|---|---:|---:|---:|---:|
| eQuilibrator | GOLD | 246 | **0.32** | 94% | 87% |
| | SILVER | 530 | 0.47 | 85% | 68% |
| | BRONZE | 18 | **3.47** | **0%** | 0% |
| dGPredictor-MS | GOLD | 184 | **0.32** | 98% | 78% |
| | SILVER | 604 | 0.55 | 82% | 64% |
| | BRONZE | 14 | **18.50** | 29% | 29% |
| Group contribution | GOLD | 0 | — | — | — |
| | SILVER | 514 | 1.28 | 69% | 41% |
| | BRONZE | 288 | **8.60** | 34% | 33% |

Monotone in every column, for every source, on data withheld from the grade.
The BRONZE row is what the asymmetric rule buys: without the corroboration
layer, GC's BRONZE median is 1.66 kcal/mol — indistinguishable from its SILVER —
and with it, 8.60.

### Database-wide result

| source | GOLD | SILVER | BRONZE | UNGRADED |
|---|---:|---:|---:|---:|
| TECRDB | 802 | 748 *(skeleton match)* | 0 | 54,452 |
| eQuilibrator | 2,443 | 13,526 | 4,125 | 35,908 |
| dGPredictor-MS | 5,808 | 12,420 | 13,185 | 24,589 |
| Group contribution | 309 | 12,998 | 14,006 | 28,689 |

Reason breakdown, e.g. dGPredictor-MS: 5,271 self-certain, 2,094 self-confident,
10,139 corroborated, 11,159 uncorroborated, 1,948 outvoted, 802 measured.

Per reaction: 33,337 have at least one source; **6,671 have at least one GOLD**,
23,090 at least one GOLD-or-SILVER, 10,247 are all-BRONZE, and 22,665 have no
thermodynamic source at all — the coverage ceiling, which no grading scheme
changes.

---

## 4. Why not the simpler alternatives

| alternative | why not |
|---|---|
| Three bins cut from ê | §5b: ρ(ê, actual) = 0.11–0.24, and Convention A GC's ê spans only 3.04–5.70, so GC is unrankable and non-gold for the wrong reason. |
| Three bins cut from raw σ | The three σ scales are not comparable (2.2× / 1.5× pessimistic, 1.6× optimistic), so grading on σ rewards pessimistic self-reporting. This is the exact failure that makes dev's `Promote_Reaction_Thermodynamics_to_Canonical.py` prefer the mis-mapped legacy dGPredictor 95.3% of the time. |
| Majority vote across sources | Discards the magnitude of disagreement that R uses, and assumes an independence that EQ and GC do not have. |
| One grade per reaction | The whole point of the correction: on rxn00001 eQuilibrator deserves GOLD and GC does not. Collapsing loses that. A reaction-level roll-up is still emitted, but as `best_grade`, derived, not primary. |
| Supervised classifier on features | 802 labels, all central metabolism, and no per-reaction audit trail. The `reason` column is worth more here than a few points of AUC. |

---

## 5. Honest limitations

1. **The anchor set is central metabolism.** TECRDB contains almost none of the
   exotic high-σ chemistry BRONZE exists to catch, so the BRONZE rows above rest
   on n = 14–288, and GOLD's guarantee is demonstrated on the easy half of the
   database. Unchanged from §5b; the cross-source layer does not remove it, it
   only adds a signal that demonstrably ranks *within* the anchor set.
2. **EQ and GC are not independent**, so R under-reports discrepancy on that
   pair. Mitigation: emit the pairwise R matrix alongside the pooled R and
   down-weight EQ–GC agreement in rule 3.
3. **p_ok leans heavily on the proxy tier.** The proxy target is a stand-in, not
   a measurement — measured against the anchor it is approximately unbiased
   (median difference ±0.01 kcal/mol) but noisier (ρ = 0.43–0.84) — and it
   supplies 63–82% of the fitting weight. Group Contribution and
   dGPredictor-ModelSEED are both calibrated against eQuilibrator, so a
   systematic error in its low-σ regime would reach two of the three curves.
   Carried per row as `n_anchor` / `n_proxy` provenance.
4. **Grades are snapshot-specific.** The Convention A rebuild changed 53% of GC
   values and doubled its σ. Every emitted table stamps the MSDB commit.
5. **GC has no GOLD path except measurement.** By construction and by evidence —
   but it means GC's 309 GOLD rows are exactly its TECRDB hits, and any future
   improvement to GC's σ reporting would need a recalibration to show up.

---

## 6. Implementation sketch

New `scripts/grade_thermo_sources.py` (~250 lines), importing `load_db`,
`fit_error_models` and the veto predicates from
`optimize_thermo_source_assignment.py` rather than re-deriving them.

```
results/thermo_grades/
  source_grades.tsv        long, the primary artifact — one row per
                           (rxn × source), 4 sources: dg, sigma, ehat, p_ok,
                           z, birge, n_src, grade, reason, provenance
  source_grades_wide.tsv   one row per rxn, one grade column per source,
                           + best_grade convenience column
  grade_calibration.json   fitted p_ok curves + the validation tables above
  grade_frontier.tsv       threshold sweep (p_ok cuts, R cuts)
```

Consumers: `load_grades()` mirroring `load_assignment()`; a `/api/grades/<rxn>`
endpoint plus a four-row grade card in `site/serve.py`, which already renders
per-reaction thermodynamics; and an optional grade floor on the direction
cascade — *"only overturn a stored bound on GOLD or SILVER thermodynamics."*

One pass over the snapshot, ~3 minutes.
