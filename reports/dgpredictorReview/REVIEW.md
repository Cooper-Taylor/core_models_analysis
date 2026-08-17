# Review — `freiburgermsu/dGPredictor` (ModelSEED fine-tuning)

Repo reviewed at `master` @ `5d656c5` ("Refresh comparisons against July retrain
and add comparative figures"), cloned to `/scratch/ctaylor/dgpredictor_repo`.
It is a fork of `maranasgroup/dGPredictor` with an added, parallel
ModelSEED-native pipeline. **No original file was modified** — the KEGG path
still works via `dg_prediction.py`; the new path is everything prefixed
`modelseed_`.

Everything numeric below was recomputed locally; scripts in `scripts/`, outputs
in `results/`, figures in `figures/` and (for the scatter) in
`core_models_analysis/reports/thermoComparison/figures/dgpredictor_base_vs_finetuned/`.

---

## 1. What is in the repo

| file | role |
|---|---|
| `retrain_modelseed.py` (949 lines) | the whole fine-tuning pipeline |
| `dg_prediction_modelseed.py` | `ModelSEEDdGPredictor` — inference by `rxn#####`, equation, or SMILES |
| `predict_all_modelseed.py` | bulk prediction over all ModelSEED reactions |
| `build_dG_comparison.py`, `make_comparison_figures.py` | old-vs-new comparison |
| `CC/` | vendored **component-contribution** (Noor/Milo) — `Compound`, pKa Legendre transform |
| `CHANGES_MODELSEED_RETRAINING.txt` | the authors' own change log; accurate |
| `REPORT.md` | unrelated: OPAM2-vs-Marvin pKa study that reuses this model |

**Not in the repo** (`.gitignore` = `*/*.pkl`, `*/*.mat`):
`model/modelseed_M12_model_BR.pkl` and `data/modelseed_training.mat`
(only the 1.3 MB `.zip` of the 1.2 GB mat ships). The fine-tuned model cannot be
loaded from a fresh clone — it has to be rebuilt.

I rebuilt it. `scripts/rebuild_finetuned_features.py` replays steps 4–5 of
`retrain_modelseed.py` sparsely and refits `BayesianRidge`; the result reproduces
the shipped `dG_model_only` for **all 32,001** reactions to a max
absolute difference of **5.8 × 10⁻⁷ kJ/mol** (`results/refit_validation.json`).
So the numbers below are the actual model, not a proxy.

---

## 2. How the fine-tuning works

It is **not** fine-tuning in the deep-learning sense. There is no warm start
from the old weights, no new experimental data, and no gradient continuation.
It is a **re-fit of the same regression on the same measurements, over a wider
structural alphabet.**

```
ModelSEED dev  ──► 45,706 compounds ──► 32,647 with SMILES
                                    ──► 27,554 with COMPLETE SMILES  (5,093 dropped for '*')
                        │
                        ├─ compound cache: pKa/pKb parsed from ModelSEED's own
                        │  columns (ChemAxon bypassed), InChI from
                        │  Unique_ModelSEED_Structures.txt, 12 special compounds
                        │  hand-defined (H+, H2, S, CO, HCO3-, Ca/K/Mg/Fe, Fd_ox/red)
                        │
                        └─ RDKit atom-environment decomposition at radius 1 and 2
                           ──► vocabulary 2,515 (r1) + 37,846 (r2) = 40,361 fragments
                               (was 1,435 + 24,882 = 26,317)

TRAINING TARGETS  ──  data/component_contribution_python.mat  (UNCHANGED)
                      train_S (673 KEGG cpds × 4,001 rxns), b (4,001 ΔG°, kJ/mol)
                        │
                        ├─ KEGG cpd ids → ModelSEED cpd ids via
                        │  Unique_ModelSEED_Compound_Aliases.txt   (67 rxns lost)
                        └─ 24 more lost to undecomposable compounds
                           ──► 3,910 training reactions

X = [Δgroups_r1 | 44 zeros | Δgroups_r2 | 44 zeros]        (3,910 × 40,449)
BayesianRidge(tol=1e-6, fit_intercept=False)  →  MSE 9.88 kJ²/mol², R² 0.9998

predict:  ΔG′° = X·coef + Σ ν_i · transform_pH7(cpd_i, pH 7.0, I 0.25 M, 298.15 K)
          →  31,924 of 55,999 ModelSEED reactions (57%)
```

What actually changed, then:

* **Compound namespace.** cpd##### native, so coverage is no longer gated on a
  compound having a KEGG id. This is the real win: **+11,352 reactions** the
  original could not touch, and no more KEGG-id inference.
* **Structural alphabet.** 26,317 → 40,361 declared fragments.
* **pKa source.** ModelSEED's stored ChemAxon values instead of live `cxcalc`,
  so the pipeline has no proprietary binary dependency.
* **Not changed:** the estimator, its hyperparameters, the feature layout
  (including the two vestigial 44-zero pads), and — critically — **the
  measurements**.

---

## 3. Is any eQuilibrator data used? — *Depends what you mean, and the answer matters*

**No eQuilibrator estimate, cache, or API call enters the pipeline.** Verified:

* `retrain_modelseed.py` and `dg_prediction_modelseed.py` import only
  `CC.compound.Compound` and `CC.thermodynamic_constants`. Neither imports
  `CompoundCacher`, which is the only thing in the repo that reads
  `CC/data_cc/equilibrator_compounds.json.gz`.
* The only other `equilibrator` references anywhere in the Python/notebook code
  are two **commented-out** lines in `db_bulk_dg_gen.py` (a MetaNetX energy table
  used for a comparison, not for training).
* No MetaNetX id crosswalk, no `equilibrator_api`, no `equilibrator_cache`.

**But the training targets are the eQuilibrator/component-contribution corpus.**
`data/component_contribution_python.mat` is a saved state of the
component-contribution model of Noor, Haraldsdóttir, Milo & Fleming (2013) — the
engine underneath eQuilibrator 2.x. It contains `dG0_cc`, `dG0_rc`, `dG0_gc`,
the 673 × 241 group matrix `G`, and the projection matrices `P_R_gc` / `P_N_gc`.
Checks I ran:

| check | result |
|---|---|
| `b` reproduced by CC's own formation energies, `Sᵀ · dG0_cc` vs `b` | **r = 0.9999**, median \|Δ\| **0.95 kJ/mol** over all 4,001 |
| training columns whose stoichiometry matches a TECRDB reaction exactly | **3,741 / 4,001** |
| remaining columns | 223 single-compound Alberty formation-energy pseudo-reactions + redox |
| `CC/data_cc/TECRDB.tsv` | 4,544 rows, the Milo-lab NIST curation shipped with eQuilibrator |

So: the fine-tuned model is **independent of eQuilibrator's predictions** and
**not independent of eQuilibrator's data**. It is trained on exactly the same
4,001 measurements eQuilibrator is trained on, reverse-transformed by the same
code. Practically, that means:

* Using it as an *independent* cross-check against eQuilibrator is not sound —
  agreement is partly guaranteed by shared measurements, and disagreement is
  purely a difference of functional form.
* Benchmarking either against TECRDB is **in-sample for both**.
* The repo's own `claude_assessment.txt` (line 101) proposed "use the
  eQuilibrator experimental training data … to augment dGPredictor's training
  set". That did not happen as a new step — it was already the training set.

---

## 4. Base vs fine-tuned scatter

Script: `core_models_analysis/scripts/plot_dgpredictor_base_vs_finetuned.py`
(house style of `plot_thermo_source_dg_scatter.py` — same reversibility-transition
colouring, same palette, same cascade). Data: ModelSEED `dev` @ `49563c6f`
(`/scratch/ctaylor/tmp/devsnap2`), where both variants sit in the same
`thermodynamics` dict, so no unit conversion or crosswalk is involved.

Coverage: original 27,715 · retrain 31,924 · **both 20,567** · original-only 7,143 ·
**retrain-only 11,352**.

| panel | n | Pearson r | Spearman ρ | median \|Δ\| | sign flips |
|---|---:|---:|---:|---:|---:|
| all co-covered | 20,567 | **0.27** | 0.44 | 7.67 kcal/mol | 22.2% |
| KEGG id vouched by a MSDB alias | 10,185 | **0.83** | 0.82 | 3.35 | 15.4% |
| KEGG id inferred (mis-mapped) | 10,382 | **0.0005** | 0.09 | 17.25 | 29.0% |

The pooled r = 0.27 is not a chemistry disagreement. Half of the co-covered set
is exactly the population where the *original* predictor was handed a KEGG
reaction id ModelSEED does not list as an alias, and it is that half that
contributes all the scatter — the mis-mapped panel shows the characteristic
vertical/horizontal stripes of a handful of reused ΔG values sprayed across
hundreds of reactions. On the vouched half the two agree at r = 0.83.
**The dominant effect of the fine-tuning is removing a data defect, not changing
the chemistry model.**

Reported σ rises from a median of **0.36 → 20.12 kcal/mol**. The original's
sub-kcal error bars were never credible; the retrain's are the ones that
correlate with real error.

Figures: `fig1_base_vs_finetuned{,_zoom}.png`,
`fig2_split_by_kegg_mask{,_zoom}.png`, `fig3_sigma.png`, `pair_stats.tsv`.
The `_zoom` renderings clip axes to ±250 kcal/mol (off-scale counts stated on
each panel); statistics are always over the full point set.

---

## 5. Inside the fine-tuned model: what fingerprints does it actually have?

Yes, you can look inside, and it is worth doing.
`results/learned_fingerprints.tsv` is the full read-out: every fragment SMILES,
its energy in kJ/mol, how many training compounds contain it, and how many
ModelSEED compounds contain it.

### 5a. The alphabet grew; the dictionary shrank

| | declared vocabulary | fragments with a non-zero coefficient | rank of the training matrix |
|---|---:|---:|---:|
| dGPredictor, original (KEGG) | 26,317 | 1,514 (247 r1 + 1,248 r2) | **431** |
| dGPredictor, ModelSEED fine-tuned | 40,361 | **1,415** (249 r1 + 1,166 r2) | **382** |
| eQuilibrator, CC 2.x (in this repo) | 241 | 136 | 129 |
| eQuilibrator, component-contribution 0.7 | 163 real + 50 placeholders | 141 | 134 |
| ModelSEED Group Contribution | 271 named groups | — (training set not distributed) | — |

A `BayesianRidge(fit_intercept=False)` coefficient for a feature column that is
identically zero across the training set is *exactly* zero — the column
contributes nothing to `Vᵀ`, nothing to `Xᵀy`, and nothing to the posterior
covariance. Since the training set was not enlarged (it shrank, 4,001 → 3,910),
**96.5% of the expanded vocabulary carries coefficient 0.00**, and the fine-tuning
*reduced* the learned vocabulary from 1,514 to 1,415 and the identifiable
dimensions from 431 to 382.

That is the central structural finding: the retrain widened what the model can
*parse*, not what it *knows*. Every genuinely novel ModelSEED substructure — the
thing the exercise was for — contributes exactly zero energy.

### 5b. The per-fingerprint energies are mostly not identifiable

1,415 non-zero coefficients, rank 382 ⇒ **1,033 directions are fixed by the L2
prior, not by data.** Roughly 73% of the "energy of a fingerprint" is a
minimum-norm artefact. Concretely visible in the read-out:

* `[Mg+2]` gets −226.83 kJ/mol at radius 1 **and** −226.83 at radius 2 — the same
  feature duplicated across the two blocks, so the ridge splits one number in
  half. Neither half means anything on its own.
* 670 of the 1,415 fingerprints occur in **exactly one** training compound.
  `[Mg+2]` (1 training compound → 11 ModelSEED compounds), `O=C=O` (1 → 1),
  `OO` (1 → 55), `O=S([O-])O` (1 → 1).

Do not quote individual fingerprint energies — same gauge trap as the
compound-level offsets in the eQ-vs-dGPMS work.

By contrast, eQuilibrator's basis is nearly fully identified: 136–141
parameters against rank 129–134 on the same ~4,000 observations.

### 5c. What happens outside the trained span

Component-contribution answers this explicitly — it splits every query into
reactant-contribution, group-contribution, and an orthogonal remainder, and
assigns the remainder `RMSE_inf`. The projection matrices for that split are
literally in the parameter file. `BayesianRidge` has no equivalent.

Measured for the fine-tuned model over 28,346 predicted reactions
(`results/rowspace_coverage.json`): the fraction of a reaction's group-change
vector lying outside the training row space has **median 0.42**, and **40% of
reactions are more than half outside**.

Where the model has no information at all it says so *loudly* in most cases —
reported σ tracks unseen-group exposure at ρ = +0.84, and the 40 reactions whose
entire group change lands on zero-coefficient fragments get median σ = 124 kJ/mol.
The one case it does **not** flag is different and larger:

**3,673 reactions (11.5% of all predictions) return ΔG_model = 0.00 ± 3.31 kJ/mol,
where 3.31 = √(1/α) is the model's smallest possible error bar.** These are
reactions whose radius-1 *and* radius-2 group-change vector is identically zero —
isomerisations, racemisations, intramolecular rearrangements. The fingerprint
basis cannot see them, and the model reports its highest confidence on exactly
those. (`results/zero_prediction_audit.json`; e.g. rxn00045, rxn00266, rxn00839,
rxn01004, rxn02243.) This is the radius-2-doesn't-cure-stereo-blindness problem,
quantified.

### 5d. Comparison with (A) Group Contribution and (B) eQuilibrator

Full crosswalk with verbatim group names in `results/vocabulary_crosswalk.md`.
The three alphabets are **not in the same space** — there is no honest
one-to-one alignment, only a structural contrast.

| property | fine-tuned dGPredictor | **(A)** ModelSEED Group Contribution | **(B)** eQuilibrator |
|---|---|---|---|
| unit | RDKit atom environment → canonical fragment SMILES | named chemical group | named chemical group |
| built by | RDKit, automatically, one entry per heavy atom per radius | MFAToolkit rules (Mavrovouniotis / Jankowski) | curated SMARTS-style rules (same lineage) |
| size | 40,361 declared / 1,415 learned | 271 | 163 real + 50 placeholders |
| median entries per compound | ~2 × #heavy atoms | 10 | 6 |
| protonation state in the label | **no** — only implicit in the pH-7 SMILES (456 learned fragments carry an explicit charge) | no | **yes** — every group is `[Hn Zq Mgm]`; 55 charged variants |
| Mg binding | no | one free-ion group `Mg` | yes, 3 Mg-bound variants |
| ring context | implicit (aromatic lowercase atoms in 432 fragments) | explicit (`RW…`, `T…`, `HeteroAromatic`) | explicit (54 `ring` / `fused rings` groups) |
| whole-molecule entries | no | yes (`H2O`, `CO2`, `urea`, `acetate`, …) | yes (50 one-hot placeholders) |
| per-molecule constant | **none** | `Origin` | `Origin [H0 Z0 Mg0]` |
| undecomposable compound | reaction is simply not predicted | `NoGroup` marker (11,488 structures) | placeholder column + `RMSE_inf` |
| out-of-span behaviour | silent, coefficient 0, no dedicated uncertainty term | n/a | explicit null-space projection + `RMSE_inf` |

(A) and (B) are cousins — both descend from the Mavrovouniotis/Jankowski group
set, and both literally contain an `Origin` group. dGPredictor is the outlier:
it replaces curated chemistry with raw substructure enumeration.

Same chemistry, three encodings — phosphate:

* **dGPredictor** (96 learned P-containing fragments): `OP`, `O=P`, `[O-]P`,
  `POP`, `COP`, `O=P(O)(O)O`, `O=P([O-])(O)O`, `O=P([O-])([O-])O`,
  `O=P([O-])([O-])OP(=O)(O)O`, …
* **ModelSEED GC** (12): `WPO3`, `WPO4nW`, `RWOPO2W`, `WCOOPO3`, `prim_phos`,
  `mid_phos`, `pyrophos`, `orthophosphate`, `formylphosphate`, `thioprim_phos`,
  `thiomid_phos`, `Itriphos`
* **eQuilibrator** (36): `-OPO3 [H0 Z-2 Mg0]`, `-OPO3 [H1 Z-1 Mg0]`,
  `-OPO3 [H2 Z0 Mg0]`, `-OPO3-OPO2- [H0 Z0 Mg1]`, `ring -OPO2-OPO2- [H1 Z-1 Mg0]`, …

The one place dGPredictor's alphabet is *richer* than either curated set is
aromatic/heteroaromatic context (432 fragments vs 5 GC groups vs 60 eQ groups) —
consistent with radius-2 resolving fused-ring environments the named vocabularies
lump together.

The protonation row is the most consequential difference. eQuilibrator carries
pKa structure **inside** the basis, so the regression can learn that a
carboxylate and a carboxylic acid are different species. dGPredictor uses one
SMILES per compound (the pH-7 major microspecies) and pushes all protonation
effects into a **post-hoc** Legendre correction (`ddG0_pH_correction`). That is
also why the OPAM2 study in `REPORT.md` could swap pKa sources without retraining:
`dG_model_only` is pKa-blind by construction.

---

## 6. Why the fine-tuned model agrees so much better with Group Contribution

Pooled on ModelSEED dev: **r(GC, original) = 0.22 → r(GC, retrain) = 0.80**
(eQuilibrator as a second reference: 0.26 → 0.86). Decomposition in
`scripts/why_gc_agreement_improved.py` / `results/why_gc_agreement_improved.json`.

**Not coverage.** Restricting both to the identical 18,413 reactions *both*
variants cover changes nothing: 0.219 vs 0.815.

**Not the retrain's 0.00 outputs.** Dropping the 2,080 reactions where the
retrain returns exactly zero leaves r at 0.814 (from 0.815).

**It is the KEGG mis-mapping, almost entirely.** Same 18,413 reactions, split on
`dgpredictor_kegg_mask.json`:

| half | n | GC vs original | GC vs retrain |
|---|---:|---|---|
| KEGG id vouched | 9,139 | r 0.817 · ρ 0.631 · med \|Δ\| 5.62 | r 0.885 · ρ **0.628** · med \|Δ\| **7.41** |
| KEGG id inferred | 9,274 | r **−0.027** · ρ 0.054 · med \|Δ\| 13.98 | r 0.788 · ρ 0.558 · med \|Δ\| 4.24 |

The mechanism is visible in the raw values: across the 9,274 mis-mapped
reactions the original carries only **956 distinct ΔG values** (10.3% — one
value reused ~10× on average, because a handful of inferred KEGG ids like R09245
and R06601 are pasted onto hundreds of reactions each), against 3,561 distinct
values from the retrain. A near-constant series cannot correlate with anything;
r ≈ 0 there is arithmetic, not chemistry. The retrain never goes through a KEGG
id at all — it decomposes the ModelSEED reaction's own stoichiometry — so the
defect simply does not exist for it.

**Second-order: the retrain and GC now read the same molecules.** GC energies
come from MFAToolkit run over ModelSEED structures; the retrain decomposes those
same structures, at the same pH-7 charge state, with the same stoichiometry. The
original decomposes *KEGG* structures for a *KEGG-written* reaction, which can
differ in protonation, hydration and stoichiometric convention even when the
mapping is correct. This is real but small — it is the residual after removing
mis-mapping, i.e. the 0.817 → 0.885 in the vouched row.

**And that residual is leverage, not agreement.** Three warnings in the vouched
half:

* Spearman ρ is **unchanged** (0.631 → 0.628). In rank terms the retrain is not
  better at all where the original was working correctly.
* Median \|Δ\| gets **worse**, 5.62 → 7.41 kcal/mol.
* Trimming the tail flips the ranking: at \|ΔG_GC\| < 100 it is 0.809 (original)
  vs 0.797 (retrain); at < 50, 0.448 vs 0.427.

Pooled r on this pair is tail-dominated because GC and the retrain are *both*
strictly additive group models over the *same* structures, so both scale with
molecule size and agree almost trivially on large reactions. Trimming the whole
like-for-like set to \|ΔG_GC\| < 50 collapses the retrain's r from 0.815 to
0.357.

**Summary.** The improvement is one real fix plus one statistical artefact: the
retrain removed a data defect that had destroyed half the original's values, and
what looks like extra agreement on the surviving half is tail leverage from two
additive models sharing an input structure file. Where the original was
answering the right reaction, the retrain is *not* closer to Group Contribution —
by typical error it is about 1.8 kcal/mol further away.

Caveat: on this snapshot GC is Convention A while both dGPredictor variants are
Convention B. Consistent A vs consistent B cancels for mass-balanced reactions
(tested previously: GC−eQ residual vs net H⁺ slope −2.67, nowhere near the
±9.539 a systematic gap would give), so the panels are drawn on raw stored
values — but it is not a zero-risk comparison.

---

## 7. Root cause of the mis-mapping, and how the retrain avoids it

Traced in `scripts/trace_kegg_mapping.py` / `results/kegg_mapping_trace.json`.

### How the original path maps a ModelSEED reaction to KEGG

Predictions are staged in
`Biochemistry/Thermodynamics/dGPredictor/json_files/reaction_*_dG.json`, keyed
`rxn##### → KEGG R-id → {dG_mean, dG_uncer}`.
`Update_Reaction_dGPredictor_Energies.py` on dev does nothing clever — it just
averages `dG_mean` over whatever KEGG ids it finds under each reaction. **The
mapping is fixed upstream, in the freiburgermsu repo**, by
`dG_prediction_modelseed_dev_branch.ipynb`: read each ModelSEED reaction's
`aliases`, pull the KEGG R-id out, look that id up in
`data/KEGG_rxn_eqn_master_branch.json`, and predict from *the KEGG equation*.

### The bug

Of 28,502 staged reactions:

| | n |
|---|---:|
| staged KEGG id **is** a ModelSEED alias of that reaction | 10,718 |
| staged KEGG id **conflicts** with a listed alias | **0** |
| reaction has **no KEGG alias at all**, yet got an id anyway | **17,784** |

Those 17,784 are not explained by any alias file on disk (dev's
`Unique_ModelSEED_Reaction_Aliases.txt`: 0; the repo's shipped copy: 0). They are
explained, at **100.0%, exact set match**, by a carry-forward: each alias-less
reaction was handed the KEGG id of the nearest *preceding* reaction in file order
that had one. The loop that did it:

```python
for i, rxn in tqdm(enumerate(json_read)):
    try:
        rxn_alias = rxn['aliases']
        for ki in rxn_alias:
            if 'KEGG' in ki:
                kegg_id_str = ki          # never initialised, never reset
        KEGG_id = kegg_id_str.replace(' ', '').split(':')[1]
        ...
    except:
        KEGG_id_ls.append('No KEGG id')
```

`kegg_id_str` lives outside the loop's control. A reaction whose alias list has
no KEGG entry falls straight through the inner `for`, `kegg_id_str` still holds
the *previous* reaction's value, no exception is raised, and the
`except: 'No KEGG id'` guard never fires. Since ModelSEED reaction ids are
sorted and provenance-clustered, alias-less reactions come in long contiguous
runs — hence R09245 on 858 reactions, R06601 on 733, R10126 on 544.

The `.py` port in the repo, `dG_prediction_modelseed_dev_branch_file_run.py:140`,
now sets `kegg_id_str = None` each iteration. **The fix exists; the staged data
on dev predates it**, and 17,271 of the affected reactions carry a stored
`dGPredictor` record there.

What it looks like in the database — six consecutive, chemically unrelated
reactions that all inherited hexaprenyl-diphosphate synthase's ΔG:

| rxn | name | stored `dGPredictor` |
|---|---|---|
| rxn13478 | Heptadecanoate transport via proton symport | −67.33 ± 2.29 `>` |
| rxn13479 | Isomerase for keto-meroacid-2 | −67.33 ± 2.29 `>` |
| rxn13480 | Isomerase for methoxy-meroacid-2 | −67.33 ± 2.29 `>` |
| rxn13481 | keto-ylation-1 for keto-meroacid-1 | −67.33 ± 2.29 `>` |
| rxn13482 | keto-ylation-1 for keto-meroacid-2 | −67.33 ± 2.29 `>` |
| rxn13483 | keto-ylation-2 for keto-meroacid-1 | −67.33 ± 2.29 `>` |

A transport reaction and two isomerases cannot share a ΔG. Note the error bar:
±2.29 kcal/mol, i.e. confidently wrong.

### How the retrain avoids it

Not by fixing the mapping — **by deleting the mapping step.** The retrain never
constructs a KEGG reaction:

1. `parse_modelseed_stoichiometry()` reads ModelSEED's own `stoichiometry`
   column into `{cpd#####: coefficient}`;
2. compounds are decomposed as `cpd#####` directly from ModelSEED structures;
3. output is staged as `modelseed_retrained_dG.json`, keyed
   `rxn##### → {dG_mean, dG_uncer}` — **one value per reaction, no KEGG layer**;
4. `Update_Reaction_dGPredictor_ModelSEED_Energies.py` reads that straight
   through, with no averaging over ids.

There is no variable that can carry the wrong reaction's identity forward,
because reaction identity is never re-derived — it is the dict key throughout.

KEGG ids do still appear once in the retrain, at a different level: the 673
*compounds* of the training set are translated `C##### → cpd#####` via
`Unique_ModelSEED_Compound_Aliases.txt`. That is a one-time, compound-level
lookup used only to build the training matrix, and it **fails closed** — 67
training reactions with an unmappable compound were dropped rather than guessed.

Its behaviour on the 17,784 affected reactions is worth being clear about: it
supplies its own value for **10,652** and **declines the other 7,132** (their
compounds have no complete structure). The block above is in the declining set —
all six now carry no `dGPredictor-ModelSEED` record at all. So the retrain does
not repair those reactions; it stops asserting a number for them, which is the
correct outcome and is why the coverage gain (+11,352) and the accuracy gain are
not the same reactions.

---

## Reproduce

```bash
cd /scratch/ctaylor/dgpredictor_finetune_review
python3 scripts/rebuild_finetuned_features.py       # rebuild + refit, ~7 s
python3 scripts/validate_refit.py                   # 6e-7 kJ/mol vs shipped
python3 scripts/identifiability_and_confidence.py
python3 scripts/rowspace_coverage.py
python3 scripts/zero_prediction_audit.py
python3 scripts/compare_group_vocabularies.py
python3 scripts/crosswalk_examples.py
python3 scripts/why_gc_agreement_improved.py            # needs results/_devsnap2_thermo.pkl
python3 scripts/trace_kegg_mapping.py                  # root-cause of the mis-mapping
python3 scripts/plot_vocabulary_comparison.py       # needs matplotlib

cd /scratch/ctaylor/core_models_analysis
python scripts/plot_dgpredictor_base_vs_finetuned.py
```

`compare_group_vocabularies.py` shells out to the `eq3` env (the only one with
`component_contribution`) to read `cc_params.npz`; everything else runs on the
base interpreter. The scatter script needs the `core_models_analysis` env
(matplotlib) and reads `/scratch/ctaylor/tmp/devsnap2`.
