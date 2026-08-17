# dGPredictor and its ModelSEED retrain — summary

Reviewed: `github.com/freiburgermsu/dGPredictor` @ `5d656c5`, against ModelSEED
`dev` @ `49563c6f` (the live dev branch as of 2026-08-14). All numbers below were
recomputed from the data, not taken from documentation. Working files and scripts:
`/scratch/ctaylor/dgpredictor_finetune_review/`.

---

## 1. What dGPredictor is

It is a **group-contribution method** — the same family as the Mavrovouniotis /
Jankowski method ModelSEED already uses, and as the component-contribution
method behind eQuilibrator. The premise is the familiar one: the free energy of a
reaction is approximately the sum of contributions from the chemical groups that
are made and broken.

The one real difference is **how the groups are defined**:

| | how groups are chosen |
|---|---|
| ModelSEED Group Contribution | ~271 named groups, defined by chemists (carboxylate, primary alcohol, …) |
| eQuilibrator | 163 named groups, each labelled with its protonation and Mg state |
| **dGPredictor** | groups are generated automatically: take every atom, walk out 1 bond (and separately 2 bonds), and record whatever substructure you find |

So dGPredictor's "groups" are not curated chemistry — they are every local atomic
environment that happens to occur in the database, written as a SMILES fragment,
e.g. `CO`, `C=O`, `O=P`, `[O-]P`, `C[O-]`. The ModelSEED version of the model has
**40,361** such fragments in its dictionary.

---

## 2. How it is trained

**Step 0 — the measurements (4,001 reactions).**
These are experimental: equilibrium constants from the NIST TECRDB, Alberty's
tabulated formation energies, and a handful of redox potentials. They are stored
as *reaction stoichiometries plus a measured ΔG* — compound IDs and coefficients
only, **no molecular structures**.

Worth knowing: this is the **same experimental corpus eQuilibrator is trained
on**. dGPredictor uses no eQuilibrator *predictions*, but it is not an
independent dataset. Benchmarking either method against TECRDB is therefore
in-sample for both.

**Step 1 — chop up every molecule.**
Software (RDKit) decomposes each compound in the database into its atom-centred
fragments. This is pure structure handling; no energies are involved. It
produces the fragment dictionary.

**Step 2 — turn each measured reaction into a "what changed" vector.**
For each of the 4,001 reactions, count the fragments on the product side minus
the fragments on the substrate side. Anything present on both sides cancels.

**Step 3 — solve for a per-fragment energy.**
Find the set of fragment energies whose sums best reproduce all 4,001 measured
values at once. (The estimator is a regularised least-squares fit; the
regularisation just means "prefer the smallest energies that fit", which keeps
the answer stable when the data under-determines it.)

**Step 4 — predict a new reaction.**
Look up the fragments of each compound, take products minus substrates, multiply
by the fitted energies, add up. Then apply the standard pH / ionic-strength
(Legendre) correction using pKa values to get ΔG′° at pH 7. That last step is
ordinary thermodynamics, not machine learning.

```
4,001 measured reactions          every compound in the database
   (IDs + stoichiometry)                      |
            |                          chop into fragments
            |                                 |
            +------------- combine ----------+
                            |
              solve: fragment energies that reproduce
                     the 4,001 measurements
                            |
                    ~1,400 fragment energies
                            |
            new reaction ---+---> sum fragments ---> + pH correction ---> ΔG'°
```

**The critical limitation.** A fragment only gets an energy if it appears in one
of the ~600 compounds involved in those 4,001 reactions. Of the 40,361 fragments
in the dictionary, **1,415 have an energy and 38,946 are exactly zero**. The
model does not warn you when a reaction contains zero-energy fragments — it
returns a number regardless. For a typical ModelSEED reaction, ~42% of its
fragment content lies outside anything the training data covered.

A second blind spot follows from step 2: if a reaction's fragments cancel
completely — isomerases, racemases, many rearrangements — the model sees nothing
and returns **0.00 kJ/mol with its smallest error bar**. This affects 3,673
reactions (11.5% of its predictions).

---

## 3. What the ModelSEED retrain changed

| | original | ModelSEED retrain |
|---|---|---|
| Experimental measurements | 4,001 | **the same 4,001, untouched** |
| Usable after ID translation | 4,001 | 3,910 |
| Molecules drawn from | KEGG structures | ModelSEED structures |
| Compound IDs | `C#####` | `cpd#####` |
| pKa source | ChemAxon, run live | ModelSEED's stored values |
| Fragment dictionary | 26,317 | 40,361 |
| Fragments with an energy | 1,514 | **1,415** |
| Reactions covered | 27,715 | 31,924 |

**This is not "fine-tuning" in the usual sense** — there is no continuation from
the previous model and no new data. It is the same regression, refit from
scratch, on the same measurements, with the molecules re-drawn from ModelSEED's
structure files instead of KEGG's.

Two consequences worth stating plainly:

- **Coverage improved, knowledge did not.** The dictionary grew by 53%, but the
  number of fragments that actually received an energy went *down* (1,514 →
  1,415), because the training set shrank slightly during ID translation. Every
  genuinely novel ModelSEED substructure is assigned an energy of zero.
- **The two models have different alphabets.** Only ~72% of their fragments
  overlap. The cause is that ~1 training compound in 7 is drawn differently in
  the two databases — different tautomer, different proton placement, or
  different stereochemistry (e.g. urate is the keto form in KEGG and the enol
  anion in ModelSEED). Thirteen compounds are outright different molecules under
  the same ID.

---

## 4. A data defect in ModelSEED dev that should be fixed

This is separate from the model itself and is the most actionable finding.

**What happened.** dGPredictor scores a *KEGG* reaction, so each ModelSEED
reaction had to be matched to a KEGG reaction first. The script that did this
(in the dGPredictor repository, not in ModelSEED) had a bug: when a ModelSEED
reaction had **no** KEGG counterpart, it silently reused the KEGG ID of the
previous reaction in the file instead of skipping it.

**Scale.** Of 28,502 reactions given a prediction, **17,784 have no KEGG ID at
all** and were handed a neighbour's. Reconstructing "take the previous
reaction's ID" reproduces the assignment for **100.0%** of them. Because
ModelSEED IDs are grouped by provenance, single IDs were reused hundreds of
times — R09245 on 858 reactions, R06601 on 733.

**What it looks like:**

| reaction | name | stored dGPredictor ΔG |
|---|---|---|
| rxn13478 | heptadecanoate transport | −67.33 ± 2.29 |
| rxn13479 | isomerase, keto-meroacid-2 | −67.33 ± 2.29 |
| rxn13480 | isomerase, methoxy-meroacid-2 | −67.33 ± 2.29 |
| rxn13481–83 | three keto-ylation steps | −67.33 ± 2.29 |

Six consecutive, chemically unrelated reactions carrying one number — farnesyl-
diphosphate synthase's. Similarly, rxn00019 (a nitroalkane oxidase) carries
RuBisCO's energy, −8.63 ± 0.04 kcal/mol.

**Where it lives in ModelSEED dev.**

| location | contents |
|---|---|
| `Biochemistry/Thermodynamics/dGPredictor/json_files/reaction_*_dG.json` | the wrong KEGG IDs (17,784 reactions) |
| `Biochemistry/reaction_*.json` → `thermodynamics.dGPredictor` | the resulting energies (**17,271** reactions, 62% of all dGPredictor records) |
| `Biochemistry/reaction_*.json` → `deltag` / `deltagerr` / `reversibility` | **7,466** reactions where that value became the database's served answer |

**Not affected:** ModelSEED's alias files are clean — they correctly list no
KEGG ID for these reactions. The bad IDs exist only inside the dGPredictor
staging files. The retrained `dGPredictor-ModelSEED` records are also unaffected,
because that pipeline never uses KEGG IDs; it reads ModelSEED's own
stoichiometry directly.

**Two aggravating factors.**

1. `Scripts/Thermodynamics/Update_Reaction_dGPredictor_Energies.py` takes the
   staged KEGG IDs at face value — it never checks them against the reaction's
   own aliases — and it is a live step in `Rerun_Thermodynamics.sh`, so the
   values are re-imported on every thermodynamics refresh.
2. `Promote_Reaction_Thermodynamics_to_Canonical.py` chooses, within a tier, the
   estimate with the **smallest stated uncertainty**. The mis-mapped values are
   the over-confident ones (±0.03–0.17), so the rule systematically prefers
   them. Its own worked example is rxn00019: it protects "a tight prediction
   (−8.6 ± 0.04)" over "a wildly-uncertain outlier (−100 ± 71)" — but the tight
   number is RuBisCO's, and the rejected one is the only estimate about the
   right reaction.

**Suggested fixes, in order of effort:**

- *Immediate:* in `Update_Reaction_dGPredictor_Energies.py`, skip any staged
  KEGG ID that is not listed in that reaction's aliases. One condition; catches
  all 17,784.
- *Proper:* regenerate the staged predictions. The corrected script already
  exists in the dGPredictor repo
  (`dG_prediction_modelseed_dev_branch_file_run.py`); only the older notebook
  had the bug.
- *Then:* re-run promotion so the 7,466 canonical values are recomputed.

---

## 5. How the three thermodynamic sources compare

| | ModelSEED Group Contribution | eQuilibrator | dGPredictor (retrain) |
|---|---|---|---|
| groups | 271, chemist-defined | 163, chemist-defined | 40,361 auto-generated |
| groups constrained by data | — | 134 of 141 used | **382 of 1,415 used** |
| protonation states | not in the group label | explicit in every label | not in the label |
| Mg binding | free ion only | explicit | absent |
| per-compound energies | yes (45,708) | yes (30,607) | not published (but computable) |
| behaviour outside training | — | flags it with a large uncertainty | silently returns zero |

The practical reading: eQuilibrator's group set is small and almost fully
determined by the data. dGPredictor's is very large and mostly *not* — about
three-quarters of its fragment energies are set by the fitting procedure's
smoothing assumption rather than by measurement. That is not fatal, because the
errors partly cancel across a balanced reaction, but it means individual
fragment energies should not be quoted, and predictions on unusual chemistry
should be treated as extrapolation.

---

## 6. Bottom line

1. The retrain's real achievement is **coverage and correct reaction identity** —
   31,924 reactions scored from ModelSEED's own stoichiometry, with no KEGG
   matching step and therefore no possibility of the mis-mapping defect.
2. It did **not** improve the underlying chemistry. Same measurements, fewer
   learned fragments, slightly worse typical error on reactions where the
   original was working correctly.
3. The KEGG mis-mapping in ModelSEED dev is a genuine data defect affecting
   17,271 stored values and 7,466 canonical ΔG / direction assignments, and it
   is re-applied on every pipeline run. This is worth fixing regardless of which
   dGPredictor variant is preferred going forward.
