# 1. The data

## 1.1 Snapshot

All reaction data is read from a read-only `git archive` of ModelSEED `dev` at
commit **49563c6f**, extracted to `/scratch/ctaylor/tmp/devsnap2`. Only
`Biochemistry/` and `Scripts/` are present in that archive.

The **cascade code** (`reversibility_heuristics.py`) is imported from the local
working checkout `/scratch/ctaylor/ModelSEEDDatabase` instead, because the
archive ships no `Libs/Python`. The two copies' `DEFAULT_HEURISTICS` lists were
compared before running and are identical in order and content:

```
atp_synthase → abc_transporter → stored_bounds → mmdeltag_band → low_energy → default
```

so the choice of checkout does not affect any direction call.

**Why not the live working tree.** It carries neither of the two things this
analysis needs:

| | live checkout (`claude-changes`) | devsnap2 (`dev` @ 49563c6f) |
|---|---|---|
| `Group contribution` | 25,826 non-sentinel, median σ 5.06 | **27,313**, median σ **10.28** (Convention A rebuild, dev `ad34d6ab`) |
| `eQuilibrator` | 19,498 | **25,028** |
| `dGPredictor` (legacy, KEGG-keyed) | 27,715 | 27,715 — **present but never read** |
| `dGPredictor-ModelSEED` (retrain) | absent | **31,924** |

The Convention A fingerprint (coverage 25,812 → 27,313, σ roughly doubling) is
how the rebuild is identified.

## 1.2 The four sources

Reaction thermodynamics live in each reaction record as
`thermodynamics[label] = [ΔG′°, σ, operator]`, in kcal/mol.

| source | short | how ΔG′° is produced | how σ is produced |
|---|---|---|---|
| **Group Contribution** | GC | sum of fitted group energies | propagated uncertainty of the fitted group energies; post-rebuild the resolver reports √(mean σᵢ² + var ΔGᵢ) |
| **eQuilibrator** | EQ | component contribution: a reactant-contribution layer anchored on measured reactions, group contribution only for the orthogonal complement. Invoked at pH 7.0, I = 0.25 M, 298.15 K, scoring ModelSEED's own stoichiometry | σ² = **ν**ᵀ**Σν** from the component-contribution covariance |
| **dGPredictor-ModelSEED** | DG | BayesianRidge on radius-1 and radius-2 atom-centred fragment count changes, retrained on ModelSEED structures; keyed directly by `rxnNNNNN` | posterior predictive standard deviation √Var[ΔG \| **x**, 𝒟] |
| **TECRDB** | TEC | NIST experimental ΔG′° = −RT ln K′ | experimental standard deviation over the contributing measurements |

The three predictors are the *predictors*; TECRDB is the *measurement* and is
treated differently everywhere.

### Coverage over the 56,002 non-EMPTY reactions

| source | reactions with a usable ΔG′° | after vetoes (§1.4) |
|---|---:|---:|
| Group Contribution | 27,313 | 27,313 |
| eQuilibrator | 25,028 | 20,059 |
| dGPredictor-ModelSEED | 31,924 | 31,413 |
| TECRDB | 1,550 | 1,550 |
| **union** | **33,337** | **33,289** |

22,665 non-EMPTY reactions have no thermodynamic source at all. That is the
coverage ceiling and no method in this folder changes it.

Sentinel values (`ΔG′° = 10000000`) are excluded throughout; in devsnap2 all
28,689 GC sentinel rows also carry operator `?`, which is how the loaders reject
them. *(That coincidence is load-bearing in
`optimize_thermo_source_assignment.load_db`, which gates on the operator rather
than on the value — worth tightening if that loader is reused on a snapshot
where the two do not agree.)*

## 1.3 The experimental reference set

TECRDB is fetched from the eQuilibrator Zenodo deposit
(doi:10.5281/zenodo.3978440): 4,544 rows of K′ with temperature and pH, keyed by
KEGG compound ids.

Matching to ModelSEED reactions is by **structure**, not by identifier. Each
compound is resolved to an RDKit structure key and a reaction becomes the
(reactant multiset, product multiset) pair with protons dropped, matched in both
directions. Two tiers:

| tier | key | reactions | meaning |
|---|---|---:|---|
| `stereo_exact` | full InChIKey of the neutralised parent | **802** | distinguishes anomers and D/L pairs |
| `skeleton` | InChIKey connectivity block (first 14 chars) | 748 | connectivity only — can conflate stereoisomers |

Only the 802 `stereo_exact` matches are used for fitting and for validation.
The skeleton tier appears in the graded output but is capped at SILVER.

Experimental scatter on the matched set: median σ 0.15 kcal/mol, p90 1.04, max
3.15. 551 of the 1,550 rest on a single measurement; both `n_measurements` and
the sd are carried through to the output so a consumer can impose its own floor.

Measured error of each predictor against this set:

| source | n | median \|error\| kcal/mol | fraction within 2 kcal/mol |
|---|---:|---:|---:|
| eQuilibrator | 794 | 0.45 | 85.6% |
| dGPredictor-ModelSEED | 802 | 0.47 | 84.8% |
| Group Contribution | 802 | 1.57 | 56.6% |

No predictor dominates on magnitude, which is why the arbitration question is
worth asking at all.

## 1.4 Known defects, and the vetoes that encode them

Four defects are established elsewhere in this analysis. Each removes a source
from a reaction outright rather than being modelled as extra uncertainty,
because in each case the value is not a noisy estimate — it is not an estimate.

| veto | reactions | why |
|---|---:|---|
| **eQuilibrator sentinel**, σ > 100 | 4,934 | eQuilibrator flags compounds it cannot estimate by inflating variance by 10⁶. Stored uncertainties are strictly bimodal — real ones cap at **65.3** kcal/mol, sentinels start at **7,504.6**, nothing in between — so the cut at 100 sits in a two-orders-of-magnitude empty gap rather than being tuned. This is the source explicitly disclaiming the reaction. |
| **eQuilibrator MetaNetX collision** | 35 | `Retrieve_eQuilibrator_Reactions_Energies.py` writes `lhs[mnx_id] = |coeff|` instead of accumulating, so two ModelSEED compounds sharing one MetaNetX id silently overwrite each other. List read from `results/eq_vs_dgpms/reconciliation.tsv`. |
| **dGPredictor-ModelSEED on quinone/quinol** | 511 | 52.8% sign disagreement with eQuilibrator on that couple, median σ 80.3 — the retrain regressed on two-electron aromatic redox. Self-flagged, but not reliably enough to leave to the σ model. |
| **legacy `dGPredictor` KEGG mis-mapping** | (n/a) | 17,271 reactions carry a value predicted from a KEGG reaction ModelSEED does not list for them. **Not applicable here**: the legacy label is never read. `dGPredictor-ModelSEED` is keyed by ModelSEED id and is structurally immune, so no mask is applied. |

## 1.5 Core-model inputs

5,683 Kegg2 core models, `data/core_models_kegg2/*.json` (symlinked to
`/scratch/ctaylor/core_models_kegg2`). Media is the standard KBase complete
media, `ModelSEEDDatabase/Media/KBaseMedia.cpd`.

A model's **unique reaction** set is its distinct `annotation['seed.reaction']`
values, normalised by `seed_annotation.seed_id` — which strips the stray
compartment-letter suffix on 17 transport annotations (`rxn11322_c` → `rxn11322`).
Exchange, sink, demand and biomass pseudo-reactions carry no SEED annotation and
are therefore excluded. A **unique compound** is the metabolite id with its
trailing `_<compartment>` removed, so a compound present in two compartments
counts once.

| | min | median | mean | max |
|---|---:|---:|---:|---:|
| unique reactions per model | 20 | 128 | 123.1 | 187 |
| unique compounds per model | 41 | 124 | 119.2 | 163 |

Combined across all 5,683 models: **239 distinct reactions**, **182 distinct
compounds**. The models are homogeneous — central metabolism — which is why the
union is so much smaller than the per-model average times the model count.

**The core set is far better measured than the database at large.** 69 of the
239 core reactions (28.9%) carry a TECRDB match, against 1,550 of 56,002 (2.8%)
overall — a 10× enrichment, because central metabolism is what NIST measured.
