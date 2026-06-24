# eQuilibrator Metabolite Identity & Static Energy Provenance

*How the reaction-direction columns are computed, how eQuilibrator resolves metabolite identity, and how the ModelSEEDDatabase static eQuilibrator energies are generated.*

Scope: explains the `Flamholz_2012` and `Beber_2022` columns of
`results/reaction_directions_literature_vs_llm.tsv`, the eQuilibrator
compound-identity mechanism behind them, and the ModelSEEDDatabase (MSDB)
scripts that produced the static `thermodynamics["eQuilibrator"]` energies.
Date: 2026-06-24.

---

## 1. The four direction sources

The columns in
[`results/reaction_directions_literature_vs_llm.tsv`](/scratch/ctaylor/core_models_analysis/results/reaction_directions_literature_vs_llm.tsv)
map to four distinct engines, assembled by two scripts
([`estimate_directions_literature.py`](/scratch/ctaylor/core_models_analysis/scripts/estimate_directions_literature.py)
for columns 1, 2, 4 and
[`estimate_directions_eq3.py`](/scratch/ctaylor/core_models_analysis/scripts/estimate_directions_eq3.py)
for column 3, merged in afterward):

| Column | Engine | Energy source | Direction rule |
|---|---|---|---|
| `Jankowski_2008` | Group contribution (MFAToolkit) | MSDB `thermodynamics["Group contribution"]` | Henry‑2007 ΔG′ feasibility window (1e‑5–0.02 M) |
| `Flamholz_2012` | eQuilibrator (bundled in MSDB) | `MetaNetX_Reaction_Energies.tbl` ln_RI column | reversibility index, \|ln γ\| < ln(1000) |
| `Beber_2022` | eQuilibrator **3.0**, run live | `equilibrator-api 0.7.0` component contribution | reversibility index, \|ln γ\| < ln(1000) |
| `LLM_Opus_4.8` | LLM | `AICurationCacheReactionDirectionality.json` | model's own `directionality` field |

**Key fact confirmed by the code:** `Flamholz_2012` and `Beber_2022` use the
*identical* direction heuristic — eQuilibrator's reversibility index with the
same `LNRI_THRESH = 6.9077553` (= ln 1000) cutoff. The docstring at
[`estimate_directions_literature.py:17-22`](/scratch/ctaylor/core_models_analysis/scripts/estimate_directions_literature.py#L17-L22)
is explicit that this is *deliberately* the same rule, "NOT the Henry‑2007
rule." So the two eQuilibrator columns differ only in the **energies** — and,
as it turns out, in **how each one resolves a metabolite's identity**.

---

## 2. How eQuilibrator 3.0 "knew what the metabolite was" (the `Beber_2022` column)

Produced live by
[`estimate_directions_eq3.py`](/scratch/ctaylor/core_models_analysis/scripts/estimate_directions_eq3.py).
Identity resolution is a two-step handshake: *we* translate each ModelSEED
compound into a cross-reference accession string, and *eQuilibrator* matches
that string against its local compound cache.

### Step 1 — build a priority-ordered accession list from the ModelSEED compound's aliases

[`estimate_directions_eq3.py:46-66`](/scratch/ctaylor/core_models_analysis/scripts/estimate_directions_eq3.py#L46-L66):

```python
def compound_accessions(c):
    accs = []
    m = _KEGG.search(text);  accs.append("kegg:" + m.group(1))             # e.g. kegg:C00003
    m = _CHEBI.search(text); accs.append("chebi:CHEBI:" + m.group(1))
    m = _BIGG.search(text);  accs.append("bigg.metabolite:" + m.group(1))
    m = _MNX.search(text);   accs.append("metanetx.chemical:" + m.group(1))  # metanetx.chemical:MNXM…
    ik = c.get("inchikey");  accs.append(ik)                                # bare InChIKey fallback
    return accs
```

For each `cpdXXXXX` we scrape its `aliases` block for an identifier in
confidence order: **KEGG → ChEBI → BiGG → MetaNetX → InChIKey**.

### Step 2 — hand each accession to eQuilibrator's `get_compound`, first hit wins

[`estimate_directions_eq3.py:92-104`](/scratch/ctaylor/core_models_analysis/scripts/estimate_directions_eq3.py#L92-L104):

```python
def resolve(cpd):
    obj = None
    for acc in cpd_accs.get(cpd, []):
        try:    obj = cc.get_compound(acc)   # cc = ComponentContribution()
        except Exception: obj = None
        if obj is not None: break            # accept the first namespace that resolves
    resolved[cpd] = obj
    return obj
```

### What `cc.get_compound(acc)` does internally

(eQuilibrator's documented mechanism, not our code.) eQuilibrator ships a local
**compound cache** (`equilibrator-cache`, a SQLite database downloaded from
Zenodo — the same reconciliation MetaNetX produces). Every row is one chemical
species carrying:

1. cross-references to *all* registries — `kegg:`, `chebi:`,
   `bigg.metabolite:`, `metanetx.chemical:`, etc.;
2. a structure (InChI / InChIKey);
3. a **precomputed group-decomposition vector** for the component-contribution
   model.

`get_compound("kegg:C00003")` looks up which cache row registers that accession
and returns its `Compound` object. That is literally how it "knows what the
metabolite is" — string-match the accession to a pre-reconciled cache row that
already has the structure attached. A bare InChIKey is matched the same way
(against the structure column) as the last-resort fallback. eQuilibrator never
sees a ModelSEED ID.

### From identity to direction

[`estimate_directions_eq3.py:127-139`](/scratch/ctaylor/core_models_analysis/scripts/estimate_directions_eq3.py#L127-L139):
the resolved `Compound` objects are assembled into an eQ `Reaction(stoich_dict)`,
checked with `r.is_balanced()`, then `cc.ln_reversibility_index(r)` returns
ln γ = (2/N)·(ΔG′m/RT) (Noor 2012). The ΔG′m underneath is eQuilibrator 3.0's
**component-contribution** estimate at pH 7, I = 0.25 M, no Mg, 298.15 K — for
each compound it pulls the group vector from the cache row and computes ΔG′f as
(reactant contribution where training data exists) + (group contribution for
the residual). Direction: `>` if `lnri < −ln1000`, `<` if `> +ln1000`, else `=`.

The **only** thing that makes this "Beber 2022" rather than "Flamholz 2012" is
the engine: `from equilibrator_api import ComponentContribution` resolving
against the **eQ 3.0 cache** (the `eq3` conda env, `equilibrator-api 0.7.0`,
`XDG_CACHE_HOME=/scratch/ctaylor/eq_cache` — see the run line at
[`estimate_directions_eq3.py:4-6`](/scratch/ctaylor/core_models_analysis/scripts/estimate_directions_eq3.py#L4-L6)).
Newer training data, updated pKa / component model, more compounds reconciled —
same reversibility-index heuristic.

---

## 3. The ModelSEEDDatabase scripts that computed the static eQuilibrator energies

Two **retrieve** scripts call eQuilibrator and write `.tbl` files; two
**update** scripts stamp those values into the JSONs. All in
`ModelSEEDDatabase/Scripts/Thermodynamics/`.

### 3a. Compound formation energies — `Retrieve_eQuilibrator_Compound_Energies.py`

[`Retrieve_eQuilibrator_Compound_Energies.py:12-25`](/scratch/ctaylor/ModelSEEDDatabase/Scripts/Thermodynamics/Retrieve_eQuilibrator_Compound_Energies.py#L12-L25)
computes a ΔG′f for **every MetaNetX compound** by treating "formation" as a
one-sided reaction:

```python
(mnx, inchikey) = line.split('\t')                                 # from Structures_in_ModelSEED_and_eQuilibrator.txt
equilibrator_reaction = Reaction.parse_formula(ccache.get_compound, ' = ' + mnx)   # " = MNXM…"
result = equilibrator_calculator.standard_dg_prime(equilibrator_reaction)
dG0_prime   = result.value.to('kilocal / mole').magnitude
uncertainty = result.error.to('kilocal / mole').magnitude
```

- **Identity is the MetaNetX accession** (`MNXMxxxxx`), fed to
  `ccache.get_compound` exactly like §2.
- The `' = ' + mnx` trick makes a degenerate "formation" reaction so
  `standard_dg_prime` returns the compound's transformed formation energy at
  pH 7 / I = 0.25 M / 298.15 K.
- Output (kcal/mol) → `Biochemistry/Thermodynamics/eQuilibrator/MetaNetX_Compound_Energies.tbl`.

### 3b. Reaction energies — `Retrieve_eQuilibrator_Reactions_Energies.py`

This produced the static reaction ΔG′° / ln_RI values. **The metabolite-identity
mechanism here is different from §2 — it is structural (InChIKey), not
alias-accession.** The map is built from
`Biochemistry/Structures/MetaNetX/Structures_in_ModelSEED_and_eQuilibrator.txt`
(a precompiled MetaNetX↔InChIKey table of structures present in *both*
databases), indexed at **three levels of strictness**
([`Retrieve_eQuilibrator_Reactions_Energies.py:33-42`](/scratch/ctaylor/ModelSEEDDatabase/Scripts/Thermodynamics/Retrieve_eQuilibrator_Reactions_Energies.py#L33-L42)):

```python
mnx_inchikey_dict[inchikey] = mnx                       # 1. full InChIKey
inchikey = "-".join(inchikey.split('-')[0:2])           # 2. proton-neutral (skeleton + stereo block)
mnx_inchikey_dict[inchikey] = mnx
inchikey = inchikey.split('-')[0]                       # 3. stereo-neutral (14-char skeleton only)
mnx_inchikey_dict[inchikey] = mnx
```

Each ModelSEED compound's InChIKey is then matched at full → proton‑neutral →
stereo‑neutral strictness
([`Retrieve_eQuilibrator_Reactions_Energies.py:60-87`](/scratch/ctaylor/ModelSEEDDatabase/Scripts/Thermodynamics/Retrieve_eQuilibrator_Reactions_Energies.py#L60-L87)),
building `seed_mnx_structural_map[cpd] = MNXM…`. The in-code comment records the
rationale and counts:

> "As per email from Elad and Moritz, we should not expect estimated energies to
> deviate between pseudoisomers (protons) and stereoisomers… 17,071 matches based
> on full inchikey / 17,863 using proton-neutral / 18,559 using stereo-neutral."

(Elad Noor and Moritz Beber — the eQuilibrator authors — sanctioning the
skeleton-level fallback.) The reaction is re-expressed entirely in MetaNetX IDs
and eQuilibrator computes it at the reaction level
([`Retrieve_eQuilibrator_Reactions_Energies.py:161-177`](/scratch/ctaylor/ModelSEEDDatabase/Scripts/Thermodynamics/Retrieve_eQuilibrator_Reactions_Energies.py#L161-L177)):

```python
equation_str = ' + '.join(f'{v} {k}' for k,v in lhs.items()) + " = " + \
               ' + '.join(f'{v} {k}' for k,v in rhs.items())
equilibrator_reaction = Reaction.parse_formula(ccache.get_compound, equation_str)
result = equilibrator_calculator.standard_dg_prime(equilibrator_reaction)   # ΔG′° (kcal/mol)
ln_RI  = equilibrator_calculator.ln_reversibility_index(equilibrator_reaction)
# output row: rxn_id \t dG0_prime \t uncertainty \t ln_RI
```

Two points worth emphasizing:

- **The reaction ΔG′° is computed at the reaction level by eQuilibrator, not
  summed from the compound `.tbl`.** It uses the full component-contribution
  covariance, so the stored `deltagerr` is eQ's own reaction uncertainty (with
  cross-compound correlations), not a naïve quadrature of per-compound errors.
  (This is exactly the correlation the stats add-on flags as missing when
  propagating per-reaction marginal σ — `standard_dg_prime_multi` is what would
  recover it.)
- Reactions are tagged `EQC` (complete — all reagents mapped) or `EQP` (partial)
  in the notes
  ([`Retrieve_eQuilibrator_Reactions_Energies.py:123-144`](/scratch/ctaylor/ModelSEEDDatabase/Scripts/Thermodynamics/Retrieve_eQuilibrator_Reactions_Energies.py#L123-L144)).

### 3c. Writing the values into the database

`Update_Compound_eQuilibrator_Energies.py` and
`Update_Reaction_eQuilibrator_Energies.py` load those two `.tbl` files and stamp
`thermodynamics["eQuilibrator"] = [dg, dge, operator]` onto each
compound/reaction JSON. The reaction updater is a **direct lookup** (no
re-summation): it reads `MetaNetX_Reaction_Energies.tbl` and writes the triple,
computing the per-source operator from that source's own (dg, dge). eQ values
**overwrite** GC where present — the README states the GC energies are written
first, "then the energies from eQuilibrator (EQ), which, in most cases, take
precedence, are used to overwrite."

### 3d. Contrast with the GC (Jankowski) path

The legacy path is structurally different: compound ΔGf comes from **MFAToolkit
group-decomposition tables** (`ModelSEED/{KEGG,MetaCyc}_{Charged,Original}_MolAnalysis.tbl`),
matched by structure-alias strings; the reaction ΔG′° is then the
**stoichiometric sum** Σ ν·ΔGf with **uncertainty in quadrature**, and it is
**all-or-nothing** (one missing reagent → the reaction gets the sentinel
`10000000`). That is the GC vs EQ split that surfaces as the
`Estimated_Reaction_Reversibility_Report_{GC,EQ}.txt` pair.

---

## 4. Bottom line: the two eQ paths differ in identity resolution *and* energy vintage

Both eQ paths ultimately hand eQuilibrator an accession string that resolves
against the same Zenodo-backed compound cache — but they reach that string by
different routes, and that difference *is* the `Flamholz_2012` / `Beber_2022`
split:

| | `Flamholz_2012` (MSDB bundled) | `Beber_2022` (live eQ 3.0) |
|---|---|---|
| Where energies came from | precomputed `MetaNetX_Reaction_Energies.tbl` | computed live, `equilibrator-api 0.7.0` |
| Metabolite → eQ identity | ModelSEED **InChIKey → MetaNetX** structural match (3-tier skeleton fallback) | ModelSEED **alias accession** (KEGG→ChEBI→BiGG→MNX→InChIKey, first hit) |
| eQ lookup call | `ccache.get_compound(MNXM…)` | `cc.get_compound(acc)` |
| eQ cache vintage | whatever generated the bundled `.tbl` (~2019) | eQ 3.0 cache at `/scratch/ctaylor/eq_cache` |
| Direction rule | reversibility index, ln 1000 cutoff | **same** reversibility index, ln 1000 cutoff |

In both cases eQuilibrator only ever sees an accession (or InChIKey) it can
match to a pre-reconciled structure carrying a group decomposition; the ΔG comes
from the component-contribution model applied to that structure.

**Provenance caveat:** the committed
[`estimate_directions_literature.py:140`](/scratch/ctaylor/core_models_analysis/scripts/estimate_directions_literature.py#L140)
writes only **4** columns (it does not itself merge `Beber_2022`). The eq3
column is emitted separately to `results/rxn_directions_eq3_2022.tsv` and was
merged into the 5-column file as a follow-on step that is not in that one
script. The data is correct; the merge just isn't a single committed function.

---

## 5. Follow-up analyses worth doing

- **Component-contribution math, one level deeper:** the reactant- vs
  group-contribution split and how the covariance feeds `deltagerr` /
  `standard_dg_prime_multi`.
- **Identity-resolution vs energy-vintage attribution:** quantify how many panel
  reactions resolve to *different* eQ compounds between the alias-route (§2) and
  the InChIKey-route (§3b) — i.e. whether the 2012/2022 columns ever disagree
  because of identity resolution rather than energy vintage.

---

## File index

- Direction-column assembly (cols 1,2,4): `core_models_analysis/scripts/estimate_directions_literature.py`
- eQ 3.0 live column (col 3): `core_models_analysis/scripts/estimate_directions_eq3.py`
- MSDB compound energies (retrieve): `ModelSEEDDatabase/Scripts/Thermodynamics/Retrieve_eQuilibrator_Compound_Energies.py`
- MSDB reaction energies (retrieve): `ModelSEEDDatabase/Scripts/Thermodynamics/Retrieve_eQuilibrator_Reactions_Energies.py`
- MSDB updaters: `ModelSEEDDatabase/Scripts/Thermodynamics/Update_Compound_eQuilibrator_Energies.py`, `Update_Reaction_eQuilibrator_Energies.py`
- Structure map: `ModelSEEDDatabase/Biochemistry/Structures/MetaNetX/Structures_in_ModelSEED_and_eQuilibrator.txt`
- Output `.tbl` files: `ModelSEEDDatabase/Biochemistry/Thermodynamics/eQuilibrator/MetaNetX_{Compound,Reaction}_Energies.tbl`
- Output TSV: `core_models_analysis/results/reaction_directions_literature_vs_llm.tsv`
