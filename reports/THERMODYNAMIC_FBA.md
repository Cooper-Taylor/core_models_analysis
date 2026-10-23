# Thermodynamic FBA (tFBA) in ModelSEEDPy — what it is, why it's broken, and what the data does to it

*Investigated 2026-10-23. Scripts in `scripts/tfba/`; regenerate the legacy database
layout with `scripts/tfba/build_legacy_msdb.sh`.*

---

## The short version

ModelSEEDPy really does contain a thermodynamic-FBA implementation. It is a faithful
version of the standard TMFA method. But:

1. **It cannot run today.** It breaks on the first compound it looks up. The cause is not a
   bug in the thermodynamics — it's that the ModelSEED *database* changed its file layout in
   May 2022 and the thermo code was never updated.
2. **It last worked against a ModelSEEDDatabase checkout from on or before 2022-05-12.**
3. **Even once you make it run, the compound ΔG numbers it reads are the wrong input.** The
   ModelSEED compound energies are stored in two different bookkeeping conventions mixed
   together. tFBA adds compound energies up, so the mix doesn't cancel — it produces garbage
   on about 30% of reactions.
4. **Net effect on core models:** as shipped, nothing (a unit bug cancels the constraint out).
   Fix the unit bug and you get a ~20% average growth loss that is *entirely an artifact* of
   the bad compound data. Use the database's reaction-level ΔG instead and the effect is zero.

---

## 1. How tFBA works, in plain terms

Ordinary FBA only knows two things: mass must balance, and fluxes have bounds. Nothing in it
prevents a reaction from running uphill energetically.

tFBA adds one idea: **give every metabolite a price.**

- Each metabolite *j* gets a number `p_j` — think of it as its energy level.
- A reaction's ΔG is just the price difference across it: `ΔG = Σ (stoichiometry × price)`.
- A reaction may only carry flux in the direction where its ΔG is negative (downhill).

That's it. The rest is plumbing to make "may only" expressible to a solver.

### The three packages

| Package | What it adds | In words |
|---|---|---|
| `RevBinPkg` | a 0/1 variable `z` per reaction | the direction switch: `z=1` allows forward flux, `z=0` allows reverse |
| `SimpleThermoPkg` | a `potential` per metabolite, a `thermo` constraint per reaction | ties the switch to the sign of ΔG: `z=1` forces `ΔG ≤ 0`, `z=0` forces `ΔG ≥ 0` |
| `FullThermoPkg` | the actual numbers | says what a metabolite's price *is* |

`FullThermoPkg` is the only one that touches real data. It sets:

```
price = ΔGf           (formation energy, from the ModelSEED database)
      + RT·ln(c)      (concentration — the solver picks c between 1 µM and 20 mM)
      + charge·F·ψ    (membrane potential, −160 mV across the cell membrane)
      + error         (slack, ±5 kJ/mol, to absorb uncertainty in ΔGf)
```

You call `FullThermoPkg`; it pulls the other two in automatically.

### Why that gives the solver freedom

A reaction isn't simply locked by its ΔG°. The solver can choose concentrations. Each
metabolite's concentration range (1 µM → 20 mM) is worth about **24.5 kJ/mol** of wiggle
room, plus **±5 kJ/mol** of error slack. So a reaction with four metabolites can shift its
ΔG by roughly ±120 kJ/mol before it is genuinely stuck.

**This is the number to keep in your head.** A reaction is only *hard-locked* when its
standard ΔG exceeds that budget. Most real reactions don't.

---

## 2. Why it isn't working right now

### The mechanism

`modelseedpy.biochem.from_local()` has **two code paths**, and it picks between them by
looking at what your ModelSEEDDatabase folder contains:

```python
def from_local(database_path):
    contents = os.listdir(f"{database_path}/Biochemistry/")
    if "compounds.tsv" in contents:
        return from_local_old(database_path)   # old path
    ...                                        # new path
```

| | old path | new path |
|---|---|---|
| triggered when | `Biochemistry/compounds.tsv` exists | it doesn't |
| returns | `ModelSEEDBiochem` | `ModelSEEDDatabase` |
| compound objects | `ModelSEEDCompound` | `ModelSEEDCompound2` |
| has `.get_seed_compound()` | **yes** | **no** |
| compound energy attribute | `.deltag` | `.delta_g` |

`FullThermoPkg` was written against the **old** path. It calls
`modelseed_api.get_seed_compound(msid)` and then reads `mscpd.deltag`. Both of those exist
only on the old classes.

Today's ModelSEEDDatabase has no `compounds.tsv` — the biochemistry is sharded into
`compound_00.tsv` … `compound_49.tsv`. So `from_local` takes the **new** path, and
`FullThermoPkg` dies immediately:

```
AttributeError: 'ModelSEEDDatabase' object has no attribute 'get_seed_compound'
```

Fix that one and you hit the second:

```
AttributeError: 'ModelSEEDCompound2' object has no attribute 'deltag'
```

Two API breaks, both from the same 2022 rewrite. Neither is in the thermodynamics.

### It's worse than the file layout — the schema changed too

Reassembling `compounds.tsv` and `reactions.tsv` from the shards is *not* enough. The old
loader also wants a `direction` column in `reactions.tsv` that the modern schema dropped:

```
KeyError: 'direction'
```

Pre-2022 `reactions.tsv` had both `reversibility` (col 9) and `direction` (col 10). The
sharded schema keeps only `reversibility`. So there are two layers of drift between the
thermo code and the current database.

---

## 3. When it last worked

| date | what happened |
|---|---|
| 2021-09-16 | last substantive change to the thermo code |
| **2022-05-12** | **MSDB commit `2f5da32e` — last commit with the legacy layout. This is the newest database the thermo code can read.** |
| 2022-05-20 | MSDB commit `8bf2281b`, *"Massive re-orientation towards breaking down biochemistry into chunks"* — deletes `compounds.tsv`, drops the `direction` column. **tFBA breaks here.** |
| 2022-09-16 | `f35ac68f` reformats the thermo files with `black` — no logic change |
| 2022-10-21 | modelseedpy adds the new loader + `from_local` dispatch |
| 2026-10-23 | still broken; `main`, `dev` and the installed 0.4.2 are byte-identical |

**So it has been broken for about three and a half years.**

### Why nobody noticed

- Both test files (`tests/fbapkg/test_full_thermo_core_model.py`,
  `test_simple_thermo_core_model.py`) have their entire contents wrapped in a triple-quoted
  string. They are collected by CI and contribute **zero** tests.
- Both in-package callers would crash *before* reaching the thermodynamics:
  `community/dfbapkg.py` calls `build_package()` with no arguments against a signature that
  requires one; `community/mscommunity.py` passes the key `modelseed_db_path` where the code
  requires `modelseed_path`.
- No GitHub issue or PR mentions `tfba`, `TMFA`, `fullthermo` or `deltaG`. Nobody has
  reported it.

### Verified: the stock code still runs against the legacy layout

I rebuilt the pre-2022 layout from the current sharded data (concatenated the shards into
`compounds.tsv`/`reactions.tsv`, re-derived the `direction` column) with
`scripts/tfba/build_legacy_msdb.sh`, and ran **unmodified** `FullThermoPkg` against it:

```
plain FBA growth : 115.0000  (optimal)
build_package OK : 7.8s   <-- STOCK CODE, NO SHIM
  potc constraints (metabolites priced) : 157
  thermo constraints (reactions coupled): 180
  logconc vars: 177   dgerr vars: 179   revbin vars: 213
```

So the code itself is intact. It's a data-compatibility failure, not rot.

---

## 4. What inputs it needs

### The call

```python
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager

pkgmgr = MSPackageManager.get_pkg_mgr(model)
pkgmgr.getpkg("FullThermoPkg").build_package({"modelseed_path": "/path/to/ModelSEEDDatabase"})
```

### Parameters

| parameter | required | default | meaning |
|---|---|---|---|
| `modelseed_path` | **yes** | — | path to a ModelSEEDDatabase checkout (must be legacy layout, see above) |
| `default_min_conc` | no | `1e-6` | lowest metabolite concentration, M |
| `default_max_conc` | no | `0.02` | highest metabolite concentration, M |
| `default_max_error` | no | `5` | ΔGf slack per metabolite, kJ/mol |
| `custom_concentrations` | no | `{}` | pin specific metabolites: `{"cpd00027_e0": [0.02, 0.02]}` |
| `custom_deltaG_error` | no | `{}` | per-metabolite error override |
| `compartment_potential` | no | `{"e0": 0, "c0": -160}` | membrane potential in **mV**; `e0` must stay 0 |
| `temperature` | no | `298` | K |
| `filter` | no | `None` | list of reaction ids to constrain; `None` = all |
| `dgbin` | no | `False` | adds binaries to maximise the number of reactions with a defined ΔG sign |
| `infeasible_model` | no | `False` | forces `dgbin` on; for diagnosing infeasibility |

### What your model must look like

- **Metabolite ids must start with `cpdXXXXX`** — the ModelSEED id is parsed off the front
  with a regex. Anything else is silently skipped (no price, no constraint).
- **Compartment ids must be `c0` / `e0`** to match `compartment_potential`. If your model
  uses `c` / `e`, the membrane-potential term silently becomes zero. (The core models use
  `c0`/`e0`, so they're fine.)
- Apply your medium to the exchange reactions **before** building the package.
- Exchange/sink/demand reactions (`EX_`, `SK_`, `DM_`) are skipped automatically.

### Gotchas in the parameter handling

- `default_deltaG_error()` returns `{"cpd00067": 0}` — intended to give the proton zero error.
  But the lookup uses the full metabolite id (`cpd00067_c0`), so **this override never fires**
  and H⁺ gets the full ±5 kJ/mol.
- The database's own per-compound `deltagerr` is **never read**. Everything gets a flat
  ±5 kJ/mol regardless of how uncertain its ΔGf actually is.
- `SimpleThermoPkg.build_package` **destroys `model.objective`** and never restores it. Save
  and restore it yourself around the call.
- It runs **two MILP solves per reaction** to compute a constant that has a closed form. On a
  213-reaction core model that's ~360 solves. See `scripts/tfba/fast_thermo.py` for a drop-in
  replacement that computes the same number directly.

---

## 5. What outputs to expect

### What gets added to the model

| kind | name | one per | units | range |
|---|---|---|---|---|
| variable | `logconc` | metabolite (not water) | ln(M) | `ln(1e-6)` … `ln(0.02)` |
| variable | `dgerr` | metabolite | kJ/mol | −5 … +5 |
| variable | `potential` | metabolite | kJ/mol | set by `potc` |
| variable | `revbin` | reaction | binary | 0 or 1 |
| constraint | `potc` | metabolite with a valid ΔGf | — | equality: price = ΔGf + RT·ln c + zFψ + err |
| constraint | `thermo` | internal reaction | — | ties ΔG sign to `revbin` |
| constraint | `revbinF` / `revbinR` | reaction | — | ties `revbin` to flux direction |

**Your model becomes a MILP.** That's the big practical consequence — it is no longer a
fast LP. On the E. coli core model (213 reactions) that's 213 binaries.

### What you read back out

```python
sol = model.optimize()
sol.objective_value                      # growth, as usual
sol.fluxes                               # fluxes, as usual

sp = pkgmgr.getpkg("SimpleThermoPkg")
ft = pkgmgr.getpkg("FullThermoPkg")

# the concentration the solver chose for a metabolite, in M
math.exp(ft.variables["logconc"]["cpd00002_c0"].primal)

# the ΔG the solver assigned to a reaction, kJ/mol
sum(s * sp.variables["potential"][m.id].primal for m, s in rxn.metabolites.items())
```

### What to expect numerically

- **Growth ≤ plain FBA growth, always.** tFBA only removes solutions.
- **Coverage is partial.** On the E. coli core model, 157 of 179 metabolites get a price;
  19 of 180 internal reactions touch a metabolite with no ΔGf.
- **A metabolite with no ΔGf does not disable the constraint — it silently makes it
  non-binding.** The price variable stays free over ±100,000 kJ/mol, so any reaction touching
  it can have whatever ΔG it needs. This is the most dangerous failure mode: you get a
  thermodynamically "constrained" model that quietly isn't.
- **The MILP may simply not solve.** The default price bounds are ±100,000 kJ/mol, which makes
  the big-M in the `thermo` constraint about 10⁶ — a hopeless LP relaxation for GLPK (the cobra
  default here; HiGHS is in optlang's list but is not wired into this cobra build). On the
  213-reaction E. coli core model neither run converged: **25 minutes** with the default bounds
  and **14 minutes** with the implied bounds written onto the price variables, both still
  running when abandoned. `model.solver.configuration.timeout` was **not honoured** by GLPK.
  `tighten_potentials()` in `scripts/tfba/tfba_experiment.py` removes no feasible point and is
  worth applying, but it is not sufficient on its own — plan on a commercial solver, a `filter`
  restricting the constraint to reactions you care about, or both.

---

## 6. The data problem (the part that actually matters)

### Two conventions, mixed

ModelSEED stores compound formation energies in **two different bookkeeping systems**:

| | Convention A (Henry 2010 / MFAToolkit) | Convention B (Alberty / eQuilibrator) |
|---|---|---|
| hydrogen atoms | left **in** the compound energy | subtracted out at compound level |
| ΔfG(H⁺) | −9.5 kcal/mol | 0 |
| water | −56.687 kcal/mol | −37.6 kcal/mol |

They are mathematically equivalent *for a balanced reaction, if you use one of them
consistently*. The difference per compound is `n_H × 9.539 kcal/mol`.

The top-level `deltag` field in ModelSEED is **a documented mix of both** — this is written
up in `Biochemistry/Thermodynamics/README.md` and
`Papers/NAR_Update_2026/data/gc_drift_source_2026-08-06.md`.

**Why this is fatal specifically for tFBA:** tFBA is the one method that adds *compound*
formation energies together. Anything that reads reaction-level ΔG is unaffected. tFBA is
exactly the consumer the mix breaks.

### The clearest example

| compound | formula | stored `deltag` | convention |
|---|---|---|---|
| cpd00101 ribose-5-phosphate | C5H9O8P | −394.12 | **A** |
| cpd00171 ribulose-5-phosphate | C5H9O8P | −294.47 | **B** |

Two isomers. Identical formula. **417 kJ/mol apart** — purely bookkeeping, not chemistry.

Consequence for ribose-5-phosphate isomerase (rxn00777), which interconverts exactly those two:

```
tFBA's ΔG (summed from compounds) : +416.9 kJ/mol
the database's own reaction ΔG    :   +4.4 kJ/mol
```

tFBA sees an impossible reaction and hard-locks it in reverse. rxn00777 appears in 395 of 400
core models sampled. The same compound also poisons rxn00770 (PRPP synthetase) and rxn01200
(transketolase) — i.e. the pentose phosphate pathway.

### How widespread

Comparing the compound-summed ΔG against the database's own reaction-level ΔG, across all
26,005 MSDB reactions that have both:

| disagreement | reactions | share |
|---|---|---|
| < 5 kJ/mol | 16,743 | 64.4% |
| 5 – 50 | 1,346 | 5.2% |
| 50 – 200 | 2,222 | 8.5% |
| 200 – 1000 | 4,247 | 16.3% |
| > 1000 | 1,447 | 5.6% |

**30.4% disagree by more than 50 kJ/mol.** Every core model sampled contains at least one.

### Plus a unit bug

`fullthermopkg.py` converts the formation energy with:

```python
constant = mscpd.deltag / calorie + ...
```

ModelSEED `deltag` is in **kcal/mol** (ATP = −548.85 kcal/mol = −2296 kJ/mol). Converting
kcal to kJ is `× 4.184`, not `÷ 4.184`. The shipped expression makes ΔGf **17.5× too small**
relative to the `RT·ln c` and error terms, which *are* correctly in kJ.

Two 2021 commits — *"remedied orders of magnitude"* and *"corrected kilo orders of
magnitude"* — edited that exact line and never touched this term.

---

## 7. Effect on core model flux

80 core models, glucose-minimal aerobic medium, 49 with non-zero baseline growth
(`scripts/tfba/lock_sweep.py`). Each row applies the direction locks implied by a different ΔG
source, then re-optimises:

| ΔG source | locks/model | models losing growth | mean growth vs FBA |
|---|---|---|---|
| compound-sum, **as shipped** (÷4.184) | 1.8 | 0/49 | **100%** |
| compound-sum, **unit-fixed** (×4.184) | 16.6 | **27/49** (2 to zero) | **79.5%** |
| reaction-level `deltag` | 2.6 | 0/49 | **100%** |
| eQuilibrator reaction values | 1.4 | 0/49 | **100%** |

Three conclusions:

1. **As shipped, tFBA does nothing.** The unit bug shrinks ΔG by 17.5× so almost nothing
   exceeds the concentration budget. The constraint is inert.
2. **Fixing the unit bug makes it worse, not better.** It unmasks the corrupted compound data:
   ~20% average growth loss and two models killed outright. All artifact.
3. **Correct thermodynamics costs these core models nothing** at the growth level. Reaction-level
   ΔG has *higher* coverage (93% vs 90%) and still produces 4.7× fewer locks — so this is real
   signal, not missing data.

Worth noting: 23 of the 27 unit-fixed locks are single-compartment reactions where the
membrane-potential term cancels exactly, so they aren't an artifact of my handling of ψ.

---

## 8. What to do instead

**Don't drive tFBA from compound formation energies.** Drive it from reaction-level ΔG:

```
ΔG_i = ΔG°_i + RT · Σ_j s_ij · ln(c_j)
```

where `ΔG°_i` comes from `reaction['deltag']` or `reaction['thermodynamics']['eQuilibrator']`.

This is standard TFA anyway, and it sidesteps the convention mix completely, because the
reaction-level values are internally consistent — whoever computed them did the sum within a
single convention.

Keep the rest of the formulation (the `revbin` coupling, the concentration variables) exactly
as ModelSEEDPy has it; that part is correct.

### If you want to use the package as-is anyway

1. Point it at a ModelSEEDDatabase checkout from **2022-05-12 or earlier**, or rebuild the
   legacy layout with `scripts/tfba/build_legacy_msdb.sh`.
2. Change `deltag / calorie` → `deltag * calorie`.
3. Save and restore `model.objective` around `build_package`.
4. Swap in `scripts/tfba/fast_thermo.py` to avoid 2 MILP solves per reaction.
5. Tighten the price bounds before solving.
6. **Accept that steps 1–5 still leave you with the 30% corrupted-compound problem.**

---

## 9. Caveats on this analysis

- The 80-model sweep fixes directions **reaction by reaction**. It ignores the coupling
  between reactions that share a metabolite price, which is where tFBA's real value lives —
  including elimination of thermodynamically infeasible internal loops / EGCs. The full MILP
  captures that; the sweep does not. So the sweep is a **lower bound** on the effect.
- The sweep omits the membrane-potential term, which matters only for the 4 transport
  reactions among the 27 locks.
- "Hard-locked" uses the concentration + error budget described in §1. Reactions below that
  threshold are left free, which is correct but conservative.

---

## Files

| file | what it does |
|---|---|
| `scripts/tfba/check_convention.py`, `check_convention2.py` | classify each compound's `deltag` as Convention A or B |
| `scripts/tfba/convention_error.py` | per-reaction error from the convention mix |
| `scripts/tfba/validate_drg.py` | compound-summed ΔG vs the DB's reaction-level ΔG, one model |
| `scripts/tfba/dbwide_drg_scan.py` | the same across all 26,005 MSDB reactions |
| `scripts/tfba/predict_locks.py` | which reactions get hard-locked, as-shipped vs unit-fixed |
| `scripts/tfba/lock_sweep.py` | the 80-model growth sweep (§7) |
| `scripts/tfba/tfba_experiment.py` | runs the real packages with API shims + the unit switch |
| `scripts/tfba/fast_thermo.py` | closed-form replacement for the per-reaction MILP loop |
| `scripts/tfba/stock_run.py` | **unmodified** FullThermoPkg against the legacy layout |
| `scripts/tfba/build_legacy_msdb.sh` | rebuilds the pre-2022 database layout from the current shards |
| `scripts/tfba/_paths.py` | `MSDB_ROOT` / `CORE_MODELS_ANALYSIS_DIR` resolution |

Results land in `results/tfba/` (gitignored). Nothing in `ModelSEEDDatabase/` or
`core_models_kegg2/` is modified — the legacy layout is symlinks plus concatenated copies.
