# ModelSEED reaction balance flags: what they are, who carries them, how fixable they are

Built 2026-10-09 against ModelSEED `dev` @ `078a395f`
(`/scratch/ctaylor/tmp/devsnap_078a395f`) — confirmed still == `origin/dev`
HEAD on this date (unchanged since the 2026-10-06 merge), so this is current,
not stale.

**Scope note**: unlike
[`reports/strictDisagreementEnzymes/`](../strictDisagreementEnzymes/STRICT_DISAGREEMENT_ENZYMES.md),
this report is **database-wide** — all 48,384 non-obsolete reactions, not a
direction-disagreement subset. It grew out of a question in that report
("ModelSEED already knows" about broken equations, via its `status` field) —
this is the full explanation of that field, plus a first pass at how tractable
each kind of flagged problem actually is.

## What `status` is

It's a fully mechanical mass/charge balance check — no literature, no
judgment, pure arithmetic. The computing function is
`BiochemPy.Reactions.balanceReaction()`
(`Libs/Python/BiochemPy/Reactions.py`, lines 316–442, in the full
ModelSEEDDatabase checkout — this function isn't present in the data-only
snapshot, only in the code checkout). For each reaction it sums
`charge × coefficient` over every reagent (net charge) and
`element-atom-count × coefficient` per chemical element (net mass, per
element), after collapsing duplicate compound/compartment reagents. Values
within `±1×10⁻⁶` are rounded to exactly 0 to absorb floating-point noise.

`status` is a single string, but it's a `|`-joined combination of codes, and
some codes carry a `:`-delimited payload. The exact parsing convention —
confirmed from `Scripts/Validation/Validate_Reactions.py`, which relies on it
— is: split on `|`; a field's *type* is whatever precedes its first `:` (or
the whole field, if there's no `:`).

| code | meaning | payload |
|---|---|---|
| `OK` | mass- and charge-balanced | none |
| `CI:N` | **Charge Imbalance** — net charge differs by `N` | `N` is the literal signed net charge difference, not a severity scale |
| `MI:Elem:N[/Elem2:N2/...]` | **Mass Imbalance** — one or more elements don't cancel | per-element signed atom-count imbalance, `/`-joined when multiple elements are off; `R` denotes a generic/undetermined substituent group |
| `MI:...\|CI:...` | both at once | `\|`-joined |
| `CPDFORMERROR` | can't even check — a reagent has no parseable formula | short-circuits before mass/charge are computed (the two massless species, photon `cpd11632` and electron `cpd12713`, are whitelisted and don't trigger this) |
| `CK` | **curator-checked override** — a human has reviewed and is accepting/flagging the computed result (e.g. `OK\|CK\|CI:1`) | not emitted by the balancer itself; woven in by curation |
| `NB` | "not balanced" — the default placeholder on a brand-new reaction before the balancer has ever run | — |
| `EMPTY` | zero reagents | — |
| `Duplicate reagents` | the same compound+compartment appears as more than one distinct reagent and can't be auto-merged (e.g. polymer-type reactions); transient, usually converted to `NB` on the next rebalance pass | — |

Two related, auto-correction provenance codes — `HB` ("hydrogen-balanced",
when the balancer silently adjusts a proton coefficient to fix a near-miss)
and `WB` (the water-coefficient equivalent) — live in a **separate** `notes`
field, not `status`.

## Database-wide breakdown

Of the 48,384 non-obsolete reactions, every flag combination that actually
occurs:

![Horizontal bar chart on a log scale of every reaction status flag combination across 48,384 non-obsolete ModelSEED reactions: OK 33,472 (69.2%), CPDFORMERROR 6,781 (14.0%), CI+MI both 3,213 (6.6%), CI only 2,824 (5.8%), MI only 2,067 (4.3%), CI+CK 21, OK+CK 5, EMPTY 1. The two CK combinations are outlined to mark them as curator-reviewed.](figures/reaction_balance_flag_counts.png)

| | n | share |
|---|---:|---:|
| `OK` | 33,472 | 69.19% |
| `CPDFORMERROR` | 6,781 | 14.01% |
| `CI` + `MI` (both) | 3,213 | 6.64% |
| `CI` only | 2,824 | 5.84% |
| `MI` only | 2,067 | 4.27% |
| `CI` + `CK` (charge-imbalanced, curator reviewed anyway) | 21 | 0.04% |
| `OK` + `CK` (balanced, curator-confirmed) | 5 | 0.01% |
| `EMPTY` | 1 | 0.00% |

Two numbers worth holding onto: **total charge-imbalanced** (any `CI`) =
6,058 (2,824 + 3,213 + 21); **total mass-imbalanced** (any `MI`) = 5,280
(2,067 + 3,213); their union is 8,125 (16.8% of the database) — the baseline
figure used throughout the strict-disagreement report. And of those 8,125,
**only 21 (0.26%) have ever been reviewed by a curator** (`CK`). The rest —
8,104 reactions — are unexamined.

*(Figure source: `scripts/plot_reaction_balance_flags.py`, registered in
`scripts/figures.tsv` as `reaction_balance_flag_counts`. Data in
`results/reaction_balance_flags/flag_counts.tsv`.)*

## Groupings and features within each flag

### `CPDFORMERROR` (6,781) — a representation choice, not a bug

6,341 of the database's 45,708 compounds have no formula at all, and the
affected reactions concentrate hard around a specific kind of compound: a
generic or protein-bound redox cofactor that genuinely has no single
molecular formula, because it's a protein, not a small molecule.

| compound | reactions affected |
|---|---:|
| `Acceptor` (fully generic, unspecified) | 967 |
| oxidized/reduced `[NADPH--hemoprotein reductase]` | 638 each |
| oxidized/reduced `[2Fe-2S]-[ferredoxin]` | 176 each |
| `Fe(II)`/`Fe(III)-[cytochrome b5]` | 145 / 144 |
| reduced electron-transferring flavoprotein (ETF) | 78 |
| oxidized/reduced `[adrenodoxin]` | 60 each |
| oxidized `[electron-transfer flavoprotein]` | 60 |
| `Fe(II)`/`Fe(III)-[cytochrome c]` | 45 each |
| oxidized/reduced flavoprotein | 40 each |
| thioredoxin-dithiol | 36 |
| generic amino-acyl-`[protein]` residues, Biomass, mRNA, Collagen, LPS, Trypsin, ... | smaller, long tail |

**4,158 distinct compounds** are implicated across the 6,781 reactions; the
**top 10 alone account for 3,092 reaction-hits (46%)**.

### `CI` only (2,824) — one real bug plus an uninvestigated remainder

The charge-imbalance magnitude is tightly clustered, not random: `|CI|=2` is
61.4%, `|CI|=1` is 15.2%, `|CI|=4` is 11.9% — almost entirely small even
integers.

**778 (27.5%) are combinatorial "expansion" reactions** — auto-generated by
pairing a generic reaction template against many specific substrates (visible
in the `abbreviation` field, e.g. `...RXN.c.maizeexp.BCAA-dehydrogenase-...`).
Of those, **614 (79%) involve H⁺** (`cpd00067`), and the substrate list is a
recognizable family: NAD(H)/NADP(H)-dependent dehydrogenase templates
expanded across specific aldehydes/alcohols (decanal, citronellal,
geranial, eicosanal, cyclohexenecarbaldehyde, ...). That is the signature of
**one systematic bug in a reaction-generation template** — a proton
coefficient that wasn't recalculated correctly when the template was expanded
per substrate — not 778 independent errors.

`CI`-only reactions are notably *under*-represented among transport
reactions (0.3% vs. 13.3% of the whole database) — transport crossing a
compartment boundary is not a driver here, which rules out the obvious
"unbalanced co-transported ion" hypothesis.

### `MI` only (2,067) — about half are mechanically triageable

Matching each reaction's exact mass residual (read straight from its
`MI:Elem:N` string) against every compound's formula in the database: **945
(45.7%) have a residual that matches a real compound's formula exactly.**

| matched compound | n |
|---|---:|
| `Deaminated-Amine-Donors` *(generic — likely coincidental)* | 292 |
| `1-Haloalkane` *(generic — likely coincidental)* | 106 |
| `L-Rhamno-1,4-lactone` | 76 |
| Ethylene | 29 |
| `2,1' dianhydride` | 26 |
| Formaldehyde | 17 |
| Sucrose | 16 |
| O₂ | 14 |
| Propene | 12 |
| H₂O₂ | 9 |
| D-Glucose | 8 |

This is the database-wide generalisation of the "dropped products that
already exist in the database" pattern the literature pass found within the
1,370 strict disagreements (oxidised glutathione, ADP, CO₂, PPi, O₂...). The
top two matches are generic, simple-formula compounds that plausibly match by
coincidence rather than real identity (the same caveat applies here as there:
**formula matching cannot distinguish isomers**, so a match is a candidate,
not a confirmation) — but entries further down (Ethylene, Formaldehyde,
Sucrose, O₂, H₂O₂, D-Glucose) are specific and chemically plausible as
genuinely dropped reagents or products.

### `CI` + `MI` both (3,213) — the hardest group, and the biggest checkable one

- **61.3% (1,968) carry no EC number at all** — mostly orphan/unannotated
  reactions, not mechanistically characterised to begin with.
- **65.7% (2,111) have three or more distinct elements imbalanced at once**
  (24.8% have exactly five). A single dropped proton or missing small
  cofactor shows up as one or two elements off; three-plus elements off
  simultaneously is the signature of a more fundamentally wrong equation —
  the wrong compound substituted, or several errors stacked — not a
  coefficient slip.
- Only 3.0% are fractional-coefficient (biomass/macromolecule-template)
  pseudo-reactions, so that's not the explanation either.
- Transport share (4.4%) is, like `CI`-only, below the database baseline
  (13.3%).

This group most resembles the 1,370 strict-disagreement reactions in kind —
it most likely needs the same case-by-case literature review, not a
mechanical fix — and at 3,213 reactions it is the single largest flagged
population after `CPDFORMERROR`.

## How much of each group reduces to a named cause

![Horizontal stacked bar chart showing, for each of the four non-OK flag groups, how many reactions match an identified concrete pattern (green) versus remain an unexplained residual (grey): CPDFORMERROR 3,092 of 6,781 (46%) match the top-10 generic cofactor placeholders; CI+MI both 1,102 of 3,213 (34%) are single/double-element imbalances; CI only 778 of 2,824 (28%) are the combinatorial-expansion template bug; MI only 945 of 2,067 (46%) have a residual matching a known compound.](figures/reaction_balance_group_tractability.png)

*(Figure source: `scripts/plot_reaction_balance_group_tractability.py`,
registered in `scripts/figures.tsv` as
`reaction_balance_group_tractability`. Data in
`results/reaction_balance_flags/group_tractability.tsv`. Note:
`CPDFORMERROR`'s green share is a **concentration** finding, not a
**fixability** one — see below; it is not directly comparable to the other
three bars, which really are "easier to fix" shares.)*

## Solutions: what it would actually take to fix each group

| group | n | what's identified | how easy to fix |
|---|---:|---|---|
| `CPDFORMERROR` | 6,781 | top 10 compounds = generic protein-bound cofactors | **Not fixable as a "bug."** A real ferredoxin or cytochrome doesn't have one formula across organisms — there's nothing correct to enter. Leave these flagged; they're an accurate reflection of generic chemistry, not an error. The one addressable slice is the 967 reactions using the fully generic `Acceptor` placeholder — identifying each one's real physiological electron acceptor is a literature-research project (same effort class as the 1,370-reaction strict-disagreement project), not a database patch. |
| `CI` only | 2,824 | 778 (27.5%) = one templating bug (NAD(P)-dehydrogenase expansions missing an H⁺ adjustment) | **Easy for that 27.5%**: find the one reaction-generation script responsible for the `.{organism}exp.` combinatorial expansion, correct its proton-coefficient logic, regenerate — likely resolves several hundred reactions in a single patch, not 778 manual edits. The remaining ~2,046 non-templated cases are an uninvestigated, probably heterogeneous population — would need the same per-reaction triage as `MI`-only below. |
| `MI` only | 2,067 | 945 (45.7%) = residual matches a known compound's formula | **Moderate.** The formula match is free and mechanical — it narrows 2,067 open-ended "what's wrong here" questions down to 945 "is it plausibly *this* compound" questions. Each candidate still needs a quick chemical-plausibility check (context: does this compound make sense in this reaction's pathway?) before being accepted, to rule out coincidental isomer matches — but that's a much smaller task than starting from an unidentified imbalance. |
| `CI` + `MI` both | 3,213 | 1,102 (34.3%) = only 1–2 elements off (simpler subset); the rest (65.7%) have 3+ elements off and 61.3% have no EC at all | **Hard.** The multi-element-imbalance majority looks like genuinely wrong equations, not coefficient bugs, and most aren't even enzyme-annotated to start research from. This is the group that most needs the strict-disagreement project's approach — individual literature review — and at 3,213 reactions (before even subtracting the "simpler" 1,102), it's the most expensive one to clear. |

**Net picture**: of the ~14,906 non-`OK`, non-`CPDFORMERROR` reactions (the
8,125 imbalanced + 6,781 uncheckable), **maybe 800–1,000 are a quick,
high-leverage fix** (the `CI`-only templating bug, fixed once at the
generator level), **roughly 1,000–2,000 more are plausible with a quick
per-reaction confirmation pass** (the `MI`-only dropped-product matches), and
**the remaining several thousand** — most of `CI+MI` both, the
non-templated `CI`-only remainder, and all of `CPDFORMERROR` apart from the
967 generic-`Acceptor` cases — **need the same kind of individual literature
review already done for the 1,370 strict disagreements**, not a mechanical
patch.

## Files

```
results/reaction_balance_flags/
  flag_counts.tsv               every status-token combination + count (figure 1 source)
  cpdformerror_compounds.tsv    every missing-formula compound, ranked by reactions affected
  cpdformerror_analysis.json    full CPDFORMERROR breakdown
  ci_only_analysis.json         charge-imbalance magnitude distribution, expansion-template overlap
  mi_only_analysis.json         mass-residual-vs-known-compound-formula matches
  both_analysis.json            EC-annotation and imbalance-multiplicity profile of the CI+MI group
  group_tractability.tsv        identified-pattern vs. residual share per group (figure 2 source)
```

## Reproducing

```bash
python3 scripts/analyze_reaction_balance_flags.py
python3 scripts/regen_figures.py reaction_balance_flag_counts reaction_balance_group_tractability
```

## Caveats

- **This is a first pass, not a verified fix list.** Every "identified
  pattern" above is a *concentration* finding (many reactions share one
  recognisable signature) — none of the candidate fixes (the templating bug,
  the dropped-product matches) have been applied or confirmed against the
  actual `balanceReaction()` code paths that would need changing.
- **Formula matching cannot distinguish isomers or confirm biological
  relevance.** A matched compound in the `MI`-only analysis is a candidate
  explanation, not a confirmed one — the same caveat the strict-disagreement
  report already carries for its own "dropped products" finding.
- **The `CI`-only templating-bug hypothesis is inferred from co-occurring
  signals** (expansion-pattern ID, H⁺ involvement, substrate family), not
  from reading the actual generator script. Confirming it would mean finding
  and reading whichever script produces `.{organism}exp.`-tagged reactions.
- **Scope is the pinned snapshot only.** As with every other report in this
  project, the unpinned local checkout (`/scratch/ctaylor/ModelSEEDDatabase`)
  has diverged from `origin/dev` and must not be mixed in with these numbers.
