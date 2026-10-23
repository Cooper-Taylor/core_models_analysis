# Direction-change target reactions: presence and FBA activity

76 ModelSEED reactions whose direction changes — 56 with an explicit change
(`= -> >` or `= -> <`) and 20 labelled "no change" in the source list but which
change because of thermodynamic uncertainty. All 76 are currently reversible
(`=`), so every change is reversible → irreversible.

Measured across all three model panels with the repo's own `cma` model-set API
(model sets, media, biomass policy, seed-key mapping). **No ModelSEEDpy fbapkg
code is involved.**

## Headline

| panel | models | growers | targets present | targets active |
|---|---:|---:|---:|---:|
| core (`core_kegg2`, kbase_complete) | 5,683 | 3,461 | **5 / 76** | **3 / 76** |
| full / GSM (`ms2_gsm`, glucose minimal) | 5,420 | 4,347 | **36 / 76** | **27 / 76** |
| full / GSM (`ms2_gsm_auxo`, auxotrophy) | 5,420 | 5,397 | **35 / 76** | **27 / 76** |
| **union** | 16,523 | | **36 / 76** | **30 / 76** |

40 of the 76 appear in no model at all. 6 more are present but never carry
flux. Split by category: explicit change 27/56 present and 22/56 active;
"no change" 9/20 present and 8/20 active.

## Activity is measured on pFBA, not plain FBA

Plain FBA picks an arbitrary vertex among alternate optima, so whether a given
reaction carries flux in it is not reproducible. pFBA is used throughout, the
same convention as `build_flux_distributions.py`, with `ACTIVE_TOL = 1e-9`.
The plain-FBA answer is recorded alongside in `active_fba`: it gives 31 rather
than 30, the extra one (rxn15230) being an alternate-optimum artefact.

## Files

| file | rows | what |
|---|---:|---|
| `direction_change_targets_by_reaction.tsv` | 76 | one row per target reaction, per-panel present/growing/active counts and flags |
| `figures/direction_change_funnel.png` | — | the three-stage drop, split by category |
| `figures/direction_change_by_reaction.png` | — | per-reaction nested counts, three facets, tall |
| `figures/direction_change_by_reaction_horizontal.png` | — | the same, wide: reactions along x, facets stacked |
| `figures/direction_change_by_reaction_horizontal_{core,gsm,gsm_auxo}.png` | — | one wide figure per model set |
| `direction_change_targets_by_model.tsv` | 154,936 | one row per (panel, model, reaction) where the reaction is present, flagged for flux |
| `presence_activity.json` | — | the raw aggregate the TSVs are built from |

`per_model_records.json` is the intermediate the by-model TSV is built from.
It is gitignored; regenerate with `--dump-records`.

### `direction_change_targets_by_reaction.tsv`

`reaction_id`, `change`, `name` (from MSDB), then a
`present_<panel>` / `growing_<panel>` / `active_<panel>` triple for each of
`core`, `gsm`, `gsm_auxo`, then `present_total`, `growing_total`,
`active_total` summed over panels, then `present_any`, `active_any` (0/1),
`n_panels_present`, `n_panels_active`.

The three counts nest: **active ≤ growing ≤ present**, checked on every row.

* `present_<panel>` — models that contain the reaction at all.
* `growing_<panel>` — of those, how many have non-zero biomass. A non-grower
  has no flux solution, so this is the denominator `active` should be read
  against, not `present`. The gap matters: rxn00011 is in 5,352 core models
  but only 3,394 of them grow.
* `active_<panel>` — of those, how many carry flux through it under pFBA.

### `direction_change_targets_by_model.tsv`

`panel`, `model_id`, `grows`, `reaction_id`, `change`, `present`,
`active_pfba`, `active_fba`. Only rows where the reaction is present in that
model are emitted.

## Per-reaction counts

Generated from `direction_change_targets_by_reaction.tsv`; that file is the
source of truth and carries all 76 reactions plus totals. Each cell is
**present / growing / active** — models containing the reaction, of those how
many grow, of those how many carry flux under pFBA.

| reaction | change | core P/G/A | gsm P/G/A | auxo P/G/A | name |
|---|---|---|---|---|---|
| `rxn00011` | = -> < | 5352/3394/2906 | 5286/4261/2921 | 5158/5138/1290 | pyruvate:thiamin diphosphate acetaldehydetrans |
| `rxn00874` | = -> < | 4110/3050/0 | 4003/3266/6 | 4003/3988/0 | butanoyl-CoA:acetyl-CoA C-butanoyltransferase |
| `rxn00441` | = -> > | 3788/3032/714 | 3787/3424/1505 | 3788/3777/3214 | R00621 |
| `rxn01961` | = -> > | 0/0/0 | 5003/4093/154 | 4999/4981/139 | Xanthosine 5'-phosphate phosphohydrolase |
| `rxn01145` | = -> > | 0/0/0 | 4994/4090/99 | 4993/4973/76 | Thymidylate 5'-phosphohydrolase |
| `rxn03909` | = -> > | 0/0/0 | 5151/4137/4137 | 4557/4543/4231 | 1-Deoxy-D-xylulose-5-phosphate pyruvate-lyase  |
| `rxn01300` | = -> > | 0/0/0 | 4558/3628/1792 | 4468/4449/2496 | ATP:L-homoserine O-phosphotransferase |
| `rxn00330` | = -> < | 2544/2371/201 | 2680/2447/2237 | 2587/2579/1254 | acetyl-CoA:glyoxylate C-acetyltransferase (thi |
| `rxn05567` | = -> > | 0/0/0 | 3660/2748/1 | 3661/3645/0 | Galactitol transport via PEP:Pyr PTS |
| `rxn05485` | *uncert.* | 0/0/0 | 3641/2729/3 | 3641/3625/1 | N-Acetyl-D-glucosamine transport via PEP:Pyr P |
| `rxn00342` | = -> > | 0/0/0 | 3531/2803/0 | 3533/3524/611 | L-Asparagine amidohydrolase |
| `rxn00950` | *uncert.* | 0/0/0 | 3039/2435/13 | 3039/3026/43 | Cystathionine L-homocysteine-lyase (deaminatin |
| `rxn02275` | *uncert.* | 0/0/0 | 2738/2299/19 | 2739/2726/44 | 4-Guanidinobutanamide amidohydrolase |
| `rxn15230` | *uncert.* | 0/0/0 | 2481/1800/0 | 2480/2468/0 | Lactose galactohydrolase |
| `rxn00985` | *uncert.* | 4186/2790/0 | 310/256/7 | 339/339/66 | ATP:propanoate phosphotransferase |
| `rxn00274` | = -> < | 0/0/0 | 2277/2032/809 | 2180/2176/680 | Acetyl-CoA:glycine C-acetyltransferase |
| `rxn00566` | = -> > | 0/0/0 | 2170/1619/121 | 2127/2114/673 | L-Cysteine hydrogen-sulfide-lyase (deaminating |
| `rxn03919` | = -> > | 0/0/0 | 1977/1851/0 | 1977/1974/0 | D-glycero-beta-D-manno-heptose 1,7-bisphosphat |
| `rxn00606` | *uncert.* | 0/0/0 | 1392/906/149 | 1349/1341/226 | alpha,alpha-Trehalose-6-phosphate phosphogluco |
| `rxn00362` | *uncert.* | 0/0/0 | 1045/1024/80 | 1075/1074/144 | Cytidine-5'-monophosphate phosphoribohydrolase |
| `rxn00456` | = -> > | 0/0/0 | 1028/880/96 | 1034/1029/230 | L-methionine methanethiol-lyase (deaminating;2 |
| `rxn01027` | = -> < | 0/0/0 | 786/689/0 | 786/783/0 | 5,6-Dihydrothymine:NAD+ oxidoreductase |
| `rxn00600` | = -> < | 0/0/0 | 684/498/498 | 627/621/621 | Succinyl-CoA:glycine C-succinyl-transferase(de |
| `rxn00202` | = -> < | 0/0/0 | 708/496/125 | 537/535/89 | acetyl-CoA:2-oxoglutarate C-acetyltransferase  |
| `rxn01630` | *uncert.* | 0/0/0 | 435/352/101 | 371/370/1 | 5-Aminopentanamide amidohydrolase |
| `rxn08571` | = -> > | 0/0/0 | 330/327/0 | 330/330/0 | D-galactose 1-phosphatase |
| `rxn17641` | = -> > | 0/0/0 | 157/155/124 | 157/156/124 | prenyltransferase |
| `rxn11859` | = -> > | 0/0/0 | 109/98/80 | 83/82/64 | R08332 |
| `rxn00680` | = -> < | 0/0/0 | 155/114/114 | 5/5/5 | propanoyl-CoA:glyoxylate C-propanoyltransferas |
| `rxn05035` | = -> > | 0/0/0 | 126/117/117 | 6/6/6 | O-phospho-L-serine:hydrogen-sulfide 2-amino-2- |
| `rxn04712` | = -> > | 0/0/0 | 23/22/0 | 23/23/0 | R06903 |
| `rxn09433` | = -> > | 0/0/0 | 6/5/5 | 9/8/8 | fatty acid synthase (n-C16:0) |
| `rxn01028` | = -> < | 0/0/0 | 6/6/0 | 6/6/0 | 5,6-Dihydrothymine:NADP+ oxidoreductase |
| `rxn01317` | *uncert.* | 0/0/0 | 2/0/0 | 4/4/4 | N-Acetylneuraminate 9-phosphate phosphohydrola |
| `rxn09431` | = -> > | 0/0/0 | 1/0/0 | 5/5/5 | fatty-acyl-CoA synthase (n-C16:0CoA) |
| `rxn01456` | = -> > | 0/0/0 | 1/1/1 | 0/0/0 | D-galacturonate:NAD+ 1-oxidoreductase |

### The 40 reactions in no model at all


`rxn00020`  `rxn00273`  `rxn00997`  `rxn00998`  `rxn01009`  `rxn01139`  `rxn01195`  `rxn01196`  `rxn01408`  `rxn01464`  `rxn01720`  `rxn01860`  `rxn01969`  `rxn02150`  `rxn02435`  `rxn02455`  `rxn02977`  `rxn03090`  `rxn03273`  `rxn03325`  `rxn03367`  `rxn03479`  `rxn03676`  `rxn03833`  `rxn03877`  `rxn03970`  `rxn03971`  `rxn04361`  `rxn04539`  `rxn04542`  `rxn04734`  `rxn10878`  `rxn13012`  `rxn13774`  `rxn15001`  `rxn15010`  `rxn15166`  `rxn15494`  `rxn19345`  `rxn19905`

## Figures

Registered in `scripts/figures.tsv` under tag `direction_change`:

```bash
python3 scripts/regen_figures.py --tag direction_change
```

* `direction_change_funnel` — all 76 → present in ≥1 model (36) → active in
  ≥1 pFBA solution (30), each stage split into the 56 explicit changes and the
  20 labelled "no change".
* `direction_change_by_reaction` — one row per reaction, three facets, tall.
  Nested bars: palest = present, mid = present *and the model grows*, solid =
  *and it carries flux*.
* `direction_change_by_reaction_horizontal` — the same data wide, with the
  reactions along the x-axis and the three panels stacked as rows.
* `direction_change_by_reaction_horizontal_core`, `..._gsm`, `..._gsm_auxo` —
  one wide figure per model set, showing only the reactions present in it
  (5, 36 and 35 respectively).

Each writes a `.tsv` of the numbers it asserts beside the image.

## Reproduce

```bash
python3 scripts/direction_change/presence_and_activity.py -j 32 --dump-records
python3 scripts/direction_change/write_tsv.py
python3 scripts/direction_change/summarise.py
```

Roughly 5 minutes for all 16,523 models at `-j 32`. Deterministic: a second
run reproduced every count.

## Validation

- The `cma` seed-key mapping (model reaction id → MSDB `rxn` id) was checked
  against a raw-JSON scan of the same 40 core models and agreed exactly.
- Counts agree with the precomputed reaction universes in
  `modelseed2_gs_models/`: 5/76 against `reaction_universe_core.txt` (233
  reactions), 36/76 against `reaction_universe.txt` (3,769).
- Reaction orientation was checked against MSDB's canonical stoichiometry in
  150 core and 60 GSM models: **0 flipped**. This matters because the `= -> >`
  labels are written in MSDB's frame.
- 0 solver errors across all 16,523 models.

## Caveat: "active" is an upper bound on impact

Because all 76 go reversible → irreversible, a reaction is only genuinely
constrained if it carries flux **in the direction being forbidden**. A reaction
already running the permitted way is unaffected. The orientation check above is
the prerequisite for that calculation, but it has not been run — so treat
30/76 as the set of candidates, not the set of reactions that will change flux.
