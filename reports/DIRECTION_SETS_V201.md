# Reaction directions in ModelSEED Biochemistry v2.0.1

Eleven direction sets, extracted from the v2.0.1 release (tagged 2026-09-15, the 2026 update) and registered so any of them can be
passed to `beginPipeline --directions <name>`.

Rebuild with `python3 scripts/build_direction_sets_v201.py`.

## What the release contains

v2.0.1 records direction two independent ways, and the paper is explicit that they are separate routes.

**Thermodynamic grading.** Each reaction is graded `gold`, `silver` or `bronze` on its collective evidence, anchored to openTECR measurements. The grade lives in a new `thermo-evidence` field beside `thermodynamics`. The recommended direction goes into `reversibility`, chosen by a fixed precedence of eQuilibrator > dGPredictor > Group contribution.

**An LLM council.** Three models predict, a fourth audits, a fifth adjudicates. It is stored as a fourth entry in `thermodynamics` under the label `LLMs`, carries no energy, and is deliberately kept out of the grading.

Two consequences shape everything below. A grade belongs to a **reaction**, not to a source, so there is no such thing as "Group contribution's gold calls" in this scheme. And because `reversibility` is populated only for graded reactions, the three grade tiers partition exactly the reactions that have a recommendation.

## The eleven sets

| set | family | calls | `>` | `<` | `=` | dropped `?` |
|---|---|---:|---:|---:|---:|---:|
| `v201_group_contribution` | source | 19246 | 10514 | 1423 | 7309 | 25625 |
| `v201_dgpredictor` | source | 12177 | 6863 | 2016 | 3298 | 17440 |
| `v201_equilibrator` | source | 21218 | 9129 | 1065 | 11024 | 3957 |
| `v201_llm_council` | source | 44850 | 38945 | 2350 | 3555 | 794 |
| `claude_opus48` | source | 51515 | 41909 | 1428 | 8178 | 0 |
| `v201_gold_only` | grade-exclusive | 3089 | 1687 | 96 | 1306 | 345 |
| `v201_silver_only` | grade-exclusive | 16150 | 6065 | 922 | 9163 | 2238 |
| `v201_bronze_only` | grade-exclusive | 8875 | 4179 | 1300 | 3396 | 2402 |
| `v201_gold` | grade-cumulative | 3089 | 1687 | 96 | 1306 | 345 |
| `v201_gold_silver` | grade-cumulative | 19239 | 7752 | 1018 | 10469 | 2583 |
| `v201_gold_silver_bronze` | grade-cumulative | 28114 | 11931 | 2318 | 13865 | 4985 |

No set contains `?`. A direction map is applied by rewriting bounds, and `?` is applied as fully reversible, so emitting it for a reaction that was graded but could not be called would *remove* a constraint in the name of evidence. Absent means no opinion, and the model keeps its own bounds. The last column counts what that dropped.

`gold + silver + bronze` (28,114 calls) reproduces the release's own `reversibility` field exactly wherever it is callable, which is the check that the extraction is faithful.

## Coverage is the first-order difference

The two LLM routes call two to four times as many reactions as any thermodynamic source:

- **Claude Opus 4.8**: 51,515
- **LLM council**: 44,850
- **eQuilibrator**: 21,218
- **Group contribution**: 19,246
- **dGPredictor**: 12,177

Both LLM sets are heavily forward-biased: the council calls 87% of its reactions forward, Claude 81%, against 43-53% for the thermodynamic sources. That is the signature to be suspicious of.

## Agreement

Percent agreement over the reactions each pair shares; `.` means no overlap.

| | Group con | dGPredict | eQuilibra | LLM counc | Claude Op | gold only | silver on | bronze on | gold | gold+silv | gold+silv |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **Group contribution** | -- | 76 | 78 | 74 | 78 | 90 | 79 | 90 | 90 | 81 | 84 |
| **dGPredictor** | 76 | -- | 73 | 63 | 65 | 51 | 90 | 90 | 51 | 76 | 81 |
| **eQuilibrator** | 78 | 73 | -- | 56 | 60 | 97 | 100 | 100 | 97 | 99 | 100 |
| **LLM council** | 74 | 63 | 56 | -- | 93 | 79 | 51 | 70 | 79 | 55 | 59 |
| **Claude Opus 4.8** | 78 | 65 | 60 | 93 | -- | 79 | 56 | 76 | 79 | 60 | 64 |
| **gold only** | 90 | 51 | 97 | 79 | 79 | -- | . | . | 100 | 100 | 100 |
| **silver only** | 79 | 90 | 100 | 51 | 56 | . | -- | . | . | 100 | 100 |
| **bronze only** | 90 | 90 | 100 | 70 | 76 | . | . | -- | . | . | 100 |
| **gold** | 90 | 51 | 97 | 79 | 79 | 100 | . | . | -- | 100 | 100 |
| **gold + silver** | 81 | 76 | 99 | 55 | 60 | 100 | 100 | . | 100 | -- | 100 |
| **gold + silver + bronze** | 84 | 81 | 100 | 59 | 64 | 100 | 100 | 100 | 100 | 100 | -- |

Three things stand out.

**The grades are mostly eQuilibrator.** The full graded set agrees with eQuilibrator on 99.52% of the 20,868 reactions they share, and with silver and bronze at 100%. That follows from the recommendation precedence putting eQuilibrator first, and it means a grade tier is closer to a confidence filter on eQuilibrator than to an independent consensus.

**The LLMs and the thermodynamics disagree.** The council agrees with eQuilibrator on only 56.22% of 17,108 shared reactions, Claude on 59.88%. They are not measuring the same thing, which is the paper's own argument for keeping them out of the grading.

**But the LLMs agree with each other**, 93.43% over 42,599 reactions. The council's members include this Claude model, so that is partly shared lineage rather than independent corroboration.

One more worth noting: dGPredictor agrees with the gold tier on only 51.08% of the 2,596 reactions they share, barely better than chance among three options.

## Effect on growth

FBA over the 100-model core panel, complete media, each set applied as a bound override. `on_disk` is the unmodified baseline.

| direction set | growers | mean flux of growers | reactions overridden per model |
|---|---:|---:|---:|
| *(baseline)* | 83/100 | 41.418 | 0.0 |
| gold | 79/100 | 34.0139 | 32.4 |
| gold + silver | 83/100 | 44.2626 | 66.2 |
| gold + silver + bronze | 83/100 | 37.1992 | 88.8 |
| eQuilibrator | 96/100 | 52.8986 | 83.9 |
| LLM council | 80/100 | 39.7834 | 93.8 |
| Claude Opus 4.8 | 59/100 | 41.5059 | 116.4 |

The incremental tiers do not behave monotonically. Gold alone costs four models their growth, because it constrains 32 reactions per model on the strength of the best evidence and nothing else loosens them. Adding silver restores all four and raises mean flux above baseline. Adding bronze keeps the grower count but pulls mean flux back down. More evidence is not uniformly more permissive.

The two LLM sets diverge sharply despite agreeing with each other 93% of the time: the council leaves 80 models growing, Claude alone only 59. Claude overrides the most reactions of any set (116 per model), and the difference is what the council's audit-and-adjudicate step removes.

## Using them

```bash
# the incremental grade story
python3 beginPipeline --models core_kegg2 --panel core_100 \
    --directions v201_gold --directions v201_gold_silver \
    --directions v201_gold_silver_bronze

# every source, comparison only, no FBA
python3 beginPipeline --models core_kegg2 --only tag:directions \
    --directions v201_equilibrator --directions v201_llm_council \
    --directions claude_opus48

# on the genome-scale models instead
python3 beginPipeline --models ms2_gsm --limit 200 \
    --media /scratch/ctaylor/modelseed2_gs_models/gmm/media/Carbon-D-Glucose.cpd \
    --directions v201_gold_silver_bronze
```

`python3 beginPipeline --list` shows all of them alongside the older maps.

## Caveats

- **Claude Opus 4.8 is not part of the release.** It is a standalone single-model run held in this repo. The council in v2.0.1 includes this model among its members, so the two sets are related and their 93% agreement should not be read as independent replication.

- **The grade tiers are not independent of the sources.** Silver and bronze agree with eQuilibrator at 100% on shared reactions. Comparing a grade tier against eQuilibrator is close to comparing eQuilibrator against itself.

- **Growth numbers above are the core panel on complete media.** The genome-scale models touch 3,769 distinct reactions against the core panel's 239, so the same direction set constrains far more of the network there.

