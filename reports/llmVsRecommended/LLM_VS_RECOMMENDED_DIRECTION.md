# Recommended reversibility vs. the LLM ensemble's direction call

Built 2026-10-07 against ModelSEED `dev` @ `078a395f`
(`/scratch/ctaylor/tmp/devsnap_078a395f`).

## Background

The current ModelSEED release carries up to **four** thermodynamic sources per
reaction, under `thermodynamics` in `reaction_*.json`:

| source | what it is |
|---|---|
| `Group contribution` | legacy estimator, ΔG′° + error + direction |
| `dGPredictor` | ML-predicted ΔG′° + error + direction |
| `eQuilibrator` | component-contribution ΔG′° + error + direction |
| `LLMs` | an **ensemble ("council") of LLMs voting a direction** — e.g. `["", "", ">"]`. It deliberately carries no ΔG or error, only the direction, because the council reasons from the reaction rather than from a number. |

Separately, the top-level `reversibility` field on each reaction is the
**recommended direction**, computed by
`Scripts/Thermodynamics/Apply_Evidence_Grades_And_Recommendation.py` with a
fixed precedence:

```
eQuilibrator  >  dGPredictor  >  Group contribution
```

— the first of those three with a callable direction wins. **`LLMs` is never
in that precedence** (by explicit design — see the script's docstring and
`Scripts/Thermodynamics/Add_LLM_Direction_Calls.py`): the council's call ships
in the database as a fourth source, but cannot itself become the
recommendation. That makes "recommended disagrees with the LLM ensemble" a
real, first-class question rather than a tautology — this report answers it.

## Coverage

| | reactions |
|---|---:|
| total reactions (dev snapshot) | 56,006 |
| obsolete (excluded throughout) | 7,622 |
| non-obsolete | 48,384 |
| → recommended reversibility is `?` (excluded) | 24,342 |
| → no `LLMs` entry at all (excluded) | 2,758 |
| → `LLMs` entry present but abstained, `?` (excluded) | 794 |
| **→ both calls directional — compared below** | **20,490** |

Most of the non-obsolete database (24,342 reactions, over half) has no
recommended direction yet, so it cannot be compared and is left out. Of the
20,490 reactions where both the recommendation and the LLM ensemble made a
real call:

- **12,024 (58.7%) agree**
- **8,466 (41.3%) disagree**

## The six mismatch groups

Excluding the three diagonal (agreement) cells, every disagreement falls into
exactly one of six groups — recommended direction vs. LLM ensemble direction:

| group key (JSON) | recommended | LLM ensemble | count |
|---|:---:|:---:|---:|
| `reversible_vs_forward` | `=` reversible | `>` forward | **6,098** |
| `reversible_vs_reverse` | `=` reversible | `<` reverse | **692** |
| `forward_vs_reversible` | `>` forward | `=` reversible | **159** |
| `forward_vs_reverse` | `>` forward | `<` reverse | **62** |
| `reverse_vs_reversible` | `<` reverse | `=` reversible | **147** |
| `reverse_vs_forward` | `<` reverse | `>` forward | **1,308** |
| **total mismatches** | | | **8,466** |

![Recommended reversibility vs. LLM ensemble direction call — a 3×3 matrix. Diagonal cells (green) are agreement: forward/forward 9,760, reversible/reversible 1,909, reverse/reverse 355. The six off-diagonal cells (blue, darker = larger) are the mismatch groups: reversible-recommended but LLM-forward 6,098, reverse-recommended but LLM-forward 1,308, reversible-recommended but LLM-reverse 692, LLM-reversible but recommended-forward 159, LLM-reversible but recommended-reverse 147, and forward-recommended but LLM-reverse 62.](figures/llm_vs_recommended_confusion.png)

*(Figure source: `scripts/plot_llm_vs_recommended_confusion.py`, registered in
`scripts/figures.tsv` as `llm_vs_recommended_confusion`. Regenerate with
`python3 scripts/regen_figures.py llm_vs_recommended_confusion`. Numbers are
also in `figures/llm_vs_recommended_confusion_values.tsv`.)*

### Reading the matrix

By far the largest mismatch group is **`reversible_vs_forward`**: 6,098
reactions ModelSEED recommends as bidirectional (`=`) that the LLM council
calls strictly forward (`>`). That alone is 72% of all mismatches. The next
largest is **`reverse_vs_forward`** (1,308): the recommendation says reverse,
the council says forward — i.e., the two don't even agree on which side is
favored, not just on whether it's reversible. The remaining four groups
(forward↔reversible, forward↔reverse) are much smaller (62–159 each) and
together account for under 5% of mismatches.

The pattern is consistent with the ensemble LLM leaning toward calling a
definite direction where the numeric cascade recommends "reversible" —
plausibly because `=` is itself often a fallback label (e.g. a small |ΔG| near
the ±threshold, or a tie among sources) rather than strong bidirectional
evidence, while the LLM council reasons qualitatively about the chemistry and
rarely abstains into "reversible."

## Reaction-level detail

Every one of the 8,466 mismatched reactions, with full identifying and
thermodynamic detail, is in:

```
results/llm_vs_recommended/mismatches_by_group.json
```

keyed by the six group names above, e.g.:

```json
{
  "reversible_vs_forward": [
    {
      "id": "rxn00004",
      "name": "4-hydroxy-4-methyl-2-oxoglutarate pyruvate-lyase (pyruvate-forming)",
      "abbreviation": "R00008",
      "equation": "(1) cpd02570[0] <=> (2) cpd00020[0]",
      "definition": "(1) Parapyruvate[0] <=> (2) Pyruvate[0]",
      "ec_numbers": ["4.1.3.17"],
      "is_transport": false,
      "recommended_reversibility": "=",
      "llm_direction": ">",
      "thermo_evidence": {"assessment": "self-confident", "cross-source": "corroborated", "grade": "silver", "source": "eQ"},
      "thermodynamics": {
        "Group contribution": [3.17, 1.55, "="],
        "LLMs": ["", "", ">"],
        "...": "eQuilibrator / dGPredictor triples when present"
      }
    },
    "... 6,097 more"
  ],
  "reversible_vs_reverse": ["..."],
  "forward_vs_reversible": ["..."],
  "forward_vs_reverse": ["..."],
  "reverse_vs_reversible": ["..."],
  "reverse_vs_forward": ["..."]
}
```

`thermo_evidence` is the grade ModelSEED assigned the *recommended* call
(gold/silver/bronze, which source it came from, and whether other sources
corroborated or disputed it) — useful for judging which mismatches are more
likely a confident recommendation overruling a shaky LLM call vs. the reverse.

The matching summary counts (the matrix plus the coverage/exclusion numbers
above) are in:

```
results/llm_vs_recommended/summary_counts.json
```

Both files are regenerated together by:

```
python3 scripts/build_llm_vs_recommended_mismatches.py
```

## Caveats

- **Obsolete reactions excluded.** This moves the total mismatch count by
  exactly 1 (8,467 → 8,466) and is not a meaningful filter — it's just hygiene.
- **"Recommended" is itself a derived field**, not a measurement — it is
  whichever of eQuilibrator / dGPredictor / Group contribution answers first.
  A `reversible_vs_forward` mismatch does not mean two independent
  measurements disagree; it means the *recommendation cascade's winner*
  disagrees with the *separate* LLM ensemble.
- **The LLM ensemble is a direction-only call** with no confidence/vote-margin
  field shipped in the database, so this report cannot (yet) distinguish a
  near-unanimous council call from a narrow one.
