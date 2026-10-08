# Enzyme functionality behind the 1,370 strict direction disagreements

Built 2026-10-08 against ModelSEED `dev` @ `078a395f`
(`/scratch/ctaylor/tmp/devsnap_078a395f`), the same snapshot as
[`reports/llmVsRecommended/LLM_VS_RECOMMENDED_DIRECTION.md`](../llmVsRecommended/LLM_VS_RECOMMENDED_DIRECTION.md).

**Status: complete — 1,370 / 1,370 reactions researched, no gaps.**

## What was asked, and what came back

The LLM-vs-recommended report splits 8,466 direction mismatches into six
groups. Two of them are *strict* disagreements, where the recommendation and
the LLM council pick opposite sides rather than merely differing over whether
a reaction is bidirectional:

| group | recommended | LLM ensemble | reactions |
|---|:---:|:---:|---:|
| `forward_vs_reverse` | `>` forward | `<` reverse | 62 |
| `reverse_vs_forward` | `<` reverse | `>` forward | 1,308 |
| **total** | | | **1,370** |

The question put to the literature was: **is each of these reactions run by a
bifunctional or polyfunctional enzyme?** — the hypothesis being that an enzyme
with two fused catalytic activities might run chemistry in both senses and so
explain why two direction calls disagree.

**The answer is no, and the margin is not close.**

| | reactions | share |
|---|---:|---:|
| monofunctional | 790 | 57.7% |
| **bifunctional** | **73** | **5.3%** |
| **polyfunctional** | **26** | **1.9%** |
| unknown (no enzyme identifiable) | 481 | 35.1% |
| **bi- or polyfunctional, combined** | **99** | **7.2%** |

![Horizontal bar chart on a log scale of enzyme functionality across the 1,370 strict direction disagreements: monofunctional 790 (57.7%), unknown 481 (35.1%), bifunctional 73 (5.3%), polyfunctional 26 (1.9%).](figures/strict_disagreement_functionality.png)

Only 50 of the 99 multifunctional calls are high-confidence, and only 35 carry
a note saying the multifunctionality bears on the *direction* at all. So the
honest headline is that **multifunctional enzymes explain a few dozen of 1,370
disagreements, not the phenomenon.**

The recurring genuine cases are worth having: lactase-phlorizin hydrolase (6
reactions, two distinct glycosidase sites on one chain), DNA polymerase I (4,
polymerase + two exonucleases), aphidicolan-16β-ol synthase (3, fused class
II/class I diterpene cyclase), class A PBPs (3, transglycosylase +
transpeptidase). The single cleanest case where bifunctionality really does
drive the conflict is **`rxn55967`**, a 1-testosterone hydratase/dehydrogenase
whose dehydrogenase half pulls the hydration forward, and **`rxn50031`**
(PaaZ), whose fused ALDH domain irreversibly consumes the hydratase domain's
product.

## What actually explains the disagreements

Having searched all 1,370 reaction-by-reaction, the dominant finding is a
different one, and it replicated in every chunk independently.

### 1. The enzymology almost always sides with the LLM council — but only in one group

| group | n | LLM right | ModelSEED right | genuinely reversible | unclear |
|---|---:|---:|---:|---:|---:|
| `forward_vs_reverse` | 62 | 31 (50%) | **12 (19%)** | 13 (21%) | 6 (10%) |
| `reverse_vs_forward` | 1,308 | **1,221 (93%)** | 17 (1%) | 41 (3%) | 29 (2%) |

**These are two different phenomena and must not be pooled.** The small group
is a genuine mixed bag where ModelSEED is often right (Renilla luciferin
sulfotransferase really does desulfate; succinate:quinone oxidoreductase
really is named misleadingly). The large group is a near-systematic failure.

### 2. More than half of these equations are broken, and ModelSEED already knows

`equation_defect`, assigned per reaction from the literature:

| | reactions |
|---|---:|
| equation sound | 504 (36.8%) |
| mass unbalanced | 430 |
| cofactors stripped | 251 |
| lumped multi-enzyme | 107 |
| source disowns the entry | 78 |
| **carrying some defect** | **866 (63.2%)** |

This is corroborated by a check that needs no literature and no model
judgment — ModelSEED's own `status` field:

| set | flagged mass/charge imbalanced |
|---|---:|
| all non-obsolete reactions (baseline) | 8,125 / 48,384 = **16.8%** |
| all 8,466 recommended-vs-LLM mismatches | 1,128 = **13.3%** |
| **the 1,370 strict disagreements** | **710 = 51.8%** |

**3.1× the database baseline.** Note that general disagreement is *not*
enriched (13.3%, below baseline) — broken stoichiometry travels specifically
with strict, opposite-direction conflicts.

![Two panels. Left: 100% stacked bars of which side the enzymology supports, per mismatch group. Right: bar chart of equation-defect categories across the researched reactions, titled with the imbalance enrichment.](figures/strict_disagreement_direction.png)

### 3. ModelSEED's recommendation contradicts its own primary source

The snapshot ships a MetaCyc dump under `Provenance/MetaCyc/`, whose equations
carry the arrow MetaCyc itself asserts. For the 709 strict reactions that
resolve there:

| | count | share |
|---|---:|---:|
| MetaCyc writes `=>` (irreversible forward), ModelSEED recommends `<` | **491** | **69.3%** |
| MetaCyc's direction agrees with the **LLM council** | 502 | **70.8%** |
| MetaCyc's direction agrees with the **recommendation** | 45 | **6.3%** |

This is the strongest result here because no LLM judgment enters it. Where the
source database has an opinion, the thermodynamic cascade overrides it and is
wrong roughly eleven times out of twelve.

Separately, of the 190 reactions carrying a stored `direction` field, **184
(97%) say `=` (reversible)** and were overridden to a strict direction.

### 4. A confident grade on a broken equation

| thermo grade | n | ModelSEED flags the equation imbalanced |
|---|---:|---:|
| **gold** | 49 | **36 (73%)** |
| silver | 422 | 188 (45%) |
| bronze | 899 | 486 (54%) |

Gold — the *highest*-confidence tier — has the *worst* imbalance rate. The
grade certifies the number's self-consistency, not the chemistry's validity.
`rxn14318` is the worked example: gold, self-certain, corroborated, on a
reaction KEGG itself annotates "unclear reaction" with no EC and no orthology.
This is the same conclusion `reports/thermoSourceMethod/` reached from a
different direction — *grade is a trust label, not a selector*.

### 5. Specific, mechanically fixable defect classes

- **Dropped products that already exist in the database.** 349 reactions
  (25%) have a mass residual matching an existing compound exactly — most
  often oxidised glutathione `cpd00111` (39 reactions, independently found by
  hand in five separate chunks), ADP `cpd00008` (16), CO₂ (13), PPi (9), O₂
  (9). These are import bugs with the counterpart sitting in the compound
  table, not research problems. *Caveat: formula matching cannot distinguish
  isomers, so entries like "D-Glucose" are candidates, not confirmations.*
- **Orphan ADP.** 17 reactions carry a residual of exactly one or three whole
  ADP molecules. **That is 17 of the 21 such cases in the entire database
  (81%)** — 28.6× enriched.
- **dGPredictor has the wrong sign on SAM methyl transfer.** Of 93
  SAM-containing reactions here, **92 (99%) get a positive ΔG**, median
  **+32.8 kJ/mol**, for chemistry that is ~−60 kJ/mol. The values repeat
  verbatim (29.04 ×11, 32.81 ×10, 32.25 ×10), so a handful of mis-signed
  class-level group estimates are flipping dozens of reactions each.
- **Free FAD/FADH₂ standing in for ETF/quinone coupling.** The single most
  reported `cofactors_stripped` pattern, hitting acyl-CoA
  dehydrogenases/oxidases across at least a dozen chunks. Removing the
  electron sink is exactly what removes the irreversibility.
- **Combinatorial expansions.** 126 reactions (9.2%) are generic-reaction
  expansions (`.metaexp.`, `.chlamyexp.`, `.maizeexp.`) that pair
  non-corresponding substrate and product instances — e.g. a tricyclic C20
  diterpene aldehyde "oxidised" to a straight-chain C20 fatty acid. Only
  **1.4× enriched**, so expansions are not themselves a cause of strict
  disagreement; they are simply reliably broken when they do appear.

## The 481 `unknown` calls

These are honest, not unresearched. They break down as: reactions KEGG or
MetaCyc explicitly disowns ("unclear reaction", "incomplete reaction",
"spontaneous", "multi-step"), entries whose source record has been *withdrawn*
upstream (several confirmed 404 at `rest.kegg.jp`), genuinely non-enzymatic
chemistry (spontaneous lactonisations, autoxidations, thermal rearrangements,
the previtamin D₃ → vitamin D₃ sigmatropic shift), lumped multi-enzyme
segments with no single catalyst, and legacy activities measured once in
1950s–70s crude extracts and never cloned. By tier: 270 `xref_only`, 179
`named`, 32 `opaque`.

## Files

```
results/strict_disagreement_enzymes/
  strict_disagreement_enzymes.json   the deliverable -- all 1,370 reactions
  summary.json                       every count and cross-check in this report
  research_queue.json                the work definition (1,370 entries)
  findings/chunk_NN.json             raw per-chunk literature findings (55)
  findings/backfill_*.json           structured verdicts for chunks 00-10
  RESEARCH_PROTOCOL.md               the classification rules given to each agent
  functionality_counts.tsv           figure source data
```

Each record in the deliverable carries:

```json
{
  "modelseed_id": "rxn02017",
  "name": "11-epi-Prostaglandin F2alpha:NADP+ 11-oxidoreductase",
  "synonyms": [], "abbreviation": "...",
  "reaction": "(1) 11-epi-PGF2a[0] => (2) H+[0] + (1) Prostaglandin D2[0]",
  "reaction_equation": "(1) cpd03548[0] => (2) cpd00067[0] + (1) cpd00527[0]",
  "ec_numbers": ["1.1.1.188"], "xrefs": {}, "pathways": [],
  "is_transport": false,
  "modelseed_status": "CI:2", "modelseed_imbalance_flags": ["CI"],
  "enzyme_name": "prostaglandin F synthase (PGFS), an aldo-keto reductase",
  "gene_symbols": ["AKR1C3"],
  "functionality": "bifunctional",
  "distinct_activities": ["PGH2 9,11-endoperoxide reductase", "PGD2 11-ketoreductase"],
  "organism_context": "...", "confidence": "high",
  "rationale": "...", "evidence_pmids": ["16475787", "3456602"],
  "direction_conflict": {
    "mismatch_group": "forward_vs_reverse",
    "modelseed_recommended": "forward", "llm_ensemble": "reverse",
    "thermo_evidence": {...}, "thermodynamics": {...},
    "enzymology_supports": "llm", "equation_defect": "none", "note": "..."
  },
  "research_tier": "named"
}
```

## Reproducing

```bash
python3 scripts/build_strict_disagreement_queue.py      # define the 1,370 + chunk them
# (run the literature subagents against results/.../chunks/, writing findings/)
python3 scripts/build_strict_disagreement_outputs.py    # merge + every cross-check
python3 scripts/regen_figures.py strict_disagreement_enzymes
```

## Method

The unit of work is deliberately the **reaction**, not the EC number. A
separate EC-keyed pipeline exists (`results/enzyme_functionality/`) and is
**not** reused here: an EC class is shared by chemistry that different proteins
run, so each reaction was researched on its own to find the enzyme that
actually catalyses *it*.

The 1,370 were split into 55 chunks of 25 and handed to `general-purpose`
subagents on **Claude Opus 5 at `xhigh` effort**, each given
[`RESEARCH_PROTOCOL.md`](../../results/strict_disagreement_enzymes/RESEARCH_PROTOCOL.md)
and the PubMed/Europe PMC MCP tools. The protocol carries two exclusions that
are easy to get wrong and were applied throughout: **moonlighting is not
bifunctionality** (a regulatory, structural or nucleic-acid-binding second role
does not count — hexokinase's sugar signalling, aconitase/IRP1, GPX4's sperm
capsule), and **a multi-subunit complex is not a fused bifunctional enzyme**
(E1/E2/E3 complexes, PdxS/PdxT, StyA/StyB). Agents were instructed to use
`unknown` honestly rather than defaulting to monofunctional.

## Caveats

- **This is a literature-search signal, not a curated-database lookup.**
  BRENDA and UniProt were not available as tools. 401 of 1,370 calls are
  `low` confidence, concentrated in secondary metabolism where the chemistry
  is clear but no protein has been cloned.
- **Chunks 00–10 (275 reactions) predate the structured fields.** They were
  dispatched before `direction_verdict` and `equation_defect` were added to
  the protocol, then backfilled by a later classification pass over the text
  those agents had already written. That backfill is text-faithful but
  noisier than natively-collected labels — sibling reactions occasionally got
  different defect labels because only one of the pair mentioned the cause.
- **`equation_defect` labelling is not perfectly consistent across chunks.**
  Some agents reserved `source_disowned` for an explicit source-side flag and
  used `none` for verified-spontaneous reactions; others did the reverse. The
  *presence* of a defect is reliable; the exact bucket is approximate.
- **MetaCyc/BioCyc and Rhea were unreachable over the network** for most of
  the run. Later chunks used the local MetaCyc dump in the snapshot, which is
  better (citable and pinned); earlier chunks fell back on chemistry plus
  literature for MetaCyc-only entries, and those carry lower confidence.
- **Several agents cross-checked `status` against
  `/scratch/ctaylor/ModelSEEDDatabase/`, a working copy at a different commit
  than the pinned snapshot.** Status strings quoted inside `rationale` text may
  therefore differ slightly from the snapshot. Every number in this report is
  computed from the pinned snapshot.
- **The LLM council is not a gold standard.** "The enzymology supports the LLM
  call" means the council happened to agree with the physiological direction,
  usually because the reaction was written in its biosynthetic sense and the
  council read the chemistry rather than a ΔG. It is evidence against the
  current recommendation cascade, not evidence that an LLM ensemble should
  replace it.
