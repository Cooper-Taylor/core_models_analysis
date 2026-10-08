# Research protocol — enzyme functionality behind strict direction disagreements

You are classifying **reactions**, one at a time. The unit of work is the
reaction, NOT the EC number. Do not classify an EC class and stamp it onto a
reaction: find the enzyme that actually catalyses *this* reaction, then judge
that enzyme. Two reactions sharing an EC may be run by different proteins.

## The question for each reaction

1. **Which enzyme (protein) catalyses this reaction?** Give its common name
   and, where you find them, gene symbol(s) and organism context.
2. **Is that enzyme mono-, bi-, or polyfunctional?**

## Categories

- `monofunctional` — one catalytic activity.
- `bifunctional` — a single polypeptide (or a fused multidomain chain)
  carrying **two distinct catalytic activities**. Textbook: PFK-2/FBPase-2;
  methylenetetrahydrofolate dehydrogenase/cyclohydrolase; DHFR-TS.
- `polyfunctional` — **three or more** fused catalytic activities. e.g. CAD
  (CPSase/ATCase/DHOase); DNA polymerase I (polymerase + two exonucleases);
  the peroxisomal trifunctional enzyme.
- `unknown` — you could not identify the catalysing enzyme, or found no
  usable evidence either way. **Use this honestly.** Do not default to
  `monofunctional` to fill the field — a real "unknown" is more useful than a
  fabricated "mono".

### Two exclusions that are easy to get wrong

- **Moonlighting is NOT bifunctional.** A second role that is regulatory,
  structural, or nucleic-acid-binding does not count. Hexokinase's sugar
  *signaling* role, aconitase/IRP1's RNA binding, GPI/neuroleukin — all stay
  `monofunctional`. Only a second genuine **catalytic** activity counts.
- **A multi-subunit complex is NOT bifunctional.** Separate polypeptides that
  associate (2-oxoglutarate dehydrogenase E1/E2/E3; biphenyl dioxygenase's
  oxygenase/ferredoxin/reductase) are not a fused bifunctional enzyme. The
  activities must sit on one chain (or a genuinely fused gene product).
  If a complex is the honest answer, say so in `rationale` and use
  `monofunctional` for the subunit that runs this reaction.

Bifunctionality that exists only in some lineages (e.g. DHFR-TS is fused in
protozoa and plants, separate in humans) counts as `bifunctional` — record the
lineage in `organism_context`.

## How to search

Use the PubMed/EuropePMC MCP tools:

- `pubmed_europepmc_search` — primary (broadest corpus).
- `pubmed_search_articles` — fallback / to confirm.
- `pubmed_fetch_articles` — only when an abstract is genuinely needed.

Budget roughly **1–3 searches per reaction**. Searching the enzyme name with
`(bifunctional OR multifunctional OR trifunctional OR "dual function")` is
usually enough to settle it; a clean monofunctional enzyme can be called in
one search, or from your own knowledge if you are confident — say so.

**Identifying the enzyme when the reaction has no name** (`research_tier` is
`xref_only` or `opaque`): the `definition` field spells the chemistry out in
plain compound names — search on substrate→product (e.g. "L-cysteate NADH
sulfolyase"). The `xrefs` block gives KEGG `R` numbers and MetaCyc reaction
ids; you may `WebFetch` a KEGG or MetaCyc entry page **solely to turn an
identifier into an enzyme name**. All bi/polyfunctionality evidence must still
come from the literature tools.

## Output

Write a JSON **array**, one object per input reaction, same order, to the
output path you were given. Every input reaction must appear exactly once.

```json
{
  "id": "rxn00837",
  "enzyme_name": "GMP reductase",
  "gene_symbols": ["guaC"],
  "functionality": "monofunctional",
  "distinct_activities": [],
  "organism_context": "",
  "confidence": "high",
  "rationale": "<=2 sentences, concrete. Name the activities if bi/poly.",
  "evidence_pmids": ["12345678"],
  "direction_verdict": "llm",
  "equation_defect": "none",
  "direction_note": ""
}
```

Field notes:

- `distinct_activities` — for `bifunctional`/`polyfunctional`, the separate
  catalytic activities by name (and EC if you know it). Empty otherwise.
- `confidence` — `high` / `medium` / `low`. Be strict: `low` if you are
  inferring the enzyme identity from chemistry rather than confirming it.
- `evidence_pmids` — PMIDs (or `PPR…`/DOIs) actually consulted. Empty list is
  acceptable for a confidently-known monofunctional enzyme; say so in
  `rationale`.
- `direction_verdict` — **required.** Which side does the enzymology actually
  support? One of:
  - `"llm"` — the LLM council's call matches the enzyme's real behaviour.
  - `"modelseed"` — ModelSEED's recommended direction is the correct one.
  - `"neither"` — the reaction is genuinely reversible in vivo, or direction is
    set by metabolite pools rather than by the enzyme.
  - `"unclear"` — not enough evidence to say.
  Judge the *enzymology*, not the thermodynamic numbers. The input gives you
  `recommended_reversibility` (ModelSEED) and `llm_direction` (the council).
- `equation_defect` — **required.** Is the ModelSEED equation itself sound
  enough for a direction to be meaningful? One of:
  - `"none"` — the equation is a fair statement of the chemistry.
  - `"cofactors_stripped"` — a redox/energy partner that sets the
    irreversibility is missing or abstracted (e.g. free FAD/FADH2 standing in
    for ETF coupling; O2 or NADPH dropped).
  - `"mass_unbalanced"` — atoms or a substrate are missing.
  - `"lumped_multienzyme"` — several enzymatic steps collapsed into one entry,
    so no single enzyme catalyses it.
  - `"source_disowned"` — KEGG/MetaCyc itself flags the entry as incomplete,
    unclear, spontaneous, or the EC has been deleted by IUBMB.
  Use `"none"` unless you have a concrete reason. If more than one applies,
  pick the one that most undermines the direction call.
- `direction_note` — **only if relevant**, one sentence: does this enzyme's
  bi/polyfunctionality, or its physiological role, bear on which *direction*
  the reaction runs? These reactions were selected because ModelSEED's
  recommended direction and an LLM council's call point opposite ways
  (`recommended_reversibility` vs `llm_direction` in the input). A
  bifunctional enzyme that runs opposite reactions on its two domains is
  exactly the kind of thing that explains such a conflict. Leave `""` when
  there is nothing real to say — do not speculate.

Return ONLY a short confirmation in your final message: how many reactions you
wrote, the count per functionality category, and anything that went wrong.
Do not paste the JSON into your reply.
