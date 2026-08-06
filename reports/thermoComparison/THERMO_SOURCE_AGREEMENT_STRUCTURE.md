# Which reactions do Group Contribution, eQuilibrator and dGPredictor agree on?

Follow-up to `THERMO_SOURCE_FBA_PIPELINE.md`, which reported the whole-database
scatter correlations: GC–eQuilibrator r = 0.91, GC–dGPredictor r = 0.08,
eQuilibrator–dGPredictor r = 0.17. This document asks *which reactions* carry
the agreement and which carry the disagreement, and what distinguishes them
chemically and structurally.

**Status: the mis-mapped dGPredictor values identified below have been filtered
out, and every figure and table in `reports/` has been regenerated against the
filtered data.** The filter is a read-time mask
(`results/thermo_agreement/dgpredictor_kegg_mask.json`, 17,271 reactions) built
by `build_dgpredictor_kegg_mask.py`. Nothing in `ModelSEEDDatabase/` is
modified — the stored records are left in place and consumers skip them.

Working set after filtering: the **8,359 ModelSEED reactions carrying a
non-sentinel ΔG′° with a defined operator from all three sources** (14,578
before). Group Contribution and eQuilibrator are untouched by the mask, so the
GC-vs-eQuilibrator comparison is unchanged at n = 18,477.

To reproduce any pre-filter number, every entry point takes `--no-dgp-mask`.

## Summary

Two separate effects were found, and they are not the same size.

1. **A reaction-identity mapping defect, which is the larger effect.**
   dGPredictor predicts from a *KEGG* reaction. Where the KEGG id it was run on
   is one ModelSEED itself lists as an alias of that reaction, dGPredictor
   tracks Group Contribution at **r = 0.74**. Where it is not, **r = −0.005**.
   Both halves are ~6–8 k reactions. The headline "dGPredictor is uncorrelated"
   was mostly this, and those 17,271 values are now masked out.
2. **A real chemical split, which is what remains after filtering.**
   **Redox chemistry** is concordant (oxidoreductases r = 0.92, NAD(P)(H) pairs
   r = 0.93, CO₂-releasing reactions r = 0.97, median |GC − dGP| ≈ 4 kcal/mol)
   and **group transfer** is not (phosphoryl transfer r = 0.29,
   phosphoanhydride change r = 0.00, acyl/CoA r = 0.08, SAM/SAH methyl transfer
   r = 0.10, amide-N transfer r = −0.21). In the group-transfer families **GC
   and eQuilibrator disagree with each other too** (ATP r = 0.10, ligases
   r = 0.06), so this is hard chemistry rather than one bad method.

### What the filter changed

| quantity | before mask | after mask |
|---|---:|---:|
| reactions with all three sources | 14,578 | 8,359 |
| r(GC, dGPredictor), all reactions | 0.103 | **0.743** |
| r(eQuilibrator, dGPredictor) | 0.161 | **0.684** |
| r(GC, eQuilibrator) | 0.938 | 0.834 |
| median \|GC − dGP\| | 8.5 | 5.5 kcal/mol |
| dGPredictor variance inside a signature class | 22.5% | **0.74%** |
| dGPredictor out-of-sample additive residual | 9.21 | **1.58** kcal/mol |
| dGPredictor per-model reaction coverage | — | 59.0% |

The last two rows matter for the "what dGPredictor is not doing wrong" section:
both the apparent non-determinism and most of the apparent non-additivity were
the mapping defect. Post-filter, dGPredictor's within-signature-class variance
(0.74%) is indistinguishable from eQuilibrator's (0.81%).

## 1. The mapping defect

`Update_Reaction_dGPredictor_Energies.py` keys predictions
`ModelSEED-rxn → KEGG-R-id → {dG_mean, dG_uncer}` and stores the mean over all
KEGG ids for that reaction. Of the 27,971 reactions with a staged prediction,
**17,359 were given a KEGG id that ModelSEED does not list as an alias** of that
reaction. (Where ModelSEED *does* list one, the staged id always matches — there
are zero disjoint cases, so the defect is confined to the no-alias group.)

Those inferred ids are heavily reused: 2,439 KEGG ids are shared by more than one
ModelSEED reaction, covering 23,516 reactions. The worst offenders:

| KEGG id | ModelSEED reactions | what KEGG says the reaction is |
|---|---:|---|
| R09245 | 861 | hexaprenyl diphosphate synthase |
| R06601 | 737 | 5-hydroxyisourate amidohydrolase |
| R10126 | 545 | citronellol dehydrogenase |
| R03877 | 220 | Mg-protoporphyrin IX chelatase |
| R10412 | 207 | farnesol dehydrogenase |

### The unvouched ids are chemically unrelated, checked independently

"Vouched" is defined off ModelSEED's alias list, so on its own it could be
measuring bookkeeping rather than correctness. It is not. Taking each staged
KEGG id's own reaction definition from `Scripts/Release/archived/kegg_reactions.txt`
and asking what fraction of the ModelSEED reaction's participants appear by name
in it (ignoring H⁺, H₂O, CO₂, Pᵢ, PPᵢ, O₂, which appear everywhere):

| subset | n | mean participant-name overlap | share with **zero** overlap |
|---|---:|---:|---:|
| vouched | 9,048 | 0.698 | 5.6% |
| unvouched | 15,178 | **0.039** | **90.2%** |

90% of the unvouched reactions share **not one named participant** with the KEGG
reaction whose ΔG they were assigned. This uses no alias information at all.

Two other checks on the flag itself: the standalone
`Biochemistry/Aliases/Unique_ModelSEED_Reaction_Aliases.txt` and the reaction
JSON `aliases` field agree exactly (14,066 reactions, identical id sets both
ways), so the split is not an artifact of reading the wrong alias source; and
the ids do not come from `linked_reaction` (the linked reactions carry no KEGG
alias either).

A worked example. 768 ModelSEED reactions have a net chemistry of exactly
ATP + H₂O → ADP + Pᵢ. Group Contribution gives every one of them −6.81
kcal/mol; eQuilibrator gives 691 of them −6.16. dGPredictor gives at least
twelve different values spanning −75.6 to +16.8:

| ModelSEED reaction | staged KEGG id | KEGG reaction | dGP (kcal/mol) | count |
|---|---|---|---:|---:|
| rxn00062 ATP phosphohydrolase | R00086 | ATP + H₂O = ADP + Pᵢ | −6.40 | 43 |
| rxn18189 Mg²⁺-ATPase | R06601 | 5-hydroxyisourate amidohydrolase | −1.85 | 487 |
| rxn05145 phosphate-transporting ATPase | R07456 | pyridoxal 5′-phosphate synthase | −75.63 | 39 |
| rxn13712 glutathione export, ABC | R09245 | hexaprenyl diphosphate synthase | −67.33 | 40 |

Only the first row is a real mapping — and it is the only one where ModelSEED
lists a KEGG alias at all. The other three reactions have **no KEGG alias in
ModelSEED**, yet a KEGG id was staged for them anyway.

Splitting the 14,578-reaction set on this one flag:

(This table is measured on the *unfiltered* data — it is the evidence for the
mask, reproducible with `--no-dgp-mask`.)

| subset | n | r(GC,eQ) | r(GC,dGP) | r(eQ,dGP) | ρ(GC,dGP) | median \|GC−dGP\| |
|---|---:|---:|---:|---:|---:|---:|
| all | 14,578 | 0.938 | **0.103** | 0.161 | 0.409 | 8.5 |
| KEGG id is a ModelSEED alias | 8,359 | 0.834 | **0.743** | 0.684 | 0.626 | 5.5 |
| KEGG id inferred | 6,219 | 0.958 | **−0.005** | 0.002 | 0.147 | 15.0 |
| …vouched *and* id unique to one reaction | 3,231 | 0.844 | **0.872** | 0.755 | 0.683 | 6.5 |

This is not a difficulty confound. GC–eQuilibrator agreement goes *down* in the
vouched subset (0.958 → 0.834) while GC–dGPredictor goes up eight-fold, and a
ΔG-matched-null control (decile-stratified resampling on the family's own GC
distribution, `family_stats.tsv`) puts the vouched family at **+0.57 excess r**
over its matched null — larger than any chemical family in the scan.

Figure: `thermoComparison/figures/thermo_agreement/fig_dg_agreement_by_kegg_trust.png`.

## 2. The chemical split, after filtering

Family scan over the 8,359 surviving reactions (`family_stats.tsv`):

| concordant family | n | r(GC,eQ) | r(GC,dGP) | median \|GC−dGP\| |
|---|---:|---:|---:|---:|
| CO₂-releasing (decarboxylase/carboxylase) | 643 | 0.943 | **0.967** | 4.5 |
| NAD(P)/NAD(P)H redox pair | 2,269 | 0.961 | **0.933** | 4.2 |
| EC 1 oxidoreductase | 2,777 | 0.943 | **0.919** | 4.8 |
| aldehyde group change | 972 | 0.921 | **0.890** | 3.4 |
| ester group change | 337 | 0.916 | **0.907** | 5.2 |
| max carbon ≤ 6 | 1,461 | 0.788 | **0.845** | 3.4 |
| exactly one group type changes | 2,277 | 0.663 | **0.822** | 4.6 |

| discordant family | n | r(GC,eQ) | r(GC,dGP) | median \|GC−dGP\| |
|---|---:|---:|---:|---:|
| phosphoanhydride change | 881 | 0.774 | **0.004** | 10.6 |
| acyl transfer (thioester/CoA) | 544 | 0.842 | **0.082** | 3.7 |
| SAM/SAH methyl transfer | 264 | 0.600 | **0.105** | 18.7 |
| amide group change | 744 | 0.908 | **0.366** | 4.5 |
| amide-N transfer (Gln/Glu) | 279 | 0.862 | **−0.208** | 2.7 |
| EC 6 ligase | 257 | **0.060** | 0.214 | 11.4 |
| phosphoryl transfer (ATP/ADP/AMP) | 809 | **0.100** | 0.291 | 10.1 |
| EC 3 hydrolase | 1,299 | **0.226** | 0.261 | 5.9 |

The line is **bond-order change at a carbon/nitrogen skeleton** (agrees) versus
**transfer of a conserved activated group between carriers** (does not):
phosphoryl, methyl, acyl, glycosyl, amide-N. Note the second table's
`r(GC,eQ)` column — on ligases, ATP reactions and hydrolases the two *additive*
sources are no better with each other than either is with dGPredictor.

Figure: `thermoComparison/figures/thermo_agreement/fig_dg_agreement_by_reaction_class.png`.

### Why: ΔG of a group transfer is a small residue of large numbers

For each reaction, using ModelSEED's own per-compound `deltag` / `deltagerr`:

    turnover  = Σ |νᵢ · ΔGf,ᵢ|          κ   = |ΔG_rxn| / turnover
    σ_rxn     = √(Σ (νᵢ · errᵢ)²)       snr = |ΔG_rxn| / σ_rxn

κ depends on where the compound-energy scale puts its zero, so `snr` is the
convention-free version and both are reported. Agreement rises monotonically
with both. By snr decile, r(GC,dGP) runs 0.45 → 0.89 from decile 1 to 8; by κ
decile, 0.16 → 0.93. Phosphoryl-transfer reactions push a median 1,427 kcal/mol
of formation energy through the reaction to net a median |ΔG| of 14 (κ = 0.002);
oxidoreductases push 906 to net 12.5 (κ = 0.030, ~15× better conditioned).

Figure: `thermoComparison/figures/thermo_agreement/fig_agreement_vs_snr.png`.
Table: `results/thermo_agreement/cancellation_by_family.tsv`.

**SAM/SAH is the one family that breaks the pattern** — high snr (39.5) yet the
worst agreement of any family (median |GC−dGP| = 18.7 kcal/mol, and GC–eQ only
0.60). Group Contribution puts SAH 182 kcal/mol above SAM (`deltag` +115.75 vs
−66.63) with a stated error of only 0.96. Whatever the right number is, the
stored uncertainty on the sulfonium centre badly understates the real spread.

### What dGPredictor is *not* doing wrong

Two mechanisms were tested. Both were largely the mapping defect in disguise,
and both shrink dramatically once it is removed.

- **Non-additivity is much milder than it first appeared.** Fitting implied
  per-compound energies by least squares over the stoichiometry matrix
  (`additivity_fit.tsv`, held-out fifth). Post-filter, Group Contribution
  reproduces out-of-sample to a median 0.003 kcal/mol (92.3% within
  1 kcal/mol), eQuilibrator to 0.010 (87.8%), dGPredictor to **1.58 (42.2%
  within 1 kcal/mol, 72.3% within 5)**. Before the filter dGPredictor's residual
  was 9.21 kcal/mol with only 13.7% within 1. It is still not a sum over
  compounds — it is a reaction-level regressor and does not have to be — but it
  behaves far more additively than the raw database suggested.
- **There is no structural null space, and dGPredictor is not non-deterministic.**
  Grouping reactions by exact atom-centred signature difference (Morgan
  radius 1, the descriptor family dGPredictor is built on) gives 3,060 classes
  post-filter. Only 0.015% of Group-Contribution variance and 0.807% of
  eQuilibrator variance lies *within* a class. **dGPredictor is now 0.742%** —
  statistically the same as eQuilibrator. Before the filter it was **22.5%**,
  which looked like a model returning different answers for the same
  transformation; it was the mapping defect measured a different way.

## 3. HDBSCAN in the 3D (GC, eQ, dGP) space

`cluster_thermo_dg_3d.py` runs HDBSCAN (`min_cluster_size=30`, `min_samples=10`,
EOM) on the standardised three-source coordinates, in raw kcal/mol and in a
signed-log transform. Post-filter, raw: **47 clusters covering 3,029 reactions,
63.8% noise** (before the filter: 77 clusters, 60.9% noise). The high noise
fraction is itself a result — most of the database is a continuous bulk, and the
clusters are archetypes sitting on top of it.

Before filtering, the clustering recovered the mapping defect without being told
about it: four clusters had dGPredictor pinned to ≤3 distinct values across the
whole cluster — one KEGG prediction broadcast over many ModelSEED reactions —
and those clusters were 0% vouched. **Post-filter no such cluster remains**, and
the point cloud collapses onto the GC = eQ = dGP diagonal.

Selected clusters (`cluster_profiles_raw.tsv`):

| cluster | n | GC | eQ | dGP | median \|GC−dGP\| | character |
|---|---:|---:|---:|---:|---:|---|
| C13 | 77 | −24.0 | −22.9 | −22.3 | **1.6** | H₂O₂ / O₂ / NH₃ oxidase; 100% O₂, 99% EC 1 |
| C35 | 76 | −5.7 | −5.8 | −6.5 | **1.1** | CoA thioester; 76% CoA, 76% thioester change |
| C36 | 80 | −6.8 | −6.3 | −6.4 | **0.4** | phosphatase; 85% Pᵢ, 74% EC 3 |
| C38 | 101 | 3.8 | 4.4 | 4.3 | 0.8 | NAD(H) ketone/aromatic redox |
| C10 | 71 | −95.8 | −106.0 | −81.7 | 14.1 | O₂ + NADPH monooxygenase |
| C23 | 282 | −14.0 | −3.2 | −3.5 | **10.5** | ATP/ADP phosphoanhydride transfer; 93% EC 2 |
| C15 | 88 | −4.4 | −1.8 | −23.1 | **18.8** | SAM/SAH methyl transfer; 98% SAM, 99% S |
| C18 | 69 | −21.6 | −2.8 | −1.1 | **19.7** | AMP/PPᵢ ligase; 61% EC 6 |

The concordant clusters are oxidase, redox, phosphatase and CoA-thioester
chemistry; the discordant ones are phosphoryl, methyl and adenylyl transfer —
the same split the correlation scan found, recovered without supervision.

Figures: `fig_dg_3d_clusters.png`, `fig_dg_3d_cluster_panels.png`.

## Reaction sets

`results/thermo_agreement/sets/` (see `SUMMARY.md` there):

| file | n | contents |
|---|---:|---|
| `set_A_three_way_concordant.tsv` | 2,722 | all three within 5 kcal/mol, vouched mapping, trivial isomerase zeros excluded |
| `set_B_group_transfer_discordant.tsv` | 385 | group-transfer chemistry, vouched, sources spread > 20 kcal/mol |
| `set_B2_…_ordinary_scale.tsv` | 352 | as B, restricted to \|ΔG\| ≤ 100 kcal/mol |
| `set_C_kegg_mismapped_withheld.tsv` | 17,271 | every reaction the mask withholds, with its staged KEGG id, that id's reuse count, and the withheld ΔG′° |
| `set_D_all_three_disagree.tsv` | 1,374 | GC and eQuilibrator themselves differ by > 15 kcal/mol |

Set C is the actionable one, and is the mask itself in human-readable form:
17,271 reactions whose stored dGPredictor ΔG′° is attributable to a KEGG
reaction that is not theirs. `results/thermo_agreement/dgpredictor_kegg_mask.tsv`
carries the same rows plus the per-reaction reason
(`unvouched_and_no_shared_participant` 14,424, `unvouched_and_unverifiable`
2,251, `unvouched` 596).

## Scripts

| script | what it does |
|---|---|
| `build_dgpredictor_kegg_mask.py` | **builds the mask.** Run first; everything else consumes it via `load_mask()`. `--lenient` rescues unvouched reactions whose participants match the KEGG definition |
| `build_thermo_agreement_features.py` | per-reaction feature table (146 cols): stoichiometry, elements, masses, cofactors, RDKit functional-group counts and deltas, EC, KEGG provenance. Needs the `eq3` env (rdkit). `--no-dgp-mask` writes the unfiltered table used by the mask-evidence figure |
| `analyze_thermo_agreement_families.py` | family scan with ΔG-matched null (`--trusted-only` is now a no-op — the input is pre-filtered) |
| `analyze_thermo_signature_nullspace.py` | additivity fit + signature-difference null space |
| `analyze_thermo_cancellation.py` | κ and snr mechanism test |
| `cluster_thermo_dg_3d.py` | HDBSCAN + cluster chemical profiling |
| `plot_thermo_agreement_findings.py`, `plot_thermo_dg_3d_clusters.py` | figures |
| `export_thermo_agreement_sets.py` | the reaction sets above |

## Caveats

- Pearson r over any reaction subset is strongly leverage-driven in this data.
  Every claim above is reported alongside Spearman ρ and median |difference|,
  and family correlations against a ΔG-matched null. Where only Pearson supports
  a claim it is not made — the apparent r = −0.70 for transport reactions, for
  instance, is 31 leverage points; inside |ΔG| < 100 it is +0.07.
- "Vouched" means the staged KEGG id is in ModelSEED's alias list. That is a
  necessary condition, not a proof the mapping is right, and absence of an alias
  is not proof it is wrong — some no-alias reactions carry legacy 1:1 KEGG
  correspondences from the original ModelSEED build (rxn00019 → R00024 and
  similar). The alias-independent name-overlap check above bounds how often that
  happens: ~10% of unvouched reactions do share participants with their assigned
  KEGG reaction. So the split is a strong statistical instrument and a good
  triage filter, but not a per-reaction verdict.
- The staged predictions themselves
  (`Biochemistry/Thermodynamics/dGPredictor/json_files/`) entered the repository
  as an opaque data upload (3870679, Vikas Upadhyay, Oct 2023). No script in the
  repository generates them, so how the KEGG ids were chosen cannot be
  reconstructed from this checkout — the evidence here is about what the mapping
  *is*, not how it was produced.
- The signature-difference classes are computed on the *ModelSEED* equation at
  Morgan radius 1, not on dGPredictor's own KEGG-side decomposition, so the
  22.5% within-class figure mixes true non-determinism with equation
  differences between the two databases.
