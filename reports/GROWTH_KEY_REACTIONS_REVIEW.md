# Reactions that control model growth — methods & literature

This note documents the analyses the explorer uses to find reactions that are
**key to keeping a model's growth high or low**, and grounds each in the
constraint-based-modeling literature. All analyses run on the 100-model
descriptive panel with the FBA objective = biomass growth, after rebinding every
reaction to the ModelSEED reversibility cascade baseline and applying the
`KBaseMedia` complete medium.

The site computes four complementary, per-reaction signals:

| Signal | What it answers | Where in the UI |
|---|---|---|
| Direction sensitivity | does *flipping a reaction's direction* move growth? | Panel Models → "Key reactions"; Analytics → "Reaction criticality" |
| Knockout essentiality | does *removing a reaction* collapse (or raise) growth? | Panel Models → "Growth control"; Analytics → "Essential reactions" |
| Flux at the optimum | how much flux does a reaction carry at max growth? | Panel Models → "Growth control" (flux-vs-essentiality scatter) |
| Reduced cost | LP marginal: growth change per unit bound relaxation | Growth-control tooltips |

---

## 1. Knockout essentiality (single-reaction deletion)

For each reaction we set its bounds to `(0, 0)`, re-optimize biomass, and record
`Δgrowth = growth_KO − growth_baseline`. A reaction is **essential** when the
knockout collapses growth (`Δ < 0`) and **growth-limiting** when the knockout
*raises* growth (`Δ > 0` — the reaction was siphoning flux away from biomass).
This is the canonical *in silico* reaction/gene-deletion analysis introduced for
genome-scale models by **Edwards & Palsson (2000)** and standardized in the FBA
methodology of **Orth, Thiele & Palsson (2010)**.

Caveat — FBA assumes the deletion mutant re-optimizes to a new biomass optimum.
For real knockouts that is often violated transiently; **Segrè, Vitkup & Church
(2002, MOMA)** argue a minimal-perturbation objective predicts mutant fluxes
better. We report the FBA-optimum knockout (fast, parameter-free, and the right
notion for "can the network still grow at all"), not MOMA.

## 2. Direction sensitivity (reversibility, not deletion)

Distinct from deletion: for each reaction we force each *alternative direction*
in `{>, <, =}` (vs the cascade baseline) and re-optimize. A reaction is "key" by
direction when some assignment causes a large `|Δgrowth|`. This is the
growth-relevant projection of the thermodynamic directionality question that the
direction variants address — directionality errors are a known driver of
infeasible or mis-predicted flux (**Henry, Broadbelt & Hatzimanikatis 2007**,
thermodynamics-based flux analysis; **Noor et al. 2014**, pathway thermodynamics
and the max-min driving force, which links a reaction's thermodynamic driving
force to whether it can carry the flux growth needs).

## 3. Flux at the optimum & the flux-vs-essentiality map

From the baseline FBA solution we read the flux each reaction carries at the
growth optimum. Plotting `|flux|` against knockout `Δgrowth` separates:
high-flux + essential = **backbone**; low-flux + essential = **bottleneck**;
high-flux + dispensable = **redundant** (alternate optima carry the load);
on the zero line = **dispensable**. Flux degeneracy here is exactly the
alternate-optima phenomenon characterized by **Mahadevan & Schilling (2003)**;
**flux variability analysis (FVA)** and the **flux-coupling analysis** of
**Burgard et al. (2004)** are the principled tools for it and are natural next
additions.

## 4. Reduced cost (LP sensitivity)

The reduced cost is the linear-program dual on a reaction's bound: the marginal
change in growth per unit relaxation of that bound at the optimum. Together with
metabolite **shadow prices** it is the sensitivity analysis of growth to network
perturbations formalized for metabolism by **Reznik, Mehta & Segrè (2013)**
(flux-imbalance analysis). We surface per-reaction reduced costs in the
growth-control tooltips; metabolite shadow prices are a planned addition.

## 5. Cross-panel aggregation

Each per-reaction signal is tallied across all 100 panel models to find
reactions that are *broadly* controlling — e.g. essential in the most models, or
key-by-direction in the most models. Broadly essential reactions are candidate
universal growth controllers; reactions essential in only a few models point to
clade-specific metabolic dependencies.

---

## References

- Edwards JS, Palsson BO (2000). *The Escherichia coli MG1655 in silico metabolic genotype: its definition, characteristics, and capabilities.* PNAS 97(10):5528–5533.
- Orth JD, Thiele I, Palsson BO (2010). *What is flux balance analysis?* Nat Biotechnol 28(3):245–248.
- Segrè D, Vitkup D, Church GM (2002). *Analysis of optimality in natural and perturbed metabolic networks (MOMA).* PNAS 99(23):15112–15117.
- Mahadevan R, Schilling CH (2003). *The effects of alternate optimal solutions in constraint-based genome-scale metabolic models (FVA).* Metab Eng 5(4):264–276.
- Burgard AP, Nikolaev EV, Schilling CH, Maranas CD (2004). *Flux coupling analysis of genome-scale metabolic network reconstructions.* Genome Res 14(2):301–312.
- Reznik E, Mehta P, Segrè D (2013). *Flux imbalance analysis and the sensitivity of cellular growth to changes in metabolite pools.* PLoS Comput Biol 9(8):e1003195.
- Henry CS, Broadbelt LJ, Hatzimanikatis V (2007). *Thermodynamics-based metabolic flux analysis.* Biophys J 92(5):1792–1805.
- Noor E, Bar-Even A, Flamholz A, Reznik E, Liebermeister W, Milo R (2014). *Pathway thermodynamics highlights kinetic obstacles in central metabolism (MDF).* PLoS Comput Biol 10(2):e1003483.
- Price ND, Reed JL, Palsson BO (2004). *Genome-scale models of microbial cells: evaluating the consequences of constraints.* Nat Rev Microbiol 2(11):886–897.
