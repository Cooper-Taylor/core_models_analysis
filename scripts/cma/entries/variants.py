"""Every reversibility variant, cascade and overlay alike.

Fourteen ``cascade`` variants are knob-sets over ``ReversibilityConfig``; their
``cfg`` callables are imported from :mod:`reversibility_lib` unchanged, so the
numbers they produce are the numbers they produced before.

Six ``overlay`` variants carry a precomputed ``{rxn: operator}`` map instead.
Those six could not be expressed at all under the old catalog, whose ``cfg``
field was a mandatory ``ReversibilityConfig`` factory. They lived only in
``thermo_variants/manifest.json`` -- and ``export_thermo_variants.py`` rebuilt
that file from the catalog on every run, so any full export deleted them.
Registering them here is what makes them survive.

``order`` replaces four separate hand-maintained ordering lists.
"""

from __future__ import annotations

import reversibility_lib as lib

from ..kinds import Variant
from ..variants import register


# ── cascade variants ────────────────────────────────────────────────────────
def _baseline_cfg():
    return lib.ReversibilityConfig()

def _v31_cfg():
    return lib.ReversibilityConfig(ln_ri_by_rxn=lib.load_ln_reversibility_index())

def _v33_cfg():
    return lib.ReversibilityConfig(per_met_conc_range=lib.BENNETT_2009_ECOLI,
                             per_met_conc=lib.BENNETT_2009_MEAN)

def _v33w_cfg():
    return lib.ReversibilityConfig(cell_min=1e-7, cell_max=1e-1)

def _v35_cfg():
    return lib.ReversibilityConfig(sigma_band_k=1.96)

def _v35w_cfg():
    return lib.ReversibilityConfig(sigma_bounds_k=1.96)

def _v36_cfg():
    return lib.ReversibilityConfig(low_energy_cpds=())

def _v37_cfg():
    return lib.ReversibilityConfig(apply_special_conc=False)

def _v310t_cfg():
    return lib.ReversibilityConfig(mm_band=1.0)

def _v310l_cfg():
    return lib.ReversibilityConfig(mm_band=4.0)

def _H4_cfg():
    return lib.ReversibilityConfig(ln_ri_by_rxn=lib.load_ln_reversibility_index(),
                             sigma_band_k=1.96,
                             per_met_conc_range=lib.BENNETT_2009_ECOLI,
                             per_met_conc=lib.BENNETT_2009_MEAN)

def _ri_gamma1_cfg():
    return lib.ReversibilityConfig(ln_ri_by_rxn=lib.load_ln_reversibility_index(),
                             ln_ri_threshold=2.302585)

def _ri_gamma2_cfg():
    return lib.ReversibilityConfig(ln_ri_by_rxn=lib.load_ln_reversibility_index(),
                             ln_ri_threshold=4.605170)

def _dgpredictor_cfg():
    return lib.ReversibilityConfig(energy_override_by_rxn=lib.load_dgpredictor_energies(),
                             energy_override_label="dGPredictor")


register(Variant(
    tag='baseline',
    kind='cascade',
    cfg=_baseline_cfg,
    title='ReversibilityConfig() default (matches MSDB)',
    apt_title='Default cascade — reproduces the fixed MSDB byte-for-byte (reference)',
    description="Reproduces the upstream ModelSEEDDatabase Estimate_Reaction_Reversibility.py cascade exactly: heuristic 1 (bounded ΔG′° check across the concentration window) → ATP-synthase → ABC-transporter → heuristic 4 (mMdeltaG ±2 kcal/mol band) → heuristic 5 (LOW_ENERGY_CPDS points rule) → '=' (reversible) default. This baseline now includes the two shadow-bug repairs H2 and H3 (adopted upstream): the phosphate accumulator, the ABC-transporter rule, the phosphate-spread term, and the CO₂/O₂/H₂ concentration overrides all fire as intended. Every other variant is a one-knob (or few-knob) diff against this reference.",
    citations=('Henry 2007',),
    section='(reference)',
    order=10,
    in_notebook=True,
))

register(Variant(
    tag='3.1',
    kind='cascade',
    cfg=_v31_cfg,
    title='Persist + use ln(reversibility_index) (Noor 2012)',
    apt_title="Use eQuilibrator's reversibility index instead of the points-and-bands heuristics",
    description="Replaces the heuristic reversibility decision with eQuilibrator's reversibility index γ, calling a reaction directional whenever |ln γ| > 6.9 — i.e. when more than a ~1000-fold concentration shift would be needed to reverse it (cfg knob: ln_ri_by_rxn). γ is the molecularity-normalized fold-change in the mass-action ratio that flips the sign of ΔG′, so it is the formal answer to the question the ad-hoc points-and-bands rules approximate.",
    citations=('Noor 2012',),
    section='§ 2.1 / 3.1',
    order=20,
    in_notebook=True,
))

register(Variant(
    tag='3.3',
    kind='cascade',
    cfg=_v33_cfg,
    title='Bennett-2009 per-metabolite concentration ranges',
    apt_title='Use measured E. coli metabolite concentrations (Bennett 2009)',
    description="Replaces the uniform [10 µM, 20 mM] concentration prior with per-metabolite ranges and geometric-mean reference points from Bennett 2009's absolute E. coli metabolomics (cfg knobs: per_met_conc_range, per_met_conc), so ATP/glutamate sit at their measured ~10/100 mM levels and rare nucleosides near 0.1 µM. The real intracellular distribution spans six orders of magnitude with 70% of metabolites below 1 mM, so the uniform prior systematically biases the reference ΔG used by the directional check.",
    citations=('Bennett 2009', 'Park 2016', 'Bar-Even 2011', 'Salvy 2019'),
    section='§ 3.3',
    order=30,
    in_notebook=True,
))

register(Variant(
    tag='3.3_wide',
    kind='cascade',
    cfg=_v33w_cfg,
    title='Wider uniform conc window [1e-7, 0.1] M',
    apt_title='Widen the uniform concentration prior to [0.1 µM, 0.1 M]',
    description='Widens the baseline three-order [10 µM, 20 mM] uniform reagent concentration window to a six-order [0.1 µM, 0.1 M] window (cfg knobs: cell_min, cell_max), covering the observed spread of E. coli metabolomics without committing to per-metabolite measurements. This serves as a measured-data-free fallback that better envelopes Bennett 2009 observations for reactions whose intermediates are unmeasured.',
    citations=('Bennett 2009', 'Park 2016'),
    section='§ 3.3 (fallback)',
    order=40,
    in_notebook=True,
))

register(Variant(
    tag='3.5',
    kind='cascade',
    cfg=_v35_cfg,
    title='Per-reaction sigma band: k=1.96 (95%) replaces ±2 kcal',
    apt_title='Use a 95% Component Contribution (CC) uncertainty band instead of fixed ±2 kcal/mol',
    description="Replaces the fixed ±2 kcal/mol reversible-band check (inherited from Henry 2007) with k·σ_rxn where k=1.96 — a 95% confidence interval drawn from the Component Contribution method's per-reaction ΔG uncertainty already on disk (cfg knob: sigma_band_k). The 2 kcal/mol constant conflates ΔG-estimation error with concentration uncertainty, while CC's per-reaction σ tightens the band for well-measured reactions and widens it for poorly-constrained ones.",
    citations=('Noor 2013', 'Gollub 2021'),
    section='§ 3.5',
    order=50,
    in_notebook=True,
))

register(Variant(
    tag='3.5_wide',
    kind='cascade',
    cfg=_v35w_cfg,
    title='Per-reaction CC bound widening: k=1.96 on stored_bounds',
    apt_title='Widen the bounded-ΔG check to a 95% CC interval (1.96σ) instead of 1σ',
    description="The directional check requires the ΔG window — computed at concentration extremes — to lie entirely on one side of zero. Baseline uses 1-σ error bars on the stored Component Contribution bounds; this variant scales them by 1.96 for a 95% confidence interval (cfg knob: sigma_bounds_k). Noor 2013's empirical coverage validation (73%/90%/95%/99% at matching CC intervals) shows that a 1-σ bound calls reactions directional more aggressively than the data support.",
    citations=('Noor 2013', 'Gollub 2021'),
    section='§ 3.5 / § 2.5',
    order=60,
    in_notebook=True,
))

register(Variant(
    tag='3.6',
    kind='cascade',
    cfg=_v36_cfg,
    title='Drop the low-energy-compounds list entirely',
    apt_title='Drop the hand-curated low-energy-compounds list (CO₂, Pi, PPi, CoA, …)',
    description="Disables the legacy points rule that forced directional calls when reactions consumed hand-curated low-energy sinks (CO₂, NH₃, ACP, Pi, PPi, CoA, dihydrolipoamide, HCO₃⁻), by emptying the list (cfg knob: low_energy_cpds=()). The frozen MFAToolkit-era compound list is superseded by principled measures like the reversibility index (§3.1) or multiTFA's P(forward). The residual phosphate-spread term in the same rule is left in place — with the H3 repair now in the baseline it is live, so this variant isolates the LOW_ENERGY_CPDS contribution specifically.",
    citations=('Noor 2012', 'Gollub 2021'),
    section='§ 3.6',
    order=70,
    in_notebook=True,
))

register(Variant(
    tag='3.7',
    kind='cascade',
    cfg=_v37_cfg,
    title='Drop the CO2 1e-4 hardcoded concentration override',
    apt_title='Disable the hardcoded 1e-4 M CO₂ concentration override',
    description='The baseline applies a 1e-4 M CO₂ override and a 1e-6 M O₂/H₂ override during the mMdeltaG walk (both now live, since the H2/H3 shadow-bug repairs are in the baseline). This variant turns the gate off (cfg knob: apply_special_conc=False) so CO₂/O₂/H₂ sit at the 1 mM default instead. eQuilibrator already models CO₂(aq)/HCO₃⁻/CO₃²⁻ speciation as a function of pH, so the override would double-count what the transform handles.',
    citations=('Alberty 2003', 'Beber 2022'),
    section='§ 3.7',
    order=80,
    in_notebook=True,
))

register(Variant(
    tag='3.10_tight',
    kind='cascade',
    cfg=_v310t_cfg,
    title='Tighten mMdeltaG band: ±1 kcal/mol',
    apt_title='Tighten the reversible band to ±1 kcal/mol',
    description="Halves the reversible-band check from |reference ΔG| ≤ 2 kcal/mol to ≤ 1 kcal/mol (cfg knob: mm_band=1.0), forcing more reactions into directional calls. The 2 kcal/mol threshold was inherited from Henry 2007 without per-reaction justification; this variant exposes the cascade's sensitivity to that hardcoded threshold by making it a configurable input.",
    citations=(),
    section='§ 3.10',
    order=90,
    in_notebook=True,
))

register(Variant(
    tag='3.10_loose',
    kind='cascade',
    cfg=_v310l_cfg,
    title='Loosen mMdeltaG band: ±4 kcal/mol',
    apt_title='Loosen the reversible band to ±4 kcal/mol',
    description="Doubles the reversible-band check from |reference ΔG| ≤ 2 kcal/mol to ≤ 4 kcal/mol (cfg knob: mm_band=4.0), pulling more reactions out of directional calls into the reversible bucket. As with 3.10_tight, the motivation is to expose the reversibility decision's sensitivity to the Henry 2007 threshold and quantify how many directional calls hinge on the 2 kcal/mol choice alone.",
    citations=(),
    section='§ 3.10',
    order=100,
    in_notebook=True,
))

register(Variant(
    tag='H4',
    kind='cascade',
    cfg=_H4_cfg,
    title='(NEW) best-evidence composite: 3.1 + 3.5 + Bennett',
    apt_title='Best-evidence stack: reversibility index + 95% CC band + Bennett concentrations',
    description="Combines three high-impact §3 changes into one 'best available evidence' configuration: eQuilibrator's ln(reversibility index), a 95% Component Contribution per-reaction σ band in place of the fixed ±2 kcal/mol, and Bennett 2009 per-metabolite E. coli concentrations (cfg knobs: ln_ri_by_rxn, sigma_band_k=1.96, per_met_conc_range, per_met_conc). The three interact non-trivially: the CC band pulls many calls into the reversible bucket, the reversibility index re-pushes a subset back to firm directions, and the measured concentrations shift the reference ΔG of those still in the band.",
    citations=('Noor 2012', 'Noor 2013', 'Gollub 2021', 'Bennett 2009', 'Park 2016', 'Bar-Even 2011', 'Salvy 2019'),
    section='§ H4',
    order=110,
    in_notebook=True,
))

register(Variant(
    tag='ri_gamma1',
    kind='cascade',
    cfg=_ri_gamma1_cfg,
    title='Reversibility index, |log10 gamma| = 1 (aggressive)',
    apt_title='eQuilibrator reversibility index at a 10-fold cutoff (Noor 2012/2013)',
    description="Same reversibility-index rule as 3.1 but with an aggressive cutoff: a reaction is called directional whenever |ln gamma| > ln(10), i.e. when only a ~10-fold concentration shift would reverse it (|log10 gamma| = 1, vs 3.1's |log10 gamma| = 3). gamma is the molecularity-normalized fold-change in the mass-action ratio that flips the sign of ΔG′ (Noor 2012). This pulls a large set of near-equilibrium reactions out of the reversible bucket into firm '>'/'<' calls (cfg: ln_ri_by_rxn, ln_ri_threshold=ln 10).",
    citations=('Noor 2012', 'Noor 2013'),
    section='§ New — reversibility-index sweep',
    order=120,
    in_notebook=True,
))

register(Variant(
    tag='ri_gamma2',
    kind='cascade',
    cfg=_ri_gamma2_cfg,
    title='Reversibility index, |log10 gamma| = 2 (mid)',
    apt_title='eQuilibrator reversibility index at a 100-fold cutoff (Noor 2012/2013)',
    description='The reversibility-index rule at |ln gamma| > ln(100) (|log10 gamma| = 2): a reaction is directional when a ~100-fold concentration shift would be needed to reverse it. Together with 3.1 (|log10 gamma| = 3) and ri_gamma1 (= 1) this spans the full Noor reversibility-index threshold sweep (cfg: ln_ri_by_rxn, ln_ri_threshold=ln 100).',
    citations=('Noor 2012', 'Noor 2013'),
    section='§ New — reversibility-index sweep',
    order=130,
    in_notebook=True,
))

register(Variant(
    tag='dgpredictor',
    kind='cascade',
    cfg=_dgpredictor_cfg,
    title='dGPredictor reaction energies (Wang 2021)',
    apt_title="Swap the ΔG′° estimator to dGPredictor's molecular-fingerprint model",
    description="Runs the full baseline cascade but replaces each reaction's stored group-contribution ΔG′° with the dGPredictor estimate (Wang, Upadhyay & Maranas 2021) for the ~27.7k reactions dGPredictor covers. dGPredictor learns standard reaction energies from automated atom-centered molecular fingerprints rather than a fixed set of manually curated functional groups, so it handles isomerases and unusual transferases the group-contribution method struggles with. This isolates the effect of the energy estimator: the heuristics are unchanged, only the ΔG′° feeding them differs (cfg: energy_override_by_rxn).",
    citations=('Wang 2021',),
    section='§ New — dGPredictor energies',
    order=140,
    in_notebook=True,
))


def _majority_vote(maps):
    """Per-reaction majority vote across contributing direction maps.

    Ties fall back to '=' (reversible), which is the cascade's own default and
    the conservative choice for FBA.
    """
    from collections import Counter

    keys = set().union(*(set(m) for m in maps)) if maps else set()
    out = {}
    for k in keys:
        votes = Counter(m[k] for m in maps if k in m)
        if not votes:
            continue
        top, n = votes.most_common(1)[0]
        out[k] = top if list(votes.values()).count(n) == 1 else "="
    return out


# ── overlay variants (precomputed direction maps) ───────────────────────────
# These six were produced by standalone builder scripts and existed only in
# thermo_variants/manifest.json. Each keeps the tag, title and description it
# already carries there, so the site and the reports do not change.

register(Variant(
    tag='ai_opus48',
    kind='overlay',
    direction_map='ai_opus48',
    title='AI-predicted directions (Claude Opus 4.8)',
    apt_title='Reaction directions predicted by Claude Opus 4.8 (via Argo), overlaid on the MSDB baseline',
    description="Reaction directions predicted by the Claude Opus 4.8 model (through the Argo gateway) from each reaction's stoichiometry. Applied as an overlay on the baseline cascade: a reaction keeps its baseline direction unless the AI expressed an opinion, so the diff and FBA effect isolate the AI's calls. forward->'>', reverse->'<', reversible->'=', uncertain->'?'.",
    citations=('Claude Opus 4.8 (Anthropic), via the Argo gateway',),
    section='AI / LLM',
    order=210,
))

register(Variant(
    tag='eq3_beber2022',
    kind='overlay',
    direction_map='eq3_beber2022',
    title='eQuilibrator 3.0 directions (Beber 2022)',
    apt_title='eQuilibrator-3.0 component-contribution directionality (Beber 2022), |log10 gamma| = 3',
    description="Reaction directions from eQuilibrator 3.0's component-contribution energies (Beber 2022) using its reversibility index at the standard |log10 gamma| = 3 cutoff (Noor 2012/2013). eQuilibrator 3.0 is the current state-of-the-art ΔG′ estimator (full covariance, Mg2+/pH corrections), so this swaps the baseline's bundled group-contribution energies for the newest data. Overlaid on baseline.",
    citations=('Beber 2022', 'Noor 2012', 'Noor 2013'),
    section='§ New — eQuilibrator 3.0',
    order=220,
))

register(Variant(
    tag='eq3_gamma1',
    kind='overlay',
    direction_map='eq3_gamma1',
    title='eQuilibrator 3.0 directions, |log10 gamma| = 1 (aggressive)',
    apt_title='eQuilibrator-3.0 energies at an aggressive 10-fold cutoff (Beber 2022 + Noor 2013)',
    description='The same eQuilibrator-3.0 component-contribution energies as eq3_beber2022, but with an aggressive reversibility-index cutoff |log10 gamma| = 1 (a reaction is directional once a ~10-fold concentration shift would reverse it). Combines the newest energy data with the most aggressive threshold, so it makes the largest direction change of the new set. Overlaid on baseline.',
    citations=('Beber 2022', 'Noor 2012', 'Noor 2013'),
    section='§ New — eQuilibrator 3.0',
    order=230,
))

register(Variant(
    tag='consensus_thermo',
    kind='derived',
    combine=_majority_vote,
    title='Consensus of GC + eQuilibrator-3.0 + dGPredictor',
    apt_title='High-confidence consensus: majority vote of three ΔG′ estimators',
    description='Majority vote across three independent reaction-ΔG′ direction estimates: the group-contribution baseline, eQuilibrator-3.0 component contribution (Beber 2022), and dGPredictor (Wang 2021). A reaction is called directional only when at least two of the three agree; single votes that disagree collapse to reversible. This keeps the baseline direction unless the two newer estimators outvote group contribution, yielding a conservative, high-confidence direction set.',
    citations=('Beber 2022', 'Wang 2021', 'Noor 2012'),
    section='§ New — consensus ensemble',
    order=240,
))

register(Variant(
    tag='kegg_implicit',
    kind='overlay',
    direction_map='kegg_implicit',
    title='KEGG implicit directions (as-built model bounds)',
    apt_title='The reaction directions the KEGG core models ship with (on-disk flux bounds), overlaid on the baseline cascade',
    description="The direction each reaction was actually built with in the 5,683 KEGG core models, read straight from their on-disk COBRA flux bounds (lower_bound/upper_bound): (0,1000)->'>', (-1000,0)->'<', (-1000,1000)->'='. Taken as the per-reaction majority across every model containing the reaction (blocked bounds ignored; models are written in MSDB-canonical orientation, verified 0 flips). Applied as an overlay on the baseline cascade: a reaction keeps its baseline direction unless the models express one, so the diff and FBA effect isolate the models' own directionality.",
    citations=('ModelSEED / KBase core models — on-disk flux bounds',),
    section='Model as-built',
    order=250,
))

register(Variant(
    tag='group_contribution',
    kind='overlay',
    direction_map='group_contribution',
    title='Group Contribution directions (Jankowski 2008)',
    apt_title='Group-contribution method reaction directionality (Jankowski 2008), MSDB default thermodynamic source',
    description="Reaction directions from the Group Contribution method (Jankowski 2008), the thermodynamic source MSDB shipped with prior to eQuilibrator/dGPredictor integration. Read directly from each reaction's bundled 'Group contribution' dG/operator triple (results/rxn_directions_group-contribution.csv, 25,812 of ~56,000 reactions with a non-sentinel estimate). Overlaid on baseline.",
    citations=('Jankowski 2008',),
    section='§ Thermodynamic source comparison',
    order=260,
))

