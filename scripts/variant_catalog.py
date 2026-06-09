"""Canonical list of ReversibilityConfig variants exercised by notebook 06.

Single source of truth shared by:
  - ``build_reversibility_notebook.py``  (the notebook generator)
  - ``export_thermo_variants.py``        (the MSDB-format report writer)
  - ``build_site_data.py``               (the website JSON builder)

Each entry has:
  ``tag``         -- short id (matches the notebook cache key suffix)
  ``title``       -- short legacy label (one-line; kept for backwards compat
                     with anything still reading the old field; new code
                     should prefer ``apt_title``)
  ``apt_title``   -- a descriptive one-line title shown to website users
  ``description`` -- 2-4 sentence technical description: what baseline does,
                     what the variant changes, and why. Self-contained for
                     a reader who has not read the heuristics review doc.
  ``citations``   -- list of citation keys appearing in §5 References of
                     Reaction_Reversibility_Heuristics_Review.md. Empty for
                     variants that are pure bug-fixes or default-tweaks.
  ``section``     -- pointer back to Reaction_Reversibility_Heuristics_Review.md
  ``cfg``         -- a callable that returns a ``ReversibilityConfig``

The callable form lets variants that need on-disk data (e.g. H4 / 3.1
loading ``ln_reversibility_index``) defer the load until they're actually
exercised.
"""

from __future__ import annotations

from typing import Callable

import reversibility_lib as lib


def _baseline_cfg() -> lib.ReversibilityConfig:
    return lib.ReversibilityConfig()


def _v31_cfg() -> lib.ReversibilityConfig:
    return lib.ReversibilityConfig(ln_ri_by_rxn=lib.load_ln_reversibility_index())


def _v33_cfg() -> lib.ReversibilityConfig:
    return lib.ReversibilityConfig(
        per_met_conc_range=lib.BENNETT_2009_ECOLI,
        per_met_conc=lib.BENNETT_2009_MEAN,
    )


def _v33w_cfg() -> lib.ReversibilityConfig:
    return lib.ReversibilityConfig(cell_min=1e-7, cell_max=1e-1)


def _v35_cfg() -> lib.ReversibilityConfig:
    return lib.ReversibilityConfig(sigma_band_k=1.96)


def _v35w_cfg() -> lib.ReversibilityConfig:
    return lib.ReversibilityConfig(sigma_bounds_k=1.96)


def _v36_cfg() -> lib.ReversibilityConfig:
    return lib.ReversibilityConfig(low_energy_cpds=())


def _v37_cfg() -> lib.ReversibilityConfig:
    return lib.ReversibilityConfig(apply_special_conc=False)


def _v310t_cfg() -> lib.ReversibilityConfig:
    return lib.ReversibilityConfig(mm_band=1.0)


def _v310l_cfg() -> lib.ReversibilityConfig:
    return lib.ReversibilityConfig(mm_band=4.0)


def _H1_cfg() -> lib.ReversibilityConfig:
    return lib.ReversibilityConfig(default_direction="?")


def _H2_cfg() -> lib.ReversibilityConfig:
    return lib.ReversibilityConfig(fix_low_local_conc=True)


def _H3_cfg() -> lib.ReversibilityConfig:
    return lib.ReversibilityConfig(fix_phosphates_shadow=True)


def _H4_cfg() -> lib.ReversibilityConfig:
    return lib.ReversibilityConfig(
        ln_ri_by_rxn=lib.load_ln_reversibility_index(),
        sigma_band_k=1.96,
        per_met_conc_range=lib.BENNETT_2009_ECOLI,
        per_met_conc=lib.BENNETT_2009_MEAN,
    )


VARIANTS: list[dict] = [
    {
        "tag": "baseline",
        "title": "ReversibilityConfig() default (matches MSDB)",
        "apt_title": "Default cascade — reproduces MSDB byte-for-byte (reference)",
        "description": (
            "Reproduces the upstream ModelSEEDDatabase "
            "Estimate_Reaction_Reversibility.py cascade exactly: heuristic 1 "
            "(bounded ΔG′° check across the concentration window) "
            "→ heuristic 4 (mMdeltaG ±2 kcal/mol band) → heuristic 5 "
            "(LOW_ENERGY_CPDS points rule) → '=' (reversible) default. "
            "Every other variant is a one-knob (or few-knob) diff against this reference."
        ),
        "citations": ["Henry 2007"],
        "section": "(reference)",
        "cfg": _baseline_cfg,
    },
    {
        "tag": "3.1",
        "title": "Persist + use ln(reversibility_index) (Noor 2012)",
        "apt_title": "Use eQuilibrator's reversibility index instead of the points-and-bands heuristics",
        "description": (
            "Replaces the heuristic reversibility decision with eQuilibrator's "
            "reversibility index γ, calling a reaction directional whenever "
            "|ln γ| > 6.9 — i.e. when more than a ~1000-fold concentration "
            "shift would be needed to reverse it (cfg knob: ln_ri_by_rxn). γ is "
            "the molecularity-normalized fold-change in the mass-action ratio that "
            "flips the sign of ΔG′, so it is the formal answer to the "
            "question the ad-hoc points-and-bands rules approximate."
        ),
        "citations": ["Noor 2012"],
        "section": "§ 2.1 / 3.1",
        "cfg": _v31_cfg,
    },
    {
        "tag": "3.3",
        "title": "Bennett-2009 per-metabolite concentration ranges",
        "apt_title": "Use measured E. coli metabolite concentrations (Bennett 2009)",
        "description": (
            "Replaces the uniform [10 µM, 20 mM] concentration prior with "
            "per-metabolite ranges and geometric-mean reference points from "
            "Bennett 2009's absolute E. coli metabolomics (cfg knobs: "
            "per_met_conc_range, per_met_conc), so ATP/glutamate sit at their "
            "measured ~10/100 mM levels and rare nucleosides near 0.1 µM. "
            "The real intracellular distribution spans six orders of magnitude "
            "with 70% of metabolites below 1 mM, so the uniform prior "
            "systematically biases the reference ΔG used by the directional check."
        ),
        "citations": ["Bennett 2009", "Park 2016", "Bar-Even 2011", "Salvy 2019"],
        "section": "§ 3.3",
        "cfg": _v33_cfg,
    },
    {
        "tag": "3.3_wide",
        "title": "Wider uniform conc window [1e-7, 0.1] M",
        "apt_title": "Widen the uniform concentration prior to [0.1 µM, 0.1 M]",
        "description": (
            "Widens the baseline three-order [10 µM, 20 mM] uniform reagent "
            "concentration window to a six-order [0.1 µM, 0.1 M] window "
            "(cfg knobs: cell_min, cell_max), covering the observed spread of "
            "E. coli metabolomics without committing to per-metabolite "
            "measurements. This serves as a measured-data-free fallback that "
            "better envelopes Bennett 2009 observations for reactions whose "
            "intermediates are unmeasured."
        ),
        "citations": ["Bennett 2009", "Park 2016"],
        "section": "§ 3.3 (fallback)",
        "cfg": _v33w_cfg,
    },
    {
        "tag": "3.5",
        "title": "Per-reaction sigma band: k=1.96 (95%) replaces ±2 kcal",
        "apt_title": "Use a 95% Component Contribution (CC) uncertainty band instead of fixed ±2 kcal/mol",
        "description": (
            "Replaces the fixed ±2 kcal/mol reversible-band check (inherited "
            "from Henry 2007) with k·σ_rxn where k=1.96 — a 95% confidence "
            "interval drawn from the Component Contribution method's per-reaction "
            "ΔG uncertainty already on disk (cfg knob: sigma_band_k). The "
            "2 kcal/mol constant conflates ΔG-estimation error with concentration "
            "uncertainty, while CC's per-reaction σ tightens the band for "
            "well-measured reactions and widens it for poorly-constrained ones."
        ),
        "citations": ["Noor 2013", "Gollub 2021"],
        "section": "§ 3.5",
        "cfg": _v35_cfg,
    },
    {
        "tag": "3.5_wide",
        "title": "Per-reaction CC bound widening: k=1.96 on stored_bounds",
        "apt_title": "Widen the bounded-ΔG check to a 95% CC interval (1.96σ) instead of 1σ",
        "description": (
            "The directional check requires the ΔG window — computed at "
            "concentration extremes — to lie entirely on one side of zero. "
            "Baseline uses 1-σ error bars on the stored Component Contribution "
            "bounds; this variant scales them by 1.96 for a 95% confidence interval "
            "(cfg knob: sigma_bounds_k). Noor 2013's empirical coverage validation "
            "(73%/90%/95%/99% at matching CC intervals) shows that a 1-σ bound "
            "calls reactions directional more aggressively than the data support."
        ),
        "citations": ["Noor 2013", "Gollub 2021"],
        "section": "§ 3.5 / § 2.5",
        "cfg": _v35w_cfg,
    },
    {
        "tag": "3.6",
        "title": "Drop the low-energy-compounds list entirely",
        "apt_title": "Drop the hand-curated low-energy-compounds list (CO₂, Pi, PPi, CoA, …)",
        "description": (
            "Disables the legacy points rule that forced directional calls when "
            "reactions consumed hand-curated low-energy sinks (CO₂, NH₃, ACP, "
            "Pi, PPi, CoA, dihydrolipoamide, HCO₃⁻), by emptying the list "
            "(cfg knob: low_energy_cpds=()). The frozen MFAToolkit-era compound list "
            "is superseded by principled measures like the reversibility index "
            "(§3.1) or multiTFA's P(forward). The residual phosphate-spread term "
            "in the same rule still runs in principle but stays inert until "
            "§H3 is stacked in."
        ),
        "citations": ["Noor 2012", "Gollub 2021"],
        "section": "§ 3.6",
        "cfg": _v36_cfg,
    },
    {
        "tag": "3.7",
        "title": "Drop the CO2 1e-4 hardcoded concentration override",
        "apt_title": "Disable the hardcoded 1e-4 M CO₂ concentration override",
        "description": (
            "Baseline declares a 1e-4 M CO₂ override (and a 1e-6 M O₂/H₂ "
            "override) but both branches are unreachable due to the shadow bug "
            "fixed in §H3 — they only fire once H3 is stacked in. This variant "
            "turns the gate off (cfg knob: apply_special_conc=False) so that even "
            "after the H3 repair, CO₂/O₂/H₂ sit at the 1 mM default. "
            "eQuilibrator already models CO₂(aq)/HCO₃⁻/CO₃²⁻ "
            "speciation as a function of pH, so the override would double-count "
            "what the transform handles."
        ),
        "citations": ["Alberty 2003", "Beber 2022"],
        "section": "§ 3.7",
        "cfg": _v37_cfg,
    },
    {
        "tag": "3.10_tight",
        "title": "Tighten mMdeltaG band: ±1 kcal/mol",
        "apt_title": "Tighten the reversible band to ±1 kcal/mol",
        "description": (
            "Halves the reversible-band check from |reference ΔG| ≤ 2 kcal/mol "
            "to ≤ 1 kcal/mol (cfg knob: mm_band=1.0), forcing more reactions "
            "into directional calls. The 2 kcal/mol threshold was inherited from "
            "Henry 2007 without per-reaction justification; this variant exposes "
            "the cascade's sensitivity to that hardcoded threshold by making it "
            "a configurable input."
        ),
        "citations": [],
        "section": "§ 3.10",
        "cfg": _v310t_cfg,
    },
    {
        "tag": "3.10_loose",
        "title": "Loosen mMdeltaG band: ±4 kcal/mol",
        "apt_title": "Loosen the reversible band to ±4 kcal/mol",
        "description": (
            "Doubles the reversible-band check from |reference ΔG| ≤ 2 kcal/mol "
            "to ≤ 4 kcal/mol (cfg knob: mm_band=4.0), pulling more reactions "
            "out of directional calls into the reversible bucket. As with "
            "3.10_tight, the motivation is to expose the reversibility decision's "
            "sensitivity to the Henry 2007 threshold and quantify how many "
            "directional calls hinge on the 2 kcal/mol choice alone."
        ),
        "citations": [],
        "section": "§ 3.10",
        "cfg": _v310l_cfg,
    },
    {
        "tag": "H1",
        "title": "(NEW) default direction = '?' for unresolved",
        "apt_title": "Distinguish 'no rule fired' from 'reversible' by returning '?'",
        "description": (
            "The baseline's final fallthrough returns '=' (reversible) whenever "
            "no earlier rule fires, conflating two genuinely different outcomes: "
            "heuristics actively agreed the reaction is reversible, versus "
            "heuristics had nothing to say. This variant returns '?' (unknown) "
            "for the bare default branch (cfg knob: default_direction='?'). "
            "6,522 of 56,012 ModelSEED Database (MSDB) reactions land on this "
            "branch, so collapsing them with confidently-reversible reactions "
            "discards information a curator would want to triage."
        ),
        "citations": [],
        "section": "§ H1",
        "cfg": _H1_cfg,
    },
    {
        "tag": "H2",
        "title": "(NEW) repair LOW_LOCAL_CONC shadow bug (O2/H2 at 1e-6 M)",
        "apt_title": "Repair the O₂/H₂ shadow bug so oxidative reactions see 1 µM, not 1 mM",
        "description": (
            "Cytoplasmic dissolved O₂ sits at 5–50 µM, not the 1 mM default; "
            "the original code intended to apply a 1 µM local concentration "
            "override for O₂ and H₂ but a variable-shadowing bug made the "
            "branch unreachable. This variant repairs the bug "
            "(cfg knob: fix_low_local_conc=True). The override shifts the reference "
            "ΔG of every O₂/H₂-coupled reaction by ~4 kcal/mol per unit "
            "stoichiometric coefficient, with sign depending on whether the gas "
            "appears as reactant or product — large enough to flip many calls in "
            "or out of the ±2 kcal reversible bucket."
        ),
        "citations": [],
        "section": "§ H2",
        "cfg": _H2_cfg,
    },
    {
        "tag": "H3",
        "title": "(NEW) repair phosphates shadow bug (ABC + low-E phosphate spread)",
        "apt_title": "Repair the phosphate-accumulator shadow bug so ABC transporters become directional",
        "description": (
            "A typo silently disables the ABC-transporter rule that forces "
            "ATP-driven uptake reactions to be directional: the accumulator loop "
            "checks the wrong field on each reagent row, so the phosphates dict "
            "is always empty and the ABC branch is dead code. Repairing it "
            "(cfg knob: fix_phosphates_shadow=True) flips 1,989 ModelSEED Database "
            "reactions — 1,209 from '=' (reversible) to '>' (forward-only) as "
            "ATP-driven uptake is correctly forced forward — and changes growth "
            "status in 21 of 100 panel models."
        ),
        "citations": [],
        "section": "§ H3",
        "cfg": _H3_cfg,
    },
    {
        "tag": "H4",
        "title": "(NEW) best-evidence composite: 3.1 + 3.5 + Bennett",
        "apt_title": "Best-evidence stack: reversibility index + 95% CC band + Bennett concentrations",
        "description": (
            "Combines three high-impact §3 changes into one 'best available "
            "evidence' configuration: eQuilibrator's ln(reversibility index), a "
            "95% Component Contribution per-reaction σ band in place of the "
            "fixed ±2 kcal/mol, and Bennett 2009 per-metabolite E. coli "
            "concentrations (cfg knobs: ln_ri_by_rxn, sigma_band_k=1.96, "
            "per_met_conc_range, per_met_conc). The three interact non-trivially: "
            "the CC band pulls many calls into the reversible bucket, the "
            "reversibility index re-pushes a subset back to firm directions, and "
            "the measured concentrations shift the reference ΔG of those still "
            "in the band."
        ),
        "citations": [
            "Noor 2012", "Noor 2013", "Gollub 2021",
            "Bennett 2009", "Park 2016", "Bar-Even 2011", "Salvy 2019",
        ],
        "section": "§ H4",
        "cfg": _H4_cfg,
    },
]


def variant_by_tag(tag: str) -> dict:
    for v in VARIANTS:
        if v["tag"] == tag:
            return v
    raise KeyError(f"unknown variant tag: {tag!r}")
