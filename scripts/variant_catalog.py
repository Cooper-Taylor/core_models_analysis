"""Canonical list of ReversibilityConfig variants exercised by notebook 06.

Single source of truth shared by:
  - ``build_reversibility_notebook.py``  (the notebook generator)
  - ``export_thermo_variants.py``        (the MSDB-format report writer)
  - ``build_site_data.py``               (the website JSON builder)

Each entry has:
  ``tag``      -- short id (matches the notebook cache key suffix)
  ``title``    -- one-line description shown to the user
  ``section``  -- pointer back to Reaction_Reversibility_Heuristics_Review.md
  ``cfg``      -- a callable that returns a ``ReversibilityConfig``

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
    {"tag": "baseline", "title": "ReversibilityConfig() default (matches MSDB)",
     "section": "(reference)", "cfg": _baseline_cfg},
    {"tag": "3.1", "title": "Persist + use ln(reversibility_index) (Noor 2012)",
     "section": "§ 2.1 / 3.1", "cfg": _v31_cfg},
    {"tag": "3.3", "title": "Bennett-2009 per-metabolite concentration ranges",
     "section": "§ 3.3", "cfg": _v33_cfg},
    {"tag": "3.3_wide", "title": "Wider uniform conc window [1e-7, 0.1] M",
     "section": "§ 3.3 (fallback)", "cfg": _v33w_cfg},
    {"tag": "3.5", "title": "Per-reaction sigma band: k=1.96 (95%) replaces ±2 kcal",
     "section": "§ 3.5", "cfg": _v35_cfg},
    {"tag": "3.5_wide", "title": "Per-reaction CC bound widening: k=1.96 on stored_bounds",
     "section": "§ 3.5 / § 2.5", "cfg": _v35w_cfg},
    {"tag": "3.6", "title": "Drop the low-energy-compounds list entirely",
     "section": "§ 3.6", "cfg": _v36_cfg},
    {"tag": "3.7", "title": "Drop the CO2 1e-4 hardcoded concentration override",
     "section": "§ 3.7", "cfg": _v37_cfg},
    {"tag": "3.10_tight", "title": "Tighten mMdeltaG band: ±1 kcal/mol",
     "section": "§ 3.10", "cfg": _v310t_cfg},
    {"tag": "3.10_loose", "title": "Loosen mMdeltaG band: ±4 kcal/mol",
     "section": "§ 3.10", "cfg": _v310l_cfg},
    {"tag": "H1", "title": "(NEW) default direction = '?' for unresolved",
     "section": "§ H1", "cfg": _H1_cfg},
    {"tag": "H2", "title": "(NEW) repair LOW_LOCAL_CONC shadow bug (O2/H2 at 1e-6 M)",
     "section": "§ H2", "cfg": _H2_cfg},
    {"tag": "H3", "title": "(NEW) repair phosphates shadow bug (ABC + low-E phosphate spread)",
     "section": "§ H3", "cfg": _H3_cfg},
    {"tag": "H4", "title": "(NEW) best-evidence composite: 3.1 + 3.5 + Bennett",
     "section": "§ H4", "cfg": _H4_cfg},
]


def variant_by_tag(tag: str) -> dict:
    for v in VARIANTS:
        if v["tag"] == tag:
            return v
    raise KeyError(f"unknown variant tag: {tag!r}")
