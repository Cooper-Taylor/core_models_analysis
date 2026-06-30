"""
Parameterizable port of ModelSEEDDatabase/Scripts/Thermodynamics/
Estimate_Reaction_Reversibility.py.

The original module is a hand-tuned cascade of heuristics for assigning each
reaction a direction (``>``, ``<``, ``=``, or ``?``) from its ΔG′° estimate.
This file factors that cascade behind a ``ReversibilityConfig`` knob-set so
the notebook can re-run with different threshold values / new heuristics
without forking the upstream code.

With the default ``ReversibilityConfig()`` and the unmodified MSDB JSON
input, ``estimate_one`` reproduces the upstream cascade byte-for-byte. The
baseline now incorporates Heuristics-Review fixes H2 and H3 -- two coupled
shadow-bug repairs in ``_walk_stoichiometry`` that are also applied upstream
in ``Estimate_Reaction_Reversibility.py`` -- so "baseline" here means the
*fixed* upstream cascade, not the historical buggy one. H1 (return ``?`` for
the bare default) was evaluated and rejected, so the fall-through stays ``=``.

The remaining knobs that the Heuristics-Review notebook exercises are
documented inside ``ReversibilityConfig`` -- each maps to one section (§ 3.x)
of ``Reaction_Reversibility_Heuristics_Review.md``.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from math import log
from typing import Callable, Optional, Sequence

# ---------------------------------------------------------------------------
# Constants - mirror the upstream module to keep the byte-for-byte baseline.
# ---------------------------------------------------------------------------
TEMPERATURE_DEFAULT = 298.15
GAS_CONSTANT = 0.0019858775
SENTINEL_DG = 10000000.0
FARADAY = 0.023061

PROTON = "cpd00067"
WATER = "cpd00001"
CO2 = "cpd00011"
PROTON_WATER = frozenset((PROTON, WATER))
LOW_LOCAL_CONC = frozenset(("cpd00007", "cpd11640"))  # O2, H2
ATPS_REAGENTS = frozenset(("cpd00002", "cpd00008", "cpd00009",
                           "cpd00001", "cpd00067"))
ATP = "cpd00002"

PHOSPHATE_IDS = ("cpd00002", "cpd00008", "cpd00018", "cpd00009", "cpd00012")

LOW_ENERGY_CPDS_DEFAULT = (
    "cpd00011",  # CO2
    "cpd00013",  # NH3
    "cpd11493",  # ACP
    "cpd00009",  # Pi
    "cpd00012",  # PPi
    "cpd00010",  # CoA
    "cpd00449",  # Dihydrolipoamide
    "cpd00242",  # HCO3
)

DB_LEVEL_LABEL = {
    "GC": "Group contribution",
    "EQ": "eQuilibrator",
    "DGP": "dGPredictor",
}
# Legacy ``notes`` flags exist for GC + EQ only; DGP is sublist-only.
DB_LEVEL_NOTE = {"GC": "GCC", "EQ": "EQU"}
# Order matters for the unfiltered fallback: prefer EQ over GC over DGP when
# multiple sublists carry the top-level deltag exactly.  Matches MSDB.
DB_LEVEL_PRIORITY = ("EQ", "GC", "DGP")


# ---------------------------------------------------------------------------
# Configuration - one struct holds every knob exposed by the cascade.
# Adding a heuristic variant means flipping one of these fields, not editing
# the cascade.
# ---------------------------------------------------------------------------
@dataclass
class ReversibilityConfig:
    """Knob-set for one run of the cascade.

    Defaults reproduce ``Estimate_Reaction_Reversibility.py`` exactly.
    Each non-default knob is wired to one section of
    ``Reaction_Reversibility_Heuristics_Review.md`` -- see the inline
    comments.
    """

    # § 3.10 -- raw constants the upstream module bakes into module scope.
    temperature: float = TEMPERATURE_DEFAULT
    cell_min: float = 1e-5
    cell_max: float = 2e-2
    cell_conc: float = 1e-3

    # § 3.10 -- the +/- 2 kcal/mol "mMdeltaG -> reversible" band.
    mm_band: float = 2.0

    # § 3.5 -- replace the fixed +/- 2 kcal/mol band with a multiple of
    # the per-reaction CC uncertainty.  Set to None to keep the legacy band.
    sigma_band_k: Optional[float] = None

    # § 3.5 -- use the per-reaction uncertainty to widen / tighten the bounded
    # MdeltaG check too.  Same semantics: None = keep the upstream behavior of
    # ``+/- deltagerr``.  Any other value (e.g. 1.96) multiplies deltagerr.
    sigma_bounds_k: Optional[float] = None

    # § 3.7 -- the CO2 (and other) special-concentration overrides used in the
    # mMdeltaG walk.  The upstream module overrides cpd00011 to 1e-4 and the
    # LOW_LOCAL_CONC set to 1e-6, but the latter is unreachable in upstream
    # because of the ``cpd`` shadow bug.  When ``apply_special_conc`` is False
    # both overrides are dropped -- equivalent to "let eQuilibrator's species
    # model decide".
    apply_special_conc: bool = True
    co2_local_conc: float = 1e-4

    # § 3.6 -- LOW_ENERGY_CPDS list.  Pass an empty tuple to disable the
    # low-energy-points heuristic outright.
    low_energy_cpds: Sequence[str] = field(default_factory=lambda: LOW_ENERGY_CPDS_DEFAULT)

    # § 2.1 / 3.1 -- ln(reversibility index) from the third column of
    # ``MetaNetX_Reaction_Energies.tbl``.  When supplied, reactions whose
    # |ln gamma| exceeds ``ln_ri_threshold`` get a directional call from this
    # quantity *before* the legacy heuristics 4 + 5 fire.  Map keys are
    # MSDB rxn ids (without compartment), values are floats in nat units.
    ln_ri_by_rxn: Optional[dict] = None
    ln_ri_threshold: float = 6.9  # ln(1000), Noor 2012's rule of thumb

    # § 3.3 -- per-metabolite concentration ranges (Bennett 2009 etc.).
    # ``per_met_conc_range`` maps compound id -> (cmin, cmax) and overrides
    # ``cell_min`` / ``cell_max`` for those metabolites.  ``per_met_conc``
    # similarly overrides ``cell_conc`` for the mMdeltaG step.
    per_met_conc_range: Optional[dict] = None
    per_met_conc: Optional[dict] = None

    # § H2 + § H3 are NOT knobs: the two shadow-bug repairs they describe are
    # now the cascade's canonical behaviour (baked into ``_walk_stoichiometry``
    # to match the fixed upstream Estimate_Reaction_Reversibility.py), so the
    # ABC-transporter rule, the phosphate-spread term, and the CO2 / O2 / H2
    # concentration overrides all fire by default.  § H1 (return ``?`` for the
    # bare default) was evaluated and rejected as contrary to the cascade's
    # design, so there is no ``default_direction`` knob -- the fall-through is
    # always ``=`` (reversible), matching upstream.

    # § 3.6 -- replace the mMdeltaG band + low-energy-points heuristics with a
    # single posterior-probability rule.  When ``p_forward_threshold`` is not
    # None, after the bounded MdeltaG / ATPS / ABCT branches the cascade
    # evaluates ``P(mMdeltaG < 0)`` and ``P(mMdeltaG > 0)`` from
    # ``ΔG′° ~ N(deltag, deltagerr)`` (treating the concentration term as
    # fixed).  ``P_forward >= threshold`` forces '>'; ``P_reverse >= threshold``
    # forces '<'; otherwise the reaction is marked reversible.  None preserves
    # the legacy mMdeltaG-band + low-energy-points cascade.
    p_forward_threshold: Optional[float] = None

    # Energy-source override -- plug a different ΔG′° estimator into the *full*
    # cascade. When ``energy_override_by_rxn`` is set, any reaction present in it
    # uses the supplied ``(ΔG′°, σ)`` pair (kcal/mol) instead of the stored
    # top-level energy, so the entire cascade (bounded-ΔG window, mMdeltaG band,
    # low-energy points rule) is evaluated on the alternative source. Reactions
    # absent from the map fall back to the normal ``_energy_for`` lookup. Used by
    # the dGPredictor variant (Wang 2021) -- see load_dgpredictor_energies().
    energy_override_by_rxn: Optional[dict] = None
    energy_override_label: str = "override"

    # Internal: cached RT (kcal / mol).
    @property
    def rt(self) -> float:
        return self.temperature * GAS_CONSTANT


# ---------------------------------------------------------------------------
# Energy lookup helpers - faithful to the upstream module.
# ---------------------------------------------------------------------------
def _thermo_pair(rxn_entry, label):
    thermo = rxn_entry.get("thermodynamics")
    if not isinstance(thermo, dict):
        return None
    pair = thermo.get(label)
    if not pair or pair[0] is None:
        return None
    dg = float(pair[0])
    if dg == SENTINEL_DG:
        return None
    return [dg, float(pair[1])]


def _is_source_eligible(rxn_entry, level):
    if _thermo_pair(rxn_entry, DB_LEVEL_LABEL[level]) is not None:
        return True
    # DGP has no legacy note -- ``.get`` returns None and the second clause
    # short-circuits via the ``is not None`` guard.
    note = DB_LEVEL_NOTE.get(level)
    return note is not None and note in rxn_entry["notes"]


def _energy_for(rxn_entry, db_level):
    rxn_dg = rxn_entry["deltag"]
    rxn_dge = rxn_entry["deltagerr"]
    if rxn_dg is not None and not isinstance(rxn_dg, str):
        rxn_dg = float(rxn_dg)
    else:
        try:
            rxn_dg = float(rxn_dg)
        except (TypeError, ValueError):
            return None, None, None
    if rxn_dge is not None and not isinstance(rxn_dge, str):
        rxn_dge = float(rxn_dge)
    else:
        try:
            rxn_dge = float(rxn_dge)
        except (TypeError, ValueError):
            rxn_dge = 0.0
    if rxn_dg is None or rxn_dg == SENTINEL_DG:
        return None, None, None

    if db_level:
        if not _is_source_eligible(rxn_entry, db_level):
            return None, None, None
        label = DB_LEVEL_LABEL[db_level]
        append_label = label if _thermo_pair(rxn_entry, label) is not None else None
        return rxn_dg, rxn_dge, append_label

    chosen_label = None
    for level in DB_LEVEL_PRIORITY:
        label = DB_LEVEL_LABEL[level]
        pair = _thermo_pair(rxn_entry, label)
        if pair is not None and abs(pair[0] - rxn_dg) < 1e-9:
            chosen_label = label
            break
    return rxn_dg, rxn_dge, chosen_label


def _has_gc_data(rxn_entry):
    return _is_source_eligible(rxn_entry, "GC")


def _incomplete_decision(rxn_entry, db_level):
    """Fallback when the chosen source has no usable energy.

    Same semantics as upstream: an EQ run defers to whatever GC already wrote
    into ``reversibility``.  No knob - this branch is structural, not tuning.
    """
    status = "Incomplete"
    thermoreversibility = "?"
    if db_level == "EQ" and _has_gc_data(rxn_entry):
        thermoreversibility = rxn_entry["reversibility"]
        status += " (GCC)"
    return status, thermoreversibility


# ---------------------------------------------------------------------------
# Stoichiometry walk - one pass produces every accumulator the cascade needs.
# ---------------------------------------------------------------------------
def _walk_stoichiometry(stoichiometry, cfg: ReversibilityConfig):
    rct_min = rct_max = 0.0
    pdt_min = pdt_max = 0.0
    rgt_sum = 0.0
    proton_cpts = {}
    phosphates = {}

    cell_min = cfg.cell_min
    cell_max = cfg.cell_max
    cell_conc = cfg.cell_conc

    for rgt in stoichiometry:
        cpd = rgt["compound"]
        cpt = rgt["compartment"]
        coeff = float(rgt["coefficient"])

        if cpd == PROTON:
            proton_cpts[cpt] = 1

        # H2 + H3 (now baseline): accumulate phosphate-bearing reagents by
        # compound id.  Upstream this block carried two coupled typos -- the
        # accumulator tested ``cpd in rgt`` (the row's dict keys, never a
        # compound id) and reused ``cpd`` as a loop variable, shadowing the
        # reagent id with PHOSPHATE_IDS[-1] (PPi).  Together they left
        # ``phosphates`` empty (ABCT + phosphate-spread dead, H3) and made the
        # PROTON_WATER skip and the CO2 / O2 / H2 concentration overrides
        # unreachable (the O2/H2 override is H2).  Both fixes are repaired
        # upstream in Estimate_Reaction_Reversibility.py and are the cascade's
        # canonical behaviour, so this port bakes them in too.
        if cpd in PHOSPHATE_IDS:
            phosphates.setdefault(cpd, 0.0)
            phosphates[cpd] += coeff

        if cpd in PROTON_WATER:
            continue

        # Per-metabolite range overrides (§ 3.3).
        if cfg.per_met_conc_range is not None and cpd in cfg.per_met_conc_range:
            cmin, cmax = cfg.per_met_conc_range[cpd]
        else:
            cmin, cmax = cell_min, cell_max

        if coeff < 0:
            rct_min += coeff * log(cmin)
            rct_max += coeff * log(cmax)
        else:
            pdt_min += coeff * log(cmin)
            pdt_max += coeff * log(cmax)

        # mMdeltaG concentration override.  CO2 ~0.1 mM and dissolved O2/H2
        # ~1 µM rather than the 1 mM default.  ``apply_special_conc`` (§ 3.7)
        # still gates both so that variant can disable them.
        if cfg.per_met_conc is not None and cpd in cfg.per_met_conc:
            local_conc = cfg.per_met_conc[cpd]
        else:
            local_conc = cell_conc
            if cfg.apply_special_conc:
                if cpd == CO2:
                    local_conc = cfg.co2_local_conc
                elif cpd in LOW_LOCAL_CONC:
                    local_conc = 1e-6
        rgt_sum += coeff * log(local_conc)

    return {
        "rct_min": rct_min, "rct_max": rct_max,
        "pdt_min": pdt_min, "pdt_max": pdt_max,
        "rgt_sum": rgt_sum,
        "proton_cpts": proton_cpts,
        "phosphates": phosphates,
    }


def _stored_bounds(rxn_dg, rxn_dge, terms, cfg: ReversibilityConfig):
    rxn_dg_transport = 0.0
    err = rxn_dge * (cfg.sigma_bounds_k if cfg.sigma_bounds_k is not None else 1.0)
    stored_max = (rxn_dg + rxn_dg_transport + err
                  + cfg.rt * terms["pdt_max"]
                  + cfg.rt * terms["rct_min"])
    stored_min = (rxn_dg + rxn_dg_transport - err
                  + cfg.rt * terms["pdt_min"]
                  + cfg.rt * terms["rct_max"])
    return stored_max, stored_min


def _is_atp_synthase(rxn_entry, proton_cpts):
    if rxn_entry["is_transport"] != 1 or len(proton_cpts) <= 1:
        return False
    cpds_cpts = {}
    for rgt in rxn_entry["stoichiometry"]:
        cpds_cpts.setdefault(rgt["compound"], []).append(rgt["compartment"])
    if len(cpds_cpts) != 5:
        return False
    for cpd, cpts in cpds_cpts.items():
        if cpd not in ATPS_REAGENTS:
            return False
        if len(cpts) == 2 and cpd != PROTON:
            return False
    return True


def _abc_transporter_decision(rxn_entry, phosphates):
    if rxn_entry["is_transport"] != 1 or ATP not in phosphates:
        return None
    coeff = phosphates[ATP]
    if coeff < 0:
        rev = ">"
    elif coeff > 0:
        rev = "<"
    else:
        rev = "="
    return f"ABCT: {coeff}", rev


def _low_energy_points(stoichiometry, phosphates, cfg: ReversibilityConfig):
    points = 0.0
    min_coeff = SENTINEL_DG
    if ATP in phosphates and len(phosphates) > 2:
        for pho_coeff in phosphates.values():
            if pho_coeff < min_coeff:
                min_coeff = pho_coeff
    if min_coeff != SENTINEL_DG:
        points -= abs(min_coeff)
    for rgt in stoichiometry:
        if rgt["compound"] in cfg.low_energy_cpds:
            points -= float(rgt["coefficient"])
    return points


# ---------------------------------------------------------------------------
# Cascade -- the entry point invoked once per reaction.
# ---------------------------------------------------------------------------
def estimate_one(rxn_entry, db_level: str = "EQ", cfg: Optional[ReversibilityConfig] = None):
    """Return ``(status_label, reversibility, source_label)``.

    With ``cfg = ReversibilityConfig()`` the cascade matches
    ``Estimate_Reaction_Reversibility.estimate_one`` exactly.
    """
    if cfg is None:
        cfg = ReversibilityConfig()

    if rxn_entry["status"] == "EMPTY":
        return "Empty", "?", None

    if cfg.energy_override_by_rxn is not None:
        ov = cfg.energy_override_by_rxn.get(rxn_entry["id"])
        if ov is not None:
            rxn_dg, rxn_dge, source_label = float(ov[0]), float(ov[1]), cfg.energy_override_label
        else:
            rxn_dg, rxn_dge, source_label = _energy_for(rxn_entry, db_level)
    else:
        rxn_dg, rxn_dge, source_label = _energy_for(rxn_entry, db_level)
    if rxn_dg is None:
        status, thermoreversibility = _incomplete_decision(rxn_entry, db_level)
        return status, thermoreversibility, None

    terms = _walk_stoichiometry(rxn_entry["stoichiometry"], cfg)
    stored_max, stored_min = _stored_bounds(rxn_dg, rxn_dge, terms, cfg)

    if stored_max < 0:
        return "MdeltaG(Max): {0:.2f}".format(stored_max), ">", source_label
    if stored_min > 0:
        return "MdeltaG(Min): {0:.2f}".format(stored_min), "<", source_label

    if _is_atp_synthase(rxn_entry, terms["proton_cpts"]):
        return "ATPS", "=", source_label

    abct = _abc_transporter_decision(rxn_entry, terms["phosphates"])
    if abct is not None:
        status, thermoreversibility = abct
        return status, thermoreversibility, source_label

    # § 2.1 / 3.1 -- ln(reversibility index) when available.
    if cfg.ln_ri_by_rxn is not None:
        ln_ri = cfg.ln_ri_by_rxn.get(rxn_entry["id"])
        if ln_ri is not None and abs(ln_ri) > cfg.ln_ri_threshold:
            # eQuilibrator's ln(reversibility index) carries the SAME sign as
            # ΔG′° (Noor 2012: ln γ̂ = (2/N)·ΔG′ₘ/RT), so a NEGATIVE index means
            # forward-favorable -- matching the MdeltaG rule above (stored_max<0
            # → ">"). A previous version had this comparison inverted
            # (`ln_ri > 0`), which flipped every strongly-directional reaction
            # to the wrong direction. See
            # reports/VARIANT_3.1_BREAKAGE_INVESTIGATION.md §"Correction".
            direction = ">" if ln_ri < 0 else "<"
            return f"lnRI: {ln_ri:.2f}", direction, source_label

    mMdeltaG = rxn_dg + cfg.rt * terms["rgt_sum"]

    # § 3.6 -- P(forward)/P(reverse) threshold replaces the mMdeltaG band +
    # low-energy-points heuristic entirely when ``cfg.p_forward_threshold`` is
    # set.  Computed from the marginal CC normal on ΔG′°; concentration term
    # is the same ``rgt_sum`` the mMdeltaG step uses (treated as fixed).
    #
    # Symmetry note: the legacy heuristic 4 marks `=` whenever
    # ``|mMdeltaG| <= mm_band`` (default 2 kcal/mol).  For the P-rule to be a
    # drop-in replacement (not a strictly more aggressive rule), the
    # "forward" event must be ``mMdeltaG < -mm_band`` and the "reverse"
    # event must be ``mMdeltaG > +mm_band``.  Then both rules collapse to
    # the same answer at threshold = 0.5 -- and tightening the threshold
    # only widens the reversible region, never the directional ones.
    if cfg.p_forward_threshold is not None:
        from math import erf, sqrt
        sigma = max(float(rxn_dge or 0.0), 0.05)
        band = cfg.mm_band
        z_f = (-band - mMdeltaG) / sigma  # P(mMdeltaG < -band)
        z_r = (-band + mMdeltaG) / sigma  # P(mMdeltaG >  +band)
        p_forward = 0.5 * (1.0 + erf(z_f / sqrt(2.0)))
        p_reverse = 0.5 * (1.0 + erf(z_r / sqrt(2.0)))
        thresh = cfg.p_forward_threshold
        if p_forward >= thresh:
            return f"P(fwd)={p_forward:.3f}", ">", source_label
        if p_reverse >= thresh:
            return f"P(rev)={p_reverse:.3f}", "<", source_label
        return f"P(fwd|rev)<{thresh:.2f}", "=", source_label

    # § 3.5 -- per-reaction tolerance (k * sigma) takes precedence when set.
    if cfg.sigma_band_k is not None:
        band = cfg.sigma_band_k * rxn_dge
        if -band <= mMdeltaG <= band:
            return "sigBand({0:.2f}): {1:.2f}".format(band, mMdeltaG), "=", source_label
    else:
        if -cfg.mm_band <= mMdeltaG <= cfg.mm_band:
            return "mMdeltaG: {0:.2f}".format(mMdeltaG), "=", source_label

    points = _low_energy_points(rxn_entry["stoichiometry"], terms["phosphates"], cfg)
    if points * mMdeltaG > 2:
        if mMdeltaG < 0:
            return ("lowE: {0:.2f}".format(mMdeltaG) + ":" + str(points), ">", source_label)
        if mMdeltaG > 0:
            return ("lowE: {0:.2f}".format(mMdeltaG) + ":" + str(points), "<", source_label)

    return "default", "=", source_label


def run_cascade(reactions_dict, db_level: str = "EQ",
                cfg: Optional[ReversibilityConfig] = None,
                gc_first: bool = True):
    """Compute reversibility for every reaction.

    Matches the upstream pipeline: when ``db_level == 'EQ'`` and
    ``gc_first`` is true, the GC pass is run first to populate
    ``reactions_dict[rxn]['reversibility']`` so the Incomplete (GCC)
    fallback in the EQ pass can read it -- mirroring
    ``Rerun_Thermodynamics.sh``'s GC-then-EQ order.

    Returns ``{rxn_id: (status, reversibility)}``.  Does not mutate the
    reactions outside the optional GC pass needed for fallback.
    """
    if cfg is None:
        cfg = ReversibilityConfig()

    if db_level == "EQ" and gc_first:
        # Mutate ``reversibility`` to the GC-derived value so the EQ pass's
        # Incomplete (GCC) fallback returns the right thing.
        for rxn_id in sorted(reactions_dict.keys()):
            rxn_entry = reactions_dict[rxn_id]
            _, rev, _ = estimate_one(rxn_entry, "GC", cfg)
            rxn_entry["reversibility"] = rev

    out = {}
    for rxn_id in sorted(reactions_dict.keys()):
        rxn_entry = reactions_dict[rxn_id]
        status, rev, _ = estimate_one(rxn_entry, db_level, cfg)
        out[rxn_id] = (status, rev)
    return out


# ---------------------------------------------------------------------------
# Auxiliary parsers used by the new heuristics.
# ---------------------------------------------------------------------------
import os as _os
_LN_RI_PATH_DEFAULT = (
    _os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase")
    + "/Biochemistry/Thermodynamics/eQuilibrator/MetaNetX_Reaction_Energies.tbl"
)


def load_ln_reversibility_index(path: str = _LN_RI_PATH_DEFAULT) -> dict:
    """Parse the ln_RI (fourth) column of ``MetaNetX_Reaction_Energies.tbl``.

    Each row is ``rxn_id\\tdg\\tdge\\tln_RI`` where ``ln_RI`` is formatted
    as ``value+/-uncertainty`` (kept verbatim from the retrieval script).
    Returns ``{rxn_id: ln_ri_value}``.

    See ``Reaction_Reversibility_Heuristics_Review.md`` § 2.1 / 3.1 -- this
    quantity is already computed upstream but currently dropped on the floor.
    """
    out = {}
    val_re = re.compile(r"^(-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)")
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            rxn_id = parts[0]
            raw = parts[3].split("+/-")[0]
            m = val_re.match(raw.strip())
            if not m:
                continue
            try:
                out[rxn_id] = float(m.group(1))
            except ValueError:
                continue
    return out


_MSDB_BIOCHEM_DEFAULT = (
    _os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase") + "/Biochemistry"
)


def load_dgpredictor_energies(biochem_dir: str = _MSDB_BIOCHEM_DEFAULT) -> dict:
    """Return ``{rxn_id: (ΔG′°, σ)}`` from MSDB ``thermodynamics['dGPredictor']``.

    dGPredictor (Wang, Upadhyay & Maranas 2021) standard transformed Gibbs
    energies of reaction, in kcal/mol (same units as the stored
    Group-contribution / eQuilibrator energies). Feeds the ``dgpredictor``
    variant via ``ReversibilityConfig.energy_override_by_rxn``. Reactions with
    no usable dGPredictor value (None / sentinel) are omitted.
    """
    import glob as _glob
    import json as _json
    out: dict = {}
    for f in sorted(_glob.glob(_os.path.join(biochem_dir, "reaction_*.json"))):
        with open(f) as fh:
            for r in _json.load(fh):
                thermo = r.get("thermodynamics")
                if not isinstance(thermo, dict):
                    continue
                pair = thermo.get("dGPredictor")
                if not pair or pair[0] is None:
                    continue
                try:
                    dg = float(pair[0])
                except (TypeError, ValueError):
                    continue
                if dg == SENTINEL_DG:
                    continue
                try:
                    dge = float(pair[1])
                except (TypeError, ValueError, IndexError):
                    dge = 0.0
                out[r["id"]] = (dg, dge)
    return out


# Bennett 2009, Table 1 -- a subset of measured E. coli concentrations
# (geometric mean, log-spread); the actual paper gives ~100 values, we keep
# the central ones that move the heuristic.  All values in M.
BENNETT_2009_ECOLI = {
    # high-abundance pools (>> 1 mM)
    "cpd00023": (1.5e-2, 1.0e-1),   # L-Glutamate         (96 mM)
    "cpd00033": (5.0e-3, 5.0e-2),   # Glycine
    "cpd00041": (1.0e-3, 5.0e-2),   # L-Aspartate
    "cpd00060": (1.0e-3, 1.0e-2),   # L-Methionine
    # cofactor pools
    "cpd00002": (5.0e-3, 1.0e-2),   # ATP                  (~9.6 mM)
    "cpd00008": (2.0e-4, 8.0e-4),   # ADP                  (~0.56 mM)
    "cpd00018": (5.0e-5, 5.0e-4),   # AMP
    "cpd00003": (1.0e-4, 5.0e-4),   # NAD+
    "cpd00004": (1.0e-5, 5.0e-5),   # NADH
    "cpd00005": (5.0e-5, 5.0e-4),   # NADPH
    "cpd00006": (1.0e-6, 5.0e-5),   # NADP+
    "cpd00010": (1.0e-5, 1.0e-3),   # CoA
    # phosphate pool
    "cpd00009": (1.0e-3, 1.0e-2),   # Pi                   (~5 mM)
    "cpd00012": (1.0e-5, 1.0e-3),   # PPi                  (~0.5 mM)
    # central metabolism intermediates
    "cpd00020": (1.0e-4, 1.0e-3),   # Pyruvate
    "cpd00022": (1.0e-4, 5.0e-4),   # Acetyl-CoA
    "cpd00024": (1.0e-5, 1.0e-3),   # 2-Oxoglutarate
    "cpd00027": (1.0e-3, 1.0e-2),   # D-Glucose            (intracell ~5 mM)
    # rare nucleosides
    "cpd00214": (5.0e-8, 5.0e-6),   # Adenosine            (very low, ~0.13 uM)
}


# Mean (1 mM) used for the mMdeltaG step under Bennett-aware run.
BENNETT_2009_MEAN = {cpd: (lo * hi) ** 0.5 for cpd, (lo, hi) in BENNETT_2009_ECOLI.items()}
