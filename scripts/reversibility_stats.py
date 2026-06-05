"""Statistical extensions to the reversibility cascade.

Implements the parts of ``Reaction_Reversibility_Heuristics_Review.md`` that
quantify *uncertainty* rather than just flipping a heuristic:

  * §3.6  ``p_direction(...)``    -- analytic P(forward)/(reverse)/(reversible)
                                     from the marginal CC normal posterior on
                                     ΔG′°.  Per-reaction summary; cheap.
  * §3.6  ``p_forward_threshold``  -- new ReversibilityConfig knob (declared in
                                     reversibility_lib): replaces the ±2 kcal
                                     mMdeltaG band + the low-energy-points
                                     heuristic with a single
                                     P(forward) ≥ threshold rule.
  * §2.5  ``sample_cascade(...)``  -- Monte Carlo cascade replay.  Per sample,
                                     each reaction's ΔG′° is resampled from
                                     ``N(μ=deltag, σ=deltagerr)``; the full
                                     cascade fires with the resampled value.
                                     Returns one rev_map per sample.
  * §2.5  ``summarize_panel_distribution`` -- per-model 5/50/95% growth-flux
                                     quantiles and P(grows) across an FBA
                                     panel that was rerun for each sample.

Independence assumption: §3.4 of the review asks for the full CC covariance
(``standard_dg_prime_multi``).  That covariance is not on disk -- only the
marginal σ (``deltagerr``) is.  The MC sampler treats reactions as
independent under their marginal normals, which is the best we can do with
the cached MSDB inputs.  This is documented in the driver output.
"""

from __future__ import annotations

import copy
import logging
import math
import multiprocessing as mp
import random
from typing import Iterable, Optional

import reversibility_lib as lib

LOG = logging.getLogger(__name__)

# Floor on σ so a deltagerr of literally 0.0 (no CC uncertainty estimate) does
# not collapse the normal into a delta.  Picked at the order of CC's tightest
# 1-σ for well-measured reactions.
_SIGMA_FLOOR = 0.05  # kcal/mol


# ---------------------------------------------------------------------------
# Analytic P(direction)
# ---------------------------------------------------------------------------
def _phi(z: float) -> float:
    """Standard normal CDF via ``math.erf`` (stdlib, no SciPy dependency)."""
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def _walk_rgt_sum(rxn_entry, cfg: lib.ReversibilityConfig) -> float:
    """Return the ``rgt_sum`` accumulator alone (mMdeltaG concentration term).

    Reuses ``lib._walk_stoichiometry`` -- this gives the same value the cascade
    uses internally, so P-based decisions are consistent with the deterministic
    cascade's `mMdeltaG` step.
    """
    terms = lib._walk_stoichiometry(rxn_entry["stoichiometry"], cfg)
    return cfg.rt * terms["rgt_sum"]


def p_direction(rxn_entry, cfg: Optional[lib.ReversibilityConfig] = None,
                threshold_kcal: float = 0.0) -> Optional[tuple]:
    """Return ``(P_forward, P_reverse, P_reversible)`` for one reaction.

    Defined on the **mMdeltaG** quantity (cascade heuristic 4's quantity), so
    the probabilities are directly comparable to the per-reaction direction
    the cascade lands on:

        mMdeltaG = ΔG′° + RT · Σ_i ν_i · ln(c_i)

    where ``Σ_i ν_i · ln(c_i)`` is computed from the same per-metabolite
    concentration tables the cascade uses (``cfg.cell_conc`` or
    ``cfg.per_met_conc``).  Under the marginal CC normal
    ``ΔG′° ~ N(μ=deltag, σ=deltagerr)`` (treating ``Σ_i …`` as fixed):

        mMdeltaG ~ N(μ + RT·Σ, σ)

    The "forward" event is mMdeltaG < -threshold, "reverse" is
    mMdeltaG > +threshold, "reversible" is everything in between.

    Returns ``None`` for reactions without a usable ΔG′° (EMPTY / Incomplete /
    SENTINEL_DG).
    """
    if cfg is None:
        cfg = lib.ReversibilityConfig()
    if rxn_entry.get("status") == "EMPTY":
        return None
    dg, dge, _ = lib._energy_for(rxn_entry, "EQ")
    if dg is None:
        dg, dge, _ = lib._energy_for(rxn_entry, "GC")
    if dg is None:
        return None
    sigma = max(float(dge or 0.0), _SIGMA_FLOOR)
    try:
        rgt_term = _walk_rgt_sum(rxn_entry, cfg)
    except Exception:
        rgt_term = 0.0
    mu = float(dg) + rgt_term
    z_f = (-threshold_kcal - mu) / sigma
    z_r = (+threshold_kcal - mu) / sigma
    p_forward = _phi(z_f)            # P(mMdeltaG < -threshold)
    p_reverse = 1.0 - _phi(z_r)      # P(mMdeltaG > +threshold)
    p_rev = max(0.0, 1.0 - p_forward - p_reverse)
    return p_forward, p_reverse, p_rev


# ---------------------------------------------------------------------------
# Monte-Carlo cascade replay
# ---------------------------------------------------------------------------
def _snapshot_dg_dge(rxns: dict) -> tuple:
    """Capture the original ``deltag`` / ``deltagerr`` per reaction.

    Returned as two parallel dicts so we can resample ``deltag`` in place
    inside the sampling loop without losing the original mean.
    """
    mu = {}
    sigma = {}
    for rid, entry in rxns.items():
        dg = entry.get("deltag")
        dge = entry.get("deltagerr")
        try:
            if dg is None:
                continue
            dg_f = float(dg)
            if dg_f == lib.SENTINEL_DG:
                continue
            dge_f = max(float(dge or 0.0), _SIGMA_FLOOR)
            mu[rid] = dg_f
            sigma[rid] = dge_f
        except (TypeError, ValueError):
            continue
    return mu, sigma


def sample_cascade(rxns: dict,
                   cfg: Optional[lib.ReversibilityConfig] = None,
                   n_samples: int = 50,
                   db_level: str = "EQ",
                   seed: int = 0,
                   verbose: bool = False) -> list:
    """Run the cascade ``n_samples`` times with resampled ΔG′°.

    Mutates ``rxns`` *in place* during the loop and restores the original
    values at the end -- avoids deep-copying the 56K-reaction dict per
    sample.  Each sample's GC pre-roll mutates ``rxn_entry['reversibility']``,
    which is overwritten on the next sample anyway, so no special handling.

    Returns ``[{rxn_id: rev}, ...]`` (one rev_map per sample).
    """
    if cfg is None:
        cfg = lib.ReversibilityConfig()
    rng = random.Random(seed)
    mu, sigma = _snapshot_dg_dge(rxns)
    original_dg = {rid: rxns[rid].get("deltag") for rid in mu}
    samples = []
    for i in range(n_samples):
        for rid, m in mu.items():
            rxns[rid]["deltag"] = rng.gauss(m, sigma[rid])
        out = lib.run_cascade(rxns, db_level=db_level, cfg=cfg,
                              gc_first=(db_level == "EQ"))
        samples.append({r: rev for r, (_, rev) in out.items()})
        if verbose:
            from collections import Counter
            c = Counter(samples[-1].values())
            LOG.info("sample %3d: %s", i, dict(c))
    # Restore originals so the caller's dict is unchanged.
    for rid, val in original_dg.items():
        rxns[rid]["deltag"] = val
    return samples


# ---------------------------------------------------------------------------
# Panel-FBA distribution + summarization
# ---------------------------------------------------------------------------
def panel_fba_for_samples(panel_ids, rev_map_samples, n_workers: int = 4,
                          verbose: bool = False) -> list:
    """Run ``gh.run_panel`` once per sample.  Returns ``[panel_fba, ...]``.

    Spawning a fresh ``multiprocessing.Pool`` per sample has a small setup
    cost but keeps memory bounded.  For N=100, this is ~150 s of overhead
    on top of the FBA solves themselves.
    """
    import growth_heuristics as gh  # local import to keep this module light
    out = []
    for i, rmap in enumerate(rev_map_samples):
        if verbose:
            LOG.info("panel-FBA sample %d/%d", i + 1, len(rev_map_samples))
        out.append(gh.run_panel(panel_ids, reversibility_map=rmap,
                                baseline_map=None, n_workers=n_workers))
    return out


def summarize_panel_distribution(panel_fba_runs: list,
                                 quantiles: tuple = (0.05, 0.5, 0.95)) -> dict:
    """For each panel model, compute flux quantiles + P(grows) across samples.

    Returns ``{model_id: {q05, q50, q95, mean_flux, p_grows, n_samples}}``.
    Quantiles are computed with the standard "nearest-rank" definition
    (no NumPy dependency).
    """
    by_model: dict = {}
    for run in panel_fba_runs:
        for r in run:
            mid = r["model_id"]
            by_model.setdefault(mid, []).append(r)

    def _quantile(sorted_xs, q):
        if not sorted_xs:
            return 0.0
        idx = max(0, min(len(sorted_xs) - 1, int(round(q * (len(sorted_xs) - 1)))))
        return sorted_xs[idx]

    out = {}
    for mid, runs in by_model.items():
        fluxes = sorted(float(r["growth_flux"]) for r in runs)
        grows = sum(1 for r in runs if r["grows"])
        stats = {
            "n_samples":     len(runs),
            "mean_flux":     sum(fluxes) / len(fluxes),
            "p_grows":       grows / len(runs),
        }
        for q in quantiles:
            stats[f"q{int(round(q * 100)):02d}"] = _quantile(fluxes, q)
        out[mid] = stats
    return out


# ---------------------------------------------------------------------------
# Convenience: aggregate per-reaction P(direction) across the database
# ---------------------------------------------------------------------------
def aggregate_p_direction(rxns: dict,
                          cfg: Optional[lib.ReversibilityConfig] = None,
                          threshold_kcal: float = 0.0,
                          restrict_to: Optional[Iterable[str]] = None) -> dict:
    """Compute ``p_direction`` for every reaction (or a subset).

    Returns ``{rxn_id: (p_fwd, p_rev, p_rev)}`` for reactions with usable ΔG′°,
    plus an aggregate count under three labels (`confident forward`,
    `confident reverse`, `ambiguous`) using a 0.95 threshold on the max prob.
    """
    if cfg is None:
        cfg = lib.ReversibilityConfig()
    iter_ids = restrict_to if restrict_to is not None else rxns
    per_rxn = {}
    n_confident_fwd = n_confident_rev = n_ambiguous = 0
    for rid in iter_ids:
        entry = rxns.get(rid)
        if entry is None:
            continue
        p = p_direction(entry, cfg, threshold_kcal)
        if p is None:
            continue
        per_rxn[rid] = p
        pf, pr, _ = p
        if pf >= 0.95:
            n_confident_fwd += 1
        elif pr >= 0.95:
            n_confident_rev += 1
        else:
            n_ambiguous += 1
    return {
        "per_rxn": per_rxn,
        "counts": {
            "confident_forward":  n_confident_fwd,
            "confident_reverse":  n_confident_rev,
            "ambiguous":          n_ambiguous,
            "total_with_dg":      len(per_rxn),
        },
    }
