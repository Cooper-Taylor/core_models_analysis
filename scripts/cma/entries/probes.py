"""Per-model probes.

Each probe supplies only the measurement; the runner supplies the prologue
(load, rebind bounds from a direction map, apply media, set the objective) that
is currently hand-written at least thirteen times across the drivers.

The seven existing probes wrap ``growth_heuristics``'s ``*_one`` functions
verbatim, so behaviour is unchanged while the drivers shed their scaffolding.
``internal_pool=True`` marks a probe that parallelizes internally -- cobra's FVA
does -- so the runner refuses to nest it inside its own pool.
"""

from __future__ import annotations

from ..kinds import Probe
from ..probes import register


def _gh():
    import growth_heuristics as gh

    return gh


register(Probe(
    name="growth",
    fn=lambda model_id, ctx: _gh().fba_one(
        model_id, reversibility_map=ctx.get("direction_map"),
        baseline_map=ctx.get("baseline_map"), ignore_bounds=ctx.get("ignore_bounds", False)),
    out_template="site/data/{ns}panel_growth.json",
    description="Biomass FBA under a direction map. The default sweep.",
))

register(Probe(
    name="rxn_pipeline",
    fn=lambda model_id, ctx: _gh().rxn_pipeline_one(
        model_id, ctx["baseline_map"], ctx["changed_dirs"]),
    out_template="site/data/{ns}panel_rxn_pipeline.json",
    description="Per-reaction decomposition of a variant's growth effect.",
))

register(Probe(
    name="key_reactions",
    fn=lambda model_id, ctx: _gh().key_reactions_one(
        model_id, ctx["baseline_map"], top_n=ctx.get("top_n", 75)),
    out_template="site/data/{ns}panel_key_reactions.json",
    description="Reactions whose direction most changes growth.",
))

register(Probe(
    name="growth_control",
    fn=lambda model_id, ctx: _gh().growth_control_one(
        model_id, ctx["baseline_map"], top_n=ctx.get("top_n", 90)),
    out_template="site/data/{ns}panel_growth_control.json",
    description="Growth-control coefficients per reaction direction.",
))

register(Probe(
    name="synthetic_lethal",
    fn=lambda model_id, ctx: _gh().synthetic_lethal_one(model_id, ctx["baseline_map"]),
    out_template="site/data/{ns}panel_synthetic_lethal.json",
    description="Direction pairs that are jointly but not singly lethal.",
))

register(Probe(
    name="fva",
    fn=lambda model_id, ctx: _gh().fva_one(
        model_id, ctx["baseline_map"], fraction=ctx.get("fraction", 0.99)),
    internal_pool=True,          # cobra's FVA runs its own pool; do not nest
    out_template="site/data/{ns}panel_fva.json",
    description="Flux variability at a fraction of optimal growth.",
))

register(Probe(
    name="flux_loops",
    fn=lambda model_id, ctx: _gh().flux_loops_one(
        model_id, reversibility_map=ctx.get("direction_map")),
    out_template="site/data/{ns}panel_flux_loops.json",
    description=("Energy-generating cycles. NOTE: relies on a small bio* reaction "
                 "being left open as the ATP probe; a model set whose BiomassPolicy "
                 "declares atp_probe=None and has only one large biomass needs an "
                 "injected probe or this silently reports zero loops."),
))
