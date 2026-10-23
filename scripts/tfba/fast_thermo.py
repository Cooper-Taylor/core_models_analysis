"""Replacement for SimpleThermoPkg.build_package that skips the per-reaction LPs.

The shipped code calls model.optimize() twice per reaction purely to learn the
range of sum(s_ij * p_j). Because the potential variables appear in no other
constraint at that point and are box-bounded at [min_potential, max_potential],
that range is just

    max = sum_j (s_j > 0 ? s_j*hi : s_j*lo)        min = -max   (symmetric bounds)

so max_energy_magnitude = max(|min|, |max|). Computing it directly is the same
number without 2*N mixed-integer solves. Equivalence is asserted against the
shipped implementation in the calling script.
"""
from modelseedpy.fbapkg.simplethermopkg import SimpleThermoPkg


def fast_build_package(self, parameters):
    self.validate_parameters(parameters, [], {
        "filter": None, "min_potential": 0, "max_potential": 1000,
        "dgbin": False, "reduced_constraints": False,
    })
    self.pkgmgr.getpkg("RevBinPkg").build_package(self.parameters["filter"])
    for metabolite in self.model.metabolites:
        self.build_variable(metabolite)

    lo = self.parameters["min_potential"]
    hi = self.parameters["max_potential"]
    for reaction in self.model.reactions:
        if reaction.id[:3] in ("EX_", "SK_", "DM_"):
            continue
        hi_sum = sum(s * (hi if s > 0 else lo) for s in reaction.metabolites.values())
        lo_sum = sum(s * (lo if s > 0 else hi) for s in reaction.metabolites.values())
        max_energy_magnitude = max(abs(hi_sum), abs(lo_sum))
        if self.parameters["filter"] is None or reaction.id in self.parameters["filter"]:
            self.build_constraint(reaction, max_energy_magnitude)

    if self.parameters["dgbin"]:
        self.optimize_dgbin()


def install():
    SimpleThermoPkg.build_package = fast_build_package
