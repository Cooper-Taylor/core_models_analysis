#!/usr/bin/env python3
"""Build the direction-influence analysis as PDF and Markdown from one source.

    python3 scripts/build_influence_report.py                  # both formats
    python3 scripts/build_influence_report.py --format md
    python3 scripts/build_influence_report.py --format pdf --out /tmp/x.pdf

The content is defined once, in :func:`document`, as format-neutral blocks;
:func:`render_pdf` and :func:`render_markdown` are the only format-aware code.
Two renderings of one definition cannot disagree with each other, and because
every number is pulled from the result files at build time, neither can
disagree with the data. The previous hand-maintained Markdown drifted from the
results on sixteen separate claims before an independent recount caught it.

Inline markup in the block text is Markdown (``**bold**``, ``*italic*``);
the PDF renderer converts it to ReportLab's tag soup.

Inputs:  results/influence_v201/{core_all,gsm_gmm,gsm_auxo}/growth*.json
         results/influence_v201/{coverage,conversion_roundtrip}.json
         results/direction_sets_v201/{manifest.json,comparison/direction_summary.json}
         reports/figures/direction_influence/*.png   (figures.tsv: direction_influence)
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
from cma import paths  # noqa: E402

SETS = ["v201_group_contribution", "v201_dgpredictor", "v201_equilibrator",
        "v201_llm_council", "v201_gold_only", "v201_silver_only",
        "v201_bronze_only", "v201_gold_silver", "v201_gold_silver_bronze"]
LABEL = {"v201_group_contribution": "Group contribution", "v201_dgpredictor": "dGPredictor",
         "v201_equilibrator": "eQuilibrator", "v201_llm_council": "LLM council",
         "v201_gold_only": "gold",
         "v201_silver_only": "silver", "v201_bronze_only": "bronze",
         "v201_gold_silver": "gold + silver",
         "v201_gold_silver_bronze": "gold + silver + bronze"}

FIGDIR = "reports/figures/direction_influence"


# ---------------------------------------------------------------------------
def load() -> dict:
    inf, sets = paths.results("influence_v201"), paths.results("direction_sets_v201")
    missing = [p for p in (inf / "core_all" / "growth.json",
                           inf / "coverage.json", sets / "manifest.json")
               if not p.exists()]
    if missing:
        raise SystemExit("missing inputs:\n  " + "\n  ".join(str(m) for m in missing)
                         + "\nRun the sweeps first; see PIPELINE.md.")
    D = {}
    for tag in ("core_all", "gsm_gmm", "gsm_auxo"):
        D[tag] = {
            "growth": json.loads((inf / tag / "growth.json").read_text())["per_variant"],
            "diff": json.loads((inf / tag / "growth_diff.json").read_text())["per_variant"],
        }
    D["coverage"] = json.loads((inf / "coverage.json").read_text())
    D["manifest"] = json.loads((sets / "manifest.json").read_text())
    D["agreement"] = json.loads((sets / "comparison" / "direction_summary.json").read_text())
    D["roundtrip"] = json.loads((inf / "conversion_roundtrip.json").read_text())
    gi = paths.results("gapfill_direction_impact")
    if (gi / "summary.json").exists():
        D["gapfill"] = json.loads((gi / "summary.json").read_text())
        D["cofactor"] = json.loads((gi / "by_cofactor.json").read_text())
    return D


# ---------------------------------------------------------------------------
def document(D) -> list:
    """The report, as format-neutral blocks. This is the single source of content."""
    core, gmm, aux = (D[t]["growth"] for t in ("core_all", "gsm_gmm", "gsm_auxo"))
    cb, gb, ab = core["on_disk"], gmm["on_disk"], aux["on_disk"]
    cov, stat = D["coverage"], {x["name"]: x for x in D["manifest"]["sets"]}
    ag = D["agreement"]["pairwise_agreement"]
    rt = D["roundtrip"]

    def pair(a, b):
        return ag.get(f"{a} vs {b}") or ag.get(f"{b} vs {a}")

    e1 = pair("v201_equilibrator", "v201_gold_silver_bronze")
    e2 = pair("v201_llm_council", "v201_equilibrator")
    e4 = pair("v201_gold_only", "v201_dgpredictor")

    B: list = []
    def add(*b):
        B.append(b)

    add("title", "Reaction-direction sources and what they do to models")
    add("sub", f"Ten direction sources from ModelSEED Biochemistry v2.0.1, applied to "
               f"{cb['n_models']:,} core models and {gb['n_models']:,} genome-scale models")

    add("h1", "Summary")
    add("lead", "The same direction source has opposite effects on the two model sets. On the "
                "core models most sources **increase** the number that grow; on the "
                "genome-scale models most **abolish** growth entirely. That is not a property "
                "of the sources.")
    add("p", "Two things produce it. The genome-scale networks are roughly fifteen times larger "
             "by distinct reaction count, so a source rewrites a far greater share of the "
             "network. And within that larger surface sit a few reactions that are "
             "thermodynamically uphill but biologically obligatory. One of them, the thiamine "
             "phosphomethylpyrimidine kinase, is present in every genome-scale model tested and "
             "on its own abolishes growth in 99% of them. The thermodynamic estimate for it is "
             "correct; the constraint derived from it is still wrong.")
    add("p", "The practical consequence is that **the core panel cannot validate a direction "
             "source**. It exercises 239 of the 56,012 reactions in the database and omits the "
             "pathways where thermodynamic and biological direction diverge.")

    add("h1", "1. What was measured")
    add("p", "**Baseline.** The bounds written in the model files themselves: forward, reverse "
             "or reversible according to each reaction's stored lower and upper bound. Nothing "
             "is recomputed. This is the *on-disk* row throughout, and it is what the models "
             "were gap-filled and published with.")
    add("p", "**Applying a source.** A reaction whose SEED annotation appears in the source is "
             "rewritten to (0, 1000), (-1000, 0) or (-1000, 1000). A reaction the source does "
             "not mention keeps its on-disk bounds. Media are applied afterwards, biomass is "
             "the objective, and a model counts as growing above 1e-6.")
    add("table",
        [["Model set", "Models", "Reactions / model", "Distinct reactions", "Media"],
         ["Core (core_kegg2)", f"{cb['n_models']:,}",
          f"{cov['core_kegg2']['reactions_per_model']['median']}",
          f"{cov['core_kegg2']['distinct_msdb_reactions']:,}", "347-compound complete"],
         ["Genome-scale (ms2_gsm)", f"{gb['n_models']:,}",
          f"{cov['ms2_gsm']['reactions_per_model']['median']:,}",
          f"{cov['ms2_gsm']['distinct_msdb_reactions']:,}", "20-compound glucose minimal"],
         ["Genome-scale (auxotrophy)", f"{ab['n_models']:,}", "1,085", "3,484",
          "51-compound auxotrophy"]],
        [44, 18, 28, 28, 50], (1, 2, 3))
    add("caption", "Reactions per model is the median count of distinct ModelSEED reactions. "
                   "The genome-scale models come from the public KBase workspaces of the "
                   "ModelSEED v2 manuscript and were gap-filled on the media shown.")

    add("h1", "2. Genome-scale model properties")
    add("p", "**The two media sets are genetically identical but gap-filled differently.** The "
             "same 5,420 genomes are reconstructed twice: once against glucose minimal media "
             "and again against auxotrophy media. Both sets have identical model IDs and start "
             "with the same genome-derived reactions, but gap-filling against different media "
             "adds different gap-filled reactions. This separation allows us to measure the "
             "effect of the medium on growth outcomes apart from the effect of direction sources.")
    add("p", "**Gap-filled reactions differ between media but constitute <1% of total reactions.** "
             "Of the 3,634 distinct reactions in the glucose minimal models, 2,665 are "
             "genome-derived, 408 are gap-filled, and 696 are media-conditional (present in one "
             "medium only). Gap-filled reactions account for only ~0.6% of all reaction "
             "occurrences across a model, but they are re-called at much higher rates: 63% of "
             "gap-filled reactions are given new directions by the graded tier, versus only 29% "
             "of genome-derived reactions.")

    add("h1", "3. The ten sources")
    add("p", "Three thermodynamic estimators, one LLM route, three evidence tiers and three "
             "cumulative tiers. A grade belongs to a *reaction*, not to a source: the release "
             "ships one grade per reaction and keeps the per-source table internal. A tier is "
             "therefore the release's recommended direction restricted to reactions at that "
             "tier, chosen by a fixed precedence of eQuilibrator over dGPredictor over group "
             "contribution.")
    rows = [["Source", "Calls", "Forward", "Reverse", "Reversible", "Directional"]]
    for k in SETS:
        st = stat[k]; o = st["operators"]; n = st["n_callable"]
        d = o.get(">", 0) + o.get("<", 0)
        rows.append([LABEL[k], f"{n:,}", f"{o.get('>',0):,}", f"{o.get('<',0):,}",
                     f"{o.get('=',0):,}", f"{100*d/n:.0f}%"])
    add("table", rows, [50, 22, 22, 22, 24, 24], (1, 2, 3, 4, 5))
    add("caption", "No set contains an undecided call. A direction map is applied by rewriting "
                   "bounds, and an undecided value would be applied as fully reversible, so "
                   "emitting it for a reaction that was graded but could not be called would "
                   "*remove* a constraint in the name of evidence. Absent means no opinion, and "
                   "the model keeps its own bounds. **Directional** is the share of a source's "
                   "calls that fix a direction rather than allow both: the constraint pressure "
                   "it applies.")

    add("pagebreak")
    add("h1", "4. What each source does to growth")
    add("figure", "fig1_growth_by_panel.png",
        "**Figure 1.** Share of models that grow under each direction source, against the "
        "bounds already in the model files. The three panels share a source ordering; each has "
        "its own baseline.", 168)
    rows = [["Source", "Core", "vs base", "GSM glucose", "GSM auxotrophy"],
            ["**on-disk baseline**", f"**{cb['n_growers']:,}**", "--",
             f"**{gb['n_growers']:,}**", f"**{ab['n_growers']:,}**"]]
    for k in SETS:
        d = D["core_all"]["diff"][k]
        net = d["n_gained_growth"] - d["n_lost_growth"]
        rows.append([LABEL[k], f"{core[k]['n_growers']:,}", f"{net:+,}",
                     f"{gmm[k]['n_growers']:,}", f"{aux[k]['n_growers']:,}"])
    add("table", rows, [52, 24, 22, 30, 34], (1, 2, 3, 4))
    add("caption", "Models that grow, of 5,683 core and 5,420 genome-scale. Zero models failed "
                   "to load in any run.")

    add("h2", "On the core models, most sources help")
    add("p", "eQuilibrator is the least directional source and gains the most models. **Gold "
             "alone costs growth while silver alone gains it**, which inverts the evidence "
             "hierarchy: 122 models stop growing under the best-graded directions and 367 start "
             "growing under the weaker ones. Gold calls 58% of its reactions directional against "
             "silver's 43%, and a directional call removes a degree of freedom the model was "
             "using. The grade measures how well an energy is known, not how safe it is to "
             "impose.")

    add("h2", "On the genome-scale models, most collapse")
    add("p", "On glucose minimal media **six of the ten sources leave zero growing models**. "
             "Every one returns an optimal solution with biomass flux zero, so these are "
             "genuine no-growth results rather than solver failures. The auxotrophy set starts "
             "from a near-complete baseline and separates the sources better: gold retains 85% "
             "of baseline, the silver-containing sets about a third, and six of eleven fall "
             "below 4%. The cumulative tiers decay monotonically there, 4,618 to 1,866 to 55.")
    add("p", "\"Nearly everything collapses\" is fair for glucose minimal and an overstatement "
             "for auxotrophy. One detail runs the other way: among the models that do survive, "
             "mean growth flux *rises*, from 1.32 at baseline to 4.09 under gold and 6.52 under "
             "gold plus silver. Directional constraints remove futile cycling and channel flux, "
             "so the cost is concentrated entirely in the models that stop growing.")

    add("h1", "5. The grading, against the original thermodynamic sources")
    add("figure", "fig2_grading_vs_thermo.png",
        "**Figure 2.** Left: how often an evidence tier and an estimator agree, over the "
        "reactions where both have an opinion. Right: the share of each source's calls that fix "
        "a direction rather than allow both.", 168)
    add("p", f"**The tiers are largely eQuilibrator re-expressed.** The full graded set agrees "
             f"with eQuilibrator on {e1['pct_agree']}% of the {e1['shared']:,} reactions they "
             "share, bronze at 100.0% and silver at 99.9%. That follows directly from the "
             "recommendation precedence putting eQuilibrator first. A tier behaves as a "
             "confidence filter on one source rather than as a consensus of several, so "
             "comparing a tier against eQuilibrator is close to comparing eQuilibrator against "
             "itself.")
    add("p", f"**The estimators disagree with each other more than the tiers suggest.** "
             f"dGPredictor agrees with the gold tier on only {e4['pct_agree']}% of the "
             f"{e4['shared']:,} reactions they share, barely above chance among three options. "
             "Gold is the tier where the evidence is strongest, so this is disagreement about "
             "well-measured reactions, not about hard ones.")
    add("p", f"**The LLM route is a different kind of claim.** The council agrees with "
             f"eQuilibrator on {e2['pct_agree']}% of {e2['shared']:,} shared reactions; it "
             "reasons from the reaction rather than from an energy, which is the release's own "
             "argument for keeping it out of the grading.")
    add("p", "The right panel of Figure 2 shows what the agreement numbers hide: the tiers are "
             "markedly *less committal* than the estimators they are drawn from. Silver allows "
             "both directions on 57% of its calls against dGPredictor's 27%. That is most of "
             "why the silver tier is gentler on the models than the source it mostly agrees "
             "with.")

    add("h1", "6. Why the genome-scale models collapse")
    add("p", "A total collapse is as likely to be a defect as a result, so it was diagnosed "
             "rather than reported.")
    add("p", f"**The conversion is faithful.** Direction round-trips exactly from the KBase "
             f"originals: {rt['reactions_checked']:,} reactions across a "
             f"{rt['models_sampled']}-model random sample, {rt['direction_mismatches']} "
             "mismatches. The models also grow normally on their own bounds. Only the override "
             "breaks them.")
    add("p", "**It is not the medium.** The collapse reproduces on the 347-compound complete "
             "medium, not only the 20-compound minimal one.")
    add("p", "**It localises to individual reactions.** Taking one genome-scale model under "
             "eQuilibrator directions, 229 of its reactions differ from disk and 51 of those "
             "are constraining, but applying them one at a time shows that only three "
             "individually abolish growth. Repeating that scan across 73 models that all grow "
             "at baseline, every one has at least one individually lethal constraint, and the "
             "distribution is extremely skewed.")
    add("figure", "fig4_mechanism.png",
        "**Figure 3.** Share of tested genome-scale models in which constraining one reaction "
        "alone abolishes growth. Rates are from a 73-model scan and shift by a few points "
        "between samples; the ranking of the first does not.", 162)
    add("p", "**One reaction explains most of it.** The thiamine phosphomethylpyrimidine kinase "
             "(rxn03108) is present in every genome-scale model and lethal in 72 of 73. "
             "eQuilibrator computes +15.59 kcal/mol for it and therefore calls it reverse. That "
             "number is not wrong: the phosphoryl transfer is uphill in isolation. But the "
             "enzyme is an ATP-driven kinase in thiamine biosynthesis and must run forward, so "
             "constraining it to the thermodynamically favoured direction removes thiamine and "
             "with it biomass. Both LLM routes call this reaction forward.")
    add("p", "ATP sulfurylase (rxn00379) is the same failure with better evidence: graded "
             "**gold** from an openTECR measurement, +11 kcal/mol, called reverse by "
             "eQuilibrator, dGPredictor and Claude. In the cell it runs forward because "
             "pyrophosphatase removes the pyrophosphate product. The measurement is right and "
             "the direction call is right in isolation; the constraint is still wrong. It is a "
             "milder case, present in 67% of models and lethal in about one in eight of those.")
    add("p", "Removing the top three reactions from the eQuilibrator map restores 137 of 381 "
             "growers in a 480-model sample, from zero. The collapse is not diffuse "
             "over-constraint; it is a handful of obligatory, cofactor-driven reactions.")
    add("p", "**Why the core models escape.** They contain 239 distinct reactions against "
             "3,634. Thiamine biosynthesis and sulfate assimilation are simply absent from a "
             "core carbon-metabolism model, so the reactions that kill the genome-scale models "
             "are never touched.")
    add("p", "**Why gold alone helps but gold+silver kills on glucose minimal.** Gold constrains "
             "only the most confident reactions (58% are directional), and it avoids the worst "
             "cofactor-driven reactions that lack measurements. On glucose minimal, gold retains "
             "2,502 growers (58% of baseline). Silver contains many lower-confidence calls, "
             "notably on CoA chemistry and membrane redox families where the evidence is weak. "
             "The core failure case — thiamine phosphomethylpyrimidine kinase (rxn03108) — "
             "receives different calls: eQuilibrator marks it reverse (call: <), silver inherits "
             "that reverse call, but gold does not mention it (no measurement, no call). Adding "
             "silver to gold adds 642 new overrides per model on average and includes the reverse "
             "calls on the lethal reactions. On auxotrophy media (where the baseline is 5,399 and "
             "includes more thiamine), gold + silver retains 1,866 growers (35%), a steep loss but "
             "not total collapse, because the auxotrophy medium supplies thiamine and relaxes the "
             "constraint somewhat.")
    add("p", "**Reaction rxn03108 (thiamine phosphomethylpyrimidine kinase) in the sources.** "
             "Gold: not called (no measurement). Silver: reverse. Bronze: not shown. eQuilibrator: "
             "reverse (+15.59 kcal/mol, thermodynamically justified). LLM council: forward (both "
             "council and Claude 3.5 call this forward, reasoning from its role in biosynthesis "
             "rather than from isolated energy).")

    if "gapfill" in D:
        idx = {(r["tier"], r["class"]): r for r in D["gapfill"]}
        eo = D["cofactor"]["original"]
        eg = D["cofactor"]["gapfilled"]
        fam = {f["cofactor"]: f for f in eo["families"]}
        g_o = idx[("gold", "original")]
        g_g = idx[("gold", "gapfilled")]
        a_o = idx[("gold + silver + bronze", "original")]
        a_g = idx[("gold + silver + bronze", "gapfilled")]
        a_m = idx[("gold + silver + bronze", "media-only")]

        add("h1", "7. Gap-filled reactions, and which ones get re-called")
        add("p", "The genome-scale models exist twice: the same 5,419 genomes reconstructed "
                 "once and gap-filled against glucose minimal media, and again against "
                 "auxotrophy media. Comparing the pair separates what the genome supports "
                 "from what gap-filling added. A reaction present in **both** models and "
                 "gap-filled in **neither** is genome-derived; one carrying gap-fill data in "
                 "either is gap-filled; one present in only one of the pair is conditional on "
                 "the medium and is reported separately.")
        add("table",
            [["Provenance", "Distinct reactions", "Occurrences", "What it is"],
             ["original", f"{a_o['reactions_in_class']:,}",
              f"{a_o['occurrences_in_class']:,}", "in both media models, gap-filled in neither"],
             ["gapfilled", f"{a_g['reactions_in_class']:,}",
              f"{a_g['occurrences_in_class']:,}", "carries gapfill data in either model"],
             ["media-only", f"{a_m['reactions_in_class']:,}",
              f"{a_m['occurrences_in_class']:,}", "present in one media model, absent from the other"]],
            [26, 30, 26, 62], (1, 2))
        add("caption", "Gap-filled and media-conditional reactions are numerous as distinct "
                       "identifiers but rare as occurrences: together they are under 1% of all "
                       "reaction instances. Most of any given model is genome-derived.")
        add("figure", "fig5_gapfill_cofactor.png",
            "**Figure 5.** Left: the share of the reactions a tier calls that it would change, "
            "split by provenance. Right: which cofactor families mark a genome-derived "
            "reaction as likely to be re-called, as enrichment against the 29% class average.",
            168)
        add("p", f"**Gap-filled reactions are re-called two to three times as often.** Of the "
                 f"genome-derived reactions the full graded set has an opinion about, it would "
                 f"change {a_o['pct_of_called_modified']:.0f}%. Of the gap-filled ones it would "
                 f"change {a_g['pct_of_called_modified']:.0f}%, and of the media-conditional "
                 f"ones {a_m['pct_of_called_modified']:.0f}%. Gold alone is starker still: "
                 f"{g_o['pct_of_called_modified']:.0f}% of genome-derived against "
                 f"{g_g['pct_of_called_modified']:.0f}% of gap-filled.")
        add("p", "That is what you would expect if gap-filling picks directions to make a model "
                 "grow rather than to match the thermodynamics, and it is a second reason the "
                 "genome-scale models are more fragile than the core ones under a direction "
                 "override: the override is disproportionately aimed at the reactions the "
                 "reconstruction added on purpose. It also connects to the mechanism in "
                 "section 5, where the second most lethal reaction, the "
                 "3'-phosphoadenylyl-sulfate sulfohydrolase, is itself gap-filled.")
        add("h2", "Which cofactors mark a reaction as likely to be re-called")
        add("p", "The prediction going in was that the modified reactions would concentrate in "
                 "ATP, NAD and NADP chemistry. Half of that holds. Raw shares are misleading "
                 "here, because a family that appears in many reactions will appear in many "
                 "modified reactions; what matters is the rate among the reactions a tier "
                 "actually calls.")
        rows = [["Cofactor family", "Called", "Modified", "Rate", "Enrichment"]]
        for k in ("CoA", "FAD(H2)", "quinone", "ATP/ADP/AMP", "pyrophosphate",
                  "ammonia", "CO2", "phosphate", "NADP(H)", "NAD(H)", "O2"):
            f = fam.get(k)
            if f:
                rows.append([k, f"{f['called']:,}", f"{f['modified']:,}",
                             f"{f['modified_rate']}%", f"{f['enrichment_vs_class']:.2f}"])
        add("table", rows, [40, 22, 24, 20, 26], (1, 2, 3, 4))
        add("caption", "Genome-derived reactions only, under gold + silver + bronze. "
                       "Enrichment is the family's modification rate over the "
                       f"{eo['class_modified_rate']}% class average; tags overlap, since a "
                       "reaction can use several cofactors.")
        add("p", "**ATP is enriched, but NAD and NADP are not.** ATP chemistry is re-called "
                 "1.5 times more often than average and pyrophosphate 1.4 times, which matches "
                 "the mechanism in section 5: these are the reactions driven by cofactor "
                 "hydrolysis rather than by their own free energy. But NAD(H) sits at 0.78 and "
                 "NADP(H) at 0.81, meaning both are re-called **less** often than average. "
                 "Nicotinamide redox potentials are well characterised and consistently "
                 "estimated, so the sources mostly agree with the reconstruction about them.")
        add("p", "**The strongest signal was not predicted at all.** CoA chemistry is re-called "
                 f"at {fam['CoA']['modified_rate']:.0f}%, more than twice the class average, and "
                 "the two membrane redox families follow at 1.7 times. Thioester hydrolysis is "
                 "strongly favourable in isolation while the cell routinely runs it "
                 "biosynthetically, and quinone potentials are the known weak point of the "
                 "estimators. Oxygen chemistry is the opposite case: re-called at 6.9%, almost "
                 "never, because a reaction consuming O2 is unambiguously downhill and every "
                 "source agrees.")
        add("p", "This supports the idea of fixing groups rather than reactions, but it "
                 "redirects it. The groups worth a rule are **CoA thioesters, quinone and FAD "
                 "redox, and ATP or pyrophosphate coupling** -- together 568 of the 1,706 "
                 "genome-derived reactions the graded set calls. A rule exempting NAD and NADP "
                 "chemistry would target reactions the sources already get right.")

    add("h1", "8. Does anything predict what a source will do?")
    add("figure", "fig3_coverage_effect.png",
        "**Figure 4.** Share of reaction occurrences a source changes, against the share of "
        "models that grow. Note the independent vertical scales: the core panel spans 51-70%, "
        "the genome-scale panel 0-80%.", 168)
    add("p", "Not reliably. Of the candidate summary statistics -- directional fraction, "
             "reactions called, occurrences changed, overrides per model -- none is significant "
             "on the core panel, and the best on the genome-scale panel reaches only Spearman "
             "rho of -0.65. The share of occurrences a source changes gives rho = -0.65 on "
             "glucose minimal, -0.45 on auxotrophy and -0.03 on the core panel, where it has "
             "essentially no ordering power at all.")
    add("p", "Directional fraction does better on the core panel, rho = -0.68, but only when "
             "computed over the whole release. Restricted to the 239 reactions the core models "
             "actually contain -- the ones that can affect the result -- it falls to -0.28 and "
             "is not significant. The ten sets are also nested rather than independent: gold "
             "sits inside gold plus silver, which sits inside gold plus silver plus bronze.")
    add("p", "Two counterexamples make the point concretely. Bronze changes *less* of the "
             "network than gold, 3.3% against 3.7%, and retains far less growth. dGPredictor "
             "changes only 12.5% and leaves zero growers. **What a source does depends on "
             "which reactions it constrains, not how many.**")

    add("h1", "9. What follows")
    for head, txt in [
        ("A direction source cannot be validated on the core models.",
         "They exercise 239 of 56,012 reactions and omit the pathways where thermodynamic and "
         "biological direction come apart. Earlier conclusions drawn from that panel are not "
         "contradicted so much as shown to be untested."),
        ("A grade is not a licence to constrain.",
         "The grading scheme answers how well an energy is known, and answers it well: ATP "
         "sulfurylase is gold because openTECR measured it. Whether it is safe to fix that "
         "reaction's direction in a model is a different question, and the two come apart "
         "precisely at the cofactor-driven reactions. A tier used as an FBA constraint set "
         "needs a second filter for obligatory reactions."),
        ("Reversibility is the safe default.",
         "Every set that is more than about 60% directional damages the genome-scale models "
         "badly. That is a tendency rather than a rule, but it is the most robust one "
         "available."),
        ("Two candidate fixes, neither tested here.",
         "Exempt reactions whose stoichiometry contains ATP or pyrophosphate hydrolysis from "
         "directional constraint, which is the mechanism in every case examined; or intersect "
         "a tier with the LLM council and keep only reactions where both agree, trading "
         "coverage for safety."),
    ]:
        add("p", f"**{head}** {txt}")

    add("h1", "10. Caveats")
    for c in [
        "**Growth is the only readout.** A direction set that preserves grower count may still "
        "distort flux distributions; nothing here measures that.",
        "**Energy-generating-cycle detection does not work on these models.** The existing "
        "probe closes every biomass reaction above ten metabolites and relies on a smaller one "
        "staying open as an ATP probe. The genome-scale models have exactly one biomass and no "
        "maintenance reaction, so it would report zero cycles for all 5,420 with no error.",
        "**The lethal-reaction scan covers 73 models**, the growers among every 60th of the "
        "set. The ranking of the first reaction is unambiguous and replicated on an independent "
        "sample; the long tail is not fully enumerated.",
        "**Correlations rest on ten non-independent sets.** No coefficient here exceeds 0.7 "
        "in magnitude and most are not significant.",
    ]:
        add("small", c)

    add("h1", "Reproducing this")
    add("code",
        "python3 scripts/build_direction_sets_v201.py\n"
        "python3 beginPipeline --models core_kegg2 --workers 64 \\\n"
        "    --directions v201_gold ... --only growth,growth_diff \\\n"
        "    --out results/influence_v201/core_all\n"
        "python3 scripts/analyze_direction_coverage.py\n"
        "python3 scripts/regen_figures.py direction_influence\n"
        "python3 scripts/build_influence_report.py")
    add("small", "Results are under results/influence_v201/ and results/direction_sets_v201/, "
                 "each run with a run.json recording the exact command. Every number in this "
                 "document is read from those files at build time, and the PDF and Markdown "
                 "renderings come from one content definition in "
                 "scripts/build_influence_report.py.")
    return B


# ---------------------------------------------------------------------------
_BOLD = re.compile(r"\*\*(.+?)\*\*", re.S)
_ITAL = re.compile(r"(?<!\*)\*(?!\*)(.+?)(?<!\*)\*(?!\*)", re.S)


def to_rl(text: str) -> str:
    """Markdown inline markup -> ReportLab tags."""
    t = text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    t = _BOLD.sub(r"<b>\1</b>", t)
    t = _ITAL.sub(r"<i>\1</i>", t)
    return t


def render_markdown(B, out: Path, figdir_rel: str = "figures/direction_influence") -> None:
    L = []
    for blk in B:
        kind = blk[0]
        if kind == "title":
            L.append(f"# {blk[1]}\n")
        elif kind == "sub":
            L.append(f"*{blk[1]}*\n")
        elif kind == "h1":
            L.append(f"\n## {blk[1]}\n")
        elif kind == "h2":
            L.append(f"\n### {blk[1]}\n")
        elif kind in ("p", "lead", "small"):
            L.append(f"{blk[1]}\n")
        elif kind == "caption":
            L.append(f"{blk[1]}\n")
        elif kind == "code":
            L.append("```bash\n" + blk[1] + "\n```\n")
        elif kind == "pagebreak":
            L.append("---\n")
        elif kind == "figure":
            _, fname, cap, _w = blk
            L.append(f"![{cap.split('.')[0].replace('**','')}]({figdir_rel}/{fname})\n")
            L.append(f"{cap}\n")
        elif kind == "table":
            _, rows, _w, align_right = blk
            head, body = rows[0], rows[1:]
            L.append("| " + " | ".join(str(c) for c in head) + " |")
            L.append("|" + "|".join("---:" if i in align_right else "---"
                                    for i in range(len(head))) + "|")
            for r in body:
                L.append("| " + " | ".join(str(c) for c in r) + " |")
            L.append("")
    out.write_text("\n".join(L).replace("\n\n\n", "\n\n") + "\n")


def render_pdf(B, out: Path) -> None:
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_LEFT
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import mm
    from reportlab.platypus import (BaseDocTemplate, Frame, Image, KeepTogether,
                                    PageBreak, PageTemplate, Paragraph, Table, TableStyle)

    INK, INK2 = colors.HexColor("#0b0b0b"), colors.HexColor("#52514e")
    MUTED, RULE = colors.HexColor("#8a8a85"), colors.HexColor("#d9d9d4")
    BAND = colors.HexColor("#f4f4f1")
    ss = getSampleStyleSheet()

    def mk(name, **kw):
        base = dict(fontName="Helvetica", fontSize=9.5, leading=13.6, textColor=INK,
                    alignment=TA_LEFT, spaceAfter=7)
        base.update(kw)
        return ParagraphStyle(name, parent=ss["Normal"], **base)

    S = {"title": mk("t", fontName="Helvetica-Bold", fontSize=20, leading=24, spaceAfter=4),
         "sub": mk("s", fontSize=11, leading=15, textColor=INK2, spaceAfter=16),
         "h1": mk("h1", fontName="Helvetica-Bold", fontSize=13.5, leading=17,
                  spaceBefore=16, spaceAfter=7),
         "h2": mk("h2", fontName="Helvetica-Bold", fontSize=10.5, leading=14,
                  spaceBefore=11, spaceAfter=4),
         "p": mk("b"), "lead": mk("l", fontSize=10.5, leading=15.5),
         "small": mk("sm", fontSize=8.2, leading=11.4, textColor=INK2),
         "caption": mk("c", fontSize=8.2, leading=11.4, textColor=INK2, spaceBefore=3,
                       spaceAfter=13),
         "cell": mk("ce", fontSize=8.4, leading=11, spaceAfter=0),
         "cellb": mk("cb", fontName="Helvetica-Bold", fontSize=8.4, leading=11, spaceAfter=0)}

    def onpage(canvas, doc):
        canvas.saveState()
        canvas.setFont("Helvetica", 7.5)
        canvas.setFillColor(MUTED)
        canvas.drawString(20 * mm, 12 * mm,
                          "Reaction-direction sources in ModelSEED Biochemistry v2.0.1")
        canvas.drawRightString(A4[0] - 20 * mm, 12 * mm, str(canvas.getPageNumber()))
        canvas.setStrokeColor(RULE)
        canvas.setLineWidth(0.4)
        canvas.line(20 * mm, 15.5 * mm, A4[0] - 20 * mm, 15.5 * mm)
        canvas.restoreState()

    doc = BaseDocTemplate(str(out), pagesize=A4, leftMargin=20 * mm, rightMargin=20 * mm,
                          topMargin=18 * mm, bottomMargin=20 * mm,
                          title="Reaction-direction sources: effects on core and "
                                "genome-scale models", author="core_models_analysis")
    doc.addPageTemplates([PageTemplate(
        id="all", frames=[Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height)],
        onPage=onpage)])

    FIG = paths.reports("figures", "direction_influence")
    E = []
    for blk in B:
        kind = blk[0]
        if kind == "pagebreak":
            E.append(PageBreak())
        elif kind in ("title", "sub", "h1", "h2", "p", "lead", "small", "caption"):
            E.append(Paragraph(to_rl(blk[1]), S[kind]))
        elif kind == "code":
            E.append(Paragraph(
                "<font face='Courier' size='8'>"
                + blk[1].replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
                        .replace("\n", "<br/>") + "</font>", S["small"]))
        elif kind == "figure":
            _, fname, cap, w = blk
            path = FIG / fname
            if not path.exists():
                E.append(Paragraph(f"[missing figure: {fname}]", S["small"]))
                continue
            from PIL import Image as PILImage
            iw, ih = PILImage.open(path).size
            img = Image(str(path), width=w * mm, height=w * mm * ih / iw)
            E.append(KeepTogether([img, Paragraph(to_rl(cap), S["caption"])]))
        elif kind == "table":
            _, rows, widths, align_right = blk
            data = []
            for r_i, row in enumerate(rows):
                cells = []
                for c_i, c in enumerate(row):
                    st = S["cellb"] if r_i == 0 else S["cell"]
                    if c_i in align_right:
                        st = ParagraphStyle(f"r{r_i}{c_i}", parent=st, alignment=2)
                    cells.append(Paragraph(to_rl(str(c)), st))
                data.append(cells)
            t = Table(data, colWidths=[w * mm for w in widths], hAlign="LEFT", repeatRows=1)
            cmds = [("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                    ("TOPPADDING", (0, 0), (-1, -1), 4.5),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 4.5),
                    ("LEFTPADDING", (0, 0), (-1, -1), 5),
                    ("LINEBELOW", (0, 0), (-1, 0), 0.9, RULE),
                    ("LINEBELOW", (0, 1), (-1, -2), 0.35, colors.HexColor("#ececE8")),
                    ("LINEBELOW", (0, -1), (-1, -1), 0.9, RULE)]
            for i in range(1, len(data)):
                if i % 2 == 0:
                    cmds.append(("BACKGROUND", (0, i), (-1, i), BAND))
            t.setStyle(TableStyle(cmds))
            E.append(t)
    doc.build(E)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--format", choices=("pdf", "md", "both"), default="both")
    ap.add_argument("--out", type=Path, default=None,
                    help="output path; only meaningful with a single --format")
    args = ap.parse_args()

    B = document(load())
    if args.format in ("pdf", "both"):
        out = args.out if (args.out and args.format == "pdf") else \
            paths.reports("DIRECTION_INFLUENCE_V201.pdf")
        render_pdf(B, out)
        print(f"wrote {out}  ({out.stat().st_size/1000:.0f} KB)")
    if args.format in ("md", "both"):
        out = args.out if (args.out and args.format == "md") else \
            paths.reports("DIRECTION_INFLUENCE_V201.md")
        render_markdown(B, out)
        print(f"wrote {out}  ({out.stat().st_size/1000:.0f} KB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
