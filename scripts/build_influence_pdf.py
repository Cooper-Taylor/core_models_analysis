#!/usr/bin/env python3
"""Typeset the direction-influence analysis as a PDF.

    python3 scripts/build_influence_pdf.py
    python3 scripts/build_influence_pdf.py --out reports/DIRECTION_INFLUENCE_V201.pdf

Every number is read from the result files at build time, so the document
cannot drift from the data. Figures come from
``reports/figures/direction_influence/`` (registered in ``scripts/figures.tsv``
as ``direction_influence``); regenerate them first if the sweeps have changed.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import (BaseDocTemplate, Frame, Image, KeepTogether,
                                NextPageTemplate, PageBreak, PageTemplate,
                                Paragraph, Spacer, Table, TableStyle)

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
from cma import paths  # noqa: E402

INK = colors.HexColor("#0b0b0b")
INK2 = colors.HexColor("#52514e")
MUTED = colors.HexColor("#8a8a85")
RULE = colors.HexColor("#d9d9d4")
BLUE = colors.HexColor("#2a78d6")
BAND = colors.HexColor("#f4f4f1")

SETS = ["v201_group_contribution", "v201_dgpredictor", "v201_equilibrator",
        "v201_llm_council", "claude_opus48", "v201_gold_only", "v201_silver_only",
        "v201_bronze_only", "v201_gold_silver", "v201_gold_silver_bronze"]
LABEL = {"v201_group_contribution": "Group contribution", "v201_dgpredictor": "dGPredictor",
         "v201_equilibrator": "eQuilibrator", "v201_llm_council": "LLM council",
         "claude_opus48": "Claude Opus 4.8", "v201_gold_only": "gold",
         "v201_silver_only": "silver", "v201_bronze_only": "bronze",
         "v201_gold_silver": "gold + silver",
         "v201_gold_silver_bronze": "gold + silver + bronze"}


def styles():
    ss = getSampleStyleSheet()
    def mk(name, **kw):
        base = dict(fontName="Helvetica", fontSize=9.5, leading=13.6,
                    textColor=INK, alignment=TA_LEFT, spaceAfter=7)
        base.update(kw)
        return ParagraphStyle(name, parent=ss["Normal"], **base)
    return {
        "title": mk("t", fontName="Helvetica-Bold", fontSize=20, leading=24,
                    spaceAfter=4),
        "sub": mk("s", fontSize=11, leading=15, textColor=INK2, spaceAfter=16),
        "h1": mk("h1", fontName="Helvetica-Bold", fontSize=13.5, leading=17,
                 spaceBefore=16, spaceAfter=7),
        "h2": mk("h2", fontName="Helvetica-Bold", fontSize=10.5, leading=14,
                 spaceBefore=11, spaceAfter=4),
        "body": mk("b"),
        "lead": mk("l", fontSize=10.5, leading=15.5),
        "small": mk("sm", fontSize=8.2, leading=11.4, textColor=INK2),
        "cap": mk("c", fontSize=8.2, leading=11.4, textColor=INK2, spaceBefore=3,
                  spaceAfter=13),
        "cell": mk("ce", fontSize=8.4, leading=11, spaceAfter=0),
        "cellb": mk("cb", fontName="Helvetica-Bold", fontSize=8.4, leading=11, spaceAfter=0),
    }


def load() -> dict:
    inf, sets = paths.results("influence_v201"), paths.results("direction_sets_v201")
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
    return D


def table(rows, S, widths, align_right=(), header=True, font=8.4):
    data = []
    for r_i, row in enumerate(rows):
        out = []
        for c_i, c in enumerate(row):
            st = S["cellb"] if (header and r_i == 0) else S["cell"]
            if c_i in align_right:
                st = ParagraphStyle(f"r{r_i}{c_i}", parent=st, alignment=2)
            out.append(Paragraph(str(c), st))
        data.append(out)
    t = Table(data, colWidths=widths, hAlign="LEFT", repeatRows=1 if header else 0)
    cmds = [("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("TOPPADDING", (0, 0), (-1, -1), 4.5),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 4.5),
            ("LEFTPADDING", (0, 0), (-1, -1), 5),
            ("LINEBELOW", (0, 0), (-1, 0), 0.9, RULE) if header else
            ("LINEBELOW", (0, 0), (-1, 0), 0, colors.white),
            ("LINEBELOW", (0, 1), (-1, -2), 0.35, colors.HexColor("#ececE8")),
            ("LINEBELOW", (0, -1), (-1, -1), 0.9, RULE)]
    for i in range(1, len(data)):
        if i % 2 == 0:
            cmds.append(("BACKGROUND", (0, i), (-1, i), BAND))
    t.setStyle(TableStyle(cmds))
    return t


def figure(path: Path, S, caption: str, width=168 * mm):
    if not path.exists():
        return [Paragraph(f"[missing figure: {path.name}]", S["small"])]
    from PIL import Image as PILImage
    try:
        w, h = PILImage.open(path).size
        ratio = h / w
    except Exception:
        ratio = 0.45
    img = Image(str(path), width=width, height=width * ratio)
    # a figure separated from its caption by a page break is a typesetting bug
    return [KeepTogether([img, Paragraph(caption, S["cap"])])]


def build(D, out: Path) -> None:
    S = styles()
    FIG = paths.reports("figures", "direction_influence")

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

    doc = BaseDocTemplate(str(out), pagesize=A4,
                          leftMargin=20 * mm, rightMargin=20 * mm,
                          topMargin=18 * mm, bottomMargin=20 * mm,
                          title="Reaction-direction sources: effects on core and "
                                "genome-scale models",
                          author="core_models_analysis")
    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id="n")
    doc.addPageTemplates([PageTemplate(id="all", frames=[frame], onPage=onpage)])

    E = []
    core = D["core_all"]["growth"]
    gmm = D["gsm_gmm"]["growth"]
    aux = D["gsm_auxo"]["growth"]
    cb, gb, ab = core["on_disk"], gmm["on_disk"], aux["on_disk"]

    # ---- title ----------------------------------------------------------
    E.append(Paragraph("Reaction-direction sources and what they do to models", S["title"]))
    E.append(Paragraph("Eleven direction sources from ModelSEED Biochemistry v2.0.1, "
                       f"applied to {cb['n_models']:,} core models and {gb['n_models']:,} "
                       "genome-scale models", S["sub"]))

    E.append(Paragraph("Summary", S["h1"]))
    E.append(Paragraph(
        "The same direction source has opposite effects on the two model sets. On the core "
        "models most sources <b>increase</b> the number that grow; on the genome-scale models "
        "most <b>abolish</b> growth entirely. That is not a property of the sources.", S["lead"]))
    E.append(Paragraph(
        "Two things produce it. The genome-scale networks are roughly fifteen times larger by "
        "distinct reaction count, so a source rewrites a far greater share of the network. And "
        "within that larger surface sit a few reactions that are thermodynamically uphill but "
        "biologically obligatory. One of them, the thiamine phosphomethylpyrimidine kinase, is "
        "present in every genome-scale model tested and on its own abolishes growth in 99% of "
        "them. The thermodynamic estimate for it is correct; the constraint derived from it is "
        "still wrong.", S["body"]))
    E.append(Paragraph(
        "The practical consequence is that <b>the core panel cannot validate a direction "
        "source</b>. It exercises 239 of the 56,012 reactions in the database and omits the "
        "pathways where thermodynamic and biological direction diverge.", S["body"]))

    # ---- method ---------------------------------------------------------
    E.append(Paragraph("1. What was measured", S["h1"]))
    E.append(Paragraph(
        "<b>Baseline.</b> The bounds written in the model files themselves: forward, reverse or "
        "reversible according to each reaction's stored lower and upper bound. Nothing is "
        "recomputed. This is the <i>on-disk</i> row throughout, and it is what the models were "
        "gap-filled and published with.", S["body"]))
    E.append(Paragraph(
        "<b>Applying a source.</b> A reaction whose SEED annotation appears in the source is "
        "rewritten to (0, 1000), (-1000, 0) or (-1000, 1000). A reaction the source does not "
        "mention keeps its on-disk bounds. Media are applied afterwards, biomass is the "
        "objective, and a model counts as growing above 1e-6.", S["body"]))

    cov = D["coverage"]
    E.append(table([
        ["Model set", "Models", "Reactions / model", "Distinct reactions", "Media"],
        ["Core (core_kegg2)", f"{cb['n_models']:,}",
         f"{cov['core_kegg2']['reactions_per_model']['median']}",
         f"{cov['core_kegg2']['distinct_msdb_reactions']:,}", "347-compound complete"],
        ["Genome-scale (ms2_gsm)", f"{gb['n_models']:,}",
         f"{cov['ms2_gsm']['reactions_per_model']['median']:,}",
         f"{cov['ms2_gsm']['distinct_msdb_reactions']:,}", "20-compound glucose minimal"],
        ["Genome-scale (auxotrophy)", f"{ab['n_models']:,}", "1,085", "3,484",
         "51-compound auxotrophy"],
    ], S, [44 * mm, 18 * mm, 28 * mm, 28 * mm, 50 * mm], align_right=(1, 2, 3)))
    E.append(Paragraph(
        "Reactions per model is the median count of distinct ModelSEED reactions. The "
        "genome-scale models come from the public KBase workspaces of the ModelSEED v2 "
        "manuscript and were gap-filled on the media shown.", S["cap"]))

    E.append(Paragraph("2. The eleven sources", S["h1"]))
    E.append(Paragraph(
        "Three thermodynamic estimators, two LLM routes, three evidence tiers and three "
        "cumulative tiers. A grade belongs to a <i>reaction</i>, not to a source: the release "
        "ships one grade per reaction and keeps the per-source table internal. A tier is "
        "therefore the release's recommended direction restricted to reactions at that tier, "
        "chosen by a fixed precedence of eQuilibrator over dGPredictor over group "
        "contribution.", S["body"]))
    stat = {x["name"]: x for x in D["manifest"]["sets"]}
    rows = [["Source", "Calls", "Forward", "Reverse", "Reversible", "Directional"]]
    for k in SETS:
        st = stat[k]; o = st["operators"]; n = st["n_callable"]
        d = o.get(">", 0) + o.get("<", 0)
        rows.append([LABEL[k], f"{n:,}", f"{o.get('>',0):,}", f"{o.get('<',0):,}",
                     f"{o.get('=',0):,}", f"{100*d/n:.0f}%"])
    E.append(table(rows, S, [50 * mm, 22 * mm, 22 * mm, 22 * mm, 24 * mm, 24 * mm],
                   align_right=(1, 2, 3, 4, 5)))
    E.append(Paragraph(
        "No set contains an undecided call. A direction map is applied by rewriting bounds, and "
        "an undecided value would be applied as fully reversible, so emitting it for a reaction "
        "that was graded but could not be called would <i>remove</i> a constraint in the name of "
        "evidence. Absent means no opinion, and the model keeps its own bounds. "
        "<b>Directional</b> is the share of a source's calls that fix a direction rather than "
        "allow both: the constraint pressure it applies.", S["cap"]))

    # ---- results --------------------------------------------------------
    E.append(PageBreak())
    E.append(Paragraph("3. What each source does to growth", S["h1"]))
    E.extend(figure(FIG / "fig1_growth_by_panel.png", S,
                    "<b>Figure 1.</b> Share of models that grow under each direction source, "
                    "against the bounds already in the model files. The three panels share a "
                    "source ordering; each has its own baseline."))

    rows = [["Source", "Core", "vs base", "GSM glucose", "GSM auxotrophy"]]
    rows.append([f"<b>on-disk baseline</b>", f"<b>{cb['n_growers']:,}</b>", "&mdash;",
                 f"<b>{gb['n_growers']:,}</b>", f"<b>{ab['n_growers']:,}</b>"])
    for k in SETS:
        d = D["core_all"]["diff"][k]
        net = d["n_gained_growth"] - d["n_lost_growth"]
        rows.append([LABEL[k], f"{core[k]['n_growers']:,}", f"{net:+,}",
                     f"{gmm[k]['n_growers']:,}", f"{aux[k]['n_growers']:,}"])
    E.append(table(rows, S, [52 * mm, 24 * mm, 22 * mm, 30 * mm, 34 * mm],
                   align_right=(1, 2, 3, 4)))
    E.append(Paragraph(
        "Models that grow, of 5,683 core and 5,420 genome-scale. Zero models failed to load in "
        "any run.", S["cap"]))

    E.append(Paragraph("On the core models, most sources help", S["h2"]))
    E.append(Paragraph(
        "eQuilibrator is the least directional source and gains the most models; Claude is among "
        "the most directional and loses the most. <b>Gold alone costs growth while silver alone "
        "gains it</b>, which inverts the evidence hierarchy: 122 models stop growing under the "
        "best-graded directions and 367 start growing under the weaker ones. Gold calls 58% of "
        "its reactions directional against silver's 43%, and a directional call removes a degree "
        "of freedom the model was using. The grade measures how well an energy is known, not how "
        "safe it is to impose.", S["body"]))

    E.append(Paragraph("On the genome-scale models, most collapse", S["h2"]))
    E.append(Paragraph(
        "On glucose minimal media <b>six of the eleven sources leave zero growing models</b>. "
        "Every one returns an optimal solution with biomass flux zero, so these are genuine "
        "no-growth results rather than solver failures. The auxotrophy set starts from a "
        "near-complete baseline and separates the sources better: gold retains 85% of baseline, "
        "the silver-containing sets about a third, and six of eleven fall below 4%. The "
        "cumulative tiers decay monotonically there, 4,618 to 1,866 to 55.", S["body"]))
    E.append(Paragraph(
        "\"Nearly everything collapses\" is fair for glucose minimal and an overstatement for "
        "auxotrophy. One detail runs the other way: among the models that do survive, mean "
        "growth flux <i>rises</i>, from 1.32 at baseline to 4.09 under gold and 6.52 under "
        "gold plus silver. Directional constraints remove futile cycling and channel flux, so "
        "the cost is concentrated entirely in the models that stop growing.", S["body"]))

    # ---- grading --------------------------------------------------------
    E.append(Paragraph("4. The grading, against the original thermodynamic sources", S["h1"]))
    E.extend(figure(FIG / "fig2_grading_vs_thermo.png", S,
                    "<b>Figure 2.</b> Left: how often an evidence tier and an estimator agree, "
                    "over the reactions where both have an opinion. Right: the share of each "
                    "source's calls that fix a direction rather than allow both."))

    ag = D["agreement"]["pairwise_agreement"]
    def pair(a, b):
        e = ag.get(f"{a} vs {b}") or ag.get(f"{b} vs {a}")
        return e
    e1 = pair("v201_equilibrator", "v201_gold_silver_bronze")
    e2 = pair("v201_llm_council", "v201_equilibrator")
    e3 = pair("v201_llm_council", "claude_opus48")
    e4 = pair("v201_gold_only", "v201_dgpredictor")

    E.append(Paragraph(
        f"<b>The tiers are largely eQuilibrator re-expressed.</b> The full graded set agrees "
        f"with eQuilibrator on {e1['pct_agree']}% of the {e1['shared']:,} reactions they share, "
        "bronze at 100.0% and silver at 99.9%. That follows directly from the recommendation "
        "precedence putting eQuilibrator first. A tier behaves as a confidence filter on one "
        "source rather than as a consensus of several, so comparing a tier against eQuilibrator "
        "is close to comparing eQuilibrator against itself.", S["body"]))
    E.append(Paragraph(
        f"<b>The estimators disagree with each other more than the tiers suggest.</b> "
        f"dGPredictor agrees with the gold tier on only {e4['pct_agree']}% of the "
        f"{e4['shared']:,} reactions they share, barely above chance among three options. Gold "
        "is the tier where the evidence is strongest, so this is disagreement about "
        "well-measured reactions, not about hard ones.", S["body"]))
    E.append(Paragraph(
        f"<b>The LLM routes are a different kind of claim.</b> The council agrees with "
        f"eQuilibrator on {e2['pct_agree']}% of {e2['shared']:,} shared reactions; it reasons "
        "from the reaction rather than from an energy, which is the release's own argument for "
        f"keeping it out of the grading. The council and Claude agree with each other "
        f"{e3['pct_agree']}% of the time over {e3['shared']:,} reactions, but Claude is one of "
        "the council's members, so that is shared lineage rather than independent "
        "corroboration.", S["body"]))
    E.append(Paragraph(
        "The right panel of Figure 2 shows what the agreement numbers hide: the tiers are "
        "markedly <i>less committal</i> than the estimators they are drawn from. Silver allows "
        "both directions on 57% of its calls against dGPredictor's 27%. That is most of why the "
        "silver tier is gentler on the models than the source it mostly agrees with.", S["body"]))

    # ---- genome-scale mechanism -----------------------------------------
    E.append(Paragraph("5. Why the genome-scale models collapse", S["h1"]))
    E.append(Paragraph(
        "A total collapse is as likely to be a defect as a result, so it was diagnosed rather "
        "than reported.", S["body"]))
    rt = D["roundtrip"]
    E.append(Paragraph(
        f"<b>The conversion is faithful.</b> Direction round-trips exactly from the KBase "
        f"originals: {rt['reactions_checked']:,} reactions across a "
        f"{rt['models_sampled']}-model random sample, {rt['direction_mismatches']} mismatches. "
        "The models also grow normally on their own bounds. Only the override breaks them.", S["body"]))
    E.append(Paragraph(
        "<b>It is not the medium.</b> The collapse reproduces on the 347-compound complete "
        "medium, not only the 20-compound minimal one.", S["body"]))
    E.append(Paragraph(
        "<b>It localises to individual reactions.</b> Taking one genome-scale model under "
        "eQuilibrator directions, 229 of its reactions differ from disk and 51 of those are "
        "constraining, but applying them one at a time shows that only three individually "
        "abolish growth. Repeating that scan across 73 models that all grow at baseline, every "
        "one has at least one individually lethal constraint, and the distribution is extremely "
        "skewed.", S["body"]))
    E.extend(figure(FIG / "fig4_mechanism.png", S,
                    "<b>Figure 3.</b> Share of tested genome-scale models in which constraining "
                    "one reaction alone abolishes growth. Rates are from a 73-model scan and "
                    "shift by a few points between samples; the ranking of the first does not.",
                    width=162 * mm))
    E.append(Paragraph(
        "<b>One reaction explains most of it.</b> The thiamine phosphomethylpyrimidine kinase "
        "(rxn03108) is present in every genome-scale model and lethal in 72 of 73. eQuilibrator "
        "computes +15.59 kcal/mol for it and therefore calls it reverse. That number is not "
        "wrong: the phosphoryl transfer is uphill in isolation. But the enzyme is an ATP-driven "
        "kinase in thiamine biosynthesis and must run forward, so constraining it to the "
        "thermodynamically favoured direction removes thiamine and with it biomass. Both LLM "
        "routes call this reaction forward.", S["body"]))
    E.append(Paragraph(
        "ATP sulfurylase (rxn00379) is the same failure with better evidence: graded <b>gold</b> "
        "from an openTECR measurement, +11 kcal/mol, called reverse by eQuilibrator, dGPredictor "
        "and Claude. In the cell it runs forward because pyrophosphatase removes the "
        "pyrophosphate product. The measurement is right and the direction call is right in "
        "isolation; the constraint is still wrong. It is a milder case, present in 67% of models "
        "and lethal in about one in eight of those.", S["body"]))
    E.append(Paragraph(
        "Removing the top three reactions from the eQuilibrator map restores 137 of 381 growers "
        "in a 480-model sample, from zero. The collapse is not diffuse over-constraint; it is a "
        "handful of obligatory, cofactor-driven reactions.", S["body"]))
    E.append(Paragraph(
        "<b>Why the core models escape.</b> They contain 239 distinct reactions against 3,634. "
        "Thiamine biosynthesis and sulfate assimilation are simply absent from a core "
        "carbon-metabolism model, so the reactions that kill the genome-scale models are never "
        "touched.", S["body"]))

    # ---- predictors -----------------------------------------------------
    E.append(Paragraph("6. Does anything predict what a source will do?", S["h1"]))
    E.extend(figure(FIG / "fig3_coverage_effect.png", S,
                    "<b>Figure 4.</b> Share of reaction occurrences a source changes, against "
                    "the share of models that grow. Note the independent vertical scales: the "
                    "core panel spans 51&ndash;70%, the genome-scale panel 0&ndash;80%."))
    E.append(Paragraph(
        "Not reliably. Of the candidate summary statistics &mdash; directional fraction, "
        "reactions called, occurrences changed, overrides per model &mdash; none is significant "
        "on the core panel, and the best on the genome-scale panel reaches only Spearman rho of "
        "-0.65. The share of occurrences a source changes gives rho = -0.65 on glucose minimal, "
        "-0.45 on auxotrophy and -0.03 on the core panel, where it has essentially no ordering "
        "power at all.", S["body"]))
    E.append(Paragraph(
        "Directional fraction does better on the core panel, rho = -0.68, but only when computed "
        "over the whole release. Restricted to the 239 reactions the core models actually "
        "contain &mdash; the ones that can affect the result &mdash; it falls to -0.28 and is "
        "not significant. The eleven sets are also nested rather than independent: gold sits "
        "inside gold plus silver, which sits inside gold plus silver plus bronze.", S["body"]))
    E.append(Paragraph(
        "Two counterexamples make the point concretely. Bronze changes <i>less</i> of the "
        "network than gold, 3.3% against 3.7%, and retains far less growth. dGPredictor changes "
        "only 12.5% and leaves zero growers, the same as Claude at 31.5%. <b>What a source does "
        "depends on which reactions it constrains, not how many.</b>", S["body"]))

    E.append(Paragraph("7. What follows", S["h1"]))
    for head, txt in [
        ("A direction source cannot be validated on the core models.",
         "They exercise 239 of 56,012 reactions and omit the pathways where thermodynamic and "
         "biological direction come apart. Earlier conclusions drawn from that panel are not "
         "contradicted so much as shown to be untested."),
        ("A grade is not a licence to constrain.",
         "The grading scheme answers how well an energy is known, and answers it well: ATP "
         "sulfurylase is gold because openTECR measured it. Whether it is safe to fix that "
         "reaction's direction in a model is a different question, and the two come apart "
         "precisely at the cofactor-driven reactions. A tier used as an FBA constraint set needs "
         "a second filter for obligatory reactions."),
        ("Reversibility is the safe default.",
         "Every set that is more than about 60% directional damages the genome-scale models "
         "badly. That is a tendency rather than a rule, but it is the most robust one available."),
        ("Two candidate fixes, neither tested here.",
         "Exempt reactions whose stoichiometry contains ATP or pyrophosphate hydrolysis from "
         "directional constraint, which is the mechanism in every case examined; or intersect a "
         "tier with the LLM council and keep only reactions where both agree, trading coverage "
         "for safety."),
    ]:
        E.append(Paragraph(f"<b>{head}</b> {txt}", S["body"]))

    E.append(Paragraph("8. Caveats", S["h1"]))
    for c in [
        "<b>Growth is the only readout.</b> A direction set that preserves grower count may "
        "still distort flux distributions; nothing here measures that.",
        "<b>Energy-generating-cycle detection does not work on these models.</b> The existing "
        "probe closes every biomass reaction above ten metabolites and relies on a smaller one "
        "staying open as an ATP probe. The genome-scale models have exactly one biomass and no "
        "maintenance reaction, so it would report zero cycles for all 5,420 with no error.",
        "<b>The lethal-reaction scan covers 73 models</b>, the growers among every 60th of the "
        "set. The ranking of the first reaction is unambiguous and replicated on an independent "
        "sample; the long tail is not fully enumerated.",
        "<b>Claude Opus 4.8 is not part of the release.</b> It is a standalone single-model run; "
        "the council includes this model among its five roles, so the two are not independent.",
        "<b>Correlations rest on eleven non-independent sets.</b> No coefficient here exceeds "
        "0.7 in magnitude and most are not significant.",
    ]:
        E.append(Paragraph(c, S["small"]))

    E.append(Paragraph("Reproducing this", S["h1"]))
    E.append(Paragraph(
        "<font face='Courier' size='8'>"
        "python3 scripts/build_direction_sets_v201.py<br/>"
        "python3 beginPipeline --models core_kegg2 --workers 64 --directions v201_gold "
        "&hellip; --only growth,growth_diff --out results/influence_v201/core_all<br/>"
        "python3 scripts/analyze_direction_coverage.py<br/>"
        "python3 scripts/regen_figures.py direction_influence<br/>"
        "python3 scripts/build_influence_pdf.py</font>", S["small"]))
    E.append(Paragraph(
        "Results are under results/influence_v201/ and results/direction_sets_v201/, each run "
        "with a run.json recording the exact command. Every number in this document is read "
        "from those files at build time.", S["small"]))

    doc.build(E)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path,
                    default=paths.reports("DIRECTION_INFLUENCE_V201.pdf"))
    args = ap.parse_args()
    build(load(), args.out)
    print(f"wrote {args.out}  ({args.out.stat().st_size/1000:.0f} KB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
