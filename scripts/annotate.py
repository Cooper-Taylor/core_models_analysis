#!/usr/bin/env python3
"""
Map cpd/rxn ModelSEED IDs to human names using ModelSEEDDatabase tsv files,
and emit a final annotated INTERPRETATION.md.
"""

import csv
import json
from pathlib import Path
from collections import Counter

ANALYSIS_DIR = Path("/scratch/ctaylor/core_models_analysis")
REPORTS_DIR = ANALYSIS_DIR / "reports"
MSDB = Path("/scratch/ctaylor/ModelSEEDDatabase/Biochemistry")

def load_lookup(prefix):
    """Load id->name from all <prefix>_NN.tsv shards."""
    lookup = {}
    for p in sorted(MSDB.glob(f"{prefix}_*.tsv")):
        if p.name.endswith(".provenance.tsv"):
            continue
        with open(p, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                lookup[row["id"]] = row.get("name", "")
    return lookup


def main():
    print("Loading compound names...")
    cpd_names = load_lookup("compound")
    print(f"  {len(cpd_names)} compounds")
    print("Loading reaction names...")
    rxn_names = load_lookup("reaction")
    print(f"  {len(rxn_names)} reactions")

    # --- compounds: top blocked precursors --------------------------------
    blocked_cpds = [
        ("cpd00002", 2222, 100.0),
        ("cpd00003", 2222, 100.0),
        ("cpd00005", 2222, 100.0),
        ("cpd00022", 2222, 100.0),
        ("cpd00024", 1850, 83.3),
        ("cpd00032", 1714, 77.1),
        ("cpd00101", 1339, 60.3),
        ("cpd00236", 1257, 56.6),
        ("cpd00079", 1179, 53.1),
        ("cpd00072", 1162, 52.3),
        ("cpd00102", 1068, 48.1),
        ("cpd00169", 1050, 47.3),
        ("cpd00061", 1029, 46.3),
        ("cpd00020", 440, 19.8),
    ]
    cpd_table = ["| compound id | name | non-growers blocked | % |",
                 "|---|---|---|---|"]
    for cpd, n, pct in blocked_cpds:
        cpd_table.append(f"| {cpd} | {cpd_names.get(cpd, '?')} | {n} | {pct:.1f}% |")

    # --- reactions: top enriched in growers / non-growers -----------------
    # Re-read the prevalence file to get full lists
    prev_md = (REPORTS_DIR / "REACTION_PREVALENCE.md").read_text().splitlines()

    def parse_section(header):
        # parse markdown table for top reactions
        items = []
        in_section = False
        for line in prev_md:
            if line.startswith("## ") and header in line:
                in_section = True
                continue
            if in_section and line.startswith("## "):
                break
            if in_section and line.startswith("| rxn"):
                parts = [p.strip() for p in line.strip("|").split("|")]
                # rid, g_count, g_pct, n_count, n_pct, delta
                items.append(parts)
        return items

    g_rows = parse_section("most enriched in GROWERS")
    n_rows = parse_section("most enriched in NON-GROWERS")

    def annotated_rxn_table(rows, title):
        out = [f"### {title}",
               "| seed.reaction | name | grower % | non-grower % | Δ |",
               "|---|---|---|---|---|"]
        for r in rows:
            rid = r[0].split("_")[0]  # rxn05466_c -> rxn05466
            name = rxn_names.get(rid, "?")
            out.append(f"| {r[0]} | {name} | {r[2]} | {r[4]} | {r[5]} |")
        return out

    # --- assemble INTERPRETATION.md ---------------------------------------
    lines = []
    w = lines.append
    w("# Core Models KEGG2 — Biological Interpretation\n")
    w(f"5,683 ModelSEED metabolic models tested for biomass production under the "
      f"ModelSEEDDatabase complete media (`KBaseMedia.cpd`, 347 compounds).\n")
    w("## Headline numbers")
    w("- **3,461 (60.9%)** produced biomass (FBA optimum > 1e-6 on `bio1`)")
    w("- **2,222 (39.1%)** solved to optimal but with zero biomass — networks with gaps")
    w("- 0 errors, 0 missing biomass reactions, 0 infeasible solves\n")

    w("## Growers are larger and gene-richer")
    w("| metric | grower median | non-grower median | gap |")
    w("|---|---|---|---|")
    w("| metabolites | 158 | 125 | +33 |")
    w("| reactions   | 175 | 118 | +57 |")
    w("| genes       | 189 | 101 | +88 |")
    w("| exchanges   | 27  | 23  | +4  |")
    w("\nThe non-growers are nearly half-size by gene count — these are draft "
      "models with too little metabolism to close the biomass equation. The "
      "exchange-count gap is small, so it's not a transport problem; it's a "
      "cytoplasmic-pathway-completeness problem.\n")

    w("## Grower biomass flux is unimodal and high")
    w("- median 52.3, mean 52.4, max 87.2")
    w("- only 33 growers fall below flux 10; almost all of the 3,461 grow vigorously\n")

    w("## What's blocking non-growers — top biomass precursors that cannot be made")
    w("\n".join(cpd_table))
    w("")
    w("**Reading the table:** four cofactors are universally blocked across all "
      "2,222 non-growers — ATP, NAD+, NADPH, and Acetyl-CoA. This is the "
      "signature of a core-energy-metabolism gap: without a closed loop for "
      "regenerating these carriers, the biomass reaction (which consumes them "
      "stoichiometrically) cannot run at any rate.\n")
    w("Beyond cofactors, the next-most-blocked precursors are TCA-cycle and "
      "central-carbon intermediates (2-oxoglutarate, oxaloacetate, "
      "3-phosphoglycerate, etc.), confirming that the missing capability is "
      "central carbon catabolism and oxidative phosphorylation, not biosynthesis "
      "of exotic biomass building blocks.\n")
    w("Per-model, non-growers are missing a median of 9 (of 14) biomass "
      "precursors — most of the biomass equation, not a single bottleneck.\n")

    w("## Which reactions correlate with growth?")
    w("Reaction prevalence compared between the two cohorts. Δ = grower% − non-grower%.\n")
    w("\n".join(annotated_rxn_table(g_rows[:20], "Top reactions enriched in GROWERS (+Δ)")))
    w("")
    w("\n".join(annotated_rxn_table(n_rows[:15], "Top reactions enriched in NON-GROWERS (−Δ)")))
    w("")
    w("**Interpretation:** the grower-enriched reactions cluster on a few axes:")
    w("- **TCA cycle**: citrate synthase (rxn00256), isocitrate dehydrogenase "
      "(rxn01387), oxalosuccinate → 2-oxoglutarate (rxn00199), 2-oxoglutarate "
      "dehydrogenase E2 (rxn01872), 2-oxoglutarate decarboxylation via TPP "
      "(rxn00441), and the glyoxylate-shunt enzymes malate synthase (rxn00330) "
      "and isocitrate lyase (rxn00336).")
    w("- **Acetyl-CoA activation**: acetate:CoA ligase (rxn00175).")
    w("- **Fatty-acid β-oxidation**: hydroxyacyl-CoA hydro-lyases / "
      "dehydrogenases (rxn02167, rxn03249/03250, rxn03240, rxn06777).")
    w("- **Respiration / electron transport**: succinate dehydrogenase "
      "(rxn10126), ferredoxin-NADP+ oxidoreductase (rxn05937), and "
      "dissimilatory nitrate reductases (rxn09001/09003, rxn14412/14414/14427).")
    w("- **Sulfur assimilation**: ATP-sulfurylase (rxn00379), needed because the "
      "biomass equation pulls sulfur-containing building blocks.")
    w("- **Dicarboxylate uptake**: rxn05561 (transports succinate/malate-like "
      "substrates into the cell).")
    w("")
    w("Non-growers are *enriched* in proton-symport transport reactions "
      "(rxn05319, rxn05466, rxn05469 pyruvate, rxn05488 acetate, rxn05559 "
      "formate) — they have the import machinery but not the downstream "
      "catabolism, so the imported carbon never reaches the biomass equation. "
      "They're also enriched in alternative/fermentative-style oxidoreductases "
      "(pyruvate:ferredoxin rxn13974, 2-oxoglutarate synthase rxn05939, "
      "CODH/ACS rxn15962), which are characteristic of obligate-anaerobe drafts "
      "where ModelSEED's default biomass — calibrated around aerobic central "
      "metabolism — doesn't close.\n")

    w("## Bottom line")
    w("Growth capability in this collection is essentially binary and tracks "
      "**model completeness**, not media composition. The 60.9% of models that "
      "grow have the central-carbon + TCA + energy backbone needed to feed the "
      "ModelSEED biomass equation; the 39.1% that don't are missing roughly "
      "two-thirds of that backbone (median 9 of 14 precursors unreachable, "
      "100% missing the four energy/redox cofactors). On complete media these "
      "non-growers would need gap-filling — not nutrient addition — to grow.\n")

    out = REPORTS_DIR / "INTERPRETATION.md"
    out.write_text("\n".join(lines))
    print(out.read_text())
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
