#!/usr/bin/env python3
"""Summarize results.csv from the growth analysis."""

import csv
from collections import Counter
from pathlib import Path
from statistics import mean, median

ANALYSIS_DIR = Path("/scratch/ctaylor/core_models_analysis")
RESULTS_CSV = ANALYSIS_DIR / "results" / "results.csv"
SUMMARY_MD = ANALYSIS_DIR / "reports" / "SUMMARY.md"
NONGROWERS_CSV = ANALYSIS_DIR / "results" / "non_growers.csv"
GROWERS_CSV = ANALYSIS_DIR / "results" / "growers.csv"

rows = []
with open(RESULTS_CSV) as f:
    for r in csv.DictReader(f):
        for k in ("n_metabolites", "n_reactions", "n_genes",
                  "n_exchanges_total", "n_exchanges_open"):
            r[k] = int(r[k])
        r["growth_flux"] = float(r["growth_flux"])
        r["grows"] = r["grows"] == "True"
        rows.append(r)

total = len(rows)
status_counts = Counter(r["status"] for r in rows)
n_errors = sum(1 for r in rows if r["error"])
n_no_biomass = sum(1 for r in rows if r["status"] == "no_biomass")
growers = [r for r in rows if r["grows"]]
non_growers = [r for r in rows if not r["grows"]]

# Split non-growers
zero_flux = [r for r in non_growers if r["status"] == "optimal" and r["growth_flux"] <= 1e-6]
infeasible = [r for r in non_growers if r["status"] not in ("optimal", "no_biomass", "")]
errored = [r for r in non_growers if r["error"]]

g_flux = [r["growth_flux"] for r in growers]

lines = []
def w(s=""): lines.append(s)

w(f"# Core Models KEGG2 — Growth Analysis Summary")
w()
w(f"- Models analyzed: **{total}**")
w(f"- Media: ModelSEEDDatabase `KBaseMedia.cpd` (complete media, 347 compounds)")
w(f"- Growth criterion: FBA on biomass reaction (`bio1`) > 1e-6")
w()
w(f"## Outcomes")
w(f"| Outcome | Count | % |")
w(f"|---|---|---|")
w(f"| **Grows (biomass flux > 1e-6)** | {len(growers)} | {100*len(growers)/total:.1f}% |")
w(f"| Optimal solve but zero biomass flux | {len(zero_flux)} | {100*len(zero_flux)/total:.1f}% |")
w(f"| Non-optimal solver status (infeasible/unbounded/etc.) | {len(infeasible)} | {100*len(infeasible)/total:.1f}% |")
w(f"| No biomass reaction found | {n_no_biomass} | {100*n_no_biomass/total:.1f}% |")
w(f"| Model load / processing error | {n_errors} | {100*n_errors/total:.1f}% |")
w()
w(f"## Solver status distribution")
for s, c in status_counts.most_common():
    w(f"- `{s or '(empty)'}` : {c}")
w()
if g_flux:
    w(f"## Biomass flux among growers (n={len(g_flux)})")
    w(f"- min: {min(g_flux):.4g}")
    w(f"- median: {median(g_flux):.4g}")
    w(f"- mean: {mean(g_flux):.4g}")
    w(f"- max: {max(g_flux):.4g}")
    w()
w(f"## Model sizes (n_reactions)")
sizes = sorted(r["n_reactions"] for r in rows)
w(f"- min: {sizes[0]}  median: {sizes[len(sizes)//2]}  max: {sizes[-1]}  mean: {mean(sizes):.1f}")
w()
w(f"## Files")
w(f"- `results.csv` — per-model row")
w(f"- `growers.csv` — subset where biomass flux > 1e-6")
w(f"- `non_growers.csv` — everything else (includes errors)")
w(f"- `failures.log` — solver-status / exception details")
w(f"- `analyze_growth.py` — analysis script")
w(f"- `summarize.py` — this summary script")

SUMMARY_MD.write_text("\n".join(lines) + "\n")
print(SUMMARY_MD.read_text())

# Write growers / non-growers split for convenience
fieldnames = list(rows[0].keys())
with open(GROWERS_CSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames); w.writeheader()
    for r in growers:
        w.writerow(r)
with open(NONGROWERS_CSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames); w.writeheader()
    for r in non_growers:
        w.writerow(r)
print(f"Wrote {GROWERS_CSV} ({len(growers)} rows)")
print(f"Wrote {NONGROWERS_CSV} ({len(non_growers)} rows)")
