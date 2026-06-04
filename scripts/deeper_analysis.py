#!/usr/bin/env python3
"""
Deeper analysis of growth results.

1. Compare grower vs non-grower model characteristics.
2. Gap analysis: for non-growers, find biomass precursors that cannot
   be produced (sink test on each precursor).
3. Reaction prevalence: which seed.reaction IDs are over-represented
   in growers vs non-growers (chi-square-ish score, but simple ratios).

Outputs are written into /scratch/ctaylor/core_models_analysis/.
"""

import csv
import json
import multiprocessing as mp
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean, median, pstdev
import os

import cobra
from cobra.io import load_json_model

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR", "/scratch/ctaylor/core_models_analysis"))
MODELS_DIR = ANALYSIS_DIR / "data" / "core_models_kegg2"
MEDIA_FILE = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase") + "/Media/KBaseMedia.cpd")
RESULTS_CSV = ANALYSIS_DIR / "results" / "results.csv"
REPORTS_DIR = ANALYSIS_DIR / "reports"
RESULTS_DIR = ANALYSIS_DIR / "results"

# ---------------------------------------------------------------------
# load growth results
def load_results():
    rows = []
    with open(RESULTS_CSV) as f:
        for r in csv.DictReader(f):
            for k in ("n_metabolites", "n_reactions", "n_genes",
                      "n_exchanges_total", "n_exchanges_open"):
                r[k] = int(r[k])
            r["growth_flux"] = float(r["growth_flux"])
            r["grows"] = r["grows"] == "True"
            rows.append(r)
    return rows


def load_media():
    with open(MEDIA_FILE) as f:
        return {line.strip() for line in f if line.strip()}


# ---------------------------------------------------------------------
# Part 1: characteristics
def part1_characteristics(rows, outpath):
    growers = [r for r in rows if r["grows"]]
    nongrowers = [r for r in rows if not r["grows"]]

    def stats(seq):
        s = sorted(seq)
        return dict(min=s[0], median=s[len(s)//2], max=s[-1],
                    mean=mean(s), sd=pstdev(s))

    keys = ("n_metabolites", "n_reactions", "n_genes",
            "n_exchanges_total", "n_exchanges_open")
    lines = ["# Part 1 — Grower vs Non-grower characteristics\n"]
    lines.append(f"Growers: {len(growers)}    Non-growers: {len(nongrowers)}\n")
    lines.append("| metric | grower median | grower mean | non-grower median | non-grower mean |")
    lines.append("|---|---|---|---|---|")
    for k in keys:
        g = stats([r[k] for r in growers])
        n = stats([r[k] for r in nongrowers])
        lines.append(f"| {k} | {g['median']:.1f} | {g['mean']:.1f} | "
                     f"{n['median']:.1f} | {n['mean']:.1f} |")
    lines.append("")

    # Flux distribution buckets among growers
    bins = [0, 1, 5, 10, 25, 50, 75, 100]
    bins_count = Counter()
    for r in growers:
        f = r["growth_flux"]
        for hi in bins[1:]:
            if f <= hi:
                bins_count[hi] += 1
                break
        else:
            bins_count["inf"] += 1
    lines.append("\n## Grower biomass flux distribution")
    lo = 0
    for hi in bins[1:]:
        lines.append(f"- ({lo}, {hi}]: {bins_count[hi]}")
        lo = hi
    if bins_count["inf"]:
        lines.append(f"- (>{bins[-1]}): {bins_count['inf']}")

    outpath.write_text("\n".join(lines) + "\n")
    print(outpath.read_text())


# ---------------------------------------------------------------------
# Part 2: gap analysis on non-growers
# For each non-grower, for each biomass precursor (negative coefficient
# in bio1), add a sink and try to maximize. If max > epsilon, precursor
# is producible; otherwise it is blocked.

MEDIA = load_media()
EPS = 1e-6


def apply_media(model: cobra.Model):
    for rxn in model.reactions:
        if not rxn.id.startswith("EX_"):
            continue
        mets = list(rxn.metabolites.keys())
        if len(mets) != 1:
            continue
        cpd_id = mets[0].id.split("_")[0]
        rxn.lower_bound = -1000.0 if cpd_id in MEDIA else 0.0
        rxn.upper_bound = 1000.0


def gap_one(path_str):
    path = Path(path_str)
    try:
        model = load_json_model(str(path))
        apply_media(model)
        if "bio1" not in model.reactions:
            return path.stem, None, "no_bio1"
        bio = model.reactions.get_by_id("bio1")
        # Precursors = metabolites consumed (negative stoichiometry)
        precursors = [m.id for m, c in bio.metabolites.items() if c < 0]
        blocked = []
        producible = []
        # Probe each precursor by maximizing a temporary demand reaction
        # done with a context manager so changes are reverted per met
        for met_id in precursors:
            with model:
                met = model.metabolites.get_by_id(met_id)
                demand = model.add_boundary(met, type="demand",
                                            reaction_id=f"_probe_{met_id}",
                                            ub=1000.0, lb=0.0)
                model.objective = demand
                sol = model.optimize()
                if sol.status == "optimal" and (sol.objective_value or 0) > EPS:
                    producible.append(met_id)
                else:
                    blocked.append(met_id)
        return path.stem, {"precursors": precursors,
                           "producible": producible,
                           "blocked": blocked}, None
    except Exception as e:
        return path.stem, None, f"{type(e).__name__}: {e}"


def part2_gap_analysis(rows, reports_dir, results_dir=None):
    if results_dir is None:
        results_dir = reports_dir
    nongrowers = [r for r in rows if not r["grows"]]
    paths = [str(MODELS_DIR / f"{r['model_id']}.json") for r in nongrowers]
    print(f"Gap analysis on {len(paths)} non-growers...")

    blocked_counter = Counter()
    per_model = {}
    n_procs = max(1, min(mp.cpu_count() - 1, 16))
    with mp.Pool(n_procs) as pool:
        for i, (mid, res, err) in enumerate(
                pool.imap_unordered(gap_one, paths, chunksize=4), 1):
            if err or res is None:
                per_model[mid] = {"error": err}
                continue
            per_model[mid] = res
            for met in res["blocked"]:
                blocked_counter[met.split("_")[0]] += 1
            if i % 200 == 0:
                print(f"  {i}/{len(paths)}")

    (results_dir / "gap_per_model.json").write_text(json.dumps(per_model, indent=2))

    # Top blocked precursors
    lines = ["# Part 2 — Gap analysis on non-growers\n"]
    lines.append(f"Non-growers probed: {len(nongrowers)}\n")
    lines.append("Each model: every biomass precursor probed with a demand reaction "
                 "under complete media. 'Blocked' = max sink flux ≤ 1e-6.\n")
    lines.append("## Most commonly blocked biomass precursors")
    lines.append("| compound | count of models | % non-growers |")
    lines.append("|---|---|---|")
    for cpd, c in blocked_counter.most_common(30):
        lines.append(f"| {cpd} | {c} | {100*c/len(nongrowers):.1f}% |")

    # Per-model: how many precursors blocked
    counts = [len(v.get("blocked", [])) for v in per_model.values() if "blocked" in v]
    if counts:
        s = sorted(counts)
        lines.append("\n## Number of blocked precursors per non-grower")
        lines.append(f"- min: {s[0]}  median: {s[len(s)//2]}  max: {s[-1]}  mean: {mean(s):.2f}")
        bins = Counter()
        for c in counts:
            bins[c] += 1
        lines.append("\nHistogram (blocked precursors → model count):")
        for c in sorted(bins):
            lines.append(f"- {c}: {bins[c]}")

    (reports_dir / "GAP_ANALYSIS.md").write_text("\n".join(lines) + "\n")
    print((reports_dir / "GAP_ANALYSIS.md").read_text())
    return per_model, blocked_counter


# ---------------------------------------------------------------------
# Part 3: reaction prevalence (growers vs non-growers)
def extract_reactions_one(path_str):
    """Return (model_id, set of seed.reaction IDs)."""
    path = Path(path_str)
    try:
        with open(path) as f:
            m = json.load(f)
        rxns = set()
        for r in m["reactions"]:
            sr = r.get("annotation", {}).get("seed.reaction")
            if sr:
                rxns.add(sr)
        return path.stem, rxns
    except Exception:
        return path.stem, set()


def part3_reaction_prevalence(rows, outpath):
    growers = {r["model_id"] for r in rows if r["grows"]}
    nongrowers = {r["model_id"] for r in rows if not r["grows"]}
    paths = [str(MODELS_DIR / f"{r['model_id']}.json") for r in rows]
    print(f"Extracting reactions from {len(paths)} models...")

    g_count = Counter()
    n_count = Counter()
    n_procs = max(1, min(mp.cpu_count() - 1, 16))
    with mp.Pool(n_procs) as pool:
        for i, (mid, rxns) in enumerate(
                pool.imap_unordered(extract_reactions_one, paths, chunksize=8), 1):
            tgt = g_count if mid in growers else n_count
            for r in rxns:
                tgt[r] += 1
            if i % 1000 == 0:
                print(f"  {i}/{len(paths)}")

    G = len(growers); N = len(nongrowers)
    all_rxns = set(g_count) | set(n_count)

    # Per reaction: fraction in growers vs non-growers, and difference
    scored = []
    for r in all_rxns:
        g_frac = g_count[r] / G if G else 0
        n_frac = n_count[r] / N if N else 0
        scored.append((r, g_count[r], g_frac, n_count[r], n_frac, g_frac - n_frac))

    # Reactions most enriched in growers (g_frac high, n_frac low)
    enriched_in_growers = sorted(scored, key=lambda x: -x[5])[:30]
    enriched_in_nongrowers = sorted(scored, key=lambda x: x[5])[:30]

    lines = ["# Part 3 — Reaction prevalence: growers vs non-growers\n"]
    lines.append(f"Growers: {G}    Non-growers: {N}\n")
    lines.append("Each row is a `seed.reaction` ID; fractions = share of models containing it.\n")

    lines.append("## Top 30 reactions most enriched in GROWERS")
    lines.append("| seed.reaction | grower count | grower % | non-grower count | non-grower % | Δ |")
    lines.append("|---|---|---|---|---|---|")
    for r, gc, gf, nc, nf, d in enriched_in_growers:
        lines.append(f"| {r} | {gc} | {100*gf:.1f}% | {nc} | {100*nf:.1f}% | {100*d:+.1f}% |")

    lines.append("\n## Top 30 reactions most enriched in NON-GROWERS")
    lines.append("| seed.reaction | grower count | grower % | non-grower count | non-grower % | Δ |")
    lines.append("|---|---|---|---|---|---|")
    for r, gc, gf, nc, nf, d in enriched_in_nongrowers:
        lines.append(f"| {r} | {gc} | {100*gf:.1f}% | {nc} | {100*nf:.1f}% | {100*d:+.1f}% |")

    outpath.write_text("\n".join(lines) + "\n")
    print(outpath.read_text())


# ---------------------------------------------------------------------
def main():
    rows = load_results()
    print(f"Loaded {len(rows)} results")
    part1_characteristics(rows, REPORTS_DIR / "CHARACTERISTICS.md")
    # part2 writes both the per-model JSON and the GAP_ANALYSIS.md;
    # split the two outputs across results/ and reports/ by passing
    # the directory and patching the writer below.
    part2_gap_analysis(rows, REPORTS_DIR, RESULTS_DIR)
    part3_reaction_prevalence(rows, REPORTS_DIR / "REACTION_PREVALENCE.md")


if __name__ == "__main__":
    main()
