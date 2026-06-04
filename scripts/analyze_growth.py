#!/usr/bin/env python3
"""
Analyze biological growth potential of core_models_kegg2 metabolic models.

For each .json model:
- Load with cobrapy
- Constrain EX_ exchange reactions so uptake (negative flux) is only allowed
  for compounds listed in the ModelSEEDDatabase KBaseMedia.cpd complete media
- Run FBA on the biomass reaction (bio1)
- Record growth flux, status, and active biomass reaction

Outputs:
- results.csv  : per-model growth summary
- failures.log : per-model errors / non-optimal solutions
"""

import csv
import json
import logging
import multiprocessing as mp
import os
import sys
import traceback
from pathlib import Path

import cobra
from cobra.io import load_json_model

# --- Paths ------------------------------------------------------------
ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR", "/scratch/ctaylor/core_models_analysis"))
MODELS_DIR = ANALYSIS_DIR / "data" / "core_models_kegg2"
MEDIA_FILE = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase") + "/Media/KBaseMedia.cpd")
RESULTS_CSV = ANALYSIS_DIR / "results" / "results.csv"
FAILURES_LOG = ANALYSIS_DIR / "logs" / "failures.log"

# --- Media ------------------------------------------------------------
def load_media_compounds(path: Path) -> set:
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}

MEDIA_COMPOUNDS = load_media_compounds(MEDIA_FILE)
GROWTH_THRESHOLD = 1e-6  # flux below this is treated as no growth


def apply_media(model: cobra.Model, media_cpds: set) -> int:
    """
    Restrict uptake: for every EX_ reaction, only allow negative (uptake) flux
    if the underlying compound id is in the media list. Secretion (positive flux)
    is always allowed.
    Returns the number of exchanges left open for uptake.
    """
    open_count = 0
    for rxn in model.reactions:
        if not rxn.id.startswith("EX_"):
            continue
        # Single-metabolite exchanges only
        mets = list(rxn.metabolites.keys())
        if len(mets) != 1:
            continue
        met = mets[0]
        # cpd id is the leading token before the compartment suffix
        # e.g. cpd00067_e0 -> cpd00067
        cpd_id = met.id.split("_")[0]
        if cpd_id in media_cpds:
            rxn.lower_bound = -1000.0
            open_count += 1
        else:
            rxn.lower_bound = 0.0
        # always allow secretion
        if rxn.upper_bound < 1000.0:
            rxn.upper_bound = 1000.0
    return open_count


def find_biomass_reaction(model: cobra.Model):
    # Prefer bio1, then anything that looks like biomass
    if "bio1" in model.reactions:
        return model.reactions.get_by_id("bio1")
    for rid in ("bio2", "biomass", "Biomass"):
        if rid in model.reactions:
            return model.reactions.get_by_id(rid)
    for r in model.reactions:
        if r.id.lower().startswith("bio") and not r.id.startswith("SK_"):
            return r
    return None


def analyze_one(path_str: str) -> dict:
    path = Path(path_str)
    res = {
        "model_id": path.stem,
        "n_metabolites": 0,
        "n_reactions": 0,
        "n_genes": 0,
        "biomass_rxn": "",
        "n_exchanges_total": 0,
        "n_exchanges_open": 0,
        "status": "",
        "growth_flux": 0.0,
        "grows": False,
        "error": "",
    }
    try:
        model = load_json_model(str(path))
        res["n_metabolites"] = len(model.metabolites)
        res["n_reactions"] = len(model.reactions)
        res["n_genes"] = len(model.genes)
        res["n_exchanges_total"] = sum(1 for r in model.reactions if r.id.startswith("EX_"))

        res["n_exchanges_open"] = apply_media(model, MEDIA_COMPOUNDS)

        bio_rxn = find_biomass_reaction(model)
        if bio_rxn is None:
            res["status"] = "no_biomass"
            return res
        res["biomass_rxn"] = bio_rxn.id

        model.objective = bio_rxn
        sol = model.optimize()
        res["status"] = sol.status
        flux = float(sol.objective_value) if sol.objective_value is not None else 0.0
        res["growth_flux"] = flux
        res["grows"] = (sol.status == "optimal") and (flux > GROWTH_THRESHOLD)
    except Exception as e:
        res["error"] = f"{type(e).__name__}: {e}"
    return res


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--single":
        # debug helper: analyze one model and print
        print(json.dumps(analyze_one(sys.argv[2]), indent=2))
        return

    logging.basicConfig(
        filename=FAILURES_LOG,
        filemode="w",
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    model_paths = sorted(str(p) for p in MODELS_DIR.glob("*.json"))
    print(f"Found {len(model_paths)} models", flush=True)
    print(f"Media has {len(MEDIA_COMPOUNDS)} compounds", flush=True)

    # Silence cobra/optlang noise from per-model solve
    cobra_logger = logging.getLogger("cobra")
    cobra_logger.setLevel(logging.ERROR)

    fields = [
        "model_id", "n_metabolites", "n_reactions", "n_genes",
        "biomass_rxn", "n_exchanges_total", "n_exchanges_open",
        "status", "growth_flux", "grows", "error",
    ]

    n_proc = max(1, min(mp.cpu_count() - 1, 16))
    print(f"Running with {n_proc} worker processes", flush=True)

    done = 0
    grew = 0
    with open(RESULTS_CSV, "w", newline="") as csvfile, mp.Pool(n_proc) as pool:
        writer = csv.DictWriter(csvfile, fieldnames=fields)
        writer.writeheader()
        for res in pool.imap_unordered(analyze_one, model_paths, chunksize=4):
            writer.writerow(res)
            done += 1
            if res["grows"]:
                grew += 1
            if res.get("error"):
                logging.warning("%s: %s", res["model_id"], res["error"])
            elif res["status"] not in ("optimal", ""):
                logging.info("%s: status=%s flux=%.4g",
                             res["model_id"], res["status"], res["growth_flux"])
            if done % 100 == 0:
                csvfile.flush()
                print(f"  progress: {done}/{len(model_paths)}  growing: {grew}",
                      flush=True)

    print(f"\nDone. {done} models analyzed; {grew} grow above threshold {GROWTH_THRESHOLD}.")
    print(f"Results: {RESULTS_CSV}")
    print(f"Failures log: {FAILURES_LOG}")


if __name__ == "__main__":
    main()
