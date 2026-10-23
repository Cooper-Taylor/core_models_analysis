"""Paths for the tFBA investigation scripts.

Follows the environment-variable convention already used by
`scripts/analyze_growth.py`: `MSDB_ROOT` for the ModelSEEDDatabase checkout and
`CORE_MODELS_ANALYSIS_DIR` for this repo, each with a sensible default.

Kept deliberately free of cobra/pandas/numpy imports so it stays cheap.
"""
import os
from pathlib import Path

REPO = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                           Path(__file__).resolve().parents[2]))
MSDB = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/ModelSEEDDatabase"))

BIOCHEM = MSDB / "Biochemistry"
MEDIA = MSDB / "Media" / "KBaseMedia.cpd"

MODELS = Path(os.environ.get("CORE_MODELS_DIR", REPO / "data" / "core_models_kegg2"))
DEFAULT_MODEL = MODELS / "GCF_000005845.2.json"

RESULTS = Path(os.environ.get("TFBA_RESULTS_DIR", REPO / "results" / "tfba"))

# Rebuilt pre-2022-05-20 ModelSEEDDatabase layout; see reports/THERMODYNAMIC_FBA.md
# and build_legacy_msdb.sh. Only needed to run the stock (unshimmed) package.
LEGACY_MSDB = Path(os.environ.get("TFBA_LEGACY_MSDB", RESULTS / "legacy_msdb"))
