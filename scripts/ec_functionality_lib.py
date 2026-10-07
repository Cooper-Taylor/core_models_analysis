"""Shared helpers for the enzyme-functionality (mono/bi/polyfunctional)
classification pipeline (scripts/build_enzyme_functionality_queue.py,
scripts/enzyme_functionality_batch.py, scripts/build_enzyme_functionality_outputs.py).

Reuses the exact "comparable reaction" filter from
build_llm_vs_recommended_mismatches.py so this pipeline's reaction universe
matches reports/llmVsRecommended/LLM_VS_RECOMMENDED_DIRECTION.md (20,490
reactions) rather than drifting from it.
"""
from __future__ import annotations

import glob
import json
import os
import re
from pathlib import Path

MSDB_ROOT = Path(os.environ.get("MSDB_ROOT", "/scratch/ctaylor/tmp/devsnap_078a395f"))
MSDB_SHA = "078a395f"
ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
BIOCHEM = MSDB_ROOT / "Biochemistry"
OUT_DIR = ANALYSIS_DIR / "results" / "enzyme_functionality"

DIRECTIONS = (">", "=", "<")

COMPLETE_EC_RE = re.compile(r"^\d+\.\d+\.\d+\.\d+$")

CATEGORIES = ("monofunctional", "bifunctional", "polyfunctional", "incomplete_ec")
CATEGORY_CODE = {"A": "monofunctional", "B": "bifunctional", "C": "polyfunctional"}
CATEGORY_RANK = {"monofunctional": 1, "bifunctional": 2, "polyfunctional": 3}


def is_complete_ec(ec: str) -> bool:
    return bool(COMPLETE_EC_RE.match(ec))


def iter_comparable_reactions():
    """Yield every reaction in the same 20,490-reaction comparable set as
    build_llm_vs_recommended_mismatches.py: non-obsolete, recommended
    reversibility and LLM-ensemble direction both present and directional.
    """
    for path in sorted(glob.glob(str(BIOCHEM / "reaction_*.json"))):
        for rxn in json.load(open(path)):
            if rxn.get("is_obsolete"):
                continue
            rec = rxn.get("reversibility")
            th = rxn.get("thermodynamics") or {}
            llm = th.get("LLMs")
            if llm is None:
                continue
            llm_dir = llm[2]
            if llm_dir == "?":
                continue
            if rec not in DIRECTIONS:
                continue
            yield rxn


def ensure_out_dir() -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    return OUT_DIR
