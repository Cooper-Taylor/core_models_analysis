"""Model sets, media and panels.

Two collections are registered:

``core_kegg2``
    the 5,683 KEGG2-derived core models this project has always used.
    ``results_namespace=None`` keeps every existing output path byte-identical.

``ms2_gsm`` / ``ms2_gsm_auxo``
    the 5,420 genome-scale models from the ModelSEED v2 manuscript
    (Faria et al. 2023), downloaded from the public KBase workspaces 155807 and
    155808 and converted to cobrapy JSON. Measured differences from the core
    set, each of which is a declared field below rather than a surprise:

    =========================  ==================  ==========================
    property                   core_kegg2          ms2_gsm
    =========================  ==================  ==========================
    reactions / model (median) 128                 1,085
    file form                  ``.json``           ``.json.gz``
    model id                   ``GCF_000005845.2`` ``GCF_000005845.2.RAST.GMM``
    biomass                    ``bio1`` + ``bio2``  ``bio1`` only, 51-61 mets
    medium                     347-compound        20-compound defined minimal
    distinct MSDB reactions    239                 3,634
    =========================  ==================  ==========================

    The biomass difference is the dangerous one. ``flux_loops_one`` closes every
    ``bio*`` reaction carrying at least ``real_biomass_min_mets`` metabolites and
    relies on a smaller one staying open as the ATP probe. The core set has
    ``bio2`` (5 metabolites) to play that part; the genome-scale set has no
    second biomass and no maintenance reaction, so everything gets closed and
    energy-generating-cycle detection returns zero for every model, with no
    error. ``atp_probe=None`` records that a probe must be *injected*.
"""

from __future__ import annotations

import gzip
import json
import os
from pathlib import Path

from .. import paths
from ..kinds import BiomassPolicy, Media, ModelSet, Panel, SeedKeyPolicy
from ..models import register, register_media, register_panel

GSM_ROOT = Path(os.environ.get("MS2_GS_MODELS_DIR", "/scratch/ctaylor/modelseed2_gs_models"))


# --- media -----------------------------------------------------------------
def _kbase_complete() -> dict:
    """The 347-compound complete medium, opened at the legacy -1000 bound."""
    p = paths.msdb("live") / "Media" / "KBaseMedia.cpd"
    return {line.strip(): -1000.0 for line in p.read_text().splitlines() if line.strip()}


def _cpd_file(path: Path, uptake: float) -> dict:
    return {c: uptake for c in path.read_text().split() if c}


register_media(Media(
    name="kbase_complete",
    loader=_kbase_complete,
    default_uptake=-1000.0,
    description="ModelSEEDDatabase/Media/KBaseMedia.cpd, 347 compounds, uptake -1000.",
))

register_media(Media(
    name="gsm_glucose_minimal",
    loader=lambda: _cpd_file(GSM_ROOT / "gmm" / "media" / "Carbon-D-Glucose.cpd", -10.0),
    default_uptake=-10.0,
    description="The 20-compound defined minimal medium the ms2_gsm models were gap-filled on.",
))

register_media(Media(
    name="gsm_auxotrophy",
    loader=lambda: _cpd_file(GSM_ROOT / "auxotrophy" / "media" / "Auxotrophy_media.cpd", -10.0),
    default_uptake=-10.0,
    description="The 51-compound auxotrophy medium of KBase workspace 155808.",
))


# --- loaders ---------------------------------------------------------------
def load_cobra_gz(path: Path):
    """cobra.io.load_json_model cannot read a gzipped file; go via the dict."""
    from cobra.io import model_from_dict

    with gzip.open(path, "rt", encoding="utf-8") as fh:
        return model_from_dict(json.load(fh))


# --- the core set ----------------------------------------------------------
register(ModelSet(
    name="core_kegg2",
    models_dir=lambda: paths.analysis_dir() / "data" / "core_models_kegg2",
    id_glob="*.json",
    loader="cobra_json",
    media="kbase_complete",
    biomass=BiomassPolicy(
        prefer=("bio1",),
        then=("bio2", "biomass", "Biomass"),
        real_biomass_min_mets=10,
        atp_probe=None,          # bio2 (5 mets) is left open and serves as the probe
    ),
    seed_key=SeedKeyPolicy(
        annotation="seed.reaction",
        strip_compartment_suffix=True,   # 17 transport reactions carry a stray _c
        allow_list_values=False,
        fallback=None,
        require_nonzero_overrides=True,
    ),
    compartments=("c0", "e0"),
    cytosol="c0",
    id_to_accession=lambda mid: mid,
    results_namespace=None,      # keeps results/ exactly as it is today
    default_panel="core_100",
    description="5,683 KEGG2-derived ModelSEED core models.",
))

register_panel(Panel(
    name="core_100", model_set="core_kegg2",
    path=lambda: paths.results("selected_ids.txt"), expected_n=100,
    description="Genome-similarity-based diverse panel (reports/DIVERSE_SELECTION.md).",
))
register_panel(Panel(
    name="core_100_tax", model_set="core_kegg2",
    path=lambda: paths.results("selected_ids_tax.txt"), expected_n=100,
    description="Taxonomy-aware panel (reports/TAXONOMY_AWARE_SELECTION.md).",
))


# --- the ModelSEED v2 genome-scale sets ------------------------------------
_GSM_BIOMASS = BiomassPolicy(
    prefer=("bio1",),
    then=(),
    real_biomass_min_mets=10,
    atp_probe=None,   # no second biomass, no maintenance reaction: inject one
)
_GSM_SEED_KEY = SeedKeyPolicy(
    annotation="seed.reaction",
    strip_compartment_suffix=True,
    allow_list_values=True,       # cheap insurance; measured str-only today
    fallback="id_prefix",         # rxn02201_c0 -> rxn02201
    require_nonzero_overrides=True,
)

for _slug, _media_name, _desc in (
    ("gmm", "gsm_glucose_minimal", "glucose minimal media (KBase ws 155807)"),
    ("auxotrophy", "gsm_auxotrophy", "auxotrophy media (KBase ws 155808)"),
):
    _name = "ms2_gsm" if _slug == "gmm" else "ms2_gsm_auxo"
    register(ModelSet(
        name=_name,
        models_dir=(lambda s=_slug: GSM_ROOT / s / "cobra"),
        id_glob="*.json.gz",
        loader="cma.entries.model_sets:load_cobra_gz",
        media=_media_name,
        biomass=_GSM_BIOMASS,
        seed_key=_GSM_SEED_KEY,
        compartments=("c0", "e0"),
        cytosol="c0",
        # GCF_000005845.2.RAST.GMM -> GCF_000005845.2, the taxonomy join key
        id_to_accession=lambda mid: mid.split(".RAST.")[0],
        results_namespace=_name,
        default_panel=None,
        description=f"5,420 ModelSEED v2 genome-scale models, {_desc}.",
    ))
