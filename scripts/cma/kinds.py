"""The plugin schemas. This file is the documentation an author reads.

Six kinds of thing can be added to this analysis: a **thermo source**, a
**heuristic variant** (or an individual **rule**), a **direction map**, a
**model set**, a **figure set**, and a **probe**. Each is a frozen dataclass
declared here and registered in ``cma/entries/``.

The design rule throughout: *every field that used to be an invisible
convention is a declared field*. A source's five different spellings, a
variant's "is this a cascade config or a precomputed map?" discriminator, a
model set's biomass-selection magic number -- all of them were previously
facts you had to know, scattered across dozens of scripts. They are data now.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

VALID_DIRECTIONS = frozenset("><=?")


# ── 1. THERMO SOURCE ────────────────────────────────────────────────────────
@dataclass(frozen=True)
class ThermoSource:
    """One estimator of reaction ΔG′° (or of directions directly).

    ``slugs`` carries every spelling in live use. Before this registry the same
    source appeared as ``gc`` / ``GC`` / ``group_contribution`` /
    ``group-contribution`` / ``Group contribution`` in different files, with no
    mapping table anywhere, so a partial edit silently dropped the source from
    one artifact.

    ``kind`` names the three genuinely incompatible ΔG mechanisms:

    ``msdb_sublist``
        read ``reaction['thermodynamics'][label]`` -> ``[dg, dge, operator]``
    ``msdb_toplevel_gated``
        read the top-level ``deltag``/``deltagerr``, with ``db_level`` only
        gating eligibility and stamping a label. Faithful to upstream
        ``top_level_energy()``; kept, not "fixed", and stamped in manifests.
    ``external_energy`` / ``external_operators`` / ``derived``
        a table off disk, or a composite of other sources.
    """

    key: str
    display: str
    slugs: dict[str, str]
    label: str | None = None
    kind: str = "msdb_sublist"
    db_level: str | None = None
    reaction_level: bool = True
    compound_level: bool = False
    incomplete_fallback: str | None = None
    energy: Callable | None = None       # lazy: () -> (rxn_entry) -> (dg, dge, label)
    operators: Callable | None = None    # lazy: () -> {rxn_id: operator}
    mask: Callable | None = None         # lazy: () -> set[rxn_id] to exclude
    color: str | None = None
    color_families: tuple[str, ...] = ()
    axis_title: str = ""
    data_ref: str = "msdb:live"
    layers: tuple[str, ...] = ("panel", "all_models")
    section: str = ""
    citations: tuple[str, ...] = ()

    def slug(self, style: str = "long") -> str:
        try:
            return self.slugs[style]
        except KeyError:
            raise KeyError(
                f"thermo source {self.key!r} has no {style!r} slug; "
                f"declared: {sorted(self.slugs)}"
            ) from None


# ── 2. HEURISTIC VARIANT (+ RULE) ───────────────────────────────────────────
@dataclass(frozen=True)
class Variant:
    """One reversibility variant.

    ``kind`` is the discriminator the old catalog lacked. Its ``cfg`` field was
    a *mandatory* zero-arg ``ReversibilityConfig`` factory, so any variant whose
    directions came from somewhere else -- an LLM, eQuilibrator 3, on-disk KEGG
    bounds, a consensus vote -- could not be registered at all. Those six lived
    only in ``thermo_variants/manifest.json``, which
    ``export_thermo_variants.py`` rebuilt from the catalog and therefore erased.
    """

    tag: str
    kind: str                     # 'cascade' | 'overlay' | 'derived'
    title: str = ""
    apt_title: str = ""
    description: str = ""
    citations: tuple[str, ...] = ()
    section: str = ""
    cfg: Callable | None = None           # required iff kind == 'cascade'
    direction_map: str | None = None      # required iff kind == 'overlay'
    overlay_on: str = "baseline"
    status_prefix: str | None = None
    combine: Callable | None = None       # required iff kind == 'derived'
    levels: tuple[str, ...] = ("EQ", "GC", "")
    order: int = 100
    in_notebook: bool = False
    in_stats: bool = False
    in_presentation: bool = False

    @property
    def key(self) -> str:
        return self.tag

    def legacy_dict(self) -> dict:
        """Exactly the seven keys the pre-registry catalog exposed."""
        return {
            "tag": self.tag,
            "title": self.title,
            "apt_title": self.apt_title,
            "description": self.description,
            "citations": list(self.citations),
            "section": self.section,
            "cfg": self.cfg,
        }


@dataclass(frozen=True)
class Rule:
    """One step of the reversibility cascade.

    Position is declared (``after`` / ``before``) instead of being typed into a
    100-line ``if``-chain. ``enabled_by`` names a ``ReversibilityConfig`` knob,
    so a rule that is inert at default settings cannot disturb baseline parity.
    """

    name: str
    fn: Callable
    after: str | None = None
    before: str | None = None
    enabled_by: str | None = None
    status_prefix: str = ""

    @property
    def key(self) -> str:
        return self.name


# ── 3. DIRECTION MAP ────────────────────────────────────────────────────────
@dataclass(frozen=True)
class DirectionMapSpec:
    """A ``{rxn_id: '>'|'<'|'='|'?'}`` table.

    ``coverage`` is declared, never inferred. ``complete`` means the map has an
    opinion (possibly ``?``) about every MSDB reaction; ``partial`` means an
    absent key is *no opinion* and the consumer must fall back. Conflating the
    two moves published growth numbers substantially, so the distinction is
    carried into every run manifest rather than normalized away.
    """

    key: str
    loader: Callable[[], dict]
    coverage: str = "partial"
    source: str | None = None
    data_ref: str = "msdb:live"
    masked: bool = False
    description: str = ""
    artifacts: dict[str, Path] = field(default_factory=dict)


# ── 4. MODEL SET ────────────────────────────────────────────────────────────
@dataclass(frozen=True)
class BiomassPolicy:
    """How to find the biomass reaction, and what counts as a real one.

    ``real_biomass_min_mets`` was a bare ``10`` inside ``flux_loops_one``: any
    ``bio*`` reaction with at least that many metabolites is closed, and a
    smaller one is left open to serve as the ATP probe. Core models have
    ``bio1`` (22 metabolites) and ``bio2`` (5), so that works by luck. Sets with
    a single large biomass and no maintenance reaction get everything closed and
    silently report zero energy-generating cycles -- hence ``atp_probe``.
    """

    prefer: tuple[str, ...] = ("bio1",)
    then: tuple[str, ...] = ("bio2", "biomass", "Biomass")
    fallback_prefix: str = "bio"
    exclude_prefix: str = "SK_"
    real_biomass_min_mets: int = 10
    atp_probe: str | None = None


@dataclass(frozen=True)
class SeedKeyPolicy:
    """How a model reaction is matched to an MSDB reaction id."""

    annotation: str = "seed.reaction"
    strip_compartment_suffix: bool = True
    allow_list_values: bool = False
    fallback: str | None = None            # 'id_prefix': rxn12345_c0 -> rxn12345
    require_nonzero_overrides: bool = True


@dataclass(frozen=True)
class Media:
    name: str
    loader: Callable[[], dict]             # {cpd_id: uptake_lower_bound}
    default_uptake: float = -1000.0
    description: str = ""

    @property
    def key(self) -> str:
        return self.name


@dataclass(frozen=True)
class ModelSet:
    """A collection of metabolic models the pipeline can sweep.

    ``results_namespace=None`` means "write to ``results/`` exactly as before".
    The core set keeps it, so adopting this registry changes no existing path.
    """

    name: str
    models_dir: Callable[[], Path]
    id_glob: str = "*.json"
    loader: str = "cobra_json"             # cobra_json | cobra_json_gz | module:fn
    media: str = "kbase_complete"
    biomass: BiomassPolicy = BiomassPolicy()
    seed_key: SeedKeyPolicy = SeedKeyPolicy()
    compartments: tuple[str, ...] = ("c0", "e0")
    cytosol: str = "c0"
    id_to_accession: Callable[[str], str] | None = None
    results_namespace: str | None = None
    default_panel: str | None = None
    description: str = ""

    @property
    def key(self) -> str:
        return self.name

    def model_ids(self) -> list[str]:
        d = self.models_dir()
        suffixes = "".join(Path(self.id_glob).suffixes)
        return sorted(p.name[: -len(suffixes)] if suffixes else p.stem
                      for p in d.glob(self.id_glob))

    def path_for(self, model_id: str) -> Path:
        suffixes = "".join(Path(self.id_glob).suffixes)
        return self.models_dir() / f"{model_id}{suffixes}"


@dataclass(frozen=True)
class Panel:
    name: str
    model_set: str
    path: Callable[[], Path]
    expected_n: int | None = None
    description: str = ""

    @property
    def key(self) -> str:
        return self.name

    def ids(self) -> list[str]:
        return self.path().read_text().split()


# ── 5. FIGURE SET ───────────────────────────────────────────────────────────
@dataclass(frozen=True)
class FigureSet:
    name: str
    script: str
    args: str = ""
    out_dir: str = ""
    tags: tuple[str, ...] = ()
    needs: tuple[str, ...] = ()
    values_sidecar: str | None = None
    palette_family: str | None = None

    @property
    def key(self) -> str:
        return self.name


# ── 6. PROBE ────────────────────────────────────────────────────────────────
@dataclass(frozen=True)
class Probe:
    """One per-model measurement run under a direction map.

    The same four-step prologue -- load, rebind bounds, apply media, maximize
    biomass -- is hand-written at least thirteen times across the drivers. A
    probe supplies only the part that differs.
    """

    name: str
    fn: Callable
    needs_media: bool = True
    needs_biomass: bool = True
    internal_pool: bool = False            # True => the runner refuses to nest it
    out_template: str = ""
    shard: bool = False
    description: str = ""

    @property
    def key(self) -> str:
        return self.name
