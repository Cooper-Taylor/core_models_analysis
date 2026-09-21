"""Tests for beginPipeline's input resolution and stage plumbing.

Run:  python3 -m pytest tests/ -q
      (needs the project interpreter; `pip install pytest` if absent)

Most of these encode defects found by an adversarial review of the first
version of this code. Each one names the wrong behaviour it guards against.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from cma import pipeline, resolve  # noqa: E402
from cma.directions import normalize_operator  # noqa: E402


# ---------------------------------------------------------------------------
# direction operators
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("raw,want", [
    (">", ">"), ("<", "<"), ("=", "="), ("?", "?"),
    ("forward", ">"), ("FWD", ">"), ("reverse", "<"), ("reversible", "="),
    # Arrow notation. The first version took the first character, so "<=>"
    # (reversible) became "<" (reverse-only) and "=>" (forward) became "="
    # -- silently inverting or loosening real constraints.
    ("<=>", "="), ("<->", "="), ("=>", ">"), ("->", ">"), ("-->", ">"),
    ("<=", "<"), ("<-", "<"), ("<--", "<"),
])
def test_operator_normalisation(raw, want):
    assert normalize_operator(raw) == want


def test_unrecognised_operator_is_rejected_not_guessed():
    """An unparseable value must raise, not silently become '?'.

    '?' is applied to the model as (-1000, 1000), i.e. fully reversible, so a
    typo used to quietly remove a constraint instead of failing.
    """
    with pytest.raises(ValueError):
        normalize_operator("banana", strict=True)
    assert normalize_operator("banana") == "?"     # lenient default still works


# ---------------------------------------------------------------------------
# model-set resolution
# ---------------------------------------------------------------------------
def test_registry_key_resolves():
    ms = resolve.resolve_model_set("core_kegg2")
    assert ms.name == "core_kegg2"


def test_directory_resolves_and_detects_loader(tmp_path):
    (tmp_path / "a.json").write_text("{}")
    (tmp_path / "b.json").write_text("{}")
    ms = resolve.resolve_model_set(str(tmp_path))
    assert ms.loader == "cobra_json"
    assert ms.id_glob == "*.json"


def test_gzip_directory_detects_gz_loader(tmp_path):
    (tmp_path / "a.json.gz").write_bytes(b"")
    ms = resolve.resolve_model_set(str(tmp_path))
    assert ms.loader == "cobra_json_gz"
    assert ms.id_glob == "*.json.gz"


def test_id_list_keeps_assembly_version(tmp_path):
    """`GCF_000005825.2` must not lose its `.2`.

    Path.suffixes treats a version number as an extension, so the first
    version resolved every id to a file that does not exist.
    """
    lst = tmp_path / "ids.txt"
    lst.write_text("GCF_000005825.2\nGCF_000005845.2\n")
    ms = resolve.resolve_model_set(str(lst))
    ids = resolve.model_ids(ms)
    assert ids == ["GCF_000005825.2", "GCF_000005845.2"]
    for i in ids:
        assert ms.path_for(i).exists(), f"{ms.path_for(i)} should exist"


def test_unknown_model_spec_raises_resolve_error():
    with pytest.raises(resolve.ResolveError):
        resolve.resolve_model_set("/definitely/not/here")


def test_mixed_format_directory_is_reported(tmp_path):
    """A directory with both .json and .json.gz must not silently drop one."""
    (tmp_path / "a.json").write_text("{}")
    (tmp_path / "b.json.gz").write_bytes(b"")
    (tmp_path / "c.json.gz").write_bytes(b"")
    with pytest.raises(resolve.ResolveError) as e:
        resolve.resolve_model_set(str(tmp_path))
    assert "json" in str(e.value)


def test_sbml_directory_raises_rather_than_failing_later(tmp_path):
    """SBML is detected but has no loader; fail at resolve, not mid-run."""
    (tmp_path / "m.xml").write_text("<sbml/>")
    with pytest.raises(resolve.ResolveError):
        resolve.resolve_model_set(str(tmp_path))


# ---------------------------------------------------------------------------
# direction-map resolution
# ---------------------------------------------------------------------------
def test_json_direction_map(tmp_path):
    p = tmp_path / "d.json"
    p.write_text(json.dumps({"rxn00001": ">", "rxn00002": "<=>"}))
    label, ops = resolve.resolve_directions(str(p))
    assert ops == {"rxn00001": ">", "rxn00002": "="}


@pytest.mark.parametrize("header", [
    "rxn_id,reversibility", "id,direction", "reaction,direction",
    "reaction_id,operator", "rxn_id,new_rev",
])
def test_csv_direction_headers(tmp_path, header):
    """Every column spelling the docs promise must actually parse."""
    p = tmp_path / "d.csv"
    p.write_text(f"{header}\nrxn00001,>\nrxn00002,=\n")
    _, ops = resolve.resolve_directions(str(p))
    assert ops == {"rxn00001": ">", "rxn00002": "="}


def test_empty_direction_table_raises(tmp_path):
    p = tmp_path / "d.json"
    p.write_text("{}")
    with pytest.raises(resolve.ResolveError):
        resolve.resolve_directions(str(p))


# ---------------------------------------------------------------------------
# heuristics
# ---------------------------------------------------------------------------
def test_inline_knobs():
    label, cfg = resolve.resolve_heuristic("sigma_band_k=1.96,mm_band=3")
    assert cfg.sigma_band_k == 1.96
    assert cfg.mm_band == 3


def test_registered_variant_tag():
    label, cfg = resolve.resolve_heuristic("baseline")
    assert label == "baseline"


def test_named_concentration_table_is_not_treated_as_a_path(tmp_path):
    p = tmp_path / "h.yaml"
    p.write_text("tag: b\nknobs:\n  per_met_conc_range: BENNETT_2009_ECOLI\n")
    _, cfg = resolve.resolve_heuristic(str(p))
    assert isinstance(cfg.per_met_conc_range, dict) and cfg.per_met_conc_range


def test_unknown_knob_raises_with_the_valid_list():
    with pytest.raises(resolve.ResolveError) as e:
        resolve.resolve_heuristic("not_a_knob=1")
    assert "sigma_band_k" in str(e.value)


# ---------------------------------------------------------------------------
# thermo
# ---------------------------------------------------------------------------
def test_thermo_sentinel_rows_dropped(tmp_path):
    """ModelSEED writes 10000000 for "no estimate"; it must never reach a statistic."""
    p = tmp_path / "t.tsv"
    p.write_text("rxn00001\t-3.5\t1.2\nrxn00002\t10000000\t10000000\n")
    _, table = resolve.resolve_thermo(str(p))
    assert table == {"rxn00001": (-3.5, 1.2)}


def test_unparseable_thermo_table_raises(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text("not a table at all\n")
    with pytest.raises(resolve.ResolveError):
        resolve.resolve_thermo(str(p))


# ---------------------------------------------------------------------------
# stage selection
# ---------------------------------------------------------------------------
def test_no_heavy_really_excludes_heavy_stages():
    """--no-heavy must not be undone by dependency expansion."""
    chosen = [s.name for s in pipeline.select(None, include_heavy=False)]
    assert "growth" not in chosen
    assert "reaction_effects" not in chosen


def test_exclude_wins_over_dependency_expansion():
    """--exclude growth must drop growth even though growth_diff needs it."""
    chosen = [s.name for s in pipeline.select(["growth_diff"], exclude=["growth"])]
    assert "growth" not in chosen


def test_only_pulls_in_dependencies():
    chosen = [s.name for s in pipeline.select(["growth_diff"])]
    assert "growth" in chosen
    assert chosen.index("growth") < chosen.index("growth_diff")


def test_unknown_stage_raises():
    from cma.registry import RegistryError

    with pytest.raises(RegistryError):
        pipeline.select(["no_such_stage"])


def test_unknown_tag_raises_rather_than_selecting_nothing():
    with pytest.raises(ValueError):
        pipeline.select(["tag:nonexistent"])


def test_tag_selection_works():
    chosen = [s.name for s in pipeline.select(["tag:inventory"])]
    assert "inputs" in chosen


# ---------------------------------------------------------------------------
# labelling
# ---------------------------------------------------------------------------
def test_duplicate_direction_labels_are_disambiguated(tmp_path):
    """Two files with the same basename must not collapse into one column."""
    a, b = tmp_path / "x" / "d.json", tmp_path / "y" / "d.json"
    for p in (a, b):
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps({"rxn00001": ">"}))
    la, _ = resolve.resolve_directions(str(a))
    lb, _ = resolve.resolve_directions(str(b))
    labels = resolve.disambiguate([la, lb], [str(a), str(b)])
    assert len(set(labels)) == 2


def test_on_disk_is_a_reserved_direction_label(tmp_path):
    """A user's file named on_disk.json must not replace the baseline column."""
    p = tmp_path / "on_disk.json"
    p.write_text(json.dumps({"rxn00001": ">"}))
    label, _ = resolve.resolve_directions(str(p))
    assert label != "on_disk"
