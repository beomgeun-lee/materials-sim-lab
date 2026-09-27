"""스키마 v0 — 수량, 성분 참조, 레시피, 결과."""

from __future__ import annotations

import datetime as dt
from pathlib import Path

import pytest
from pydantic import ValidationError

from msl.schema.provenance import DataProvenance, RunProvenance, canonical_hash
from msl.schema.quantity import Dimension, Quantity
from msl.schema.recipe import Recipe, RecipeError, load_recipe
from msl.schema.refs import SubstanceRef, cas_checksum_ok
from msl.schema.result import AssayResult, Fidelity, ResultValue, Status, ValueKind

EXAMPLES = Path(__file__).resolve().parents[1] / "examples" / "recipes"

# ── 수량 ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    ("text", "dim", "si"),
    [
        ("1,000 °C", Dimension.TEMPERATURE, 1273.15),
        ("1273 K", Dimension.TEMPERATURE, 1273.0),
        ("1 atm", Dimension.PRESSURE, 101_325.0),
        ("0.1 mol/L", Dimension.CONCENTRATION, 100.0),
        ("30 vol%", Dimension.VOLUME_FRACTION, 0.3),
        ("250 mg", Dimension.MASS, 2.5e-4),
        ("1e-3 mol", Dimension.AMOUNT, 1e-3),
    ],
)
def test_quantity_parse(text: str, dim: Dimension, si: float) -> None:
    q = Quantity.model_validate(text)
    assert q.dimension is dim
    assert q.to_si() == pytest.approx(si)


@pytest.mark.parametrize("text", ["1000 furlongs", "°C", "abc K", ""])
def test_quantity_rejects(text: str) -> None:
    with pytest.raises(ValidationError):
        Quantity.model_validate(text)


# ── 성분 참조 ─────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "ref",
    [
        "element:Fe",
        "formula:Fe2O3",
        "formula:Ca(OH)2",
        "formula:CuSO4·5H2O",
        "formula:Li0.5CoO2",
        "mineral:monazite-(Ce)",
        "cas:7732-18-5",
        "mp:mp-19770",
        "cod:9000660",
        "inchikey:XLYOFNOQVPJJNP-UHFFFAOYSA-N",
        "name:표백제",
        "material:e-glass-fiber",
    ],
)
def test_ref_accepts(ref: str) -> None:
    assert str(SubstanceRef.model_validate(ref)) == ref


@pytest.mark.parametrize(
    "ref",
    [
        "Fe",  # 네임스페이스 없음
        "element:Xx",
        "formula:Fe2Q3",  # Q 는 원소가 아님
        "formula:Ca(OH2",  # 괄호 짝
        "cas:7732-18-4",  # 검사 숫자 불일치
        "mp:19770",
        "planet:earth",
    ],
)
def test_ref_rejects(ref: str) -> None:
    with pytest.raises(ValidationError):
        SubstanceRef.model_validate(ref)


def test_cas_checksum() -> None:
    assert cas_checksum_ok("7647-01-0")  # HCl
    assert cas_checksum_ok("7681-52-9")  # NaOCl
    assert not cas_checksum_ok("7681-52-8")


# ── 레시피 ────────────────────────────────────────────────────────────────


def test_examples_are_valid() -> None:
    paths = sorted(EXAMPLES.glob("*.yaml"))
    assert len(paths) >= 5  # 데모 레시피 (계획서 1절 5종 + 실행 데모)
    for path in paths:
        load_recipe(path)


def _recipe(**overrides: object) -> dict[str, object]:
    base: dict[str, object] = {
        "id": "rcp-test",
        "name": "시험",
        "components": [
            {"ref": "formula:MgO", "amount": "1 mol"},
            {"ref": "formula:Al2O3", "amount": "1 mol"},
        ],
    }
    return base | overrides


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"id": "test"}, "String should match pattern"),
        ({"components": [{"ref": "formula:MgO"}]}, "모든 성분에 양이 필요"),
        ({"components": [{"ref": "formula:MgO", "amount": "-1 mol"}]}, "0보다 커야"),
        ({"components": [{"ref": "formula:MgO", "amount": "1 atm"}]}, "쓸 수 없는 단위"),
        (
            {"components": [{"ref": "formula:MgO", "amount": "1 mol"}, {"ref": "formula:MgO", "amount": "2 mol"}]},
            "두 번",
        ),
        (
            {"components": [{"ref": "element:Cu", "amount": "60 at%"}, {"ref": "element:Ni", "amount": "30 at%"}]},
            "100%",
        ),
        (
            {"components": [{"ref": "element:Cu", "amount": "70 at%"}, {"ref": "element:Ni", "amount": "30 wt%"}]},
            "같은 종류의 분율",
        ),
        ({"conditions": {"T": "-300 °C"}}, "0 K 보다"),
        ({"conditions": {"T": "1 atm"}}, "온도에 쓸 수 없는"),
        ({"assays": ["S0", "B1"]}, "시험 코드 형식"),
        ({"mode": "system"}, "원소만"),
        ({"sweep": {"param": "ratio", "start": 0.0, "stop": 1.5, "steps": 5}}, "0~1"),
        ({"sweep": {"param": "T", "start": 0.1, "stop": 0.9, "steps": 5}}, "단위가 있는 수량"),
        ({"extra_field": 1}, "Extra inputs"),
    ],
)
def test_recipe_rejects(overrides: dict[str, object], message: str) -> None:
    with pytest.raises(ValidationError, match=message):
        Recipe.model_validate(_recipe(**overrides))


def test_recipe_unknown_assay_against_catalog(tmp_path: Path) -> None:
    path = tmp_path / "r.yaml"
    path.write_text("id: rcp-x\nname: x\ncomponents: [{ref: 'element:Fe', amount: 1 mol}]\nassays: [A99]\n")
    with pytest.raises(RecipeError, match="A99"):
        load_recipe(path, known_assays={"S0", "A1"})


# ── 결과·프로비넌스 ───────────────────────────────────────────────────────


def _run() -> RunProvenance:
    now = dt.datetime(2026, 9, 26, tzinfo=dt.UTC)
    return RunProvenance(
        input_hash="0" * 64, msl_version="0.0.1", engine_versions={}, started_at=now, finished_at=now
    )


def _result(**overrides: object) -> dict[str, object]:
    base: dict[str, object] = {
        "assay": "A2",
        "recipe_id": "rcp-test",
        "status": Status.OK,
        "fidelity": Fidelity.L0,
        "engine": "pymatgen",
        "engine_version": "2026.9.24",
        "conditions_basis": "0 K DFT",
        "values": [ResultValue(name="반응에너지", kind=ValueKind.ENERGY, value=-0.1, unit="eV/atom")],
        "energy_reference": "MP2020",
        "provenance": _run(),
    }
    return base | overrides


def test_result_requires_energy_reference() -> None:
    AssayResult.model_validate(_result())
    with pytest.raises(ValidationError, match="energy_reference"):
        AssayResult.model_validate(_result(energy_reference=None))


def test_result_l2_requires_weights() -> None:
    with pytest.raises(ValidationError, match="model_weights"):
        AssayResult.model_validate(_result(fidelity=Fidelity.L2))
    AssayResult.model_validate(_result(fidelity=Fidelity.L2, model_weights="orb-v3@sha256:abc"))


def test_result_ok_needs_values() -> None:
    with pytest.raises(ValidationError, match="값이 하나 이상"):
        AssayResult.model_validate(_result(values=[], energy_reference=None))
    AssayResult.model_validate(_result(status=Status.NOT_APPLICABLE, values=[], energy_reference=None))


def test_canonical_hash_ignores_key_order() -> None:
    assert canonical_hash({"a": 1, "b": [1, 2]}) == canonical_hash({"b": [1, 2], "a": 1})
    assert canonical_hash({"a": 1}) != canonical_hash({"a": 2})
    r1 = Recipe.model_validate(_recipe())
    r2 = Recipe.model_validate(_recipe())
    assert canonical_hash(r1) == canonical_hash(r2)


def test_data_provenance_method() -> None:
    fields = {
        "source": "materials-project",
        "source_id": "mp-19770",
        "source_version": "v2026.04.13",
        "retrieved_at": dt.datetime(2026, 9, 26, tzinfo=dt.UTC),
        "license": "CC-BY-4.0",
    }
    DataProvenance(**fields, method="DFT-PBE+U", correction_scheme="MP2020")
    with pytest.raises(ValidationError):
        DataProvenance(**fields, method="guess")
