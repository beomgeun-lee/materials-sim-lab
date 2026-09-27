"""2단계 — 질량수지, A10 함량 기준, 검증 채점, 리포트 (네트워크 없이)."""

from __future__ import annotations

import datetime as dt

import pytest

from msl.assays.regulation import _applies, _classes, _judge
from msl.recipe.balance import balance
from msl.resolve import Resolved
from msl.schema.provenance import RunProvenance
from msl.schema.quantity import Quantity
from msl.schema.recipe import State
from msl.schema.refs import SubstanceRef, to_pymatgen_formula
from msl.schema.result import AssayResult, Fidelity, ResultValue, Status, ValueKind
from msl.validation import _check, score, CaseResult


def comp(ref: str, formula: str | None, amount: str, state: State | None = None) -> Resolved:
    return Resolved(ref=SubstanceRef.model_validate(ref), name=ref, formula=formula, state=state,
                    amount=Quantity.model_validate(amount))


# ── 질량수지 ──────────────────────────────────────────────────────────────


def test_balance_absolute_with_water_volume() -> None:
    b = balance([comp("cas:7647-01-0", "HCl", "0.1 mol", State.AQUEOUS), comp("cas:7732-18-5", "H2O", "1 L", State.LIQUID)])
    assert b.basis == "absolute"
    assert b.mass_fraction("cas:7647-01-0") == pytest.approx(3.646 / (3.646 + 1000), rel=1e-3)
    assert b.elements["Cl"] == pytest.approx(0.1)


def test_balance_concentration_uses_water_volume() -> None:
    b = balance([comp("cas:7647-14-5", "NaCl", "0.5 mol/L", State.AQUEOUS), comp("cas:7732-18-5", "H2O", "2 L", State.LIQUID)])
    assert b.elements["Na"] == pytest.approx(1.0)


def test_balance_fractions_relative() -> None:
    b = balance([comp("element:Cu", "Cu", "70 at%"), comp("element:Ni", "Ni", "30 at%")])
    assert b.basis == "relative"
    cu = next(r for r in b.rows if r["formula"] == "Cu")
    assert cu["mass_fraction"] == pytest.approx(0.7 * 63.546 / (0.7 * 63.546 + 0.3 * 58.693), rel=1e-4)


def test_balance_warns_on_unconvertible() -> None:
    b = balance([comp("cas:64-17-5", "C2H6O", "10 vol%"), comp("cas:7732-18-5", "H2O", "90 vol%")])
    assert b.basis == "none" and len(b.warnings) == 2


def test_hydrate_formula() -> None:
    assert to_pymatgen_formula("CuSO4·5H2O") == "CuSO4(H2O)5"
    assert to_pymatgen_formula("Na2CO3·H2O") == "Na2CO3(H2O)"


# ── A10 함량 기준 ─────────────────────────────────────────────────────────


def test_a10_threshold_parsing_and_judgement() -> None:
    (acc, human) = _classes({"typeList": [
        {"sbstnClsfTypeNm": "사고대비물질", "unqNo": "42", "contInfo": "염화 수소 및 이를 10% 이상 함유한 혼합물"},
        {"sbstnClsfTypeNm": "인체등유해성물질", "unqNo": "97-1-203", "contInfo": "인체급성유해성 : 10%, 생태유해성 : 2.5%"},
        {"sbstnClsfTypeNm": "기존화학물질", "unqNo": "KE-20189"}]})
    assert acc["threshold"] == 10 and human["threshold"] == 2.5  # 가장 낮은 % 가 기준
    assert not _applies(acc, 0.36) and _applies(acc, 36.0) and _applies(acc, None)  # 농도 모르면 보수적으로 해당
    assert "미만" in _judge(acc, 0.36) and "해당" in _judge(acc, 36.0)


# ── 검증 채점 ─────────────────────────────────────────────────────────────


def _res(**kw: object) -> AssayResult:
    now = dt.datetime(2026, 9, 27, tzinfo=dt.UTC)
    base = dict(assay="A2", recipe_id="rcp-t", status=Status.OK, fidelity=Fidelity.L0, engine="pymatgen",
                engine_version="x", conditions_basis="0 K", energy_reference="MP2020",
                values=[ResultValue(name="pH", kind=ValueKind.PROPERTY, value=8.2)],
                provenance=RunProvenance(input_hash="0" * 64, msl_version="t", engine_versions={}, started_at=now, finished_at=now))
    return AssayResult.model_validate(base | kw)


def test_check_value_between_and_products() -> None:
    assert _check(_res(), {"value_between": {"name": "pH", "range": [8.1, 8.3]}})[0]
    res = _res(data={"curves": [{"label": "1300 K", "kinks": [
        {"x": 0, "e": 0, "rxn": "ZnO -> ZnO"}, {"x": 0.3, "e": -55.8, "rxn": "0.5 ZnO + 0.5 Fe2O3 -> 0.5 Zn(FeO2)2"},
        {"x": 1, "e": 0, "rxn": "Fe2O3 -> Fe2O3"}]}]})
    assert _check(res, {"product_contains": "ZnFe2O4"})[0]  # 표기 순서가 달라도 조성으로 비교
    assert _check(res, {"curve_min_below": -1})[0] and not _check(res, {"curve_min_above": -1})[0]


def test_score_requires_full_s0_recall() -> None:
    ok = [CaseResult("a", "S0", True, hazard_expected=True)] + [CaseResult(f"b{i}", "A4", True) for i in range(9)]
    assert score(ok)["passed"]
    miss = [CaseResult("a", "S0", False, hazard_expected=True)] + ok[1:]
    assert not score(miss)["passed"]


# ── 리포트 ────────────────────────────────────────────────────────────────


def test_markdown_report_renders() -> None:
    from msl.report import to_markdown

    body = {"recipe": {"name": "시험", "id": "rcp-t"}, "elapsed": 0.1,
            "components": [{"ref": "formula:MgO", "name": "MgO", "formula": "MgO", "amount": "1 mol", "error": None, "via": "화학식"}],
            "balance": {"basis": "absolute", "rows": [{"formula": "MgO", "moles": 1.0, "mass_g": 40.3, "mass_fraction": 1.0}]},
            "results": [_res().model_dump(mode="json")], "assay_names": {"A2": "고상·계면 반응"}, "source_names": {}}
    md = to_markdown(body)
    assert "# 시험" in md and "## 질량수지" in md and "A2 고상·계면 반응" in md
