"""웹 작업대 — 성분 추정, 폼 ↔ 레시피 ↔ YAML, 오류 문구, 저장, API 함수 (서버·httpx 없이 함수 직접 호출)."""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from fastapi import HTTPException
from pymatgen.core import Composition

from msl.registry.load import load_registry
from msl.resolve import korean_aliases, molecular_formula
from msl.schema.recipe import load_recipe
from msl.web import workbench as wb

REG = load_registry()
EXAMPLES = sorted((Path(__file__).resolve().parents[1] / "examples" / "recipes").glob("*.yaml"))


# ── 성분 추정 ─────────────────────────────────────────────────────────────


@pytest.mark.parametrize(("text", "ref"), [
    ("염산", "name:염산"), ("7647-01-0", "cas:7647-01-0"), ("Fe", "element:Fe"), ("Fe2O3", "formula:Fe2O3"),
    ("CuSO4·5H2O", "formula:CuSO4·5H2O"), ("cas:64-17-5", "cas:64-17-5"), ("CAS: 64-17-5", "cas:64-17-5"),
    ("KE-20189", "ke:KE-20189"), ("sodium chloride", "name:sodium chloride"),
])
def test_guess_ref(text: str, ref: str) -> None:
    assert wb.guess_ref(text) == ref


def test_guess_ref_mineral() -> None:
    from msl.resolve import mineral_db

    if not mineral_db():
        pytest.skip("IMA 광물 목록 미적재")
    assert wb.guess_ref("hematite") == "mineral:hematite"


# ── 폼 ↔ 레시피 ↔ YAML ────────────────────────────────────────────────────


@pytest.mark.parametrize("path", EXAMPLES, ids=[p.stem for p in EXAMPLES])
def test_form_roundtrip(path: Path) -> None:
    """예제 레시피 → 폼 → 검증 → YAML → 다시 읽기 가 같은 레시피."""
    r = load_recipe(path, known_assays=set(REG.assays))
    again, errors = wb.validate_recipe(wb.to_form(r), REG)
    assert not errors and again.model_dump() == r.model_dump()
    data, errors = wb.parse_yaml(wb.dump_yaml(r))
    assert not errors and wb.validate_recipe(data, REG)[0].model_dump() == r.model_dump()


def test_blank_fields_are_dropped() -> None:
    form = {"id": "rcp-t", "name": "t", "components": [{"ref": "element:Fe", "amount": "1 g", "state": "", "label": ""}],
            "conditions": {"T": "", "pH": None}, "notes": "", "assays": "auto"}
    r, errors = wb.validate_recipe(form, REG)
    assert not errors and r.components[0].state is None and r.conditions.T is None and r.notes is None
    _, errors = wb.validate_recipe({**form, "components": [{"ref": "element:Fe", "amount": ""}]}, REG)
    assert errors == ["(레시피): mixture 모드는 모든 성분에 양이 필요함: ['element:Fe']"]


def test_errors_in_korean() -> None:
    _, e1 = wb.validate_recipe({"id": "bad id", "name": "x", "components": []}, REG)
    assert e1 == ["id: 형식이 맞지 않음 — rcp- 로 시작하고 영소문자·숫자·하이픈만 (예: rcp-my-test)",
                  "성분: 꼭 적어야 하는 항목이 없음 (성분은 하나 이상)"]
    _, e2 = wb.validate_recipe({"id": "rcp-t", "name": "x", "components": [{"ref": "element:Fe", "amount": "1 g"}], "conditions": {"pH": 20}}, REG)
    assert e2 == ["조건 · pH: 값이 허용 범위를 벗어남 (최대 16)"]
    _, e3 = wb.validate_recipe({"id": "rcp-t", "name": "x", "components": [{"ref": "element:Fe", "amount": "1 g"}], "assays": ["A99"]}, REG)
    assert "A99" in e3[0]


# ── 저장 ──────────────────────────────────────────────────────────────────


def test_save_text(tmp_path: Path) -> None:
    path, err = wb.save_text(tmp_path, "rcp-my-test", "a", overwrite=False)
    assert err is None and path.read_text() == "a"
    assert wb.save_text(tmp_path, "rcp-my-test", "b", overwrite=False)[1].startswith("이미 있는 파일")
    assert wb.save_text(tmp_path, "rcp-my-test", "b", overwrite=True)[0].read_text() == "b"
    assert wb.save_text(tmp_path, "../evil", "x", overwrite=True)[0] is None
    listed = [(tmp_path / "rcp-my-test.yaml", "user")]
    assert wb.find_file(listed, "../../etc/passwd") is None  # 목록에 있는 파일만 돌려준다 — 경로로 빠져나갈 수 없음
    assert wb.find_file(listed, "../rcp-my-test") == listed[0][0]
    assert wb.find_file(listed, "rcp my test") is None


# ── API 함수 ──────────────────────────────────────────────────────────────


@pytest.fixture
def api(tmp_path, monkeypatch):
    from msl.web import app as web

    monkeypatch.setattr(wb, "USER_RECIPES", tmp_path / "recipes")
    monkeypatch.setattr(wb, "USER_SPACES", tmp_path / "spaces")
    return web


FORM = {"id": "rcp-my-mgo-al2o3", "name": "스피넬 (웹)", "components": [
    {"ref": "formula:MgO", "amount": "1 mol", "state": "solid"}, {"ref": "formula:Al2O3", "amount": "1 mol", "state": "solid"}],
    "conditions": {"T": "1400 °C"}, "assays": "auto"}


def test_api_check_and_save(api) -> None:
    res = api.recipe_check(api.RecipeRequest(recipe=FORM))
    assert res["ok"] and [r["code"] for r in res["route"]][:3] == ["S0", "A1", "A2"]
    assert res["balance"]["basis"] == "absolute" and "components:" in res["yaml"]
    bad = api.recipe_check(api.RecipeRequest(recipe={**FORM, "id": "x"}))
    assert not bad["ok"] and bad["errors"][0].startswith("id:")

    saved = api.recipe_save(api.RecipeRequest(recipe=FORM))
    assert (wb.USER_RECIPES / saved["file"]).exists()
    with pytest.raises(HTTPException) as dup:
        api.recipe_save(api.RecipeRequest(recipe=FORM))
    assert dup.value.status_code == 409
    api.recipe_save(api.RecipeRequest(recipe={**FORM, "name": "바뀐 이름"}, overwrite=True))
    assert "바뀐 이름" in (wb.USER_RECIPES / saved["file"]).read_text(encoding="utf-8")
    with pytest.raises(HTTPException) as ex:
        api.recipe_save(api.RecipeRequest(recipe={**FORM, "id": "rcp-iron-water"}))
    assert ex.value.status_code == 409 and "예제" in ex.value.detail
    assert any(r["source"] == "user" and r["id"] == FORM["id"] for r in api.recipes())


def test_api_space_check_and_save(api) -> None:
    text = """id: spc-my-t
name: t
base: {assays: [A2]}
generator: {kind: ratio, components: ["formula:MgO", "formula:Al2O3"], step: 0.25}
collect: [{assay: A2, value: 최저 반응에너지}]
"""
    res = api.space_check(api.YamlRequest(yaml=text))
    assert res["ok"] and res["n"] == 3
    assert not api.space_check(api.YamlRequest(yaml=text.replace("step: 0.25", "step: 0.3")))["ok"]
    assert api.space_save(api.YamlRequest(yaml=text))["file"] == "spc-my-t.yaml"
    with pytest.raises(HTTPException):
        api.space_save(api.YamlRequest(yaml=text.replace("spc-my-t", "spc-cu-ni-composition")))


# ── 국문 관용명 · 분자식 ──────────────────────────────────────────────────


def test_korean_alias_table() -> None:
    rows = korean_aliases()
    assert rows["염산"]["cas"] == "7647-01-0" and rows["가성소다"]["cas"] == "1310-73-2"
    for name, row in rows.items():
        assert re.fullmatch(r"\d{2,7}-\d{2}-\d", row["cas"]), name
        Composition(row["formula"])
        assert row.get("state") in (None, "solid", "liquid", "gas", "aqueous"), name


@pytest.mark.parametrize(("formula", "shown"), [
    ("C2H4O2", "C2H4O2"), ("C6H12O6", "C6H12O6"), ("C2H6O", "C2H6O"), ("CH4", "CH4"), ("H2O2", "H2O2"),
    ("NaHCO3", "NaHCO3"), ("Fe2O3", "Fe2O3"), ("ClH", "HCl"), ("C12H22O11", "C12H22O11"),
])
def test_molecular_formula(formula: str, shown: str) -> None:
    """약분하면 분자가 바뀌는 화학식(초산·포도당)은 그대로, 유기물은 힐 표기."""
    assert molecular_formula(Composition(formula)) == shown
