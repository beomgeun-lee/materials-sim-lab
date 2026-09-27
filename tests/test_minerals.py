"""IMA 광물 목록 파서 · COD 화학식 · 광물명 매칭 · mineral: 해석 — 네트워크 없이."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from msl.connectors import cod, ima_cnmnc
from msl.connectors.ima_cnmnc import FormulaError, build_index, chem, match_name, name_key, parse_formula

FIXTURE = Path(__file__).parent / "fixtures" / "ima_cnmnc_excerpt.tsv"


@pytest.fixture(scope="module")
def rows() -> dict[str, dict]:
    words = ima_cnmnc.read_tsv(FIXTURE.read_text(encoding="utf-8"))
    return {r["name"]: r for r in ima_cnmnc.parse_words(words)}


# ── PDF 표 복원 (pdftotext -tsv 발췌) ────────────────────────────────────


def test_excerpt_row_count(rows: dict[str, dict]) -> None:
    assert len(rows) == 14  # 머리말·표 머리는 행이 아니다


def test_single_line_row(rows: dict[str, dict]) -> None:
    r = rows["Abellaite"]
    assert r == {"name": "Abellaite", "formula_ima": "NaPb2(CO3)2(OH)", "status": "A",
                 "ima_number": "2014-111", "year": 2014, "country": "Spain"}


def test_wrapped_cells(rows: dict[str, dict]) -> None:
    # 나라 칸이 이름 줄 위아래로 감긴 행
    assert rows["Abelloemringerite"]["country"] == "Republic of the Congo"
    assert rows["Actinolite"]["country"] == "Germany / Austria"
    assert rows["Actinolite"]["year"] == 2012 and rows["Actinolite"]["ima_number"] is None  # '2012 s.p.'
    # 두 줄 화학식: 둘째 줄의 아래첨자가 다음 행 이름에 더 가까워도 제 행으로 간다
    assert rows["Carlfrancisite"]["formula_ima"] == \
        "Mn2+3(Mn2+,Mg,Fe3+,Al)42(As3+O3)2(As5+O4)4[(Si,As5+)O4]8(OH)42"
    assert rows["Carlfriesite"]["formula_ima"] == "CaTe4+2Te6+O8"


def test_superscripts_rejoined(rows: dict[str, dict]) -> None:
    # 따로 찍힌 위첨자(산화수)를 제자리에 넣는다 (-layout 텍스트에서는 윗줄로 떨어져 나간다)
    assert rows["Adanite"]["formula_ima"] == "Pb2(Te4+O3)(SO4)"
    assert rows["Abhurite"]["formula_ima"] == "Sn2+21O6(OH)14Cl16"
    # 이름 칸의 첨자도 이름의 일부
    assert rows["Julgoldite-(Fe2+)"]["formula_ima"] == "Ca2(Fe2+Fe3+2)(Si2O7)(SiO4)(OH)2(H2O)"
    assert "Julgoldite-(Fe3+)" in rows


def test_year_cell_split_around_status(rows: dict[str, dict]) -> None:
    r = rows["Cryptomelane"]  # 연도 칸이 '1982 s.p.' / '?' 두 줄이라 상태 기호 줄에는 연도가 없다
    assert (r["status"], r["year"], r["ima_number"]) == ("A", 1982, None)
    assert r["formula_ima"] == "K(Mn4+7Mn3+)O16"


# ── IMA 화학식 ────────────────────────────────────────────────────────────


@pytest.mark.parametrize("formula, expected", [
    ("Cu2+Mn3+6O8(SiO4)", {"Cu": 1, "Mn": 6, "O": 12, "Si": 1}),  # Mn3+6 = Mn(3+) × 6
    ("Sn2+21O6(OH)14Cl16", {"Sn": 21, "O": 20, "H": 14, "Cl": 16}),  # Sn2+21 = Sn(2+) × 21
    ("Cu+Cu2+2(OH)3Cl2", {"Cu": 3, "O": 3, "H": 3, "Cl": 2}),  # 숫자 없는 전하
    ("Fe3+", {"Fe": 1}),
    ("K(UO2)(AsO4)·3H2O", {"K": 1, "U": 1, "As": 1, "O": 9, "H": 6}),  # 수화물
    ("Fe2+Fe3+2Al3(PO4)4(OH)5(H2O)4∙2H2O", {"Fe": 3, "Al": 3, "P": 4, "O": 27, "H": 17}),  # ∙ 변형
    ("Mg2[B6O7(OH)6]2·9H2O", {"Mg": 2, "B": 12, "O": 35, "H": 30}),
    ("☐Al3Al6(Si5AlO18)(BO3)3(OH)3O", {"Al": 10, "Si": 5, "B": 3, "O": 31, "H": 3}),  # 공공 □
    ("Ba2Sr(Nb4.5Fe3+0.5)O15", {"Ba": 2, "Sr": 1, "Nb": 4.5, "Fe": 0.5, "O": 15}),
    ("K2PbCl4 (?)", {"K": 2, "Pb": 1, "Cl": 4}),  # 의문 표시
])
def test_parse_formula(formula: str, expected: dict[str, float]) -> None:
    assert parse_formula(formula) == pytest.approx(expected)


@pytest.mark.parametrize("formula", [
    "(K,Na)3(Fe3+,Ti,Al,Mg)5O2(AsO4)5",  # 치환 쉼표
    "Mo3O8·nH2O",  # 변수 n
    "Ag1-xSbx (x ~ 0.09-0.16)",  # 변수 x
    "☐Ca2(Mg4.5-2.5Fe2+0.5-2.5)Si8O22(OH)2",  # 범위
    "(Na0.5REE0.5)TiO3",  # 자리표시
    "Ca3Fe2[(AlO3(OH)]3",  # 괄호 짝 오류 (원문 오타)
])
def test_unparseable_formula(formula: str) -> None:
    with pytest.raises(FormulaError):
        parse_formula(formula)
    assert chem(formula)["formula_reduced"] is None


def test_chem_reduced_and_elements() -> None:
    assert chem("Ca(CO3)") == {"formula_reduced": "CaCO3", "elements": "C,Ca,O"}
    assert chem("SiO2")["formula_reduced"] == "SiO2"
    # 치환식은 축약식이 없어도 원소 집합은 안다. 자리표시(REE)가 있으면 원소 집합도 모른다
    assert chem("(Mg,Fe)2SiO4") == {"formula_reduced": None, "elements": "Fe,Mg,O,Si"}
    assert chem("(Na0.5REE0.5)TiO3")["elements"] is None


# ── COD ───────────────────────────────────────────────────────────────────


def test_cod_formula() -> None:
    assert cod.parse_formula("- Ca C O3 -") == {"Ca": 1, "C": 1, "O": 3}
    assert cod.parse_formula("- C3 D3 O7 Sr -") == {"C": 3, "H": 3, "O": 7, "Sr": 1}  # 중수소 → H
    assert cod.parse_formula("- F0.5 Fe6.1 Si8 -") == pytest.approx({"F": 0.5, "Fe": 6.1, "Si": 8})
    assert cod.reduced_formula(cod.parse_formula("- Ca4 Mg4 O24 Si8 -")) == "CaMg(SiO3)2"
    assert cod.parse_formula(None) == {}
    with pytest.raises(ValueError):
        cod.parse_formula("- Ca C O3 x -")


def test_cod_row() -> None:
    entry = {"id": "1010928", "attributes": {
        "_cod_mineral": "Calcite ", "_cod_calcformula": "- C Ca O3 -", "_cod_sg": "R -3 c :H", "_cod_sgnumber": None,
        "_cod_a": 4.99, "_cod_flags": "has coordinates,has Fobs", "_cod_duplicateof": None, "_cod_year": "1999"}}
    r = cod._row(entry)
    assert (r["cod_id"], r["mineral"], r["formula_reduced"], r["elements"]) == (1010928, "Calcite", "CaCO3", "C,Ca,O")
    assert r["sg_number"] == 167 and r["has_coordinates"] and r["year"] == 1999 and r["duplicate_of"] is None


# ── 광물명 정규화·매칭 ────────────────────────────────────────────────────


def test_name_key() -> None:
    assert name_key("Abenakiite-(Ce)") == name_key("abenakiite (Ce)") == name_key("ABENAKIITE-CE") == "abenakiitece"
    assert name_key("Alumoåkermanite") == "alumoakermanite"
    assert name_key("Bøggildite") == "boggildite"
    assert name_key("Zvěstovite-(Zn)") == "zvestovitezn"


IMA = ["Calcite", "Quartz", "Muscovite", "Baryte", "Wüstite", "Celsian", "Epidote", "Abenakiite-(Ce)",
       "Monazite-(Ce)", "Monazite-(La)", "Chevkinite-(Ce)", "Bastnäsite-(La)", "Diopside"]


@pytest.mark.parametrize("cod_name, expected", [
    ("Calcite", "Calcite"),
    ("calcite", "Calcite"),
    ("Calcite, magnesian", "Calcite"),  # 쉼표 뒤 수식어
    ("Quartz low", "Quartz"),
    ("Muscovite-2M1", "Muscovite"),  # 폴리타입
    ("Muscovite 3T", "Muscovite"),
    ("Diopside (deuterated)", "Diopside"),  # 괄호 주석
    ("Diopside-ferrian", "Diopside"),
    ("chromian epidote", "Epidote"),
    ("Barite", "Baryte"),  # 표기 동의어
    ("Wuestite", "Wüstite"),  # ue → ü
    ("Bastnaesite (La)", "Bastnäsite-(La)"),
    ("Celsian", "Celsian"),  # '-ian' 이지만 수식어가 아니라 광물명
    ("Abenakiite", "Abenakiite-(Ce)"),  # 뿌리가 한 종뿐
    ("Monazite", None),  # 뿌리가 여러 종 → 고르지 않는다
    ("Chevkinite-(Nd)", None),  # 원소 접미가 다르면 다른 종
    ("Olivine", None),  # 광물군 이름
    (None, None),
])
def test_match_name(cod_name: str | None, expected: str | None) -> None:
    assert match_name(cod_name, build_index(IMA)) == expected


def _cod_rows(**cols) -> pd.DataFrame:
    n = len(next(iter(cols.values())))
    base = {"cod_id": list(range(1, n + 1)), "elements": ["C,Ca,O"] * n, "formula_reduced": ["CaCO3"] * n,
            "sg_number": [167] * n, "has_coordinates": [True] * n, "duplicate_of": [None] * n,
            "cod_status": [None] * n, "cell_temp": [None] * n, "cell_pressure": [None] * n,
            "r_obs": [None] * n, "year": [2000] * n}
    return pd.DataFrame(base | cols)


def test_pick_best_order() -> None:
    pick = lambda df: ima_cnmnc.pick_best(df, "C,Ca,O", "CaCO3").cod_id
    assert pick(_cod_rows(has_coordinates=[False, True], r_obs=[0.01, 0.2])) == 2  # 좌표 있음이 먼저
    assert pick(_cod_rows(duplicate_of=[2, None], r_obs=[0.01, 0.2])) == 2  # 중복 아님
    assert pick(_cod_rows(cod_status=["retracted", None], r_obs=[0.01, 0.2])) == 2  # 철회되지 않음
    assert pick(_cod_rows(elements=["C,Ca,Mg,O", "C,Ca,O"], formula_reduced=["Ca0.9Mg0.1CO3", "CaCO3"],
                          r_obs=[0.01, 0.2])) == 2  # 원소 구성 일치
    assert pick(_cod_rows(formula_reduced=["Ca0.98C1.02O3", "CaCO3"], r_obs=[0.01, 0.2])) == 2  # 축약식 일치
    assert pick(_cod_rows(elements=["C,Ca,O", "C,Ca,H,O"], formula_reduced=["CaC0.9O3", "CaCO3"])) == 2  # H 무시
    assert pick(_cod_rows(cell_temp=[1073.0, 295.0], r_obs=[0.01, 0.2])) == 2  # 상온
    assert pick(_cod_rows(cell_pressure=[5.0e6, None], r_obs=[0.01, 0.2])) == 2  # 상압
    assert pick(_cod_rows(r_obs=[0.08, 0.03])) == 2  # R 낮음
    assert pick(_cod_rows(r_obs=[None, 0.05])) == 2  # R 없음은 뒤로
    assert pick(_cod_rows(year=[1990, 2010])) == 2  # 나머지가 같으면 최신


def test_pick_best_prefers_modal_space_group() -> None:
    # 'Spinel' 이름이 붙은 고압상(Pbnm)보다 대다수(Fd-3m) 구조를 고른다
    df = _cod_rows(sg_number=[227, 227, 62], year=[1980, 1990, 2007], elements=["Al,Mg,O"] * 3,
                   formula_reduced=["MgAl2O4"] * 3)
    assert ima_cnmnc.pick_best(df, "Al,Mg,O", "MgAl2O4").cod_id == 2


def test_pick_best_substitution_formula() -> None:
    # IMA 치환식 (Mg,Fe)2SiO4: COD 원소가 목록 안에 있으면 맞는 것으로 본다
    df = _cod_rows(elements=["Mg,Ni,O,Si", "Mg,O,Si"], formula_reduced=["x", "Mg2SiO4"], r_obs=[0.01, 0.2])
    assert ima_cnmnc.pick_best(df, "Fe,Mg,O,Si", None).cod_id == 2


# ── mineral: 해석 (DB 먼저, 없으면 PubChem) ───────────────────────────────


def _component(ref: str):
    from msl.schema.recipe import Component

    return Component.model_validate({"ref": ref, "amount": "1 mol", "state": "solid"})


def test_resolve_mineral_from_db(monkeypatch: pytest.MonkeyPatch) -> None:
    import msl.resolve as res

    db = {"calcite": {"name": "Calcite", "formula": "CaCO3", "version": "2026-09", "best_cod_id": 1010928,
                      "best_formula": "CaCO3"},
          "augite": {"name": "Augite", "formula": None, "version": "2026-09", "best_cod_id": 9000123,
                     "best_formula": "CaMg(SiO3)2"}}
    monkeypatch.setattr(res, "mineral_db", lambda: db)
    monkeypatch.setattr(res.pubchem, "lookup", lambda key: None)  # 네트워크 없음
    r = res.resolve(_component("mineral:Calcite"))
    assert (r.name, r.formula, r.via, r.error) == ("Calcite", "CaCO3", "IMA 2026-09 · COD 1010928", None)
    r = res.resolve(_component("mineral:augite"))  # IMA 치환식 → 대표 COD 구조의 화학식
    assert r.formula == "CaMg(SiO3)2" and "화학식은 COD 구조" in r.via


def test_resolve_mineral_falls_back_to_pubchem(monkeypatch: pytest.MonkeyPatch) -> None:
    import msl.resolve as res

    monkeypatch.setattr(res, "mineral_db", dict)
    monkeypatch.setattr(res.pubchem, "lookup", lambda key: {"CID": 516889, "MolecularFormula": "CCaO3",
                                                            "Title": "Calcium Carbonate"})
    monkeypatch.setattr(res.pubchem, "reactive_groups", lambda cid: [])
    monkeypatch.setattr(res.pubchem, "cas_number", lambda cid: "471-34-1")
    r = res.resolve(_component("mineral:calcite"))
    assert (r.formula, r.cid, r.cas) == ("CaCO3", 516889, "471-34-1")
    assert r.via.startswith("PubChem")
