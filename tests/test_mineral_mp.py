"""IMA 광물 ↔ MP phases_calc 매칭 (mineral_mp) · 공간군 정규화 · mineral: 해석의 MP ID — 네트워크 없이."""

from __future__ import annotations

import duckdb
import pandas as pd
import pytest

from msl.connectors import cod, ima_cnmnc
from msl.connectors.ima_cnmnc import (LEVEL_FORMULA, LEVEL_SG, LEVEL_SUPER, match_mp_frame, mp_formula, sg_key,
                                      sg_tier)

# ── 공간군 ────────────────────────────────────────────────────────────────


def test_setting_variants_share_number() -> None:
    # H-M 기호의 설정 차이는 번호로 비교하면 사라진다 (aragonite 는 COD 에 Pmcn, MP 에 Pnma 로 올라 있다)
    assert cod._sg_from_symbol("P m c n") == cod._sg_from_symbol("P n m a") == cod._sg_from_symbol("P b n m") == 62
    assert cod._sg_from_symbol("R -3 c :H") == cod._sg_from_symbol("R -3 c :R") == 167


@pytest.mark.parametrize("symbol, number", [
    ("P 21 21 21 (a+1/4,b-1/4,c)", 19),  # 원점 이동
    ("P 4 21 2 (a-1/4,b+1/4,c)", 90),  # 원점 이동을 떼도 기호 넷은 그대로
    ("R -3 c RS", 167),
    ("P m c n S1", 62),
    ("P 4/n n c Z1", 126),
    ("C c c b :1", 68),
    ("I 41/a m d 1", 141),  # 기호 넷 뒤의 원점 선택 숫자
    ("?", None),
])
def test_cod_symbol_tags_stripped(symbol: str, number: int | None) -> None:
    assert cod._sg_from_symbol(symbol) == number


def test_sg_key_enantiomorphs() -> None:
    assert sg_key(152) == sg_key(154) == 152  # 석영 P3₁21 / P3₂21
    assert sg_key(92.0) == sg_key(96) == 92  # 크리스토발석, pandas 실수 열도 받는다
    assert sg_key(62) == 62
    assert sg_key(None) is None and sg_key(float("nan")) is None and sg_key(pd.NA) is None


def test_sg_tier() -> None:
    assert sg_tier(62, 62) == 0
    assert sg_tier(154, 152) == 0  # 거울상 쌍
    assert sg_tier(194, 186) == 1  # MP P6₃/mmc 은 COD P6₃mc 의 바로 위 군 (흑연)
    assert sg_tier(15, 9) == 1  # C2/c ⊃ Cc
    assert sg_tier(186, 194) == 2  # MP 쪽이 더 낮은 대칭이면 근거가 아니다
    assert sg_tier(62, 167) == 2
    assert sg_tier(62, None) == 2 and sg_tier(None, 62) == 2


# ── 비교할 화학식 ──────────────────────────────────────────────────────────


def test_mp_formula_prefers_ima() -> None:
    assert mp_formula("CaCO3", "C,Ca,O", "CaCO3") == ("CaCO3", "ima")
    assert mp_formula("Mg(HO)2", "H,Mg,O", "MgO2") == ("Mg(HO)2", "ima")  # IMA 식이 있으면 COD 식은 안 본다


def test_mp_formula_cod_fallback_guarded() -> None:
    # 치환식 (Mg,Fe)2SiO4: IMA 축약식이 없으면 대표 COD 식
    assert mp_formula(None, "Fe,Mg,O,Si", "Mg2SiO4") == ("Mg2SiO4", "cod")
    assert mp_formula(None, "H,Mg,O", "MgO2") == (None, None)  # X선 구조라 H 가 빠진 식 → 과산화물에 붙지 않게
    assert mp_formula(None, "Fe,Mg,O,Si", "Ni2SiO4") == (None, None)  # IMA 에 없는 원소
    assert mp_formula(None, None, "CeCO3F") == ("CeCO3F", "cod")  # REE 자리표시라 원소 목록이 없으면 COD 식
    assert mp_formula(None, None, None) == (None, None)


# ── 매칭 규칙 ─────────────────────────────────────────────────────────────


def _mp(*rows) -> pd.DataFrame:
    cols = ["material_id", "formula_pretty", "spacegroup_symbol", "spacegroup_number", "energy_above_hull_eV",
            "deprecated"]
    return pd.DataFrame(rows, columns=cols)


def _minerals(*rows) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=["name", "formula_reduced", "elements", "best_formula", "cod_sg_number"])


MP = _mp(
    ("mp-cal", "CaCO3", "R-3c", 167, 0.0, False),
    ("mp-ara", "CaCO3", "Pnma", 62, 0.001, False),
    ("mp-ara2", "CaCO3", "Pnma", 62, 0.013, False),
    ("mp-tri", "CaCO3", "P-1", 2, 0.003, False),
    ("mp-old", "CaCO3", "R-3c", 167, -0.01, True),  # deprecated: 더 낮아도 후보가 아니다
    ("mp-qz", "SiO2", "P3_121", 152, 0.0, False),
    ("mp-cri", "SiO2", "P4_32_12", 96, 0.005, False),
    ("mp-sti", "SiO2", "P4_2/mnm", 136, 0.2, False),
    ("mp-c12", "C", "C2/m", 12, 0.0, False),  # 가장 낮지만 흑연이 아니다
    ("mp-gra", "C", "P6_3/mmc", 194, 0.003, False),
    ("mp-dia", "C", "Fd-3m", 227, 0.11, False),
    ("mp-fo", "Mg2SiO4", "Pnma", 62, 0.0, False),
)

MINERALS = _minerals(
    ("Calcite", "CaCO3", "C,Ca,O", "CaCO3", 167),
    ("Aragonite", "CaCO3", "C,Ca,O", "CaCO3", 62),  # COD 는 Pmcn → 번호 62
    ("Vaterite", "CaCO3", "C,Ca,O", None, None),  # COD 구조 없음 → 다형이 있어 대표를 두지 않는다
    ("Quartz", "SiO2", "O,Si", "SiO2", 154),  # COD P3₂21, MP P3₁21
    ("Cristobalite", "SiO2", "O,Si", "SiO2", 92),
    ("Stishovite", "SiO2", "O,Si", "SiO2", 136),
    ("Graphite", "C", "C", "C", 186),  # COD 대표가 P6₃mc
    ("Diamond", "C", "C", "C", 227),
    ("Forsterite", None, "Fe,Mg,O,Si", "Mg2SiO4", 62),  # 치환식 → COD 식
    ("Halite", "NaCl", "Cl,Na", "NaCl", 225),  # MP 에 없음
)


@pytest.fixture(scope="module")
def matched() -> pd.DataFrame:
    return match_mp_frame(MINERALS, MP).set_index("name")


@pytest.mark.parametrize("name, best, level, sg", [
    ("Calcite", "mp-cal", LEVEL_SG, "R-3c"),
    ("Aragonite", "mp-ara", LEVEL_SG, "Pnma"),
    ("Quartz", "mp-qz", LEVEL_SG, "P3_121"),
    ("Cristobalite", "mp-cri", LEVEL_SG, "P4_32_12"),
    ("Stishovite", "mp-sti", LEVEL_SG, "P4_2/mnm"),  # E_hull 이 높아도 공간군이 맞는 쪽
    ("Graphite", "mp-gra", LEVEL_SUPER, "P6_3/mmc"),  # 가장 낮은 C2/m 이 아니라 상위군
    ("Diamond", "mp-dia", LEVEL_SG, "Fd-3m"),
    ("Forsterite", "mp-fo", LEVEL_SG, "Pnma"),
])
def test_polymorphs_get_distinct_ids(matched: pd.DataFrame, name: str, best: str, level: str, sg: str) -> None:
    r = matched.loc[name]
    assert (r.best_mp_id, r.match_level, r.mp_sg) == (best, level, sg)
    assert not r.polymorph_ambiguous and pd.isna(r.mp_shared_with)


def test_formula_only_polymorph_left_open(matched: pd.DataFrame) -> None:
    r = matched.loc["Vaterite"]
    assert r.match_level == LEVEL_FORMULA and r.polymorph_ambiguous
    assert pd.isna(r.best_mp_id) and pd.isna(r.mp_e_above_hull_summary_eV)  # 빈 값은 pandas 결측
    assert r.mp_ids == ["mp-cal", "mp-ara", "mp-tri", "mp-ara2"]  # 후보는 남긴다 (E_hull 순)


def test_unmatched_and_basis(matched: pd.DataFrame) -> None:
    r = matched.loc["Halite"]
    assert r.n_mp == 0 and r.mp_ids == [] and pd.isna(r.best_mp_id) and pd.isna(r.match_level)
    assert pd.isna(r.polymorph_ambiguous)
    assert matched.loc["Forsterite", "formula_basis"] == "cod" and matched.loc["Calcite", "formula_basis"] == "ima"


def test_deprecated_excluded(matched: pd.DataFrame) -> None:
    assert "mp-old" not in matched.loc["Calcite", "mp_ids"]
    assert matched.loc["Calcite", "n_mp"] == 4


def test_candidate_order(matched: pd.DataFrame) -> None:
    # 같은 공간군(E_hull 순) → 상위군 → 나머지(E_hull 순)
    assert matched.loc["Aragonite", "mp_ids"] == ["mp-ara", "mp-ara2", "mp-cal", "mp-tri"]
    assert matched.loc["Graphite", "mp_ids"] == ["mp-gra", "mp-c12", "mp-dia"]
    # E_hull 이 없으면 뒤로, 같으면 material_id 순
    mp = _mp(("mp-b", "ZnS", "F-43m", 216, 0.0, False), ("mp-a", "ZnS", "F-43m", 216, 0.0, False),
             ("mp-0", "ZnS", "F-43m", 216, None, False))
    df = match_mp_frame(_minerals(("Sphalerite", "ZnS", "S,Zn", "ZnS", 216)), mp)
    assert df.loc[0, "mp_ids"] == ["mp-a", "mp-b", "mp-0"] and df.loc[0, "best_mp_id"] == "mp-a"


def test_supergroup_skips_id_claimed_by_sibling() -> None:
    # R3c(161) 구조의 가상 광물: 상위군 R-3c 후보 중 calcite 가 공간군까지 맞춰 잡은 항목은 건너뛴다
    mp = pd.concat([MP, _mp(("mp-cal2", "CaCO3", "R-3c", 167, 0.05, False))])
    minerals = pd.concat([MINERALS, _minerals(("Foo", "CaCO3", "C,Ca,O", "CaCO3", 161))])
    df = match_mp_frame(minerals, mp).set_index("name")
    assert (df.loc["Foo", "best_mp_id"], df.loc["Foo", "match_level"]) == ("mp-cal2", LEVEL_SUPER)
    assert df.loc["Calcite", "best_mp_id"] == "mp-cal"
    # 상위군 후보가 모두 잡혀 있으면 화학식 단계로 내려가고, 다형이 있으니 대표를 두지 않는다
    df = match_mp_frame(minerals, MP).set_index("name")
    assert pd.isna(df.loc["Foo", "best_mp_id"]) and df.loc["Foo", "polymorph_ambiguous"]


def test_formula_only_without_polymorph_keeps_ground_state() -> None:
    df = match_mp_frame(_minerals(("Forsterite", "Mg2SiO4", "Mg,O,Si", None, None)), MP)
    assert (df.loc[0, "best_mp_id"], df.loc[0, "match_level"], df.loc[0, "polymorph_ambiguous"]) == \
        ("mp-fo", LEVEL_FORMULA, False)


def test_same_formula_and_spacegroup_flagged_shared() -> None:
    # 규회석·브레이석처럼 축약식·공간군이 모두 같으면 공간군으로 못 가른다 → 서로를 표시
    mp = _mp(("mp-w", "CaSiO3", "P-1", 2, 0.0, False), ("mp-w2", "CaSiO3", "P-1", 2, 0.01, False))
    df = match_mp_frame(_minerals(("Wollastonite", "CaSiO3", "Ca,O,Si", "CaSiO3", 2),
                                  ("Breyite", "CaSiO3", "Ca,O,Si", "CaSiO3", 2)), mp).set_index("name")
    assert df.loc["Wollastonite", "mp_shared_with"] == "Breyite"
    assert df.loc["Breyite", "mp_shared_with"] == "Wollastonite"


# ── DB 적재 (in-memory duckdb) ─────────────────────────────────────────────


def _fake_db() -> duckdb.DuckDBPyConnection:
    con = duckdb.connect()
    minerals = MINERALS[["name", "formula_reduced", "elements"]].assign(source_version="2026-09")
    structures = MINERALS[["name", "best_formula"]].assign(best_cod_id=range(1, len(MINERALS) + 1))
    exp = pd.DataFrame({"cod_id": structures["best_cod_id"], "sg_number": MINERALS["cod_sg_number"].astype("Int64"),
                        "source_version": "2026-09-27"})
    phases = MP.assign(source_version="v2026.04.13", correction_scheme="MP summary 혼합")
    for name, df in {"minerals": minerals, "mineral_structures": structures, "structures_exp": exp,
                     "phases_calc": phases}.items():
        con.register(f"_{name}", df)
        con.execute(f"CREATE VIEW {name} AS SELECT * FROM _{name}")
    return con


def test_match_mp_writes_with_provenance(monkeypatch: pytest.MonkeyPatch) -> None:
    written: dict[str, pd.DataFrame] = {}
    monkeypatch.setattr(ima_cnmnc, "connect", lambda **kw: _fake_db())
    monkeypatch.setattr(ima_cnmnc, "write_table", lambda t, df: written.setdefault(t, df) is not None and {"open": 1})
    ima_cnmnc.match_mp()
    df = written["mineral_mp"].set_index("name")
    assert df.loc["Aragonite", "best_mp_id"] == "mp-ara" and df.loc["Aragonite", "cod_sg_number"] == 62
    assert df.loc["Graphite", "best_mp_id"] == "mp-gra"
    assert set(df["source"]) == {"ima-cnmnc"} and set(df["license"]) == {"CC-BY-SA-3.0"}
    assert set(df["method"]) == {"compiled"}
    assert set(df["source_version"]) == {"2026-09+cod-2026-09-27+mp-v2026.04.13"}
    assert set(df["correction_scheme"]) == {"mp_e_above_hull_summary_eV: MP summary 혼합"}


@pytest.mark.parametrize("have, calls", [
    ({"minerals", "structures_exp", "phases_calc"}, ["match", "match_mp"]),  # COD 매칭 먼저 (대표 공간군을 쓴다)
    ({"minerals", "phases_calc"}, ["match_mp"]),
    ({"minerals", "structures_exp"}, ["match"]),
    ({"structures_exp", "phases_calc"}, []),  # 광물 목록이 없으면 아무것도 안 한다
])
def test_refresh_matches(monkeypatch: pytest.MonkeyPatch, have: set[str], calls: list[str]) -> None:
    seen: list[str] = []
    monkeypatch.setattr(ima_cnmnc, "tables", lambda: sorted(have))
    monkeypatch.setattr(ima_cnmnc, "match", lambda: seen.append("match") or {"open": 1})
    monkeypatch.setattr(ima_cnmnc, "match_mp", lambda: seen.append("match_mp") or {"open": 1})
    out = ima_cnmnc.refresh_matches()
    assert seen == calls
    assert list(out) == [{"match": "mineral_structures", "match_mp": "mineral_mp"}[c] for c in calls]


# ── COD commonname 확장 ───────────────────────────────────────────────────


def test_cod_row_keeps_commonname() -> None:
    r = cod._row({"id": "2300702", "attributes": {"_cod_mineral": None, "_cod_commonname": " Diamond ",
                                                  "_cod_calcformula": "- C -", "_cod_flags": "has coordinates"}})
    assert (r["mineral"], r["commonname"], r["formula_reduced"]) == (None, "Diamond", "C")


@pytest.mark.parametrize("cod_els, ima_els, exact, ok", [
    ("C,Ca,O", "C,Ca,O", True, True),
    ("C,Ca,H,O", "C,Ca,O", True, True),  # H 는 비교하지 않는다
    ("Ge,O", "O,Si", True, False),  # GeO2 'quartz' 같은 합성 유사체
    ("As,Fe,O", "Fe", True, False),  # 'Iron (II, III) arsenate' 가 괄호를 떼면 'Iron' 이 되는 경우
    ("Mg,O,Si", "Fe,Mg,O,Si", False, True),  # 치환식: 목록 안이면 된다
    ("Mg,O,Si", "Fe,Mg,O,Si", True, False),  # 축약식이 있으면 같아야 한다
    ("C,Ca,O", None, False, False),  # IMA 원소를 모르면 받지 않는다
    (None, "C,Ca,O", True, False),
])
def test_common_fits(cod_els, ima_els, exact: bool, ok: bool) -> None:
    assert ima_cnmnc.common_fits(cod_els, ima_els, exact) is ok


def _cod_db(rows: list[dict]) -> duckdb.DuckDBPyConnection:
    base = {"mineral": None, "commonname": None, "sg": "R -3 c :H", "sg_number": 167, "has_coordinates": True,
            "duplicate_of": None, "cod_status": None, "cell_temp": None, "cell_pressure": None, "r_obs": None,
            "year": 2000, "source_version": "2026-09-27"}
    con = duckdb.connect()
    minerals = pd.DataFrame({"name": ["Calcite", "Quartz", "Stepanovite"],
                             "formula_reduced": ["CaCO3", "SiO2", "NaMgFeH18(C2O7)3"],
                             "elements": ["C,Ca,O", "O,Si", "C,Fe,H,Mg,Na,O"]})
    exp = pd.DataFrame([base | r for r in rows]).astype(
        {"duplicate_of": "Int64", "r_obs": "Float64", "cod_status": "str"})
    for name, df in {"minerals": minerals, "structures_exp": exp}.items():
        con.register(f"_{name}", df)
        con.execute(f"CREATE VIEW {name} AS SELECT * FROM _{name}")
    return con


def test_match_uses_commonname_but_prefers_mineral_entries(monkeypatch: pytest.MonkeyPatch) -> None:
    calcite = {"elements": "C,Ca,O", "formula_reduced": "CaCO3"}
    con = _cod_db([
        {"cod_id": 9000001, "mineral": "Calcite", **calcite},  # AMCSD 광물 항목: R 값 없음
        {"cod_id": 2000002, "commonname": "calcite", "r_obs": 0.02, **calcite},  # R 값이 있어도 대표가 아니다
        {"cod_id": 2000003, "commonname": "Quartz", "elements": "Ge,O", "formula_reduced": "GeO2"},  # 조성이 다름
        {"cod_id": 2000004, "commonname": "Stepanovite", "elements": "C,Fe,H,Mg,Na,O",
         "formula_reduced": "NaMgFeH18(C2O7)3", "sg": "P 3 c 1", "sg_number": 158},  # 광물명 항목이 없는 종
        {"cod_id": 2000005, "commonname": "Sodium chloride", "elements": "Cl,Na", "formula_reduced": "NaCl"},
    ])
    written: dict[str, pd.DataFrame] = {}
    monkeypatch.setattr(ima_cnmnc, "connect", lambda **kw: con)
    monkeypatch.setattr(ima_cnmnc, "write_table", lambda t, df: written.setdefault(t, df) is not None and {"open": 1})
    ima_cnmnc.match()
    df = written["mineral_structures"].set_index("name")
    assert df.loc["Calcite", "cod_ids"] == [2000002, 9000001] and df.loc["Calcite", "n_cod_commonname"] == 1
    assert df.loc["Calcite", "best_cod_id"] == 9000001
    assert df.loc["Quartz", "n_cod"] == 0
    assert (df.loc["Stepanovite", "best_cod_id"], df.loc["Stepanovite", "best_sg"]) == (2000004, "P 3 c 1")


# ── mineral: 해석 → MP ID ─────────────────────────────────────────────────


def _component(ref: str):
    from msl.schema.recipe import Component

    return Component.model_validate({"ref": ref, "amount": "1 mol", "state": "solid"})


def test_resolve_mineral_carries_mp_id(monkeypatch: pytest.MonkeyPatch) -> None:
    import msl.resolve as res

    db = {"aragonite": {"name": "Aragonite", "formula": "CaCO3", "version": "2026-09", "best_cod_id": 9014355,
                        "best_formula": "CaCO3", "best_mp_id": "mp-aaaaagvy", "mp_match": LEVEL_SG},
          "graphite": {"name": "Graphite", "formula": "C", "version": "2026-09", "best_cod_id": 9012230,
                       "best_formula": "C", "best_mp_id": "mp-aaaaaabw", "mp_match": LEVEL_SUPER},
          "vaterite": {"name": "Vaterite", "formula": "CaCO3", "version": "2026-09", "best_cod_id": None,
                       "best_formula": None, "best_mp_id": None, "mp_match": LEVEL_FORMULA}}
    monkeypatch.setattr(res, "mineral_db", lambda: db)
    monkeypatch.setattr(res.pubchem, "lookup", lambda key: None)  # 네트워크 없음
    r = res.resolve(_component("mineral:aragonite"))
    assert (r.formula, r.cod_id, r.mp_id) == ("CaCO3", 9014355, "mp-aaaaagvy")
    assert r.via == "IMA 2026-09 · COD 9014355 · MP mp-aaaaagvy"
    r = res.resolve(_component("mineral:Graphite"))
    assert r.mp_id == "mp-aaaaaabw" and r.via.endswith("(공간군은 상위군 일치)")
    r = res.resolve(_component("mineral:vaterite"))  # 다형을 못 가렸으면 MP ID 없음
    assert (r.formula, r.mp_id, r.cod_id, r.via) == ("CaCO3", None, None, "IMA 2026-09")


def test_formula_ref_has_no_mp_id(monkeypatch: pytest.MonkeyPatch) -> None:
    import msl.resolve as res

    monkeypatch.setattr(res.pubchem, "lookup", lambda key: None)
    r = res.resolve(_component("formula:CaCO3"))  # 화학식에는 다형 정보가 없다
    assert r.formula == "CaCO3" and r.mp_id is None and r.cod_id is None
