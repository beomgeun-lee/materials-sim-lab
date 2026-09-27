"""A9 공급·경제성 — 집중도·판정 계산과 조회 (가짜 DataFrame·메모리 DuckDB, 네트워크·실데이터 없음)."""

from __future__ import annotations

import duckdb
import pandas as pd
import pytest

from msl.assays import supply
from msl.assays.base import Context
from msl.connectors import usgs_mcs, world_bank_pink_sheet
from msl.registry.load import load_registry
from msl.resolve import Resolved
from msl.schema.recipe import Recipe
from msl.schema.refs import SubstanceRef
from msl.schema.result import Status

W, O = supply.WORLD, supply.OTHERS


def prod(rows: list[tuple[str, int, float | None]], unit: str = "metric tons") -> pd.DataFrame:
    return pd.DataFrame(rows, columns=["country", "year", "value"]).assign(unit=unit)


# ── 순수 계산 ──────────────────────────────────────────────────────────────


def test_hhi_monopoly_and_even_split() -> None:
    assert supply.hhi([1.0]) == pytest.approx(10_000)
    assert supply.hhi([0.25] * 4) == pytest.approx(2_500)
    assert supply.hhi([0.5, 0.3, 0.2]) == pytest.approx(3_800)


def test_concentration_uses_reported_world_total_and_latest_year() -> None:
    df = prod([("A", 2024, 10), ("B", 2024, 90), (W, 2024, 100),
               ("A", 2025, 60), ("B", 2025, 20), (O, 2025, 20), (W, 2025, 100)])
    c = supply.concentration(df)
    assert c is not None and c.year == 2025 and c.world_reported
    assert (c.top_country, c.top_share) == ("A", pytest.approx(0.6))
    assert c.hhi == pytest.approx(60**2 + 20**2)  # 'Other countries' 는 HHI 에서 빠진다 (하한)
    assert c.n_countries == 2


def test_concentration_falls_back_to_sum_without_world_total() -> None:
    c = supply.concentration(prod([("A", 2025, 30), ("B", 2025, 10), (O, 2025, 10)]))
    assert c is not None and not c.world_reported
    assert c.world == 50 and c.top_share == pytest.approx(0.6)


def test_concentration_skips_withheld_and_empty() -> None:
    assert supply.concentration(prod([("A", 2025, None), (W, 2025, 100)])) is None
    c = supply.concentration(prod([("US", 2025, None), ("A", 2025, 80), (W, 2025, 100)]))
    assert c is not None and c.top_country == "A" and c.top_share == pytest.approx(0.8)


@pytest.mark.parametrize(("critical", "top", "hhi", "status"), [
    (False, 0.3, 1_500, Status.OK),
    (True, 0.3, 1_500, Status.WARNING),
    (False, 0.51, 1_500, Status.WARNING),
    (False, 0.50, 2_500, Status.OK),  # 경계값은 초과가 아니므로 보통
    (False, 0.4, 2_501, Status.WARNING),
    (False, None, None, Status.OK),
])
def test_judge_thresholds(critical: bool, top: float | None, hhi: float | None, status: Status) -> None:
    assert supply.judge(critical, top, hhi)[0] is status


def test_judge_lists_reasons() -> None:
    _, why = supply.judge(True, 0.67, 4_559)
    assert why == ["핵심광물", "1위 점유율 67%", "HHI 4,559"]


def test_pick_nir_latest_year_prefers_named_row_then_total() -> None:
    df = pd.DataFrame({
        "measure_detail": ["NIR: Ferrosilicon", "NIR: Silicon metal", "NIR: Total", "NIR: Total"],
        "year": [2025, 2025, 2025, 2024], "value": [50.0, 50.0, 50.0, 40.0],
        "value_flag": ["<", ">", ">", None], "value_raw": ["<50", ">50", ">50", "40"],
    })
    assert supply.pick_nir(df, "silicon metal")["text"] == ">50%"
    assert supply.pick_nir(df)["year"] == 2025 and "Total" in df.loc[2, "measure_detail"]
    exporter = pd.DataFrame({"measure_detail": ["NIR"], "year": [2025], "value": [None], "value_flag": ["E"],
                             "value_raw": ["E"]})
    assert supply.pick_nir(exporter)["text"] == "순수출"


def test_price_change_against_same_month_last_year() -> None:
    df = pd.DataFrame({"month": ["2025-08", "2026-07", "2026-08"], "price": [100.0, 140.0, 150.0], "unit": "$/mt"})
    p = supply.price_change(df)
    assert p["month"] == "2026-08" and p["price"] == 150 and p["change"] == pytest.approx(0.5)
    assert supply.price_change(df.iloc[1:])["change"] is None


# ── 조회 (메모리 DuckDB) ───────────────────────────────────────────────────


def stats_rows() -> pd.DataFrame:
    base = {"section": "World Mine Production and Reserves", "measure": "production", "unit": "metric tons",
            "value_flag": None, "value_raw": None, "is_critical_2025": True}
    rows = [
        {**base, "commodity": "Nickel", "country": "Indonesia", "measure_detail": "Mine production", "year": 2025, "value": 70.0},
        {**base, "commodity": "Nickel", "country": "Canada", "measure_detail": "Mine production", "year": 2025, "value": 20.0},
        {**base, "commodity": "Nickel", "country": O, "measure_detail": "Mine production", "year": 2025, "value": 10.0},
        {**base, "commodity": "Nickel", "country": W, "measure_detail": "Mine production: rounded", "year": 2025, "value": 100.0},
        # 다른 항목의 세계 합계는 섞이면 안 된다
        {**base, "commodity": "Nickel", "country": W, "measure_detail": "Mine production: other, rounded", "year": 2025, "value": 999.0},
        {**base, "commodity": "Lime", "country": "China", "measure_detail": "Production", "year": 2025, "value": 30.0,
         "is_critical_2025": False},
        {**base, "commodity": "Lime", "country": "India", "measure_detail": "Production", "year": 2025, "value": 30.0,
         "is_critical_2025": False},
        {**base, "commodity": "Lime", "country": "Japan", "measure_detail": "Production", "year": 2025, "value": 40.0,
         "is_critical_2025": False},
        {**base, "section": "Salient Statistics—United States", "measure": "net import reliance", "unit": "percent",
         "commodity": "Nickel", "country": "United States", "value": 41.0, "value_raw": "41", "year": 2025,
         "measure_detail": "Net import reliance as a percentage of apparent consumption"},
    ]
    return pd.DataFrame(rows)


CMAP = {"exclude": ["O"], "elements": {
    "Ni": {"name": "니켈", "critical": "Nickel", "world_bank": "Nickel",
           "stages": [{"commodity": "Nickel", "label": "광산", "production": "Mine production"}]},
    "Ca": {"name": "칼슘", "stages": [{"commodity": "Lime", "label": "석회", "production": "Production"}]},
}}


@pytest.fixture
def con() -> duckdb.DuckDBPyConnection:
    c = duckdb.connect()
    c.register("commodity_stats", stats_rows())
    c.register("critical_minerals", pd.DataFrame({"mineral": ["Nickel ", "Copper"]}))
    c.register("commodity_prices", pd.DataFrame({"commodity": "Nickel", "month": ["2025-08", "2026-08"],
                                                 "price": [15_000.0, 16_500.0], "unit": "$/mt"}))
    return c


def test_assess_reads_tables_and_judges(con: duckdb.DuckDBPyConnection) -> None:
    items = {it["element"]: it for it in supply.assess(
        con, ["Ni", "Ca", "Xx"], CMAP, {"commodity_stats", "critical_minerals", "commodity_prices"})}
    ni = items["Ni"]
    c = ni["stages"][0]["conc"]
    assert c.world == 100 and c.top_country == "Indonesia" and c.top_share == pytest.approx(0.7)
    assert c.hhi == pytest.approx(70**2 + 20**2)
    assert ni["critical"] and ni["status"] is Status.WARNING
    assert ni["stages"][0]["nir"]["text"] == "41%"
    assert ni["price"]["change"] == pytest.approx(0.1)
    ca = items["Ca"]
    assert not ca["critical"] and ca["status"] is Status.WARNING  # 1위 40% 이지만 HHI 3,400 > 2,500
    assert ca["stages"][0]["conc"].hhi == pytest.approx(3_400)
    assert not items["Xx"]["mapped"]


def test_assess_hhi_over_threshold_warns_even_if_not_critical(con: duckdb.DuckDBPyConnection) -> None:
    (ca,) = supply.assess(con, ["Ca"], CMAP, {"commodity_stats"})
    assert ca["why"] == ["HHI 3,400"]


def test_assess_falls_back_to_mcs_critical_flag(con: duckdb.DuckDBPyConnection) -> None:
    (ni,) = supply.assess(con, ["Ni"], CMAP, {"commodity_stats"})  # critical_minerals 없음
    assert ni["critical"] and ni["price"] is None


# ── 시험 진입점 ────────────────────────────────────────────────────────────


def ctx_of(formulas: dict[str, str | None]) -> Context:
    comps = [Resolved(ref=SubstanceRef.model_validate(ref), name=ref, formula=f) for ref, f in formulas.items()]
    recipe = Recipe.model_validate({"id": "rcp-t", "name": "t",
                                    "components": [{"ref": r, "amount": "1 mol"} for r in formulas]})
    return Context(recipe=recipe, comps=comps, registry=load_registry())


def test_elements_of_skips_excluded_and_formula_less(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(supply, "commodity_map", lambda: CMAP)
    ctx = ctx_of({"formula:NiO": "NiO", "formula:CaO": "CaO", "material:epoxy": None})
    assert supply.elements_of(ctx) == (["Ni", "Ca"], ["material:epoxy"])


def test_a9_pending_without_mcs(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(supply.db, "tables", lambda: [])
    out = supply.a9(ctx_of({"formula:NiO": "NiO"}))
    assert out.status is Status.PENDING and "usgs-mcs" in out.summary


def test_real_commodity_map_is_consistent() -> None:
    cmap = supply.commodity_map()
    assert {"H", "C", "N", "O"} <= set(cmap["exclude"])
    for el, spec in cmap["elements"].items():
        assert el not in cmap["exclude"]
        for s in spec["stages"]:
            assert s["commodity"] and s.get("label"), el


# ── 커넥터 정규화 ──────────────────────────────────────────────────────────


@pytest.mark.parametrize(("raw", "value", "flag"), [
    ("3,900,000", 3_900_000.0, None), ("—", 0.0, "zero"), ("W", None, "W"), ("E", None, "E"),
    (">140,000,000", 140_000_000.0, ">"), ("<1", 1.0, "<"), ("500–16,000", 8_250.0, "range"),
    ("330 - 390", 360.0, "range"), ("Large", None, "text"), ("", None, "NA"), ("432.3", 432.3, None),
])
def test_mcs_parse_value(raw: str, value: float | None, flag: str | None) -> None:
    assert usgs_mcs.parse_value(raw) == (value, flag)


def test_mcs_normalize_canonical_commodity_names() -> None:
    raw = pd.DataFrame({
        "MCS chapter": ["TITANIUM"] * 3, "Section": ["World Production"] * 3,
        "Commodity": ["TiO2 Pigment", "Tio2 Pigment", "TiO2 Pigment"], "Country": ["China"] * 3,
        "Statistics": ["Production"] * 3, "Statistics_detail": ["Pigment"] * 3, "Unit": ["metric tons"] * 3,
        "Year": ["2024", "2025", "2021–24"], "Value": ["10", "W", "5"], "Notes": ["", "n", ""],
        "Is critical mineral 2025": ["No", "Yes", ""], "Other notes": [""] * 3,
    })
    df = usgs_mcs.normalize(raw)
    assert df["commodity"].unique().tolist() == ["TiO2 Pigment"]
    assert df["year"].tolist()[:2] == [2024, 2025] and pd.isna(df["year"].iloc[2])
    assert df["measure"].iloc[0] == "production" and df["record_id"].iloc[0].endswith(":2")
    assert df["is_critical_2025"].tolist()[:2] == [False, True] and pd.isna(df["is_critical_2025"].iloc[2])


def test_pink_sheet_normalize_long_format() -> None:
    rows = [
        ("World Bank Commodity Price Data (The Pink Sheet)", None, None),
        ("Updated on September 02, 2026", None, None),
        (None, "Copper", "Potassium chloride **"),
        (None, "($/mt)", "($/mt)"),
        ("2026M07", 13_543, "…"),
        ("2026M08", 14_326, 386.9),
    ]
    assert world_bank_pink_sheet._updated(rows) == "2026-09-02"
    df = world_bank_pink_sheet.normalize(rows)
    assert len(df) == 3  # '…' 칸은 적재하지 않는다
    k = df[df["commodity"] == "Potassium chloride"].iloc[0]
    assert (k["month"], k["price"], k["unit"], k["series"]) == ("2026-08", 386.9, "$/mt", "Potassium chloride **")
    assert df["record_id"].is_unique


def test_mp_row_license_tags_gnome_as_noncommercial() -> None:
    from msl.connectors import materials_project as mp

    doc = {"material_id": "mp-x", "symmetry": {"symbol": "Fm-3m", "number": 225}, "database_IDs": {"icsd": ["icsd-1"]},
           "builder_meta": {"license": "BY-C", "batch_id": "mp_2018"}}
    ok = mp._row(doc, {mp.GNOME_BATCH})
    assert ok["license"] == "CC-BY-4.0" and ok["has_icsd"] and ok["spacegroup_number"] == 225
    by_nc = {**doc, "builder_meta": {"license": "BY-NC", "batch_id": "other"}}
    gnome = {**doc, "builder_meta": {"license": None, "batch_id": mp.GNOME_BATCH}, "database_IDs": {}}
    assert mp._row(by_nc, set())["license"] == mp._row(gnome, {mp.GNOME_BATCH})["license"] == "CC-BY-NC-4.0"
    assert not mp._row(gnome, {mp.GNOME_BATCH})["has_icsd"]
