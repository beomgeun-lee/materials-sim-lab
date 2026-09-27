"""A9 공급·경제성 — USGS MCS 2026 생산 집중도·미국 순수입 의존도, USGS 2025 핵심광물, World Bank 월간 가격 (L0).

원소 → 공급 품목은 kb/commodity_map.yaml 규칙으로 정한다. 조회는 DuckDB 적재 테이블
(commodity_stats, critical_minerals, commodity_prices)만 쓴다 — 네트워크 호출 없음.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from typing import Any, Iterable

import duckdb
import pandas as pd
import yaml

from msl import db
from msl.assays.base import Context, Outcome, not_applicable, pending, val
from msl.env import KB_DIR
from msl.schema.result import Fidelity, Status, ValueKind

ENGINE = "msl-supply-risk"
MAP_FILE = KB_DIR / "commodity_map.yaml"
HHI_WARN = 2500.0  # 고집중 시장 기준 (HHI, 0~10,000)
TOP_WARN = 0.50  # 1위 국가 점유율
WORLD, OTHERS, US = "World total", "Other countries", "United States"
COLUMNS = ["원소", "품목 (USGS)", "단계", "기준 연도", "세계 생산", "1위 국가", "1위 점유율", "HHI",
           "미국 순수입 의존도", "핵심광물 (USGS 2025)", "가격 (World Bank)", "12개월 변화", "판정"]


@cache
def commodity_map() -> dict[str, Any]:
    return yaml.safe_load(MAP_FILE.read_text(encoding="utf-8"))


# ── 계산 (순수 함수) ────────────────────────────────────────────────────────


def hhi(shares: Iterable[float]) -> float:
    """허핀달-허쉬만 지수. shares 는 0~1 비율, 결과는 0~10,000."""
    return float(sum((100 * s) ** 2 for s in shares))


@dataclass
class Concentration:
    year: int
    unit: str
    world: float
    world_reported: bool  # False 면 공표 국가·기타 합으로 대신함
    top_country: str
    top_share: float
    hhi: float
    n_countries: int


def concentration(rows: pd.DataFrame) -> Concentration | None:
    """세계 생산 표 한 항목(열: country, year, value, unit)에서 최근 연도의 집중도를 구한다.

    'Other countries' 는 국가가 아니라 HHI·1위 계산에서 빼고 세계 합계에만 들어간다 (HHI 는 하한).
    """
    rows = rows.dropna(subset=["value", "year"])
    if rows.empty:
        return None
    year = int(rows["year"].max())
    r = rows[rows["year"] == year]
    countries = r[~r["country"].isin([WORLD, OTHERS]) & (r["value"] > 0)]
    if countries.empty:
        return None
    reported = r.loc[(r["country"] == WORLD) & (r["value"] > 0), "value"]
    world = float(reported.iloc[0]) if len(reported) else float(r.loc[r["country"] != WORLD, "value"].sum())
    shares = (countries.groupby("country")["value"].sum() / world).clip(upper=1.0)
    unit = str(r["unit"].dropna().iloc[0]) if "unit" in r and r["unit"].notna().any() else ""
    return Concentration(year=year, unit=unit, world=world, world_reported=bool(len(reported)),
                         top_country=str(shares.idxmax()), top_share=float(shares.max()), hhi=hhi(shares),
                         n_countries=len(shares))


def judge(critical: bool, top_share: float | None, hhi_value: float | None) -> tuple[Status, list[str]]:
    """핵심광물이거나 1위 점유율 > 50% 이거나 HHI > 2,500 이면 주의."""
    why = []
    if critical:
        why.append("핵심광물")
    if top_share is not None and top_share > TOP_WARN:
        why.append(f"1위 점유율 {top_share:.0%}")
    if hhi_value is not None and hhi_value > HHI_WARN:
        why.append(f"HHI {hhi_value:,.0f}")
    return (Status.WARNING if why else Status.OK), why


def pick_nir(rows: pd.DataFrame, prefer: str | None = None) -> dict[str, Any] | None:
    """미국 순수입 의존도 행(열: measure_detail, year, value, value_flag, value_raw)에서 최근 연도 하나를 고른다."""
    rows = rows.dropna(subset=["year"])
    if rows.empty:
        return None
    r = rows[rows["year"] == rows["year"].max()]
    detail = r["measure_detail"].str.lower()
    for mask in ([detail.str.contains(prefer.lower(), regex=False)] if prefer else []) + [detail.str.contains("total")]:
        if mask.any():
            r = r[mask]
            break
    row = r.iloc[0]
    flag = row.get("value_flag") if isinstance(row.get("value_flag"), str) else None
    value = None if pd.isna(row["value"]) else float(row["value"])
    text = "순수출" if flag == "E" else f"{row['value_raw']}%"
    return {"year": int(row["year"]), "value": value, "flag": flag, "text": text}


def price_change(rows: pd.DataFrame) -> dict[str, Any] | None:
    """월별 가격(열: month 'YYYY-MM', price, unit) → 최신월 가격과 12개월 전 대비 변화율."""
    rows = rows.dropna(subset=["price"]).sort_values("month")
    if rows.empty:
        return None
    last = rows.iloc[-1]
    y, m = int(last["month"][:4]), int(last["month"][5:7])
    ago = rows[rows["month"] == f"{y - 1:04d}-{m:02d}"]
    change = float(last["price"] / ago["price"].iloc[0] - 1) if len(ago) and ago["price"].iloc[0] else None
    return {"month": last["month"], "price": float(last["price"]), "unit": last["unit"], "change": change}


# ── 조회 ──────────────────────────────────────────────────────────────────


def _production(con: duckdb.DuckDBPyConnection, commodity: str, detail: str) -> pd.DataFrame:
    """국가 행은 detail 과 같은 항목, 세계 합계 행은 detail + ': rounded' / ', rounded' 도 받는다."""
    return con.execute(
        """SELECT country, year, value, unit FROM commodity_stats
           WHERE commodity = $c AND measure = 'production' AND section LIKE 'World%'
             AND (lower(measure_detail) = lower($d)
                  OR (country = $w AND lower(measure_detail) IN (lower($d) || ': rounded', lower($d) || ', rounded')))""",
        {"c": commodity, "d": detail, "w": WORLD}).df()


def _nir(con: duckdb.DuckDBPyConnection, commodity: str) -> pd.DataFrame:
    return con.execute(
        """SELECT measure_detail, year, value, value_flag, value_raw FROM commodity_stats
           WHERE commodity = ? AND measure = 'net import reliance' AND country = ? ORDER BY measure_detail""",
        [commodity, US]).df()


def _critical(con: duckdb.DuckDBPyConnection, tables: set[str], name: str | None, commodities: list[str]) -> bool:
    if "critical_minerals" in tables:
        return bool(name) and con.execute(
            "SELECT count(*) FROM critical_minerals WHERE lower(trim(mineral)) = lower(?)", [name]).fetchone()[0] > 0
    # 목록 테이블이 없으면 MCS 의 'Is critical mineral 2025' 열로 대신한다
    return bool(commodities) and con.execute(
        "SELECT coalesce(bool_or(is_critical_2025), false) FROM commodity_stats WHERE commodity IN (SELECT unnest(?))",
        [commodities]).fetchone()[0]


def _prices(con: duckdb.DuckDBPyConnection, commodity: str) -> pd.DataFrame:
    return con.execute(
        "SELECT month, price, unit FROM commodity_prices WHERE commodity = ? ORDER BY month DESC LIMIT 13",
        [commodity]).df()


def assess(con: duckdb.DuckDBPyConnection, elements: list[str], cmap: dict[str, Any],
           tables: set[str]) -> list[dict[str, Any]]:
    """원소마다 판정과 단계별 수치. tables 는 con 에 있는 테이블 이름."""
    out = []
    for el in elements:
        spec = cmap.get("elements", {}).get(el)
        if spec is None:
            out.append({"element": el, "mapped": False, "status": Status.OK, "why": [], "stages": []})
            continue
        stages = spec.get("stages") or []
        critical = _critical(con, tables, spec.get("critical"), [s["commodity"] for s in stages])
        price = None
        if spec.get("world_bank") and "commodity_prices" in tables:
            price = price_change(_prices(con, spec["world_bank"]))
        rows = []
        for s in stages:
            conc = concentration(_production(con, s["commodity"], s["production"])) if s.get("production") else None
            rows.append({"commodity": s["commodity"], "label": s.get("label", ""), "conc": conc,
                         "nir": pick_nir(_nir(con, s["commodity"]), s.get("nir"))})
        concs = [r["conc"] for r in rows if r["conc"]]
        top = max(concs, key=lambda c: c.top_share, default=None)
        worst_hhi = max((c.hhi for c in concs), default=None)
        status, why = judge(critical, top.top_share if top else None, worst_hhi)
        out.append({"element": el, "name": spec.get("name", el), "mapped": True, "critical": critical,
                    "critical_name": spec.get("critical"), "world_bank": spec.get("world_bank"), "price": price,
                    "stages": rows, "status": status, "why": why})
    return out


# ── 시험 ──────────────────────────────────────────────────────────────────


def elements_of(ctx: Context) -> tuple[list[str], list[str]]:
    """레시피 성분의 원소(첫 등장 순, 제외 원소 빼고)와 화학식이 없어 건너뛴 성분 이름."""
    exclude = set(commodity_map().get("exclude", []))
    seen: list[str] = []
    skipped = []
    for c in ctx.comps:
        comp = c.composition
        if comp is None:
            skipped.append(c.name)
            continue
        for el in comp.elements:
            sym = el.symbol
            if sym not in exclude and sym not in seen:
                seen.append(sym)
    return seen, skipped


def _fmt(x: float) -> str:
    return f"{x:,.0f}" if x >= 100 else f"{x:,.3g}"


def _table_rows(items: list[dict[str, Any]]) -> list[list[Any]]:
    rows = []
    for it in items:
        if not it["mapped"]:
            rows.append([it["element"], "매핑 없음", "—", "—", "—", "—", "—", "—", "—", "—", "—", "—", "—"])
            continue
        verdict = "주의 (" + ", ".join(it["why"]) + ")" if it["why"] else "보통"
        crit = f"해당 ({it['critical_name']})" if it["critical"] else "비해당"
        p = it["price"]
        price = f"{p['price']:,.2f} {p['unit']} ({p['month']}, {it['world_bank']})" if p else "—"
        change = f"{p['change']:+.1%}" if p and p["change"] is not None else "—"
        stages = it["stages"] or [{"commodity": "—", "label": "—", "conc": None, "nir": None}]
        for i, s in enumerate(stages):
            c, n = s["conc"], s["nir"]
            rows.append([
                f"{it['element']} ({it['name']})" if i == 0 else "", s["commodity"], s["label"],
                c.year if c else (n["year"] if n else "—"),
                f"{_fmt(c.world)} {c.unit}" + ("" if c.world_reported else " (국가 합)") if c else "—",
                c.top_country if c else "—", f"{c.top_share:.0%}" if c else "—", f"{c.hhi:,.0f}" if c else "—",
                n["text"] if n else "—",
                crit if i == 0 else "", price if i == 0 else "", change if i == 0 else "", verdict if i == 0 else "",
            ])
    return rows


def _values(items: list[dict[str, Any]]) -> list:
    values = []
    for it in items:
        if not it["mapped"]:
            continue
        el = it["element"]
        values.append(val(f"{el} · 공급 위험", "주의 — " + ", ".join(it["why"]) if it["why"] else "보통",
                          kind=ValueKind.FLAG))
        for s in it["stages"]:
            if c := s["conc"]:
                values.append(val(f"{el} {s['label']} · 1위 {c.top_country} 점유율", round(100 * c.top_share, 1), "%"))
                values.append(val(f"{el} {s['label']} · HHI", round(c.hhi)))
        if n := next((s["nir"] for s in it["stages"] if s["nir"]), None):
            values.append(val(f"{el} · 미국 순수입 의존도 ({n['year']})", n["text"], kind=ValueKind.TEXT))
        if p := it["price"]:
            values.append(val(f"{el} · 가격 {p['month']}", p["price"], p["unit"]))
            if p["change"] is not None:
                values.append(val(f"{el} · 가격 12개월 변화", round(100 * p["change"], 1), "%"))
    return values


def _version(con: duckdb.DuckDBPyConnection, table: str) -> str | None:
    try:
        return str(con.execute(f"SELECT any_value(source_version) FROM {table}").fetchone()[0])
    except duckdb.Error:
        return None


def a9(ctx: Context) -> Outcome:
    tables = set(db.tables())
    if "commodity_stats" not in tables:
        return pending("USGS MCS 자료가 적재되지 않음 — msl db load usgs-mcs 필요", ENGINE)
    elements, skipped = elements_of(ctx)
    if not elements:
        return not_applicable("평가할 원소 없음 (H·C·N·O 등 제외 원소만 있거나 화학식이 없는 성분)", ENGINE)
    con = db.connect()
    try:
        items = assess(con, elements, commodity_map(), tables)
        versions = {t: _version(con, t) for t in ("commodity_stats", "critical_minerals", "commodity_prices") if t in tables}
    finally:
        con.close()

    flagged = [it for it in items if it["why"]]
    unmapped = [it["element"] for it in items if not it["mapped"]]
    warnings = []
    if unmapped:
        warnings.append(f"공급 품목 매핑이 없는 원소: {', '.join(unmapped)} (kb/commodity_map.yaml)")
    if skipped:
        warnings.append(f"화학식이 없어 평가에서 뺀 성분: {', '.join(skipped)}")
    if "commodity_prices" not in tables:
        warnings.append("가격 자료 없음 — msl db load world-bank-pink-sheet 필요")
    summary = (f"{len(items)}개 원소 중 {len(flagged)}개 공급 주의 — "
               + "; ".join(f"{it['element']}({', '.join(it['why'])})" for it in flagged)
               if flagged else f"{len(items)}개 원소 모두 핵심광물 아님·생산 집중 낮음")
    sources = [("usgs-mcs", versions.get("commodity_stats"))]
    if "critical_minerals" in versions:
        sources.append(("usgs-critical-minerals", versions["critical_minerals"]))
    if "commodity_prices" in versions:
        sources.append(("world-bank-pink-sheet", versions["commodity_prices"]))
    return Outcome(
        status=Status.WARNING if flagged else Status.OK, fidelity=Fidelity.L0, engine=ENGINE, engine_version="v0",
        conditions_basis="USGS MCS 2026 세계 생산(최근 추정 연도)·미국 통계, World Bank 월평균 명목 가격 최신월",
        values=_values(items), summary=summary, sources=sources, warnings=warnings,
        caveats=[
            f"판정 기준: 핵심광물(USGS 2025) 또는 1위 국가 점유율 > {TOP_WARN:.0%} 또는 HHI > {HHI_WARN:,.0f} 이면 주의",
            "HHI 는 국가별로 공표된 생산량으로 계산 — 'Other countries' 묶음과 비공개(W) 값은 빠지므로 하한값",
            "원소 → 품목은 kb/commodity_map.yaml 규칙. 화합물 원료(예: CaCO3)도 원소 품목(석회)으로 평가하며 실제 조달 원료와 다를 수 있음",
            "순수입 의존도는 미국 기준 — 국내 수급(한국광해광업공단)은 아직 연결되지 않음",
            "가격은 World Bank 가 다루는 품목만 (금속 10종·비료 등)",
        ],
        data={"table": {"columns": COLUMNS, "rows": _table_rows(items)}},
    )
