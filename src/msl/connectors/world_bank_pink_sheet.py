"""World Bank Commodity Price Data — Pink Sheet (CC-BY-4.0). 월간 명목 달러 가격, 1960년~.

파일 주소가 발간마다 바뀌므로 commodity-markets 페이지에서 CMO-Historical-Data-Monthly.xlsx 링크를 찾는다.
버전은 시트 머리의 'Updated on …' 날짜(YYYY-MM-DD)이고, 원본은 data/raw/world-bank-pink-sheet/{날짜}/ 에 둔다.
'…' 등 숫자가 아닌 칸(자료 없음)은 적재하지 않는다.
"""

from __future__ import annotations

import datetime as dt
import re
import shutil
import ssl
import urllib.request
from pathlib import Path

import certifi
import openpyxl
import pandas as pd

from msl.db import RAW, download, raw_dir, record_checksums, with_provenance, write_table

SOURCE = "world-bank-pink-sheet"
TABLES = ["commodity_prices"]
LICENSE = "CC-BY-4.0"
PAGE = "https://www.worldbank.org/en/research/commodity-markets"
FILE = "CMO-Historical-Data-Monthly.xlsx"
SHEET = "Monthly Prices"
UA = {"User-Agent": "materials-sim-lab/0.0.1 (+research)"}


def _file_url() -> str:
    ctx = ssl.create_default_context(cafile=certifi.where())
    with urllib.request.urlopen(urllib.request.Request(PAGE, headers=UA), timeout=60, context=ctx) as resp:
        html = resp.read().decode("utf-8", "replace")
    urls = re.findall(rf"https://thedocs\.worldbank\.org/[^\"' ]+/{FILE}", html)
    if not urls:
        raise RuntimeError(f"{PAGE} 에서 {FILE} 링크를 찾지 못함")
    return urls[0]


def _rows(path: Path) -> list[tuple]:
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        return list(wb[SHEET].iter_rows(values_only=True))
    finally:
        wb.close()


def _updated(rows: list[tuple]) -> str:
    for row in rows[:6]:
        if isinstance(row[0], str) and row[0].startswith("Updated on"):
            return dt.datetime.strptime(row[0].removeprefix("Updated on").strip(), "%B %d, %Y").date().isoformat()
    raise ValueError(f"{SHEET}: 'Updated on' 날짜 줄이 없음")


def fetch() -> Path:
    url = _file_url()
    tmp = RAW / SOURCE / "_incoming" / FILE
    download(url, tmp, overwrite=True)
    d = raw_dir(SOURCE, _updated(_rows(tmp)))
    shutil.move(tmp, d / FILE)
    shutil.rmtree(tmp.parent, ignore_errors=True)
    (d / "SOURCE_URL").write_text(url + "\n", encoding="utf-8")
    record_checksums(d)
    return d


def latest_version() -> str | None:
    versions = sorted(p.name for p in (RAW / SOURCE).glob("????-??-??") if (p / FILE).exists())
    return versions[-1] if versions else None


def normalize(rows: list[tuple]) -> pd.DataFrame:
    """'Monthly Prices' 시트(머리 5번째 줄 = 품목, 6번째 = 단위, 이후 '1960M01' 형식 월) → 긴 형식."""
    head = next(i for i, r in enumerate(rows) if r[0] is None and isinstance(r[1], str) and r[1].strip())
    names, units = rows[head], rows[head + 1]
    out = []
    for r in rows[head + 2:]:
        if not (isinstance(r[0], str) and re.fullmatch(r"\d{4}M\d{2}", r[0])):
            continue
        month = f"{r[0][:4]}-{r[0][5:]}"
        for name, unit, v in zip(names[1:], units[1:], r[1:]):
            if not name or not isinstance(v, (int, float)):
                continue
            series = str(name).strip()
            out.append({
                "commodity": series.removesuffix("**").strip(), "series": series, "month": month,
                "price": float(v), "unit": str(unit or "").strip().strip("()"),
            })
    df = pd.DataFrame(out)
    df["record_id"] = df["commodity"] + "|" + df["month"]
    return df


def load() -> dict[str, dict[str, int]]:
    version = latest_version()
    if version is None:
        version = fetch().name
    df = normalize(_rows(raw_dir(SOURCE, version) / FILE))
    df = with_provenance(df, source=SOURCE, version=version, license=LICENSE, method="compiled", id_column="record_id")
    return {"commodity_prices": write_table("commodity_prices", df)}
