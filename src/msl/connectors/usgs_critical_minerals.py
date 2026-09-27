"""USGS 2025 핵심광물 목록 (퍼블릭 도메인) — 60종 (연방관보 2025-11-07).

목록 자체는 기계 판독 파일이 따로 없어서, 같은 USGS 가 MCS 2026 데이터 릴리스에 실은
표 T6(MCS2026_T6_Critical_Minerals_End_Use.csv — '2025 Final List' 60종과 주 용도)을 원본으로 삼는다.
설명 페이지(about-2025-list-critical-minerals)도 함께 떠 둔다.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd

from msl.connectors import usgs_mcs
from msl.db import download, raw_dir, record_checksums, with_provenance, write_table

SOURCE = "usgs-critical-minerals"
TABLES = ["critical_minerals"]
VERSION = "2025"
LICENSE = "public-domain"
LIST_NAME = "USGS 2025"
TABLE_FILE = "MCS2026_T6_Critical_Minerals_End_Use.csv"
PAGE = "https://www.usgs.gov/programs/mineral-resources-program/science/about-2025-list-critical-minerals"


def fetch() -> Path:
    d = raw_dir(SOURCE, VERSION)
    mcs = raw_dir(usgs_mcs.SOURCE, usgs_mcs.VERSION) / TABLE_FILE
    if not mcs.exists():
        usgs_mcs.fetch()
    shutil.copy2(mcs, d / TABLE_FILE)
    download(PAGE, d / "about-2025-list-critical-minerals.html")
    record_checksums(d)
    return d


def load() -> dict[str, dict[str, int]]:
    src = raw_dir(SOURCE, VERSION) / TABLE_FILE
    if not src.exists():
        fetch()
    raw = pd.read_csv(src, encoding="utf-8-sig", dtype=str, keep_default_na=False)
    df = pd.DataFrame({
        "mineral": raw["Critical_Mineral"].str.strip(),
        "list": LIST_NAME,
        "primary_applications": raw["Primary_Applications"].str.strip(),
        "note": raw["Category_Note"].str.strip().replace({"": None}),
    })
    df = with_provenance(df, source=SOURCE, version=VERSION, license=LICENSE, method="compiled", id_column="mineral")
    return {"critical_minerals": write_table("critical_minerals", df)}
