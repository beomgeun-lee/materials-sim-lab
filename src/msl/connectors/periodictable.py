"""periodictable (퍼블릭 도메인) — 배포용 원소 레이어의 기본 소스 (결정 D10).

파이썬 패키지에 들어 있는 값을 CSV 로 떠서 원본 스냅샷으로 삼는다.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import periodictable as pt

from msl.db import raw_dir, record_checksums, with_provenance, write_table

SOURCE = "periodictable"
TABLES = ["elements"]
VERSION = pt.__version__
LICENSE = "public-domain"


def fetch() -> Path:
    d = raw_dir(SOURCE, VERSION)
    rows = []
    for el in pt.elements:
        if el.number == 0:  # 중성자(n) 제외
            continue
        rows.append({
            "Z": el.number, "symbol": el.symbol, "name": el.name,
            "atomic_mass": el.mass,
            "density_g_cm3": el.density if el.density else None,
        })
    pd.DataFrame(rows).to_csv(d / "elements.csv", index=False)
    record_checksums(d)
    return d


def load() -> dict[str, dict[str, int]]:
    src = raw_dir(SOURCE, VERSION) / "elements.csv"
    if not src.exists():
        fetch()
    df = pd.read_csv(src)
    df = with_provenance(df, source=SOURCE, version=VERSION, license=LICENSE, method="compiled", id_column="symbol")
    return {"elements": write_table("elements", df)}
