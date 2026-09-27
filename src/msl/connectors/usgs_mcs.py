"""USGS Mineral Commodity Summaries 2026 데이터 릴리스 (CC0-1.0) — 품목별 미국·세계 통계 (계획서 4.4절).

DOI 10.5066/P1WKQ63T (ScienceBase 696a75d5d4be0228872d3bf8, 2026-05-27 개정).
ScienceBase 는 Cloudflare 브라우저 검사로 스크립트 요청을 403 으로 막는다. fetch() 는 먼저 직접 받기를 시도하고,
막히면 브라우저로 받아 둘 파일 목록을 알려 준다. 어느 쪽이든 ScienceBase 가 공개한 MD5 와 대조한다.

CSV 는 cp1252 인코딩(Section 의 em dash 등)이다. 값 기호는 메타데이터 정의를 따른다:
W 비공개, NA 자료 없음, XX 세계 합계 산출 불가, E 순수출, s 극소, — 0. '>n'·'<n'·'a–b' 범위는 숫자를 살리고 value_flag 에 남긴다.
"""

from __future__ import annotations

import hashlib
import re
import urllib.error
from pathlib import Path

import pandas as pd

from msl.db import download, raw_dir, record_checksums, with_provenance, write_table

SOURCE = "usgs-mcs"
TABLES = ["commodity_stats"]
VERSION = "2026"
LICENSE = "CC0-1.0"
MAIN = "MCS2026_Commodities_Data.csv"
_SB = "https://www.sciencebase.gov/catalog/file/get/"
# 파일명 → (ScienceBase 다운로드 URL, ScienceBase 공개 MD5 — None 이면 공개값 없음)
FILES = {
    MAIN: (_SB + "696a75d5d4be0228872d3bf8?f=__disk__eb%2F4d%2Fd1%2Feb4dd129d4f08cfba010a8b11cbbd2d88bbefdcb",
           "36185ff3742087e1dd90c52fe634fe12"),
    "MCS2026_Commodities_Data.xml": (
        _SB + "69837e43b66b01367d7ec7c7?f=__disk__59%2F9f%2F78%2F599f78dcdae6cfc820ec8410e8a0e684c28be9d2", None),
    "MCS2026_T6_Critical_Minerals_End_Use.csv": (
        _SB + "696a75d5d4be0228872d3bf8?f=__disk__41%2F88%2F81%2F4188812fc589685e01e9e9761155ee89e97e6f8a",
        "802cee69ef5991c8fc01a1c659199af1"),
    "MCS2026_T7_Critical_Minerals_Salient.csv": (
        _SB + "696a75d5d4be0228872d3bf8?f=__disk__72%2Fec%2F18%2F72ec18177f3ba33eaa6e6cd09ce433456837b1c6",
        "e185664f1e5d77b7cd3cc3487d76e5cb"),
    "MCS2026_Fig2_Net_Import_Reliance.csv": (
        _SB + "696a75d5d4be0228872d3bf8?f=__disk__90%2F85%2Fd4%2F9085d4823bf86edb0f54dd63eb8023c11d097c61",
        "ac4593283060f39d39caa2af2bcd37c9"),
    "MCS2026_Fig3_Major_Import_Sources.csv": (
        _SB + "696a75d5d4be0228872d3bf8?f=__disk__01%2F9e%2F39%2F019e39aab1e394ef119ffe69b225b951face974f",
        "013c72116dcbc3e492f17bf5fcef965f"),
}
FLAGS = {"W": "W", "NA": "NA", "XX": "XX", "E": "E", "s": "s"}
_NUM = r"\d[\d,]*(?:\.\d+)?"


class ManualDownloadRequired(RuntimeError):
    """ScienceBase 가 스크립트 요청을 막아 사람이 브라우저로 받아야 할 때."""


def fetch() -> Path:
    d = raw_dir(SOURCE, VERSION)
    blocked = []
    for name, (url, _) in FILES.items():
        try:
            download(url, d / name)
        except urllib.error.HTTPError as exc:
            if exc.code != 403:
                raise
            blocked.append(name)
    if blocked:
        listing = "\n".join(f"  {n}  {FILES[n][0]}" for n in blocked)
        raise ManualDownloadRequired(
            f"ScienceBase 가 Cloudflare 검사로 막음 — 브라우저로 받아 {d}/ 에 두고 다시 실행:\n{listing}")
    for name, (_, md5) in FILES.items():
        got = hashlib.md5((d / name).read_bytes()).hexdigest()
        if md5 and got != md5:
            raise ValueError(f"{name}: MD5 {got} ≠ ScienceBase {md5}")
    record_checksums(d)
    return d


def parse_value(raw: str) -> tuple[float | None, str | None]:
    """MCS 값 문자열 → (숫자, 기호). 기호는 None(보통 숫자)·W·NA·XX·E·s·zero·>·<·range·text."""
    v = (raw or "").strip()
    if not v:
        return None, "NA"
    if v in FLAGS:
        return None, FLAGS[v]
    if v == "—":
        return 0.0, "zero"
    if re.fullmatch(_NUM, v):
        return float(v.replace(",", "")), None
    if m := re.fullmatch(rf"([<>])\s*({_NUM})", v):
        return float(m.group(2).replace(",", "")), m.group(1)
    if m := re.fullmatch(rf"({_NUM})\s*[–-]\s*({_NUM})", v):
        lo, hi = (float(x.replace(",", "")) for x in m.groups())
        return (lo + hi) / 2, "range"
    return None, "text"


def _canonical(names: pd.Series) -> pd.Series:
    """대소문자만 다른 품목명(TiO2 Pigment / Tio2 Pigment 등)을 가장 흔한 표기로 모은다."""
    key = names.str.strip().str.lower()
    best = names.str.strip().groupby(key).agg(lambda s: s.value_counts().index[0])
    return key.map(best)


def normalize(raw: pd.DataFrame) -> pd.DataFrame:
    parsed = raw["Value"].map(parse_value)
    year = pd.to_numeric(raw["Year"].str.strip(), errors="coerce").astype("Int64")
    crit = raw["Is critical mineral 2025"].str.strip().map({"Yes": True, "No": False})
    return pd.DataFrame({
        "record_id": [f"{MAIN}:{i + 2}" for i in range(len(raw))],  # CSV 줄 번호 (헤더 = 1)
        "chapter": raw["MCS chapter"].str.strip(),
        "section": raw["Section"].str.strip(),
        "commodity": _canonical(raw["Commodity"]),
        "country": raw["Country"].str.strip(),
        "measure": raw["Statistics"].str.strip().str.lower(),
        "measure_detail": raw["Statistics_detail"].str.strip(),
        "year": year,
        "period": raw["Year"].str.strip(),
        "value": [p[0] for p in parsed],
        "value_flag": [p[1] for p in parsed],
        "value_raw": raw["Value"].str.strip(),
        "unit": raw["Unit"].str.strip(),
        "is_critical_2025": crit.astype("boolean"),
        "notes": raw["Notes"].str.strip().replace({"": None}),
        "other_notes": raw["Other notes"].str.strip().replace({"": None}),
    })


def read_raw(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, encoding="cp1252", dtype=str, keep_default_na=False)


def load() -> dict[str, dict[str, int]]:
    src = raw_dir(SOURCE, VERSION) / MAIN
    if not src.exists():
        fetch()
    df = normalize(read_raw(src))
    df = with_provenance(df, source=SOURCE, version=VERSION, license=LICENSE, method="compiled", id_column="record_id")
    return {"commodity_stats": write_table("commodity_stats", df)}
