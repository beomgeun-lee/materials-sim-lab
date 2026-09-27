"""COD (CC0) — 광물명이 달린 실험 결정 구조의 메타데이터 (계획서 4.1절 X 레이어, 결정 D4).

CIF 전량 미러는 범위 밖이다. 광물명(`_cod_mineral`)이 있는 항목과, 광물명 없이 commonname 만 있는 무기물 항목을
OPTIMADE 로 받는다.
- REST 검색(/cod/result)에는 광물 필터가 없어 전량(53만 건)을 받거나 이름별로 6천 번 질의해야 한다.
- OPTIMADE 는 `_cod_mineral IS KNOWN` 필터와 `response_fields` 로 필요한 열만 준다 (페이지 상한 100건).
- commonname 만 있는 항목(5.3만 건)은 대부분 유기물이라 C·H 를 함께 가진 것은 이름에 'ite' 가 있을 때만 받는다
  (말라카이트·웨웰라이트 같은 탄산염·유기 광물). 이름이 IMA 광물명과 맞는지는 ima_cnmnc.match() 가 본다.
응답 JSON 을 받은 그대로 data/raw/cod/{조회일}/{page|common}-{offset}.json 에 쌓는다. 끊기면 다시 돌려 이어받는다.
서버 예의: 요청은 하나씩, 사이에 1초 이상 쉰다 (User-Agent 는 msl.db.download 가 붙인다).
"""

from __future__ import annotations

import datetime as dt
import json
import re
import time
import urllib.error
import urllib.parse
import warnings
from functools import cache
from pathlib import Path

import pandas as pd
from pymatgen.core import Composition

from msl.db import download, now, raw_dir, record_checksums, with_provenance, write_table

SOURCE = "cod"
TABLES = ["structures_exp"]
LICENSE = "CC0-1.0"
BASE = "https://www.crystallography.net/cod/optimade/v1/structures"
FILTER = "_cod_mineral IS KNOWN"
COMMON_FILTER = ('_cod_mineral IS UNKNOWN AND _cod_commonname IS KNOWN AND '
                 '(NOT (elements HAS ALL "C","H") OR _cod_commonname CONTAINS "ite")')
QUERIES = {"page": FILTER, "common": COMMON_FILTER}  # 원본 파일 이름 앞머리 → 필터
FIELDS = [
    "_cod_mineral", "_cod_commonname", "_cod_calcformula", "_cod_sg", "_cod_sgnumber",
    "_cod_a", "_cod_b", "_cod_c", "_cod_alpha", "_cod_beta", "_cod_gamma", "_cod_vol",
    "_cod_celltemp", "_cod_cellpressure", "_cod_year", "_cod_doi", "_cod_robs", "_cod_rall",
    "_cod_flags", "_cod_duplicateof", "_cod_status", "_cod_method", "_cod_svnrevision",
    "chemical_formula_reduced", "elements",
]
PAGE = 100  # 서버 상한 (page_limit > 100 이면 403)
PAUSE = 1.0  # 요청 사이 최소 간격(초)


def _snapshots() -> list[Path]:
    return sorted(p for p in raw_dir(SOURCE, "_").parent.glob("*") if p.is_dir() and list(p.glob("page-*.json")))


def _url(offset: int, flt: str = FILTER) -> str:
    query = {"filter": flt, "page_limit": PAGE, "page_offset": offset, "response_fields": ",".join(FIELDS)}
    return f"{BASE}?{urllib.parse.urlencode(query)}"


def _get(url: str, dest: Path) -> Path:
    for attempt in range(4):
        try:
            return download(url, dest)
        except (urllib.error.URLError, TimeoutError):
            if attempt == 3:
                raise
            time.sleep(10 * (attempt + 1))
    return dest


def _fetch_query(d: Path, prefix: str, flt: str) -> dict:
    offset = 0
    while True:
        dest = d / f"{prefix}-{offset:06d}.json"
        if not dest.exists():
            _get(_url(offset, flt), dest)
            time.sleep(PAUSE)
        meta = json.loads(dest.read_text(encoding="utf-8"))["meta"]
        if not meta.get("more_data_available"):
            break
        offset += PAGE
    return {"filter": flt, "data_returned": meta.get("data_returned"), "data_available": meta.get("data_available")}


def fetch() -> Path:
    d = raw_dir(SOURCE, dt.date.today().isoformat())
    started = now()
    queries = {prefix: _fetch_query(d, prefix, flt) for prefix, flt in QUERIES.items()}
    (d / "query.json").write_text(json.dumps({
        "endpoint": BASE, "queries": queries, "response_fields": FIELDS, "page_limit": PAGE,
        "started_at": started, "retrieved_at": now(),
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    record_checksums(d)
    return d


# ── 정규화 ────────────────────────────────────────────────────────────────

_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*\.?\d*)")


def parse_formula(text: str | None) -> dict[str, float]:
    """COD 화학식 '- Ca C O3 -' → {'Ca': 1, 'C': 1, 'O': 3}. 중수소 D·삼중수소 T 는 H 로 센다."""
    counts: dict[str, float] = {}
    if not text:
        return counts
    for tok in text.strip().strip("-").split():
        m = _TOKEN.fullmatch(tok)
        if not m:
            raise ValueError(f"COD 화학식 토큰을 못 읽음: {tok!r} ({text!r})")
        el = "H" if m.group(1) in ("D", "T") else m.group(1)
        counts[el] = counts.get(el, 0.0) + (float(m.group(2)) if m.group(2) else 1.0)
    return counts


def reduced_formula(counts: dict[str, float]) -> str | None:
    if not counts:
        return None
    try:
        with warnings.catch_warnings():  # 비활성 기체의 전기음성도 경고
            warnings.simplefilter("ignore")
            return Composition(counts).reduced_formula
    except Exception:  # 원소가 아닌 기호 등
        return None


def _sg_plain(symbol: str) -> str:
    """번호와 상관없는 꼬리표를 뗀다: 원점 이동 '(a,b+1/2,c)', 설정 ':1'·':H', 'S'·'S1'·'Z1'·'RS',
    기호 넷(격자+세 방향) 뒤에 붙은 원점 선택 숫자 ('I 41/a m d 1')."""
    s = re.sub(r"\s*\([^)]*\)\s*$", "", symbol.strip())
    s = re.sub(r"\s*:\s*\S+$", "", s)
    s = re.sub(r"\s+(S\d?|Z\d?|RS)$", "", s)
    parts = s.split()
    return " ".join(parts[:4]) if len(parts) > 4 and parts[-1].isdigit() else s


@cache
def _sg_from_symbol(symbol: str) -> int | None:
    from pymatgen.symmetry.groups import SpaceGroup

    for s in dict.fromkeys((symbol, _sg_plain(symbol))):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                return SpaceGroup(s).int_number
        except Exception:
            continue
    return None


def _sg_number(attrs: dict) -> int | None:
    """공간군 번호. 비어 있으면(대부분) H-M 기호에서 되짚는다 — 못 읽으면 None."""
    if attrs.get("_cod_sgnumber") is not None:
        return int(attrs["_cod_sgnumber"])
    symbol = (attrs.get("_cod_sg") or "").strip()
    return _sg_from_symbol(symbol) if symbol else None


def _int(v) -> int | None:
    try:
        return int(v) if v not in (None, "") else None
    except (TypeError, ValueError):
        return None


def _row(entry: dict) -> dict:
    a = entry["attributes"]
    try:
        counts = parse_formula(a.get("_cod_calcformula"))
    except ValueError:
        counts = {}
    elements = sorted(counts) if counts else sorted(a.get("elements") or [])
    return {
        "cod_id": int(entry["id"]),
        "mineral": (a.get("_cod_mineral") or "").strip() or None,
        "commonname": (a.get("_cod_commonname") or "").strip() or None,
        "formula_cod": a.get("_cod_calcformula"),
        "formula_reduced": reduced_formula(counts),
        "elements": ",".join(elements) or None,
        "sg": a.get("_cod_sg"),
        "sg_number": _sg_number(a),
        "a": a.get("_cod_a"), "b": a.get("_cod_b"), "c": a.get("_cod_c"),
        "alpha": a.get("_cod_alpha"), "beta": a.get("_cod_beta"), "gamma": a.get("_cod_gamma"),
        "volume": a.get("_cod_vol"),
        "cell_temp": a.get("_cod_celltemp"),  # K
        "cell_pressure": a.get("_cod_cellpressure"),  # kPa
        "year": _int(a.get("_cod_year")),
        "doi": a.get("_cod_doi"),
        "r_obs": a.get("_cod_robs"),
        "r_all": a.get("_cod_rall"),
        "has_coordinates": "has coordinates" in (a.get("_cod_flags") or ""),
        "duplicate_of": _int(a.get("_cod_duplicateof")),
        "cod_status": a.get("_cod_status"),  # errors / warnings / retracted …
        "svn_revision": _int(a.get("_cod_svnrevision")),
    }


def load() -> dict[str, dict[str, int]]:
    snaps = _snapshots()
    d = snaps[-1] if snaps else fetch()
    entries: dict[str, dict] = {}
    files = [f for prefix in QUERIES for f in sorted(d.glob(f"{prefix}-*.json"))]  # 옛 스냅샷에는 common- 이 없다
    for f in files:
        for e in json.loads(f.read_text(encoding="utf-8"))["data"]:
            entries[e["id"]] = e  # 조회 중 DB 가 갱신돼 페이지 경계에서 겹친 항목은 하나로
    df = pd.DataFrame([_row(e) for e in entries.values()]).sort_values("cod_id")
    q = d / "query.json"
    retrieved = json.loads(q.read_text(encoding="utf-8"))["retrieved_at"] if q.exists() else None
    df = with_provenance(df, source=SOURCE, version=d.name, license=LICENSE, method="experimental",
                         id_column="cod_id", retrieved_at=retrieved)
    out = {"structures_exp": write_table("structures_exp", df)}
    from msl.connectors import ima_cnmnc  # 광물 목록이 이미 있으면 매칭(COD·MP)도 새로 만든다

    return out | ima_cnmnc.refresh_matches()
