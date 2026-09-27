"""CALPHAD TDB 커넥터 공통 — TDB 파일 하나를 tdb_systems 한 행으로 색인한다 (A6 합금 상평형).

행 = 파일 하나. system 은 파일이 다루는 원소(알파벳순, '-' 연결), binaries 는 G·L 파라미터에
두 원소가 함께 나오는 쌍(= 평가된 2원계)이다. A6 은 레시피 원소의 모든 쌍이 binaries 에 있는
파일만 쓴다 (없는 쌍을 이상용액으로 계산하면 그럴듯한 오답이 나온다).
"""

from __future__ import annotations

import re
import warnings
from functools import lru_cache
from itertools import combinations
from pathlib import Path
from typing import Any

import pandas as pd

from msl.db import with_provenance, write_table
from msl.schema.refs import ELEMENTS

TABLE = "tdb_systems"
_ELEMENTS = {e.upper(): e for e in ELEMENTS}
_ENERGY_PARAMS = {"G", "L"}


def read_text(path: Path) -> str:
    """TDB 는 인코딩 표기가 없다 — UTF-8 로 읽고 안 되면 latin-1 (MatCalc·SGTE 파일의 ö·é 등)."""
    raw = path.read_bytes()
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("latin-1")


@lru_cache(maxsize=16)
def read_database(path: str) -> Any:
    from pycalphad import Database

    with warnings.catch_warnings():  # TYPE_DEFINITION 누락 같은 무해한 경고가 수천 줄 쌓이지 않게
        warnings.simplefilter("ignore", UserWarning)
        return Database.from_string(read_text(Path(path)), fmt="tdb")


def elements_of(dbf: Any) -> list[str]:
    """실제 원소만 (VA·/-·ZE 같은 가상 성분 제외), 원소 기호 표기."""
    return sorted(_ELEMENTS[str(e).upper()] for e in dbf.elements if str(e).upper() in _ELEMENTS)


def binaries(dbf: Any) -> list[str]:
    """G·L 파라미터에 함께 나오는 원소 쌍 ('Cu-Ni' 형식, 정렬)."""
    pairs: set[tuple[str, str]] = set()
    for p in dbf._parameters.all():
        if p["parameter_type"] not in _ENERGY_PARAMS:
            continue
        els = {_ELEMENTS[str(el).upper()] for subl in p["constituent_array"] for sp in subl
               for el in getattr(sp, "constituents", {}) if str(el).upper() in _ELEMENTS}
        pairs.update(combinations(sorted(els), 2))
    return sorted(f"{a}-{b}" for a, b in pairs)


def header_reference(text: str, lines: int = 3) -> str:
    """파일 첫 주석 블록에서 의미 있는 줄 몇 개 (출처·논문 표기)."""
    out = []
    for raw in text.splitlines():
        s = raw.strip()
        if not s:
            continue
        if not s.startswith("$"):
            break
        s = s.lstrip("$").strip()
        if s and not re.fullmatch(r"[-=*#$ ]*", s):
            out.append(s)
        if len(out) == lines:
            break
    return " / ".join(out)[:300]


def index_file(path: Path, rel: str, **extra: Any) -> dict[str, Any]:
    """TDB 한 파일 → tdb_systems 행. pycalphad 가 못 읽으면 parse_ok=False 로 남긴다 (A6 은 건너뜀)."""
    text = read_text(path)
    row: dict[str, Any] = {"file": rel, "reference": header_reference(text), **extra}
    try:
        dbf = read_database(str(path))
    except Exception as exc:  # 비표준 명령(C_S·STATUS 등) — 행은 남기고 이유를 적는다
        row.update(parse_ok=False, parse_error=f"{type(exc).__name__}: {str(exc).splitlines()[0][:160]}",
                   phases="", binaries="")
        row.setdefault("system", "")
        return row
    els = elements_of(dbf)
    row.update(system="-".join(els), phases=",".join(sorted(dbf.phases)),
               binaries=",".join(binaries(dbf)), parse_ok=True, parse_error=None)
    read_database.cache_clear()  # 적재 중 수천 개를 열어도 메모리가 쌓이지 않게
    return row


def write_rows(rows: list[dict[str, Any]], *, source: str, version: str, license: str) -> dict[str, int]:
    df = pd.DataFrame(rows)
    df["n_elements"] = df["system"].map(lambda s: len(s.split("-")) if s else 0)
    df = with_provenance(df, source=source, version=version, license=license, method="CALPHAD", id_column="file")
    return write_table(TABLE, df)
