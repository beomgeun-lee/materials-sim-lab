"""COST 507 경합금 열역학 DB (Round II, 1999) — 29원소, 평가된 2원계·3원계 (Al·Mg·Cu·Zn·Si 등).

OpenCalphad 배포본(원본)은 문법 오류가 있어 pycalphad 가 못 읽는다 (Y 원자량 '8.89059+01', 'REF: 0',
중복 FUNCTION R 등). pycalphad 팀이 고쳐 테스트 DB 로 동봉한 사본을 계산용으로 쓰고, 원본 zip 도 함께 둔다.
파일 안 문구는 "It may be used freely but on your own risk" 뿐이라 상업 이용·재배포 명시 없음 → unknown (restricted).
"""

from __future__ import annotations

import shutil
import zipfile
from importlib.metadata import version
from pathlib import Path

from msl.connectors._tdb import index_file, write_rows
from msl.db import download, raw_dir, record_checksums

SOURCE = "cost507"
TABLES = ["tdb_systems"]
VERSION = f"round-ii-1999+pycalphad-{version('pycalphad')}"
LICENSE = "unknown"
URL = "https://www.opencalphad.com/databases/COST507.zip"
USED = "COST507.pycalphad.tdb"
REFERENCE = ("COST 507 Thermochemical database for light metal alloys, Round II (1999) — "
             "I. Ansara, A.T. Dinsdale, M.H. Rand (eds.), EUR 18499 (pycalphad 동봉 수정본)")


def bundled() -> Path:
    import pycalphad

    return Path(pycalphad.__file__).parent / "tests" / "databases" / "COST507.tdb"


def fetch() -> Path:
    d = raw_dir(SOURCE, VERSION)
    zpath = download(URL, d / "COST507.zip")
    with zipfile.ZipFile(zpath) as zf:
        (d / "COST507.tdb").write_bytes(zf.read("COST507.tdb"))
    shutil.copy2(bundled(), d / USED)
    record_checksums(d)
    return d


def load() -> dict[str, dict[str, int]]:
    d = raw_dir(SOURCE, VERSION)
    if not (d / USED).exists():
        fetch()
    row = index_file(d / USED, USED, quality="evaluated", dataset="COST507", year=1999, variant=None, preferred=True)
    row["reference"] = REFERENCE
    return {"tdb_systems": write_rows([row], source=SOURCE, version=VERSION, license=LICENSE)}
