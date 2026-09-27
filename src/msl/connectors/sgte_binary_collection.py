"""SGTE 2원계 데이터셋 모음 (Hallstedt 2025, Calphad 89, 102833 부록) — A6 합금 상평형의 기본 TDB.

부록 zip(mmc3) 에 평가된 2원계 TDB 2,773개가 들어 있다. 한 계에 여러 평가가 있으면
Landolt-Börnstein 수록본(-LB)을, 없으면 가장 최근 평가를 preferred 로 표시한다.
라이선스: 논문은 CC BY 4.0 이지만 부록 파일에 그대로 적용되는지 확인하지 못해 unknown (restricted).
"""

from __future__ import annotations

import re
import shutil
import subprocess
import urllib.error
import zipfile
from pathlib import Path

from msl.connectors._tdb import index_file, write_rows
from msl.db import download, raw_dir, record_checksums

SOURCE = "sgte-binary-collection"
TABLES = ["tdb_systems"]
VERSION = "calphad-89-102833"
LICENSE = "unknown"
URL = "https://ars.els-cdn.com/content/image/1-s2.0-S0364591625000367-mmc3.zip"
REFERENCE = "B. Hallstedt, The SGTE collection of binary datasets, Calphad 89 (2025) 102833"
# 파일 이름: <원소쌍>-<연도 2자리><저자 3자>[-변형…]  예) CuNi-92Mey-LB, AlSi-24Zob-3g
_NAME = re.compile(r"^(?P<pair>(?:[A-Z][a-z]?){2})-(?P<yy>\d{2})(?P<author>[A-Za-z]+)(?:-(?P<variant>.+))?$")


def _download(dest: Path) -> Path:
    """Elsevier CDN 은 urllib 요청을 403 으로 막는다 (같은 User-Agent 의 curl 은 통과) → curl 로 한 번 더."""
    try:
        return download(URL, dest)
    except urllib.error.HTTPError as exc:
        if exc.code != 403 or shutil.which("curl") is None:
            raise
    dest.with_suffix(dest.suffix + ".part").unlink(missing_ok=True)
    subprocess.run(["curl", "-fsSL", "-A", "materials-sim-lab/0.0.1 (+research)", "-o", str(dest), URL], check=True)
    return dest


def fetch() -> Path:
    d = raw_dir(SOURCE, VERSION)
    zpath = _download(d / "mmc3.zip")
    out = d / "datasets"
    out.mkdir(exist_ok=True)
    with zipfile.ZipFile(zpath) as zf:
        for info in zf.infolist():
            if info.filename.endswith(".tdb"):
                (out / Path(info.filename).name).write_bytes(zf.read(info))
    record_checksums(d)
    return d


def parse_name(stem: str) -> dict:
    m = _NAME.match(stem)
    if not m:
        return {"dataset": stem, "year": None, "variant": None, "system_hint": ""}
    yy = int(m["yy"])
    pair = sorted(re.findall(r"[A-Z][a-z]?", m["pair"]))
    return {"dataset": stem, "year": (1900 if yy >= 30 else 2000) + yy, "variant": m["variant"],
            "system_hint": "-".join(pair)}


def mark_preferred(rows: list[dict]) -> None:
    """계마다 하나: -LB 수록본 → 변형 없는 최근 평가 → 아무 최근 평가 (pycalphad 로 읽히는 것 중)."""
    by_system: dict[str, list[dict]] = {}
    for r in rows:
        r["preferred"] = False
        if r["parse_ok"]:
            by_system.setdefault(r["system"], []).append(r)
    for group in by_system.values():
        lb = [r for r in group if (r["variant"] or "").split("-")[0] == "LB"]
        plain = [r for r in group if not r["variant"]]
        pick = lb or plain or group
        max(pick, key=lambda r: (r["year"] or 0, r["dataset"]))["preferred"] = True


def load() -> dict[str, dict[str, int]]:
    d = raw_dir(SOURCE, VERSION)
    if not (d / "datasets").exists():
        fetch()
    rows = []
    for f in sorted((d / "datasets").glob("*.tdb")):
        meta = parse_name(f.stem)
        hint = meta.pop("system_hint")
        row = index_file(f, f"datasets/{f.name}", quality="evaluated", **meta)
        row["system"] = row["system"] or hint
        rows.append(row)
    mark_preferred(rows)
    for r in rows:
        r["reference"] = f"{r['reference']} — {REFERENCE}" if r["reference"] else REFERENCE
    return {"tdb_systems": write_rows(rows, source=SOURCE, version=VERSION, license=LICENSE)}
