"""PHREEQC 동봉 열역학 DB (USGS 공개) — 수용액 평형(A4)이 쓰는 광물·기체 상과 마스터 화학종.

phreeqc.dat·pitzer.dat 는 phreeqpython 동봉본을, llnl.dat·minteq.v4.dat·sit.dat 는 USGS phreeqc3
저장소에서 받는다. Thermoddem 파일은 BRGM 약관 문제로 받지 않는다 (sources.yaml license_note).
"""

from __future__ import annotations

import re
import shutil
from importlib.metadata import version
from pathlib import Path

import pandas as pd

from msl.db import download, raw_dir, record_checksums, with_provenance, write_table

SOURCE = "phreeqc-db"
TABLES = ["aqueous_phases", "aqueous_masters"]
VERSION = f"phreeqpython-{version('phreeqpython')}+phreeqc3-master"
LICENSE = "USGS-public"
REMOTE = "https://raw.githubusercontent.com/usgs-coupled/phreeqc3/master/database/{name}"
BUNDLED = ("phreeqc.dat", "pitzer.dat")
FETCHED = ("llnl.dat", "minteq.v4.dat", "sit.dat")


def fetch() -> Path:
    import phreeqpython

    d = raw_dir(SOURCE, VERSION)
    for name in BUNDLED:
        shutil.copy2(Path(phreeqpython.__file__).parent / "database" / name, d / name)
    for name in FETCHED:
        download(REMOTE.format(name=name), d / name)
    record_checksums(d)
    return d


def _block(text: str, keyword: str) -> str:
    """KEYWORD 블록 본문 (다음 대문자 키워드 줄 전까지)."""
    m = re.search(rf"^{keyword}\s*$(.*?)(?=^[A-Z_]{{4,}}\s*$)", text, re.M | re.S)
    return m.group(1) if m else ""


def _num(tokens: list[str]) -> float | None:
    try:
        return float(tokens[0]) if tokens else None
    except ValueError:
        return None


def parse_phases(text: str) -> list[dict]:
    rows, current = [], None
    for raw in _block(text, "PHASES").splitlines():
        line = raw.split("#", 1)[0].split(";", 1)[0].rstrip()  # ';' 뒤는 같은 줄의 다음 문장
        if not line.strip():
            continue
        if not raw[:1].isspace():  # 상 이름 줄
            current = {"phase": line.split()[0], "reaction": None, "log_k": None, "delta_h": None}
            rows.append(current)
        elif current is not None:
            s = line.strip()
            if current["reaction"] is None and "=" in s:
                current["reaction"] = s
            elif re.match(r"^-?log_?k\b", s, re.I):
                current["log_k"] = _num(s.split()[1:2])
            elif re.match(r"^-?delta_?h\b", s, re.I):
                current["delta_h"] = _num(s.split()[1:2])
    for r in rows:
        r["formula"] = r["reaction"].split("=")[0].strip().split(" + ")[0] if r["reaction"] else None
        r["is_gas"] = "(g)" in r["phase"]
    return [r for r in rows if r["reaction"]]


def parse_masters(text: str) -> list[dict]:
    rows = []
    for raw in _block(text, "SOLUTION_MASTER_SPECIES").splitlines():
        parts = raw.split("#", 1)[0].split()
        if len(parts) >= 2 and re.match(r"^[A-Z][a-z]?(\(-?\d+\))?$", parts[0]):
            rows.append({"master": parts[0], "species": parts[1]})
    return rows


def load() -> dict[str, dict[str, int]]:
    d = raw_dir(SOURCE, VERSION)
    if not all((d / n).exists() for n in BUNDLED + FETCHED):
        fetch()
    phases, masters = [], []
    for name in BUNDLED + FETCHED:
        text = (d / name).read_text(encoding="latin-1")
        phases += [{"database": name, **r} for r in parse_phases(text)]
        masters += [{"database": name, **r} for r in parse_masters(text)]
    ph = pd.DataFrame(phases)
    ph["key"] = ph["database"] + ":" + ph["phase"]
    ms = pd.DataFrame(masters)
    ms["key"] = ms["database"] + ":" + ms["master"]
    kw = dict(source=SOURCE, version=VERSION, license=LICENSE, method="compiled", id_column="key")
    return {"aqueous_phases": write_table("aqueous_phases", with_provenance(ph, **kw)),
            "aqueous_masters": write_table("aqueous_masters", with_provenance(ms, **kw))}
