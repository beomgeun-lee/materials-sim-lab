"""CAMEO Chemicals 반응성 그룹 — PubChem PUG-View 'Reactive Group' 주석 경유 (CID ↔ 그룹).

CAMEO 사이트는 robots.txt 로 수집을 막고 벌크를 제공하지 않으므로 PubChem 이 재배포하는 주석을 쓴다.
라이선스는 미확인(unknown) → restricted 파티션 (계획서 10절 #1, NOAA 문의 대상).
"""

from __future__ import annotations

import datetime as dt
import json
import time
from pathlib import Path

import pandas as pd

from msl.db import download, raw_dir, record_checksums, with_provenance, write_table

SOURCE = "cameo"
TABLES = ["reactive_groups"]
URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/annotations/heading/Reactive%20Group/JSON?page={page}"
LICENSE = "unknown"


def _snapshots() -> list[Path]:
    return sorted(p for p in raw_dir(SOURCE, "_").parent.glob("*") if p.is_dir() and list(p.glob("page_*.json")))


def fetch() -> Path:
    d = raw_dir(SOURCE, dt.date.today().isoformat())
    first = download(URL.format(page=1), d / "page_001.json")
    total = json.loads(first.read_text(encoding="utf-8"))["Annotations"].get("TotalPages", 1)
    for page in range(2, total + 1):
        time.sleep(0.4)  # PubChem 권장 한도 안쪽
        download(URL.format(page=page), d / f"page_{page:03d}.json")
    record_checksums(d)
    return d


def load() -> dict[str, dict[str, int]]:
    snaps = _snapshots()
    d = snaps[-1] if snaps else fetch()
    rows = []
    for f in sorted(d.glob("page_*.json")):
        for ann in json.loads(f.read_text(encoding="utf-8"))["Annotations"]["Annotation"]:
            cids = (ann.get("LinkedRecords") or {}).get("CID") or [None]
            groups = sorted({s["String"] for item in ann.get("Data", [])
                             for s in item.get("Value", {}).get("StringWithMarkup", []) if s.get("String")})
            for cid in cids:
                for g in groups:
                    rows.append({"cid": cid, "name": ann.get("Name"), "reactive_group": g,
                                 "cameo_id": ann.get("SourceID")})
    df = pd.DataFrame(rows).dropna(subset=["cid"])
    df["cid"] = df["cid"].astype(int)
    df["key"] = df["cid"].astype(str) + ":" + df["reactive_group"]
    df = with_provenance(df, source=SOURCE, version=d.name, license=LICENSE,
                         method="compiled", id_column="key")
    return {"reactive_groups": write_table("reactive_groups", df)}
