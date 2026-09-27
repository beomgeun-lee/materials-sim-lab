"""MP 수용액 이온 기준 데이터 (MPContribs 'ion_ref_data', CC BY 4.0) — A5 Pourbaix 의 이온 에너지.

pymatgen/mp-api 가 Pourbaix 도표를 만들 때 쓰는 실험 이온 생성 자유에너지 362건
(NBS Technical Note 270, Pourbaix 도감 등에서 모은 값, Persson 2012 체계).
mp-api 는 이 표를 MPContribs 클라이언트로 받는데, 그 클라이언트의 의존성(pint>=0.25)이
mendeleev(pint<0.25)와 충돌해 공개 REST 로 직접 받아 스냅샷으로 둔다.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from msl.db import raw_dir, record_checksums, with_provenance, write_table

SOURCE = "mp-ion-ref-data"
TABLES = ["aqueous_ions"]
VERSION = "2022-07-28"  # 기여 항목의 마지막 수정일 (fetch 가 확인)
LICENSE = "CC-BY-4.0"
API = "https://contribs-api.materialsproject.org"
FILE = "ion_ref_data.json"


def _get(path: str, params: dict[str, str | int]) -> dict:
    import requests  # urllib 기본 User-Agent 는 403 으로 막힌다

    r = requests.get(f"{API}{path}", params=params, timeout=60,
                     headers={"User-Agent": "materials-sim-lab/0.0.1 (+research)", "Accept": "application/json"})
    r.raise_for_status()
    return r.json()


def fetch() -> Path:
    d = raw_dir(SOURCE, VERSION)
    project = _get("/projects/ion_ref_data", {"_fields": "name,title,license,is_public,references,description,columns"})
    fields = ",".join(["identifier", "formula", "last_modified", *(c["path"] for c in project["columns"])])
    rows: list[dict] = []
    for _ in range(100):  # _page 는 무시되고 _skip 만 동작한다
        chunk = _get("/contributions/", {"project": "ion_ref_data", "_fields": fields, "_limit": 100, "_skip": len(rows)})
        rows += chunk["data"]
        if not chunk.get("has_more") or not chunk["data"]:
            break
    else:
        raise RuntimeError("ion_ref_data 페이지가 끝나지 않음")
    if len({r["identifier"] for r in rows}) != len(rows):
        raise RuntimeError("ion_ref_data 중복 항목 — 페이지 이동이 바뀌었는지 확인")
    latest = max(r["last_modified"] for r in rows)[:10]
    if latest != VERSION:
        raise RuntimeError(f"ion_ref_data 가 갱신됨 ({latest}) — VERSION 을 올리고 다시 받을 것")
    rows.sort(key=lambda r: r["identifier"])
    (d / "project.json").write_text(json.dumps(project, ensure_ascii=False, indent=1), encoding="utf-8")
    (d / FILE).write_text(json.dumps(rows, ensure_ascii=False, indent=1), encoding="utf-8")
    record_checksums(d)
    return d


def records() -> list[dict]:
    """mp-api get_ion_reference_data() 와 같은 형식의 목록 (identifier, formula, data{...})."""
    src = raw_dir(SOURCE, VERSION) / FILE
    if not src.exists():
        fetch()
    return json.loads(src.read_text(encoding="utf-8"))


def load() -> dict[str, dict[str, int]]:
    rows = [{
        "identifier": r["identifier"],
        "formula": r["formula"],
        "charge": r["data"]["charge"]["value"],
        "dgf_kj_mol": r["data"]["ΔGᶠ"]["value"],
        "major_element": r["data"]["MajElements"],
        "ref_solid": r["data"]["RefSolid"],
        "ref_solid_dgf_kj_mol": r["data"]["ΔGᶠRefSolid"]["value"],
        "reference": r["data"]["reference"],
    } for r in records()]
    df = with_provenance(pd.DataFrame(rows), source=SOURCE, version=VERSION, license=LICENSE,
                         method="experimental", id_column="identifier")
    return {"aqueous_ions": write_table("aqueous_ions", df)}
