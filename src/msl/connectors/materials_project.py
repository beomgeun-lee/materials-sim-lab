"""Materials Project (CC-BY-4.0) — summary 전량을 구조 없이 받아 phases_calc 로 적재 (계획서 4.4절).

API 를 한 번에 한 요청씩(1,000건 단위) 훑는다. deprecated=False 가 API 기본값이라 True 도 따로 받는다.
응답 문서는 받은 그대로 summary.jsonl.gz 에 한 줄씩 쓴다.

GNoME 유래 배치(gnome_r2scan_statics, CC BY-NC 4.0)는 이용 동의가 없는 키에는 서버가 돌려주지 않는다.
그래도 행마다 builder_meta.license·batch_id 를 보고, 섞이면 CC-BY-NC-4.0 을 달아 restricted 로 보낸다.

thermo 값(formation_energy_per_atom, energy_above_hull)은 summary 기본인 GGA_GGA+U_R2SCAN 혼합 체계다
(thermo 엔드포인트와 대조해 확인). A1 이 쓰는 GGA_GGA+U(MP2020) hull 과 섞으면 안 된다.
"""

from __future__ import annotations

import gzip
import json
import time
from pathlib import Path
from typing import Any

import pandas as pd

from msl.db import now, raw_dir, record_checksums, with_provenance, write_table
from msl.engines.mp import DB_VERSION, _rester

SOURCE = "materials-project"
TABLES = ["phases_calc"]
VERSION = DB_VERSION
LICENSE = "CC-BY-4.0"
NC_LICENSE = "CC-BY-NC-4.0"
GNOME_BATCH = "gnome_r2scan_statics"
METHOD = "DFT-PBE+U"
CORRECTION = "MP GGA_GGA+U_R2SCAN 혼합 (summary 기본 thermo) — A1 의 GGA_GGA+U(MP2020) 와 다른 기준"
FIELDS = [
    "material_id", "formula_pretty", "chemsys", "nelements", "nsites", "symmetry", "density", "volume",
    "formation_energy_per_atom", "energy_above_hull", "energy_per_atom", "uncorrected_energy_per_atom",
    "is_stable", "band_gap", "is_metal", "theoretical", "database_IDs", "deprecated", "builder_meta",
]
PAGE = 1000
RAW_FILE = "summary.jsonl.gz"
MANIFEST = "manifest.json"


def _get(session: Any, url: str, params: dict[str, Any], retries: int = 4) -> dict[str, Any]:
    for attempt in range(retries):
        try:
            resp = session.get(url, params=params, timeout=120)
            resp.raise_for_status()
            return resp.json()
        except Exception:
            if attempt == retries - 1:
                raise
            time.sleep(5 * (attempt + 1))
    raise RuntimeError("unreachable")


def fetch() -> Path:
    d = raw_dir(SOURCE, VERSION)
    t0 = time.monotonic()
    with _rester() as mpr:
        api = mpr.materials.summary
        live = str(api.db_version)
        if f"v{live}" != VERSION:
            raise RuntimeError(f"MP 라이브 DB {live} 가 스냅샷 버전 {VERSION} 과 다름 — DB_VERSION 을 먼저 올릴 것")
        url, session = api.endpoint, api.session
        gnome_visible = _get(session, url, {"batch_id": GNOME_BATCH, "_fields": "material_id", "_limit": 1})["meta"]["total_doc"]
        counts: dict[str, int] = {}
        meta: dict[str, Any] = {}
        tmp = d / (RAW_FILE + ".part")
        with gzip.open(tmp, "wt", encoding="utf-8") as fh:
            for deprecated in (False, True):
                skip, total = 0, None
                while total is None or skip < total:
                    page = _get(session, url, {"_fields": ",".join(FIELDS), "_limit": PAGE, "_skip": skip,
                                               "deprecated": deprecated})
                    meta = page["meta"]
                    total = meta["total_doc"]
                    if not page["data"]:
                        break
                    for doc in page["data"]:
                        fh.write(json.dumps(doc, ensure_ascii=False, default=str) + "\n")
                    skip += len(page["data"])
                counts[f"deprecated={deprecated}"] = total or 0
        tmp.rename(d / RAW_FILE)
    manifest = {
        "db_version": live, "api_version": meta.get("api_version"), "retrieved_at": now(),
        "endpoint": url, "fields": FIELDS, "page_size": PAGE, "total_doc": counts,
        "access_controlled_batch_ids": list(api.access_controlled_batch_ids),
        "gnome_docs_visible_to_key": gnome_visible, "elapsed_s": round(time.monotonic() - t0, 1),
    }
    (d / MANIFEST).write_text(json.dumps(manifest, ensure_ascii=False, indent=1), encoding="utf-8")
    record_checksums(d)
    return d


def _row(doc: dict[str, Any], controlled: set[str]) -> dict[str, Any]:
    sym = doc.get("symmetry") or {}
    bm = doc.get("builder_meta") or {}
    icsd = (doc.get("database_IDs") or {}).get("icsd") or []
    nc = bm.get("license") == "BY-NC" or bm.get("batch_id") in controlled
    return {
        "material_id": doc["material_id"], "formula_pretty": doc.get("formula_pretty"), "chemsys": doc.get("chemsys"),
        "nelements": doc.get("nelements"), "nsites": doc.get("nsites"),
        "spacegroup_symbol": sym.get("symbol"), "spacegroup_number": sym.get("number"),
        "crystal_system": sym.get("crystal_system"),
        "density_g_cm3": doc.get("density"), "volume_A3": doc.get("volume"),
        "formation_energy_per_atom_eV": doc.get("formation_energy_per_atom"),
        "energy_above_hull_eV": doc.get("energy_above_hull"),
        "energy_per_atom_eV": doc.get("energy_per_atom"),
        "uncorrected_energy_per_atom_eV": doc.get("uncorrected_energy_per_atom"),
        "is_stable": doc.get("is_stable"), "band_gap_eV": doc.get("band_gap"), "is_metal": doc.get("is_metal"),
        "theoretical": doc.get("theoretical"), "has_icsd": bool(icsd), "icsd_ids": ",".join(map(str, icsd)) or None,
        "deprecated": doc.get("deprecated"), "batch_id": bm.get("batch_id"), "mp_license": bm.get("license"),
        "license": NC_LICENSE if nc else LICENSE,
    }


def load() -> dict[str, dict[str, int]]:
    d = raw_dir(SOURCE, VERSION)
    if not (d / RAW_FILE).exists():
        fetch()
    manifest = json.loads((d / MANIFEST).read_text(encoding="utf-8"))
    controlled = set(manifest.get("access_controlled_batch_ids") or [GNOME_BATCH])
    with gzip.open(d / RAW_FILE, "rt", encoding="utf-8") as fh:
        rows = [_row(json.loads(line), controlled) for line in fh]
    df = pd.DataFrame(rows).drop_duplicates("material_id")
    lic = df.pop("license")
    df = with_provenance(df, source=SOURCE, version=VERSION, license=LICENSE, method=METHOD, id_column="material_id",
                         retrieved_at=manifest["retrieved_at"], correction_scheme=CORRECTION)
    df["license"] = lic.values
    out = {"phases_calc": write_table("phases_calc", df)}
    from msl.connectors import ima_cnmnc  # 광물 목록이 이미 있으면 광물 ↔ MP 매칭도 새로 만든다

    if ima_cnmnc.minerals_loaded():
        out["mineral_mp"] = ima_cnmnc.match_mp()
    return out
