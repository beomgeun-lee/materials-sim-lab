"""Materials Project 어댑터 — 화학계 엔트리와 요약 물성 (mp-api, 결과는 디스크 캐시)."""

from __future__ import annotations

import hashlib
import json
import logging
import pickle
import os
import warnings
from importlib.metadata import version

from pymatgen.entries.computed_entries import ComputedStructureEntry

from msl.env import CACHE_DIR, load_dotenv

CACHE = CACHE_DIR / "mp"
DB_VERSION = "v2026.04.13"
THERMO_TYPE = "GGA_GGA+U"  # MP2020 보정 체계 — 다른 thermo type 과 한 hull 에 섞지 않는다
ENERGY_REFERENCE = "MP2020 (GGA/GGA+U)"


class MPUnavailable(RuntimeError):
    """API 키가 없거나 MP 에 접근할 수 없을 때."""


def engine_version() -> str:
    return version("mp-api")


def _rester():
    load_dotenv()
    if not os.environ.get("MP_API_KEY"):
        raise MPUnavailable("MP_API_KEY 가 .env 에 없음")
    logging.getLogger("mp_api").setLevel(logging.ERROR)
    warnings.filterwarnings("ignore", module="mp_api")
    from mp_api.client import MPRester

    return MPRester(mute_progress_bars=True)


def mp_int(mid: str | None) -> int | None:
    """MP ID 를 정수로. 새 형식(mp-aaaaafwb, AlphaID)과 옛 형식(mp-3953-GGA+U)이 같은 정수를 가리킨다."""
    if not mid:
        return None
    from emmet.core.mpid import AlphaID

    tail = str(mid).split("-GGA")[0].split("-r2SCAN")[0].split("-")[-1]
    try:
        return int(tail) if tail.isdigit() else int(AlphaID(tail))
    except Exception:
        return None


def entries_in_chemsys(elements: list[str]) -> list[ComputedStructureEntry]:
    """화학계의 모든 엔트리(MP2020 보정 적용). 같은 계는 캐시에서 읽는다."""
    chemsys = "-".join(sorted(set(elements)))
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"entries_{chemsys}_{THERMO_TYPE}.pkl"  # 엔트리 data 에 Element 키가 있어 JSON 대신 pickle (로컬 캐시 전용)
    if path.exists():
        return pickle.loads(path.read_bytes())
    with _rester() as mpr:
        entries = mpr.get_entries_in_chemsys(
            sorted(set(elements)), additional_criteria={"thermo_types": [THERMO_TYPE]}
        )
    path.write_bytes(pickle.dumps(entries))
    return entries


def summaries(material_ids: list[str], fields: list[str]) -> dict[str, dict]:
    """material_id → 요약 물성. 같은 요청은 캐시에서 읽는다."""
    if not material_ids:
        return {}
    CACHE.mkdir(parents=True, exist_ok=True)
    key = "-".join(sorted(material_ids))[:180] + "_" + "-".join(sorted(fields))
    path = CACHE / f"summary_{hashlib.sha256(key.encode()).hexdigest()[:16]}.json"
    if path.exists():
        cached = json.loads(path.read_text(encoding="utf-8"))
        if cached.get("key") == key:
            return cached["docs"]
    with _rester() as mpr:
        docs = mpr.materials.summary.search(material_ids=material_ids, fields=["material_id", *fields])
    out = {str(d.material_id): {f: getattr(d, f, None) for f in fields} for d in docs}
    path.write_text(json.dumps({"key": key, "docs": out}, default=str), encoding="utf-8")
    return out
