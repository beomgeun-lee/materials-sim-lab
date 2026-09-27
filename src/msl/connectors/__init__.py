"""데이터 소스 커넥터 — 소스마다 모듈 하나 (계획서 1단계).

각 모듈은 다음을 제공한다.
    SOURCE: str            kb/sources.yaml 의 id
    TABLES: list[str]      적재하는 테이블
    fetch() -> Path        원본 스냅샷을 data/raw/{SOURCE}/{version}/ 에 받고 SHA256SUMS 기록
    load() -> dict         스냅샷을 정규화해 msl.db.write_table 로 적재. {테이블: {파티션: 행 수}}

모듈 이름은 소스 id 의 하이픈을 밑줄로 바꾼 것이다 (예: nasa-cea-thermo → nasa_cea_thermo).
"""

from __future__ import annotations

import importlib
import pkgutil
from types import ModuleType


def available() -> dict[str, ModuleType]:
    """connectors 패키지의 모든 커넥터 모듈 (SOURCE 가 있는 것)."""
    out: dict[str, ModuleType] = {}
    for info in pkgutil.iter_modules(__path__):
        if info.name.startswith("_"):
            continue
        mod = importlib.import_module(f"{__name__}.{info.name}")
        if hasattr(mod, "SOURCE"):
            out[mod.SOURCE] = mod
    return dict(sorted(out.items()))
