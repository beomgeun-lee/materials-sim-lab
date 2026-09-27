"""프로비넌스 — 데이터 레코드의 출처와 시험 실행 이력 (계획서 4.2·3.5절)."""

from __future__ import annotations

import datetime as dt
import hashlib
import json
from typing import Any

from pydantic import BaseModel, ConfigDict, Field

# experimental / estimated / compiled / CALPHAD / DFT-<범함수> / ML-<모델>
METHOD_PATTERN = r"^(experimental|estimated|compiled|CALPHAD|DFT-[A-Za-z0-9+\-]+|ML-[A-Za-z0-9.+\-]+)$"


class DataProvenance(BaseModel):
    """DB 에 적재되는 모든 레코드가 갖는 출처 필드."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    source: str  # sources.yaml id
    source_id: str  # 원 소스에서의 식별자 (mp-…, COD ID, CID 등)
    source_version: str
    retrieved_at: dt.datetime
    license: str  # licenses.yaml id
    method: str = Field(pattern=METHOD_PATTERN)
    correction_scheme: str | None = None  # MP2020, OQMD-fit-mu 등


class RunProvenance(BaseModel):
    """시험 한 번 실행의 이력. input_hash 가 같으면 결과를 캐시에서 재사용한다."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    input_hash: str = Field(pattern=r"^[0-9a-f]{64}$")
    msl_version: str
    engine_versions: dict[str, str]
    parameters: dict[str, Any] = Field(default_factory=dict)
    started_at: dt.datetime
    finished_at: dt.datetime


def canonical_hash(obj: Any) -> str:
    """입력을 정규화한 JSON 의 SHA-256. 키 순서·공백과 무관하게 같은 입력이면 같은 해시."""
    if isinstance(obj, BaseModel):
        obj = obj.model_dump(mode="json")
    payload = json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":"), default=str)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()
