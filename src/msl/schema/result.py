"""시험 결과 스키마 — 계획서 3.5절 '결과 신뢰 규칙'.

모든 결과는 충실도·엔진·에너지 기준·조건 기준·출처·한계를 반드시 갖는다.
서로 다른 에너지 기준의 값이 한 판정에 섞이지 않게 하는 것이 목적이다.
"""

from __future__ import annotations

from enum import StrEnum
from typing import Any

from pydantic import BaseModel, ConfigDict, Field, model_validator

from msl.schema.provenance import RunProvenance


class Fidelity(StrEnum):
    """충실도 계층 — 계획서 3.3절. T 는 평형 트랙."""

    L0 = "L0"  # 조회·규칙
    L1 = "L1"  # ML 대리모델
    L2 = "L2"  # uMLIP
    L3 = "L3"  # DFT·상용 CALPHAD
    T = "T"  # 평형 솔버


class Status(StrEnum):
    OK = "ok"
    WARNING = "warning"  # 결과는 있으나 주의 필요 (S0 '주의' 등)
    FAILED = "failed"  # 실행 실패
    NOT_APPLICABLE = "not-applicable"  # 이 조합에는 해당 없음
    PENDING = "pending"  # 아직 실행할 수 없음 (미구현 단계·키 대기·데이터 확인 중)


class ValueKind(StrEnum):
    ENERGY = "energy"  # 형성에너지·반응에너지·E_hull 등 → energy_reference 필수
    PROPERTY = "property"
    FLAG = "flag"  # 판정 (적합/주의/부적합 등)
    TEXT = "text"


class ResultValue(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    name: str
    kind: ValueKind
    value: float | str | bool
    unit: str | None = None
    uncertainty: float | None = Field(default=None, ge=0)  # ± 값 (같은 단위)


class SourceCitation(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    source: str  # sources.yaml id
    version: str | None = None
    record_ids: list[str] = Field(default_factory=list)
    license: str  # licenses.yaml id


class AssayResult(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    assay: str = Field(pattern=r"^(S0|A[1-9][0-9]*)$")
    recipe_id: str
    status: Status
    fidelity: Fidelity
    engine: str  # engines.yaml id
    engine_version: str
    model_weights: str | None = None  # 이름@sha256:… (uMLIP 등)
    energy_reference: str | None = None  # MP2020 / OQMD / self-consistent:<모델> …
    conditions_basis: str  # "0 K DFT", "298.15 K 표준상태", "입력 조건" …
    values: list[ResultValue] = Field(default_factory=list)
    sources: list[SourceCitation] = Field(default_factory=list)
    caveats: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    data: dict[str, Any] = Field(default_factory=dict)  # 표·곡선 등 화면용 구조화 결과
    provenance: RunProvenance

    @model_validator(mode="after")
    def _trust_rules(self) -> AssayResult:
        if any(v.kind is ValueKind.ENERGY for v in self.values) and not self.energy_reference:
            raise ValueError("에너지 값이 있으면 energy_reference 가 필요함")
        if self.fidelity is Fidelity.L2 and self.model_weights is None:
            raise ValueError("L2(uMLIP) 결과는 model_weights 가 필요함")
        if self.status in (Status.OK, Status.WARNING) and not self.values:
            raise ValueError("ok/warning 결과에는 값이 하나 이상 필요함")
        return self
