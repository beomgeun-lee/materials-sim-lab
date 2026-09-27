"""시험(assay) 공통 틀 — 계획서 3.2절 플러그인 인터페이스의 v0 구현."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable

from msl.registry.load import Registry
from msl.resolve import Resolved
from msl.schema.recipe import Recipe
from msl.schema.result import Fidelity, ResultValue, Status, ValueKind


@dataclass
class Context:
    recipe: Recipe
    comps: list[Resolved]
    registry: Registry
    shared: dict[str, Any] = field(default_factory=dict)  # 시험 사이에 넘기는 계산 결과 (상태도 등)

    @property
    def T(self) -> float | None:
        t = self.recipe.conditions.T
        return t.to_si() if t else None

    @property
    def P(self) -> float:
        p = self.recipe.conditions.P
        return p.to_si() if p else 101_325.0


@dataclass
class Outcome:
    """시험 함수가 돌려주는 값. 실행기가 프로비넌스·출처 라이선스를 붙여 AssayResult 로 만든다."""

    status: Status
    fidelity: Fidelity
    engine: str  # engines.yaml id
    engine_version: str
    conditions_basis: str
    values: list[ResultValue] = field(default_factory=list)
    sources: list[tuple[str, str | None]] = field(default_factory=list)  # (sources.yaml id, 버전)
    energy_reference: str | None = None
    model_weights: str | None = None
    caveats: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    data: dict[str, Any] = field(default_factory=dict)
    summary: str = ""  # 한 줄 요약 (리포트 카드 첫 줄)


AssayFn = Callable[[Context], Outcome]


def val(name: str, value: float | str | bool, unit: str | None = None, *,
        kind: ValueKind = ValueKind.PROPERTY, unc: float | None = None) -> ResultValue:
    if isinstance(value, float):
        value = round(value, 6)
    return ResultValue(name=name, kind=kind, value=value, unit=unit, uncertainty=unc)


def not_applicable(reason: str, engine: str, fidelity: Fidelity = Fidelity.L0) -> Outcome:
    return Outcome(status=Status.NOT_APPLICABLE, fidelity=fidelity, engine=engine, engine_version="—",
                   conditions_basis="—", summary=reason)


def pending(reason: str, engine: str, fidelity: Fidelity = Fidelity.L0) -> Outcome:
    return Outcome(status=Status.PENDING, fidelity=fidelity, engine=engine, engine_version="—",
                   conditions_basis="—", summary=reason)
