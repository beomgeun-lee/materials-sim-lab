"""조합 레시피 스키마 v0 — 계획서 3.1절.

레시피 = 성분 목록(무엇을 얼마나, 어떤 상태로) + 조건(T, P, 분위기, pH, 시간) + 적용할 시험.
두 가지 모드가 있다.
- mixture : 구체적인 양을 섞는다 (석회석 1 mol + 석영 1 mol)
- system  : 원소 계(chemical system)를 탐색한다 (Li–Co–O 에 안정한 화합물이 있는가)
"""

from __future__ import annotations

import re
from enum import StrEnum
from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, ValidationError, field_validator, model_validator

from msl.schema.quantity import FRACTION_DIMENSIONS, Dimension, Quantity
from msl.schema.refs import Namespace, SubstanceRef

ASSAY_CODE_RE = re.compile(r"^(S0|A[1-9][0-9]*)$")
AMOUNT_DIMENSIONS = frozenset(
    {Dimension.AMOUNT, Dimension.MASS, Dimension.VOLUME, Dimension.CONCENTRATION, Dimension.MOLALITY}
    | FRACTION_DIMENSIONS
)
FRACTION_SUM_TOLERANCE = 0.5  # %


class RecipeError(ValueError):
    """레시피 파일을 읽거나 검증하다 실패했을 때."""


class State(StrEnum):
    SOLID = "solid"
    LIQUID = "liquid"
    GAS = "gas"
    AQUEOUS = "aqueous"


class Mode(StrEnum):
    MIXTURE = "mixture"
    SYSTEM = "system"


class _Model(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)


def _require(q: Quantity | None, dims: frozenset[Dimension] | Dimension, what: str) -> None:
    if q is None:
        return
    allowed = {dims} if isinstance(dims, Dimension) else dims
    if q.dimension not in allowed:
        raise ValueError(f"{what}에 쓸 수 없는 단위: {q} (차원 {q.dimension.value})")


class Component(_Model):
    ref: SubstanceRef
    amount: Quantity | None = None
    state: State | None = None  # 비우면 해석기가 표준상태로 추정한다 (2단계)
    label: str | None = None

    @model_validator(mode="after")
    def _check_amount(self) -> Component:
        _require(self.amount, AMOUNT_DIMENSIONS, "성분 양")
        if self.amount is not None and self.amount.value <= 0:
            raise ValueError(f"성분 양은 0보다 커야 함: {self.amount}")
        return self


class Conditions(_Model):
    T: Quantity | None = None
    P: Quantity | None = None
    atmosphere: str | None = None
    pH: float | None = Field(default=None, ge=-2, le=16)
    time: Quantity | None = None

    @model_validator(mode="after")
    def _check(self) -> Conditions:
        _require(self.T, Dimension.TEMPERATURE, "온도")
        _require(self.P, Dimension.PRESSURE, "압력")
        _require(self.time, Dimension.TIME, "시간")
        if self.T is not None and self.T.to_si() <= 0:
            raise ValueError(f"온도는 0 K 보다 높아야 함: {self.T}")
        if self.P is not None and self.P.to_si() <= 0:
            raise ValueError(f"압력은 0 보다 커야 함: {self.P}")
        return self


class Sweep(_Model):
    """한 변수를 구간에서 훑는다. ratio 는 두 성분의 혼합비 x (첫 성분 비율, 0~1)."""

    param: Literal["T", "P", "ratio"]
    start: Quantity | float
    stop: Quantity | float
    steps: int = Field(ge=2, le=10_000)

    @model_validator(mode="after")
    def _check(self) -> Sweep:
        if self.param == "ratio":
            for v in (self.start, self.stop):
                if not isinstance(v, float) or not 0.0 <= v <= 1.0:
                    raise ValueError("ratio 스윕의 start/stop 은 0~1 사이 숫자여야 함")
        else:
            dim = Dimension.TEMPERATURE if self.param == "T" else Dimension.PRESSURE
            for v in (self.start, self.stop):
                if not isinstance(v, Quantity):
                    raise ValueError(f"{self.param} 스윕의 start/stop 은 단위가 있는 수량이어야 함")
                _require(v, dim, f"{self.param} 스윕")
        return self


class Recipe(_Model):
    id: str = Field(pattern=r"^rcp-[a-z0-9][a-z0-9\-]*$")
    name: str
    mode: Mode = Mode.MIXTURE
    components: list[Component] = Field(min_length=1)
    conditions: Conditions = Conditions()
    sweep: Sweep | None = None
    assays: Literal["auto"] | list[str] = "auto"
    notes: str | None = None

    @field_validator("assays")
    @classmethod
    def _check_assay_codes(cls, v: Literal["auto"] | list[str]) -> Literal["auto"] | list[str]:
        if v == "auto":
            return v
        bad = [c for c in v if not ASSAY_CODE_RE.match(c)]
        if bad:
            raise ValueError(f"시험 코드 형식이 아님: {bad} (예: S0, A1, A2)")
        if len(set(v)) != len(v):
            raise ValueError("시험 코드가 중복됨")
        return v

    @model_validator(mode="after")
    def _check_components(self) -> Recipe:
        refs = [str(c.ref) for c in self.components]
        dupes = sorted({r for r in refs if refs.count(r) > 1})
        if dupes:
            raise ValueError(f"같은 성분이 두 번 나옴: {dupes}")

        if self.mode is Mode.SYSTEM:
            if len(self.components) < 2:
                raise ValueError("system 모드는 원소가 2개 이상이어야 함")
            for c in self.components:
                if c.ref.namespace is not Namespace.ELEMENT:
                    raise ValueError(f"system 모드 성분은 원소만 가능: {c.ref}")
                if c.amount is not None:
                    raise ValueError(f"system 모드 성분에는 양을 적지 않음: {c.ref}")
        else:
            missing = [str(c.ref) for c in self.components if c.amount is None]
            if missing:
                raise ValueError(f"mixture 모드는 모든 성분에 양이 필요함: {missing}")
            self._check_fractions()

        if self.sweep is not None and self.sweep.param == "ratio":
            if self.mode is not Mode.MIXTURE or len(self.components) != 2:
                raise ValueError("ratio 스윕은 성분이 2개인 mixture 에서만 가능함")
        return self

    def _check_fractions(self) -> None:
        amounts = [c.amount for c in self.components if c.amount is not None]
        fractions = [q for q in amounts if q.dimension in FRACTION_DIMENSIONS]
        if not fractions:
            return
        if len(fractions) != len(amounts) or len({q.dimension for q in fractions}) > 1:
            raise ValueError("분율(wt%/vol%/mol%)을 쓰면 모든 성분이 같은 종류의 분율이어야 함")
        total = sum(q.value for q in fractions)
        if abs(total - 100.0) > FRACTION_SUM_TOLERANCE:
            raise ValueError(f"분율 합이 100%가 아님: {total:g}%")


def load_recipe(path: Path, known_assays: set[str] | None = None) -> Recipe:
    """YAML 레시피를 읽어 검증한다. known_assays 를 주면 시험 코드가 카탈로그에 있는지도 본다."""
    try:
        data = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise RecipeError(f"읽기 실패: {exc}") from exc
    try:
        recipe = Recipe.model_validate(data)
    except ValidationError as exc:
        lines = [f"{'.'.join(map(str, e['loc'])) or '(레시피)'}: {e['msg']}" for e in exc.errors()]
        raise RecipeError("\n".join(lines)) from exc
    if known_assays is not None and recipe.assays != "auto":
        unknown = [c for c in recipe.assays if c not in known_assays]
        if unknown:
            raise RecipeError(f"assays: 시험 카탈로그에 없는 코드 {unknown}")
    return recipe
