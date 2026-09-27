"""조합 공간 생성기 — 레시피 틀 하나에서 여러 레시피를 만들고, 시험 결과를 한 표로 모은다 (계획서 6장 3단계).

생성 방식 (generator.kind)
    ratio    2~4 성분의 혼합비 격자 (분율 step 간격, 합 1). total 을 주면 그 절대량(예: 1 mol)을 나눈다
    amount   한 성분의 양을 start~stop 사이 steps 개로 (예: 적정 곡선). 0 이면 그 성분을 뺀다
    subsets  pool 에서 size 개씩 뽑는 모든 조합 (+ require 는 항상 포함). base.mode: system 이면 원소계 탐색

base 에는 레시피 필드(conditions, assays, mode, sweep)와 고정 성분(components, 예: 물 1 L)을 적는다.
collect 는 결과 표의 열이다: {assay, value(값 이름 또는 앞부분), unit} 이나 {assay, field(verdict·summary·status 등)}.
"""

from __future__ import annotations

import itertools
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Annotated, Any, Callable, Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, ValidationError, model_validator

from msl.registry.load import Registry
from msl.schema.quantity import Quantity
from msl.schema.recipe import Recipe, RecipeError, State


class _Model(BaseModel):
    model_config = ConfigDict(extra="forbid")


class ComponentSpec(_Model):
    ref: str
    state: State | None = None
    amount: str | None = None
    label: str | None = None  # 표·축에 쓸 이름 (없으면 참조 키)

    @model_validator(mode="before")
    @classmethod
    def _from_str(cls, data: Any) -> Any:
        return {"ref": data} if isinstance(data, str) else data

    @property
    def short(self) -> str:
        """표·축 이름: label, 없으면 참조 키 (element:Cu → Cu, formula:MgO → MgO)."""
        return self.label or self.ref.split(":", 1)[-1]

    def component(self, amount: str | None = None) -> dict[str, Any]:
        out: dict[str, Any] = {"ref": self.ref}
        if amount or self.amount:
            out["amount"] = amount or self.amount
        if self.state:
            out["state"] = self.state.value
        if self.label:
            out["label"] = self.label
        return out


class RatioGen(_Model):
    kind: Literal["ratio"]
    components: list[ComponentSpec] = Field(min_length=2, max_length=4)
    step: float = Field(gt=0, le=0.5)
    unit: str = "mol%"  # total 이 없을 때 성분 양의 단위 (mol%·at%·wt%·vol%)
    total: str | None = None  # 예: "1 mol", "10 g" — 주면 분율 × total
    include_pure: bool = False  # 한 성분만 100% 인 끝점도 넣을지

    def fractions(self) -> list[tuple[float, ...]]:
        n = round(1 / self.step)
        if abs(n * self.step - 1) > 1e-9:
            raise ValueError(f"step {self.step} 은 1 을 나누어떨어지게 해야 함 (예: 0.1, 0.05, 0.25)")
        k = len(self.components)
        out = []
        for cut in itertools.combinations(range(n + k - 1), k - 1):  # 별과 막대: n 을 k 개 정수로
            parts = [b - a - 1 for a, b in zip((-1, *cut), (*cut, n + k - 1))]
            if self.include_pure or all(p > 0 for p in parts):
                out.append(tuple(p / n for p in parts))
        return out


class AmountGen(_Model):
    kind: Literal["amount"]
    component: ComponentSpec
    start: Quantity
    stop: Quantity
    steps: int = Field(ge=2, le=500)

    @model_validator(mode="after")
    def _same_unit(self) -> AmountGen:
        if self.start.unit != self.stop.unit:
            raise ValueError(f"start·stop 단위가 같아야 함: {self.start.unit} ≠ {self.stop.unit}")
        return self


class SubsetGen(_Model):
    kind: Literal["subsets"]
    pool: list[ComponentSpec] = Field(min_length=2)
    size: int | list[int] = 2
    require: list[ComponentSpec] = Field(default_factory=list)

    @property
    def sizes(self) -> list[int]:
        return [self.size] if isinstance(self.size, int) else self.size


Generator = Annotated[RatioGen | AmountGen | SubsetGen, Field(discriminator="kind")]


class Collect(_Model):
    assay: str
    value: str | None = None  # ResultValue 이름 (같거나 앞부분이 같으면)
    unit: str | None = None  # 같은 이름이 단위별로 여럿일 때 (예: K / °C)
    field: str | None = None  # data 키 (verdict, summary …) 또는 status
    label: str | None = None
    chart: str | None = None  # 같은 이름끼리 한 차트에 겹쳐 그린다 (예: 고상선·액상선)

    @model_validator(mode="after")
    def _one(self) -> Collect:
        if (self.value is None) == (self.field is None):
            raise ValueError("collect 는 value 와 field 중 하나만")
        return self

    @property
    def title(self) -> str:
        return self.label or f"{self.assay} {self.value or self.field}" + (f" ({self.unit})" if self.unit else "")


class Space(_Model):
    id: str = Field(pattern=r"^spc-[a-z0-9][a-z0-9\-]*$")
    name: str
    description: str | None = None
    base: dict[str, Any] = Field(default_factory=dict)
    generator: Generator
    collect: list[Collect] = Field(min_length=1)
    max_variants: int = Field(default=300, ge=1, le=5000)


@dataclass
class Variant:
    recipe: Recipe
    label: str
    params: dict[str, Any]


def load_space(path: Path) -> Space:
    try:
        return Space.model_validate(yaml.safe_load(Path(path).read_text(encoding="utf-8")))
    except (OSError, yaml.YAMLError) as exc:
        raise RecipeError(f"읽기 실패: {exc}") from exc
    except ValidationError as exc:
        raise RecipeError("\n".join(f"{'.'.join(map(str, e['loc'])) or '(공간)'}: {e['msg']}" for e in exc.errors())) from exc


def _fmt(x: float) -> str:
    return f"{x:.6g}"


def expand(space: Space) -> list[Variant]:
    """공간의 모든 레시피. 레시피 검증 오류는 어느 변형에서 났는지와 함께 RecipeError 로 알린다."""
    base = dict(space.base)
    fixed = [ComponentSpec.model_validate(c).component() for c in base.pop("components", [])]
    g = space.generator
    raw: list[tuple[list[dict], str, dict[str, Any]]] = []
    if isinstance(g, RatioGen):
        for fr in g.fractions():
            comps, params = [], {}
            for c, x in zip(g.components, fr):
                params[c.short] = round(x, 6)
                if x <= 0:
                    continue
                if g.total:
                    tot = Quantity.model_validate(g.total)
                    comps.append(c.component(f"{_fmt(tot.value * x)} {tot.unit}"))
                else:
                    comps.append(c.component(f"{_fmt(x * 100)} {g.unit}"))
            label = " · ".join(f"{k} {_fmt(v * 100)}" for k, v in params.items()) + ("" if g.total else f" {g.unit}")
            raw.append((comps, label, params))
    elif isinstance(g, AmountGen):
        a, b = g.start.value, g.stop.value
        for i in range(g.steps):
            v = a + (b - a) * i / (g.steps - 1)
            comps = [g.component.component(f"{_fmt(v)} {g.start.unit}")] if v > 0 else []
            raw.append((comps, f"{g.component.short} {_fmt(v)} {g.start.unit}", {g.component.short: round(v, 9)}))
    else:
        for k in g.sizes:
            for combo in itertools.combinations(g.pool, k):
                members = [*combo, *g.require]
                raw.append(([c.component() for c in members], " + ".join(c.short for c in combo),
                            {"members": [c.short for c in combo]}))
    if len(raw) > space.max_variants:
        raise RecipeError(f"변형 {len(raw)}개가 max_variants {space.max_variants} 를 넘음 — step·size 를 줄이거나 한도를 올릴 것")
    stem = space.id.removeprefix("spc-")
    out = []
    for i, (comps, label, params) in enumerate(raw):
        data = {**base, "id": f"rcp-{stem}-{i:03d}", "name": f"{space.name} · {label}", "components": comps + fixed}
        try:
            out.append(Variant(Recipe.model_validate(data), label, params))
        except ValidationError as exc:
            msg = "; ".join(f"{'.'.join(map(str, e['loc']))}: {e['msg']}" for e in exc.errors())
            raise RecipeError(f"변형 {i} ({label}): {msg}") from exc
    return out


def pick(results: list, c: Collect) -> Any:
    """한 레시피의 시험 결과에서 collect 가 가리키는 값."""
    r = next((x for x in results if x.assay == c.assay), None)
    if r is None:
        return None
    if c.field == "status":
        return r.status.value
    if c.field:
        return (r.data or {}).get(c.field)
    for v in r.values:
        if (v.name == c.value or v.name.startswith(c.value)) and (c.unit is None or v.unit == c.unit):
            return v.value
    return None


@dataclass
class SpaceResult:
    space: Space
    rows: list[dict[str, Any]] = field(default_factory=list)  # {label, params, values{title: v}, status{assay: s}}
    elapsed: float = 0.0

    @property
    def columns(self) -> list[str]:
        return [c.title for c in self.space.collect]

    def axis(self) -> tuple[str, list[float]] | None:
        """1차원 공간이면 x 축 (이름, 값). 이원 혼합비는 첫 성분 분율, amount 는 양."""
        g = self.space.generator
        if isinstance(g, RatioGen) and len(g.components) == 2:
            k = g.components[0].short
            return f"{k} 분율", [r["params"][k] for r in self.rows]
        if isinstance(g, AmountGen):
            k = g.component.short
            return f"{k} ({g.start.unit})", [r["params"][k] for r in self.rows]
        return None

    def matrix(self) -> dict[str, Any] | None:
        """2개씩 뽑은 subsets 이면 쌍별 행렬 (첫 collect 기준)."""
        g = self.space.generator
        if not (isinstance(g, SubsetGen) and g.sizes == [2]):
            return None
        names = [c.short for c in g.pool]
        title = self.columns[0]
        cells = {tuple(r["params"]["members"]): r["values"].get(title) for r in self.rows}
        return {"names": names, "title": title,
                "cells": [[cells.get((a, b), cells.get((b, a))) if a != b else None for b in names] for a in names]}

    def to_json(self) -> dict[str, Any]:
        ax = self.axis()
        return {"space": self.space.model_dump(mode="json"), "columns": self.columns, "rows": self.rows,
                "elapsed": round(self.elapsed, 2), "axis": {"name": ax[0], "values": ax[1]} if ax else None,
                "matrix": self.matrix()}


def run_space(space: Space, registry: Registry, use_cache: bool = True,
              progress: Callable[[int, int, Variant], None] | None = None) -> SpaceResult:
    from msl.runtime.runner import run_recipe

    variants = expand(space)
    out = SpaceResult(space)
    t0 = time.monotonic()
    for i, v in enumerate(variants):
        if progress:
            progress(i, len(variants), v)
        rep = run_recipe(v.recipe, use_cache=use_cache, registry=registry)
        out.rows.append({
            "label": v.label, "params": v.params, "recipe_id": v.recipe.id,
            "values": {c.title: pick(rep.results, c) for c in space.collect},
            "status": {r.assay: r.status.value for r in rep.results},
            "summary": {r.assay: (r.data or {}).get("summary", "") for r in rep.results},
        })
    out.elapsed = time.monotonic() - t0
    return out
