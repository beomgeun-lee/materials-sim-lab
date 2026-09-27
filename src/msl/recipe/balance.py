"""질량수지 — 성분 양을 mol·g 으로 바꾸고 원소 합계와 질량 분율을 낸다 (계획서 6장 2단계 '레시피 검증').

바꿀 수 있는 양
    mol·질량        → 절대량
    물의 부피(L)     → 밀도 1 g/mL 로 가정해 질량
    몰농도(mol/L)   → 레시피의 물 부피를 용매로 보고 mol
    mol%·at%·wt%    → 상대량 (전체 합 기준)
바꿀 수 없는 양(vol%, 물이 아닌 액체의 부피, 화학식 없는 소재)은 경고와 함께 뺀다.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from pymatgen.core import Composition

from msl.resolve import Resolved
from msl.schema.quantity import Dimension


@dataclass
class Balance:
    basis: str  # absolute / relative / none
    rows: list[dict[str, Any]] = field(default_factory=list)
    elements: dict[str, float] = field(default_factory=dict)  # 원소 → mol (relative 면 상대값)
    warnings: list[str] = field(default_factory=list)

    def mass_fraction(self, ref: str) -> float | None:
        row = next((r for r in self.rows if r["ref"] == ref), None)
        return row.get("mass_fraction") if row else None

    def to_json(self) -> dict[str, Any]:
        return {"basis": self.basis, "rows": self.rows,
                "elements": {k: round(v, 6) for k, v in self.elements.items()}, "warnings": self.warnings}


def _water_liters(comps: list[Resolved]) -> float | None:
    for c in comps:
        if c.is_water and c.amount is not None and c.amount.dimension is Dimension.VOLUME:
            return c.amount.to_si() * 1000
    return None


def balance(comps: list[Resolved]) -> Balance:
    water_l = _water_liters(comps)
    absolute: dict[str, tuple[float, Composition]] = {}
    relative: dict[str, tuple[float, Composition]] = {}
    warnings: list[str] = []
    for c in comps:
        ref = str(c.ref)
        if c.formula is None or c.amount is None:
            warnings.append(f"{ref}: 화학식 또는 양이 없어 질량수지에서 뺌")
            continue
        comp = Composition(c.formula)
        molar_mass = comp.weight  # g/mol
        a, dim = c.amount, c.amount.dimension
        if dim is Dimension.AMOUNT:
            absolute[ref] = (a.to_si(), comp)
        elif dim is Dimension.MASS:
            absolute[ref] = (a.to_si() * 1000 / molar_mass, comp)
        elif dim is Dimension.VOLUME and c.is_water:
            absolute[ref] = (a.to_si() * 1e6 / molar_mass, comp)  # m³ → g (밀도 1 g/mL 가정)
        elif dim is Dimension.CONCENTRATION and water_l:
            absolute[ref] = (a.to_si() / 1000 * water_l, comp)  # mol/m³ → mol/L × L
        elif dim is Dimension.MOLE_FRACTION:
            relative[ref] = (a.to_si(), comp)
        elif dim is Dimension.MASS_FRACTION:
            relative[ref] = (a.to_si() / molar_mass, comp)
        else:
            warnings.append(f"{ref}: {a} 는 mol 로 바꿀 수 없어 질량수지에서 뺌 (밀도·용매 정보 필요)")
    if absolute and relative:
        warnings.append("절대량(mol·g)과 분율이 섞여 있어 분율 성분은 빼고 계산")
        relative = {}
    chosen, basis = (absolute, "absolute") if absolute else (relative, "relative") if relative else ({}, "none")
    total_mass = sum(n * comp.weight for n, comp in chosen.values())
    rows, elements = [], {}
    for ref, (n, comp) in chosen.items():
        mass = n * comp.weight
        rows.append({"ref": ref, "formula": comp.reduced_formula, "moles": round(n, 8), "mass_g": round(mass, 6),
                     "mass_fraction": round(mass / total_mass, 8) if total_mass else None})
        for el, amt in comp.get_el_amt_dict().items():
            elements[el] = elements.get(el, 0.0) + n * amt
    if basis == "relative":
        rows = [{**r, "mass_g": None} for r in rows]
    return Balance(basis=basis, rows=rows, elements=elements, warnings=warnings)
