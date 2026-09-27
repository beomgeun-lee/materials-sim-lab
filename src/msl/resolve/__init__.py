"""성분 참조 해석 — `namespace:key` → 화학식·이름·PubChem CID·반응성 그룹 (계획서 3.1절).

v0 범위:
- element, formula : pymatgen 으로 화학식을 정규화하고 PubChem 에서 CID 를 찾는다
- cas, cid, name   : PubChem 조회
- mineral          : PubChem 이름 조회 (IMA 목록 ↔ 결정 구조 매칭은 1단계)
- material         : kb/materials.yaml 의 사용자 정의 소재
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cache
from typing import Any

import yaml
from pymatgen.core import Composition, Element

from msl.env import KB_DIR
from msl.resolve import pubchem
from msl.schema.quantity import Quantity
from msl.schema.recipe import Component, State
from msl.schema.refs import Namespace, SubstanceRef


@dataclass
class Resolved:
    ref: SubstanceRef
    name: str
    formula: str | None  # pymatgen 축약 화학식 (소재는 None)
    cid: int | None = None
    cas: str | None = None  # CAS 등록번호 (cas: 참조거나 PubChem 동의어에서)
    groups: list[str] = field(default_factory=list)  # CAMEO 반응성 그룹
    state: State | None = None
    amount: Quantity | None = None
    via: str = ""  # 어떻게 해석했는지 (리포트 표시용)
    props: dict[str, Any] = field(default_factory=dict)  # 사용자 정의 소재 물성
    error: str | None = None

    @property
    def composition(self) -> Composition | None:
        return Composition(self.formula) if self.formula else None

    @property
    def is_water(self) -> bool:
        return self.formula == "H2O"

    def moles(self) -> float | None:
        """물질량(mol). 몰·질량 단위일 때만 계산한다."""
        if self.amount is None:
            return None
        dim = self.amount.dimension.value
        if dim == "amount":
            return self.amount.to_si()
        if dim == "mass" and self.composition is not None:
            return self.amount.to_si() * 1000 / self.composition.weight
        return None


@cache
def materials() -> dict[str, dict[str, Any]]:
    data = yaml.safe_load((KB_DIR / "materials.yaml").read_text(encoding="utf-8")) or []
    return {m["id"]: m for m in data}


def _from_pubchem(props: dict[str, Any] | None) -> tuple[int | None, str | None, str | None]:
    if not props:
        return None, None, None
    cid = int(props["CID"])
    formula = props.get("MolecularFormula")
    try:
        formula = Composition(formula).reduced_formula if formula else None
    except Exception:  # 유기 복합 화학식 등 pymatgen 이 못 읽는 경우
        pass
    return cid, formula, props.get("Title")


def resolve(component: Component) -> Resolved:
    ref = component.ref
    ns, key = ref.namespace, ref.key
    base = {"ref": ref, "state": component.state, "amount": component.amount}
    try:
        if ns is Namespace.MATERIAL:
            m = materials().get(key)
            if m is None:
                return Resolved(**base, name=key, formula=None, error=f"kb/materials.yaml 에 {key!r} 없음")
            return Resolved(**base, name=m["name"], formula=None, via="kb/materials.yaml", props=m)
        if ns is Namespace.ELEMENT:
            name = Element(key).long_name
            cid, _, _ = _from_pubchem(pubchem.lookup(name))
            return _with_groups(Resolved(**base, name=name, formula=key, cid=cid, via="원소 기호"))
        if ns is Namespace.FORMULA:
            formula = Composition(key).reduced_formula
            cid, _, title = _from_pubchem(pubchem.lookup(key))
            return _with_groups(Resolved(**base, name=title or formula, formula=formula, cid=cid, via="화학식"))
        if ns is Namespace.CID:
            cid, formula, title = _from_pubchem(pubchem.by_cid(int(key)))
        elif ns in (Namespace.CAS, Namespace.NAME, Namespace.MINERAL):
            cid, formula, title = _from_pubchem(pubchem.lookup(key))
        else:
            return Resolved(**base, name=key, formula=None, error=f"{ns.value}: 해석기 v0 미지원 (1단계)")
        if cid is None:
            return Resolved(**base, name=key, formula=None, error=f"PubChem 에서 {key!r} 를 찾지 못함")
        via = "PubChem 이름 조회 (IMA 매칭은 1단계)" if ns is Namespace.MINERAL else f"PubChem ({ns.value})"
        name = f"{key} ({title})" if ns is Namespace.MINERAL and title else (title or key)
        return _with_groups(Resolved(**base, name=name, formula=formula, cid=cid, via=via))
    except Exception as exc:  # 네트워크 오류 등 — 리포트에 그대로 보여 준다
        return Resolved(**base, name=key, formula=None, error=f"해석 실패: {exc}")


def _with_groups(r: Resolved) -> Resolved:
    if r.ref.namespace is Namespace.CAS:
        r.cas = r.ref.key
    if r.cid is not None:
        r.groups = pubchem.reactive_groups(r.cid)
        if r.cas is None:
            r.cas = pubchem.cas_number(r.cid)
    return r
