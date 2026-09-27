"""성분 참조 해석 — `namespace:key` → 화학식·이름·PubChem CID·반응성 그룹 (계획서 3.1절).

v0 범위:
- element, formula : pymatgen 으로 화학식을 정규화하고 PubChem 에서 CID 를 찾는다
- cas, cid, name   : PubChem 조회. 한글 이름은 KOSHA 국문명 검색으로 CAS 를 찾은 뒤 PubChem (data.go.kr 키 필요)
- ke               : KOSHA 에서 KE 번호 → CAS → PubChem
- mineral          : 적재된 IMA 목록(minerals)·COD 매칭(mineral_structures)에서 먼저 찾고, 없으면 PubChem 이름 조회
- material         : kb/materials.yaml 의 사용자 정의 소재
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from functools import cache
from typing import Any

import yaml
from pymatgen.core import Composition, Element

from msl.env import KB_DIR
from msl.resolve import pubchem
from msl.schema.quantity import Quantity
from msl.schema.recipe import Component, State
from msl.schema.refs import Namespace, SubstanceRef, to_pymatgen_formula

HANGUL = re.compile(r"[가-힣]")


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
            formula = Composition(to_pymatgen_formula(key)).reduced_formula
            cid, _, title = _from_pubchem(pubchem.lookup(key))
            return _with_groups(Resolved(**base, name=title or formula, formula=formula, cid=cid, via="화학식"))
        if ns is Namespace.MINERAL and (hit := mineral_db().get(_mineral_key(key))) and \
                (hit["formula"] or hit["best_formula"]):
            return _from_mineral_db(base, key, hit)
        if ns is Namespace.KE or (ns is Namespace.NAME and HANGUL.search(key)):
            return _via_korean_registry(base, ns, key)
        if ns is Namespace.CID:
            cid, formula, title = _from_pubchem(pubchem.by_cid(int(key)))
        elif ns in (Namespace.CAS, Namespace.NAME, Namespace.MINERAL):
            cid, formula, title = _from_pubchem(pubchem.lookup(key))
        else:
            return Resolved(**base, name=key, formula=None, error=f"{ns.value}: 해석기 v0 미지원 (1단계)")
        if cid is None:
            return Resolved(**base, name=key, formula=None, error=f"PubChem 에서 {key!r} 를 찾지 못함")
        via = "PubChem 이름 조회 (IMA 목록에 없음)" if ns is Namespace.MINERAL else f"PubChem ({ns.value})"
        name = f"{key} ({title})" if ns is Namespace.MINERAL and title else (title or key)
        return _with_groups(Resolved(**base, name=name, formula=formula, cid=cid, via=via))
    except Exception as exc:  # 네트워크 오류 등 — 리포트에 그대로 보여 준다
        return Resolved(**base, name=key, formula=None, error=f"해석 실패: {exc}")


def _mineral_key(name: str) -> str:
    from msl.connectors.ima_cnmnc import name_key

    return name_key(name)


@cache
def mineral_db() -> dict[str, dict[str, Any]]:
    """적재된 IMA 광물 목록(+COD 매칭) — {광물명 비교 키: 행}. DB 가 비었거나 못 읽으면 빈 dict (PubChem 폴백)."""
    try:
        from msl.connectors.ima_cnmnc import build_index
        from msl.db import connect, tables

        have = set(tables())
        if "minerals" not in have:
            return {}
        cols = "s.best_cod_id, s.best_formula" if "mineral_structures" in have else "NULL, NULL"
        join = "LEFT JOIN mineral_structures s USING (name)" if "mineral_structures" in have else ""
        rows = connect().execute(
            f"SELECT m.name, m.formula_reduced, m.source_version, {cols} FROM minerals m {join}").fetchall()
    except Exception:
        return {}
    by_name = {r[0]: {"name": r[0], "formula": r[1], "version": r[2],
                      "best_cod_id": None if r[3] is None else int(r[3]), "best_formula": r[4]} for r in rows}
    return {k: by_name[n] for k, n in build_index(list(by_name)).items()}


def _from_mineral_db(base: dict[str, Any], key: str, hit: dict[str, Any]) -> Resolved:
    """IMA 정본명·화학식(치환식이라 축약식이 없으면 대표 COD 구조의 화학식). CID·CAS·반응성 그룹은 PubChem 이름 조회로 채운다."""
    formula = hit["formula"] or hit["best_formula"]
    via = f"IMA {hit['version']}" + (f" · COD {hit['best_cod_id']}" if hit["best_cod_id"] else "")
    if not hit["formula"] and formula:
        via += " (화학식은 COD 구조)"
    try:
        cid, _, _ = _from_pubchem(pubchem.lookup(key))
    except Exception:  # PubChem 이 안 돼도 광물 해석은 유지한다
        cid = None
    r = Resolved(**base, name=hit["name"], formula=formula, cid=cid, via=via)
    try:
        return _with_groups(r)
    except Exception:
        return r


def _via_korean_registry(base: dict[str, Any], ns: Namespace, key: str) -> Resolved:
    """국문명·KE 번호 → (KOSHA MSDS 목록) → CAS → PubChem. 이름은 정확히 일치할 때만 채택한다."""
    from msl.engines import datagokr as dg

    try:
        cas = dg.cas_from_ke(key) if ns is Namespace.KE else dg.cas_from_korean_name(key)
    except dg.DataGoKrUnavailable as exc:
        return Resolved(**base, name=key, formula=None, error=f"국문명·KE 해석에는 data.go.kr 키가 필요: {exc}")
    if cas is None:
        return Resolved(**base, name=key, formula=None, error=f"KOSHA 목록에서 {key!r} 와 정확히 같은 이름을 찾지 못함")
    cid, formula, title = _from_pubchem(pubchem.lookup(cas))
    if cid is None:
        return Resolved(**base, name=key, formula=None, cas=cas, error=f"CAS {cas} 를 PubChem 에서 찾지 못함")
    r = Resolved(**base, name=f"{key} ({title})" if title else key, formula=formula, cid=cid, cas=cas,
                 via=f"KOSHA {'KE' if ns is Namespace.KE else '국문명'} → CAS {cas} → PubChem")
    return _with_groups(r)


def _with_groups(r: Resolved) -> Resolved:
    if r.ref.namespace is Namespace.CAS:
        r.cas = r.ref.key
    if r.cid is not None:
        r.groups = pubchem.reactive_groups(r.cid)
        if r.cas is None:
            r.cas = pubchem.cas_number(r.cid)
    return r
