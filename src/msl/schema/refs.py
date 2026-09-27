"""성분 참조(`namespace:key`) 형식 검증.

레시피의 각 성분은 `element:Fe`, `mineral:calcite`, `cas:7647-01-0` 처럼 참조한다.
여기서는 형식만 검사한다. 실제 엔티티로 해석하는 일은 해석기(resolve, 1단계)가 맡는다.
"""

from __future__ import annotations

import re
from enum import StrEnum
from typing import Any

from pydantic import BaseModel, ConfigDict, model_serializer, model_validator

ELEMENTS: tuple[str, ...] = (
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
    "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga",
    "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
    "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os",
    "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa",
    "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg",
    "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
)  # fmt: skip
_ELEMENT_SET = frozenset(ELEMENTS)


class Namespace(StrEnum):
    ELEMENT = "element"  # element:Fe
    FORMULA = "formula"  # formula:Fe2O3
    MP = "mp"  # mp:mp-19770 (Materials Project ID)
    COD = "cod"  # cod:9000660 (COD ID)
    MINERAL = "mineral"  # mineral:hematite (IMA 정본명)
    CAS = "cas"  # cas:7647-01-0
    CID = "cid"  # cid:313 (PubChem)
    INCHIKEY = "inchikey"
    KE = "ke"  # ke:KE-20189 (국내 기존화학물질 고유번호)
    NAME = "name"  # name:표백제 — 해석기가 후보를 제시한다
    MATERIAL = "material"  # material:epoxy-generic — 사용자·상용 등급 소재


_FORMULA_RE = re.compile(r"^(?:\d*(?:\.\d+)?)(?:[A-Z][a-z]?|\d+(?:\.\d+)?|\.\d+|[()\[\]])+$")
_KEY_RULES: dict[Namespace, re.Pattern[str]] = {
    Namespace.MP: re.compile(r"^mp-[0-9a-z]+$"),
    Namespace.COD: re.compile(r"^\d{7}$"),
    Namespace.MINERAL: re.compile(r"^[A-Za-z][A-Za-z0-9\-() ]*$"),
    Namespace.CAS: re.compile(r"^\d{2,7}-\d{2}-\d$"),
    Namespace.CID: re.compile(r"^[1-9]\d*$"),
    Namespace.INCHIKEY: re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$"),
    Namespace.KE: re.compile(r"^KE-?\d{3,}$"),
    Namespace.MATERIAL: re.compile(r"^[a-z0-9][a-z0-9\-]*$"),
}


def cas_checksum_ok(cas: str) -> bool:
    """CAS 등록번호 검사 숫자 확인. 마지막 숫자 = 나머지 숫자를 오른쪽부터 1,2,3…배 한 합 mod 10."""
    digits = cas.replace("-", "")
    body, check = digits[:-1], int(digits[-1])
    return sum(int(d) * i for i, d in enumerate(reversed(body), start=1)) % 10 == check


def formula_elements(formula: str) -> list[str]:
    """화학식에 나오는 원소 기호 목록 (수화물 `·` 구분 허용)."""
    return re.findall(r"[A-Z][a-z]?", formula)


def _check_formula(formula: str) -> None:
    for part in formula.split("·"):
        if not part or not _FORMULA_RE.match(part):
            raise ValueError(f"화학식 형식이 아님: {formula!r}")
    unknown = [s for s in formula_elements(formula) if s not in _ELEMENT_SET]
    if unknown:
        raise ValueError(f"화학식에 없는 원소 기호: {unknown} ({formula!r})")
    if formula.count("(") != formula.count(")") or formula.count("[") != formula.count("]"):
        raise ValueError(f"화학식 괄호 짝이 맞지 않음: {formula!r}")


class SubstanceRef(BaseModel):
    """`namespace:key` 형식의 성분 참조. 직렬화하면 다시 문자열이 된다."""

    model_config = ConfigDict(frozen=True)

    namespace: Namespace
    key: str

    @model_validator(mode="before")
    @classmethod
    def _parse(cls, data: Any) -> Any:
        if isinstance(data, str):
            ns, sep, key = data.partition(":")
            if not sep or not key.strip():
                raise ValueError(f"성분 참조는 'namespace:key' 형식이어야 함: {data!r}")
            return {"namespace": ns.strip(), "key": key.strip()}
        return data

    @model_validator(mode="after")
    def _check_key(self) -> SubstanceRef:
        ns, key = self.namespace, self.key
        if ns is Namespace.ELEMENT:
            if key not in _ELEMENT_SET:
                raise ValueError(f"원소 기호가 아님: {key!r}")
        elif ns is Namespace.FORMULA:
            _check_formula(key)
        elif ns in _KEY_RULES:
            if not _KEY_RULES[ns].match(key):
                raise ValueError(f"{ns.value} 참조 형식이 아님: {key!r}")
            if ns is Namespace.CAS and not cas_checksum_ok(key):
                raise ValueError(f"CAS 번호 검사 숫자 불일치: {key!r}")
        return self

    @model_serializer
    def _serialize(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return f"{self.namespace.value}:{self.key}"


def to_pymatgen_formula(formula: str) -> str:
    """수화물 표기 `CuSO4·5H2O` 를 pymatgen 이 읽는 `CuSO4(H2O)5` 로 바꾼다."""
    head, *rest = formula.split("·")
    out = head
    for part in rest:
        m = re.match(r"^(\d+(?:\.\d+)?)?(.+)$", part)
        n, body = (m.group(1) or "1"), m.group(2)
        out += f"({body}){n}" if n != "1" else f"({body})"
    return out
