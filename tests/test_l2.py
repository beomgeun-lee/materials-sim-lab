"""L2 — 치환 배치, 캐시 키, MP2020 보정 엔트리, 허용 가중치 (무거운 이완은 l2 설치·캐시가 있을 때만)."""

from __future__ import annotations

import importlib.util

import numpy as np
import pytest
from pymatgen.core import Composition, Lattice, Structure

from msl import l2
from msl.engines import umlip

HAVE_L2 = importlib.util.find_spec("mace") is not None and importlib.util.find_spec("orb_models") is not None
needs_l2 = pytest.mark.skipif(not HAVE_L2, reason="uMLIP 없음 — uv sync --extra l2")


def rocksalt() -> Structure:
    return Structure(Lattice.cubic(4.2), ["Ni", "O"], [[0, 0, 0], [0.5, 0.5, 0.5]]).get_primitive_structure()


def test_orderings_rocksalt() -> None:
    t = Composition("Li6MnNi3O10")
    prim = rocksalt()
    out = l2.orderings(prim, t, 4, 40)
    assert len(out) == 4 and all(len(s) == 20 for s in out)
    assert all(s.composition.reduced_composition.almost_equals(t.reduced_composition) for s in out)
    o_sites = [[tuple(np.round(site.frac_coords, 3)) for site in s if str(site.specie) == "O"] for s in out]
    assert all(x == o_sites[0] for x in o_sites)  # 음이온 자리는 그대로
    assert len({tuple(str(site.specie) for site in s) for s in out}) > 1  # 양이온 배치는 서로 다름
    again = l2.orderings(prim, t, 4, 40)
    assert [tuple(str(x.specie) for x in s) for s in again] == [tuple(str(x.specie) for x in s) for s in out]  # 시드 고정
    assert l2.orderings(prim, t, 4, 10) == []  # 20원자가 필요한데 상한 10
    assert l2._key("a", out[0]) == l2._key("a", out[0].copy()) != l2._key("b", out[0])


def test_mp_entry_corrections() -> None:
    s = Structure(Lattice.cubic(4.3), ["Fe", "O"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    e = umlip.mp_entry(s, -15.0)
    assert e.parameters["run_type"] == "GGA+U" and e.parameters["hubbards"]["Fe"] == 5.3
    names = {a.name for a in e.energy_adjustments}
    assert any("oxide" in n for n in names) and any("Fe" in n for n in names)
    mgo = umlip.mp_entry(Structure(Lattice.cubic(4.2), ["Mg", "O"], [[0, 0, 0], [0.5, 0.5, 0.5]]), -12.0)
    assert mgo.parameters["run_type"] == "GGA"


def test_only_allowed_weights() -> None:
    with pytest.raises(ValueError, match="허용되지 않은"):
        umlip.calculator("medium-omat-0")  # ASL(비상업) 가중치 — 허용 목록 밖


@needs_l2
def test_single_point_matches_mp_for_gga() -> None:
    """GGA 물질(MgO 암염)은 MP 원시 DFT 에너지와 원자당 0.02 eV 안 — LASPH 영향이 작은 s·p 원소."""
    s = Structure(Lattice.cubic(4.25), ["Mg"] * 4 + ["O"] * 4,
                  [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5], [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]])
    for name in umlip.MODELS:
        r = umlip.relax(s, name)
        assert r.converged and -6.02 < r.energy / len(s) < -5.95  # MP mp-1265 원시 약 −5.98 eV/atom


def test_l2_web_saved_results(tmp_path, monkeypatch) -> None:
    """저장된 L2 결과는 계산 없이 바로 돌려준다 · 잘못된 화학식은 422."""
    from fastapi import HTTPException

    from msl.web import app as web

    monkeypatch.setattr(l2, "RESULTS", tmp_path)
    fake = {"formula": "MgAl2O4", "known": "mp-3536-GGA", "mp_ehull": 0.0, "n_structures": 1, "ehull_mean": 0.0,
            "ehull_spread": 0.0, "models": {"mace-mpa-0": {"ehull": 0.0, "decomposition": {"MgAl2O4": 1.0}}}, "seconds": 1.0}
    l2.save(fake)
    assert l2.saved("Mg2Al4O8")["formula"] == "MgAl2O4"  # 약분 조성으로 찾는다
    r = web.l2_start(web.L2Request(formula="MgAl2O4"))
    assert r["status"] == "done" and r["cached"] and r["result"]["ehull_mean"] == 0.0
    assert set(web.l2_saved(web.L2SavedRequest(formulas=["MgAl2O4", "NaCl"]))) == {"MgAl2O4"}
    with pytest.raises(HTTPException):
        web.l2_start(web.L2Request(formula="Xx9"))
    with pytest.raises(HTTPException):
        web.l2_status("없는작업")


@pytest.mark.parametrize(("formula", "expect"), [
    ("Li6MnNi3O10", {"Li": 1, "Mn": 3.5, "Ni": 3.5, "O": -2}),  # Ni 는 ICSD 에 +1~+4 → 가변 (pymatgen 추정은 Mn⁷⁺ 를 고름)
    ("LiCoO2", {"Li": 1, "Co": 3, "O": -2}), ("Mg2SiO4", {"Mg": 2, "Si": 4, "O": -2}), ("Li2MnO3", {"Li": 1, "Mn": 4, "O": -2}),
])
def test_charges(formula: str, expect: dict) -> None:
    q = l2.charges(Composition(formula))
    assert q == pytest.approx(expect)
    assert sum(Composition(formula)[e] * v for e, v in q.items()) == pytest.approx(0)


def test_ewald_orderings_rocksalt() -> None:
    t = Composition("Li6MnNi3O10")
    out = l2.ewald_orderings(rocksalt(), t, 4, 40)
    assert 1 <= len(out) <= 4 and all(len(s) == 20 for s in out)
    assert all(s.composition.reduced_composition.almost_equals(t.reduced_composition) for s in out)
    assert all(type(site.specie).__name__ == "Element" for s in out for site in s)  # 산화수 표시는 지운다
    li_sets = {frozenset(tuple(np.round(site.frac_coords, 3)) for site in s if str(site.specie) == "Li") for s in out}
    assert len(li_sets) == len(out)  # 서로 다른 배치 (Ewald 변환은 원자 순서를 원소별로 정렬해 돌려준다)


@pytest.mark.skipif(not (l2.mp.CACHE / "entries_Li-Mn-Ni-O_GGA_GGA+U.pkl").exists(), reason="MP Li-Mn-Ni-O 캐시 없음")
def test_candidates_mix() -> None:
    labels = [lb for lb, _ in l2.candidates(Composition("Li6MnNi3O10"), use_prototypes=False)]
    assert any("Ewald" in lb for lb in labels) and any("무작위" in lb for lb in labels)
    assert not any("Li6MnNi3O10" in lb.split("(")[0] for lb in labels)


@pytest.mark.skipif(not (l2.mp.CACHE / "proto_AB2C4_O.pkl").exists(), reason="원형 캐시 없음 (MP 조회 필요)")
def test_prototypes_spinel() -> None:
    """LiMn2O4 를 모른다고 치면 다른 화학계의 스피넬(AB2O4)을 Li·Mn 으로 장식한다."""
    got = l2.prototypes(Composition("LiMn2O4"), exclude_chemsys="Li-Mn-O")
    assert got and all(s.composition.reduced_composition.almost_equals(Composition("LiMn2O4").reduced_composition) for _, s in got)
    assert all("Li-Mn-O" not in lb for lb, _ in got)
