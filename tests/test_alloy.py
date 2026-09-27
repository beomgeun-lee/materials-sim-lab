"""A6 합금 상평형 — pycalphad 동봉 COST507.tdb 로 네트워크 없이 검사 (실제 data/ 는 건드리지 않는다)."""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest

import msl.db as db
from msl.assays.alloy import Calc, a6, describe, mole_fractions, transitions
from msl.assays.base import Context
from msl.connectors._tdb import index_file, read_database, write_rows
from msl.connectors.cost507 import bundled
from msl.connectors.sgte_binary_collection import mark_preferred, parse_name
from msl.registry.load import load_registry
from msl.resolve import Resolved
from msl.schema.recipe import Recipe
from msl.schema.result import Status

COST507 = bundled()


def _recipe(comp: dict[str, str], sweep: tuple[str, str, int] | None = ("1300 K", "1700 K", 9)) -> Recipe:
    data = {"id": "rcp-test-alloy", "name": "test",
            "components": [{"ref": f"element:{el}", "amount": amt} for el, amt in comp.items()]}
    if sweep:
        data["sweep"] = {"param": "T", "start": sweep[0], "stop": sweep[1], "steps": sweep[2]}
    return Recipe.model_validate(data)


def _context(recipe: Recipe) -> Context:
    comps = [Resolved(ref=c.ref, name=c.ref.key, formula=c.ref.key, amount=c.amount) for c in recipe.components]
    return Context(recipe=recipe, comps=comps, registry=load_registry())


@pytest.fixture
def tdb_data(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """임시 data/ 에 COST507 한 행짜리 tdb_systems 를 만든다 (license unknown → restricted)."""
    monkeypatch.setattr(db, "DATA_DIR", tmp_path)
    folder = tmp_path / "raw" / "cost507" / "test"
    folder.mkdir(parents=True)
    shutil.copy2(COST507, folder / "COST507.tdb")
    row = index_file(folder / "COST507.tdb", "COST507.tdb", quality="evaluated", preferred=True)
    assert write_rows([row], source="cost507", version="test", license="unknown") == {"restricted": 1}
    return tmp_path


def test_index_finds_assessed_pairs() -> None:
    row = index_file(COST507, "COST507.tdb")
    pairs = row["binaries"].split(",")
    assert row["parse_ok"] and {"Cu", "Ni", "Al", "Zn"} <= set(row["system"].split("-"))
    assert {"Cu-Ni", "Al-Cu", "Al-Zn", "Cu-Zn"} <= set(pairs)
    assert "Fe-Ni" not in pairs  # COST507 에 없는 2원계 — A6 이 이상용액으로 계산하지 않게 걸러야 함


def test_cu_ni_solidus_liquidus() -> None:
    """an Mey 1992 Cu–Ni, Ni 30 at%: 평형 고상선 ≈ 1467 K, 액상선 ≈ 1513 K (C71500 1170–1240 °C 와 ±30 K)."""
    calc = Calc(read_database(str(COST507)), ["Cu", "Ni"], {"Cu": 0.7, "Ni": 0.3}, 101_325.0)
    T = np.linspace(1300, 1700, 9)
    fracs = calc.at(T)
    assert all(abs(sum(f.values()) - 1) < 1e-6 for f in fracs)
    assert set(fracs[0]) == {"FCC_A1"} and set(fracs[-1]) == {"LIQUID"}
    tr = transitions(calc, T, fracs)
    assert tr["solidus"] == pytest.approx(1467, abs=3) and tr["liquidus"] == pytest.approx(1513, abs=3)
    assert abs(tr["solidus"] - (1170 + 273.15)) < 30 and abs(tr["liquidus"] - (1240 + 273.15)) < 30
    assert tr["liquidus"] - tr["solidus"] > 30 and set(tr["below"]) == {"FCC_A1"}


def test_room_temperature_miscibility_gap() -> None:
    calc = Calc(read_database(str(COST507)), ["Cu", "Ni"], {"Cu": 0.7, "Ni": 0.3}, 101_325.0)
    (room,) = calc.at(np.array([298.15]))
    assert describe(room) == "FCC_A1 ×2 (상분리)"


def test_mole_fractions_from_wt() -> None:
    x = mole_fractions(_context(_recipe({"Cu": "70 wt%", "Ni": "30 wt%"})).comps)
    assert x["Ni"] == pytest.approx(0.3166, abs=1e-3)


def test_a6_end_to_end(tdb_data: Path) -> None:
    out = a6(_context(_recipe({"Cu": "70 at%", "Ni": "30 at%"})))
    assert out.status is Status.OK and out.sources == [("cost507", "test")]
    vals = {(v.name, v.unit): v.value for v in out.values}
    assert vals[("고상선 (평형, 액상 첫 출현)", "K")] == pytest.approx(1467, abs=3)
    assert vals[("액상선 (평형, 고상 소멸)", "°C")] == pytest.approx(1240, abs=3)
    assert "연구용 TDB" in out.caveats[0]
    pf = out.data["phase_fractions"]
    assert len(pf["T"]) == 9 and set(pf["phases"]) == {"FCC_A1", "LIQUID"}
    assert out.data["table"]["columns"][:3] == ["T (K)", "T (°C)", "안정 상"]


def test_a6_not_applicable_without_tdb(tdb_data: Path) -> None:
    out = a6(_context(_recipe({"Au": "50 at%", "Pt": "50 at%"})))
    assert out.status is Status.NOT_APPLICABLE and "공개 TDB 에 Au–Pt 없음" in out.summary


def test_a6_rejects_unassessed_pair(tdb_data: Path) -> None:
    out = a6(_context(_recipe({"Fe": "50 at%", "Ni": "50 at%"})))
    assert out.status is Status.NOT_APPLICABLE and "평가되지 않은 원소 쌍" in out.summary


def test_a6_pending_when_not_loaded(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(db, "DATA_DIR", tmp_path)
    assert a6(_context(_recipe({"Cu": "70 at%", "Ni": "30 at%"}))).status is Status.PENDING


def test_sgte_names_and_preferred() -> None:
    assert parse_name("CuNi-92Mey-LB") == {"dataset": "CuNi-92Mey-LB", "year": 1992, "variant": "LB", "system_hint": "Cu-Ni"}
    assert parse_name("AlSi-24Zob-3g")["year"] == 2024
    rows = [dict(parse_name(n), system="Cu-Ni", parse_ok=True) for n in ("CuNi-87Jan", "CuNi-92Mey-LB", "CuNi-07Tur")]
    rows += [dict(parse_name(n), system="Al-Cu", parse_ok=True) for n in ("AlCu-02Mie", "AlCu-15Lia", "AlCu-21Kro-mod")]
    mark_preferred(rows)
    assert {r["dataset"] for r in rows if r["preferred"]} == {"CuNi-92Mey-LB", "AlCu-15Lia"}
