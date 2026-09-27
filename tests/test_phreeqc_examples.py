"""PHREEQC 공식 예제 재현 — A4 경로(phreeqpython IPhreeqc)가 공식 phreeqc 3.8.6 과 같은 답을 내는가.

기준값은 USGS 배포본 소스로 빌드한 공식 실행 파일의 결과(kb/validation/phreeqc_examples_3.8.6.json).
예제 입력은 tests/fixtures/phreeqc (배포본에서 수정 없이 복사, NOTICE 동봉). 실행 파일·네트워크 없이 돈다.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from msl.bench.phreeqc import bundled_db, compare, parse, passes, reference, run_ours

ROOT = Path(__file__).resolve().parents[1]
REF = reference(ROOT / "kb" / "validation" / "phreeqc_examples_3.8.6.json")
FIX = ROOT / "tests" / "fixtures" / "phreeqc"


def test_reference_summary() -> None:
    """비교 가능한 공식 예제는 모두 재현, 비교 불가는 USER_GRAPH 전용 3개뿐."""
    rows = REF["examples"]
    assert len(rows) == 31
    assert all(r["reproduced"] for r in rows.values() if r["comparable"])
    assert sorted(k for k, r in rows.items() if not r["comparable"]) == ["ex19", "ex19b", "ex2b"]


@pytest.mark.parametrize("name", ["ex1", "ex3", "ex4", "ex5", "ex9"])
def test_official_example_reproduced(name: str) -> None:
    row = REF["examples"][name]
    ex = FIX / name
    assert hashlib.sha256(ex.read_bytes()).hexdigest() == row["input_sha256"]  # 배포본 원본 그대로
    got = parse(run_ours(ex, bundled_db())["out"])
    c = compare(row["official"], got)
    assert passes(c), c


def test_seawater_speciation_values() -> None:
    """예제 1 (해수 종분화, A4 의 동봉 phreeqc.dat): pH 8.220, pe 8.451, 이온 세기 0.6747 mol/kgw, 방해석 SI 0.75 (배포본 3.8.6 DB 는 0.6704·0.79)."""
    got = parse(run_ours(FIX / "ex1", bundled_db())["out"])
    sol = got["solutions"][0]
    assert (sol["pH"], sol["pe"]) == (8.22, 8.451) and sol["I"] == pytest.approx(0.6747, abs=1e-4)
    assert got["si"][0]["Calcite"] == pytest.approx(0.75, abs=0.005)  # 3.8.6 DB 는 0.79
