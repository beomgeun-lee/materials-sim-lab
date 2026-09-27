"""저장 계층 — 프로비넌스 검사, 라이선스 파티션, 조회 (실제 data/ 는 건드리지 않는다)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

import msl.db as db
from msl.connectors import available
from msl.connectors.phreeqc_db import parse_masters, parse_phases


@pytest.fixture
def tmp_data(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    monkeypatch.setattr(db, "DATA_DIR", tmp_path)
    return tmp_path


def _rows(license_id: str, n: int = 3) -> pd.DataFrame:
    df = pd.DataFrame({"symbol": [f"X{i}" for i in range(n)], "value": range(n)})
    return db.with_provenance(df, source="periodictable", version="t", license=license_id,
                              method="compiled", id_column="symbol")


def test_open_license_goes_to_open(tmp_data: Path) -> None:
    assert db.write_table("elements", _rows("CC0-1.0")) == {"open": 3}
    assert (tmp_data / "open" / "elements" / "periodictable.parquet").exists()


@pytest.mark.parametrize("license_id", ["unknown", "CC-BY-NC-4.0", "NIST-SRD", "no-such-license"])
def test_non_open_license_goes_to_restricted(tmp_data: Path, license_id: str) -> None:
    assert db.write_table("elements", _rows(license_id)) == {"restricted": 3}


def test_mixed_licenses_are_split(tmp_data: Path) -> None:
    """MP 처럼 한 소스 안에 BY 와 BY-NC 가 섞이면 행 단위로 나뉘어야 한다."""
    df = pd.concat([_rows("CC-BY-4.0", 2), _rows("CC-BY-NC-4.0", 1)])
    assert db.write_table("phases", df) == {"open": 2, "restricted": 1}
    open_only = db.connect(include_restricted=False).sql("SELECT count(*) FROM phases").fetchone()[0]
    everything = db.connect().sql("SELECT count(*) FROM phases").fetchone()[0]
    assert (open_only, everything) == (2, 3)


def test_missing_provenance_is_rejected(tmp_data: Path) -> None:
    with pytest.raises(db.ProvenanceError, match="프로비넌스 열 없음"):
        db.write_table("elements", pd.DataFrame({"symbol": ["Fe"]}))
    bad = _rows("CC0-1.0")
    bad.loc[0, "source_version"] = ""
    with pytest.raises(db.ProvenanceError, match="빈 행"):
        db.write_table("elements", bad)


def test_status_reports_rows(tmp_data: Path) -> None:
    db.write_table("elements", _rows("CC0-1.0", 5))
    (row,) = db.status()
    assert (row.table, row.source, row.partition, row.rows) == ("elements", "periodictable", "open", 5)


def test_checksums(tmp_path: Path) -> None:
    (tmp_path / "a.txt").write_text("hello")
    text = db.record_checksums(tmp_path).read_text()
    assert "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824  a.txt" in text


def test_connectors_follow_contract() -> None:
    reg = __import__("msl.registry.load", fromlist=["load_registry"]).load_registry()
    for sid, mod in available().items():
        assert sid in reg.sources, f"{sid}: sources.yaml 에 없는 소스"
        assert callable(mod.fetch) and callable(mod.load) and mod.TABLES


PHREEQC_SAMPLE = """
SOLUTION_MASTER_SPECIES
Ca       Ca+2           0.0     Ca              40.08
C(4)     HCO3-          1.0     HCO3            12.0111
PHASES
Calcite
        CaCO3 = CO3-2 + Ca+2
        -log_k  -8.48;  -delta_h -2.297 kcal
        -analytic -171.9065 -0.077993 2839.319 71.595
CO2(g)
        CO2 = CO2
        log_k -1.468
EXCHANGE_MASTER_SPECIES
"""


def test_phreeqc_parsers() -> None:
    phases = {p["phase"]: p for p in parse_phases(PHREEQC_SAMPLE)}
    assert phases["Calcite"]["formula"] == "CaCO3" and phases["Calcite"]["log_k"] == -8.48
    assert phases["CO2(g)"]["is_gas"] and phases["CO2(g)"]["log_k"] == -1.468
    assert [m["master"] for m in parse_masters(PHREEQC_SAMPLE)] == ["Ca", "C(4)"]
