"""NASA 열역학 데이터 (Apache-2.0) — 평형 트랙(A2 교차검증·A3)이 쓰는 화학종 목록.

원본은 NASA CEA 의 data/thermo.inp. 엔진은 이를 변환한 Cantera 번들 nasa_gas.yaml·nasa_condensed.yaml 을
쓰므로, 적재 테이블도 엔진이 실제로 보는 화학종 기준으로 만든다 (두 원본 모두 스냅샷).
"""

from __future__ import annotations

import shutil
from importlib.metadata import version
from pathlib import Path

import cantera as ct
import pandas as pd
from pymatgen.core import Composition

from msl.db import download, raw_dir, record_checksums, with_provenance, write_table

SOURCE = "nasa-cea-thermo"
TABLES = ["thermo_species"]
VERSION = f"cantera-{version('cantera')}"
LICENSE = "Apache-2.0"
THERMO_INP = "https://raw.githubusercontent.com/nasa/cea/main/data/thermo.inp"


def fetch() -> Path:
    d = raw_dir(SOURCE, VERSION)
    data = Path(ct.__file__).parent / "data"
    for name in ("nasa_gas.yaml", "nasa_condensed.yaml"):
        shutil.copy2(data / name, d / name)
    download(THERMO_INP, d / "thermo.inp")
    record_checksums(d)
    return d


def _rows(path: Path, phase: str) -> list[dict]:
    rows = []
    for sp in ct.Species.list_from_file(str(path)):
        comp = {el: n for el, n in sp.composition.items() if el != "E"}
        try:
            formula = Composition(comp).formula.replace(" ", "") if comp else None
        except Exception:
            formula = None
        T = 298.15
        valid = sp.thermo.min_temp <= T <= sp.thermo.max_temp
        rows.append({
            "species": sp.name, "phase": phase, "formula": formula, "charge": sp.charge,
            "t_min_k": sp.thermo.min_temp, "t_max_k": sp.thermo.max_temp,
            "h298_kj_mol": sp.thermo.h(T) / 1e6 if valid else None,  # J/kmol → kJ/mol
            "s298_j_mol_k": sp.thermo.s(T) / 1e3 if valid else None,
        })
    return rows


def load() -> dict[str, dict[str, int]]:
    d = raw_dir(SOURCE, VERSION)
    if not (d / "nasa_gas.yaml").exists():
        fetch()
    df = pd.DataFrame(_rows(d / "nasa_gas.yaml", "gas") + _rows(d / "nasa_condensed.yaml", "condensed"))
    df["key"] = df["phase"] + ":" + df["species"]
    df = with_provenance(df, source=SOURCE, version=VERSION, license=LICENSE, method="compiled", id_column="key")
    return {"thermo_species": write_table("thermo_species", df)}
