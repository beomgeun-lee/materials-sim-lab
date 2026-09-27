"""mendeleev — 연구용 원소 레이어 (속성이 가장 풍부, 수록값 라이선스 미확인 → restricted, D10)."""

from __future__ import annotations

from importlib.metadata import version
from pathlib import Path

import pandas as pd

from msl.db import raw_dir, record_checksums, with_provenance, write_table

SOURCE = "mendeleev"
TABLES = ["elements"]
VERSION = version("mendeleev")
LICENSE = "unknown"
KEEP = {  # mendeleev 열 → 우리 열 (수치 위주로 추림)
    "atomic_number": "Z", "symbol": "symbol", "name": "name", "atomic_weight": "atomic_mass",
    "block": "block", "period": "period", "group_id": "group", "density": "density_g_cm3",
    "melting_point": "melting_point_k", "boiling_point": "boiling_point_k",
    "en_pauling": "electronegativity_pauling", "electron_affinity": "electron_affinity_ev",
    "covalent_radius_pyykko": "covalent_radius_pm", "vdw_radius": "vdw_radius_pm",
    "atomic_radius": "atomic_radius_pm", "dipole_polarizability": "dipole_polarizability_au",
    "thermal_conductivity": "thermal_conductivity_w_mk", "fusion_heat": "fusion_heat_kj_mol",
    "evaporation_heat": "evaporation_heat_kj_mol", "lattice_structure": "lattice_structure",
    "lattice_constant": "lattice_constant_a", "abundance_crust": "abundance_crust_mg_kg",
    "is_radioactive": "is_radioactive",
}


def fetch() -> Path:
    """패키지에 든 SQLite 를 그대로 읽는다. mendeleev.fetch 의 조회 함수는 pandas 3 와 호환되지 않아
    (1.3.0, 2026-09-27 확인) 쓰지 않는다."""
    import sqlite3

    import mendeleev as md

    d = raw_dir(SOURCE, VERSION)
    with sqlite3.connect(Path(md.__file__).parent / "elements.db") as con:
        pd.read_sql_query("SELECT * FROM elements", con).to_csv(d / "elements.csv", index=False)
        pd.read_sql_query(
            "SELECT atomic_number AS Z, ionization_energy AS ionization_energy_ev FROM ionizationenergies "
            "WHERE ion_charge = 0 ORDER BY atomic_number", con,
        ).to_csv(d / "ionization_energies.csv", index=False)
    record_checksums(d)
    return d


def load() -> dict[str, dict[str, int]]:
    d = raw_dir(SOURCE, VERSION)
    if not (d / "ionization_energies.csv").exists():
        fetch()
    raw = pd.read_csv(d / "elements.csv")
    df = raw[[c for c in KEEP if c in raw.columns]].rename(columns=KEEP)
    df = df.merge(pd.read_csv(d / "ionization_energies.csv"), on="Z", how="left")
    df = with_provenance(df, source=SOURCE, version=VERSION, license=LICENSE, method="compiled", id_column="symbol")
    return {"elements": write_table("elements", df)}
