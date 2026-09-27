"""PubChem 주기율표 CSV (NLM 생성분 퍼블릭 도메인) — 배포용 원소 레이어 교차검증 소스 (D10)."""

from __future__ import annotations

import datetime as dt
from pathlib import Path

import pandas as pd

from msl.db import download, raw_dir, record_checksums, with_provenance, write_table

SOURCE = "pubchem-periodic-table"
TABLES = ["elements"]
URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/periodictable/CSV"
LICENSE = "public-domain"
COLUMNS = {  # PubChem 열 → 우리 열
    "AtomicNumber": "Z", "Symbol": "symbol", "Name": "name", "AtomicMass": "atomic_mass",
    "ElectronConfiguration": "electron_configuration", "Electronegativity": "electronegativity_pauling",
    "AtomicRadius": "atomic_radius_pm", "IonizationEnergy": "ionization_energy_ev",
    "ElectronAffinity": "electron_affinity_ev", "OxidationStates": "oxidation_states",
    "StandardState": "standard_state", "MeltingPoint": "melting_point_k", "BoilingPoint": "boiling_point_k",
    "Density": "density_g_cm3", "GroupBlock": "group_block", "YearDiscovered": "year_discovered",
}


def _version() -> str:
    return dt.date.today().isoformat()


def _latest() -> Path | None:
    snaps = sorted((raw_dir(SOURCE, "_").parent).glob("*/periodictable.csv"))
    return snaps[-1] if snaps else None


def fetch() -> Path:
    d = raw_dir(SOURCE, _version())
    download(URL, d / "periodictable.csv")
    record_checksums(d)
    return d


def load() -> dict[str, dict[str, int]]:
    src = _latest() or fetch() / "periodictable.csv"
    df = pd.read_csv(src).rename(columns=COLUMNS)[list(COLUMNS.values())]
    df = with_provenance(df, source=SOURCE, version=src.parent.name, license=LICENSE,
                         method="compiled", id_column="symbol")
    return {"elements": write_table("elements", df)}
