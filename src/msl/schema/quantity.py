"""수량 파싱과 단위 정규화.

레시피에는 사람이 쓰는 문자열("1,000 °C", "0.1 mol/L", "30 vol%")이 들어온다.
v0은 조합 시험에 필요한 단위만 허용하고, 그 밖의 단위는 거부한다.
"""

from __future__ import annotations

import re
from enum import StrEnum
from typing import Any

from pydantic import BaseModel, ConfigDict, model_validator


class Dimension(StrEnum):
    TEMPERATURE = "temperature"
    PRESSURE = "pressure"
    AMOUNT = "amount"  # 물질량
    MASS = "mass"
    VOLUME = "volume"
    CONCENTRATION = "concentration"  # 몰농도
    MOLALITY = "molality"
    MASS_FRACTION = "mass_fraction"
    VOLUME_FRACTION = "volume_fraction"
    MOLE_FRACTION = "mole_fraction"
    TIME = "time"


D = Dimension

# 단위 → (차원, SI 기준 환산 계수). 온도는 오프셋이 있어 to_si()에서 따로 처리한다.
UNITS: dict[str, tuple[Dimension, float]] = {
    "K": (D.TEMPERATURE, 1.0),
    "°C": (D.TEMPERATURE, 1.0),
    "degC": (D.TEMPERATURE, 1.0),
    "Pa": (D.PRESSURE, 1.0),
    "kPa": (D.PRESSURE, 1e3),
    "MPa": (D.PRESSURE, 1e6),
    "GPa": (D.PRESSURE, 1e9),
    "bar": (D.PRESSURE, 1e5),
    "atm": (D.PRESSURE, 101_325.0),
    "mol": (D.AMOUNT, 1.0),
    "mmol": (D.AMOUNT, 1e-3),
    "kmol": (D.AMOUNT, 1e3),
    "kg": (D.MASS, 1.0),
    "g": (D.MASS, 1e-3),
    "mg": (D.MASS, 1e-6),
    "t": (D.MASS, 1e3),
    "m3": (D.VOLUME, 1.0),
    "L": (D.VOLUME, 1e-3),
    "mL": (D.VOLUME, 1e-6),
    "mol/L": (D.CONCENTRATION, 1e3),
    "M": (D.CONCENTRATION, 1e3),
    "mmol/L": (D.CONCENTRATION, 1.0),
    "mM": (D.CONCENTRATION, 1.0),
    "mol/kg": (D.MOLALITY, 1.0),
    "wt%": (D.MASS_FRACTION, 1e-2),
    "vol%": (D.VOLUME_FRACTION, 1e-2),
    "mol%": (D.MOLE_FRACTION, 1e-2),
    "at%": (D.MOLE_FRACTION, 1e-2),
    "s": (D.TIME, 1.0),
    "min": (D.TIME, 60.0),
    "h": (D.TIME, 3600.0),
}

FRACTION_DIMENSIONS = frozenset({D.MASS_FRACTION, D.VOLUME_FRACTION, D.MOLE_FRACTION})

_QUANTITY_RE = re.compile(r"^\s*([+-]?(?:\d[\d,]*)?\.?\d+(?:[eE][+-]?\d+)?)\s*(\S+)\s*$")


class Quantity(BaseModel):
    """값 + 단위. `"1,000 °C"` 같은 문자열이나 `{value, unit}` 매핑으로 만든다."""

    model_config = ConfigDict(frozen=True)

    value: float
    unit: str

    @model_validator(mode="before")
    @classmethod
    def _parse(cls, data: Any) -> Any:
        if isinstance(data, str):
            m = _QUANTITY_RE.match(data)
            if not m:
                raise ValueError(f"수량 형식이 아님: {data!r} (예: '1273 K', '0.1 mol/L', '30 vol%')")
            return {"value": float(m.group(1).replace(",", "")), "unit": m.group(2)}
        return data

    @model_validator(mode="after")
    def _known_unit(self) -> Quantity:
        if self.unit not in UNITS:
            raise ValueError(f"지원하지 않는 단위: {self.unit!r} (허용: {', '.join(UNITS)})")
        return self

    @property
    def dimension(self) -> Dimension:
        return UNITS[self.unit][0]

    def to_si(self) -> float:
        """SI 기준값. 온도는 K, 압력은 Pa, 분율은 0~1."""
        if self.unit in ("°C", "degC"):
            return self.value + 273.15
        return self.value * UNITS[self.unit][1]

    def __str__(self) -> str:
        return f"{self.value:g} {self.unit}"
