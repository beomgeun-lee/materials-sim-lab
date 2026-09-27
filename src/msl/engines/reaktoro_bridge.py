"""Reaktoro 연결 — conda(micromamba) 전용 환경의 Python 을 별도 프로세스로 불러 수계 평형을 푼다 (A4 이중화, D22).

Reaktoro(LGPL-2.1)는 conda-forge 로만 배포돼 프로젝트 가상환경에 넣지 않는다. 환경 위치는 MSL_REAKTORO_PYTHON 또는
~/micromamba/envs/reaktoro/bin/python. 없으면 available() 이 False 이고 A4 는 교차검증을 건너뛴다.

DB 는 A4 와 같은 phreeqc.dat 만 쓴다 (USGS 공공 — 제품 경로 허용). Reaktoro 동봉 SUPCRTBL 은 CC BY-NC-ND 라 제외(D2·D22),
SUPCRT98/07 은 라이선스 확인 전이라 연결하지 않는다.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
from functools import cache
from pathlib import Path
from typing import Any

WORKER = r"""
import json, sys, reaktoro as rkt
spec = json.load(sys.stdin)
db = rkt.PhreeqcDatabase.fromFile(spec["db_path"])
aq = rkt.AqueousPhase(rkt.speciate(" ".join(spec["elements"])))
aq.set(rkt.ActivityModelPhreeqc(db))
phases = [aq] + ([rkt.MineralPhases(" ".join(spec["minerals"]))] if spec["minerals"] else [])
system = rkt.ChemicalSystem(db, *phases)
st = rkt.ChemicalState(system)
st.temperature(spec["T_K"], "K"); st.pressure(spec["P_bar"], "bar"); st.set("H2O", spec["water_kg"], "kg")
for name, n in spec["species"].items():
    st.set(name, n, "mol")
for name, n in spec["minerals"].items():
    st.set(name, n, "mol")
res = rkt.EquilibriumSolver(system).solve(st)
p = rkt.AqueousProps(st)
out = {"ok": bool(res.succeeded()), "pH": float(p.pH()), "I": float(p.ionicStrength()),
       "elements": {e: float(p.elementMolality(e)) for e in spec["elements"] if e not in ("H", "O")},
       "si": {m: float(p.saturationIndex(m)) for m in spec["minerals"]}}
json.dump(out, sys.stdout)
"""


def python() -> Path:
    return Path(os.environ.get("MSL_REAKTORO_PYTHON", Path.home() / "micromamba/envs/reaktoro/bin/python"))


@cache
def available() -> bool:
    return os.environ.get("MSL_REAKTORO", "1") != "0" and python().exists()


@cache
def version() -> str:
    r = subprocess.run([str(python()), "-c", "import reaktoro; print(reaktoro.__version__)"], capture_output=True, text=True, timeout=60)
    return r.stdout.strip() or "?"


@cache
def _masters(db_path: str) -> dict[str, str]:
    """phreeqc.dat SOLUTION_MASTER_SPECIES: 마스터 이름(Na, C(4) …) → 화학종 (Na+, CO3-2 …)."""
    txt = Path(db_path).read_text(encoding="latin-1")
    block = txt.split("SOLUTION_MASTER_SPECIES", 1)[1].split("SOLUTION_SPECIES", 1)[0]
    out = {}
    for line in block.splitlines():
        parts = line.split("#", 1)[0].split()
        if len(parts) >= 2:
            out[parts[0]] = parts[1]
    return out


def _charge(species: str) -> int:
    m = re.search(r"([+-])(\d*)$", species)
    return 0 if not m else (1 if m.group(1) == "+" else -1) * int(m.group(2) or 1)


def _element(master: str) -> str:
    return re.match(r"[A-Z][a-z]?", master).group(0)


def from_a4(totals: dict[str, float], minerals: list[str], db_path: str, liters: float, T_K: float, P_bar: float) -> dict[str, Any]:
    """A4 입력(PHREEQC 마스터 이름별 mol/kgw, 평형 광물) → Reaktoro 입력. 전하는 H⁺/OH⁻ 로 맞춘다 (PHREEQC 'pH charge' 와 같은 뜻).
    광물은 PHREEQC equalize 기본값과 같이 10 mol 을 넣어 과잉 고체로 둔다. 광물 원소도 원소 목록에 넣어야 한다."""
    masters = _masters(db_path)
    species: dict[str, float] = {}
    els = {"H", "O"}
    for name, molal in totals.items():
        sp = masters.get(name)
        if sp is None:
            raise KeyError(f"마스터 화학종을 모름: {name}")
        species[sp] = species.get(sp, 0.0) + molal * liters
        els.add(_element(name))
    q = sum(_charge(s) * n for s, n in species.items())
    if abs(q) > 1e-12:
        comp = "OH-" if q > 0 else "H+"
        species[comp] = species.get(comp, 0.0) + abs(q)
    for m in minerals:
        els.update(re.findall(r"[A-Z][a-z]?", _mineral_formula(db_path, m)))
    return {"db_path": db_path, "elements": sorted(els), "species": species, "minerals": {m: 10.0 for m in minerals},
            "T_K": T_K, "P_bar": P_bar, "water_kg": liters}


@cache
def _mineral_formula(db_path: str, name: str) -> str:
    txt = Path(db_path).read_text(encoding="latin-1").split("\nPHASES", 1)[1]
    m = re.search(rf"^{re.escape(name)}\s*\n\s+([^=\n]+?)\s*=", txt, re.M)
    if m is None:
        raise KeyError(f"광물을 모름: {name}")
    return m.group(1).split("+")[0]


def solve(spec: dict[str, Any], timeout: int = 120) -> dict[str, Any]:
    r = subprocess.run([str(python()), "-c", WORKER], input=json.dumps(spec), capture_output=True, text=True, timeout=timeout)
    if r.returncode != 0:
        raise RuntimeError(r.stderr.strip().splitlines()[-1] if r.stderr.strip() else "Reaktoro 실패")
    return json.loads(r.stdout)
