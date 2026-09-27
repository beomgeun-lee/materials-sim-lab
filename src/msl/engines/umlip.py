"""범용 원자간 퍼텐셜(uMLIP) — 구조 이완·에너지 (L2, 계획서 4단계, D5 에 따라 CPU 기준).

허용 가중치만 불러온다: 같은 패키지로 비상업(ASL 등) 가중치도 로드되므로 이름 허용 목록으로 막는다 (engines.yaml).
에너지는 MPtrj 류로 학습한 모델이라 MP 원시 DFT(GGA/GGA+U) 에너지와 같은 기준이다 → MP2020 보정을 붙여
MP 상태도(A1 과 같은 GGA_GGA+U 엔트리)와 바로 비교한다.
"""

from __future__ import annotations

import time
import warnings
from dataclasses import dataclass
from functools import cache
from importlib.metadata import version
from typing import Any

import numpy as np
from pymatgen.core import Structure
from pymatgen.entries.computed_entries import ComputedStructureEntry

MODELS = {  # 이름 → (패키지, 가중치, engines.yaml id). 모두 허용 라이선스(MIT·Apache-2.0)
    "mace-mpa-0": ("mace-torch", "medium-mpa-0", "mace-mpa-0"),
    "orb-v3": ("orb-models", "orb_v3_conservative_inf_mpa", "orb-v3"),
}
MP2020_U = {"Co": 3.32, "Cr": 3.7, "Fe": 5.3, "Mn": 3.9, "Mo": 4.38, "Ni": 6.2, "V": 3.25, "W": 6.2}  # MP GGA+U 값 (O·F 화합물)


class UmlipUnavailable(RuntimeError):
    """선택 설치 묶음 l2 가 없을 때 — `uv sync --extra l2`."""


@cache
def calculator(name: str, device: str = "cpu"):
    """ASE 계산기. 허용 목록 밖의 이름은 거부한다."""
    if name not in MODELS:
        raise ValueError(f"허용되지 않은 uMLIP: {name} (허용: {', '.join(MODELS)})")
    warnings.filterwarnings("ignore")
    import os

    import certifi

    os.environ.setdefault("SSL_CERT_FILE", certifi.where())  # python.org 판 Python 의 urllib 은 인증서 묶음이 없다 (가중치 내려받기)
    try:
        if name == "mace-mpa-0":
            from mace.calculators import mace_mp

            return mace_mp(model=MODELS[name][1], device=device, default_dtype="float64")
        from orb_models.forcefield import pretrained
        from orb_models.forcefield.inference.calculator import ORBCalculator

        orbff, adapter = pretrained.orb_v3_conservative_inf_mpa(  # 0.7.0: (모델, 원자 어댑터)
            device=device, precision="float32-high", compile=False)  # CPU 에서 compile 하면 새 크기마다 20~35초씩 다시 컴파일
        return ORBCalculator(orbff, atoms_adapter=adapter, device=device)
    except ImportError as exc:
        raise UmlipUnavailable(f"uMLIP 패키지가 없음 — `uv sync --extra l2` ({exc})") from exc


def engine_version(name: str) -> str:
    pkg, weights, _ = MODELS[name]
    return f"{pkg} {version(pkg)} / {weights}"


@dataclass
class Relaxed:
    structure: Structure
    energy: float  # eV (셀 전체)
    steps: int
    converged: bool
    seconds: float
    fmax: float  # 마지막 최대 힘 (eV/Å)


def single_point(structure: Structure, name: str) -> tuple[float, float]:
    """(에너지 eV, 걸린 초)."""
    from pymatgen.io.ase import AseAtomsAdaptor

    atoms = AseAtomsAdaptor.get_atoms(structure)
    atoms.calc = calculator(name)
    t = time.perf_counter()
    e = float(atoms.get_potential_energy())
    return e, time.perf_counter() - t


def relax(structure: Structure, name: str, fmax: float = 0.05, steps: int = 300, cell: bool = True) -> Relaxed:
    """원자 위치(+격자) 이완 — ASE FIRE, 격자는 FrechetCellFilter."""
    from ase.filters import FrechetCellFilter
    from ase.optimize import FIRE
    from pymatgen.io.ase import AseAtomsAdaptor

    atoms = AseAtomsAdaptor.get_atoms(structure)
    atoms.calc = calculator(name)
    target = FrechetCellFilter(atoms) if cell else atoms
    opt = FIRE(target, logfile=None)
    t = time.perf_counter()
    converged = bool(opt.run(fmax=fmax, steps=steps))
    secs = time.perf_counter() - t
    forces = atoms.get_forces()
    return Relaxed(AseAtomsAdaptor.get_structure(atoms), float(atoms.get_potential_energy()), int(opt.get_number_of_steps()),
                   converged, secs, float(np.linalg.norm(forces, axis=1).max()))


def mp_entry(structure: Structure, energy: float, entry_id: str = "umlip") -> ComputedStructureEntry | None:
    """uMLIP 에너지 → MP2020 보정을 붙인 엔트리 (MP 상태도와 같은 기준). 보정할 수 없으면 None."""
    from pymatgen.entries.compatibility import MaterialsProject2020Compatibility

    els = {str(e) for e in structure.composition.elements}
    is_u = bool(els & {"O", "F"}) and bool(els & set(MP2020_U))
    params: dict[str, Any] = {"run_type": "GGA+U" if is_u else "GGA", "is_hubbard": is_u,
                              "hubbards": {e: MP2020_U.get(e, 0.0) for e in els} if is_u else {}}
    entry = ComputedStructureEntry(structure, energy, parameters=params, entry_id=entry_id)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return MaterialsProject2020Compatibility(check_potcar=False).process_entry(entry)
