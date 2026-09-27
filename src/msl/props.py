"""L2 물성 — uMLIP 으로 탄성률·포논 안정성 (계획서 4단계 A7 L2).

- 탄성: 이완한 구조에 6개 방향 변형(±0.5·±1%)을 주고 원자만 다시 이완해 응력을 잰 뒤, 기울기로 탄성 텐서 C(GPa) 를 맞춘다.
  체적·전단 탄성률은 Voigt–Reuss–Hill 평균 (pymatgen ElasticTensor). 응력은 ASE 관례(인장 +).
- 포논: 한 변 ≥ 10 Å 초격자에서 원자 변위 0.01 Å 의 힘으로 phonopy 힘 상수를 만들고 격자(mesh) 진동수를 본다.
  허수 진동수(음수로 표기)가 −0.3 THz 보다 작으면 동역학적으로 불안정 (Γ 근처 음향 모드의 수치 잡음은 봐준다).
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np
from pymatgen.core import Structure

from msl.engines import umlip

EV_A3_TO_GPA = 160.21766
STRAINS = (-0.01, -0.005, 0.005, 0.01)
IMAG_TOL_THZ = -0.3


def _atoms(structure: Structure, model: str):
    from pymatgen.io.ase import AseAtomsAdaptor

    a = AseAtomsAdaptor.get_atoms(structure)
    a.calc = umlip.calculator(model)
    return a


def _stress_after_ion_relax(structure: Structure, model: str) -> np.ndarray:
    """셀 고정, 원자만 이완한 뒤 응력 (Voigt xx,yy,zz,yz,xz,xy, GPa)."""
    from ase.optimize import FIRE

    a = _atoms(structure, model)
    if len(a) > 1:
        FIRE(a, logfile=None).run(fmax=0.01, steps=300)
    return np.asarray(a.get_stress(voigt=True)) * EV_A3_TO_GPA


def elastic(structure: Structure, model: str, relax: bool = True) -> dict[str, Any]:
    from pymatgen.analysis.elasticity import ElasticTensor

    t0 = time.perf_counter()
    s0 = umlip.relax(structure, model, fmax=0.005, steps=500).structure if relax else structure
    C = np.zeros((6, 6))
    for j in range(6):
        rows = []
        for d in STRAINS:
            eps = np.zeros((3, 3))
            if j < 3:
                eps[j, j] = d
            else:  # 공학 전단 변형 γ = d → 텐서 성분 d/2
                a, b = {3: (1, 2), 4: (0, 2), 5: (0, 1)}[j]
                eps[a, b] = eps[b, a] = d / 2
            s = s0.copy()
            s.lattice = s.lattice.__class__(s0.lattice.matrix @ (np.eye(3) + eps))
            rows.append(_stress_after_ion_relax(s, model))
        slope = np.polyfit(np.array(STRAINS), np.array(rows), 1)[0]  # (6,) 각 응력 성분의 기울기
        C[:, j] = slope
    C = 0.5 * (C + C.T)
    et = ElasticTensor.from_voigt(C)
    eig = np.linalg.eigvalsh(C)
    return {"C_GPa": np.round(C, 1).tolist(), "K_vrh": round(float(et.k_vrh), 1), "G_vrh": round(float(et.g_vrh), 1),
            "mechanically_stable": bool(eig.min() > 0), "min_eigen_GPa": round(float(eig.min()), 1),
            "seconds": round(time.perf_counter() - t0, 1), "engine": umlip.engine_version(model)}


def phonon(structure: Structure, model: str, min_length: float = 10.0, mesh: int = 12, relax: bool = True) -> dict[str, Any]:
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms

    t0 = time.perf_counter()
    s0 = umlip.relax(structure, model, fmax=0.005, steps=500).structure if relax else structure
    unit = PhonopyAtoms(symbols=[str(x.specie) for x in s0], cell=s0.lattice.matrix, scaled_positions=s0.frac_coords)
    sc = np.diag([max(1, int(np.ceil(min_length / n))) for n in s0.lattice.abc])
    ph = Phonopy(unit, supercell_matrix=sc)
    ph.generate_displacements(distance=0.01)
    calc = umlip.calculator(model)
    from ase import Atoms

    forces = []
    for cell in ph.supercells_with_displacements:
        a = Atoms(cell.symbols, cell=cell.cell, scaled_positions=cell.scaled_positions, pbc=True)
        a.calc = calc
        forces.append(a.get_forces())
    ph.forces = np.array(forces)
    ph.produce_force_constants()
    ph.symmetrize_force_constants()  # 음향 합 규칙·대칭 — Γ 근처 가짜 허수 모드를 없앤다
    ph.run_mesh([mesh] * 3, is_gamma_center=True)
    freqs = ph.get_mesh_dict()["frequencies"]  # THz, 허수는 음수
    fmin = float(freqs.min())
    return {"min_freq_THz": round(fmin, 3), "max_freq_THz": round(float(freqs.max()), 2),
            "imaginary_fraction": round(float((freqs < IMAG_TOL_THZ).mean()), 4), "dynamically_stable": fmin >= IMAG_TOL_THZ,
            "supercell": [int(x) for x in np.diag(sc)], "n_displacements": len(forces), "n_atoms_supercell": len(forces[0]),
            "seconds": round(time.perf_counter() - t0, 1), "engine": umlip.engine_version(model)}
