"""A8 블렌드·복합 물성 — 혼합법칙과 Voigt·Reuss·Hashin–Shtrikman 경계 (자체 구현, L0)."""

from __future__ import annotations

from msl.assays.base import Context, Outcome, not_applicable, val
from msl.schema.quantity import Dimension
from msl.schema.result import Fidelity, Status


def _kg(E: float, nu: float) -> tuple[float, float]:
    """영률·푸아송비 → 체적탄성률 K, 전단탄성률 G."""
    return E / (3 * (1 - 2 * nu)), E / (2 * (1 + nu))


def _e(K: float, G: float) -> float:
    return 9 * K * G / (3 * K + G)


def _hs(phi: list[float], K: list[float], G: list[float]) -> tuple[float, float]:
    """다상 Hashin–Shtrikman 경계로 구한 영률 (하한, 상한)."""
    def bound(Kr: float, Gr: float) -> float:
        k = 1 / sum(f / (k_ + 4 / 3 * Gr) for f, k_ in zip(phi, K)) - 4 / 3 * Gr
        z = Gr / 6 * (9 * Kr + 8 * Gr) / (Kr + 2 * Gr)
        g = 1 / sum(f / (g_ + z) for f, g_ in zip(phi, G)) - z
        return _e(k, g)

    return bound(min(K), min(G)), bound(max(K), max(G))


def a8(ctx: Context) -> Outcome:
    comps = ctx.comps
    if not all(c.props for c in comps):
        return not_applicable("사용자 정의 소재(material:)끼리의 블렌드에만 적용", "msl-mixing-rules")
    need = ("density", "youngs_modulus", "poisson_ratio")
    lacking = [c.name for c in comps if any(c.props.get(k) is None for k in need)]
    if lacking:
        return not_applicable(f"물성 부족 ({', '.join(need)}): {', '.join(lacking)}", "msl-mixing-rules")
    dims = {c.amount.dimension if c.amount else None for c in comps}
    if dims == {Dimension.VOLUME_FRACTION}:
        phi = [c.amount.to_si() for c in comps]
    elif dims == {Dimension.MASS_FRACTION}:
        vol = [c.amount.to_si() / c.props["density"] for c in comps]
        phi = [v / sum(vol) for v in vol]
    else:
        return not_applicable("블렌드 성분 양은 모두 vol% 또는 wt% 로 적어야 함", "msl-mixing-rules")

    rho = [c.props["density"] for c in comps]
    E = [c.props["youngs_modulus"] for c in comps]
    nu = [c.props["poisson_ratio"] for c in comps]
    K, G = zip(*(_kg(e, n) for e, n in zip(E, nu)))
    density = sum(f * r for f, r in zip(phi, rho))
    voigt = sum(f * e for f, e in zip(phi, E))
    reuss = 1 / sum(f / e for f, e in zip(phi, E))
    hs_lo, hs_hi = _hs(phi, list(K), list(G))
    return Outcome(
        status=Status.OK, fidelity=Fidelity.L0, engine="msl-mixing-rules", engine_version="v0",
        conditions_basis="선형 탄성, 등방성 성분, 완전 접착 계면",
        values=[
            val("밀도 (혼합법칙)", round(density, 3), "g/cm³"),
            val("영률 상한 · Voigt (섬유 방향)", round(voigt, 2), "GPa"),
            val("영률 하한 · Reuss (섬유 수직)", round(reuss, 2), "GPa"),
            val("영률 · Hashin–Shtrikman 하한", round(hs_lo, 2), "GPa"),
            val("영률 · Hashin–Shtrikman 상한", round(hs_hi, 2), "GPa"),
        ],
        summary=f"영률은 배열에 따라 {reuss:.1f}~{voigt:.1f} GPa, 무작위 등방 배열이면 {hs_lo:.1f}~{hs_hi:.1f} GPa",
        caveats=["미세구조(섬유 길이·배향·계면)를 넣지 않으면 경계값만 제시", "입력 물성은 kb/materials.yaml 대표값"],
        data={"bounds": {
            "unit": "GPa", "voigt": voigt, "reuss": reuss, "hs": [hs_lo, hs_hi],
            "components": [{"name": c.name, "value": e, "fraction": f} for c, e, f in zip(comps, E, phi)],
        }},
    )
