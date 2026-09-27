"""레지스트리 — 실제 kb/ 통과 여부와 라이선스 정책 규칙."""

from __future__ import annotations

import datetime as dt

import pytest

from msl.registry.load import Registry, check_registry, load_registry
from msl.registry.models import (
    EngineEntry,
    Isolation,
    LicenseEntry,
    Partition,
    SourceEntry,
)

TODAY = dt.date(2026, 9, 26)


def test_real_kb_passes() -> None:
    reg = load_registry()
    assert reg.licenses and reg.assays
    issues = check_registry(reg)
    assert issues == [], "\n".join(issues)


def test_real_kb_mvp_top10_is_open_except_known() -> None:
    """MVP 적재 소스 중 restricted 파티션은 라이선스 확인 대기 중인 것만이어야 한다."""
    reg = load_registry()
    restricted = sorted(
        s.id for s in reg.sources.values() if s.mvp_rank and reg.source_partition(s) is Partition.RESTRICTED
    )
    # KOSHA MSDS: 포털·기관 문구 상충, CAMEO·CIAAW: 재사용 조건 미확인 (계획서 10절)
    # mendeleev: 코드는 MIT 지만 수록값 출처(CRC·위키백과 등)가 섞여 있음 → 연구용 (D10)
    assert set(restricted) <= {"kosha-msds", "cameo", "ciaaw", "mendeleev"}, restricted


def test_real_kb_open_element_layer_exists() -> None:
    """D10: 원소 레이어는 배포 빌드용 open 소스가 따로 있어야 한다."""
    reg = load_registry()
    open_elements = {
        s.id for s in reg.sources.values() if s.layer.value == "E" and reg.source_partition(s) is Partition.OPEN
    }
    assert {"periodictable", "pubchem-periodic-table"} <= open_elements, open_elements
    assert reg.source_partition(reg.sources["wikidata"]) is Partition.OPEN


@pytest.mark.parametrize(
    ("commercial", "redistribution", "nd", "expected"),
    [
        (True, True, False, Partition.OPEN),
        (True, True, True, Partition.RESTRICTED),
        (False, True, False, Partition.RESTRICTED),
        (None, None, False, Partition.RESTRICTED),
    ],
)
def test_partition_rule(commercial: bool | None, redistribution: bool | None, nd: bool, expected: Partition) -> None:
    lic = LicenseEntry(
        id="x", name="x", commercial_use=commercial, redistribution=redistribution, no_derivatives=nd
    )
    assert lic.partition is expected


def _mini_registry(**engine: object) -> Registry:
    licenses = {
        "MIT": LicenseEntry(id="MIT", name="MIT", commercial_use=True, redistribution=True),
        "GPL-3.0": LicenseEntry(id="GPL-3.0", name="GPL", commercial_use=True, redistribution=True, copyleft="strong"),
        "AGPL-3.0": LicenseEntry(
            id="AGPL-3.0", name="AGPL", commercial_use=True, redistribution=True, copyleft="network"
        ),
        "ASL": LicenseEntry(id="ASL", name="ASL", commercial_use=False, redistribution=None),
        "unknown": LicenseEntry(id="unknown", name="미확인", commercial_use=None, redistribution=None),
    }
    fields: dict[str, object] = {
        "id": "eng",
        "name": "eng",
        "kind": "library",
        "category": "시험",
        "code_license": "MIT",
        "isolation": "in-process",
        "allowed_in": ["research", "product"],
        "install": "pip",
        "fit": "상",
        "checked": TODAY,
    } | engine
    return Registry(licenses=licenses, engines={"eng": EngineEntry.model_validate(fields)})


@pytest.mark.parametrize(
    ("engine", "message"),
    [
        ({"code_license": "GPL-3.0"}, "GPL"),
        ({"code_license": "AGPL-3.0", "isolation": "excluded"}, "AGPL"),
        ({"kind": "model", "weights_license": "ASL"}, "비상업"),
        ({"kind": "model"}, "weights_license"),
        ({"code_license": "unknown"}, "확인되지 않은"),
        ({"kind": "commercial"}, "reference-only"),
        ({"code_license": "WTFPL"}, "없는 라이선스"),
        ({"assays": ["A1"]}, "없는 시험 코드"),
    ],
)
def test_engine_policy_violations(engine: dict[str, object], message: str) -> None:
    issues = check_registry(_mini_registry(**engine))
    assert any(message in i for i in issues), issues


@pytest.mark.parametrize(
    "engine",
    [
        {"code_license": "GPL-3.0", "isolation": Isolation.WORKER},
        {"code_license": "AGPL-3.0", "isolation": "excluded", "allowed_in": ["research"]},
        {"kind": "model", "weights_license": "ASL", "allowed_in": ["research"]},
    ],
)
def test_engine_policy_ok(engine: dict[str, object]) -> None:
    assert check_registry(_mini_registry(**engine)) == []


def test_source_mvp_requires_high_priority() -> None:
    reg = _mini_registry()
    reg.sources["s"] = SourceEntry.model_validate(
        {
            "id": "s",
            "name": "s",
            "layer": "E",
            "operator": "x",
            "url": "https://example.org",
            "data": ["x"],
            "access": ["bulk"],
            "auth": "none",
            "license": "MIT",
            "priority": "중",
            "mvp_rank": 1,
            "checked": TODAY,
        }
    )
    issues = check_registry(reg)
    assert any("우선순위 '상'" in i for i in issues)
    assert any("Top 10 순위 누락" in i for i in issues)
