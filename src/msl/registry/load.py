"""kb/*.yaml 로드와 교차 검사.

항목 하나하나의 형식은 pydantic 모델이 검사하고, 여기서는 파일 사이의 규칙을 검사한다.
- 참조 무결성: 라이선스·엔진·소스·시험 id 가 실제로 있는가
- 라이선스 정책: GPL 은 격리, AGPL·비상업은 제품 제외
- 계획서와의 일치: MVP 적재 Top 10 순위가 모두 채워졌는가
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, TypeVar

import yaml
from pydantic import BaseModel, ValidationError

from msl.env import KB_DIR
from msl.registry.models import (
    AssayEntry,
    EngineEntry,
    EngineKind,
    Isolation,
    LicenseEntry,
    Partition,
    SourceEntry,
    SourceStatus,
    Priority,
    Use,
)

M = TypeVar("M", bound=BaseModel)



@dataclass
class Registry:
    licenses: dict[str, LicenseEntry] = field(default_factory=dict)
    sources: dict[str, SourceEntry] = field(default_factory=dict)
    engines: dict[str, EngineEntry] = field(default_factory=dict)
    assays: dict[str, AssayEntry] = field(default_factory=dict)
    errors: list[str] = field(default_factory=list)  # 형식 오류 (로드 단계)

    def source_partition(self, source: SourceEntry) -> Partition:
        lic = self.licenses.get(source.license)
        return lic.partition if lic else Partition.RESTRICTED


def _load_list(path: Path, model: type[M], key: str, errors: list[str]) -> dict[str, M]:
    if not path.exists():
        return {}
    raw: Any = yaml.safe_load(path.read_text(encoding="utf-8")) or []
    if not isinstance(raw, list):
        errors.append(f"{path.name}: 최상위는 목록이어야 함")
        return {}
    out: dict[str, M] = {}
    for i, item in enumerate(raw):
        label = item.get(key, f"#{i}") if isinstance(item, dict) else f"#{i}"
        try:
            entry = model.model_validate(item)
        except ValidationError as exc:
            for err in exc.errors():
                loc = ".".join(str(p) for p in err["loc"])
                errors.append(f"{path.name}[{label}].{loc}: {err['msg']}")
            continue
        entry_key = getattr(entry, key)
        if entry_key in out:
            errors.append(f"{path.name}[{entry_key}]: id 중복")
        out[entry_key] = entry
    return out


def load_registry(kb_dir: Path = KB_DIR) -> Registry:
    reg = Registry()
    reg.licenses = _load_list(kb_dir / "licenses.yaml", LicenseEntry, "id", reg.errors)
    reg.sources = _load_list(kb_dir / "sources.yaml", SourceEntry, "id", reg.errors)
    reg.engines = _load_list(kb_dir / "engines.yaml", EngineEntry, "id", reg.errors)
    reg.assays = _load_list(kb_dir / "assays.yaml", AssayEntry, "code", reg.errors)
    return reg


def check_registry(reg: Registry) -> list[str]:
    """교차 검사. 문제 목록을 돌려준다 (빈 목록이면 통과)."""
    issues = list(reg.errors)
    lic = reg.licenses

    for s in reg.sources.values():
        where = f"sources.yaml[{s.id}]"
        if s.license not in lic:
            issues.append(f"{where}: 없는 라이선스 id {s.license!r}")
        for code in s.uses:
            if code not in reg.assays:
                issues.append(f"{where}.uses: 없는 시험 코드 {code!r}")
        if s.priority is Priority.EXCLUDED and s.status is not SourceStatus.EXCLUDED:
            issues.append(f"{where}: 우선순위 '제외'면 status 도 excluded 여야 함")
        if s.mvp_rank is not None and s.priority is not Priority.HIGH:
            issues.append(f"{where}: MVP 적재 소스는 우선순위 '상'이어야 함")

    ranks = {s.mvp_rank for s in reg.sources.values() if s.mvp_rank is not None}
    missing = sorted(set(range(1, 11)) - ranks)
    if reg.sources and missing:
        issues.append(f"sources.yaml: MVP 적재 Top 10 순위 누락 {missing}")

    for e in reg.engines.values():
        where = f"engines.yaml[{e.id}]"
        licenses = [("code_license", e.code_license)]
        if e.weights_license is not None:
            licenses.append(("weights_license", e.weights_license))
        elif e.kind is EngineKind.MODEL:
            issues.append(f"{where}: 모델(kind=model)은 weights_license 가 필요함")
        for fld, lid in licenses:
            entry = lic.get(lid)
            if entry is None:
                issues.append(f"{where}.{fld}: 없는 라이선스 id {lid!r}")
                continue
            if entry.copyleft == "strong" and e.isolation is Isolation.IN_PROCESS:
                issues.append(f"{where}: GPL 계열({lid})은 in-process 로 쓸 수 없음 → worker")
            if entry.copyleft == "network" and Use.PRODUCT in e.allowed_in:
                issues.append(f"{where}: AGPL({lid})은 제품(product)에서 쓸 수 없음")
            if entry.commercial_use is False and Use.PRODUCT in e.allowed_in:
                issues.append(f"{where}: 비상업 라이선스({lid})는 제품(product)에서 쓸 수 없음")
            if entry.commercial_use is None and Use.PRODUCT in e.allowed_in:
                issues.append(f"{where}: 상업 이용이 확인되지 않은 라이선스({lid})는 제품(product)에서 쓸 수 없음 (D2)")
        if e.kind is EngineKind.COMMERCIAL and e.isolation is not Isolation.REFERENCE_ONLY:
            issues.append(f"{where}: 상용 제품은 isolation=reference-only 여야 함")
        for code in e.assays:
            if code not in reg.assays:
                issues.append(f"{where}.assays: 없는 시험 코드 {code!r}")

    for a in reg.assays.values():
        where = f"assays.yaml[{a.code}]"
        for eid in a.engines:
            if eid not in reg.engines:
                issues.append(f"{where}.engines: 없는 엔진 id {eid!r}")
            elif reg.engines[eid].isolation is Isolation.EXCLUDED:
                issues.append(f"{where}.engines: 제외된 엔진 {eid!r} 을 참조함")
        for sid in a.sources:
            if sid not in reg.sources:
                issues.append(f"{where}.sources: 없는 소스 id {sid!r}")
    return issues


def summarize(reg: Registry) -> dict[str, Counter[str]]:
    """CLI 요약용 집계."""
    return {
        "소스·레이어": Counter(s.layer.value for s in reg.sources.values()),
        "소스·우선순위": Counter(s.priority.value for s in reg.sources.values()),
        "소스·파티션": Counter(reg.source_partition(s).value for s in reg.sources.values()),
        "소스·동일조건(SA·GPL)": Counter(
            "적용" for s in reg.sources.values() if s.license in reg.licenses and reg.licenses[s.license].same_license_required
        ),
        "엔진·격리": Counter(e.isolation.value for e in reg.engines.values()),
        "엔진·종류": Counter(e.kind.value for e in reg.engines.values()),
    }
