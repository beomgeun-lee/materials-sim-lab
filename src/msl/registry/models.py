"""지식 레지스트리(kb/*.yaml) 항목 모델.

- licenses.yaml : 라이선스 정책 표
- sources.yaml  : 데이터 소스 (계획서 4장)
- engines.yaml  : 계산 엔진·모델 가중치 (계획서 3·5장)
- assays.yaml   : 시험 카탈로그 (계획서 3.2절)
"""

from __future__ import annotations

import datetime as dt
from enum import StrEnum
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from msl.schema.result import Fidelity

SLUG = r"^[a-z0-9][a-z0-9\-]*$"
ASSAY_CODE = r"^(S0|A[1-9][0-9]*)$"


class _Entry(BaseModel):
    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)


# ── 라이선스 ──────────────────────────────────────────────────────────────


class Partition(StrEnum):
    OPEN = "open"  # 배포 빌드에 넣을 수 있다
    RESTRICTED = "restricted"  # 연구용 파티션에만 둔다


class LicenseEntry(_Entry):
    id: str
    name: str
    commercial_use: bool | None
    redistribution: bool | None
    share_alike: bool = False
    no_derivatives: bool = False
    copyleft: Literal["none", "weak", "strong", "network"] = "none"
    url: str | None = None
    note: str | None = None

    @property
    def partition(self) -> Partition:
        if self.commercial_use is True and self.redistribution is True and not self.no_derivatives:
            return Partition.OPEN
        return Partition.RESTRICTED

    @property
    def same_license_required(self) -> bool:
        """파생물에 같은 라이선스를 적용해야 하는가. 데이터에 붙은 GPL 도 SA 와 같이 취급한다."""
        return self.share_alike or self.copyleft in ("strong", "network")


# ── 데이터 소스 ───────────────────────────────────────────────────────────


class Layer(StrEnum):
    """계획서 4.1절 데이터 레이어."""

    E = "E"  # 원소
    C = "C"  # 계산 결정 (ML 학습 데이터 포함)
    X = "X"  # 실험 결정 구조
    M = "M"  # 광물
    Q = "Q"  # 분자·안전 (식별자, GHS, 반응성)
    T = "T"  # 열역학
    G = "G"  # 수급·경제
    K = "K"  # 규제 (국내 중심)


class Priority(StrEnum):
    HIGH = "상"
    MID = "중"
    LOW = "하"
    EXCLUDED = "제외"


class Access(StrEnum):
    REST_API = "rest-api"
    BULK = "bulk"  # 파일·덤프·S3·rsync·FTP 벌크
    OPTIMADE = "optimade"
    PYTHON_PACKAGE = "python-package"
    SPARQL = "sparql"
    WEB = "web"  # 웹 조회만 가능
    DESKTOP_APP = "desktop-app"
    PDF = "pdf"


class Auth(StrEnum):
    NONE = "none"
    API_KEY = "api-key"
    REGISTRATION = "registration"  # 회원가입·로그인
    APPLICATION = "application"  # 신청·승인 필요
    PAID = "paid"


class SourceStatus(StrEnum):
    PLANNED = "planned"
    FETCHED = "fetched"  # 원본 스냅샷 확보
    LOADED = "loaded"  # DB 적재 완료
    EXCLUDED = "excluded"


class SourceEntry(_Entry):
    id: str = Field(pattern=SLUG)
    name: str
    layer: Layer
    operator: str
    url: str
    data: list[str] = Field(min_length=1)
    access: list[Access] = Field(min_length=1)
    auth: Auth
    license: str  # licenses.yaml id
    license_note: str | None = None
    version: str | None = None
    size: str | None = None
    update: str | None = None
    python_client: str | None = None
    priority: Priority
    mvp_rank: int | None = Field(default=None, ge=1, le=10)  # 계획서 4장 MVP 적재 Top 10 순위
    status: SourceStatus = SourceStatus.PLANNED
    country: str | None = Field(default=None, pattern=r"^[A-Z]{2}$|^INT$")
    uses: list[str] = Field(default_factory=list)  # 이 소스를 쓰는 시험 코드
    caveats: str | None = None
    refs: list[str] = Field(default_factory=list)
    checked: dt.date


# ── 엔진·모델 ─────────────────────────────────────────────────────────────


class EngineKind(StrEnum):
    LIBRARY = "library"  # 파이썬 라이브러리·계산 코드
    MODEL = "model"  # 사전학습 모델 가중치 (uMLIP 등)
    SERVICE = "service"  # 호스팅 API
    RULES = "rules"  # 규칙 체계 (CAMEO 반응성 등)
    IN_HOUSE = "in-house"  # 우리가 직접 구현
    COMMERCIAL = "commercial"  # 상용 제품 (교차검증 기준)


class Isolation(StrEnum):
    IN_PROCESS = "in-process"  # 코어에서 import
    WORKER = "worker"  # 격리 워커 컨테이너 (GPL 등)
    EXCLUDED = "excluded"  # 쓰지 않음
    REFERENCE_ONLY = "reference-only"  # 교차검증·참고만


class Use(StrEnum):
    RESEARCH = "research"
    PRODUCT = "product"


class Install(StrEnum):
    PIP = "pip"
    CONDA = "conda"
    SOURCE = "source"
    BINARY = "binary"
    COMMERCIAL = "commercial"
    SERVICE = "service"
    NONE = "none"  # 자체 구현


class Compute(StrEnum):
    LAPTOP = "laptop"
    GPU = "gpu"
    HPC = "hpc"


class Fit(StrEnum):
    HIGH = "상"
    MID = "중"
    LOW = "하"
    REFERENCE = "기준용"


class EngineEntry(_Entry):
    id: str = Field(pattern=SLUG)
    name: str
    kind: EngineKind
    category: str
    assays: list[str] = Field(default_factory=list)
    fidelity: list[Fidelity] = Field(default_factory=list)
    version: str | None = None
    released: dt.date | None = None
    code_license: str  # licenses.yaml id
    weights_license: str | None = None  # 모델 가중치 라이선스 (kind=model 이면 필수)
    isolation: Isolation
    allowed_in: list[Use] = Field(min_length=1)
    install: Install
    python_requires: str | None = None
    compute: list[Compute] = Field(default_factory=list)
    fit: Fit
    role: str | None = None
    caveats: str | None = None
    refs: list[str] = Field(default_factory=list)
    checked: dt.date


# ── 시험 카탈로그 ─────────────────────────────────────────────────────────


class AssayEntry(_Entry):
    code: str = Field(pattern=ASSAY_CODE)
    name: str
    question: str
    applies_to: list[str] = Field(min_length=1)
    fidelity: list[Fidelity] = Field(min_length=1)
    engines: list[str] = Field(default_factory=list)  # engines.yaml id
    sources: list[str] = Field(default_factory=list)  # sources.yaml id
    stage: int = Field(ge=0, le=6)  # 도입 단계 (계획서 6장)
    note: str | None = None
