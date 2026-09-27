"""데이터 저장 계층 — 원본 스냅샷, 프로비넌스 검사, 라이선스 파티션, DuckDB 조회 (계획서 4.2·4.3절).

배치
    data/raw/{source}/{version}/...         원본 그대로 + SHA256SUMS
    data/open/{table}/{source}.parquet      재배포 가능 (배포 빌드에 들어감)
    data/restricted/{table}/{source}.parquet 비상업·SA·허가 필요·미확인 (연구용)

규칙
    - 모든 행은 PROVENANCE 열을 가진다. 비어 있으면 적재를 거부한다.
    - 행의 license 로 파티션을 정한다 (kb/licenses.yaml 정책). 사람이 고르지 않는다.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import shutil
import urllib.request
from dataclasses import dataclass
from functools import cache
from pathlib import Path

import duckdb
import pandas as pd

from msl.env import DATA_DIR
from msl.registry.load import load_registry
from msl.registry.models import Partition

RAW = DATA_DIR / "raw"
PROVENANCE = ["source", "source_id", "source_version", "retrieved_at", "license", "method", "correction_scheme"]
REQUIRED = [c for c in PROVENANCE if c != "correction_scheme"]


class ProvenanceError(ValueError):
    """프로비넌스 열이 없거나 비어 있을 때."""


# ── 원본 스냅샷 ───────────────────────────────────────────────────────────


def raw_dir(source: str, version: str) -> Path:
    d = RAW / source / version.replace("/", "-")
    d.mkdir(parents=True, exist_ok=True)
    return d


def record_checksums(folder: Path) -> Path:
    """폴더 안 파일들의 SHA256 을 SHA256SUMS 에 기록한다 (shasum -c 로 검증 가능)."""
    lines = []
    for f in sorted(p for p in folder.rglob("*") if p.is_file() and p.name != "SHA256SUMS"):
        h = hashlib.sha256(f.read_bytes()).hexdigest()
        lines.append(f"{h}  {f.relative_to(folder)}")
    out = folder / "SHA256SUMS"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def download(url: str, dest: Path, *, overwrite: bool = False, headers: dict[str, str] | None = None) -> Path:
    """URL 을 dest 로 내려받는다. 이미 있으면 건너뛴다 (스냅샷은 버전 폴더로 구분)."""
    if dest.exists() and not overwrite:
        return dest
    import ssl

    import certifi

    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "materials-sim-lab/0.0.1 (+research)", **(headers or {})})
    ctx = ssl.create_default_context(cafile=certifi.where())
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(req, timeout=120, context=ctx) as resp, tmp.open("wb") as fh:
        shutil.copyfileobj(resp, fh)
    tmp.rename(dest)
    return dest


# ── 적재 ──────────────────────────────────────────────────────────────────


def now() -> str:
    return dt.datetime.now(dt.UTC).isoformat(timespec="seconds")


def with_provenance(df: pd.DataFrame, *, source: str, version: str, license: str, method: str,
                    id_column: str, retrieved_at: str | None = None,
                    correction_scheme: str | None = None) -> pd.DataFrame:
    """행마다 같은 프로비넌스를 붙인다. source_id 는 id_column 값."""
    out = df.copy()
    out["source"] = source
    out["source_id"] = out[id_column].astype(str)
    out["source_version"] = version
    out["retrieved_at"] = retrieved_at or now()
    out["license"] = license
    out["method"] = method
    out["correction_scheme"] = correction_scheme
    return out


@cache
def _partition_of(license_id: str) -> Partition:
    lic = load_registry().licenses.get(license_id)
    return lic.partition if lic else Partition.RESTRICTED


def write_table(table: str, df: pd.DataFrame) -> dict[str, int]:
    """프로비넌스를 검사하고 라이선스별로 파티션을 나눠 Parquet 으로 쓴다. 반환: {파티션: 행 수}."""
    missing = [c for c in PROVENANCE if c not in df.columns]
    if missing:
        raise ProvenanceError(f"{table}: 프로비넌스 열 없음 {missing}")
    empty = [c for c in REQUIRED if df[c].isna().any() or (df[c].astype(str).str.strip() == "").any()]
    if empty:
        raise ProvenanceError(f"{table}: 프로비넌스 값이 빈 행이 있음 {empty}")
    sources = df["source"].unique()
    if len(sources) != 1:
        raise ProvenanceError(f"{table}: 한 번에 한 소스만 적재 ({list(sources)})")
    source = str(sources[0])
    counts: dict[str, int] = {}
    for part in Partition:
        path = DATA_DIR / part.value / table / f"{source}.parquet"
        keep = {lid for lid in df["license"].unique() if _partition_of(str(lid)) is part}
        rows = df[df["license"].isin(keep)]
        if rows.empty:
            path.unlink(missing_ok=True)
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        rows.to_parquet(path, index=False)
        counts[part.value] = len(rows)
    return counts


# ── 조회 ──────────────────────────────────────────────────────────────────


def tables() -> list[str]:
    names = {p.name for part in Partition for p in (DATA_DIR / part.value).glob("*") if p.is_dir()}
    return sorted(names)


def connect(include_restricted: bool = True) -> duckdb.DuckDBPyConnection:
    """테이블마다 뷰를 만든 DuckDB 연결. include_restricted=False 면 배포용 파티션만 보인다."""
    con = duckdb.connect()
    parts = [Partition.OPEN] + ([Partition.RESTRICTED] if include_restricted else [])
    for t in tables():
        files = [str(f) for part in parts for f in (DATA_DIR / part.value / t).glob("*.parquet")]
        if files:
            listing = ", ".join(f"'{f}'" for f in files)
            con.execute(f"CREATE VIEW {t} AS SELECT * FROM read_parquet([{listing}], union_by_name=true)")
    return con


@dataclass
class TableStatus:
    table: str
    source: str
    partition: str
    rows: int
    version: str
    retrieved_at: str


def status() -> list[TableStatus]:
    out = []
    for part in Partition:
        for f in sorted((DATA_DIR / part.value).glob("*/*.parquet")):
            meta = duckdb.sql(
                f"SELECT count(*), any_value(source_version), max(retrieved_at) FROM read_parquet('{f}')"
            ).fetchone()
            out.append(TableStatus(f.parent.name, f.stem, part.value, meta[0], str(meta[1]), str(meta[2])))
    return out
