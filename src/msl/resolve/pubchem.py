"""PubChem 조회 — 식별자·화학식·CAMEO 반응성 그룹 (키 불필요, 결과는 디스크 캐시)."""

from __future__ import annotations

import json
import re
import ssl
import time
import urllib.error
import urllib.parse
import urllib.request
from functools import cache
from pathlib import Path
from typing import Any

import certifi

from msl.env import CACHE_DIR

BASE = "https://pubchem.ncbi.nlm.nih.gov/rest"
CACHE = CACHE_DIR / "pubchem"
_MIN_INTERVAL = 0.25  # PubChem 권장 한도(초당 5건) 안쪽
_last_call = 0.0
# requests(urllib3)로 보내면 PubChem 이 503 을 돌려주는 현상이 있어(2026-09-26 확인) 표준 urllib 을 쓴다
_SSL = ssl.create_default_context(cafile=certifi.where())
_HEADERS = {"User-Agent": "materials-sim-lab/0.0.1 (+research)"}


def _fetch(url: str) -> tuple[int, bytes]:
    req = urllib.request.Request(url, headers=_HEADERS)
    try:
        with urllib.request.urlopen(req, timeout=20, context=_SSL) as resp:
            return resp.status, resp.read()
    except urllib.error.HTTPError as exc:
        return exc.code, exc.read()


def _get(url: str) -> dict[str, Any] | None:
    """GET + 캐시. 404(없음)도 캐시해 같은 질의를 반복하지 않는다."""
    global _last_call
    CACHE.mkdir(parents=True, exist_ok=True)
    key = Path(CACHE / (urllib.parse.quote(url.removeprefix(BASE), safe="") + ".json"))
    if key.exists():
        cached = json.loads(key.read_text(encoding="utf-8"))
        return cached or None
    for attempt in range(4):  # 503 ServerBusy·429 는 잠시 뒤 다시 시도한다
        wait = _MIN_INTERVAL - (time.monotonic() - _last_call)
        if wait > 0:
            time.sleep(wait)
        _last_call = time.monotonic()
        status, body = _fetch(url)
        if status not in (429, 503) or attempt == 3:
            break
        time.sleep(1.5 * 2**attempt)
    if status == 404:
        key.write_text("{}", encoding="utf-8")
        return None
    if status != 200:
        raise RuntimeError(f"PubChem HTTP {status}: {body[:120]!r}")
    data = json.loads(body)
    key.write_text(json.dumps(data, ensure_ascii=False), encoding="utf-8")
    return data


def lookup(name_or_cas: str) -> dict[str, Any] | None:
    """이름·CAS 번호로 첫 번째 CID 와 화학식·대표명을 찾는다."""
    q = urllib.parse.quote(name_or_cas)
    data = _get(f"{BASE}/pug/compound/name/{q}/property/MolecularFormula,Title/JSON")
    props = (data or {}).get("PropertyTable", {}).get("Properties", [])
    return props[0] if props else None


def by_cid(cid: int) -> dict[str, Any] | None:
    data = _get(f"{BASE}/pug/compound/cid/{cid}/property/MolecularFormula,Title/JSON")
    props = (data or {}).get("PropertyTable", {}).get("Properties", [])
    return props[0] if props else None


@cache
def _local_groups() -> dict[int, list[str]]:
    """적재된 reactive_groups 테이블 (msl db load cameo). 없으면 빈 dict."""
    try:
        from msl.db import connect, tables

        if "reactive_groups" not in tables():
            return {}
        rows = connect().sql("SELECT cid, list(DISTINCT reactive_group ORDER BY reactive_group) FROM reactive_groups GROUP BY cid").fetchall()
        return {int(cid): groups for cid, groups in rows}
    except Exception:
        return {}


def reactive_groups(cid: int) -> list[str]:
    """CAMEO Chemicals 가 부여한 반응성 그룹. 적재된 로컬 테이블을 먼저 보고, 없으면 PubChem 에 묻는다."""
    local = _local_groups()
    if cid in local:
        return local[cid]
    data = _get(f"{BASE}/pug_view/data/compound/{cid}/JSON?heading=Reactive+Group")
    found: set[str] = set()

    def walk(node: Any) -> None:
        if isinstance(node, dict):
            for item in node.get("StringWithMarkup", []):
                if item.get("String"):
                    found.add(item["String"])
            for value in node.values():
                walk(value)
        elif isinstance(node, list):
            for value in node:
                walk(value)

    walk((data or {}).get("Record", {}))
    return sorted(found)


def cas_number(cid: int) -> str | None:
    """PubChem 동의어에서 검사 숫자가 맞는 첫 CAS 번호."""
    from msl.schema.refs import cas_checksum_ok

    data = _get(f"{BASE}/pug/compound/cid/{cid}/synonyms/JSON")
    info = (data or {}).get("InformationList", {}).get("Information", [])
    for syn in (info[0].get("Synonym", []) if info else []):
        if re.fullmatch(r"\d{2,7}-\d{2}-\d", syn) and cas_checksum_ok(syn):
            return syn
    return None
