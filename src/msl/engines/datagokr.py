"""공공데이터포털(data.go.kr) 국내 화학물질 API — 한국환경공단 화학물질·유독물 GHS, KOSHA MSDS.

키는 .env 의 DATA_GO_KR_SERVICE_KEY (포털 '일반 인증키'). Encoding 키(%xx 포함)는 그대로,
Decoding 키는 URL 인코딩해서 붙인다. 응답은 data/cache/datagokr/ 에 캐시해 일일 호출 한도
(개발계정: 한국환경공단 10,000건, KOSHA 2,000건)를 아낀다.
"""

from __future__ import annotations

import hashlib
import json
import os
import ssl
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from typing import Any

import certifi

from msl.env import CACHE_DIR, load_dotenv

BASE = "https://apis.data.go.kr"
CACHE = CACHE_DIR / "datagokr"
_SSL = ssl.create_default_context(cafile=certifi.where())

KECO_CHEM = "B552584/kecoapi/ncissbstn/chemSbstnList"
KECO_GHS = "B552584/kecoapi/ncisghs/ghsList"
KOSHA_LIST = "B552468/msdschem1/getChemList001"
KOSHA_DETAIL = "B552468/msdschem1/getChemDetail{n:02d}1"  # 항목 번호 + '1' (예: 15번 법적 규제 → getChemDetail151)
SEARCH_BY_CAS = {"keco": "2", "kosha": "1"}  # searchGubun / searchCnd 에서 CAS 검색 코드 (2026-09-27 확인)


class DataGoKrUnavailable(RuntimeError):
    """키가 없거나, 등록되지 않았거나, 서비스 오류일 때."""


def _key() -> str:
    load_dotenv()
    key = os.environ.get("DATA_GO_KR_SERVICE_KEY", "")
    if not key:
        raise DataGoKrUnavailable("DATA_GO_KR_SERVICE_KEY 가 .env 에 없음")
    return key if "%" in key else urllib.parse.quote(key, safe="")


def _get(path: str, params: dict[str, str]) -> str:
    query = urllib.parse.urlencode(params)
    cache_id = hashlib.sha256(f"{path}?{query}".encode()).hexdigest()[:24]
    CACHE.mkdir(parents=True, exist_ok=True)
    cached = CACHE / f"{cache_id}.txt"
    if cached.exists():
        return cached.read_text(encoding="utf-8")
    url = f"{BASE}/{path}?serviceKey={_key()}&{query}"
    req = urllib.request.Request(url, headers={"User-Agent": "materials-sim-lab/0.0.1 (+research)"})
    try:
        with urllib.request.urlopen(req, timeout=20, context=_SSL) as resp:
            body = resp.read().decode("utf-8")
    except urllib.error.HTTPError as exc:
        raise DataGoKrUnavailable(f"HTTP {exc.code} ({path})") from exc
    if "SERVICE_KEY_IS_NOT_REGISTERED" in body or "SERVICE_ACCESS_DENIED" in body:
        raise DataGoKrUnavailable("인증키가 이 API 에 등록되지 않음 — 활용신청·동기화 확인")
    if "LIMITED_NUMBER_OF_SERVICE_REQUESTS" in body:
        raise DataGoKrUnavailable("일일 호출 한도 초과")
    cached.write_text(body, encoding="utf-8")
    return body


def _keco(path: str, cas: str) -> list[dict[str, Any]]:
    body = _get(path, {"pageNo": "1", "numOfRows": "10", "returnType": "JSON",
                       "searchGubun": SEARCH_BY_CAS["keco"], "searchNm": cas})
    data = json.loads(body)
    header = data.get("header", {})
    if str(header.get("resultCode")) != "200":
        raise DataGoKrUnavailable(f"한국환경공단 응답 {header.get('resultCode')}: {header.get('resultMsg')}")
    items = data.get("body", {}).get("items") or []
    return [i for i in items if i.get("casNo") == cas]


def keco_substance(cas: str) -> dict[str, Any] | None:
    """화학물질 정보 — 규제 분류(typeList) 포함. 같은 CAS 가 여러 건이면 첫 건(대표 물질)."""
    items = _keco(KECO_CHEM, cas)
    return items[0] if items else None


def keco_ghs(cas: str) -> dict[str, Any] | None:
    """유독물 GHS — 신호어, 그림문자, UN 번호, 유해성 분류(H 코드)."""
    items = _keco(KECO_GHS, cas)
    return items[0] if items else None


def _xml_items(body: str) -> list[dict[str, str]]:
    root = ET.fromstring(body)
    code = root.findtext("./header/resultCode")
    if code not in ("00", "0"):
        raise DataGoKrUnavailable(f"KOSHA 응답 {code}: {root.findtext('./header/resultMsg')}")
    return [{child.tag: (child.text or "").strip() for child in item} for item in root.iter("item")]


def kosha_chem(cas: str) -> dict[str, str] | None:
    body = _get(KOSHA_LIST, {"pageNo": "1", "numOfRows": "5", "searchWrd": cas, "searchCnd": SEARCH_BY_CAS["kosha"]})
    items = [i for i in _xml_items(body) if i.get("casNo") == cas]
    return items[0] if items else None


def kosha_section(chem_id: str, n: int) -> list[tuple[str, str]]:
    """MSDS n번 항목의 (소항목 이름, 내용) 목록. 내용의 '|' 구분은 줄바꿈 목록이다."""
    body = _get(KOSHA_DETAIL.format(n=n), {"chemId": chem_id})
    return [(i.get("msdsItemNameKor", ""), i.get("itemDetail", "")) for i in _xml_items(body) if i.get("itemDetail")]
