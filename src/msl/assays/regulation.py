"""A10 국내 규제 — 한국환경공단 화학물질·유독물 GHS, KOSHA MSDS 법적 규제 (data.go.kr, L0)."""

from __future__ import annotations

from typing import Any

from msl.assays.base import Context, Outcome, not_applicable, pending, val
from msl.engines import datagokr as dg
from msl.schema.result import Fidelity, Status, ValueKind

# 화학물질관리법 등에서 취급 의무가 생기는 분류 — 결과를 '주의'로 표시한다
KEY_CLASSES = ("유독물질", "허가물질", "제한물질", "금지물질", "사고대비물질", "인체등유해성물질",
               "생태등유해성물질", "인체만성유해성물질", "중점관리물질")
ENGINE = "msl-kr-regulation"
EMPTY = {"자료없음", "해당없음", "해당 없음", "-", ""}


def _classes(sub: dict[str, Any] | None) -> list[str]:
    out = []
    for t in (sub or {}).get("typeList", []) or []:
        name = (t.get("sbstnClsfTypeNm") or "").strip()
        if not name or name in ("기존화학물질", "등록대상기존화학물질"):
            continue
        detail = (t.get("contInfo") or "").strip()
        out.append(f"{name}({t.get('unqNo', '').strip()})" + (f" — {detail}" if detail else ""))
    return out


def a10(ctx: Context) -> Outcome:
    comps = [c for c in ctx.comps if c.cas]
    if not comps:
        return not_applicable("CAS 번호가 있는 성분이 없음 (사용자 정의 소재·미해석 성분)", ENGINE)
    rows: list[list[Any]] = []
    kosha: list[dict[str, Any]] = []
    values = []
    flagged = 0
    try:
        for c in comps:
            sub, ghs = dg.keco_substance(c.cas), dg.keco_ghs(c.cas)
            classes = _classes(sub)
            key = [x for x in classes if x.startswith(KEY_CLASSES)]
            flagged += bool(key)
            h_codes = sorted({h.get("hrmDngrCd") for h in (ghs or {}).get("hrmflnList", []) or [] if h.get("hrmDngrCd")})
            name = (sub or {}).get("sbstnNmKor") or c.name
            rows.append([
                name, c.cas, (sub or {}).get("korexst") or "—",
                "; ".join(classes) or ("목록에 없음" if sub is None else "규제 분류 없음 (기존화학물질)"),
                (ghs or {}).get("sfsgwd") or "—",
                ((ghs or {}).get("pctgrmCd") or "—").replace("^", " "),
                " ".join(h_codes) or "—",
            ])
            values.append(val(f"{name} · 국내 규제", ", ".join(x.split("(")[0] for x in key) or "핵심 규제 해당 없음",
                              kind=ValueKind.FLAG))
            ms = dg.kosha_chem(c.cas)
            if ms and ms.get("chemId"):
                kosha.append({"name": ms.get("chemNameKor") or name, "chem_id": ms["chemId"],
                              "updated": ms.get("lastDate"),
                              "items": [(n, t) for n, t in dg.kosha_section(ms["chemId"], 15)
                                        if t.strip() not in EMPTY]})
    except dg.DataGoKrUnavailable as exc:
        return pending(f"data.go.kr 조회 불가: {exc}", ENGINE)

    status = Status.WARNING if flagged else Status.OK
    summary = (f"{flagged}개 성분이 국내 관리 대상 — 취급·보관·신고 의무 확인 필요" if flagged
               else "핵심 규제(유독·사고대비·허가·제한·금지 등) 해당 성분 없음")
    return Outcome(
        status=status, fidelity=Fidelity.L0, engine=ENGINE, engine_version="v0",
        conditions_basis="조회 시점의 고시 목록 (응답은 data/cache/datagokr/ 캐시)",
        values=values, summary=summary,
        sources=[("keco-chem", None), ("keco-ghs", None), ("kosha-msds", None)],
        caveats=[
            "성분별 판정 — 혼합물 함량 기준(예: '10% 이상 함유 혼합물')과 레시피 농도는 아직 비교하지 않음 (v0)",
            "KOSHA MSDS 는 상업 이용 조건 확인 전까지 연구용 (계획서 10절 #2)",
            "법령·고시 개정으로 바뀔 수 있음 — 실제 취급 전 최신 고시 확인",
        ],
        data={"table": {"columns": ["성분", "CAS", "KE 번호", "국내 규제 분류", "GHS 신호어", "그림문자", "H 코드"], "rows": rows},
              "kosha": kosha},
    )
