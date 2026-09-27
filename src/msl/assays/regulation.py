"""A10 국내 규제 — 한국환경공단 화학물질·유독물 GHS, KOSHA MSDS 법적 규제 (data.go.kr, L0).\n\nv1: 분류별 함량 기준(%)을 레시피 질량 분율과 비교한다 (msl.recipe.balance).\n"""

from __future__ import annotations

import re
from typing import Any

from msl.assays.base import Context, Outcome, not_applicable, pending, val
from msl.engines import datagokr as dg
from msl.schema.result import Fidelity, Status, ValueKind

# 화학물질관리법 등에서 취급 의무가 생기는 분류 — 결과를 '주의'로 표시한다
KEY_CLASSES = ("유독물질", "허가물질", "제한물질", "금지물질", "사고대비물질", "인체등유해성물질",
               "생태등유해성물질", "인체만성유해성물질", "중점관리물질")
ENGINE = "msl-kr-regulation"
EMPTY = {"자료없음", "해당없음", "해당 없음", "-", ""}


def _classes(sub: dict[str, Any] | None) -> list[dict[str, Any]]:
    """규제 분류 목록. 함량 기준(예: '10% 이상 함유한 혼합물', '인체급성유해성 : 10%')이 있으면 가장 낮은 % 를 threshold 로."""
    out = []
    for t in (sub or {}).get("typeList", []) or []:
        name = (t.get("sbstnClsfTypeNm") or "").strip()
        if not name or name in ("기존화학물질", "등록대상기존화학물질"):
            continue
        detail = (t.get("contInfo") or "").strip()
        pcts = [float(x) for x in re.findall(r"(\d+(?:\.\d+)?)\s*%", detail)]
        out.append({"name": name, "unq": (t.get("unqNo") or "").strip(), "detail": detail,
                    "threshold": min(pcts) if pcts else None, "key": name.startswith(KEY_CLASSES)})
    return out


def _judge(cls: dict[str, Any], wt_pct: float | None) -> str:
    """함량 기준 판정 문구."""
    if cls["threshold"] is None:
        return "함량 기준 없음"
    if wt_pct is None:
        return f"기준 {cls['threshold']:g}% — 레시피 농도 계산 불가"
    return (f"기준 {cls['threshold']:g}% 이상 → 이 레시피 {wt_pct:.3g}% 로 해당" if wt_pct >= cls["threshold"]
            else f"기준 {cls['threshold']:g}% 미만 (이 레시피 {wt_pct:.3g}%)")


def _applies(cls: dict[str, Any], wt_pct: float | None) -> bool:
    """핵심 분류가 이 레시피 농도에서 해당하는가. 기준이 없거나 농도를 모르면 보수적으로 해당으로 본다."""
    return cls["key"] and (cls["threshold"] is None or wt_pct is None or wt_pct >= cls["threshold"])


def a10(ctx: Context) -> Outcome:
    comps = [c for c in ctx.comps if c.cas]
    if not comps:
        return not_applicable("CAS 번호가 있는 성분이 없음 (사용자 정의 소재·미해석 성분)", ENGINE)
    bal = ctx.shared.get("balance")
    rows: list[list[Any]] = []
    kosha: list[dict[str, Any]] = []
    values = []
    flagged = 0
    try:
        for c in comps:
            sub, ghs = dg.keco_substance(c.cas), dg.keco_ghs(c.cas)
            classes = _classes(sub)
            frac = bal.mass_fraction(str(c.ref)) if bal else None
            wt_pct = frac * 100 if frac is not None else None
            key = [x for x in classes if x["key"]]
            applied = [x for x in key if _applies(x, wt_pct)]
            below = [x for x in key if not _applies(x, wt_pct)]
            flagged += bool(applied)
            h_codes = sorted({h.get("hrmDngrCd") for h in (ghs or {}).get("hrmflnList", []) or [] if h.get("hrmDngrCd")})
            name = (sub or {}).get("sbstnNmKor") or c.name
            rows.append([
                name, c.cas, (sub or {}).get("korexst") or "—",
                "; ".join(f"{x['name']}({x['unq']})" + (f" — {_judge(x, wt_pct)}" if x["key"] else "") for x in classes)
                or ("목록에 없음" if sub is None else "규제 분류 없음 (기존화학물질)"),
                (ghs or {}).get("sfsgwd") or "—",
                ((ghs or {}).get("pctgrmCd") or "—").replace("^", " "),
                " ".join(h_codes) or "—",
            ])
            verdict = ", ".join(x["name"] for x in applied) or (
                "성분은 " + ", ".join(x["name"] for x in below) + " — 이 농도에서는 기준 미만" if below else "핵심 규제 해당 없음")
            values.append(val(f"{name} · 국내 규제", verdict, kind=ValueKind.FLAG))
            if wt_pct is not None:
                values.append(val(f"{name} · 레시피 내 함량", round(wt_pct, 4), "wt%"))
            ms = dg.kosha_chem(c.cas)
            if ms and ms.get("chemId"):
                kosha.append({"name": ms.get("chemNameKor") or name, "chem_id": ms["chemId"],
                              "updated": ms.get("lastDate"),
                              "items": [(n, t) for n, t in dg.kosha_section(ms["chemId"], 15)
                                        if t.strip() not in EMPTY]})
    except dg.DataGoKrUnavailable as exc:
        return pending(f"data.go.kr 조회 불가: {exc}", ENGINE)

    status = Status.WARNING if flagged else Status.OK
    summary = (f"{flagged}개 성분이 이 농도에서 국내 관리 대상 — 취급·보관·신고 의무 확인 필요" if flagged
               else "이 레시피 농도에서 핵심 규제(유독·사고대비·허가·제한·금지 등) 해당 성분 없음")
    return Outcome(
        status=status, fidelity=Fidelity.L0, engine=ENGINE, engine_version="v0",
        conditions_basis="조회 시점의 고시 목록 (응답은 data/cache/datagokr/ 캐시)",
        values=values, summary=summary,
        sources=[("keco-chem", None), ("keco-ghs", None), ("kosha-msds", None)],
        caveats=[
            "함량 기준은 레시피 질량 분율(wt%)과 비교 — 물 부피는 밀도 1 g/mL 로 가정, 농도를 모르면 보수적으로 '해당'",
            "KOSHA MSDS 는 상업 이용 조건 확인 전까지 연구용 (계획서 10절 #2)",
            "법령·고시 개정으로 바뀔 수 있음 — 실제 취급 전 최신 고시 확인",
        ],
        data={"table": {"columns": ["성분", "CAS", "KE 번호", "국내 규제 분류", "GHS 신호어", "그림문자", "H 코드"], "rows": rows},
              "kosha": kosha},
    )
