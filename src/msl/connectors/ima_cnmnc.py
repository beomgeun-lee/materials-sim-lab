"""IMA-CNMNC 광물 목록 (CC BY-SA 3.0) — 광물 정본명·승인 화학식·상태 (계획서 4.1절 M 레이어, 결정 D4).

원본은 PDF 뿐이다. `pdftotext -tsv` 로 낱말마다 좌표를 받아 표를 복원한다.
`-layout` 텍스트는 따로 찍힌 위첨자(산화수)를 윗줄·아랫줄로 떼어 내고 두 줄 화학식을 옆 행과 섞어서 쓰지 않는다.
- 이름 칸은 쪽의 왼쪽 끝에서 시작한다. 한 행의 다른 칸(화학식·나라)은 이름 줄을 가운데 두고 위아래로 감긴다.
- 쪽마다 열 경계를 데이터에서 잰다 (표 머리는 첫 쪽에만 있다).
- 위·아래 첨자는 작은 글씨 낱말이다. 가까운 줄에 붙여 x 순서로 이으면 원문 표기 'Cu2+Mn3+6O8' 이 된다.

매칭 테이블 mineral_structures(match())는 COD structures_exp 가, mineral_mp(match_mp())는 MP phases_calc 가 적재된 뒤
만든다 (refresh_matches() 가 적재된 것만 골라 다시 만든다). IMA 목록의 파생물이므로 SA 조건을 그대로 이어
CC-BY-SA-3.0 으로 둔다 (MP 의 CC-BY-4.0 보다 SA 쪽이 더 제한적이다).
"""

from __future__ import annotations

import datetime as dt
import re
import shutil
import statistics
import subprocess
import unicodedata
import warnings
from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import cache
from pathlib import Path

import pandas as pd
from pymatgen.core import Composition, Element

from msl.db import connect, download, raw_dir, record_checksums, tables, with_provenance, write_table

SOURCE = "ima-cnmnc"
TABLES = ["minerals", "mineral_structures", "mineral_mp"]
VERSION = "2026-09"
LICENSE = "CC-BY-SA-3.0"
URL = "https://cnmnc.units.it/files/editor/IMA_Master_List_(2026-09).pdf"
PDF = f"IMA_list_of_minerals_{VERSION}.pdf"
TOTAL = 6239  # PDF 머리말: "the 6239 currently valid species"


def fetch() -> Path:
    d = raw_dir(SOURCE, VERSION)
    pdf = d / PDF
    if not pdf.exists():
        download(URL, pdf)
        record_checksums(d)
    return pdf


# ── PDF → 표 행 ───────────────────────────────────────────────────────────


@dataclass(frozen=True)
class Word:
    page: int
    x: float
    y: float  # 위쪽 끝
    w: float
    h: float
    text: str

    @property
    def right(self) -> float:
        return self.x + self.w

    @property
    def cy(self) -> float:
        return self.y + self.h / 2


def read_tsv(text: str) -> list[Word]:
    """`pdftotext -tsv` 출력에서 낱말(level 5)만 읽는다."""
    words = []
    for line in text.splitlines():
        f = line.split("\t")
        if len(f) >= 12 and f[0] == "5":
            words.append(Word(int(f[1]), float(f[6]), float(f[7]), float(f[8]), float(f[9]), f[11]))
    return words


def pdf_words(pdf: Path) -> list[Word]:
    exe = shutil.which("pdftotext") or "/opt/homebrew/bin/pdftotext"
    out = subprocess.run([exe, "-tsv", str(pdf), "-"], capture_output=True, text=True, check=True).stdout
    return read_tsv(out)


STATUS = re.compile(r"(A|G|Rd|Rn|Q)\??")
YEAR = re.compile(r"\d{4}(-(\d{3}[a-z]?|xxx))?|\?")
IMA_NO = re.compile(r"\d{4}-(\d{3}[a-z]?|xxx)")


def _small(w: Word, body_h: float) -> bool:
    return w.h < 0.8 * body_h  # 첨자 (위첨자+아래첨자가 겹친 '3+6' 은 본문 높이보다 크다)


def _lines(words: list[Word], body_h: float) -> list[list[Word]]:
    """본문 글씨로 줄을 묶고, 첨자는 세로 중심이 가장 가까운 줄에 붙인다."""
    body = sorted((w for w in words if not _small(w, body_h)), key=lambda w: w.y)
    lines: list[list[Word]] = []
    for w in body:
        if lines and abs(lines[-1][0].y - w.y) < 2:
            lines[-1].append(w)
        else:
            lines.append([w])
    for w in words:
        if _small(w, body_h) and lines:
            min(lines, key=lambda ln: abs(w.cy - _center(ln, body_h))).append(w)
    return [sorted(ln, key=lambda w: w.x) for ln in lines]


def _center(line: list[Word], body_h: float) -> float:
    return statistics.fmean(w.cy for w in line if not _small(w, body_h))


def _text(line: list[Word]) -> str:
    """한 줄을 x 순서로 잇는다. 낱말 사이가 띄어쓰기만큼(1.8pt 초과) 벌어진 곳만 띄운다."""
    out, right = "", None
    for w in line:
        if right is not None and w.x - right > 1.8:
            out += " "
        out += w.text
        right = w.right if right is None else max(right, w.right)
    return out


def _mode(values: list[float]) -> float:
    return Counter(round(v) for v in values).most_common(1)[0][0]


def _parse_page(ws: list[Word]) -> list[dict]:
    body_h = Counter(round(w.h, 2) for w in ws).most_common(1)[0][0]
    body = [w for w in ws if abs(w.h - body_h) < 0.5]
    rows_by_y: dict[int, list[Word]] = defaultdict(list)
    for w in body:
        rows_by_y[round(w.y)].append(w)

    def same_line(a: Word) -> list[Word]:
        near = [w for k in (round(a.y) - 1, round(a.y), round(a.y) + 1) for w in rows_by_y.get(k, [])]
        return sorted((w for w in near if abs(w.y - a.y) < 1.5 and w is not a), key=lambda w: w.x)

    x0 = min(w.x for w in body)
    names = [(w, same_line(w)) for w in sorted(body, key=lambda w: w.y) if w.x < x0 + 1]
    # 상태 열: '상태 기호 + 연도' 가 붙어 나오는 줄에서 가운데 위치를 잰다 → 연도가 빈 행도 잡는다
    anchors = [ln[i] for _, ln in names for i in range(len(ln) - 1)
               if STATUS.fullmatch(ln[i].text) and YEAR.fullmatch(ln[i + 1].text)]
    if not anchors:
        return []
    status_cx = statistics.median(w.x + w.w / 2 for w in anchors)
    rows = []
    for n, ln in names:
        st = next((w for w in ln if STATUS.fullmatch(w.text) and abs(w.x + w.w / 2 - status_cx) < 6), None)
        if st is not None:
            rows.append((n, ln, st))
    if not rows:
        return []

    # 열 경계: 이름 | 화학식 | 상태·연도 | 나라 | 참고문헌. 화학식 열 = 이름 뒤 처음으로 크게 벌어진 곳
    firsts = []
    for n, ln, st in rows:
        right = n.right
        for w in ln:
            if w.x - right > 8:
                if w is not st:
                    firsts.append(w.x)
                break
            right = max(right, w.right)
    fx = _mode(firsts) - 1 if firsts else max(n.right for n, _, _ in rows) + 1
    sx = min(st.x for *_, st in rows) - 2
    after_year = []
    for _, ln, st in rows:
        rest = [w for w in ln if w.x > st.right]
        while rest and (YEAR.fullmatch(rest[0].text) or rest[0].text == "s.p."):
            rest.pop(0)
        if rest:
            after_year.append(rest[0].x)
    cx = _mode(after_year) - 1 if after_year else float("inf")
    starts: Counter[int] = Counter()  # 참고문헌 열 = 나라 열 오른쪽에서 줄이 시작되는 x 중 가장 왼쪽
    for ln in _lines([w for w in body if w.x > cx + 5], body_h):
        prev = None
        for w in ln:
            if prev is None or w.x - prev.right > 8:
                starts[round(w.x)] += 1
            prev = w
    rx = min((x for x, c in starts.items() if c >= 2), default=float("inf")) - 1

    centers = [n.cy for n, _, _ in rows]
    cells: list[dict[str, list[list[Word]]]] = [defaultdict(list) for _ in rows]
    for col, lo, hi in (("name", -1.0, fx), ("formula", fx, sx), ("status", sx, cx), ("country", cx, rx)):
        for line in _lines([w for w in ws if lo <= w.x < hi], body_h):
            cy = _center(line, body_h)
            i = min(range(len(centers)), key=lambda k: abs(centers[k] - cy))
            if abs(centers[i] - cy) < 2.5 * body_h:
                cells[i][col].append(line)

    out = []
    for (n, _, st), cell in zip(rows, cells):
        texts = {col: [_text(line) for line in sorted(lines, key=lambda ln: _center(ln, body_h))]
                 for col, lines in cell.items()}
        formula = ""
        for s in texts.get("formula", []):
            formula += (" " if formula and " " in s else "") + s  # 감긴 화학식은 붙이고, 설명(공백 포함)은 띄운다
        tokens = [t for s in texts.get("status", []) for t in s.split()]
        tokens.remove(st.text)  # 연도 칸이 두 줄이면 상태 기호 위아래로 나뉜다 ('1982 s.p.' / 'A' / '?')
        years = [t for t in tokens if re.match(r"\d{4}", t)]
        out.append({
            "name": _text(next(line for line in cell["name"] if n in line)),
            "formula_ima": formula or None,
            "status": st.text,
            "ima_number": next((t for t in tokens if IMA_NO.fullmatch(t)), None),
            "year": int(years[0][:4]) if years else None,
            "country": " ".join(texts.get("country", [])) or None,
        })
    return out


def parse_words(words: list[Word]) -> list[dict]:
    """낱말 좌표 → 광물 행. 표 머리('Name … CNMMN/CNMNC') 앞(머리말 쪽)은 건너뛴다."""
    header = next((w for w in words if w.text == "Name" and any(
        v.text == "CNMMN/CNMNC" and v.page == w.page and abs(v.y - w.y) < 1 for v in words)), None)
    if header is not None:
        words = [w for w in words if (w.page, w.y) > (header.page, header.y + 8)]
    pages: dict[int, list[Word]] = defaultdict(list)
    for w in words:
        pages[w.page].append(w)
    return [row for _, ws in sorted(pages.items()) for row in _parse_page(ws)]


# ── 화학식 ────────────────────────────────────────────────────────────────

SYMBOLS = set(Element.__members__) - {"D", "T"}
# 가운뎃점 변형·빼기 기호, 라틴 글자와 모양이 같은 키릴 글자 ('ОН' 이 섞인 행이 있다)
_NORMAL = str.maketrans({"∙": "·", "•": "·", "⋅": "·", "−": "-", "–": "-",
                         "О": "O", "Н": "H", "С": "C", "Р": "P", "В": "B", "К": "K", "М": "M", "Т": "T"})
_OPEN = {"(": ")", "[": "]", "{": "}"}
_NUM = re.compile(r"\d+/\d+|\d+(?:\.\d+)?|\.\d+")


class FormulaError(ValueError):
    """축약식으로 바꿀 수 없는 표기 (치환 쉼표, 변수 x·n, REE 같은 자리표시 등)."""


def _symbol(s: str, i: int) -> tuple[str, int]:
    if i + 1 < len(s) and s[i + 1].islower() and s[i:i + 2] in SYMBOLS:
        return s[i:i + 2], i + 2
    if s[i] in SYMBOLS:
        return s[i], i + 1
    raise FormulaError(f"원소 기호가 아님: {s[i:i + 3]!r}")


def _number(s: str, i: int) -> tuple[float, int]:
    m = _NUM.match(s, i)
    if not m:
        return 1.0, i
    num, _, den = m.group().partition("/")
    return float(num) / float(den or 1), m.end()


def _charge(s: str, i: int) -> int:
    """원소 기호·닫는 괄호 바로 뒤의 산화수 표기('3+', '2-', '+')를 건너뛴다 (Mn3+6 = Mn³⁺ × 6)."""
    if i + 1 < len(s) and s[i].isdigit() and s[i + 1] in "+-":
        return i + 2
    if i < len(s) and s[i] in "+-":
        return i + 1
    return i


def _group(s: str, i: int, close: str | None) -> tuple[dict[str, float], int]:
    counts: dict[str, float] = {}
    while i < len(s):
        c = s[i]
        if c in _OPEN:
            sub, i = _group(s, i + 1, _OPEN[c])
            n, i = _number(s, _charge(s, i))
            for el, k in sub.items():
                counts[el] = counts.get(el, 0.0) + k * n
        elif c in ")]}":
            if c != close:
                raise FormulaError(f"괄호 짝이 안 맞음: {s!r}")
            return counts, i + 1
        elif c.isupper():
            el, i = _symbol(s, i)
            n, i = _number(s, _charge(s, i))
            counts[el] = counts.get(el, 0.0) + n
        else:
            raise FormulaError(f"읽을 수 없는 문자 {c!r}: {s!r}")
    if close is not None:
        raise FormulaError(f"괄호가 닫히지 않음: {s!r}")
    return counts, i


def parse_formula(text: str) -> dict[str, float]:
    """IMA 승인 화학식(평문) → 원소별 개수. 못 읽으면 FormulaError.

    공공(□)은 원자가 없으므로 빼고, 자리 합계 표기(Σ4)·근삿값 표시(~)·의문 표시 '(?)' 는 떼고 읽는다.
    """
    s = re.sub(r"\s+", "", text.translate(_NORMAL))
    s = re.sub(r"\(\?\)$", "", s)
    s = re.sub(r"[☐□](\d+(?:\.\d+)?)?|Σ\d+(?:\.\d+)?|~", "", s)
    if not s:
        raise FormulaError("빈 화학식")
    counts: dict[str, float] = {}
    for k, part in enumerate(s.split("·")):  # 수화물: ·3H2O, ·0.5H2O
        n, j = (1.0, 0) if k == 0 else _number(part, 0)
        sub, _ = _group(part[j:], 0, None)
        if not sub:
            raise FormulaError(f"빈 성분: {text!r}")
        for el, c in sub.items():
            counts[el] = counts.get(el, 0.0) + c * n
    return counts


def formula_elements(text: str | None) -> list[str] | None:
    """화학식에 나오는 원소 (치환식도 포함). REE·Ln 같은 자리표시가 있으면 원소 집합을 확정할 수 없어 None."""
    if not text:
        return None
    found: set[str] = set()
    for m in re.finditer(r"[A-Z][a-z]?", text.translate(_NORMAL)):
        tok = m.group()
        if tok in SYMBOLS:
            found.add(tok)
        elif tok[0] in SYMBOLS and len(tok) == 2:  # 'Sbx' 처럼 변수가 붙은 경우
            found.add(tok[0])
        else:
            return None
    return sorted(found) or None


def chem(text: str | None) -> dict:
    """formula_reduced(pymatgen 축약식, 못 읽으면 None)·elements(쉼표 목록)."""
    reduced = None
    if text:
        try:
            reduced = Composition(parse_formula(text)).reduced_formula
        except (FormulaError, ValueError):
            reduced = None
    els = formula_elements(text)
    return {"formula_reduced": reduced, "elements": ",".join(els) if els else None}


# ── 적재 ──────────────────────────────────────────────────────────────────


def minerals_loaded() -> bool:
    return "minerals" in tables()


def load() -> dict[str, dict[str, int]]:
    pdf = fetch()
    rows = parse_words(pdf_words(pdf))
    df = pd.DataFrame([{**r, **chem(r["formula_ima"])} for r in rows])
    df = df.drop_duplicates("name")
    if len(df) != TOTAL:  # 새 판을 받으면 VERSION·TOTAL 을 함께 고친다
        warnings.warn(f"IMA 목록 {len(df)}종 복원 — PDF 머리말의 {TOTAL}종과 다름", stacklevel=2)
    retrieved = dt.datetime.fromtimestamp(pdf.stat().st_mtime, dt.UTC).isoformat(timespec="seconds")
    df = with_provenance(df, source=SOURCE, version=VERSION, license=LICENSE, method="compiled",
                         id_column="name", retrieved_at=retrieved)
    return {"minerals": write_table("minerals", df)} | refresh_matches()


def refresh_matches() -> dict[str, dict[str, int]]:
    """광물 목록에서 파생되는 매칭 테이블을 지금 적재된 입력으로 다시 만든다.

    IMA·COD·MP 어느 것을 먼저 적재해도 마지막 적재 뒤에는 최신이 된다. mineral_mp 는 대표 COD 구조의 공간군을
    쓰므로 mineral_structures 뒤에 만든다.
    """
    have = set(tables())
    out: dict[str, dict[str, int]] = {}
    if "minerals" not in have:
        return out
    if "structures_exp" in have:
        out["mineral_structures"] = match()
    if "phases_calc" in have:
        out["mineral_mp"] = match_mp()
    return out


# ── 광물명 ↔ COD 구조 매칭 ────────────────────────────────────────────────

_LETTERS = str.maketrans({"ø": "o", "æ": "ae", "œ": "oe", "ł": "l", "ß": "ss", "đ": "d", "ı": "i", "þ": "th"})


def name_key(name: str) -> str:
    """광물명 비교 키: 발음 구별 기호 제거, 소문자, 글자·숫자만 (Abenakiite-(Ce) ≡ abenakiite (Ce) ≡ ABENAKIITE-CE)."""
    s = unicodedata.normalize("NFKD", name.lower().translate(_LETTERS))
    s = "".join(c for c in s if not unicodedata.combining(c))
    return re.sub(r"[^a-z0-9]", "", s)


# COD 광물명에 붙는 수식어 (광물종이 아니라 상태·변종·시료 표시)
_MODIFIER = re.compile(r"(low|high|alpha|beta|gamma|delta|natural|synthetic|deuterated|dehydrated|hydrated|"
                       r"proto|ordered|disordered|ht|lt|hp|[a-z]+(ian|oan))", re.IGNORECASE)
_POLYTYPE = re.compile(r"[\s_-]+\d+[A-Z]{1,2}\d*[a-z]?$")  # Muscovite-2M1, Muscovite 3T, Lizardite-2H1
# 표기만 다르거나 IMA 가 이름을 바꾼 경우 (비교 키 → 비교 키). 뜻이 하나로 정해지는 것만 둔다
SYNONYMS = {"barite": "baryte", "sulfur": "sulphur", "aluminum": "aluminium", "wuestite": "wustite",
            "rhodocrosite": "rhodochrosite", "mayenite": "chlormayenite",
            "ferrocolumbite": "columbitefe", "manganocolumbite": "columbitemn",
            "ferrotantalite": "tantalitefe", "manganotantalite": "tantalitemn"}


def build_index(names) -> dict[str, str]:
    """IMA 이름 목록 → {비교 키: IMA 이름}. 이름 뿌리(접미 '-(Ce)', '-Ca' 앞)가 한 종에만 쓰이면 뿌리 키도 넣는다
    (COD 'Loparite' → IMA 'Loparite-(Ce)'). 뿌리가 여러 종에 걸치면(Monazite-(Ce/La/Nd…)) 넣지 않는다."""
    index = {name_key(n): n for n in names}
    roots: dict[str, list[str]] = defaultdict(list)
    for n in names:
        root = re.sub(r"-(\([^)]*\)|[A-Z][a-z]?)$", "", n)
        if root != n:
            roots[name_key(root)].append(n)
    return index | {k: v[0] for k, v in roots.items() if len(v) == 1 and k not in index}


def _variants(cod_name: str, index: dict[str, str]) -> list[str]:
    """COD 광물명에서 IMA 이름 후보를 만든다 (앞의 것부터 시도)."""
    s = cod_name.strip()
    base = s.split(",")[0]  # 'Calcite, magnesian'
    # 괄호 주석('(deuterated)', 잘린 '(F')은 떼되, 종을 가르는 원소 접미('(Ce)', '(La)')는 남긴다
    plain = re.sub(r"\s*\((?!\s*(?:[A-Z][a-z]?\s*)+\))[^)]*(\)|$)", "", base).strip()
    nopoly = _POLYTYPE.sub("", plain)
    words = [w for w in re.split(r"[\s_-]+", nopoly) if w]
    kept = [w for w in words if not _MODIFIER.fullmatch(w) or name_key(w) in index]  # Celsian 은 수식어가 아니다
    return list(dict.fromkeys(v for v in (s, base, plain, nopoly, " ".join(kept)) if v))


def match_name(cod_name: str | None, index: dict[str, str]) -> str | None:
    """COD 광물명 → IMA 정본명 (index 는 build_index 결과). 못 찾으면 None.

    순서: 그대로 → 쉼표 뒤 수식어 제거 → 괄호 주석 제거 → 폴리타입 기호 제거 → 수식어 낱말 제거.
    단계마다 표기 동의어(barite → baryte)와 독일어식 모음(ue → u)도 본다.
    """
    if not cod_name:
        return None
    for v in _variants(cod_name, index):
        k = name_key(v)
        for cand in (k, SYNONYMS.get(k, k), k.replace("ue", "u").replace("oe", "o").replace("ae", "a")):
            if cand in index:
                return index[cand]
    return None


AMBIENT_T = (280.0, 310.0)  # K
AMBIENT_P = 200.0  # kPa 이하 (상압 101.325 kPa)


def _els(csv: str | None) -> frozenset[str]:
    """쉼표 원소 목록 → 집합 (H 제외: X선 구조는 H 위치를 자주 빼서 COD 화학식에 H 가 없는 경우가 많다)."""
    return frozenset(e for e in csv.split(",") if e) - {"H"} if isinstance(csv, str) else frozenset()


def _rank(row, ima_els: frozenset[str] | None, ima_formula: str | None, modal_sg: float | None) -> tuple:
    t, p = row.cell_temp, row.cell_pressure
    ambient = (pd.isna(t) or AMBIENT_T[0] <= t <= AMBIENT_T[1]) and (pd.isna(p) or p <= AMBIENT_P)
    cod = _els(row.elements)
    if ima_formula:  # 0 축약식 같음 · 1 원소 집합만 같음 · 2 다름
        chem_fit = 0 if row.formula_reduced == ima_formula else 1 if cod == ima_els else 2
    else:  # 치환식: COD 원소가 IMA 원소 목록 안에 있으면 맞는 것으로 본다
        chem_fit = 0 if ima_els is None or cod <= ima_els else 2
    r = row.r_obs if not pd.isna(row.r_obs) else float("inf")
    return (not row.has_coordinates, not pd.isna(row.duplicate_of), row.cod_status in ("errors", "retracted"),
            chem_fit, not ambient, modal_sg is not None and row.sg_number != modal_sg, r,
            -(row.year if not pd.isna(row.year) else 0), row.cod_id)


def pick_best(cod: pd.DataFrame, ima_elements: str | None, ima_formula: str | None):
    """한 광물의 COD 구조 중 대표 하나 (itertuples 행).

    좌표 있음 → 중복 아님 → 오류·철회 표시 없음 → 화학식이 IMA 와 맞음(축약식 같음 > 원소 집합 같음, H 제외)
    → 상온·상압 → 가장 흔한 공간군 → R_obs 낮음 → 최신 → COD 번호 순.
    가장 흔한 공간군을 보는 것은 AMCSD 가 출발 물질 이름을 붙인 고압상(예: 'Spinel' 이름의 CaFe2O4형)을 거르기 위해서다.
    """
    ima_els = _els(ima_elements) if isinstance(ima_elements, str) else None
    ima_formula = ima_formula if isinstance(ima_formula, str) else None
    usable = cod[cod["has_coordinates"] & cod["duplicate_of"].isna()]
    sgs = (usable if len(usable) else cod)["sg_number"].dropna()
    modal_sg = float(sgs.mode().iloc[0]) if len(sgs) else None
    return min(cod.itertuples(), key=lambda r: _rank(r, ima_els, ima_formula, modal_sg))


def common_fits(cod_elements: str | None, ima_elements: str | None, exact: bool) -> bool:
    """commonname 으로만 이름이 맞은 COD 항목을 받을지. commonname 은 화학 관용명이라 같은 이름의 합성 유사체
    (GeO₂ 'quartz')나 조성이 다른 것이 섞이므로 원소 구성이 IMA 와 맞을 때만 받는다 (H 제외).
    축약식이 있는 광물은 원소 집합이 같아야 하고, 치환식 광물은 COD 원소가 IMA 원소 목록 안에 있어야 한다.
    IMA 원소 목록을 모르면(REE 같은 자리표시) 받지 않는다."""
    if not isinstance(ima_elements, str):
        return False
    cod, ima = _els(cod_elements), _els(ima_elements)
    return bool(cod) and (cod == ima if exact else cod <= ima)


def match() -> dict[str, int]:
    """minerals × structures_exp → mineral_structures (IMA 광물마다 COD 구조 목록과 대표 구조).

    COD 이름은 광물명(mineral)을 먼저 보고, 없으면 commonname 을 본다 (commonname 항목은 common_fits 를 통과할 때만).
    commonname 항목은 R 값이 있어 AMCSD 광물 항목(R 값 없음)보다 앞서기 쉬우므로, 대표 구조는 광물명 항목이 하나도
    없는 광물에서만 commonname 항목 중에 고른다.
    formula_check: IMA 축약식과 대표 COD 구조의 원소 집합이 같은가. X선 구조는 H 위치를 자주 빼므로 H 는 비교에서 뺀다.
    IMA 화학식을 축약식으로 못 바꾼 광물(치환식 등)은 None.
    """
    con = connect()
    minerals = con.execute("SELECT name, formula_reduced, elements FROM minerals").df()
    common = "commonname" in {r[0] for r in con.execute("DESCRIBE structures_exp").fetchall()}  # 옛 적재본에는 없다
    cod = con.execute("SELECT * FROM structures_exp WHERE mineral IS NOT NULL"
                      + (" OR commonname IS NOT NULL" if common else "")).df()
    cod_version = str(cod["source_version"].iloc[0]) if len(cod) else "none"
    index = build_index(list(minerals["name"]))
    cod["from_common"] = cod["mineral"].isna()
    cod["cod_name"] = cod["mineral"].fillna(cod["commonname"]) if common else cod["mineral"]
    names = {m: match_name(m, index) for m in cod["cod_name"].unique()}
    cod["ima_name"] = cod["cod_name"].map(names)
    ima = minerals.set_index("name")
    rejected = [r.from_common and isinstance(r.ima_name, str) and not common_fits(
        r.elements, ima.at[r.ima_name, "elements"], isinstance(ima.at[r.ima_name, "formula_reduced"], str))
        for r in cod.itertuples()]
    cod.loc[rejected, "ima_name"] = None
    groups = {k: g for k, g in cod.dropna(subset=["ima_name"]).groupby("ima_name")}

    rows = []
    for m in minerals.itertuples():
        g = groups.get(m.name)
        base = {"name": m.name, "n_cod": 0, "n_cod_commonname": 0, "cod_ids": [], "cod_names": None,
                "best_cod_id": None, "best_formula": None, "best_sg": None, "formula_check": None}
        if g is None:
            rows.append(base)
            continue
        exact = isinstance(m.formula_reduced, str)
        named = g[~g["from_common"]]  # 대표는 광물명 항목에서 고른다. commonname 항목은 광물명 항목이 없을 때만
        best = pick_best(named if len(named) else g, m.elements, m.formula_reduced)
        rows.append({**base, "n_cod": len(g), "n_cod_commonname": int(g["from_common"].sum()),
                     "cod_ids": sorted(int(i) for i in g["cod_id"]),
                     "cod_names": "; ".join(sorted(g["cod_name"].unique())),
                     "best_cod_id": int(best.cod_id), "best_formula": best.formula_reduced, "best_sg": best.sg,
                     "formula_check": (_els(best.elements) == _els(m.elements)) if exact and isinstance(m.elements, str)
                     else None})
    df = pd.DataFrame(rows)
    df["best_cod_id"] = df["best_cod_id"].astype("Int64")
    df["formula_check"] = df["formula_check"].astype("boolean")
    df = with_provenance(df, source=SOURCE, version=f"{VERSION}+cod-{cod_version}", license=LICENSE,
                         method="compiled", id_column="name")
    return write_table("mineral_structures", df)


# ── 광물 ↔ MP 계산 구조 매칭 ──────────────────────────────────────────────

LEVEL_SG = "formula+spacegroup"  # 축약식 + 대표 COD 구조와 같은 공간군
LEVEL_SUPER = "formula+supergroup"  # 축약식 + MP 공간군이 대표 COD 공간군의 바로 위 군
LEVEL_FORMULA = "formula"  # 축약식만
# 거울상 쌍(좌·우형). 같은 구조라 DFT 에너지가 같고, MP 와 COD 가 서로 다른 쪽을 올린 경우가 많다 (석영 P3₁21 ↔ P3₂21)
ENANTIOMORPHS = {76: 78, 91: 95, 92: 96, 144: 145, 151: 153, 152: 154, 169: 170, 171: 172, 178: 179, 180: 181,
                 212: 213}
_ENANTIO = {n: min(a, b) for a, b in ENANTIOMORPHS.items() for n in (a, b)}


def sg_key(number) -> int | None:
    """공간군 비교 키. 번호로 비교하므로 H-M 기호의 설정·표기 차이(Pnma = Pmcn = Pbnm = 62번)는 따로 맞출 필요가
    없고, 거울상 쌍은 한 번호로 묶는다. 번호가 없으면 None."""
    if number is None or pd.isna(number):
        return None
    return _ENANTIO.get(int(number), int(number))


@cache
def _maximal_subgroups() -> dict[int, frozenset[int]]:
    from pymatgen.symmetry.groups import SYMM_DATA

    return {int(k): frozenset(int(v) for v in vs) for k, vs in SYMM_DATA["maximal_subgroups"].items()}


def sg_tier(mp_sg, cod_sg) -> int:
    """0 같은 공간군 · 1 MP 공간군이 COD 공간군의 바로 위 군 (COD 군이 MP 군의 극대 부분군) · 2 그 밖·모름.

    1 은 옛 정밀화가 대칭을 낮춰 기술한 구조(흑연 P6₃mc ↔ P6₃/mmc, NiAs 형 황화물)나 DFT 이완이 대칭을 올린 경우다.
    """
    a, b = sg_key(mp_sg), sg_key(cod_sg)
    if a is None or b is None:
        return 2
    if a == b:
        return 0
    return 1 if int(cod_sg) in _maximal_subgroups().get(int(mp_sg), ()) else 2


@cache
def _reduced(formula: str) -> str | None:
    try:
        with warnings.catch_warnings():  # 비활성 기체(MP 의 Ar·He·Ne 단체)의 전기음성도 경고
            warnings.simplefilter("ignore")
            return Composition(formula).reduced_formula
    except Exception:  # 원소가 아닌 기호 등
        return None


def mp_formula(ima_formula, ima_elements, cod_formula) -> tuple[str | None, str | None]:
    """MP 와 비교할 축약식과 그 출처 ('ima' · 'cod').

    IMA 축약식이 있으면 그것만 쓴다 (COD 식은 X선 구조라 H 가 빠진 경우가 많아 Mg(OH)2 가 MgO2 에 붙는다).
    치환식이라 IMA 축약식이 없으면 대표 COD 구조의 식을 쓰되, COD 원소가 IMA 원소 목록 안에 있고 H 유무가 같을 때만.
    """
    if isinstance(ima_formula, str):
        return _reduced(ima_formula), "ima"
    if not isinstance(cod_formula, str) or (f := _reduced(cod_formula)) is None:
        return None, None
    if isinstance(ima_elements, str):
        ima_els, cod_els = set(ima_elements.split(",")), {str(e) for e in Composition(f).elements}
        if not cod_els <= ima_els or ("H" in ima_els) != ("H" in cod_els):
            return None, None
    return f, "cod"


def match_mp_frame(minerals: pd.DataFrame, mp: pd.DataFrame) -> pd.DataFrame:
    """광물마다 MP 후보 목록과 대표 하나 (match_mp 의 계산 부분 — DB·네트워크 없음).

    minerals: name, formula_reduced, elements, best_formula(대표 COD 식), cod_sg_number(대표 COD 공간군 번호)
    mp: material_id, formula_pretty, spacegroup_symbol, spacegroup_number, energy_above_hull_eV, deprecated

    후보는 축약식이 같고 deprecated 가 아닌 MP 항목이다. 순서: 대표 COD 구조와 같은 공간군 → 바로 위 군 → 나머지,
    그 안에서 energy_above_hull 낮은 순 → material_id.
    대표(best_mp_id):
    - 같은 공간군 후보가 있으면 그 첫째 (formula+spacegroup).
    - 없으면, 같은 축약식의 다른 광물이 공간군까지 맞춰 대표로 잡은 항목은 건너뛴다 (vaterite 가 calcite 의 항목을
      갖지 않게).
    - 공간군 근거 없이 축약식만 맞았는데 같은 축약식의 IMA 광물이 둘 이상이면(다형) 대표를 두지 않는다
      (polymorph_ambiguous). 이때도 후보 목록 mp_ids 는 남긴다.
    mp_shared_with 는 대표가 같은 다른 광물이다. 축약식·공간군이 모두 같은 다른 구조(규회석·브레이석)는 공간군으로
    가를 수 없어서 생긴다.
    """
    live = mp[~mp["deprecated"].fillna(False).astype(bool)]
    live = live.assign(key=live["formula_pretty"].map(lambda f: _reduced(f) if isinstance(f, str) else None),
                       e_sort=live["energy_above_hull_eV"].fillna(float("inf")))
    by_formula = {f: g for f, g in live.groupby("key")}
    basis = [mp_formula(m.formula_reduced, m.elements, m.best_formula) for m in minerals.itertuples()]
    n_same = Counter(f for f, _ in basis if f)

    ranked: dict[str, list] = {}
    for m, (f, _) in zip(minerals.itertuples(), basis):
        g = by_formula.get(f)
        tiers = [] if g is None else [sg_tier(sg, m.cod_sg_number) for sg in g["spacegroup_number"]]
        ranked[m.name] = [] if g is None else sorted(
            zip(tiers, g.itertuples()), key=lambda c: (c[0], c[1].e_sort, c[1].material_id))
    claimed = {c[0][1].material_id for c in ranked.values() if c and c[0][0] == 0}

    rows = []
    for m, (f, src) in zip(minerals.itertuples(), basis):
        cands = ranked[m.name]
        best, level, ambiguous = None, None, None
        if cands:
            tier, best = cands[0] if cands[0][0] == 0 else \
                next((c for c in cands if c[1].material_id not in claimed), (2, None))
            level, ambiguous = (LEVEL_SG, LEVEL_SUPER, LEVEL_FORMULA)[tier], False
            if tier == 2 and n_same[f] > 1:
                best, ambiguous = None, True
        rows.append({
            "name": m.name, "formula_basis": src, "match_formula": f,
            "cod_sg_number": None if pd.isna(m.cod_sg_number) else int(m.cod_sg_number),
            "n_mp": len(cands), "mp_ids": [c[1].material_id for c in cands],
            "best_mp_id": None if best is None else best.material_id, "match_level": level,
            "polymorph_ambiguous": ambiguous,
            "mp_sg": None if best is None else best.spacegroup_symbol,
            "mp_sg_number": None if best is None else best.spacegroup_number,
            "mp_e_above_hull_summary_eV": None if best is None else best.energy_above_hull_eV,
        })
    df = pd.DataFrame(rows)
    owners = df.dropna(subset=["best_mp_id"]).groupby("best_mp_id")["name"].agg(list)
    df["mp_shared_with"] = ["; ".join(n for n in owners[b] if n != name) or None if isinstance(b, str) else None
                            for name, b in zip(df["name"], df["best_mp_id"])]
    for col in ("cod_sg_number", "mp_sg_number"):
        df[col] = df[col].astype("Int64")
    df["polymorph_ambiguous"] = df["polymorph_ambiguous"].astype("boolean")
    df["mp_e_above_hull_summary_eV"] = df["mp_e_above_hull_summary_eV"].astype("Float64")
    return df


def match_mp() -> dict[str, int]:
    """minerals (+ mineral_structures · structures_exp) × phases_calc → mineral_mp.

    배포용(open) 입력만 읽는다. 파생 테이블이 SA 로 open 파티션에 가므로 NC 항목(GNoME)이 섞이면 안 된다.
    mp_e_above_hull_summary_eV 는 MP summary 기준(GGA/GGA+U/r²SCAN 혼합)이라 조회용이다. 안정성 판정(A1·A2)은
    best_mp_id 로 GGA/GGA+U 엔트리를 찾아 쓴다 (결정 D12). COD 매칭이 없으면 축약식만으로 맞춘다.
    """
    con = connect(include_restricted=False)
    have = {r[0] for r in con.execute("SELECT table_name FROM information_schema.tables").fetchall()}
    with_cod = {"mineral_structures", "structures_exp"} <= have
    minerals = con.execute(
        "SELECT m.name, m.formula_reduced, m.elements, s.best_formula, e.sg_number AS cod_sg_number FROM minerals m "
        "LEFT JOIN mineral_structures s USING (name) LEFT JOIN structures_exp e ON e.cod_id = s.best_cod_id"
        if with_cod else
        "SELECT name, formula_reduced, elements, NULL AS best_formula, NULL AS cod_sg_number FROM minerals").df()
    mp = con.execute("SELECT material_id, formula_pretty, spacegroup_symbol, spacegroup_number, energy_above_hull_eV, "
                     "deprecated FROM phases_calc").df()
    mp_version, scheme = con.execute(
        "SELECT any_value(source_version), any_value(correction_scheme) FROM phases_calc").fetchone()
    version = VERSION
    if with_cod:
        version += f"+cod-{con.execute('SELECT any_value(source_version) FROM structures_exp').fetchone()[0]}"
    df = with_provenance(match_mp_frame(minerals, mp), source=SOURCE, version=f"{version}+mp-{mp_version}",
                         license=LICENSE, method="compiled", id_column="name",
                         correction_scheme=f"mp_e_above_hull_summary_eV: {scheme}")
    return write_table("mineral_mp", df)
