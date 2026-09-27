"""PHREEQC 공식 예제 재현 — 계획서 6장 3단계 완료 기준 (A4 엔진 경로 검증).

USGS 배포본(phreeqc-3.8.6, 퍼블릭 도메인)의 examples/ex* 입력을
(가) 같은 배포본 소스로 빌드한 공식 phreeqc 실행 파일과 (나) msl A4 가 쓰는 경로(phreeqpython 의 IPhreeqc)로
같은 데이터베이스에서 돌려, 출력의 용액 기술(pH·pe·이온 세기)과 포화지수(SI)를 비교한다.
배포본에는 기대 출력(.out)이 없어서 공식 실행 파일 결과를 기준으로 삼고, 그 값을 kb/validation 에 남긴다
(테스트는 실행 파일 없이 이 기준값과 비교한다).
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

SOL_RE = re.compile(r"-+Description of solution-+\n(.*?)(?:\n-{5,}|\Z)", re.S)
NUM = r"([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)"
KEEP_VALUES = 12  # 공식 결과 값을 통째로 남길 최대 용액 수
SI_BLOCK_RE = re.compile(r"-+Saturation indices-+\n\s*\n\s*Phase\s+SI\*\*.*?\n\s*\n(.*?)(?:\n\s*\n)", re.S)


def parse(out: str) -> dict[str, list]:
    """출력 문자열에서 용액 기술(pH·pe·이온 세기)과 포화지수 표를 순서대로 뽑는다."""
    sols = []
    for block in SOL_RE.findall(out):
        row = {}
        for key, pat in (("pH", rf"^\s*pH\s*=\s*{NUM}"), ("pe", rf"^\s*pe\s*=\s*{NUM}"),
                         ("I", rf"Ionic strength \(mol/kgw\)\s*=\s*{NUM}")):
            m = re.search(pat, block, re.M)
            if m:
                row[key] = float(m.group(1))
        if row:
            sols.append(row)
    sis = []
    for block in SI_BLOCK_RE.findall(out):
        table = {}
        for line in block.splitlines():
            parts = line.split()
            if len(parts) >= 2:
                try:
                    table[parts[0]] = float(parts[1])
                except ValueError:
                    continue
        sis.append(table)
    return {"solutions": sols, "si": sis}


def _new_files(work: Path, before: set[str]) -> dict[str, str]:
    """실행이 새로 만든 텍스트 파일 (SELECTED_OUTPUT -file 등)."""
    skip = {"out.txt", "phreeqc.log"}
    return {f.name: f.read_text(encoding="latin-1") for f in sorted(work.iterdir())
            if f.is_file() and f.name not in before and f.name not in skip}


def run_official(binary: Path, example: Path, db: Path, timeout: int = 300) -> dict[str, Any]:
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        for f in example.parent.iterdir():  # 예제가 같은 폴더의 파일을 읽기도 한다
            if f.is_file():
                shutil.copy2(f, work / f.name)
        before = {f.name for f in work.iterdir()}
        r = subprocess.run([str(binary), example.name, "out.txt", str(db)], cwd=work, timeout=timeout,
                           capture_output=True, check=False)
        out = work / "out.txt"
        return {"out": out.read_text(encoding="latin-1") if out.exists() else "", "files": _new_files(work, before),
                "error": r.stderr.decode("latin-1")[-2000:] if r.returncode else ""}


_OURS = r"""
import ctypes, json, os, shutil, sys, tempfile
from phreeqpython.viphreeqc import VIPhreeqc
ex, db = sys.argv[1], sys.argv[2]
tmp = tempfile.mkdtemp()
for f in os.listdir(os.path.dirname(ex)):
    p = os.path.join(os.path.dirname(ex), f)
    if os.path.isfile(p):
        shutil.copy2(p, tmp)
before = set(os.listdir(tmp))
os.chdir(tmp)
ip = VIPhreeqc()
ip.load_database(db)
d = ip.dll
for fn, args, res in (("SetOutputStringOn", [ctypes.c_int, ctypes.c_int], ctypes.c_int),
                      ("SetSelectedOutputFileOn", [ctypes.c_int, ctypes.c_int], ctypes.c_int),
                      ("GetOutputString", [ctypes.c_int], ctypes.c_char_p),
                      ("GetErrorString", [ctypes.c_int], ctypes.c_char_p),
                      ("RunString", [ctypes.c_int, ctypes.c_char_p], ctypes.c_int)):
    getattr(d, fn).argtypes, getattr(d, fn).restype = args, res
d.SetOutputStringOn(ip.id_, 1)
d.SetSelectedOutputFileOn(ip.id_, 1)
n_err = d.RunString(ip.id_, open(ex, encoding="latin-1").read().encode("latin-1"))
files = {f: open(f, encoding="latin-1").read() for f in sorted(os.listdir(".")) if f not in before and os.path.isfile(f)}
json.dump({"out": d.GetOutputString(ip.id_).decode("latin-1"), "files": files,
           "error": d.GetErrorString(ip.id_).decode("latin-1")[-2000:] if n_err else ""}, sys.stdout)
"""


def run_ours(example: Path, db: Path, timeout: int = 300) -> dict[str, Any]:
    """msl A4 경로: phreeqpython 이 싣는 IPhreeqc 라이브러리 (별도 프로세스 — 시간 제한)."""
    r = subprocess.run([sys.executable, "-c", _OURS, str(example), str(db)], capture_output=True, timeout=timeout,
                       check=False)
    try:
        return json.loads(r.stdout.decode("latin-1"))
    except json.JSONDecodeError:
        return {"out": "", "files": {}, "error": r.stderr.decode("latin-1")[-2000:]}


def _numbers(text: str) -> list[float | str]:
    out: list[float | str] = []
    for tok in text.split():
        try:
            out.append(float(tok))
        except ValueError:
            out.append(tok)
    return out


def compare_files(ref: dict[str, str], got: dict[str, str]) -> dict[str, Any]:
    """SELECTED_OUTPUT 파일: 같은 이름끼리 토큰 단위로 비교 (숫자는 상대오차)."""
    names = sorted(ref.keys() & got.keys())
    worst, cells, same_shape = 0.0, 0, True
    for n in names:
        a, b = _numbers(ref[n]), _numbers(got[n])
        same_shape &= len(a) == len(b)
        for x, y in zip(a, b):
            if isinstance(x, float) and isinstance(y, float):
                cells += 1
                worst = max(worst, abs(x - y) / max(abs(x), 1e-12) if abs(x) > 1e-12 else abs(y))
            elif x != y:
                same_shape = False
    return {"files": [sorted(ref), sorted(got)], "cells": cells, "max_rel": worst, "same_shape": same_shape}


def compare(ref: dict[str, list], got: dict[str, list]) -> dict[str, Any]:
    n = min(len(ref["solutions"]), len(got["solutions"]))
    d_ph = max((abs(a["pH"] - b["pH"]) for a, b in zip(ref["solutions"], got["solutions"]) if "pH" in a and "pH" in b), default=0.0)
    d_pe = max((abs(a["pe"] - b["pe"]) for a, b in zip(ref["solutions"], got["solutions"]) if "pe" in a and "pe" in b), default=0.0)
    d_i = max((abs(a["I"] - b["I"]) / max(abs(a["I"]), 1e-12) for a, b in zip(ref["solutions"], got["solutions"])
               if "I" in a and "I" in b), default=0.0)
    d_si = max((abs(ta[k] - tb[k]) for ta, tb in zip(ref["si"], got["si"]) for k in ta.keys() & tb.keys()
                if abs(ta[k]) <= 10 and abs(tb[k]) <= 10), default=0.0)  # |SI|>10 은 산화환원 기체 등 수치 잡음
    return {"solutions": [len(ref["solutions"]), len(got["solutions"])], "si_tables": [len(ref["si"]), len(got["si"])],
            "compared": n, "max_dpH": round(d_ph, 4), "max_dpe": round(d_pe, 4), "max_rel_dI": round(d_i, 5),
            "max_dSI": round(d_si, 4)}


def passes(c: dict[str, Any]) -> bool:
    """출력 표기 자릿수(pH·pe 소수 3자리, SI 2자리) 안에서 같으면 재현."""
    return (c["compared"] > 0 and c["solutions"][0] == c["solutions"][1] and c["max_dpH"] <= 0.002
            and c["max_dpe"] <= 0.002 and c["max_rel_dI"] <= 1e-3 and c["max_dSI"] <= 0.011)


def examples(dist: Path) -> list[Path]:
    key = lambda p: (int(re.match(r"ex(\d+)", p.name).group(1)), p.name)
    return sorted((p for p in (dist / "examples").iterdir() if re.fullmatch(r"ex\d+[a-z]?", p.name)), key=key)


def reference(ref_json: Path) -> dict[str, Any]:
    return json.loads(ref_json.read_text(encoding="utf-8"))


def bundled_db() -> Path:
    """msl A4 가 쓰는 데이터베이스 — phreeqpython 동봉 phreeqc.dat."""
    import phreeqpython

    return Path(phreeqpython.__file__).parent / "database" / "phreeqc.dat"


def database_for(name: str, dist: Path) -> tuple[Path, str]:
    """배포본 CMakeLists.txt 의 예제별 데이터베이스. phreeqc.dat 예제는 A4 와 같은 동봉 파일을 쓴다."""
    if name.startswith("ex15"):
        return dist / "examples" / "ex15.dat", "ex15.dat (배포본)"
    if name.startswith("ex17"):
        return dist / "database" / "pitzer.dat", "pitzer.dat (배포본)"
    if name.startswith("ex20"):
        return dist / "database" / "iso.dat", "iso.dat (배포본)"
    return bundled_db(), "phreeqc.dat (phreeqpython 동봉 = A4)"


def bench(dist: Path, binary: Path, only: list[str] | None = None,
          progress=None) -> dict[str, Any]:
    import hashlib

    rows = {}
    for ex in examples(dist):
        if only and ex.name not in only:
            continue
        if progress:
            progress(ex.name)
        db, db_label = database_for(ex.name, dist)
        try:
            off, ours = run_official(binary, ex, db), run_ours(ex, db)
        except subprocess.TimeoutExpired:
            rows[ex.name] = {"db": db_label, "error": "시간 초과"}
            continue
        ref, got = parse(off["out"]), parse(ours["out"])
        c = compare(ref, got)
        fc = compare_files(off["files"], ours["files"])
        text_ok = passes(c) if ref["solutions"] else None  # 용액 기술을 찍지 않는 예제는 파일로만 판정
        file_ok = (fc["cells"] > 0 and fc["same_shape"] and fc["max_rel"] <= 1e-6) if off["files"] else None
        verdicts = [v for v in (text_ok, file_ok) if v is not None]
        keep = len(ref["solutions"]) <= KEEP_VALUES  # 기준값 파일 크기 — 큰 예제는 비교 요약만
        row = {"db": db_label, "input_sha256": hashlib.sha256(ex.read_bytes()).hexdigest(), "official": ref if keep else None,
               "compare": c, "files": fc, "reproduced": bool(verdicts) and all(verdicts),
               "comparable": bool(verdicts), "error_ours": ours["error"], "error_official": off["error"]}
        if db == bundled_db():  # 참고: 배포본 3.8.6 phreeqc.dat 과의 데이터베이스 판 차이
            row["db_version_effect"] = compare(parse(run_official(binary, ex, dist / "database" / "phreeqc.dat")["out"]), ref)
        rows[ex.name] = row
    return {"phreeqc": dist.name, "binary": "공식 phreeqc (배포본 소스 빌드)", "engine": "phreeqpython IPhreeqc", "examples": rows}
