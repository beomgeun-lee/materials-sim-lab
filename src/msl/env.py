"""프로젝트 경로와 .env 로딩."""

from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
KB_DIR = ROOT / "kb"
DATA_DIR = ROOT / "data"
CACHE_DIR = DATA_DIR / "cache"


def load_dotenv(path: Path = ROOT / ".env") -> None:
    """KEY=VALUE 형식의 .env 를 환경변수로 올린다. 이미 설정된 값은 덮어쓰지 않는다."""
    if not path.exists():
        return
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, value = line.partition("=")
        value = value.strip().strip('"').strip("'")
        if value:
            os.environ.setdefault(key.strip(), value)


def has_key(name: str) -> bool:
    load_dotenv()
    return bool(os.environ.get(name))
