"""materials-sim-lab — 원소·광물·화학물질 조합 시뮬레이션 프레임워크."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("msl")
except PackageNotFoundError:  # 설치하지 않고 소스에서 실행할 때
    __version__ = "0.0.0+local"
