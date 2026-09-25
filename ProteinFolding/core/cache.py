"""结果缓存:同 序列+参数+管线版本 的模拟直接复用。"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from .config import PIPELINE_VERSION, SimParams

RESULT_FILES = ("native.pdb", "linear.pdb", "topology.pdb", "trajectory.dcd", "meta.json")


def data_root() -> Path:
    return Path(__file__).resolve().parents[1] / "data"


def cache_root() -> Path:
    return data_root() / "cache"


def cache_key(sequence: str, params: SimParams) -> str:
    payload = json.dumps(
        {"seq": sequence, "params": params.to_meta(), "v": PIPELINE_VERSION},
        sort_keys=True,
        ensure_ascii=False,
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]


def find_cached(key: str) -> Path | None:
    """缓存命中要求全部产物齐备,缺失视为未命中。"""
    d = cache_root() / key
    if all((d / f).exists() for f in RESULT_FILES):
        return d
    return None


def new_cache_dir(key: str) -> Path:
    d = cache_root() / key
    d.mkdir(parents=True, exist_ok=True)
    return d
