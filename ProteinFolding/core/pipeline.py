"""端到端管线编排:预测 → 线性化 → 模拟,带缓存与进度回调。"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from .cache import cache_key, find_cached, new_cache_dir
from .config import SimParams
from .linear import make_linear_pdb
from .pdbfix import ensure_sidechains
from .predictors.base import predict_with_fallback
from .simulation import run_targeted_md

logger = logging.getLogger(__name__)


def run_pipeline(
    sequence: str,
    params: SimParams,
    platform_pref: str = "auto",
    force_refresh: bool = False,
    progress: Any = None,
    out_dir: Path | None = None,
    native_pdb_text: str | None = None,
) -> dict[str, Any]:
    """执行完整管线,返回 {key, cached, dir, meta, predictor}。

    native_pdb_text:直接给定天然态 PDB 时跳过结构预测(折叠已知结构的
    蛋白,或测试/复现用途)。此时缓存 key 以该 PDB 内容为基准。

    progress(stage, frac, message) 的 stage 与分值约定见 app 前端:
    predict 0~0.30 / linear 0.30~0.45 / minimize 0.45~0.50 /
    guide 0.50~0.90 / relax 0.90~1.00。
    """
    from .sequence import clean_and_check

    def report(stage: str, frac: float, message: str = "") -> None:
        if progress:
            progress(stage, frac, message)

    check = clean_and_check(sequence)
    if native_pdb_text is None:
        if not check.ok:
            raise ValueError(";".join(check.errors))
    sequence = check.sequence

    if native_pdb_text is None:
        key = cache_key(sequence, params)
    else:
        key = cache_key(native_pdb_text, params)
    if not force_refresh:
        cached = find_cached(key)
        if cached:
            meta = (cached / "meta.json").read_text(encoding="utf-8")
            report("done", 1.0, "命中缓存")
            return {
                "key": key,
                "cached": True,
                "dir": str(cached),
                "meta": meta,
                "predictor": None,
            }

    work = new_cache_dir(key) if out_dir is None else Path(out_dir)
    work.mkdir(parents=True, exist_ok=True)

    # ---- 1. 结构预测 ----
    if native_pdb_text is not None:
        predictor_name = "given-pdb"
        native_pdb = native_pdb_text
        report("predict", 0.30, "使用给定的天然态 PDB")
    else:
        report("predict", 0.02, "结构预测中")
        native_pdb, predictor = predict_with_fallback(sequence)
        predictor_name = predictor.name
        report("predict", 0.30, f"结构预测完成({predictor.name})")
    # ESM3 可能输出主链-only(缺失原子 inf 标记),用 PDBFixer 重建侧链
    native_pdb = ensure_sidechains(native_pdb)
    (work / "native.pdb").write_text(native_pdb, encoding="ascii", newline="\n")

    # ---- 2. 线性链 ----
    report("linear", 0.32, "生成全延伸线性链")
    linear_pdb, _ = make_linear_pdb(native_pdb)
    (work / "linear.pdb").write_text(linear_pdb, encoding="ascii", newline="\n")
    report("linear", 0.45, "线性链就绪")

    # ---- 3. 模拟 ----
    meta = run_targeted_md(
        linear_pdb,
        native_pdb,
        params,
        work,
        platform_pref=platform_pref,
        progress=progress,
    )
    return {
        "key": key,
        "cached": False,
        "dir": str(work),
        "meta": meta,
        "predictor": predictor_name,
    }
