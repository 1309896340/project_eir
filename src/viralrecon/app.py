"""ViralRecon 演示服务入口:FastAPI 后端 + 原生前端(无外部 JS 依赖)。

运行:
    uv run python -m viralrecon.app    # 或:uv run viralrecon-app
浏览器自动打开 http://127.0.0.1:7862 。计算为纯模拟、秒级完成,
POST /api/run 同步返回完整报告。

流程概念与工具映射(演示为教学模拟,生产请用真实 nf-core/viralrecon):
    模拟 FASTQ       <- Illumina 风格 reads(质量衰减 + 低质量 reads)
    质控统计         <- FastQC / MultiQC
    理想化比对+pileup<- bowtie2 + samtools mpileup
    变异检测         <- ivar variants(depth/qual/频率过滤)
    共有序列         <- ivar consensus(多数碱基,min_frac 阈值)
"""

from __future__ import annotations

import logging
import os
import threading
import webbrowser
from pathlib import Path
from typing import Any

from fastapi import FastAPI, HTTPException
from fastapi.staticfiles import StaticFiles

from .core.pipeline import run_pipeline
from .core.simulate import SimParams
from .core.variants import VariantFilters

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s: %(message)s")
logger = logging.getLogger("viralrecon.app")

ROOT = Path(__file__).resolve().parent
FRONTEND = ROOT / "frontend"

app = FastAPI(title="ViralRecon Demo", version="1.0.0")

# ---------------- 参数预设与校验 ----------------

PRESETS: dict[str, SimParams] = {
    "quick": SimParams(          # 最小规模,毫秒级
        genome_len=2000, coverage=60, read_len=100,
        n_snps=10, n_indels=3, coverage_mode="amplicon", seed=42,
    ),
    "standard": SimParams(       # 默认教学档
        genome_len=5000, coverage=100, read_len=150,
        n_snps=20, n_indels=6, coverage_mode="amplicon", seed=42,
    ),
    "coronavirus": SimParams(    # 冠状病毒基因组规模(~30kb)
        genome_len=30000, coverage=120, read_len=150,
        n_snps=40, n_indels=12, coverage_mode="amplicon", seed=42,
    ),
}
PRESET_LABELS = {
    "quick": "快速预览(2kb · 60x)",
    "standard": "标准教学(5kb · 100x)",
    "coronavirus": "冠状病毒规模(30kb · 120x)",
}

# 参数覆盖白名单与取值范围(clamp 到边界)
_PARAM_RANGES: dict[str, tuple[float, float]] = {
    "genome_len": (500, 30000),
    "coverage": (5, 200),
    "read_len": (50, 250),
    "error_rate_scale": (0.0, 20.0),
    "n_snps": (0, 100),
    "n_indels": (0, 50),
    "indel_max_len": (1, 10),
    "low_quality_frac": (0.0, 0.5),
    "q_start": (20, 41),
    "q_end": (10, 41),
    "seed": (0, 2**31 - 1),
}
_INT_PARAMS = {
    "genome_len", "coverage", "read_len", "n_snps", "n_indels",
    "indel_max_len", "seed",
}
# 总碱基预算,防止参数组合导致同步计算过久
_BASE_BUDGET = 4_500_000


def resolve_params(preset: str, overrides: dict[str, Any],
                    seed: int | None) -> SimParams:
    """从预设出发应用覆盖参数,类型/范围校验,超预算直接 400。"""
    base = PRESETS.get(preset)
    if base is None:
        raise ValueError(f"未知预设:{preset}(可选 {'/'.join(PRESETS)})")
    values: dict[str, Any] = {
        f.name: getattr(base, f.name)
        for f in base.__dataclass_fields__.values()
    }
    for key, raw in overrides.items():
        if key == "coverage_mode":
            if raw not in ("uniform", "amplicon"):
                raise ValueError("coverage_mode 仅支持 uniform/amplicon")
            values["coverage_mode"] = str(raw)
            continue
        if key not in _PARAM_RANGES:
            raise ValueError(f"未知参数:{key}")
        try:
            val = float(raw)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"参数 {key} 不是数值:{raw!r}") from exc
        if val != val:  # NaN
            raise ValueError(f"参数 {key} 非法(NaN)")
        lo, hi = _PARAM_RANGES[key]
        val = min(hi, max(lo, val))
        values[key] = int(val) if key in _INT_PARAMS else val
    if seed is not None:
        values["seed"] = int(min(2**31 - 1, max(0, seed)))
    params = SimParams(**values)
    budget = params.genome_len * params.coverage
    if budget > _BASE_BUDGET:
        raise ValueError(
            f"计算预算超限:基因组长度×覆盖度 = {budget:,} 碱基"
            f"(上限 {_BASE_BUDGET:,});请调小 genome_len 或 coverage"
        )
    return params


# ---------------- 流程元信息(前端流程条) ----------------

PIPELINE_STEPS = [
    {"key": "fastq", "label": "模拟 FASTQ", "tool": "Illumina 风格 reads",
     "desc": "按覆盖度采样 reads,质量随 cycle 衰减,含少量低质量 reads"},
    {"key": "qc", "label": "质控统计", "tool": "FastQC / MultiQC",
     "desc": "per-cycle 质量分布、GC 含量、碱基总量"},
    {"key": "align", "label": "比对 + pileup", "tool": "bowtie2 · samtools",
     "desc": "reads 映射回参考基因组,统计逐位点覆盖深度(演示为理想化比对)"},
    {"key": "variants", "label": "变异检测", "tool": "ivar variants",
     "desc": "SNP/Indel 位点与频率,depth/qual/AF 三重过滤"},
    {"key": "consensus", "label": "共有序列", "tool": "ivar consensus",
     "desc": "多数碱基调用,低置信位写 N,高频缺失移除"},
]


def _public_params(p: SimParams) -> dict[str, Any]:
    return {f.name: getattr(p, f.name) for f in p.__dataclass_fields__.values()}


@app.get("/api/presets")
def presets() -> dict[str, Any]:
    return {
        "presets": [
            {"key": k, "label": PRESET_LABELS[k],
             "params": _public_params(PRESETS[k])}
            for k in PRESETS
        ],
        "param_ranges": {k: list(v) for k, v in _PARAM_RANGES.items()},
        "base_budget": _BASE_BUDGET,
    }


@app.get("/api/pipeline")
def pipeline_steps() -> list[dict[str, str]]:
    return PIPELINE_STEPS


@app.post("/api/run")
def run(body: dict[str, Any]) -> dict[str, Any]:
    preset = str(body.get("preset", "standard"))
    overrides = body.get("params") or {}
    if not isinstance(overrides, dict):
        raise HTTPException(status_code=400, detail="params 须为对象")
    seed = body.get("seed")
    try:
        params = resolve_params(preset, overrides, seed)
        filters = VariantFilters()
        report = run_pipeline(params, filters)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    logger.info("完成模拟: %s (seed=%s, %.2fs)", preset, params.seed,
                float(report["timings"]["total_s"]))
    return report


app.mount("/", StaticFiles(directory=FRONTEND, html=True), name="frontend")


def main() -> None:
    import uvicorn

    url = "http://127.0.0.1:7862"
    # 调试器(VS Code F5)下由 serverReadyAction 开浏览器,避免开两个标签页
    if not os.environ.get("VIRALRECON_NO_BROWSER"):
        threading.Timer(1.2, lambda: webbrowser.open(url)).start()
    logger.info("ViralRecon 演示服务启动: %s", url)
    uvicorn.run(app, host="127.0.0.1", port=7862, log_level="warning")


if __name__ == "__main__":
    main()
