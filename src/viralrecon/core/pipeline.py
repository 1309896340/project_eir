"""端到端管线:合成数据 -> reads 模拟 -> QC/pileup -> 变异 -> consensus。

对应 nf-core/viralrecon 的主流程(工具映射见 app.py 的流程条):
FASTQ 模拟(FastQC 质控口径) -> 理想化比对(bowtie2 的简化) ->
ivar variants 风格变异检测 -> ivar consensus 风格共识调用。
全部计算为纯 Python/numpy,标准参数下亚秒级完成。
"""

from __future__ import annotations

import random
import time
from bisect import bisect_right
from typing import Any

import numpy as np

from .consensus import build_consensus, to_fasta
from .pileup import Pileup
from .qc import QcCollector
from .simulate import (
    BASE_IDX,
    SimParams,
    build_sample,
    sample_read_starts,
)
from .variants import VariantFilters, call_variants, score_against_truth

# 深度图下发前端的下采样点数上限
_MAX_DEPTH_POINTS = 4000
# Phred Q=2..41 的测序错误概率查表(p_err = 10^(-Q/10))
_ERR_TABLE = {q: 10 ** (-q / 10) for q in range(2, 42)}


def run_pipeline(params: SimParams,
                 filters: VariantFilters | None = None,
                 consensus_min_frac: float = 0.75) -> dict[str, Any]:
    """跑完整演示流程,返回 JSON 友好的报告字典。"""
    filters = filters or VariantFilters()
    t0 = time.perf_counter()

    rng = random.Random(params.seed)
    build = build_sample(params, rng)
    starts = sample_read_starts(params, len(build.sample_genome), rng)
    t_sim_start = time.perf_counter()

    pileup = Pileup(params.genome_len)
    qc = QcCollector(params.read_len)

    # ---- 单遍 reads 模拟:质量 -> 测序错误 -> pileup/QC 流式收集 ----
    genome = build.sample_genome
    ref_of = build.ref_of_sample
    # 插入碱基标记(构建后预处理,避免热循环里比较相邻坐标)
    is_ins = [False] * len(ref_of)
    for i in range(1, len(ref_of)):
        if ref_of[i] == ref_of[i - 1]:
            is_ins[i] = True
    seam_positions = [s for s, _ in build.del_seams]

    L = params.read_len
    q0, q1 = params.q_start, params.q_end
    q_slope = (q1 - q0) / max(1, L - 1)
    err_scale = params.error_rate_scale
    lowq_frac = params.low_quality_frac
    counts = pileup.counts
    qual_sums = pileup.qual_sums
    rand = rng.random
    gauss = rng.gauss
    err_table = _ERR_TABLE
    base_idx = BASE_IDX

    for s in starts:
        end = s + L
        seg = genome[s:end]
        bias = -12.0 if rand() < lowq_frac else 0.0
        quals: list[int] = []
        bases: list[str] = []
        for c, ch in enumerate(seg):
            q = q0 + q_slope * c + gauss(0.0, 1.2) + bias
            qi = 2 if q < 2 else (41 if q > 41 else int(q))
            quals.append(qi)
            if rand() < err_scale * err_table[qi]:
                # 测序错误:随机替换为其他三种碱基之一
                ch = "ACGT"[(base_idx[ch] + 1 + (1 if rand() < 0.5 else 0)) % 4]
            bases.append(ch)

        gc = 0
        ins_run = -1
        for c, ch in enumerate(bases):
            if ch == "G" or ch == "C":
                gc += 1
            sp = s + c
            if is_ins[sp]:
                if ins_run < 0:
                    ins_run = c
                continue
            if ins_run >= 0:
                pileup.add_insertion(
                    ref_of[s + ins_run], "".join(bases[ins_run:c])
                )
                ins_run = -1
            k = ref_of[sp] * 4 + base_idx[ch]
            counts[k] += 1
            qual_sums[k] += quals[c]
        if ins_run >= 0:  # read 末端恰为插入串
            pileup.add_insertion(
                ref_of[s + ins_run], "".join(bases[ins_run:])
            )

        # 跨过删除接缝的 reads:对 del 参考区间每位计缺失
        i = bisect_right(seam_positions, s)
        while i < len(seam_positions) and seam_positions[i] < end:
            v = build.del_seams[i][1]
            pileup.add_del_span(v.pos, len(v.ref))
            i += 1

        qc.add_read(quals, gc, len(bases))
    t_sim_end = time.perf_counter()

    # ---- 覆盖度统计 + 深度图下采样 ----
    depth = pileup.depth_array()
    n_pts = min(_MAX_DEPTH_POINTS, len(depth))
    step = -(-len(depth) // n_pts)  # ceil
    trimmed = len(depth) // step * step
    depth_ds = (
        depth[:trimmed].reshape(-1, step).mean(axis=1)
        if trimmed
        else depth.astype(np.float64)
    )
    coverage = {
        "mean": float(depth.mean()),
        "median": float(np.median(depth)),
        "breadth": float((depth > 0).mean()),          # 至少 1x 覆盖比例
        "pct_ge_10x": float((depth >= 10).mean()),
        "max": int(depth.max()),
        "depth": [round(float(x), 1) for x in depth_ds],
    }

    # ---- 变异检测 + 真值对照 ----
    called = call_variants(pileup, build.reference, filters)
    truth_score = score_against_truth(called, build.truth)
    variants = {
        **truth_score,
        "rows": [
            {
                "pos": c.pos, "ref": c.ref, "alt": c.alt, "kind": c.kind,
                "dp": c.depth, "af": round(c.af, 4),
                "qual": round(c.qual, 1), "filter": c.filter_status,
            }
            for c in called
        ],
        "truth_detail": [
            {"pos": v.pos + 1, "kind": v.kind, "ref": v.ref, "alt": v.alt}
            for v in build.truth
        ],
    }

    # ---- consensus ----
    cons = build_consensus(pileup, build.reference, consensus_min_frac)
    consensus = {
        "length": cons.length,
        "n_count": cons.n_count,
        "diffs_vs_ref": cons.diffs_vs_ref,
        "fasta": to_fasta("simulated_consensus", cons.seq),
    }

    return {
        "params": _params_dict(params),
        "timings": {
            "total_s": round(time.perf_counter() - t0, 3),
            "simulate_s": round(t_sim_end - t_sim_start, 3),
        },
        "qc": qc.report(),
        "coverage": coverage,
        "variants": variants,
        "consensus": consensus,
    }


def _params_dict(p: SimParams) -> dict[str, Any]:
    return {
        "genome_len": p.genome_len,
        "coverage": p.coverage,
        "read_len": p.read_len,
        "error_rate_scale": p.error_rate_scale,
        "n_snps": p.n_snps,
        "n_indels": p.n_indels,
        "low_quality_frac": p.low_quality_frac,
        "coverage_mode": p.coverage_mode,
        "q_start": p.q_start,
        "q_end": p.q_end,
        "seed": p.seed,
    }
