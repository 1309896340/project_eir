"""命令行端到端:不启动浏览器,直接跑完 预测→线性化→模拟 并落盘。

示例:
    uv run python -m ProteinFolding.demo --sample ins --preset standard
    uv run python -m ProteinFolding.demo --native-pdb path/to/structure.pdb
    uv run python -m ProteinFolding.demo --sequence-file protein.fa --preset fast
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from .core.config import PRESET_KEYS, resolve_params
from .core.pipeline import run_pipeline
from .core.samples import SAMPLES_BY_KEY
from .core.sequence import clean_and_check

logger = logging.getLogger("proteinfolding.demo")

STAGE_LABELS = {
    "predict": "结构预测", "linear": "线性链构建", "minimize": "能量最小化",
    "guide": "引导折叠", "relax": "终点松弛", "done": "完成",
}


def main() -> int:
    parser = argparse.ArgumentParser(description="蛋白质折叠动画命令行管线")
    src = parser.add_mutually_exclusive_group(required=False)
    src.add_argument("--sequence", help="氨基酸序列(裸序列或 FASTA)")
    src.add_argument("--sequence-file", help="序列文件路径")
    src.add_argument("--sample", choices=sorted(SAMPLES_BY_KEY), help="内置样例")
    parser.add_argument(
        "--native-pdb", type=Path, default=None,
        help="跳过预测,直接给定天然态 PDB(可与序列来源组合,也可单独使用)",
    )
    parser.add_argument("--preset", default="standard", choices=PRESET_KEYS)
    parser.add_argument("--platform", default="auto", help="auto/CPU/CUDA/OpenCL")
    parser.add_argument("--force", action="store_true", help="忽略缓存强制重算")
    parser.add_argument("--params", default="{}", help="高级参数 JSON,如 '{\"k_max\": 3000}'")
    parser.add_argument("--out", type=Path, default=None, help="输出目录(默认走缓存)")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s: %(message)s")

    if args.sequence_file:
        sequence = Path(args.sequence_file).read_text(encoding="utf-8")
    elif args.sample:
        sequence = SAMPLES_BY_KEY[args.sample].sequence
    elif args.sequence:
        sequence = args.sequence
    else:
        sequence = ""

    try:
        params = resolve_params(args.preset, json.loads(args.params))
    except (ValueError, json.JSONDecodeError) as exc:
        logger.error("参数错误: %s", exc)
        return 2

    native_pdb_text = args.native_pdb.read_text(encoding="utf-8") if args.native_pdb else None
    check = clean_and_check(sequence)
    if native_pdb_text is None and not check.ok:
        logger.error("序列不合法: %s", ";".join(check.errors))
        return 2

    def progress(stage: str, frac: float, message: str) -> None:
        bar = "#" * int(frac * 30)
        print(f"\r[{bar:<30}] {int(frac*100):3d}% {STAGE_LABELS.get(stage, stage)} · {message}",
              end="", flush=True)
        if frac >= 1.0:
            print()

    result = run_pipeline(
        sequence=sequence,
        params=params,
        platform_pref=args.platform,
        force_refresh=args.force,
        progress=progress,
        out_dir=args.out,
        native_pdb_text=native_pdb_text,
    )

    meta = result["meta"] if isinstance(result["meta"], dict) else json.loads(result["meta"])
    print(f"\n输出目录 : {result['dir']}")
    print(f"缓存     : {'命中' if result['cached'] else '新计算'}")
    print(f"平台     : {meta['platform']}  原子数: {meta['n_atoms']}  帧数: {len(meta['rmsd_A'])}")
    print(f"RMSD     : {meta['rmsd_start_A']} Å → {meta['rmsd_end_A']} Å")
    print("产物     : native.pdb / linear.pdb / topology.pdb / trajectory.dcd / meta.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
