"""模拟参数与预设档位。

参数含义见 docs/物理方法说明.md。所有参数参与缓存 key 计算,
同一序列 + 同一参数 + 同一管线版本命中同一条缓存。
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from typing import Any

# 管线行为变化时递增,避免旧缓存被误复用
PIPELINE_VERSION = "2"

PRESET_KEYS = ("fast", "standard", "high")
PRESET_LABELS = {
    "fast": "快速预览",
    "standard": "标准",
    "high": "高质量",
}


@dataclass(frozen=True)
class SimParams:
    # 引导阶段(有外力)步数
    guide_steps: int = 3500
    # 终点松弛阶段(撤力自由 MD)步数
    relax_steps: int = 1500
    # 轨迹总帧数(两阶段合并,匀帧导出)
    n_frames: int = 250
    # 朗之万恒温器温度
    temperature_K: float = 300.0
    # 简谐引导力最大力常数 kJ/mol/nm^2
    k_max: float = 5000.0
    # 积分步长 fs(配合 H 键约束)
    timestep_fs: float = 2.0
    # 朗之万摩擦系数(1/ps);偏高的摩擦帮助耗散引导拉拽中的碰撞动能
    friction_per_ps: float = 2.0
    # 起点能量最小化迭代数
    minimization_iters: int = 500
    # 随机种子(None = 每次随机)
    seed: int | None = 42

    def with_overrides(self, overrides: dict[str, Any]) -> "SimParams":
        """高级面板的逐项覆盖;未知键与非法值直接拒绝。"""
        allowed: dict[str, type] = {
            f: type(getattr(self, f)) for f in SimParams.__dataclass_fields__
        }
        clean: dict[str, Any] = {}
        for k, v in (overrides or {}).items():
            if k not in allowed:
                raise ValueError(f"未知模拟参数: {k}")
            if v is None:
                continue
            if k == "seed":
                clean[k] = v
                continue
            want = allowed[k]
            # JSON 里的数字一律是 int;数值字段放宽 int → float(布尔除外)
            if isinstance(v, bool) or not isinstance(v, want):
                if want is float and isinstance(v, int):
                    clean[k] = float(v)
                    continue
                raise ValueError(f"参数 {k} 需要类型 {want.__name__}, 收到 {type(v).__name__}")
            clean[k] = v
        return replace(self, **clean)

    def to_meta(self) -> dict[str, Any]:
        return asdict(self)


PRESETS: dict[str, SimParams] = {
    "fast": SimParams(guide_steps=2000, relax_steps=800, n_frames=200, k_max=6000.0),
    "standard": SimParams(),
    "high": SimParams(guide_steps=5000, relax_steps=2500, n_frames=300, k_max=4000.0),
}


def resolve_params(preset: str, overrides: dict[str, Any] | None = None) -> SimParams:
    if preset not in PRESETS:
        raise ValueError(f"未知预设档位: {preset!r}, 可选 {PRESET_KEYS}")
    base = PRESETS[preset]
    return base.with_overrides(overrides or {})
