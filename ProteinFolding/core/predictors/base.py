"""结构预测器抽象接口。

约定:predict_pdb 输入无空白的标准氨基酸序列,返回单链 PDB 文本(重原子,
坐标为天然态预测结果)。预测器失败抛 PredictorUnavailable,由调用方决定
是否回退到下一个预测器。
"""

from __future__ import annotations

from abc import ABC, abstractmethod


class PredictorUnavailable(RuntimeError):
    """当前预测器不可用(缺权重、缺 key 等),可尝试回退。"""


class PredictionError(RuntimeError):
    """预测器可用但本次预测失败。"""


class StructurePredictor(ABC):
    name: str = "base"
    description: str = ""

    @abstractmethod
    def is_available(self) -> bool: ...

    @abstractmethod
    def predict_pdb(self, sequence: str) -> str: ...


def get_predictor_chain() -> list[StructurePredictor]:
    """按优先级返回预测器链:本地 ESM3 优先,Forge 远程兜底。"""
    from .esm3_local import ESM3LocalPredictor
    from .forge_remote import ForgeRemotePredictor

    return [ESM3LocalPredictor(), ForgeRemotePredictor()]


def predict_with_fallback(sequence: str) -> tuple[str, StructurePredictor]:
    """依次尝试预测器链,返回 (pdb 文本, 成功使用的预测器)。"""
    errors: list[str] = []
    for predictor in get_predictor_chain():
        if not predictor.is_available():
            errors.append(f"{predictor.name}: 不可用({predictor.description})")
            continue
        try:
            return predictor.predict_pdb(sequence), predictor
        except Exception as exc:  # noqa: BLE001 - 任一失败都记录并尝试下一个
            errors.append(f"{predictor.name}: 预测失败({exc})")
    raise PredictorUnavailable("所有预测器均失败:\n" + "\n".join(errors))
