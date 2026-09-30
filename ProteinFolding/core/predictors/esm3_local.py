"""本地 ESM3 (esm3-sm-open-v1) 结构预测器。

权重由用户手动下载到 ProteinFolding/data/esm3/(目录结构与 HF 仓库一致),
运行时通过补丁把 esm SDK 的 data_root 指向该目录,避免触发
huggingface_hub 联网下载(见 CodonTransformer 的同类先例)。
"""

from __future__ import annotations

import logging
import threading
from pathlib import Path
from typing import cast

from .base import PredictionError, StructurePredictor

logger = logging.getLogger(__name__)

REQUIRED_FILES = (
    "config.json",
    "data/weights/esm3_sm_open_v1.pth",
    "data/weights/esm3_structure_encoder_v0.pth",
    "data/weights/esm3_structure_decoder_v0.pth",
    "data/weights/esm3_function_decoder_v0.pth",
)

_model = None
_model_lock = threading.Lock()
_model_device = ""


def weights_dir() -> Path:
    """ESM3 权重目录(与 HF 仓库 biohub/esm3-sm-open-v1 布局一致)。"""
    return Path(__file__).resolve().parents[2] / "data" / "esm3"


def _patch_data_root(target: Path) -> None:
    """把 esm SDK 的权重根目录重定向到本地,阻断联网下载。

    esm.pretrained 以 `from ... import data_root` 方式按名导入,
    因此两个绑定都要补丁。
    """
    import esm.pretrained as esm_pretrained
    import esm.utils.constants.esm3 as esm3_constants

    def local_data_root(model: str, _target: Path = target) -> Path:
        return _target

    esm3_constants.data_root = local_data_root
    esm_pretrained.data_root = local_data_root
    logger.info("esm data_root -> %s", target)


def _check_weights() -> tuple[bool, str]:
    d = weights_dir()
    missing = [f for f in REQUIRED_FILES if not (d / f).exists()]
    if missing:
        return False, "缺少权重文件: " + ", ".join(missing)
    return True, ""


def _load_model():
    global _model, _model_device
    with _model_lock:
        if _model is not None:
            return _model, _model_device
        import torch

        device = "cuda" if torch.cuda.is_available() else "cpu"
        ok, reason = _check_weights()
        if not ok:
            raise RuntimeError(reason)
        _patch_data_root(weights_dir())
        from esm.models.esm3 import ESM3

        logger.info("加载本地 ESM3 权重 (device=%s) ...", device)
        _model = ESM3.from_pretrained("esm3-sm-open-v1", device=torch.device(device))
        _model.eval()
        _model_device = device
        return _model, _model_device


class ESM3LocalPredictor(StructurePredictor):
    name = "esm3-local"
    description = "本地 esm3-sm-open-v1 权重未就绪"

    def is_available(self) -> bool:
        ok, _ = _check_weights()
        return ok

    def predict_pdb(self, sequence: str) -> str:
        import torch
        from esm.sdk.api import ESMProtein
        from esm.sdk.api import GenerationConfig

        model, device = _load_model()
        protein = ESMProtein(sequence=sequence)
        # esm 3.2.3 本地 generate 走 GenerationConfig(SamplingConfig 是 Forge
        # 远端专用,本地路径会因缺 condition_on_coordinates_only 属性崩溃);
        # track="structure" 即折叠,temperature=0 确定性解码
        config = GenerationConfig(track="structure", temperature=0.0)
        logger.info("ESM3 结构预测开始 (%d aa, device=%s)", len(sequence), device)
        try:
            with torch.no_grad():
                # SDK 注解只写了 GenerationConfig,但 ESMProtein 输入在
                # 运行时分发到 iterative_sampling_raw(官方示例的标准用法)
                output = model.generate(
                    protein,
                    config,
                )
        except RuntimeError as exc:
            raise PredictionError(f"ESM3 生成失败: {exc}") from exc
        result = cast(ESMProtein, output)
        if result.coordinates is None:
            raise PredictionError("ESM3 未返回结构坐标")
        pdb_text = result.to_pdb_string()
        logger.info("ESM3 结构预测完成 (%d 字符 PDB)", len(pdb_text))
        return pdb_text

    @staticmethod
    def loaded_device() -> str:
        return _model_device
