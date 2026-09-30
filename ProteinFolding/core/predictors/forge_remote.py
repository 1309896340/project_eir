"""Forge 远程结构预测回退(EvolutionaryScale 官方 API)。

使用前提:在 https://forge.evolutionaryscale.ai 注册并获取 API key,
以环境变量 FORGE_API_KEY 提供。序列会上传至第三方服务,仅在本地权重
不可用或本地预测失败时启用。
"""

from __future__ import annotations

import logging
import os

from .base import PredictionError, StructurePredictor

logger = logging.getLogger(__name__)

FORGE_MODEL = "esm3-sm-open-v1"


class ForgeRemotePredictor(StructurePredictor):
    name = "forge-remote"
    description = "未设置 FORGE_API_KEY 环境变量"

    def __init__(self) -> None:
        self._token = os.environ.get("FORGE_API_KEY", "").strip()

    def is_available(self) -> bool:
        return bool(self._token)

    def predict_pdb(self, sequence: str) -> str:
        from esm.sdk.api import ESMProteinError
        from esm.sdk.forge import SequenceStructureForgeInferenceClient

        client = SequenceStructureForgeInferenceClient(
            model=FORGE_MODEL, token=self._token
        )
        logger.info("Forge fold 请求 (%d aa)", len(sequence))
        result = client.fold(sequence)
        if isinstance(result, ESMProteinError):
            raise PredictionError(f"Forge fold 失败: {result}")
        if getattr(result, "coordinates", None) is None:
            raise PredictionError("Forge fold 未返回结构坐标")
        return result.to_pdb_string()
