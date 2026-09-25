"""预测结构的侧链/缺失原子修复。

ESM3 的结构输出可能只含主链 N/CA/C(缺失原子以 inf 标记,经 to_pdb_string
后被丢弃)。amber14 力场模板要求完整侧链,因此用 PDBFixer 按残基模板
重建缺失的重原子(氢交给 OpenMM Modeller,不在此处理)。
"""

from __future__ import annotations

import io
import logging

logger = logging.getLogger(__name__)


def ensure_sidechains(pdb_text: str) -> str:
    """重建 PDB 中缺失的侧链重原子;输入原子齐全时原样返回(幂等)。"""
    from pdbfixer import PDBFixer

    fixer = PDBFixer(pdbfile=io.StringIO(pdb_text))
    # 不允许插入缺失残基(链编号连续,无需补链);只补残基内缺失原子
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    n_missing = sum(len(v) for v in fixer.missingAtoms.values())
    if n_missing == 0 and not any(fixer.missingTerminals.values()):
        return pdb_text
    fixer.addMissingAtoms()
    out = io.StringIO()
    import openmm.app as app

    app.PDBFile.writeFile(fixer.topology, fixer.positions, out, keepIds=True)
    logger.info("PDBFixer 重建了 %d 个缺失重原子", n_missing)
    return out.getvalue()
