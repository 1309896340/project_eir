"""全延伸线性肽链生成。

给定天然态 PDB(单链),重建一条 φ=-139°/ψ=+135°(β 折叠式全延伸)的
线性主链;每个残基的 CB 与侧链原子通过 Kabsch 叠合(N,CA,C 三点局部
参考系)从天然构象搬运过来,保留天然侧链内角与手性。O 按理想肽平面
几何放置,C 端缺 OXT 时补齐(amber14 力场模板必需)。

输出重原子 PDB,后续由 OpenMM Modeller 加氢并构建模拟拓扑。
"""

from __future__ import annotations

import io
import logging
from typing import Any

import numpy as np

logger = logging.getLogger(__name__)

# float64 多维数组别名(供严格类型检查用)
FloatArray = np.ndarray[tuple[int, ...], np.dtype[np.float64]]

# 全延伸(β-strand 型)主链内坐标;单位 Å / 度
PHI_EXT = -139.0
PSI_EXT = 135.0
OMEGA_EXT = 180.0

BOND_N_CA = 1.458
BOND_CA_C = 1.525
BOND_C_N = 1.329
ANGLE_N_CA_C = 111.2
ANGLE_CA_C_N = 116.2
ANGLE_C_N_CA = 121.7

STANDARD_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

# 解析结果:(残基名, [(原子名, 坐标)])
ParsedResidue = tuple[str, list[tuple[str, FloatArray]]]


def _place(
    a: FloatArray,
    b: FloatArray,
    c: FloatArray,
    bond: float,
    angle_deg: float,
    dihedral_deg: float,
) -> FloatArray:
    """NeRF 内坐标定位:给定 A,B,C 与 |CD|、∠BCD、二面角 ABCD,求 D。"""
    ang = np.radians(angle_deg)
    tor = np.radians(dihedral_deg)
    bc = c - b
    bc /= np.linalg.norm(bc)
    ab = b - a
    nv = np.cross(ab, bc)
    nv /= np.linalg.norm(nv)
    m = np.array([bc, np.cross(nv, bc), nv]).T
    d2 = np.array(
        [
            -bond * np.cos(ang),
            bond * np.sin(ang) * np.cos(tor),
            bond * np.sin(ang) * np.sin(tor),
        ]
    )
    return c + m @ d2


def _kabsch(P: FloatArray, Q: FloatArray) -> tuple[FloatArray, FloatArray]:
    """求把 P 叠合到 Q 的旋转矩阵与平移(det=+1,保手性)。"""
    p_mean = P.mean(axis=0)
    q_mean = Q.mean(axis=0)
    cov = (P - p_mean).T @ (Q - q_mean)
    u, _, vt = np.linalg.svd(cov)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return rot, q_mean - rot @ p_mean


def _parse_native(pdb_text: str) -> list[ParsedResidue]:
    """解析天然态 PDB,取最长标准残基链的重原子坐标。"""
    from Bio.PDB import PDBParser  # pyright: ignore[reportPrivateImportUsage]

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("native", io.StringIO(pdb_text))
    assert structure is not None, "PDB 解析失败"
    model = next(structure.get_models())

    best_residues: list[Any] = []
    for chain in model:
        residues = [
            r
            for r in chain
            if r.get_id()[0] == " " and r.get_resname() in STANDARD_RESIDUES
        ]
        if len(residues) > len(best_residues):
            best_residues = residues
    if not best_residues:
        raise ValueError("PDB 中未找到标准氨基酸残基")

    parsed: list[ParsedResidue] = []
    for res in best_residues:
        atoms: list[tuple[str, FloatArray]] = []
        for atom in res:
            if atom.element == "H":
                continue
            atoms.append((atom.get_name(), atom.get_coord().astype(float)))
        parsed.append((res.get_resname(), atoms))
    return parsed


def build_linear(native_pdb_text: str) -> tuple[str, list[str]]:
    """返回 (线性链 PDB 文本, 序列列表)。"""
    native = _parse_native(native_pdb_text)
    n_res = len(native)
    if n_res < 2:
        raise ValueError("残基数过少,无法构建线性链")

    # ---- 重建延伸主链 N/CA/C ----
    backbone: list[dict[str, FloatArray]] = [dict() for _ in range(n_res)]
    n0: FloatArray = np.array([0.0, 0.0, 0.0])
    ca0: FloatArray = np.array([BOND_N_CA, 0.0, 0.0])
    ca_n_dir = (n0 - ca0) / np.linalg.norm(n0 - ca0)
    ang = np.radians(ANGLE_N_CA_C)
    ca_c_dir = np.array(
        [
            ca_n_dir[0] * np.cos(ang) - ca_n_dir[1] * np.sin(ang),
            ca_n_dir[0] * np.sin(ang) + ca_n_dir[1] * np.cos(ang),
            0.0,
        ]
    )
    c0 = ca0 + BOND_CA_C * ca_c_dir
    backbone[0]["N"] = n0
    backbone[0]["CA"] = ca0
    backbone[0]["C"] = c0

    for i in range(1, n_res):
        n_i = _place(
            backbone[i - 1]["N"], backbone[i - 1]["CA"], backbone[i - 1]["C"],
            BOND_C_N, ANGLE_CA_C_N, PSI_EXT,
        )
        ca_i = _place(
            backbone[i - 1]["CA"], backbone[i - 1]["C"], n_i,
            BOND_N_CA, ANGLE_C_N_CA, OMEGA_EXT,
        )
        c_i = _place(
            backbone[i - 1]["C"], n_i, ca_i,
            BOND_CA_C, ANGLE_N_CA_C, PHI_EXT,
        )
        backbone[i]["N"] = n_i
        backbone[i]["CA"] = ca_i
        backbone[i]["C"] = c_i

    # ---- 残基内原子:CB/侧链按 (N,CA,C) 局部参考系从天然构象搬运 ----
    # O 不参与叠合:肽基碳为 sp2 平面,O 在 (CA, C, N(i+1)) 平面内、
    # 绕 C-N(i+1) 肽键轴与 CA 反式(ω=180)。绕 CA-C 轴放置会与 N(i+1) 冲突。
    out_residues: list[ParsedResidue] = []
    for i, (resname, atoms) in enumerate(native):
        pos = {name: xyz for name, xyz in atoms}
        if i < n_res - 1:
            o_i = _place(
                backbone[i]["CA"], backbone[i + 1]["N"], backbone[i]["C"],
                1.231, 123.0, 180.0,
            )
        else:
            # 末残基无 N(i+1),退回 N-CA-C-O = 180(无冲突对象,羧基另补 OXT)
            o_i = _place(
                backbone[i]["N"], backbone[i]["CA"], backbone[i]["C"],
                1.231, 120.5, 180.0,
            )
        linear_res: list[tuple[str, FloatArray]] = [
            ("N", backbone[i]["N"]), ("CA", backbone[i]["CA"]),
            ("C", backbone[i]["C"]), ("O", o_i),
        ]
        fixed = [name for name, _ in atoms]
        moving = [name for name in fixed if name not in ("N", "CA", "C", "O", "OXT")]
        if moving:
            P = np.array([pos[name] for name in moving])
            ref_nat = np.array([pos[name] for name in ("N", "CA", "C")])
            ref_lin = np.array(
                [backbone[i]["N"], backbone[i]["CA"], backbone[i]["C"]]
            )
            rot, trans = _kabsch(ref_nat, ref_lin)
            moved = (rot @ P.T).T + trans
            linear_res.extend(zip(moving, moved))
        out_residues.append((resname, linear_res))

    _ensure_c_terminal_oxt(out_residues)

    pdb_text = _write_pdb(out_residues)
    sequence = [resname for resname, _ in out_residues]
    return pdb_text, sequence


def _ensure_c_terminal_oxt(residues: list[ParsedResidue]) -> None:
    """末残基缺 OXT 时按标准羧基几何补一个(amber14 C 端模板必需)。"""
    atoms = dict(residues[-1][1])
    if "OXT" in atoms:
        return
    if not {"N", "CA", "C", "O"} <= atoms.keys():
        raise ValueError("末残基骨架原子不全,无法补 OXT")
    # 二面角 N-CA-C-OXT = 0(与 O 成对:O 在 180,两者夹角约 126°)
    oxt = _place(atoms["N"], atoms["CA"], atoms["C"], 1.25, 116.6, 0.0)
    residues[-1][1].append(("OXT", oxt))


def make_linear_pdb(native_pdb_text: str) -> tuple[str, list[str]]:
    """对外入口:天然态 PDB → 线性态 PDB。"""
    linear_text, sequence = build_linear(native_pdb_text)
    logger.info("线性链构建完成:%d 残基", len(sequence))
    return linear_text, sequence


# ---------------- 构象形变(morph)目标 ----------------


def _dihedral(p0: FloatArray, p1: FloatArray, p2: FloatArray, p3: FloatArray) -> float:
    """四个点的二面角(度,弧度域 [-180,180])。"""
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    x = np.dot(v, w)
    y = np.dot(np.cross(b1n, v), w)
    return float(np.degrees(np.arctan2(y, x)))


def _nearest_angle_lerp(a: float, b: float, t: float) -> float:
    """角度插值,走最短弧(处理 ±180 环绕)。"""
    d = (b - a + 180.0) % 360.0 - 180.0
    return a + d * t


class MorphTarget:
    """内坐标构象形变目标:s=0 全延伸,s=1 天然态二面角。

    拉力目标不沿笛卡尔直线走——那会让链段穿越自身、被键合几何卡死——
    而是按 (φ,ψ,ω) 从全延伸值插值到天然值,用 NeRF 实时重建主链坐标,
    侧链按局部参考系叠合。任意中间态都是键长键角理想的可行构象,
    原子跟随的是一条真实可折叠的路径。
    """

    def __init__(self, native_pdb_text: str):
        self._parsed = _parse_native(native_pdb_text)
        self.n_res = len(self._parsed)
        self.phi_nat, self.psi_nat, self.omega_nat = self._native_dihedrals()

        # 预存每个残基的侧链原子(相对 N,CA,C 的叠合源)
        self._moving: list[tuple[list[str], FloatArray]] = []
        for _resname, atoms in self._parsed:
            pos = {name: xyz for name, xyz in atoms}
            moving = [
                n for n, _ in atoms if n not in ("N", "CA", "C", "O", "OXT")
            ]
            xyz = np.array([pos[n] for n in moving]) if moving else np.zeros((0, 3))
            self._moving.append((moving, xyz))
        self._native_pos = [
            {name: xyz for name, xyz in atoms} for _, atoms in self._parsed
        ]
        self.sequence = [resname for resname, _ in self._parsed]

    def _native_dihedrals(self) -> tuple[list[float], list[float], list[float]]:
        """从天然态坐标提取 φ/ψ/ω;缺失(端点)用延伸值兜底。"""
        phis: list[float] = []
        psis: list[float] = []
        omegas: list[float] = []
        for i in range(self.n_res):
            cur = {name: xyz for name, xyz in self._parsed[i][1]}
            nxt = (
                {name: xyz for name, xyz in self._parsed[i + 1][1]}
                if i + 1 < self.n_res
                else None
            )
            phis_i, psi_i, omega_i = PHI_EXT, PSI_EXT, OMEGA_EXT
            if i > 0:
                prev = {name: xyz for name, xyz in self._parsed[i - 1][1]}
                if "C" in prev:
                    phis_i = _dihedral(prev["C"], cur["N"], cur["CA"], cur["C"])
            if nxt is not None and "O" in cur and "N" in nxt:
                psi_i = _dihedral(cur["N"], cur["CA"], cur["C"], nxt["N"])
                omega_i = _dihedral(cur["CA"], cur["C"], nxt["N"], nxt["CA"])
            phis.append(phis_i)
            psis.append(psi_i)
            omegas.append(omega_i)
        return phis, psis, omegas

    def coords_at(self, s: float) -> list[list[tuple[str, FloatArray]]]:
        """s ∈ [0,1] 的完整重原子坐标(与 make_linear_pdb 输出同构)。"""
        s = min(max(s, 0.0), 1.0)
        n_res = self.n_res
        phi = [
            _nearest_angle_lerp(PHI_EXT, self.phi_nat[i], s) for i in range(n_res)
        ]
        psi = [
            _nearest_angle_lerp(PSI_EXT, self.psi_nat[i], s) for i in range(n_res)
        ]
        omega = [
            _nearest_angle_lerp(OMEGA_EXT, self.omega_nat[i], s) for i in range(n_res)
        ]

        backbone: list[dict[str, FloatArray]] = [dict() for _ in range(n_res)]
        n0 = np.array([0.0, 0.0, 0.0])
        ca0 = np.array([BOND_N_CA, 0.0, 0.0])
        ca_n_dir = (n0 - ca0) / np.linalg.norm(n0 - ca0)
        ang = np.radians(ANGLE_N_CA_C)
        ca_c_dir = np.array(
            [
                ca_n_dir[0] * np.cos(ang) - ca_n_dir[1] * np.sin(ang),
                ca_n_dir[0] * np.sin(ang) + ca_n_dir[1] * np.cos(ang),
                0.0,
            ]
        )
        backbone[0]["N"] = n0
        backbone[0]["CA"] = ca0
        backbone[0]["C"] = ca0 + BOND_CA_C * ca_c_dir

        for i in range(1, n_res):
            n_i = _place(
                backbone[i - 1]["N"], backbone[i - 1]["CA"], backbone[i - 1]["C"],
                BOND_C_N, ANGLE_CA_C_N, psi[i - 1],
            )
            ca_i = _place(
                backbone[i - 1]["CA"], backbone[i - 1]["C"], n_i,
                BOND_N_CA, ANGLE_C_N_CA, omega[i - 1],
            )
            c_i = _place(
                backbone[i - 1]["C"], n_i, ca_i,
                BOND_CA_C, ANGLE_N_CA_C, phi[i],
            )
            backbone[i]["N"] = n_i
            backbone[i]["CA"] = ca_i
            backbone[i]["C"] = c_i

        out: list[list[tuple[str, FloatArray]]] = []
        for i in range(n_res):
            cur = backbone[i]
            if i < n_res - 1:
                o_i = _place(cur["CA"], backbone[i + 1]["N"], cur["C"], 1.231, 123.0, 180.0)
            else:
                o_i = _place(cur["N"], cur["CA"], cur["C"], 1.231, 120.5, 180.0)
            res: list[tuple[str, FloatArray]] = [
                ("N", cur["N"]), ("CA", cur["CA"]), ("C", cur["C"]), ("O", o_i),
            ]
            moving, native_xyz = self._moving[i]
            if moving:
                nat = self._native_pos[i]
                ref_nat = np.array([nat["N"], nat["CA"], nat["C"]])
                ref_cur = np.array([cur["N"], cur["CA"], cur["C"]])
                rot, trans = _kabsch(ref_nat, ref_cur)
                moved = (rot @ native_xyz.T).T + trans
                res.extend(zip(moving, moved))
            if i == n_res - 1 and "OXT" in self._native_pos[i]:
                # 末残基 OXT:天然态有就叠合,没有按理想几何补
                res.append(("OXT", self._oxt(cur)))
            elif i == n_res - 1:
                res.append(("OXT", self._oxt(cur)))
            out.append(res)
        return out

    @staticmethod
    def _oxt(cur: dict[str, FloatArray]) -> FloatArray:
        return _place(cur["N"], cur["CA"], cur["C"], 1.25, 116.6, 0.0)


def _write_pdb(residues: list[ParsedResidue]) -> str:
    lines = []
    serial = 1
    for res_idx, (resname, atoms) in enumerate(residues, start=1):
        for atom_name, xyz in atoms:
            elem = _element_of(atom_name)
            name_field = f" {atom_name:<3}" if len(atom_name) < 4 else atom_name
            lines.append(
                f"ATOM  {serial:5d} {name_field}{' ':1s}{resname:>3s} A{res_idx:4d}    "
                f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1.00:6.2f}{0.00:6.2f}"
                f"          {elem:>2s}  "
            )
            serial += 1
    lines.append("TER")
    lines.append("END")
    return "\n".join(lines) + "\n"


def _element_of(atom_name: str) -> str:
    """标准残基原子名首字母即元素(C/N/O/S/H);H 已在上游剔除。"""
    return atom_name.strip()[:1]
