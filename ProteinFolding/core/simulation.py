"""OpenMM targeted MD:线性链 → 天然态的物理引导折叠动画。

流程:线性链加氢建拓扑 → AMBER14SB + GBn2 隐式溶剂 → 能量最小化 →
引导阶段(重原子受简谐外力,目标坐标沿线性插值从线性态连续形变到
天然态) → 终点松弛阶段(撤除外力自由 MD)。全程导出 DCD 轨迹与逐帧 RMSD。

数值稳定性说明:不采用"固定天然态目标 + 力常数爬升"的经典 TMD——
全延伸直棒链被反向折叠时,两段链会中途对撞导致 vdW 爆炸(NaN)。
本实现固定中等力常数、让目标本身按平滑时间表推进,原子沿近似路径行进,
数值上稳定得多。动画为物理引导示意,非真实折叠动力学。

openmm.unit 的模块属性(nanometer 等)为运行时动态创建,静态分析
不可见,故统一在导入处加 pyright 忽略。
"""

from __future__ import annotations

import io
import json
import logging
import time
from pathlib import Path
from typing import Any

import numpy as np

from .config import SimParams
from openmm.unit import femtoseconds, kelvin, nanometer, picosecond  # pyright: ignore[reportAttributeAccessIssue]

logger = logging.getLogger(__name__)

FloatArray = np.ndarray[tuple[int, ...], np.dtype[np.float64]]

FORCEFIELD_FILES = ("amber14-all.xml", "implicit/gbn2.xml")
NONBONDED_CUTOFF_NM = 1.6

# 自动平台优先级:OpenCL 在 Intel 核显上常慢于 CPU,不参与自动选择
AUTO_PLATFORM_ORDER = ("CUDA", "CPU")

# NaN 自动降档重试:力常数减半、步数加倍,最多重试次数
MAX_RETRY = 2


def available_platforms() -> list[str]:
    import openmm as mm

    return [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]


def choose_platform(preference: str = "auto") -> str:
    """preference ∈ {auto, CPU, CUDA, OpenCL, Reference}。"""
    have = available_platforms()
    if preference != "auto":
        if preference not in have:
            raise RuntimeError(
                f"平台 {preference} 不可用;当前可用: {have}。"
                "若需要 CUDA,请先安装 openmm[cuda13](见环境安装交接文档)。"
            )
        return preference
    for name in AUTO_PLATFORM_ORDER:
        if name in have:
            return name
    return have[0]


class _PullForce:
    """全部重原子的简谐引导力,目标坐标沿内坐标构象形变路径推进。

    目标不沿笛卡尔直线插值——那会让反向折叠的链段穿越自身、被键合几何
    卡死或顶爆——而是 MorphTarget 按 (φ,ψ,ω) 插值实时重建的可行构象:
    原子跟随一条键长键角始终理想的折叠路径。
    """

    def __init__(self, system: Any, topology: Any, positions: Any, morph: Any):
        import openmm as mm

        force = mm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        force.addGlobalParameter("k", 0.0)
        for p in ("x0", "y0", "z0"):
            force.addPerParticleParameter(p)

        self._morph = morph
        res_list = list(topology.residues())
        pos_nm = positions.value_in_unit(nanometer)
        assert len(pos_nm) == topology.getNumAtoms()

        # (残基序号, 原子名) → 拓扑原子序;同时保留初始坐标供 set_progress(0) 对齐
        atom_index: dict[tuple[int, str], int] = {}
        atoms_list = list(topology.atoms())
        for atom in atoms_list:
            if atom.element is None or atom.element.symbol == "H":
                continue
            key = (int(res_list[atom.residue.index].id), atom.name)
            atom_index[key] = atom.index

        coords = morph.coords_at(0.0)
        pulled: list[int] = []
        for res_idx, res_atoms in enumerate(coords, start=1):
            for name, _xyz in res_atoms:
                idx = atom_index.get((res_idx, name))
                if idx is None:
                    raise ValueError(
                        f"线性链与形变目标不一致:缺原子 ({res_idx}, {name})"
                    )
                force.addParticle(idx, list(pos_nm[idx]))
                pulled.append(idx)
        if not pulled:
            raise ValueError("未匹配到任何可引导原子:线性链与天然态拓扑不一致")
        system.addForce(force)
        self._force = force
        self.context: Any = None
        self.indices = pulled
        self._current_target: FloatArray = np.array(
            [pos_nm[i] for i in pulled], dtype=float
        )
        # 天然态终态目标(预计算,RMSD 监测用)
        final = morph.coords_at(1.0)
        self._final_target: FloatArray = np.array(
            [[xyz[0] * 0.1, xyz[1] * 0.1, xyz[2] * 0.1] for res in final for _, xyz in res]
        )
        # 形变旅程(平均位移,Å):用于按漂移速度自动预算引导步数
        self.journey_A = float(
            np.linalg.norm(self._final_target - self._current_target, axis=1).mean() * 10.0
        )
        logger.info("引导力覆盖重原子 %d 个", len(pulled))

    def bind(self, context: Any) -> None:
        self.context = context

    def set_progress(self, s: float) -> None:
        """s ∈ [0,1]:0 = 全延伸构象,1 = 天然态构象。"""
        assert self.context is not None, "先调用 bind()"
        s = min(max(s, 0.0), 1.0)
        coords = self._morph.coords_at(s)
        flat: list[list[float]] = []
        for res_atoms in coords:
            for _, xyz in res_atoms:
                flat.append([xyz[0] * 0.1, xyz[1] * 0.1, xyz[2] * 0.1])  # Å → nm
        self._current_target = np.array(flat)
        for entry, atom_index in enumerate(self.indices):
            self._force.setParticleParameters(
                entry, atom_index, list(self._current_target[entry])
            )
        self._force.updateParametersInContext(self.context)

    def current_target_nm(self) -> FloatArray:
        """当前引导目标坐标(nm,与拓扑原子一一对应的 pulled 子集)。"""
        return self._current_target.copy()

    def set_k(self, k: float) -> None:
        assert self.context is not None, "先调用 bind()"
        self.context.setParameter("k", k)

    def current_rmsd_to_target(self, positions_nm: FloatArray) -> float:
        """相对天然态终点的重原子 RMSD(Å);引导过程中用于展示收敛进度。"""
        diff = positions_nm[self.indices] - self._final_target
        return float(np.sqrt((diff**2).sum(axis=1).mean())) * 10.0  # nm → Å

    def gap_to_current_target(self, positions_nm: FloatArray) -> float:
        """相对当前引导目标的重原子 RMSD(nm);自适应步调的滞后测量。"""
        diff = positions_nm[self.indices] - self._current_target
        return float(np.sqrt((diff**2).sum(axis=1).mean()))


class _HTransport:
    """几何形变阶段的全原子坐标生成:重原子来自 MorphTarget(顺序一致),
    氢原子跟随其成键重原子、保持原始偏移向量(键长恒定,数值上不可能爆)。
    setPositions 直接驱动,无力传输。"""

    def __init__(self, topology: Any, positions: Any, morph: Any):
        res_list = list(topology.residues())
        pos_nm = positions.value_in_unit(nanometer)
        self._full0 = np.array(pos_nm, dtype=float)

        # (残基序号, 原子名) → 拓扑原子序
        atom_index: dict[tuple[int, str], int] = {}
        for atom in topology.atoms():
            if atom.element is None:
                continue
            atom_index[(int(res_list[atom.residue.index].id), atom.name)] = atom.index

        # 每个 H 的成键重原子伙伴与原始偏移向量
        self._h_follow: list[tuple[int, int, FloatArray]] = []
        for bond in topology.bonds():
            a1, a2 = bond
            h, partner = None, None
            if a1.element is not None and a1.element.symbol == "H":
                h, partner = a1.index, a2.index
            elif a2.element is not None and a2.element.symbol == "H":
                h, partner = a2.index, a1.index
            if h is None or partner is None:
                continue
            offset = self._full0[h] - self._full0[partner]
            self._h_follow.append((h, partner, offset))

        # 1-2 / 1-3 键合近邻集合(松弛时豁免;这些原子对本就该近)
        n = topology.getNumAtoms()
        excluded = np.zeros((n, n), dtype=bool)
        adj: dict[int, set[int]] = {i: set() for i in range(n)}
        for bond in topology.bonds():
            a, b = bond[0].index, bond[1].index
            excluded[a, b] = excluded[b, a] = True
            adj[a].add(b)
            adj[b].add(a)
        for i in range(n):
            for j in adj[i]:
                for k2 in adj[j]:
                    if k2 != i:
                        excluded[i, k2] = True
        self._excluded = excluded
        self._atom_index = atom_index
        self._morph = morph

    def _relax_overlaps(self, full: FloatArray, min_nm: float = 0.24, max_iter: int = 60) -> None:
        """把非键合且距离小于 min_nm 的原子对推开(纯几何,保键合几何)。

        形变中段链段互相穿越不可避免;全强度 LJ/GB 在近重合时发散,
        幽灵缩放也压不住 r^-12。几何上直接把这类原子对分摊推开,
        使后续短程 MD 数值稳定。每迭代位移封顶,防止推挤发散。
        """
        ex = self._excluded
        sq = (full * full).sum(axis=1)
        for _ in range(max_iter):
            d2 = sq[:, None] + sq[None, :] - 2.0 * (full @ full.T)
            bad = (d2 < min_nm * min_nm) & (~ex) & (d2 > 1e-18)
            if not bad.any():
                return
            ii, jj = np.nonzero(bad)
            diff = full[jj] - full[ii]
            dist = np.sqrt(np.maximum(d2[ii, jj], 0.0))[:, None]
            unit = diff / np.maximum(dist, 1e-9)
            push = (min_nm - dist) * 0.5
            disp = np.zeros_like(full)
            np.add.at(disp, ii, -unit * push)
            np.add.at(disp, jj, unit * push)
            # 单原子总位移封顶,抑制多对推挤的发散
            norm = np.linalg.norm(disp, axis=1, keepdims=True)
            scale = np.minimum(1.0, 0.05 / np.maximum(norm, 1e-12))
            candidate = full + disp * scale
            if not np.isfinite(candidate).all():
                return
            full = candidate
            sq = (full * full).sum(axis=1)

    def full_at(self, s: float) -> FloatArray:
        """s ∈ [0,1] 的全原子坐标(nm)。重原子取形变目标,氢按偏移随迁,
        随后做最小间距松弛防止穿越重合。"""
        coords = self._morph.coords_at(s)
        full = self._full0.copy()
        for res_idx, res_atoms in enumerate(coords, start=1):
            for name, xyz in res_atoms:
                full[self._atom_index[(res_idx, name)]] = np.array(
                    [xyz[0] * 0.1, xyz[1] * 0.1, xyz[2] * 0.1]
                )
        for h, partner, offset in self._h_follow:
            full[h] = full[partner] + offset
        self._relax_overlaps(full)
        return full


def _build_systems(
    linear_pdb_text: str, native_pdb_text: str
) -> tuple[Any, Any, Any, _PullForce, _HTransport]:
    """加氢、建双系统;返回 (modeller, morph_system, relax_system, pull, transport)。

    morph_system:不含 GB(GBn2 的 Born 半径链在原子近重合时产生 NaN,
    IEEE 的 0×NaN=NaN 使"零电荷"也无效)、非键电荷与 epsilon 全零——
    形变阶段只有键合力 + 温和引导约束,不存在数值爆散通道。
    relax_system:完整 AMBER14SB + GBn2,用于终点松弛的物理收敛。
    """
    import openmm as mm
    import openmm.app as app

    from .linear import MorphTarget

    pdb = app.PDBFile(io.StringIO(linear_pdb_text))
    forcefield = app.ForceField(*FORCEFIELD_FILES)
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield)

    morph_system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.CutoffNonPeriodic,  # pyright: ignore[reportArgumentType]
        nonbondedCutoff=NONBONDED_CUTOFF_NM * nanometer,
        constraints=app.HBonds,
    )
    # 移除 GB,非键电荷与 epsilon 全零
    for i in range(morph_system.getNumForces() - 1, -1, -1):
        force = morph_system.getForce(i)
        if isinstance(force, mm.CustomGBForce):
            morph_system.removeForce(i)
        elif isinstance(force, mm.NonbondedForce):
            for p in range(force.getNumParticles()):
                q, sigma, epsilon = force.getParticleParameters(p)
                force.setParticleParameters(
                    p, q * 0.0, sigma, epsilon * 0.0
                )

    morph = MorphTarget(native_pdb_text)
    pull = _PullForce(morph_system, modeller.topology, modeller.positions, morph)
    transport = _HTransport(modeller.topology, modeller.positions, morph)

    relax_system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.CutoffNonPeriodic,  # pyright: ignore[reportArgumentType]
        nonbondedCutoff=NONBONDED_CUTOFF_NM * nanometer,
        constraints=app.HBonds,
    )
    return modeller, morph_system, relax_system, pull, transport


def _simulate_once(
    modeller: Any,
    morph_system: Any,
    relax_system: Any,
    pull: _PullForce,
    transport: _HTransport,
    params: SimParams,
    out_dir: Path,
    platform_name: str,
    progress: Any,
) -> dict[str, Any]:
    """单次模拟尝试(参数已定),返回 meta dict。"""
    import openmm as mm
    import openmm.app as app
    from openmm.app import DCDFile

    def report(stage: str, frac: float, message: str = "") -> None:
        if progress:
            progress(stage, frac, message)

    seed = (
        params.seed
        if params.seed is not None
        else int(np.random.default_rng().integers(0, 2**31 - 1))
    )
    platform = mm.Platform.getPlatformByName(platform_name)
    # 形变阶段用 Verlet + 零速度:几何关键帧插值与高频 setPositions 组合下,
    # CPU 平台的 Langevin 随机速度 + 约束投影会偶发数值发散(实测);Verlet 无随机项,稳定。
    integrator = mm.VerletIntegrator(
        params.timestep_fs * femtoseconds  # pyright: ignore[reportOperatorIssue]
    )
    simulation = app.Simulation(modeller.topology, morph_system, integrator, platform)
    simulation.context.setPositions(modeller.positions)
    pull.bind(simulation.context)

    # 终点松弛用完整 GBn2 物理的独立系统
    relax_integrator = mm.LangevinMiddleIntegrator(
        params.temperature_K * kelvin,  # pyright: ignore[reportOperatorIssue]
        params.friction_per_ps / picosecond,  # pyright: ignore[reportOperatorIssue]
        params.timestep_fs * femtoseconds,  # pyright: ignore[reportOperatorIssue]
    )
    relax_integrator.setRandomNumberSeed(seed + 1)
    relax_simulation = app.Simulation(modeller.topology, relax_system, relax_integrator, platform)

    with open(out_dir / "topology.pdb", "w", encoding="ascii", newline="\n") as fh:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, fh, keepIds=True)

    total_steps = params.guide_steps + params.relax_steps
    interval = max(total_steps // max(params.n_frames, 1), 1)

    dcd_path = out_dir / "trajectory.dcd"
    with open(dcd_path, "wb") as dcd_fh:
        dcd = DCDFile(
            dcd_fh,
            modeller.topology,
            dt=params.timestep_fs * femtoseconds,  # pyright: ignore[reportOperatorIssue]
            interval=interval,
        )

        report("minimize", 0.0, "能量最小化中")
        t0 = time.perf_counter()
        simulation.minimizeEnergy(maxIterations=params.minimization_iters)
        minimize_s = time.perf_counter() - t0
        logger.info("能量最小化完成 (%.1fs)", minimize_s)


        rmsd_curve: list[float] = []

        def record(stage: str, frac: float, ctx: Any = None) -> None:
            source = ctx if ctx is not None else simulation
            pos = source.context.getState(getPositions=True).getPositions(asNumpy=True)
            dcd.writeModel(pos)
            rmsd = pull.current_rmsd_to_target(pos.value_in_unit(nanometer))
            rmsd_curve.append(round(rmsd, 2))
            report(stage, frac, f"第 {len(rmsd_curve)} 帧 RMSD={rmsd:.1f} Å")

        sim_t0 = time.perf_counter()

        # ---- 引导阶段:纯几何关键帧形变(无积分,不存在数值爆散路径) ----
        # 逐关键帧生成全原子坐标(MorphTarget 内坐标插值 + 最小间距松弛)
        # 并直接写入 context。物理收敛全部由末端松弛阶段承担。
        n_keyframes = max(60, params.n_frames // 2)
        sim_t0 = time.perf_counter()
        for kf in range(1, n_keyframes + 1):
            s = kf / n_keyframes
            simulation.context.setPositions(transport.full_at(s) * nanometer)
            record("guide", 0.45 + 0.45 * (kf / n_keyframes))
        guide_s = time.perf_counter() - sim_t0

        # ---- 终点松弛阶段:切换到完整 GBn2 物理系统,最小化 + 自由 MD ----
        relax_frames = max(params.n_frames - len(rmsd_curve), 10)
        interval = max(params.relax_steps // relax_frames, 1)
        final_nm = simulation.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(nanometer)
        relax_simulation.context.setPositions(final_nm * nanometer)
        relax_simulation.minimizeEnergy(maxIterations=params.minimization_iters)
        relax_t0 = time.perf_counter()
        done = 0
        while done < params.relax_steps:
            batch = min(interval, params.relax_steps - done)
            relax_simulation.step(batch)
            done += batch
            if done % interval == 0 or done >= params.relax_steps:
                record("relax", 0.90 + 0.10 * (done / params.relax_steps), relax_simulation)
        relax_s = time.perf_counter() - relax_t0
    sim_total_s = time.perf_counter() - sim_t0

    morph_steps = n_keyframes
    return {
        "platform": platform_name,
        "platforms_available": available_platforms(),
        "n_atoms": modeller.topology.getNumAtoms(),
        "total_steps": morph_steps + params.relax_steps,
        "morph_keyframes": n_keyframes,
        "relax_steps": params.relax_steps,
        "minimize_seconds": round(minimize_s, 2),
        "guide_seconds": round(guide_s, 2),
        "relax_seconds": round(relax_s, 2),
        "simulate_seconds": round(sim_total_s, 2),
        "seed": seed,
        "rmsd_A": rmsd_curve,
        "rmsd_start_A": rmsd_curve[0] if rmsd_curve else None,
        "rmsd_end_A": rmsd_curve[-1] if rmsd_curve else None,
    }


def run_targeted_md(
    linear_pdb_text: str,
    native_pdb_text: str,
    params: SimParams,
    out_dir: Path,
    platform_pref: str = "auto",
    progress: Any = None,
) -> dict[str, Any]:
    """执行完整模拟,out_dir 下产出 topology.pdb / trajectory.dcd / meta.json。"""
    import openmm as mm
    from dataclasses import replace as dc_replace

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    modeller, morph_system, relax_system, pull, transport = _build_systems(linear_pdb_text, native_pdb_text)
    platform_name = choose_platform(platform_pref)

    # 松弛步数按体系大小兜底(形变阶段为几何关键帧,与旅程长度无关)
    n_res = modeller.topology.getNumResidues()
    effective = dc_replace(
        params,
        relax_steps=max(params.relax_steps, n_res * 20),
    )
    if effective.relax_steps != params.relax_steps:
        logger.info(
            "按体系大小(残基 %d)调整松弛步数:%d→%d",
            n_res, params.relax_steps, effective.relax_steps,
        )

    # NaN/数值发散时自动降档重试(几何形变下基本不会触发,纯保险)
    attempt = effective
    retried = 0
    while True:
        try:
            meta = _simulate_once(
                modeller, morph_system, relax_system, pull, transport, attempt,
                out_dir, platform_name, progress,
            )
            break
        except mm.OpenMMException as exc:
            if retried >= MAX_RETRY or "NaN" not in str(exc):
                raise
            retried += 1
            attempt = _fallback_params(effective, retried)
            logger.warning(
                "模拟数值发散(第 %d 次),降档重试:k_max=%s, guide_steps=%s",
                retried, attempt.k_max, attempt.guide_steps,
            )

    meta["retries"] = retried
    meta["params"] = attempt.to_meta()
    meta["params_requested"] = params.to_meta()
    meta["n_residues"] = n_res
    meta["journey_A"] = round(pull.journey_A, 1)
    (out_dir / "meta.json").write_text(
        json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    if progress:
        progress("done", 1.0, "模拟完成")
    logger.info(
        "模拟完成:platform=%s, %d 帧, 引导 %.1fs + 松弛 %.1fs, RMSD %.1f→%.1f Å",
        platform_name,
        len(meta["rmsd_A"]),
        meta["guide_seconds"],
        meta["relax_seconds"],
        meta["rmsd_start_A"] or -1,
        meta["rmsd_end_A"] or -1,
    )
    return meta


def _fallback_params(params: SimParams, retry: int) -> SimParams:
    """第 retry 次重试的降档参数。"""
    from dataclasses import replace

    return replace(
        params,
        k_max=params.k_max / (2.0**retry),
        guide_steps=params.guide_steps * (2**retry),
        relax_steps=params.relax_steps,
    )
