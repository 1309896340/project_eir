"""ProteinFolding 本地服务入口:FastAPI 后端 + Mol* 前端。

运行:
    uv run python -m ProteinFolding.app
浏览器自动打开 http://127.0.0.1:7861 。单任务串行执行(MD 为重计算),
新任务排队;结果按 序列+参数 缓存,重复提交瞬时返回。
"""

from __future__ import annotations

import logging
import queue
import threading
import uuid
import webbrowser
from pathlib import Path
from typing import Any

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

from .core.cache import cache_key, find_cached
from .core.config import PRESET_KEYS, PRESET_LABELS, resolve_params
from .core.pipeline import run_pipeline
from .core.predictors.base import get_predictor_chain
from .core.samples import SAMPLES
from .core.sequence import MIN_LENGTH, MAX_LENGTH, clean_and_check
from .core.simulation import available_platforms, choose_platform

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s: %(message)s")
logger = logging.getLogger("proteinfolding.app")

ROOT = Path(__file__).resolve().parent
FRONTEND = ROOT / "frontend"

app = FastAPI(title="ProteinFolding", version="1.0.0")

# ---------------- 任务管理(单工作线程 + 队列) ----------------

JOBS: dict[str, dict[str, Any]] = {}
_JOB_LOCK = threading.Lock()
_TASKQ: queue.Queue[str] = queue.Queue(maxsize=8)


def _worker() -> None:
    while True:
        job_id = _TASKQ.get()
        job = JOBS.get(job_id)
        if job is None:
            continue
        # 闭包内窄化失效,先绑定到非 Optional 局部变量
        job_dict: dict[str, Any] = job
        with _JOB_LOCK:
            job_dict["status"] = "running"
            job_dict["stage"] = "predict"
            job_dict["frac"] = 0.0
            job_dict["message"] = "排队结束,开始执行"

        def progress(stage: str, frac: float, message: str = "") -> None:
            with _JOB_LOCK:
                job_dict["stage"] = stage
                job_dict["frac"] = frac
                job_dict["message"] = message

        try:
            result = run_pipeline(
                sequence=job["sequence"],
                params=job["params"],
                platform_pref=job["platform"],
                force_refresh=job["force"],
                progress=progress,
                native_pdb_text=job["native_pdb"],
            )
            with _JOB_LOCK:
                job["status"] = "done"
                job["result"] = {
                    "key": result["key"],
                    "cached": result["cached"],
                    "predictor": result["predictor"],
                    "meta": result["meta"],
                }
                job["frac"] = 1.0
                job["message"] = "缓存命中" if result["cached"] else "完成"
        except Exception as exc:  # noqa: BLE001 - 面向前端转成错误信息
            logger.exception("任务 %s 失败", job_id)
            with _JOB_LOCK:
                job["status"] = "error"
                job["message"] = str(exc)


_threading_started = False


def _ensure_worker() -> None:
    global _threading_started
    if not _threading_started:
        threading.Thread(target=_worker, daemon=True, name="pf-worker").start()
        _threading_started = True


# ---------------- API ----------------


@app.get("/api/system")
def system_info() -> dict[str, Any]:
    from .core.predictors.esm3_local import weights_dir

    chain = get_predictor_chain()
    return {
        "platforms": available_platforms(),
        "platform_active": _safe_platform(),
        "length_range": [MIN_LENGTH, MAX_LENGTH],
        "predictors": [
            {"name": p.name, "available": p.is_available(), "description": p.description}
            for p in chain
        ],
        "weights_dir": str(weights_dir()),
    }


def _safe_platform() -> str:
    try:
        return choose_platform("auto")
    except Exception:  # noqa: BLE001
        return "unknown"


@app.get("/api/samples")
def samples() -> list[dict[str, Any]]:
    return [{"key": s.key, "label": s.label, "sequence": s.sequence} for s in SAMPLES]


@app.get("/api/presets")
def presets() -> dict[str, str]:
    return {k: PRESET_LABELS[k] for k in PRESET_KEYS}


@app.post("/api/validate")
def validate(body: dict[str, Any]) -> dict[str, Any]:
    check = clean_and_check(str(body.get("sequence", "")))
    return {
        "ok": check.ok,
        "length": len(check.sequence),
        "errors": check.errors,
        "warnings": check.warnings,
    }


@app.post("/api/jobs")
def submit_job(body: dict[str, Any]) -> dict[str, Any]:
    _ensure_worker()
    sequence = str(body.get("sequence", ""))
    check = clean_and_check(sequence)
    if not check.ok:
        raise HTTPException(status_code=400, detail=";".join(check.errors))

    native_pdb = body.get("native_pdb")  # 可选:直接给定天然态 PDB(调试/已知结构)
    if native_pdb is not None:
        native_pdb = str(native_pdb)

    try:
        params = resolve_params(
            str(body.get("preset", "standard")), body.get("params") or {}
        )
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc

    platform = str(body.get("platform", "auto"))
    force = bool(body.get("force", False))

    # 缓存直返:相同 序列+参数 直接复用,不占队列
    key = cache_key(native_pdb if native_pdb else check.sequence, params)
    if not force:
        cached = find_cached(key)
        if cached:
            job_id = uuid.uuid4().hex[:8]
            with _JOB_LOCK:
                JOBS[job_id] = {
                    "status": "done",
                    "stage": "done",
                    "frac": 1.0,
                    "message": "缓存命中",
                    "sequence": check.sequence,
                    "params": params,
                    "platform": platform,
                    "force": force,
                    "native_pdb": native_pdb,
                    "result": {"key": key, "cached": True, "predictor": None, "meta": None},
                }
            return {"job_id": job_id}

    job_id = uuid.uuid4().hex[:8]
    with _JOB_LOCK:
        JOBS[job_id] = {
            "status": "queued",
            "stage": "queued",
            "frac": 0.0,
            "message": "排队中",
            "sequence": check.sequence,
            "params": params,
            "platform": platform,
            "force": force,
            "native_pdb": native_pdb,
            "result": None,
        }
    try:
        _TASKQ.put_nowait(job_id)
    except queue.Full:
        with _JOB_LOCK:
            del JOBS[job_id]
        raise HTTPException(status_code=429, detail="任务队列已满,请稍后再试") from None
    return {"job_id": job_id}


@app.get("/api/jobs/{job_id}")
def job_status(job_id: str) -> JSONResponse:
    job = JOBS.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="任务不存在")
    with _JOB_LOCK:
        payload = {
            "status": job["status"],
            "stage": job["stage"],
            "frac": job["frac"],
            "message": job["message"],
            "result": job["result"],
        }
    return JSONResponse(payload)


RESULT_FILES = ("native.pdb", "linear.pdb", "topology.pdb", "trajectory.dcd", "meta.json")


@app.get("/api/results/{key}/{filename}")
def result_file(key: str, filename: str) -> FileResponse:
    if filename not in RESULT_FILES:
        raise HTTPException(status_code=404, detail="未知文件")
    # 防路径穿越
    if "/" in key or "\\" in key or ".." in key:
        raise HTTPException(status_code=400, detail="非法 key")
    path = ROOT / "data" / "cache" / key / filename
    if not path.exists():
        raise HTTPException(status_code=404, detail="文件不存在")
    media = "chemical/x-pdb" if filename.endswith(".pdb") else "application/octet-stream"
    return FileResponse(path, media_type=media, filename=filename)


# ---------------- 前端静态资源 ----------------

app.mount("/", StaticFiles(directory=FRONTEND, html=True), name="frontend")


def main() -> None:
    import os

    import uvicorn

    _ensure_worker()
    url = "http://127.0.0.1:7861"
    # 调试器(VS Code F5)下由 serverReadyAction 开浏览器,避免开两个标签页
    if not os.environ.get("PROTEINFOLDING_NO_BROWSER"):
        threading.Timer(1.2, lambda: webbrowser.open(url)).start()
    logger.info("ProteinFolding 服务启动: %s", url)
    uvicorn.run(app, host="127.0.0.1", port=7861, log_level="warning")


if __name__ == "__main__":
    main()
