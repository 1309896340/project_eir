/* ProteinFolding 前端逻辑:Mol* 查看器 + 任务轮询 + 播放器 + webm 导出。 */
"use strict";

const $ = (id) => document.getElementById(id);
const state = {
  viewer: null,
  plugin: null,
  frameCell: null,      // 带有 modelIndex 参数的状态节点
  frameCount: 0,
  currentFrame: 0,
  playing: false,
  lastTick: 0,
  acc: 0,
  key: null,
  meta: null,
  exporting: false,     // 静默导出进行中
  exportAbort: false,   // 用户请求取消
  frameCenters: [],     // 每帧几何中心缓存(轨迹坐标)
  gridReprRefs: null,   // 网格结构的 representation refs(样式切换时跳过)
  gridStructureRef: null, // 网格结构在 hierarchy 中的 cell ref(中心计算时跳过)
};

// ---------------- 初始化 ----------------

async function initViewer() {
  state.viewer = await molstar.Viewer.create("viewport", {
    layoutIsExpanded: false,
    layoutShowControls: false,
    layoutShowLog: false,
    layoutShowSequence: true,
    layoutShowLeftPanel: false,
    viewportShowAnimate: false,
    viewportShowControls: false,
    viewportShowSelectionMode: false,
    pluginSpec: { animations: [] },
  });
  state.plugin = state.viewer.plugin;
}

async function loadSystemInfo() {
  try {
    const sys = await (await fetch("/api/system")).json();
    const pred = sys.predictors.map((p) =>
      `${p.name}${p.available ? "" : "(未就绪)"}`
    ).join(" / ");
    $("sysinfo").innerHTML =
      `平台 <b>${sys.platform_active}</b> · 预测器 <b>${pred}</b> · 长度 ${sys.length_range[0]}~${sys.length_range[1]} aa`;
    const sel = $("p-platform");
    for (const name of sys.platforms) {
      const opt = document.createElement("option");
      opt.value = name; opt.textContent = `强制 ${name}`;
      sel.appendChild(opt);
    }
  } catch (e) {
    $("sysinfo").textContent = "系统信息加载失败: " + e.message;
  }
}

async function loadSamplesAndPresets() {
  const samples = await (await fetch("/api/samples")).json();
  for (const s of samples) {
    const btn = document.createElement("button");
    btn.className = "btn-small";
    btn.textContent = s.label;
    btn.onclick = () => { $("sequence").value = s.sequence; updateCount(); };
    $("sample-btns").appendChild(btn);
  }
  const presets = await (await fetch("/api/presets")).json();
  const seg = $("preset-seg");
  Object.entries(presets).forEach(([key, label], idx) => {
    const lab = document.createElement("label");
    lab.textContent = label;
    const input = document.createElement("input");
    input.type = "radio"; input.name = "preset"; input.value = key;
    if (idx === 1) { input.checked = true; lab.classList.add("on"); }
    input.onchange = () => {
      seg.querySelectorAll("label").forEach((l) => l.classList.remove("on"));
      lab.classList.add("on");
    };
    lab.prepend(input);
    seg.appendChild(lab);
  });
}

// ---------------- 任务提交与轮询 ----------------

function updateCount() {
  const n = $("sequence").value.replace(/[^A-Za-z]/g, "").length;
  $("seq-count").textContent = n + " aa";
}

function setError(msg) {
  const box = $("error-box");
  if (!msg) { box.hidden = true; return; }
  box.hidden = false;
  box.textContent = msg;
}

async function submitJob() {
  setError("");
  const sequence = $("sequence").value;
  const preset = document.querySelector('input[name="preset"]:checked').value;
  const body = { sequence, preset, platform: $("p-platform").value };

  const adv = {};
  const map = { "p-guide": "guide_steps", "p-relax": "relax_steps", "p-frames": "n_frames", "p-kmax": "k_max", "p-temp": "temperature_K" };
  for (const [elId, name] of Object.entries(map)) {
    const raw = $(elId).value.trim();
    if (raw !== "") adv[name] = Number(raw);
  }
  if (Object.keys(adv).length) body.params = adv;
  if ($("p-force").checked) body.force = true;

  $("run-btn").disabled = true;
  $("progress-wrap").hidden = false;
  setProgress(0, "提交任务…");
  $("result-card").hidden = true;
  $("player").hidden = true;

  try {
    const res = await fetch("/api/jobs", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
    });
    if (!res.ok) {
      const err = await res.json();
      throw new Error(err.detail || res.statusText);
    }
    const { job_id } = await res.json();
    pollJob(job_id);
  } catch (e) {
    setError(String(e.message || e));
    $("run-btn").disabled = false;
    $("progress-wrap").hidden = true;
  }
}

function setProgress(frac, msg) {
  $("progress-bar").style.width = `${Math.round(frac * 100)}%`;
  $("progress-msg").textContent = msg;
}

async function pollJob(jobId) {
  for (;;) {
    let job;
    try {
      job = await (await fetch(`/api/jobs/${jobId}`)).json();
    } catch (e) {
      setError("轮询失败: " + e.message);
      $("run-btn").disabled = false;
      return;
    }
    setProgress(job.frac || 0, `${stageLabel(job.stage)} · ${job.message || ""}`);
    if (job.status === "done") {
      $("run-btn").disabled = false;
      $("progress-wrap").hidden = true;
      await showResult(job.result);
      return;
    }
    if (job.status === "error") {
      setError("任务失败: " + (job.message || "未知错误"));
      $("run-btn").disabled = false;
      $("progress-wrap").hidden = true;
      return;
    }
    await new Promise((r) => setTimeout(r, 600));
  }
}

function stageLabel(stage) {
  return ({
    queued: "排队中", predict: "结构预测", linear: "线性链构建",
    minimize: "能量最小化", guide: "引导折叠", relax: "终点松弛", done: "完成",
  })[stage] || stage;
}

// ---------------- 结果加载与播放器 ----------------

async function showResult(result) {
  state.key = result.key;
  state.meta = result.meta;
  $("result-card").hidden = false;
  $("hint").style.display = "none";

  for (const [id, file] of [
    ["dl-native", "native.pdb"], ["dl-linear", "linear.pdb"],
    ["dl-dcd", "trajectory.dcd"], ["dl-meta", "meta.json"],
  ]) {
    const a = $(id); a.href = `/api/results/${result.key}/${file}`;
  }

  if (result.meta) renderMeta(result.meta);

  const url = `/api/results/${result.key}`;
  await loadTrajectoryIntoViewer(`${url}/topology.pdb`, `${url}/trajectory.dcd`);
  $("player").hidden = false;
}

function renderMeta(meta) {
  const p = meta.params || {};
  $("result-meta").innerHTML =
    `平台 <b>${meta.platform}</b> · 原子数 <b>${meta.n_atoms}</b> · 帧数 <b>${(meta.rmsd_A || []).length}</b><br>` +
    `RMSD(相对天然态) <b>${meta.rmsd_start_A}</b> Å → <b>${meta.rmsd_end_A}</b> Å · ` +
    `模拟用时 <b>${meta.simulate_seconds}s</b>(引导 ${meta.guide_seconds}s + 松弛 ${meta.relax_seconds}s)`;
  drawRmsdChart(meta.rmsd_A || []);
}

function drawRmsdChart(rmsd) {
  const canvas = $("rmsd-chart");
  const ctx = canvas.getContext("2d");
  const W = canvas.width, H = canvas.height;
  ctx.clearRect(0, 0, W, H);
  if (!rmsd.length) return;
  const maxV = Math.max(...rmsd), minV = Math.min(...rmsd);
  const pad = 6;
  ctx.strokeStyle = "#4da3ff";
  ctx.lineWidth = 1.6;
  ctx.beginPath();
  rmsd.forEach((v, i) => {
    const x = pad + (i / (rmsd.length - 1)) * (W - 2 * pad);
    const y = H - pad - ((v - minV) / Math.max(maxV - minV, 1e-6)) * (H - 2 * pad);
    i === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
  });
  ctx.stroke();
  ctx.fillStyle = "#7d8ca1";
  ctx.font = "10px sans-serif";
  ctx.fillText(`RMSD ${minV.toFixed(1)} → ${maxV.toFixed(1)} Å`, pad, 11);
}

async function loadTrajectoryIntoViewer(topologyUrl, dcdUrl) {
  const plugin = state.plugin;
  state.frameCell = null;
  // 清空既有状态(molstar 5.x 的 PluginContext.clear)
  await plugin.clear();
  const res = await state.viewer.loadTrajectory({
    model: { kind: "model-url", url: topologyUrl, format: "pdb" },
    coordinates: { kind: "coordinates-url", url: dcdUrl, format: "dcd", isBinary: true },
    preset: "default",
  });
  state.viewer.handleResize();

  // 帧切换:定位坐标轨迹派生的 Model 节点(modelIndex 参数),帧数取自轨迹
  state.frameCell = null;
  state.frameCount = 0;
  plugin.state.data.cells.forEach((cell) => {
    if (!cell.obj || !cell.params || !cell.params.values) return;
    const type = cell.obj.type.name;
    if (type === "Trajectory" && "modelRef" in cell.params.values) {
      const d = cell.obj.data;
      state.frameCount = Math.max(state.frameCount, d.frameCount ?? d.length ?? 0);
    }
    if (type === "Model" && typeof cell.params.values.modelIndex === "number") {
      state.frameCell = cell;  // 同名参数的 Model 节点取最后创建的(坐标轨迹派生)
    }
  });
  if (!state.frameCell || !state.frameCount) {
    console.warn("未找到 modelIndex 参数节点,播放控制不可用");
    showPlayerError("Mol* 轨迹节点未找到,播放控制不可用");
    return;
  }
  state.currentFrame = 0;
  $("frame-slider").max = state.frameCount - 1;
  $("frame-slider").value = 0;
  updateFrameLabel();
  state.playing = false;
  $("play-btn").textContent = "▶";
  // 原点网格参照(每次加载轨迹后重建——开头 clear() 会连网格一起清掉)
  state.gridStructureRef = null;
  state.gridReprRefs = null;
  try {
    await loadOriginGrid();
  } catch (e) {
    console.warn("网格参照加载失败:", e);
  }
  // 相机适配结构并自适应画布尺寸,pivot 落在第 0 帧几何中心
  setTimeout(() => {
    try {
      state.viewer.handleResize();
      if (state.plugin.managers.camera) state.plugin.managers.camera.reset();
    } catch (e) { /* 相机复位失败不影响功能 */ }
    recenterCameraOnFrame();
  }, 300);
}

async function setFrame(i) {
  if (!state.frameCell) return;
  i = Math.max(0, Math.min(i, state.frameCount - 1));
  state.currentFrame = i;
  const cell = state.frameCell;
  const params = { ...cell.params.values, modelIndex: i };
  const d = state.plugin.state.data;
  const b = d.build().to(cell).update(params);
  await state.plugin.runTask(d.updateTree(b));
  updateFrameLabel();
  $("frame-slider").value = i;
  recenterCameraOnFrame();  // 相机 pivot 跟随当前帧几何中心
}

function updateFrameLabel() {
  $("frame-label").textContent = `${state.currentFrame + 1} / ${state.frameCount}`;
}

function showPlayerError(msg) {
  $("frame-label").textContent = "播放不可用";
  console.error(msg);
}

// 播放循环:rAF + 时间累积(setFrame 为异步,防重入避免状态堆积)
let advancing = false;
function tick(t) {
  requestAnimationFrame(tick);
  if (!state.playing || !state.frameCell || advancing) {
    state.lastTick = t;
    return;
  }
  const fps = Number($("speed-sel").value);
  state.acc += (t - state.lastTick) / 1000;
  state.lastTick = t;
  const frameDur = 1 / fps;
  if (state.acc < frameDur) return;
  state.acc = state.acc % frameDur;
  const next = state.currentFrame + 1 >= state.frameCount ? 0 : state.currentFrame + 1;
  advancing = true;
  setFrame(next).finally(() => { advancing = false; });
}
requestAnimationFrame(tick);

$("play-btn").onclick = () => {
  if (state.exporting) return;  // 导出期间禁止播放
  state.playing = !state.playing;
  $("play-btn").textContent = state.playing ? "⏸" : "▶";
};
$("frame-slider").oninput = (e) => {
  if (state.exporting) return;  // 导出期间禁止拖动
  state.playing = false;
  $("play-btn").textContent = "▶";
  setFrame(Number(e.target.value));
};

// 显示样式切换(网格参照的 point 表示保持不变)
$("style-sel").onchange = async (e) => {
  const kind = e.target.value;
  const plugin = state.plugin;
  const reprCells = [];
  plugin.state.data.cells.forEach((cell) => {
    if (cell.obj && cell.obj.type && cell.obj.type.name === "representation") reprCells.push(cell);
  });
  for (const cell of reprCells) {
    if (state.gridReprRefs && state.gridReprRefs.has(cell.transform.ref)) continue;
    const params = JSON.parse(JSON.stringify(cell.transform.params));
    if (params && params.type) {
      params.type.name = kind;
      params.type.params = params.type.params || {};
      const b = plugin.state.data.build().to(cell).update(params);
      await plugin.runTask(plugin.state.data.updateTree(b));
    }
  }
};

// ---------------- 逐帧几何中心与相机跟随 ----------------

function currentStructureCenter() {
  // 取蛋白结构(排除网格参照物)的包围球中心,即当前帧的几何中心。
  // 判据:原子数最多的结构即蛋白(网格仅数百个参考点);gridStructureRef
  // 的 hierarchy 识别存在异步时序问题,只作日志参考。
  const structs = state.plugin.managers.structure.hierarchy.current.structures;
  let target = null, best = -1;
  for (const s of structs) {
    const d = s.cell.obj && s.cell.obj.data;
    const n = d && typeof d.elementCount === "number" ? d.elementCount : 0;
    if (n > best) { best = n; target = s; }
  }
  if (!target) return null;
  const data = target.cell.obj && target.cell.obj.data;
  if (!data || !data.boundary) return null;
  const c = data.boundary.sphere.center;
  return [c[0], c[1], c[2]];
}

function recenterCameraOnFrame() {
  // 把相机 target 平移到当前帧几何中心(position 同步平移,视角/缩放不变):
  // 之后鼠标旋转始终围绕当前帧中心,折叠漂移不会被视锥体裁剪
  const center = currentStructureCenter();
  if (!center) return;
  state.frameCenters[state.currentFrame] = center;
  const c3d = state.plugin.canvas3d;
  if (!c3d) return;
  const cam = c3d.camera;
  const snap = cam.getSnapshot();
  const [tx, ty, tz] = snap.target;
  const dx = center[0] - tx, dy = center[1] - ty, dz = center[2] - tz;
  if (Math.abs(dx) + Math.abs(dy) + Math.abs(dz) < 1e-4) return;
  const [px, py, pz] = snap.position;
  snap.target = center.slice();
  snap.position = [px + dx, py + dy, pz + dz];
  cam.setState(snap);
  c3d.requestDraw();
}

// ---------------- 原点网格参照 ----------------

function buildGridPdb(spacing = 10, half = 100) {
  // xz 平面(y=0)点阵网格 + 原点标记(N 原子,默认元素着色为蓝)
  const lines = [];
  let serial = 1;
  const het = (elem, resname, x, y, z) =>
    `HETATM${String(serial++).padStart(5)}  ${elem.padEnd(2)}  ${resname} A   1    ` +
    `${x.toFixed(3).padStart(8)}${y.toFixed(3).padStart(8)}${z.toFixed(3).padStart(8)}` +
    `  1.00  0.00          ${elem.padStart(2)}`;
  for (let x = -half; x <= half + 1e-6; x += spacing) {
    for (let z = -half; z <= half + 1e-6; z += spacing) {
      lines.push(het("C", "GRD", x, 0, z));
    }
  }
  lines.push(het("N", "ORG", 0, 0, 0));  // 原点标记
  lines.push("END");
  return lines.join("\n");
}

async function loadOriginGrid() {
  // 网格作为独立结构加载,用 point 表示渲染;以加载前后差集识别网格子树
  // (Mol* 的 cell ref 为随机 ID,无法按前缀归属),记录 refs 供
  // 样式切换与中心计算排除
  const plugin = state.plugin;
  const structsBefore = new Set(
    plugin.managers.structure.hierarchy.current.structures.map((s) => s.cell.ref)
  );
  const cellsBefore = new Set();
  plugin.state.data.cells.forEach((_, key) => cellsBefore.add(String(key)));

  await state.viewer.loadStructureFromData(buildGridPdb(), "pdb");

  // hierarchy 行为是异步刷新的,轮询等待网格结构出现(最多 ~2s)
  let gridStruct = null;
  for (let i = 0; i < 20 && !gridStruct; i++) {
    const cur = plugin.managers.structure.hierarchy.current.structures;
    gridStruct = cur.find((s) => !structsBefore.has(s.cell.ref));
    if (!gridStruct) await new Promise((r) => setTimeout(r, 100));
  }
  state.gridStructureRef = gridStruct ? gridStruct.cell.ref : null;

  state.gridReprRefs = new Set();
  const newReprs = [];
  plugin.state.data.cells.forEach((cell, key) => {
    if (cellsBefore.has(String(key))) return;
    const t = cell.obj && cell.obj.type && cell.obj.type.name;
    if (t === "Structure 3D") newReprs.push(cell);
  });
  for (const cell of newReprs) {
    state.gridReprRefs.add(cell.transform.ref);
    const params = JSON.parse(JSON.stringify(cell.transform.params || {}));
    params.type = { name: "point", params: {} };
    const b = plugin.state.data.build().to(cell).update(params);
    await plugin.runTask(plugin.state.data.updateTree(b));
  }
}


// 重置视图:复位相机并重新居中到当前帧
$("reset-btn").onclick = () => {
  try {
    state.viewer.handleResize();
    if (state.plugin.managers.camera) state.plugin.managers.camera.reset();
  } catch (e) { /* 相机复位失败不影响后续居中 */ }
  setTimeout(recenterCameraOnFrame, 250);  // 等 reset 过渡完成后落回当前帧中心
};


$("export-btn").onclick = () => {
  if (state.exporting) { state.exportAbort = true; return; }
  exportVideo();
};

function getRenderCanvas() {
  // Mol* 5.11 的 canvas 挂在 canvas3d.webgl.gl.canvas 上;DOM 查询作回退
  const c3d = state.plugin && state.plugin.canvas3d;
  const viaPlugin = c3d && c3d.webgl && c3d.webgl && c3d.webgl.gl && c3d.webgl.gl.canvas;
  return viaPlugin instanceof HTMLCanvasElement
    ? viaPlugin
    : document.querySelector("#viewport canvas");
}

const sleepMs = (ms) => new Promise((r) => setTimeout(r, ms));
// 三重 rAF:等待 Mol* requestDraw 之后的合成完成,确保快照的是新帧
const nextPaint = () =>
  new Promise((r) =>
    requestAnimationFrame(() =>
      requestAnimationFrame(() => requestAnimationFrame(r))
    )
  );

// 静默导出:不播放预览,逐帧 set -> 渲染 -> requestFrame 推入流。
// captureStream(0) 手动帧模式保证"生成多少帧导出多少帧,不多不少";
// 节拍按所选帧率补偿,渲染慢于帧间隔时以渲染耗时为准(帧数仍精确)。
async function exportVideo() {
  const canvas = getRenderCanvas();
  if (!canvas || !canvas.captureStream || typeof MediaRecorder === "undefined") {
    alert("当前浏览器不支持画布录制(MediaRecorder/captureStream)");
    return;
  }
  if (!state.frameCell || !state.frameCount) {
    alert("尚未加载动画结果,无法导出");
    return;
  }
  const fps = Number($("speed-sel").value);
  const mime = MediaRecorder.isTypeSupported("video/webm;codecs=vp9")
    ? "video/webm;codecs=vp9" : "video/webm";
  const stream = canvas.captureStream(0);  // 手动帧模式
  const track = stream.getVideoTracks()[0];
  if (!track || typeof track.requestFrame !== "function") {
    alert("当前浏览器不支持手动推帧(requestFrame),无法精确导出");
    stream.getTracks().forEach((t) => t.stop());
    return;
  }
  const rec = new MediaRecorder(stream, { mimeType: mime, videoBitsPerSecond: 12_000_000 });
  const chunks = [];
  rec.ondataavailable = (e) => { if (e.data.size) chunks.push(e.data); };

  const N = state.frameCount;
  const startFrame = state.currentFrame;
  state.exporting = true;
  state.exportAbort = false;
  state.playing = false;
  $("play-btn").textContent = "▶";
  const btn = $("export-btn");
  const overlay = $("export-overlay"), bar = $("export-bar"), msg = $("export-msg");
  btn.classList.add("recording");
  btn.textContent = "✕ 取消导出";
  overlay.hidden = false;
  bar.style.width = "0%";
  msg.textContent = `准备导出 ${N} 帧 @ ${fps} 帧/秒`;

  const finish = (ok) => {
    state.exporting = false;
    overlay.hidden = true;
    btn.classList.remove("recording");
    btn.textContent = "● 导出视频";
    if (ok) {
      rec.onstop = () => {
        const blob = new Blob(chunks, { type: "video/webm" });
        const url = URL.createObjectURL(blob);
        const a = document.createElement("a");
        a.href = url;
        a.download = `folding_${state.key || "traj"}.webm`;
        a.click();
        setTimeout(() => URL.revokeObjectURL(url), 5000);
      };
      rec.stop();
    } else {
      rec.onstop = null;
      chunks.length = 0;
      try { rec.stop(); } catch (e) { /* 忽略中止异常 */ }
    }
    setFrame(startFrame);
  };

  rec.start();
  const interval = 1000 / fps;
  const t0 = performance.now();
  try {
    for (let i = 0; i < N; i++) {
      if (state.exportAbort) { finish(false); return; }
      await setFrame(i);
      state.plugin.canvas3d.requestDraw();
      await nextPaint();
      track.requestFrame();  // 恰好一帧
      bar.style.width = `${Math.round(((i + 1) / N) * 100)}%`;
      msg.textContent = `导出中:第 ${i + 1} / ${N} 帧`;
      const wait = t0 + (i + 1) * interval - performance.now();
      if (wait > 0) await sleepMs(wait);
    }
    await sleepMs(Math.max(interval, 150));  // 收尾,确保末帧完整入流
    finish(true);
  } catch (e) {
    console.error("导出失败:", e);
    finish(false);
  }
}

// ---------------- 启动 ----------------

$("sequence").addEventListener("input", updateCount);
$("fasta-file").addEventListener("change", (e) => {
  const file = e.target.files[0];
  if (!file) return;
  file.text().then((t) => { $("sequence").value = t; updateCount(); });
});
$("run-btn").onclick = submitJob;

(async () => {
  updateCount();  // 输入框默认填有 INS 测试样例,初始化字数
  await Promise.all([loadSystemInfo(), loadSamplesAndPresets()]);
  try {
    await initViewer();
  } catch (e) {
    console.error(e);
    $("hint").textContent = "Mol* 初始化失败: " + e.message;
  }
  // 支持 ?key=<cache-key> 直达加载既有结果
  const urlKey = new URLSearchParams(location.search).get("key");
  if (urlKey && /^[0-9a-f]{16}$/.test(urlKey)) {
    try {
      const meta = await (await fetch(`/api/results/${urlKey}/meta.json`)).json();
      await showResult({ key: urlKey, meta });
    } catch (e) {
      setError(`加载结果 ${urlKey} 失败: ` + e.message);
    }
  }
})();
