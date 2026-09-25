/* ProteinFolding 前端逻辑:Mol* 查看器 + 任务轮询 + 播放器 + webm 导出。 */
"use strict";

const $ = (id) => document.getElementById(id);
const state = {
  viewer: null,
  plugin: null,
  frameCell: null,      // 带有 frameIndex 参数的状态节点
  frameCount: 0,
  currentFrame: 0,
  playing: false,
  lastTick: 0,
  acc: 0,
  key: null,
  meta: null,
  recording: null,      // MediaRecorder 会话
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
  // 相机适配结构并自适应画布尺寸
  setTimeout(() => {
    try {
      state.viewer.handleResize();
      if (state.plugin.managers.camera) state.plugin.managers.camera.reset();
    } catch (e) { /* 相机复位失败不影响功能 */ }
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
  if (!state.playing) $("frame-slider").value = i;
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
  state.playing = !state.playing;
  $("play-btn").textContent = state.playing ? "⏸" : "▶";
};
$("frame-slider").oninput = (e) => {
  state.playing = false;
  $("play-btn").textContent = "▶";
  setFrame(Number(e.target.value));
};

// 显示样式切换
$("style-sel").onchange = async (e) => {
  const kind = e.target.value;
  const plugin = state.plugin;
  const reprCells = [];
  plugin.state.data.cells.forEach((cell) => {
    if (cell.obj && cell.obj.type && cell.obj.type.name === "representation") reprCells.push(cell);
  });
  for (const cell of reprCells) {
    const params = JSON.parse(JSON.stringify(cell.transform.params));
    if (params && params.type) {
      params.type.name = kind;
      params.type.params = params.type.params || {};
      await plugin.state.data.update(cell.transform.ref, params);
    }
  }
};

// ---------------- webm 视频导出(手动触发) ----------------

$("export-btn").onclick = () => {
  if (state.recording) { stopRecording(); return; }
  startRecording();
};

function startRecording() {
  const canvas = state.plugin.canvas3d && state.plugin.canvas3d.canvas;
  if (!canvas || !canvas.captureStream) {
    alert("当前浏览器不支持画布录制(MediaRecorder/captureStream)");
    return;
  }
  const stream = canvas.captureStream(Number($("speed-sel").value));
  const mime = MediaRecorder.isTypeSupported("video/webm;codecs=vp9")
    ? "video/webm;codecs=vp9" : "video/webm";
  const recorder = new MediaRecorder(stream, { mimeType: mime, videoBitsPerSecond: 12_000_000 });
  const chunks = [];
  recorder.ondataavailable = (e) => { if (e.data.size) chunks.push(e.data); };
  recorder.onstop = () => {
    const blob = new Blob(chunks, { type: "video/webm" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `folding_${state.key || "traj"}.webm`;
    a.click();
    setTimeout(() => URL.revokeObjectURL(url), 5000);
  };
  state.recording = recorder;
  recorder.start();
  $("export-btn").classList.add("recording");
  $("export-btn").textContent = "■ 停止录制";
  // 从头播放一遍
  state.playing = false;
  setFrame(0).then(() => { state.playing = true; $("play-btn").textContent = "⏸"; });
}

function stopRecording() {
  if (!state.recording) return;
  state.recording.stop();
  state.recording = null;
  $("export-btn").classList.remove("recording");
  $("export-btn").textContent = "● 导出视频";
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
