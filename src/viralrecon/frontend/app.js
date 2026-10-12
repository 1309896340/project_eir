/* ViralRecon 演示台前端逻辑:参数表单 -> POST /api/run -> MultiQC 风格报告渲染 */
"use strict";

const $ = (sel) => document.querySelector(sel);

const state = {
  preset: "standard",
  presets: [],
  ranges: {},
  report: null,
};

/* ---------------- 初始化 ---------------- */

async function init() {
  const [presetData, steps] = await Promise.all([
    fetch("/api/presets").then((r) => r.json()),
    fetch("/api/pipeline").then((r) => r.json()),
  ]);
  state.presets = presetData.presets;
  state.ranges = presetData.param_ranges;
  renderPresetChips();
  renderAdvancedForm();
  renderPipelineStrip(steps);
  $("#run-btn").addEventListener("click", runAnalysis);
  $("#variant-filter").addEventListener("change", () => renderVariantTable());
  $("#download-tsv").addEventListener("click", downloadTsv);
  $("#download-fasta").addEventListener("click", downloadFasta);
}

function renderPresetChips() {
  const row = $("#preset-row");
  row.innerHTML = "";
  for (const p of state.presets) {
    const chip = document.createElement("div");
    chip.className = "preset-chip" + (p.key === state.preset ? " active" : "");
    chip.innerHTML = `<div class="name">${p.label}</div>
      <div class="desc">${describeParams(p.params)}</div>`;
    chip.addEventListener("click", () => {
      state.preset = p.key;
      document.querySelectorAll(".preset-chip").forEach((c) => c.classList.remove("active"));
      chip.classList.add("active");
    });
    row.appendChild(chip);
  }
}

function describeParams(p) {
  const kb = p.genome_len >= 1000 ? `${(p.genome_len / 1000).toFixed(0)}kb` : `${p.genome_len}bp`;
  return `${kb} 基因组 · ${p.coverage}x · read ${p.read_len}bp · ${p.n_snps} SNP / ${p.n_indels} indel`;
}

/* 高级参数表单:字段名 -> 中文标签/单位 */
const ADV_FIELDS = [
  ["genome_len", "基因组长度 (bp)"],
  ["coverage", "覆盖深度 (x)"],
  ["read_len", "read 长度 (bp)"],
  ["n_snps", "注入 SNP 数"],
  ["n_indels", "注入 Indel 数"],
  ["error_rate_scale", "错误率缩放"],
  ["low_quality_frac", "低质量 reads 比例"],
  ["q_start", "起始质量 Q"],
  ["q_end", "末端质量 Q"],
];

function renderAdvancedForm() {
  const grid = $("#adv-grid");
  grid.innerHTML = "";
  for (const [key, label] of ADV_FIELDS) {
    const range = state.ranges[key];
    const wrap = document.createElement("label");
    wrap.innerHTML = `${label}
      <input type="number" data-param="${key}"
        min="${range ? range[0] : 0}" max="${range ? range[1] : 100}"
        step="${key.includes("frac") || key.includes("scale") ? "0.05" : "1"}"
        placeholder="默认">`;
    grid.appendChild(wrap);
  }
  const modeWrap = document.createElement("label");
  modeWrap.innerHTML = `覆盖模式
    <select data-param="coverage_mode">
      <option value="">默认 (amplicon)</option>
      <option value="amplicon">amplicon 扩增子波浪</option>
      <option value="uniform">uniform 均匀覆盖</option>
    </select>`;
  grid.appendChild(modeWrap);
}

function renderPipelineStrip(steps) {
  const strip = $("#pipeline-strip");
  strip.innerHTML = "";
  steps.forEach((s, i) => {
    if (i > 0) {
      const arrow = document.createElement("div");
      arrow.className = "step-arrow";
      arrow.textContent = "➜";
      strip.appendChild(arrow);
    }
    const el = document.createElement("div");
    el.className = "step";
    el.innerHTML = `<div class="step-label">${i + 1}. ${s.label}</div>
      <div class="step-tool">${s.tool}</div>
      <div class="step-desc">${s.desc}</div>`;
    strip.appendChild(el);
  });
}

/* ---------------- 运行与报告渲染 ---------------- */

async function runAnalysis() {
  const params = {};
  document.querySelectorAll("[data-param]").forEach((el) => {
    if (el.value !== "") params[el.dataset.param] = el.value;
  });
  const seedRaw = $("#seed-input").value;
  const body = {
    preset: state.preset,
    params,
    seed: seedRaw === "" ? null : Number(seedRaw),
  };

  const btn = $("#run-btn");
  btn.disabled = true;
  btn.textContent = "⏳ 计算中…";
  $("#report").hidden = false;
  $("#report").classList.add("loading");
  clearError();
  try {
    const res = await fetch("/api/run", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
    });
    if (!res.ok) {
      const err = await res.json().catch(() => ({}));
      throw new Error(err.detail || `HTTP ${res.status}`);
    }
    state.report = await res.json();
    renderReport(state.report);
  } catch (e) {
    showError(`运行失败:${e.message}`);
  } finally {
    btn.disabled = false;
    btn.textContent = "▶ 生成分析报告";
    $("#report").classList.remove("loading");
  }
}

function renderReport(r) {
  $("#report-meta").textContent =
    `· ${r.params.genome_len}bp · ${r.params.coverage}x · seed=${r.params.seed} · ${r.timings.total_s}s`;

  renderSummaryCards(r);
  renderTruthCards(r.variants);
  drawQualityChart(r.qc);
  drawGcChart(r.qc);
  drawDepthChart(r.coverage, r.params.genome_len);
  renderVariantTable();
  renderConsensusCards(r.consensus);
}

function renderSummaryCards(r) {
  const v = r.variants;
  const passN = v.rows.filter((x) => x.filter === "PASS").length;
  const items = [
    { val: fmtInt(r.qc.n_reads), key: "reads 数" },
    { val: fmtInt(r.qc.total_bases), key: "总碱基" },
    { val: r.coverage.mean.toFixed(0) + "x", key: "平均深度" },
    { val: (r.coverage.breadth * 100).toFixed(1) + "%", key: "覆盖广度(≥1x)" },
    { val: (r.coverage.pct_ge_10x * 100).toFixed(1) + "%", key: "≥10x 比例", cls: r.coverage.pct_ge_10x > 0.95 ? "good" : "warn" },
    { val: String(passN), key: "PASS 变异", cls: "good" },
    { val: `${v.tp}/${v.truth_total}`, key: "检出/注入(真值)" },
    { val: (r.qc.mean_gc * 100).toFixed(1) + "%", key: "平均 GC" },
    { val: r.timings.total_s + "s", key: "计算耗时" },
  ];
  $("#summary-cards").innerHTML = items.map((c) =>
    `<div class="card ${c.cls || ""}"><div class="val">${c.val}</div><div class="key">${c.key}</div></div>`
  ).join("");
}

function renderTruthCards(v) {
  const items = [
    { val: String(v.truth_total), key: "注入真值变异" },
    { val: String(v.tp), key: "TP 正确检出", cls: "good" },
    { val: String(v.fn), key: "FN 漏检", cls: v.fn > 0 ? "warn" : "good" },
    { val: String(v.fp), key: "FP 假阳性", cls: v.fp > 0 ? "bad" : "good" },
  ];
  $("#truth-cards").innerHTML = items.map((c) =>
    `<div class="card ${c.cls || ""}"><div class="val">${c.val}</div><div class="key">${c.key}</div></div>`
  ).join("");
}

function renderConsensusCards(c) {
  const nPct = (c.n_count / Math.max(1, c.length) * 100).toFixed(1);
  const items = [
    { val: fmtInt(c.length) + "bp", key: "共识长度" },
    { val: String(c.diffs_vs_ref), key: "与参考差异数" },
    { val: nPct + "%", key: "N 比例", cls: parseFloat(nPct) < 2 ? "good" : "warn" },
  ];
  $("#consensus-cards").innerHTML = items.map((c2) =>
    `<div class="card ${c2.cls || ""}"><div class="val">${c2.val}</div><div class="key">${c2.key}</div></div>`
  ).join("");
}

/* ---------------- 变异表 ---------------- */

function renderVariantTable() {
  const filter = $("#variant-filter").value;
  const rows = state.report.variants.rows.filter((r) =>
    filter === "all" ? true :
    filter === "PASS" ? r.filter === "PASS" : r.filter !== "PASS"
  );
  const tbody = $("#variant-table tbody");
  tbody.innerHTML = rows.map((r) => `<tr>
      <td>${r.pos}</td>
      <td class="kind-${r.kind}">${r.kind.toUpperCase()}</td>
      <td>${r.ref || "—"}</td>
      <td>${r.alt || "—"}</td>
      <td>${r.dp}</td>
      <td>${(r.af * 100).toFixed(1)}%</td>
      <td>${r.qual.toFixed(1)}</td>
      <td class="filter-${r.filter}">${r.filter}</td>
    </tr>`).join("");
  if (!rows.length) {
    tbody.innerHTML = `<tr><td colspan="8" style="text-align:center;color:#64748b;padding:20px">
      (无${filter === "PASS" ? "通过过滤的" : filter === "filtered" ? "被过滤的" : ""}变异)</td></tr>`;
  }
}

/* ---------------- 下载 ---------------- */

function downloadBlob(filename, text, mime) {
  const blob = new Blob([text], { type: mime });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = filename;
  a.click();
  URL.revokeObjectURL(a.href);
}

function downloadTsv() {
  if (!state.report) return;
  const header = "POS\tREF\tALT\tTYPE\tDP\tAF\tQUAL\tFILTER";
  const lines = state.report.variants.rows.map((r) =>
    [r.pos, r.ref, r.alt, r.kind, r.dp, r.af.toFixed(4), r.qual, r.filter].join("\t")
  );
  downloadBlob("variants.tsv", [header, ...lines].join("\n"), "text/tab-separated-values");
}

function downloadFasta() {
  if (!state.report) return;
  downloadBlob("consensus.fasta", state.report.consensus.fasta, "text/plain");
}

/* ---------------- Canvas 图表 ---------------- */

function setupCanvas(canvas) {
  const dpr = window.devicePixelRatio || 1;
  const w = canvas.clientWidth || canvas.parentElement.clientWidth;
  const height = parseInt(canvas.getAttribute("height") || "240", 10);
  canvas.width = w * dpr;
  canvas.height = height * dpr;
  canvas.style.height = height + "px";
  const ctx = canvas.getContext("2d");
  ctx.scale(dpr, dpr);
  return { ctx, w, h: height };
}

function drawQualityChart(qc) {
  const canvas = $("#quality-chart");
  const { ctx, w, h } = setupCanvas(canvas);
  const pad = { l: 42, r: 12, t: 14, b: 30 };
  const n = qc.cycle_q50.length;
  const x = (i) => pad.l + (i / Math.max(1, n - 1)) * (w - pad.l - pad.r);
  const y = (q) => pad.t + (1 - q / 42) * (h - pad.t - pad.b);

  ctx.clearRect(0, 0, w, h);
  // FastQC 风格背景色带(Q 0-10 红 / 10-28 橙黄 / 28-41 绿)
  const bands = [[0, 10, "#f6d5d5"], [10, 28, "#fdf3d0"], [28, 42, "#d9efd9"]];
  for (const [q0, q1, color] of bands) {
    ctx.fillStyle = color;
    ctx.fillRect(pad.l, y(q1), w - pad.l - pad.r, y(q0) - y(q1));
  }
  // 轴
  ctx.strokeStyle = "#b9c6d4";
  ctx.fillStyle = "#64748b";
  ctx.font = "11px sans-serif";
  for (const q of [0, 10, 20, 28, 34, 40]) {
    ctx.fillText(String(q), 16, y(q) + 4);
    ctx.beginPath();
    ctx.moveTo(pad.l, y(q));
    ctx.lineTo(w - pad.r, y(q));
    ctx.globalAlpha = 0.35;
    ctx.stroke();
    ctx.globalAlpha = 1;
  }
  ctx.fillText("cycle (read 内位置) →", w / 2 - 60, h - 8);
  ctx.save();
  ctx.translate(11, h / 2 + 24);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText("Phred Q", 0, 0);
  ctx.restore();
  // 四分位带 + 中位线
  ctx.fillStyle = "rgba(15, 111, 184, 0.16)";
  ctx.beginPath();
  qc.cycle_q25.forEach((q, i) => i === 0 ? ctx.moveTo(x(i), y(q)) : ctx.lineTo(x(i), y(q)));
  for (let i = n - 1; i >= 0; i--) ctx.lineTo(x(i), y(qc.cycle_q75[i]));
  ctx.closePath();
  ctx.fill();
  ctx.strokeStyle = "#0f6fb8";
  ctx.lineWidth = 2;
  ctx.beginPath();
  qc.cycle_q50.forEach((q, i) => i === 0 ? ctx.moveTo(x(i), y(q)) : ctx.lineTo(x(i), y(q)));
  ctx.stroke();
  // Q20 参考线
  ctx.strokeStyle = "#b91c1c";
  ctx.setLineDash([5, 4]);
  ctx.beginPath();
  ctx.moveTo(pad.l, y(20));
  ctx.lineTo(w - pad.r, y(20));
  ctx.stroke();
  ctx.setLineDash([]);
  ctx.fillStyle = "#b91c1c";
  ctx.fillText("Q20", w - pad.r - 30, y(20) - 4);
}

function drawGcChart(qc) {
  const canvas = $("#gc-chart");
  const { ctx, w, h } = setupCanvas(canvas);
  const pad = { l: 42, r: 12, t: 14, b: 30 };
  const hist = qc.gc_hist;
  const maxN = Math.max(...hist, 1);
  const x = (i) => pad.l + (i / 100) * (w - pad.l - pad.r);
  const y = (n) => pad.t + (1 - n / maxN) * (h - pad.t - pad.b);

  ctx.clearRect(0, 0, w, h);
  ctx.fillStyle = "rgba(15, 111, 184, 0.55)";
  const bw = (w - pad.l - pad.r) / 101;
  hist.forEach((n, i) => {
    ctx.fillRect(x(i) - bw / 2, y(n), bw, Math.max(0, h - pad.b - y(n)));
  });
  // 轴
  ctx.strokeStyle = "#b9c6d4";
  ctx.fillStyle = "#64748b";
  ctx.font = "11px sans-serif";
  for (const gc of [0, 20, 40, 60, 80, 100]) {
    ctx.fillText(String(gc), x(gc) - 8, h - 12);
  }
  ctx.fillText("GC% (per read) →", w / 2 - 50, h - 1);
  ctx.save();
  ctx.translate(11, h / 2 + 20);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText("reads 数", 0, 0);
  ctx.restore();
}

function drawDepthChart(cov, genomeLen) {
  const canvas = $("#depth-chart");
  const { ctx, w, h } = setupCanvas(canvas);
  const pad = { l: 46, r: 12, t: 14, b: 28 };
  const depth = cov.depth;
  const maxD = Math.max(...depth, 1) * 1.08;
  const x = (i) => pad.l + (i / Math.max(1, depth.length - 1)) * (w - pad.l - pad.r);
  const y = (d) => pad.t + (1 - d / maxD) * (h - pad.t - pad.b);

  ctx.clearRect(0, 0, w, h);
  // 面积
  const grad = ctx.createLinearGradient(0, pad.t, 0, h - pad.b);
  grad.addColorStop(0, "rgba(23, 130, 88, 0.35)");
  grad.addColorStop(1, "rgba(23, 130, 88, 0.05)");
  ctx.beginPath();
  depth.forEach((d, i) => i === 0 ? ctx.moveTo(x(i), y(d)) : ctx.lineTo(x(i), y(d)));
  ctx.lineTo(x(depth.length - 1), h - pad.b);
  ctx.lineTo(x(0), h - pad.b);
  ctx.closePath();
  ctx.fillStyle = grad;
  ctx.fill();
  ctx.strokeStyle = "#177a58";
  ctx.lineWidth = 1.2;
  ctx.beginPath();
  depth.forEach((d, i) => i === 0 ? ctx.moveTo(x(i), y(d)) : ctx.lineTo(x(i), y(d)));
  ctx.stroke();
  // 平均深度参考线
  ctx.strokeStyle = "#0f6fb8";
  ctx.setLineDash([5, 4]);
  ctx.beginPath();
  ctx.moveTo(pad.l, y(cov.mean));
  ctx.lineTo(w - pad.r, y(cov.mean));
  ctx.stroke();
  ctx.setLineDash([]);
  ctx.fillStyle = "#0f6fb8";
  ctx.font = "11px sans-serif";
  ctx.fillText(`平均 ${cov.mean.toFixed(0)}x`, w - pad.r - 80, y(cov.mean) - 4);
  // 轴
  ctx.fillStyle = "#64748b";
  ctx.fillText("基因组位置 (5'→3') →", w / 2 - 70, h - 8);
  ctx.save();
  ctx.translate(12, h / 2 + 26);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText("深度 (x)", 0, 0);
  ctx.restore();
  for (let k = 0; k <= 4; k++) {
    const d = (maxD * k) / 4;
    ctx.fillText(Math.round(d), 20, y(d) + 4);
  }
  $("#depth-hint").textContent =
    `覆盖模式呈 ${depth.length} 个采样点;峰值 ${cov.max}x,中位 ${cov.median.toFixed(0)}x,` +
    `≥10x 占比 ${(cov.pct_ge_10x * 100).toFixed(1)}%(amplicon 模式的波浪形态来自扩增子效率差异)。`;
}

/* ---------------- 工具与错误 ---------------- */

function fmtInt(n) { return Number(n).toLocaleString("en-US"); }

function showError(msg) {
  clearError();
  const div = document.createElement("div");
  div.className = "error-banner";
  div.textContent = msg;
  $("#report").prepend(div);
}
function clearError() {
  document.querySelectorAll(".error-banner").forEach((e) => e.remove());
}

init().catch((e) => showError(`初始化失败:${e.message}`));
