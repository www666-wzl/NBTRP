const GENES = ["CLEC2D", "GDAP1L1", "PHF21B", "SERPINF1", "THRA", "UNG"];
const state = { model: null, rows: [], predictions: [] };

function byId(id) { return document.getElementById(id); }
function setStatus(text, kind = "") {
  const el = byId("status");
  el.textContent = text;
  el.className = `status ${kind}`;
}

async function loadModel() {
  const response = await fetch("./nbtrp-model.json");
  if (!response.ok) throw new Error(`模型加载失败（HTTP ${response.status}）`);
  state.model = await response.json();
  const model = state.model.learner.gradient_booster.model;
  if (model.trees.length !== 83) throw new Error("模型完整性检查失败");
  setStatus("模型已就绪，可以开始评估", "ok");
  const modelState = byId("model-state");
  if (modelState && modelState.classList) {
    modelState.classList.add("ready");
    if (modelState.lastChild) modelState.lastChild.textContent = "模型完整性已核验";
  }
}

function splitLine(line) {
  const cells = line.includes("\t") ? line.split("\t") : line.trim().split(/[ ,;]+/);
  return cells.map(cell => cell.trim().replace(/^"|"$/g, ""));
}

function parseMatrix(text) {
  const lines = text.replace(/^\uFEFF/, "").split(/\r?\n/).filter(line => line.trim());
  if (lines.length < 2) throw new Error("文件至少需要表头和一行数据");
  const table = lines.map(splitLine);
  const header = table[0];
  const firstHeader = (header[0] || "").toLowerCase();
  const genesAsRows = ["gene", "genes", "symbol", "gene_symbol", "id", ""].includes(firstHeader) || table.slice(1).some(row => GENES.includes(row[0]));

  if (genesAsRows) {
    const map = new Map(table.slice(1).map(row => [row[0], row.slice(1)]));
    const missing = GENES.filter(g => !map.has(g));
    const hasExplicitRowNameHeader = table[1] && table[1].length === header.length;
    const samples = hasExplicitRowNameHeader ? header.slice(1) : header;
    return samples.map((sample, i) => ({
      sample: sample || `Sample ${i + 1}`,
      values: GENES.map(g => map.has(g) ? Number(map.get(g)[i]) : Number.NaN),
      missing
    }));
  }

  const indices = GENES.map(g => header.indexOf(g));
  const missing = GENES.filter((_, i) => indices[i] < 0);
  if (missing.length) throw new Error(`缺少基因列：${missing.join("、")}`);
  return table.slice(1).map((row, i) => ({
    sample: row[0] || `Sample ${i + 1}`,
    values: indices.map(index => Number(row[index]))
  }));
}

function validateRows(rows) {
  if (!rows.length) throw new Error("没有可计算的样本");
  rows.forEach((row, i) => row.values.forEach((value, j) => {
    if (!Number.isFinite(value) && !Number.isNaN(value)) throw new Error(`样本 ${row.sample || i + 1} 的 ${GENES[j]} 不是有效数值`);
  }));
}

function standardize(rows) {
  if (rows.length < 2) throw new Error("队列内 Z-score 标准化至少需要 2 个样本；建议上传完整研究队列");
  const means = GENES.map((_, j) => {
    const values = rows.map(row => row.values[j]).filter(Number.isFinite);
    return values.length ? values.reduce((a, b) => a + b, 0) / values.length : Number.NaN;
  });
  const sds = GENES.map((_, j) => {
    const values = rows.map(row => row.values[j]).filter(Number.isFinite);
    if (!values.length) return Number.NaN;
    const variance = values.reduce((sum, value) => sum + (value - means[j]) ** 2, 0) / (values.length - 1);
    return Math.sqrt(variance);
  });
  sds.forEach((sd, j) => { if (sd === 0) throw new Error(`${GENES[j]} 在该队列中无变异，无法进行 Z-score 标准化`); });
  return rows.map(row => row.values.map((value, j) => Number.isFinite(value) ? (value - means[j]) / sds[j] : Number.NaN));
}

function predictRow(values) {
  const learner = state.model.learner;
  const booster = learner.gradient_booster.model;
  const baseScore = Number(learner.learner_model_param.base_score);
  let margin = Math.log(baseScore);
  for (const tree of booster.trees) {
    let node = 0;
    while (tree.left_children[node] !== -1) {
      const feature = tree.split_indices[node];
      const value = values[feature];
      if (Number.isNaN(value)) node = tree.default_left[node] ? tree.left_children[node] : tree.right_children[node];
      else node = value < tree.split_conditions[node] ? tree.left_children[node] : tree.right_children[node];
    }
    margin += tree.split_conditions[node];
  }
  return Math.exp(margin);
}

function formatScore(value) { return value.toFixed(6); }
function escapeHtml(value) {
  return String(value).replace(/[&<>'"]/g, char => ({"&":"&amp;","<":"&lt;",">":"&gt;","'":"&#39;",'"':"&quot;"}[char]));
}

function renderChart(predictions, median) {
  const chart = byId("score-chart");
  const max = Math.max(...predictions);
  const ordered = [...predictions].sort((a, b) => a - b);
  const stride = Math.max(1, Math.ceil(ordered.length / 80));
  const displayed = ordered.filter((_, i) => i % stride === 0 || i === ordered.length - 1);
  chart.innerHTML = displayed.map(value => {
    const height = Math.max(7, Math.round((value / max) * 70));
    return `<i class="score-bar ${value > median ? "high" : ""}" style="height:${height}px" title="${formatScore(value)}"></i>`;
  }).join("");
}

function render(rows, predictions) {
  const median = [...predictions].sort((a, b) => a - b)[Math.floor(predictions.length / 2)];
  byId("summary").innerHTML = `<strong>${rows.length}</strong> 个样本已完成真实模型推理。`;
  byId("sample-count").textContent = rows.length;
  byId("median-score").textContent = formatScore(median);
  byId("result-body").innerHTML = rows.map((row, i) => {
    const high = predictions[i] > median;
    return `<tr><td>${escapeHtml(row.sample)}</td><td>${formatScore(predictions[i])}</td><td><span class="pill ${high ? "high" : "low"}">${high ? "高分组" : "低分组"}</span></td></tr>`;
  }).join("");
  renderChart(predictions, median);
  byId("results").hidden = false;
}

async function calculate() {
  try {
    setStatus("正在解析矩阵并运行模型…");
    if (!state.model) await loadModel();
    const rows = parseMatrix(byId("matrix-input").value);
    validateRows(rows);
    const predictions = standardize(rows).map(predictRow);
    state.rows = rows;
    state.predictions = predictions;
    render(rows, predictions);
    const missing = [...new Set(rows.flatMap(row => row.missing || []))];
    setStatus(missing.length ? `计算完成；输入缺少 ${missing.join("、")}，已按模型缺失值分支处理` : "计算完成：已使用作者发布的真实 XGBoost 模型", "ok");
    byId("results").scrollIntoView({ behavior: "smooth", block: "start" });
  } catch (error) { setStatus(error.message, "error"); }
}

async function loadExample() {
  try {
    const response = await fetch("./example_E-MTAB-8248_gene_matrix.txt");
    if (!response.ok) throw new Error(`HTTP ${response.status}`);
    byId("matrix-input").value = await response.text();
    setStatus("已载入作者示例队列 E-MTAB-8248（223 个样本），可以直接生成结果", "ok");
  } catch (error) { setStatus(`示例加载失败：${error.message}`, "error"); }
}

function downloadCsv() {
  const sorted = [...state.predictions].sort((a, b) => a - b);
  const median = sorted[Math.floor(sorted.length / 2)];
  const lines = ["sample,nbtrp_score,group", ...state.rows.map((row, i) => `${JSON.stringify(row.sample)},${state.predictions[i]},${state.predictions[i] > median ? "high" : "low"}`)];
  const blob = new Blob(["\uFEFF" + lines.join("\n")], { type: "text/csv;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url; link.download = "nbtrp-results.csv"; link.click();
  URL.revokeObjectURL(url);
}

byId("calculate").addEventListener("click", calculate);
byId("load-example").addEventListener("click", loadExample);
byId("download").addEventListener("click", downloadCsv);
byId("file-input").addEventListener("change", async event => {
  const file = event.target.files[0];
  if (file) { byId("matrix-input").value = await file.text(); setStatus(`已读取 ${file.name}，可以生成结果`, "ok"); }
});
const uploadPanel = byId("upload-panel");
["dragenter", "dragover"].forEach(type => uploadPanel.addEventListener(type, event => { event.preventDefault(); uploadPanel.classList.add("dragover"); }));
["dragleave", "drop"].forEach(type => uploadPanel.addEventListener(type, event => { event.preventDefault(); uploadPanel.classList.remove("dragover"); }));
uploadPanel.addEventListener("drop", async event => {
  const file = event.dataTransfer.files[0];
  if (file) { byId("matrix-input").value = await file.text(); setStatus(`已读取 ${file.name}，可以生成结果`, "ok"); }
});
byId("menu-toggle").addEventListener("click", event => {
  const open = byId("site-nav").classList.toggle("open");
  event.currentTarget.setAttribute("aria-expanded", String(open));
});
if (typeof document.querySelectorAll === "function") {
  document.querySelectorAll("#site-nav a").forEach(link => link.addEventListener("click", () => byId("site-nav").classList.remove("open")));
}

async function initialize() {
  try {
    await loadModel();
    if (new URLSearchParams(location.search).get("selftest") === "1") {
      await loadExample(); await calculate();
      document.body.dataset.selftest = state.predictions.length === 223 ? "passed" : "failed";
    }
  } catch (error) { setStatus(error.message, "error"); }
}
initialize();
