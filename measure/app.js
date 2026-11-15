/* Histology Measure, draft 0.1
 * Static, no build step, no server. Images are read from a local folder and never leave the machine.
 * Geometry is stored in IMAGE PIXELS; the calibration is stored separately (see measure/README.md).
 */
(function () {
"use strict";
const TOOL_VERSION = "0.1.0-draft";
const SCHEMA_VERSION = "0.1";
const IMG_EXT = /\.(png|jpe?g|webp|gif|bmp)$/i;
const KINDS = {
  periosteum: { label: "Periosteum", color: "var(--c-periosteum)" },
  medullary:  { label: "Medullary margin", color: "var(--c-medullary)" },
  mark:       { label: "Growth mark", color: "var(--c-mark)" },
  efs_inner:  { label: "EFS inner boundary", color: "var(--c-efs_inner)" },
};
const MARK_TYPES = ["LAG", "annulus", "double_multiple", "uncertain"];
const SCORES = [
  ["tissue_dominant", "Dominant tissue", ["", "woven", "fibrolamellar", "parallel_fibred", "lamellar", "mixed", "unknown"]],
  ["vascular_organisation", "Vascular organisation", ["", "avascular", "longitudinal", "circumferential_laminar", "plexiform", "reticular", "radial", "mixed", "unknown"]],
  ["growth_mark_type", "Growth marks", ["", "LAG", "annulus", "both", "none_visible", "unknown"]],
  ["efs", "EFS / OCL", ["", "present", "absent", "unknown"]],
  ["efs_traceable", "EFS traceable", ["", "yes", "no_present_not_traceable", "not_applicable"]],
  ["secondary_remodelling", "Secondary remodelling", ["", "none", "scattered", "extensive", "unknown"]],
  ["medullary_bone", "Medullary bone", ["", "present", "absent", "unknown"]],
  ["kastschenko_line", "Kastschenko's line", ["", "present", "absent", "unknown"]],
  ["image_quality", "Image quality", ["", "good", "fair", "poor"]],
];
const IDENTITY = ["source_id", "citation", "taxon", "specimen", "element", "view", "light", "quality_tier", "scale_value"];

const $ = (id) => document.getElementById(id);
const svg = $("svg"), layerImg = $("layerImg"), layerShapes = $("layerShapes"), layerVerts = $("layerVerts");
const NS = "http://www.w3.org/2000/svg";
const S = { tasks: [], idx: -1, files: new Map(), dirHandle: null, doc: null, tool: "pan", active: null,
            view: { x: 0, y: 0, w: 1000, h: 800 }, img: null, undo: [], scaleClick: null, drag: null,
            space: false, lastActivity: Date.now(), selCurve: null, current: null, filter: { contrast: 100, brightness: 100 } };

// ---------- storage (per-viewer convenience only; the exported file is the record) ----------
const store = {
  get(k) { try { const v = localStorage.getItem(k); return v ? JSON.parse(v) : null; } catch (e) { return null; } },
  set(k, v) { try { localStorage.setItem(k, JSON.stringify(v)); } catch (e) { /* storage unavailable */ } },
};
const keyDoc = (fid) => "hm:doc:" + fid, keyStatus = (fid) => "hm:status:" + fid;

// ---------- helpers ----------
const el = (tag, attrs = {}, parent) => { const e = document.createElementNS(NS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]); if (parent) parent.appendChild(e); return e; };
const h = (tag, props = {}, ...kids) => { const e = document.createElement(tag);
  for (const k in props) { if (k === "class") e.className = props[k]; else if (k.startsWith("on")) e.addEventListener(k.slice(2), props[k]);
    else if (k in e && k !== "list") e[k] = props[k]; else e.setAttribute(k, props[k]); }
  for (const c of kids) e.append(c); return e; };
const dist = (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1]);
const fmt = (v, d = 2) => (v == null || !isFinite(v)) ? "–" : (Math.abs(v) >= 100 ? v.toFixed(0) : v.toFixed(d));
const uid = () => Math.random().toString(36).slice(2, 9);
function status(msg) { $("status").textContent = msg; }
function hint(msg) { $("hint").textContent = msg || ""; }
function touch() { S.lastActivity = Date.now(); }

function pxPerUnit() { const c = S.doc && S.doc.calibration;
  return c && c.p1 && c.p2 && c.length > 0 ? dist(c.p1, c.p2) / c.length : null; }
function unit() { return (S.doc && S.doc.calibration && S.doc.calibration.unit) || "px"; }
function conv(px) { const k = pxPerUnit(); return k ? px / k : px; }

// ---------- folder and manifest ----------
async function openFolder() {
  if (window.showDirectoryPicker) {
    try {
      const dir = await window.showDirectoryPicker({ mode: "readwrite" });
      const files = new Map();
      for await (const [name, entry] of dir.entries()) if (entry.kind === "file") files.set(name, await entry.getFile());
      S.dirHandle = dir; await loadFiles(files); return;
    } catch (e) { if (e.name === "AbortError") return; console.warn("Directory picker failed, using file input", e); }
  }
  $("dirInput").click();
}
$("dirInput").addEventListener("change", async (ev) => {
  const files = new Map();
  for (const f of ev.target.files) { const parts = (f.webkitRelativePath || f.name).split("/");
    if (parts.length <= 2) files.set(f.name, f); }       // top level of the chosen folder only
  S.dirHandle = null; await loadFiles(files); ev.target.value = "";
});

async function loadFiles(files) {
  saveDoc();
  S.files = files; let tasks = [];
  const mf = files.get("manifest.json");
  if (mf) {
    try { const m = JSON.parse(await mf.text()); tasks = (m.tasks || []).map((t) => Object.assign({ batch_id: m.batch_id || null }, t)); }
    catch (e) { status("manifest.json could not be read: " + e.message); }
  }
  if (!tasks.length) {
    tasks = [...files.keys()].filter((n) => IMG_EXT.test(n)).sort().map((n) => ({ figure_id: n.replace(/\.[^.]+$/, ""), file: n, _no_manifest: true }));
    if (!mf) status("No manifest.json: tasks built from image file names; identity fields are editable.");
  }
  for (const t of tasks) t._missing = !files.has(t.file);
  S.tasks = tasks; S.idx = -1; renderTasks();
  const first = tasks.findIndex((t) => !t._missing && store.get(keyStatus(t.figure_id)) !== "done");
  if (tasks.length) openTask(first >= 0 ? first : 0);
}

function renderTasks() {
  const ul = $("taskList"); ul.innerHTML = "";
  const done = S.tasks.filter((t) => store.get(keyStatus(t.figure_id)) === "done").length;
  $("taskCount").textContent = S.tasks.length ? `${done} / ${S.tasks.length} done` : "";
  if (!S.tasks.length) { ul.append(h("li", { class: "muted" }, "Open a folder or a demo image.")); return; }
  S.tasks.forEach((t, i) => {
    const st = store.get(keyStatus(t.figure_id));
    const li = h("li", { class: i === S.idx ? "current" : "", onclick: () => openTask(i), title: t.file },
      h("span", { class: "dot " + (st === "done" ? "done" : st === "prog" ? "prog" : "") }),
      h("span", {}, t.figure_id + (t._missing ? " (file missing)" : "")));
    ul.append(li);
  });
}

async function sha256(file) {
  if (!(window.crypto && crypto.subtle)) return null;
  const buf = await file.arrayBuffer(); const d = await crypto.subtle.digest("SHA-256", buf);
  return [...new Uint8Array(d)].map((b) => b.toString(16).padStart(2, "0")).join("");
}

async function openTask(i) {
  if (i < 0 || i >= S.tasks.length) return;
  saveDoc();
  const t = S.tasks[i];
  if (t._missing) { status(`Image file not in folder: ${t.file}`); return; }
  S.idx = i; S.active = null; S.undo = []; S.scaleClick = null; S.selCurve = null; S.current = null; closeTag();
  const file = S.files.get(t.file);
  if (S.img && S.img.url) URL.revokeObjectURL(S.img.url);
  const url = URL.createObjectURL(file);
  const dims = await new Promise((res, rej) => { const im = new Image(); im.onload = () => res([im.naturalWidth, im.naturalHeight]);
    im.onerror = () => rej(new Error("cannot display this image format")); im.src = url; }).catch((e) => { status(`${t.file}: ${e.message}`); return null; });
  if (!dims) return;
  S.img = { url, w: dims[0], h: dims[1] };
  const saved = store.get(keyDoc(t.figure_id));
  S.doc = saved && saved.schema_version === SCHEMA_VERSION ? saved : newDoc(t);
  S.doc.image.px_width = dims[0]; S.doc.image.px_height = dims[1];
  const hash = await sha256(file);
  S.doc.image.sha256 = hash;
  S.doc.image.sha256_matches_manifest = (t.sha256 && hash) ? t.sha256 === hash : null;
  if (S.doc.image.sha256_matches_manifest === false) status("Warning: image hash differs from the manifest.");
  else status(`${t.file} · ${dims[0]} × ${dims[1]} px`);
  layerImg.innerHTML = "";
  el("image", { href: url, x: 0, y: 0, width: dims[0], height: dims[1], id: "theImage" }, layerImg);
  applyFilter(); renderTasks(); renderMeta(); renderScores(); fit(); renderAll(); setTool("pan");
}

function newDoc(t) {
  const task = {}; for (const k in t) if (!k.startsWith("_")) task[k] = t[k];
  return { schema_version: SCHEMA_VERSION, tool_version: TOOL_VERSION, task, identity_edits: {},
           image: { file: t.file, sha256: null, px_width: null, px_height: null },
           measurer: $("measurer").value.trim() || null, repeat: !!(t.assignment && t.assignment.repeat),
           started_at: new Date().toISOString(), time_spent_s: 0,
           calibration: null, centre: null, flags: {}, curves: [], scores: {}, notes: "" };
}
function saveDoc() {
  if (!S.doc) return;
  S.doc.notes = $("notes").value;
  store.set(keyDoc(S.doc.task.figure_id), S.doc);
  if (store.get(keyStatus(S.doc.task.figure_id)) !== "done" && (S.doc.curves.length || S.doc.calibration))
    store.set(keyStatus(S.doc.task.figure_id), "prog");
}
setInterval(() => { if (S.doc && document.visibilityState === "visible" && Date.now() - S.lastActivity < 120000) S.doc.time_spent_s += 5; }, 5000);
setInterval(saveDoc, 15000);
window.addEventListener("beforeunload", saveDoc);

function renderMeta() {
  const dl = $("meta"); dl.innerHTML = ""; const t = S.doc.task;
  for (const k of IDENTITY) {
    const v = t[k];
    dl.append(h("dt", {}, k.replace(/_/g, " ")));
    if (v != null && v !== "") dl.append(h("dd", {}, String(v)));
    else { const inp = h("input", { value: S.doc.identity_edits[k] || "", placeholder: "not in manifest",
        oninput: (e) => { S.doc.identity_edits[k] = e.target.value; touch(); } });
      dl.append(h("dd", {}, inp)); }
  }
  $("caption").textContent = t.caption || "";
  $("notes").value = S.doc.notes || "";
  document.querySelectorAll("[data-flag]").forEach((cb) => { cb.checked = !!S.doc.flags[cb.dataset.flag]; });
}

// ---------- view ----------
function svgSize() { const r = svg.getBoundingClientRect(); return [Math.max(r.width, 1), Math.max(r.height, 1)]; }
function setView() { const v = S.view; svg.setAttribute("viewBox", `${v.x} ${v.y} ${v.w} ${v.h}`);
  const k = upp(); const im = $("theImage"); if (im) im.style.imageRendering = k < 0.5 ? "pixelated" : "auto"; renderAll(); }
function upp() { const [cw] = svgSize(); return S.view.w / cw; }       // image units per screen pixel
function fit() {
  if (!S.img) return; const [cw, ch] = svgSize(); const pad = 1.04;
  let w = S.img.w * pad, hh = S.img.h * pad;
  if (w / hh < cw / ch) w = hh * cw / ch; else hh = w * ch / cw;
  S.view = { x: S.img.w / 2 - w / 2, y: S.img.h / 2 - hh / 2, w, h: hh }; setView();
}
window.addEventListener("resize", () => { const [cw, ch] = svgSize(); const cx = S.view.x + S.view.w / 2, cy = S.view.y + S.view.h / 2;
  S.view.h = S.view.w * ch / cw; S.view.x = cx - S.view.w / 2; S.view.y = cy - S.view.h / 2; setView(); });
function toImg(ev) { const p = svg.createSVGPoint(); p.x = ev.clientX; p.y = ev.clientY;
  const q = p.matrixTransform(svg.getScreenCTM().inverse()); return [q.x, q.y]; }
function zoomAt(p, f) { const v = S.view; const nw = v.w * f; if (nw < 5 || nw > 200000) return;
  S.view = { x: p[0] - (p[0] - v.x) * f, y: p[1] - (p[1] - v.y) * f, w: nw, h: v.h * f }; setView(); }
function applyFilter() { const im = $("theImage"); if (im) im.style.filter = `contrast(${S.filter.contrast}%) brightness(${S.filter.brightness}%)`; }

// ---------- undo ----------
function pushUndo() { if (!S.doc) return; S.undo.push(JSON.stringify({ c: S.doc.calibration, z: S.doc.centre, v: S.doc.curves }));
  if (S.undo.length > 200) S.undo.shift(); touch(); }
function undo() { const u = S.undo.pop(); if (!u) return; const o = JSON.parse(u);
  S.doc.calibration = o.c; S.doc.centre = o.z; S.doc.curves = o.v;
  if (S.active && !S.doc.curves.find((c) => c.id === S.active.curveId)) S.active = null; renderAll(); renderCurveList(); }

// ---------- tools ----------
function setTool(t) {
  S.tool = t; S.scaleClick = null;
  if (t !== "trace") S.active = null;
  document.querySelectorAll("[data-tool]").forEach((b) => b.classList.toggle("active", b.dataset.tool === t));
  svg.classList.toggle("pan", t === "pan");
  updateHint(); renderAll(); renderCurveList();
}
function updateHint() {
  if (!S.doc) return hint("Open a folder of images, or try the demo image.");
  if (S.tool === "scale") return hint(S.scaleClick ? "Click the OTHER end of the scale bar." : "Click one end of the scale bar.");
  if (S.tool === "centre") return hint("Click the centre of the section.");
  if (S.tool === "trace" && S.active) { const c = curve(S.active.curveId);
    const seg = S.active.seg == null ? c.segments.length + 1 : S.active.seg + 1;
    return hint(`Tracing ${curveName(c)} · segment ${seg}: click along it. Enter = gap (new segment) · T = tag · Esc = stop.`); }
  hint("");
}
document.querySelectorAll("[data-tool]").forEach((b) => b.addEventListener("click", () => setTool(b.dataset.tool)));
document.querySelectorAll("[data-flag]").forEach((cb) => cb.addEventListener("change", () => { if (!S.doc) return;
  S.doc.flags[cb.dataset.flag] = cb.checked; touch(); renderValidation(); }));

const curve = (id) => S.doc.curves.find((c) => c.id === id);
function curveName(c) { return c.kind === "mark" ? `Mark ${c.mark_index == null ? "?" : c.mark_index}${c.index_uncertain ? "?" : ""}` : KINDS[c.kind].label; }
function markDefaults(c) {
  const marks = S.doc.curves.filter((x) => x.kind === "mark" && x.id !== c.id);
  Object.assign(c, { mark_index: marks.reduce((m, x) => Math.max(m, x.mark_index || 0), 0) + 1,
    index_uncertain: false, mark_type: "LAG", confidence: 3, linked_to: null, in_efs: false });
}
function setKind(c, kind) {
  if (c.kind === kind) return; pushUndo(); c.kind = kind;
  if (kind === "mark") markDefaults(c);
  else for (const k of ["mark_index", "index_uncertain", "mark_type", "confidence", "linked_to", "in_efs"]) delete c[k];
}
// Start = a new curve, tagged as the next growth mark by default (retag with Tag)
function newCurve(kind = "mark") {
  if (!S.doc) return; pushUndo();
  const c = { id: uid(), kind, segments: [] };
  if (kind === "mark") { S.doc.curves.push(c); markDefaults(c); } else S.doc.curves.push(c);
  startTrace(c.id);
}
function startTrace(id) { S.tool = "trace"; S.active = { curveId: id, seg: null }; S.current = id;
  document.querySelectorAll("[data-tool]").forEach((b) => b.classList.remove("active")); svg.classList.remove("pan");
  updateHint(); renderAll(); renderCurveList(); svg.focus(); }
function endSegment() { if (S.active) { S.active.seg = null; updateHint(); renderAll(); renderCurveList(); } }
function stopTrace() { endSegment(); setTool("pan"); }
$("btnStart").addEventListener("click", () => newCurve());
$("btnStop").addEventListener("click", stopTrace);
$("btnTag").addEventListener("click", (e) => { e.stopPropagation(); toggleTag($("btnTag")); });

// ---------- tag mini-menu ----------
function tagTarget() { const id = S.active ? S.active.curveId : S.current; return id && S.doc ? curve(id) : null; }
function closeTag() { $("tagMenu").hidden = true; }
function toggleTag(anchor) { const m = $("tagMenu"); if (!m.hidden) return closeTag(); openTag(anchor); }
function openTag(anchor) {
  const c = tagTarget(); if (!c) { status("Start a curve, or select one in the list, before tagging."); return; }
  renderTag(c); const m = $("tagMenu"); m.hidden = false;
  const r = (anchor || $("btnTag")).getBoundingClientRect(), mw = m.offsetWidth, mh = m.offsetHeight;
  m.style.left = Math.max(8, Math.min(r.left, window.innerWidth - mw - 8)) + "px";
  m.style.top = (r.bottom + 4 + mh > window.innerHeight ? Math.max(8, r.top - mh - 4) : r.bottom + 4) + "px";
}
function renderTag(c) {
  const m = $("tagMenu"); m.innerHTML = "";
  const re = () => { renderTag(c); renderAll(); renderCurveList(); updateHint(); renderValidation(); };
  m.append(h("div", { class: "tm-title" }, "Tag: " + curveName(c)));
  const kinds = h("div", { class: "seg" });
  for (const k in KINDS) kinds.append(h("button", { class: c.kind === k ? "on" : "", onclick: () => { setKind(c, k); re(); } },
    h("span", { class: "sw", style: `background:${KINDS[k].color}` }), KINDS[k].label));
  m.append(kinds);
  if (c.kind === "mark") {
    const marks = S.doc.curves.filter((x) => x.kind === "mark" && x.id !== c.id);
    m.append(h("div", { class: "tm-row" }, h("label", {}, "Index"),
      h("input", { type: "number", min: 1, step: 1, value: c.mark_index ?? "", title: "1 = innermost preserved mark",
        onchange: (e) => { pushUndo(); c.mark_index = e.target.value === "" ? null : parseInt(e.target.value, 10); re(); } }),
      h("label", { class: "chk" }, h("input", { type: "checkbox", checked: c.index_uncertain,
        onchange: (e) => { c.index_uncertain = e.target.checked; re(); } }), "uncertain")));
    const ty = h("div", { class: "seg" });
    for (const t of MARK_TYPES) ty.append(h("button", { class: c.mark_type === t ? "on" : "", onclick: () => { c.mark_type = t; re(); } }, t.replace("_", "/")));
    m.append(h("div", { class: "tm-row" }, h("label", {}, "Type"), ty));
    const cf = h("div", { class: "seg" });
    for (const v of [3, 2, 1]) cf.append(h("button", { class: c.confidence === v ? "on" : "", title: v === 3 ? "sure" : v === 1 ? "doubtful" : "",
      onclick: () => { c.confidence = v; re(); } }, String(v)));
    m.append(h("div", { class: "tm-row" }, h("label", {}, "Confidence"), cf));
    m.append(h("div", { class: "tm-row" }, h("label", {}, ""), h("label", { class: "chk" }, h("input", { type: "checkbox", checked: c.in_efs,
      onchange: (e) => { c.in_efs = e.target.checked; re(); } }), "inside EFS")));
    if (c.mark_type === "double_multiple") {
      const sel = h("select", { onchange: (e) => { c.linked_to = e.target.value || null; re(); } });
      sel.append(h("option", { value: "" }, "—"));
      for (const x of marks) sel.append(h("option", { value: x.id, selected: c.linked_to === x.id }, curveName(x)));
      m.append(h("div", { class: "tm-row" }, h("label", {}, "Paired with"), sel));
    }
  }
  m.append(h("div", { class: "tm-foot" }, h("button", { class: "primary", onclick: closeTag }, "Done")));
}
document.addEventListener("pointerdown", (e) => { const m = $("tagMenu"); if (!m.hidden && !m.contains(e.target) && e.target !== $("btnTag")) closeTag(); });

// ---------- pan pad ----------
function panBy(fx, fy) { S.view.x += S.view.w * fx; S.view.y += S.view.h * fy; setView(); }
document.querySelectorAll("[data-pan]").forEach((b) => b.addEventListener("click", () => { if (!S.doc) return;
  const [fx, fy] = b.dataset.pan.split(",").map(Number); panBy(fx * 0.2, fy * 0.2); }));
$("padIn").addEventListener("click", () => S.doc && zoomAt([S.view.x + S.view.w / 2, S.view.y + S.view.h / 2], 1 / 1.25));
$("padOut").addEventListener("click", () => S.doc && zoomAt([S.view.x + S.view.w / 2, S.view.y + S.view.h / 2], 1.25));
$("padFit").addEventListener("click", fit);

function addVertex(p) {
  const c = curve(S.active.curveId); pushUndo();
  if (S.active.seg == null) { c.segments.push([]); S.active.seg = c.segments.length - 1; }
  c.segments[S.active.seg].push([+p[0].toFixed(2), +p[1].toFixed(2)]);
  updateHint(); renderAll(); renderCurveList();
}
function removeLastVertex() {
  if (!S.active) return; const c = curve(S.active.curveId); if (!c || !c.segments.length) return;
  pushUndo(); const si = S.active.seg == null ? c.segments.length - 1 : S.active.seg; const seg = c.segments[si];
  seg.pop(); if (!seg.length) { c.segments.splice(si, 1); S.active.seg = null; } else S.active.seg = si;
  updateHint(); renderAll(); renderCurveList();
}

// ---------- pointer ----------
svg.addEventListener("contextmenu", (e) => e.preventDefault());
svg.addEventListener("pointerdown", (ev) => {
  if (!S.doc) return; touch(); svg.focus();
  const tgt = ev.target; const p = toImg(ev);
  if (tgt.classList.contains("vert")) {
    const c = curve(tgt.dataset.c), si = +tgt.dataset.s, vi = +tgt.dataset.v;
    if (ev.button === 2) { pushUndo(); c.segments[si].splice(vi, 1); if (!c.segments[si].length) c.segments.splice(si, 1);
      if (S.active && S.active.curveId === c.id) S.active.seg = null; renderAll(); renderCurveList(); return; }
    pushUndo(); S.drag = { type: "vert", c, si, vi }; svg.setPointerCapture(ev.pointerId); return;
  }
  if (tgt.dataset && tgt.dataset.handle) { pushUndo(); S.drag = { type: "handle", which: tgt.dataset.handle }; svg.setPointerCapture(ev.pointerId); return; }
  if (ev.button === 1 || ev.button === 2 || S.space || S.tool === "pan") {
    S.drag = { type: "pan", sx: ev.clientX, sy: ev.clientY, v: Object.assign({}, S.view) };
    svg.classList.add("panning"); svg.setPointerCapture(ev.pointerId); return;
  }
  if (S.tool === "scale") {
    if (!S.scaleClick) { S.scaleClick = p; updateHint(); renderAll(); }
    else { pushUndo(); const prev = S.doc.calibration || {};
      const guess = parseScale(S.doc.task.scale_value);
      S.doc.calibration = { p1: rnd(S.scaleClick), p2: rnd(p), length: prev.length || guess.length || null, unit: prev.unit || guess.unit || "mm" };
      S.scaleClick = null; setTool("pan"); renderCalib(); }
    return;
  }
  if (S.tool === "centre") { pushUndo(); S.doc.centre = rnd(p); setTool("pan"); renderCurveList(); return; }
  if (S.tool === "trace" && S.active) addVertex(p);
});
svg.addEventListener("pointermove", (ev) => {
  if (!S.doc) return; const p = toImg(ev);
  const k = pxPerUnit();
  status(`x ${p[0].toFixed(1)}, y ${p[1].toFixed(1)} px` + (k ? ` · ${(p[0] / k).toFixed(3)}, ${(p[1] / k).toFixed(3)} ${unit()}` : "") +
         ` · zoom ${(1 / upp()).toFixed(2)}×`);
  if (!S.drag) return; touch();
  if (S.drag.type === "pan") { const f = upp(); const d = S.drag;
    S.view.x = d.v.x - (ev.clientX - d.sx) * f; S.view.y = d.v.y - (ev.clientY - d.sy) * f; setView(); }
  else if (S.drag.type === "vert") { S.drag.c.segments[S.drag.si][S.drag.vi] = rnd(p); renderAll(); }
  else if (S.drag.type === "handle") { const w = S.drag.which;
    if (w === "centre") S.doc.centre = rnd(p); else S.doc.calibration[w] = rnd(p); renderAll(); }
});
const endDrag = () => { if (S.drag && S.drag.type !== "pan") { renderCurveList(); renderCalib(); } S.drag = null; svg.classList.remove("panning"); };
svg.addEventListener("pointerup", endDrag); svg.addEventListener("pointercancel", endDrag);
svg.addEventListener("wheel", (ev) => { if (!S.doc) return; ev.preventDefault();
  zoomAt(toImg(ev), Math.exp(ev.deltaY * (ev.ctrlKey ? 0.01 : 0.0015))); }, { passive: false });
const rnd = (p) => [+p[0].toFixed(2), +p[1].toFixed(2)];

// ---------- keys ----------
document.addEventListener("keydown", (ev) => {
  const typing = /^(INPUT|TEXTAREA|SELECT)$/.test(document.activeElement.tagName);
  if ((ev.metaKey || ev.ctrlKey) && ev.key.toLowerCase() === "z" && !typing) { ev.preventDefault(); undo(); return; }
  if (typing) return;
  if (ev.key === " ") { S.space = true; svg.classList.add("pan"); ev.preventDefault(); }
  else if (ev.key === "Enter") endSegment();
  else if (ev.key === "Escape") { if (!$("tagMenu").hidden) closeTag(); else stopTrace(); }
  else if (ev.key === "t" || ev.key === "T") openTag($("btnTag"));
  else if (ev.key.startsWith("Arrow")) { ev.preventDefault(); if (!S.doc) return; const st = ev.shiftKey ? 0.5 : 0.1;
    panBy(ev.key === "ArrowLeft" ? -st : ev.key === "ArrowRight" ? st : 0, ev.key === "ArrowUp" ? -st : ev.key === "ArrowDown" ? st : 0); }
  else if (ev.key === "Backspace" || ev.key === "Delete") { ev.preventDefault(); removeLastVertex(); }
  else if (ev.key === "0") fit();
  else if (ev.key === "+" || ev.key === "=") zoomAt([S.view.x + S.view.w / 2, S.view.y + S.view.h / 2], 1 / 1.25);
  else if (ev.key === "-") zoomAt([S.view.x + S.view.w / 2, S.view.y + S.view.h / 2], 1.25);
});
document.addEventListener("keyup", (ev) => { if (ev.key === " ") { S.space = false; svg.classList.toggle("pan", S.tool === "pan"); } });

// ---------- drawing ----------
function renderAll() {
  layerShapes.innerHTML = ""; layerVerts.innerHTML = "";
  if (!S.doc) return;
  const k = upp(), sw = 2 * k, r = 4 * k;
  for (const c of S.doc.curves) {
    const isAct = S.active && S.active.curveId === c.id, isSel = S.selCurve === c.id || S.current === c.id;
    const color = KINDS[c.kind].color;
    c.segments.forEach((seg, si) => {
      if (!seg.length) return;
      el("polyline", { points: seg.map((q) => q.join(",")).join(" "), fill: "none", stroke: color,
        "stroke-width": (isAct || isSel ? 2.2 : 1.4) * sw, "stroke-linejoin": "round", "stroke-linecap": "round",
        "stroke-dasharray": c.kind === "mark" && c.mark_type === "uncertain" ? `${6 * k} ${4 * k}` : "", opacity: 0.95 }, layerShapes);
      if (S.tool === "pan" || isAct || isSel) seg.forEach((q, vi) => {
        const last = isAct && S.active.seg === si && vi === seg.length - 1;
        el("circle", { cx: q[0], cy: q[1], r: last ? r * 1.4 : r, fill: last ? "var(--sel)" : color, stroke: "#000", "stroke-width": 0.8 * k,
          class: "vert", "data-c": c.id, "data-s": si, "data-v": vi }, layerVerts); });
    });
    if (c.kind === "mark" && c.segments.length && c.segments[0].length) {
      const q = c.segments[0][0];
      const t = el("text", { x: q[0] + 6 * k, y: q[1] - 6 * k, fill: color, "font-size": 13 * k, "font-weight": 700,
        stroke: "#000", "stroke-width": 3 * k, "paint-order": "stroke" }, layerShapes); t.textContent = c.mark_index == null ? "?" : c.mark_index;
    }
  }
  const cal = S.doc.calibration;
  if (cal && cal.p1 && cal.p2) {
    el("line", { x1: cal.p1[0], y1: cal.p1[1], x2: cal.p2[0], y2: cal.p2[1], stroke: "var(--c-scale)", "stroke-width": 2 * k }, layerShapes);
    for (const w of ["p1", "p2"]) el("circle", { cx: cal[w][0], cy: cal[w][1], r: r * 1.2, fill: "var(--c-scale)", stroke: "#000",
      "stroke-width": 0.8 * k, "data-handle": w, class: "handle", style: "cursor:move" }, layerVerts);
  }
  if (S.scaleClick) el("circle", { cx: S.scaleClick[0], cy: S.scaleClick[1], r: r * 1.2, fill: "var(--c-scale)" }, layerVerts);
  const z = S.doc.centre;
  if (z) { const L = 14 * k;
    el("path", { d: `M${z[0] - L},${z[1]}H${z[0] + L}M${z[0]},${z[1] - L}V${z[1] + L}`, stroke: "var(--c-centre)", "stroke-width": 2 * k }, layerShapes);
    el("circle", { cx: z[0], cy: z[1], r: r * 1.2, fill: "transparent", stroke: "var(--c-centre)", "stroke-width": 1.5 * k,
      "data-handle": "centre", style: "cursor:move" }, layerVerts); }
  renderValidation();
}

// ---------- calibration panel ----------
function parseScale(s) { if (!s) return {}; const m = String(s).match(/([\d.]+)\s*(µm|um|μm|mm|cm|m)\b/i);
  return m ? { length: +m[1], unit: m[2].replace(/[µμ]m|um/i, "um").toLowerCase() } : {}; }
function renderCalib() {
  const box = $("calibBox"); box.innerHTML = ""; if (!S.doc) return;
  const c = S.doc.calibration;
  if (!c) { box.append(h("span", { class: "muted" }, "No scale yet: choose 'Scale bar' and click both ends.")); return; }
  const len = h("input", { type: "number", step: "any", min: "0", value: c.length ?? "", style: "width:70px",
    oninput: (e) => { c.length = parseFloat(e.target.value) || null; renderCurveList(); renderValidation(); } });
  const un = h("select", { onchange: (e) => { c.unit = e.target.value; renderCurveList(); } });
  for (const u of ["um", "mm", "cm"]) un.append(h("option", { value: u, selected: c.unit === u }, u));
  const k = pxPerUnit();
  box.append(h("div", {}, "Scale bar = ", len, " ", un),
    h("div", { class: "muted" }, `${fmt(dist(c.p1, c.p2), 1)} px` + (k ? ` · ${fmt(k, 2)} px per ${c.unit}` : " · enter the bar length")));
  if (S.doc.task.scale_value) box.append(h("div", { class: "muted" }, `Printed: ${S.doc.task.scale_value}`));
  const ctl = h("div", { class: "muted" }, "Contrast ",
    h("input", { type: "range", min: 50, max: 250, value: S.filter.contrast, oninput: (e) => { S.filter.contrast = +e.target.value; applyFilter(); } }),
    " Brightness ", h("input", { type: "range", min: 50, max: 200, value: S.filter.brightness, oninput: (e) => { S.filter.brightness = +e.target.value; applyFilter(); } }));
  box.append(ctl);
}

// ---------- curve list and live (preview) statistics ----------
function resample(seg, step) { const out = [];
  for (let i = 0; i < seg.length - 1; i++) { const a = seg[i], b = seg[i + 1], L = dist(a, b), n = Math.max(1, Math.ceil(L / step));
    for (let j = 0; j < n; j++) out.push([a[0] + (b[0] - a[0]) * j / n, a[1] + (b[1] - a[1]) * j / n]); }
  if (seg.length) out.push(seg[seg.length - 1]); return out; }
function stats(c) {
  let L = 0, nv = 0; for (const s of c.segments) { nv += s.length; for (let i = 1; i < s.length; i++) L += dist(s[i - 1], s[i]); }
  const o = { length: conv(L), nv, nseg: c.segments.length };
  const z = S.doc.centre;
  if (z && nv > 1) { const pts = c.segments.flatMap((s) => resample(s, 2)); const rs = pts.map((q) => dist(q, z));
    const bins = new Set(pts.map((q) => Math.floor(((Math.atan2(q[1] - z[1], q[0] - z[0]) * 180 / Math.PI) + 360) % 360 / 5)));
    o.rmean = conv(rs.reduce((a, b) => a + b, 0) / rs.length); o.rmin = conv(Math.min(...rs)); o.rmax = conv(Math.max(...rs)); o.span = bins.size * 5; }
  return o;
}
function curveDetail(c) {
  if (c.kind !== "mark") return "";
  return [c.mark_type.replace("_", "/"), "conf " + c.confidence, c.in_efs ? "in EFS" : null,
          c.linked_to ? "paired with " + curveName(curve(c.linked_to) || { kind: "mark", mark_index: "?" }) : null].filter(Boolean).join(" · ");
}
function renderCurveList() {
  const ul = $("curveList"); ul.innerHTML = ""; if (!S.doc) return;
  const order = ["periosteum", "efs_inner", "mark", "medullary"];
  const cs = [...S.doc.curves].sort((a, b) => order.indexOf(a.kind) - order.indexOf(b.kind) || (b.mark_index || 0) - (a.mark_index || 0));
  if (!cs.length) { ul.append(h("li", { class: "muted" }, "Nothing traced yet. Press Start.")); }
  for (const c of cs) {
    const st = stats(c), u = pxPerUnit() ? unit() : "px";
    const isAct = S.active && S.active.curveId === c.id;
    const li = h("li", { class: (isAct ? "activeCurve " : "") + (S.current === c.id ? "currentCurve" : ""), title: "Click to select",
      onclick: () => { S.current = c.id; renderAll(); renderCurveList(); },
      onmouseenter: () => { S.selCurve = c.id; renderAll(); }, onmouseleave: () => { S.selCurve = null; renderAll(); } });
    const stop = (f) => (e) => { e.stopPropagation(); f(e); };
    li.append(h("div", { class: "row" }, h("span", { class: "sw", style: `background:${KINDS[c.kind].color}` }),
      h("span", { class: "kind" }, curveName(c)), h("span", { class: "detail" }, curveDetail(c)),
      h("span", { class: "rowbtns" },
        h("button", { title: isAct ? "Stop tracing" : "Continue tracing this curve", onclick: stop(() => isAct ? stopTrace() : startTrace(c.id)) }, isAct ? "■" : "▶"),
        h("button", { title: "Tag", onclick: stop((e) => { S.current = c.id; if (!isAct) S.active = null; openTag(e.currentTarget); }) }, "✎"),
        h("button", { title: "Delete curve", onclick: stop(() => { pushUndo(); S.doc.curves = S.doc.curves.filter((x) => x.id !== c.id);
          if (isAct) S.active = null; if (S.current === c.id) S.current = null; closeTag(); renderAll(); renderCurveList(); }) }, "✕"))));
    li.append(h("div", { class: "stats" }, `${st.nseg} seg · ${st.nv} pts · ${fmt(st.length)} ${u}` +
      (st.rmean != null ? ` · r ${fmt(st.rmean)} (${fmt(st.rmin)}–${fmt(st.rmax)}) ${u} · ${st.span}°` : "")));
    ul.append(li);
  }
}

// ---------- scores ----------
function renderScores() {
  const f = $("scores"); f.innerHTML = "";
  for (const [k, lab, opts] of SCORES) {
    const sel = h("select", { onchange: (e) => { S.doc.scores[k] = e.target.value || null; touch(); renderValidation(); } });
    for (const o of opts) sel.append(h("option", { value: o, selected: (S.doc.scores[k] || "") === o }, o ? o.replace(/_/g, " ") : "—"));
    f.append(h("label", {}, lab), sel);
  }
  renderCalib(); renderCurveList();
}
$("notes").addEventListener("input", () => { if (S.doc) { S.doc.notes = $("notes").value; touch(); } });
$("measurer").addEventListener("input", () => { store.set("hm:measurer", $("measurer").value.trim()); if (S.doc) S.doc.measurer = $("measurer").value.trim() || null; renderValidation(); });

// ---------- validation and export ----------
function validate() {
  const d = S.doc, err = [], warn = [];
  if (!$("measurer").value.trim()) err.push("Enter a measurer ID (top bar).");
  const k = pxPerUnit();
  if (!k && !d.flags.no_scale) err.push(d.calibration ? "Enter the scale-bar length." : "Set the scale bar, or tick 'No scale'.");
  if (!d.centre && !d.flags.no_centre) warn.push("No centre: click it, or tick 'No centre' if the section is incomplete.");
  const marks = d.curves.filter((c) => c.kind === "mark");
  if (!marks.length) warn.push("No growth marks traced.");
  for (const m of marks) if (m.mark_index == null && !m.index_uncertain) err.push("A mark has no index (enter one or tick '?').");
  const idx = marks.filter((m) => !m.index_uncertain && m.mark_index != null).map((m) => m.mark_index);
  const dup = idx.filter((v, i) => idx.indexOf(v) !== i);
  if (dup.length) warn.push(`Mark index used twice: ${[...new Set(dup)].join(", ")} (fine for segments of one mark only if entered as one curve).`);
  for (const c of d.curves) { if (!c.segments.length) warn.push(`${curveName(c)} has no points.`);
    if (c.segments.some((s) => s.length < 2)) warn.push(`${curveName(c)} has a segment with a single point.`); }
  if (!d.curves.some((c) => c.kind === "periosteum")) warn.push("Periosteum not traced.");
  if (d.scores.efs === "present" && d.scores.efs_traceable === "yes" && !d.curves.some((c) => c.kind === "efs_inner")) warn.push("EFS marked traceable but its inner boundary is not traced.");
  return { err, warn };
}
function renderValidation() {
  const box = $("validation"); box.innerHTML = ""; if (!S.doc) return; const v = validate();
  if (!v.err.length && !v.warn.length) { box.append(h("div", { class: "ok" }, "Ready to export.")); return; }
  for (const e of v.err) box.append(h("div", { class: "bad" }, "✕ " + e));
  for (const w of v.warn) box.append(h("div", { class: "warn" }, "! " + w));
}
async function exportDoc() {
  if (!S.doc) return; saveDoc();
  const v = validate(); if (v.err.length) { renderValidation(); status("Fix the items marked ✕ before exporting."); return; }
  const d = JSON.parse(JSON.stringify(S.doc));
  d.measurer = $("measurer").value.trim(); d.exported_at = new Date().toISOString(); d.validation_warnings = v.warn;
  d.calibration_px_per_unit = pxPerUnit();
  const stamp = d.exported_at.replace(/[:.]/g, "").replace("T", "_").slice(0, 17);   // YYYY-MM-DD_HHMMSS (UTC)
  const name = `${d.task.figure_id}__${d.measurer.replace(/[^\w-]+/g, "")}__${stamp}.json`;
  const text = JSON.stringify(d, null, 1);
  let where = null;
  if (S.dirHandle) { try { const ex = await S.dirHandle.getDirectoryHandle("exports", { create: true });
      const fh = await ex.getFileHandle(name, { create: true }); const w = await fh.createWritable(); await w.write(text); await w.close(); where = "exports/" + name; }
    catch (e) { console.warn("Could not write into the folder; downloading instead", e); } }
  if (!where) { const a = document.createElement("a"); a.href = URL.createObjectURL(new Blob([text], { type: "application/json" }));
    a.download = name; document.body.append(a); a.click(); a.remove(); where = "Downloads/" + name; }
  store.set(keyStatus(d.task.figure_id), "done"); renderTasks(); status(`Exported ${where}`);
}
$("btnExport").addEventListener("click", exportDoc);
$("btnNext").addEventListener("click", () => { if (S.idx + 1 < S.tasks.length) openTask(S.idx + 1); });
$("btnUndo").addEventListener("click", undo);
$("btnFit").addEventListener("click", fit);
$("btnOpen").addEventListener("click", openFolder);

// ---------- demo image (synthetic section with drift and a broken sector) ----------
async function demo() {
  const W = 1400, H = 1000, cv = document.createElement("canvas"); cv.width = W; cv.height = H; const g = cv.getContext("2d");
  g.fillStyle = "#e9e4da"; g.fillRect(0, 0, W, H);
  const C = [660, 500], rot = 0.3, rings = 9, out = [420, 350], med = [150, 120];
  const shape = (f) => { const dx = (1 - f) * -60, dy = (1 - f) * 25;                  // centre drifts as the bone grows
    return { cx: C[0] + dx, cy: C[1] + dy, rx: med[0] + (out[0] - med[0]) * f, ry: med[1] + (out[1] - med[1]) * f }; };
  const ell = (e) => { g.beginPath(); g.ellipse(e.cx, e.cy, e.rx, e.ry, rot, 0, Math.PI * 2); };
  ell(shape(1)); g.fillStyle = "#c9a77c"; g.fill();
  for (let i = 0; i < 2500; i++) { const a = Math.random() * 6.283, f = 0.05 + Math.random() * 0.93, e = shape(f);
    const x = e.cx + e.rx * Math.cos(a) * Math.cos(rot) - e.ry * Math.sin(a) * Math.sin(rot), y = e.cy + e.rx * Math.cos(a) * Math.sin(rot) + e.ry * Math.sin(a) * Math.cos(rot);
    g.fillStyle = "rgba(70,45,25,.55)"; g.beginPath(); g.ellipse(x, y, 1.5 + Math.random() * 2.5, 1 + Math.random() * 1.5, a, 0, 6.283); g.fill(); }
  const fs = [...Array(rings)].map((_, i) => { const t = i + 1; return 1 - Math.exp(-0.32 * (t + 1.2)); });
  const fmax = fs[fs.length - 1]; g.strokeStyle = "#5b3a1e"; g.lineWidth = 2.2;
  for (const f of fs) { ell(shape(f / fmax * 0.97)); g.stroke(); }
  ell(shape(0)); g.fillStyle = "#4a3b30"; g.fill();
  g.save(); g.translate(C[0], C[1]); g.beginPath(); g.moveTo(0, 0);           // broken sector
  for (let a = -0.9; a <= 0.55; a += 0.05) g.lineTo(Math.cos(a) * 700 + (Math.random() - .5) * 25, Math.sin(a) * 700 + (Math.random() - .5) * 25);
  g.closePath(); g.fillStyle = "#e9e4da"; g.fill(); g.restore();
  g.fillStyle = "#111"; g.fillRect(1120, 930, 200, 8); g.font = "22px sans-serif"; g.fillText("1 mm", 1190, 922);
  const blob = await new Promise((r) => cv.toBlob(r, "image/png"));
  const file = new File([blob], "demo_section.png", { type: "image/png" });
  saveDoc(); S.dirHandle = null; S.files = new Map([["demo_section.png", file]]);
  S.tasks = [{ figure_id: "demo_section", file: "demo_section.png", source_id: "synthetic", citation: "Synthetic demo (generated in the browser)",
    taxon: "Synthetica demonstrata", specimen: "DEMO-1", element: "femur", view: "partial_section", light: "normal", quality_tier: "B",
    scale_value: "1 mm", caption: "Synthetic mid-diaphyseal section: 9 growth marks, centre drifting with growth, one broken sector. Try: scale bar → centre → periosteum → marks (inner = 1)." }];
  S.idx = -1; renderTasks(); openTask(0);
}
$("btnDemo").addEventListener("click", demo);

// ---------- start ----------
$("ver").textContent = "v" + TOOL_VERSION;
$("measurer").value = store.get("hm:measurer") || "";
renderValidation(); updateHint();
})();
