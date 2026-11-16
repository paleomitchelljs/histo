/* Histology Measure 0.3
 * Static, no build step, no server. Images are read from a local folder and never leave the machine.
 * Geometry is stored in IMAGE PIXELS. Derived measurements (lengths, radii, coverage, zone order, lacuna density) are
 * computed at export and written under `derived`; the screen shows only what the measurer needs to act on.
 */
(function () {
"use strict";
const TOOL_VERSION = "0.3.0";
const SCHEMA_VERSION = "0.3";
const IMG_EXT = /\.(png|jpe?g|webp|gif|bmp)$/i;

// Curve kinds. `boundary` kinds bound tissue zones; count_area only bounds a lacuna count.
const KINDS = {
  periosteum: { label: "Periosteum", short: "Peri", color: "var(--c-periosteum)", boundary: true },
  medullary:  { label: "Medullary margin", short: "Med", color: "var(--c-medullary)", boundary: true },
  mark:       { label: "Growth mark", short: "", color: "var(--c-mark)", boundary: true },
  label:      { label: "Fluorochrome label", short: "L", color: "var(--c-label)", boundary: true },
  efs_inner:  { label: "EFS inner boundary", short: "EFS", color: "var(--c-efs_inner)", boundary: true },
  count_area: { label: "Count area", short: "Area", color: "var(--c-count)", boundary: false, closed: true },
  ignore:     { label: "Ignore box", short: "Ign", color: "var(--c-ignore)", boundary: false, closed: true },
};
const TAGGABLE = ["mark", "label", "periosteum", "medullary", "efs_inner"];
const MARK_TYPES = ["LAG", "annulus", "uncertain"];
const DYES = ["calcein", "alizarin", "tetracycline", "xylenol", "other"];
const MARK_FIELDS = ["mark_index", "index_uncertain", "mark_type", "confidence", "in_efs", "branch_of"];
const LABEL_FIELDS = ["label_index", "dye", "label_date", "age_days"];

// Zone descriptors: components, not composite tissue names.
const ZONE_FIELDS = [
  ["matrix", "Matrix", ["woven", "parallel_fibred", "lamellar", "mixed", "unclear"]],
  ["vasc_density", "Vascular density", ["none", "low", "moderate", "high", "unclear"]],
  ["vasc_orient", "Canal orientation", ["longitudinal", "circumferential", "radial", "reticular", "plexiform", "mixed", "unclear"]],
  ["remodelling", "Secondary osteons", ["none", "scattered", "dense", "unclear"]],
];
const SECTION_FIELDS = [
  ["position", "Position", ["midshaft", "proximal_diaphysis", "distal_diaphysis", "metaphysis", "unknown"]],
  ["plane", "Plane", ["transverse", "oblique", "longitudinal"]],
  ["min_circumference", "Min. circ.", ["yes", "no", "unknown"]],
];
const FLAGS = [["no_scale", "No scale"], ["no_centre", "No centre"], ["outer_cortex_eroded", "Outer eroded"],
               ["deformation_visible", "Deformed"], ["cracks_infill", "Cracks"]];
const SCORES = [
  ["efs", "EFS / OCL", ["present", "absent", "unclear"]],
  ["efs_traceable", "EFS traceable", ["yes", "no", "n/a"]],
  ["medullary_bone", "Medullary bone", ["present", "absent", "unclear"]],
  ["kastschenko_line", "Kastschenko line", ["present", "absent", "unclear"]],
  ["image_quality", "Image quality", ["good", "fair", "poor"]],
];
const IDENTITY = ["source_id", "citation", "taxon", "specimen", "element", "view", "light", "quality_tier", "scale_value"];

const $ = (id) => document.getElementById(id);
const svg = $("svg"), layerImg = $("layerImg"), layerShapes = $("layerShapes"), layerVerts = $("layerVerts");
const NS = "http://www.w3.org/2000/svg";
const S = { tasks: [], idx: -1, files: new Map(), dirHandle: null, doc: null, tool: "pan", active: null,
            view: { x: 0, y: 0, w: 1000, h: 800 }, img: null, undo: [], scaleClick: null, drag: null,
            space: false, lastActivity: Date.now(), hover: null, current: null, filter: { contrast: 100, brightness: 100 },
            popFor: null };

// ---------- storage (per-viewer convenience only; the exported file is the record) ----------
// Extension points (features.js): open/loaded image, drawing, label-tool clicks, export.
const hooks = { open: [], loaded: [], render: [], label: [], beforeExport: [], exportFiles: [], derive: [] };
const store = {
  get(k) { try { const v = localStorage.getItem(k); return v ? JSON.parse(v) : null; } catch (e) { return null; } },
  set(k, v) { try { localStorage.setItem(k, JSON.stringify(v)); return true; } catch (e) { return false; } },
  del(k) { try { localStorage.removeItem(k); } catch (e) { /* storage unavailable */ } },
};
const keyDoc = (fid) => "hm:doc:" + fid, keyStatus = (fid) => "hm:status:" + fid;

// ---------- helpers ----------
const el = (tag, attrs = {}, parent) => { const e = document.createElementNS(NS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]); if (parent) parent.appendChild(e); return e; };
const h = (tag, props = {}, ...kids) => { const e = document.createElement(tag);
  for (const k in props) { if (k === "class") e.className = props[k]; else if (k.startsWith("on")) e.addEventListener(k.slice(2), props[k]);
    else if (k in e && k !== "list") e[k] = props[k]; else e.setAttribute(k, props[k]); }
  for (const c of kids) if (c != null) e.append(c); return e; };
const dist = (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1]);
const uid = () => Math.random().toString(36).slice(2, 9);
const nice = (s) => String(s).replace(/_/g, " ");
const rnd = (p) => [+p[0].toFixed(2), +p[1].toFixed(2)];
let statusTimer = null;
function status(msg, sticky) { $("status").textContent = msg; clearTimeout(statusTimer);
  if (!sticky) statusTimer = setTimeout(() => { $("status").textContent = ""; }, 4000); }
function hint(msg) { $("hint").textContent = msg || ""; }
function touch() { S.lastActivity = Date.now(); }

function pxPerUnit(doc = S.doc) { const c = doc && doc.calibration;
  return c && c.p1 && c.p2 && c.length > 0 ? dist(c.p1, c.p2) / c.length : null; }

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
    catch (e) { status("manifest.json could not be read: " + e.message, true); }
  }
  if (!tasks.length) {
    tasks = [...files.keys()].filter((n) => IMG_EXT.test(n)).sort().map((n) => ({ figure_id: n.replace(/\.[^.]+$/, ""), file: n, _no_manifest: true }));
    if (!mf) status("No manifest: one task per image.");
  }
  for (const t of tasks) t._missing = !files.has(t.file);
  S.tasks = tasks; S.idx = -1; renderTasks();
  const first = tasks.findIndex((t) => !t._missing && store.get(keyStatus(t.figure_id)) !== "done");
  if (tasks.length) openTask(first >= 0 ? first : 0);
}

function renderTasks() {
  const ul = $("taskList"); ul.innerHTML = "";
  const done = S.tasks.filter((t) => store.get(keyStatus(t.figure_id)) === "done").length;
  $("taskCount").textContent = S.tasks.length ? `${done}/${S.tasks.length}` : "";
  if (!S.tasks.length) { ul.append(h("li", { class: "muted" }, "Open a folder.")); return; }
  S.tasks.forEach((t, i) => {
    const st = store.get(keyStatus(t.figure_id));
    ul.append(h("li", { class: i === S.idx ? "current" : "", onclick: () => openTask(i), title: t.file },
      h("span", { class: "dot " + (st === "done" ? "done" : st === "prog" ? "prog" : "") }),
      h("span", {}, t.figure_id + (t._missing ? " (missing)" : ""))));
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
  if (t._missing) { status(`Not in folder: ${t.file}`); return; }
  S.idx = i; S.active = null; S.undo = []; S.scaleClick = null; S.hover = null; S.current = null; closePop();
  const file = S.files.get(t.file);
  if (S.img && S.img.url) URL.revokeObjectURL(S.img.url);
  const url = URL.createObjectURL(file);
  const dims = await new Promise((res, rej) => { const im = new Image(); im.onload = () => res([im.naturalWidth, im.naturalHeight]);
    im.onerror = () => rej(new Error("cannot display this image format")); im.src = url; }).catch((e) => { status(`${t.file}: ${e.message}`, true); return null; });
  if (!dims) return;
  S.img = { url, w: dims[0], h: dims[1] };
  S.doc = migrate(store.get(keyDoc(t.figure_id))) || newDoc(t);
  S.doc.image.px_width = dims[0]; S.doc.image.px_height = dims[1];
  const hash = await sha256(file);
  S.doc.image.sha256 = hash;
  S.doc.image.sha256_matches_manifest = (t.sha256 && hash) ? t.sha256 === hash : null;
  if (S.doc.image.sha256_matches_manifest === false) status("Image differs from the manifest hash.", true);
  layerImg.innerHTML = "";
  el("image", { href: url, x: 0, y: 0, width: dims[0], height: dims[1], id: "theImage" }, layerImg);
  applyFilter(); renderTasks(); renderMeta(); renderPanels(); fit(); setTool("pan");
  for (const f of hooks.open) f();
}

function newDoc(t) {
  const task = {}; for (const k in t) if (!k.startsWith("_")) task[k] = t[k];
  return { schema_version: SCHEMA_VERSION, tool_version: TOOL_VERSION, task, identity_edits: {},
           image: { file: t.file, sha256: null, px_width: null, px_height: null },
           measurer: $("measurer").value.trim() || null, repeat: !!(t.assignment && t.assignment.repeat),
           started_at: new Date().toISOString(), time_spent_s: 0,
           calibration: null, centre: null, section: {}, flags: {}, curves: [], lacunae: [], scores: {}, notes: "",
           features: { settings: {}, objects: [] } };
}
// 0.1 -> 0.2: double/multiple marks with a pair become branches; "section oblique" becomes the plane.
// 0.2 -> 0.3: adds `features` (labelled canals, osteons etc.) and the `ignore` curve kind.
function migrate(d) {
  if (!d) return null;
  if (d.schema_version === SCHEMA_VERSION) return d;
  if (d.schema_version === "0.2") { d.features = d.features || { settings: {}, objects: [] }; d.schema_version = SCHEMA_VERSION; d.tool_version = TOOL_VERSION; return d; }
  if (d.schema_version !== "0.1") return null;
  for (const c of d.curves || []) if (c.kind === "mark") {
    if (c.mark_type === "double_multiple") { c.mark_type = "LAG"; c.branch_of = c.linked_to || null; }
    delete c.linked_to; if (!("branch_of" in c)) c.branch_of = null;
  }
  d.section = d.section || {};
  if (d.flags && d.flags.section_oblique) { d.section.plane = "oblique"; delete d.flags.section_oblique; }
  d.lacunae = d.lacunae || [];
  d.migrated_from = "0.1"; d.schema_version = SCHEMA_VERSION; d.tool_version = TOOL_VERSION;
  d.features = { settings: {}, objects: [] };
  return d;
}
function saveDoc() {
  if (!S.doc) return;
  S.doc.notes = $("notes").value;
  let ok = store.set(keyDoc(S.doc.task.figure_id), S.doc);
  if (!ok) {   // storage full (labelled features are large): drop drafts of images already exported, then retry
    for (const t of S.tasks) if (t.figure_id !== S.doc.task.figure_id && store.get(keyStatus(t.figure_id)) === "done") store.del(keyDoc(t.figure_id));
    ok = store.set(keyDoc(S.doc.task.figure_id), S.doc);
  }
  if (!ok && !S.warnedStorage) {
    S.warnedStorage = true; status("Browser storage is full: autosave is off for this image. Export to keep your work.", true); }
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
    dl.append(h("dt", {}, nice(k)));
    if (v != null && v !== "") dl.append(h("dd", {}, String(v)));
    else dl.append(h("dd", {}, h("input", { value: S.doc.identity_edits[k] || "", placeholder: "—",
      oninput: (e) => { S.doc.identity_edits[k] = e.target.value; touch(); } })));
  }
  $("caption").textContent = t.caption || "";
  $("captionBox").hidden = !t.caption;
  $("notes").value = S.doc.notes || "";
}

// ---------- view ----------
function svgSize() { const r = svg.getBoundingClientRect(); return [Math.max(r.width, 1), Math.max(r.height, 1)]; }
function setView() { const v = S.view; svg.setAttribute("viewBox", `${v.x} ${v.y} ${v.w} ${v.h}`);
  const im = $("theImage"); if (im) im.style.imageRendering = upp() < 0.5 ? "pixelated" : "auto"; renderAll(); }
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
const zoomMid = (f) => S.doc && zoomAt([S.view.x + S.view.w / 2, S.view.y + S.view.h / 2], f);
function applyFilter() { const im = $("theImage"); if (im) im.style.filter = `contrast(${S.filter.contrast}%) brightness(${S.filter.brightness}%)`; }
$("contrast").addEventListener("input", (e) => { S.filter.contrast = +e.target.value; applyFilter(); });
$("brightness").addEventListener("input", (e) => { S.filter.brightness = +e.target.value; applyFilter(); });

// ---------- undo ----------
function snapshot() { const d = S.doc; return JSON.stringify({ c: d.calibration, z: d.centre, v: d.curves, l: d.lacunae }); }
function pushUndo() { if (!S.doc) return; S.undo.push(snapshot()); if (S.undo.length > 200) S.undo.shift(); touch(); }
function undo() { const u = S.undo.pop(); if (!u) return; const o = JSON.parse(u);
  Object.assign(S.doc, { calibration: o.c, centre: o.z, curves: o.v, lacunae: o.l });
  if (S.active && !curve(S.active.curveId)) S.active = null; renderAll(); renderPanels(); }

// ---------- tools ----------
function setTool(t) {
  S.tool = t; S.scaleClick = null;
  if (t !== "trace") S.active = null;
  document.querySelectorAll("[data-tool]").forEach((b) => b.classList.toggle("active", b.dataset.tool === t));
  svg.classList.toggle("pan", t === "pan");
  updateHint(); renderAll(); renderCurveList();
}
function updateHint() {
  if (!S.doc) return hint("Open a folder, or try the demo.");
  if (S.tool === "scale") return hint(S.scaleClick ? "Other end of the scale bar" : "One end of the scale bar");
  if (S.tool === "centre") return hint("Click the centre");
  if (S.tool === "count") return hint("Click lacunae · Esc to stop");
  if (S.tool === "ignore") return hint("Drag a box over panels, photos or labels to leave out · Esc to stop");
  if (S.tool === "label") return hint("Click a feature to label it · hover + 1–6 / X · Esc to stop");
  if (S.tool === "trace" && S.active) return hint(`Tracing ${curveName(curve(S.active.curveId))} · Enter = gap · Esc = stop`);
  hint("");
}
document.querySelectorAll("[data-tool]").forEach((b) => b.addEventListener("click", () => S.doc && setTool(b.dataset.tool)));

const curve = (id) => S.doc.curves.find((c) => c.id === id);
function curveName(c) {
  if (!c) return "?";
  if (c.kind === "mark") { const p = c.branch_of && curve(c.branch_of);
    return p ? `Mark ${p.mark_index ?? "?"} branch` : `Mark ${c.mark_index ?? "?"}${c.index_uncertain ? "?" : ""}`; }
  if (c.kind === "label") return `Label ${c.label_index ?? "?"}`;
  return KINDS[c.kind].label;
}
function shortName(c) {
  if (!c) return "·";
  if (c.kind === "mark") { const p = c.branch_of && curve(c.branch_of);
    return p ? `${p.mark_index ?? "?"}b` : `${c.mark_index ?? "?"}${c.index_uncertain ? "?" : ""}`; }
  if (c.kind === "label") return `L${c.label_index ?? "?"}`;
  return KINDS[c.kind].short;
}
const nextIndex = (kind, field, self) => S.doc.curves.filter((x) => x.kind === kind && x.id !== self.id && !x.branch_of)
  .reduce((m, x) => Math.max(m, x[field] || 0), 0) + 1;
function setDefaults(c) {
  if (c.kind === "mark") Object.assign(c, { mark_index: nextIndex("mark", "mark_index", c), index_uncertain: false,
    mark_type: "LAG", confidence: 3, in_efs: false, branch_of: null });
  if (c.kind === "label") Object.assign(c, { label_index: nextIndex("label", "label_index", c), dye: null, label_date: "", age_days: null });
}
function setKind(c, kind) {
  if (c.kind === kind) return; pushUndo();
  for (const k of [...MARK_FIELDS, ...LABEL_FIELDS]) delete c[k];
  for (const x of S.doc.curves) if (x.branch_of === c.id) x.branch_of = null;
  c.kind = kind; setDefaults(c);
}
function newCurve(kind = "mark") {
  if (!S.doc) return; pushUndo();
  const c = { id: uid(), kind, segments: [] }; setDefaults(c); S.doc.curves.push(c);
  startTrace(c.id);
}
function startTrace(id) { closePop(); S.tool = "trace"; S.active = { curveId: id, seg: null }; S.current = id;
  document.querySelectorAll("[data-tool]").forEach((b) => b.classList.remove("active")); svg.classList.remove("pan");
  updateHint(); renderAll(); renderPanels(); svg.focus(); }
function endSegment() { if (S.active) { S.active.seg = null; updateHint(); renderAll(); } }
function stopTrace() { endSegment(); setTool("pan"); renderPanels(); }
$("btnStart").addEventListener("click", () => newCurve());
$("btnStop").addEventListener("click", stopTrace);
$("btnTag").addEventListener("click", (e) => { e.stopPropagation(); openTag($("btnTag")); });
$("btnArea").addEventListener("click", () => newCurve("count_area"));

// ---------- popover (tags, zones) ----------
function closePop() { $("pop").hidden = true; S.popFor = null; }
function showPop(anchor, render, key) {
  if (!$("pop").hidden && S.popFor === key) return closePop();
  S.popFor = key; render(); const m = $("pop"); m.hidden = false;
  const r = anchor.getBoundingClientRect(), mw = m.offsetWidth, mh = m.offsetHeight;
  m.style.left = Math.max(8, Math.min(r.left, window.innerWidth - mw - 8)) + "px";
  m.style.top = (r.bottom + 4 + mh > window.innerHeight ? Math.max(8, r.top - mh - 4) : r.bottom + 4) + "px";
}
const seg = (opts, cur, onPick, fmt = nice) => { const d = h("div", { class: "seg" });
  for (const o of opts) d.append(h("button", { class: cur === o ? "on" : "", onclick: () => onPick(cur === o ? null : o) }, fmt(o)));
  return d; };
const row = (label, ...kids) => h("div", { class: "tm-row" }, h("label", {}, label), ...kids);

function tagTarget() { const id = S.active ? S.active.curveId : S.current; return id && S.doc ? curve(id) : null; }
function openTag(anchor) {
  const c = tagTarget(); if (!c) { status("Select or start a curve first."); return; }
  showPop(anchor, () => renderTag(c), "tag:" + c.id);
}
function renderTag(c) {
  const m = $("pop"); m.innerHTML = "";
  const re = () => { renderTag(c); renderAll(); renderPanels(); updateHint(); };
  m.append(h("div", { class: "tm-title" }, curveName(c)));
  if (c.kind !== "count_area") {
    const kinds = h("div", { class: "seg" });
    for (const k of TAGGABLE) kinds.append(h("button", { class: c.kind === k ? "on" : "", onclick: () => { setKind(c, k); re(); } },
      h("span", { class: "sw", style: `background:${KINDS[k].color}` }), KINDS[k].label));
    m.append(kinds);
  }
  if (c.kind === "mark") {
    const others = S.doc.curves.filter((x) => x.kind === "mark" && x.id !== c.id && !x.branch_of);
    if (!c.branch_of) m.append(row("Index",
      h("input", { type: "number", min: 1, step: 1, value: c.mark_index ?? "", title: "1 = innermost preserved mark",
        onchange: (e) => { pushUndo(); c.mark_index = e.target.value === "" ? null : parseInt(e.target.value, 10); re(); } }),
      h("label", { class: "chk" }, h("input", { type: "checkbox", checked: c.index_uncertain,
        onchange: (e) => { c.index_uncertain = e.target.checked; re(); } }), "?")));
    m.append(row("Type", seg(MARK_TYPES, c.mark_type, (v) => { c.mark_type = v || "LAG"; re(); }, String)));
    m.append(row("Confidence", seg([3, 2, 1], c.confidence, (v) => { c.confidence = v || 3; re(); }, String)));
    const sel = h("select", { onchange: (e) => { pushUndo(); c.branch_of = e.target.value || null;
      if (c.branch_of) { c.mark_index = curve(c.branch_of).mark_index; c.index_uncertain = false; } re(); } },
      h("option", { value: "" }, "no (main trace)"));
    for (const x of others) sel.append(h("option", { value: x.id, selected: c.branch_of === x.id }, curveName(x)));
    m.append(row("Branch of", sel));
    m.append(row("", h("label", { class: "chk" }, h("input", { type: "checkbox", checked: c.in_efs,
      onchange: (e) => { c.in_efs = e.target.checked; re(); } }), "inside EFS")));
  }
  if (c.kind === "label") {
    m.append(row("Index", h("input", { type: "number", min: 1, step: 1, value: c.label_index ?? "",
      onchange: (e) => { c.label_index = e.target.value === "" ? null : parseInt(e.target.value, 10); re(); } })));
    m.append(row("Dye", seg(DYES, c.dye, (v) => { c.dye = v; re(); })));
    m.append(row("Date", h("input", { type: "date", value: c.label_date || "", onchange: (e) => { c.label_date = e.target.value; re(); } })));
    m.append(row("Age (days)", h("input", { type: "number", min: 0, step: 1, value: c.age_days ?? "", title: "Animal's age at injection",
      onchange: (e) => { c.age_days = e.target.value === "" ? null : +e.target.value; re(); } })));
  }
  m.append(h("div", { class: "tm-foot" }, h("button", { class: "primary", onclick: closePop }, "Done")));
}
document.addEventListener("pointerdown", (e) => {
  const m = $("pop"); if (!m.hidden && !m.contains(e.target) && !e.target.closest("[data-pop]") && e.target !== $("btnTag")) closePop();
  const hp = $("help"); if (!hp.hidden && !hp.contains(e.target) && e.target !== $("btnHelp")) hp.hidden = true;
});
$("btnHelp").addEventListener("click", () => { const hp = $("help"); hp.hidden = !hp.hidden; });

// ---------- geometry ----------
function resample(seg, step) { const out = [];
  for (let i = 0; i < seg.length - 1; i++) { const a = seg[i], b = seg[i + 1], L = dist(a, b), n = Math.max(1, Math.ceil(L / step));
    for (let j = 0; j < n; j++) out.push([a[0] + (b[0] - a[0]) * j / n, a[1] + (b[1] - a[1]) * j / n]); }
  if (seg.length) out.push(seg[seg.length - 1]); return out; }
const allPts = (c) => c.segments.flat();
function sample(c, max = 60) { const p = c.segments.flatMap((s) => resample(s, 2)); if (p.length <= max) return p;
  const k = p.length / max; return [...Array(max)].map((_, i) => p[Math.floor(i * k)]); }
function distToPolyline(p, segs) { let best = Infinity;
  for (const s of segs) for (let i = 0; i < s.length - 1; i++) { const a = s[i], b = s[i + 1];
    const dx = b[0] - a[0], dy = b[1] - a[1], L2 = dx * dx + dy * dy;
    const t = L2 ? Math.max(0, Math.min(1, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / L2)) : 0;
    best = Math.min(best, Math.hypot(p[0] - a[0] - t * dx, p[1] - a[1] - t * dy)); }
  return best; }
const angle = (q, z) => ((Math.atan2(q[1] - z[1], q[0] - z[0]) * 180 / Math.PI) + 360) % 360;
// Angular coverage (5-degree bins) and the arc occupied, as [start, end] degrees clockwise in image coordinates.
function angular(c, z) {
  const pts = c.segments.flatMap((s) => resample(s, 2)); if (!z || pts.length < 2) return null;
  const bins = new Array(72).fill(false); for (const q of pts) bins[Math.floor(angle(q, z) / 5) % 72] = true;
  const n = bins.filter(Boolean).length; if (n === 72) return { coverage_deg: 360, arc_deg: [0, 360] };
  let gapStart = 0, gapLen = 0;                       // the largest empty run marks where the arc is not
  for (let i = 0; i < 72; i++) if (!bins[i]) { let L = 0; while (L < 72 && !bins[(i + L) % 72]) L++; if (L > gapLen) { gapLen = L; gapStart = i; } }
  return { coverage_deg: n * 5, arc_deg: [((gapStart + gapLen) % 72) * 5, (gapStart * 5) % 360 || 360] };
}
function polyArea(pts) { let a = 0; for (let i = 0; i < pts.length; i++) { const p = pts[i], q = pts[(i + 1) % pts.length]; a += p[0] * q[1] - q[0] * p[1]; } return Math.abs(a) / 2; }
function inPoly(p, pts) { let inside = false;
  for (let i = 0, j = pts.length - 1; i < pts.length; j = i++) { const a = pts[i], b = pts[j];
    if ((a[1] > p[1]) !== (b[1] > p[1]) && p[0] < (b[0] - a[0]) * (p[1] - a[1]) / (b[1] - a[1]) + a[0]) inside = !inside; }
  return inside; }

// Boundaries in radial order, innermost first. Order by: distance from the centre if set; otherwise distance from the
// periosteum (works for wedges with no centre); otherwise by kind and index.
function orderedBoundaries(doc = S.doc) {
  const bs = doc.curves.filter((c) => KINDS[c.kind].boundary && !c.branch_of && c.segments.some((s) => s.length));
  const z = doc.centre, peri = doc.curves.find((c) => c.kind === "periosteum" && c.segments.some((s) => s.length > 1));
  let key, method;
  if (z) { method = "centre"; key = (c) => { const p = sample(c); return p.reduce((a, q) => a + dist(q, z), 0) / p.length; }; }
  else if (peri) { method = "periosteum_distance"; key = (c) => { if (c === peri) return Infinity; const p = sample(c);
    return -p.reduce((a, q) => a + distToPolyline(q, peri.segments), 0) / p.length; }; }
  else { method = "index"; const rank = { medullary: 0, mark: 1, label: 1, efs_inner: 2, periosteum: 3 };
    key = (c) => rank[c.kind] * 1000 + (c.mark_index || c.label_index || 0); }
  const k = new Map(bs.map((c) => [c.id, key(c)]));
  return { method, list: bs.sort((a, b) => k.get(a.id) - k.get(b.id)) };
}
// Zones: the tissue just inside each boundary. Descriptors are stored on the outer boundary (`zone`), so they stay
// with the curve when traces are added or reordered.
function zones(doc = S.doc) {
  const { method, list } = orderedBoundaries(doc), out = [];
  list.forEach((c, i) => { if (c.kind === "medullary") return; out.push({ inner: list[i - 1] || null, outer: c }); });
  return { method, zones: out };
}

// ---------- pan pad ----------
function panBy(fx, fy) { S.view.x += S.view.w * fx; S.view.y += S.view.h * fy; setView(); }
document.querySelectorAll("[data-pan]").forEach((b) => b.addEventListener("click", () => { if (!S.doc) return;
  const [fx, fy] = b.dataset.pan.split(",").map(Number); panBy(fx * 0.2, fy * 0.2); }));
$("padIn").addEventListener("click", () => zoomMid(1 / 1.25));
$("padOut").addEventListener("click", () => zoomMid(1.25));
$("padFit").addEventListener("click", fit);

function addVertex(p) {
  const c = curve(S.active.curveId); pushUndo();
  if (S.active.seg == null) { c.segments.push([]); S.active.seg = c.segments.length - 1; }
  c.segments[S.active.seg].push(rnd(p));
  renderAll(); renderPanels();
}
function removeLastVertex() {
  if (!S.active) return; const c = curve(S.active.curveId); if (!c || !c.segments.length) return;
  pushUndo(); const si = S.active.seg == null ? c.segments.length - 1 : S.active.seg; const sg = c.segments[si];
  sg.pop(); if (!sg.length) { c.segments.splice(si, 1); S.active.seg = null; } else S.active.seg = si;
  renderAll(); renderPanels();
}

// ---------- pointer ----------
svg.addEventListener("contextmenu", (e) => e.preventDefault());
svg.addEventListener("pointerdown", (ev) => {
  if (!S.doc) return; touch(); svg.focus();
  const tgt = ev.target; const p = toImg(ev);
  if (tgt.classList.contains("lac") && ev.button === 2) { pushUndo(); S.doc.lacunae.splice(+tgt.dataset.i, 1); renderAll(); renderCount(); return; }
  if (tgt.classList.contains("vert")) {
    const c = curve(tgt.dataset.c), si = +tgt.dataset.s, vi = +tgt.dataset.v;
    if (ev.button === 2) { pushUndo(); c.segments[si].splice(vi, 1); if (!c.segments[si].length) c.segments.splice(si, 1);
      if (S.active && S.active.curveId === c.id) S.active.seg = null; renderAll(); renderPanels(); return; }
    if (S.tool === "pan" || (S.active && S.active.curveId === c.id && ev.altKey)) {
      pushUndo(); S.drag = { type: "vert", c, si, vi }; svg.setPointerCapture(ev.pointerId); return; }
  }
  if (tgt.dataset && tgt.dataset.handle && S.tool === "pan") { pushUndo(); S.drag = { type: "handle", which: tgt.dataset.handle }; svg.setPointerCapture(ev.pointerId); return; }
  if (ev.button === 1 || ev.button === 2 || S.space || S.tool === "pan") {
    S.drag = { type: "pan", sx: ev.clientX, sy: ev.clientY, v: Object.assign({}, S.view) };
    svg.classList.add("panning"); svg.setPointerCapture(ev.pointerId); return;
  }
  if (S.tool === "scale") {
    if (!S.scaleClick) { S.scaleClick = p; updateHint(); renderAll(); }
    else { pushUndo(); const prev = S.doc.calibration || {}, guess = parseScale(S.doc.task.scale_value);
      S.doc.calibration = { p1: rnd(S.scaleClick), p2: rnd(p), length: prev.length || guess.length || null, unit: prev.unit || guess.unit || "mm" };
      S.scaleClick = null; setTool("pan"); renderCalib(); if (!S.doc.calibration.length) $("calLen").focus(); }
    return;
  }
  if (S.tool === "centre") { pushUndo(); S.doc.centre = rnd(p); setTool("pan"); renderPanels(); return; }
  if (S.tool === "count") { pushUndo(); S.doc.lacunae.push(rnd(p)); renderAll(); renderCount(); return; }
  if (S.tool === "ignore") { S.drag = { type: "box", a: p, b: p }; svg.setPointerCapture(ev.pointerId); return; }
  if (S.tool === "label") { for (const f of hooks.label) f(p, ev); return; }
  if (S.tool === "trace" && S.active) addVertex(p);
});
svg.addEventListener("pointermove", (ev) => {
  if (!S.doc || !S.drag) return; touch(); const p = toImg(ev);
  if (S.drag.type === "pan") { const f = upp(); const d = S.drag;
    S.view.x = d.v.x - (ev.clientX - d.sx) * f; S.view.y = d.v.y - (ev.clientY - d.sy) * f; setView(); }
  else if (S.drag.type === "vert") { S.drag.c.segments[S.drag.si][S.drag.vi] = rnd(p); renderAll(); }
  else if (S.drag.type === "box") { S.drag.b = p; renderAll(); }
  else if (S.drag.type === "handle") { const w = S.drag.which;
    if (w === "centre") S.doc.centre = rnd(p); else S.doc.calibration[w] = rnd(p); renderAll(); }
});
const endDrag = () => {
  const d = S.drag;
  if (d && d.type === "box") { const [x0, x1] = [d.a[0], d.b[0]].sort((a, b) => a - b), [y0, y1] = [d.a[1], d.b[1]].sort((a, b) => a - b);
    if (x1 - x0 > 3 && y1 - y0 > 3) { pushUndo(); S.doc.curves.push({ id: uid(), kind: "ignore", segments: [[rnd([x0, y0]), rnd([x1, y0]), rnd([x1, y1]), rnd([x0, y1])]] }); } }
  if (d && d.type !== "pan") renderPanels(); S.drag = null; svg.classList.remove("panning"); renderAll(); };
svg.addEventListener("pointerup", endDrag); svg.addEventListener("pointercancel", endDrag);
svg.addEventListener("wheel", (ev) => { if (!S.doc) return; ev.preventDefault();
  zoomAt(toImg(ev), Math.exp(ev.deltaY * (ev.ctrlKey ? 0.01 : 0.0015))); }, { passive: false });

// ---------- keys ----------
document.addEventListener("keydown", (ev) => {
  const typing = /^(INPUT|TEXTAREA|SELECT)$/.test(document.activeElement.tagName);
  if ((ev.metaKey || ev.ctrlKey) && ev.key.toLowerCase() === "z" && !typing) { ev.preventDefault(); undo(); return; }
  if (typing || ev.metaKey || ev.ctrlKey) return;
  if (ev.key === " ") { S.space = true; svg.classList.add("pan"); ev.preventDefault(); }
  else if (ev.key === "Enter") endSegment();
  else if (ev.key === "Escape") { if (!$("pop").hidden) closePop(); else if (!$("help").hidden) $("help").hidden = true; else stopTrace(); }
  else if (ev.key === "t" || ev.key === "T") openTag($("btnTag"));
  else if (ev.key === "n" || ev.key === "N") newCurve();
  else if (ev.key === "?") $("help").hidden = !$("help").hidden;
  else if (ev.key.startsWith("Arrow")) { ev.preventDefault(); if (!S.doc) return; const st = ev.shiftKey ? 0.5 : 0.1;
    panBy(ev.key === "ArrowLeft" ? -st : ev.key === "ArrowRight" ? st : 0, ev.key === "ArrowUp" ? -st : ev.key === "ArrowDown" ? st : 0); }
  else if (ev.key === "Backspace" || ev.key === "Delete") { ev.preventDefault(); removeLastVertex(); }
  else if (ev.key === "0") fit();
  else if (ev.key === "+" || ev.key === "=") zoomMid(1 / 1.25);
  else if (ev.key === "-") zoomMid(1.25);
});
document.addEventListener("keyup", (ev) => { if (ev.key === " ") { S.space = false; svg.classList.toggle("pan", S.tool === "pan"); } });

// ---------- drawing ----------
function renderAll() {
  layerShapes.innerHTML = ""; layerVerts.innerHTML = "";
  if (!S.doc) return;
  const k = upp(), sw = 2 * k, r = 4 * k;
  const hiZone = S.hover && S.hover.startsWith("zone:") ? S.hover.slice(5).split("|") : [];
  for (const c of S.doc.curves) {
    const isAct = S.active && S.active.curveId === c.id;
    const isSel = S.hover === c.id || S.current === c.id || hiZone.includes(c.id);
    const color = KINDS[c.kind].color, closed = KINDS[c.kind].closed;
    c.segments.forEach((sg, si) => {
      if (!sg.length) return;
      const pts = sg.map((q) => q.join(",")).join(" ");
      el(closed ? "polygon" : "polyline", { points: pts, fill: c.kind === "ignore" ? "rgba(0,0,0,.45)" : closed ? "rgba(148,163,184,.12)" : "none", stroke: color,
        "stroke-width": (isAct || isSel ? 2.2 : 1.4) * sw, "stroke-linejoin": "round", "stroke-linecap": "round",
        "stroke-dasharray": (c.kind === "mark" && (c.mark_type === "uncertain" || c.branch_of)) || closed ? `${6 * k} ${4 * k}` : "",
        opacity: 0.95 }, layerShapes);
      if (S.tool === "pan" || isAct || isSel) sg.forEach((q, vi) => {
        const last = isAct && S.active.seg === si && vi === sg.length - 1;
        el("circle", { cx: q[0], cy: q[1], r: last ? r * 1.4 : r, fill: last ? "var(--sel)" : color, stroke: "#000", "stroke-width": 0.8 * k,
          class: "vert", "data-c": c.id, "data-s": si, "data-v": vi }, layerVerts); });
    });
    if ((c.kind === "mark" || c.kind === "label") && c.segments.length && c.segments[0].length) {
      const q = c.segments[0][0];
      const t = el("text", { x: q[0] + 6 * k, y: q[1] - 6 * k, fill: color, "font-size": 13 * k, "font-weight": 700,
        stroke: "#000", "stroke-width": 3 * k, "paint-order": "stroke" }, layerShapes); t.textContent = shortName(c);
    }
  }
  S.doc.lacunae.forEach((q, i) => el("circle", { cx: q[0], cy: q[1], r: r * 0.9, fill: "var(--c-lac)", stroke: "#000", "stroke-width": 0.6 * k,
    class: "lac", "data-i": i }, layerVerts));
  const cal = S.doc.calibration;
  if (cal && cal.p1 && cal.p2) {
    el("line", { x1: cal.p1[0], y1: cal.p1[1], x2: cal.p2[0], y2: cal.p2[1], stroke: "var(--c-scale)", "stroke-width": 2 * k }, layerShapes);
    for (const w of ["p1", "p2"]) el("circle", { cx: cal[w][0], cy: cal[w][1], r: r * 1.2, fill: "var(--c-scale)", stroke: "#000",
      "stroke-width": 0.8 * k, "data-handle": w, class: "handle" }, layerVerts);
  }
  if (S.drag && S.drag.type === "box") { const { a, b } = S.drag;
    el("rect", { x: Math.min(a[0], b[0]), y: Math.min(a[1], b[1]), width: Math.abs(b[0] - a[0]), height: Math.abs(b[1] - a[1]),
      fill: "rgba(0,0,0,.35)", stroke: "var(--c-ignore)", "stroke-width": 1.5 * k, "stroke-dasharray": `${6 * k} ${4 * k}` }, layerShapes); }
  if (S.scaleClick) el("circle", { cx: S.scaleClick[0], cy: S.scaleClick[1], r: r * 1.2, fill: "var(--c-scale)" }, layerVerts);
  const z = S.doc.centre;
  if (z) { const L = 14 * k;
    el("path", { d: `M${z[0] - L},${z[1]}H${z[0] + L}M${z[0]},${z[1] - L}V${z[1] + L}`, stroke: "var(--c-centre)", "stroke-width": 2 * k }, layerShapes);
    el("circle", { cx: z[0], cy: z[1], r: r * 1.2, fill: "transparent", stroke: "var(--c-centre)", "stroke-width": 1.5 * k,
      "data-handle": "centre", class: "handle" }, layerVerts); }
  for (const f of hooks.render) f();
}

// ---------- panels ----------
function renderPanels() { if (!S.doc) return; renderCalib(); renderSection(); renderCurveList(); renderZones(); renderCount(); renderScores(); renderValidation(); }

function parseScale(s) { if (!s) return {}; const m = String(s).match(/([\d.]+)\s*(µm|um|μm|mm|cm|m)\b/i);
  return m ? { length: +m[1], unit: m[2].replace(/[µμ]m|um/i, "um").toLowerCase() } : {}; }
function renderCalib() {
  const box = $("calibBox"); box.innerHTML = ""; const c = S.doc.calibration;
  if (!c) { box.append(h("span", { class: "muted" }, S.doc.flags.no_scale ? "No scale" : "Scale not set")); return; }
  const un = h("select", { onchange: (e) => { c.unit = e.target.value; renderValidation(); } });
  for (const u of ["um", "mm", "cm"]) un.append(h("option", { value: u, selected: c.unit === u }, u === "um" ? "µm" : u));
  box.append(h("span", {}, "Bar "), h("input", { id: "calLen", type: "number", step: "any", min: "0", value: c.length ?? "", class: "num",
    oninput: (e) => { c.length = parseFloat(e.target.value) || null; renderValidation(); } }), un,
    S.doc.task.scale_value ? h("span", { class: "muted", title: "Printed scale value" }, ` (${S.doc.task.scale_value})`) : null);
}
function renderSection() {
  const g = $("section"); g.innerHTML = ""; const s = S.doc.section;
  for (const [k, lab, opts] of SECTION_FIELDS) {
    const sel = h("select", { onchange: (e) => { s[k] = e.target.value || null; touch(); renderValidation(); } }, h("option", { value: "" }, "—"));
    for (const o of opts) sel.append(h("option", { value: o, selected: s[k] === o }, nice(o)));
    g.append(h("label", {}, lab), sel);
  }
  const f = $("flags"); f.innerHTML = "";
  for (const [k, lab] of FLAGS) f.append(h("button", { class: "chip" + (S.doc.flags[k] ? " on" : ""),
    onclick: () => { S.doc.flags[k] = !S.doc.flags[k]; touch(); renderSection(); renderCalib(); renderValidation(); } }, lab));
}
function curveDetail(c) {
  if (c.kind === "mark") return [c.mark_type !== "LAG" ? c.mark_type : null, c.confidence < 3 ? "conf " + c.confidence : null, c.in_efs ? "EFS" : null].filter(Boolean).join(" · ");
  if (c.kind === "label") return [c.dye, c.label_date || (c.age_days != null ? c.age_days + " d" : null)].filter(Boolean).join(" · ");
  if (c.kind === "count_area") { const n = S.doc.lacunae.filter((q) => c.segments[0] && c.segments[0].length > 2 && inPoly(q, c.segments[0])).length; return n + " lacunae"; }
  return "";
}
function renderCurveList() {
  const ul = $("curveList"); ul.innerHTML = "";
  const order = ["periosteum", "efs_inner", "label", "mark", "medullary", "count_area", "ignore"];
  const cs = [...S.doc.curves].sort((a, b) => order.indexOf(a.kind) - order.indexOf(b.kind) || (b.mark_index || b.label_index || 0) - (a.mark_index || a.label_index || 0));
  if (!cs.length) { ul.append(h("li", { class: "muted empty" }, "Press New (N) and click along a curve.")); return; }
  for (const c of cs) {
    const isAct = S.active && S.active.curveId === c.id, empty = !c.segments.some((s) => s.length);
    const stop = (f) => (e) => { e.stopPropagation(); f(e); };
    ul.append(h("li", { class: (isAct ? "activeCurve " : "") + (S.current === c.id ? "currentCurve" : ""),
        onclick: () => { S.current = c.id; renderAll(); renderCurveList(); },
        onmouseenter: () => { S.hover = c.id; renderAll(); }, onmouseleave: () => { S.hover = null; renderAll(); } },
      h("span", { class: "sw", style: `background:${KINDS[c.kind].color}` }),
      h("span", { class: "kind" + (c.branch_of ? " branch" : "") }, curveName(c)),
      h("span", { class: "detail" }, empty ? "no points" : curveDetail(c)),
      h("span", { class: "rowbtns" },
        h("button", { title: isAct ? "Stop" : "Continue tracing", onclick: stop(() => isAct ? stopTrace() : startTrace(c.id)) }, isAct ? "■" : "▶"),
        h("button", { title: "Tag", "data-pop": 1, onclick: stop((e) => { S.current = c.id; if (!isAct) S.active = null; openTag(e.currentTarget); }) }, "✎"),
        h("button", { title: "Delete", onclick: stop(() => { pushUndo(); S.doc.curves = S.doc.curves.filter((x) => x.id !== c.id);
          for (const x of S.doc.curves) if (x.branch_of === c.id) x.branch_of = null;
          if (isAct) S.active = null; if (S.current === c.id) S.current = null; closePop(); renderAll(); renderPanels(); }) }, "✕"))));
  }
}
function zoneSummary(z) { const d = z || {};
  const parts = ZONE_FIELDS.map(([k]) => d[k]).filter(Boolean).map(nice); if (d.varies) parts.push("varies");
  return parts.length ? parts.join(" · ") : "—"; }
function renderZones() {
  const ul = $("zoneList"); ul.innerHTML = "";
  const { zones: zs } = zones();
  if (!zs.length) { ul.append(h("li", { class: "muted empty" }, "Appear once boundaries are traced.")); return; }
  for (const z of [...zs].reverse()) {
    const key = (z.inner ? z.inner.id : "") + "|" + z.outer.id;
    ul.append(h("li", { "data-pop": 1, class: z.outer.zone && Object.values(z.outer.zone).some(Boolean) ? "" : "unset",
        onclick: (e) => showPop(e.currentTarget, () => renderZonePop(z), "zone:" + key),
        onmouseenter: () => { S.hover = "zone:" + key; renderAll(); }, onmouseleave: () => { S.hover = null; renderAll(); } },
      h("span", { class: "zname" }, `${shortName(z.inner)} → ${shortName(z.outer)}`),
      h("span", { class: "detail" }, zoneSummary(z.outer.zone))));
  }
}
function renderZonePop(z) {
  const m = $("pop"); m.innerHTML = ""; const c = z.outer; c.zone = c.zone || {};
  const re = () => { renderZonePop(z); renderZones(); touch(); };
  m.append(h("div", { class: "tm-title" }, `Zone ${shortName(z.inner)} → ${shortName(z.outer)}`));
  for (const [k, lab, opts] of ZONE_FIELDS) m.append(h("div", { class: "tm-block" }, h("label", {}, lab), seg(opts, c.zone[k], (v) => { c.zone[k] = v; re(); })));
  m.append(h("label", { class: "chk" }, h("input", { type: "checkbox", checked: !!c.zone.varies,
    onchange: (e) => { c.zone.varies = e.target.checked; re(); } }), "varies around the section (say how in Notes)"));
  const inner = z.inner && z.inner.zone;
  m.append(h("div", { class: "tm-foot" },
    h("button", { disabled: !inner, title: "Copy the descriptors of the zone just inside", onclick: () => { c.zone = Object.assign({}, inner); re(); } }, "Same as inner"),
    h("button", { class: "primary", onclick: closePop }, "Done")));
}
function renderCount() {
  const n = S.doc.lacunae.length, areas = S.doc.curves.filter((c) => c.kind === "count_area").length;
  $("countInfo").textContent = n || areas ? `${n} counted${areas ? "" : " · no area"}` : "";
}
function renderScores() {
  const f = $("scores"); f.innerHTML = "";
  for (const [k, lab, opts] of SCORES) {
    const sel = h("select", { onchange: (e) => { S.doc.scores[k] = e.target.value || null; touch(); renderValidation(); } }, h("option", { value: "" }, "—"));
    for (const o of opts) sel.append(h("option", { value: o, selected: S.doc.scores[k] === o }, nice(o)));
    f.append(h("label", {}, lab), sel);
  }
}
$("notes").addEventListener("input", () => { if (S.doc) { S.doc.notes = $("notes").value; touch(); } });
$("measurer").addEventListener("input", () => { store.set("hm:measurer", $("measurer").value.trim()); if (S.doc) S.doc.measurer = $("measurer").value.trim() || null; renderValidation(); });

// ---------- validation ----------
function validate() {
  const d = S.doc, err = [], warn = [];
  if (!$("measurer").value.trim()) err.push("Measurer ID");
  if (!pxPerUnit() && !d.flags.no_scale) err.push(d.calibration ? "Scale-bar length" : "Scale bar (or No scale)");
  if (!d.centre && !d.flags.no_centre) warn.push("Centre (or No centre)");
  if (!d.section.position) warn.push("Section position");
  const marks = d.curves.filter((c) => c.kind === "mark" && !c.branch_of);
  if (!marks.length) warn.push("No growth marks");
  if (marks.some((m) => m.mark_index == null && !m.index_uncertain)) err.push("Mark without index");
  const idx = marks.filter((m) => !m.index_uncertain && m.mark_index != null).map((m) => m.mark_index);
  const dup = [...new Set(idx.filter((v, i) => idx.indexOf(v) !== i))];
  if (dup.length) warn.push(`Index used twice: ${dup.join(", ")}`);
  if (d.curves.some((c) => !c.segments.some((s) => s.length))) warn.push("Curve with no points");
  if (d.curves.some((c) => c.segments.some((s) => s.length === 1))) warn.push("Single-point segment");
  if (!d.curves.some((c) => c.kind === "periosteum")) warn.push("Periosteum not traced");
  if (d.scores.efs === "present" && d.scores.efs_traceable === "yes" && !d.curves.some((c) => c.kind === "efs_inner")) warn.push("EFS inner boundary not traced");
  if (d.curves.some((c) => c.kind === "label" && !c.label_date && c.age_days == null)) warn.push("Label without date or age");
  if (d.lacunae.length && !d.curves.some((c) => c.kind === "count_area")) warn.push("Lacunae counted without an area");
  return { err, warn };
}
function renderValidation() {
  const box = $("validation"); box.innerHTML = ""; if (!S.doc) return; const v = validate();
  if (!v.err.length && !v.warn.length) { box.append(h("div", { class: "ok" }, "Ready")); return; }
  if (v.err.length) box.append(h("div", { class: "bad" }, "Needed: " + v.err.join(" · ")));
  if (v.warn.length) box.append(h("div", { class: "warn" }, "Check: " + v.warn.join(" · ")));
}

// ---------- export: raw geometry plus derived measurements ----------
function derive(d) {
  const k = pxPerUnit(d), u = k ? d.calibration.unit : "px", cv = (px) => (k ? px / k : px);
  const z = d.centre, num = (v) => (v == null || !isFinite(v) ? null : +v.toFixed(4));
  const curves = d.curves.map((c) => {
    const closed = KINDS[c.kind].closed;
    let L = 0, nv = 0; for (const s of c.segments) { nv += s.length; for (let i = 1; i < s.length; i++) L += dist(s[i - 1], s[i]);
      if (closed && s.length > 2) L += dist(s[s.length - 1], s[0]); }
    const o = { id: c.id, kind: c.kind, name: curveName(c), n_segments: c.segments.length, n_points: nv, [closed ? "perimeter" : "length"]: num(cv(L)) };
    if (z && nv > 1 && !closed) { const pts = c.segments.flatMap((s) => resample(s, 2)), rs = pts.map((q) => dist(q, z));
      Object.assign(o, { radius_mean: num(cv(rs.reduce((a, b) => a + b, 0) / rs.length)), radius_min: num(cv(Math.min(...rs))),
        radius_max: num(cv(Math.max(...rs))) }, angular(c, z)); }
    if (KINDS[c.kind].closed && c.segments[0] && c.segments[0].length > 2) o.area = num(cv(cv(polyArea(c.segments[0]))));
    return o;
  });
  const zz = zones(d);
  const zoneRows = zz.zones.map((x, i) => Object.assign({ zone: i + 1, inner_id: x.inner ? x.inner.id : null, inner: x.inner ? curveName(x.inner) : null,
    outer_id: x.outer.id, outer: curveName(x.outer) }, x.outer.zone || {}));
  const counts = d.curves.filter((c) => c.kind === "count_area" && c.segments[0] && c.segments[0].length > 2).map((c) => {
    const n = d.lacunae.filter((q) => inPoly(q, c.segments[0])).length, a = cv(cv(polyArea(c.segments[0])));
    return { area_id: c.id, n_lacunae: n, area: num(a), density_per_unit2: num(a ? n / a : null) };
  });
  return { note: "Computed at export as a convenience. Analysis should recompute from the raw points.",
    unit: u, px_per_unit: num(k), reference_for_zone_order: zz.method, curves, zones: zoneRows, lacuna_counts: counts,
    lacunae_outside_areas: d.lacunae.filter((q) => !d.curves.some((c) => c.kind === "count_area" && c.segments[0] && c.segments[0].length > 2 && inPoly(q, c.segments[0]))).length };
}
function deriveAll(d) { const out = derive(d); for (const f of hooks.derive) f(d, out); return out; }
async function writeOut(name, data, type) {
  if (S.dirHandle) { try { const ex = await S.dirHandle.getDirectoryHandle("exports", { create: true });
      const fh = await ex.getFileHandle(name, { create: true }); const w = await fh.createWritable(); await w.write(data); await w.close(); return "exports/" + name; }
    catch (e) { console.warn("Could not write into the folder; downloading instead", e); } }
  const a = document.createElement("a"); a.href = URL.createObjectURL(data instanceof Blob ? data : new Blob([data], { type }));
  a.download = name; document.body.append(a); a.click(); a.remove(); return "Downloads/" + name;
}
async function exportDoc() {
  if (!S.doc) return; saveDoc();
  const v = validate(); if (v.err.length) { renderValidation(); status("Needed before export: " + v.err.join(", ")); return; }
  const at = new Date().toISOString(), measurer = $("measurer").value.trim();
  const stamp = at.replace(/[:.]/g, "").replace("T", "_").slice(0, 17);   // YYYY-MM-DD_HHMMSS (UTC)
  const base = `${S.doc.task.figure_id}__${measurer.replace(/[^\w-]+/g, "")}__${stamp}`;
  for (const f of hooks.beforeExport) f(S.doc, base);
  const d = JSON.parse(JSON.stringify(S.doc));
  d.measurer = measurer; d.exported_at = at; d.validation_warnings = v.warn;
  d.derived = deriveAll(S.doc);
  const where = await writeOut(base + ".json", JSON.stringify(d, null, 1), "application/json");
  for (const f of hooks.exportFiles) { const x = await f(S.doc, base); if (x) await writeOut(x.name, x.blob); }
  store.set(keyStatus(d.task.figure_id), "done"); renderTasks(); status(`Saved ${where}`);
}
// Reopen an earlier export of this image (and its mask, if chosen too) to review or relabel it. The next export is a
// new file that records which export it was based on.
async function loadExport(files) {
  const js = files.find((f) => /\.json$/i.test(f.name)), mask = files.find((f) => /\.png$/i.test(f.name));
  if (!js) { status("Choose an export .json (and its .mask.png, if any)."); return; }
  let d; try { d = JSON.parse(await js.text()); } catch (e) { status("Not readable: " + e.message, true); return; }
  const fid = d.task && d.task.figure_id, i = S.tasks.findIndex((t) => t.figure_id === fid);
  if (i < 0) { status(`Open the folder that holds ${fid || "this image"} first.`, true); return; }
  if (i !== S.idx) await openTask(i);
  d = migrate(d); if (!d) { status("Unknown export version.", true); return; }
  d.based_on = { file: js.name, measurer: d.measurer || null, exported_at: d.exported_at || null };
  for (const k of ["derived", "exported_at", "validation_warnings"]) delete d[k];
  d.measurer = $("measurer").value.trim() || null; d.started_at = new Date().toISOString(); d.time_spent_s = 0;
  d.image = Object.assign({}, d.image, S.doc.image);
  S.doc = d; S.undo = []; S.current = null; S.active = null; closePop();
  renderMeta(); renderPanels(); renderAll(); saveDoc();
  for (const f of hooks.loaded) f(mask || null);
  status(`Loaded ${js.name}${mask ? " with its mask" : ""}. Export saves a new file.`);
}
$("btnLoad").addEventListener("click", () => { if (!S.doc) { status("Open the image folder first."); return; } $("loadInput").click(); });
$("loadInput").addEventListener("change", async (e) => { await loadExport([...e.target.files]); e.target.value = ""; });
$("btnExport").addEventListener("click", exportDoc);
$("btnNext").addEventListener("click", () => { if (S.idx + 1 < S.tasks.length) openTask(S.idx + 1); });
$("btnUndo").addEventListener("click", undo);
$("btnOpen").addEventListener("click", openFolder);

// ---------- demo image (synthetic section with drift, a broken sector and one branching mark) ----------
async function demo() {
  const W = 1400, H = 1000, cv = document.createElement("canvas"); cv.width = W; cv.height = H; const g = cv.getContext("2d");
  g.fillStyle = "#e9e4da"; g.fillRect(0, 0, W, H);
  const C = [660, 500], rot = 0.3, rings = 9, out = [420, 350], med = [150, 120];
  const shape = (f) => { const dx = (1 - f) * -60, dy = (1 - f) * 25;                  // centre drifts as the bone grows
    return { cx: C[0] + dx, cy: C[1] + dy, rx: med[0] + (out[0] - med[0]) * f, ry: med[1] + (out[1] - med[1]) * f }; };
  const ell = (e, a0 = 0, a1 = Math.PI * 2) => { g.beginPath(); g.ellipse(e.cx, e.cy, e.rx, e.ry, rot, a0, a1); };
  ell(shape(1)); g.fillStyle = "#c9a77c"; g.fill();
  for (let i = 0; i < 2500; i++) { const a = Math.random() * 6.283, f = 0.05 + Math.random() * 0.93, e = shape(f);
    const x = e.cx + e.rx * Math.cos(a) * Math.cos(rot) - e.ry * Math.sin(a) * Math.sin(rot), y = e.cy + e.rx * Math.cos(a) * Math.sin(rot) + e.ry * Math.sin(a) * Math.cos(rot);
    g.fillStyle = "rgba(70,45,25,.55)"; g.beginPath(); g.ellipse(x, y, 1.5 + Math.random() * 2.5, 1 + Math.random() * 1.5, a, 0, 6.283); g.fill(); }
  const fs = [...Array(rings)].map((_, i) => 1 - Math.exp(-0.32 * (i + 2.2)));
  const fmax = fs[fs.length - 1]; g.strokeStyle = "#5b3a1e"; g.lineWidth = 2.2;
  fs.forEach((f, i) => { ell(shape(f / fmax * 0.97)); g.stroke();
    if (i === 5) { ell(shape(f / fmax * 0.97 + 0.018), 1.6, 3.3); g.stroke(); } });   // mark 6 branches over part of the cortex
  ell(shape(0)); g.fillStyle = "#4a3b30"; g.fill();
  g.save(); g.translate(C[0], C[1]); g.beginPath(); g.moveTo(0, 0);           // broken sector
  for (let a = -0.9; a <= 0.55; a += 0.05) g.lineTo(Math.cos(a) * 700 + (Math.random() - .5) * 25, Math.sin(a) * 700 + (Math.random() - .5) * 25);
  g.closePath(); g.fillStyle = "#e9e4da"; g.fill(); g.restore();
  g.fillStyle = "#111"; g.fillRect(1120, 930, 200, 8); g.font = "22px sans-serif"; g.fillText("1 mm", 1190, 922);
  const blob = await new Promise((r) => cv.toBlob(r, "image/png"));
  const file = new File([blob], "demo_section.png", { type: "image/png" });
  saveDoc(); S.dirHandle = null; S.files = new Map([["demo_section.png", file]]);
  S.tasks = [{ figure_id: "demo_section", file: "demo_section.png", source_id: "synthetic", citation: "Synthetic demo",
    taxon: "Synthetica demonstrata", specimen: "DEMO-1", element: "femur", view: "partial_section", light: "normal", quality_tier: "B",
    scale_value: "1 mm", caption: "Synthetic midshaft section: 9 growth marks, mark 6 branching over part of the cortex, drifting centre, one broken sector." }];
  S.idx = -1; renderTasks(); openTask(0);
}
$("btnDemo").addEventListener("click", demo);

// ---------- API for features.js ----------
window.HM = { S, $, h, el, svg, layerVerts, hooks, status, touch, uid, seg, row, showPop, closePop, pxPerUnit, KINDS, curveName,
  resample, angle, angular, dist, inPoly, toImg, upp, setView, zoomAt, setTool, renderAll, renderPanels };

// ---------- start ----------
$("ver").textContent = TOOL_VERSION;
$("measurer").value = store.get("hm:measurer") || "";
updateHint();
})();
