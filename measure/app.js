/* Histology Measure 0.4
 * Static, no build step, no server. Images are read from a local folder and never leave the machine.
 * Geometry is stored in IMAGE PIXELS. Derived measurements are computed at export and written under `derived`.
 * The screen is a sequence of steps (Image, Scale, Outline, Marks, Tissue, Whole bone, More, Export). Each step can be
 * skipped with a reason, and the reason is exported: the model treats "absent" and "not observable" differently.
 * Tissue questions and their interpretation rules live in guide.js (window.HG).
 */
(function () {
"use strict";
const TOOL_VERSION = "0.4.0";
const SCHEMA_VERSION = "0.4";
const IMG_EXT = /\.(png|jpe?g|webp|gif|bmp)$/i;
const HG = window.HG;

// Curve kinds. `boundary` kinds bound tissue zones; `outer` kinds are the outer edge (always outermost).
const KINDS = {
  periosteum:    { label: "Periosteum", short: "Peri", color: "var(--c-periosteum)", boundary: true, outer: true },
  broken_edge:   { label: "Broken edge", short: "Brk", color: "var(--c-broken)", boundary: true, outer: true },
  medullary:     { label: "Medullary margin", short: "Med", color: "var(--c-medullary)", boundary: true },
  remodel_front: { label: "Remodelling front", short: "Rem", color: "var(--c-remodel)", boundary: true },
  mark:          { label: "Growth mark", short: "", color: "var(--c-mark)", boundary: true },
  label:         { label: "Fluorochrome label", short: "L", color: "var(--c-label)", boundary: true },
  efs_inner:     { label: "EFS inner boundary", short: "EFS", color: "var(--c-efs_inner)", boundary: true },
  transect:      { label: "Transect", short: "T", color: "var(--c-transect)", boundary: false },
  count_area:    { label: "Count area", short: "Area", color: "var(--c-count)", boundary: false, closed: true },
  ignore:        { label: "Ignore box", short: "Ign", color: "var(--c-ignore)", boundary: false, closed: true },
};
const OUTER = new Set(["periosteum", "broken_edge"]);
const TAGGABLE = ["mark", "label", "efs_inner", "periosteum", "broken_edge", "medullary", "remodel_front"];
const MARK_TYPES = [["LAG", "LAG (thin line)"], ["annulus", "Annulus (band)"], ["uncertain", "Not sure"]];
const DYES = ["calcein", "alizarin", "tetracycline", "xylenol", "other"];
const MARK_FIELDS = ["mark_index", "index_auto", "index_uncertain", "mark_type", "confidence", "in_efs", "branch_of"];
const LABEL_FIELDS = ["label_index", "dye", "label_date", "age_days"];
const VIEWS = [["whole", "A whole cross-section"], ["part", "Part of a section (wedge or fragment)"],
               ["detail", "A close-up of the tissue"], ["unusable", "Nothing usable"]];
const LIGHTS = [["plain", "Plain light"], ["xpl", "Crossed polars"], ["lambda", "Polars + lambda plate"],
                ["fluorescence", "Fluorescence"], ["ct", "CT / synchrotron"], ["other", "Other"]];
const POLARIZED = new Set(["xpl", "lambda"]);
const SECTION_FIELDS = [
  ["position", "Position on the shaft", ["midshaft", "proximal_diaphysis", "distal_diaphysis", "metaphysis", "unknown"]],
  ["plane", "Plane", ["transverse", "oblique", "longitudinal"]],
  ["min_circumference", "Minimum circumference", ["yes", "no", "unknown"]],
];
const FLAGS = [["deformation_visible", "Crushed or deformed"], ["cracks_infill", "Cracks or infill"]];
const IDENTITY = ["source_id", "citation", "taxon", "specimen", "element", "view", "light", "quality_tier", "scale_value"];
const XING_KINDS = [["mark", "Growth mark"], ["efs_inner", "EFS inner boundary"], ["label", "Fluorochrome label"]];
const ENDS = { outer: [["periosteum", "Periosteum"], ["broken_edge", "Broken edge"], ["image_edge", "Image edge"], ["unclear", "Unclear"]],
               inner: [["medullary", "Medullary margin"], ["remodel_front", "Remodelling front"], ["centre", "Centre"], ["image_edge", "Image edge"], ["unclear", "Unclear"]] };

const STEPS = [
  { id: "image", rail: "Image", title: "Image" },
  { id: "scale", rail: "Scale", title: "Scale",
    reasons: ["no scale bar in the figure", "scale bar cut off or unreadable", "other"] },
  { id: "outline", rail: "Outline", title: "Outline",
    reasons: ["edges not in the image", "too damaged to trace", "other"] },
  { id: "marks", rail: "Marks", title: "Growth marks",
    reasons: ["marks too faint to place", "image too coarse", "other"] },
  { id: "tissue", rail: "Tissue", title: "Tissue",
    reasons: ["image too coarse to see canals", "tissue obscured (stain, damage)", "other"] },
  { id: "bone", rail: "Bone", title: "Whole bone" },
  { id: "more", rail: "More", title: "More (optional)" },
  { id: "export", rail: "Export", title: "Export" },
];

const $ = (id) => document.getElementById(id);
const svg = $("svg"), layerImg = $("layerImg"), layerShapes = $("layerShapes"), layerVerts = $("layerVerts");
const NS = "http://www.w3.org/2000/svg";
const S = { tasks: [], idx: -1, files: new Map(), dirHandle: null, doc: null, tool: "pan", active: null,
            view: { x: 0, y: 0, w: 1000, h: 800 }, img: null, comp: null, blend: 0, undo: [], scaleClick: null, drag: null,
            space: false, lastActivity: Date.now(), hover: null, current: null, filter: { contrast: 100, brightness: 100 },
            popFor: null, step: "image", openZone: null, tStart: null, activeT: null, align: null };

// ---------- storage (per-viewer convenience only; the exported file is the record) ----------
// Extension points (features.js): open/loaded image, drawing, label-tool clicks, export.
const hooks = { open: [], loaded: [], render: [], label: [], beforeExport: [], exportFiles: [], derive: [] };
const store = {
  get(k) { try { const v = localStorage.getItem(k); return v ? JSON.parse(v) : null; } catch (e) { return null; } },
  set(k, v) { try { localStorage.setItem(k, JSON.stringify(v)); return true; } catch (e) { return false; } },
  del(k) { try { localStorage.removeItem(k); } catch (e) { /* storage unavailable */ } },
};
const keyDoc = (fid) => "hm:doc:" + fid, keyStatus = (fid) => "hm:status:" + fid, keyComp = (file) => "hm:comp:" + file;

// ---------- helpers ----------
const el = (tag, attrs = {}, parent) => { const e = document.createElementNS(NS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]); if (parent) parent.appendChild(e); return e; };
const h = (tag, props = {}, ...kids) => { const e = document.createElement(tag);
  for (const k in props) { if (props[k] == null) continue; if (k === "class") e.className = props[k]; else if (k.startsWith("on")) e.addEventListener(k.slice(2), props[k]);
    else if (k in e && k !== "list") e[k] = props[k]; else e.setAttribute(k, props[k]); }
  for (const c of kids) if (c != null && c !== false) e.append(c); return e; };
const dist = (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1]);
const uid = () => Math.random().toString(36).slice(2, 9);
const nice = (s) => String(s).replace(/_/g, " ");
const rnd = (p) => [+p[0].toFixed(2), +p[1].toFixed(2)];
const hasPts = (c) => c.segments.some((s) => s.length);
let statusTimer = null;
function status(msg, sticky) { $("status").textContent = msg; clearTimeout(statusTimer);
  if (!sticky) statusTimer = setTimeout(() => { $("status").textContent = ""; }, 5000); }
function hint(msg) { $("hint").textContent = msg || ""; }
function touch() { S.lastActivity = Date.now(); }

function pxPerUnit(doc = S.doc) { const c = doc && doc.calibration;
  return c && c.p1 && c.p2 && c.length > 0 ? dist(c.p1, c.p2) / c.length : null; }
const zoneDone = (zone, ctx) => !!zone && (!!zone.legacy || HG.complete(zone.obs, ctx));
const isPolarized = (d = S.doc) => POLARIZED.has(d.setup.light) || !!(d.setup.companion && POLARIZED.has(d.setup.companion.light));

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
  const first = tasks.findIndex((t) => !t._missing && store.get(keyStatus(t.figure_id)) !== "done" && !isCompanionTask(t));
  if (tasks.length) openTask(first >= 0 ? first : 0);
}
const isCompanionTask = (t) => { const o = store.get(keyComp(t.file)); return o && o !== t.figure_id; };

function renderTasks() {
  const ul = $("taskList"); ul.innerHTML = "";
  const done = S.tasks.filter((t) => store.get(keyStatus(t.figure_id)) === "done").length;
  $("taskCount").textContent = S.tasks.length ? `${done}/${S.tasks.length}` : "";
  if (!S.tasks.length) { ul.append(h("li", { class: "muted" }, "Open a folder.")); return; }
  S.tasks.forEach((t, i) => {
    const st = store.get(keyStatus(t.figure_id)), comp = isCompanionTask(t);
    ul.append(h("li", { class: (i === S.idx ? "current" : "") + (comp ? " companion" : ""), onclick: () => openTask(i),
        title: comp ? `${t.file}: other-light image of ${store.get(keyComp(t.file))}` : t.file },
      h("span", { class: "dot " + (st === "done" ? "done" : st === "prog" ? "prog" : "") }),
      h("span", {}, t.figure_id + (t._missing ? " (missing)" : "") + (comp ? " ◑" : ""))));
  });
}

async function sha256(file) {
  if (!(window.crypto && crypto.subtle)) return null;
  const buf = await file.arrayBuffer(); const d = await crypto.subtle.digest("SHA-256", buf);
  return [...new Uint8Array(d)].map((b) => b.toString(16).padStart(2, "0")).join("");
}
const imgDims = (url) => new Promise((res, rej) => { const im = new Image(); im.onload = () => res([im.naturalWidth, im.naturalHeight]);
  im.onerror = () => rej(new Error("cannot display this image format")); im.src = url; });

async function openTask(i) {
  if (i < 0 || i >= S.tasks.length) return;
  saveDoc();
  const t = S.tasks[i];
  if (t._missing) { status(`Not in folder: ${t.file}`); return; }
  S.idx = i; S.active = null; S.undo = []; S.scaleClick = null; S.hover = null; S.current = null; S.tStart = null; S.activeT = null;
  S.align = null; S.openZone = null; closePop();
  const file = S.files.get(t.file);
  if (S.img && S.img.url) URL.revokeObjectURL(S.img.url);
  const url = URL.createObjectURL(file);
  const dims = await imgDims(url).catch((e) => { status(`${t.file}: ${e.message}`, true); return null; });
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
  await loadCompanion();
  $("welcome").hidden = true; document.body.classList.remove("nodoc");
  applyFilter(); renderTasks(); fit(); setTool("pan");
  S.step = firstOpenStep();
  renderPanels();
  for (const f of hooks.open) f();
}

const guessView = (v) => ({ full_section: "whole", partial_section: "part", detail: "detail" })[v] || null;
function guessLight(s) { s = String(s || "");
  if (/lambda|gypsum/i.test(s)) return "lambda"; if (/pol|xpl|cross/i.test(s)) return "xpl"; if (/fluor|uv/i.test(s)) return "fluorescence";
  if (/\bct\b|synchro|tomog/i.test(s)) return "ct"; if (/normal|plain|ppl|bright/i.test(s)) return "plain"; return null; }
function newDoc(t) {
  const task = {}; for (const k in t) if (!k.startsWith("_")) task[k] = t[k];
  const d = { schema_version: SCHEMA_VERSION, tool_version: TOOL_VERSION, task, identity_edits: {},
           image: { file: t.file, sha256: null, px_width: null, px_height: null },
           measurer: $("measurer").value.trim() || null, repeat: !!(t.assignment && t.assignment.repeat),
           started_at: new Date().toISOString(), time_spent_s: 0,
           setup: { view: guessView(t.view), light: guessLight(t.light), companion: null, marks_mode: "trace" },
           calibration: null, centre: null, centre_method: null, centre_fit: null, section: {}, flags: {}, curves: [], lacunae: [],
           counts: { visible_marks: null, visible_at_least: false, efs_lines_min: null }, field: null, scores: {}, notes: "", steps: {},
           features: { settings: {}, objects: [] } };
  if (t.companion_file) d.setup.companion = { file: t.companion_file, light: guessLight(t.companion_light), sha256: null, transform: null, transform_method: null };
  return d;
}
// 0.1 -> 0.2: double/multiple marks with a pair become branches; "section oblique" becomes the plane.
// 0.2 -> 0.3: adds `features` and the `ignore` curve kind.
// 0.3 -> 0.4: adds setup, counts, steps, field, centre_method; existing mark indices are kept as typed (index_auto false).
function migrate(d) {
  if (!d) return null;
  if (d.schema_version === SCHEMA_VERSION) return d;
  if (d.schema_version === "0.1") {
    for (const c of d.curves || []) if (c.kind === "mark") {
      if (c.mark_type === "double_multiple") { c.mark_type = "LAG"; c.branch_of = c.linked_to || null; }
      delete c.linked_to; if (!("branch_of" in c)) c.branch_of = null;
    }
    d.section = d.section || {};
    if (d.flags && d.flags.section_oblique) { d.section.plane = "oblique"; delete d.flags.section_oblique; }
    d.lacunae = d.lacunae || []; d.migrated_from = "0.1"; d.schema_version = "0.3";
  }
  if (d.schema_version === "0.2") d.schema_version = "0.3";
  if (d.schema_version !== "0.3") return null;
  d.features = d.features || { settings: {}, objects: [] };
  const t = d.task || {};
  d.setup = { view: guessView(t.view), light: guessLight(t.light), companion: null, marks_mode: "trace" };
  if (d.flags && d.flags.no_centre === undefined) d.flags.no_centre = false;
  d.centre_method = d.centre ? "clicked" : null; d.centre_fit = null;
  d.counts = { visible_marks: null, visible_at_least: false, efs_lines_min: null };
  d.field = null; d.steps = {};
  if (d.flags && d.flags.no_scale) d.steps.scale = { status: "skipped", reason: "no scale (0.3 flag)" };
  for (const c of d.curves || []) { if (c.kind === "mark") c.index_auto = false;
    if (c.zone && !c.zone.obs && ["matrix", "vasc_density", "vasc_orient", "remodelling"].some((f) => c.zone[f])) c.zone.legacy = true; }
  d.migrated_from = d.migrated_from || "0.3"; d.schema_version = SCHEMA_VERSION; d.tool_version = TOOL_VERSION;
  return d;
}
function saveDoc() {
  if (!S.doc) return;
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

// ---------- other-light image of the same section ----------
// Drawn over the main image through a similarity transform [a, b, tx, ty]: main = [a -b; b a] * companion + t.
async function loadCompanion() {
  if (S.comp && S.comp.url) URL.revokeObjectURL(S.comp.url);
  S.comp = null; const old = $("theImage2"); if (old) old.remove();
  const c = S.doc.setup.companion; $("padSwap").hidden = true;
  if (!c) return;
  const f = S.files.get(c.file); if (!f) { status(`Other-light image ${c.file} is not in this folder.`); return; }
  const url = URL.createObjectURL(f), dims = await imgDims(url).catch(() => null);
  if (!dims) { status(`${c.file}: cannot display this image.`); return; }
  S.comp = { url, w: dims[0], h: dims[1] };
  c.px_width = dims[0]; c.px_height = dims[1];
  if (!c.sha256) c.sha256 = await sha256(f);
  if (!c.transform) {
    const s = S.img.w / dims[0];
    c.transform = [s, 0, 0, 0]; c.transform_method = dims[0] === S.img.w && dims[1] === S.img.h ? "same_size" : "scaled_to_width";
  }
  el("image", { href: url, x: 0, y: 0, width: dims[0], height: dims[1], id: "theImage2" }, layerImg);
  store.set(keyComp(c.file), S.doc.task.figure_id);
  $("padSwap").hidden = false; applyCompanion();
}
function applyCompanion() {
  const im2 = $("theImage2"), im = $("theImage"); if (!im) return;
  if (!im2 || !S.doc.setup.companion) { im.style.opacity = 1; return; }
  const raw = S.align && S.align.pts.length % 2 === 1;               // clicking the landmark on the other image
  const [a, b, tx, ty] = S.doc.setup.companion.transform;
  im2.setAttribute("transform", raw ? "" : `matrix(${a} ${b} ${-b} ${a} ${tx} ${ty})`);
  im2.style.opacity = raw ? 1 : S.blend; im.style.opacity = raw ? 0 : 1;
  im2.style.filter = im.style.filter;
}
function swapLight() { if (!S.comp) return; S.blend = S.blend > 0.5 ? 0 : 1; applyCompanion(); renderPanels();
  status(S.blend ? `Showing ${lightName(S.doc.setup.companion.light)} (P switches back)` : `Showing ${lightName(S.doc.setup.light)}`); }
const lightName = (v) => (LIGHTS.find(([k]) => k === v) || [0, "the other image"])[1].toLowerCase();
function startAlign() { if (!S.comp) return; S.align = { pts: [] }; S.blend = 0; setTool("align"); applyCompanion(); }
function alignClick(p) {
  const A = S.align; A.pts.push(p);
  if (A.pts.length === 4) {
    const [p1, q1, p2, q2] = A.pts, dq = [q2[0] - q1[0], q2[1] - q1[1]], dp = [p2[0] - p1[0], p2[1] - p1[1]], n = dq[0] ** 2 + dq[1] ** 2;
    if (n < 4) { status("The two landmarks are too close together. Try again."); S.align = { pts: [] }; applyCompanion(); updateHint(); return; }
    const a = (dp[0] * dq[0] + dp[1] * dq[1]) / n, b = (dp[1] * dq[0] - dp[0] * dq[1]) / n;
    const c = S.doc.setup.companion; pushUndo();
    c.transform = [+a.toFixed(6), +b.toFixed(6), +(p1[0] - (a * q1[0] - b * q1[1])).toFixed(2), +(p1[1] - (b * q1[0] + a * q1[1])).toFixed(2)];
    c.transform_method = "two_landmarks"; c.landmarks = { main: [rnd(p1), rnd(p2)], companion: [rnd(q1), rnd(q2)] };
    S.align = null; S.blend = 0.5; setTool("pan"); applyCompanion(); renderPanels();
    status("Aligned. Check the overlap with the slider; P switches between the two.");
    return;
  }
  applyCompanion(); updateHint();
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
function onResize() { const [cw, ch] = svgSize(); const cx = S.view.x + S.view.w / 2, cy = S.view.y + S.view.h / 2;
  S.view.h = S.view.w * ch / cw; S.view.x = cx - S.view.w / 2; S.view.y = cy - S.view.h / 2; setView(); }
window.addEventListener("resize", onResize);
function toImg(ev) { const p = svg.createSVGPoint(); p.x = ev.clientX; p.y = ev.clientY;
  const q = p.matrixTransform(svg.getScreenCTM().inverse()); return [q.x, q.y]; }
function zoomAt(p, f) { const v = S.view; const nw = v.w * f; if (nw < 5 || nw > 200000) return;
  S.view = { x: p[0] - (p[0] - v.x) * f, y: p[1] - (p[1] - v.y) * f, w: nw, h: v.h * f }; setView(); }
const zoomMid = (f) => S.doc && zoomAt([S.view.x + S.view.w / 2, S.view.y + S.view.h / 2], f);
function applyFilter() { const f = `contrast(${S.filter.contrast}%) brightness(${S.filter.brightness}%)`;
  for (const id of ["theImage", "theImage2"]) { const im = $(id); if (im) im.style.filter = f; } }
$("contrast").addEventListener("input", (e) => { S.filter.contrast = +e.target.value; applyFilter(); });
$("brightness").addEventListener("input", (e) => { S.filter.brightness = +e.target.value; applyFilter(); });

// ---------- undo ----------
function snapshot() { const d = S.doc; return JSON.stringify({ c: d.calibration, z: d.centre, zm: d.centre_method, zf: d.centre_fit, v: d.curves, l: d.lacunae,
  t: d.setup.companion && d.setup.companion.transform }); }
function pushUndo() { if (!S.doc) return; S.undo.push(snapshot()); if (S.undo.length > 200) S.undo.shift(); touch(); }
function undo() { const u = S.undo.pop(); if (!u) return; const o = JSON.parse(u);
  Object.assign(S.doc, { calibration: o.c, centre: o.z, centre_method: o.zm, centre_fit: o.zf, curves: o.v, lacunae: o.l });
  if (S.doc.setup.companion && o.t) S.doc.setup.companion.transform = o.t;
  if (S.active && !curve(S.active.curveId)) S.active = null; if (S.activeT && !curve(S.activeT)) { S.activeT = null; setTool("pan"); }
  applyCompanion(); renderAll(); renderPanels(); }

// ---------- tools ----------
function setTool(t) {
  S.tool = t; S.scaleClick = null; S.tStart = null;
  if (t !== "trace") S.active = null;
  if (t !== "cross") S.activeT = null;
  if (t !== "align" && S.align) { S.align = null; applyCompanion(); }
  document.querySelectorAll("[data-tool]").forEach((b) => b.classList.toggle("active", b.dataset.tool === t));
  svg.classList.toggle("pan", t === "pan");
  updateHint(); renderAll();
}
function updateHint() {
  if (!S.doc) return hint("");
  if (S.tool === "scale") return hint(`Scale · point ${S.scaleClick ? 2 : 1} of 2`);
  if (S.tool === "centre") return hint("Centre");
  if (S.tool === "count") return hint("Lacunae");
  if (S.tool === "ignore") return hint("Ignore box");
  if (S.tool === "label") return hint("Label features");
  if (S.tool === "transect") return hint(`Transect · ${S.tStart ? "inner" : "outer"} end`);
  if (S.tool === "cross") return hint("Transect · crossings");
  if (S.tool === "align") { const n = S.align.pts.length; return hint(`Align · landmark ${n < 2 ? 1 : 2} · ${n % 2 ? "other-light image" : "main image"}`); }
  if (S.tool === "trace" && S.active) return hint(curveName(curve(S.active.curveId)));
  hint("");
}
document.querySelectorAll("[data-tool]").forEach((b) => b.addEventListener("click", () => S.doc && setTool(S.tool === b.dataset.tool ? "pan" : b.dataset.tool)));

const curve = (id) => S.doc.curves.find((c) => c.id === id);
function curveName(c) {
  if (!c) return "?";
  if (c.kind === "mark") { const p = c.branch_of && curve(c.branch_of);
    return p ? `Mark ${p.mark_index ?? "?"} branch` : `Mark ${c.mark_index ?? "?"}${c.index_uncertain ? "?" : ""}`; }
  if (c.kind === "label") return `Label ${c.label_index ?? "?"}`;
  if (c.kind === "transect") return `Transect ${S.doc.curves.filter((x) => x.kind === "transect").indexOf(c) + 1}`;
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
  if (c.kind === "mark") Object.assign(c, { mark_index: nextIndex("mark", "mark_index", c), index_auto: true, index_uncertain: false,
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
// Outline kinds: one curve per kind; tracing again adds a new piece (segment) to it.
function traceKind(kind) {
  const c = S.doc.curves.find((x) => x.kind === kind);
  if (c) { startTrace(c.id); status(`Adding a new piece to the ${KINDS[kind].label.toLowerCase()}`); } else newCurve(kind);
}
function startTrace(id) { closePop(); setTool("trace"); S.active = { curveId: id, seg: null }; S.current = id;
  svg.classList.remove("pan"); updateHint(); renderAll(); renderPanels(); svg.focus(); }
function endSegment() { if (S.active) { S.active.seg = null; updateHint(); renderAll(); } }
function stopTrace() { endSegment(); setTool("pan"); renderPanels(); }

// ---------- popover ----------
function closePop() { $("pop").hidden = true; S.popFor = null; }
function showPop(anchor, render, key) {
  if (!$("pop").hidden && S.popFor === key) return closePop();
  S.popFor = key; render(); const m = $("pop"); m.hidden = false;
  const r = anchor.getBoundingClientRect(), mw = m.offsetWidth, mh = m.offsetHeight;
  m.style.left = Math.max(8, Math.min(r.left - mw - 8 > 8 ? r.left - mw - 8 : r.left, window.innerWidth - mw - 8)) + "px";
  m.style.top = Math.max(8, Math.min(r.top, window.innerHeight - mh - 8)) + "px";
}
const seg = (opts, cur, onPick, fmt = nice) => { const d = h("div", { class: "seg" });
  for (const o of opts) { const [v, lab] = Array.isArray(o) ? o : [o, fmt(o)];
    d.append(h("button", { class: cur === v ? "on" : "", onclick: () => onPick(cur === v ? null : v) }, lab)); }
  return d; };
const row = (label, ...kids) => h("div", { class: "tm-row" }, h("label", {}, label), ...kids);

function tagTarget() { const id = S.active ? S.active.curveId : S.current; return id && S.doc ? curve(id) : null; }
function openTag(anchor) {
  const c = tagTarget(); if (!c) { status("Select a line in the list first."); return; }
  if (c.kind === "transect") return showPop(anchor, () => renderTransectPop(c), "tr:" + c.id);
  showPop(anchor, () => renderTag(c), "tag:" + c.id);
}
function renderTag(c) {
  const m = $("pop"); m.innerHTML = "";
  const re = () => { renderTag(c); renderAll(); renderPanels(); updateHint(); };
  m.append(h("div", { class: "tm-title" }, curveName(c)));
  if (TAGGABLE.includes(c.kind)) {
    const sel = h("select", { onchange: (e) => { setKind(c, e.target.value); re(); } });
    for (const k of TAGGABLE) sel.append(h("option", { value: k, selected: c.kind === k }, KINDS[k].label));
    m.append(row("This line is", sel));
  }
  if (c.kind === "mark") {
    const others = S.doc.curves.filter((x) => x.kind === "mark" && x.id !== c.id && !x.branch_of);
    if (!c.branch_of) m.append(row("Number",
      h("input", { type: "number", min: 1, step: 1, value: c.mark_index ?? "", disabled: c.index_auto, title: "1 = innermost mark you can see",
        onchange: (e) => { pushUndo(); c.mark_index = e.target.value === "" ? null : parseInt(e.target.value, 10); re(); } }),
      h("label", { class: "chk", title: "Numbered from the inside out by position" }, h("input", { type: "checkbox", checked: c.index_auto !== false,
        onchange: (e) => { c.index_auto = e.target.checked; re(); } }), "auto"),
      h("label", { class: "chk" }, h("input", { type: "checkbox", checked: c.index_uncertain,
        onchange: (e) => { c.index_uncertain = e.target.checked; re(); } }), "unsure")));
    m.append(row("Type", seg(MARK_TYPES, c.mark_type, (v) => { c.mark_type = v || "LAG"; re(); })));
    m.append(row("Confidence", seg([[3, "Sure"], [2, "Fairly"], [1, "Guess"]], c.confidence, (v) => { c.confidence = v || 3; re(); })));
    const sel = h("select", { onchange: (e) => { pushUndo(); c.branch_of = e.target.value || null;
      if (c.branch_of) { c.mark_index = curve(c.branch_of).mark_index; c.index_uncertain = false; } re(); } },
      h("option", { value: "" }, "no: its own mark"));
    for (const x of others) sel.append(h("option", { value: x.id, selected: c.branch_of === x.id }, curveName(x)));
    m.append(row("Splits off", sel));
    m.append(row("", h("label", { class: "chk" }, h("input", { type: "checkbox", checked: c.in_efs,
      onchange: (e) => { c.in_efs = e.target.checked; re(); } }), "inside the EFS")));
  }
  if (c.kind === "label") {
    m.append(row("Number", h("input", { type: "number", min: 1, step: 1, value: c.label_index ?? "",
      onchange: (e) => { c.label_index = e.target.value === "" ? null : parseInt(e.target.value, 10); re(); } })));
    m.append(row("Dye", seg(DYES, c.dye, (v) => { c.dye = v; re(); })));
    m.append(row("Date", h("input", { type: "date", value: c.label_date || "", onchange: (e) => { c.label_date = e.target.value; re(); } })));
    m.append(row("Age (days)", h("input", { type: "number", min: 0, step: 1, value: c.age_days ?? "", title: "Animal's age at injection",
      onchange: (e) => { c.age_days = e.target.value === "" ? null : +e.target.value; re(); } })));
  }
  m.append(h("div", { class: "tm-foot" }, h("button", { class: "primary", onclick: closePop }, "Done")));
}
function renderTransectPop(c) {
  const m = $("pop"); m.innerHTML = ""; const re = () => { renderTransectPop(c); renderAll(); renderPanels(); };
  m.append(h("div", { class: "tm-title" }, curveName(c)));
  for (const end of ["outer", "inner"]) {
    const sel = h("select", { onchange: (e) => { c.ends[end] = e.target.value; re(); } });
    for (const [v, lab] of ENDS[end]) sel.append(h("option", { value: v, selected: c.ends[end] === v }, lab));
    m.append(row(end === "outer" ? "Outer end at" : "Inner end at", sel));
  }
  const xs = crossings(c);
  if (!xs.length) m.append(h("p", { class: "muted" }, "No crossings yet. Click marks along the line."));
  for (const x of xs) {
    const sel = h("select", { onchange: (e) => { x.kind = e.target.value; re(); } });
    for (const [v, lab] of XING_KINDS) sel.append(h("option", { value: v, selected: x.kind === v }, lab));
    m.append(h("div", { class: "tm-row" }, h("label", {}, x.kind === "mark" ? `Mark ${x._idx ?? "?"}` : "·"),
      h("span", { class: "row" }, sel,
        x.kind === "mark" ? h("input", { type: "number", class: "num", min: 1, value: x.index_auto === false ? x.mark_index ?? "" : "", placeholder: "auto",
          title: "Leave empty to number automatically", onchange: (e) => { const v = e.target.value; x.index_auto = v === ""; x.mark_index = v === "" ? null : parseInt(v, 10); re(); } }) : null,
        h("button", { class: "icon", title: "Delete", onclick: () => { pushUndo(); c.crossings.splice(c.crossings.indexOf(x), 1); re(); } }, "✕"))));
  }
  m.append(h("div", { class: "tm-foot" },
    h("button", { onclick: () => { closePop(); S.activeT = c.id; S.tool = "cross"; updateHint(); renderAll(); } }, "Add crossings"),
    h("button", { class: "primary", onclick: closePop }, "Done")));
}
document.addEventListener("pointerdown", (e) => {
  const m = $("pop"); if (!m.hidden && !m.contains(e.target) && !e.target.closest("[data-pop]")) closePop();
  const hp = $("help"); if (!hp.hidden && !hp.contains(e.target) && e.target !== $("btnHelp")) hp.hidden = true;
});
$("btnHelp").addEventListener("click", () => { const hp = $("help"); hp.hidden = !hp.hidden; });

// ---------- geometry ----------
function resample(sg, step) { const out = [];
  for (let i = 0; i < sg.length - 1; i++) { const a = sg[i], b = sg[i + 1], L = dist(a, b), n = Math.max(1, Math.ceil(L / step));
    for (let j = 0; j < n; j++) out.push([a[0] + (b[0] - a[0]) * j / n, a[1] + (b[1] - a[1]) * j / n]); }
  if (sg.length) out.push(sg[sg.length - 1]); return out; }
function sample(c, max = 60) { const p = c.segments.flatMap((s) => resample(s, 2)); if (p.length <= max) return p;
  const k = p.length / max; return [...Array(max)].map((_, i) => p[Math.floor(i * k)]); }
function distToPolyline(p, segs) { let best = Infinity;
  for (const s of segs) for (let i = 0; i < s.length - 1; i++) { const a = s[i], b = s[i + 1];
    const dx = b[0] - a[0], dy = b[1] - a[1], L2 = dx * dx + dy * dy;
    const t = L2 ? Math.max(0, Math.min(1, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / L2)) : 0;
    best = Math.min(best, Math.hypot(p[0] - a[0] - t * dx, p[1] - a[1] - t * dy)); }
  return best; }
const angle = (q, z) => ((Math.atan2(q[1] - z[1], q[0] - z[0]) * 180 / Math.PI) + 360) % 360;
function coverageOf(pts, z) { const bins = new Array(72).fill(false); for (const q of pts) bins[Math.floor(angle(q, z) / 5) % 72] = true; return bins; }
// Angular coverage (5-degree bins) and the arc occupied, as [start, end] degrees clockwise in image coordinates.
function angular(c, z) {
  const pts = c.segments.flatMap((s) => resample(s, 2)); if (!z || pts.length < 2) return null;
  const bins = coverageOf(pts, z), n = bins.filter(Boolean).length; if (n === 72) return { coverage_deg: 360, arc_deg: [0, 360] };
  let gapStart = 0, gapLen = 0;                       // the largest empty run marks where the arc is not
  for (let i = 0; i < 72; i++) if (!bins[i]) { let L = 0; while (L < 72 && !bins[(i + L) % 72]) L++; if (L > gapLen) { gapLen = L; gapStart = i; } }
  return { coverage_deg: n * 5, arc_deg: [((gapStart + gapLen) % 72) * 5, (gapStart * 5) % 360 || 360] };
}
function polyArea(pts) { let a = 0; for (let i = 0; i < pts.length; i++) { const p = pts[i], q = pts[(i + 1) % pts.length]; a += p[0] * q[1] - q[0] * p[1]; } return Math.abs(a) / 2; }
function inPoly(p, pts) { let inside = false;
  for (let i = 0, j = pts.length - 1; i < pts.length; j = i++) { const a = pts[i], b = pts[j];
    if ((a[1] > p[1]) !== (b[1] > p[1]) && p[0] < (b[0] - a[0]) * (p[1] - a[1]) / (b[1] - a[1]) + a[0]) inside = !inside; }
  return inside; }

// Boundaries in radial order, innermost first. The medullary margin is always innermost and the outer edge
// (periosteum, broken edge) always outermost; the rest are ordered by distance from the centre if set, otherwise
// by distance below the outer edge (wedges with no centre), otherwise by kind and index.
function orderedBoundaries(doc = S.doc) {
  const bs = doc.curves.filter((c) => KINDS[c.kind].boundary && !c.branch_of && hasPts(c));
  const z = doc.centre, outerSegs = doc.curves.filter((c) => OUTER.has(c.kind)).flatMap((c) => c.segments).filter((s) => s.length > 1);
  const fixed = (c) => (c.kind === "medullary" ? -Infinity : OUTER.has(c.kind) ? Infinity : null);
  let key, method;
  if (z) { method = "centre"; key = (c) => { const p = sample(c); return p.reduce((a, q) => a + dist(q, z), 0) / p.length; }; }
  else if (outerSegs.length) { method = "periosteum_distance"; key = (c) => { const p = sample(c);
    return -p.reduce((a, q) => a + distToPolyline(q, outerSegs), 0) / p.length; }; }
  else { method = "index"; const rank = { remodel_front: 0, mark: 1, label: 1, efs_inner: 2 };
    key = (c) => rank[c.kind] * 1000 + (c.mark_index || c.label_index || 0); }
  const k = new Map(bs.map((c) => [c.id, fixed(c) ?? key(c)]));
  return { method, list: bs.sort((a, b) => k.get(a.id) - k.get(b.id) || (a.kind === "periosteum" ? -1 : b.kind === "periosteum" ? 1 : 0)) };
}
// Zones: the tissue just inside each boundary. Descriptors are stored on the zone's outer boundary (`zone`). The outer
// edge counts once: its zone is stored on the periosteum if traced, otherwise on the broken edge.
function zones(doc = S.doc) {
  const { method, list } = orderedBoundaries(doc);
  const outers = list.filter((c) => OUTER.has(c.kind)), rep = outers[0];
  const L = list.filter((c) => !OUTER.has(c.kind)); if (rep) L.push(rep);
  const out = [];
  L.forEach((c, i) => { if (c.kind === "medullary") return; out.push({ inner: L[i - 1] || null, outer: c }); });
  return { method, zones: out };
}
// Marks numbered by radial order (1 = innermost), except those the measurer numbered by hand. Branches take their
// main mark's number.
function autoIndex(d = S.doc) {
  const { method, list } = orderedBoundaries(d);
  if (method !== "index") list.filter((c) => c.kind === "mark").forEach((c, i) => { if (c.index_auto !== false) c.mark_index = i + 1; });
  for (const c of d.curves) if (c.kind === "mark" && c.branch_of) { const p = d.curves.find((x) => x.id === c.branch_of); if (p) c.mark_index = p.mark_index; }
}

// Centre from the outer edge: centroid of the outline if it closes (>= 330 degrees), else a least-squares circle
// fitted to the arc (refused under 90 degrees).
function circleFit(pts) {
  const n = pts.length, mx = pts.reduce((a, p) => a + p[0], 0) / n, my = pts.reduce((a, p) => a + p[1], 0) / n;
  let suu = 0, svv = 0, suv = 0, suuu = 0, svvv = 0, suvv = 0, svuu = 0;
  for (const p of pts) { const u = p[0] - mx, v = p[1] - my; suu += u * u; svv += v * v; suv += u * v; suuu += u * u * u; svvv += v * v * v; suvv += u * v * v; svuu += v * u * u; }
  const det = suu * svv - suv * suv; if (Math.abs(det) < 1e-9) return null;
  const b1 = (suuu + suvv) / 2, b2 = (svvv + svuu) / 2, uc = (b1 * svv - b2 * suv) / det, vc = (suu * b2 - suv * b1) / det;
  const c = [mx + uc, my + vc], r = Math.sqrt(uc * uc + vc * vc + (suu + svv) / n);
  return { c, r, rms: Math.sqrt(pts.reduce((a, p) => a + (dist(p, c) - r) ** 2, 0) / n) };
}
// Uses the periosteum alone when it spans 90 degrees or more: broken edges often run radially (cracks) and would bias it.
function fitCentre(d = S.doc) {
  const ptsOf = (kinds) => d.curves.filter((c) => kinds.includes(c.kind)).flatMap((c) => c.segments.flatMap((s) => resample(s, 3)));
  let pts = ptsOf(["periosteum"]), cf = pts.length >= 8 ? circleFit(pts) : null;
  if (!cf || coverageOf(pts, cf.c).filter(Boolean).length * 5 < 90) { pts = ptsOf(["periosteum", "broken_edge"]); cf = pts.length >= 8 ? circleFit(pts) : null; }
  if (pts.length < 8) return { err: "Trace the outer edge (periosteum or broken edge) first." };
  if (!cf) return { err: "Could not fit a circle to the outer edge." };
  const cov = coverageOf(pts, cf.c).filter(Boolean).length * 5;
  if (cov < 90) return { err: `The outer edge covers only ${cov}° around its centre: too little to fit one. Click the centre, or choose No centre.` };
  if (cov >= 330) {
    const ring = [...pts].sort((a, b) => angle(a, cf.c) - angle(b, cf.c)); let A = 0, cx = 0, cy = 0;
    for (let i = 0; i < ring.length; i++) { const p = ring[i], q = ring[(i + 1) % ring.length], x = p[0] * q[1] - q[0] * p[1]; A += x; cx += (p[0] + q[0]) * x; cy += (p[1] + q[1]) * x; }
    if (Math.abs(A) > 1) return { c: [cx / (3 * A), cy / (3 * A)], fit: { method: "outline_centroid", arc_deg: cov, rms_px: +cf.rms.toFixed(2) } };
  }
  return { c: cf.c, fit: { method: "circle_fit", arc_deg: cov, rms_px: +cf.rms.toFixed(2) } };
}

// Transect crossings, outer end first, with mark numbers (1 = innermost mark crossing on this transect).
function crossings(c) {
  const xs = c.crossings || []; xs.sort((a, b) => a.t - b.t);
  const marks = xs.filter((x) => x.kind === "mark"); let n = 0;
  for (let i = marks.length - 1; i >= 0; i--) { n++; marks[i]._idx = marks[i].index_auto === false ? marks[i].mark_index : n; }
  return xs;
}
function project(p, a, b) { const dx = b[0] - a[0], dy = b[1] - a[1], L2 = dx * dx + dy * dy;
  const t = L2 ? Math.max(0, Math.min(1, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / L2)) : 0; return { t, p: [a[0] + t * dx, a[1] + t * dy] }; }
function nearestKind(p, kinds) { const tol = 15 * upp(); let best = null, bd = tol;
  for (const c of S.doc.curves) if (kinds.includes(c.kind)) { const dd = distToPolyline(p, c.segments); if (dd < bd) { bd = dd; best = c.kind; } }
  return best; }

// ---------- pan pad ----------
function panBy(fx, fy) { S.view.x += S.view.w * fx; S.view.y += S.view.h * fy; setView(); }
$("padIn").addEventListener("click", () => zoomMid(1 / 1.25));
$("padOut").addEventListener("click", () => zoomMid(1.25));
$("padFit").addEventListener("click", fit);
$("padSwap").addEventListener("click", swapLight);

function addVertex(p) {
  const c = curve(S.active.curveId); pushUndo();
  if (S.active.seg == null) { c.segments.push([]); S.active.seg = c.segments.length - 1; }
  c.segments[S.active.seg].push(rnd(p));
  renderAll(); renderPanels();
}
function removeLastVertex() {
  if (S.tool === "cross" && S.activeT) { const c = curve(S.activeT); if (c && c.crossings.length) { pushUndo(); c.crossings.pop(); renderAll(); renderPanels(); } return; }
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
  if (tgt.classList.contains("lac") && ev.button === 2) { pushUndo(); S.doc.lacunae.splice(+tgt.dataset.i, 1); renderAll(); renderPanels(); return; }
  if (tgt.classList.contains("xing") && ev.button === 2) { const c = curve(tgt.dataset.c); pushUndo();
    c.crossings.splice(+tgt.dataset.i, 1); renderAll(); renderPanels(); return; }
  if (tgt.classList.contains("vert")) {
    const c = curve(tgt.dataset.c), si = +tgt.dataset.s, vi = +tgt.dataset.v;
    if (ev.button === 2) { pushUndo(); c.segments[si].splice(vi, 1); if (!c.segments[si].length) c.segments.splice(si, 1);
      if (S.active && S.active.curveId === c.id) S.active.seg = null; renderAll(); renderPanels(); return; }
    if (S.tool === "pan" || (S.active && S.active.curveId === c.id && ev.altKey)) {
      pushUndo(); S.drag = { type: "vert", c, si, vi }; svg.setPointerCapture(ev.pointerId); return; }
  }
  if (tgt.classList.contains("line") && S.tool === "pan" && ev.button === 0 && !S.space) { S.current = tgt.dataset.c; renderAll(); renderPanels(); return; }
  if (tgt.dataset && tgt.dataset.handle && S.tool === "pan") { pushUndo(); S.drag = { type: "handle", which: tgt.dataset.handle }; svg.setPointerCapture(ev.pointerId); return; }
  if (ev.button === 1 || ev.button === 2 || S.space || S.tool === "pan") {
    S.drag = { type: "pan", sx: ev.clientX, sy: ev.clientY, v: Object.assign({}, S.view) };
    svg.classList.add("panning"); svg.setPointerCapture(ev.pointerId); return;
  }
  if (S.tool === "scale") {
    if (!S.scaleClick) { S.scaleClick = p; updateHint(); renderAll(); }
    else { pushUndo(); const prev = S.doc.calibration || {}, guess = parseScale(S.doc.task.scale_value);
      S.doc.calibration = { p1: rnd(S.scaleClick), p2: rnd(p), length: prev.length || guess.length || null, unit: prev.unit || guess.unit || "mm" };
      S.scaleClick = null; delete S.doc.steps.scale; S.doc.flags.no_scale = false; setTool("pan"); renderPanels();
      const inp = $("calLen"); if (inp && !S.doc.calibration.length) inp.focus(); }
    return;
  }
  if (S.tool === "centre") { pushUndo(); S.doc.centre = rnd(p); S.doc.centre_method = "clicked"; S.doc.centre_fit = null; S.doc.flags.no_centre = false;
    setTool("pan"); renderPanels(); return; }
  if (S.tool === "count") { pushUndo(); S.doc.lacunae.push(rnd(p)); renderAll(); renderPanels(); return; }
  if (S.tool === "ignore") { S.drag = { type: "box", a: p, b: p }; svg.setPointerCapture(ev.pointerId); return; }
  if (S.tool === "label") { for (const f of hooks.label) f(p, ev); return; }
  if (S.tool === "align") { alignClick(rnd(p)); return; }
  if (S.tool === "transect") {
    if (!S.tStart) { S.tStart = rnd(p); updateHint(); renderAll(); return; }
    pushUndo(); const a = S.tStart, b = rnd(p);
    const c = { id: uid(), kind: "transect", segments: [[a, b]], crossings: [],
      ends: { outer: nearestKind(a, ["periosteum", "broken_edge"]) || "unclear",
              inner: nearestKind(b, ["medullary", "remodel_front"]) || (S.doc.centre && dist(b, S.doc.centre) < 20 * upp() ? "centre" : "unclear") } };
    S.doc.curves.push(c); S.tStart = null; S.tool = "cross"; S.activeT = c.id; S.current = c.id; updateHint(); renderAll(); renderPanels(); return;
  }
  if (S.tool === "cross" && S.activeT) {
    const c = curve(S.activeT), [a, b] = c.segments[0], pr = project(p, a, b); pushUndo();
    c.crossings.push({ p: rnd(pr.p), t: +pr.t.toFixed(5), kind: "mark", index_auto: true, mark_index: null }); renderAll(); renderPanels(); return;
  }
  if (S.tool === "trace" && S.active) addVertex(p);
});
svg.addEventListener("pointermove", (ev) => {
  if (!S.doc) return;
  if (S.tool === "transect" && S.tStart) { S.rubber = toImg(ev); renderAll(); }
  if (!S.drag) return; touch(); const p = toImg(ev);
  if (S.drag.type === "pan") { const f = upp(); const d = S.drag;
    S.view.x = d.v.x - (ev.clientX - d.sx) * f; S.view.y = d.v.y - (ev.clientY - d.sy) * f; setView(); }
  else if (S.drag.type === "vert") { S.drag.c.segments[S.drag.si][S.drag.vi] = rnd(p); renderAll(); }
  else if (S.drag.type === "box") { S.drag.b = p; renderAll(); }
  else if (S.drag.type === "handle") { const w = S.drag.which;
    if (w === "centre") { S.doc.centre = rnd(p); S.doc.centre_method = "clicked"; S.doc.centre_fit = null; } else S.doc.calibration[w] = rnd(p); renderAll(); }
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
  else if (ev.key === "?") $("help").hidden = !$("help").hidden;
  else if (!S.doc) return;
  else if (ev.key === "t" || ev.key === "T") openTag($("cards"));
  else if (ev.key === "n" || ev.key === "N") { goStep("marks"); newCurve("mark"); }
  else if (ev.key === "p" || ev.key === "P") swapLight();
  else if (ev.key === "[") stepBy(-1);
  else if (ev.key === "]") stepBy(1);
  else if (ev.key.startsWith("Arrow")) { ev.preventDefault(); const st = ev.shiftKey ? 0.5 : 0.1;
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
  const label = (q, txt, color) => { const t = el("text", { x: q[0] + 6 * k, y: q[1] - 6 * k, fill: color, "font-size": 13 * k, "font-weight": 700,
    stroke: "#000", "stroke-width": 3 * k, "paint-order": "stroke" }, layerShapes); t.textContent = txt; };
  for (const c of S.doc.curves) {
    const isAct = (S.active && S.active.curveId === c.id) || S.activeT === c.id;
    const isSel = S.hover === c.id || S.current === c.id || hiZone.includes(c.id);
    const color = KINDS[c.kind].color, closed = KINDS[c.kind].closed;
    if (c.kind === "transect") { drawTransect(c, k, isAct || isSel); continue; }
    c.segments.forEach((sg, si) => {
      if (!sg.length) return;
      const pts = sg.map((q) => q.join(",")).join(" ");
      el(closed ? "polygon" : "polyline", { points: pts, fill: c.kind === "ignore" ? "rgba(0,0,0,.45)" : closed ? "rgba(148,163,184,.12)" : "none", stroke: color,
        "stroke-width": (isAct || isSel ? 2.2 : 1.4) * sw, "stroke-linejoin": "round", "stroke-linecap": "round",
        "stroke-dasharray": (c.kind === "mark" && (c.mark_type === "uncertain" || c.branch_of)) || closed || c.kind === "broken_edge" ? `${6 * k} ${4 * k}` : "",
        opacity: 0.95, "pointer-events": "none" }, layerShapes);
      if (S.tool === "pan") el(closed ? "polygon" : "polyline", { points: pts, fill: "none", stroke: "transparent", "stroke-width": 12 * k,
        class: "line", "data-c": c.id, style: "cursor:pointer" }, layerShapes);         // wide invisible hit area: click to select
      if (isAct || isSel) sg.forEach((q, vi) => {
        const last = isAct && S.active && S.active.seg === si && vi === sg.length - 1;
        el("circle", { cx: q[0], cy: q[1], r: last ? r * 1.4 : r, fill: last ? "var(--sel)" : color, stroke: "#000", "stroke-width": 0.8 * k,
          class: "vert", "data-c": c.id, "data-s": si, "data-v": vi }, layerVerts); });
    });
    if ((c.kind === "mark" || c.kind === "label") && c.segments.length && c.segments[0].length) label(c.segments[0][0], shortName(c), color);
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
  if (S.tStart) { el("circle", { cx: S.tStart[0], cy: S.tStart[1], r: r * 1.2, fill: "var(--c-transect)" }, layerVerts);
    if (S.rubber) el("line", { x1: S.tStart[0], y1: S.tStart[1], x2: S.rubber[0], y2: S.rubber[1], stroke: "var(--c-transect)", "stroke-width": 1.5 * k,
      "stroke-dasharray": `${4 * k} ${3 * k}`, "pointer-events": "none" }, layerShapes); }
  if (S.align) S.align.pts.forEach((q, i) => { const main = i % 2 === 0, show = (S.align.pts.length % 2 === 1) !== main;
    if (show) el("circle", { cx: q[0], cy: q[1], r: r * 1.6, fill: "none", stroke: "var(--c-transect)", "stroke-width": 2 * k }, layerVerts); });
  const z = S.doc.centre;
  if (z) { const L = 14 * k;
    el("path", { d: `M${z[0] - L},${z[1]}H${z[0] + L}M${z[0]},${z[1] - L}V${z[1] + L}`, stroke: "var(--c-centre)", "stroke-width": 2 * k,
      "stroke-dasharray": S.doc.centre_method === "fitted" ? `${3 * k} ${2 * k}` : "" }, layerShapes);
    el("circle", { cx: z[0], cy: z[1], r: r * 1.2, fill: "transparent", stroke: "var(--c-centre)", "stroke-width": 1.5 * k,
      "data-handle": "centre", class: "handle" }, layerVerts); }
  for (const f of hooks.render) f();
}
function drawTransect(c, k, hi) {
  const [a, b] = c.segments[0], color = KINDS.transect.color, L = dist(a, b) || 1, nx = -(b[1] - a[1]) / L, ny = (b[0] - a[0]) / L;
  el("line", { x1: a[0], y1: a[1], x2: b[0], y2: b[1], stroke: color, "stroke-width": (hi ? 2.4 : 1.6) * k, opacity: 0.9 }, layerShapes);
  for (const q of [a, b]) el("rect", { x: q[0] - 3 * k, y: q[1] - 3 * k, width: 6 * k, height: 6 * k, fill: color, stroke: "#000", "stroke-width": 0.6 * k }, layerShapes);
  crossings(c).forEach((x, i) => {
    const s = 9 * k, col = x.kind === "mark" ? "var(--c-mark)" : KINDS[x.kind].color;
    el("line", { x1: x.p[0] - nx * s, y1: x.p[1] - ny * s, x2: x.p[0] + nx * s, y2: x.p[1] + ny * s, stroke: col, "stroke-width": 3 * k,
      class: "xing", "data-c": c.id, "data-i": c.crossings.indexOf(x), style: "cursor:pointer" }, layerVerts);
    if (x.kind === "mark") { const t = el("text", { x: x.p[0] + nx * (s + 4 * k), y: x.p[1] + ny * (s + 4 * k) + 4 * k, fill: col, "font-size": 12 * k,
      "font-weight": 700, stroke: "#000", "stroke-width": 3 * k, "paint-order": "stroke", "pointer-events": "none" }, layerShapes); t.textContent = x._idx ?? "?"; }
  });
}

// ---------- steps ----------
function applies(id, d = S.doc) {
  const v = d.setup.view;
  if (v === "unusable") return id === "image" || id === "export";
  if (v === "detail" && (id === "outline" || id === "marks")) return false;
  return true;
}
const skipped = (id, d = S.doc) => d.steps[id] && d.steps[id].status === "skipped";
function stepState(id, d = S.doc) {
  if (!applies(id, d)) return "na";
  if (skipped(id, d)) return "skipped";
  const done = {
    image: () => d.setup.view && d.setup.light,
    scale: () => pxPerUnit(d),
    outline: () => d.curves.some((c) => OUTER.has(c.kind) && hasPts(c)) && (d.centre || d.flags.no_centre),
    marks: () => d.curves.some((c) => c.kind === "mark" && hasPts(c)) || d.curves.some((c) => c.kind === "transect" && c.crossings.length) || d.counts.visible_marks != null,
    tissue: () => { const zs = zones(d).zones; const ctx = { polarized: isPolarized(d) };
      return zs.length ? zs.every((z) => zoneDone(z.outer.zone, ctx)) : zoneDone(d.field, ctx); },
    bone: () => d.scores.efs != null,
    more: () => false,
    export: () => store.get(keyStatus(d.task.figure_id)) === "done",
  }[id]();
  return done ? "done" : "todo";
}
function firstOpenStep() { const s = STEPS.find((x) => x.id !== "more" && stepState(x.id) === "todo"); return s ? s.id : "export"; }
function goStep(id) { if (S.step === id) return; S.step = id; closePop(); if (["trace", "transect", "cross", "scale", "centre"].includes(S.tool)) setTool("pan");
  renderPanels(); $("cards").scrollTop = 0; }
function stepBy(n) { const ids = STEPS.map((s) => s.id).filter((id) => applies(id)); let i = ids.indexOf(S.step) + n;
  if (i >= 0 && i < ids.length) goStep(ids[i]); }
$("stepBack").addEventListener("click", () => stepBy(-1));
$("stepNext").addEventListener("click", () => stepBy(1));

function renderRail() {
  const nav = $("rail"); nav.innerHTML = "";
  STEPS.forEach((s, i) => { const st = stepState(s.id);
    nav.append(h("button", { class: [S.step === s.id ? "on" : "", st].join(" "), title: s.title + (st === "skipped" ? ` (skipped: ${S.doc.steps[s.id].reason})` : st === "na" ? " (not needed for this image)" : ""),
      onclick: () => goStep(s.id) }, h("span", { class: "pip" }, st === "done" ? "✓" : st === "skipped" ? "–" : String(i + 1)), s.rail)); });
  const ids = STEPS.map((s) => s.id).filter((id) => applies(id)), i = ids.indexOf(S.step);
  $("stepBack").disabled = i <= 0; $("stepNext").disabled = i >= ids.length - 1;
  $("stepCount").textContent = `${i + 1} of ${ids.length}`;
}
function renderIdent() {
  const t = S.doc.task, ie = S.doc.identity_edits, g = (k) => t[k] || ie[k];
  $("ident").innerHTML = "";
  $("ident").append(h("b", {}, g("taxon") || "Unknown taxon"), " ", h("span", {}, [g("specimen"), g("element")].filter(Boolean).join(" · ")),
    h("span", { class: "muted" }, g("citation") || t.figure_id));
  const n = S.tasks.length, done = S.tasks.filter((x) => store.get(keyStatus(x.figure_id)) === "done").length;
  $("taskName").innerHTML = ""; $("taskName").append(t.figure_id, h("span", { class: "muted" }, `  ${S.idx + 1}/${n}` + (done ? ` · ${done} done` : "")));
}
function renderPanels() {
  if (!S.doc) return;
  autoIndex();
  renderRail(); renderIdent();
  document.querySelectorAll(".card").forEach((c) => c.classList.toggle("on", c.dataset.step === S.step));
  const body = $("body-" + S.step); body.innerHTML = "";
  const s = STEPS.find((x) => x.id === S.step);
  body.append(h("h3", {}, s.title));
  if (!applies(S.step)) body.append(h("div", { class: "box" }, S.doc.setup.view === "unusable" ? "Not needed: the image was marked as unusable." : "Not needed for a close-up."));
  else STEP_BODY[S.step](body);
  if (s.reasons && applies(S.step)) body.append(skipBar(s));
  updateHint();
}
function skipBar(s) {
  const st = S.doc.steps[s.id];
  if (st && st.status === "skipped") return h("div", { class: "skipbar" }, `Skipped: ${st.reason}. `,
    h("button", { class: "link", onclick: () => { delete S.doc.steps[s.id]; if (s.id === "scale") S.doc.flags.no_scale = false; touch(); renderPanels(); } }, "Undo"));
  const sel = h("select", { onchange: (e) => { if (!e.target.value) return; S.doc.steps[s.id] = { status: "skipped", reason: e.target.value };
    if (s.id === "scale") S.doc.flags.no_scale = true;
    if (s.id === "outline" && !S.doc.centre) S.doc.flags.no_centre = true;
    touch(); stepBy(1); renderPanels(); } }, h("option", { value: "" }, "reason…"));
  for (const r of s.reasons) sel.append(h("option", { value: r }, r));
  return h("div", { class: "skipbar" }, "Not done: ", sel);
}

// One answer from a list of buttons, with optional pictures and a "Can't tell" button.
function choices(list, cur, onPick, { unclear = true, unclearLabel = "Can't tell" } = {}) {
  const pics = list.some((o) => o[2]);
  const d = h("div", { class: "opts" + (pics ? " pics" : "") });
  for (const [v, label, icon] of list) { const b = h("button", { class: cur === v ? "on" : "", onclick: () => onPick(cur === v ? null : v) });
    if (icon && HG.ICONS[icon]) b.innerHTML = HG.ICONS[icon]; b.append(h("span", {}, label)); d.append(b); }
  if (unclear) d.append(h("button", { class: "unclear" + (cur === HG.UNCLEAR ? " on" : ""), onclick: () => onPick(cur === HG.UNCLEAR ? null : HG.UNCLEAR) }, unclearLabel));
  return d;
}
const question = (text, tip, cur, ...kids) => h("div", { class: "q" + (cur != null ? " done" : "") }, h("div", { class: "qt" }, text), tip ? h("div", { class: "tip" }, tip) : null, ...kids);
function kindButton(kind, onclick, label) { const has = S.doc.curves.some((c) => c.kind === kind && hasPts(c));
  return h("button", { class: "btn-kind" + (S.active && curve(S.active.curveId) && curve(S.active.curveId).kind === kind ? " active" : ""), 
    onclick: onclick || (() => traceKind(kind)) }, h("span", { class: "sw", style: `background:${KINDS[kind].color}` }), label || (has ? `${KINDS[kind].label} +` : KINDS[kind].label)); }

const STEP_BODY = {
  image(b) {
    const d = S.doc, su = d.setup, re = () => { touch(); renderPanels(); };
    if (d.task.caption) b.append(h("details", { class: "fold" }, h("summary", {}, "Figure caption"), h("div", { class: "small" }, d.task.caption)));
    b.append(question("What does this image show?", null, su.view,
      choices(VIEWS, su.view, (v) => { su.view = v; if (v !== "unusable") delete d.steps.image; re(); }, { unclear: false })));
    if (su.view === "unusable") b.append(h("div", { class: "row" }, "Why? ",
      h("input", { value: (d.steps.image && d.steps.image.reason) || "", placeholder: "e.g. drawing, no bone, too small",
        onchange: (e) => { d.steps.image = { status: "skipped", reason: e.target.value || "unusable" }; touch(); } })));
    b.append(question("How was it lit?", null,
      su.light, choices(LIGHTS, su.light, (v) => { su.light = v; re(); }, { unclear: false })));
    // other-light image
    const others = S.tasks.filter((t) => t.file !== d.image.file && !t._missing);
    const comp = su.companion;
    const blk = h("div", { class: "q" + (comp ? " done" : "") }, h("div", { class: "qt" }, "Same section, other light"));
    if (!comp) {
      const sel = h("select", { onchange: async (e) => { if (!e.target.value) return; pushUndo();
        su.companion = { file: e.target.value, light: su.light === "plain" ? "xpl" : su.light && POLARIZED.has(su.light) ? "plain" : null, sha256: null, transform: null, transform_method: null };
        await loadCompanion(); S.blend = 0; renderTasks(); re(); } }, h("option", { value: "" }, others.length ? "choose an image…" : "no other images in this folder"));
      for (const t of others) sel.append(h("option", { value: t.file }, t.figure_id));
      blk.append(sel);
    } else {
      const aligned = comp.transform_method === "two_landmarks" || comp.transform_method === "same_size";
      blk.append(h("div", { class: "row" }, h("b", {}, comp.file), h("button", { class: "link", onclick: () => { pushUndo(); store.del(keyComp(comp.file)); su.companion = null; loadCompanion(); renderTasks(); re(); } }, "remove")),
        h("div", { class: "row" }, "Its light: ", seg(LIGHTS.slice(0, 4), comp.light, (v) => { comp.light = v; re(); })),
        h("div", { class: "row" }, h("button", { class: S.tool === "align" ? "active" : "", onclick: startAlign }, "Align (2 landmarks)"),
          h("span", { class: aligned ? "ok small" : "warn small" }, comp.transform_method === "two_landmarks" ? "aligned" : comp.transform_method === "same_size" ? "same size: assumed aligned" : "not aligned yet")),
        h("div", { class: "row" }, "Show: ", h("input", { type: "range", min: 0, max: 100, value: Math.round(S.blend * 100), style: "flex:1",
          oninput: (e) => { S.blend = e.target.value / 100; applyCompanion(); } }), h("kbd", {}, "P")));
    }
    b.append(blk);
    const meta = h("dl", { class: "meta" }), t = d.task; let missing = 0;
    for (const k of IDENTITY) { const v = t[k]; meta.append(h("dt", {}, nice(k)));
      if (v != null && v !== "") meta.append(h("dd", {}, String(v)));
      else { if (["taxon", "specimen", "element"].includes(k) && !d.identity_edits[k]) missing++;
        meta.append(h("dd", {}, h("input", { value: d.identity_edits[k] || "", placeholder: "—", onchange: (e) => { d.identity_edits[k] = e.target.value; touch(); renderIdent(); } }))); } }
    b.append(h("details", { class: "fold", open: missing > 0 }, h("summary", {}, missing ? `Specimen details (${missing} missing)` : "Specimen details"), meta));
  },
  scale(b) {
    const d = S.doc, c = d.calibration;
    b.append(h("div", { class: "tools" }, h("button", { class: S.tool === "scale" ? "active primary" : "primary", onclick: () => setTool(S.tool === "scale" ? "pan" : "scale") },
      c ? "Reset scale bar" : "Scale bar")));
    if (c) {
      const un = h("select", { onchange: (e) => { c.unit = e.target.value; touch(); renderPanels(); } });
      for (const u of ["um", "mm", "cm"]) un.append(h("option", { value: u, selected: c.unit === u }, u === "um" ? "µm" : u));
      b.append(h("div", { class: "row" }, "Length ", h("input", { id: "calLen", type: "number", step: "any", min: "0", value: c.length ?? "", class: "num",
        onchange: (e) => { c.length = parseFloat(e.target.value) || null; touch(); renderPanels(); } }), un));
      if (d.task.scale_value) b.append(h("div", { class: "muted small" }, `The figure says: ${d.task.scale_value}`));
      const k = pxPerUnit(); if (k) b.append(h("div", { class: "box ok" }, `${(k).toFixed(2)} px per ${c.unit === "um" ? "µm" : c.unit}`));
    }
  },
  outline(b) {
    const d = S.doc, re = () => { touch(); renderPanels(); };
    b.append(h("div", { class: "sub" }, "Outer edge"),
      h("div", { class: "tools" }, kindButton("periosteum"), kindButton("broken_edge")));
    b.append(h("div", { class: "sub" }, "Inner edge"),
      h("div", { class: "tools" }, kindButton("medullary"), kindButton("remodel_front")));
    b.append(curveList(["periosteum", "broken_edge", "medullary", "remodel_front"]));
    b.append(h("div", { class: "sub" }, "Centre"));
    const fitted = d.centre_method === "fitted", f = d.centre_fit;
    b.append(h("div", { class: "tools" },
      h("button", { class: d.centre_method === "clicked" ? "active" : "", onclick: () => setTool("centre") }, "Click it"),
      h("button", { class: fitted ? "active" : "", onclick: () => { const r = fitCentre(); if (r.err) { status(r.err, true); return; }
        pushUndo(); d.centre = rnd(r.c); d.centre_method = "fitted"; d.centre_fit = r.fit; d.flags.no_centre = false; re(); } }, "Fit from outline"),
      h("button", { class: d.flags.no_centre ? "active" : "", onclick: () => { pushUndo(); d.centre = null; d.centre_method = null; d.centre_fit = null; d.flags.no_centre = true; re(); } }, "No centre")));
    b.append(h("div", { class: "small muted" }, d.centre ? (fitted ? `Fitted: ${f.method === "outline_centroid" ? "outline centroid" : "circle"}, ${f.arc_deg}°` : "Clicked") : d.flags.no_centre ? "None" : "Not set"));
    b.append(h("div", { class: "sub" }, "Section"));
    const g = h("div", { class: "grid2" });
    for (const [k, lab, opts] of SECTION_FIELDS) {
      const sel = h("select", { onchange: (e) => { d.section[k] = e.target.value || null; re(); } }, h("option", { value: "" }, "—"));
      for (const o of opts) sel.append(h("option", { value: o, selected: d.section[k] === o }, nice(o)));
      g.append(h("label", {}, lab), sel);
    }
    b.append(g, h("div", { class: "chips" }, ...FLAGS.map(([k, lab]) => h("button", { class: "chip" + (d.flags[k] ? " on" : ""), onclick: () => { d.flags[k] = !d.flags[k]; re(); } }, lab))));
  },
  marks(b) {
    const d = S.doc, su = d.setup, re = () => { touch(); renderPanels(); };
    b.append(h("div", { class: "row" }, seg([["trace", "Trace"], ["transect", "Transect"]], su.marks_mode || "trace", (v) => { su.marks_mode = v || "trace"; re(); })));
    if ((su.marks_mode || "trace") === "trace") {
      b.append(h("div", { class: "tools" }, kindButton("mark", () => newCurve("mark"), "Growth mark"), h("kbd", {}, "N"),
          su.light === "fluorescence" || (su.companion && su.companion.light === "fluorescence") ? kindButton("label", () => newCurve("label"), "Fluorochrome label") : null,
          kindButton("efs_inner", () => traceKind("efs_inner"), "EFS inner boundary")));
      b.append(curveList(["mark", "label", "efs_inner"]));
    } else {
      b.append(h("div", { class: "tools" }, h("button", { class: S.tool === "transect" ? "active primary" : "primary", onclick: () => setTool("transect") }, "Add transect"),
          S.tool === "cross" ? h("button", { onclick: () => setTool("pan") }, "Done") : null));
      b.append(curveList(["transect"]));
    }
    b.append(h("div", { class: "sub" }, "Count"));
    const ct = d.counts;
    b.append(question("Total growth marks counted (traced and untraced)", null, ct.visible_marks,
      h("div", { class: "row" }, h("input", { type: "number", min: 0, class: "num", value: ct.visible_marks ?? "",
        onchange: (e) => { ct.visible_marks = e.target.value === "" ? null : parseInt(e.target.value, 10); re(); } }),
        h("label", { class: "chk" }, h("input", { type: "checkbox", checked: ct.visible_at_least, onchange: (e) => { ct.visible_at_least = e.target.checked; re(); } }), "at least"),
        h("button", { class: ct.visible_marks === 0 ? "active" : "", onclick: () => { ct.visible_marks = 0; ct.visible_at_least = false; re(); } }, "None visible"))));
  },
  tissue(b) {
    const d = S.doc, ctx = { polarized: isPolarized(d) };
    const zs = zones().zones;
    if (!zs.length) {
      if (!d.field) d.field = { obs: {} };
      b.append(h("div", { class: "sub" }, "Whole field"));
      b.append(zoneQuestions(d.field, ctx, null));
      return;
    }
    const list = [...zs].reverse();
    if (!S.openZone || !list.some((z) => z.outer.id === S.openZone)) { const first = list.find((z) => !zoneDone(z.outer.zone, ctx)); S.openZone = (first || list[0]).outer.id; }
    list.forEach((z, i) => {
      const key = (z.inner ? z.inner.id : "") + "|" + z.outer.id, zd = z.outer.zone, ok = zoneDone(zd, ctx), open = S.openZone === z.outer.id;
      const name = `${z.inner ? shortName(z.inner) : "inner"} → ${shortName(z.outer)}`;
      const card = h("div", { class: "zonecard" + (open ? " open" : "") + (ok ? "" : " unset"),
          onmouseenter: () => { S.hover = "zone:" + key; renderAll(); }, onmouseleave: () => { S.hover = null; renderAll(); } },
        h("div", { class: "zhead", onclick: () => { S.openZone = open ? null : z.outer.id; renderPanels(); } },
          h("span", { class: "zname" }, (ok ? "✓ " : "") + name), h("span", { class: "detail" }, ok ? "described" : zd && zd.obs && Object.keys(zd.obs).length ? "in progress" : "not described")));
      if (open) { z.outer.zone = z.outer.zone || { obs: {} }; z.outer.zone.obs = z.outer.zone.obs || {};
        const innerZ = list[i + 1] && list[i + 1].outer.zone;
        card.append(h("div", { class: "zbody" }, zoneQuestions(z.outer.zone, ctx, { innerZ, all: list, next: list[i + 1] }))); }
      b.append(card);
    });
  },
  bone(b) {
    const d = S.doc, sc = d.scores, re = () => { touch(); renderPanels(); };
    const ask = (q) => {
      const cur = sc[q.key];
      const kids = [choices(q.opts.map(([v, l, ic]) => [v, l, ic]), cur, (v) => { sc[q.key] = v; if (v !== HG.UNCLEAR) delete sc[q.key + "_reason"]; re(); })];
      if (cur === HG.UNCLEAR && q.reasons) { const sel = h("select", { onchange: (e) => { sc[q.key + "_reason"] = e.target.value || null; re(); } }, h("option", { value: "" }, "why?"));
        for (const r of q.reasons) sel.append(h("option", { value: r, selected: sc[q.key + "_reason"] === r }, r)); kids.push(h("div", { class: "row small" }, "Because: ", sel)); }
      const tip = q.tip;
      if (q.key === "efs" && cur === "present") kids.push(h("div", { class: "row" }, kindButton("efs_inner", () => traceKind("efs_inner"), "EFS inner boundary")),
        h("div", { class: "row small" }, "Unresolved lines in EFS, at least ",
          h("input", { type: "number", min: 0, class: "num", value: d.counts.efs_lines_min ?? "", onchange: (e) => { d.counts.efs_lines_min = e.target.value === "" ? null : parseInt(e.target.value, 10); touch(); } })));
      return question(q.q, tip, cur, ...kids);
    };
    for (const q of HG.BONE_Q.filter((x) => !x.more)) b.append(ask(q));
    b.append(h("details", { class: "fold" }, h("summary", {}, "More questions"), ...HG.BONE_Q.filter((x) => x.more).map(ask)));
    b.append(h("div", { class: "sub" }, "Notes"), h("textarea", { rows: 3, placeholder: "", value: d.notes || "",
      oninput: (e) => { d.notes = e.target.value; touch(); } }));
  },
  more(b) {
    const n = S.doc.lacunae.length, areas = S.doc.curves.filter((c) => c.kind === "count_area").length;
    $("countInfo").textContent = n || areas ? `${n} counted${areas ? "" : " · no area"}` : "";
  },
  export(b) {
    const v = validate();
    if (!v.err.length && !v.warn.length) b.append(h("div", { class: "box ok" }, "Complete"));
    const ul = h("ul", { class: "check" });
    for (const [cls, items] of [["bad", v.err], ["warn", v.warn]]) for (const it of items)
      ul.append(h("li", { class: cls }, h("span", {}, it.msg), it.step ? h("button", { class: "link", onclick: () => goStep(it.step) }, "go") : null));
    b.append(ul);
    b.append(h("div", { class: "tools" }, h("button", { class: "primary", disabled: !!v.err.length, onclick: exportDoc }, "Export"),
      h("button", { disabled: S.idx + 1 >= S.tasks.length, onclick: nextTask }, "Next image ›")));
    if (stepState("export") === "done") b.append(h("p", { class: "small ok" }, "Exported"));
  },
};

function zoneQuestions(zone, ctx, nav) {
  zone.obs = zone.obs || {}; const o = zone.obs;
  const wrap = h("div", {});
  const re = () => { delete zone.legacy; Object.assign(zone, HG.interpretZone(o, ctx)); touch(); renderPanels(); };
  if (zone.legacy && !Object.keys(o).length) wrap.append(h("div", { class: "box" }, "0.3 description kept"));
  for (const q of HG.ZONE_Q) {
    if (q.when && !q.when(o, ctx)) { if (o[q.key] != null && q.key !== "matrix_xpl" && q.key !== "lacuna_shape") delete o[q.key]; continue; }
    wrap.append(question(q.q, q.tip, o[q.key], choices(q.opts, o[q.key], (v) => { if (v == null) delete o[q.key]; else o[q.key] = v; re(); })));
  }
  const t = HG.interpretZone(o, ctx);
  if (!(zone.legacy && !Object.keys(o).length)) Object.assign(zone, t);
  wrap.append(h("label", { class: "chk small" }, h("input", { type: "checkbox", checked: !!zone.varies, onchange: (e) => { zone.varies = e.target.checked; touch(); } }),
    "Varies around the section"));
  if (nav) {
    const copy = (src, dst) => { dst.obs = JSON.parse(JSON.stringify(src.obs || {})); Object.assign(dst, HG.interpretZone(dst.obs, ctx)); };
    wrap.append(h("div", { class: "tools" },
      h("button", { disabled: !(nav.innerZ && nav.innerZ.obs && Object.keys(nav.innerZ.obs).length), 
        onclick: () => { copy(nav.innerZ, zone); re(); } }, "Same as inner"),
      h("button", { disabled: !HG.complete(o, ctx),
        onclick: () => { let n = 0; for (const z of nav.all) { const c = z.outer; if (zoneDone(c.zone, ctx)) continue; c.zone = c.zone || {}; copy(zone, c.zone); n++; }
          status(`Copied to ${n} zone${n === 1 ? "" : "s"}. Check each one.`); renderPanels(); } }, "Copy to undescribed zones"),
      nav.next ? h("button", { class: "primary", onclick: () => { S.openZone = nav.next.outer.id; renderPanels(); } }, "Next zone ›") : null));
  }
  return wrap;
}

function curveList(kinds) {
  const ul = h("ul", { class: "curves" });
  const order = ["periosteum", "broken_edge", "efs_inner", "label", "mark", "remodel_front", "medullary", "transect", "count_area", "ignore"];
  const cs = S.doc.curves.filter((c) => kinds.includes(c.kind))
    .sort((a, b) => order.indexOf(a.kind) - order.indexOf(b.kind) || (b.mark_index || b.label_index || 0) - (a.mark_index || a.label_index || 0));
  if (!cs.length) { ul.append(h("li", { class: "empty" }, "Nothing traced yet.")); return ul; }
  for (const c of cs) {
    const isAct = (S.active && S.active.curveId === c.id) || S.activeT === c.id, empty = c.kind !== "transect" && !hasPts(c);
    const stop = (f) => (e) => { e.stopPropagation(); f(e); };
    const cont = () => { if (c.kind === "transect") { S.activeT = c.id; S.tool = "cross"; S.current = c.id; updateHint(); renderAll(); renderPanels(); } else startTrace(c.id); };
    ul.append(h("li", { class: (isAct ? "activeCurve " : "") + (S.current === c.id ? "currentCurve" : ""),
        onclick: () => { S.current = c.id; renderAll(); renderPanels(); },
        onmouseenter: () => { S.hover = c.id; renderAll(); }, onmouseleave: () => { S.hover = null; renderAll(); } },
      h("span", { class: "sw", style: `background:${KINDS[c.kind].color}` }),
      h("span", { class: "kind" + (c.branch_of ? " branch" : "") }, curveName(c)),
      h("span", { class: "detail" }, empty ? "no points" : curveDetail(c)),
      h("span", { class: "rowbtns" },
        h("button", { title: isAct ? "Stop" : "Continue", onclick: stop(() => isAct ? stopTrace() : cont()) }, isAct ? "■" : "▶"),
        h("button", { title: "Edit (T)", "data-pop": 1, onclick: stop((e) => { S.current = c.id; if (!isAct) S.active = null; openTag(e.currentTarget); }) }, "✎"),
        h("button", { title: "Delete", onclick: stop(() => { pushUndo(); S.doc.curves = S.doc.curves.filter((x) => x.id !== c.id);
          for (const x of S.doc.curves) if (x.branch_of === c.id) x.branch_of = null;
          if (isAct) setTool("pan"); if (S.current === c.id) S.current = null; closePop(); renderAll(); renderPanels(); }) }, "✕"))));
  }
  return ul;
}
function curveDetail(c) {
  if (c.kind === "mark") return [c.mark_type !== "LAG" ? c.mark_type : null, c.confidence < 3 ? ["", "guess", "fairly sure"][c.confidence] : null,
    c.in_efs ? "in EFS" : null, c.index_auto === false ? "numbered by hand" : null].filter(Boolean).join(" · ");
  if (c.kind === "label") return [c.dye, c.label_date || (c.age_days != null ? c.age_days + " d" : null)].filter(Boolean).join(" · ");
  if (c.kind === "transect") { const n = (c.crossings || []).filter((x) => x.kind === "mark").length; return `${n} mark${n === 1 ? "" : "s"} · ${nice(c.ends.outer)} → ${nice(c.ends.inner)}`; }
  if (c.kind === "count_area") { const n = S.doc.lacunae.filter((q) => c.segments[0] && c.segments[0].length > 2 && inPoly(q, c.segments[0])).length; return n + " lacunae"; }
  const n = c.segments.filter((s) => s.length > 1).length; return n > 1 ? `${n} pieces` : "";
}
function parseScale(s) { if (!s) return {}; const m = String(s).match(/([\d.]+)\s*(µm|um|μm|mm|cm|m)\b/i);
  return m ? { length: +m[1], unit: m[2].replace(/[µμ]m|um/i, "um").toLowerCase() } : {}; }

$("measurer").addEventListener("input", () => { store.set("hm:measurer", $("measurer").value.trim()); if (S.doc) { S.doc.measurer = $("measurer").value.trim() || null; if (S.step === "export") renderPanels(); } });

// ---------- validation ----------
function validate() {
  const d = S.doc, err = [], warn = [], E = (msg, step) => err.push({ msg, step }), W = (msg, step) => warn.push({ msg, step });
  if (!$("measurer").value.trim()) E("Measurer ID");
  if (!d.setup.view) E("Image type", "image");
  if (d.setup.view === "unusable") return { err, warn };
  if (!d.setup.light) W("Light", "image");
  if (!pxPerUnit() && !skipped("scale")) E(d.calibration ? "Scale-bar length" : "Scale bar", "scale");
  const c = d.setup.companion;
  if (c && c.transform_method === "scaled_to_width") W("Other-light image not aligned", "image");
  if (c && !c.light) W("Other-light image: light", "image");
  const marks = d.curves.filter((x) => x.kind === "mark" && !x.branch_of);
  if (applies("outline")) {
    if (!skipped("outline")) {
      if (!d.curves.some((x) => OUTER.has(x.kind) && hasPts(x))) W("Outer edge", "outline");
      if (!d.centre && !d.flags.no_centre) W("Centre", "outline");
      if (!d.section.position) W("Section position", "outline");
    }
    if (!skipped("marks") && stepState("marks") !== "done") W("Growth marks or count", "marks");
    if (marks.some((m) => m.mark_index == null && !m.index_uncertain)) E("Mark without number", "marks");
    const idx = marks.filter((m) => !m.index_uncertain && m.mark_index != null).map((m) => m.mark_index);
    const dup = [...new Set(idx.filter((v, i) => idx.indexOf(v) !== i))];
    if (dup.length) W(`Duplicate mark number: ${dup.join(", ")}`, "marks");
    const tr = marks.length + 0, cnt = d.counts.visible_marks;
    if (cnt != null && !d.counts.visible_at_least && cnt < tr) W(`Count ${cnt} < marks traced ${tr}`, "marks");
  }
  if (d.curves.some((x) => x.kind !== "transect" && !hasPts(x))) W("Line with no points");
  if (d.curves.some((x) => x.kind !== "transect" && x.segments.some((s) => s.length === 1))) W("Single-point segment");
  if (d.curves.some((x) => x.kind === "transect" && !x.crossings.length)) W("Transect without crossings", "marks");
  if (!skipped("tissue") && stepState("tissue") !== "done") W("Tissue", "tissue");
  if (d.scores.efs == null) W("EFS", "bone");
  if (d.scores.efs === "present" && !d.curves.some((x) => x.kind === "efs_inner")) W("EFS inner boundary", "bone");
  for (const q of HG.BONE_Q) if (d.scores[q.key] === HG.UNCLEAR && q.reasons && !d.scores[q.key + "_reason"]) W(`Reason missing: ${q.key.replace(/_/g, " ")}`, "bone");
  if (d.curves.some((x) => x.kind === "label" && !x.label_date && x.age_days == null)) W("Label date or age", "marks");
  if (d.lacunae.length && !d.curves.some((x) => x.kind === "count_area")) W("Lacuna count area", "more");
  return { err, warn };
}

// ---------- export: raw geometry plus derived measurements ----------
function derive(d) {
  const k = pxPerUnit(d), u = k ? d.calibration.unit : "px", cv = (px) => (k ? px / k : px);
  const z = d.centre, num = (v) => (v == null || !isFinite(v) ? null : +v.toFixed(4));
  const curves = d.curves.filter((c) => c.kind !== "transect").map((c) => {
    const closed = KINDS[c.kind].closed;
    let L = 0, nv = 0; for (const s of c.segments) { nv += s.length; for (let i = 1; i < s.length; i++) L += dist(s[i - 1], s[i]);
      if (closed && s.length > 2) L += dist(s[s.length - 1], s[0]); }
    const o = { id: c.id, kind: c.kind, name: curveName(c), n_segments: c.segments.length, n_points: nv, [closed ? "perimeter" : "length"]: num(cv(L)) };
    if (z && nv > 1 && !closed) { const pts = c.segments.flatMap((s) => resample(s, 2)), rs = pts.map((q) => dist(q, z));
      Object.assign(o, { radius_mean: num(cv(rs.reduce((a, b) => a + b, 0) / rs.length)), radius_min: num(cv(Math.min(...rs))),
        radius_max: num(cv(Math.max(...rs))) }, angular(c, z)); }
    if (closed && c.segments[0] && c.segments[0].length > 2) o.area = num(cv(cv(polyArea(c.segments[0]))));
    return o;
  });
  const transects = d.curves.filter((c) => c.kind === "transect").map((c) => { const a = c.segments[0][0];
    return { id: c.id, name: curveName(c), ends: c.ends, length: num(cv(dist(c.segments[0][0], c.segments[0][1]))),
      crossings: crossings(c).map((x, i) => ({ order_from_outer: i + 1, kind: x.kind, mark_index: x.kind === "mark" ? x._idx : null,
        dist_from_outer: num(cv(dist(x.p, a))), radius: z ? num(cv(dist(x.p, z))) : null })) }; });
  const zz = zones(d), ctx = { polarized: isPolarized(d) };
  const zoneRows = zz.zones.map((x, i) => Object.assign({ zone: i + 1, inner_id: x.inner ? x.inner.id : null, inner: x.inner ? curveName(x.inner) : null,
    outer_id: x.outer.id, outer: curveName(x.outer), complete: zoneDone(x.outer.zone, ctx) }, x.outer.zone || {}));
  const counts = d.curves.filter((c) => c.kind === "count_area" && c.segments[0] && c.segments[0].length > 2).map((c) => {
    const n = d.lacunae.filter((q) => inPoly(q, c.segments[0])).length, a = cv(cv(polyArea(c.segments[0])));
    return { area_id: c.id, n_lacunae: n, area: num(a), density_per_unit2: num(a ? n / a : null) };
  });
  const cover = (kind) => { if (!z) return null; const pts = d.curves.filter((c) => c.kind === kind).flatMap((c) => c.segments.flatMap((s) => resample(s, 2)));
    return pts.length > 1 ? coverageOf(pts, z).filter(Boolean).length * 5 : 0; };
  return { note: "Computed at export as a convenience. Analysis should recompute from the raw points.",
    unit: u, px_per_unit: num(k), reference_for_zone_order: zz.method, curves, transects, zones: zoneRows, lacuna_counts: counts,
    outer_edge: { periosteum_coverage_deg: cover("periosteum"), broken_edge_coverage_deg: cover("broken_edge"),
      surface_missing_deg: z && d.setup.view === "whole" ? 360 - cover("periosteum") : null },
    marks_traced: d.curves.filter((c) => c.kind === "mark" && !c.branch_of && hasPts(c)).length,
    step_states: Object.fromEntries(STEPS.map((s) => [s.id, stepState(s.id, d)])),
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
  if (!S.doc) return; autoIndex(); saveDoc();
  const v = validate(); if (v.err.length) { goStep("export"); status("Needed before export: " + v.err.map((x) => x.msg).join(", ")); return; }
  const at = new Date().toISOString(), measurer = $("measurer").value.trim();
  const stamp = at.replace(/[:.]/g, "").replace("T", "_").slice(0, 17);   // YYYY-MM-DD_HHMMSS (UTC)
  const base = `${S.doc.task.figure_id}__${measurer.replace(/[^\w-]+/g, "")}__${stamp}`;
  for (const f of hooks.beforeExport) f(S.doc, base);
  for (const c of S.doc.curves) if (c.kind === "transect") crossings(c).forEach((x) => { if (x.kind === "mark") x.mark_index = x._idx; });
  const d = JSON.parse(JSON.stringify(S.doc, (key, val) => (key === "_idx" ? undefined : val)));
  d.measurer = measurer; d.exported_at = at; d.validation_warnings = v.warn.map((x) => x.msg); d.rules_version = HG.RULES_VERSION;
  d.derived = deriveAll(S.doc);
  const where = await writeOut(base + ".json", JSON.stringify(d, null, 1), "application/json");
  for (const f of hooks.exportFiles) { const x = await f(S.doc, base); if (x) await writeOut(x.name, x.blob); }
  store.set(keyStatus(d.task.figure_id), "done"); renderTasks(); renderPanels(); status(`Saved ${where}`);
}
function nextTask() { for (let i = S.idx + 1; i < S.tasks.length; i++) if (!S.tasks[i]._missing && !isCompanionTask(S.tasks[i])) return openTask(i); status("This was the last image."); }
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
  await loadCompanion(); renderPanels(); renderAll(); saveDoc();
  for (const f of hooks.loaded) f(mask || null);
  status(`Loaded ${js.name}${mask ? " with its mask" : ""}. Export saves a new file.`);
}
$("btnLoad").addEventListener("click", () => { if (!S.doc) { status("Open the image folder first."); return; } $("loadInput").click(); });
$("loadInput").addEventListener("change", async (e) => { await loadExport([...e.target.files]); e.target.value = ""; });
$("btnExport").addEventListener("click", () => { if (S.doc) { goStep("export"); exportDoc(); } });
$("btnNext").addEventListener("click", nextTask);
$("btnPrev").addEventListener("click", () => { for (let i = S.idx - 1; i >= 0; i--) if (!S.tasks[i]._missing && !isCompanionTask(S.tasks[i])) return openTask(i); });
$("btnUndo").addEventListener("click", undo);
for (const id of ["btnOpen", "btnOpen2"]) $(id).addEventListener("click", openFolder);
$("btnDrawer").addEventListener("click", () => { $("drawer").hidden = !$("drawer").hidden; onResize(); });

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
    taxon: "Synthetica demonstrata", specimen: "DEMO-1", element: "femur", view: "full_section", light: "normal", quality_tier: "B",
    scale_value: "1 mm", caption: "Synthetic midshaft section: 9 growth marks, mark 6 branching over part of the cortex, drifting centre, one broken sector." }];
  S.idx = -1; renderTasks(); openTask(0);
}
for (const id of ["btnDemo", "btnDemo2"]) $(id).addEventListener("click", demo);

// ---------- API for features.js and scripted tests ----------
window.HM = { S, $, h, el, svg, layerVerts, hooks, status, touch, uid, seg, row, showPop, closePop, pxPerUnit, KINDS, curveName,
  resample, angle, angular, dist, inPoly, toImg, upp, setView, zoomAt, setTool, renderAll, renderPanels,
  zones, orderedBoundaries, fitCentre, autoIndex, validate, derive: deriveAll, migrate, goStep, openTask, loadFiles };

// ---------- start ----------
$("ver").textContent = TOOL_VERSION;
document.body.classList.add("nodoc");
$("measurer").value = store.get("hm:measurer") || "";
})();
